#!/usr/bin/env python3

from __future__ import annotations

import argparse
import csv
import gzip
import math
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import TextIO

from pycirclize import Circos


@dataclass(frozen=True)
class Gene:
    chrom: str
    start: int
    end: int
    name: str


@dataclass
class Loop:
    loop_id: str
    chrom_a: str
    start_a: int
    end_a: int
    chrom_b: str
    start_b: int
    end_b: int
    score: float
    gene_a: Gene | None = None
    gene_b: Gene | None = None


def open_text(path: str | Path) -> TextIO:
    path = Path(path)
    if path.suffix == ".gz":
        return gzip.open(path, "rt")
    return path.open("r")


def read_chrom_sizes(path: str | Path) -> dict[str, int]:
    chrom_sizes: dict[str, int] = {}

    with open_text(path) as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue

            chrom, size, *_ = line.rstrip().split("\t")
            chrom_sizes[chrom] = int(size)

    if not chrom_sizes:
        raise ValueError("chrom.sizesから染色体を読み込めませんでした")

    return chrom_sizes


def parse_attributes(text: str) -> dict[str, str]:
    """GTF形式とGFF3形式の属性を簡易的に解釈する。"""

    attributes: dict[str, str] = {}

    for item in text.strip().rstrip(";").split(";"):
        item = item.strip()
        if not item:
            continue

        if "=" in item:
            # GFF3: ID=gene1;Name=abc1
            key, value = item.split("=", 1)
        else:
            # GTF: gene_id "gene1"; gene_name "abc1";
            parts = item.split(None, 1)
            if len(parts) != 2:
                continue
            key, value = parts

        attributes[key.strip()] = value.strip().strip('"')

    return attributes


def read_genes(
    annotation_path: str | Path,
    allowed_chroms: set[str],
) -> dict[str, list[Gene]]:
    genes: dict[str, list[Gene]] = defaultdict(list)

    with open_text(annotation_path) as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue

            fields = line.rstrip().split("\t")
            if len(fields) < 9:
                continue

            chrom, _, feature_type, start, end, _, _, _, attr_text = fields

            if chrom not in allowed_chroms:
                continue
            if feature_type not in {"gene", "pseudogene"}:
                continue

            attrs = parse_attributes(attr_text)

            gene_name = (
                attrs.get("gene_name")
                or attrs.get("Name")
                or attrs.get("gene")
                or attrs.get("gene_id")
                or attrs.get("ID")
            )

            if not gene_name:
                continue

            # GTF/GFF3は1-based closedなので、0-based half-openへ変換
            genes[chrom].append(
                Gene(
                    chrom=chrom,
                    start=int(start) - 1,
                    end=int(end),
                    name=gene_name,
                )
            )

    for chrom in genes:
        genes[chrom].sort(key=lambda gene: gene.start)

    return genes


def read_cloops2_loops(
    path: str | Path,
    rank_field: str,
    chrom_sizes: dict[str, int],
) -> list[Loop]:
    loops: list[Loop] = []

    with open_text(path) as handle:
        header_line = handle.readline().rstrip()
        header = header_line.split("\t")
        header[0] = header[0].lstrip("#")

        column = {name: index for index, name in enumerate(header)}

        required = {
            "loopId",
            "chrA",
            "startA",
            "endA",
            "chrB",
            "startB",
            "endB",
            rank_field,
        }

        missing = required - column.keys()
        if missing:
            raise ValueError(
                "cLoops2 loopファイルに必要な列がありません: "
                + ", ".join(sorted(missing))
            )

        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue

            fields = line.rstrip().split("\t")

            if "significant" in column:
                significant = fields[column["significant"]].lower()
                if significant not in {"1", "true"}:
                    continue

            chrom_a = fields[column["chrA"]]
            chrom_b = fields[column["chrB"]]

            if chrom_a not in chrom_sizes or chrom_b not in chrom_sizes:
                continue

            start_a = int(fields[column["startA"]])
            end_a = int(fields[column["endA"]])
            start_b = int(fields[column["startB"]])
            end_b = int(fields[column["endB"]])

            # 描画範囲からはみ出さないように補正
            start_a = max(0, min(start_a, chrom_sizes[chrom_a] - 1))
            end_a = max(start_a + 1, min(end_a, chrom_sizes[chrom_a]))
            start_b = max(0, min(start_b, chrom_sizes[chrom_b] - 1))
            end_b = max(start_b + 1, min(end_b, chrom_sizes[chrom_b]))

            loops.append(
                Loop(
                    loop_id=fields[column["loopId"]],
                    chrom_a=chrom_a,
                    start_a=start_a,
                    end_a=end_a,
                    chrom_b=chrom_b,
                    start_b=start_b,
                    end_b=end_b,
                    score=float(fields[column[rank_field]]),
                )
            )

    return loops


def interval_distance(
    start1: int,
    end1: int,
    start2: int,
    end2: int,
) -> int:
    if end1 <= start2:
        return start2 - end1
    if end2 <= start1:
        return start1 - end2
    return 0


def nearest_gene(
    chrom: str,
    start: int,
    end: int,
    genes: dict[str, list[Gene]],
    max_distance: int,
) -> Gene | None:
    candidates = genes.get(chrom, [])
    if not candidates:
        return None

    gene = min(
        candidates,
        key=lambda item: interval_distance(
            start,
            end,
            item.start,
            item.end,
        ),
    )

    distance = interval_distance(start, end, gene.start, gene.end)
    if distance > max_distance:
        return None

    return gene


def select_top_loops(
    loops: list[Loop],
    top_percent: float,
    max_links: int,
) -> list[Loop]:
    if not loops:
        return []

    loops = sorted(loops, key=lambda item: item.score, reverse=True)

    number = max(1, math.ceil(len(loops) * top_percent / 100.0))

    if max_links > 0:
        number = min(number, max_links)

    return loops[:number]


def write_selected_loops(path: Path, loops: list[Loop]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(
            [
                "loopId",
                "chrA",
                "startA",
                "endA",
                "geneA",
                "chrB",
                "startB",
                "endB",
                "geneB",
                "score",
            ]
        )

        for loop in loops:
            writer.writerow(
                [
                    loop.loop_id,
                    loop.chrom_a,
                    loop.start_a,
                    loop.end_a,
                    loop.gene_a.name if loop.gene_a else ".",
                    loop.chrom_b,
                    loop.start_b,
                    loop.end_b,
                    loop.gene_b.name if loop.gene_b else ".",
                    loop.score,
                ]
            )


def draw_circos(
    chrom_sizes: dict[str, int],
    loops: list[Loop],
    output: Path,
    top_percent: float,
    rank_field: str,
    tick_interval: int,
) -> None:
    circos = Circos(chrom_sizes, space=4)

    chromosome_colors = [
        "#4C78A8",
        "#F58518",
        "#54A24B",
        "#E45756",
        "#72B7B2",
    ]

    chromosome_tracks = {}

    for index, sector in enumerate(circos.sectors):
        color = chromosome_colors[index % len(chromosome_colors)]

        sector.text(
            sector.name,
            r=112,
            size=13,
            weight="bold",
        )

        track = sector.add_track((92, 100))
        chromosome_tracks[sector.name] = track

        track.axis(
            fc=color,
            ec="black",
            lw=0.6,
            alpha=0.85,
        )

        track.xticks_by_interval(
            tick_interval,
            label_size=6,
            label_orientation="vertical",
            show_endlabel=False,
            label_formatter=lambda value: f"{value / 1_000_000:.1f} Mb",
            line_kws={"color": "black", "lw": 0.4},
        )

    scores = [loop.score for loop in loops]
    score_min = min(scores)
    score_max = max(scores)

    def normalized_score(score: float) -> float:
        if score_max == score_min:
            return 1.0
        return (score - score_min) / (score_max - score_min)

    annotated_genes: dict[tuple[str, int, int, str], Gene] = {}

    for loop in loops:
        normalized = normalized_score(loop.score)

        # 高スコアほど濃く表示
        alpha = 0.20 + 0.65 * normalized
        color = "#3B4CC0" if loop.chrom_a == loop.chrom_b else "#B40426"

        circos.link(
            (
                loop.chrom_a,
                loop.start_a,
                loop.end_a,
            ),
            (
                loop.chrom_b,
                loop.start_b,
                loop.end_b,
            ),
            r1=91,
            r2=91,
            color=color,
            alpha=alpha,
            ec="none",
            allow_twist=False,
        )

        if loop.gene_a:
            gene = loop.gene_a
            annotated_genes[
                (gene.chrom, gene.start, gene.end, gene.name)
            ] = gene

        if loop.gene_b:
            gene = loop.gene_b
            annotated_genes[
                (gene.chrom, gene.start, gene.end, gene.name)
            ] = gene

    # 上位loopに関連するgeneだけを外周に明記
    for gene in annotated_genes.values():
        position = (gene.start + gene.end) / 2
        track = chromosome_tracks[gene.chrom]

        track.annotate(
            position,
            gene.name,
            min_r=100,
            max_r=108,
            label_size=6,
            shorten=None,
            line_kws={
                "color": "#555555",
                "lw": 0.4,
            },
            text_kws={
                "color": "#202020",
            },
        )

    figure = circos.plotfig(
        dpi=300,
        figsize=(11, 11),
    )

    figure.suptitle(
        f"Top {top_percent:g}% loops ranked by {rank_field}\n"
        f"n = {len(loops)}",
        fontsize=14,
        y=0.98,
    )

    figure.savefig(
        output,
        dpi=300,
        bbox_inches="tight",
    )


def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "cLoops2 loop出力から上位loopを抽出し、"
            "最寄りgene名付きCircos図を作成する"
        )
    )

    parser.add_argument("--chrom-sizes", required=True)
    parser.add_argument("--loops", required=True)
    parser.add_argument("--annotation", required=True)
    parser.add_argument("--output", required=True)

    parser.add_argument(
        "--top-percent",
        type=float,
        default=2.0,
    )
    parser.add_argument(
        "--max-links",
        type=int,
        default=50,
    )
    parser.add_argument(
        "--rank-field",
        default="PETs",
        choices=["PETs", "enrichmentScore", "P2LL"],
    )
    parser.add_argument(
        "--max-gene-distance",
        type=int,
        default=2000,
    )
    parser.add_argument(
        "--tick-interval",
        type=int,
        default=1_000_000,
    )

    args = parser.parse_args()

    chrom_sizes = read_chrom_sizes(args.chrom_sizes)

    all_loops = read_cloops2_loops(
        args.loops,
        args.rank_field,
        chrom_sizes,
    )

    selected_loops = select_top_loops(
        all_loops,
        args.top_percent,
        args.max_links,
    )

    if not selected_loops:
        raise SystemExit("描画対象となるloopがありません")

    genes = read_genes(
        args.annotation,
        set(chrom_sizes),
    )

    for loop in selected_loops:
        loop.gene_a = nearest_gene(
            loop.chrom_a,
            loop.start_a,
            loop.end_a,
            genes,
            args.max_gene_distance,
        )
        loop.gene_b = nearest_gene(
            loop.chrom_b,
            loop.start_b,
            loop.end_b,
            genes,
            args.max_gene_distance,
        )

    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)

    selected_path = output.with_suffix(".selected.tsv")
    write_selected_loops(selected_path, selected_loops)

    draw_circos(
        chrom_sizes,
        selected_loops,
        output,
        args.top_percent,
        args.rank_field,
        args.tick_interval,
    )

    print(f"Total significant loops: {len(all_loops)}")
    print(f"Selected loops: {len(selected_loops)}")
    print(f"Figure: {output}")
    print(f"Selected table: {selected_path}")


if __name__ == "__main__":
    main()