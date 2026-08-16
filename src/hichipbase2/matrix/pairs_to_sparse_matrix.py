from __future__ import annotations

import argparse
import gzip
from collections import Counter
from pathlib import Path


def read_chromsizes(
    path: Path,
    binsize: int,
) -> tuple[
    dict[str, int],
    dict[str, int],
    list[tuple[str, int, int]],
]:
    chromosome_lengths: dict[str, int] = {}
    chromosome_offsets: dict[str, int] = {}
    bins: list[tuple[str, int, int]] = []

    global_bin_offset = 0

    with path.open("rt", encoding="utf-8") as handle:
        for line in handle:
            chromosome, length_text = line.rstrip().split("\t")[:2]
            length = int(length_text)

            chromosome_lengths[chromosome] = length
            chromosome_offsets[chromosome] = global_bin_offset

            for start in range(0, length, binsize):
                end = min(start + binsize, length)
                bins.append((chromosome, start, end))
                global_bin_offset += 1

    return chromosome_lengths, chromosome_offsets, bins


def parse_columns(line: str) -> list[str]:
    return (
        line.removeprefix("#columns:")
        .strip()
        .split()
    )


def make_matrix(
    pairs_path: Path,
    chromsizes_path: Path,
    binsize: int,
    bins_output: Path,
    pixels_output: Path,
) -> None:
    chromosome_lengths, chromosome_offsets, bins = read_chromsizes(
        chromsizes_path,
        binsize,
    )

    counts: Counter[tuple[int, int]] = Counter()
    columns: list[str] | None = None
    total_pairs = 0

    with gzip.open(
        pairs_path,
        "rt",
        encoding="utf-8",
    ) as handle:
        for line in handle:
            if line.startswith("#columns:"):
                columns = parse_columns(line)
                continue

            if line.startswith("#"):
                continue

            values = line.rstrip("\n").split("\t")

            if columns is None:
                # Standard .pairs first columns.
                columns = [
                    "readID",
                    "chrom1",
                    "pos1",
                    "chrom2",
                    "pos2",
                    "strand1",
                    "strand2",
                    "pair_type",
                ]

            record = dict(zip(columns, values))

            chrom1 = record["chrom1"]
            chrom2 = record["chrom2"]
            pos1 = int(record["pos1"])
            pos2 = int(record["pos2"])

            if chrom1 not in chromosome_lengths:
                raise ValueError(f"Unknown chromosome: {chrom1}")

            if chrom2 not in chromosome_lengths:
                raise ValueError(f"Unknown chromosome: {chrom2}")

            if not 1 <= pos1 <= chromosome_lengths[chrom1]:
                raise ValueError(
                    f"Position outside chromosome: {chrom1}:{pos1}"
                )

            if not 1 <= pos2 <= chromosome_lengths[chrom2]:
                raise ValueError(
                    f"Position outside chromosome: {chrom2}:{pos2}"
                )

            bin1 = (
                chromosome_offsets[chrom1]
                + (pos1 - 1) // binsize
            )
            bin2 = (
                chromosome_offsets[chrom2]
                + (pos2 - 1) // binsize
            )

            if bin1 > bin2:
                bin1, bin2 = bin2, bin1

            counts[(bin1, bin2)] += 1
            total_pairs += 1

    with bins_output.open("wt", encoding="utf-8") as output:
        for chromosome, start, end in bins:
            output.write(
                f"{chromosome}\t{start}\t{end}\n"
            )

    with gzip.open(
        pixels_output,
        "wt",
        encoding="utf-8",
    ) as output:
        for (bin1, bin2), count in sorted(counts.items()):
            output.write(
                f"{bin1}\t{bin2}\t{count}\n"
            )

    print(f"Bin size: {binsize:,} bp")
    print(f"Bins: {len(bins):,}")
    print(f"Input pairs: {total_pairs:,}")
    print(f"Nonzero pixels: {len(counts):,}")
    print(f"Matrix sum: {sum(counts.values()):,}")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--pairs", required=True, type=Path)
    parser.add_argument("--chromsizes", required=True, type=Path)
    parser.add_argument("--binsize", required=True, type=int)
    parser.add_argument("--bins-output", required=True, type=Path)
    parser.add_argument("--pixels-output", required=True, type=Path)

    arguments = parser.parse_args()

    make_matrix(
        pairs_path=arguments.pairs,
        chromsizes_path=arguments.chromsizes,
        binsize=arguments.binsize,
        bins_output=arguments.bins_output,
        pixels_output=arguments.pixels_output,
    )


if __name__ == "__main__":
    main()