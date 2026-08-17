"""Filter a pairs file against a BED blacklist and a chromosome allow-list."""

from __future__ import annotations

import argparse
import bisect
import gzip
import io
from collections import defaultdict
from collections.abc import Iterable, Iterator
from contextlib import contextmanager
from itertools import chain
from pathlib import Path
from typing import TextIO

PAIRTOOLS_REQUIRED_COLUMNS = {
    "chrom1",
    "pos1",
    "chrom2",
    "pos2",
    "strand1",
    "strand2",
    "pair_type",
}


@contextmanager
def open_text(path: str | Path, mode: str) -> Iterator[TextIO]:
    """Open plain text, gzip input, or pairtools-compatible BGZF output."""

    path = Path(path)
    if path.suffix == ".gz":
        if mode == "r":
            # gzip can read both ordinary gzip and BGZF streams.
            with gzip.open(path, "rt", encoding="utf-8", newline="") as handle:
                yield handle
        else:
            try:
                import pysam
            except ImportError as error:
                raise RuntimeError(
                    "pysam is required to write pairtools-compatible .gz output"
                ) from error
            with (
                pysam.BGZFile(str(path), "w") as compressed,
                io.TextIOWrapper(compressed, encoding="utf-8", newline="") as handle,
            ):
                yield handle
    else:
        with path.open(mode, encoding="utf-8", newline="") as handle:
            yield handle


def load_chromosomes(path: str | Path) -> dict[str, int]:
    """Load chromosome names and lengths, preserving input order."""

    chromosomes: dict[str, int] = {}
    with open_text(path, "r") as handle:
        for line_number, line in enumerate(handle, start=1):
            if not line.strip() or line.startswith("#"):
                continue

            fields = line.split()
            if len(fields) < 2:
                raise ValueError(
                    f"{path}:{line_number}: expected chromosome and length"
                )
            chrom = fields[0]
            try:
                length = int(fields[1])
            except ValueError as error:
                raise ValueError(
                    f"{path}:{line_number}: invalid chromosome length {fields[1]!r}"
                ) from error
            if length <= 0:
                raise ValueError(
                    f"{path}:{line_number}: chromosome length must be positive"
                )
            if chrom in chromosomes:
                raise ValueError(
                    f"{path}:{line_number}: duplicate chromosome {chrom!r}"
                )
            chromosomes[chrom] = length

    if not chromosomes:
        raise ValueError(f"{path}: no chromosomes found")
    return chromosomes


def load_blacklist(path: str | Path) -> dict[str, tuple[list[int], list[int]]]:
    """Load BED intervals and merge overlapping or directly adjacent entries."""

    raw_intervals: defaultdict[str, list[tuple[int, int]]] = defaultdict(list)

    with open_text(path, "r") as handle:
        for line_number, line in enumerate(handle, start=1):
            if not line.strip() or line.startswith("#"):
                continue

            fields = line.split()
            if len(fields) < 3:
                raise ValueError(f"{path}:{line_number}: expected a BED3 record")
            chrom = fields[0]
            try:
                start = int(fields[1])
                end = int(fields[2])
            except ValueError as error:
                raise ValueError(
                    f"{path}:{line_number}: BED coordinates must be integers"
                ) from error
            if start < 0:
                raise ValueError(
                    f"{path}:{line_number}: BED start must be non-negative"
                )
            if start >= end:
                continue
            raw_intervals[chrom].append((start, end))

    merged: dict[str, tuple[list[int], list[int]]] = {}
    for chrom, intervals in raw_intervals.items():
        intervals.sort()
        chrom_merged: list[list[int]] = []

        for start, end in intervals:
            if not chrom_merged or start > chrom_merged[-1][1]:
                chrom_merged.append([start, end])
            else:
                chrom_merged[-1][1] = max(chrom_merged[-1][1], end)

        starts = [interval[0] for interval in chrom_merged]
        ends = [interval[1] for interval in chrom_merged]
        merged[chrom] = (starts, ends)

    return merged


def overlaps_blacklist(
    blacklist: dict[str, tuple[list[int], list[int]]], chrom: str, pos_1based: int
) -> bool:
    """Return whether a 1-based pairs position overlaps a BED interval."""

    if chrom not in blacklist:
        return False

    # pairs coordinates are 1-based; BED intervals are 0-based and half-open.
    pos_0based = pos_1based - 1
    starts, ends = blacklist[chrom]
    index = bisect.bisect_right(starts, pos_0based) - 1
    return index >= 0 and pos_0based < ends[index]


def read_pairs_header(src: TextIO) -> tuple[list[str], str | None]:
    """Consume the leading pairs header and return it plus the first body line."""

    header: list[str] = []
    for line in src:
        if line.startswith("#"):
            header.append(line.rstrip("\r\n"))
        else:
            return header, line
    return header, None


def parse_columns(header: Iterable[str]) -> list[str]:
    """Read and validate the pairtools column declaration."""

    declarations = [
        line.split(":", 1)[1].strip().split()
        for line in header
        if line.startswith("#columns:")
    ]
    if len(declarations) != 1:
        raise ValueError("input must contain exactly one '#columns:' header line")

    columns = declarations[0]
    if len(columns) != len(set(columns)):
        raise ValueError("the '#columns:' header contains duplicate column names")
    missing = sorted(PAIRTOOLS_REQUIRED_COLUMNS.difference(columns))
    if missing:
        raise ValueError(
            "input is not compatible with pairtools stats/scaling; "
            f"missing columns: {', '.join(missing)}"
        )
    return columns


def output_header(header: Iterable[str], chrom_sizes: dict[str, int]) -> list[str]:
    """Replace chromosome metadata with the chromosome allow-list."""

    retained = [
        line
        for line in header
        if not line.startswith("#chromsize:") and not line.startswith("#chromosomes:")
    ]
    chrom_lines = [
        f"#chromsize: {chrom} {length}" for chrom, length in chrom_sizes.items()
    ]

    # pairtools expects chromosome metadata before the column declaration.
    columns_index = next(
        (index for index, line in enumerate(retained) if line.startswith("#columns:")),
        len(retained),
    )
    retained[columns_index:columns_index] = chrom_lines
    return retained


def filter_pairs(
    input_path: str | Path,
    output_path: str | Path,
    blacklist_path: str | Path,
    chrom_sizes_path: str | Path,
    stats_path: str | Path,
) -> dict[str, int]:
    """Filter pairs and return filtering counts."""

    input_path = Path(input_path)
    output_path = Path(output_path)
    blacklist_path = Path(blacklist_path)
    chrom_sizes_path = Path(chrom_sizes_path)
    stats_path = Path(stats_path)
    if output_path.resolve() == stats_path.resolve():
        raise ValueError("output and stats paths must be different")
    read_paths = {
        "input": input_path.resolve(),
        "blacklist": blacklist_path.resolve(),
        "chrom-sizes": chrom_sizes_path.resolve(),
    }
    for output_name, path in {
        "output": output_path.resolve(),
        "stats": stats_path.resolve(),
    }.items():
        for input_name, input_path_resolved in read_paths.items():
            if path == input_path_resolved:
                raise ValueError(
                    f"{output_name} and {input_name} paths must be different"
                )

    chrom_sizes = load_chromosomes(chrom_sizes_path)
    blacklist = load_blacklist(blacklist_path)
    counts = {
        "input_pairs": 0,
        "kept_pairs": 0,
        "blacklisted_pairs": 0,
        "excluded_chromosome_pairs": 0,
        "malformed_pairs": 0,
    }

    with open_text(input_path, "r") as src:
        header, first_body_line = read_pairs_header(src)
        columns = parse_columns(header)
        column_indexes = {column: index for index, column in enumerate(columns)}
        prepared_header = output_header(header, chrom_sizes)

        with open_text(output_path, "w") as dst:
            for line in prepared_header:
                dst.write(f"{line}\n")

            body: Iterable[str]
            if first_body_line is None:
                body = src
            else:
                body = chain((first_body_line,), src)

            for line in body:
                counts["input_pairs"] += 1
                fields = line.rstrip("\r\n").split("\t")
                if len(fields) != len(columns):
                    counts["malformed_pairs"] += 1
                    continue

                chrom1 = fields[column_indexes["chrom1"]]
                chrom2 = fields[column_indexes["chrom2"]]
                if chrom1 not in chrom_sizes or chrom2 not in chrom_sizes:
                    counts["excluded_chromosome_pairs"] += 1
                    continue

                try:
                    pos1 = int(fields[column_indexes["pos1"]])
                    pos2 = int(fields[column_indexes["pos2"]])
                except ValueError:
                    counts["malformed_pairs"] += 1
                    continue
                if not (1 <= pos1 <= chrom_sizes[chrom1]) or not (
                    1 <= pos2 <= chrom_sizes[chrom2]
                ):
                    counts["malformed_pairs"] += 1
                    continue

                if overlaps_blacklist(blacklist, chrom1, pos1) or overlaps_blacklist(
                    blacklist, chrom2, pos2
                ):
                    counts["blacklisted_pairs"] += 1
                    continue

                dst.write(line)
                counts["kept_pairs"] += 1

    with open_text(stats_path, "w") as handle:
        handle.write("category\tcount\n")
        for category, count in counts.items():
            handle.write(f"{category}\t{count}\n")

    return counts


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--blacklist", required=True)
    parser.add_argument("--chrom-sizes", required=True)
    parser.add_argument("--stats", required=True)
    args = parser.parse_args()

    try:
        filter_pairs(
            args.input,
            args.output,
            args.blacklist,
            args.chrom_sizes,
            args.stats,
        )
    except (OSError, RuntimeError, ValueError) as error:
        parser.error(str(error))


if __name__ == "__main__":
    main()
