"""Build an annotation-derived S. pombe blacklist BED."""

from __future__ import annotations

import argparse
from collections import Counter, defaultdict
from collections.abc import Iterable, Sequence
from dataclasses import dataclass
from pathlib import Path
from urllib.parse import unquote

from Bio import SeqIO

PRIMARY_CHROMOSOMES = ("I", "II", "III")
DEFAULT_SUBTELOMERE_SIZE = 50_000
RDNA_CHROMOSOME = "III"
REPEAT_KEYWORDS = ("duplicat", "centromer", "telomer", "subtelomer")

FEATURE_REASONS = {
    "regional_centromere": "centromere",
    "gap": "assembly_gap",
    "mating_type_region": "mating_type_region",
    "LTR_retrotransposon": "ltr_retrotransposon",
    "long_terminal_repeat": "long_terminal_repeat",
    "low_complexity_region": "low_complexity",
    "rRNA": "rrna",
    "nuclear_mt_pseudogene": "numt",
}


@dataclass(frozen=True)
class BlacklistProfile:
    """Rules used to select and synthesize blacklist intervals."""

    feature_types: frozenset[str]
    annotated_repeat_regions: bool
    all_repeat_regions: bool
    subtelomere_chromosomes: tuple[str, ...]
    terminal_rdna: bool = True
    non_primary_contigs: bool = True


CORE_FEATURE_TYPES = frozenset({"regional_centromere", "gap"})
BALANCED_FEATURE_TYPES = CORE_FEATURE_TYPES | frozenset(
    {
        "mating_type_region",
        "LTR_retrotransposon",
        "long_terminal_repeat",
        "low_complexity_region",
        "rRNA",
        "nuclear_mt_pseudogene",
    }
)

PROFILES = {
    "core": BlacklistProfile(
        feature_types=CORE_FEATURE_TYPES,
        annotated_repeat_regions=False,
        all_repeat_regions=False,
        subtelomere_chromosomes=("I", "II"),
    ),
    "balanced": BlacklistProfile(
        feature_types=BALANCED_FEATURE_TYPES,
        annotated_repeat_regions=True,
        all_repeat_regions=False,
        subtelomere_chromosomes=("I", "II"),
    ),
    "comprehensive": BlacklistProfile(
        feature_types=BALANCED_FEATURE_TYPES,
        annotated_repeat_regions=False,
        all_repeat_regions=True,
        subtelomere_chromosomes=PRIMARY_CHROMOSOMES,
    ),
}


@dataclass(frozen=True)
class GffFeature:
    """A validated GFF3 feature using BED-style coordinates."""

    seqid: str
    feature_type: str
    start: int
    end: int
    attributes: dict[str, str]


@dataclass(frozen=True)
class BlacklistInterval:
    """A BED interval and its biological or technical reasons."""

    chrom: str
    start: int
    end: int
    reasons: frozenset[str]


def load_sequence_lengths(fasta_path: Path) -> dict[str, int]:
    """Read FASTA identifiers and lengths while preserving sequence order."""

    lengths: dict[str, int] = {}

    with fasta_path.open("rt", encoding="utf-8") as fasta_input:
        for record in SeqIO.parse(fasta_input, "fasta"):
            if record.id in lengths:
                raise ValueError(
                    f"{fasta_path}: duplicate FASTA sequence {record.id!r}"
                )
            lengths[record.id] = len(record.seq)

    if not lengths:
        raise ValueError(f"{fasta_path}: no FASTA sequences found")

    missing = [chrom for chrom in PRIMARY_CHROMOSOMES if chrom not in lengths]
    if missing:
        raise ValueError(
            f"{fasta_path}: missing primary chromosome(s): {', '.join(missing)}"
        )

    return lengths


def parse_attributes(text: str) -> dict[str, str]:
    """Parse the GFF3 attributes column."""

    if text == ".":
        return {}

    attributes: dict[str, str] = {}

    for item in text.split(";"):
        if not item:
            continue
        if "=" not in item:
            raise ValueError(f"invalid GFF3 attribute {item!r}")
        key, value = item.split("=", maxsplit=1)
        if not key:
            raise ValueError(f"invalid GFF3 attribute {item!r}")
        attributes[key] = unquote(value)

    return attributes


def load_gff_features(
    gff_path: Path,
    sequence_lengths: dict[str, int],
) -> list[GffFeature]:
    """Load GFF3 features and validate them against the FASTA."""

    features: list[GffFeature] = []

    with gff_path.open("rt", encoding="utf-8") as gff_input:
        for line_number, line in enumerate(gff_input, start=1):
            if not line.strip() or line.startswith("#"):
                continue

            fields = line.rstrip("\r\n").split("\t")
            if len(fields) != 9:
                raise ValueError(f"{gff_path}:{line_number}: expected 9 GFF3 columns")

            seqid, _, feature_type, start_text, end_text, _, _, _, attr_text = fields

            if seqid not in sequence_lengths:
                raise ValueError(
                    f"{gff_path}:{line_number}: sequence {seqid!r} is absent "
                    "from the FASTA"
                )

            try:
                start_1based = int(start_text)
                end_1based = int(end_text)
            except ValueError as error:
                raise ValueError(
                    f"{gff_path}:{line_number}: invalid GFF3 coordinates"
                ) from error

            if not (1 <= start_1based <= end_1based <= sequence_lengths[seqid]):
                raise ValueError(
                    f"{gff_path}:{line_number}: coordinates "
                    f"{start_1based}-{end_1based} are outside {seqid} "
                    f"(length {sequence_lengths[seqid]})"
                )

            try:
                attributes = parse_attributes(attr_text)
            except ValueError as error:
                raise ValueError(f"{gff_path}:{line_number}: {error}") from error

            features.append(
                GffFeature(
                    seqid=seqid,
                    feature_type=feature_type,
                    start=start_1based - 1,
                    end=end_1based,
                    attributes=attributes,
                )
            )

    if not features:
        raise ValueError(f"{gff_path}: no GFF3 features found")

    return features


def repeat_region_is_annotated(feature: GffFeature) -> bool:
    """Return whether a repeat annotation describes a mapping-risk region."""

    annotation = " ".join(feature.attributes.values()).lower()
    return any(keyword in annotation for keyword in REPEAT_KEYWORDS)


def feature_reason(
    feature: GffFeature,
    profile: BlacklistProfile,
) -> str | None:
    """Return the blacklist reason for a selected GFF3 feature."""

    if feature.feature_type in profile.feature_types:
        return FEATURE_REASONS[feature.feature_type]

    if feature.feature_type != "repeat_region":
        return None

    if profile.all_repeat_regions:
        return "repeat_region"
    if profile.annotated_repeat_regions and repeat_region_is_annotated(feature):
        return "annotated_repeat"

    return None


def interval(
    chrom: str,
    start: int,
    end: int,
    reason: str,
) -> BlacklistInterval:
    """Construct one blacklist interval."""

    return BlacklistInterval(chrom, start, end, frozenset({reason}))


def terminal_rdna_intervals(
    features: Iterable[GffFeature],
    chromosome_length: int,
    terminal_window: int,
) -> list[BlacklistInterval]:
    """Infer chromosome III terminal rDNA arrays from terminal gene groups."""

    gene_groups = [
        feature
        for feature in features
        if feature.seqid == RDNA_CHROMOSOME and feature.feature_type == "gene_group"
    ]
    left_groups = [feature for feature in gene_groups if feature.end <= terminal_window]
    right_groups = [
        feature
        for feature in gene_groups
        if feature.start >= chromosome_length - terminal_window
    ]

    if not left_groups or not right_groups:
        raise ValueError("could not identify terminal chromosome III rDNA gene groups")

    return [
        interval(
            RDNA_CHROMOSOME,
            0,
            max(feature.end for feature in left_groups),
            "rdna_array",
        ),
        interval(
            RDNA_CHROMOSOME,
            min(feature.start for feature in right_groups),
            chromosome_length,
            "rdna_array",
        ),
    ]


def collect_intervals(
    sequence_lengths: dict[str, int],
    features: Sequence[GffFeature],
    profile: BlacklistProfile,
    subtelomere_size: int = DEFAULT_SUBTELOMERE_SIZE,
) -> list[BlacklistInterval]:
    """Collect GFF-derived and synthetic intervals for a profile."""

    if subtelomere_size <= 0:
        raise ValueError("subtelomere size must be positive")

    primary = frozenset(PRIMARY_CHROMOSOMES)
    intervals: list[BlacklistInterval] = []

    for feature in features:
        if feature.seqid not in primary:
            continue

        reason = feature_reason(feature, profile)
        if reason is not None:
            intervals.append(
                interval(feature.seqid, feature.start, feature.end, reason)
            )

    for chrom in profile.subtelomere_chromosomes:
        chrom_length = sequence_lengths[chrom]
        mask_size = min(subtelomere_size, chrom_length)
        intervals.extend(
            [
                interval(chrom, 0, mask_size, "subtelomere"),
                interval(
                    chrom,
                    chrom_length - mask_size,
                    chrom_length,
                    "subtelomere",
                ),
            ]
        )

    if profile.terminal_rdna:
        intervals.extend(
            terminal_rdna_intervals(
                features,
                sequence_lengths[RDNA_CHROMOSOME],
                terminal_window=min(
                    subtelomere_size,
                    sequence_lengths[RDNA_CHROMOSOME],
                ),
            )
        )

    if profile.non_primary_contigs:
        for chrom, chrom_length in sequence_lengths.items():
            if chrom not in primary:
                intervals.append(interval(chrom, 0, chrom_length, "non_primary_contig"))

    return intervals


def merge_intervals(
    intervals: Iterable[BlacklistInterval],
    sequence_lengths: dict[str, int],
) -> list[BlacklistInterval]:
    """Merge overlapping or adjacent intervals while retaining all reasons."""

    by_chrom: defaultdict[str, list[BlacklistInterval]] = defaultdict(list)

    for item in intervals:
        if item.chrom not in sequence_lengths:
            raise ValueError(f"unknown interval chromosome {item.chrom!r}")
        if not (0 <= item.start < item.end <= sequence_lengths[item.chrom]):
            raise ValueError(f"invalid interval {item.chrom}:{item.start}-{item.end}")
        by_chrom[item.chrom].append(item)

    merged: list[BlacklistInterval] = []

    for chrom in sequence_lengths:
        chrom_intervals = sorted(
            by_chrom.get(chrom, []),
            key=lambda item: (item.start, item.end),
        )

        for item in chrom_intervals:
            if not merged or merged[-1].chrom != chrom or item.start > merged[-1].end:
                merged.append(item)
                continue

            previous = merged[-1]
            merged[-1] = BlacklistInterval(
                chrom=chrom,
                start=previous.start,
                end=max(previous.end, item.end),
                reasons=previous.reasons | item.reasons,
            )

    return merged


def write_bed(output_path: Path, intervals: Iterable[BlacklistInterval]) -> None:
    """Write sorted blacklist intervals as BED4."""

    output_path.parent.mkdir(parents=True, exist_ok=True)

    with output_path.open("wt", encoding="utf-8") as output:
        for item in intervals:
            reasons = ",".join(sorted(item.reasons))
            output.write(f"{item.chrom}\t{item.start}\t{item.end}\t{reasons}\n")


def print_summary(
    profile_name: str,
    raw_intervals: Sequence[BlacklistInterval],
    merged_intervals: Sequence[BlacklistInterval],
    sequence_lengths: dict[str, int],
) -> None:
    """Print profile, reason counts, and merged coverage."""

    reason_counts: Counter[str] = Counter()
    for item in raw_intervals:
        reason_counts.update(item.reasons)

    print(f"Profile: {profile_name}")
    print("Input intervals by reason:")
    for reason in sorted(reason_counts):
        print(f"  {reason}: {reason_counts[reason]:,}")

    print("Merged blacklist coverage:")
    for chrom in sequence_lengths:
        chrom_intervals = [item for item in merged_intervals if item.chrom == chrom]
        if not chrom_intervals:
            continue
        masked_bases = sum(item.end - item.start for item in chrom_intervals)
        print(f"  {chrom}: intervals={len(chrom_intervals):,}, bases={masked_bases:,}")


def generate_blacklist(
    fasta_path: Path,
    gff_path: Path,
    output_path: Path,
    profile_name: str = "balanced",
) -> list[BlacklistInterval]:
    """Generate a BED4 blacklist and return its merged intervals."""

    if profile_name not in PROFILES:
        raise ValueError(f"unknown blacklist profile {profile_name!r}")
    if output_path.resolve() in {fasta_path.resolve(), gff_path.resolve()}:
        raise ValueError("output path must differ from FASTA and GFF paths")

    sequence_lengths = load_sequence_lengths(fasta_path)
    features = load_gff_features(gff_path, sequence_lengths)
    raw_intervals = collect_intervals(
        sequence_lengths,
        features,
        PROFILES[profile_name],
    )
    merged_intervals = merge_intervals(raw_intervals, sequence_lengths)
    write_bed(output_path, merged_intervals)
    print_summary(
        profile_name,
        raw_intervals,
        merged_intervals,
        sequence_lengths,
    )
    return merged_intervals


def build_parser() -> argparse.ArgumentParser:
    """Build the command-line parser."""

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--fasta",
        required=True,
        type=Path,
        help="S. pombe reference FASTA",
    )
    parser.add_argument(
        "--gff",
        required=True,
        type=Path,
        help="PomBase GFF3 annotations",
    )
    parser.add_argument(
        "--output",
        required=True,
        type=Path,
        help="Output blacklist BED4",
    )
    parser.add_argument(
        "--profile",
        choices=tuple(PROFILES),
        default="balanced",
        help="Blacklist rule profile (default: balanced)",
    )
    return parser


def main(argv: Sequence[str] | None = None) -> None:
    """Run the blacklist generator command-line interface."""

    parser = build_parser()
    arguments = parser.parse_args(argv)

    try:
        generate_blacklist(
            fasta_path=arguments.fasta,
            gff_path=arguments.gff,
            output_path=arguments.output,
            profile_name=arguments.profile,
        )
    except (OSError, ValueError) as error:
        parser.error(str(error))


if __name__ == "__main__":
    main()
