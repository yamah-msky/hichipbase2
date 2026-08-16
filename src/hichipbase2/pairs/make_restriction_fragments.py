from __future__ import annotations

import argparse
import re
from collections import defaultdict
from collections.abc import Sequence
from dataclasses import dataclass
from itertools import pairwise
from pathlib import Path
from typing import Final, cast

from Bio import Restriction, SeqIO
from Bio.Restriction.Restriction import AbstractCut
from Bio.Seq import Seq

IUPAC_REGEX: Final[dict[str, str]] = {
    "A": "A",
    "C": "C",
    "G": "G",
    "T": "T",
    "R": "[AG]",
    "Y": "[CT]",
    "S": "[GC]",
    "W": "[AT]",
    "K": "[GT]",
    "M": "[AC]",
    "B": "[CGT]",
    "D": "[AGT]",
    "H": "[ACT]",
    "V": "[ACG]",
    "N": "[ACGT]",
}

IUPAC_COMPLEMENT: Final[dict[str, str]] = {
    "A": "T",
    "C": "G",
    "G": "C",
    "T": "A",
    "R": "Y",
    "Y": "R",
    "S": "S",
    "W": "W",
    "K": "M",
    "M": "K",
    "B": "V",
    "D": "H",
    "H": "D",
    "V": "B",
    "N": "N",
}


@dataclass(frozen=True)
class EnzymeSpec:
    """A restriction enzyme definition used for an in silico digest."""

    name: str
    site: str
    cut_offset: int
    biopython_enzyme: type[AbstractCut] | None = None


def recognition_pattern(site: str) -> re.Pattern[str]:
    """Convert an IUPAC recognition sequence into an overlapping regex."""

    try:
        regex = "".join(IUPAC_REGEX[base] for base in site.upper())
    except KeyError as error:
        unsupported = error.args[0]
        raise ValueError(f"Unsupported IUPAC character: {unsupported}") from error

    return re.compile(f"(?=({regex}))")


def reverse_complement_site(site: str) -> str:
    """Return the reverse complement of an IUPAC recognition sequence."""

    try:
        return "".join(IUPAC_COMPLEMENT[base] for base in reversed(site.upper()))
    except KeyError as error:
        unsupported = error.args[0]
        raise ValueError(f"Unsupported IUPAC character: {unsupported}") from error


def parse_custom_enzyme(definition: str) -> EnzymeSpec:
    """Parse a custom enzyme written as ``NAME=LEFT^RIGHT``."""

    if definition.count("=") != 1:
        raise ValueError(
            "Custom enzyme must have the form NAME=LEFT^RIGHT: "
            f"{definition!r}"
        )

    name, marked_site = definition.split("=", maxsplit=1)
    name = name.strip()
    marked_site = marked_site.strip().upper()

    if not name:
        raise ValueError("Custom enzyme name must not be empty")
    if marked_site.count("^") != 1:
        raise ValueError(
            "Custom enzyme site must contain exactly one '^': "
            f"{definition!r}"
        )

    left, right = marked_site.split("^", maxsplit=1)
    site = left + right

    if not site:
        raise ValueError("Custom enzyme recognition sequence must not be empty")

    recognition_pattern(site)

    return EnzymeSpec(
        name=name,
        site=site,
        cut_offset=len(left),
    )


def resolve_builtin_enzyme(name: str) -> EnzymeSpec:
    """Resolve one enzyme name using Biopython's restriction database."""

    enzyme_name = name.strip()

    try:
        enzyme = Restriction.AllEnzymes.get(enzyme_name)
    except ValueError:
        enzyme = None

    if enzyme is None:
        raise ValueError(f"Unknown restriction enzyme: {enzyme_name!r}")
    if enzyme.fst5 is None:
        raise ValueError(
            f"Restriction enzyme {enzyme_name!r} has no known cleavage position"
        )

    typed_enzyme = cast("type[AbstractCut]", enzyme)

    return EnzymeSpec(
        name=enzyme_name,
        site=str(enzyme.site),
        cut_offset=int(enzyme.fst5),
        biopython_enzyme=typed_enzyme,
    )


def resolve_enzymes(
    builtin_names: Sequence[str],
    custom_definitions: Sequence[str],
) -> list[EnzymeSpec]:
    """Resolve and validate all requested built-in and custom enzymes."""

    enzymes = [resolve_builtin_enzyme(name) for name in builtin_names]
    enzymes.extend(parse_custom_enzyme(item) for item in custom_definitions)

    if not enzymes:
        raise ValueError(
            "Specify at least one enzyme with --enzymes or --custom-enzyme"
        )

    seen_names: set[str] = set()

    for enzyme in enzymes:
        if enzyme.name in seen_names:
            raise ValueError(f"Duplicate restriction enzyme name: {enzyme.name!r}")
        seen_names.add(enzyme.name)

    return enzymes


def find_custom_cut_positions(sequence: str, enzyme: EnzymeSpec) -> set[int]:
    """Find zero-based cut boundaries for a custom enzyme on both strands."""

    positions = {
        match.start() + enzyme.cut_offset
        for match in recognition_pattern(enzyme.site).finditer(sequence)
    }
    reverse_site = reverse_complement_site(enzyme.site)

    if reverse_site != enzyme.site:
        reverse_offset = len(enzyme.site) - enzyme.cut_offset
        positions.update(
            match.start() + reverse_offset
            for match in recognition_pattern(reverse_site).finditer(sequence)
        )

    return positions


def find_cut_sites(
    sequence: str,
    enzymes: Sequence[EnzymeSpec],
) -> dict[int, set[str]]:
    """
    Return zero-based cut boundaries and the enzymes responsible for each cut.

    Biopython reports cut positions as one-based boundaries, so registered
    enzyme positions are converted before being used as BED boundaries.
    """

    normalized_sequence = sequence.upper()
    biopython_sequence = Seq(normalized_sequence)
    cut_sources: dict[int, set[str]] = defaultdict(set)

    for enzyme in enzymes:
        if enzyme.biopython_enzyme is None:
            cut_positions = find_custom_cut_positions(normalized_sequence, enzyme)
        else:
            cut_positions = {
                position - 1
                for position in enzyme.biopython_enzyme.search(
                    biopython_sequence,
                    linear=True,
                )
            }

        for cut_position in cut_positions:
            cut_sources[cut_position].add(enzyme.name)

    return cut_sources


def digest_genome(
    fasta_path: Path,
    bed_path: Path,
    cutsite_path: Path,
    enzymes: Sequence[EnzymeSpec],
) -> None:
    """Digest every FASTA record and write BED fragments and cut annotations."""

    if not enzymes:
        raise ValueError("At least one restriction enzyme is required")

    total_fragments = 0
    total_enzyme_sites = {enzyme.name: 0 for enzyme in enzymes}

    with (
        fasta_path.open("rt", encoding="utf-8") as fasta_input,
        bed_path.open("wt", encoding="utf-8") as bed_output,
        cutsite_path.open("wt", encoding="utf-8") as cutsite_output,
    ):
        cutsite_output.write("chrom\tcut_position\tenzymes\n")

        for record in SeqIO.parse(fasta_input, "fasta"):
            chromosome = record.id
            sequence = str(record.seq).upper()
            chromosome_length = len(sequence)
            cut_sources = find_cut_sites(sequence, enzymes)

            for enzymes_at_site in cut_sources.values():
                for enzyme_name in enzymes_at_site:
                    total_enzyme_sites[enzyme_name] += 1

            sorted_cuts = sorted({0, chromosome_length, *cut_sources})
            chromosome_fragment_count = 0

            for start, end in pairwise(sorted_cuts):
                if start >= end:
                    continue

                bed_output.write(f"{chromosome}\t{start}\t{end}\n")
                chromosome_fragment_count += 1
                total_fragments += 1

            for cut_position in sorted(cut_sources):
                enzyme_names = ",".join(sorted(cut_sources[cut_position]))
                cutsite_output.write(
                    f"{chromosome}\t{cut_position}\t{enzyme_names}\n"
                )

            print(
                f"{chromosome}: "
                f"length={chromosome_length:,}, "
                f"cut_sites={len(cut_sources):,}, "
                f"fragments={chromosome_fragment_count:,}"
            )

    print()
    print("Recognition-site counts:")

    for enzyme in enzymes:
        print(
            f"  {enzyme.name}: "
            f"site={enzyme.site}, "
            f"offset={enzyme.cut_offset}, "
            f"count={total_enzyme_sites[enzyme.name]:,}"
        )

    print(f"Total fragments: {total_fragments:,}")
    print(f"BED: {bed_path}")
    print(f"Cut-site table: {cutsite_path}")


def build_parser() -> argparse.ArgumentParser:
    """Build the command-line argument parser."""

    parser = argparse.ArgumentParser(
        description=(
            "Digest a reference genome in silico with one or more restriction "
            "enzymes."
        )
    )
    parser.add_argument(
        "--fasta",
        required=True,
        type=Path,
        help="Reference genome FASTA",
    )
    parser.add_argument(
        "--bed",
        required=True,
        type=Path,
        help="Output restriction-fragment BED",
    )
    parser.add_argument(
        "--cuts",
        required=True,
        type=Path,
        help="Output cut-site annotation TSV",
    )
    parser.add_argument(
        "--enzymes",
        nargs="+",
        default=[],
        metavar="NAME",
        help="Biopython restriction enzyme name(s), for example MboI HinfI",
    )
    parser.add_argument(
        "--custom-enzyme",
        action="append",
        default=[],
        metavar="NAME=LEFT^RIGHT",
        help=(
            "Custom IUPAC recognition site and cut boundary; may be specified "
            "more than once"
        ),
    )

    return parser


def main(argv: Sequence[str] | None = None) -> None:
    """Run the restriction-fragment command-line interface."""

    parser = build_parser()
    arguments = parser.parse_args(argv)

    try:
        enzymes = resolve_enzymes(arguments.enzymes, arguments.custom_enzyme)
    except ValueError as error:
        parser.error(str(error))

    digest_genome(
        fasta_path=arguments.fasta,
        bed_path=arguments.bed,
        cutsite_path=arguments.cuts,
        enzymes=enzymes,
    )


if __name__ == "__main__":
    main()
