from __future__ import annotations

import argparse
import re
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
    """A built-in or custom restriction enzyme definition."""

    name: str
    site: str
    cut_offset: int
    biopython_enzyme: type[AbstractCut] | None = None


def recognition_pattern(site: str) -> re.Pattern[str]:
    """Convert an IUPAC recognition sequence to an overlapping regex."""

    try:
        regex = "".join(IUPAC_REGEX[base] for base in site.upper())
    except KeyError as error:
        raise ValueError(
            f"Unsupported IUPAC character: {error.args[0]}"
        ) from error

    return re.compile(f"(?=({regex}))")


def reverse_complement_site(site: str) -> str:
    """Return the reverse complement of an IUPAC recognition sequence."""

    try:
        return "".join(
            IUPAC_COMPLEMENT[base] for base in reversed(site.upper())
        )
    except KeyError as error:
        raise ValueError(
            f"Unsupported IUPAC character: {error.args[0]}"
        ) from error


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

    return EnzymeSpec(name=name, site=site, cut_offset=len(left))


def resolve_builtin_enzyme(name: str) -> EnzymeSpec:
    """Resolve an enzyme name using Biopython's restriction database."""

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
    """Resolve and validate all requested enzymes."""

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
    """Find custom-enzyme cut boundaries on both strands."""

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


def find_cut_positions(
    sequence: str,
    enzymes: Sequence[EnzymeSpec],
) -> set[int]:
    """Return all zero-based cut boundaries in a linear sequence."""

    normalized_sequence = sequence.upper()
    biopython_sequence = Seq(normalized_sequence)
    positions: set[int] = set()

    for enzyme in enzymes:
        if enzyme.biopython_enzyme is None:
            positions.update(find_custom_cut_positions(normalized_sequence, enzyme))
        else:
            positions.update(
                position - 1
                for position in enzyme.biopython_enzyme.search(
                    biopython_sequence,
                    linear=True,
                )
            )

    return positions


def make_multidigest_fragments(
    fasta_path: Path,
    chrom_sizes_path: Path,
    output_path: Path,
    enzymes: Sequence[EnzymeSpec],
) -> None:
    """Digest chromosomes and write their restriction fragments as BED."""

    if not enzymes:
        raise ValueError("At least one restriction enzyme is required")

    with fasta_path.open("rt", encoding="utf-8") as fasta_input:
        records = SeqIO.to_dict(SeqIO.parse(fasta_input, "fasta"))

    output_path.parent.mkdir(parents=True, exist_ok=True)

    with (
        chrom_sizes_path.open("rt", encoding="utf-8") as sizes_input,
        output_path.open("wt", encoding="utf-8") as output,
    ):
        for line in sizes_input:
            if not line.strip() or line.startswith("#"):
                continue

            fields = line.split()

            if len(fields) < 2:
                raise ValueError(f"Invalid chrom.sizes line: {line.rstrip()!r}")

            chromosome, size_text = fields[:2]

            try:
                expected_size = int(size_text)
            except ValueError as error:
                raise ValueError(
                    f"Invalid chromosome size for {chromosome!r}: {size_text!r}"
                ) from error

            if chromosome not in records:
                raise ValueError(f"{chromosome} is absent from {fasta_path}")

            sequence = records[chromosome].seq

            if len(sequence) != expected_size:
                raise ValueError(
                    f"{chromosome}: FASTA length={len(sequence)}, "
                    f"chrom.sizes length={expected_size}"
                )

            cut_positions = find_cut_positions(str(sequence), enzymes)
            boundaries = sorted(
                {
                    0,
                    expected_size,
                    *(
                        position
                        for position in cut_positions
                        if 0 <= position <= expected_size
                    ),
                }
            )

            for start, end in pairwise(boundaries):
                if start < end:
                    output.write(f"{chromosome}\t{start}\t{end}\n")


def build_parser() -> argparse.ArgumentParser:
    """Build the command-line argument parser."""

    parser = argparse.ArgumentParser(
        description=(
            "Digest a reference genome with multiple restriction enzymes and "
            "write BED fragments."
        )
    )
    parser.add_argument(
        "--fasta",
        required=True,
        type=Path,
        help="Reference genome FASTA",
    )
    parser.add_argument(
        "--chrom-sizes",
        required=True,
        type=Path,
        help="Chromosome sizes file",
    )
    parser.add_argument(
        "--output",
        required=True,
        type=Path,
        help="Output restriction-fragment BED",
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
    """Run the multidigest-fragment command-line interface."""

    parser = build_parser()
    arguments = parser.parse_args(argv)

    try:
        enzymes = resolve_enzymes(arguments.enzymes, arguments.custom_enzyme)
    except ValueError as error:
        parser.error(str(error))

    make_multidigest_fragments(
        fasta_path=arguments.fasta,
        chrom_sizes_path=arguments.chrom_sizes,
        output_path=arguments.output,
        enzymes=enzymes,
    )


if __name__ == "__main__":
    main()
