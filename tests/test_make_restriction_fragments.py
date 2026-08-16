from __future__ import annotations

import io
import tempfile
import unittest
from contextlib import redirect_stderr, redirect_stdout
from pathlib import Path

from hichipbase2.pairs.make_restriction_fragments import (
    digest_genome,
    find_cut_sites,
    parse_custom_enzyme,
    resolve_enzymes,
)


class RestrictionEnzymeTest(unittest.TestCase):
    def test_registered_type_iis_enzyme_uses_biopython_cut_position(self) -> None:
        enzymes = resolve_enzymes(["BsaI"], [])

        self.assertEqual(
            find_cut_sites("AAGGTCTCAAAAATT", enzymes),
            {9: {"BsaI"}},
        )

    def test_custom_non_palindromic_site_is_found_on_both_strands(self) -> None:
        enzyme = parse_custom_enzyme("Custom=A^CG")

        self.assertEqual(
            find_cut_sites("TTACGAATTTCGTAA", [enzyme]),
            {3: {"Custom"}, 12: {"Custom"}},
        )

    def test_overlapping_custom_sites_are_all_found(self) -> None:
        enzyme = parse_custom_enzyme("Overlap=^AAA")

        self.assertEqual(
            find_cut_sites("AAAA", [enzyme]),
            {0: {"Overlap"}, 1: {"Overlap"}},
        )

    def test_enzyme_names_are_combined_at_the_same_cut(self) -> None:
        enzymes = resolve_enzymes(["MboI"], ["Custom=^GATC"])

        self.assertEqual(
            find_cut_sites("AAGATCTT", enzymes),
            {2: {"Custom", "MboI"}},
        )

    def test_invalid_enzyme_definitions_are_rejected(self) -> None:
        invalid_custom_definitions = [
            "",
            "NoEquals",
            "=G^ATC",
            "Empty=^",
            "NoCaret=GATC",
            "TwoCarets=G^AT^C",
            "Unsupported=GA^UC",
        ]

        for definition in invalid_custom_definitions:
            with self.subTest(definition=definition), self.assertRaises(ValueError):
                parse_custom_enzyme(definition)

        with self.assertRaisesRegex(ValueError, "at least one enzyme"):
            resolve_enzymes([], [])
        with self.assertRaisesRegex(ValueError, "Unknown restriction enzyme"):
            resolve_enzymes(["NotAnEnzyme"], [])
        with self.assertRaisesRegex(ValueError, "no known cleavage position"):
            resolve_enzymes(["Asu14238IV"], [])
        with self.assertRaisesRegex(ValueError, "Duplicate restriction enzyme"):
            resolve_enzymes(["MboI"], ["MboI=^GATC"])


class DigestGenomeTest(unittest.TestCase):
    def test_digest_writes_fragments_cuts_and_selected_enzyme_counts(self) -> None:
        enzymes = resolve_enzymes(["MboI"], ["Custom=A^ATT"])

        with tempfile.TemporaryDirectory() as temporary_directory:
            directory = Path(temporary_directory)
            fasta_path = directory / "reference.fa"
            bed_path = directory / "fragments.bed"
            cuts_path = directory / "cuts.tsv"
            fasta_path.write_text(">chr1\nAAGATCAATT\n", encoding="utf-8")
            stdout = io.StringIO()

            with redirect_stdout(stdout):
                digest_genome(fasta_path, bed_path, cuts_path, enzymes)

            self.assertEqual(
                bed_path.read_text(encoding="utf-8"),
                "chr1\t0\t2\nchr1\t2\t7\nchr1\t7\t10\n",
            )
            self.assertEqual(
                cuts_path.read_text(encoding="utf-8"),
                (
                    "chrom\tcut_position\tenzymes\n"
                    "chr1\t2\tMboI\n"
                    "chr1\t7\tCustom\n"
                ),
            )
            self.assertIn("MboI: site=GATC, offset=0, count=1", stdout.getvalue())
            self.assertIn(
                "Custom: site=AATT, offset=1, count=1",
                stdout.getvalue(),
            )
            self.assertNotIn("MluCI", stdout.getvalue())

    def test_empty_enzyme_list_is_rejected(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory:
            directory = Path(temporary_directory)

            with self.assertRaisesRegex(ValueError, "At least one"):
                digest_genome(
                    directory / "reference.fa",
                    directory / "fragments.bed",
                    directory / "cuts.tsv",
                    [],
                )


class CommandLineTest(unittest.TestCase):
    def test_cli_accepts_builtin_and_custom_enzymes(self) -> None:
        from hichipbase2.pairs.make_restriction_fragments import main

        with tempfile.TemporaryDirectory() as temporary_directory:
            directory = Path(temporary_directory)
            fasta_path = directory / "reference.fa"
            bed_path = directory / "fragments.bed"
            cuts_path = directory / "cuts.tsv"
            fasta_path.write_text(">chr1\nAAGATCAATT\n", encoding="utf-8")

            with redirect_stdout(io.StringIO()):
                main(
                    [
                        "--fasta",
                        str(fasta_path),
                        "--bed",
                        str(bed_path),
                        "--cuts",
                        str(cuts_path),
                        "--enzymes",
                        "MboI",
                        "--custom-enzyme",
                        "Custom=A^ATT",
                    ]
                )

            self.assertEqual(
                bed_path.read_text(encoding="utf-8"),
                "chr1\t0\t2\nchr1\t2\t7\nchr1\t7\t10\n",
            )

    def test_cli_requires_an_enzyme(self) -> None:
        from hichipbase2.pairs.make_restriction_fragments import main

        stderr = io.StringIO()

        with redirect_stderr(stderr), self.assertRaises(SystemExit) as error:
            main(
                [
                    "--fasta",
                    "reference.fa",
                    "--bed",
                    "fragments.bed",
                    "--cuts",
                    "cuts.tsv",
                ]
            )

        self.assertEqual(error.exception.code, 2)
        self.assertIn("Specify at least one enzyme", stderr.getvalue())


if __name__ == "__main__":
    unittest.main()
