from __future__ import annotations

import io
import tempfile
import unittest
from contextlib import redirect_stderr
from pathlib import Path

from hichipbase2.digest.make_multidigest_fragments import (
    find_cut_positions,
    main,
    make_multidigest_fragments,
    parse_custom_enzyme,
    resolve_enzymes,
)


class CustomEnzymeTest(unittest.TestCase):
    def test_non_palindromic_site_is_found_on_both_strands(self) -> None:
        enzyme = parse_custom_enzyme("Custom=A^CG")

        self.assertEqual(
            find_cut_positions("TTACGAATTTCGTAA", [enzyme]),
            {3, 12},
        )

    def test_overlapping_sites_are_found(self) -> None:
        enzyme = parse_custom_enzyme("Overlap=^AAA")

        self.assertEqual(find_cut_positions("AAAA", [enzyme]), {0, 1})

    def test_invalid_definitions_and_duplicate_names_are_rejected(self) -> None:
        invalid_definitions = [
            "",
            "NoEquals",
            "=G^ATC",
            "Empty=^",
            "NoCaret=GATC",
            "TwoCarets=G^AT^C",
            "Unsupported=GA^UC",
        ]

        for definition in invalid_definitions:
            with self.subTest(definition=definition), self.assertRaises(ValueError):
                parse_custom_enzyme(definition)

        with self.assertRaisesRegex(ValueError, "Duplicate restriction enzyme"):
            resolve_enzymes(["MboI"], ["MboI=^GATC"])


class MultidigestTest(unittest.TestCase):
    def test_registered_and_custom_enzymes_write_expected_bed(self) -> None:
        enzymes = resolve_enzymes(["MboI"], ["Custom=A^ATT"])

        with tempfile.TemporaryDirectory() as temporary_directory:
            directory = Path(temporary_directory)
            fasta_path = directory / "reference.fa"
            sizes_path = directory / "chrom.sizes"
            output_path = directory / "nested" / "fragments.bed"
            fasta_path.write_text(">chr1\nAAGATCAATT\n", encoding="utf-8")
            sizes_path.write_text(
                "# test chromosomes\n\nchr1\t10\n",
                encoding="utf-8",
            )

            make_multidigest_fragments(
                fasta_path,
                sizes_path,
                output_path,
                enzymes,
            )

            self.assertEqual(
                output_path.read_text(encoding="utf-8"),
                "chr1\t0\t2\nchr1\t2\t7\nchr1\t7\t10\n",
            )

    def test_same_cut_from_multiple_enzymes_does_not_create_empty_fragment(
        self,
    ) -> None:
        enzymes = resolve_enzymes(["MboI"], ["Custom=^GATC"])

        with tempfile.TemporaryDirectory() as temporary_directory:
            directory = Path(temporary_directory)
            fasta_path = directory / "reference.fa"
            sizes_path = directory / "chrom.sizes"
            output_path = directory / "fragments.bed"
            fasta_path.write_text(">chr1\nAAGATCTT\n", encoding="utf-8")
            sizes_path.write_text("chr1\t8\n", encoding="utf-8")

            make_multidigest_fragments(
                fasta_path,
                sizes_path,
                output_path,
                enzymes,
            )

            self.assertEqual(
                output_path.read_text(encoding="utf-8"),
                "chr1\t0\t2\nchr1\t2\t8\n",
            )

    def test_missing_chromosome_and_size_mismatch_are_rejected(self) -> None:
        enzymes = resolve_enzymes(["MboI"], [])

        with tempfile.TemporaryDirectory() as temporary_directory:
            directory = Path(temporary_directory)
            fasta_path = directory / "reference.fa"
            sizes_path = directory / "chrom.sizes"
            output_path = directory / "fragments.bed"
            fasta_path.write_text(">chr1\nGATC\n", encoding="utf-8")

            sizes_path.write_text("chr2\t4\n", encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "chr2 is absent"):
                make_multidigest_fragments(
                    fasta_path,
                    sizes_path,
                    output_path,
                    enzymes,
                )

            sizes_path.write_text("chr1\t5\n", encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "FASTA length=4"):
                make_multidigest_fragments(
                    fasta_path,
                    sizes_path,
                    output_path,
                    enzymes,
                )


class CommandLineTest(unittest.TestCase):
    def test_cli_accepts_registered_and_custom_enzymes(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory:
            directory = Path(temporary_directory)
            fasta_path = directory / "reference.fa"
            sizes_path = directory / "chrom.sizes"
            output_path = directory / "fragments.bed"
            fasta_path.write_text(">chr1\nAAGATCAATT\n", encoding="utf-8")
            sizes_path.write_text("chr1\t10\n", encoding="utf-8")

            main(
                [
                    "--fasta",
                    str(fasta_path),
                    "--chrom-sizes",
                    str(sizes_path),
                    "--output",
                    str(output_path),
                    "--enzymes",
                    "MboI",
                    "--custom-enzyme",
                    "Custom=A^ATT",
                ]
            )

            self.assertEqual(
                output_path.read_text(encoding="utf-8"),
                "chr1\t0\t2\nchr1\t2\t7\nchr1\t7\t10\n",
            )

    def test_cli_rejects_missing_and_unknown_enzymes(self) -> None:
        required_arguments = [
            "--fasta",
            "reference.fa",
            "--chrom-sizes",
            "chrom.sizes",
            "--output",
            "fragments.bed",
        ]

        for enzyme_arguments, message in (
            ([], "Specify at least one enzyme"),
            (["--enzymes", "NotAnEnzyme"], "Unknown restriction enzyme"),
        ):
            stderr = io.StringIO()

            with (
                self.subTest(enzyme_arguments=enzyme_arguments),
                redirect_stderr(stderr),
                self.assertRaises(SystemExit) as error,
            ):
                main([*required_arguments, *enzyme_arguments])

            self.assertEqual(error.exception.code, 2)
            self.assertIn(message, stderr.getvalue())


if __name__ == "__main__":
    unittest.main()
