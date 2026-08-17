from __future__ import annotations

import io
import tempfile
import unittest
from contextlib import redirect_stdout
from pathlib import Path

from hichipbase2.blacklist.make_spombe_blacklist import (
    PROFILES,
    BlacklistInterval,
    build_parser,
    collect_intervals,
    feature_reason,
    generate_blacklist,
    load_gff_features,
    load_sequence_lengths,
    merge_intervals,
)


def write_reference(directory: Path) -> tuple[Path, Path]:
    fasta_path = directory / "genome.fa"
    gff_path = directory / "genes.gff3"
    fasta_path.write_text(
        "".join(
            f">{chrom}\n{'A' * length}\n"
            for chrom, length in (
                ("I", 200_000),
                ("II", 200_000),
                ("III", 200_000),
                ("mitochondrial", 100),
            )
        ),
        encoding="utf-8",
    )
    gff_path.write_text(
        """##gff-version 3
I\tPomBase\tregional_centromere\t80001\t90000\t.\t+\t.\tID=cen1
I\tPomBase\tlong_terminal_repeat\t70001\t70300\t.\t+\t.\tID=ltr1
I\tPomBase\trepeat_region\t100001\t100100\t.\t+\t.\tID=dup;note=duplicated region
I\tPomBase\trepeat_region\t110001\t110100\t.\t+\t.\tID=simple;note=(at)50
II\tPomBase\tgap\t60001\t60100\t.\t+\t.\tID=gap1
II\tPomBase\tmating_type_region\t100001\t101000\t.\t+\t.\tID=mat
II\tPomBase\tnuclear_mt_pseudogene\t110001\t110100\t.\t+\t.\tID=numt
III\tPomBase\tgene_group\t1\t1000\t.\t-\t.\tID=left-rdna
III\tPomBase\tgene_group\t195001\t200000\t.\t+\t.\tID=right-rdna
III\tPomBase\trRNA\t70001\t70200\t.\t+\t.\tID=rrna
mitochondrial\tPomBase\trRNA\t1\t100\t.\t+\t.\tID=mt-rrna
""",
        encoding="utf-8",
    )
    return fasta_path, gff_path


class ProfileTest(unittest.TestCase):
    def test_balanced_is_the_cli_default(self) -> None:
        parser = build_parser()
        arguments = parser.parse_args(
            ["--fasta", "genome.fa", "--gff", "genes.gff3", "--output", "x.bed"]
        )

        self.assertEqual(arguments.profile, "balanced")
        self.assertEqual(set(PROFILES), {"core", "balanced", "comprehensive"})

    def test_profiles_select_expected_features_and_terminal_regions(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory:
            directory = Path(temporary_directory)
            fasta_path, gff_path = write_reference(directory)
            lengths = load_sequence_lengths(fasta_path)
            features = load_gff_features(gff_path, lengths)

            core = collect_intervals(lengths, features, PROFILES["core"])
            balanced = collect_intervals(lengths, features, PROFILES["balanced"])
            comprehensive = collect_intervals(
                lengths,
                features,
                PROFILES["comprehensive"],
            )

            core_reasons = {reason for item in core for reason in item.reasons}
            balanced_reasons = {reason for item in balanced for reason in item.reasons}
            comprehensive_reasons = {
                reason for item in comprehensive for reason in item.reasons
            }

            self.assertNotIn("long_terminal_repeat", core_reasons)
            self.assertNotIn("rrna", core_reasons)
            self.assertIn("long_terminal_repeat", balanced_reasons)
            self.assertIn("annotated_repeat", balanced_reasons)
            self.assertNotIn("repeat_region", balanced_reasons)
            self.assertIn("repeat_region", comprehensive_reasons)
            self.assertIn(
                BlacklistInterval(
                    "III",
                    0,
                    50_000,
                    frozenset({"subtelomere"}),
                ),
                comprehensive,
            )

    def test_balanced_repeat_selection_uses_annotation_keywords(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory:
            directory = Path(temporary_directory)
            fasta_path, gff_path = write_reference(directory)
            lengths = load_sequence_lengths(fasta_path)
            features = load_gff_features(gff_path, lengths)
            repeat_features = [
                feature
                for feature in features
                if feature.feature_type == "repeat_region"
            ]

            self.assertEqual(
                feature_reason(repeat_features[0], PROFILES["balanced"]),
                "annotated_repeat",
            )
            self.assertIsNone(feature_reason(repeat_features[1], PROFILES["balanced"]))


class BlacklistGenerationTest(unittest.TestCase):
    def test_balanced_generation_writes_merged_bed4_and_summary(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory:
            directory = Path(temporary_directory)
            fasta_path, gff_path = write_reference(directory)
            output_path = directory / "results" / "blacklist.bed"
            stdout = io.StringIO()

            with redirect_stdout(stdout):
                intervals = generate_blacklist(
                    fasta_path,
                    gff_path,
                    output_path,
                )

            bed = output_path.read_text(encoding="utf-8")
            self.assertIn("I\t0\t50000\tsubtelomere\n", bed)
            self.assertIn("III\t0\t1000\trdna_array\n", bed)
            self.assertIn(
                "mitochondrial\t0\t100\tnon_primary_contig\n",
                bed,
            )
            self.assertTrue(intervals)
            self.assertIn("Profile: balanced", stdout.getvalue())
            self.assertIn("Merged blacklist coverage:", stdout.getvalue())

    def test_merge_combines_adjacent_intervals_and_reasons(self) -> None:
        merged = merge_intervals(
            [
                BlacklistInterval("I", 10, 20, frozenset({"first"})),
                BlacklistInterval("I", 20, 30, frozenset({"second"})),
            ],
            {"I": 100},
        )

        self.assertEqual(
            merged,
            [
                BlacklistInterval(
                    "I",
                    10,
                    30,
                    frozenset({"first", "second"}),
                )
            ],
        )

    def test_gff_sequence_and_coordinate_mismatches_are_rejected(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory:
            directory = Path(temporary_directory)
            fasta_path, gff_path = write_reference(directory)
            lengths = load_sequence_lengths(fasta_path)

            gff_path.write_text(
                "unknown\tPomBase\tgap\t1\t10\t.\t+\t.\tID=gap\n",
                encoding="utf-8",
            )
            with self.assertRaisesRegex(ValueError, "absent from the FASTA"):
                load_gff_features(gff_path, lengths)

            gff_path.write_text(
                "I\tPomBase\tgap\t1\t200001\t.\t+\t.\tID=gap\n",
                encoding="utf-8",
            )
            with self.assertRaisesRegex(ValueError, "outside I"):
                load_gff_features(gff_path, lengths)


if __name__ == "__main__":
    unittest.main()
