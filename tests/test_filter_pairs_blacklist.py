from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

from hichipbase2.blacklist.filter_pairs_blacklist import filter_pairs


class FilterPairsBoundaryTest(unittest.TestCase):
    def test_bed4_boundaries_are_applied_to_one_based_pairs_positions(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory:
            directory = Path(temporary_directory)
            input_path = directory / "input.pairs"
            output_path = directory / "output.pairs"
            blacklist_path = directory / "blacklist.bed"
            sizes_path = directory / "chrom.sizes"
            stats_path = directory / "stats.tsv"
            input_path.write_text(
                """## pairs format v1.0
#chromsize: I 100
#columns: readID chrom1 pos1 chrom2 pos2 strand1 strand2 pair_type
before\tI\t10\tI\t90\t+\t-\tUU
first\tI\t11\tI\t90\t+\t-\tUU
last\tI\t20\tI\t90\t+\t-\tUU
after\tI\t21\tI\t90\t+\t-\tUU
""",
                encoding="utf-8",
            )
            blacklist_path.write_text(
                "I\t10\t20\ttest_reason\n",
                encoding="utf-8",
            )
            sizes_path.write_text("I\t100\n", encoding="utf-8")

            counts = filter_pairs(
                input_path,
                output_path,
                blacklist_path,
                sizes_path,
                stats_path,
            )

            output = output_path.read_text(encoding="utf-8")
            self.assertIn("before\tI\t10", output)
            self.assertNotIn("first\tI\t11", output)
            self.assertNotIn("last\tI\t20", output)
            self.assertIn("after\tI\t21", output)
            self.assertEqual(counts["input_pairs"], 4)
            self.assertEqual(counts["kept_pairs"], 2)
            self.assertEqual(counts["blacklisted_pairs"], 2)


if __name__ == "__main__":
    unittest.main()
