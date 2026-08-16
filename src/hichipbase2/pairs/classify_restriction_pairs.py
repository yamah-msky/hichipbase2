from __future__ import annotations

import argparse
import gzip
from collections import Counter
from contextlib import ExitStack
from pathlib import Path
from typing import TextIO

from hichipbase2.util import is_gzip

CLASS_NAMES = (
    "dangling_end",
    "self_circle",
    "re_ligation",
    "valid_interaction",
    "unclassified",
)

def open_text(path: Path, mode: str) -> TextIO:
    if is_gzip(path):
        return gzip.open(path, mode, encoding="utf-8")
    return path.open(mode, encoding="utf-8")

def endpoint(
    record: dict[str, str],
    side: int,
) -> dict[str, str | int ]:
    return {
        "chromosome": record[f"chrom{side}"],
        "position": int(record[f"pos{side}"]),
        "strand": record[f"strand{side}"],
        "fragment": int(record[f"rfrag{side}"]),
        "fragment_start": int(record[f"rfrag_start{side}"]),
        "fragment_end": int(record[f"rfrag_end{side}"]),
    }

def classify_pair(record: dict[str, str]) -> str:
    first = endpoint(record, 1)
    second = endpoint(record, 2)

    if first["fragment"] < 0 or second["fragment"] < 0:
        return "unclassified"
    
    if first["chromosome"] != second["chromosome"]:
        return "valid_interaction"
    
    if first["position"] > second["position"]:
        first, second = second, first
    
    same_fragment = (
        first["fragment"] == second["fragment"]
        and first["fragment_start"] == second["fragment_start"]
        and first["fragment_end"] == second["fragment_end"]
    )

    if same_fragment:
        orientation = (
            str(first["strand"]),
            str(second["strand"]),
        )

        if orientation == ("+", "-"):
            return "dangling_end"
        
        if orientation == ("-", "+"):
            return "self_circle"
        
        return "unclassified"

    directly_adjacent = (
        first["fragment_end"] == second["fragment_start"]
        or second["fragment_end"] == first["fragment_start"]
    )

    if directly_adjacent:
        return "re_ligation"
    
    return "valid_interaction"

def classify_file(
    input_path: Path,
    output_prefix: Path,
) -> None:
    output_prefix.parent.mkdir(parents=True, exist_ok=True)

    counts: Counter[str] = Counter()
    columns: list[str] | None = None
    
    with ExitStack() as stack:
        input_handle = stack.enter_context(
            open_text(input_path, "rt")
        )

        combined_output = stack.enter_context(
            gzip.open(
                f"{output_prefix}.classified.pairs.gz",
                "wt",
                encoding="utf-8",
            )
        )

        class_outputs = {
            class_name: stack.enter_context(
                gzip.open(
                    f"{output_prefix}.{class_name}.pairs.gz",
                    "wt",
                    encoding="utf-8",
                )
            )
            for class_name in CLASS_NAMES
        }

        for line in input_handle:
            if line.startswith("#columns:"):
                columns = (
                    line.removeprefix("#columns:")
                    .strip()
                    .split()
                )

                new_header = (
                    "#columns: "
                    + " ".join(columns)
                    + " interaction_type\n"
                )

                combined_output.write(new_header)

                for output_handle in class_outputs.values():
                    output_handle.write(new_header)
                
                continue
            
            if line.startswith("#"):
                combined_output.write(line)

                for output_handle in class_outputs.values():
                    output_handle.write(line)
                
                continue
            
            if columns is None:
                raise ValueError(
                    "Input file has no #columns header"
                )

            values = line.rstrip("\n").split("\t")

            if len(values) != len(columns):
                raise ValueError(
                    "Column count mismatch: "
                    f"expected {len(columns)}, "
                    f"observed {len(values)}"
                )

            record = dict(zip(columns, values))
            interaction_type = classify_pair(record)


            output_line = (
                line.rstrip("\n")
                + f"\t{interaction_type}\n"
            )

            combined_output.write(output_line)
            class_outputs[interaction_type].write(output_line)

            counts[interaction_type] += 1

    total = sum(counts.values())

    stats_path = Path(
        f"{output_prefix}.classification.tsv"
    )

    with stats_path.open("wt", encoding="utf-8") as stats:
        stats.write(
            "interaction_type\ncount\tfraction\n"
        )

        for class_name in CLASS_NAMES:
            count = counts[class_name]
            fraction = count / total if total else 0.0

            stats.write(
                f"{class_name}\t"
                f"{count}\t"
                f"{fraction:.8f}\n"
            )

    print(f"Total classified pairs: {total:,}")

    for class_name in CLASS_NAMES:
        count = counts[class_name]
        fraction = count / total if total else 0.0

        print(
            f"{class_name}: "
            f"{count:,} "
            f"({fraction:.2%})"
        )

def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Classify restriction-fragment-annotated "
            "Hi-C/HiChIP pairs."
        )
    )

    parser.add_argument(
        "--input",
        required=True,
        type=Path,
    )
    parser.add_argument(
        "--output-prefix",
        required=True,
        type=Path,
    )

    arguments = parser.parse_args()

    classify_file(
        input_path=arguments.input,
        output_prefix=arguments.output_prefix,
    )

if __name__ == "__main__":
    main()