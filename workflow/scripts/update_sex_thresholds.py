#!/usr/bin/env python3

"""Update sex-QC F thresholds in the workflow YAML config.

This intentionally edits only the female_max_f and male_min_f entries inside
sex_qc, preserving the rest of the config file and its comments/formatting.
"""

from __future__ import annotations

import argparse
import math
import re
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", required=True, help="Workflow YAML config")
    parser.add_argument("--female-max-f", required=True, type=float)
    parser.add_argument("--male-min-f", required=True, type=float)
    return parser.parse_args()


def main() -> None:
    args = parse_args()

    if not math.isfinite(args.female_max_f) or not math.isfinite(args.male_min_f):
        raise SystemExit("Sex-QC thresholds must be finite numbers.")

    if args.female_max_f >= args.male_min_f:
        raise SystemExit("female_max_f must be lower than male_min_f.")

    path = Path(args.config)
    if not path.is_file():
        raise SystemExit(f"Config file does not exist: {path}")

    lines = path.read_text().splitlines(keepends=True)

    in_sex_qc = False
    found_section = False
    found_female = False
    found_male = False
    updated: list[str] = []

    for line in lines:
        if re.match(r"^sex_qc:\s*(?:#.*)?$", line.rstrip("\n")):
            in_sex_qc = True
            found_section = True
            updated.append(line)
            continue

        if in_sex_qc and re.match(r"^[^\s#][^:]*:\s*", line):
            in_sex_qc = False

        if in_sex_qc and re.match(r"^\s+female_max_f\s*:", line):
            newline = "\n" if line.endswith("\n") else ""
            updated.append(f"  female_max_f: {args.female_max_f:g}{newline}")
            found_female = True
            continue

        if in_sex_qc and re.match(r"^\s+male_min_f\s*:", line):
            newline = "\n" if line.endswith("\n") else ""
            updated.append(f"  male_min_f: {args.male_min_f:g}{newline}")
            found_male = True
            continue

        updated.append(line)

    if not found_section:
        raise SystemExit("Config is missing the 'sex_qc:' section.")
    if not found_female:
        raise SystemExit("Config is missing 'sex_qc/female_max_f'.")
    if not found_male:
        raise SystemExit("Config is missing 'sex_qc/male_min_f'.")

    path.write_text("".join(updated))

    print("Updated sex-QC thresholds:")
    print(f"  female_max_f: {args.female_max_f:g}")
    print(f"  male_min_f:   {args.male_min_f:g}")


if __name__ == "__main__":
    main()
