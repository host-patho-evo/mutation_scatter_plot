#!/usr/bin/env python3
# vim: set fileencoding=utf-8 ts=4 sw=4 expandtab :
import sys
import os
import re
import argparse


def main():
    parser = argparse.ArgumentParser(description="High-performance sequence regex string replacer.")
    default_tsv = os.path.join(os.path.dirname(os.path.abspath(__file__)), "s_protein_indel_misalignments.tsv")
    parser.add_argument("--rules", default=default_tsv,
                        help="Path to TSV rules. Defaults to s_protein_indel_misalignments.tsv alongside this script")
    args = parser.parse_args()

    # Pre-compile all regex rules safely
    compiled_rules = []
    try:
        with open(args.rules, "r", encoding="utf-8") as f:
            _ = next(f, None)
            for line_idx, line in enumerate(f, start=2):
                clean_line = line.rstrip('\r\n')
                if not clean_line:
                    continue
                parts = clean_line.split('\t')
                if len(parts) != 2:
                    print(f"Warning: Skipping rule on line {line_idx} (requires 2 columns)", file=sys.stderr)
                    continue

                original, resulting = parts
                # Using re.IGNORECASE to inherently match the original sed /i flag efficiency
                compiled_rules.append((re.compile(original, flags=re.IGNORECASE), resulting))
    except Exception as e:  # pylint: disable=broad-exception-caught
        print(f"Error loading rules from {args.rules}: {e}", file=sys.stderr)
        sys.exit(1)

    if not compiled_rules:
        print("Error: No regex rules loaded.", file=sys.stderr)
        sys.exit(1)

    # Fast, buffered C-level stdin line iteration
    # Because sed processed line-by-line sequentially, we strictly replicate
    # that mathematical guarantee sequentially across all rules per line.
    for line in sys.stdin:
        for pattern, replacement in compiled_rules:
            line = pattern.sub(replacement, line)
        sys.stdout.write(line)


if __name__ == "__main__":
    main()
