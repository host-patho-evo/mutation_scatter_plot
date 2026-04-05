#!/usr/bin/env python3
# vim: set fileencoding=utf-8 ts=4 sw=4 expandtab :
import sys
import os
import re
import argparse
import multiprocessing
from itertools import islice

# Global variable for workers to store compiled rules without IPC overhead
_COMPILED_RULES = []


def _init_worker(rules_path):
    """
    Initializes each worker process by dynamically compiling the TSV regex rules.
    This entirely prevents cross-process thread-safety and IPC pickling overhead.
    """
    global _COMPILED_RULES  # pylint: disable=global-statement
    _COMPILED_RULES = []
    try:
        with open(rules_path, "r", encoding="utf-8") as f:
            _ = next(f, None)
            for line_idx, line in enumerate(f, start=2):
                clean_line = line.rstrip('\r\n')
                if not clean_line:
                    continue
                parts = clean_line.split('\t')
                if len(parts) == 2:
                    original, resulting = parts
                    # We implicitly deploy robust case-insensitive alignment safely simulating sed /i
                    _COMPILED_RULES.append((re.compile(original, flags=re.IGNORECASE), resulting))
                else:
                    print(f"Warning: Skipping rule on line {line_idx} (requires 2 columns)", file=sys.stderr)
    except Exception as e:  # pylint: disable=broad-exception-caught
        print(f"Worker Error: Failed to load rules from {rules_path}: {e}", file=sys.stderr)


def _process_chunk(chunk):
    """
    Processes a rigid multi-line chunk sequentially matching bash strictness structurally.
    """
    results = []
    for line in chunk:
        for pattern, replacement in _COMPILED_RULES:
            line = pattern.sub(replacement, line)
        results.append(line)
    return results


def get_chunks(iterable, size=10000):
    """
    Safely yields chunks of sizes dynamically buffered structurally bypassing RAM explosions.
    """
    it = iter(iterable)
    while True:
        chunk = list(islice(it, size))
        if not chunk:
            break
        yield chunk


def main():
    parser = argparse.ArgumentParser(description="High-performance multi-core sequence regex string replacer.")
    default_tsv = os.path.join(
        os.path.dirname(os.path.abspath(__file__)),
        "manual_SARS-CoV-2_InDel_realignment_rules.tsv"
    )
    parser.add_argument(
        "--rules", default=default_tsv,
        help="Path to TSV rules. Defaults to manual_SARS-CoV-2_InDel_realignment_rules.tsv alongside this script"
    )
    parser.add_argument(
        "--cores", type=int, default=1,
        help="Number of multi-processing CPU cores to structurally allocate (default: 1)"
    )
    args = parser.parse_args()

    # Pre-flight check exclusively to ensure the system is structurally sound before deploying multiprocessing pools
    if not os.path.isfile(args.rules):
        print(f"Error: Rules TSV file not found at {args.rules}", file=sys.stderr)
        sys.exit(1)

    # Secure evaluation of execution logic intelligently
    # Single-core fallback is kept structurally efficient without multiprocessing mapping matrices physically
    if args.cores <= 1:
        _init_worker(args.rules)
        if not _COMPILED_RULES:
            print("Error: No regex rules loaded.", file=sys.stderr)
            sys.exit(1)

        for line in sys.stdin:
            for pattern, replacement in _COMPILED_RULES:
                line = pattern.sub(replacement, line)
            sys.stdout.write(line)
    else:
        # Utilizing a multiprocessing configuration correctly
        # `imap` intrinsically guarantees mathematical sequential result parity identically matched to original stdout
        # chunksize=1 here means 1 outer chunk (which contains 10000 lines internally) natively dispatched
        with multiprocessing.Pool(
            processes=args.cores,
            initializer=_init_worker,
            initargs=(args.rules,)
        ) as pool:
            for processed_chunk in pool.imap(_process_chunk, get_chunks(sys.stdin, size=10000), chunksize=1):
                sys.stdout.writelines(processed_chunk)


if __name__ == "__main__":
    main()
