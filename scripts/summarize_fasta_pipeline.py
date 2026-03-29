#!/usr/bin/env python3
"""Summarize FASTA pipeline: record counts and NNNNx count sums across pipeline stages.

For each file found under <search_path> matching <filename_prefix>*.fasta{,.old,.ori,.orig},
computes:
  - Number of FASTA records  (grep -c '^>')
  - Sum of NNNNx count prefixes in the FASTA IDs
      (non-NNNNx IDs, e.g. GISAID accessions, contribute 0 to the sum)

Then prints a table with per-step deltas so you can track how many
records and sequence counts are gained or lost at each pipeline stage.

Usage:
    summarize_fasta_pipeline.py <search_path> <filename_prefix>

Example:
    summarize_fasta_pipeline.py . spikenuc1207.native2ascii.no_junk
    summarize_fasta_pipeline.py .. spikenuc1207.native2ascii.no_junk
"""

import glob
import os
import subprocess
import sys

FASTA_SUFFIXES = ('.fasta', '.fasta.old', '.fasta.ori', '.fasta.orig')

# ── helpers ───────────────────────────────────────────────────────────────────

def _count_records(path: str) -> int:
    """Number of '>' header lines in a FASTA file."""
    result = subprocess.run(
        ['grep', '-c', '^>', path],
        capture_output=True, text=True,
    )
    # grep exits 1 when there are zero matches; that is still valid.
    text = result.stdout.strip()
    return int(text) if text else 0


def _sum_nnnx_counts(path: str) -> int:
    """Sum of the leading NNNNx integer prefixes across all FASTA header IDs.

    Equivalent to:
        grep '^>' $file | cut -c 2- | awk '{print $1}' | sed -e 's/x.*//' \\
            | awk '{SUM += $1} END {print SUM}'

    Non-NNNNx IDs (e.g. GISAID accessions that contain no 'x') contribute 0.
    """
    cmd = (
        "grep '^>' " + _shell_quote(path) + r" | cut -c 2-"
        r" | awk '{print $1}'"
        r" | sed -e 's/x.*//'"
        r" | awk '{SUM += $1} END {print SUM+0}'"
    )
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    text = result.stdout.strip()
    return int(text) if text else 0


def _shell_quote(s: str) -> str:
    """Minimal single-quote escaping for shell embedding."""
    return "'" + s.replace("'", "'\\''") + "'"


def _delta_str(current: int, previous: int) -> str:
    diff = current - previous
    sign = '+' if diff >= 0 else ''
    return f"{sign}{diff:,}"


def _pct_str(current: int, reference: int) -> str:
    if reference == 0:
        return "n/a"
    return f"{current / reference * 100:.2f}%"


# ── main ──────────────────────────────────────────────────────────────────────

def main() -> None:
    if len(sys.argv) < 3:
        print(__doc__, file=sys.stderr)
        sys.exit(1)

    search_path = sys.argv[1]
    prefix = sys.argv[2]

    # Collect matching files (any depth under search_path).
    found: set[str] = set()
    for suffix in FASTA_SUFFIXES:
        pattern = os.path.join(search_path, '**', prefix + '*' + suffix)
        found.update(glob.glob(pattern, recursive=True))
        # Also search non-recursively at the top level.
        pattern_flat = os.path.join(search_path, prefix + '*' + suffix)
        found.update(glob.glob(pattern_flat))

    if not found:
        print(f"No files found under '{search_path}' matching '{prefix}*.fasta{{,.old,.ori,.orig}}'",
              file=sys.stderr)
        sys.exit(1)

    files = sorted(found)
    print(f"Found {len(files):,} file(s).\n", file=sys.stderr)

    # ── gather data ───────────────────────────────────────────────────────────
    rows: list[tuple[str, int, int]] = []  # (display_name, n_records, n_sum)
    for f in files:
        display = os.path.relpath(f, search_path)
        print(f"  scanning {display} …", file=sys.stderr, flush=True)
        n_rec = _count_records(f)
        n_sum = _sum_nnnx_counts(f)
        rows.append((display, n_rec, n_sum))

    # ── compute column widths ─────────────────────────────────────────────────
    col_file = max(len(r[0]) for r in rows)
    col_file = max(col_file, len("File"))

    # ── print table ───────────────────────────────────────────────────────────
    SEP = "  "
    W_NUM = 14   # width for numeric columns
    W_DELTA = 16 # width for delta columns

    header = (
        f"{'File':<{col_file}}{SEP}"
        f"{'Records':>{W_NUM}}{SEP}"
        f"{'ΔRecords':>{W_DELTA}}{SEP}"
        f"{'Sum of NNNNx':>{W_NUM}}{SEP}"
        f"{'ΔSum':>{W_DELTA}}"
    )
    rule = '-' * len(header)
    print()
    print(header)
    print(rule)

    prev_rec = prev_sum = None
    for display, n_rec, n_sum in rows:
        d_rec = _delta_str(n_rec, prev_rec) if prev_rec is not None else '—'
        d_sum = _delta_str(n_sum, prev_sum) if prev_sum is not None else '—'
        print(
            f"{display:<{col_file}}{SEP}"
            f"{n_rec:>{W_NUM},}{SEP}"
            f"{d_rec:>{W_DELTA}}{SEP}"
            f"{n_sum:>{W_NUM},}{SEP}"
            f"{d_sum:>{W_DELTA}}"
        )
        prev_rec, prev_sum = n_rec, n_sum

    print(rule)

    # ── summary: total lost vs first file ─────────────────────────────────────
    if len(rows) >= 2:
        first_rec, first_sum = rows[0][1], rows[0][2]
        last_rec,  last_sum  = rows[-1][1], rows[-1][2]
        print()
        print(f"Overall change  (last vs first):")
        print(f"  Records : {_delta_str(last_rec, first_rec):>+16,}  ({_pct_str(last_rec, first_rec)} of first)")
        print(f"  Sum     : {_delta_str(last_sum, first_sum):>+16,}  ({_pct_str(last_sum, first_sum)} of first)")


if __name__ == '__main__':
    main()
