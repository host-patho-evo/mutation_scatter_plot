#!/usr/bin/env python3
"""Thread-scaling benchmark for calculate_codon_frequencies.

Run this on the target machine (192-core Xeon) to measure actual wall-clock
scaling across thread counts.  Prints a live result line after EVERY timed run.

Key features
------------
* Inspects the alignment at startup: sequence count and length.
* --row-sizes lets you benchmark on escalating subsets (100 → 1000 → all).
  Subsets are written to a temp file — no modification of the original.
* --min-col / --max-col limit to a codon range (natural 1-based codon numbers).
  Example: --min-col 4 --max-col 999  tests codons 4–999 only.
* Live per-run output with running min/median/max and speedup.

Usage (from project root):
    PYTHONPATH=src python scripts/bench_thread_scaling.py \\
        --reference tests/inputs/MN908947.3_S_full.fasta \\
        --alignment tests/inputs/test2_full.fasta \\
        [--threads 1 4 8 16 32 64 128 192] \\
        [--runs 3] \\
        [--row-sizes 500 5000 50000] \\
        [--min-col 1] [--max-col 999]

Tip: use --runs 1 for a quick single-pass scan.
"""

import argparse
import os
import subprocess
import sys
import tempfile
import time


# ── File inspection ──────────────────────────────────────────────────────────

def inspect_fasta(path):
    """Return (n_sequences, seq_len) by streaming through the file once.

    Counts header lines ('>') for sequence count.
    Reads only the very first sequence line to obtain length (all sequences in
    a padded alignment have the same length, so this is sufficient).
    """
    n_seqs = 0
    seq_len = None
    with open(path, "rb") as fh:
        for raw in fh:
            if raw.startswith(b">"):
                n_seqs += 1
            elif seq_len is None:
                seq_len = len(raw.rstrip())
    return n_seqs, seq_len


def write_subset(src_path, dst_path, n_seqs):
    """Write the first *n_seqs* sequences from src_path into dst_path."""
    written = 0
    with open(src_path, "rb") as src, open(dst_path, "wb") as dst:
        buf = []
        in_seq = False
        for raw in src:
            if raw.startswith(b">"):
                if written >= n_seqs:
                    break
                written += 1
                in_seq = True
                buf.append(raw)
            elif in_seq:
                buf.append(raw)
                # Flush every 1000 lines to avoid huge in-memory buffers
                if len(buf) >= 1000:
                    dst.writelines(buf)
                    buf = []
        if buf:
            dst.writelines(buf)


# ── Benchmark helpers ────────────────────────────────────────────────────────

def build_auto_row_sizes(total_seqs):
    """Geometric series from 100 up to total_seqs, always including the last."""
    sizes = []
    n = 100
    while n < total_seqs:
        sizes.append(n)
        n = int(n * 3.16)  # ≈ one order of magnitude per 2 steps
    sizes.append(total_seqs)
    return sizes


def run_one(ref, aln_path, threads, extra_args, env):
    with tempfile.TemporaryDirectory() as tmpdir:
        prefix = os.path.join(tmpdir, "bench")
        cmd = [
            sys.executable, "-m",
            "mutation_scatter_plot.calculate_codon_frequencies.cli",
            "--reference-infile", ref,
            "--alignment-file", aln_path,
            "--outfile-prefix", prefix,
            "--padded-reference",
            "--x-after-count",
            "--overwrite",
            f"--threads={threads}",
        ] + extra_args
        t0 = time.perf_counter()
        r = subprocess.run(cmd, env=env, capture_output=True, check=False)
        elapsed = time.perf_counter() - t0
        if r.returncode != 0:
            print(f"  FAILED threads={threads}:\n{r.stderr.decode()[:500]}", file=sys.stderr)
            return None
        return elapsed


def _stats(times):
    s = sorted(times)
    return s[0], s[len(s) // 2], s[-1]


# ── Main ─────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--reference", required=True, help="Reference FASTA file")
    parser.add_argument("--alignment", required=True, help="Alignment FASTA file")
    parser.add_argument("--threads", nargs="+", type=int,
                        default=[1, 2, 4, 8, 16, 32, 64, 128],
                        help="Thread counts to benchmark")
    parser.add_argument("--runs", type=int, default=3,
                        help="Timed runs per (thread, row-size) combination.  "
                             "Use --runs 1 for a quick single-pass scan.")
    parser.add_argument("--row-sizes", nargs="+", type=int, default=None,
                        metavar="N",
                        help="Run on subsets of N sequences.  "
                             "Omit to use --row-auto or the full file.")
    parser.add_argument("--row-auto", action="store_true",
                        help="Auto-generate a geometric row-size ladder up to --max-rows "
                             "(or total sequences if --max-rows 0).")
    parser.add_argument("--max-rows", type=int, default=100_000,
                        metavar="N",
                        help="Cap for --row-auto ladder (default: 100000).  "
                             "Use 0 to go all the way to total sequences.")
    parser.add_argument("--min-col", type=int, default=1, metavar="CODON",
                        help="First codon to include (1-based natural numbering, default=1).")
    parser.add_argument("--max-col", type=int, default=0, metavar="CODON",
                        help="Last codon to include (1-based natural numbering, default=all).  "
                             "Examples: --max-col 9  --max-col 999")
    parser.add_argument("--aa-start", type=int, default=0, help="--aa_start value")
    parser.add_argument("--extra", nargs="*", default=[], help="Extra CLI args to pass through")
    args = parser.parse_args()

    env = os.environ.copy()
    if "PYTHONPATH" not in env:
        env["PYTHONPATH"] = os.path.join(
            os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "src")

    # Capture current git SHA for output file naming and header display.
    try:
        import subprocess as _sp
        git_sha = _sp.check_output(
            ["git", "rev-parse", "--short", "HEAD"],
            stderr=_sp.DEVNULL, text=True).strip()
    except Exception:
        git_sha = "unknown"

    extra = list(args.extra or [])
    if args.aa_start:
        extra += [f"--aa_start={args.aa_start}"]
    # Convert natural codon numbering → 1-based nucleotide positions for the CLI.
    # CLI internally does: _min_start = user_val - 1  (0-indexed)
    #                      _max_stop  = user_val + 1  (exclusive)
    # For codon C (1-based): nt start = (C-1)*3+1, exclusive bound uses (C-1)*3.
    if args.min_col > 1:
        extra += [f"--min_start={(args.min_col - 1) * 3 + 1}"]
    if args.max_col > 0:
        extra += [f"--max_stop={(args.max_col - 1) * 3}"]

    # ── Inspect the alignment file ───────────────────────────────────────────
    print("# Inspecting alignment file … ", end="", flush=True)
    t_inspect = time.perf_counter()
    total_seqs, seq_len = inspect_fasta(args.alignment)
    inspect_s = time.perf_counter() - t_inspect
    aln_size_mb = os.path.getsize(args.alignment) / 1024 / 1024
    ncpu = os.popen("nproc").read().strip()

    print(f"done in {inspect_s:.1f}s")
    print(f"# Alignment : {os.path.basename(args.alignment)} ({aln_size_mb:.1f} MB)")
    print(f"#   sequences : {total_seqs:,}")
    print(f"#   seq length : {seq_len:,} nt  ({seq_len // 3:,} codons)")
    if args.max_col > 0:
        col_range = f"codons {args.min_col}–{args.max_col}"
    elif args.min_col > 1:
        col_range = f"codons {args.min_col}–end"
    else:
        col_range = f"all {seq_len // 3:,} codons"
    print(f"#   col range  : {col_range}")
    print(f"# Reference : {os.path.basename(args.reference)}")
    print(f"# CPU cores : {ncpu} logical")
    print(f"# Git SHA   : {git_sha}")
    print(f"# Runs/combo: {args.runs}")

    # ── Determine row-size ladder ────────────────────────────────────────────
    if args.row_sizes:
        row_sizes = [min(n, total_seqs) for n in sorted(set(args.row_sizes))]
    elif args.row_auto:
        cap = total_seqs if args.max_rows == 0 else min(args.max_rows, total_seqs)
        row_sizes = build_auto_row_sizes(cap)
        if cap < total_seqs:
            print(f"# (--row-auto capped at {cap:,}; use --max-rows 0 to go up to {total_seqs:,})")
    else:
        row_sizes = [total_seqs]  # full file only

    print(f"# Row sizes : {row_sizes}")
    print()

    # ── Column header ────────────────────────────────────────────────────────
    hdr = (f"{'rows':>10}  {'threads':>7}  {'run':>3}  "
           f"{'this_s':>7}  {'sf_min':>7}  {'sf_med':>7}  {'sf_max':>7}  "
           f"{'speedup':>8}  time")
    print(hdr)
    print("-" * len(hdr))
    sys.stdout.flush()

    all_results = []  # (row_size, threads, mn, med, mx, speedup)

    with tempfile.TemporaryDirectory(dir=os.environ.get("TMPDIR")) as subset_dir:
        for row_size in row_sizes:
            # Write the subset FASTA once per row_size (reused across all thread counts)
            if row_size < total_seqs:
                subset_path = os.path.join(subset_dir, f"subset_{row_size}.fasta")
                if not os.path.exists(subset_path):
                    print(f"# Creating subset: {row_size:,} sequences … ",
                          end="", flush=True)
                    t0 = time.perf_counter()
                    write_subset(args.alignment, subset_path, row_size)
                    print(f"done in {time.perf_counter() - t0:.1f}s")
            else:
                subset_path = args.alignment

            baseline_median = None
            for t in sorted(args.threads):
                times = []
                for i in range(args.runs):
                    ts = time.strftime("%H:%M:%S")
                    elapsed = run_one(args.reference, subset_path, t, extra, env)
                    if elapsed is None:
                        break
                    times.append(elapsed)
                    mn, med, mx = _stats(times)
                    sp_str = "--"
                    if baseline_median is not None and med > 0:
                        sp_str = f"{baseline_median / med:.2f}x"
                    print(f"{row_size:>10,}  {t:>7}  {i+1:>3}  "
                          f"{elapsed:>7.2f}  {mn:>7.2f}  {med:>7.2f}  {mx:>7.2f}  "
                          f"{sp_str:>8}  {ts}")
                    sys.stdout.flush()

                if not times:
                    print(f"{row_size:>10,}  {t:>7}  ---  FAILED")
                    sys.stdout.flush()
                    continue

                mn, median, mx = _stats(times)
                if baseline_median is None:
                    baseline_median = median
                speedup = baseline_median / median if median > 0 else float("inf")
                all_results.append((row_size, t, mn, median, mx, speedup))

            print(f"# rows={row_size:,} done")
            print()
            sys.stdout.flush()

    # ── Final summary table ──────────────────────────────────────────────────
    print()
    print(f"{'rows':>10}  {'threads':>7}  {'min_s':>7}  {'median_s':>8}  "
          f"{'max_s':>7}  {'speedup':>8}")
    print("-" * 58)
    for row_size, t, mn, med, mx, sp in all_results:
        print(f"{row_size:>10,}  {t:>7}  {mn:>7.2f}  {med:>8.2f}  {mx:>7.2f}  {sp:>7.2f}x")
    sys.stdout.flush()

    tsv_path = f"bench_scaling_{os.path.basename(args.alignment)}_{git_sha}.tsv"
    with open(tsv_path, "w", encoding="utf-8") as f:
        f.write("rows\tthreads\tmin_s\tmedian_s\tmax_s\tspeedup\n")
        for row in all_results:
            f.write("\t".join(str(x) for x in row) + "\n")
    print(f"# TSV written to {tsv_path}")


if __name__ == "__main__":
    main()
