#!/usr/bin/env python3
"""Thread-scaling benchmark for calculate_codon_frequencies.

Run this on the target machine (192-core Xeon) to measure actual wall-clock
scaling across thread counts.  Prints a live result line after EVERY timed
run so you can monitor progress even when each run takes minutes.

Usage (from project root):
    PYTHONPATH=src python scripts/bench_thread_scaling.py \\
        --reference tests/inputs/MN908947.3_S_full.fasta \\
        --alignment tests/inputs/test2_full.fasta \\
        [--threads 1 2 4 8 16 32 64 128 192] \\
        [--runs 3]

For a quick single-pass scan use --runs 1.
"""

import argparse
import os
import subprocess
import sys
import tempfile
import time


def run_one(ref, aln, threads, extra_args, env):
    with tempfile.TemporaryDirectory() as tmpdir:
        prefix = os.path.join(tmpdir, "bench")
        cmd = [
            sys.executable, "-m",
            "mutation_scatter_plot.calculate_codon_frequencies.cli",
            "--reference-infile", ref,
            "--alignment-file", aln,
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
            print(f"  FAILED (threads={threads}):\n{r.stderr.decode()[:500]}", file=sys.stderr)
            return None
        return elapsed


def _stats(times):
    """Return (min, median, max) from an unsorted list."""
    s = sorted(times)
    return s[0], s[len(s) // 2], s[-1]


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--reference", required=True, help="Reference FASTA file")
    parser.add_argument("--alignment", required=True, help="Alignment FASTA file")
    parser.add_argument("--threads", nargs="+", type=int,
                        default=[1, 2, 4, 8, 16, 32, 64, 128],
                        help="Thread counts to benchmark")
    parser.add_argument("--runs", type=int, default=3,
                        help="Timed runs per thread count (uses median).  "
                             "Use --runs 1 for a quick single-pass scan.")
    parser.add_argument("--aa-start", type=int, default=0, help="--aa_start value")
    parser.add_argument("--extra", nargs="*", default=[], help="Extra CLI args to pass through")
    args = parser.parse_args()

    env = os.environ.copy()
    if "PYTHONPATH" not in env:
        env["PYTHONPATH"] = os.path.join(
            os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "src")

    extra = args.extra or []
    if args.aa_start:
        extra += [f"--aa_start={args.aa_start}"]

    aln_size_mb = os.path.getsize(args.alignment) / 1024 / 1024
    ncpu = os.popen("nproc").read().strip()
    print(f"# Benchmark : {os.path.basename(args.alignment)} ({aln_size_mb:.1f} MB)")
    print(f"# Reference : {os.path.basename(args.reference)}")
    print(f"# Runs/count: {args.runs}")
    print(f"# CPU cores : {ncpu} logical")
    print()

    # Live header — printed once and flushed
    hdr = (f"{'threads':>8}  {'run':>3}  {'this_s':>8}  "
           f"{'so_far_min':>10}  {'so_far_med':>10}  {'so_far_max':>10}  {'speedup':>8}")
    print(hdr)
    print("-" * len(hdr))
    sys.stdout.flush()

    baseline_median = None
    results = []

    for t in sorted(args.threads):
        times = []
        for i in range(args.runs):
            t_start_wall = time.strftime("%H:%M:%S")
            elapsed = run_one(args.reference, args.alignment, t, extra, env)
            if elapsed is None:
                break
            times.append(elapsed)
            mn, med, mx = _stats(times)
            sp_str = "--"
            if baseline_median is not None and med > 0:
                sp_str = f"{baseline_median / med:.2f}x"
            # One live line per run — always flushed immediately
            print(f"{t:>8}  {i+1:>3}  {elapsed:>8.2f}  "
                  f"{mn:>10.2f}  {med:>10.2f}  {mx:>10.2f}  {sp_str:>8}"
                  f"  [{t_start_wall}]")
            sys.stdout.flush()

        if not times:
            print(f"{t:>8}  ---  FAILED")
            sys.stdout.flush()
            continue

        mn, median, mx = _stats(times)
        if baseline_median is None:
            baseline_median = median
        speedup = baseline_median / median if median > 0 else float("inf")
        results.append((t, mn, median, mx, speedup))
        print(f"#   t={t:>3} DONE  min={mn:.2f}s  median={median:.2f}s  "
              f"max={mx:.2f}s  speedup={speedup:.2f}x")
        print()
        sys.stdout.flush()

    # ── Final summary table ──────────────────────────────────────────────────
    print()
    print(f"{'threads':>8}  {'min_s':>8}  {'median_s':>8}  {'max_s':>8}  {'speedup':>8}")
    print("-" * 50)
    for t, mn, med, mx, sp in results:
        print(f"{t:>8}  {mn:>8.2f}  {med:>8.2f}  {mx:>8.2f}  {sp:>7.2f}x")
    sys.stdout.flush()

    tsv_path = f"bench_scaling_{os.path.basename(args.alignment)}.tsv"
    with open(tsv_path, "w", encoding="utf-8") as f:
        f.write("threads\tmin_s\tmedian_s\tmax_s\tspeedup\n")
        for row in results:
            f.write("\t".join(str(x) for x in row) + "\n")
    print(f"# TSV written to {tsv_path}")


if __name__ == "__main__":
    main()
