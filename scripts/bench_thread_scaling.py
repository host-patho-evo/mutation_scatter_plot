#!/usr/bin/env python3
"""Thread-scaling benchmark for calculate_codon_frequencies.

Run this on the target machine (192-core Xeon) to measure actual wall-clock
scaling across thread counts. Outputs a TSV-formatted table.

Usage (from project root):
    PYTHONPATH=src python scripts/bench_thread_scaling.py \
        --reference tests/inputs/MN908947.3_S_full.fasta \
        --alignment tests/inputs/test2_full.fasta \
        [--threads 1 2 4 8 16 32 64 128 192] \
        [--runs 5]

Or with a large custom file:
    PYTHONPATH=src python scripts/bench_thread_scaling.py \
        --reference /path/to/reference.fasta \
        --alignment /path/to/large_alignment.fasta \
        --threads 1 2 4 8 16 32 64 128 192 \
        --runs 3
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
            print(f"  FAILED (threads={threads}):\n{r.stderr.decode()[:400]}", file=sys.stderr)
            return None
        return elapsed


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--reference", required=True, help="Reference FASTA file")
    parser.add_argument("--alignment", required=True, help="Alignment FASTA file")
    parser.add_argument("--threads", nargs="+", type=int,
                        default=[1, 2, 4, 8, 16, 32, 64, 128],
                        help="Thread counts to benchmark")
    parser.add_argument("--runs", type=int, default=3,
                        help="Number of timed runs per thread count (uses median)")
    parser.add_argument("--aa-start", type=int, default=0, help="--aa_start value")
    parser.add_argument("--extra", nargs="*", default=[], help="Extra CLI args to pass through")
    args = parser.parse_args()

    env = os.environ.copy()
    if "PYTHONPATH" not in env:
        env["PYTHONPATH"] = os.path.join(os.path.dirname(os.path.dirname(
            os.path.abspath(__file__))), "src")

    extra = args.extra or []
    if args.aa_start:
        extra += [f"--aa_start={args.aa_start}"]

    aln_size_mb = os.path.getsize(args.alignment) / 1024 / 1024
    print(f"# Benchmark: {os.path.basename(args.alignment)} ({aln_size_mb:.1f} MB)")
    print(f"# Reference: {os.path.basename(args.reference)}")
    print(f"# Runs per thread count: {args.runs}")
    print(f"# CPU info: {os.popen('nproc').read().strip()} logical cores available")
    print()
    print(f"{'threads':>8}  {'min_s':>8}  {'median_s':>8}  {'max_s':>8}  {'speedup':>8}")
    print("-" * 50)

    baseline_median = None
    results = []
    for t in sorted(args.threads):
        times = []
        for i in range(args.runs):
            elapsed = run_one(args.reference, args.alignment, t, extra, env)
            if elapsed is None:
                break
            times.append(elapsed)
            print(f"  threads={t} run {i+1}: {elapsed:.3f}s", file=sys.stderr)

        if not times:
            print(f"{t:>8}  {'FAILED':>8}")
            continue

        times.sort()
        median = times[len(times) // 2]
        if baseline_median is None:
            baseline_median = median
        speedup = baseline_median / median if median > 0 else float("inf")
        results.append((t, min(times), median, max(times), speedup))
        print(f"{t:>8}  {min(times):>8.3f}  {median:>8.3f}  {max(times):>8.3f}  {speedup:>7.2f}x")

    print()
    print("# Saved results above. Copy this table into docs/profiling_report.md")

    # Also write a TSV for easy import
    tsv_path = f"bench_scaling_{os.path.basename(args.alignment)}.tsv"
    with open(tsv_path, "w", encoding="utf-8") as f:
        f.write("threads\tmin_s\tmedian_s\tmax_s\tspeedup\n")
        for row in results:
            f.write("\t".join(str(x) for x in row) + "\n")
    print(f"# TSV written to {tsv_path}")


if __name__ == "__main__":
    main()
