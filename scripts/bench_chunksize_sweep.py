#!/usr/bin/env python3
"""Chunksize sweep for Pool.starmap/imap on calculate_codon_frequencies.

Monkey-patches Pool.starmap to use imap with a given chunksize, then
benchmarks parse_alignment in-process to find the optimal chunksize
for the target machine and workload.

Run on 192-core machine from project root:
    PYTHONPATH=src python scripts/bench_chunksize_sweep.py \
        --reference tests/inputs/MN908947.3_S_full.fasta \
        --alignment tests/inputs/test2_full.fasta \
        --threads 4 \
        --chunksizes 1 2 4 8 16 32 64 128 256
"""

import argparse
import multiprocessing
import multiprocessing.pool
import os
import sys
import tempfile
import time
import types

# Allow running from project root without installing
sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(
    os.path.abspath(__file__))), "src"))

from Bio import SeqIO  # noqa: E402  (after sys.path fixup)
from mutation_scatter_plot.calculate_codon_frequencies import (  # noqa: E402
    get_codons, parse_alignment,
)
from mutation_scatter_plot import alt_translate  # noqa: E402


def make_options(x_after_count=True, minimum_aln_length=50,
                 left=0, right=0, discard_leading=0, discard_trailing=0, debug=0):
    """Build a minimal myoptions namespace matching what cli.py passes."""
    o = types.SimpleNamespace()
    o.debug = debug
    o.left_reference_offset = left
    o.right_reference_offset = right
    o.minimum_aln_length = minimum_aln_length
    o.discard_this_many_leading_nucs = discard_leading
    o.discard_this_many_trailing_nucs = discard_trailing
    o.x_after_count = x_after_count
    return o


def bench_one(ref_seq, prot_seq, codons, alignment_file, n_workers, myoptions):
    """Run parse_alignment once and return the elapsed wall-clock seconds."""
    with tempfile.TemporaryDirectory() as tmpdir:
        with (open(os.path.join(tmpdir, "out.tsv"), 'w', encoding='utf-8') as outf,
              open(os.path.join(tmpdir, "unch.tsv"), 'w', encoding='utf-8') as unchf,
              open(os.path.join(tmpdir, "count.txt"), 'w', encoding='utf-8') as cntf):
            t0 = time.perf_counter()
            parse_alignment(
                myoptions, alignment_file, ref_seq, prot_seq, codons,
                outf, unchf, cntf, 0, 0, 0, threads=n_workers,
            )
            return time.perf_counter() - t0


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--reference", required=True)
    parser.add_argument("--alignment", required=True)
    parser.add_argument("--threads", type=int, default=None,
                        help="Worker count (default: all available cores)")
    parser.add_argument("--chunksizes", nargs="+", type=int,
                        default=[1, 4, 8, 16, 32, 64, 128, 256, 512],
                        help="chunksize values to sweep")
    parser.add_argument("--runs", type=int, default=3)
    parser.add_argument("--x-after-count", action="store_true", default=True,
                        help="Pass --x-after-count to parse_alignment")
    args = parser.parse_args()

    with open(args.reference, encoding="utf-8") as f:
        for rec in SeqIO.parse(f, "fasta"):
            ref_seq = str(rec.seq)
            break
    prot_seq = alt_translate(ref_seq)
    codons = get_codons(ref_seq)
    myoptions = make_options(x_after_count=args.x_after_count)

    n_workers = args.threads or multiprocessing.cpu_count()
    aln_size_kb = os.path.getsize(args.alignment) / 1024
    total_cores = multiprocessing.cpu_count()

    print(f"# Chunksize sweep for {os.path.basename(args.alignment)} ({aln_size_kb:.0f} KB)")
    print(f"# Workers: {n_workers}  |  Total cores: {total_cores}  |  Runs: {args.runs}")
    print()
    print(f"{'chunksize':>10}  {'min_s':>8}  {'median_s':>8}  {'max_s':>8}")
    print("-" * 45)

    orig_starmap = multiprocessing.pool.Pool.starmap

    # Include auto chunksize (None) as the baseline
    sweep = [(None, 'auto')] + [(cs, str(cs)) for cs in sorted(args.chunksizes)]

    results = []
    for cs, cs_label in sweep:
        # Override chunksize by patching Pool.starmap to pass our value
        def _patched_starmap(self, func, iterable, chunksize=None, _cs=cs):
            return orig_starmap(self, func, iterable, chunksize=_cs)

        multiprocessing.pool.Pool.starmap = _patched_starmap

        times = []
        for i in range(args.runs):
            elapsed = bench_one(ref_seq, prot_seq, codons,
                                args.alignment, n_workers, myoptions)
            times.append(elapsed)
            print(f"  chunksize={cs_label} run {i+1}: {elapsed:.3f}s", file=sys.stderr)

        multiprocessing.pool.Pool.starmap = orig_starmap

        times.sort()
        median = times[len(times) // 2]
        results.append((cs_label, min(times), median, max(times)))
        print(f"{cs_label:>10}  {min(times):>8.3f}  {median:>8.3f}  {max(times):>8.3f}")

    tsv_path = f"bench_chunksize_{os.path.basename(args.alignment)}.tsv"
    with open(tsv_path, "w", encoding="utf-8") as f:
        f.write("chunksize\tmin_s\tmedian_s\tmax_s\n")
        for row in results:
            f.write("\t".join(str(x) for x in row) + "\n")
    print(f"\n# TSV written to {tsv_path}")
    best = min((r for r in results if r[0] != 'auto'), key=lambda r: r[2])
    print(f"# Best explicit chunksize (by median): {best[0]} ({best[2]:.3f}s)")


if __name__ == "__main__":
    main()
