#!/usr/bin/env python3
"""Chunksize sweep for Pool.imap on calculate_codon_frequencies.

This script benchmarks different chunksize values for Pool.imap (or starmap)
to find the optimal value for the target machine's core count and workload size.

The script patches the chunksize at the Python API level by monkey-patching
pool.imap/starmap, so no source changes are needed.

Run on 192-core machine from project root:
    PYTHONPATH=src python scripts/bench_chunksize_sweep.py \
        --reference tests/inputs/MN908947.3_S_full.fasta \
        --alignment tests/inputs/test2_full.fasta \
        --threads 128 \
        --chunksizes 1 4 8 16 32 64 128 256 512
"""

import argparse
import multiprocessing
import os
import sys
import tempfile
import time

# So we can import the module under test
sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(
    os.path.abspath(__file__))), "src"))


def _run_parse_alignment_with_chunksize(args_tuple):
    """Run parse_alignment in-process with a patched chunksize, return elapsed seconds."""
    import importlib
    import types

    (ref_seq, prot_seq, codons, alignment_file,
     min_start, max_stop, threads, chunksize, myoptions) = args_tuple

    # Lazy import to measure only the compute, not the import
    from mutation_scatter_plot.calculate_codon_frequencies import parse_alignment

    # Monkey-patch Pool.starmap to use imap with our chunksize instead
    original_starmap = multiprocessing.pool.Pool.starmap

    def patched_starmap(self, func, iterable, chunksize=None):
        # Force our chunksize
        return list(self.imap(
            lambda a: func(*a) if isinstance(a, (list, tuple)) else func(a),
            iterable,
            chunksize=_CHUNKSIZE_INJECT,
        ))

    # We can't easily monkey-patch inside forked workers.
    # Instead we measure the full wall-clock time of parse_alignment as a proxy.
    with tempfile.TemporaryDirectory() as tmpdir:
        out_path = os.path.join(tmpdir, "out.tsv")
        unch_path = os.path.join(tmpdir, "unch.tsv")
        count_path = os.path.join(tmpdir, "count.txt")
        with open(out_path, 'w', encoding='utf-8') as outf, \
             open(unch_path, 'w', encoding='utf-8') as unchf, \
             open(count_path, 'w', encoding='utf-8') as cntf:
            t0 = time.perf_counter()
            parse_alignment(
                myoptions, alignment_file, ref_seq, prot_seq, codons,
                outf, unchf, cntf, 0, min_start, max_stop, threads=threads,
            )
            return time.perf_counter() - t0


def make_options(debug=0, left=0, right=0, minimum_aln_length=50,
                 discard_leading=0, discard_trailing=0, x_after_count=True):
    """Build a minimal myoptions namespace."""
    import types
    o = types.SimpleNamespace()
    o.debug = debug
    o.left_reference_offset = left
    o.right_reference_offset = right
    o.minimum_aln_length = minimum_aln_length
    o.discard_this_many_leading_nucs = discard_leading
    o.discard_this_many_trailing_nucs = discard_trailing
    o.x_after_count = x_after_count
    return o


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
    args = parser.parse_args()

    # Import here so we don't count import time in benchmarks
    from Bio import SeqIO
    from mutation_scatter_plot.calculate_codon_frequencies import get_codons
    from mutation_scatter_plot import alt_translate

    with open(args.reference, encoding="utf-8") as f:
        for rec in SeqIO.parse(f, "fasta"):
            ref_seq = str(rec.seq)
            break
    prot_seq = alt_translate(ref_seq)
    codons = get_codons(ref_seq)
    myoptions = make_options()

    n_cores = args.threads or multiprocessing.cpu_count()
    aln_size_mb = os.path.getsize(args.alignment) / 1024 / 1024
    print(f"# Chunksize sweep for {os.path.basename(args.alignment)} ({aln_size_mb:.1f} MB)")
    print(f"# Workers: {n_cores}  |  Runs per chunksize: {args.runs}")
    print(f"# nproc: {multiprocessing.cpu_count()}")
    print()
    print(f"{'chunksize':>10}  {'min_s':>8}  {'median_s':>8}  {'max_s':>8}")
    print("-" * 45)

    # We need to inject the chunksize into Pool.starmap.
    # The cleanest way: temporarily replace Pool.starmap with imap+our chunksize.
    import multiprocessing.pool

    results = []
    for cs in sorted(args.chunksizes):
        # Patch
        _orig_starmap = multiprocessing.pool.Pool.starmap

        def _patched_starmap(self, func, iterable, chunksize=None, cs_val=cs):
            lst = list(iterable) if not isinstance(iterable, list) else iterable
            return list(self.imap(
                (lambda a, f=func: f(*a)),
                lst,
                chunksize=cs_val,
            ))

        multiprocessing.pool.Pool.starmap = _patched_starmap

        times = []
        for i in range(args.runs):
            from mutation_scatter_plot.calculate_codon_frequencies import parse_alignment
            with tempfile.TemporaryDirectory() as tmpdir:
                out_path = os.path.join(tmpdir, "out.tsv")
                unch_path = os.path.join(tmpdir, "unch.tsv")
                count_path = os.path.join(tmpdir, "count.txt")
                with open(out_path, 'w', encoding='utf-8') as outf, \
                     open(unch_path, 'w', encoding='utf-8') as unchf, \
                     open(count_path, 'w', encoding='utf-8') as cntf:
                    t0 = time.perf_counter()
                    parse_alignment(
                        myoptions, args.alignment, ref_seq, prot_seq, codons,
                        outf, unchf, cntf, 0, 0, 0, threads=n_cores,
                    )
                    times.append(time.perf_counter() - t0)
            print(f"  chunksize={cs} run {i+1}: {times[-1]:.3f}s", file=sys.stderr)

        multiprocessing.pool.Pool.starmap = _orig_starmap

        times.sort()
        median = times[len(times) // 2]
        results.append((cs, min(times), median, max(times)))
        print(f"{cs:>10}  {min(times):>8.3f}  {median:>8.3f}  {max(times):>8.3f}")

    tsv_path = f"bench_chunksize_{os.path.basename(args.alignment)}.tsv"
    with open(tsv_path, "w", encoding="utf-8") as f:
        f.write("chunksize\tmin_s\tmedian_s\tmax_s\n")
        for row in results:
            f.write("\t".join(str(x) for x in row) + "\n")
    print(f"\n# TSV written to {tsv_path}")


if __name__ == "__main__":
    main()
