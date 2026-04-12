#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# vim: set fileencoding=utf-8 ts=4 sw=4 expandtab :
"""Chunksize sweep for Pool.starmap/imap on calculate_codon_frequencies.

Monkey-patches Pool.starmap to use imap with a given chunksize, then
benchmarks parse_alignment in-process to find the optimal chunksize
for the target machine and workload.

Run on 192-core machine from project root:
    PYTHONPATH=src python scripts/bench_chunksize_sweep.py \
        --reference tests/inputs/MN908947.3_S.3873.fasta \
        --alignment tests/inputs/test2_full.fasta \
        --threads 4 \
        --chunksizes 1 2 4 8 16 32 64 128 256 \
        [--max-rows 50000]

Tip: use --max-rows to benchmark on a manageable subset of a very large
alignment file (e.g. millions of sequences) without changing the original.
"""

import argparse
import multiprocessing
import multiprocessing.pool
import os
import sys
import tempfile
import time
import types

from Bio import SeqIO

from mutation_scatter_plot import alt_translate
from mutation_scatter_plot.calculate_codon_frequencies import (get_codons,
                                                               parse_alignment)

# Allow running from project root without installing
sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(
    os.path.abspath(__file__))), "src"))

# Bio and mutation_scatter_plot must follow sys.path fixup above.
# pylint: disable=wrong-import-position

# ── Subset helper (mirrors bench_thread_scaling.py) ──────────────────────────


def _count_seqs(path):
    """Count FASTA records by counting '>' header lines."""
    n = 0
    with open(path, "rb") as fh:
        for raw in fh:
            if raw.startswith(b">"):
                n += 1
    return n


def _write_subset(src_path, dst_path, n_seqs):
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
                if len(buf) >= 1000:
                    dst.writelines(buf)
                    buf = []
        if buf:
            dst.writelines(buf)


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
    """Parse CLI arguments and run the chunksize sweep benchmark."""
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
    parser.add_argument("--max-rows", type=int, default=0, metavar="N",
                        help="Benchmark on only the first N sequences of the alignment "
                             "(written to a temp file; original is not modified). "
                             "Use 0 (default) to use the full file.")
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

    # Capture current git SHA for output file naming.
    try:
        import subprocess as _sp
        git_sha = _sp.check_output(
            ["git", "rev-parse", "--short", "HEAD"],
            stderr=_sp.DEVNULL, text=True).strip()
    except Exception:
        git_sha = "unknown"

    n_workers = args.threads or multiprocessing.cpu_count()
    aln_size_kb = os.path.getsize(args.alignment) / 1024
    total_cores = multiprocessing.cpu_count()

    # ── Optionally subset the alignment ─────────────────────────────────────
    _subset_tmpdir = None
    alignment_file = args.alignment
    if args.max_rows > 0:
        print(f"# Counting sequences in {os.path.basename(args.alignment)} … ",
              end="", flush=True)
        total_seqs = _count_seqs(args.alignment)
        actual_rows = min(args.max_rows, total_seqs)
        print(f"{total_seqs:,} found")
        if actual_rows < total_seqs:
            _tmpdir = (os.environ.get("TMPDIR")
                       or ("/scratch.ssd/mmokrejs"
                           if os.path.isdir("/scratch.ssd/mmokrejs") else None))
            _subset_tmpdir = tempfile.mkdtemp(dir=_tmpdir)
            alignment_file = os.path.join(
                _subset_tmpdir,
                f"subset_{actual_rows}_{os.path.basename(args.alignment)}")
            print(f"# Writing subset of {actual_rows:,} sequences … ",
                  end="", flush=True)
            _t0 = time.perf_counter()
            _write_subset(args.alignment, alignment_file, actual_rows)
            print(f"done in {time.perf_counter() - _t0:.1f}s")
            aln_size_kb = os.path.getsize(alignment_file) / 1024
        else:
            print(f"# --max-rows {args.max_rows} >= total sequences ({total_seqs:,}); "
                  f"using full file")

    try:
        _run_sweep(args, alignment_file, n_workers, aln_size_kb, total_cores,
                   myoptions, ref_seq, prot_seq, codons, git_sha)
    finally:
        if _subset_tmpdir:
            import shutil
            shutil.rmtree(_subset_tmpdir, ignore_errors=True)


def _run_sweep(args, alignment_file, n_workers, aln_size_kb, total_cores,
               myoptions, ref_seq, prot_seq, codons, git_sha):
    """Inner sweep loop, separated so temp-dir cleanup always runs."""
    print(f"# Chunksize sweep for {os.path.basename(alignment_file)} ({aln_size_kb:.0f} KB)")
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
        def _patched_starmap(self, func, iterable, _chunksize=None, _cs=cs):
            return orig_starmap(self, func, iterable, chunksize=_cs)

        multiprocessing.pool.Pool.starmap = _patched_starmap

        times = []
        for i in range(args.runs):
            ts = time.strftime("%H:%M:%S")
            elapsed = bench_one(ref_seq, prot_seq, codons,
                                alignment_file, n_workers, myoptions)
            times.append(elapsed)
            mn, med, mx = sorted(times)[0], sorted(times)[len(times)//2], sorted(times)[-1]
            print(f"  [{ts}] chunksize={cs_label:>6}  run {i+1}/{args.runs}: "
                  f"{elapsed:.3f}s  (so far: min={mn:.3f} med={med:.3f} max={mx:.3f})",
                  flush=True)

        multiprocessing.pool.Pool.starmap = orig_starmap

        times.sort()
        median = times[len(times) // 2]
        results.append((cs_label, min(times), median, max(times)))
        print(f"{cs_label:>10}  {min(times):>8.3f}  {median:>8.3f}  {max(times):>8.3f}",
              flush=True)

    tsv_path = f"bench_chunksize_{os.path.basename(alignment_file)}_{git_sha}.tsv"
    with open(tsv_path, "w", encoding="utf-8") as f:
        f.write("chunksize\tmin_s\tmedian_s\tmax_s\n")
        for row in results:
            f.write("\t".join(str(x) for x in row) + "\n")
    print(f"\n# TSV written to {tsv_path}")
    best = min((r for r in results if r[0] != 'auto'), key=lambda r: r[2])
    print(f"# Best explicit chunksize (by median): {best[0]} ({best[2]:.3f}s)")


if __name__ == "__main__":
    main()
