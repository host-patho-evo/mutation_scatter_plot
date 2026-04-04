#!/usr/bin/env bash
# vim: set fileencoding=utf-8 ts=4 sw=4 expandtab :
# Line-level profiling for calculate_codon_frequencies using line_profiler.
# Decorates parse_alignment (and optionally _process_one_site) at runtime
# without modifying source, using line_profiler's explicit API.
#
# Prerequisites:  pip install line_profiler
#
# Run from project root:
#   bash scripts/profile_line.sh <alignment.fasta> [threads] [max_rows]
#
# Outputs (all names include the short git commit hash for traceability):
#   $TMPDIR/profile_line_<SHA>_tN.lprof   - binary profile
#   $TMPDIR/profile_line_<SHA>_tN.txt     - human-readable line-by-line report (ms units)
#
# Examples:
#   REFERENCE=.../MN908947.3_S.fasta bash scripts/profile_line.sh \
#       .../spikenuc1207...exactly_3822.fasta 4 100000

set -euo pipefail

ALIGNMENT="${1:-tests/inputs/test2_full.fasta}"
REFERENCE="${REFERENCE:-tests/inputs/MN908947.3_S_full.fasta}"
THREADS="${2:-1}"
MAX_ROWS="${3:-0}"

# Prefer fast local SSD for temp files; fall back to /tmp.
if [[ -z "${TMPDIR:-}" && -d /scratch.ssd/mmokrejs ]]; then
    export TMPDIR=/scratch.ssd/mmokrejs
fi
# Embed git SHA so output files are unambiguously tied to the code revision.
GIT_SHA="$(git rev-parse --short HEAD 2>/dev/null || echo unknown)"
OUTBASE="${TMPDIR:-/tmp}/profile_line_${GIT_SHA}_t${THREADS}"

export PYTHONPATH="${PYTHONPATH:-$(pwd)/src}"

# ── Optional subsetting ───────────────────────────────────────────────────────
SUBSET_FILE=""
if [[ "$MAX_ROWS" -gt 0 ]]; then
    SUBSET_FILE="$(mktemp "${TMPDIR:-/tmp}/profile_subset_XXXXXX.fasta")"
    trap 'rm -f "$SUBSET_FILE"' EXIT
    echo "=== Subsetting: writing first ${MAX_ROWS} sequences to ${SUBSET_FILE} ==="
    python3 - "$ALIGNMENT" "$SUBSET_FILE" "$MAX_ROWS" <<'PYEOF'
import sys
src_path, dst_path, n_seqs = sys.argv[1], sys.argv[2], int(sys.argv[3])
written = 0
with open(src_path, "rb") as src, open(dst_path, "wb") as dst:
    buf = []
    for raw in src:
        if raw.startswith(b">"):
            if written >= n_seqs:
                break
            written += 1
            buf.append(raw)
        else:
            buf.append(raw)
        if len(buf) >= 1000:
            dst.writelines(buf)
            buf = []
    if buf:
        dst.writelines(buf)
print(f"Written {written} sequences.")
PYEOF
    ALIGNMENT="$SUBSET_FILE"
fi

echo "=== line_profiler: calculate_codon_frequencies ==="
echo "  alignment : $ALIGNMENT  ($(du -h "$ALIGNMENT" | cut -f1))"
echo "  reference : $REFERENCE"
echo "  threads   : $THREADS"
echo "  max_rows  : ${MAX_ROWS:-all}"
echo "  lprof out : ${OUTBASE}.lprof"
echo ""

# Run via line_profiler's explicit API (no @profile decorator needed in source).
python3 - \
    "$REFERENCE" "$ALIGNMENT" "$OUTBASE" "$THREADS" \
    <<'PYEOF'
import sys, line_profiler

reference_file  = sys.argv[1]
alignment_file  = sys.argv[2]
outbase         = sys.argv[3]
threads         = int(sys.argv[4]) if int(sys.argv[4]) > 0 else None

import mutation_scatter_plot.calculate_codon_frequencies as _mod

lp = line_profiler.LineProfiler()
# Profile the two most interesting functions.
lp.add_function(_mod.parse_alignment)
lp.add_function(_mod._process_one_site)
lp.add_function(_mod.fast_fasta_iter)

# Minimal CLI setup to call parse_alignment directly.
import argparse, io
from mutation_scatter_plot.calculate_codon_frequencies import open_file
from mutation_scatter_plot.calculate_codon_frequencies.cli import build_option_parser
from Bio import SeqIO

myoptions = build_option_parser().parse_args([
    "--reference-infile", reference_file,
    "--alignment-file",   alignment_file,
    "--outfile-prefix",   outbase + "_out",
    "--padded-reference",
    "--x-after-count",
    "--overwrite",
    f"--threads={threads or 0}",
    f"--chunksize=256",
])

# Read reference (mirrors cli.main logic).
_ref_records = list(SeqIO.parse(myoptions.reference_infilename, "fasta"))
_padded_ref_dna   = str(_ref_records[0].seq)
_ref_protein      = str(_ref_records[0].seq.translate(table=1))[:-1]
from mutation_scatter_plot.calculate_codon_frequencies import get_codons
_ref_as_codons    = get_codons(_padded_ref_dna)

import os
_outfile  = open_file(outbase + "_out.tsv",             overwrite=True, encoding="utf-8")
_unchanged = open_file(outbase + "_out.unchanged.tsv",  overwrite=True, encoding="utf-8")
_count_fh  = open_file(outbase + "_out.count.txt",      overwrite=True, encoding="utf-8")

_aa_start  = (myoptions.aa_start  - 1) if myoptions.aa_start  else 0
_min_start = (myoptions.min_start - 1) if myoptions.min_start else 0
_max_stop  = (myoptions.max_stop  + 1) if myoptions.max_stop  else 0
_threads   = myoptions.threads if myoptions.threads > 0 else None
_chunksize = myoptions.chunksize if myoptions.chunksize > 0 else None

# Wrap parse_alignment with line profiler.
wrapped = lp(_mod.parse_alignment)
wrapped(
    myoptions, alignment_file,
    _padded_ref_dna, _ref_protein, _ref_as_codons,
    _outfile, _unchanged, _count_fh,
    _aa_start, _min_start, _max_stop,
    threads=_threads, chunksize=_chunksize,
)

_outfile.close(); _unchanged.close(); _count_fh.close()

lp.dump_stats(outbase + ".lprof")
print(f"\nProfile saved: {outbase}.lprof")

# Print the report.
with open(outbase + ".txt", "w") as f:
    lp.print_stats(stream=f, output_unit=1e-3)  # ms units
lp.print_stats(output_unit=1e-3)
PYEOF

echo ""
echo "Done. Retrieve files:"
echo "  ${OUTBASE}.txt     (human-readable line timings, ms)"
echo "  ${OUTBASE}.lprof   (binary, view with: python -m line_profiler ${OUTBASE}.lprof)"
