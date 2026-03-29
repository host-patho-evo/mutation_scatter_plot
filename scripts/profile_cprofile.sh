#!/usr/bin/env bash
# cProfile run for calculate_codon_frequencies on the target machine.
# Run from project root:  bash scripts/profile_cprofile.sh <alignment.fasta> [threads] [max_rows]
#
# Outputs:
#   /tmp/profile_codon_tN.stats              - cProfile binary (open with pstats or snakeviz)
#   /tmp/profile_codon_tN_top50_*.txt        - pre-rendered top-50 hotspots
#
# Usage examples:
#   bash scripts/profile_cprofile.sh tests/inputs/test2_full.fasta
#   bash scripts/profile_cprofile.sh /data/large_alignment.fasta 1
#   bash scripts/profile_cprofile.sh /data/large_alignment.fasta 64
#   bash scripts/profile_cprofile.sh /data/large_alignment.fasta 4 100000
#                                                                    ^-- profile on first 100k seqs only

set -euo pipefail

ALIGNMENT="${1:-tests/inputs/test2_full.fasta}"
REFERENCE="${REFERENCE:-tests/inputs/MN908947.3_S_full.fasta}"
THREADS="${2:-1}"
MAX_ROWS="${3:-0}"

# Prefer fast local SSD for temp files; fall back to /tmp.
if [[ -z "${TMPDIR:-}" && -d /scratch.ssd/mmokrejs ]]; then
    export TMPDIR=/scratch.ssd/mmokrejs
fi
OUTBASE="${TMPDIR:-/tmp}/profile_codon_t${THREADS}"

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
print(f"Written {written} sequences.")
PYEOF
    ALIGNMENT="$SUBSET_FILE"
fi

echo "=== cProfile: calculate_codon_frequencies ==="
echo "  alignment : $ALIGNMENT  ($(du -h "$ALIGNMENT" | cut -f1))"
echo "  reference : $REFERENCE"
echo "  threads   : $THREADS"
echo "  max_rows  : ${MAX_ROWS:-all}"
echo "  stats out : ${OUTBASE}.stats"
echo ""

# Run with cProfile
python -m cProfile -o "${OUTBASE}.stats" \
    -m mutation_scatter_plot.calculate_codon_frequencies.cli \
    --reference-infile="$REFERENCE" \
    --alignment-file="$ALIGNMENT" \
    --outfile-prefix="${OUTBASE}_out" \
    --padded-reference \
    --x-after-count \
    --overwrite \
    --threads="$THREADS"

# Dump top-50 hotspots by tottime and cumtime
python3 - "${OUTBASE}" <<'EOF'
import pstats, io, sys
OUTBASE = sys.argv[1]

for sort_key in ('tottime', 'cumulative'):
    s = io.StringIO()
    p = pstats.Stats(f'{OUTBASE}.stats', stream=s)
    p.sort_stats(sort_key)
    p.print_stats(50)
    txt = s.getvalue()
    with open(f'{OUTBASE}_top50_{sort_key}.txt', 'w') as f:
        f.write(txt)
    print(f"--- Top 50 by {sort_key} ---")
    lines = txt.splitlines()
    for line in lines[:60]:
        print(line)
    print("  ...")
    print()
EOF

echo ""
echo "Done. Retrieve files:"
echo "  ${OUTBASE}_top50_tottime.txt"
echo "  ${OUTBASE}_top50_cumulative.txt"
echo "  ${OUTBASE}.stats  (for snakeviz)"
