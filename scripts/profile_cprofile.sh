#!/usr/bin/env bash
# cProfile run for calculate_codon_frequencies on the target machine.
# Run from project root:  bash scripts/profile_cprofile.sh <alignment.fasta>
#
# Outputs:
#   /tmp/profile_codon_N.stats   - cProfile binary (open with pstats or snakeviz)
#   /tmp/profile_codon_N.txt     - pre-rendered top-50 hotspots
#
# Usage examples:
#   bash scripts/profile_cprofile.sh tests/inputs/test2_full.fasta
#   bash scripts/profile_cprofile.sh /data/large_alignment.fasta --threads 1
#   bash scripts/profile_cprofile.sh /data/large_alignment.fasta --threads 64

set -euo pipefail

ALIGNMENT="${1:-tests/inputs/test2_full.fasta}"
REFERENCE="${REFERENCE:-tests/inputs/MN908947.3_S_full.fasta}"
THREADS="${2:-1}"
OUTBASE="/tmp/profile_codon_t${THREADS}"

export PYTHONPATH="${PYTHONPATH:-$(pwd)/src}"

echo "=== cProfile: calculate_codon_frequencies ==="
echo "  alignment : $ALIGNMENT  ($(du -h "$ALIGNMENT" | cut -f1))"
echo "  reference : $REFERENCE"
echo "  threads   : $THREADS"
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
python - <<EOF
import pstats, io

for sort_key in ('tottime', 'cumulative'):
    s = io.StringIO()
    p = pstats.Stats('${OUTBASE}.stats', stream=s)
    p.sort_stats(sort_key)
    p.print_stats(50)
    txt = s.getvalue()
    with open('${OUTBASE}_top50_${sort_key}.txt', 'w') as f:
        f.write(txt)
    print(f"--- Top 50 by {sort_key} ---")
    # Print just the header + first 20 lines of stats table
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
