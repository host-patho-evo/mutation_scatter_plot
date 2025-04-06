#! /usr/bin/env sh

# This script has been used to traverse per-sample TSV files and
# to output observed frequencies for a single codon position.
# This allows to easily compare frequencies across samples.

START=`for f in *.frequencies.tsv ; do head -n 1 "$f" | awk -F '\t' '{print $1}'; done | sort -n | head -n 1`

