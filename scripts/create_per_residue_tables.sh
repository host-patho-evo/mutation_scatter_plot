#! /usr/bin/env bash

rm -f data/outputs/per_residue_tables/*.tsv

cd data/intermediates
START=`for f in *.WTref.gofasta.frequencies.tsv; do head -n 1 "$f" | awk -F '\t' '{print $1}'; done | sort -n | head -n 1`
END=`for f in *.WTref.gofasta.frequencies.tsv ; do tail -n 1 "$f" | awk -F '\t' '{print $1}'; done | sort -n | tail -n 1`

for (( i=START ; i <= END; i++ )); do for f in *WTref.frequencies.tsv; do awk '{if ($1 == '"$i"') print FILENAME, $0}' "$f" >> ../outputs/per_residue_tables/"$i".tsv; done; done

cd ../outputs/per_residue_tables
printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\n' "filename" "position" "original_aa" "mutant_aa" "frequency" "original_codon" "mutant_codon" > header.csv || exit 255
for f in *.tsv; do
  # Alpha-7th-library-LSS__7-BR-CC.BRPA.WTref.frequencies.tsv 340	E	K	0.000424	GAA	AAA
  cat header.csv "$f" > "$f".$$
  mv "$f".$$ "$f"
done
rm -f header.csv
