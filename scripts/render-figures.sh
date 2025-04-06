#! /usr/bin/env sh

# This is a main wrapper used to re-generate figures for all samples.

mutation_scatter_plot.py --xmin 340 --xmax 516 --tsv WT-7th-library-LSS__7-WU-FF.WUPA.WTref.frequencies.tsv --outfile WT-7th-library-LSS__7-WU-FF.WUPA.WTref.aa.frequencies.png --aminoacids
mutation_scatter_plot.py --xmin 340 --xmax 516 --tsv WT-7th-library-LSS__7-WU-FF.WUPA.WTref.frequencies.tsv --outfile WT-7th-library-LSS__7-WU-FF.WUPA.WTref.codon.frequencies.png
mutation_scatter_plot.py --xmin 340 --xmax 517 --tsv Alfa-7th-library-LSS__7-BR-CC.BRPA.WTref.frequencies.tsv --outfile Alfa-7th-library-LSS__7-BR-CC.BRPA.WTref.aa.frequencies.png --aminoacids
mutation_scatter_plot.py --xmin 340 --xmax 517 --tsv Alfa-7th-library-LSS__7-BR-CC.BRPA.WTref.frequencies.tsv --outfile Alfa-7th-library-LSS__7-BR-CC.BRPA.WTref.codon.frequencies.png
mutation_scatter_plot.py --xmin 340 --xmax 515 --tsv Beta-7th-library-LSS__7-SA-AA.SAPA.WTref.frequencies.tsv --outfile Beta-7th-library-LSS__7-SA-AA.SAPA.WTref.aa.frequencies.png --aminoacids
mutation_scatter_plot.py --xmin 340 --xmax 515 --tsv Beta-7th-library-LSS__7-SA-AA.SAPA.WTref.frequencies.tsv --outfile Beta-7th-library-LSS__7-SA-AA.SAPA.WTref.codon.frequencies.png
mutation_scatter_plot.py --xmin 340 --xmax 518 --tsv RBD-v48-7th-library-LSS__7-RY-BB.RYPA.WTref.frequencies.tsv --outfile RBD-v48-7th-library-LSS__7-RY-BB.RYPA.WTref.aa.frequencies.png --aminoacids
mutation_scatter_plot.py --xmin 340 --xmax 518 --tsv RBD-v48-7th-library-LSS__7-RY-BB.RYPA.WTref.frequencies.tsv --outfile RBD-v48-7th-library-LSS__7-RY-BB.RYPA.WTref.codon.frequencies.png
mutation_scatter_plot.py --xmin 338 --xmax 515 --tsv BA1-7th-library-LSS__7-OM-EE.OMPA.WTref.frequencies.tsv --outfile BA1-7th-library-LSS__7-OM-EE.OMPA.WTref.aa.frequencies.png --aminoacids
mutation_scatter_plot.py --xmin 338 --xmax 515 --tsv BA1-7th-library-LSS__7-OM-EE.OMPA.WTref.frequencies.tsv --outfile BA1-7th-library-LSS__7-OM-EE.OMPA.WTref.codon.frequencies.png
mutation_scatter_plot.py --xmin 430 --xmax 528 --tsv WT-unsorted__D5.WT1.WTref.frequencies.tsv --outfile WT-unsorted__D5.WT1.WTref.aa.frequencies.png --aminoacids
mutation_scatter_plot.py --xmin 430 --xmax 528 --tsv WT-unsorted__D5.WT1.WTref.frequencies.tsv --outfile WT-unsorted__D5.WT1.WTref.codon.frequencies.png
mutation_scatter_plot.py --xmin 430 --xmax 528 --tsv I358F-unsorted__E5.I3.WTref.frequencies.tsv --outfile I358F-unsorted__E5.I3.WTref.aa.frequencies.png --aminoacids
mutation_scatter_plot.py --xmin 430 --xmax 528 --tsv I358F-unsorted__E5.I3.WTref.frequencies.tsv --outfile I358F-unsorted__E5.I3.WTref.codon.frequencies.png
mutation_scatter_plot.py --xmin 430 --xmax 528 --tsv BA1-unsorted__F5.B1.WTref.frequencies.tsv --outfile BA1-unsorted__F5.B1.WTref.aa.frequencies.png --aminoacids
mutation_scatter_plot.py --xmin 430 --xmax 528 --tsv BA1-unsorted__F5.B1.WTref.frequencies.tsv --outfile BA1-unsorted__F5.B1.WTref.codon.frequencies.png
mutation_scatter_plot.py --xmin 430 --xmax 528 --tsv BA2-unsorted__G5.B2.WTref.frequencies.tsv --outfile BA2-unsorted__G5.B2.WTref.aa.frequencies.png --aminoacids
mutation_scatter_plot.py --xmin 430 --xmax 528 --tsv BA2-unsorted__G5.B2.WTref.frequencies.tsv --outfile BA2-unsorted__G5.B2.WTref.codon.frequencies.png
mutation_scatter_plot.py --xmin 430 --xmax 528 --tsv WT-2nd-round-of-sort__H5.WT2.WTref.frequencies.tsv --outfile WT-2nd-round-of-sort__H5.WT2.WTref.aa.frequencies.png --aminoacids
mutation_scatter_plot.py --xmin 430 --xmax 528 --tsv WT-2nd-round-of-sort__H5.WT2.WTref.frequencies.tsv --outfile WT-2nd-round-of-sort__H5.WT2.WTref.codon.frequencies.png
mutation_scatter_plot.py --xmin 430 --xmax 528 --tsv I358F-2nd-round-of-sort__A6.I3.WTref.frequencies.tsv --outfile I358F-2nd-round-of-sort__A6.I3.WTref.aa.frequencies.png --aminoacids
mutation_scatter_plot.py --xmin 430 --xmax 528 --tsv I358F-2nd-round-of-sort__A6.I3.WTref.frequencies.tsv --outfile I358F-2nd-round-of-sort__A6.I3.WTref.codon.frequencies.png
mutation_scatter_plot.py --xmin 430 --xmax 528 --tsv BA1-2nd-round-of-sort__B6.BM.WTref.frequencies.tsv --outfile BA1-2nd-round-of-sort__B6.BM.WTref.aa.frequencies.png --aminoacids
mutation_scatter_plot.py --xmin 430 --xmax 528 --tsv BA1-2nd-round-of-sort__B6.BM.WTref.frequencies.tsv --outfile BA1-2nd-round-of-sort__B6.BM.WTref.codon.frequencies.png
mutation_scatter_plot.py --xmin 430 --xmax 528 --tsv BA2-2nd-round-of-sort__C6.BA2.WTref.frequencies.tsv --outfile BA2-2nd-round-of-sort__C6.BA2.WTref.aa.frequencies.png --aminoacids
mutation_scatter_plot.py --xmin 430 --xmax 528 --tsv BA2-2nd-round-of-sort__C6.BA2.WTref.frequencies.tsv --outfile BA2-2nd-round-of-sort__C6.BA2.WTref.codon.frequencies.png
mutation_scatter_plot.py --xmin 430 --xmax 528 --tsv WT-4th-round-of-sort__D6.WT4.WTref.frequencies.tsv --outfile WT-4th-round-of-sort__D6.WT4.WTref.aa.frequencies.png --aminoacids
mutation_scatter_plot.py --xmin 430 --xmax 528 --tsv WT-4th-round-of-sort__D6.WT4.WTref.frequencies.tsv --outfile WT-4th-round-of-sort__D6.WT4.WTref.codon.frequencies.png
mutation_scatter_plot.py --xmin 430 --xmax 528 --tsv I358F-4th-round-of-sort__E6.I34.WTref.frequencies.tsv --outfile I358F-4th-round-of-sort__E6.I34.WTref.aa.frequencies.png --aminoacids
mutation_scatter_plot.py --xmin 430 --xmax 528 --tsv I358F-4th-round-of-sort__E6.I34.WTref.frequencies.tsv --outfile I358F-4th-round-of-sort__E6.I34.WTref.codon.frequencies.png
mutation_scatter_plot.py --xmin 430 --xmax 528 --tsv BA1-4th-round-of-sort__F6.BA1.WTref.frequencies.tsv --outfile BA1-4th-round-of-sort__F6.BA1.WTref.aa.frequencies.png --aminoacids
mutation_scatter_plot.py --xmin 430 --xmax 528 --tsv BA1-4th-round-of-sort__F6.BA1.WTref.frequencies.tsv --outfile BA1-4th-round-of-sort__F6.BA1.WTref.codon.frequencies.png
mutation_scatter_plot.py --xmin 430 --xmax 528 --tsv BA2-4th-round-of-sort__G6.BA2.WTref.frequencies.tsv --outfile BA2-4th-round-of-sort__G6.BA2.WTref.aa.frequencies.png --aminoacids
mutation_scatter_plot.py --xmin 430 --xmax 528 --tsv BA2-4th-round-of-sort__G6.BA2.WTref.frequencies.tsv --outfile BA2-4th-round-of-sort__G6.BA2.WTref.codon.frequencies.png

# move the resulting files to their respective folders
d=`cwd`
rm -f *.aa.frequencies.legend.png # zap empty figures
mkdir -p "$d"data/outputs "$d"data/outputs/codon "$d"data/outputs/codon/legend "$d"data/outputs/codon/jpg "$d"data/outputs/codon/png "$d"data/outputs/codon/pdf "$d"data/outputs/codon/html "$d"data/outputs/codon/color "$d"data/outputs/aa "$d"data/outputs/aa/jpg "$d"data/outputs/aa/png "$d"data/outputs/aa/pdf "$d"data/outputs/aa/html "$d"data/outputs/aa/color/ 
mv *.codon.frequencies.legend.* "$d"/data/outputs/codon/legend/
for ext in jpg png pdf html; do
    mv *.codon.frequencies."$ext" "$d"/data/outputs/codon/"$ext"/
    mv *.aa.frequencies."$ext" "$d"/data/outputs/aa/"$ext"
done

mv *.codon.frequencies.colors.tsv "$d"/data/outputs/codon/color/
mv *.aa.frequencies.colors.tsv "$d"/data/outputs/aa/color/
mv *.frequencies.actually_rendered.tsv "$d"/data/outputs/

