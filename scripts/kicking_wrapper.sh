#! /bin/sh

# Kick out sequences not being 3822 nt long after padding. In this way we move
# sequences with INSertions occurring in the sample reads into a separate file.
# This allows to use pair-wise alignment results and use them as a multiple-sequence
# alignment. Alternativley, the reference sequence would eed to be inflated by padding
# dashes to create space for the INSertion in the sample.

find . -name \*fastp.amplicons.clean.counts.fasta | sort | while read f; do bash kick.sh $f; done

# For eventual analysis, use NCBI blastn to align the discarded sample reads to
# their proper reference.

find . -name \*.YEASTref.scores_above_84.fastp.amplicons.counts.kicked_out.filterbyname.fasta | sort | uniq | grep -v results | grep -v 'MN908947.3_or_OR995580_or_T3_or_T5' | while read f; do
  if [ -s "$f" ]; then
    d=`dirname $f`;
    p=`basename $f .fasta`;
    fileref=`echo $p | cut -c 3-4` # cut out from 2-T5 the 'T5' portion
    echo "${d}/${p} $fileref"
    blastn -task blastn -reward 2 -max_hsps 1 -num_alignments 1 -word_size 4 -num_threads 16 -dust no -evalue 1e-25 -db ${fileref}.fasta -query $f > $d/$p.blastn.txt
  fi
done

find . -name \*.WTref.scores_above_84.fastp.amplicons.counts.kicked_out.filterbyname.fasta | sort | uniq | grep -v results | grep -v 'MN908947.3_or_OR995580_or_T3_or_T5' | while read f; do
  if [ -s "$f" ]; then
    d=`dirname $f`;
    p=`basename $f .fasta`;
    fileref=MN908947.3_S
    echo "${d}/${p} $fileref"
    blastn -task blastn -reward 2 -max_hsps 1 -num_alignments 1 -word_size 4 -num_threads 16 -dust no -evalue 1e-25 -db ${fileref}.fasta -query $f > $d/$p.blastn.txt
  fi
done
