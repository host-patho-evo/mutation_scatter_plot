#! /bin/sh

# Kick out sequences not being exactly 3822 nucleotides long and place them into a separate file.
# The input file contains unique sequences sorted in a top-down ordering with FASTA ID containing
# an incidence count, followed by a 'x.' and then by an SHA256 checksum. It happens there are
# typically singleton sequences with incidence 1x and it is ipossible to track back the original
# sequence before read trimming, etc. Therefore, the sequence checksum is very handy.
#
# Notably it uses grep --no-group-separator to discard the '--' which is the delimiter added
# by 'grep -B' to the output stream.
# It will strip away trailing .fasta and add .filtered.fasta and .kicked_out.fasta.lsta
#
# The script uses filterbyname.sh from BBmap bundle to extract the offending sequences into
# a separate file. See https://sourceforge.net/projects/bbmap/
#
# Usage: kick.sh somefile.fastp.amplicons.clean.counts.fasta

p=`echo $1 | sed -e 's/.fasta//'`;
f="$1"
rm -f "$p".kicked_out.fasta "$p".kicked_out.fasta.lst
grep -v '^>' "$f" | awk '{ if ( length!=3822 ) print length,$1 }' | sort | uniq | sort | while read total_len sequence; do
  echo $sequence | grep --no-group-separator -B 1 -F -f - "$f" >> "$p".kicked_out.fasta;
done;
grep '^>' "$p".kicked_out.fasta | cut -c 2- | awk '{print $1}' | sort -nr > "$p".kicked_out.fasta.lst
original=`echo $1 | sed -e 's/.clean.counts.fasta//'`
filterbyname.sh include=f names="$p".kicked_out.fasta.lst in=$1 out="$p".filtered.fasta ignorejunk=t overwrite=t fastawrap=0 >"$p".filtered.log 2>&1
echo "Info: Wrote $p.filtered.fasta with entries long just 3822 nt"
filterbyname.sh include=t names="$p".kicked_out.fasta.lst in="$original".counts.fasta out="$original".counts.kicked_out.filterbyname.fasta overwrite=t fastawrap=0 >"$original".counts.kicked_out.filterbyname.log 2>&1
echo "Info: Wrote $original.counts.kicked_out.filterbyname.fasta"
