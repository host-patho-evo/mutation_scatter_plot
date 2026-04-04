#! /usr/bin/env python3

# >A00808:1538:HWL7GDRX3:2:2101:7491:1125 minus
# AGAAAGTACTACTAATCTGTATGGTAGGTGACCAACACCATAAGTGGGTCGGAAACCATATGATTGTAAAGGAGAGTAACAATTAGGACCTGCAACACAAATACAAGGTTTGATACCGGCCAGATAGATATCAGTTGAAATATATCTCTCAAAAGGTTTCAGTTTAGACTTCCTAAACAATATATACCAGTAATCATAAATACCACTATGCTTAGAATCAAGCATGTTAGAATTCCAAGCTATATCGCAGCCTTTTAAATCATCTGGTAATTTTTAATTATAATCAGCAATATTTCCAGTTTGCCCTGGAGCGATTTGGCTGACTTCATTACCTTTAATTACAAATGAATCTGCATAGACATGAGTAAAGCAGTGATCATTTATTTTAGTAGGAGACACTCC
# >A00808:1538:HWL7GDRX3:2:2101:7419:1157 minus
# AGAAAGTACTACTACTCTGTATGGTTGGTGACCAACACCATAAGTGGGTCGGAAACCATATGATTGTAAAGGAGAGTAACAATTAGGACCTGCAACACCATTACAAGGTTTGTTACCGGCCTGATAGATTTCAGTTGAAATATCTCTCTCAAAAGGGTTGAGTTTAGACTTCCTAAACAATCTATACCAGTAATCATAATTACCACTATGCTTAGAATCAAGCTTGTTAGAATTCCAAGCTATAACGCAGCCTGTAAAATCATCTGGTAATTTATAATTATAATCAGCAATATTTCCAGTTTGCCCTGGAGCGATTTGGCTGACTTCATTACCTTTAATTACAAATGAATCTGCATAGACATTAGTAAAGCAGAGATCATTTAATTTAGTAGGAGACACTCC
# >A00808:1538:HWL7GDRX3:2:2101:21296:1204 plus
# GGAGTGTCTCCTACTAAATTAAATGATCTCTGCTTTACTAATGTCTATGCAGATTCATTTGTAATTAAAGGTAATGAAGTCAGCCAAATCGCTCCAGGGCAAACTGGAAATATTGCTGATTATAATTATAAATTACCAGATGATTTTACAGGCTGCGTTATAGCTTGGAATTCTAACAAGCTTGATTCTAAGCATAGTGGTAATTATGATTACTGGTATAGATTGTTTAGGAAGTCTAAACTCAAACCTTTTGAGAGAGATATTTCAACTGAAATCTATCAGGCCGGTAACAAACCTTGTAAAGG---TAAAGGTCCTAATTGTTACTTTCCTTTACAATCATATGGTTTCCGACCCACTTATGGTGTTGGTCACCAACCATACAGAGTAGTAGTACTTTCT
# >A00808:1538:HWL7GDRX3:2:2101:24831:1219 minus
# AGAAAGTACTACTACTCTGTATGGTTGGTGACCAACACCATAAGTGGGTCGGAAACCATATGATTGTAAAGGAAAGTAACAATTAGGACCTTTA---CCTTTACAAGGTTTGTTACCGGCCTGATAGATTTCAGTTGAAATATCTCTCTCAAAAGGTTTGAGTTTAGACTTCCTAAACAATCTATACCAGTAATCATAATTACCACTATGCTTAGAATCAAGCTTGTTAGAATTCCAAGCTATAACGCAGCCTGTAAAATCATCTGGTAATTTATAATTATAATCAGCAATATTTCCAGTTTGCCCTGGAGCGATTTGGCTGACTTCATTACCTTTAATTACAAATGAATCTGCATAGACATTAGTAAAGCAGAGATCATTTAATTTAGTAGGAGACACTCC
# >A00808:1538:HWL7GDRX3:2:2101:25021:1235 minus
# AGAAAGTACTACTACTCTGTATGGTTGGTGACCAACACCATAAGTGGGTCGGAAACCATATGATTGTAAAGGAAAGTAACAATTAGGACCTTTA---CCTTTACAAGGTTTGTTACCGGCCTGATAGATTTCAGTTGAAATATCTCTCTCAAAAGGTTTGAGTTTAGACTTCCTAAACAATCTATACCAGTAATCATAATTACCACTATGCTTAGAATCAAGCTTGTTAGAATTCCAAGCTATAACGCAGCCTGTAAAATCATCTGGTAATTTATAATTATAATCAGCAATATTTCCAGTTTGCCCTGGAGCGATTTGGCTGACTTCATTACCTTTAATTACAAATGAATCTGCATAGACATTAGTAAAGCAGAGATCATTTAATTTAGTAGGAGACACTCC

"""
Parse TSV input stream with BLASTN results and figure out position of false,
frame-breaking InDels.


Examples of erroneous DELetions:

>MN908947.3_S_protein  3822nt 1274codons
Length=3822

 Score = 404 bits (447),  Expect = 3e-116
 Identities = 227/228 (99%), Gaps = 1/228 (0%)
 Strand=Plus/Plus

Query  1     AATCTTGATTCTAAGGTTGGTGGTAATTATAATTACCTGTATAGATTGTTTAGGAAGTCT  60
             ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
Sbjct  1318  AATCTTGATTCTAAGGTTGGTGGTAATTATAATTACCTGTATAGATTGTTTAGGAAGTCT  1377

Query  61    AATCTCAAACCTTT-GAGAGAGATATTTCAACTGAAATCTATCAGGCCGGTAGCACACCT  119
             |||||||||||||| |||||||||||||||||||||||||||||||||||||||||||||
Sbjct  1378  AATCTCAAACCTTTTGAGAGAGATATTTCAACTGAAATCTATCAGGCCGGTAGCACACCT  1437

Query  120   TGTAATGGTGTTGAAGGTTTTAATTGTTACTTTCCTTTACAATCATATGGTTTCCAACCC  179
             ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
Sbjct  1438  TGTAATGGTGTTGAAGGTTTTAATTGTTACTTTCCTTTACAATCATATGGTTTCCAACCC  1497

Query  180   ACTAATGGTGTTGGTTACCAACCATACAGAGTAGTAGTACTTTCTTTT  227
             ||||||||||||||||||||||||||||||||||||||||||||||||
Sbjct  1498  ACTAATGGTGTTGGTTACCAACCATACAGAGTAGTAGTACTTTCTTTT  1545


>MN908947.3_S_protein  3822nt 1274codons
Length=3822

 Score = 404 bits (447),  Expect = 3e-116
 Identities = 227/228 (99%), Gaps = 1/228 (0%)
 Strand=Plus/Plus

Query  1     AATCTTGATTCTAAGGTTGGTGGTAATTATAATTACCTGTATAGATTGTTTAGGAAGTCT  60
             ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
Sbjct  1318  AATCTTGATTCTAAGGTTGGTGGTAATTATAATTACCTGTATAGATTGTTTAGGAAGTCT  1377

Query  61    AATCTCAAACCTTTTGAGAGAGATATTTCAACTGAAATCTATCAGGCCGGTAGCACACCT  120
             ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
Sbjct  1378  AATCTCAAACCTTTTGAGAGAGATATTTCAACTGAAATCTATCAGGCCGGTAGCACACCT  1437

Query  121   TGTAATGGTGTTGAAGGTTTTAATTGTTA-TTTCCTTTACAATCATATGGTTTCCAACCC  179
             ||||||||||||||||||||||||||||| ||||||||||||||||||||||||||||||
Sbjct  1438  TGTAATGGTGTTGAAGGTTTTAATTGTTACTTTCCTTTACAATCATATGGTTTCCAACCC  1497

Query  180   ACTAATGGTGTTGGTTACCAACCATACAGAGTAGTAGTACTTTCTTTT  227
             ||||||||||||||||||||||||||||||||||||||||||||||||
Sbjct  1498  ACTAATGGTGTTGGTTACCAACCATACAGAGTAGTAGTACTTTCTTTT  1545


>MN908947.3_S_protein  3822nt 1274codons
Length=3822

 Score = 343 bits (379),  Expect = 9e-98
 Identities = 217/235 (92%), Gaps = 4/235 (2%)
 Strand=Plus/Plus

Query  1     AACAATCTTGATTCTAAGGTTGGTGGTAATTACAATTACCTGTTTAGATTGTTTAGGAAG  60
             |||||||||||||||||||||||||||||||| |||||||||| ||||||||||||||||
Sbjct  1315  AACAATCTTGATTCTAAGGTTGGTGGTAATTATAATTACCTGTATAGATTGTTTAGGAAG  1374

Query  61    TCTAATCTCAAACCTTTTGAGAAAGATATTTCAACTGAAATCTATCAGGCCGGCAGCACA  120
             |||||||||||||||||||||| |||||||||||||||||||||||||||||| ||||||
Sbjct  1375  TCTAATCTCAAACCTTTTGAGAGAGATATTTCAACTGAAATCTATCAGGCCGGTAGCACA  1434

Query  121   CCTTGTAATGATGTTGAAGGTTTTAATTGTTATTTTCCTCTACAATCATATGGTTTCCAT  180
             |||||||||| ||||||||||||||||||||| |||||| |||||||||||||||||||
Sbjct  1435  CCTTGTAATGGTGTTGAAGGTTTTAATTGTTACTTTCCTTTACAATCATATGGTTTCCAA  1494

Query  181   CGTACCAATGGTGTTGGTACGATGGAAACCATACAGAGTAGTAGTACTTTCTTTT  235
             |  || ||||||||||||    |   |||||||||||||||||||||||||||||
Sbjct  1495  CCCACTAATGGTGTTGGT----TACCAACCATACAGAGTAGTAGTACTTTCTTTT  1545



blastn -task blastn -reward 2 -max_hsps 1 -num_alignments 1 -word_size 4 -num_threads 4 -dust no -evalue 1e-5 -db T3.fasta -query wrong4.fa

>T3
Length=3822

 Score = 326 bits (361),  Expect = 6e-93
 Identities = 206/221 (93%), Gaps = 4/221 (2%)
 Strand=Plus/Plus

Query  1     ATTGCCTGGAACAGCAACAAACTGGACAGCAAGGTGGGAGGCAACTACAACTACCTCTAC  60
             ||||||||||||||||||||||||||||||||||||   |||||||||||||||||||||
Sbjct  1300  ATTGCCTGGAACAGCAACAAACTGGACAGCAAGGTGTCTGGCAACTACAACTACCTCTAC  1359

Query  61    AGACTGTTCAGGAAGAGCAACCTGAAACCATTTGAGAGGGACATCAGCACAGAGATTTAC  120
             ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
Sbjct  1360  AGACTGTTCAGGAAGAGCAACCTGAAACCATTTGAGAGGGACATCAGCACAGAGATTTAC  1419

Query  121   CAGGCAGGCAACAAACCATGTAATGGAGTGGCCGGCTTTAACTGTTACCATCCACTC--G  178
             ||||| |||||||||||||||||||||||||||||||| |||||||||  |||||||  |
Sbjct  1420  CAGGCTGGCAACAAACCATGTAATGGAGTGGCCGGCTTCAACTGTTACTTTCCACTCAAG  1479

Query  179   TCCCATATGCCTTCCGTCCAACCTACGGAGTTGGTCATCAA  219
             |||  ||| | |||||||||||||||||||| || ||||||
Sbjct  1480  TCC--TATTCTTTCCGTCCAACCTACGGAGTGGGCCATCAA  1518



blastn -task blastn -reward 1 -max_hsps 1 -num_alignments 1 -word_size 4 -num_threads 4 -dust no -evalue 1e-5 -db T3.fasta -query wrong4.fa

>T3
Length=3822

 Score = 307 bits (155),  Expect = 2e-87
 Identities = 203/219 (93%), Gaps = 0/219 (0%)
 Strand=Plus/Plus

Query  1     ATTGCCTGGAACAGCAACAAACTGGACAGCAAGGTGGGAGGCAACTACAACTACCTCTAC  60
             ||||||||||||||||||||||||||||||||||||   |||||||||||||||||||||
Sbjct  1300  ATTGCCTGGAACAGCAACAAACTGGACAGCAAGGTGTCTGGCAACTACAACTACCTCTAC  1359

Query  61    AGACTGTTCAGGAAGAGCAACCTGAAACCATTTGAGAGGGACATCAGCACAGAGATTTAC  120
             ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
Sbjct  1360  AGACTGTTCAGGAAGAGCAACCTGAAACCATTTGAGAGGGACATCAGCACAGAGATTTAC  1419

Query  121   CAGGCAGGCAACAAACCATGTAATGGAGTGGCCGGCTTTAACTGTTACCATCCACTCGTC  180
             ||||| |||||||||||||||||||||||||||||||| |||||||||  |||||||
Sbjct  1420  CAGGCTGGCAACAAACCATGTAATGGAGTGGCCGGCTTCAACTGTTACTTTCCACTCAAG  1479

Query  181   CCATATGCCTTCCGTCCAACCTACGGAGTTGGTCATCAA  219
              | ||| | |||||||||||||||||||| || ||||||
Sbjct  1480  TCCTATTCTTTCCGTCCAACCTACGGAGTGGGCCATCAA  1518


>T3
Length=3822

 Score = 310 bits (343),  Expect = 5e-88
 Identities = 199/216 (92%), Gaps = 1/216 (0%)
 Strand=Plus/Plus

Query  1     GCCTGGAACAGCAACAACCTGGACAGCAAGAAGGGAGGCAACTACAACTACCTCTACAGA  60
             ||||||||||||||||| ||||||||||||  |   ||||||||||||||||||||||||
Sbjct  1303  GCCTGGAACAGCAACAAACTGGACAGCAAGGTGTCTGGCAACTACAACTACCTCTACAGA  1362

Query  61    CTGTTCAGGAAGAGCAAACTGAAACCATTTGAGAGGGACATCAGCACAGAGATTTACCAG  120
             ||||||||||||||||| ||||||||||||||||||||||||||||||||||||||||||
Sbjct  1363  CTGTTCAGGAAGAGCAACCTGAAACCATTTGAGAGGGACATCAGCACAGAGATTTACCAG  1422

Query  121   GCAGGCAACAAACCATGTAATGGAGTGGCCGGCTTAAACTGTTACCATCCACTCCAGT-A  179
             || |||||||||||||||||||||||||||||||| |||||||||  ||||||| |||
Sbjct  1423  GCTGGCAACAAACCATGTAATGGAGTGGCCGGCTTCAACTGTTACTTTCCACTCAAGTCC  1482

Query  180   TATGGCTTCCGTCCAACCTACGGAGTGGGCCATCAA  215
             |||   ||||||||||||||||||||||||||||||
Sbjct  1483  TATTCTTTCCGTCCAACCTACGGAGTGGGCCATCAA  1518


>MN908947.3_S_protein  3822nt 1274codons
Length=3822

 Score = 530 bits (587),  Expect = 4e-154
 Identities = 297/298 (99%), Gaps = 1/298 (0%)
 Strand=Plus/Plus

Query  1     ACAGGCTGCGTTATAGCTTGGAATTCTAACAATCTTGATTCTAAGGTTGGTGGGTAATTA  60
             ||||||||||||||||||||||||||||||||||||||||||||||||||||| ||||||
Sbjct  1288  ACAGGCTGCGTTATAGCTTGGAATTCTAACAATCTTGATTCTAAGGTTGGTGG-TAATTA  1346

Query  61    TAATTACCTGTATAGATTGTTTAGGAAGTCTAATCTCAAACCTTTTGAGAGAGATATTTC  120
             ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
Sbjct  1347  TAATTACCTGTATAGATTGTTTAGGAAGTCTAATCTCAAACCTTTTGAGAGAGATATTTC  1406

Query  121   AACTGAAATCTATCAGGCCGGTAGCACACCTTGTAATGGTGTTGAAGGTTTTAATTGTTA  180
             ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
Sbjct  1407  AACTGAAATCTATCAGGCCGGTAGCACACCTTGTAATGGTGTTGAAGGTTTTAATTGTTA  1466

Query  181   CTTTCCTTTACAATCATATGGTTTCCAACCCACTAATGGTGTTGGTTACCAACCATACAG  240
             ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
Sbjct  1467  CTTTCCTTTACAATCATATGGTTTCCAACCCACTAATGGTGTTGGTTACCAACCATACAG  1526

Query  241   AGTAGTAGTACTTTCTTTTGAACTTCTACATGCACCAGCAACTGTTTGTGGACCTAAA  298
             ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
Sbjct  1527  AGTAGTAGTACTTTCTTTTGAACTTCTACATGCACCAGCAACTGTTTGTGGACCTAAA  1584


>MN908947.3_S_protein  3822nt 1274codons
Length=3822

 Score = 292 bits (323),  Expect = 2e-82
 Identities = 222/261 (85%), Gaps = 1/261 (0%)
 Strand=Plus/Plus

Query  1     GATTTTAGGGTTGGTGGTAATTTTTATTTAATGTATAGATTGGTTTGGGAATTTAAGCTT  60
             |||| || |||||||||||||| | |||   ||||||||||| || || | | ||| || 
Sbjct  1324  GATTCTAAGGTTGGTGGTAATTATAATTACCTGTATAGATTGTTTAGGAAGTCTAATCTC  1383

Query  61    AAACCTTTTTAGAGGGGTATTTTAACAGAAATTTAACAAGGCGGTAAAAAACCTTGTAAT  120
             ||||||||| |||| | ||||| ||| ||||| || || | |||||  | ||||||||||
Sbjct  1384  AAACCTTTTGAGAGAGATATTTCAACTGAAATCTATCAGGCCGGTAGCACACCTTGTAAT  1443

Query  121   TGTGTAGAGGGTTTTAAT-GTTACTTTTCCTTACAATTATATTGTTTTCGACCCACCTAT  179
              |||| || ||||||||| |||||||| | ||||||| |||| |||| | ||||||  ||
Sbjct  1444  GGTGTTGAAGGTTTTAATTGTTACTTTCCTTTACAATCATATGGTTTCCAACCCACTAAT  1503

Query  180   GGTGTTGGTTACCAACCATACAGAGTAGTAGTACTTTCTTTTGAACTTCTACATGGACCA  239
             ||||||||||||||||||||||||||||||||||||||||||||||||||||||| ||||
Sbjct  1504  GGTGTTGGTTACCAACCATACAGAGTAGTAGTACTTTCTTTTGAACTTCTACATGCACCA  1563

Query  240   GCAACTGTTTGTGGACCTAAA  260
             |||||||||||||||||||||
Sbjct  1564  GCAACTGTTTGTGGACCTAAA  1584


>MN908947.3_S_protein  3822nt 1274codons
Length=3822

 Score = 267 bits (295),  Expect = 5e-75
 Identities = 202/237 (85%), Gaps = 1/237 (0%)
 Strand=Plus/Plus

Query  1     AATTACCGGTATAGATTGTTTAGCAAGTATAAACTCAAACCTTTAGATAGAGAATTTTCA  60
             ||||||| ||||||||||||||| |||| ||| ||||||||||| || |||||  |||||
Sbjct  1348  AATTACCTGTATAGATTGTTTAGGAAGTCTAATCTCAAACCTTTTGAGAGAGATATTTCA  1407

Query  61    ACTGAAATATATCGGGCAGGTAGAACAC-TTGTAATGGTGTTGAAGGTTTTAATTGTTAG  119
             |||||||| |||| ||| ||||| |||| ||||||||||||||||||||||||||||||
Sbjct  1408  ACTGAAATCTATCAGGCCGGTAGCACACCTTGTAATGGTGTTGAAGGTTTTAATTGTTAC  1467

Query  120   TTTCATTTACGATCATATGGTTTCCGACACACTTATGGTGCTGGTTACCAACCATAAACA  179
             |||| ||||| |||||||||||||| || |||| |||||| ||||||||||||||| | |
Sbjct  1468  TTTCCTTTACAATCATATGGTTTCCAACCCACTAATGGTGTTGGTTACCAACCATACAGA  1527

Query  180   CTACAACTAATATCTCTTCAAAAAATACATGCACCAGCTACTGTTTGTGGACCTAAA  236
              ||  | || | ||| || ||    ||||||||||||| ||||||||||||||||||
Sbjct  1528  GTAGTAGTACTTTCTTTTGAACTTCTACATGCACCAGCAACTGTTTGTGGACCTAAA  1584



Show lines not having length 3822:

grep -v '^>' C2/2-T3BA1-C1.IT3N__scores_above_84.yeast.fastp.amplicons.counts.clean.fasta | awk '{ if ( length!=3822 ) print }' |  grep -v 3822 | sort | uniq | sort | while read s; do grep -B 1 -F -f - C2/2-T3BA1-C1.IT3N__scores_above_84.yeast.fastp.amplicons.counts.clean.fasta; done


Example
gzip -dc WW_waste-water_samples__WW-P1.parental@JZ244_54_F.CATGCTCAATAGC.HWL7GDRX3.2.fastp.fastq.gz | fastq_to_fasta | blastn -task blastn -reward 2 -max_hsps 1 -num_alignments 1 -word_size 4 -num_threads 64 -dust no -evalue 1e-5 -outfmt "6 qacc sstart send sstrand evalue bitscore score length pident nident mismatch positive gapopen gaps ppos qseq sseq" -db MN908947.3_S.fasta -query - 2>/dev/null | awk '{print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17}' | drop_erroneous_insertions.py --infile=- outfile=- | reversecomplement_reads_on_minus.py --infile=- | translate.py --infile=- | grep -v '^>' | sort | uniq -c | sort -nr

echo "1x.d2b2b8089265daaac94d3c3553ddf45d00841b13da817b2f2f703cd834d05437 1288 1584 plus ACTGGCTGCGTTATAGCTTGGAATTCTAACAAGCTTGATTCTAAGGTTGGTGGTAATTATAATTACCTGTATAGATTGTTTAGGAAGTCTAATCTCAAACCTTTTGAGAGAGATATTTCAACTGAAATCTATCAGGCCGGTAACAAACCTTGTAATGGTGTTGCAGGTTTTAATTGTTACTTTCCTTTACGATCATATGGTTTCCGACCCACTTATGGTGTTGGTCACCAACCATACAGA---GTAGTACTTTCATTTGAACTTCTACATGCACCAGGAACTGTTTGTGGACCTAcAA ACTGGCTGCGTTATAGCTTGGAATTCTAACAAGCTTGATTCTAAGGTTGGTGGTAATTATAATTACCTGTATAGATTGTTTAGGAAGTCTAATCTCAAACCTTTTGAGAGAGATATTTCAACTGAAATCTATCAGGCCGGTAACAAACCTTGTAATGGTGTTGCAGGTTTTAATTGTTACTTTCCTTTACGATCATATGGTTTCCGACCCACTTATGGTGTTGGTCACCAACCATACAGANNNGTAGTACTTTCATTTGAACTTCTACATGCACCAGGAACTGTTTGTGGACCTA-AA" | drop_erroneous_insertions.py --infile=- --outfile=- --drop-misinsertions
>1x.d2b2b8089265daaac94d3c3553ddf45d00841b13da817b2f2f703cd834d05437 1288 1584 plus
ACTGGCTGCGTTATAGCTTGGAATTCTAACAAGCTTGATTCTAAGGTTGGTGGTAATTATAATTACCTGTATAGATTGTTTAGGAAGTCTAATCTCAAACCTTTTGAGAGAGATATTTCAACTGAAATCTATCAGGCCGGTAACAAACCTTGTAATGGTGTTGCAGGTTTTAATTGTTACTTTCCTTTACGATCATATGGTTTCCGACCCACTTATGGTGTTGGTCACCAACCATACAGA---GTAGTACTTTCATTTGAACTTCTACATGCACCAGGAACTGTTTGTGGACCTAAA
"""

import sys
import os
import re
from optparse import OptionParser

VERSION = "202601061600"

myparser = OptionParser(version="%s version %s" % ('%prog', VERSION))
myparser.add_option("--infile", action="store", type="string", dest="infile", default='',
    help="Input FASTA file with start and stop coordinates and a word plus or minus as the last word of the header, use minus (dash) for stdin")
myparser.add_option("--outfile", action="store", type="string", dest="outfile", default='',
    help="Output filename with words plus and minus removed from FASTA header lines, use minus (dash) for stdout")
myparser.add_option("--infileformat", action="store", type="string", dest="infileformat", default='fasta',
    help="Input file format (Default: fasta)")
myparser.add_option("--outfileformat", action="store", type="string", dest="outfileformat", default='fasta-2line',
    help="Output file format (Default: fasta-2line)")
myparser.add_option("--drop-misinsertions", action="store_true", dest="drop_misinsertions", default=False,
    help="Drop evidently mis-INSerted single- or double- nucleotides breaking the reading frame which must be a sequencing error. It drops the nucleotides from sample sequence if there is a dash in the reference.")
myparser.add_option("--debug", action="store", type="int", dest="debug", default=0,
    help="Enable debug output")
(myoptions, myargs) = myparser.parse_args()

r=re.compile('-')

def find_insertions_in_sseq_and_drop_them_from_qseq(qseq, sseq):
    """Replace instances of single or double dashes and replace them with nothing.
    We want to look for INSertions relative to reference sequence which are
    not adjacent 3 INSertions.

    Maybe one could guess that there is an nucleotide over-call on at one of the INSerted
    nucleotides or immediately before the INSertion or immediately after the INSertion.
    For simplicity, we retain the first three nucleotides and drop the trailing one or two.

    We need to pass the sseq downstream in the output so that the number of INSertions
    in the reference can be used for subtraction.
    """
    _protect_multiples_of_3 = sseq.replace('---','###')
    
    if myoptions.debug:
        print(f"BEFORE: {qseq}")
        
    _qseq = "".join(q for q, p in zip(qseq, _protect_multiples_of_3) if p != '-')
    
    if myoptions.debug:
        print(f"AFTER:  {_qseq}")
        
    return _qseq

if __name__ == "__main__":
    if myoptions.infile == '-':
        _infile = sys.stdin
    elif os.path.exists(myoptions.infile):
        _infile = myoptions.infile
    else:
        raise RuntimeError("File %s does not exist" % myoptions.infile)

    if not myoptions.outfile:
        _outfile = sys.stdout
    elif myoptions.outfile == '-':
        _outfile = sys.stdout
    else:
        _outfile = myoptions.outfile

    for _line in _infile:
        myid, checksum, sstart, send, sstrand, evalue, bitscore, score, length, pident, nident, mismatch, positive, gapopen, gaps, ppos, qseq, sseq = 18 * ['']
        _fasta_header_items = _line.split()
        # qacc sstart send sstrand qseq sseq
        # qacc sstart send sstrand evalue bitscore score length pident nident mismatch positive gapopen gaps ppos qseq sseq
        if len(_fasta_header_items) == 17:
            myid, sstart, send, sstrand, evalue, bitscore, score, length, pident, nident, mismatch, positive, gapopen, gaps, ppos, qseq, sseq = _fasta_header_items
            if 'x.' in myid:
                myid, checksum = myid.split('.')
            else:
                checksum = ''
        elif len(_fasta_header_items) == 18:
            myid, checksum, sstart, send, sstrand, evalue, bitscore, score, length, pident, nident, mismatch, positive, gapopen, gaps, ppos, qseq, sseq = _fasta_header_items
        elif len(_fasta_header_items) == 6:
            myid, sstart, send, sstrand, qseq, sseq = _fasta_header_items
            if 'x.' in myid:
                myid, checksum = myid.split('.')
            else:
                checksum = ''
        elif len(_fasta_header_items) == 7:
            myid, checksum, sstart, send, sstrand, qseq, sseq = _fasta_header_items
        else:
            # >A00877:1511:HWLFTDRX3:2:2101:30572:1141 1541 1141 minus
            raise ValueError("Cannot parse line %s" % _line)

        _number_of_insertions_in_sseq = sseq.lstrip('-').rstrip('-').count('-') # we need the number of INSertions in the aligned sseq because that will make a difference when just subtracting coordinates
        if myoptions.drop_misinsertions and _number_of_insertions_in_sseq:
            # some INSertion in sample relative to reference
            if myoptions.debug: print("some INSertion in sample")
            _qseq = find_insertions_in_sseq_and_drop_them_from_qseq(qseq, sseq)
        else:
            _qseq = qseq

        if checksum:
            output_str = ">%s.%s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s%s" % (myid, checksum, sstart, send, sstrand, evalue, bitscore, score, length, pident, nident, mismatch, positive, gapopen, gaps, ppos, sseq, os.linesep, _qseq)
        else:
            output_str = ">%s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s%s" % (myid, sstart, send, sstrand, evalue, bitscore, score, length, pident, nident, mismatch, positive, gapopen, gaps, ppos, sseq, os.linesep, _qseq)
            
        sys.stdout.write(output_str + '\n')
