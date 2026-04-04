#! /usr/bin/env python

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

#>1x.d2b2b8089265daaac94d3c3553ddf45d00841b13da817b2f2f703cd834d05437 1288 1584 plus
#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ACTGGCTGCGTTATAGCTTGGAATTCTAACAAGCTTGATTCTAAGGTTGGTGGTAATTATAATTACCTGTATAGATTGTTTAGGAAGTCTAATCTCAAACCTTTTGAGAGAGATATTTCAACTGAAATCTATCAGGCCGGTAACAAACCTTGTAATGGTGTTGCAGGTTTTAATTGTTACTTTCCTTTACGATCATATGGTTTCCGACCCACTTATGGTGTTGGTCACCAACCATACAGA---GTAGTACTTTCATTTGAACTTCTACATGCACCAGGAACTGTTTGTGGACCTAAA--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#
# >1x.4a00c0c4ef25b507488b88a4a3ff5e9ab5c5abcb523664fffb12da2eaae85309 1 2066 plus 3363 3729 2078 95.813 1991 51 1991 6 36 95.81
# ATGTTTGTTTTTCTTGTTTTATTGCCACTAGTCTCTAGTCAGTGTGTCATGCCGCTGTTTAATCTTMTAACTACAACTCAAT---------CATACACTAATTCTTTCACACGTGGTGTTTATTACCCTGACAAAGTTTTCAGATCCTCAGTTTTACATTTAACTCAGGACTTGTTCTTACCTTTCTTTTCCAATGTTACTTGGTTCCATGCTA------TCTCTGGGACCAATGGTACTAAGAGGTTTGATAACCCTGTCCTACCATTTAATGATGGTGTTTATTTTGCTTCCACTGAGAAGTCTAACATAATAAGAGGCTGGATTTTTGGTACTACTTTAGATTCGAAGACCCAGTCCCTACTTATTGTTAATAACGCTACTAATGTTTTTATTAAAGTCTGTGAATTTCAATTTTGTAATGATCCATTTTTGGATGTTTA---CCACAAAAACAACAAAAGTTGGATGGAAAGTGAGTCAGGAGTTTATTCTAGTGCGAATAATTGCACTTTTGAATATGTCTCTCAGCCTTTTCTTATGGACCTTGAAGGAAAACAGGGTAATTTCAAAAATCTTAGGGAATTTGTGTTTAAGAATATTGATGGTTATTTTAAAATATATTCTAAGCACACGCCTATTA---TAGGGCGTGATTTCCCTCAGGGTTTTTCGGCTTTAGAACCATTGGTAGATTTGCCAATAGGTATTAACATCACTAGGTTTCAAACTTTACTTGCTTTAAATAGAAGTTATTTGACTCCTGGTGATTCTTCTTCAGGTTGGACAGCTGGTGCTGCAGATTATTATGTGGGTTATCTTCAACCTAGGACTTTTCTATTAAAATATAATGAAAATGGAACCATTACAGATGCTGTAGACTGTGCACTTGACCCTCTCTCAGAAACAAAGTGTACGTTGAAATCCTTCACTGTAGAAAAAGGAATCTATCAAACTTCTAACTTTAGAGTCCAACCAACAGAATCTATTGTTAGATTTCCTAATGTTACAAACTTGTGCCCTTTTCATGAAGTTTTTAACGCCACCAGATTTGCATCTGTTTATGCTTGGAACAGGACGAGAATCAGCAACTGTGTTGCTGATTATTCTGTCCTATATAATTTCGCACCATTTTTCGCTTTTAAGTGTTATGGAGTGTCTCCTACTAAATTAAATGATCTCTGCTTTACTAATGTCTATGCAGATTCATTTGTAATTAAAGGTAATGAAGTCAGCCAAATCGCTCCAGGGCAAACTGGAAATATTGCTGATTATAATTATAAATTACCAGATGATTTTACAGGCTGCGTTATAGCTTGGAATTCTAACAAGCTTGATTCTAAGCATAGTGGTAATTATGATTACTGGTATAGATCGTTTAGGAAGTCTAAACTCAAACCTTTTGAGAGAGATATTTCAACTGAAATCTATCAGGCCGGTAACAAACCTTGTAAAGGT---AAAGGTCCTAATTGTTACTTTCCTTTACAATCATATGGTTTCCGACCCACTTATGGTGTTGGTCACCAACCATACAGAGTAGTAGTACTTTCTTTTGAACTTCTACATGCACCAGCAACTGTTTGTGGACCTAAAAAGTCTACTAATTTGGTTAAAAACAAATGTGTCAATTTCAACTTCAATGGTTTAACAGGCACAGGTGTTCTTACTAAGTCTAACAAAAAGTTTCTGCCTTTCCAACAATTTGGCAGAGACATTGTTGACACTACTGATGCTGTCCGTGATCCACAGACACTTGAGATTCTTGACATTACACCATGTTCTTTTGGTGGTGTCAGTGTTATAACACCAGGAACAAATACTTCTAACCAGGTTGCTGTTCTTTATCAGGGTGTTAACTGCACAGAAGTCTCTGTTGCTATTCATGCAGATCAACTTACTCCTACTTGGCGTGTTTATTCTACAGGTTCTAATGTTTTTCAAACACGTGCAGGCTGTTTAATAGGGGCTGAATATGTCAACAACTCATATGAGTGTGACATACCCATTGGTGCAGGTATATGCGCTAGTTATCAGACTCAGACTAAGTCTCGTCGGCGGGCACGTAGTGTAGCTAG

# BUG: It breaks on the following sequence with 12nt insertion and results in a shorter line only 3810 nt long instead of 3822 nt.
# >536x.a3f44694f8b950e98fb28d7b79e21e7ab4fbc8d6bb8b7e8f16f588184dd48b06 1 2352 plus 3877 4299 2364 96.320 2277 51 2277 6 36 96.32
# ATGTTTGTTTTTCTTGTTTTATTGCCACTAGTCTCTAGTCAGTGTGTCATGCCGCTGTTTAATCTTATAACTACAACTCAA---------TCATACACTAATTCTTTCACACGTGGTGTTTATTACCCTGACAAAGTTTTCAGATCCTCAGTTTTACATTTAACTCAGGACTTGTTCTTACCTTTCTTTTCCAATGTTACTTGGTTCCATGCTATC------TCTGGGACCAATGGTACTAAGAGGTTTGATAACCCTGTCCTACCATTTAATGATGGTGTTTATTTTGCTTCCACTGAGAAGTCTAACATAATAAGAGGCTGGATTTTTGGTACTACTTTAGATTCGAAGACCCAGTCCCTACTTATTGTTAATAACGCTACTAATGTTTTTATTAAAGTCTGTGAATTTCAATTTTGTAATGATCCATTTTTGGATGTTTAC---CACAAAAACAACAAAAGTTGGATGGAAAGTGAGTCAGGAGTTTATTCTAGTGCGAATAATTGCACTTTTGAATATGTCTCTCAGCCTTTTCTTATGGACCTTGAAGGAAAACAGGGTAATTTCAAAAATCTTAGGGAATTTGTGTTTAAGAATATTGATGGTTATTTTAAAATATATTCTAAGCACACGCCTATT---ATAGGGCGTGATTTCCCTCAGGGTTTTTCGGCTTTAGAACCATTGGTAGATTTGCCAATAGGTATTAACATCACTAGGTTTCAAACTTTACTTGCTTTAAATAGAAGTTATTTGACTCCTGGTGATTCTTCTTCAGGTTGGACAGCTGGTGCTGCAGATTATTATGTGGGTTATCTTCAACCTAGGACTTTTCTATTAAAATATAATGAAAATGGAACCATTACAGATGCTGTAGACTGTGCACTTGACCCTCTCTCAGAAACAAAGTGTACGTTGAAATCCTTCACTGTAGAAAAAGGAATCTATCAAACTTCTAACTTTAGAGTCCAACCAACAGAATCTATTGTTAGATTTCCTAATGTTACAAACTTGTGCCCTTTTCATGAAGTTTTTAACGCCACCAGATTTGCATCTGTTTATGCTTGGAACAGGACGAGAATCAGCAACTGTGTTGCTGATTATTCTGTCCTATATAATTTCGCACCATTTTTCGCTTTTAAGTGTTATGGAGTGTCTCCTACTAAATTAAATGATCTCTGCTTTACTAATGTCTATGCAGATTCATTTGTAATTAAAGGTAATGAAGTCAGCCAAATCGCTCCAGGGCAAACTGGAAATATTGCTGATTATAATTATAAATTACCAGATGATTTTACAGGCTGCGTTATAGCTTGGAATTCTAACAAGCTTGATTCTAAGCATAGTGGTAATTATGATTACTGGTATAGATCGTTTAGGAAGTCTAAACTCAAACCTTTTGAGAGAGATATTTCAACTGAAATCTATCAGGCCGGTAACAAACCTTGTAAAGGT---AAAGGTCCTAATTGTTACTTTCCTTTACAATCATATGGTTTCCGACCCACTTATGGTGTTGGTCACCAACCATACAGAGTAGTAGTACTTTCTTTTGAACTTCTACATGCACCAGCAACTGTTTGTGGACCTAAAAAGTCTACTAATTTGGTTAAAAACAAATGTGTCAATTTCAACTTCAATGGTTTAACAGGCACAGGTGTTCTTACTAAGTCTAACAAAAAGTTTCTGCCTTTCCAACAATTTGGCAGAGACATTGTTGACACTACTGATGCTGTCCGTGATCCACAGACACTTGAGATTCTTGACATTACACCATGTTCTTTTGGTGGTGTCAGTGTTATAACACCAGGAACAAATACTTCTAACCAGGTTGCTGTTCTTTATCAGGGTGTTAACTGCACAGAAGTCTCTGTTGCTATTCATGCAGATCAACTTACTCCTACTTGGCGTGTTTATTCTACAGGTTCTAATGTTTTTCAAACACGTGCAGGCTGTTTAATAGGGGCTGAATATGTCAACAACTCATATGAGTGTGACATACCCATTGGTGCAGGTATATGCGCTAGTTATCAGACTCAGACTAAGTCTCGTCGGCGGGCACGTAGTGTAGCTAGTCAATCCATCATTGCCTACACTATGTCACTTGGTGCAGAAAATTCAGTTGCTTACTCTAATAACTCTATTGCCATACCCACAAATTTTACTATTAGTGTTACCACAGAAATTCTACCAGTGTCTATGACCAAGACATCAGTAGATTGTACAATGTACATTTGTGGTGATTCAACTGAATGCAGCAATCTTTTGTTGCAATATGGCAGTTTTTGTACACAATTAAAACGTGCTTTAACTGGAATAGCTGTTGAACAAGACAAAAACACCCAAGAAGTTTTTGCACAA------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

"""
This script can be used to flip FASTA/FASTQ input stream into a plus stranded FASTA
and optionally, can prepend and append dashes to make the resulting FASTA aligned.
For this it requires a reference sequence.

One can use UNIX awk command to parse the blastn outpu stream and print it out as
a multiFASTA stream, directly consumed by this utility.

It also adjusts the start and stop alignment positions in the FASTA header if they
meet expectations and are e.g. from blastn TSV output.

Example:
gzip -dc HWL7GDRX3.fastp.fastq.gz | fastq_to_fasta | blastn -task blastn -reward 2 -max_hsps 1 -num_alignments 1 -word_size 4 -num_threads 64 -dust no -evalue 1e-10 -outfmt "6 qacc sstart send sstrand evalue bitscore score length pident nident mismatch positive gapopen gaps ppos qseq sseq" -db "$fileref" -query - 2>/dev/null | awk '{print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17}' | drop_erroneous_insertions.py --drop-misinsertions --infile=- outfile=- | reversecomplement_reads_on_minus.py --reference="$fileref" --infile=- --min_start=1 --max_stop=1833 --x-after-count > "$prefix".aln

Check the lengths of the entries with:
grep -v '^>' "$prefix".aln | while read line; do echo "$line" | wc -m ; done | sort | uniq -c | sort -nr

Check for lines longer or shorter than expected:
grep -v '^>' "$prefix".aln | awk '{ if ( length!=3822 ) print }' |  grep -v 3822 | sort | uniq | sort | while read s; do grep -B 1 -F -f - "$prefix".aln; done

"""

import sys
import os
from optparse import OptionParser

import textwrap
from Bio import SeqIO

VERSION = "202604042205"

_COMPLEMENT_TAB = str.maketrans(
    'ACGTRYSWKMBDHVNacgtryswkmbdhvn',
    'TGCAYRSWMKVHDBNtgcayrswmkvhdbn'
)

def string_reverse_complement(seq: str) -> str:
    """Natively calculate sequence reverse complement at C string speed without BioPython allocation overhead."""
    return str(seq).translate(_COMPLEMENT_TAB)[::-1]

myparser = OptionParser(version="%s version %s" % ('%prog', VERSION))
myparser.add_option("--infile", action="store", type="string", dest="infile", default='',
    help="Input FASTA file with word plus or minus as the second word of the header, use minus (dash) for stdin")
myparser.add_option("--outfile", action="store", type="string", dest="outfile", default='',
    help="Output filename with words plus and minus removed from FASTA header lines, use minus (dash) for stdout")
myparser.add_option("--infileformat", action="store", type="string", dest="infileformat", default='fasta',
    help="Input file format (Default: fasta)")
myparser.add_option("--outfileformat", action="store", type="string", dest="outfileformat", default='fasta-2line',
    help="Output file format (Default: fasta-2line)")
myparser.add_option("--reference", action="store", type="string", dest="referencefilename", default='',
    help="Reference filename to infer length of the reference sequence from. It is necessary to calculate lengths of paddings around the query sequence to fit width of the multiple sequence alignment. (Default: '')")
myparser.add_option("--min_start", action="store", type="int", dest="min_start", default=0,
    help="Minimum start position of the amplicon region to be included in output. The length of the sequence if filled with padding dashes as necessary.")
myparser.add_option("--max_stop", action="store", type="int", dest="max_stop", default=0,
    help="Specify width of the alignment to be output, including paddings. Maximum stop position of the amplicon region to be included in output. The length of the sequence if filled with padding dashes as necessary. If not specified, no trailing dashes will be added after the HSP match.")
myparser.add_option("--aln_start", action="store", type="int", dest="aln_start", default=0,
    help="Start coordinate of the FASTA sequence when aligned to a reference")
myparser.add_option("--aln_stop", action="store", type="int", dest="aln_stop", default=0,
    help="End coordinate of the FASTA sequence when aligned to a reference")
myparser.add_option("--min_count", action="store", type="int", dest="min_count", default=0,
    help="Minimum read count")
myparser.add_option("--x-after-count", action="store_true", dest="x_after_count", default=False,
    help="The FASTA file ID contains the count value followed by lowercase 'x', then followed by a dot and then by a checksum")
myparser.add_option("--debug", action="store", type="int", dest="debug", default=0,
    help="Set debug level")
(myoptions, myargs) = myparser.parse_args()


def shorten_sequence(sequence, min_start, number_of_leading_dashes, number_of_trailing_dashes, max_stop, aln_stop, qseq, aln_start_qseq, aln_stop_qseq):
    """
    A00808:1538:HWL7GDRX3:2:2101:21377:18975 1141 1543
    A00808:1538:HWL7GDRX3:2:2101:28230:1658 1141 1371

    $ blastn -task blastn -reward 2 -max_hsps 1 -num_alignments 1 -word_size 4 -num_threads 4 -dust no -evalue 1e-5 -outfmt '6 qacc sstart send sstrand qseq sseq' -db MN908947.3_S.fasta -query LH00211:37:222VFLLT1:1:1101:4920:2218:ACTATTTTC.fasta
    LH00211:37:222VFLLT1:1:1101:4920:2218:ACTATTTTC	1315	1545	plus	AACAATCTTGATTCTAAGGTTGGTGGTAATTACAATTACCTGTTTAGATTGTTTAGGAAGTCTAATCTCAAACCTTTTGAGAAAGATATTTCAACTGAAATCTATCAGGCCGGCAGCACACCTTGTAATGATGTTGAAGGTTTTAATTGTTATTTTCCTCTACAATCATATGGTTTCCATCGTACCAATGGTGTTGGTACGATGGAAACCATACAGAGTAGTAGTACTTTCTTTT	AACAATCTTGATTCTAAGGTTGGTGGTAATTATAATTACCTGTATAGATTGTTTAGGAAGTCTAATCTCAAACCTTTTGAGAGAGATATTTCAACTGAAATCTATCAGGCCGGTAGCACACCTTGTAATGGTGTTGAAGGTTTTAATTGTTACTTTCCTTTACAATCATATGGTTTCCAACCCACTAATGGTGTTGGT----TACCAACCATACAGAGTAGTAGTACTTTCTTTT
    $

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


    >MN908947.3_S_protein  3822nt 1274codons
    Length=3822
    
     Score = 376 bits (416),  Expect = 4e-108
     Identities = 220/228 (96%), Gaps = 0/228 (0%)
     Strand=Plus/Plus
    
    Query  1     AATCTTGATTCTAAGGTTGGTGGTAATTACAATTACCTGTTTAGATTGTTTAGGAAGTCT  60
                 ||||||||||||||||||||||||||||| |||||||||| |||||||||||||||||||
    Sbjct  1318  AATCTTGATTCTAAGGTTGGTGGTAATTATAATTACCTGTATAGATTGTTTAGGAAGTCT  1377
    
    Query  61    AATCTCAAACCTTTTGAGAGAGATATTTCAACTGAAATCTTTCAGGCCGGCAGCACACCT  120
                 |||||||||||||||||||||||||||||||||||||||| ||||||||| |||||||||
    Sbjct  1378  AATCTCAAACCTTTTGAGAGAGATATTTCAACTGAAATCTATCAGGCCGGTAGCACACCT  1437
    
    Query  121   TGTAATGGTGTTGAAGGTTTTAATTGTTATTTTCCTGTACAATCATATGGTTTCCAACCT  180
                 ||||||||||||||||||||||||||||| |||||| |||||||||||||||||||||| 
    Sbjct  1438  TGTAATGGTGTTGAAGGTTTTAATTGTTACTTTCCTTTACAATCATATGGTTTCCAACCC  1497
    
    Query  181   ACCAATGGTGTTGGTTACCAACCATACAGAGTAGTAGTACTTTCTTTT  228
                 || |||||||||||||||||||||||||||||||||||||||||||||
    Sbjct  1498  ACTAATGGTGTTGGTTACCAACCATACAGAGTAGTAGTACTTTCTTTT  1545


    >MN908947.3_S_protein  3822nt 1274codons
    Length=3822
    
     Score = 366 bits (405),  Expect = 7e-105
     Identities = 218/228 (96%), Gaps = 4/228 (2%)
     Strand=Plus/Plus
    
    Query  1     AATCTTGATTCTAAGGTTGGTGGTAATTACAATTACCTGTTTAGATTGTTTAGGAAGTCT  60
                 ||||||||||||||||||||||||||||| |||||||||| |||||||||||||||||||
    Sbjct  1318  AATCTTGATTCTAAGGTTGGTGGTAATTATAATTACCTGTATAGATTGTTTAGGAAGTCT  1377
    
    Query  61    AATCTCAAACCTTTTGAGAAAGATATTTCAACTGAAATCTTTCAGGCCGGCAGCACACCT  120
                 ||||||||||||||||||| |||||||||||||||||||| ||||||||| |||||||||
    Sbjct  1378  AATCTCAAACCTTTTGAGAGAGATATTTCAACTGAAATCTATCAGGCCGGTAGCACACCT  1437
    
    Query  121   TGTAATGGTGTTGAAGGTTTTAATTGTTACTTTCCTTTACAATCATATGGTTTCCAACCC  180
                 ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
    Sbjct  1438  TGTAATGGTGTTGAAGGTTTTAATTGTTACTTTCCTTTACAATCATATGGTTTCCAACCC  1497
    
    Query  181   ACTAATGGTGTTGG----AAACCATACAGAGTAGTAGTACTTTCTTTT  224
                 ||||||||||||||     |||||||||||||||||||||||||||||
    Sbjct  1498  ACTAATGGTGTTGGTTACCAACCATACAGAGTAGTAGTACTTTCTTTT  1545

    $


$ blastn -task blastn -reward 2 -max_hsps 1 -word_size 4 -num_threads 64 -dust no -evalue 1e-5 -db MN908947.3_S.fasta -query spikenuc0719.unidecode.native2ascii.unidecode.no_junk.2010-12.counts.fasta

Query= 1x.94eb9ddd80265473c73ea4f1477bf485d2c4015928e803eef13cbd1ab87e3668

Length=3756
                                                                      Score     E
Sequences producing significant alignments:                          (Bits)  Value

MN908947.3_S_protein 3822nt 1274codons                                3664    0.0  


>MN908947.3_S_protein  3822nt 1274codons
Length=3822

 Score = 3664 bits (4063),  Expect = 0.0
 Identities = 2553/2899 (88%), Gaps = 35/2899 (1%)
 Strand=Plus/Plus

Query  889   AAAGGCATTTACCAAACTTCTAATTTTAGGACTTCTCCTACTACACAAG---TTGTTAGG  945
             ||||| || || ||||||||||| ||||| | | |  | |  ||| ||    ||||||| 
Sbjct  928   AAAGGAATCTATCAAACTTCTAACTTTAG-AGTCCAACCA--ACAGAATCTATTGTTAGA  984

Query  946   TTTCCTAATATTACAAATTTATGCCCTTTTGGTGAAGTTTTTAACGCCACCACTTTCGCT  1005
             ||||||||||||||||| || |||||||||||||||||||||||||||||||  || || 
Sbjct  985   TTTCCTAATATTACAAACTTGTGCCCTTTTGGTGAAGTTTTTAACGCCACCAGATTTGCA  1044

Query  1006  TCAGTTTATGCATGGAACAGAAGAAGAATCAGCAATTGTGTTGCAGATTATTCTGTACTA  1065
             || |||||||| |||||||| |  ||||||||||| |||||||| ||||||||||| |||
Sbjct  1045  TCTGTTTATGCTTGGAACAGGAAGAGAATCAGCAACTGTGTTGCTGATTATTCTGTCCTA  1104

Query  1066  TATAACACAACCTCATTTTCAACTTTTAAATGTTATGGGGTTTCACCCACTAAATTAAAT  1125
             |||||  |  | |||||||| |||||||| |||||||| || || || ||||||||||||
Sbjct  1105  TATAATTCCGCATCATTTTCCACTTTTAAGTGTTATGGAGTGTCTCCTACTAAATTAAAT  1164

Query  1126  GATCTCTGTTTTACTAATGTTTATGCAGACTCGTTTGTTGTTAGAGGTGATGAAGTTAGG  1185
             |||||||| ||||||||||| |||||||| || |||||  |||||||||||||||| || 
Sbjct  1165  GATCTCTGCTTTACTAATGTCTATGCAGATTCATTTGTAATTAGAGGTGATGAAGTCAGA  1224

Query  1186  CAAATAGCTCCAGGTCAGACTGGCAAAATTGCTGACTATAATTATAAACTCCCAGATGAT  1245
             ||||| |||||||| || ||||| || |||||||| |||||||||||| | |||||||||
Sbjct  1225  CAAATCGCTCCAGGGCAAACTGGAAAGATTGCTGATTATAATTATAAATTACCAGATGAT  1284

Query  1246  TTTATGGGTTGTGTAATAGCATGGAATTCTATAAGTTTAGATGCT------GGTGGTTCT  1299
             ||||  || || || ||||| ||||||||||  | | | ||| ||      ||||||  |
Sbjct  1285  TTTACAGGCTGCGTTATAGCTTGGAATTCTAACAATCTTGATTCTAAGGTTGGTGGTAAT  1344

Query  1300  TATTATTA------TAGACTTTTTAGAAAGTCTGTTCTTAAGCCTTTTGAAAGAGATATA  1353
             ||| ||||      |||| | ||||| ||||||  ||| || |||||||| |||||||| 
Sbjct  1345  TATAATTACCTGTATAGATTGTTTAGGAAGTCTAATCTCAAACCTTTTGAGAGAGATATT  1404

Query  1354  TCTACTCAACTTTATCAAGCAGGTGATAAACCTTGCTCA----GTTGAGGGTCCTGATTG  1409
             || ||| || | ||||| || |||   | |||||| | |    ||||| |||  | ||||
Sbjct  1405  TCAACTGAAATCTATCAGGCCGGTAGCACACCTTG-TAATGGTGTTGAAGGTTTTAATTG  1463

Query  1410  TTACTACCCGTTGCAGTCCTATTATTTTCAGAGCACTAATGGTGTTGGTTACCAACCTTA  1469
             |||||  || || || || |||  ||| ||   |||||||||||||||||||||||| ||
Sbjct  1464  TTACTTTCCTTTACAATCATATGGTTTCCAACCCACTAATGGTGTTGGTTACCAACCATA  1523

Query  1470  TAGGGTAGTAGTGCTATCTTTTGAGCTTCTAAATGCACCAGCAACCGTCTGTGGGCCTAA  1529
              || |||||||| || |||||||| |||||| ||||||||||||| || ||||| |||||
Sbjct  1524  CAGAGTAGTAGTACTTTCTTTTGAACTTCTACATGCACCAGCAACTGTTTGTGGACCTAA  1583

Query  1530  AAAATCCACCCACTTGGTTGTAAACAAATGTGTCAATTTTAACTTTAATGGTTTAACTGG  1589
             ||| || ||  | ||||||  |||||||||||||||||| ||||| ||||||||||| ||
Sbjct  1584  AAAGTCTACTAATTTGGTTAAAAACAAATGTGTCAATTTCAACTTCAATGGTTTAACAGG  1643

Query  1590  TACAGGTGTTCTTACGTCTTCTACAAAGAAATTCTTGCCTTTCCAGCAATTTGGTAGAGA  1649
              ||||||||||||||    ||||  || || ||  |||||||||| |||||||| |||||
Sbjct  1644  CACAGGTGTTCTTACTGAGTCTAACAAAAAGTTTCTGCCTTTCCAACAATTTGGCAGAGA  1703

Query  1650  TGTAGCTGATACTACAAATGCTGTTCGAGACCCACAGACACTTGAGGTGTTAGACATAAC  1709
               | ||||| |||||  ||||||| || || ||||||||||||||| |  | ||||| ||
Sbjct  1704  CATTGCTGACACTACTGATGCTGTCCGTGATCCACAGACACTTGAGATTCTTGACATTAC  1763

Query  1710  ACCATGTTCTTTTGGTGGTGTTAGTGTCATAACGCCTGGAACTAATACCTCTACTCAAGT  1769
             ||||||||||||||||||||| ||||| ||||| || ||||| ||||| ||||  || ||
Sbjct  1764  ACCATGTTCTTTTGGTGGTGTCAGTGTTATAACACCAGGAACAAATACTTCTAACCAGGT  1823

Query  1770  AGCTGTTCTTTATCAAGATGTCAATTGTACTGACGTCCCGTCTGCAATCCATGCTGACCA  1829
              |||||||||||||| ||||| || || || || |||||   ||| || ||||| || ||
Sbjct  1824  TGCTGTTCTTTATCAGGATGTTAACTGCACAGAAGTCCCTGTTGCTATTCATGCAGATCA  1883

Query  1830  GCTTTCATCCACTTGGCGTGTTTATTCTACAGGTCCTAATGTTTTTCAAACACGCGCAGG  1889
              ||| |  | |||||||||||||||||||||||| ||||||||||||||||||| |||||
Sbjct  1884  ACTTACTCCTACTTGGCGTGTTTATTCTACAGGTTCTAATGTTTTTCAAACACGTGCAGG  1943

Query  1890  CTGTTTAATAGGGGCTGAACATGTCAACAATTCATATGATTGTGACATACCCATTGGTGC  1949
             |||||||||||||||||||||||||||||| |||||||| ||||||||||||||||||||
Sbjct  1944  CTGTTTAATAGGGGCTGAACATGTCAACAACTCATATGAGTGTGACATACCCATTGGTGC  2003

Query  1950  AGGTATATGCGCCAGTTACCAGACTCAAACTAATTC------------ACGTAGTGTAAC  1997
             |||||||||||| ||||| |||||||| ||||||||            |||||||||| |
Sbjct  2004  AGGTATATGCGCTAGTTATCAGACTCAGACTAATTCTCCTCGGCGGGCACGTAGTGTAGC  2063

Query  1998  CAGTCAATCCATTATTGCTTACACTATGTCACTTGGTGCAGAAAATTCAGTTGCTTATTC  2057
              ||||||||||| ||||| |||||||||||||||||||||||||||||||||||||| ||
Sbjct  2064  TAGTCAATCCATCATTGCCTACACTATGTCACTTGGTGCAGAAAATTCAGTTGCTTACTC  2123

Query  2058  TAATAACTCTATTGCCATACCTACAAATTTCACCATTAGTGTCACCACAGAAATTCTACC  2117
             ||||||||||||||||||||| |||||||| || |||||||| |||||||||||||||||
Sbjct  2124  TAATAACTCTATTGCCATACCCACAAATTTTACTATTAGTGTTACCACAGAAATTCTACC  2183

Query  2118  TGTGTCTATGACCAAGACGTCGGTAGATTGTACAATGTACATTTGTGGTGACTCAACAGA  2177
              ||||||||||||||||| || ||||||||||||||||||||||||||||| ||||| ||
Sbjct  2184  AGTGTCTATGACCAAGACATCAGTAGATTGTACAATGTACATTTGTGGTGATTCAACTGA  2243

Query  2178  GTGCAGCAATCTTTTGTTGCAATATGGCAGCTTTTGCTCACAACTAAATCGTGCTTTAAC  2237
              ||||||||||||||||||||||||||||| |||||  ||||| |||| |||||||||||
Sbjct  2244  ATGCAGCAATCTTTTGTTGCAATATGGCAGTTTTTGTACACAATTAAACCGTGCTTTAAC  2303

Query  2238  TGGAATAGCTGTTGAACAAGATAAAAACACTCAAGAGGTTTTTGCACAAGTCAAACAACT  2297
             ||||||||||||||||||||| |||||||| ||||| ||||||||||||||||||||| |
Sbjct  2304  TGGAATAGCTGTTGAACAAGACAAAAACACCCAAGAAGTTTTTGCACAAGTCAAACAAAT  2363

Query  2298  TTATAAAACACCACCAATCAAAGACTTTGGTGGTTTCAATTTTTCACAAATATTACCAGA  2357
             ||| |||||||||||||| ||||| ||||||||||| |||||||||||||||||||||||
Sbjct  2364  TTACAAAACACCACCAATTAAAGATTTTGGTGGTTTTAATTTTTCACAAATATTACCAGA  2423

Query  2358  TCCATCAAAACCAAGCAAGCGGTCATTTATTGAGGACTTGCTTTTCAATAAAGTGACACT  2417
             ||||||||||||||||||| ||||||||||||| ||  | |||||||| |||||||||||
Sbjct  2424  TCCATCAAAACCAAGCAAGAGGTCATTTATTGAAGATCTACTTTTCAACAAAGTGACACT  2483

Query  2418  TGCTGATGCTGGCTTCATCAAACAATATGGTGACTGCCTTGGTGATATTGCTGCTAGAGA  2477
             ||| ||||||||||||||||||||||||||||| ||||||||||||||||||||||||||
Sbjct  2484  TGCAGATGCTGGCTTCATCAAACAATATGGTGATTGCCTTGGTGATATTGCTGCTAGAGA  2543

Query  2478  TCTCATTTGTGCACAAAAGTTTAACGGCCTTACTGTTCTGCCACCTTTACTCACAGATGA  2537
              |||||||||||||||||||||||||||||||||||| |||||||||| |||||||||||
Sbjct  2544  CCTCATTTGTGCACAAAAGTTTAACGGCCTTACTGTTTTGCCACCTTTGCTCACAGATGA  2603

Query  2538  AATGATTGCTCAATACACCTCTGCACTATTAGCAGGTACAATAACTTCAGGTTGGACCTT  2597
             |||||||||||||||||| |||||||| ||||| |||||||| ||||| |||||||||||
Sbjct  2604  AATGATTGCTCAATACACTTCTGCACTGTTAGCGGGTACAATCACTTCTGGTTGGACCTT  2663

Query  2598  TGGTGCAGGTGCTGCGTTACAAATACCATTTGCAATGCAAATGGCTTATAGGTTTAATGG  2657
             ||||||||||||||| ||||||||||||||||| ||||||||||||||||||||||||||
Sbjct  2664  TGGTGCAGGTGCTGCATTACAAATACCATTTGCTATGCAAATGGCTTATAGGTTTAATGG  2723

Query  2658  TATTGGAGTTACACAGAATGTACTTTATGAGAATCAGAAATTGATTGCCAATCAATTCAA  2717
             ||||||||||||||||||||| || |||||||| || |||||||||||||| ||||| ||
Sbjct  2724  TATTGGAGTTACACAGAATGTTCTCTATGAGAACCAAAAATTGATTGCCAACCAATTTAA  2783

Query  2718  CAGTGCCATTGGCAAAATCCAAGACTCGCTTTCTTCAACTGCTAGTGCACTTGGAAAACT  2777
              ||||| ||||||||||| |||||||| |||||||| || || |||||||||||||||||
Sbjct  2784  TAGTGCTATTGGCAAAATTCAAGACTCACTTTCTTCCACAGCAAGTGCACTTGGAAAACT  2843

Query  2778  TCAAGATGTTGTCAACCAAAATGCACAGGCTCTAAATACACTTGTAAAACAACTTAGTTC  2837
             ||||||||| ||||||||||||||||| ||| |||| || ||||| ||||||||||| ||
Sbjct  2844  TCAAGATGTGGTCAACCAAAATGCACAAGCTTTAAACACGCTTGTTAAACAACTTAGCTC  2903

Query  2838  TAACTTTGGAGCTATTTCAAGTGTGCTAAATGATATTCTTTCACGTCTCGACAAAGTTGA  2897
              || ||||| || |||||||||||  |||||||||| ||||||||||| |||||||||||
Sbjct  2904  CAATTTTGGTGCAATTTCAAGTGTTTTAAATGATATCCTTTCACGTCTTGACAAAGTTGA  2963

Query  2898  GGCTGAAGTTCAAATTGATAGGCTGATCACAGGTAGACTCCAAAGTCTGCAGACTTATGT  2957
             ||||||||| |||||||||||| |||||||||| ||||| |||||| ||||||| |||||
Sbjct  2964  GGCTGAAGTGCAAATTGATAGGTTGATCACAGGCAGACTTCAAAGTTTGCAGACATATGT  3023

Query  2958  GACTCAACAACTAATCAGAGCCGCAGAGATCAGAGCTTCTGCTAATCTTGCTGCAACTAA  3017
             |||||||||| |||| ||||| ||||| |||||||||||||||||||||||||| |||||
Sbjct  3024  GACTCAACAATTAATTAGAGCTGCAGAAATCAGAGCTTCTGCTAATCTTGCTGCTACTAA  3083

Query  3018  AATGTCAGAGTGTGTTCTTGGACAATCTAAAAGAGTTGACTTTTGTGGAAAAGGCTATCA  3077
             ||||||||||||||| ||||||||||| ||||||||||| ||||||||||| ||||||||
Sbjct  3084  AATGTCAGAGTGTGTACTTGGACAATCAAAAAGAGTTGATTTTTGTGGAAAGGGCTATCA  3143

Query  3078  TTTAATGTCTTTTCCTCAGTCAGCACCTCATGGTGTGGTCTTTTTGCATGTGACCTACGT  3137
             | | ||||| || ||||||||||||||||||||||| ||||| ||||||||||| || ||
Sbjct  3144  TCTTATGTCCTTCCCTCAGTCAGCACCTCATGGTGTAGTCTTCTTGCATGTGACTTATGT  3203

Query  3138  CCCTGCACAACAAAAGAACTTCACAACTGCTCCTGCCATTTGTCATGATGGAAAAGCACA  3197
             |||||||||| |||||||||||||||||||||||||||||||||||||||||||||||||
Sbjct  3204  CCCTGCACAAGAAAAGAACTTCACAACTGCTCCTGCCATTTGTCATGATGGAAAAGCACA  3263

Query  3198  CTTCCCTCGCGAAGGTGTCTTTGTGTCAAATGGTACGCATTGGTTTGTAACACAAAGGAA  3257
             ||| ||||| |||||||||||||| |||||||| || || ||||||||||||||||||||
Sbjct  3264  CTTTCCTCGTGAAGGTGTCTTTGTTTCAAATGGCACACACTGGTTTGTAACACAAAGGAA  3323

Query  3258  TTTCTATGAACCACAAATTATTACTACAGACAACACATTTGTCTCTGGCCACTGTGATGT  3317
             ||| |||||||||||||| ||||||||||||||||||||||| |||||  ||||||||||
Sbjct  3324  TTTTTATGAACCACAAATCATTACTACAGACAACACATTTGTGTCTGGTAACTGTGATGT  3383

Query  3318  TGTAATTGGAATTGTCAACAACACAGTTTATGATCCTTTGCAACCAGAACTTGATTCATT  3377
             |||||| |||||||||||||||||||||||||||||||||||||| ||| | || |||||
Sbjct  3384  TGTAATAGGAATTGTCAACAACACAGTTTATGATCCTTTGCAACCTGAATTAGACTCATT  3443

Query  3378  CAAAGAGGAGCTAGACAAGTATTTTAAAAATCATACATCGCCAGATGTTGATTTAGGCGA  3437
             ||| |||||| |||| || |||||||| ||||||||||| ||||||||||||||||| ||
Sbjct  3444  CAAGGAGGAGTTAGATAAATATTTTAAGAATCATACATCACCAGATGTTGATTTAGGTGA  3503

Query  3438  CATCTCTGGCATTAATGCTTCAGTTGTCAATATTCAAAAAGAAATTGACCGCCTCAATGA  3497
             ||||||||||||||||||||||||||| || |||||||||||||||||||||||||||||
Sbjct  3504  CATCTCTGGCATTAATGCTTCAGTTGTAAACATTCAAAAAGAAATTGACCGCCTCAATGA  3563

Query  3498  GGTTGCCAAGAATTTAAATGAATCTCTCATCGATCTCCAAGAACTTGGAAAGTATGAGCA  3557
             ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
Sbjct  3564  GGTTGCCAAGAATTTAAATGAATCTCTCATCGATCTCCAAGAACTTGGAAAGTATGAGCA  3623

Query  3558  GTACATAAAATGGCCATGGTACATTTGGCTAGGTTTTATAGCTGGCTTGATTGCCATAGT  3617
             ||| ||||||||||||||||||||||||||||||||||||||||||||||||||||||||
Sbjct  3624  GTATATAAAATGGCCATGGTACATTTGGCTAGGTTTTATAGCTGGCTTGATTGCCATAGT  3683

Query  3618  AATGGTGACAATTATGCTTTGCTGTATGACCAGTTGCTGCAGTTGTCTCAAGGGCTGTTG  3677
             ||||||||||||||||||||||||||||||||||||||| ||||||||||||||||||||
Sbjct  3684  AATGGTGACAATTATGCTTTGCTGTATGACCAGTTGCTGTAGTTGTCTCAAGGGCTGTTG  3743

Query  3678  TTCTTGTGGATCCTGCTGCAAATTTGATGAAGACGACTCTGAGCCAGTGCTCAAAGGAGT  3737
             ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
Sbjct  3744  TTCTTGTGGATCCTGCTGCAAATTTGATGAAGACGACTCTGAGCCAGTGCTCAAAGGAGT  3803

Query  3738  CAAATTACATTACACATAA  3756
             |||||||||||||||||||
Sbjct  3804  CAAATTACATTACACATAA  3822



    $ blastn -task blastn -reward 2 -max_hsps 1 -num_alignments 1 -word_size 4 -num_threads 4 -dust no -evalue 1e-5 -outfmt '6 qacc sstart send sstrand qseq sseq' -db MN908947.3_S.fasta -query LH00211:37:222VFLLT1:1:1101:27573:3724:CAATGCCCG.fasta
    LH00211:37:222VFLLT1:1:1101:27573:3724:CAATGCCCG	1291	1569	plus	GGCTGCGTTATAGCTTGGAATTCTAACAATCTTGATTCTAAGGTTGGTGGTAATTACAATTACCTGTTTAGATTGTTTAGGAAGTCTAATCTCAAACCTTTTGAGAAAGATATTTCAACTGAAATCTTTCAGGCCGCCAGCACACCTTGTAATAATGTTGAAGGTTTTAATTGTTATTTTCCTCTACAATCATATGGTTTCCATCCTTCCAATGGTTGTGGTTACCAACCATACAGAGTAGTAGTACTTTCTTTTGAACTTCTACATGCACCAGCAACT	GGCTGCGTTATAGCTTGGAATTCTAACAATCTTGATTCTAAGGTTGGTGGTAATTATAATTACCTGTATAGATTGTTTAGGAAGTCTAATCTCAAACCTTTTGAGAGAGATATTTCAACTGAAATCTATCAGGCCGGTAGCACACCTTGTAATGGTGTTGAAGGTTTTAATTGTTACTTTCCTTTACAATCATATGGTTTCCAACCCACTAATGGTGTTGGTTACCAACCATACAGAGTAGTAGTACTTTCTTTTGAACTTCTACATGCACCAGCAACT

    Check the length of the paddded sequences and report these with different length:
    find . -name \\*.fastp.amplicons.counts.clean.fasta  | while read f; do echo $f; grep -v '^>' "$f" | awk '{ if ( length!=3822 ) print }' |  grep -v 3822 | sort | uniq | sort | while read s; do grep -B 1 -F -f - "$f"; done; done
    """

    if max_stop or min_start:
        _required_alignment_width = max_stop - min_start + 1
    else:
        _required_alignment_width = 0

    # determine how many nucleotides should should be stripped away from the beginning of the padded query/sample sequence, after figuring out if it is in plus or minus orientation
    if myoptions.debug: print("Info: entered shorten_sequence() with %s" % sequence)
    if aln_start_qseq < aln_stop_qseq:
        # is on plus
        _unmatched_nucleotides_at_the_left_query_end = aln_start_qseq - 1 # length of the leading unmatched query sequence
    else:
        # is on minus
        _unmatched_nucleotides_at_the_left_query_end = len(qseq) - aln_stop_qseq # length of the trailig unmatched query sequence

    if max_stop:
        if max_stop < aln_stop:
            if myoptions.debug: print("Info: Seems we need to drop trailing %d nucleotides at the right end" % number_of_trailing_dashes)
            _sequence_to_return = number_of_leading_dashes * '-' + sequence + number_of_trailing_dashes * '-'
        elif max_stop > aln_stop:
            if myoptions.debug: print("Info1: seems we need extra %d dashes at the right end for extra sequence" % number_of_trailing_dashes)
            _sequence_to_return = number_of_leading_dashes * '-' + sequence + number_of_trailing_dashes * '-'
        else:
            if myoptions.debug: print("Info2: seems we need extra %d dashes at the right end for extra sequence" % number_of_trailing_dashes)
            _sequence_to_return = number_of_leading_dashes * '-' + sequence + number_of_trailing_dashes * '-'
    else:

        if myoptions.debug: print("Info: max_stop was not defined, will not append trailing dashes")
        _sequence_to_return = number_of_leading_dashes * '-' + sequence + number_of_trailing_dashes * '-'

    # slice the sequence before returning to match requested width as specified at the command line
    if max_stop:
        if min_start:
            return _sequence_to_return[min_start:max_stop], min_start, max_stop
        else:
            return _sequence_to_return[:max_stop], min_start, max_stop
    elif min_start:
        return _sequence_to_return[min_start:], min_start, max_stop
    else:
        return _sequence_to_return, min_start, max_stop


def parse_input(infile, reference_sequence, infileformat, min_start=0, max_stop=0):
    """The function assumes nucleotide sequence starts in frame +1. If not, it prepends
    dashes in front so that the number of dashes can be divided by three without remainder.
    """

    if min_start:
        min_start -= 1 # treat real numbering like 0 for 0-based pythonic counting, undo off-by-one error when --min_start=1 is provided

    if max_stop or min_start:
        _required_alignment_width = max_stop - min_start
    else:
        _required_alignment_width = 0

    for record in SeqIO.parse(infile, infileformat):
        _fasta_description_items = record.description.split()
        _min_start = 0
        # >A00808:1538:HWL7GDRX3:2:2101:7491:1125 1542 1141 minus
        _description = record.description.replace(' plus','').replace(' minus','')
        if myoptions.x_after_count:
            _myid = _fasta_description_items[0]
            if 'x.' in _myid:
                _id, _checksum = _myid.split('x.')
            else:
                _id = _myid.rstrip('x')
                _checksum = ''
            try:
                _some_count = int(_id) # cut out the numeric value
            except:
                _some_count = 1
                if myoptions.debug: print("Error: parsing the number failed, _id=%s is probably not an integer" % str(_id))
        else:
            _some_count = 1
            if myoptions.debug: print("Debug: You did not provide --x-after-count, assuming count 1")
        if myoptions.debug: print("Debug: %" % str(_fasta_description_items))
        # >60x.83b5909c7ce41adb2921d29154dfe53e5521ae9f769c33e187db0ecd5f02cfc5 1 3822 plus 0.0 6634 7356 3822 98.482 3764 46 3764 2 12 98.48
        if int(_some_count) >= myoptions.min_count and len(_fasta_description_items) > 1:
            if _fasta_description_items[3] == 'minus':
                # >A00808:1538:HWL7GDRX3:2:2101:7491:1125 1542 1141 minus
                _fasta_description_items = _description.split()
                if len(_fasta_description_items) < 15:
                    raise ValueError(f"Error: Header missing upstream blastn fields. Expected at least 15 items, found {len(_fasta_description_items)}. Ensure awk pipeline passes qseq/sseq columns natively.")
                    
                if len(_fasta_description_items) == 4:
                    _aln_start = int(_fasta_description_items[3])
                    _aln_stop = int(_fasta_description_items[2])
                elif len(_fasta_description_items) == 3:
                    _aln_start = int(_fasta_description_items[2])
                    _aln_stop = int(_fasta_description_items[1])
                elif len(_fasta_description_items) > 4:
                    _aln_start = int(_fasta_description_items[1])
                    _aln_stop = int(_fasta_description_items[2])
                else:
                    _aln_start = 0
                    _aln_stop = 0
                _aln_start_qseq = int(_fasta_description_items[9])
                _aln_stop_qseq = int(_fasta_description_items[10])
                _qseq = record.seq # original query sequence matching on minus
                _sseq = _fasta_description_items[14]
                
                # Perform blazing fast C-level native string reverse complement
                _native_rc_seq = string_reverse_complement(_qseq)

                _number_of_leading_dashes = _aln_stop - 1
                # compare reference_sequence length and end of the minus-oriented alignment to figure out how many leading dashes should be prepended to reverse-complemented query/sample sequence
                if max_stop:
                    _number_of_trailing_dashes = min(len(reference_sequence),max_stop) - _aln_start
                else:
                    _number_of_trailing_dashes = len(reference_sequence) - _aln_start
                if myoptions.debug: print("Debug: BEFORE: id=%s, len=%d padded query minus-oriented sequence:\n%s" % (record.id, len(_native_rc_seq), _native_rc_seq))
                _sequence, _min_start, _max_stop = shorten_sequence(_native_rc_seq, min_start, _number_of_leading_dashes, _number_of_trailing_dashes, max_stop, _aln_stop, _qseq, _aln_start_qseq, _aln_stop_qseq)
                if myoptions.debug: print("Debug: AFTER: id=%s, len=%d padded query minus-oriented sequence:\n%s" % (record.id, len(_sequence), _sequence))
                if _required_alignment_width and len(_sequence) != _required_alignment_width:
                    raise ValueError("Length %d of the padded sequence %s does not match required alignment width %d" % (len(_sequence), _sequence, _required_alignment_width))
                _aln_start += _min_start # increase the aln start position
                _aln_stop += _max_stop # decrease the aln stop position by addition of a negative number
                _reordered_header_items = ' '.join([_fasta_description_items[0], str(_aln_start), str(_aln_stop), 'plus'] + list(_fasta_description_items[4:]))
                
                yield _reordered_header_items, _sequence
            elif _fasta_description_items[3] == 'plus':
                record.description = _description
                _fasta_description_items = _description.split()
                _aln_start_qseq = int(_fasta_description_items[9])
                _aln_stop_qseq = int(_fasta_description_items[10])
                _sseq = _fasta_description_items[14]
                _qseq = record.seq # original query sequence matching on plus
                if len(_fasta_description_items) < 15:
                    raise ValueError(f"Error: Header missing upstream blastn fields. Expected at least 15 items, found {len(_fasta_description_items)}. Ensure awk pipeline passes qseq/sseq columns natively.")
                    
                if len(_fasta_description_items) == 4:
                    _aln_start = int(_fasta_description_items[2])
                    _aln_stop = int(_fasta_description_items[3])
                elif len(_fasta_description_items) == 3:
                    _aln_start = int(_fasta_description_items[1])
                    _aln_stop = int(_fasta_description_items[2])
                elif len(_fasta_description_items) > 4:
                    _aln_start = int(_fasta_description_items[1])
                    _aln_stop = int(_fasta_description_items[2])
                else:
                    _aln_start = 0
                    _aln_stop = 0
                _number_of_leading_dashes = _aln_start - 1
                if max_stop:
                    _number_of_trailing_dashes = min(len(reference_sequence),max_stop) - _aln_stop
                else:
                    _number_of_trailing_dashes = len(reference_sequence) - _aln_stop
                if myoptions.debug: print("Debug: BEFORE: id=%s, len=%d padded query plus-oriented sequence\n%s" % (record.id, len(record.seq), record.seq))
                _sequence, _min_start, _max_stop = shorten_sequence(record.seq, min_start, _number_of_leading_dashes, _number_of_trailing_dashes, max_stop, _aln_stop, _qseq, _aln_start_qseq, _aln_stop_qseq)
                if myoptions.debug: print("Debug: AFTER:  id=%s, len=%d padded query plus-oriented sequence\n%s" % (record.id, len(_sequence), _sequence))
                if _required_alignment_width and len(_sequence) != _required_alignment_width:
                    raise ValueError("Length %d of the padded sequence %s does not match required alignment width %d" % (len(_sequence), _sequence, _required_alignment_width))
                
                _aln_start += _min_start # increase the aln start position
                _aln_stop += _max_stop # decrease the aln stop position by addition of a negative number
                _reordered_header_items = ' '.join([_fasta_description_items[0], str(_aln_start), str(_aln_stop), 'plus'] + list(_fasta_description_items[4:]))
                yield _reordered_header_items, _sequence
            elif myoptions.aln_start or myoptions.aln_stop: # untested code below
                if myoptions.debug: print("Debug: --aln_start or --aln_stop was set")
                record.description = _description
                _fasta_description_items = _description.split()
                if len(_fasta_description_items) < 15:
                    raise ValueError(f"Error: Header missing upstream blastn fields. Expected at least 15 items, found {len(_fasta_description_items)}. Ensure awk pipeline passes qseq/sseq columns natively.")
                    
                _aln_start_qseq = int(_fasta_description_items[9])
                _aln_stop_qseq = int(_fasta_description_items[10])
                _sseq = _fasta_description_items[14]
                _qseq = record.seq # original query sequence
                _aln_start = myoptions.aln_start
                _aln_stop = myoptions.aln_stop
                _number_of_leading_dashes = _aln_start - 1
                if max_stop:
                    _number_of_trailing_dashes = min(len(reference_sequence),max_stop) - _aln_stop
                else:
                    _number_of_trailing_dashes = len(reference_sequence) - _aln_stop
                _sequence, _min_start, _max_stop = shorten_sequence(record.seq, min_start, _number_of_leading_dashes, _number_of_trailing_dashes, max_stop, _aln_stop, _qseq, _aln_start_qseq, _aln_stop_qseq)
                if myoptions.debug: print("Debug: id=%s, len=%d padded query sequence\n%s" % (record.id, len(_sequence), _sequence))
                if _required_alignment_width and len(_sequence) != _required_alignment_width:
                    raise ValueError("Length %d of the padded sequence %s does not match required alignment width %d" % (len(_sequence), _sequence, _required_alignment_width))
                
                _aln_start += _min_start # increase the aln start position
                _aln_stop += _max_stop # decrease the aln stop position by addition of a negative number
                _reordered_header_items = ' '.join([_fasta_description_items[0], str(_aln_start), str(_aln_stop), 'plus'] + list(_fasta_description_items[4:]))
                yield _reordered_header_items, _sequence
            else:
                # >60x.83b5909c7ce41adb2921d29154dfe53e5521ae9f769c33e187db0ecd5f02cfc5 1 3822 plus 0.0 6634 7356 3822 98.482 3764 46 3764 2 12 98.48
                raise RuntimeError("Either provide FASTA input with alignment start and stop coordinates or provide --min_start and --max_stop coordinates which will ensure the result is pre-padded and post-padded with dashes to be of requested width and the entries must specify plus or minus orientation as the last word in the description line.")
        else:
            if myoptions.debug: print("Debug: incidence %s is =< %d or only one item is in FASTA description line" % (_some_count, myoptions.min_count))


if __name__ == "__main__":
    if myoptions.infile == '-':
        _infile = sys.stdin
        if myoptions.debug: print(f"Info: Parsing input file {sys.stdin}")
    elif os.path.exists(myoptions.infile):
        _infile = str(myoptions.infile)
        if myoptions.debug: print(f"Info: Parsing input file {myoptions.infile}")
    else:
        raise RuntimeError("File %s does not exist" % str(myoptions.infile))

    if not myoptions.referencefilename:
        raise RuntimeError("File with reference sequence is required to provide calculate length of the alignment on the fly")
    else:
        _reference_sequence = None
        _reference_iterator = SeqIO.parse(myoptions.referencefilename, "fasta")
        for _reference in _reference_iterator:
            if _reference_sequence:
                raise RuntimeError("There should be only one one reference sequence in the %s file" % myoptions.referencefilename)
            _reference_sequence = _reference.seq
        if myoptions.debug: print("Debug: Reference sequence %s has length %d" % (_reference_sequence, len(_reference_sequence)))
    if not myoptions.outfile or myoptions.outfile == '-':
        _outfile = sys.stdout
#    elif os.path.exists(myoptions.outfile):
#        raise RuntimeError("File %s already exists" % myoptions.outfile)
    else:
        _outfile = myoptions.outfile

    _sequences = parse_input(_infile, _reference_sequence, myoptions.infileformat, min_start=myoptions.min_start, max_stop=myoptions.max_stop)
    
    if myoptions.outfile == '-' or not myoptions.outfile:
        out_f = sys.stdout
    else:
        out_f = open(myoptions.outfile, 'w')
        
    try:
        # Avoid extremely slow Bio.SeqIO.write looping, write directly at C speed bindings
        for header, seq in _sequences:
            out_f.write(f">{header}\n")
            if myoptions.outfileformat == 'fasta':
                # standard wrapped format
                out_f.write('\n'.join(textwrap.wrap(str(seq), width=80)) + '\n')
            else:
                # default fasta-2line format
                out_f.write(f"{seq}\n")
    finally:
        if out_f is not sys.stdout:
            out_f.close()
