#! /usr/bin/env python3

"""
This script can parse a multiple-sequence alignment file and compare the nucleotides
or aminoacid residues of every row with a reference sequence. 

Parse an external reference file and compare to it a multi-FASTA file with sequences
aligned to it. Report only nucleotides or aminoacid residues changed when compared to the
reference. User must provide both input files in same alphabet, that is either both in DNA
or in protein sequence type.

The functions discard leading and trailing X aminoacids from the aligned sequence and output
a dot (unchanged nucleotide or aminoacid) to silence sequencing noise.


Usage example: alignment2dots.py --reference-infile=MN908947.3_S.fasta --alignment-file=E6/I358F-4th-round-of-sort__E6.I34.WTref.scores_above_84.fastp.amplicons.clean.counts.fasta --print-fasta-ids --aln_start=1288 --aln_stop=1584 --relative_threshold=0.001

find . -name *.WTref.scores_above_84.fastp.amplicons.clean.counts.fasta | while read f; do echo $f; alignment2dots.py --reference-infile=MN908947.3_S.fasta --alignment-file=$f --print-fasta-ids --aln_start=1288 --aln_stop=1584 --threshold=30; done
"""

VERSION = 202509181335

import os
import sys
from optparse import OptionParser
from Bio import SeqIO

myparser = OptionParser()
myparser.add_option("--reference-infile", action="store", type="string", dest="reference_infilename", default=None, metavar="FILE",
    help="FASTA formatted input file with reference padded sequence or not")
myparser.add_option("--alignment-file", action="store", type="string", dest="alignment_file", default=None, metavar="FILE",
    help="FASTA formatted input file with padded sequence or not")
myparser.add_option("--aminoacids", action="store_true", dest="aminoacids", default=False,
    help="FASTA formatted input files are protein sequences instead of DNA [default: DNA]")
myparser.add_option("--outfilename", action="store", type="string", dest="outfilename", default=None, metavar="FILE",
    help="Output filename. If the filename ends with .tsv the output lines will be split using TABs and original FASTA ID will be on the same line with the split sequence")
myparser.add_option("--aln_start", action="store", type="int", dest="aln_start", default=1,
    help="First nucleotide of the ORF region of interest to be sliced out from the input sequences")
myparser.add_option("--aln_stop", action="store", type="int", dest="aln_stop", default=0,
    help="Last nucleotide of the last codon of interest to be sliced out from the input sequences")
myparser.add_option("--print-fasta-ids", action="store_true", dest="print_ids", default=False,
    help="The FASTA file ID should be printed for each aligned entry")
myparser.add_option("--threshold", action="store", type="int", dest="threshold", default=0,
    help="Set the minimum absolute count of different characters in sequence to be output. Set this to 1 or higher if you want to see at aleast 2 aa residues being changed and in turn, do not want to see just dashes in some cases for synonymous changes [default: 0].")
myparser.add_option("--relative_threshold", action="store", type="float", dest="relative_threshold", default=0,
    help="Set the minimum relative incidence threshold of the whole sequence hash. Maybe you want something like 0.001 [default: 0].")
myparser.add_option("--top_n", action="store", type="int", dest="top_n", default=0,
    help="Write only first N entries containing some difference [default: all]")
myparser.add_option("--debug", action="store", type="int", dest="debug", default=0,
    help="Set debug level to some real number")

(myoptions, myargs) = myparser.parse_args()


def replace_unchanged_chars(reference_infilename, alignment_file, outfilename, aln_start, aln_stop, top_n=myoptions.top_n, aminoacids=myoptions.aminoacids, split_using_tabs=False):
    """Parse an external reference file and compare to it a multi-FASTA file with sequences
    aligned to it. Report only nucleotides or aminoacid residues changed when compared to the
    reference. User must provide both input files in same alphabet, that is either both in DNA
    or in protein sequence type.

    The functions discard leading and trailing X aminoacids from the aligned sequence and output
    a dot (unchanged nucleotide or aminoacid) to silence sequencing noise.
    """

    _refseq_record_handle = SeqIO.parse(reference_infilename,'fasta')
    _output_counter = 0
    if aminoacids:
        _unknown_char = 'X'
    else:
        _unknown_char = 'N'
    for _refseq_record in _refseq_record_handle:
        _refseq_sequence = _refseq_record.seq.upper()
        print("Info: reference sequence %s '%s' is %s characters long" % (_refseq_record.id, _refseq_sequence, len(_refseq_sequence)))
        _aligned_records_handle = SeqIO.parse(alignment_file, 'fasta')
        _highest_count = 0
        for _seq_record in _aligned_records_handle:
            if not top_n or _output_counter < top_n:
                _count = 0
                if not _highest_count:
                    if 'x.' in _seq_record.id:
                        _highest_count = int(_seq_record.id.split('x.')[0])
                        _count = _highest_count
                    elif 'x' in _seq_record.id:
                        _highest_count = int(_seq_record.id.split('x')[0])
                        _count = _highest_count
                    else:
                        _highest_count = 0
                        _count = 0
                else:
                    if 'x.' in _seq_record.id:
                        _count = int(_seq_record.id.split('x.')[0])
                    elif 'x' in _seq_record.id:
                        _count = int(_seq_record.id.split('x')[0])
                if not myoptions.relative_threshold or float(_count) / float(_highest_count) > myoptions.relative_threshold:
                    _aligned_sequence = _seq_record.seq.upper()
                    _output = ''
                    for _i in range(aln_start, aln_stop or len(_refseq_sequence) - 1):
                        _char1 = _refseq_sequence[_i]
                        try:
                            _char2 = _aligned_sequence[_i]
                        except IndexError:
                            raise IndexError("Cannot slice aligned sequence %s '%s' at position %s" % (_seq_record.id, _aligned_sequence, _i))
                        if (_i == 1 and _char2 == _unknown_char) or (_i > 0 and (_char2 == _unknown_char and (not _aligned_sequence[:_i-1].replace('-','') or not _aligned_sequence[_i+1:].replace('-','') ))):
                            # output '.' instead of 'X' or 'N' when this is at the very beginning of the aligned segment or at its very end
                            _output += '.'
                        elif _char1 != _char2:
                            _output += _char2
                        else:
                            _output += '.'
                    if not myoptions.threshold or len(_output.replace('.','').replace('-','')) > myoptions.threshold:
                        if split_using_tabs:
                            # split the _output using TABs
                            _tabbed_output = '\t'.join(list(_output))
                            if myoptions.print_ids:
                                outfilename.write(f"{_seq_record.id}\t{_tabbed_output}{os.linesep}")
                            else:
                                outfilename.write(f"{_tabbed_output}{os.linesep}")
                        else:
                            if myoptions.print_ids:
                                outfilename.write(f"{_seq_record.id}{os.linesep}")
                            outfilename.write(f"{_output}{os.linesep}")
                        _output_counter += 1
                elif myoptions.debug:
                    print("Debug: Below threshold: %f is < %f" % (float(_count) / float(_highest_count), myoptions.relative_threshold))


def main():
    if myoptions.aln_start:
        _aln_start = myoptions.aln_start - 1
    else:
        _aln_start = 0
    if myoptions.aln_stop:
        _aln_stop = myoptions.aln_stop - 1
    else:
        _aln_stop = 0

    if myoptions.outfilename:
        if not os.path.exists(myoptions.outfilename):
            _outfileh = open(myoptions.outfilename,'w')
        else:
            raise RuntimeError("File %s already exists, exitting." % myoptions.outfilename)
    else:
        _outfileh = sys.stdout
    if myoptions.outfilename and myoptions.outfilename.endswith('.tsv'):
        replace_unchanged_chars(myoptions.reference_infilename, myoptions.alignment_file, _outfileh, _aln_start, _aln_stop, split_using_tabs=True)
    else:
        replace_unchanged_chars(myoptions.reference_infilename, myoptions.alignment_file, _outfileh, _aln_start, _aln_stop, split_using_tabs=False)
    _outfileh.close()

if __name__ == "__main__":
    main()

# vim:ts=4:sw=4:expandtab:smartindent
