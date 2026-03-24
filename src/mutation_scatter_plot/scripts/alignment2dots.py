#! /usr/bin/env python3

"""Parse a multiple-sequence alignment and replace unchanged positions with dots.

Compare the nucleotides or amino acid residues of every row with a reference
sequence. Parse an external reference file and compare to it a multi-FASTA
file with sequences aligned to it. Report only nucleotides or amino acid
residues changed when compared to the reference. User must provide both input
files in the same alphabet, that is either both in DNA or in protein sequence
type.

The functions discard leading and trailing X amino acids from the aligned
sequence and output a dot (unchanged nucleotide or amino acid) to silence
sequencing noise.

Usage example::

    alignment2dots --reference-infile=MN908947.3_S.fasta \\
        --alignment-file=E6/sample.counts.fasta \\
        --print-fasta-ids --aln_start=1288 --aln_stop=1584 \\
        --relative_threshold=0.001

    find . -name '*.counts.fasta' | while read f; do
        alignment2dots --reference-infile=MN908947.3_S.fasta \\
            --alignment-file="$f" --print-fasta-ids \\
            --aln_start=1288 --aln_stop=1584 --threshold=30
    done
"""

import os
import sys
import argparse
from decimal import Decimal
from Bio import SeqIO

VERSION = "0.3"


class NoWrapFormatter(argparse.RawDescriptionHelpFormatter, argparse.ArgumentDefaultsHelpFormatter):
    """Help formatter that does not wrap long lines, preserving URLs."""

    def _split_lines(self, text, width):
        return text.splitlines()


def build_option_parser():
    """Build and return the command-line argument parser."""
    myparser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=NoWrapFormatter,
    )
    myparser.add_argument("--version", action="version", version=f"%(prog)s {VERSION}")
    myparser.add_argument(
        "--reference-infile", action="store", type=str,
        dest="reference_infilename", default=None, metavar="FILE",
        help="FASTA formatted input file with reference padded sequence or not",
    )
    myparser.add_argument(
        "--alignment-file", action="store", type=str,
        dest="alignment_file", default=None, metavar="FILE",
        help="FASTA formatted input file with padded sequence or not",
    )
    myparser.add_argument(
        "--aminoacids", action="store_true", dest="aminoacids", default=False,
        help="FASTA formatted input files are protein sequences instead of DNA"
             " [default: DNA]",
    )
    myparser.add_argument(
        "--outfilename", action="store", type=str, dest="outfilename",
        default=None, metavar="FILE",
        help="Output filename. If the filename ends with .tsv the output lines"
             " will be split using TABs and original FASTA ID will be on the"
             " same line with the split sequence",
    )
    myparser.add_argument(
        "--aln_start", action="store", type=int, dest="aln_start", default=1,
        help="First nucleotide of the ORF region of interest to be sliced out"
             " from the input sequences",
    )
    myparser.add_argument(
        "--aln_stop", action="store", type=int, dest="aln_stop", default=0,
        help="Last nucleotide of the last codon of interest to be sliced out"
             " from the input sequences",
    )
    myparser.add_argument(
        "--print-fasta-ids", action="store_true", dest="print_ids",
        default=False,
        help="The FASTA file ID should be printed for each aligned entry",
    )
    myparser.add_argument(
        "--threshold", action="store", type=int, dest="threshold", default=0,
        help="Set the minimum absolute count of different characters in"
             " sequence to be output. Set this to 1 or higher if you want to"
             " see at least 2 aa residues being changed and in turn, do not"
             " want to see just dashes in some cases for synonymous changes"
             " [default: 0].",
    )
    myparser.add_argument(
        "--relative_threshold", action="store", type=float,
        dest="relative_threshold", default=0,
        help="Set the minimum relative incidence threshold of the whole"
             " sequence hash. Maybe you want something like 0.001 [default: 0].",
    )
    myparser.add_argument(
        "--top_n", action="store", type=int, dest="top_n", default=0,
        help="Write only first N entries containing some difference [default: all]",
    )
    myparser.add_argument(
        "--debug", action="store", type=int, dest="debug", default=0,
        help="Set debug level to some real number",
    )
    return myparser


# pylint: disable=too-many-arguments,too-many-positional-arguments,too-many-locals
# pylint: disable=too-many-branches,too-many-statements,too-many-nested-blocks
def replace_unchanged_chars(
        myoptions, reference_infilename, alignment_file, outfilename,
        aln_start, aln_stop, top_n=0, aminoacids=False,
        split_using_tabs=False):
    """Parse an external reference file and compare to it a multi-FASTA file with sequences
    aligned to it. Report only nucleotides or aminoacid residues changed when compared to the
    reference. User must provide both input files in same alphabet, that is either both in DNA
    or in protein sequence type.

    The functions discard leading and trailing X aminoacids from the aligned sequence and output
    a dot (unchanged nucleotide or aminoacid) to silence sequencing noise.
    """
    _refseq_record_handle = SeqIO.parse(reference_infilename, 'fasta')
    _output_counter = 0
    _unknown_char = 'X' if aminoacids else 'N'

    for _refseq_record in _refseq_record_handle:
        _refseq_sequence = _refseq_record.seq.upper()
        print(
            f"Info: reference sequence {_refseq_record.id}"
            f" '{_refseq_sequence}' is {len(_refseq_sequence)} characters long"
        )
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
                if (not myoptions.relative_threshold
                        or Decimal(_count) / Decimal(_highest_count)
                        > Decimal(str(myoptions.relative_threshold))):
                    _aligned_sequence = _seq_record.seq.upper()
                    _output = ''
                    for _i in range(aln_start, aln_stop or len(_refseq_sequence) - 1):
                        _char1 = _refseq_sequence[_i]
                        try:
                            _char2 = _aligned_sequence[_i]
                        except IndexError as exc:
                            raise IndexError(
                                f"Cannot slice aligned sequence {_seq_record.id}"
                                f" '{_aligned_sequence}' at position {_i}"
                            ) from exc
                        _is_leading_or_trailing_unknown = (
                            _char2 == _unknown_char
                            and (_i == 1
                                 or not _aligned_sequence[:_i - 1].replace('-', '')
                                 or not _aligned_sequence[_i + 1:].replace('-', ''))
                        )  # pylint: disable=too-many-boolean-expressions
                        if _is_leading_or_trailing_unknown:
                            _output += '.'
                        elif _char1 != _char2:
                            _output += _char2
                        else:
                            _output += '.'
                    if (not myoptions.threshold
                            or len(_output.replace('.', '').replace('-', ''))
                            > myoptions.threshold):
                        if split_using_tabs:
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
                    print(
                        f"Debug: Below threshold:"
                        f" {Decimal(_count) / Decimal(_highest_count):.6f}"
                        f" is < {myoptions.relative_threshold}"
                    )


# pylint: enable=too-many-arguments,too-many-positional-arguments,too-many-locals
# pylint: enable=too-many-branches,too-many-statements,too-many-nested-blocks
def main():
    """Entry point: open output file and call replace_unchanged_chars."""
    myparser = build_option_parser()
    myoptions, _ = myparser.parse_args()

    _aln_start = (myoptions.aln_start - 1) if myoptions.aln_start else 0
    _aln_stop  = (myoptions.aln_stop  - 1) if myoptions.aln_stop  else 0

    if not myoptions.reference_infilename:
        raise RuntimeError("Please specify --reference-infile")
    if not os.path.exists(myoptions.reference_infilename):
        raise RuntimeError(
            f"Reference file not found: {myoptions.reference_infilename}"
        )
    if os.path.getsize(myoptions.reference_infilename) == 0:
        raise RuntimeError(
            f"Reference file is empty: {myoptions.reference_infilename}"
        )
    if not myoptions.alignment_file:
        raise RuntimeError("Please specify --alignment-file")
    if not os.path.exists(myoptions.alignment_file):
        raise RuntimeError(
            f"Alignment file not found: {myoptions.alignment_file}"
        )
    if os.path.getsize(myoptions.alignment_file) == 0:
        raise RuntimeError(
            f"Alignment file is empty: {myoptions.alignment_file}"
        )

    _split_tabs = (
        myoptions.outfilename is not None
        and myoptions.outfilename.endswith('.tsv')
    )
    if myoptions.outfilename:
        if os.path.exists(myoptions.outfilename):
            raise RuntimeError(
                f"File {myoptions.outfilename} already exists, exitting."
            )
        with open(myoptions.outfilename, 'w', encoding='utf-8') as _outfileh:
            replace_unchanged_chars(
                myoptions, myoptions.reference_infilename,
                myoptions.alignment_file, _outfileh,
                _aln_start, _aln_stop,
                top_n=myoptions.top_n,
                aminoacids=myoptions.aminoacids,
                split_using_tabs=_split_tabs,
            )
    else:
        replace_unchanged_chars(
            myoptions, myoptions.reference_infilename,
            myoptions.alignment_file, sys.stdout,
            _aln_start, _aln_stop,
            top_n=myoptions.top_n,
            aminoacids=myoptions.aminoacids,
            split_using_tabs=False,
        )


if __name__ == "__main__":
    main()

# vim:ts=4:sw=4:expandtab:smartindent
