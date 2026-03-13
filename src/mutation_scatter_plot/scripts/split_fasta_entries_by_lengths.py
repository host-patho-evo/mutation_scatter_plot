#! /usr/bin/env python3

"""Split FASTA entries into separate files by sequence length.

Kick out sequences not being exactly full_length nucleotides long and place
them into a separate file.  The input file contains unique sequences sorted
in a top-down ordering with FASTA ID containing an incidence count, followed
by a 'x.' and then by an SHA256 checksum. It happens there are typically
singleton sequences with incidence 1x and it is impossible to track back the
original sequence before read trimming, etc. Therefore, the sequence checksum
is very handy.

Usage: split_fasta_entries_by_lengths
           --infile=spikenuc1207.native2ascii.no_junk.clean.counts.fasta
           --outfile-prefix=spikenuc1207.native2ascii.no_junk.clean.counts
           --full-length=3822
"""

import os
from optparse import OptionParser

from Bio import SeqIO

VERSION = "0.3"


def build_option_parser():
    """Build and return the command-line option parser."""
    myparser = OptionParser(version=f"%prog version {VERSION}")
    myparser.add_option(
        "--infile", action="store", type="string", dest="infile", default='',
        help="Input FASTA/Q file path.",
    )
    myparser.add_option(
        "--outfile-prefix", action="store", type="string",
        dest="outfile_prefix", default='',
        help="Output file path prefix.",
    )
    myparser.add_option(
        "--full-length", action="store", type="int", dest="full_length",
        default=0,
        help="Full length required for perfect alignment [0]",
    )
    myparser.add_option(
        "--format", action="store", type="string", dest="format",
        default="fasta",
        help="Input file format.",
    )
    myparser.add_option(
        "--debug", action="store", type="int", dest="debug", default=0,
        help="Set debug to some value",
    )
    return myparser


def main():
    """Entry point: parse arguments and split FASTA file by sequence length."""
    myparser = build_option_parser()
    myoptions, _ = myparser.parse_args()

    if not myoptions.infile:
        raise RuntimeError("Please specify --infile")
    if not os.path.exists(myoptions.infile):
        raise RuntimeError(f"Input file not found: {myoptions.infile}")
    if os.path.getsize(myoptions.infile) == 0:
        raise RuntimeError(f"Input file is empty: {myoptions.infile}")

    _exact_length_name = (
        myoptions.outfile_prefix
        + ".exactly_" + str(myoptions.full_length) + ".fasta"
    )
    _shorter_length_name = (
        myoptions.outfile_prefix
        + ".shorter_" + str(myoptions.full_length) + ".fasta"
    )
    _longer_length_name = (
        myoptions.outfile_prefix
        + ".longer_" + str(myoptions.full_length) + ".fasta"
    )

    _exact_length_cnt   = 0
    _shorter_length_cnt = 0
    _longer_length_cnt  = 0

    with open(_exact_length_name,   'w', encoding='utf-8') as _exact_length, \
         open(_shorter_length_name, 'w', encoding='utf-8') as _shorter_length, \
         open(_longer_length_name,  'w', encoding='utf-8') as _longer_length:

        for _record in SeqIO.parse(myoptions.infile, myoptions.format):
            if len(_record.seq) < myoptions.full_length:
                SeqIO.write(_record, _shorter_length, "fasta-2line")
                _shorter_length_cnt += 1
            elif len(_record.seq) > myoptions.full_length:
                SeqIO.write(_record, _longer_length, "fasta-2line")
                _longer_length_cnt += 1
            else:
                SeqIO.write(_record, _exact_length, "fasta-2line")
                _exact_length_cnt += 1

    print(f"Info: Wrote {_exact_length_cnt} entries into {_exact_length_name}")
    print(f"Info: Wrote {_shorter_length_cnt} entries into {_shorter_length_name}")
    print(f"Info: Wrote {_longer_length_cnt} entries into {_longer_length_name}")


if __name__ == "__main__":
    main()
