#! /usr/bin/env python3

# https://creativecommons.org/licenses/by/4.0/legalcode.en

"""Cut out a region of interest out of 2-line FASTA alignment ALN file padded with dashes
for gaps and count incidence of those words. You can highlight multiple words in the resulting
UNIX terminal on each line by REGEXP using grep and multiple patterns:

count_motifs_in_sequences --infile=$somefile --start-position=498 --end-position=501 | \
    grep -e QPTN -e QPTY -e HPTN -e HPTY -e HPTT -e RPTY

"""

import os
import sys
from collections import Counter, OrderedDict
from optparse import OptionParser

from Bio import SeqIO

VERSION = "0.3"


def build_option_parser():
    """Build and return the command-line option parser."""
    myparser = OptionParser(version=f"%prog version {VERSION}")
    myparser.add_option(
        "--infilename", action="store", type="string", dest="infilename",
        default='stdin.fasta',
        help="Input FASTA/Q file path.",
    )
    myparser.add_option(
        "--infile-format", action="store", type="string",
        dest="infile_format", default='fasta',
        help="Input FASTA/Q file format [default: fasta]. Outfile is always fasta.",
    )
    myparser.add_option(
        "--motif", action="store", type="string", dest="motif", default='',
        help="Protein or nucleotide motif to count. No REGEXPs allowed (yet)",
    )
    myparser.add_option(
        "--start-position", action="store", type="int", dest="startpos",
        default=0,
        help="Exact start position of the query",
    )
    myparser.add_option(
        "--end-position", action="store", type="int", dest="endpos", default=0,
        help="Exact end position of the query",
    )
    myparser.add_option(
        "--debug", action="store", type="int", dest="debug", default=0,
        help="Set debug to some value",
    )
    return myparser


def main():  # pylint: disable=too-many-locals,too-many-branches
    """Entry point: count motif occurrences in a FASTA alignment file."""
    myparser = build_option_parser()
    myoptions, _ = myparser.parse_args()

    if not myoptions.infilename:
        raise RuntimeError("Please provide input filename via --infilename")
    if not myoptions.infilename.endswith('.counts.fasta'):
        sys.stderr.write(
            "Warning: Ideally input filename must end with .prot.counts.fasta or "
            ".nuc.counts.fasta and have FASTA id as >12345x after unification by "
            f"'sort | uniq -c | sort -nr' beforehand.{os.linesep}"
        )
    if not os.path.exists(myoptions.infilename):
        raise RuntimeError(
            f"File {myoptions.infilename} does not exist, please check --infilename"
        )

    _motif_counter = Counter()
    _total = 0
    _motif = myoptions.motif.upper()
    _motif_length = len(_motif)

    for seq_record in SeqIO.parse(myoptions.infilename, "fasta"):
        try:
            _count = int(seq_record.id.replace('x', ''))
        except ValueError:
            _count = 1
        _sequence = str(seq_record.seq).upper()
        _substring = _sequence[myoptions.startpos - 1:myoptions.endpos]
        _left_side = _sequence[:myoptions.endpos].replace('-', '')
        _right_side = _sequence[myoptions.startpos - 1:].replace('-', '')
        if len(_left_side) < _motif_length or len(_right_side) < _motif_length:
            pass
        else:
            if _substring:
                _motif_counter[_substring] += _count
            _total += _count

    if _motif not in _motif_counter.keys():
        _motif_counter[_motif] = 0

    _motif_counter_ordered = OrderedDict(_motif_counter.most_common())
    _percentages = OrderedDict()
    for _key, _value in _motif_counter_ordered.items():
        if _total:
            _percentages[_key] = round((_value / _total) * 100, 2)
        else:
            _percentages[_key] = 0

    if _motif not in _percentages.keys():
        _percentages[_motif] = 0

    print(
        f"Motif: {'#filename'} {'motif'} {'motif_count'} {'total'} "
        f"{'frequency'} {'collection'}"
    )
    print(
        f"Motif: {myoptions.infilename} {_motif} {_motif_counter[_motif]} "
        f"{_total} {_percentages[_motif]:.2f} {_percentages}"
    )


if __name__ == "__main__":
    main()
