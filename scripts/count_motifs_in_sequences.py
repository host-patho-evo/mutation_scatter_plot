#! /usr/bin/env python3

# Copyright (c) 2025 Charles University in Prague - First Faculty of Medicine
#
# Permission is hereby granted, free of charge, to any person obtaining
# a copy of this software and associated documentation files (the
# "Software"), to deal in the Software without restriction, including
# without limitation the rights to use, copy, modify, merge, publish,
# distribute, sublicense, and/or sell copies of the Software, and to
# permit persons to whom the Software is furnished to do so, subject to
# the following conditions:
#
# The above copyright notice and this permission notice shall be
# included in all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
# EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
# MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
# IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY
# CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
# TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
# SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
#
# Except as contained in this notice, the name of Charles University in
# Prague - First Faculty of Medicine shall not be used in advertising or
# otherwise to promote the sale, use or other dealings in this Software
# without prior written authorization from Charles University in Prague
# - First Faculty of Medicine.

"""Cut out a region of interest out of 2-line FASTA alignment ALN file padded with dashes
for gaps and count incidence of those words. You can highlight multiple words in the resulting
UNIX terminal on each line by REGEXP using grep and multiple patterns:

count_motifs_in_sequences.py --infile=$somefile --start-position=498 --end-position=501 | \
    grep -e QPTN -e QPTY -e HPTN -e HPTY -e HPTT -e RPTY

"""

version = "202504011000"

from optparse import OptionParser
import subprocess, shlex
import os, sys, io
import Bio
from Bio import SeqIO
import re

myparser = OptionParser(version="%s version %s" % ('%prog', version))
myparser.add_option("--infilename", action="store", type="string", dest="infilename", default='stdin.fasta',
    help="Input FASTA/Q file path.")
myparser.add_option("--infile-format", action="store", type="string", dest="infile_format", default='fasta',
    help="Input FASTA/Q file format [default: fasta]. Outfile  is always fasta.")
myparser.add_option("--motif", action="store", type="string", dest="motif", default='',
    help="Protein or nucleotide motif to count. No REGEXPs allowed (yet)")
myparser.add_option("--start-position", action="store", type="int", dest="startpos", default=0,
    help="Exact start position of the query")
myparser.add_option("--end-position", action="store", type="int", dest="endpos", default=0,
    help="Exact end position of the query")
myparser.add_option("--debug", action="store", type="int", dest="debug", default=0,
    help="Set debug to some value")
(myoptions, myargs) = myparser.parse_args()

from collections import defaultdict, Counter, OrderedDict

if not myoptions.infilename:
    raise RuntimeError("Please provide input filename via --infilename")
elif not myoptions.infilename.endswith('.counts.fasta') and not myoptions.infilename.endswith('.counts.fasta'):
    sys.stderr.write("Warning: Ideally input filename must end with .prot.counts.fasta or .nuc.counts.fasta and have FASTA id as >12345x after unification by 'sort | uniq -c | sort -nr' beforehand.%s" % os.linesep)
elif not os.path.exists(myoptions.infilename):
    raise RuntimeError("File %s does not exist, please check --infilename" % myoptions.infilename)

_both_approaches = False

# _counter = defaultdict
# sequence_dict = defaultdict(lambda: [None, 0])  # {sequence: [header, count]}
_motif_counter = Counter()
_counter = 0
_total = 0
_motif = myoptions.motif.upper()
_engine = re.compile(_motif)
_motif_length = len(_motif)

for seq_record in SeqIO.parse(myoptions.infilename, "fasta"):
    # lets assume the file has FASTA entries like:
    # >12345x
    try:
        _count = int(seq_record.id.replace('x',''))
    except ValueError:
        # probably has a different file format, do not die and count just as a single match
        _count = 1
    _sequence = str(seq_record.seq).upper()
    _substring = _sequence[myoptions.startpos - 1:myoptions.endpos]
    _left_side = _sequence[:myoptions.endpos].replace('-','')
    _right_side = _sequence[myoptions.startpos - 1:].replace('-','')
    if len(_left_side) < _motif_length or len(_right_side) < _motif_length:
        # do not add to the total count cases when the motif is spanning either end of the sequence
        pass
    else:
        if _substring:
            _motif_counter[_substring] += _count # increase a counter per the respective motif
        _total += _count

if _motif not in _motif_counter.keys():
    # prevent KeyError if the motif was not found at all
    _motif_counter[_motif] = 0

_motif_counter_ordered = OrderedDict(_motif_counter.most_common())
_percentages = OrderedDict()
for _key, _value in _motif_counter_ordered.items():
    if _total:
        _percentages[_key] = round((_value / _total) * 100, 2) # round to two decimal positions
    else:
        _percentages[_key] = 0

if _motif not in _percentages.keys():
    # prevent KeyError if the motif was not found at all
    _percentages[_motif] = 0

print("Motif: %s %s %s %s %s %s" % ('#filename', 'motif', 'motif_count', 'total', 'frequency', 'collection'))
print("Motif: %s %s %s %s %s %s" % (myoptions.infilename, _motif, _motif_counter[_motif], _total, '{0:.2f}'.format(_percentages[_motif]), str(_percentages)))
