#! /usr/bin/env python3

# Kick out sequences not being exactly 3822 nucleotides long and place them into a separate file.
# The input file contains unique sequences sorted in a top-down ordering with FASTA ID containing
# an incidence count, followed by a 'x.' and then by an SHA256 checksum. It happens there are
# typically singleton sequences with incidence 1x and it is impossible to track back the original
# sequence before read trimming, etc. Therefore, the sequence checksum is very handy.
#
# Usage: kick.py --infile=spikenuc1207.native2ascii.no_junk.clean.counts.fasta --outfile-prefix=spikenuc1207.native2ascii.no_junk.clean.counts --full-length=3822

import os, sys 
import subprocess
import shutil
import io, hashlib
#import tempfile
from Bio import SeqIO

from optparse import OptionParser
import gzip

VERSION = 202512161435

myparser = OptionParser(version="%s version %s" % ('%prog', VERSION))
myparser.add_option("--infile", action="store", type="string", dest="infile", default='',
    help="Input FASTA/Q file path.")
myparser.add_option("--outfile-prefix", action="store", type="string", dest="outfile_prefix", default='',
    help="Output file path. [infile.replace('.fast[aq]', '_with_lineages.fasta')]")
myparser.add_option("--full-length", action="store", type="int", dest="full_length", default=0,
    help="Full length required for perfect alignment [0]")
myparser.add_option("--format", action="store", type="string", dest="format", default="fasta",
    help="Input file format.")
myparser.add_option("--debug", action="store", type="int", dest="debug", default=0,
    help="Set debug to some value")
(myoptions, myargs) = myparser.parse_args()

_exact_length_name = myoptions.outfile_prefix + ".exactly_" + str(myoptions.full_length) + ".fasta"
_shorter_length_name = myoptions.outfile_prefix + ".shorter_" + str(myoptions.full_length) + ".fasta"
_longer_length_name = myoptions.outfile_prefix + ".longer_" + str(myoptions.full_length) + ".fasta"

_exact_length = open(_exact_length_name, 'w')
_shorter_length = open(_shorter_length_name, 'w')
_longer_length = open(_longer_length_name, 'w')

_exact_length_cnt = 0
_shorter_length_cnt = 0
_longer_length_cnt = 0

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

print("Info: Wrote %d entries into %s" % (_exact_length_cnt, _exact_length_name))
print("Info: Wrote %d entries into %s" % (_shorter_length_cnt, _shorter_length_name))
print("Info: Wrote %d entries into %s" % (_longer_length_cnt, _longer_length_name))

    
