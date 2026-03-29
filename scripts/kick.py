#! /usr/bin/env python3

# Kick out sequences not being exactly for example 3822 nucleotides long and place them into a separate file.
# The input file contains unique sequences sorted in a top-down ordering with FASTA ID containing
# an incidence count, followed by a 'x.' and then by an SHA256 checksum. It happens there are
# typically singleton sequences with incidence 1x and it is impossible to track back the original
# sequence before read trimming, etc. Therefore, the sequence checksum is very handy.
#
# The input sequences are padded and can start with many '-'.
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

VERSION = 202603292130

myparser = OptionParser(version="%s version %s" % ('%prog', VERSION))
myparser.add_option("--infile", action="store", type="string", dest="infile", default='',
    help="Input FASTA/Q file path.")
myparser.add_option("--outfile-prefix", action="store", type="string", dest="outfile_prefix", default='',
    help="Output file path prefix. Suffixes .exactly_N.fasta, .shorter_N.fasta, .longer_N.fasta are appended.")
myparser.add_option("--full-length", action="store", type="int", dest="full_length", default=0,
    help="Full length required for perfect alignment [0]")
myparser.add_option("--format", action="store", type="string", dest="format", default="fasta",
    help="Input file format. Output format is fasta-2line")
myparser.add_option("--overwrite", action="store_true", dest="overwrite", default=False,
    help="Overwrite output files if they already exist. By default the script skips if "
         "all outputs are up-to-date and errors if any output is stale.")
myparser.add_option("--debug", action="store", type="int", dest="debug", default=0,
    help="Set debug to some value")
(myoptions, myargs) = myparser.parse_args()

if not myoptions.infile:
    myparser.error("--infile is required")
if not os.path.exists(myoptions.infile):
    myparser.error("File does not exist: %s" % myoptions.infile)
if not myoptions.full_length:
    myparser.error("--full-length is required")

_exact_length_name   = myoptions.outfile_prefix + ".exactly_"  + str(myoptions.full_length) + ".fasta"
_shorter_length_name = myoptions.outfile_prefix + ".shorter_" + str(myoptions.full_length) + ".fasta"
_longer_length_name  = myoptions.outfile_prefix + ".longer_"  + str(myoptions.full_length) + ".fasta"

# ── timestamp-aware output guard ──────────────────────────────────────────────
# Make-style: skip if all outputs exist and are newer than the input.
_outputs     = [_exact_length_name, _shorter_length_name, _longer_length_name]
_input_mtime = os.path.getmtime(myoptions.infile)
_all_exist   = all(os.path.exists(p) for p in _outputs)

if _all_exist:
    _min_out_mtime = min(os.path.getmtime(p) for p in _outputs)
    if _min_out_mtime > _input_mtime and not myoptions.overwrite:
        print("Info: all outputs are up-to-date (newer than %s), skipping." % myoptions.infile,
              file=sys.stderr)
        sys.exit(0)
    elif not myoptions.overwrite:
        _stale = next(p for p in _outputs if os.path.getmtime(p) <= _input_mtime)
        raise RuntimeError(
            "Output is stale (older than input): %s\n"
            "Use --overwrite to regenerate." % _stale
        )
else:
    for _outpath in _outputs:
        if os.path.exists(_outpath) and not myoptions.overwrite:
            raise RuntimeError(
                "Output file already exists: %s\n"
                "Use --overwrite to replace it." % _outpath
            )

_exact_length   = open(_exact_length_name,   'w')
_shorter_length = open(_shorter_length_name, 'w')
_longer_length  = open(_longer_length_name,  'w')

_exact_length_cnt   = 0
_shorter_length_cnt = 0
_longer_length_cnt  = 0

for _record in SeqIO.parse(myoptions.infile, myoptions.format):
    _record_length = len(_record.seq)
    if myoptions.debug:
        print("Info: Record %s has length %d of padded sequence %s" % (_record.id, _record_length, _record.seq))
    if _record_length < myoptions.full_length:
        SeqIO.write(_record, _shorter_length, "fasta-2line")
        _shorter_length_cnt += 1
    elif _record_length > myoptions.full_length:
        SeqIO.write(_record, _longer_length, "fasta-2line")
        _longer_length_cnt += 1
    else:
        SeqIO.write(_record, _exact_length, "fasta-2line")
        _exact_length_cnt += 1

print("Info: Wrote %d entries into %s" % (_exact_length_cnt,   _exact_length_name))
print("Info: Wrote %d entries into %s" % (_shorter_length_cnt, _shorter_length_name))
print("Info: Wrote %d entries into %s" % (_longer_length_cnt,  _longer_length_name))
