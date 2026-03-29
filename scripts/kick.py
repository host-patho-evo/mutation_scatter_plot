#!/usr/bin/env python3
"""Kick out sequences that are not exactly the required length.

Reads a deduplicated FASTA file (NNNNx.sha256 IDs) and splits records into
three output files based on whether each sequence's padded length is:
  - exactly *full_length* nucleotides  → {prefix}.exactly_N.fasta
  - shorter than *full_length*         → {prefix}.shorter_N.fasta
  - longer than *full_length*          → {prefix}.longer_N.fasta

Usage:
    kick.py --infile=spikenuc1207.native2ascii.no_junk.clean.counts.fasta \\
            --outfile-prefix=spikenuc1207.native2ascii.no_junk.clean.counts \\
            --full-length=3822
"""

import os
import sys
from optparse import OptionParser

from Bio import SeqIO

VERSION = 202603292130

myparser = OptionParser(version=f'%prog version {VERSION}')
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
    myparser.error(f"File does not exist: {myoptions.infile}")
if not myoptions.full_length:
    myparser.error("--full-length is required")

_exact_length_name   = f"{myoptions.outfile_prefix}.exactly_{myoptions.full_length}.fasta"
_shorter_length_name = f"{myoptions.outfile_prefix}.shorter_{myoptions.full_length}.fasta"
_longer_length_name  = f"{myoptions.outfile_prefix}.longer_{myoptions.full_length}.fasta"

# ── timestamp-aware output guard ──────────────────────────────────────────────
# Make-style: skip if all outputs exist and are newer than the input.
_outputs     = [_exact_length_name, _shorter_length_name, _longer_length_name]
_input_mtime = os.path.getmtime(myoptions.infile)
_all_exist   = all(os.path.exists(p) for p in _outputs)

if _all_exist:
    _min_out_mtime = min(os.path.getmtime(p) for p in _outputs)
    if _min_out_mtime > _input_mtime and not myoptions.overwrite:
        print(f"Info: all outputs are up-to-date (newer than {myoptions.infile}), skipping.",
              file=sys.stderr)
        sys.exit(0)
    elif not myoptions.overwrite:
        _stale = next(p for p in _outputs if os.path.getmtime(p) <= _input_mtime)
        raise RuntimeError(
            f"Output is stale (older than input): {_stale}\n"
            "Use --overwrite to regenerate."
        )
else:
    for _outpath in _outputs:
        if os.path.exists(_outpath) and not myoptions.overwrite:
            raise RuntimeError(
                f"Output file already exists: {_outpath}\n"
                "Use --overwrite to replace it."
            )

_exact_length_cnt   = 0
_shorter_length_cnt = 0
_longer_length_cnt  = 0

with (open(_exact_length_name,   'w', encoding='utf-8') as _exact_length,
      open(_shorter_length_name, 'w', encoding='utf-8') as _shorter_length,
      open(_longer_length_name,  'w', encoding='utf-8') as _longer_length):
    for _record in SeqIO.parse(myoptions.infile, myoptions.format):
        _record_length = len(_record.seq)
        if myoptions.debug:
            print(f"Info: Record {_record.id} has length {_record_length} "
                  f"of padded sequence {_record.seq}")
        if _record_length < myoptions.full_length:
            SeqIO.write(_record, _shorter_length, "fasta-2line")
            _shorter_length_cnt += 1
        elif _record_length > myoptions.full_length:
            SeqIO.write(_record, _longer_length, "fasta-2line")
            _longer_length_cnt += 1
        else:
            SeqIO.write(_record, _exact_length, "fasta-2line")
            _exact_length_cnt += 1

print(f"Info: Wrote {_exact_length_cnt} entries into {_exact_length_name}")
print(f"Info: Wrote {_shorter_length_cnt} entries into {_shorter_length_name}")
print(f"Info: Wrote {_longer_length_cnt} entries into {_longer_length_name}")
