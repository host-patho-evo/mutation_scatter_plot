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

import argparse
import os
import sys

from Bio import SeqIO

VERSION = 202603292130

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--infile", required=True, help="Input FASTA/Q file path.")
parser.add_argument("--outfile-prefix", dest="outfile_prefix", default="",
    help="Output file path prefix. Suffixes .exactly_N.fasta etc. are appended.")
parser.add_argument("--full-length", dest="full_length", type=int, required=True,
    help="Full length required for perfect alignment.")
parser.add_argument("--format", dest="format", default="fasta",
    help="Input file format. Output format is fasta-2line.")
parser.add_argument("--overwrite", action="store_true",
    help="Overwrite output files if they already exist.")
parser.add_argument("--debug", type=int, default=0, help="Debug verbosity level.")
parser.add_argument("--version", action="version", version=f"%(prog)s {VERSION}")


def main():
    """Parse arguments and split FASTA into exact/shorter/longer output files."""
    print(
        f"kick.py  version {VERSION}"
        f"  invoked: {' '.join(sys.argv)}",
        file=sys.stderr,
    )
    myoptions = parser.parse_args()

    if not os.path.exists(myoptions.infile):
        parser.error(f"File does not exist: {myoptions.infile}")

    exact_name   = f"{myoptions.outfile_prefix}.exactly_{myoptions.full_length}.fasta"
    shorter_name = f"{myoptions.outfile_prefix}.shorter_{myoptions.full_length}.fasta"
    longer_name  = f"{myoptions.outfile_prefix}.longer_{myoptions.full_length}.fasta"

    outputs     = [exact_name, shorter_name, longer_name]
    input_mtime = os.path.getmtime(myoptions.infile)
    all_exist   = all(os.path.exists(p) for p in outputs)

    if all_exist:
        min_out_mtime = min(os.path.getmtime(p) for p in outputs)
        if min_out_mtime > input_mtime and not myoptions.overwrite:
            print(f"Info: all outputs are up-to-date (newer than {myoptions.infile}), skipping.",
                  file=sys.stderr)
            sys.exit(0)
        elif not myoptions.overwrite:
            stale = next(p for p in outputs if os.path.getmtime(p) <= input_mtime)
            raise RuntimeError(
                f"Output is stale (older than input): {stale}\n"
                "Use --overwrite to regenerate."
            )
    else:
        for outpath in outputs:
            if os.path.exists(outpath) and not myoptions.overwrite:
                raise RuntimeError(
                    f"Output file already exists: {outpath}\n"
                    "Use --overwrite to replace it."
                )

    exact_cnt = shorter_cnt = longer_cnt = 0

    with (open(exact_name,   'w', encoding='utf-8') as exact_fh,
          open(shorter_name, 'w', encoding='utf-8') as shorter_fh,
          open(longer_name,  'w', encoding='utf-8') as longer_fh):
        for record in SeqIO.parse(myoptions.infile, myoptions.format):
            record_len = len(record.seq)
            if myoptions.debug:
                print(f"Info: Record {record.id} has length {record_len} "
                      f"of padded sequence {record.seq}")
            if record_len < myoptions.full_length:
                SeqIO.write(record, shorter_fh, "fasta-2line")
                shorter_cnt += 1
            elif record_len > myoptions.full_length:
                SeqIO.write(record, longer_fh, "fasta-2line")
                longer_cnt += 1
            else:
                SeqIO.write(record, exact_fh, "fasta-2line")
                exact_cnt += 1

    print(f"Info: Wrote {exact_cnt} entries into {exact_name}")
    print(f"Info: Wrote {shorter_cnt} entries into {shorter_name}")
    print(f"Info: Wrote {longer_cnt} entries into {longer_name}")


if __name__ == "__main__":
    main()
