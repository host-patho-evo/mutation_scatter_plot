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

import argparse
import datetime
import os
import sys

from Bio import SeqIO

VERSION = "0.3"


def _get_git_version() -> str:
    """Return ``git describe --always --dirty --tags`` output, or ``'unknown'``.

    Mirrors the three-tier logic in ``mutation_scatter_plot/__init__.py``.
    Falls back gracefully so the script always starts even without git.
    """
    import subprocess as _sp  # local import — subprocess not needed at module level
    _here = os.path.dirname(os.path.abspath(__file__))
    try:
        result = _sp.run(
            ["git", "describe", "--always", "--dirty", "--tags"],
            capture_output=True, text=True, check=True,
            cwd=_here,
        )
        ver = result.stdout.strip()
        if ver:
            return ver
    except Exception:  # pylint: disable=broad-except
        pass
    env_ver = os.environ.get("GIT_COMMIT", "").strip()
    return env_ver[:12] if env_ver else "unknown"


_GIT_VERSION: str = _get_git_version()


def build_option_parser():
    """Build and return the command-line argument parser."""
    myparser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    myparser.add_argument("--version", action="version", version=f"%(prog)s {VERSION}")
    myparser.add_argument(
        "--infile", action="store", type=str, dest="infile", default='',
        help="Input FASTA/Q file path.",
    )
    myparser.add_argument(
        "--outfile-prefix", action="store", type=str,
        dest="outfile_prefix", default='',
        help="Output file path prefix. Suffixes .exactly_N.fasta etc. are appended.",
    )
    myparser.add_argument(
        "--full-length", action="store", type=int, dest="full_length",
        default=0,
        help="Full length required for perfect alignment [0]",
    )
    myparser.add_argument(
        "--format", action="store", type=str, dest="format",
        default="fasta",
        help="Input file format. Output format is always fasta-2line.",
    )
    myparser.add_argument(
        "--overwrite", action="store_true",
        help="Overwrite output files if they already exist (bypasses staleness check).",
    )
    myparser.add_argument(
        "--debug", action="store", type=int, dest="debug", default=0,
        help="Set debug to some value",
    )
    myparser.add_argument(
        "--convert-to-upper", action="store_true", dest="convert_to_upper",
        help="Convert nucleotide strings to uppercase before writing the outputs.",
    )
    return myparser


def main():
    """Entry point: parse arguments and split FASTA file by sequence length."""
    _start_ts = datetime.datetime.now().strftime('[%Y-%m-%d %H:%M:%S]')
    print(
        f"{_start_ts} split_fasta_entries_by_lengths"
        f"  version {VERSION}  git:{_GIT_VERSION}"
        f"  invoked: {' '.join(sys.argv)}",
        file=sys.stderr,
    )
    myparser = build_option_parser()
    myoptions = myparser.parse_args()

    if not myoptions.infile:
        raise RuntimeError("Please specify --infile")
    if myoptions.full_length <= 0:
        myparser.error("--full-length must be a positive integer")
    if not os.path.exists(myoptions.infile):
        raise RuntimeError(f"Input file not found: {myoptions.infile}")
    if os.path.getsize(myoptions.infile) == 0:
        raise RuntimeError(f"Input file is empty: {myoptions.infile}")

    _exact_length_name = (
        f"{myoptions.outfile_prefix}.exactly_{myoptions.full_length!s}.fasta"
    )
    _shorter_length_name = (
        f"{myoptions.outfile_prefix}.shorter_{myoptions.full_length!s}.fasta"
    )
    _longer_length_name = (
        f"{myoptions.outfile_prefix}.longer_{myoptions.full_length!s}.fasta"
    )

    outputs = [_exact_length_name, _shorter_length_name, _longer_length_name]
    input_mtime = os.path.getmtime(myoptions.infile)
    all_exist = all(os.path.exists(p) for p in outputs)

    if all_exist:
        min_out_mtime = min(os.path.getmtime(p) for p in outputs)
        if min_out_mtime > input_mtime and not myoptions.overwrite:
            print(
                f"Info: all outputs are up-to-date (newer than {myoptions.infile}), skipping.",
                file=sys.stderr,
            )
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

    _exact_length_cnt   = 0
    _shorter_length_cnt = 0
    _longer_length_cnt  = 0

    with (open(_exact_length_name,   'w', encoding='utf-8') as _exact_length,
          open(_shorter_length_name, 'w', encoding='utf-8') as _shorter_length,
          open(_longer_length_name,  'w', encoding='utf-8') as _longer_length):

        for _record in SeqIO.parse(myoptions.infile, myoptions.format):
            if myoptions.convert_to_upper:
                _record.seq = _record.seq.upper()
                
            if myoptions.debug:
                print(
                    f"Info: Record {_record.id} has length {len(_record.seq)} "
                    f"of padded sequence {_record.seq}"
                )
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
