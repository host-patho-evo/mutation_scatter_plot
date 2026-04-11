#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# vim: set fileencoding=utf-8 ts=4 sw=4 expandtab :
"""Split GISAID FASTA files into per-month (YYYY-MM) subsets.

The GISAID header format places the collection date either in the 3rd or
4th ``|``-delimited column.  Both variants are handled transparently.

Example headers::

    >Spike|hCoV-19/Wuhan/WIV04/2019|2019-12-30|EPI_ISL_402124|...
    >Spike|hCoV-19/France/BRE-IPP16678/2021|EPI_ISL_3219515|2021-07-13|...

Output
------
For each distinct YYYY-MM found in the input, a file
``{outfile_prefix}.{YYYY-MM}.fasta`` is written containing all sequences
deposited in that month.

Usage::

    split_GISAID_sequences_by_month.py \\
        --infilename spikenuc1207.fasta \\
        --outfile-prefix spikenuc1207 --debug 1
"""

import argparse
import os
import subprocess
from collections import defaultdict

VERSION = "202604111600"


def _get_git_version() -> str:
    """Return ``git describe --always --dirty --tags``, or ``'unknown'``."""
    _here = os.path.dirname(os.path.abspath(__file__))
    try:
        result = subprocess.run(
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


# ---------------------------------------------------------------------------
# Header parsing
# ---------------------------------------------------------------------------

def sanitize_year_and_month(year, month, original_date=''):
    """Validate and zero-pad year/month strings."""
    _ctx = f" (from '{original_date}')" if original_date else ''

    if len(month) < 2:
        _month = f"0{month}"
    else:
        _month = month

    if _month == '00':
        print(f"Info: Weird month in {year}-{month}{_ctx}")
    else:
        try:
            if int(_month.lstrip('0')) > 12:
                raise ValueError(
                    f"Weird month {_month}{_ctx}"
                )
        except ValueError as exc:
            _month_names = {
                'Jan': '01', 'Feb': '02', 'Mar': '03',
                'Apr': '04', 'May': '05', 'Jun': '06',
                'Jul': '07', 'Aug': '08', 'Sep': '09',
                'Oct': '10', 'Nov': '11', 'Dec': '12',
            }
            _stripped = _month.lstrip('0')
            if _stripped not in _month_names:
                raise ValueError(
                    f"Weird month {_month}{_ctx}"
                ) from exc
            _month = _month_names[_stripped]

    if len(year) < 4:
        raise ValueError(f"Weird year {year}{_ctx}")

    return year, _month


def _normalize_date_parts(year, month, day):
    """Swap year/day when the date is in DD-Mon-YYYY order."""
    if len(year) <= 2 and len(day) == 4:
        return day, month, year
    return year, month, day


def _parse_header_date(line):
    """Parse (year, month) from a GISAID FASTA header line.

    Handles both column orderings:
      ``>Spike|virusname|DATE|EPI_ID|...``
      ``>Spike|virusname|EPI_ID|DATE|...``

    Also handles DD-Mon-YYYY date formats (e.g. ``01-Jul-2021``).
    """
    try:
        _virusname, _date, _epi_id = line.split('|')[1:4]
    except ValueError as exc:
        raise ValueError(
            "Cannot split entries from line '%s'" % str(line)
        ) from exc

    try:
        _year, _month, _day = _date.split('-')
    except ValueError:
        # Date is in the 4th column instead of the 3rd
        try:
            _virusname, _epi_id, _date = line.split('|')[1:4]
        except ValueError as exc:
            raise ValueError(
                "Cannot split date from line '%s'" % str(line)
            ) from exc
        _year, _month, _day = _date.split('-')

    _year, _month, _day = _normalize_date_parts(_year, _month, _day)

    return sanitize_year_and_month(_year, _month, original_date=_date)


# ---------------------------------------------------------------------------
# Core logic
# ---------------------------------------------------------------------------

def get_list_of_yyyymm_present(infilename):
    """Pre-parse input file to collect distinct YYYY-MM dates.

    We cannot use ``Bio.SeqIO`` because the GISAID "FASTA" format
    contains white spaces inside the description line.
    """
    _yyyymm_dictionary = defaultdict(list)
    with open(infilename) as _infile:
        for _line in _infile:
            if _line[0] == '>':
                _year, _month = _parse_header_date(_line)
                _yyyymm = f"{_year}-{_month}"
                if _yyyymm not in _yyyymm_dictionary:
                    _yyyymm_dictionary[_yyyymm] = []
    return _yyyymm_dictionary


def split_sequences_by_month(
    infilename, outfile_prefix, yyyymm_dictionary, debug=0
):
    """Write per-month FASTA files."""
    for yyyymm in sorted(yyyymm_dictionary.keys()):
        with (
            open(infilename) as _infile,
            open(
                f"{outfile_prefix}.{yyyymm}.fasta", 'w'
            ) as _outfileh,
        ):
            _output_following_lines = False
            for _line in _infile:
                if _line[0] == '>':
                    _output_following_lines = False
                    _year, _month = _parse_header_date(_line)
                    _yyyymm = f"{_year}-{_month}"
                    if _yyyymm == yyyymm:
                        _outfileh.write(_line)
                        _output_following_lines = True
                elif _output_following_lines:
                    _outfileh.write(_line)
        if debug:
            print(
                f"Debug: Closed {outfile_prefix}.{yyyymm}.fasta"
            )


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main():
    """Entry point."""
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "--infilename", required=True,
        help="Input FASTA file path.")
    parser.add_argument(
        "--outfile-prefix", dest="outfile_prefix", default='',
        help="Output FASTA file path prefix. "
             "Appended by '.YYYY-MM.fasta'.")
    parser.add_argument(
        "--debug", type=int, default=0,
        help="Set debug to some value [default: 0].")
    parser.add_argument(
        "--version", action="version",
        version=f"%(prog)s {VERSION}  git:{_GIT_VERSION}")
    args = parser.parse_args()

    if not os.path.exists(args.infilename):
        parser.error(
            f"File {args.infilename} does not exist, "
            f"please check --infilename"
        )

    if not args.outfile_prefix:
        for _ext in ('.fastq.gz', '.fastq', '.fasta.gz', '.fasta'):
            if args.infilename.endswith(_ext):
                outfile_prefix = args.infilename[:-len(_ext)]
                break
        else:
            outfile_prefix = args.infilename
    else:
        outfile_prefix = args.outfile_prefix

    print(
        f"split_GISAID_sequences_by_month"
        f"  version {VERSION}  git:{_GIT_VERSION}"
    )

    yyyymm_dictionary = get_list_of_yyyymm_present(args.infilename)
    if args.debug:
        print(
            f"Debug: Entries are from the following months: "
            f"{sorted(yyyymm_dictionary.keys())}"
        )
    split_sequences_by_month(
        args.infilename, outfile_prefix, yyyymm_dictionary,
        debug=args.debug,
    )


if __name__ == "__main__":
    main()
