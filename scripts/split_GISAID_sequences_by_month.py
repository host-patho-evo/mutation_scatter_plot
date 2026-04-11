#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# vim: set fileencoding=utf-8 ts=4 sw=4 expandtab :
"""Split GISAID FASTA files into per-month (YYYY-MM) subsets.

The GISAID header format places the collection date in varying
``|``-delimited columns.  All columns are scanned to locate a
date-like value.  The script also handles dates without a day
component (``YYYY-MM``) and the ``DD-Mon-YYYY`` format.

Example headers (all supported)::

    >Spike|hCoV-19/Wuhan/WIV04/2019|2019-12-30|EPI_ISL_402124|...
    >Spike|hCoV-19/France/BRE-IPP16678/2021|EPI_ISL_3219515|2021-07-13|...
    >Spike|hCoV-19/example/2021|EPI_ISL_999|2021-07|NorthAmerica
    >Spike|hCoV-19/example/2021|EPI_ISL_999|01-Jul-2021|NorthAmerica

Some virus names contain stray ``|`` characters, e.g.::

    >Spike|hCoV-19/USA/LA-OD-|O-4336284453/2023|EPI_ISL_16576342|2023-01-10|NorthAmerica

These are detected and sed edit recipes are written to a
``.stray_pipes.sed`` file for fixing the input.

Output
------
For each distinct YYYY-MM found in the input, a file
``{outfile_prefix}.{YYYY-MM}.fasta`` is written containing all sequences
deposited in that month.

Usage::

    split_GISAID_sequences_by_month.py \\\\
        --infilename spikenuc1207.fasta \\\\
        --outfile-prefix spikenuc1207 --debug 1
"""

import argparse
import re
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

def sanitize_year_and_month(
    year, month, original_date='', header_line='',
):
    """Validate and zero-pad year/month strings.

    Handles numeric months (``'1'`` → ``'01'``, ``'07'``) as well as
    abbreviated month names (``'Jul'`` → ``'07'``).

    Examples::

        >>> sanitize_year_and_month('2023', '7')
        ('2023', '07')
        >>> sanitize_year_and_month('2021', 'Jul')
        ('2021', '07')
    """
    _hdr = header_line.rstrip('\n')
    _ctx = ''
    if original_date:
        _ctx += f" (date='{original_date}'"
        if _hdr:
            _ctx += f", line='{_hdr}'"
        _ctx += ')'
    elif _hdr:
        _ctx += f" (line='{_hdr}')"

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
    """Swap year/day when the date is in DD-Mon-YYYY order.

    Example::

        >>> _normalize_date_parts('01', 'Jul', '2021')
        ('2021', 'Jul', '01')
        >>> _normalize_date_parts('2021', '07', '01')
        ('2021', '07', '01')
    """
    if len(year) <= 2 and len(day) == 4:
        return day, month, year
    return year, month, day


def _looks_like_date(text):
    """Return True if *text* could be a date (YYYY-MM-DD, YYYY-MM, DD-Mon-YYYY)."""
    parts = text.split('-')
    if len(parts) < 2:
        return False
    # YYYY-MM or YYYY-MM-DD
    if len(parts[0]) == 4 and parts[0].isdigit():
        return True
    # DD-Mon-YYYY
    if len(parts[0]) <= 2 and parts[0].isdigit() and len(parts) == 3:
        return True
    return False


def _parse_date_column(date_str):
    """Split a date string into (year, month, day) parts.

    Returns ``(year, month, day)``; *day* may be ``''`` for YYYY-MM.
    """
    parts = date_str.split('-')
    if len(parts) == 3:
        year, month, day = parts
    elif len(parts) == 2:
        year, month = parts
        day = ''
    else:
        raise ValueError(f"Unexpected date format '{date_str}'")
    if day:
        year, month, day = _normalize_date_parts(year, month, day)
    return year, month, day


def _is_epi_id(text):
    """Return True if *text* looks like a GISAID accession."""
    return bool(re.match(r'EPI_ISL_\d+$', text.strip()))


def _parse_header_date(line):
    """Parse (year, month) from a GISAID FASTA header line.

    Scans all ``|``-delimited columns (after the first) for a value
    that looks like a date, to handle varying GISAID column orderings.

    Supported date formats:

    - ``YYYY-MM-DD`` (e.g. ``2023-01-10``)
    - ``YYYY-MM``    (e.g. ``2021-07``, no day component)
    - ``DD-Mon-YYYY`` (e.g. ``01-Jul-2021``)

    Examples::

        >>> _parse_header_date(
        ...     '>Spike|hCoV-19/Wuhan/WIV04/2019|2019-12-30|EPI_ISL_402124|Asia')
        ('2019', '12')
        >>> _parse_header_date(
        ...     '>Spike|hCoV-19/example/2021|EPI_ISL_999|01-Jul-2021|Region')
        ('2021', '07')
        >>> _parse_header_date(
        ...     '>Spike|hCoV-19/USA/LA-OD-|O-4336284453/2023|EPI_ISL_16576342|2023-01-10|NA')
        ('2023', '01')
    """
    columns = line.split('|')[1:]   # skip the gene/sequence name
    if not columns:
        raise ValueError(
            "No pipe-delimited columns in line '%s'"
            % line.rstrip('\n')
        )

    for col in columns:
        col = col.strip()
        if _looks_like_date(col):
            _year, _month, _day = _parse_date_column(col)
            return sanitize_year_and_month(
                _year, _month, original_date=col,
                header_line=line,
            )

    raise ValueError(
        "No date column found in line '%s'" % line.rstrip('\n')
    )


def _detect_stray_pipes(line):
    """Detect stray ``|`` inside the virus name.

    For example, the header::

        >Spike|hCoV-19/USA/LA-OD-|O-4336284453/2023|EPI_ISL_16576342|2023-01-10|NA

    has a stray ``|`` splitting ``hCoV-19/USA/LA-OD-O-4336284453/2023``.

    Returns ``None`` if the header looks normal, or a dict with:

    - ``original``: the original header (stripped)
    - ``fixed``:    the corrected header with stray pipes removed
    - ``virusname_original``: the broken virus name fragments
    - ``virusname_fixed``:    the reconstructed virus name
    - ``extra_pipes``: number of stray pipes found
    """
    raw = line.rstrip('\n')
    all_cols = raw.split('|')

    # Identify which columns are "known" (date or EPI_ISL).
    # Everything between the gene (col 0) and the first known
    # column is the virus name (possibly split by stray pipes).
    first_known_idx = None
    for idx, col in enumerate(all_cols[1:], start=1):
        stripped = col.strip()
        if _looks_like_date(stripped) or _is_epi_id(stripped):
            first_known_idx = idx
            break

    # Virus name should be at index 1; if the first known column
    # is at index 2, that's normal (gene|virusname|date_or_epi|...).
    if first_known_idx is None or first_known_idx <= 2:
        return None

    # Columns 1..first_known_idx-1 are fragments of the virus name
    virusname_fragments = all_cols[1:first_known_idx]
    virusname_fixed = ''.join(virusname_fragments)
    fixed_cols = (
        [all_cols[0]]
        + [virusname_fixed]
        + all_cols[first_known_idx:]
    )

    return {
        'original': raw,
        'fixed': '|'.join(fixed_cols),
        'virusname_original': '|'.join(virusname_fragments),
        'virusname_fixed': virusname_fixed,
        'extra_pipes': first_known_idx - 2,
    }


# ---------------------------------------------------------------------------
# Core logic
# ---------------------------------------------------------------------------

def get_list_of_yyyymm_present(infilename):
    """Pre-parse input file to collect distinct YYYY-MM dates.

    We cannot use ``Bio.SeqIO`` because the GISAID "FASTA" format
    contains white spaces inside the description line.

    Returns
    -------
    yyyymm_dictionary : dict
        Mapping of YYYY-MM strings to empty lists.
    stray_pipe_records : list of dict
        Records describing headers with stray ``|`` in the virus
        name, suitable for generating edit recipes.
    """
    _yyyymm_dictionary = defaultdict(list)
    _stray_pipe_records = []
    with open(infilename) as _infile:
        for _lineno, _line in enumerate(_infile, start=1):
            if _line[0] == '>':
                _year, _month = _parse_header_date(_line)
                _yyyymm = f"{_year}-{_month}"
                if _yyyymm not in _yyyymm_dictionary:
                    _yyyymm_dictionary[_yyyymm] = []
                _stray = _detect_stray_pipes(_line)
                if _stray is not None:
                    _stray['lineno'] = _lineno
                    _stray_pipe_records.append(_stray)
    return _yyyymm_dictionary, _stray_pipe_records


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

    yyyymm_dictionary, stray_pipe_records = \
        get_list_of_yyyymm_present(args.infilename)

    if stray_pipe_records:
        _recipe_path = f"{outfile_prefix}.stray_pipes.sed"
        with open(_recipe_path, 'w') as _fh:
            _fh.write(
                "# sed recipes to fix stray pipe characters "
                "in virus names\n"
            )
            _fh.write(
                f"# Generated from: {args.infilename}\n"
            )
            _fh.write(
                f"# Total headers with stray pipes: "
                f"{len(stray_pipe_records)}\n\n"
            )
            for rec in stray_pipe_records:
                # Escape sed special chars in the original
                _orig_esc = rec['original'].replace(
                    '/', '\\/'
                ).replace('&', '\\&')
                _fixed_esc = rec['fixed'].replace(
                    '/', '\\/'
                ).replace('&', '\\&')
                _fh.write(
                    f"# Line {rec['lineno']}: "
                    f"{rec['extra_pipes']} stray pipe(s) "
                    f"in virus name "
                    f"'{rec['virusname_original']}'\n"
                )
                _fh.write(
                    f"s/{_orig_esc}/{_fixed_esc}/\n\n"
                )
        print(
            f"Warning: {len(stray_pipe_records)} header(s) with "
            f"stray '|' in virus name. "
            f"Edit recipes written to: {_recipe_path}"
        )
        print(
            f"  Applying: sed -i -f {_recipe_path} "
            f"{args.infilename}"
        )
        subprocess.run(
            ['sed', '-i', '-f', _recipe_path, args.infilename],
            check=True,
        )
        # Re-parse after fixing the input file
        yyyymm_dictionary, stray_pipe_records2 = \
            get_list_of_yyyymm_present(args.infilename)
        if stray_pipe_records2:
            print(
                f"Warning: {len(stray_pipe_records2)} header(s) "
                f"still have stray pipes after fix"
            )

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
