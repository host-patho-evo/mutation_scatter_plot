#! /usr/bin/env python3
"""Split GISAID FASTA files by month (YYYY-MM).

The GISAID header format is broken and the date appears either
in the 3rd or 4th column. Great design!

Example headers::

    >Spike|hCoV-19/Wuhan/WIV04/2019|2019-12-30|EPI_ISL_402124|...
    >Spike|hCoV-19/France/BRE-IPP16678/2021|EPI_ISL_3219515|2021-07-13|...

Usage::

    split_GISAID_sequences_by_month.py \\
        --infilename=spikenuc1207.fasta \\
        --outfile-prefix=spikenuc1207 --debug=1
"""

import os
from collections import defaultdict
from optparse import OptionParser  # noqa: E501  # legacy, kept for compat

VERSION = "202601021910"

myparser = OptionParser(version="%s version %s" % ('%prog', VERSION))
myparser.add_option(
    "--infilename", action="store", type="string",
    dest="infilename", default='stdin.fasta',
    help="Input FASTA/Q file path.")
myparser.add_option(
    "--infile-format", action="store", type="string",
    dest="infile_format", default='fasta',
    help="Input FASTA/Q file format [default: fasta]. "
         "Outfile is always fasta.")
myparser.add_option(
    "--outfile-prefix", action="store", type="string",
    dest="outfile_prefix", default='',
    help="Output FASTA file path prefix. "
         "It will be appended by '.counts.fast[aq]'.")
myparser.add_option(
    "--debug", action="store", type="int",
    dest="debug", default=0,
    help="Set debug to some value")
(myoptions, myargs) = myparser.parse_args()


def sanitize_year_and_month(year, month):
    """Validate and zero-pad year/month strings."""
    if len(month) < 2:
        _month = f"0{month}"
    else:
        _month = month

    if _month == '00':
        # raise ValueError(f"Weird month {_month}")
        print(f"Info: Weird month in {year}-{month}")
    elif int(_month.lstrip('0')) > 12:
        raise ValueError(f"Weird month {_month}")

    if len(year) < 4:
        raise ValueError(f"Weird year {year}")
    else:
        _year = year

    return _year, _month


def _parse_header_date(line):
    """Parse YYYY, MM from a GISAID header line.

    Handles both column orderings:
      >Spike|virusname|DATE|EPI_ID|...
      >Spike|virusname|EPI_ID|DATE|...
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
        # Date is in 4th column instead of 3rd
        # >Spike|hCoV-19/.../2021|EPI_ISL_...|2021-07-13|...
        try:
            _virusname, _epi_id, _date = line.split('|')[1:4]
        except ValueError as exc:
            raise ValueError(
                "Cannot split date from line '%s'" % str(line)
            ) from exc
        _year, _month, _day = _date.split('-')

    return sanitize_year_and_month(_year, _month)


def get_list_of_yyyymm_present(infilename, infile_format):
    """Pre-parse input file to collect YYYY-MM dates.

    We cannot use Bio.SeqIO because the GISAID "FASTA" format
    contains white spaces inside the description line.
    """
    _yyyymm_dictionary = defaultdict(list)
    with open(infilename) as _infile:
        for _line in _infile:
            if _line[0] == '>':
                _year, _month = _parse_header_date(_line)
                _yyyymm = f"{_year}-{_month}"
                # We only need to track that this month exists;
                # the EPI IDs are not used downstream.
                if _yyyymm not in _yyyymm_dictionary:
                    _yyyymm_dictionary[_yyyymm] = []
    if myoptions.debug:
        print(
            f"Debug: Entries are from the following months: "
            f"{sorted(_yyyymm_dictionary.keys())}"
        )
    return _yyyymm_dictionary


def split_sequences_by_month(infilename, outfile_prefix, yyyymm_dictionary):
    """Write per-month FASTA files."""
    for yyyymm in sorted(yyyymm_dictionary.keys()):
        with (
            open(infilename) as _infile,
            open(f"{outfile_prefix}.{yyyymm}.fasta", 'w') as _outfileh,
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
        if myoptions.debug:
            print(f"Debug: Closed {outfile_prefix}.{yyyymm}.fasta")


if not myoptions.infilename:
    raise RuntimeError("Please provide input filename via --infilename")
if not os.path.exists(myoptions.infilename):
    raise RuntimeError(
        "File %s does not exist, please check --infilename"
        % myoptions.infilename
    )

if not myoptions.outfile_prefix:
    for _ext in ('.fastq.gz', '.fastq', '.fasta.gz', '.fasta'):
        if myoptions.infilename.endswith(_ext):
            _outfile_prefix = myoptions.infilename[:-len(_ext)]
            break
    else:
        _outfile_prefix = myoptions.infilename
else:
    _outfile_prefix = myoptions.outfile_prefix

# Create a list of YYYYMM entries parsed from the input file
_yyyymm_dictionary = get_list_of_yyyymm_present(
    myoptions.infilename, infile_format=myoptions.infile_format
)
split_sequences_by_month(
    myoptions.infilename, _outfile_prefix, _yyyymm_dictionary
)
