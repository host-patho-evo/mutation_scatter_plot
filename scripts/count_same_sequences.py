#! /usr/bin/env python3

VERSION = "202603291800"

from optparse import OptionParser
import subprocess, shlex
import os, sys, io
import hashlib

myparser = OptionParser(version="%s version %s" % ('%prog', VERSION))
myparser.add_option("--infilename", action="store", type="string", dest="infilename", default='stdin.fasta',
    help="Input FASTA/Q file path.")
myparser.add_option("--infile-format", action="store", type="string", dest="infile_format", default='fasta',
    help="Input FASTA/Q file format [default: fasta]. Outfile  is always fasta.")
myparser.add_option("--outfile-prefix", action="store", type="string", dest="outfile_prefix", default='',
    help="Output FASTA file path prefix. It will be appended by '.counts.fast[aq]'.")
myparser.add_option("--sort-bucket-size", action="store", type="string", dest="sort_bucket_size", default='30%',
    help="Size of the sort buffer bucket size to be used by UNIX sort commain in TMPDIR [default: 30% of RAM].")
myparser.add_option("--min-count", action="store", type="int", dest="min_count", default=0,
    help="Write out sequences above some minimum incidence into {outfile_prefix}.min_count_{min_count}.fast[aq]")
myparser.add_option("--top-n", action="store", type="int", dest="top_n", default=0,
    help="Write out the top n-most sequences into {outfile_prefix}.top_{top_n}_counts.fastq or {outfile_prefix}.counts.fast[aq]")
myparser.add_option("--mapping-outfile", action="store", type="string", dest="mapping_outfile", default='',
    help=(
        "Optional TSV file mapping sha256 -> original FASTA IDs. "
        "Columns (tab-separated, no header): sha256hex, count, id_1, id_2, ..."
        "Allows tracing any deduplicated record back to all original sequences "
        "that share its sequence content."
    ))
myparser.add_option("--debug", action="store", type="int", dest="debug", default=0,
    help="Set debug to some value")
(myoptions, myargs) = myparser.parse_args()

from collections import defaultdict

"""
Read a FASTQ file and return a FASTA file with sequence counts. It truncates the output file if it already existed.
It also discards any description text from input FASTA data due to 'awk 'NR %% 2 == 0' infile | sort -S %s | uniq -c | sort -S %s -nr' .

Usage:
find . -name *.fastp.amplicons.fasta | sort | uniq | grep -v results | while read f; do
  d=`dirname $f`;
  p=`basename $f .fasta`;
  echo "$d"/"$p";
  if [ ! -e "$d"/"$p".counts.fasta ]; then
    count_same_sequences.py --infilename="$f" --infile-format=fasta --outfile-prefix="$d"/"$p";
  fi;
done

TODO:
Merge paired-end mates together and rerun this script
"""

def read_and_count_sequences(infilename, outfileh, infile_format, top_n=0, min_count=0, sort_bucket_size='10G'):
    """Count unique sequences already padded with dashes. Do not remove the dashes
    otherwise slicing later fails on shorter sequences, as it uses the indexes from
    the padded alignment.

    This runs the command
    reformat.sh fastawrap=0 in=%s out=stdout.fasta | awk 'NR % 2 == 0' | sort -S $sort_bucket_size | uniq -c | sort -nr | head -n 100
    and parse output with counts.

    Unix sort respects TMPDIR variable to place the temporary files i there, instead of cwd. To make counting faster
    we allow to specify percentage of total RAM to be used '50%' or specify size as '10G' from the [KMGTP].
    """

    if top_n:
        _cmd = "cat %s | reformat.sh fastawrap=0 in=stdin.%s out=stdout.fasta ignorejunk=t simd=f | awk 'NR %% 2 == 0' | sort -S %s | uniq -c | sort -S %s -nr | head -n %d" % (infilename, infile_format, sort_bucket_size, sort_bucket_size, top_n)
    elif min_count:
        _cmd = "cat %s | reformat.sh fastawrap=0 in=stdin.%s out=stdout.fasta ignorejunk=t simd=f | awk 'NR %% 2 == 0' | sort -S %s | uniq -c | sort -S %s -nr | awk '{if ($1 >= %d) print}'" % (infilename, infile_format, sort_bucket_size, sort_bucket_size, min_count)
    else:
        _cmd = "cat %s | reformat.sh fastawrap=0 in=stdin.%s out=stdout.fasta ignorejunk=t simd=f | awk 'NR %% 2 == 0' | sort -S %s | uniq -c | sort -S %s -nr" % (infilename, infile_format, sort_bucket_size, sort_bucket_size)
    #_cmdlist = shlex.split(_cmd)
    print("Info: %s" % str(_cmd))
    _stdout, _stderr = subprocess.Popen(_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, shell=True).communicate()
    #
    # check if java binary is available in PATH
    # if _stderr and _stderr[0]:
    #     raise RuntimeError("%s failed with:\n%s" % (_cmd, str(_stderr)))
    sys.stdout.flush()
    sys.stderr.flush()

    # print("Result: %s" % str(_stdout))
    for _line in io.StringIO(_stdout).readlines():
        if _line:
            # print("_line: %s" % _line)
            try:
                _count, _sequence = _line.split()
            except ValueError:
                try:
                    _count, _sequence, _hash = _line.split()
                except ValueError:
                    # do not break on some entries from GISAID
                    # pick the first item and the very last item
                    # 1  Cytologiczn Szpital Specjalistyczny im Edmunda Biernackiego w Mielcu|1. Tricity SARS-CoV-2 sequencing ATGTTT...
                    _myvalues = _line.split()
                    _count, _sequence = _myvalues[0], _myvalues[-1]
            _m = hashlib.sha256() # always create a new object because .update() would just append to it new string incrementally
            _m.update(_sequence.encode()) # calculate the hash without a trailing newline
            _hash = _m.hexdigest()
            # print(">%sx\n%s" % (_count, _sequence))
            outfileh.write(">%sx.%s%s%s%s" % (_count, _hash, os.linesep, _sequence, os.linesep)) # use a dot to merge the hash to the count followed by 'x' so that NCBI blastn does not discard the hash from results
            del(_m)


def build_sha256_id_mapping(infilename, mapping_outfile):
    """Build a sha256 -> original FASTA IDs translation table.

    The external sort | uniq -c pipeline in read_and_count_sequences() discards
    FASTA headers, so this function performs a separate single-pass scan of the
    original file to reconstruct the mapping.

    For each record the same SHA-256 is computed as in read_and_count_sequences()
    (raw sequence bytes, no case conversion, no trailing newline) so the hashes
    match exactly.

    Output TSV columns (tab-separated, no header line):
        sha256hex   count   id_1    id_2    ...

    Args:
        infilename:     Path to the original (non-deduplicated) FASTA file.
        mapping_outfile: Path for the output TSV.
    """
    # sha256 -> [count, [id1, id2, ...]]
    _mapping = {}
    _n_in = 0

    _name = None
    _seq_parts = []

    def _flush(name, seq_parts):
        _seq = "".join(seq_parts).rstrip("\r\n")
        _digest = hashlib.sha256(_seq.encode()).hexdigest()
        if _digest in _mapping:
            _mapping[_digest][0] += 1
            _mapping[_digest][1].append(name)
        else:
            _mapping[_digest] = [1, [name]]

    with open(infilename, "r", encoding="utf-8", errors="replace") as _fh:
        for _line in _fh:
            _line = _line.rstrip("\r\n")
            if not _line:
                continue
            if _line[0] == ">":
                if _name is not None:
                    _flush(_name, _seq_parts)
                    _n_in += 1
                    if myoptions.debug and _n_in % 500_000 == 0:
                        print("Info: mapping pass: processed %d records, %d unique" % (_n_in, len(_mapping)), file=sys.stderr)
                # first word after '>' is the ID
                _parts = _line[1:].split()
                _name = _parts[0] if _parts else ""
                _seq_parts = []
            else:
                _seq_parts.append(_line)
        if _name is not None:
            _flush(_name, _seq_parts)
            _n_in += 1

    print("Info: mapping pass: %d total records -> %d unique sequences" % (_n_in, len(_mapping)), file=sys.stderr)

    with open(mapping_outfile, "w", encoding="utf-8") as _out:
        for _digest, (_count, _ids) in _mapping.items():
            _out.write("\t".join([_digest, str(_count)] + _ids) + "\n")

    print("Info: wrote sha256->IDs mapping (%d entries) to %s" % (len(_mapping), mapping_outfile), file=sys.stderr)


if not myoptions.infilename:
    raise RuntimeError("Please provide input filename via --infilename")
elif not os.path.exists(myoptions.infilename):
    raise RuntimeError("File %s does not exist, please check --infilename" % myoptions.infilename)

if not myoptions.outfile_prefix:
    if myoptions.infilename.endswith('.fastq.gz'):
        _outfile_prefix = myoptions.infilename.replace('.fastq.gz', '')
    elif myoptions.infilename.endswith('.fastq'):
        _outfile_prefix = myoptions.infilename.replace('.fastq', '')
    elif myoptions.infilename.endswith('.fasta.gz'):
        _outfile_prefix = myoptions.infilename.replace('.fasta.gz', '')
    elif myoptions.infilename.endswith('.fasta'):
        _outfile_prefix = myoptions.infilename.replace('.fasta', '')
    else:
        _outfile_prefix = myoptions.infilename
else:
    _outfile_prefix = myoptions.outfile_prefix

_outfileh = open(_outfile_prefix + '.counts.fasta', 'w') # truncate the output file

if myoptions.min_count > 0:
    read_and_count_sequences(myoptions.infilename, _outfileh, infile_format=myoptions.infile_format, min_count=myoptions.min_count, sort_bucket_size=myoptions.sort_bucket_size)
elif myoptions.top_n > 0:
    read_and_count_sequences(myoptions.infilename, _outfileh, infile_format=myoptions.infile_format, top_n=myoptions.top_n, sort_bucket_size=myoptions.sort_bucket_size)
else:
    read_and_count_sequences(myoptions.infilename, _outfileh, infile_format=myoptions.infile_format, sort_bucket_size=myoptions.sort_bucket_size)

_outfileh.close()

if myoptions.mapping_outfile:
    build_sha256_id_mapping(myoptions.infilename, myoptions.mapping_outfile)
