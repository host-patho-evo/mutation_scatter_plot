#!/usr/bin/env python3
"""Count identical sequences in a FASTA/Q file and write a deduplicated output.

Each unique sequence gets a FASTA ID of the form::

    {count}x.{sha256hex}

where *count* is the number of occurrences and *sha256hex* is the 64-character
hexadecimal SHA-256 of the uppercase sequence (alignment dashes preserved,
CR/LF stripped) — identical to what ``reformat.sh fastawrap=0`` produces.

In addition to the deduplicated FASTA (``{prefix}.counts.fasta``) the script
produces a TSV mapping file (``{prefix}.sha256_to_ids.tsv``) that maps every
sha256 to the list of original FASTA IDs sharing that sequence.  This enables
later traceability (see ``create_list_of_discarded_sequences.py``).

Behaviour
---------
* Output guard: if both output files already exist and are **newer** than the
  input the script prints ``Info: up-to-date, skipping`` and exits 0 without
  doing any work.  Pass ``--overwrite`` to force regeneration.
* If an output exists but is **older** than the input the script raises a
  ``RuntimeError`` — the file is stale; add ``--overwrite`` to regenerate.
* Values containing ``=`` passed to ``--outfile-prefix`` or
  ``--mapping-outfile`` are rejected with a helpful suggestion, catching the
  common typo ``--outfile-prefix=infilename=filename_prefix``.

Example usage::

    count_same_sequences.py \\
        --infilename=filename_prefix.fasta \\
        --outfile-prefix=filename_prefix \\
        --mapping-outfile=filename_prefix.sha256_to_ids.tsv

    # Batch loop (skips already-finished files automatically):
    find . -name '*.fasta' | while read f; do
      p="${f%.fasta}"
      count_same_sequences.py --infilename="$f" --outfile-prefix="$p"
    done
"""

import argparse
import hashlib
import io
import os
import re
import subprocess
import sys

VERSION = "202603292130"


def _decode_fasta_line(raw: bytes) -> str:
    r"""Decode one raw FASTA byte line to a clean Unicode str.

    Handles GISAID FASTA files that mix UTF-8 and Latin-1 *within the same
    line* (e.g. a header containing both Polish UTF-8 multi-byte chars like
    ``ś`` (\xC5\x9B) and raw Latin-1 bytes like ``é`` (\xe9)).

    Strategy:
      1. Decode as UTF-8 with ``errors='surrogateescape'``: valid multi-byte
         sequences decode normally; every byte that is not part of a valid
         UTF-8 sequence becomes a surrogate code point (U+DC80..U+DCFF).
      2. Map each surrogate U+DC80+b back to chr(b), i.e. the Latin-1
         interpretation of the original byte.  This is lossless — all bytes
         end up correctly decoded.
      3. Convert literal ``\uXXXX`` escape sequences (produced by some GISAID
         export pipelines) to real Unicode codepoints.
    """
    s = raw.decode("utf-8", errors="surrogateescape")
    if any(0xDC80 <= ord(ch) <= 0xDCFF for ch in s):
        s = "".join(
            chr(ord(ch) - 0xDC00) if 0xDC80 <= ord(ch) <= 0xDCFF else ch
            for ch in s
        )
    # Convert literal \uXXXX escape sequences to real codepoints.
    s = re.sub(r'\\u([0-9a-fA-F]{4})',
               lambda m: chr(int(m.group(1), 16)), s)
    return s


def read_and_count_sequences(infilename, outfileh, infile_format,
                             top_n=0, min_count=0, sort_bucket_size='10G'):
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
        cmd = (
            f"cat {infilename} | reformat.sh fastawrap=0 in=stdin.{infile_format}"
            f" out=stdout.fasta ignorejunk=t simd=f"
            f" | awk 'NR %% 2 == 0'"
            f" | sort -S {sort_bucket_size}"
            f" | uniq -c | sort -S {sort_bucket_size} -nr | head -n {top_n}"
        )
    elif min_count:
        cmd = (
            f"cat {infilename} | reformat.sh fastawrap=0 in=stdin.{infile_format}"
            f" out=stdout.fasta ignorejunk=t simd=f"
            f" | awk 'NR %% 2 == 0'"
            f" | sort -S {sort_bucket_size}"
            f" | uniq -c | sort -S {sort_bucket_size} -nr"
            f" | awk '{{if ($1 >= {min_count}) print}}'"
        )
    else:
        cmd = (
            f"cat {infilename} | reformat.sh fastawrap=0 in=stdin.{infile_format}"
            f" out=stdout.fasta ignorejunk=t simd=f"
            f" | awk 'NR %% 2 == 0'"
            f" | sort -S {sort_bucket_size}"
            f" | uniq -c | sort -S {sort_bucket_size} -nr"
        )
    print(f"Info: {cmd}")
    with subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                          text=True, shell=True) as proc:
        stdout, _stderr = proc.communicate()
    sys.stdout.flush()
    sys.stderr.flush()

    for line in io.StringIO(stdout).readlines():
        if line:
            try:
                count, sequence = line.split()
            except ValueError:
                try:
                    count, sequence, _hash = line.split()
                except ValueError:
                    # do not break on some entries from GISAID
                    # pick the first item and the very last item
                    # 1  Cytologiczn Szpital Specjalistyczny im Edmunda Biernackiego w Mielcu|1. Tricity SARS-CoV-2 sequencing ATGTTT...
                    myvalues = line.split()
                    count, sequence = myvalues[0], myvalues[-1]
            digest = hashlib.sha256(sequence.encode()).hexdigest()
            outfileh.write(f">{count}x.{digest}{os.linesep}{sequence}{os.linesep}")


def build_sha256_id_mapping(infilename, mapping_outfile, debug=0):
    """Build a sha256 -> original FASTA IDs translation table.

    The external sort | uniq -c pipeline in read_and_count_sequences() discards
    FASTA headers, so this function performs a separate single-pass scan of the
    original file to reconstruct the mapping.

    For each record the same SHA-256 is computed as in read_and_count_sequences():
    uppercase sequence, dashes preserved (alignment padding must not be removed),
    no embedded newlines — identical to what reformat.sh + sort|uniq produces.

    Output TSV columns (tab-separated, no header line):
        sha256hex   count   id_1    id_2    ...

    Args:
        infilename:      Path to the original (non-deduplicated) FASTA file.
        mapping_outfile: Path for the output TSV.
        debug:           Debug verbosity level.
    """
    mapping = {}
    n_in = 0

    name = None
    seq_parts = []

    def _flush(flush_name, flush_parts):
        # Normalise to match what reformat.sh + sort|uniq does in
        # read_and_count_sequences(): uppercase only.  Dashes are intentionally
        # KEPT here, because read_and_count_sequences() also keeps them
        # (removing them would break subsequence slicing on padded alignments,
        # per the docstring of that function).  Stripping dashes would produce
        # a different sha256 than the one embedded in the counts FASTA ID,
        # making TSV lookups fail for any alignment-padded input file.
        seq = "".join(flush_parts).replace("\r", "").replace("\n", "").upper()
        digest = hashlib.sha256(seq.encode()).hexdigest()
        if digest in mapping:
            mapping[digest][0] += 1
            mapping[digest][1].append(flush_name)
        else:
            mapping[digest] = [1, [flush_name]]

    with open(infilename, "rb") as fh:
        for raw in fh:
            line = _decode_fasta_line(raw).rstrip("\r\n")
            if not line:
                continue
            if line[0] == ">":
                if name is not None:
                    _flush(name, seq_parts)
                    n_in += 1
                    if debug and n_in % 500_000 == 0:
                        print(f"Info: mapping pass: processed {n_in} records,"
                              f" {len(mapping)} unique", file=sys.stderr)
                parts = line[1:].split()
                name = parts[0] if parts else ""
                seq_parts = []
            else:
                seq_parts.append(line)
        if name is not None:
            _flush(name, seq_parts)
            n_in += 1

    print(f"Info: mapping pass: {n_in} total records -> {len(mapping)} unique sequences",
          file=sys.stderr)

    with open(mapping_outfile, "w", encoding="utf-8") as out:
        for digest, (count, ids) in mapping.items():
            out.write("\t".join([digest, str(count)] + ids) + "\n")

    print(f"Info: wrote sha256->IDs mapping ({len(mapping)} entries) to {mapping_outfile}",
          file=sys.stderr)


def main():
    """Parse CLI arguments and run the deduplication pipeline."""
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--infilename", required=True,
                        help="Input FASTA/Q file path.")
    parser.add_argument("--infile-format", dest="infile_format", default="fasta",
                        help="Input FASTA/Q file format [fasta].")
    parser.add_argument("--outfile-prefix", dest="outfile_prefix", default="",
                        help="Output path prefix; '.counts.fasta' is appended.")
    parser.add_argument("--sort-bucket-size", dest="sort_bucket_size", default="30%",
                        help="Sort buffer size for UNIX sort [30%%].")
    parser.add_argument("--min-count", dest="min_count", type=int, default=0,
                        help="Only write sequences with at least this many occurrences.")
    parser.add_argument("--top-n", dest="top_n", type=int, default=0,
                        help="Only write the top N most frequent sequences.")
    parser.add_argument("--mapping-outfile", dest="mapping_outfile", default="",
                        help=(
                            "TSV file mapping sha256 -> original FASTA IDs. "
                            "Columns: sha256hex, count, id_1, id_2, ... "
                            "Defaults to --infilename with FASTA extension removed "
                            "and '.sha256_to_ids.tsv' appended."
                        ))
    parser.add_argument("--overwrite", action="store_true",
                        help="Overwrite output files if they already exist.")
    parser.add_argument("--debug", type=int, default=0,
                        help="Debug verbosity level.")
    parser.add_argument("--version", action="version", version=f"%(prog)s {VERSION}")
    myoptions = parser.parse_args()

    if not os.path.exists(myoptions.infilename):
        parser.error(f"File does not exist: {myoptions.infilename}")

    if myoptions.outfile_prefix and '=' in myoptions.outfile_prefix:
        parser.error(
            f"--outfile-prefix value '{myoptions.outfile_prefix}' contains '=' — "
            "looks like you accidentally included the argument name in the value.\n"
            f"Did you mean: --outfile-prefix={myoptions.outfile_prefix.split('=', 1)[-1]}"
        )
    if myoptions.mapping_outfile and '=' in myoptions.mapping_outfile:
        parser.error(
            f"--mapping-outfile value '{myoptions.mapping_outfile}' contains '=' — "
            "looks like you accidentally included the argument name in the value.\n"
            f"Did you mean: --mapping-outfile={myoptions.mapping_outfile.split('=', 1)[-1]}"
        )

    if not myoptions.outfile_prefix:
        for ext in ('.fastq.gz', '.fastq', '.fasta.gz', '.fasta'):
            if myoptions.infilename.endswith(ext):
                outfile_prefix = myoptions.infilename[:-len(ext)]
                break
        else:
            outfile_prefix = myoptions.infilename
    else:
        outfile_prefix = myoptions.outfile_prefix

    mapping_outfile = myoptions.mapping_outfile or (outfile_prefix + '.sha256_to_ids.tsv')

    counts_outfile = outfile_prefix + '.counts.fasta'
    outputs = [counts_outfile, mapping_outfile]
    input_mtime = os.path.getmtime(myoptions.infilename)

    all_exist = all(os.path.exists(p) for p in outputs)
    if all_exist:
        min_out_mtime = min(os.path.getmtime(p) for p in outputs)
        if min_out_mtime > input_mtime and not myoptions.overwrite:
            print(
                f"Info: outputs are up-to-date (all newer than {myoptions.infilename}), skipping.",
                file=sys.stderr,
            )
            sys.exit(0)
        elif not myoptions.overwrite:
            raise RuntimeError(
                "Output file is stale (older than input): "
                + next(p for p in outputs if os.path.getmtime(p) <= input_mtime)
                + "\nUse --overwrite to regenerate it."
            )
    else:
        for outpath in outputs:
            if os.path.exists(outpath) and not myoptions.overwrite:
                raise RuntimeError(
                    f"Output file already exists: {outpath}\n"
                    "Use --overwrite to replace it."
                )

    with open(counts_outfile, 'w', encoding='utf-8') as outfileh:
        if myoptions.min_count > 0:
            read_and_count_sequences(myoptions.infilename, outfileh,
                                     infile_format=myoptions.infile_format,
                                     min_count=myoptions.min_count,
                                     sort_bucket_size=myoptions.sort_bucket_size)
        elif myoptions.top_n > 0:
            read_and_count_sequences(myoptions.infilename, outfileh,
                                     infile_format=myoptions.infile_format,
                                     top_n=myoptions.top_n,
                                     sort_bucket_size=myoptions.sort_bucket_size)
        else:
            read_and_count_sequences(myoptions.infilename, outfileh,
                                     infile_format=myoptions.infile_format,
                                     sort_bucket_size=myoptions.sort_bucket_size)

    build_sha256_id_mapping(myoptions.infilename, mapping_outfile,
                            debug=myoptions.debug)


if __name__ == "__main__":
    main()
