#!/usr/bin/env python3
"""Expand a FASTA file with deduplicated NNNNx.sha256 IDs back to original FASTA IDs.

Given a FASTA file produced by count_same_sequences.py (whose records have IDs
of the form ``{count}x.{sha256hex}``), this tool reports the original FASTA IDs
that were compacted into each entry.  SHA-256 must be present in every input
record ID; records without it are skipped with a warning.

Usage examples::

    # Mode 1 – list discarded deduplicated IDs
    create_list_of_discarded_sequences.py \\
        --infilename=spikenuc1207.counts.clean.longer_3822.fasta

    # Mode 2 – list original IDs for discarded records via mapping TSV
    create_list_of_discarded_sequences.py \\
        --infilename=spikenuc1207.counts.clean.longer_3822.fasta \\
        --mapping-outfile=spikenuc1207.sha256_to_ids.tsv

    # Mode 3 – full original header lines by re-scanning the source FASTA
    create_list_of_discarded_sequences.py \\
        --infilename=spikenuc1207.counts.clean.longer_3822.fasta \\
        --original-infilename=spikenuc1207.native2ascii.no_junk.fasta

    # Mode 2+3 – original IDs from TSV, full headers from source FASTA
    create_list_of_discarded_sequences.py \\
        --infilename=spikenuc1207.counts.clean.longer_3822.fasta \\
        --mapping-outfile=spikenuc1207.sha256_to_ids.tsv \\
        --original-infilename=spikenuc1207.native2ascii.no_junk.fasta \\
        --outfile=discarded_ids.txt
"""

import hashlib
import os
import sys
from optparse import OptionParser

VERSION = "202603291900"

myparser = OptionParser(version="%s version %s" % ('%prog', VERSION))
myparser.add_option(
    "--infilename", action="store", type="string", dest="infilename", default="",
    help="FASTA file with deduplicated NNNNx.sha256 IDs (e.g. longer_3822.fasta).",
)
myparser.add_option(
    "--mapping-outfile", action="store", type="string", dest="mapping_outfile", default="",
    help=(
        "TSV mapping file produced by count_same_sequences.py --mapping-outfile. "
        "Columns: sha256hex, count, id_1, id_2, ..."
    ),
)
myparser.add_option(
    "--original-infilename", action="store", type="string", dest="original_infilename", default="",
    help=(
        "Original (pre-compaction) FASTA file. When supplied the tool re-scans "
        "it and outputs the full original header lines (everything after '>') "
        "for matching records."
    ),
)
myparser.add_option(
    "--outfile", action="store", type="string", dest="outfile", default="",
    help="Output file path. Defaults to stdout.",
)
myparser.add_option(
    "--debug", action="store", type="int", dest="debug", default=0,
    help="Debug verbosity level [0].",
)
(myoptions, _myargs) = myparser.parse_args()

if not myoptions.infilename:
    myparser.error("--infilename is required")
if not os.path.exists(myoptions.infilename):
    myparser.error("File does not exist: %s" % myoptions.infilename)
if myoptions.mapping_outfile and not os.path.exists(myoptions.mapping_outfile):
    myparser.error("Mapping TSV does not exist: %s" % myoptions.mapping_outfile)
if myoptions.original_infilename and not os.path.exists(myoptions.original_infilename):
    myparser.error("Original FASTA does not exist: %s" % myoptions.original_infilename)

# ── Step 1: collect sha256 hashes from the input FASTA ──────────────────────
# Each record ID has the form  NNNNx.SHA256HEX  (possibly with extra fields).
# We extract the sha256 as the part after the first 'x.'.

def _extract_sha256(record_id):
    """Return the sha256hex from an ID of the form NNNNx.SHA256HEX, or None."""
    xdot = record_id.find('x.')
    if xdot >= 0:
        candidate = record_id[xdot + 2:].split()[0]
        if len(candidate) == 64:
            return candidate
    return None

sha256_to_dedup_id = {}   # sha256 -> the NNNNx.sha256 ID
_skipped = 0

with open(myoptions.infilename, "r", encoding="utf-8", errors="replace") as _fh:
    for _line in _fh:
        _line = _line.rstrip("\r\n")
        if not _line or _line[0] != ">":
            continue
        _parts = _line[1:].split()
        _name = _parts[0] if _parts else ""
        _sha = _extract_sha256(_name)
        if _sha:
            sha256_to_dedup_id[_sha] = _name
        else:
            print("Warning: no sha256 found in ID '%s', skipping" % _name, file=sys.stderr)
            _skipped += 1

if _skipped:
    print("Warning: %d record(s) skipped (no sha256 in ID)" % _skipped, file=sys.stderr)

_target_sha256s = set(sha256_to_dedup_id.keys())
print("Info: %d records in %s" % (len(_target_sha256s), myoptions.infilename), file=sys.stderr)

# ── Step 2: resolve original IDs ────────────────────────────────────────────
# Priority: original-infilename (full headers) > mapping-outfile (IDs only)
# > fallback (deduplicated IDs from the input file).

# lines_to_emit: list of strings to write (without trailing newline)
lines_to_emit = []

if myoptions.original_infilename:
    # Scan the original FASTA; for every record whose sha256 matches one in
    # the target set, emit its full original header (everything after '>').
    _n_scanned = 0
    _n_matched = 0
    _cur_id = None
    _cur_header = None
    _cur_seq_parts = []

    def _flush_orig(cur_id, cur_header, cur_seq_parts):
        global _n_matched
        _seq = "".join(cur_seq_parts)
        _digest = hashlib.sha256(_seq.encode()).hexdigest()
        if _digest in _target_sha256s:
            lines_to_emit.append(cur_header)
            _n_matched += 1

    with open(myoptions.original_infilename, "r", encoding="utf-8", errors="replace") as _fh:
        for _line in _fh:
            _line = _line.rstrip("\r\n")
            if not _line:
                continue
            if _line[0] == ">":
                if _cur_id is not None:
                    _flush_orig(_cur_id, _cur_header, _cur_seq_parts)
                    _n_scanned += 1
                    if myoptions.debug and _n_scanned % 500_000 == 0:
                        print("Info: scanned %d original records, %d matched" % (_n_scanned, _n_matched), file=sys.stderr)
                _cur_header = _line[1:]   # everything after '>'
                _parts = _cur_header.split()
                _cur_id = _parts[0] if _parts else ""
                _cur_seq_parts = []
            else:
                _cur_seq_parts.append(_line)
        if _cur_id is not None:
            _flush_orig(_cur_id, _cur_header, _cur_seq_parts)
            _n_scanned += 1

    print("Info: scanned %d original records, %d matched target sha256s" % (_n_scanned, _n_matched), file=sys.stderr)

elif myoptions.mapping_outfile:
    # Load the TSV and collect original IDs for matching sha256s.
    _n_matched = 0
    with open(myoptions.mapping_outfile, "r", encoding="utf-8") as _fh:
        for _line in _fh:
            _fields = _line.rstrip("\n").split("\t")
            if len(_fields) < 3:
                continue
            _digest, _count = _fields[0], _fields[1]
            if _digest in _target_sha256s:
                _orig_ids = _fields[2:]
                for _orig_id in _orig_ids:
                    lines_to_emit.append(_orig_id)
                _n_matched += len(_orig_ids)
    print("Info: found %d original IDs for %d target sha256s from mapping TSV" % (_n_matched, len(_target_sha256s)), file=sys.stderr)

else:
    # Fallback: print the deduplicated NNNNx.sha256 IDs from --infilename.
    lines_to_emit = list(sha256_to_dedup_id.values())
    print("Info: no mapping TSV or original FASTA supplied — outputting deduplicated IDs", file=sys.stderr)

# ── Step 3: write output ─────────────────────────────────────────────────────
_out = open(myoptions.outfile, "w", encoding="utf-8") if myoptions.outfile else sys.stdout
for _line in lines_to_emit:
    _out.write(_line + "\n")
if myoptions.outfile:
    _out.close()
    print("Info: wrote %d lines to %s" % (len(lines_to_emit), myoptions.outfile), file=sys.stderr)
else:
    print("Info: wrote %d lines to stdout" % len(lines_to_emit), file=sys.stderr)
