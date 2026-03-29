#!/usr/bin/env python3
"""Deduplicate a FASTA file by sequence content using SHA-256.

For each unique sequence the script emits one output record whose FASTA ID
has the form::

    {count}x.{sha256hex}  [original fields from the first-seen header]

where *count* is the number of input records that share the same sequence and
*sha256hex* is the lowercase hex-digest of the raw (ASCII) sequence.

Optionally a translation table (``--mapping-outfile``) is written so that
every sha256 can be traced back to the full list of original FASTA IDs.
The TSV has the following columns (tab-separated, no header)::

    sha256hex   count   id_1    id_2    …

Usage example::

    count_same_sequences.py \\
        --infile=spikenuc1207.native2ascii.no_junk.fasta \\
        --outfile=spikenuc1207.native2ascii.no_junk.counts.clean.fasta \\
        --mapping-outfile=spikenuc1207.sha256_to_ids.tsv
"""

import hashlib
import sys
from optparse import OptionParser

VERSION = 202603291600

myparser = OptionParser(version="%s version %s" % ('%prog', VERSION))
myparser.add_option(
    "--infile", action="store", type="string", dest="infile", default="",
    help="Input FASTA file (uncompressed)."
)
myparser.add_option(
    "--outfile", action="store", type="string", dest="outfile", default="",
    help="Output deduplicated FASTA file."
)
myparser.add_option(
    "--mapping-outfile", action="store", type="string", dest="mapping_outfile",
    default="",
    help=(
        "Optional TSV file mapping sha256 → original FASTA IDs. "
        "Columns: sha256hex, count, id_1, id_2, … (tab-separated, no header)."
    ),
)
myparser.add_option(
    "--format", action="store", type="string", dest="format", default="fasta",
    help="Input file format understood by BioPython SeqIO [fasta]."
)
myparser.add_option(
    "--debug", action="store", type="int", dest="debug", default=0,
    help="Debug verbosity level [0]."
)
(myoptions, _myargs) = myparser.parse_args()

if not myoptions.infile:
    myparser.error("--infile is required")
if not myoptions.outfile:
    myparser.error("--outfile is required")

# ── Pass 1: group records by SHA-256 of their sequence ──────────────────────
# We keep, per unique sequence:
#   counts[sha256]  = int           (how many records share this sequence)
#   ids[sha256]     = list[str]     (original FASTA IDs, in order seen)
#   first_seq[sha256] = str         (the sequence string, uppercase)
#   first_desc[sha256] = str        (everything after the first word of the
#                                    header, for the representative record)

counts = {}        # sha256 → int
ids = {}           # sha256 → list[str]
first_seq = {}     # sha256 → str
first_desc = {}    # sha256 → str (the extra header fields after the ID word)

_n_in = 0

# Use Bio.SeqIO for robust multi-line FASTA support
from Bio import SeqIO  # pylint: disable=import-outside-toplevel

with open(myoptions.infile, "r", encoding="utf-8") as _fh:
    for _record in SeqIO.parse(_fh, myoptions.format):
        _n_in += 1
        _seq_upper = str(_record.seq).upper().replace("\r", "")
        _digest = hashlib.sha256(_seq_upper.encode("ascii")).hexdigest()
        _orig_id = _record.id
        # _record.description includes the id + the rest of the header line;
        # strip the leading id word to get the extra fields.
        _desc_extra = _record.description[len(_orig_id):].strip()

        if _digest in counts:
            counts[_digest] += 1
            ids[_digest].append(_orig_id)
        else:
            counts[_digest] = 1
            ids[_digest] = [_orig_id]
            first_seq[_digest] = _seq_upper
            first_desc[_digest] = _desc_extra

        if myoptions.debug and _n_in % 100_000 == 0:
            print(
                f"Info: processed {_n_in:,} records, "
                f"{len(counts):,} unique so far",
                file=sys.stderr,
            )

print(
    f"Info: read {_n_in:,} records → {len(counts):,} unique sequences",
    file=sys.stderr,
)

# ── Pass 2: write deduplicated FASTA and optional mapping TSV ───────────────
_n_out = 0
_mapping_fh = (
    open(myoptions.mapping_outfile, "w", encoding="utf-8")
    if myoptions.mapping_outfile
    else None
)

with open(myoptions.outfile, "w", encoding="utf-8") as _out_fh:
    for _digest, _count in counts.items():
        _n_out += 1
        # Build new FASTA ID: "{count}x.{sha256hex}"
        _new_id = f"{_count}x.{_digest}"
        _desc = first_desc[_digest]
        _header = f">{_new_id} {_desc}" if _desc else f">{_new_id}"
        _out_fh.write(f"{_header}\n{first_seq[_digest]}\n")

        if _mapping_fh is not None:
            # sha256 \t count \t id1 \t id2 \t …
            _mapping_fh.write(
                "\t".join([_digest, str(_count)] + ids[_digest]) + "\n"
            )

if _mapping_fh is not None:
    _mapping_fh.close()
    print(f"Info: wrote mapping to {myoptions.mapping_outfile}", file=sys.stderr)

print(
    f"Info: wrote {_n_out:,} unique records to {myoptions.outfile}",
    file=sys.stderr,
)
