#!/usr/bin/env python3
"""List original FASTA IDs corresponding to (or absent from) a deduplicated FASTA.

The input ``--infilename`` is a deduplicated counts FASTA produced by
``count_same_sequences.py``.  Record IDs may be in the modern format
``{count}x.{sha256hex}`` or the legacy format ``{count}x`` (no sha256 in ID).
In both cases sha256 is derived from sequence content when not present in the ID.

**Default mode** (``--infilename`` = a file of *discarded* sequences):
  Scan ``--original-infilename``; emit records whose sha256 **IS** present in
  ``--infilename``.  I.e. "which original IDs were compacted into the discarded
  entries?"

**Inverted mode** (``--inverted``, ``--infilename`` = the *kept* sequences):
  Scan ``--original-infilename``; emit records whose sha256 is **NOT** present
  in ``--infilename``.  I.e. "which original IDs were dropped / discarded?"

The mapping TSV (``--mapping-outfile``) provides a fast path for the default
mode when IDs already carry sha256 *and* the TSV was pre-built by
``count_same_sequences.py --mapping-outfile``.

Usage examples::

    # Default mode – expand sha256-format discarded file to original IDs (via TSV)
    create_list_of_discarded_sequences.py \\
        --infilename=spikenuc1207.counts.clean.longer_3822.fasta \\
        --mapping-outfile=spikenuc1207.sha256_to_ids.tsv

    # Default mode – expand any format discarded file by re-scanning original FASTA
    create_list_of_discarded_sequences.py \\
        --infilename=spikenuc1207.counts.clean.longer_3822.fasta \\
        --original-infilename=spikenuc1207.native2ascii.no_junk.fasta

    # Inverted mode – given the KEPT file, list what was discarded from original
    create_list_of_discarded_sequences.py \\
        --infilename=spikenuc1207.counts.clean.filtered.fasta.old \\
        --original-infilename=spikenuc1207.native2ascii.no_junk.fasta \\
        --inverted \\
        --outfile=discarded_original_ids.txt
"""

import hashlib
import os
import sys
from optparse import OptionParser

VERSION = "202603292000"

myparser = OptionParser(version="%s version %s" % ('%prog', VERSION))
myparser.add_option(
    "--infilename", action="store", type="string", dest="infilename", default="",
    help=(
        "Deduplicated counts FASTA (NNNNx.sha256 or legacy NNNNx format). "
        "In default mode this is the *discarded* file; in --inverted mode it "
        "is the *kept* file."
    ),
)
myparser.add_option(
    "--mapping-outfile", action="store", type="string", dest="mapping_outfile", default="",
    help=(
        "TSV mapping file produced by count_same_sequences.py --mapping-outfile. "
        "Columns: sha256hex, count, id_1, id_2, ... "
        "Fast path for default mode when IDs carry sha256."
    ),
)
myparser.add_option(
    "--original-infilename", action="store", type="string", dest="original_infilename", default="",
    help=(
        "Original (pre-compaction) FASTA file. Required for --inverted mode and "
        "for default mode when --mapping-outfile is not supplied or IDs lack sha256."
    ),
)
myparser.add_option(
    "--inverted", action="store_true", dest="inverted", default=False,
    help=(
        "Invert match: emit original records whose sha256 is NOT present in "
        "--infilename (i.e. --infilename is the kept file, output is discarded). "
        "Requires --original-infilename."
    ),
)
myparser.add_option(
    "--outfile", action="store", type="string", dest="outfile", default="",
    help="Output file path. Defaults to {infilename_stem}.discarded_original_ids.txt.",
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
if myoptions.inverted and not myoptions.original_infilename:
    myparser.error("--inverted requires --original-infilename")

# Derive the stem from --infilename by stripping backup suffix + FASTA extension.
# Used for auto-detecting the mapping TSV and defaulting --outfile.
_infile_stem = myoptions.infilename
for _ext in ('.old', '.ori', '.orig', '.bak', '.backup'):
    if _infile_stem.endswith(_ext):
        _infile_stem = _infile_stem[:-len(_ext)]
        break
for _ext in ('.fasta.gz', '.fastq.gz', '.fasta', '.fastq', '.fa', '.fq'):
    if _infile_stem.endswith(_ext):
        _infile_stem = _infile_stem[:-len(_ext)]
        break

# Auto-detect mapping TSV when not provided.
if not myoptions.mapping_outfile:
    _guessed_mapping = _infile_stem + '.sha256_to_ids.tsv'
    if os.path.exists(_guessed_mapping):
        myoptions.mapping_outfile = _guessed_mapping
        print("Info: auto-detected mapping TSV: %s" % _guessed_mapping, file=sys.stderr)

# Default --outfile to {stem}.discarded_original_ids.txt when not given.
if not myoptions.outfile:
    myoptions.outfile = _infile_stem + '.discarded_original_ids.txt'
    print("Info: output will be written to %s" % myoptions.outfile, file=sys.stderr)



# ── helpers ──────────────────────────────────────────────────────────────────

def _extract_sha256(record_id):
    """Return the sha256hex from an ID of the form NNNNx.SHA256HEX, or None."""
    xdot = record_id.find('x.')
    if xdot >= 0:
        candidate = record_id[xdot + 2:].split()[0]
        if len(candidate) == 64:
            return candidate
    return None


def _iter_fasta(path):
    """Yield (name, full_header, seq) triples from a FASTA file.

    *name*        – first word after '>'
    *full_header* – everything after '>' (for output purposes)
    *seq*         – joined sequence string, \\r stripped
    """
    name = full_header = None
    parts = []
    with open(path, "r", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            line = line.rstrip("\r\n")
            if not line:
                continue
            if line[0] == ">":
                if name is not None:
                    yield name, full_header, "".join(parts)
                full_header = line[1:]
                ws = full_header.split()
                name = ws[0] if ws else ""
                parts = []
            else:
                parts.append(line.replace("\r", ""))
    if name is not None:
        yield name, full_header, "".join(parts)


# ── Step 1: build sha256 set from --infilename ────────────────────────────────
# Always read sequences so we can compute sha256 when the ID lacks it.

infile_sha256s = {}    # sha256 -> dedup_id (NNNNx or NNNNx.sha256)
_ids_computed = 0      # how many sha256s were computed (not extracted from ID)

for _name, _header, _seq in _iter_fasta(myoptions.infilename):
    _sha = _extract_sha256(_name)
    if _sha is None:
        # Legacy NNNNx format: compute sha256 from sequence content.
        _sha = hashlib.sha256(_seq.encode()).hexdigest()
        _ids_computed += 1
    infile_sha256s[_sha] = _name

if _ids_computed:
    print(
        "Info: %d/%d infilename IDs lacked sha256 — computed from sequence content."
        % (_ids_computed, len(infile_sha256s)),
        file=sys.stderr,
    )
    if not myoptions.mapping_outfile:
        print(
            "Tip: to avoid this, regenerate --infilename in modern format by running:\n"
            "  count_same_sequences.py --infilename=<original_pre-compaction.fasta>"
            " --outfile-prefix=<prefix>\n"
            "This produces NNNNx.sha256 IDs and a .sha256_to_ids.tsv mapping.",
            file=sys.stderr,
        )
print("Info: %d unique sequences in %s" % (len(infile_sha256s), myoptions.infilename), file=sys.stderr)

_target_sha256s = set(infile_sha256s.keys())

# ── Step 2: resolve / filter IDs ─────────────────────────────────────────────

lines_to_emit = []

if myoptions.inverted:
    # Emit originals NOT in infile (infile = kept, output = discarded).
    _n_scanned = _n_emitted = 0
    for _name, _header, _seq in _iter_fasta(myoptions.original_infilename):
        _sha = hashlib.sha256(_seq.encode()).hexdigest()
        _n_scanned += 1
        if _sha not in _target_sha256s:
            lines_to_emit.append(_header)
            _n_emitted += 1
        if myoptions.debug and _n_scanned % 500_000 == 0:
            print("Info: scanned %d, emitted %d" % (_n_scanned, _n_emitted), file=sys.stderr)
    print(
        "Info: scanned %d original records, %d not in infilename (discarded)"
        % (_n_scanned, _n_emitted),
        file=sys.stderr,
    )

elif myoptions.original_infilename:
    # Default mode: emit originals whose sha256 IS in infile.
    _n_scanned = _n_matched = 0
    for _name, _header, _seq in _iter_fasta(myoptions.original_infilename):
        _sha = hashlib.sha256(_seq.encode()).hexdigest()
        _n_scanned += 1
        if _sha in _target_sha256s:
            lines_to_emit.append(_header)
            _n_matched += 1
        if myoptions.debug and _n_scanned % 500_000 == 0:
            print("Info: scanned %d, matched %d" % (_n_scanned, _n_matched), file=sys.stderr)
    print(
        "Info: scanned %d original records, %d matched infilename sha256s"
        % (_n_scanned, _n_matched),
        file=sys.stderr,
    )

elif myoptions.mapping_outfile and not _ids_computed:
    # Fast path: all IDs had sha256 and we have a pre-built TSV.
    _n_matched = 0
    with open(myoptions.mapping_outfile, "r", encoding="utf-8") as _fh:
        for _line in _fh:
            _fields = _line.rstrip("\n").split("\t")
            if len(_fields) < 3:
                continue
            _digest = _fields[0]
            if _digest in _target_sha256s:
                for _orig_id in _fields[2:]:
                    lines_to_emit.append(_orig_id)
                _n_matched += len(_fields) - 2
    print(
        "Info: found %d original IDs via mapping TSV" % _n_matched,
        file=sys.stderr,
    )

else:
    # Fallback: just output the deduplicated IDs from --infilename.
    lines_to_emit = list(infile_sha256s.values())
    print("Info: outputting deduplicated IDs from --infilename", file=sys.stderr)

# ── Step 3: write output ──────────────────────────────────────────────────────

def _line_count(line):
    """Extract the integer count prefix from a NNNNx or NNNNx.sha256 first word.
    Returns 1 for plain IDs without such a prefix (e.g. original GISAID IDs).
    """
    first_word = line.split()[0] if line.split() else ""
    x_pos = first_word.find('x')
    if x_pos > 0 and first_word[:x_pos].isdigit():
        return int(first_word[:x_pos])
    return 1

_total_count = sum(_line_count(l) for l in lines_to_emit)

_out = open(myoptions.outfile, "w", encoding="utf-8") if myoptions.outfile else sys.stdout
for _line in lines_to_emit:
    _out.write(_line + "\n")
if myoptions.outfile:
    _out.close()
    print(
        "Info: wrote %d lines (total sequence count: %d) to %s"
        % (len(lines_to_emit), _total_count, myoptions.outfile),
        file=sys.stderr,
    )
else:
    print(
        "Info: wrote %d lines (total sequence count: %d) to stdout"
        % (len(lines_to_emit), _total_count),
        file=sys.stderr,
    )
