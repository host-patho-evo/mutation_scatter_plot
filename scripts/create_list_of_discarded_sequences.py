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

Behaviour
---------
* Output guard: if the output file exists and is *newer* than all input files
  the script prints ``Info: up-to-date, skipping`` and exits 0.  Pass
  ``--overwrite`` to force regeneration even when the output is fresh.
* If the output exists and is *older* than any input the script raises a
  ``RuntimeError`` — add ``--overwrite`` to regenerate.
* Passing ``--outfile=/dev/null`` bypasses the timestamp guard entirely (used
  internally by ``summarize_fasta_pipeline.py`` for inline stats).

Usage examples::

    # Default mode – expand sha256-format discarded file to original IDs (via TSV)
    create_list_of_discarded_sequences.py \\
        --infilename=filename_prefix.counts.clean.longer_3822.fasta \\
        --mapping-outfile=filename_prefix.sha256_to_ids.tsv

    # Default mode – expand any format discarded file by re-scanning original FASTA
    create_list_of_discarded_sequences.py \\
        --infilename=filename_prefix.counts.clean.longer_3822.fasta \\
        --original-infilename=filename_prefix.fasta

    # Inverted mode – given the KEPT file, list what was discarded from original
    create_list_of_discarded_sequences.py \\
        --infilename=filename_prefix.counts.clean.exactly_3822.fasta \\
        --original-infilename=filename_prefix.fasta \\
        --inverted \\
        --outfile=filename_prefix.counts.clean.exactly_3822.discarded_original_ids.txt
"""

import argparse
import hashlib
import os
import sys

VERSION = "202603292130"

_parser = argparse.ArgumentParser(description=__doc__,
                                   formatter_class=argparse.RawDescriptionHelpFormatter)
_parser.add_argument(
    "--infilename", required=True,
    help=(
        "Deduplicated counts FASTA (NNNNx.sha256 or legacy NNNNx format). "
        "In default mode this is the *discarded* file; in --inverted mode it "
        "is the *kept* file."
    ),
)
_parser.add_argument(
    "--mapping-outfile", dest="mapping_outfile", default="",
    help=(
        "TSV mapping file produced by count_same_sequences.py --mapping-outfile. "
        "Columns: sha256hex, count, id_1, id_2, ... "
        "Fast path for default mode when IDs carry sha256."
    ),
)
_parser.add_argument(
    "--original-infilename", dest="original_infilename", default="",
    help=(
        "Original (pre-compaction) FASTA file. Required for --inverted mode and "
        "for default mode when --mapping-outfile is not supplied or IDs lack sha256."
    ),
)
_parser.add_argument(
    "--inverted", action="store_true",
    help=(
        "Invert match: emit original records whose sha256 is NOT present in "
        "--infilename (i.e. --infilename is the kept file, output is discarded). "
        "Requires --original-infilename."
    ),
)
_parser.add_argument(
    "--outfile", default="",
    help="Output file path. Defaults to {infilename_stem}.discarded_original_ids.txt.",
)
_parser.add_argument(
    "--debug", type=int, default=0,
    help="Debug verbosity level [0].",
)
_parser.add_argument(
    "--overwrite", action="store_true",
    help="Overwrite the output file if it already exists.",
)
_parser.add_argument("--version", action="version", version=f"%(prog)s {VERSION}")


# ── helpers ──────────────────────────────────────────────────────────────────

def _extract_sha256(record_id):
    """Return the sha256hex from an ID of the form NNNNx.SHA256HEX, or None."""
    xdot = record_id.find('x.')
    if xdot >= 0:
        candidate = record_id[xdot + 2:].split()[0]
        if len(candidate) == 64:
            return candidate
    return None


def _decode_fasta_line(raw: bytes) -> str:
    """Decode a FASTA line to str, trying UTF-8 first then falling back to
    Latin-1. Handles GISAID headers that mix UTF-8 and Latin-1/Latin-2
    encoded characters in sample descriptions."""
    try:
        return raw.decode("utf-8")
    except UnicodeDecodeError:
        return raw.decode("latin-1")


def _iter_fasta(path):
    """Yield (name, full_header, sequence) for each record in a FASTA file.

    Opened in binary mode with UTF-8 → Latin-1 fallback per line so that
    non-ASCII characters in GISAID sample descriptions are converted to proper
    Unicode rather than replaced with U+FFFD.
    """
    name = full_header = None
    parts = []
    with open(path, "rb") as fh:
        for raw in fh:
            line = _decode_fasta_line(raw).rstrip("\r\n")
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    yield name, full_header, "".join(parts)
                header = line[1:]
                toks = header.split()
                name = toks[0] if toks else ""
                full_header = header
                parts = []
            else:
                parts.append(line)
        if name is not None:
            yield name, full_header, "".join(parts)


def _line_count(line):
    """Extract the integer count prefix from a NNNNx or NNNNx.sha256 first word.
    Returns 1 for plain IDs without such a prefix (e.g. original GISAID IDs).
    """
    first_word = line.split()[0] if line.split() else ""
    x_pos = first_word.find('x')
    if x_pos > 0 and first_word[:x_pos].isdigit():
        return int(first_word[:x_pos])
    return 1


def main():
    """Parse arguments and emit original FASTA IDs matching a deduplicated FASTA."""
    myoptions = _parser.parse_args()

    if not os.path.exists(myoptions.infilename):
        _parser.error(f"File does not exist: {myoptions.infilename}")
    if myoptions.mapping_outfile and not os.path.exists(myoptions.mapping_outfile):
        _parser.error(f"Mapping TSV does not exist: {myoptions.mapping_outfile}")
    if myoptions.original_infilename and not os.path.exists(myoptions.original_infilename):
        _parser.error(f"Original FASTA does not exist: {myoptions.original_infilename}")
    if myoptions.inverted and not myoptions.original_infilename and not myoptions.mapping_outfile:
        _parser.error("--inverted requires --original-infilename or --mapping-outfile")

    # Derive the stem from --infilename by stripping backup suffix + FASTA extension.
    infile_stem = myoptions.infilename
    for ext in ('.old', '.ori', '.orig', '.bak', '.backup'):
        if infile_stem.endswith(ext):
            infile_stem = infile_stem[:-len(ext)]
            break
    for ext in ('.fasta.gz', '.fastq.gz', '.fasta', '.fastq', '.fa', '.fq'):
        if infile_stem.endswith(ext):
            infile_stem = infile_stem[:-len(ext)]
            break

    # Auto-detect mapping TSV when not provided.
    if not myoptions.mapping_outfile:
        guessed_mapping = infile_stem + '.sha256_to_ids.tsv'
        if os.path.exists(guessed_mapping):
            myoptions.mapping_outfile = guessed_mapping
            print(f"Info: auto-detected mapping TSV: {guessed_mapping}", file=sys.stderr)

    # Default --outfile.
    if not myoptions.outfile:
        myoptions.outfile = infile_stem + '.discarded_original_ids.txt'
        print(f"Info: output will be written to {myoptions.outfile}", file=sys.stderr)

    # ── timestamp-aware output guard ─────────────────────────────────────────
    input_files = [myoptions.infilename]
    if myoptions.original_infilename:
        input_files.append(myoptions.original_infilename)
    if myoptions.mapping_outfile:
        input_files.append(myoptions.mapping_outfile)
    input_mtime_max = max(os.path.getmtime(f) for f in input_files)

    if myoptions.outfile != '/dev/null' and os.path.exists(myoptions.outfile):
        out_mtime = os.path.getmtime(myoptions.outfile)
        if out_mtime > input_mtime_max and not myoptions.overwrite:
            print(f"Info: output is up-to-date, skipping: {myoptions.outfile}",
                  file=sys.stderr)
            sys.exit(0)
        if out_mtime <= input_mtime_max and not myoptions.overwrite:
            raise RuntimeError(
                f"Output is stale (older than one or more inputs): {myoptions.outfile}\n"
                "Use --overwrite to regenerate it."
            )
        # else: --overwrite set, proceed

    # ── Step 1: build sha256 set from --infilename ────────────────────────────
    infile_sha256s = {}    # sha256 -> dedup_id (NNNNx or NNNNx.sha256)
    ids_computed = 0
    infile_total_count = 0

    for rec_name, _header, rec_seq in _iter_fasta(myoptions.infilename):
        sha = _extract_sha256(rec_name)
        if sha is None:
            sha = hashlib.sha256(rec_seq.replace('-', '').upper().encode()).hexdigest()
            ids_computed += 1
        infile_sha256s[sha] = rec_name
        x_pos = rec_name.find('x')
        if x_pos > 0 and rec_name[:x_pos].isdigit():
            infile_total_count += int(rec_name[:x_pos])
        else:
            infile_total_count += 1

    if ids_computed:
        print(
            f"Info: {ids_computed:,} input records ({len(infile_sha256s):,} unique sequences)"
            " in --infilename lacked sha256 — computed from sequence content.",
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
    print(
        f"Info: {len(infile_sha256s):,} unique sequences"
        f" (total count: {infile_total_count:,}) in {myoptions.infilename}",
        file=sys.stderr,
    )

    target_sha256s = set(infile_sha256s.keys())

    # ── Step 2: resolve / filter IDs ─────────────────────────────────────────
    lines_to_emit = []

    if myoptions.inverted:
        if myoptions.original_infilename:
            n_scanned = n_emitted = 0
            for rec_name, rec_header, rec_seq in _iter_fasta(myoptions.original_infilename):
                sha = hashlib.sha256(rec_seq.upper().encode()).hexdigest()
                n_scanned += 1
                if sha not in target_sha256s:
                    lines_to_emit.append(rec_header)
                    n_emitted += 1
                if myoptions.debug and n_scanned % 500_000 == 0:
                    print(f"Info: scanned {n_scanned:,}, emitted {n_emitted:,}",
                          file=sys.stderr)
            print(
                f"Info: scanned {n_scanned:,} original records,"
                f" {n_emitted:,} not in infilename (discarded)",
                file=sys.stderr,
            )
        else:
            n_emitted = n_emitted_total = 0
            with open(myoptions.mapping_outfile, "r", encoding="utf-8") as fh:
                for tsv_line in fh:
                    fields = tsv_line.rstrip("\n").split("\t")
                    if len(fields) < 3:
                        continue
                    digest = fields[0]
                    if digest not in target_sha256s:
                        for orig_id in fields[2:]:
                            lines_to_emit.append(orig_id)
                        n_emitted += len(fields) - 2
                        try:
                            n_emitted_total += int(fields[1])
                        except (ValueError, IndexError):
                            n_emitted_total += len(fields) - 2
            print(
                f"Info: found {n_emitted:,} discarded original IDs"
                f" (total sequence count: {n_emitted_total:,}) via mapping TSV",
                file=sys.stderr,
            )

    elif myoptions.original_infilename:
        n_scanned = n_matched = 0
        for rec_name, rec_header, rec_seq in _iter_fasta(myoptions.original_infilename):
            sha = hashlib.sha256(rec_seq.upper().encode()).hexdigest()
            n_scanned += 1
            if sha in target_sha256s:
                lines_to_emit.append(rec_header)
                n_matched += 1
            if myoptions.debug and n_scanned % 500_000 == 0:
                print(f"Info: scanned {n_scanned:,}, matched {n_matched:,}", file=sys.stderr)
        print(
            f"Info: scanned {n_scanned:,} original records,"
            f" {n_matched:,} matched infilename sha256s",
            file=sys.stderr,
        )

    elif myoptions.mapping_outfile and not ids_computed:
        n_matched = n_matched_total = 0
        with open(myoptions.mapping_outfile, "r", encoding="utf-8") as fh:
            for tsv_line in fh:
                fields = tsv_line.rstrip("\n").split("\t")
                if len(fields) < 3:
                    continue
                digest = fields[0]
                if digest in target_sha256s:
                    for orig_id in fields[2:]:
                        lines_to_emit.append(orig_id)
                    n_matched += len(fields) - 2
                    try:
                        n_matched_total += int(fields[1])
                    except (ValueError, IndexError):
                        n_matched_total += len(fields) - 2
        print(
            f"Info: found {n_matched:,} original IDs"
            f" (total sequence count: {n_matched_total:,}) via mapping TSV",
            file=sys.stderr,
        )

    else:
        lines_to_emit = list(infile_sha256s.values())
        print("Info: outputting deduplicated IDs from --infilename", file=sys.stderr)

    # ── Step 3: write output ──────────────────────────────────────────────────
    total_count = sum(_line_count(ln) for ln in lines_to_emit)

    if myoptions.outfile:
        with open(myoptions.outfile, "w", encoding="utf-8") as out:
            for line in lines_to_emit:
                out.write(line + "\n")
        print(
            f"Info: wrote {len(lines_to_emit):,} FASTA ID tags"
            f" (total sum of the count values in the FASTA IDs: {total_count:,})"
            f" to {myoptions.outfile}",
            file=sys.stderr,
        )
    else:
        for line in lines_to_emit:
            sys.stdout.write(line + "\n")
        print(
            f"Info: wrote {len(lines_to_emit):,} FASTA ID tags"
            f" (total sum of the count values in the FASTA IDs: {total_count:,})"
            " to stdout",
            file=sys.stderr,
        )


if __name__ == "__main__":
    main()
