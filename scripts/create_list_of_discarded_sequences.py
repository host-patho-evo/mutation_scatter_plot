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

VERSION = "202603311728"

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
    help=(
        "Output path for sha256-hash entries (NNNNx.sha256hex lines from --infilename). "
        "Defaults to {infilename_stem}.discarded_sha256_hashes.txt. "
        "A companion .discarded_original_ids.txt is always written alongside it "
        "when --mapping-outfile or --original-infilename is available."
    ),
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

    # Default --outfile: now the sha256-hashes file.
    if not myoptions.outfile:
        myoptions.outfile = infile_stem + '.discarded_sha256_hashes.txt'
        print(f"Info: sha256 hashes will be written to {myoptions.outfile}", file=sys.stderr)

    # Companion original-IDs file (always alongside the sha256 hashes file).
    if myoptions.outfile == '/dev/null':
        _original_ids_outfile = '/dev/null'
    else:
        _original_ids_outfile = infile_stem + '.discarded_original_ids.txt'
        print(f"Info: original IDs will be written to {_original_ids_outfile}", file=sys.stderr)

    # ── timestamp-aware output guard ────────────────────────────────────────────
    input_files = [myoptions.infilename]
    if myoptions.original_infilename:
        input_files.append(myoptions.original_infilename)
    if myoptions.mapping_outfile:
        input_files.append(myoptions.mapping_outfile)
    input_mtime_max = max(os.path.getmtime(f) for f in input_files)

    for _outpath in (myoptions.outfile, _original_ids_outfile):
        if _outpath == '/dev/null' or not os.path.exists(_outpath):
            continue
        _out_mtime = os.path.getmtime(_outpath)
        if _out_mtime > input_mtime_max and not myoptions.overwrite:
            print(f"Info: output is up-to-date, skipping: {_outpath}", file=sys.stderr)
            sys.exit(0)
        if _out_mtime <= input_mtime_max and not myoptions.overwrite:
            raise RuntimeError(
                f"Output is stale (older than one or more inputs): {_outpath}\n"
                "Use --overwrite to regenerate it."
            )

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
                "Tip: to avoid this, regenerate --infilename with modern sha256 IDs:\n"
                "  count_same_sequences.py --infilename=<original.fasta>"
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

    # ── Step 2: build sha256_lines and original_id_lines ──────────────────────────
    #
    # sha256_lines    → *.discarded_sha256_hashes.txt
    #   One NNNNx.sha256hex entry per unique sequence.  In default (non-inverted)
    #   mode these come directly from infile_sha256s.values().  In inverted mode
    #   they are reconstructed from the mapping TSV or omitted (FASTA scan gives
    #   no count prefix).
    #
    # original_id_lines → *.discarded_original_ids.txt
    #   One original FASTA ID per individual sequence (count expanded).  Total
    #   lines must equal the sum of the count prefixes; discrepancies are reported.
    #   Generated only when --mapping-outfile or --original-infilename is given.

    sha256_lines = []        # NNNNx.sha256hex entries for the sha256 hashes file
    original_id_lines = []   # expanded original FASTA IDs for the original IDs file
    expected_original_count = infile_total_count  # sum of counts from --infilename
    actual_original_count = 0

    if myoptions.inverted:
        # In inverted mode --infilename is the *kept* set; we want the discarded.
        if myoptions.mapping_outfile:
            n_sha = 0
            with open(myoptions.mapping_outfile, "r", encoding="utf-8") as fh:
                for tsv_line in fh:
                    fields = tsv_line.rstrip("\n").split("\t")
                    if len(fields) < 3:
                        continue
                    digest = fields[0]
                    if digest not in target_sha256s:
                        try:
                            count = int(fields[1])
                        except (ValueError, IndexError):
                            count = len(fields) - 2
                        orig_ids = fields[2:]
                        sha256_lines.append(f"{count}x.{digest}")
                        original_id_lines.extend(orig_ids)
                        actual_original_count += len(orig_ids)
                        if len(orig_ids) != count:
                            print(
                                f"Warning: sha256 {digest[:16]}...: count {count:,} in mapping"
                                f" but {len(orig_ids):,} IDs stored —"
                                f" discrepancy of {abs(count - len(orig_ids)):,}",
                                file=sys.stderr,
                            )
                        n_sha += 1
            # expected_original_count is unknwon for inverted mode without pre-scan;
            # set it to actual so the final check is meaningful only if a FASTA
            # scan path was used.
            expected_original_count = actual_original_count
            print(
                f"Info: {n_sha:,} discarded sha256 entries"
                f" ({actual_original_count:,} original IDs) via mapping TSV (inverted)",
                file=sys.stderr,
            )
        elif myoptions.original_infilename:
            n_scanned = 0
            for _rec_name, rec_header, rec_seq in _iter_fasta(myoptions.original_infilename):
                sha = hashlib.sha256(rec_seq.upper().encode()).hexdigest()
                n_scanned += 1
                if sha not in target_sha256s:
                    original_id_lines.append(rec_header)
                    actual_original_count += 1
                if myoptions.debug and n_scanned % 500_000 == 0:
                    print(f"Info: scanned {n_scanned:,}, emitted {actual_original_count:,}",
                          file=sys.stderr)
            expected_original_count = actual_original_count  # no count prefix available
            print(
                f"Info: scanned {n_scanned:,} original records,"
                f" {actual_original_count:,} not in infilename (discarded, inverted)",
                file=sys.stderr,
            )
        else:
            print(
                "Warning: --inverted mode without --mapping-outfile or"
                " --original-infilename — no output can be generated",
                file=sys.stderr,
            )

    else:
        # Default (non-inverted) mode: sha256_lines = dedup IDs from --infilename.
        sha256_lines = list(infile_sha256s.values())

        if myoptions.mapping_outfile and not ids_computed:
            # Fast path: expand via pre-built mapping TSV.
            with open(myoptions.mapping_outfile, "r", encoding="utf-8") as fh:
                for tsv_line in fh:
                    fields = tsv_line.rstrip("\n").split("\t")
                    if len(fields) < 3:
                        continue
                    digest = fields[0]
                    if digest in target_sha256s:
                        try:
                            count = int(fields[1])
                        except (ValueError, IndexError):
                            count = len(fields) - 2
                        orig_ids = fields[2:]
                        original_id_lines.extend(orig_ids)
                        actual_original_count += len(orig_ids)
                        if len(orig_ids) != count:
                            print(
                                f"Warning: sha256 {digest[:16]}...: expected {count:,} IDs"
                                f" but mapping has {len(orig_ids):,} —"
                                f" discrepancy of {abs(count - len(orig_ids)):,}",
                                file=sys.stderr,
                            )
            print(
                f"Info: found {actual_original_count:,} original IDs"
                f" (expected {expected_original_count:,}) via mapping TSV",
                file=sys.stderr,
            )

        elif myoptions.original_infilename:
            # Slow path: scan original FASTA and match by sha256.
            n_scanned = 0
            sha256_hit_counts: dict = {}  # sha256 → how many times seen in original
            for _rec_name, rec_header, rec_seq in _iter_fasta(myoptions.original_infilename):
                sha = hashlib.sha256(rec_seq.upper().encode()).hexdigest()
                n_scanned += 1
                if sha in target_sha256s:
                    original_id_lines.append(rec_header)
                    actual_original_count += 1
                    sha256_hit_counts[sha] = sha256_hit_counts.get(sha, 0) + 1
                if myoptions.debug and n_scanned % 500_000 == 0:
                    print(f"Info: scanned {n_scanned:,}, matched {actual_original_count:,}",
                          file=sys.stderr)
            # Per-sha256 discrepancy check.
            for sha, rec_name in infile_sha256s.items():
                x_pos = rec_name.find('x')
                exp = int(rec_name[:x_pos]) if (x_pos > 0 and rec_name[:x_pos].isdigit()) else 1
                got = sha256_hit_counts.get(sha, 0)
                if got != exp:
                    print(
                        f"Warning: {rec_name}: expected {exp:,} occurrences in"
                        f" {myoptions.original_infilename} but found {got:,}"
                        f" — discrepancy of {abs(exp - got):,}",
                        file=sys.stderr,
                    )
            print(
                f"Info: scanned {n_scanned:,} original records,"
                f" {actual_original_count:,} matched (expected {expected_original_count:,})",
                file=sys.stderr,
            )

        else:
            print(
                "Info: no --mapping-outfile or --original-infilename provided;"
                f" {_original_ids_outfile} will not be generated",
                file=sys.stderr,
            )

    # Overall discrepancy report for the original IDs file.
    if original_id_lines:
        if actual_original_count != expected_original_count:
            print(
                f"Warning: {actual_original_count:,} original IDs found but expected"
                f" {expected_original_count:,} (sum of count prefixes in --infilename);"
                f" total discrepancy of {abs(expected_original_count - actual_original_count):,}",
                file=sys.stderr,
            )
        else:
            print(
                f"Info: original ID count ({actual_original_count:,}) matches"
                f" expected total ({expected_original_count:,}) — no discrepancy",
                file=sys.stderr,
            )

    # ── Step 3: write sha256 hashes file ──────────────────────────────────────────
    total_count = sum(_line_count(ln) for ln in sha256_lines)

    if myoptions.outfile == '/dev/null':
        print(
            f"Info: {len(sha256_lines):,} sha256 hash entries"
            f" (total count: {total_count:,}) — not written (stats-only mode)",
            file=sys.stderr,
        )
    elif myoptions.outfile:
        with open(myoptions.outfile, "w", encoding="utf-8") as out:
            for line in sha256_lines:
                out.write(line + "\n")
        print(
            f"Info: wrote {len(sha256_lines):,} sha256 hash entries"
            f" (total count: {total_count:,}) to {myoptions.outfile}",
            file=sys.stderr,
        )
    else:
        for line in sha256_lines:
            sys.stdout.write(line + "\n")
        print(
            f"Info: wrote {len(sha256_lines):,} sha256 hash entries"
            f" (total count: {total_count:,}) to stdout",
            file=sys.stderr,
        )

    # ── Step 4: write original IDs file ─────────────────────────────────────────
    if original_id_lines and _original_ids_outfile != '/dev/null':
        with open(_original_ids_outfile, "w", encoding="utf-8") as out:
            for line in original_id_lines:
                out.write(line + "\n")
        print(
            f"Info: wrote {len(original_id_lines):,} original FASTA IDs"
            f" to {_original_ids_outfile}",
            file=sys.stderr,
        )
    elif not original_id_lines and _original_ids_outfile != '/dev/null':
        print(
            f"Info: no original IDs resolved — {_original_ids_outfile} not written",
            file=sys.stderr,
        )


if __name__ == "__main__":
    main()
