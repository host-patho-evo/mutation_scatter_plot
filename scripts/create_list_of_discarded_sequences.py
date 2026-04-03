#!/usr/bin/env python3
"""List original FASTA IDs corresponding to (or absent from) a deduplicated FASTA.

The input ``--infilename`` is a deduplicated counts FASTA produced by
``count_same_sequences.py``.  Record IDs may be in the modern format
``{count}x.{sha256hex}`` or the legacy format ``{count}x`` (no sha256 in ID).
In both cases sha256 is derived from sequence content when not present in the ID.

**Default mode** (``--infilename`` = a file of *discarded* sequences):
  Scan ``--original-infilename`` or ``--mapping-outfile``; emit records whose
  sha256 **IS** present in ``--infilename``.  I.e. "which original IDs were
  compacted into the discarded entries?"

**Inverted mode** (``--inverted``, ``--infilename`` = the *kept* sequences):
  Scan ``--original-infilename`` or ``--mapping-outfile``; emit records whose
  sha256 is **NOT** present in ``--infilename``.  I.e. "which original IDs
  were dropped / discarded?"

Output files
------------
Two files are always written (unless ``--outfile=/dev/null``):

``{stem}.discarded_sha256_hashes.txt``  (controlled by ``--outfile``)
    One ``NNNNx.sha256hex`` entry per *unique* discarded sequence.  The sum
    of all ``NNNNx`` prefixes equals the total number of individual original
    sequences that were discarded.  Used as the fast-read cache by
    ``summarize_fasta_pipeline.py``.

``{stem}.discarded_original_ids.txt``  (auto-generated companion)
    One original FASTA header per *individual* discarded sequence
    (count-expanded).  ``wc -l`` of this file should equal the sum of all
    ``NNNNx`` prefixes in the sha256-hashes file when there are no
    discrepancies.  Generated only when ``--mapping-outfile`` or
    ``--original-infilename`` is provided (auto-detected from the stem when
    omitted).

The mapping TSV (``--mapping-outfile``, typically ``*.sha256_to_ids.tsv``)
provides a fast path without needing the original FASTA: it maps each sha256
to the list of original FASTA IDs and is auto-detected when not provided.

Behaviour
---------
* Output guard: if the output files exist and are *newer* than all inputs the
  script prints ``Info: up-to-date, skipping`` and exits 0.  Pass
  ``--overwrite`` to force regeneration.
* If an output is *older* than any input the script raises ``RuntimeError``
  — add ``--overwrite`` to regenerate.
* ``--outfile=/dev/null`` bypasses the timestamp guard entirely (stats-only
  mode used internally by ``summarize_fasta_pipeline.py``).
* Legacy IDs (no sha256 in ID): sha256 is computed from sequence content and
  a tip to regenerate with modern NNNNx.sha256hex IDs is printed.

Usage examples::

    # Default mode – sha256-format discarded file; TSV auto-detected
    create_list_of_discarded_sequences.py \\
        --infilename=prefix.counts.clean.longer_3822.fasta
    # Writes:
    #   prefix.counts.clean.longer_3822.discarded_sha256_hashes.txt
    #   prefix.counts.clean.longer_3822.discarded_original_ids.txt

    # Default mode – explicit TSV
    create_list_of_discarded_sequences.py \\
        --infilename=prefix.counts.clean.longer_3822.fasta \\
        --mapping-outfile=prefix.sha256_to_ids.tsv

    # Default mode – re-scan from original FASTA (slow)
    create_list_of_discarded_sequences.py \\
        --infilename=prefix.counts.clean.longer_3822.fasta \\
        --original-infilename=prefix.fasta

    # Inverted mode – given the KEPT file, list what was discarded
    create_list_of_discarded_sequences.py \\
        --infilename=prefix.counts.clean.exactly_3822.fasta \\
        --mapping-outfile=prefix.sha256_to_ids.tsv \\
        --inverted
    # Writes:
    #   prefix.counts.clean.exactly_3822.discarded_sha256_hashes.txt
    #   prefix.counts.clean.exactly_3822.discarded_original_ids.txt
    #
    # Verify no discrepancy:
    #   wc -l prefix.counts.clean.exactly_3822.discarded_original_ids.txt
    #   awk '{print $1}' prefix.counts.clean.exactly_3822.discarded_sha256_hashes.txt \\
    #       | sed 's/x.*//' | awk '{s+=$1} END {print s}'
"""

import argparse
import datetime
import hashlib
import os
import re
import subprocess
import sys

VERSION = "202603311815"


def _get_git_version() -> str:
    """Return ``git describe --always --dirty --tags`` output, or ``'unknown'``."""
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
        "--infilename (i.e. --infilename is the kept file, output is what was "
        "discarded).  Requires --original-infilename or --mapping-outfile."
    ),
)
_parser.add_argument(
    "--output-context", dest="output_context",
    choices=["discarded", "effectively_used"], default="discarded",
    help=(
        "Label used in the output filenames to indicate what the file contains. "
        "'discarded' (default): output files are named *.discarded_sha256_hashes.txt "
        "and *.discarded_original_ids.txt. "
        "'effectively_used': output files are named *.effectively_used_sha256_hashes.txt "
        "and *.effectively_used_original_ids.txt. "
        "Set to 'effectively_used' when --infilename is the *kept* file and "
        "--inverted is NOT set (i.e. you are expanding what was kept, not discarded)."
    ),
)
_parser.add_argument(
    "--outfile", default="",
    help=(
        "Output path for sha256-hash entries (NNNNx.sha256hex lines). "
        "Defaults to {infilename_stem}.{output_context}_sha256_hashes.txt "
        "where output_context is 'discarded' or 'effectively_used'. "
        "A companion {output_context}_original_ids.txt is always written alongside "
        "when --mapping-outfile or --original-infilename is available."
    ),
)
_parser.add_argument(
    "--full-fasta-header", dest="full_fasta_header", action="store_true",
    help=(
        "When reading the mapping TSV, prefer *.sha256_to_descr_lines.tsv "
        "(full FASTA header lines, everything after '>') over "
        "*.sha256_to_ids.tsv (first-word ID only).  The descr TSV is "
        "auto-detected from the mapping TSV path when not provided explicitly. "
        "Embedded '\\t' two-character sequences in description fields are "
        "unescaped back to real TAB characters on read."
    ),
)
_parser.add_argument(
    "--outfile-prefix", dest="outfile_prefix", default="",
    help=(
        "Stem (basename) used for all output filenames instead of the one derived "
        "from --infilename.  Combined with --path when provided. "
        "E.g. '--path=/data --outfile-prefix=spikenuc.native2ascii' produces "
        "/data/spikenuc.native2ascii.discarded_sha256_hashes.txt. "
        "Follows the same convention as count_same_sequences.py."
    ),
)
_parser.add_argument(
    "--path", dest="path", default="",
    help=(
        "Directory in which to look for auto-detected ancillary files "
        "(*.sha256_to_ids.tsv) and where output files are written when "
        "--outfile-prefix is a bare filename (no directory part). "
        "Has no effect when --outfile-prefix already contains a directory "
        "component or when --outfile is given explicitly."
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



def _sort_inplace(path: str, sort_args: list) -> None:
    """Sort *path* in-place using the system 'sort' command.

    Uses ``sort -o path path`` which GNU sort guarantees to be safe
    (reads input fully before overwriting the output).
    """
    result = subprocess.run(
        ['sort'] + sort_args + ['-o', path, path],
        capture_output=True, text=True, check=False,
    )
    if result.returncode != 0:
        print(
            f"Warning: sort failed for {path}: {result.stderr.strip()}",
            file=sys.stderr,
        )


def _extract_sha256(record_id):
    """Return the sha256hex from an ID of the form NNNNx.SHA256HEX, or None."""
    xdot = record_id.find('x.')
    if xdot >= 0:
        candidate = record_id[xdot + 2:].split()[0]
        if len(candidate) == 64:
            return candidate
    return None


def _unescape_unicode(s: str) -> str:
    r"""Convert literal \uXXXX escape sequences to real Unicode characters.

    Some GISAID exports encode accented characters as the literal six-character
    sequence \u00e9 rather than the proper UTF-8 glyph.  This function is
    applied after byte-decoding so both forms end up as the same Unicode string.
    The fast path (no '\\u' substring) adds negligible overhead for normal lines.
    """
    if '\\u' not in s:
        return s
    return re.sub(r'\\u([0-9a-fA-F]{4})',
                  lambda m: chr(int(m.group(1), 16)), s)


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
    # Convert surrogates from non-UTF-8 bytes to their Latin-1 equivalents.
    if any(0xDC80 <= ord(ch) <= 0xDCFF for ch in s):
        s = "".join(
            chr(ord(ch) - 0xDC00) if 0xDC80 <= ord(ch) <= 0xDCFF else ch
            for ch in s
        )
    return _unescape_unicode(s)


def _iter_fasta(path):
    """Yield (name, full_header, sequence) for each record in a FASTA file.

    Opened in binary mode with UTF-8 → Latin-1 fallback per line so that
    non-ASCII characters in GISAID sample descriptions are converted to proper
    Unicode rather than replaced with U+FFFD.
    """
    name = full_header = None
    parts = []
    skip_seq = False
    with open(path, "rb") as fh:
        for raw in fh:
            if not raw.rstrip(b"\r\n"):
                continue
            if raw.startswith(b">"):
                if name is not None:
                    yield name, full_header, "".join(parts)
                line = _decode_fasta_line(raw).rstrip("\r\n")
                header = line[1:]
                toks = header.split()
                name = toks[0] if toks else ""
                full_header = header
                parts = []
                skip_seq = _extract_sha256(name) is not None
            else:
                if not skip_seq:
                    line = _decode_fasta_line(raw).rstrip("\r\n")
                    if line:
                        parts.append(line)
                elif not parts:
                    line = raw.decode("ascii", errors="ignore").strip()
                    if line:
                        parts.append(line[0])
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
    _start_ts = datetime.datetime.now().strftime('[%Y-%m-%d %H:%M:%S]')
    print(
        f"{_start_ts} create_list_of_discarded_sequences.py"
        f"  version {VERSION}  git:{_GIT_VERSION}"
        f"  invoked: {' '.join(sys.argv)}",
        file=sys.stderr,
    )
    myoptions = _parser.parse_args()

    if not os.path.exists(myoptions.infilename):
        _parser.error(f"File does not exist: {myoptions.infilename}")
    if myoptions.mapping_outfile and not os.path.exists(myoptions.mapping_outfile):
        _parser.error(f"Mapping TSV does not exist: {myoptions.mapping_outfile}")
    if myoptions.original_infilename and not os.path.exists(myoptions.original_infilename):
        _parser.error(f"Original FASTA does not exist: {myoptions.original_infilename}")
    if myoptions.inverted and not myoptions.original_infilename and not myoptions.mapping_outfile:
        _parser.error("--inverted requires --original-infilename or --mapping-outfile")

    # Guard against common typo --outfile-prefix=--infilename=...
    if myoptions.outfile_prefix and '=' in myoptions.outfile_prefix:
        _parser.error(
            f"--outfile-prefix value '{myoptions.outfile_prefix}' contains '=' — "
            "this looks like a double-flag typo. "
            f"Did you mean: --outfile-prefix={myoptions.outfile_prefix.split('=', 1)[-1]}"
        )

    # Derive the stem:
    #   1. --outfile-prefix given: use it directly (may already include directory).
    #   2. --path given: join it with the basename derived from --infilename.
    #   3. Neither: derive from --infilename as-is (strips FASTA extensions).
    if myoptions.outfile_prefix:
        if myoptions.path and not os.path.isabs(myoptions.outfile_prefix) \
                and not os.path.dirname(myoptions.outfile_prefix):
            # Bare prefix + explicit path -> combine.
            infile_stem = os.path.join(myoptions.path, myoptions.outfile_prefix)
        else:
            infile_stem = myoptions.outfile_prefix
    else:
        # Strip backup + FASTA extensions from --infilename.
        raw = myoptions.infilename
        for ext in ('.old', '.ori', '.orig', '.bak', '.backup'):
            if raw.endswith(ext):
                raw = raw[:-len(ext)]
                break
        for ext in ('.fasta.gz', '.fastq.gz', '.fasta', '.fastq', '.fa', '.fq'):
            if raw.endswith(ext):
                raw = raw[:-len(ext)]
                break
        if myoptions.path and not os.path.dirname(raw):
            # Bare filename stem + explicit path -> combine.
            infile_stem = os.path.join(myoptions.path, os.path.basename(raw))
        else:
            infile_stem = raw

    # Auto-detect mapping TSV when not provided.
    if not myoptions.mapping_outfile:
        guessed_mapping = infile_stem + '.sha256_to_ids.tsv'
        if os.path.exists(guessed_mapping):
            myoptions.mapping_outfile = guessed_mapping
            print(f"Info: auto-detected mapping TSV: {guessed_mapping}", file=sys.stderr)

    # When --full-fasta-header, prefer *.sha256_to_descr_lines.tsv over the
    # ID-only TSV.  Auto-detect it from the mapping TSV path when possible.
    if myoptions.full_fasta_header and myoptions.mapping_outfile:
        if myoptions.mapping_outfile.endswith('.sha256_to_descr_lines.tsv'):
            # summarize_fasta_pipeline.py already provided the descr TSV directly;
            # nothing to switch — suppress the spurious "not found" warning.
            pass
        else:
            descr_candidate = myoptions.mapping_outfile.replace(
                '.sha256_to_ids.tsv', '.sha256_to_descr_lines.tsv'
            )
            if descr_candidate != myoptions.mapping_outfile and os.path.exists(descr_candidate):
                print(
                    f"Info: --full-fasta-header: switching to descr TSV: {descr_candidate}",
                    file=sys.stderr,
                )
                myoptions.mapping_outfile = descr_candidate
            else:
                print(
                    f"Warning: --full-fasta-header: descr TSV not found ({descr_candidate}); "
                    "falling back to ID-only TSV.  Run summarize_fasta_pipeline.py with "
                    "--full-fasta-header to build *.sha256_to_descr_lines.tsv files.",
                    file=sys.stderr,
                )

    # Default --outfile: context-aware sha256-hashes filename.
    _ctx = myoptions.output_context  # 'discarded' or 'effectively_used'
    if not myoptions.outfile:
        myoptions.outfile = infile_stem + f'.{_ctx}_sha256_hashes.txt'
        print(f"Info: sha256 hashes will be written to {myoptions.outfile}", file=sys.stderr)

    # Companion original-IDs file (always alongside the sha256 hashes file).
    if myoptions.outfile == '/dev/null':
        _original_ids_outfile = '/dev/null'
    else:
        _original_ids_outfile = infile_stem + f'.{_ctx}_original_ids.txt'
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
            sha = hashlib.sha256(rec_seq.replace('\r', '').replace('\n', '').replace('-', '').upper().encode()).hexdigest()
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
            _is_descr_tsv = myoptions.mapping_outfile.endswith('.sha256_to_descr_lines.tsv')
            _id_groups: list[tuple[int, list[str]]] = []  # (nnnx_count, orig_ids) — sorted later
            with open(myoptions.mapping_outfile, "r", encoding="utf-8") as fh:
                for tsv_line in fh:
                    fields = tsv_line.rstrip("\n").split("\t")
                    if len(fields) < 3:
                        continue
                    digest = fields[0]
                    if digest not in target_sha256s:
                        try:
                            tsv_count = int(fields[1])
                        except (ValueError, IndexError):
                            tsv_count = len(fields) - 2
                        orig_ids = fields[2:]
                        # Unescape '\t' -> TAB in description fields
                        # (only the descr TSV uses this encoding).
                        if _is_descr_tsv:
                            orig_ids = [s.replace('\\t', '\t') for s in orig_ids]
                        # Use the NNNNx count embedded in each ID (e.g. "17342x.sha256")
                        # rather than fields[1] (how many times the sha256 appears in
                        # the parent FASTA file).  For a deduplicated counts FASTA,
                        # fields[1] is always 1, but the IDs carry the true multiplicity.
                        # For non-deduplicated FASTAs (plain GISAID IDs), _line_count
                        # returns 1 per ID, so the sum equals len(orig_ids) = tsv_count.
                        nnnx_count = sum(_line_count(oid) for oid in orig_ids)
                        sha256_lines.append(f"{nnnx_count}x.{digest}")
                        _id_groups.append((nnnx_count, orig_ids))
                        actual_original_count += nnnx_count
                        if len(orig_ids) != tsv_count:
                            _dir = "fewer" if len(orig_ids) < tsv_count else "more"
                            _tsv_mtime = datetime.datetime.fromtimestamp(
                                os.path.getmtime(myoptions.mapping_outfile)
                            ).strftime('%Y-%m-%d %H:%M:%S')
                            print(
                                f"Warning: sha256 {digest[:16]}...: expected {tsv_count:,} IDs"
                                f" but mapping file {myoptions.mapping_outfile}"
                                f" created on {_tsv_mtime} has {len(orig_ids):,}"
                                f" ({_dir} than expected by {abs(tsv_count - len(orig_ids)):,})",
                                file=sys.stderr,
                            )
                        n_sha += 1
            # Sort groups by NNNNx count descending (mirrors sort -rn on the sha256_hashes
            # file), then flatten into original_id_lines.
            _id_groups.sort(key=lambda g: g[0], reverse=True)
            original_id_lines.extend(oid for _, ids in _id_groups for oid in ids)
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
                sha = hashlib.sha256(rec_seq.replace('\r', '').replace('\n', '').upper().encode()).hexdigest()
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
            _id_groups2: list[tuple[int, list[str]]] = []  # (nnnx_count, orig_ids) — sorted later
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
                        # Use the NNNNx count from the dedup ID in infile_sha256s as the
                        # group sort key (same integer used in sha256_hashes sort -rn).
                        nnnx = _line_count(infile_sha256s.get(digest, ''))
                        _id_groups2.append((nnnx, orig_ids))
                        actual_original_count += len(orig_ids)
                        if len(orig_ids) != count:
                            _dir = "fewer" if len(orig_ids) < count else "more"
                            _tsv_mtime = datetime.datetime.fromtimestamp(
                                os.path.getmtime(myoptions.mapping_outfile)
                            ).strftime('%Y-%m-%d %H:%M:%S')
                            print(
                                f"Warning: sha256 {digest[:16]}...: expected {count:,} IDs"
                                f" but mapping file {myoptions.mapping_outfile}"
                                f" created on {_tsv_mtime} has {len(orig_ids):,}"
                                f" ({_dir} than expected by {abs(count - len(orig_ids)):,})",
                                file=sys.stderr,
                            )
            # Sort groups by NNNNx count descending then flatten.
            _id_groups2.sort(key=lambda g: g[0], reverse=True)
            original_id_lines.extend(oid for _, ids in _id_groups2 for oid in ids)
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
                sha = _extract_sha256(_rec_name)
                if not sha:
                    sha = hashlib.sha256(rec_seq.replace('\r', '').replace('\n', '').upper().encode()).hexdigest()
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
                    _dir = "fewer" if got < exp else "more"
                    _orig_mtime = datetime.datetime.fromtimestamp(
                        os.path.getmtime(myoptions.original_infilename)
                    ).strftime('%Y-%m-%d %H:%M:%S')
                    print(
                        f"Warning: {rec_name}: expected {exp:,} occurrences"
                        f" but original FASTA file {myoptions.original_infilename}"
                        f" created on {_orig_mtime} has {got:,}"
                        f" ({_dir} than expected by {abs(exp - got):,})",
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

    # ── Step 3: write and sort sha256 hashes file ────────────────────────────
    # sha256 hashes file: sort -rn so lines are ordered by the leading
    # integer count (NNNNx.) descending — most-duplicated sequences first.
    # original IDs file: sort -V (version/natural sort) so embedded numbers
    # in FASTA IDs are compared by value, not lexicographically.
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
        # Sort by leading count descending: most-duplicated sequences first.
        _sort_inplace(myoptions.outfile, ['-rn'])
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
        # Data is already ordered by NNNNx group descending (mirrors sort -rn
        # on the sha256_hashes file); no external sort needed for TSV paths.
        # FASTA-scan paths (rare fallback) retain encounter order.
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
