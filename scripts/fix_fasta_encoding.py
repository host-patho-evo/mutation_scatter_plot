#!/usr/bin/env python3
"""Normalize the character encoding of FASTA files in-place.

Reads each input file with a two-stage decode:
  1. Bytes → str: UTF-8 attempted first; on UnicodeDecodeError falls back to
     Latin-1 (ISO-8859-1).  This handles GISAID FASTA files that mix UTF-8
     and Latin-1 encoded characters in sample description lines.
  2. Literal escape conversion: ``\\uXXXX`` sequences (produced by some GISAID
     export or native2ascii pipelines) are replaced with the real Unicode
     codepoint (e.g. ``\\u00e9`` → ``é``).

The cleaned file is written as strict UTF-8 and the original is renamed to
``{infile}.orig`` as a backup.  If ``{infile}.orig`` already exists, the file
is skipped to prevent overwriting a previous backup.

The equivalent manual pipeline this replaces::

    # 1. Convert \\uXXXX escapes to real Unicode (Java native2ascii):
    native2ascii -encoding UTF-8 -reverse input.fasta output.fasta

    # 2. Optionally transliterate remaining non-ASCII to ASCII (unidecode):
    unidecode output.fasta > output.fasta.tmp && mv output.fasta.tmp output.fasta

This script only performs step 1 (lossless).  Step 2 (lossy ASCII
transliteration via the ``unidecode`` CLI/library) is not applied here because
it discards information; re-run ``unidecode`` manually if ASCII-only output is
required.

Usage examples::

    # Fix a single file:
    fix_fasta_encoding.py spikenuc1207.native2ascii.no_junk.fasta

    # Fix all matching FASTA files under a directory:
    fix_fasta_encoding.py /data/seqs/*.fasta

    # Dry-run (show what would be done without modifying files):
    fix_fasta_encoding.py --dry-run spikenuc1207.native2ascii.no_junk.fasta

    # Force re-processing even when .orig backup already exists:
    fix_fasta_encoding.py --overwrite spikenuc1207.native2ascii.no_junk.fasta
"""

import argparse
import os
import re
import sys

VERSION = "202603312000"

_parser = argparse.ArgumentParser(
    description=__doc__,
    formatter_class=argparse.RawDescriptionHelpFormatter,
)
_parser.add_argument(
    "infiles", nargs="+", metavar="FASTA",
    help="One or more FASTA files to normalize in-place.",
)
_parser.add_argument(
    "--dry-run", action="store_true",
    help=(
        "Show what would be done without modifying any files.  "
        "Prints a summary of changed lines per file."
    ),
)
_parser.add_argument(
    "--overwrite", action="store_true",
    help=(
        "Re-process files even when a .orig backup already exists.  "
        "The existing .orig is overwritten."
    ),
)
_parser.add_argument(
    "--stats-only", action="store_true",
    help=(
        "Like --dry-run but suppresses per-line output; only prints the "
        "per-file summary line."
    ),
)
_parser.add_argument("--version", action="version", version=f"%(prog)s {VERSION}")


# ── helpers ──────────────────────────────────────────────────────────────────

def _unescape_unicode(s: str) -> str:
    r"""Convert literal \uXXXX escape sequences to real Unicode characters.

    Fast path: if '\\u' is not in the string the function returns immediately
    with negligible overhead.
    """
    if "\\u" not in s:
        return s
    return re.sub(
        r"\\u([0-9a-fA-F]{4})",
        lambda m: chr(int(m.group(1), 16)),
        s,
    )


def _decode_line(raw: bytes) -> str:
    """Decode one raw FASTA byte line to a clean Unicode str.

    Steps (same as summarize_fasta_pipeline.py / create_list_of_discarded_sequences.py):
      1. Try UTF-8; on failure fall back to Latin-1.
      2. Convert literal \\uXXXX escape sequences to real Unicode codepoints.
    """
    try:
        line = raw.decode("utf-8")
    except UnicodeDecodeError:
        line = raw.decode("latin-1")
    return _unescape_unicode(line)


def _process_file(path: str, dry_run: bool, overwrite: bool,
                  stats_only: bool) -> int:
    """Normalize *path* in-place.  Returns the number of lines changed."""
    orig_path = path + ".orig"

    if os.path.exists(orig_path) and not overwrite:
        print(
            f"  Skipping {path}: backup {os.path.basename(orig_path)} already "
            "exists (use --overwrite to force).",
            file=sys.stderr,
        )
        return 0

    # Read all lines, decoding as we go.
    changed = 0
    out_lines: list[str] = []
    with open(path, "rb") as fh:
        for raw in fh:
            # Preserve the original line ending (strip nothing here).
            original_text = raw.decode("latin-1")   # lossless reference decode
            clean = _decode_line(raw)
            out_lines.append(clean if clean.endswith("\n") else clean.rstrip("\r\n") + "\n")
            if clean.rstrip("\r\n") != original_text.rstrip("\r\n"):
                changed += 1
                if not stats_only and not dry_run:
                    pass  # changes written silently; flag them in dry-run below
                elif not stats_only:
                    # dry-run: show diff
                    print(
                        f"  ~ {original_text.rstrip()!r}\n"
                        f"  + {clean.rstrip()!r}",
                    )

    if changed == 0:
        print(f"  {path}: no changes needed.")
        return 0

    action = "Would write" if dry_run else "Writing"
    print(f"  {action} {path}: {changed:,} line(s) changed.")

    if dry_run:
        return changed

    # Write clean UTF-8 to a temp path, then rename into place.
    tmp_path = path + ".encoding_fix_tmp"
    with open(tmp_path, "w", encoding="utf-8") as fh:
        fh.writelines(out_lines)

    # Backup original (or overwrite existing backup if --overwrite).
    os.rename(path, orig_path)
    print(f"  Renamed {os.path.basename(path)} -> {os.path.basename(orig_path)}")
    os.rename(tmp_path, path)
    print(f"  Wrote clean UTF-8 -> {os.path.basename(path)}")
    return changed


# ── entry point ──────────────────────────────────────────────────────────────

def main() -> None:
    opts = _parser.parse_args()

    total_files = 0
    total_changed = 0
    for path in opts.infiles:
        if not os.path.isfile(path):
            print(f"Warning: not a file, skipping: {path}", file=sys.stderr)
            continue
        print(f"{path}:", file=sys.stderr)
        n = _process_file(
            path,
            dry_run=opts.dry_run,
            overwrite=opts.overwrite,
            stats_only=opts.stats_only,
        )
        total_files += 1
        total_changed += n

    print(
        f"\nDone: {total_files} file(s) processed, "
        f"{total_changed:,} total line(s) changed.",
        file=sys.stderr,
    )


if __name__ == "__main__":
    main()
