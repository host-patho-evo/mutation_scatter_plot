#!/usr/bin/env python3
"""Summarize FASTA pipeline: record counts and NNNNx count sums across pipeline stages.

For each file found under <search_path> matching <filename_prefix>*.fasta{,.old,.ori,.orig},
computes:
  - Number of FASTA records  (grep -c '^>')
  - Sum of NNNNx count prefixes in the FASTA IDs
      (non-NNNNx IDs, e.g. GISAID accessions, contribute 0 to the sum)

Parent-child relationships between pipeline files are inferred automatically
from filename suffixes: file B is a child of file A if strip_fasta_suffix(B)
starts with strip_fasta_suffix(A) + ".".  For each parent->child pair,
pre-computed ancillary files are used when they exist and are up-to-date
(newer than both FASTA files), otherwise create_list_of_discarded_sequences.py
is invoked in --inverted mode.

Per-step output files written alongside each child FASTA
---------------------------------------------------------
{child_stem}.discarded_sha256_hashes.txt
    One NNNNx.sha256hex entry per unique discarded sequence.  The sum of all
    NNNNx prefixes equals the total number of individual sequences discarded
    at this pipeline step.  Used as the tier-1 speed cache by this script.

{child_stem}.discarded_original_ids.txt
    One original FASTA header per individual discarded sequence (count-expanded).
    'wc -l' of this file should equal the sum of all NNNNx prefixes in the
    .discarded_sha256_hashes.txt.  Any discrepancy is reported as a Warning
    by create_list_of_discarded_sequences.py.  Generated automatically whenever
    a mapping TSV or original FASTA is available.

Auto-migration of legacy files
-------------------------------
Older runs wrote NNNNx.sha256hex content into *.discarded_original_ids.txt
(before the sha256-hashes / original-IDs naming split).  On first run,
summarize_fasta_pipeline.py detects such files by peeking at the first line
and renames them automatically to *.discarded_sha256_hashes.txt.  Genuine
*.discarded_original_ids.txt files (containing plain FASTA headers) are never
renamed.

Cache priority (fastest first):
  1. <child_stem>.discarded_sha256_hashes.txt  newer than parent + child FASTA
  2. <parent_stem>.sha256_to_ids.tsv           newer than parent FASTA
  3. Full FASTA scan via --original-infilename (slow fallback)

Usage:
    summarize_fasta_pipeline.py <search_path> <filename_prefix> [options]

Options:
    --no-discard-stats   Skip the per-step discard-statistics calls entirely.
    --verbose            Print progress messages during TSV generation and
                         discard-stats computation (the ↳, [auto-generating
                         TSV], and Info: lines).  By default only the table
                         itself is shown.
    --disable-discarded-original-ids-file
                         Do not write per-step .discarded_sha256_hashes.txt
                         files (or the companion .discarded_original_ids.txt).
                         By default both are always written alongside each
                         pipeline step; .discarded_sha256_hashes.txt serves as
                         a fast-read cache of NNNNx.sha256hex entries for
                         subsequent runs, and .discarded_original_ids.txt
                         provides the expanded original FASTA IDs for review.
    --add-missing-checksums-to-fasta-files
                         If a pipeline FASTA file has legacy NNNNx IDs (no
                         sha256 hex in the ID), rewrite it in-place with
                         NNNNx.sha256hex IDs and rename the original to
                         .fasta.orig.  This makes all subsequent analyses
                         faster: sha256 can be read from the ID directly
                         instead of being recomputed from the sequence.
                         Skipped if a .fasta.orig backup already exists.
                         Only applies to deduplicated files whose IDs have
                         an NNNNx count prefix.
    --full-fasta-header
                         When building *.sha256_to_ids.tsv files, also write
                         a companion *.sha256_to_descr_lines.tsv that stores
                         the full FASTA header line (everything after '>'),
                         not just the first whitespace-delimited word.  This
                         second file is built in the same single FASTA scan
                         so no extra I/O is needed.  When present,
                         *.sha256_to_descr_lines.tsv is used as the mapping
                         source for *.discarded_original_ids.txt, giving
                         the full description in those files.

Examples:
    # Show full pipeline table with per-step discard statistics
    summarize_fasta_pipeline.py . spikenuc1207.native2ascii.no_junk

    # Fast table only (skip discard stats)
    summarize_fasta_pipeline.py .. spikenuc1207.native2ascii.no_junk --no-discard-stats

    # Upgrade legacy NNNNx IDs to NNNNx.sha256hex in all matching FASTA files
    summarize_fasta_pipeline.py . spikenuc1207.native2ascii.no_junk \
        --add-missing-checksums-to-fasta-files

    # Per step the following files are written (example for one step):
    #   spikenuc1207...counts.clean.exactly_3822.discarded_sha256_hashes.txt
    #   spikenuc1207...counts.clean.exactly_3822.discarded_original_ids.txt
    #
    # Verify no discrepancy for that step:
    #   wc -l *exactly_3822.discarded_original_ids.txt
    #   awk '{print $1}' *exactly_3822.discarded_sha256_hashes.txt \
    #       | sed 's/x.*//' | awk '{s+=$1} END {print s}'
"""

import datetime
import glob
import hashlib
import os
import subprocess
import sys

FASTA_SUFFIXES = ('.fasta.orig', '.fasta.ori', '.fasta.old', '.fasta')

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
DISCARD_SCRIPT = os.path.join(SCRIPT_DIR, 'create_list_of_discarded_sequences.py')

# ── helpers ───────────────────────────────────────────────────────────────────

def _strip_fasta_suffix(path: str) -> str:
    """Return path with the trailing FASTA suffix removed (longest match first)."""
    for sfx in FASTA_SUFFIXES:
        if path.endswith(sfx):
            return path[:-len(sfx)]
    return path


def _mtime(path: str) -> float:
    """Return mtime of path, or 0.0 if it does not exist."""
    try:
        return os.path.getmtime(path)
    except OSError:
        return 0.0


def _count_records(path: str) -> int:
    """Number of '>' header lines in a FASTA file."""
    result = subprocess.run(['grep', '-c', '^>', path], capture_output=True, text=True, check=False)
    text = result.stdout.strip()
    return int(text) if text else 0


def _sum_nnnx_counts(path: str) -> int:
    """Sum of the leading NNNNx integer prefixes across all FASTA header IDs."""
    cmd = (
        "grep '^>' " + _shell_quote(path) + r" | cut -c 2-"
        r" | awk '{print $1}'"
        r" | sed -e 's/x.*//'"
        r" | awk '{SUM += $1} END {print SUM+0}'"
    )
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True, check=False)
    text = result.stdout.strip()
    return int(text) if text else 0


def _shell_quote(s: str) -> str:
    return "'" + s.replace("'", "'\\''") + "'"


def _delta_str(current: int, previous: int) -> str:
    diff = current - previous
    sign = '+' if diff >= 0 else ''
    return f"{sign}{diff:,}"


def _pct_str(current: int, reference: int) -> str:
    if reference == 0:
        return "n/a"
    return f"{current / reference * 100:.2f}%"


# ── ancillary-file helpers ─────────────────────────────────────────────────────

def _is_sha256_hashes_content(path: str) -> bool:
    """Return True if the first non-empty line of *path* looks like a
    NNNNx.sha256hex entry (i.e. the file contains sha256-hash lines, not
    expanded original FASTA headers).

    Detection heuristic: the first word of the first non-empty line must match
    ``\\d+x.[0-9a-f]{64}`` — the NNNNx.sha256hex format produced by
    create_list_of_discarded_sequences.py.  Original FASTA headers (GISAID
    accessions, plain descriptions) never match this pattern.
    """
    try:
        with open(path, 'r', encoding='utf-8', errors='replace') as fh:
            for raw_line in fh:
                line = raw_line.strip()
                if not line:
                    continue
                first_word = line.split()[0]
                xdot = first_word.find('x.')
                if xdot > 0 and first_word[:xdot].isdigit():
                    suffix = first_word[xdot + 2:]
                    if len(suffix) == 64 and all(c in '0123456789abcdefABCDEF' for c in suffix):
                        return True
                return False  # first non-empty line does not match
    except OSError:
        pass
    return False


def _migrate_legacy_discarded_txt(stem: str) -> None:
    """Rename a legacy *.discarded_original_ids.txt to *.discarded_sha256_hashes.txt
    when it contains NNNNx.sha256hex content (the old naming convention).

    The rename is skipped when:
      - the new *.discarded_sha256_hashes.txt already exists (nothing to do), or
      - the file content looks like expanded original FASTA headers (real
        original-IDs file — must not be renamed).
    """
    old_path = stem + '.discarded_original_ids.txt'
    new_path = stem + '.discarded_sha256_hashes.txt'
    if os.path.exists(new_path):
        return  # already migrated
    if not os.path.exists(old_path):
        return  # nothing to migrate
    if not _is_sha256_hashes_content(old_path):
        return  # genuine original-IDs file — leave it alone
    os.rename(old_path, new_path)
    print(
        f"  [migrate] renamed {os.path.basename(old_path)}"
        f" -> {os.path.basename(new_path)}"
        f" (detected sha256-hash content in legacy file)",
        flush=True,
    )


def _fresh_discarded_txt(child_path: str, parent_path: str) -> str | None:
    """Return path to a valid .discarded_sha256_hashes.txt if it exists and is
    newer than BOTH the child and parent FASTA files, else None.

    Automatically migrates a legacy .discarded_original_ids.txt that contains
    NNNNx.sha256hex content to the new filename before the freshness check.
    """
    stem = _strip_fasta_suffix(child_path)
    _migrate_legacy_discarded_txt(stem)
    candidate = stem + '.discarded_sha256_hashes.txt'
    if not os.path.exists(candidate):
        return None
    txt_mtime = _mtime(candidate)
    if txt_mtime > _mtime(child_path) and txt_mtime > _mtime(parent_path):
        return candidate
    return None  # stale


def _fresh_tsv(parent_path: str) -> str | None:
    """Return path to a valid .sha256_to_ids.tsv if it exists and is newer
    than the parent FASTA file, else None."""
    candidate = _strip_fasta_suffix(parent_path) + '.sha256_to_ids.tsv'
    if not os.path.exists(candidate):
        return None
    if _mtime(candidate) > _mtime(parent_path):
        return candidate
    return None  # stale


def _fresh_descr_tsv(parent_path: str) -> str | None:
    """Return path to a valid .sha256_to_descr_lines.tsv if it exists and is
    newer than the parent FASTA file, else None.

    This file has the same tab-separated format as .sha256_to_ids.tsv but
    stores the full FASTA header line (everything after '>') instead of just
    the first whitespace-delimited word.
    """
    candidate = _strip_fasta_suffix(parent_path) + '.sha256_to_descr_lines.tsv'
    if not os.path.exists(candidate):
        return None
    if _mtime(candidate) > _mtime(parent_path):
        return candidate
    return None  # stale


def _decode_fasta_line(raw: bytes) -> str:
    """Decode a raw FASTA byte line to str (UTF-8 with Latin-1 fallback)."""
    try:
        return raw.decode("utf-8")
    except UnicodeDecodeError:
        return raw.decode("latin-1")


def _extract_sha256_from_id(record_id: str) -> str | None:
    """Return the 64-char hex sha256 embedded in a NNNNx.sha256hex ID, or None."""
    xdot = record_id.find('x.')
    if xdot >= 0:
        candidate = record_id[xdot + 2:].split()[0]
        if len(candidate) == 64 and all(c in '0123456789abcdefABCDEF' for c in candidate):
            return candidate.lower()
    return None


def _build_tsv(fasta_path: str, tsv_path: str,
               descr_tsv_path: str = "",
               verbose: bool = True) -> dict:
    """Single-pass scan of *fasta_path* to build a sha256->IDs mapping TSV.

    For each record the sha256 is computed over the uppercase, dash-stripped
    sequence — identical to the normalisation used by count_same_sequences.py.

    Output TSV columns (tab-separated, no header):
        sha256hex   count   id_1   id_2   ...

    When *descr_tsv_path* is non-empty, a second TSV is written in the same
    scan using the full FASTA header line (everything after '>') instead of
    just the first word.  Both files are produced from a single I/O pass.

    Returns a dict mapping each FASTA ID (first word) to its sha256hex string.
    """
    id_map: dict = {}    # sha256 -> [count, [first-word-id, ...]]
    descr_map: dict = {} # sha256 -> [count, [full-header, ...]]
    name = None
    full_header = None
    seq_parts: list = []

    def _flush(flush_name: str, flush_full: str, flush_parts: list) -> None:
        seq = "".join(flush_parts).rstrip("\r\n").replace("-", "").upper()
        digest = hashlib.sha256(seq.encode()).hexdigest()
        if digest in id_map:
            id_map[digest][0] += 1
            id_map[digest][1].append(flush_name)
            if descr_tsv_path:
                descr_map[digest][0] += 1
                descr_map[digest][1].append(flush_full)
        else:
            id_map[digest] = [1, [flush_name]]
            if descr_tsv_path:
                descr_map[digest] = [1, [flush_full]]

    n_in = 0
    with open(fasta_path, "rb") as fh:
        for raw in fh:
            line = _decode_fasta_line(raw).rstrip("\r\n")
            if not line:
                continue
            if line[0] == ">":
                if name is not None:
                    _flush(name, full_header, seq_parts)
                    n_in += 1
                hdr = line[1:]
                toks = hdr.split()
                name = toks[0] if toks else ""
                # Escape any embedded TAB characters as the two-character
                # literal '\t' so the TSV field delimiter is unambiguous.
                # Readers of *.sha256_to_descr_lines.tsv must unescape
                # '\\t' back to '\t' when consuming description fields.
                # (Actual TABs in FASTA headers are extremely rare, but
                # this encoding is lossless unlike replacing with a space.)
                full_header = hdr.replace("\t", "\\t")
                seq_parts = []
            else:
                seq_parts.append(line)
        if name is not None:
            _flush(name, full_header, seq_parts)
            n_in += 1

    with open(tsv_path, "w", encoding="utf-8") as out:
        for digest, (count, ids) in id_map.items():
            out.write("\t".join([digest, str(count)] + ids) + "\n")

    if descr_tsv_path:
        with open(descr_tsv_path, "w", encoding="utf-8") as out:
            for digest, (count, descrs) in descr_map.items():
                out.write("\t".join([digest, str(count)] + descrs) + "\n")

    if verbose:
        extra = f" + {os.path.basename(descr_tsv_path)}" if descr_tsv_path else ""
        print(f"    Info: built TSV with {len(id_map):,} unique sequences"
              f" from {n_in:,} records -> {tsv_path}{extra}", flush=True)

    # Return inverted mapping: id -> sha256 (one entry per original record).
    return {orig_id: sha for sha, (_, ids) in id_map.items() for orig_id in ids}


def _has_legacy_ids(id_to_sha: dict) -> bool:
    """Return True if ANY ID in the file lacks an embedded sha256.

    An ID is 'modern' if it matches NNNNx.sha256hex (64-char hex after 'x.').
    Only files composed entirely of NNNNx-prefix IDs (deduplicated pipeline
    output) are candidates for enrichment.
    """
    return any(_extract_sha256_from_id(rec_id) is None for rec_id in id_to_sha)


def _enrich_fasta(fasta_path: str, id_to_sha: dict) -> str | None:
    """Rewrite *fasta_path* with sha256hex appended to each NNNNx ID.

    The original file is renamed to ``{fasta_path}.orig`` before the new
    version is written.  If a ``.orig`` backup already exists the enrichment
    is skipped to avoid overwriting a previous backup.

    Only IDs that lack an embedded sha256 are modified; any ID that already
    carries sha256 is left unchanged.

    Returns the path of the renamed original on success, or None if skipped.
    """
    orig_path = fasta_path + '.orig'
    if os.path.exists(orig_path):
        print(f"  [add-missing-checksums] skipping {os.path.basename(fasta_path)}"
              f" — backup {os.path.basename(orig_path)} already exists", flush=True)
        return None
    os.rename(fasta_path, orig_path)
    print(f"  [add-missing-checksums] renamed {os.path.basename(fasta_path)}"
          f" -> {os.path.basename(orig_path)}", flush=True)

    enriched = 0
    with (open(orig_path, "rb") as src,
          open(fasta_path, "w", encoding="utf-8") as dst):
        for raw in src:
            line = _decode_fasta_line(raw).rstrip("\r\n")
            if not line:
                dst.write("\n")
                continue
            if line[0] == ">":
                toks = line[1:].split()
                rec_id = toks[0] if toks else ""
                if _extract_sha256_from_id(rec_id) is None and rec_id in id_to_sha:
                    # Keep the description (everything after the ID) intact.
                    rest = line[1 + len(rec_id):]
                    dst.write(f">{rec_id}.{id_to_sha[rec_id]}{rest}\n")
                    enriched += 1
                else:
                    dst.write(line + "\n")
            else:
                dst.write(line + "\n")

    print(f"  [add-missing-checksums] wrote {enriched:,} enriched IDs"
          f" -> {os.path.basename(fasta_path)}", flush=True)
    return orig_path


def _ensure_tsv(parent_path: str, add_checksums: bool = False,
                full_fasta_header: bool = False,
                verbose: bool = True) -> tuple[str | None, str | None]:
    """Return (ids_tsv_path, descr_tsv_path) for *parent_path*.

    ids_tsv_path  : path to *.sha256_to_ids.tsv (first-word IDs), or None.
    descr_tsv_path: path to *.sha256_to_descr_lines.tsv (full headers), or
                    None when *full_fasta_header* is False or generation failed.

    Both TSVs are produced from a single FASTA scan when they need to be
    (re)built.  Any embedded TABs in full headers are replaced with spaces.

    When *add_checksums* is True and the file has legacy NNNNx IDs (no sha256
    embedded), the FASTA is rewritten in-place with NNNNx.sha256hex IDs and
    the original is renamed to .fasta.orig.
    """
    stem = _strip_fasta_suffix(parent_path)
    candidate     = stem + '.sha256_to_ids.tsv'
    descr_candidate = stem + '.sha256_to_descr_lines.tsv' if full_fasta_header else ""

    ids_tsv   = _fresh_tsv(parent_path)
    descr_tsv = _fresh_descr_tsv(parent_path) if full_fasta_header else None

    # If both TSVs are fresh and no checksums requested, return immediately.
    if ids_tsv and (not full_fasta_header or descr_tsv) and not add_checksums:
        return ids_tsv, descr_tsv

    # If add_checksums: check whether the FASTA already has modern IDs.
    if ids_tsv and add_checksums:
        try:
            with open(parent_path, "rb") as fh:
                for raw in fh:
                    line = _decode_fasta_line(raw)
                    if line.startswith(">"):
                        first_id = line[1:].split()[0]
                        if _extract_sha256_from_id(first_id) is not None:
                            # Already modern; only rebuild descr_tsv if missing.
                            if full_fasta_header and not descr_tsv:
                                break  # fall through to build
                            return ids_tsv, descr_tsv
                        break  # legacy IDs — fall through
        except OSError:
            return ids_tsv, descr_tsv

    if verbose:
        print(f"  [auto-generating TSV] scanning {os.path.basename(parent_path)} \u2026",
              flush=True)
    try:
        id_to_sha = _build_tsv(parent_path, candidate,
                               descr_tsv_path=descr_candidate,
                               verbose=verbose)
    except OSError as exc:
        print(f"    [TSV generation failed: {exc}]", flush=True)
        return None, None

    if add_checksums and _has_legacy_ids(id_to_sha):
        _enrich_fasta(parent_path, id_to_sha)
        if verbose:
            print("  [auto-generating TSV] rebuilding TSV from enriched FASTA \u2026",
                  flush=True)
        try:
            _build_tsv(parent_path, candidate,
                       descr_tsv_path=descr_candidate,
                       verbose=verbose)
        except OSError as exc:
            print(f"    [TSV rebuild failed: {exc}]", flush=True)

    return (candidate,
            descr_candidate if full_fasta_header else None)


def _read_discarded_txt_stats(txt_path: str) -> tuple[int, int]:
    """Return (n_ids, nnnx_sum) from a .discarded_sha256_hashes.txt.

    Each line is a NNNNx.sha256hex entry (one per unique discarded sequence).
    - n_ids   = number of lines (unique sequences)
    - nnnx_sum = sum of the NNNNx count prefixes (original-sequence total)
    """
    n_ids = int(subprocess.run(
        f"wc -l < {_shell_quote(txt_path)}",
        shell=True, capture_output=True, text=True, check=False,
    ).stdout.strip() or 0)
    nnnx_sum = int(subprocess.run(
        f"awk '{{print $1}}' {_shell_quote(txt_path)}"
        r" | sed -e 's/x.*//'"
        r" | awk '{SUM += $1} END {print SUM+0}'",
        shell=True, capture_output=True, text=True, check=False,
    ).stdout.strip() or 0)
    return n_ids, nnnx_sum


# ── per-pair discard stats ────────────────────────────────────────────────────

def _compute_discard_stats(parent_path: str, child_path: str,
                           add_checksums: bool = False,
                           full_fasta_header: bool = False,
                           save_discard_list: bool = True,
                           verbose: bool = True
                           ) -> tuple[int, int] | tuple[None, None]:
    """Compute discarded-ID stats for a parent->child pipeline pair.

    Cache priority (fastest first — checked by timestamp):
      1. Existing .discarded_sha256_hashes.txt newer than both FASTAs -> read directly.
         Exception: if the companion .discarded_original_ids.txt is missing and
         save_discard_list is True, falls through to tiers 2-4 so the companion is
         generated without touching any FASTA files.
      2. Existing .sha256_to_ids.tsv newer than parent FASTA         -> TSV scan.
      3. Auto-generate .sha256_to_ids.tsv (in-process FASTA scan)   -> then TSV scan.
      4. Fallback: full FASTA scan via --original-infilename (slow).

    Returns (n_discard_ids, n_discard_nnnx_sum), or (None, None) on failure.
    When *save_discard_list* is True (default), both output files are written:
      {child_stem}.discarded_sha256_hashes.txt  (NNNNx.sha256hex entries; fast cache)
      {child_stem}.discarded_original_ids.txt   (expanded original FASTA headers)
    """
    parent_base = os.path.basename(parent_path)
    child_base  = os.path.basename(child_path)
    child_stem  = _strip_fasta_suffix(child_path)

    # ── tier 1: fresh .discarded_sha256_hashes.txt ───────────────────────────
    txt = _fresh_discarded_txt(child_path, parent_path)
    if txt:
        original_ids_path = child_stem + '.discarded_original_ids.txt'
        companion_missing = save_discard_list and not os.path.exists(original_ids_path)
        if not companion_missing:
            if verbose:
                print(f"  \u21b3 [cached TXT] {os.path.basename(txt)}", flush=True)
            return _read_discarded_txt_stats(txt)
        # Companion .discarded_original_ids.txt is missing: fall through to
        # regenerate both files.  No FASTA timestamp change needed.
        if verbose:
            print(
                f"  \u21b3 [cached TXT] {os.path.basename(txt)}"
                f" but companion {os.path.basename(original_ids_path)} missing"
                " \u2014 regenerating both",
                flush=True,
            )

    # ── tiers 2–4: invoke create_list_of_discarded_sequences.py ────────
    if not os.path.exists(DISCARD_SCRIPT):
        print(f"  [discard stats] script not found: {DISCARD_SCRIPT}", flush=True)
        return None, None

    tsv, descr_tsv = _ensure_tsv(parent_path, add_checksums=add_checksums,
                                  full_fasta_header=full_fasta_header,
                                  verbose=verbose)
    # Prefer the descr TSV (full headers) when the caller requested it
    # and it was successfully built; fall back to ID-only TSV otherwise.
    mapping_tsv = descr_tsv if (full_fasta_header and descr_tsv) else tsv
    if mapping_tsv:
        source_arg   = f'--mapping-outfile={mapping_tsv}'
        source_label = f"TSV: {os.path.basename(mapping_tsv)}"
    else:
        source_arg   = f'--original-infilename={parent_path}'
        source_label = f"FASTA scan: {parent_base}"

    if save_discard_list:
        # Write the sha256-hashes file (NNNNx.sha256hex entries); this is what
        # _read_discarded_txt_stats() reads for the pipeline summary table.
        # create_list_of_discarded_sequences.py also auto-generates the companion
        # .discarded_original_ids.txt (expanded original FASTA headers) alongside.
        outfile = child_stem + '.discarded_sha256_hashes.txt'
        outfile_args = [f'--outfile={outfile}', '--overwrite']
    else:
        outfile      = None
        outfile_args = ['--outfile=/dev/null']

    cmd = [
        sys.executable, DISCARD_SCRIPT,
        f'--infilename={child_path}',
        source_arg,
        '--inverted',
        '--output-context=discarded',
        *(['--full-fasta-header'] if full_fasta_header else []),
        *outfile_args,
    ]
    if verbose:
        print(f"  \u21b3 [{source_label}] -> {child_base}", flush=True)
    result = subprocess.run(cmd, capture_output=True, text=True, check=False)
    # Suppress Info: lines — the key numbers will appear in the table instead.
    for line in result.stderr.splitlines():
        if line.startswith(('Warning:', 'Error:')):
            print(f"    {line}", flush=True)
    if result.returncode != 0:
        print(f"    [exit code {result.returncode}]", flush=True)
        return None, None

    if outfile and os.path.exists(outfile):
        return _read_discarded_txt_stats(outfile)
    # Fallback: check if the sha256-hashes file was written with a non-default name.
    sha_file = child_stem + '.discarded_sha256_hashes.txt'
    if os.path.exists(sha_file):
        return _read_discarded_txt_stats(sha_file)
    return None, None


# ── main ──────────────────────────────────────────────────────────────────────

def main() -> None:
    """Entry point: parse args, scan files, print table, optionally show discard stats."""
    args = sys.argv[1:]
    if len(args) < 2 or '--help' in args or '-h' in args:
        print(__doc__, file=sys.stderr)
        sys.exit(0 if '--help' in args or '-h' in args else 1)

    search_path = args[0]
    prefix      = args[1]
    do_discard         = '--no-discard-stats'                      not in args
    verbose            = '--verbose'                               in args
    save_discard_list  = '--disable-discarded-original-ids-file'   not in args
    add_checksums      = '--add-missing-checksums-to-fasta-files'  in args
    full_fasta_header  = '--full-fasta-header'                     in args

    # ── collect matching files ────────────────────────────────────────────────
    found: set[str] = set()
    for suffix in FASTA_SUFFIXES:
        for pat in (
            os.path.join(search_path, prefix + '*' + suffix),
            os.path.join(search_path, '**', prefix + '*' + suffix),
        ):
            found.update(glob.glob(pat, recursive=True))

    if not found:
        print(f"No files found under '{search_path}' matching '{prefix}*.fasta{{,.old,.ori,.orig}}'",
              file=sys.stderr)
        sys.exit(1)

    # Sort by stem length (shorter = earlier in pipeline), then alphabetically.
    files = sorted(found, key=lambda p: (len(_strip_fasta_suffix(os.path.basename(p))),
                                          os.path.basename(p)))
    print(f"Found {len(files):,} file(s).\n", file=sys.stderr)

    # ── infer parent->child pairs from naming convention ──────────────────────
    # base_B.startswith(base_A + ".") defines the relationship.
    bases = [_strip_fasta_suffix(os.path.basename(f)) for f in files]

    def _direct_parent_idx(child_idx: int) -> int | None:
        child_base = bases[child_idx]
        best = None
        for i, b in enumerate(bases):
            if i != child_idx and child_base.startswith(b + '.'):
                if best is None or len(b) > len(bases[best]):
                    best = i
        return best

    parent_map: dict[int, int] = {}
    for i in range(len(files)):
        p = _direct_parent_idx(i)
        if p is not None:
            parent_map[i] = p

    # ── gather per-file data ─────────────────────────────────────────────────
    rows: list[tuple[str, str, int, int]] = []
    for f in files:
        display = os.path.relpath(f, search_path)
        print(f"  scanning {display} \u2026", file=sys.stderr, flush=True)
        mtime_s = datetime.datetime.fromtimestamp(os.path.getmtime(f)).strftime('%Y-%m-%d %H:%M')
        rows.append((display, mtime_s, _count_records(f), _sum_nnnx_counts(f)))

    # ── phase 2: compute discard stats for all pairs ─────────────────────────
    # Runs before table printing so the numbers can appear as proper columns.
    discard_data: dict[int, tuple[int, int]] = {}  # child_idx -> (n_ids, nnnx_sum)
    if do_discard:
        print(file=sys.stderr)
        for i, f_child in enumerate(files):
            p = parent_map.get(i)
            if p is not None:
                n_d, s_d = _compute_discard_stats(
                    files[p], f_child,
                    add_checksums=add_checksums,
                    full_fasta_header=full_fasta_header,
                    save_discard_list=save_discard_list,
                    verbose=verbose,
                )
                if n_d is not None:
                    discard_data[i] = (n_d, s_d)

    # ── print table ──────────────────────────────────────────────────────────
    col_file = max(max(len(r[0]) for r in rows), len("File"))
    sep, w_num, w_delta, w_ts = "  ", 14, 16, 16
    w_disc1 = len("'Discarded original FASTA IDs'")   # 30
    w_disc2 = len("'Sum of discarded sequences'")       # 28

    _hdr_disc1  = "'Discarded original FASTA IDs'"
    _hdr_disc2  = "'Sum of discarded sequences'"
    _hdr_nnnx   = "'Sum of NNNNx'"
    _hdr_drec   = '\u0394Records'
    _hdr_dsum   = '\u0394SumToParent'
    disc_header = (
        f"{sep}{_hdr_disc1:>{w_disc1}}{sep}{_hdr_disc2:>{w_disc2}}"
        if do_discard else ""
    )
    header = (
        f"{'File':<{col_file}}{sep}"
        f"{'Modified':<{w_ts}}{sep}"
        f"{'Records':>{w_num}}{sep}{_hdr_drec:>{w_delta}}{sep}"
        f"{_hdr_nnnx:>{w_num}}{sep}{_hdr_dsum:>{w_delta}}"
        + disc_header
    )
    rule = '-' * len(header)
    print()
    print(header)
    print(rule)

    for i, (display, mtime_s, n_rec, n_sum) in enumerate(rows):
        p = parent_map.get(i)
        if p is not None:
            prev_rec, prev_sum = rows[p][2], rows[p][3]
            d_rec = _delta_str(n_rec, prev_rec)
            d_sum = _delta_str(n_sum, prev_sum)
            parent_label = f"vs {os.path.basename(files[p])}"
        else:
            d_rec = d_sum = '\u2014'
            parent_label = ''

        if do_discard and p is not None:
            _em = '\u2014'  # em dash — pre-assigned to avoid backslash inside f-string {} (Python < 3.12)
            if i in discard_data:
                n_d, s_d = discard_data[i]
                disc_cols = f"{sep}{n_d:>{w_disc1},}{sep}{s_d:>{w_disc2},}"
            else:
                disc_cols = f"{sep}{_em:>{w_disc1}}{sep}{_em:>{w_disc2}}"
        else:
            disc_cols = ''

        print(
            f"{display:<{col_file}}{sep}"
            f"{mtime_s:<{w_ts}}{sep}"
            f"{n_rec:>{w_num},}{sep}{d_rec:>{w_delta}}{sep}"
            f"{n_sum:>{w_num},}{sep}{d_sum:>{w_delta}}"
            + disc_cols
            + (f"  ({parent_label})" if parent_label else '')
        )



    print(rule)

    # ── overall summary ───────────────────────────────────────────────────────
    if len(rows) >= 2:
        first_rec, first_sum = rows[0][2], rows[0][3]
        last_rec,  last_sum  = rows[-1][2], rows[-1][3]
        print()
        print("Overall change  (last vs first):")
        print(f"  Records : {_delta_str(last_rec, first_rec):>16}  ({_pct_str(last_rec, first_rec)} of first)")
        print(f"  Sum     : {_delta_str(last_sum, first_sum):>16}  ({_pct_str(last_sum, first_sum)} of first)")


if __name__ == '__main__':
    main()
