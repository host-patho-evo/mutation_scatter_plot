#!/usr/bin/env python3
# pylint: disable=too-many-lines
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
    --jobs N             Number of parallel workers for Phase 0 (sha256
                         identity check) and Phase 1 (per-file gather).
                         Default is 1 (sequential, unchanged behaviour).
                         On HPC / parallel-filesystem nodes higher values
                         (e.g. --jobs 4) typically halve wall time because
                         file scans are I/O-bound and release the GIL.
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

    --verify-sha256
                         For each NNNNx.sha256hex file, recompute sha256 from
                         the unpadded sequence (dashes stripped, uppercased) and
                         compare it against the value embedded in the ID.
                         Mismatches indicate the sequence was modified after
                         deduplication (e.g. truncated by a downstream tool).
                         Adds four diagnostic columns to the table:
                           Seq clipped(dup)   records whose new sha256 is in the
                                              direct parent's sha256 set (became
                                              a known duplicate after modification)
                           NNNNx clipped(dup) sum of NNNNx prefixes for those
                           Seq clipped(new)   records whose new sha256 is NOT in
                                              the parent (genuinely novel result)
                           NNNNx clipped(new) sum of NNNNx prefixes for those
    --classify-mismatches
                         When --verify-sha256 detects mismatches, perform a
                         second pass over each file's direct parent FASTA to
                         compare the depadded post- and pre-alignment sequences
                         and classify each mismatch as:
                           Altered due to end-clipping
                                       right-, left-, or both-ends clipped
                                       (internal content intact — expected
                                       alignment trimming)
                           Altered inside the sequence
                                       current is not a substring of original
                                       (internal bases were changed, not just
                                       end-trimmed — investigate further)
                         Adds two columns to the table.  Gated behind this flag
                         because it requires an extra scan of each parent FASTA.

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
import re
import subprocess
import sys
import threading
from concurrent.futures import ThreadPoolExecutor, as_completed

FASTA_SUFFIXES = ('.fasta.orig', '.fasta.ori', '.fasta.old', '.fasta')

# Protein FASTA suffixes — these files have NNNNx.dna_sha256 IDs but their
# sequences are amino acid, so sha256(sequence) ≠ sha256(DNA).  They are
# included in the pipeline table with an extra 'Prot unique' column.
PROT_SUFFIXES = ('.prot.fasta', '.prot', '.faa')

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
DISCARD_SCRIPT = os.path.join(SCRIPT_DIR, 'create_list_of_discarded_sequences.py')

VERSION = "202604030900"

# ── helpers ───────────────────────────────────────────────────────────────────

def _ts() -> str:
    """Return a wall-clock timestamp prefix for stderr progress lines.

    Example output: ``'[09:34:51] '``
    """
    return datetime.datetime.now().strftime('[%H:%M:%S] ')


def _get_git_version() -> str:
    """Return a live git version string (``git describe --always --dirty --tags``).

    Mirrors the three-tier logic in ``mutation_scatter_plot/__init__.py``:

    * **Tier 1** — ``git describe`` in the script's own directory (works in any
      live checkout, including editable installs).
    * **Tier 2** — ``GIT_COMMIT`` environment variable (CI/CD fallback).
    * **Tier 3** — ``"unknown"`` if neither is available.

    Falls back gracefully so the script always starts even without git.
    """
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
    if env_ver:
        return env_ver[:12]
    return "unknown"


_GIT_VERSION: str = _get_git_version()

# Protects all stderr writes in parallel Phase 0 / Phase 1 worker threads.
_stderr_lock = threading.Lock()


def _fmt_size(n: int) -> str:
    """Format *n* bytes as a human-readable string (e.g. '31.8 GB')."""
    for unit in ('B', 'KB', 'MB', 'GB', 'TB'):
        if n < 1024:
            return f"{n:.1f} {unit}"
        n //= 1024
    return f"{n:.1f} PB"


def _strip_fasta_suffix(path: str) -> str:
    """Return path with the trailing FASTA suffix removed (longest match first)."""
    for sfx in FASTA_SUFFIXES:
        if path.endswith(sfx):
            return path[:-len(sfx)]
    return path


def _is_prot_file(path: str) -> bool:
    """Return True when *path* is a protein-sequence FASTA (.prot.fasta/.prot/.faa)."""
    return any(path.endswith(sfx) for sfx in PROT_SUFFIXES)


def _fasta_sha256(path: str) -> str:
    """Return sha256 hex digest of the raw byte content of *path*.

    Reads the file in 4 MB chunks so arbitrarily large files (tens of GB)
    can be checked without loading them into RAM.  Used for content-identity
    checks before the expensive per-file gather phase.
    """
    h = hashlib.sha256()
    with open(path, 'rb') as fh:
        while chunk := fh.read(4 << 20):
            h.update(chunk)
    return h.hexdigest()


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
    """Sum of the leading NNNNx integer prefixes across all FASTA header IDs.

    Peeks at the first FASTA header only.  If it has no NNNNx prefix (e.g.
    a GISAID accession like 'Spike|hCoV-19/…'), the sum must be 0 for every
    record — no further file reading is needed.
    """
    # ── fast path: check first header ────────────────────────────────────────
    try:
        with open(path, "rb") as fh:
            for raw in fh:
                if raw[:1] != b'>':
                    continue
                first_id = _decode_fasta_line(raw)[1:].split()[0] if raw else ""
                xpos = first_id.find('x.')
                if xpos > 0:
                    try:
                        int(first_id[:xpos])  # valid NNNNx prefix
                    except ValueError:
                        return 0  # not a number before x.
                    break  # NNNNx file — fall through to full scan
                return 0  # no 'x.' at all — plain GISAID/legacy ID
    except OSError:
        return 0
    # ── full scan: file has NNNNx IDs ────────────────────────────────────────
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


def _unescape_unicode(s: str) -> str:
    r"""Convert literal \uXXXX escape sequences to real Unicode characters.

    Some GISAID exports encode accented characters as the literal six-character
    sequence \u00e9 rather than the proper UTF-8 glyph.  Applied after
    byte-decoding so both forms end up as the same Unicode string.
    """
    if '\\u' not in s:
        return s
    return re.sub(r'\\u([0-9a-fA-F]{4})',
                  lambda m: chr(int(m.group(1), 16)), s)


def _decode_fasta_line(raw: bytes) -> str:
    r"""Decode one raw FASTA byte line to a clean Unicode str.

    Handles GISAID FASTA files that mix UTF-8 and Latin-1 *within the same
    line*.  Uses surrogateescape so each non-UTF-8 byte is individually
    mapped to its Latin-1 equivalent rather than triggering a whole-line
    fallback that would corrupt valid UTF-8 multi-byte sequences.
    Also converts literal ``\uXXXX`` escape sequences to real codepoints.
    """
    s = raw.decode("utf-8", errors="surrogateescape")
    if any(0xDC80 <= ord(ch) <= 0xDCFF for ch in s):
        s = "".join(
            chr(ord(ch) - 0xDC00) if 0xDC80 <= ord(ch) <= 0xDCFF else ch
            for ch in s
        )
    return _unescape_unicode(s)


def _extract_sha256_from_id(record_id: str) -> str | None:
    """Return the 64-char hex sha256 embedded in a NNNNx.sha256hex ID, or None."""
    xdot = record_id.find('x.')
    if xdot >= 0:
        candidate = record_id[xdot + 2:].split()[0]
        if len(candidate) == 64 and all(c in '0123456789abcdefABCDEF' for c in candidate):
            return candidate.lower()
    return None


def _collect_sha256_set(fasta_path: str) -> tuple[set[str], int]:
    """Collect sha256 values associated with *fasta_path*.

    Priority:
      1. Read the first column of the cached *.sha256_to_ids.tsv — up to 60×
         faster than reading the full FASTA because the TSV contains only
         sha256 hex strings and IDs, not gigabytes of sequence data.
      2. Fall back to a header-only FASTA scan when no valid TSV cache exists.

    The sha256 set is used to build the universe of all known sequence
    fingerprints so that --verify-sha256 can classify mismatches as
    →existing (sequence changed to a known sha256) or →novel (new sha256).

    Returns:
        sha256_set : set of lowercase 64-char sha256hex strings
        n_legacy   : number of records with no embedded sha256 (GISAID IDs)
                     — meaningful only from the FASTA fallback path; from
                     the TSV path every first-column entry IS a sha256 so
                     n_legacy is always 0 there.
    """
    # ── TSV fast path ─────────────────────────────────────────────────────
    tsv_candidate = _strip_fasta_suffix(fasta_path) + '.sha256_to_ids.tsv'
    if (os.path.exists(tsv_candidate)
            and os.path.getmtime(tsv_candidate) >= os.path.getmtime(fasta_path)):
        sha256_set: set[str] = set()
        try:
            with open(tsv_candidate, "rb") as fh:
                for raw in fh:
                    # TSV first column is sha256hex; split on first tab only.
                    sha_bytes = raw.split(b'\t', 1)[0].strip()
                    if len(sha_bytes) == 64:
                        sha256_set.add(sha_bytes.decode('ascii').lower())
        except OSError:
            pass
        return sha256_set, 0  # every TSV row is a sha256 — no legacy count

    # ── FASTA fallback: header-only scan ───────────────────────────────────
    sha256_set = set()
    n_legacy = 0
    try:
        with open(fasta_path, "rb") as fh:
            for raw in fh:
                if raw[:1] != b'>':
                    continue
                line = _decode_fasta_line(raw).rstrip("\r\n")
                toks = line[1:].split()
                sha = _extract_sha256_from_id(toks[0]) if toks else None
                if n_legacy == 0 and sha is None:
                    # First header has no embedded sha256 → GISAID/legacy file.
                    # All records use the same ID format, so no sha256 exists
                    # anywhere in this file.  Return immediately using wc
                    # just for the n_legacy count.
                    return set(), _count_records(fasta_path)
                if sha is not None:
                    sha256_set.add(sha)
                else:
                    n_legacy += 1
    except OSError:
        pass
    return sha256_set, n_legacy


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
        # Prefer the sha256 already embedded in the ID (NNNNx.sha256hex format).
        # Recomputing from the current sequence is WRONG for FASTAs whose sequences
        # were transformed after deduplication (e.g. blastn alignment, reverse-
        # complementing, indel correction) while preserving the original NNNNx.sha256
        # ID from count_same_sequences.py.  Recomputing creates a different sha256
        # that no longer matches the ID, causing false duplicates in the TSV and
        # corrupting all downstream discard statistics.
        # Fall back to sequence-based computation only for legacy NNNNx IDs (no sha256)
        # or plain GISAID accession IDs.
        digest = _extract_sha256_from_id(flush_name)
        if digest is None:
            # Sequences in counts.fasta and downstream files may carry
            # alignment-padding dashes, but sha256 is always of the unpadded
            # biological sequence (dashes stripped) to be consistent with what
            # count_same_sequences.py embeds in the NNNNx.sha256hex ID.
            seq = "".join(flush_parts).replace("\r", "").replace("\n", "").replace("-", "").upper()
            if not seq:
                # Empty sequence — sha256("") would be wrong and misleading.
                # Warn and skip rather than polluting the TSV.
                print(f"  Warning: empty sequence for record '{flush_name}' in"
                      f" {fasta_path} — skipped in TSV",
                      file=sys.stderr, flush=True)
                return
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


def _count_prot_unique(path: str) -> tuple[int, int] | None:
    """Group protein sequences by sha256; count unique groups and sum NNNNx per group.

    For each record in a protein FASTA (NNNNx.dna_sha256 IDs), computes
    sha256(depadded_uppercase_protein_seq) and accumulates the NNNNx prefix
    into a per-group running total.

    Returns (n_unique_proteins, total_nnnx_sum) where:
      n_unique_proteins – number of distinct protein sequences
      total_nnnx_sum    – sum of all NNNNx counts (same as _sum_nnnx_counts,
                          but cross-checks that translate.py dropped nothing)

    Returns None on I/O error.
    """
    groups: dict[str, int] = {}  # protein_sha256 → cumulative NNNNx
    name: str | None = None
    parts: list[str] = []
    cur_nnnx: int = 1

    def _chk(chk_name: str, chk_parts: list[str]) -> None:
        seq = (''.join(chk_parts)
               .replace('\r', '').replace('\n', '').replace('-', '').upper())
        if not seq:
            return
        psha = hashlib.sha256(seq.encode()).hexdigest()
        groups[psha] = groups.get(psha, 0) + cur_nnnx

    try:
        with open(path, 'rb') as fh:
            for raw in fh:
                line = _decode_fasta_line(raw).rstrip('\r\n')
                if not line:
                    continue
                if line[0] == '>':
                    if name is not None:
                        _chk(name, parts)
                    tok = line[1:].split()[0] if line[1:].split() else ''
                    name = tok
                    # Extract NNNNx count from ID (format: NNNNx.sha256hex or NNNNx)
                    xpos = tok.find('x.')
                    if xpos < 0:
                        xpos = tok.find('x')
                    try:
                        cur_nnnx = int(tok[:xpos]) if xpos > 0 else 1
                    except ValueError:
                        cur_nnnx = 1
                    parts = []
                else:
                    parts.append(line)
            if name is not None:
                _chk(name, parts)
    except OSError:
        return None
    return len(groups), sum(groups.values())


def _verify_sha256(
        fasta_path: str,
        known_sha256s: set[str] | None = None,
) -> tuple[int, int, int, int] | None:
    """Check ID-embedded sha256 against current sequence content.

    For each NNNNx.sha256hex record the sequence is read, dashes stripped,
    uppercased, and sha256 recomputed.  If the result differs from the sha256
    in the ID the sequence was modified (truncated, corrected, re-aligned …)
    after the deduplication step that created the ID.

    Mismatching records are split into two groups:
      * "clipped→dup": the new computed sha256 IS in *known_sha256s* (the
        sequence was end-clipped/modified to become identical to a sequence
        that already existed in the immediate parent file).
      * "padded→new": the new computed sha256 is NOT in *known_sha256s*
        (modified sequence is novel relative to its direct parent).

    Args:
        fasta_path:    FASTA file to verify.
        known_sha256s: sha256 set of the direct parent file in the pipeline
                       (from _collect_sha256_set).  A mismatch is clipped(dup)
                       when the new sha256 was already present in the parent;
                       padded(new) when it was not.  Pass None to treat every
                       mismatch as padded(new).

    Returns:
        (n_existing, sum_nnnx_existing, n_novel, sum_nnnx_novel) or None when
        the file has no ID-embedded sha256 (GISAID / legacy) or is empty.

    Returns:
        (n_clipped_dup, sum_nnnx_clipped_dup, n_clipped_new, sum_nnnx_clipped_new,
         altered_id_sha256s)
        where *altered_id_sha256s* is the set of id-embedded sha256 strings for every
        record whose current sequence sha256 differs from the one in its ID.
        Returns None when the file has no ID-embedded sha256 or is empty.
    """
    n_existing = 0
    sum_existing = 0
    n_novel = 0
    sum_novel = 0
    altered_id_shas: set[str] = set()   # id-embedded sha256 of every mismatch
    mismatch_seqs:   dict[str, str] = {}  # id-embedded sha256 → depadded current seq
    is_nnnx_file: bool | None = None
    name: str | None = None
    seq_parts: list = []

    def _chk(flush_name: str, flush_parts: list) -> None:
        nonlocal n_existing, sum_existing, n_novel, sum_novel
        id_sha = _extract_sha256_from_id(flush_name)
        if id_sha is None:
            return
        seq = "".join(flush_parts).replace("\r", "").replace("\n", "").replace("-", "").upper()
        if not seq:
            # Empty sequence — would hash to sha256("") and be falsely
            # flagged as →novel.  Warn and skip so it doesn't inflate counts.
            print(f"  Warning: empty sequence for record '{flush_name}' in"
                  f" {fasta_path} — skipped in sha256 verification",
                  file=sys.stderr, flush=True)
            return
        new_sha = hashlib.sha256(seq.encode()).hexdigest()
        if new_sha == id_sha:
            return  # unchanged
        altered_id_shas.add(id_sha)
        mismatch_seqs[id_sha] = seq  # store depadded current seq for classification
        xpos = flush_name.find('x.')
        try:
            nnnx = int(flush_name[:xpos]) if xpos > 0 else 1
        except ValueError:
            nnnx = 1
        if known_sha256s is not None and new_sha in known_sha256s:
            n_existing  += 1
            sum_existing += nnnx
        else:
            n_novel  += 1
            sum_novel += nnnx

    try:
        with open(fasta_path, "rb") as fh:
            for raw in fh:
                line = _decode_fasta_line(raw).rstrip("\r\n")
                if not line:
                    continue
                if line[0] == ">":
                    if name is not None and is_nnnx_file:
                        _chk(name, seq_parts)
                    hdr  = line[1:]
                    toks = hdr.split()
                    name = toks[0] if toks else ""
                    seq_parts = []
                    if is_nnnx_file is None:
                        if _extract_sha256_from_id(name) is None:
                            return None  # GISAID / legacy file — skip
                        is_nnnx_file = True
                else:
                    seq_parts.append(line)
            if name is not None and is_nnnx_file:
                _chk(name, seq_parts)
    except OSError:
        return None

    if is_nnnx_file is None:
        return None  # empty file
    return n_existing, sum_existing, n_novel, sum_novel, altered_id_shas, mismatch_seqs


def _classify_mismatches(
        mismatch_seqs: dict[str, str],
        parent_path: str,
) -> dict[str, int]:
    """Classify mismatch sequences against original sequences in *parent_path*.

    *mismatch_seqs* maps id-embedded sha256 → depadded post-alignment sequence.
    Scans *parent_path* (the direct-parent FASTA with NNNNx.sha256hex IDs) and
    for each record whose sha256 matches a key in *mismatch_seqs*, compares
    the original (depadded) and current sequences:

      right_clipped     – current is a prefix of original (right end trimmed)
      left_clipped      – current is a suffix of original (left end trimmed)
      both_ends_clipped – current is an interior substring of original
      other             – current is not a substring of original at all
                          (internal bases were altered, not just ends clipped)

    Early-exits as soon as all sha256s in *mismatch_seqs* are resolved, so for
    rare mismatches (e.g. 130 out of 800 K records) only a small fraction of
    the parent FASTA is actually read.
    """
    cats: dict[str, int] = {
        'right_clipped': 0, 'left_clipped': 0, 'both_ends_clipped': 0, 'other': 0,
    }
    if not mismatch_seqs:
        return cats

    remaining = dict(mismatch_seqs)  # sha256 → current depadded seq; shrinks as resolved
    name: str | None = None
    parts: list = []

    def _resolve(orig_name: str, orig_parts: list) -> None:
        sha = _extract_sha256_from_id(orig_name)
        if sha not in remaining:
            return
        orig_seq = (''.join(orig_parts)
                    .replace('\r', '').replace('\n', '').replace('-', '').upper())
        curr_seq = remaining.pop(sha)
        if orig_seq.startswith(curr_seq):
            cats['right_clipped'] += 1
        elif orig_seq.endswith(curr_seq):
            cats['left_clipped'] += 1
        elif curr_seq in orig_seq:
            cats['both_ends_clipped'] += 1
        else:
            cats['other'] += 1

    try:
        with open(parent_path, 'rb') as fh:
            for raw in fh:
                line = _decode_fasta_line(raw).rstrip('\r\n')
                if not line:
                    continue
                if line[0] == '>':
                    if name is not None:
                        _resolve(name, parts)
                    name = line[1:].split()[0] if line[1:].split() else ''
                    parts = []
                    if not remaining:
                        break  # all resolved — stop early
                else:
                    parts.append(line)
            if name is not None and remaining:
                _resolve(name, parts)
    except OSError as exc:
        print(f'  Warning: _classify_mismatches: cannot read {parent_path}: {exc}',
              file=sys.stderr)
    return cats


def _write_original_descr_lines(
        root_descr_tsv: str,
        files: list[str],
        parent_map: dict[int, int],
        sha256_sets: list[tuple[set[str], int]],
        verify_data: list[tuple | None],
) -> None:
    """Stream *root_descr_tsv* once, writing per-step traceability TSV files.

    For each child file (those with a direct parent in the pipeline) three
    output files are written alongside the FASTA:

      *.sha256_to_original_descr_lines_of_survived.tsv
          sha256s present in this child file whose sequence is UNCHANGED
          (altered sequences are excluded — they appear only in _of_survived_altered).
      *.sha256_to_original_descr_lines_of_discarded.tsv
          sha256s in the parent that are absent from this child.
      *.sha256_to_original_descr_lines_of_survived_altered.tsv
          sha256s whose embedded sequence sha256 changed (the pipeline step
          modified the sequence).  Mutually exclusive with _of_survived.

    Output columns: sha256hex<TAB>original_gisaid_header
    (the count column from the root TSV is omitted; columns 3 and onwards
    of root_descr_tsv are included verbatim as the header field).

    All output files are opened lazily and written in one streaming pass.
    """
    # Build per-file category sets ─────────────────────────────────────────────
    # Each element: (survivor_set | None, discarded_set | None, altered_set | None, stem)
    file_cats: list[tuple] = []
    for i, f in enumerate(files):
        p = parent_map.get(i)
        child_sha256s, _ = sha256_sets[i]
        stem = _strip_fasta_suffix(f)
        surv_set = child_sha256s if (p is not None and child_sha256s) else None
        disc_set: set[str] | None = None
        if p is not None:
            parent_sha256s, _ = sha256_sets[p]
            diff = parent_sha256s - child_sha256s
            disc_set = diff if diff else None
        vd = verify_data[i]
        alt_set: set[str] | None = (
            vd[4] if (vd is not None and len(vd) > 4 and vd[4]) else None
        )
        file_cats.append((surv_set, disc_set, alt_set, stem))

    # Open output files lazily ─────────────────────────────────────────────────
    open_fhs: dict[str, object] = {}

    def _fh(path: str) -> object:
        if path not in open_fhs:
            open_fhs[path] = open(path, 'w', encoding='utf-8')  # type: ignore[assignment]  # pylint: disable=consider-using-with
        return open_fhs[path]

    # Stream root TSV, dispatch each sha256 to matching output files ────────────
    n_written = 0
    try:
        print(f"\n  streaming {os.path.basename(root_descr_tsv)}"
              " for traceability files …", file=sys.stderr, flush=True)
        with open(root_descr_tsv, 'r', encoding='utf-8', errors='replace') as rf:
            for line in rf:
                parts = line.rstrip('\n').split('\t', 2)
                if len(parts) < 3:
                    continue
                sha256, _, header = parts
                # Output: sha256 TAB original_header (cols 3+ from root TSV)
                out_line = sha256 + '\t' + header + '\n'
                matched = False
                for surv_set, disc_set, alt_set, stem in file_cats:
                    # Survivors: present in child AND not altered (mutually exclusive
                    # with _altered — an altered record goes only into _altered).
                    if surv_set and sha256 in surv_set and not (alt_set and sha256 in alt_set):
                        _fh(stem + '.sha256_to_original_descr_lines_of_survived.tsv').write(out_line)  # type: ignore[union-attr]
                        matched = True
                    if disc_set and sha256 in disc_set:
                        _fh(stem + '.sha256_to_original_descr_lines_of_discarded.tsv').write(out_line)  # type: ignore[union-attr]
                        matched = True
                    if alt_set and sha256 in alt_set:
                        _fh(stem + '.sha256_to_original_descr_lines_of_survived_altered.tsv').write(out_line)  # type: ignore[union-attr]
                        matched = True
                if matched:
                    n_written += 1
    except OSError as exc:
        print(f"  Warning: could not read root TSV {root_descr_tsv}: {exc}",
              file=sys.stderr)
    finally:
        for fh in open_fhs.values():
            fh.close()  # type: ignore[union-attr]
        if open_fhs:
            print(f"  Info: wrote original GISAID headers for {n_written:,} unique"
                  f" sha256s into {len(open_fhs):,} traceability TSV file(s).",
                  file=sys.stderr, flush=True)


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


# ── Phase 1 helper ────────────────────────────────────────────────────────────

def _scan_primary_file(
        f: str,
        idx: int,
        n_total: int,
        search_path: str,
) -> dict:
    """Scan one primary FASTA and return all stats as a plain dict.

    Designed to run either sequentially or concurrently (ThreadPoolExecutor).
    All stderr output is collected in ``result['lines']`` and must be printed
    by the caller via ``_emit_lines()`` so that concurrent scans each print
    their block atomically without interleaving.

    Returned keys:
        'idx'        – original file index (int)
        'row'        – (display, mtime_s, n_rec, n_sum)
        'sha256_set' – (set[str], int)  sha256 set + n_legacy_ids count
        'prot_unique'– None | (int, int)  distinct proteins + NNNNx sum
        'lines'      – list[str]  stderr progress lines (print atomically)
    """
    lines: list[str] = []
    display  = os.path.relpath(f, search_path)
    mtime_s  = datetime.datetime.fromtimestamp(
        os.path.getmtime(f)).strftime('%Y-%m-%d %H:%M')
    tag      = f"  [{idx + 1}/{n_total}]"
    sz_str   = _fmt_size(os.path.getsize(f))
    is_prot  = _is_prot_file(f)
    kind     = 'protein FASTA' if is_prot else 'FASTA'
    lines.append(f"{_ts()}{tag} {display}  ({sz_str}, {kind})")
    lines.append("        counting records \u2026")
    n_rec = _count_records(f)
    lines.append("        summing NNNNx counts \u2026")
    n_sum = _sum_nnnx_counts(f)
    lines.append("        collecting sha256 IDs \u2026")
    sha256_set_result = _collect_sha256_set(f)
    sha_set, n_legacy = sha256_set_result
    if is_prot:
        lines.append("        counting unique protein sequences \u2026")
    prot_u = _count_prot_unique(f) if is_prot else None
    lines.append(
        f"{_ts()}        done: {n_rec:,} records, NNNNx sum={n_sum:,}, "
        f"{len(sha_set):,} unique sha256s"
        + (f", {n_legacy:,} legacy IDs" if n_legacy else "") + "."
    )
    # ── internal consistency: sha256 set size == record count ─────────────────
    # For NNNNx files every record has a unique sha256 ID, so the set
    # cardinality must equal the record count.  A mismatch signals
    # duplicate sha256 IDs or a stale / mismatched TSV.
    if n_legacy == 0 and len(sha_set) != n_rec:
        lines.append(
            f"        INTERNAL CHECK FAILED: sha256 set size {len(sha_set):,} "
            f"\u2260 record count {n_rec:,} (duplicate sha256 IDs or stale TSV?)"
        )
    return {
        'idx':         idx,
        'row':         (display, mtime_s, n_rec, n_sum),
        'sha256_set':  sha256_set_result,
        'prot_unique': prot_u,
        'lines':       lines,
    }


def _emit_lines(lines: list[str]) -> None:
    """Print a block of stderr lines atomically (thread-safe via _stderr_lock)."""
    with _stderr_lock:
        for ln in lines:
            print(ln, file=sys.stderr, flush=True)


# ── main ──────────────────────────────────────────────────────────────────────

def main() -> None:
    """Entry point: parse args, scan files, print table, optionally show discard stats."""
    _start_ts = datetime.datetime.now().strftime('[%Y-%m-%d %H:%M:%S]')
    print(
        f"{_start_ts} summarize_fasta_pipeline.py"
        f"  version {VERSION}  git:{_GIT_VERSION}"
        f"  invoked: {' '.join(sys.argv)}",
        file=sys.stderr,
    )
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
    # verify_sha256 is ON by default; pass --no-verify-sha256 to disable.
    # The old --verify-sha256 flag is accepted silently for backwards compat.
    verify_sha256        = '--no-verify-sha256'             not in args
    write_original_descr = '--write-original-descr-lines'  in args
    # Classify mismatches as clip-only (right/left/both-ends) vs internal change.
    # Requires an extra scan of each parent FASTA; opt-in only.
    classify_mismatches  = '--classify-mismatches'          in args
    # --jobs N: parallel workers for Phase 0 (sha256) and Phase 1 (scan).
    _jobs_str = (
        next((a.split('=', 1)[1] for a in args if a.startswith('--jobs=')), None)
        or next((args[i + 1] for i, a in enumerate(args)
                 if a == '--jobs' and i + 1 < len(args)), None)
    )
    try:
        jobs = max(1, int(_jobs_str)) if _jobs_str else 1
    except ValueError:
        jobs = 1

    # ── collect matching files ────────────────────────────────────────────────
    found: set[str] = set()
    for suffix in FASTA_SUFFIXES:
        for pat in (
            os.path.join(search_path, prefix + '*' + suffix),
            os.path.join(search_path, '**', prefix + '*' + suffix),
        ):
            found.update(glob.glob(pat, recursive=True))
    # Also pick up protein FASTA files (.prot.fasta already matched above;
    # .prot and .faa have no .fasta extension so need separate globs).
    for suffix in ('.prot', '.faa'):
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
    print(f"{_ts()}Found {len(files):,} file(s).\n", file=sys.stderr)

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

    # ── detect identical files (size + sha256) to skip redundant scans ──────
    # Two files with the same raw content produce identical records, sha256 IDs,
    # and prot-unique counts.  We process the first (primary) and copy results.
    # Only triggered for size-tied pairs — computing sha256 of a 33 GB file
    # takes ~3 min but saves ~38 min of Python-level _verify_sha256 per twin.
    content_twin: dict[int, int] = {}   # duplicate file-index → primary file-index
    _size_groups: dict[int, list[int]] = {}  # size → [file indices]
    for idx, f in enumerate(files):
        try:
            sz = os.path.getsize(f)
        except OSError:
            continue
        _size_groups.setdefault(sz, []).append(idx)
    _size_tied = {sz: idxs for sz, idxs in _size_groups.items()
                  if len(idxs) >= 2 and sz > 0}
    if _size_tied:
        print(
            f"\n{_ts()}── Phase 0: Identity check (size + sha256) "
            f"─────────────────────────────────────────────────────────",
            file=sys.stderr,
        )
        print(
            f"  {len(_size_tied)} size group(s) with \u22652 files will be sha256-checked;"
            f" identical files share scan results.",
            file=sys.stderr,
        )
    def _sha256_worker(i: int) -> tuple[int, str]:
        """Compute sha256 for files[i]; print start message thread-safely."""
        _name = os.path.basename(files[i])
        with _stderr_lock:
            print(f"{_ts()}    computing sha256: {_name} \u2026",
                  file=sys.stderr, flush=True)
        return i, _fasta_sha256(files[i])

    for sz, idxs in _size_tied.items():
        print(
            f"\n{_ts()}  Size group {_fmt_size(sz)} \u2014 {len(idxs)} file(s):",
            file=sys.stderr,
        )
        # Compute sha256 for every file in this size group.
        # Runs in parallel when --jobs > 1 (I/O-bound; GIL released during read).
        max_w0 = min(jobs, len(idxs))
        idx_to_sha: dict[int, str] = {}
        if max_w0 == 1:
            for _i in idxs:
                _ii, _sha = _sha256_worker(_i)
                idx_to_sha[_ii] = _sha
        else:
            with ThreadPoolExecutor(max_workers=max_w0) as _exc:
                for _ii, _sha in _exc.map(_sha256_worker, idxs):
                    idx_to_sha[_ii] = _sha

        # Determine primary vs twin in sorted index order (reproducible).
        sha_to_primary: dict[str, int] = {}
        for idx in sorted(idxs):
            name = os.path.basename(files[idx])
            sha  = idx_to_sha[idx]
            if sha in sha_to_primary:
                pri_name = os.path.basename(files[sha_to_primary[sha]])
                content_twin[idx] = sha_to_primary[sha]
                print(
                    f"{_ts()}    \u2192 {name} \u2261 {pri_name}"
                    f" (sha256={sha[:16]}\u2026) \u2014 scan results will be reused.",
                    file=sys.stderr, flush=True,
                )
            else:
                sha_to_primary[sha] = idx
                print(f"{_ts()}    \u2192 {name}: sha256={sha[:16]}\u2026 (primary)",
                      file=sys.stderr, flush=True)

    # ── gather per-file data ─────────────────────────────────────────────────
    print(
        f"\n{_ts()}── Phase 1: Gather per-file statistics "
        f"──────────────────────────────────────────────────────────────",
        file=sys.stderr,
    )
    n_twins = len(content_twin)
    n_primaries = len(files) - n_twins
    print(
        f"  {len(files)} file(s) total: {n_primaries} to scan, "
        f"{n_twins} to reuse from content-identical primary.",
        file=sys.stderr,
    )
    # Pre-size all accumulator lists so index-based (thread-safe) writes are
    # possible in the parallel path (no .append() ordering constraints).
    rows:          list = [None] * len(files)   # (display, mtime_s, n_rec, n_sum)
    sha256_sets:   list = [None] * len(files)   # (set[str], int) sha256 set + n_legacy
    verify_data:   list = [None] * len(files)   # filled in Phase 2
    classify_data: list = [None] * len(files)   # filled after Phase 2
    prot_unique:   list = [None] * len(files)   # distinct protein seq count (prot files only)

    primary_indices = [i for i in range(len(files)) if i not in content_twin]

    if jobs <= 1 or len(primary_indices) <= 1:
        # Sequential (original behaviour — predictable ordered output).
        for idx in primary_indices:
            result = _scan_primary_file(files[idx], idx, len(files), search_path)
            _emit_lines(result['lines'])
            rows[idx]        = result['row']
            sha256_sets[idx] = result['sha256_set']
            prot_unique[idx] = result['prot_unique']
    else:
        # Parallel: all primaries submitted together to the thread pool.
        # Each file’s progress block is printed atomically as it completes
        # (order non-deterministic, but per-file output never interleaves).
        max_w = min(jobs, len(primary_indices))
        print(
            f"  Parallel Phase 1: scanning {len(primary_indices)} primary file(s) "
            f"with {max_w} worker(s).",
            file=sys.stderr,
        )
        with ThreadPoolExecutor(max_workers=max_w) as exc:
            futs = {
                exc.submit(
                    _scan_primary_file, files[idx], idx, len(files), search_path
                ): idx
                for idx in primary_indices
            }
            for fut in as_completed(futs):
                result = fut.result()
                _emit_lines(result['lines'])
                rows[result['idx']]        = result['row']
                sha256_sets[result['idx']] = result['sha256_set']
                prot_unique[result['idx']] = result['prot_unique']

    # Resolve twins (instant — copy results from already-scanned primary).
    for idx, f in enumerate(files):
        pri = content_twin.get(idx)
        if pri is None:
            continue
        display = os.path.relpath(f, search_path)
        mtime_s = datetime.datetime.fromtimestamp(
            os.path.getmtime(f)).strftime('%Y-%m-%d %H:%M')
        _, _, n_rec, n_sum = rows[pri]
        pri_display = os.path.relpath(files[pri], search_path)
        rows[idx]        = (display, mtime_s, n_rec, n_sum)
        sha256_sets[idx] = sha256_sets[pri]
        prot_unique[idx] = prot_unique[pri]
        print(
            f"{_ts()}  [{idx+1}/{len(files)}] {display}\n"
            f"        ↳ reusing results from [{pri+1}] {pri_display}"
            f" ({n_rec:,} records, NNNNx sum={n_sum:,}).",
            file=sys.stderr, flush=True,
        )

    # sha256_sets is still needed in the table-printing loop for the
    # 'Novel sha256s' column.  It is freed after that loop completes.

    if verify_sha256:
        print(
            f"\n{_ts()}── Phase 2: sha256 integrity verification "
            f"────────────────────────────────────────────────────────",
            file=sys.stderr,
        )
        n_verify = sum(1 for i, f in enumerate(files)
                       if not _is_prot_file(f)
                       and not (content_twin.get(i) is not None
                                and (parent_map.get(i) == content_twin[i]
                                     or content_twin.get(parent_map.get(i)) == parent_map.get(content_twin[i])
                                     or parent_map.get(i) == parent_map.get(content_twin[i]))))
        print(
            f"  Verifying sha256 integrity for {n_verify} file(s) "
            f"(protein files and content twins with equivalent parent skipped).",
            file=sys.stderr,
        )
        for idx, f in enumerate(files):
            if _is_prot_file(f):
                display = os.path.relpath(f, search_path)
                print(
                    f"  [{idx+1}/{len(files)}] {display}: skipped — protein file "
                    f"(ID contains DNA sha256, not protein sha256).",
                    file=sys.stderr, flush=True,
                )
                continue
            pri = content_twin.get(idx)
            if pri is not None:
                p = parent_map.get(idx)
                p_pri = parent_map.get(pri)
                if (p == pri
                        or content_twin.get(p) == p_pri
                        or p == p_pri):
                    display = os.path.relpath(f, search_path)
                    pri_display = os.path.relpath(files[pri], search_path)
                    print(
                        f"  [{idx+1}/{len(files)}] {display}: skipped — content-identical "
                        f"to [{pri+1}] {pri_display} and parent is also equivalent; "
                        f"no mismatches possible.",
                        file=sys.stderr, flush=True,
                    )
                    verify_data[idx] = (0, 0, 0, 0, set(), {})
                    continue
                # Different parent context — fall through to normal verify below.
            display = os.path.relpath(f, search_path)
            p = parent_map.get(idx)
            parent_sha256s = sha256_sets[p][0] if p is not None else None
            parent_display = os.path.relpath(files[p], search_path) if p is not None else '(no parent)'
            print(
                f"{_ts()}  [{idx+1}/{len(files)}] {display}: verifying sha256 IDs "
                f"against parent [{p+1 if p is not None else '?'}] {parent_display} …",
                file=sys.stderr, flush=True,
            )
            vd = _verify_sha256(f, parent_sha256s)
            verify_data[idx] = vd
            if vd is not None and (vd[0] > 0 or vd[2] > 0):
                print(
                    f"{_ts()}    Warning: {display}:"
                    + (f" {vd[0]:,} record(s) sha256→existing" if vd[0] else "")
                    + (f" {vd[2]:,} record(s) sha256→novel" if vd[2] else "")
                    + f" (NNNNx: {vd[1]+vd[3]:,} total)",
                    file=sys.stderr,
                )
            else:
                print(
                    f"{_ts()}    OK: all sha256 IDs match their sequences and parent set.",
                    file=sys.stderr, flush=True,
                )
            if classify_mismatches and vd is not None and p is not None:
                mismatch_seqs_vd: dict[str, str] = vd[5] if len(vd) > 5 else {}
                if mismatch_seqs_vd:
                    print(
                        f"    classifying {len(mismatch_seqs_vd):,} mismatch(es) "
                        f"(end-clipping vs internal change) …",
                        file=sys.stderr, flush=True,
                    )
                    classify_data[idx] = _classify_mismatches(
                        mismatch_seqs_vd, files[p]
                    )

    # ── phase 2: compute discard stats for all pairs ─────────────────────────
    # Runs before table printing so the numbers can appear as proper columns.
    discard_data: dict[int, tuple[int, int]] = {}  # child_idx -> (n_ids, nnnx_sum)
    if do_discard:
        print(
            f"\n── Phase 3: Discard statistics (parent→child pairs) "
            f"─────────────────────────────────────────────────────",
            file=sys.stderr,
        )
        print(
            f"  Computing discarded-sequence stats for {len(parent_map)} parent→child pair(s).",
            file=sys.stderr,
        )
        for i, f_child in enumerate(files):
            p = parent_map.get(i)
            if p is not None:
                p_display = os.path.relpath(files[p], search_path)
                c_display = os.path.relpath(f_child, search_path)
                print(
                    f"  [{i+1}/{len(files)}] {p_display} → {c_display} …",
                    file=sys.stderr, flush=True,
                )
                n_d, s_d = _compute_discard_stats(
                    files[p], f_child,
                    add_checksums=add_checksums,
                    full_fasta_header=full_fasta_header,
                    save_discard_list=save_discard_list,
                    verbose=verbose,
                )
                if n_d is not None:
                    discard_data[i] = (n_d, s_d)
                    print(
                        f"    {n_d:,} discarded original ID(s), NNNNx sum={s_d:,}.",
                        file=sys.stderr, flush=True,
                    )

    # ── print table ──────────────────────────────────────────────────────────
    col_file = max(max(len(r[0]) for r in rows), len("File"))
    sep, w_num, w_delta, w_ts = "  ", 14, 16, 16
    w_disc1  = len("'Discarded original FASTA IDs'")   # 30
    w_disc2  = len("'Sum of discarded sequences'")       # 28
    w_novel  = len("'Novel sha256s'")                   # 15
    w_chg1   = max(len("'Seq clipped(dup)'"), w_num)
    w_chg2   = max(len("'NNNNx clipped(dup)'"), w_num)
    w_chg3   = max(len("'Seq clipped(new)'"), w_num)
    w_chg4   = max(len("'NNNNx clipped(new)'"), w_num)
    w_clip   = max(len("'Altered due to end-clipping'"), w_num)
    w_itrn   = max(len("'Altered inside the sequence'"), w_num)
    w_prot   = max(len("'Prot unique'"), w_num)
    w_protn  = max(len("'Prot NNNNx sum'"), w_num)

    _hdr_disc1 = "'Discarded original FASTA IDs'"
    _hdr_disc2 = "'Sum of discarded sequences'"
    _hdr_nnnx  = "'Sum of NNNNx'"
    _hdr_drec  = '\u0394Records'
    _hdr_dsum  = '\u0394SumToParent'
    _hdr_novel = "'Novel sha256s'"
    _hdr_chg1  = "'Seq clipped(dup)'"
    _hdr_chg2  = "'NNNNx clipped(dup)'"
    _hdr_chg3  = "'Seq clipped(new)'"
    _hdr_chg4  = "'NNNNx clipped(new)'"
    _hdr_clip  = "'Altered due to end-clipping'"
    _hdr_itrn  = "'Altered inside the sequence'"
    _hdr_prot  = "'Prot unique'"
    _hdr_protn = "'Prot NNNNx sum'"
    verify_cols_hdr = (
        f"{sep}{_hdr_chg1:>{w_chg1}}{sep}{_hdr_chg2:>{w_chg2}}"
        f"{sep}{_hdr_chg3:>{w_chg3}}{sep}{_hdr_chg4:>{w_chg4}}"
        + (f"{sep}{_hdr_clip:>{w_clip}}{sep}{_hdr_itrn:>{w_itrn}}"
           if classify_mismatches else "")
        if verify_sha256 else ""
    )
    disc_header = (
        f"{sep}{_hdr_disc1:>{w_disc1}}{sep}{_hdr_disc2:>{w_disc2}}"
        if do_discard else ""
    )
    header = (
        f"{'File':<{col_file}}{sep}"
        f"{'Modified':<{w_ts}}{sep}"
        f"{'Records':>{w_num}}{sep}{_hdr_drec:>{w_delta}}{sep}"
        f"{_hdr_nnnx:>{w_num}}{sep}{_hdr_dsum:>{w_delta}}{sep}"
        f"{_hdr_novel:>{w_novel}}{sep}{_hdr_prot:>{w_prot}}{sep}{_hdr_protn:>{w_protn}}"
        + verify_cols_hdr
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

        # ── novel sha256s column ──────────────────────────────────────────────
        _em = '\u2014'  # em dash — pre-assigned to avoid backslash in f-string (Python < 3.12)
        if p is not None:
            child_sha_set, child_legacy = sha256_sets[i]
            parent_sha_set, parent_legacy = sha256_sets[p]
            parent_total = rows[p][2]  # record count of parent
            if parent_legacy == parent_total:
                # Parent has only plain (non-NNNNx) IDs — sha256 set is empty,
                # comparison is not meaningful.
                novel_col = f"{sep}{_em:>{w_novel}}"
            else:
                n_novel = len(child_sha_set - parent_sha_set)
                novel_str = f"{n_novel:,}"
                if child_legacy > 0:
                    novel_str += "+"  # '+' = some records in child also lack embedded sha256
                novel_col = f"{sep}{novel_str:>{w_novel}}"
        else:
            novel_col = f"{sep}{_em:>{w_novel}}"

        # ── protein-unique column ─────────────────────────────────────────────
        pu_result = prot_unique[i]  # None | (n_unique, nnnx_sum)
        if pu_result is None:
            prot_col = f"{sep}{_em:>{w_prot}}{sep}{_em:>{w_protn}}"
        else:
            n_pu, nnnx_pu = pu_result
            prot_col = f"{sep}{n_pu:>{w_prot},}{sep}{nnnx_pu:>{w_protn},}"

        # ── sha256 verification columns ───────────────────────────────────────
        if verify_sha256:
            vd = verify_data[i]
            if vd is None:
                verify_cols = (
                    f"{sep}{_em:>{w_chg1}}{sep}{_em:>{w_chg2}}"
                    f"{sep}{_em:>{w_chg3}}{sep}{_em:>{w_chg4}}"
                    + (f"{sep}{_em:>{w_clip}}{sep}{_em:>{w_itrn}}"
                       if classify_mismatches else "")
                )
            else:
                n_ex, s_ex, n_nv, s_nv, *_ = vd
                if classify_mismatches:
                    cd = classify_data[i]
                    if cd is None:
                        clip_str = _em
                        itrn_str = _em
                    else:
                        n_clip = cd['right_clipped'] + cd['left_clipped'] + cd['both_ends_clipped']
                        clip_str = f"{n_clip:,}"
                        itrn_str = f"{cd['other']:,}"
                    classify_cols = f"{sep}{clip_str:>{w_clip}}{sep}{itrn_str:>{w_itrn}}"
                else:
                    classify_cols = ""
                verify_cols = (
                    f"{sep}{n_ex:>{w_chg1},}{sep}{s_ex:>{w_chg2},}"
                    f"{sep}{n_nv:>{w_chg3},}{sep}{s_nv:>{w_chg4},}"
                    + classify_cols
                )
        else:
            verify_cols = ""

        if do_discard and p is not None:
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
            + novel_col
            + prot_col
            + verify_cols
            + disc_cols
            + (f"  ({parent_label})" if parent_label else '')
        )



    print(rule)

    # ── write original-header traceability TSV files ───────────────────────────
    if write_original_descr:
        root_stem = _strip_fasta_suffix(files[0])
        root_descr_tsv = root_stem + '.sha256_to_descr_lines.tsv'
        if not os.path.exists(root_descr_tsv):
            root_descr_tsv = root_stem + '.sha256_to_ids.tsv'
        if os.path.exists(root_descr_tsv):
            _write_original_descr_lines(
                root_descr_tsv, files, parent_map, sha256_sets, verify_data,
            )
        else:
            print(f"  Warning: --write-original-descr-lines: root mapping TSV not found"
                  f" ({root_stem}.sha256_to_descr_lines.tsv). Run with"
                  " --full-fasta-header first to build it.", file=sys.stderr)

    del sha256_sets  # no longer needed; free RAM

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
