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

Cache priority (fastest first):
  1. <child_base>.discarded_original_ids.txt  newer than parent + child FASTA
  2. <parent_base>.sha256_to_ids.tsv          newer than parent FASTA
  3. Full FASTA scan via --original-infilename (slow fallback)

Usage:
    summarize_fasta_pipeline.py <search_path> <filename_prefix> [options]

Options:
    --no-discard-stats   Skip the per-step discard-statistics calls entirely.
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

Example:
    summarize_fasta_pipeline.py . spikenuc1207.native2ascii.no_junk
    summarize_fasta_pipeline.py .. spikenuc1207.native2ascii.no_junk --no-discard-stats
"""

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

def _fresh_discarded_txt(child_path: str, parent_path: str) -> str | None:
    """Return path to a valid .discarded_original_ids.txt if it exists and is
    newer than BOTH the child and parent FASTA files, else None."""
    candidate = _strip_fasta_suffix(child_path) + '.discarded_original_ids.txt'
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


def _build_tsv(fasta_path: str, tsv_path: str) -> dict:
    """Single-pass scan of *fasta_path* to build a sha256→IDs mapping TSV.

    For each record the sha256 is computed over the uppercase, dash-stripped
    sequence — identical to the normalisation used by count_same_sequences.py.

    Output TSV columns (tab-separated, no header):
        sha256hex   count   id_1   id_2   ...

    Returns a dict mapping each original ID to its sha256hex string, which
    callers can use to rewrite the FASTA with sha256 embedded in the IDs.
    """
    mapping: dict = {}   # sha256 -> [count, [id, ...]]
    name = None
    seq_parts: list = []

    def _flush(flush_name: str, flush_parts: list) -> None:
        seq = "".join(flush_parts).rstrip("\r\n").replace("-", "").upper()
        digest = hashlib.sha256(seq.encode()).hexdigest()
        if digest in mapping:
            mapping[digest][0] += 1
            mapping[digest][1].append(flush_name)
        else:
            mapping[digest] = [1, [flush_name]]

    n_in = 0
    with open(fasta_path, "rb") as fh:
        for raw in fh:
            line = _decode_fasta_line(raw).rstrip("\r\n")
            if not line:
                continue
            if line[0] == ">":
                if name is not None:
                    _flush(name, seq_parts)
                    n_in += 1
                toks = line[1:].split()
                name = toks[0] if toks else ""
                seq_parts = []
            else:
                seq_parts.append(line)
        if name is not None:
            _flush(name, seq_parts)
            n_in += 1

    with open(tsv_path, "w", encoding="utf-8") as out:
        for digest, (count, ids) in mapping.items():
            out.write("\t".join([digest, str(count)] + ids) + "\n")

    print(f"    Info: built TSV with {len(mapping):,} unique sequences"
          f" from {n_in:,} records \u2192 {tsv_path}", flush=True)

    # Return inverted mapping: id -> sha256 (one entry per original record).
    return {orig_id: sha for sha, (_, ids) in mapping.items() for orig_id in ids}


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


def _ensure_tsv(parent_path: str, add_checksums: bool = False) -> str | None:
    """Return a fresh .sha256_to_ids.tsv for *parent_path*, generating it by
    scanning the FASTA in-process if one does not already exist or is stale.

    When *add_checksums* is True and the file has legacy NNNNx IDs (no sha256
    embedded), the FASTA is rewritten in-place with NNNNx.sha256hex IDs and
    the original is renamed to .fasta.orig.  This makes all future analyses
    faster since sha256 can be read from the ID instead of being recomputed.

    Returns the TSV path on success, or None if generation failed.
    """
    tsv = _fresh_tsv(parent_path)
    if tsv and not add_checksums:
        return tsv
    # If add_checksums: always rebuild the TSV after enrichment so the TSV
    # reflects the new IDs (which now embed sha256 directly).
    if tsv and add_checksums:
        # Check quickly whether the file already has modern IDs.
        try:
            with open(parent_path, "rb") as fh:
                for raw in fh:
                    line = _decode_fasta_line(raw)
                    if line.startswith(">"):
                        first_id = line[1:].split()[0]
                        if _extract_sha256_from_id(first_id) is not None:
                            return tsv  # already modern format
                        break  # first record lacks sha256 — fall through
        except OSError:
            return tsv  # can't read, just use existing TSV

    candidate = _strip_fasta_suffix(parent_path) + '.sha256_to_ids.tsv'
    print(f"  [auto-generate TSV] scanning {os.path.basename(parent_path)} \u2026",
          flush=True)
    try:
        id_to_sha = _build_tsv(parent_path, candidate)
    except OSError as exc:
        print(f"    [TSV generation failed: {exc}]", flush=True)
        return None

    if add_checksums and _has_legacy_ids(id_to_sha):
        _enrich_fasta(parent_path, id_to_sha)
        # Rebuild the TSV against the new file so IDs match the enriched form.
        print("  [auto-generate TSV] rebuilding TSV from enriched FASTA …", flush=True)
        try:
            _build_tsv(parent_path, candidate)
        except OSError as exc:
            print(f"    [TSV rebuild failed: {exc}]", flush=True)

    return candidate


def _read_discarded_txt_stats(txt_path: str) -> tuple[int, int]:
    """Return (n_ids, nnnx_sum) from a .discarded_original_ids.txt."""
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

def _run_discard_stats(parent_path: str, child_path: str,
                       add_checksums: bool = False) -> None:
    """Show discarded-ID stats for a parent->child pipeline pair.

    Cache priority (fastest first — checked by timestamp):
      1. Existing .discarded_original_ids.txt newer than both FASTAs -> read directly.
      2. Existing .sha256_to_ids.tsv newer than parent FASTA         -> TSV scan.
      3. Fallback: full FASTA scan via --original-infilename (slow).
    """
    child_base  = os.path.basename(child_path)
    parent_base = os.path.basename(parent_path)

    # ── path 1: fresh .discarded_original_ids.txt ──────────────────────────
    txt = _fresh_discarded_txt(child_path, parent_path)
    if txt:
        n_ids, nnnx_sum = _read_discarded_txt_stats(txt)
        print(
            f"  \u21b3 [cached TXT] {os.path.basename(txt)}: "
            f"discarded IDs: {n_ids:,}  NNNNx sum: {nnnx_sum:,}",
            flush=True,
        )
        return

    # ── paths 2 & 3: invoke create_list_of_discarded_sequences.py ──────────
    if not os.path.exists(DISCARD_SCRIPT):
        print(f"  [discard stats] script not found: {DISCARD_SCRIPT}", flush=True)
        return

    # Try to get (or auto-generate) a fresh TSV for the parent — faster than
    # a full FASTA scan and the result is then cached for future calls.
    tsv = _ensure_tsv(parent_path, add_checksums=add_checksums)
    if tsv:
        source_arg   = f'--mapping-outfile={tsv}'
        source_label = f"TSV: {os.path.basename(tsv)}"
    else:
        source_arg   = f'--original-infilename={parent_path}'
        source_label = f"FASTA scan: {parent_base}"

    cmd = [
        sys.executable, DISCARD_SCRIPT,
        f'--infilename={child_path}',
        source_arg,
        '--inverted',
        '--outfile=/dev/null',
    ]
    print(f"  \u21b3 [{source_label}] -> {child_base}", flush=True)
    result = subprocess.run(cmd, capture_output=True, text=True, check=False)
    for line in result.stderr.splitlines():
        if line.startswith(('Info:', 'Warning:', 'Error:')):
            print(f"    {line}", flush=True)
    if result.returncode != 0:
        print(f"    [exit code {result.returncode}]", flush=True)


# ── main ──────────────────────────────────────────────────────────────────────

def main() -> None:
    """Entry point: parse args, scan files, print table, optionally show discard stats."""
    args = sys.argv[1:]
    if len(args) < 2 or '--help' in args or '-h' in args:
        print(__doc__, file=sys.stderr)
        sys.exit(0 if '--help' in args or '-h' in args else 1)

    search_path = args[0]
    prefix      = args[1]
    do_discard     = '--no-discard-stats' not in args
    add_checksums  = '--add-missing-checksums-to-fasta-files' in args

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

    # ── gather per-file data ──────────────────────────────────────────────────
    rows: list[tuple[str, int, int]] = []
    for f in files:
        display = os.path.relpath(f, search_path)
        print(f"  scanning {display} …", file=sys.stderr, flush=True)
        rows.append((display, _count_records(f), _sum_nnnx_counts(f)))

    # ── print table ───────────────────────────────────────────────────────────
    col_file = max(max(len(r[0]) for r in rows), len("File"))
    sep, w_num, w_delta = "  ", 14, 16

    header = (
        f"{'File':<{col_file}}{sep}"
        f"{'Records':>{w_num}}{sep}{'ΔRecords':>{w_delta}}{sep}"
        f"{'Sum of NNNNx':>{w_num}}{sep}{'ΔSum':>{w_delta}}"
    )
    rule = '-' * len(header)
    print()
    print(header)
    print(rule)

    for i, (display, n_rec, n_sum) in enumerate(rows):
        p = parent_map.get(i)
        if p is not None:
            prev_rec, prev_sum = rows[p][1], rows[p][2]
            d_rec = _delta_str(n_rec, prev_rec)
            d_sum = _delta_str(n_sum, prev_sum)
            parent_label = f"vs {os.path.basename(files[p])}"
        else:
            d_rec = d_sum = '—'
            parent_label = ''

        print(
            f"{display:<{col_file}}{sep}"
            f"{n_rec:>{w_num},}{sep}{d_rec:>{w_delta}}{sep}"
            f"{n_sum:>{w_num},}{sep}{d_sum:>{w_delta}}"
            + (f"  ({parent_label})" if parent_label else '')
        )

        if do_discard and p is not None:
            _run_discard_stats(files[p], files[i], add_checksums=add_checksums)

    print(rule)

    # ── overall summary ───────────────────────────────────────────────────────
    if len(rows) >= 2:
        first_rec, first_sum = rows[0][1], rows[0][2]
        last_rec,  last_sum  = rows[-1][1], rows[-1][2]
        print()
        print("Overall change  (last vs first):")
        print(f"  Records : {_delta_str(last_rec, first_rec):>16}  ({_pct_str(last_rec, first_rec)} of first)")
        print(f"  Sum     : {_delta_str(last_sum, first_sum):>16}  ({_pct_str(last_sum, first_sum)} of first)")


if __name__ == '__main__':
    main()
