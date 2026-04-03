#!/usr/bin/env python3
"""Diagnostic: classify sequences changed by alignment.

For each record in an aligned FASTA with NNNNx.sha256hex IDs, compare:
  sha256(depadded_aligned_seq.upper())  vs  sha256 stored in the ID

A mismatch means the alignment trimmed or modified the sequence so that
its depadded form no longer matches the original pre-alignment sequence.

When ``--original-fasta`` is provided (the *.counts.fasta that was fed into
the alignment step), mismatching records are further classified by comparing
the post-alignment depadded sequence against the pre-alignment original
(also depadded and uppercased):

  right_clipped     : current seq is a PREFIX of the original  (right end trimmed)
  left_clipped      : current seq is a SUFFIX of the original  (left end trimmed)
  both_ends_clipped : current seq is an INTERIOR substring of the original
                      (both ends trimmed, internal content intact)
  other             : current is not a substring of original at all
                      — internal bases were modified, not just ends clipped

Two-pass strategy to stay memory-efficient on large files:
  Pass 1: scan aligned FASTA → collect {sha256_in_id: depadded_current_seq}
           only for mismatching records (typically very few).
  Pass 2: scan original FASTA → for each matching sha256, compare sequences.

Usage::

    # Basic: count mismatches only
    check_alignment_trimming.py \\
        --infilename=spikenuc1207.no_junk.counts.3822.clean.fasta

    # With classification: also explain what changed
    check_alignment_trimming.py \\
        --infilename=spikenuc1207.no_junk.counts.3822.clean.fasta \\
        --original-fasta=spikenuc1207.no_junk.counts.fasta
"""
import hashlib
import sys
import argparse

VERSION = "202604022300"

parser = argparse.ArgumentParser(
    description=__doc__,
    formatter_class=argparse.RawDescriptionHelpFormatter,
)
parser.add_argument(
    "--infilename", required=True,
    help="Aligned FASTA with NNNNx.sha256hex IDs (post-alignment file to audit).",
)
parser.add_argument(
    "--original-fasta", dest="original_fasta", default="",
    help=(
        "Pre-alignment *.counts.fasta (NNNNx.sha256hex IDs).  When provided, "
        "each sha256 mismatch is classified as right_clipped / left_clipped / "
        "both_ends_clipped / other by comparing the depadded post-alignment "
        "sequence against the original depadded sequence from this file."
    ),
)
parser.add_argument("--debug", action="store_true",
                    help="Print per-record mismatch details.")
parser.add_argument("--version", action="version", version=f"%(prog)s {VERSION}")
args = parser.parse_args()


def _extract_sha256(record_id: str) -> str | None:
    """Return 64-char hex sha256 from a NNNNx.sha256hex ID, or None."""
    xdot = record_id.find('x.')
    if xdot >= 0:
        candidate = record_id[xdot + 2:].split()[0]
        if len(candidate) == 64 and all(c in '0123456789abcdefABCDEF' for c in candidate):
            return candidate.lower()
    return None


def _depad(seq: str) -> str:
    """Strip dashes and uppercase — the canonical pre-hash normalisation."""
    return seq.replace('-', '').upper()


def _iter_fasta(path: str):
    """Yield (name, depadded_seq) for each record."""
    name = None
    parts: list[str] = []
    with open(path, 'r', encoding='utf-8', errors='replace') as fh:
        for line in fh:
            line = line.rstrip('\r\n')
            if not line:
                continue
            if line.startswith('>'):
                if name is not None:
                    yield name, _depad(''.join(parts))
                name = line[1:].split()[0]
                parts = []
            else:
                parts.append(line)
    if name is not None:
        yield name, _depad(''.join(parts))


def _classify(current: str, original: str) -> str:
    """Classify how *current* (post-alignment, depadded) relates to *original*.

    Categories (mutually exclusive, checked in order):
      right_clipped     – current is a non-empty prefix of original
      left_clipped      – current is a non-empty suffix of original
      both_ends_clipped – current is an interior substring of original
                          (neither a prefix nor a suffix)
      other             – current is not a substring of original at all;
                          internal bases were altered, not just ends clipped
    """
    if not current or not original:
        return 'other'
    if original.startswith(current):
        return 'right_clipped'
    if original.endswith(current):
        return 'left_clipped'
    if current in original:
        return 'both_ends_clipped'
    return 'other'


def main() -> None:
    """Scan aligned FASTA and report (and optionally classify) sha256 mismatches."""
    print(
        f"check_alignment_trimming.py  version {VERSION}"
        f"  invoked: {' '.join(sys.argv)}",
        file=sys.stderr,
    )

    # ── Pass 1: detect mismatches ───────────────────────────────────────────
    # mismatch_seqs: sha256_in_id → depadded post-alignment sequence
    mismatch_seqs: dict[str, str] = {}
    n_total = n_with_sha = n_mismatch = 0

    for rec_id, dep_seq in _iter_fasta(args.infilename):
        n_total += 1
        sha_in_id = _extract_sha256(rec_id)
        if sha_in_id is None:
            continue
        n_with_sha += 1
        sha_computed = hashlib.sha256(dep_seq.encode()).hexdigest()
        if sha_computed != sha_in_id:
            n_mismatch += 1
            mismatch_seqs[sha_in_id] = dep_seq
            if args.debug:
                print(
                    f"MISMATCH id={rec_id}  id_sha={sha_in_id}"
                    f"  seq_sha={sha_computed}  depadded_len={len(dep_seq)}",
                    file=sys.stderr,
                )
        if n_total % 500_000 == 0:
            print(f"Info: {n_total:,} records processed, "
                  f"{n_mismatch:,} mismatches so far", file=sys.stderr)

    # ── Pass 2 (optional): classify mismatches against the original FASTA ──
    categories: dict[str, int] = {
        'right_clipped': 0,
        'left_clipped': 0,
        'both_ends_clipped': 0,
        'other': 0,
    }
    n_classified = 0

    if args.original_fasta and mismatch_seqs:
        print(f"\nInfo: classifying {len(mismatch_seqs):,} mismatches against "
              f"{args.original_fasta} …", file=sys.stderr)
        remaining = dict(mismatch_seqs)  # copy; remove as we find matches
        for orig_id, orig_seq in _iter_fasta(args.original_fasta):
            sha_orig = _extract_sha256(orig_id)
            if sha_orig is None:
                # Legacy NNNNx ID — compute sha256 from sequence content.
                sha_orig = hashlib.sha256(orig_seq.encode()).hexdigest()
            if sha_orig not in remaining:
                continue
            curr_seq = remaining.pop(sha_orig)
            cat = _classify(curr_seq, orig_seq)
            categories[cat] += 1
            n_classified += 1
            if args.debug:
                print(
                    f"  sha={sha_orig[:16]}…  cat={cat}"
                    f"  orig_len={len(orig_seq)}  curr_len={len(curr_seq)}",
                    file=sys.stderr,
                )
            if not remaining:
                break  # all mismatches resolved; no need to read further

        if remaining:
            print(
                f"  Warning: {len(remaining):,} mismatching sha256(s) not found "
                f"in {args.original_fasta} — original sequences missing or "
                "from a different pipeline step.",
                file=sys.stderr,
            )

    # ── Summary ─────────────────────────────────────────────────────────────
    print("\nSummary")
    print(f"  Total records           : {n_total:,}")
    print(f"  Records with sha256 ID  : {n_with_sha:,}")
    print(f"  sha256 mismatches       : {n_mismatch:,}"
          "  (sequences trimmed/changed by alignment)")
    if n_with_sha:
        print(f"  Match rate              : "
              f"{(n_with_sha - n_mismatch) / n_with_sha * 100:.4f}%")

    if args.original_fasta:
        print(f"\nClassification vs {args.original_fasta} ({n_classified:,} resolved):")
        width = max(len(k) for k in categories) + 2
        for cat, count in categories.items():
            print(f"  {cat:<{width}}: {count:,}")
        if n_mismatch > n_classified:
            print(f"  (unresolved            : {n_mismatch - n_classified:,}"
                  "  — sha256 absent from original FASTA)")


if __name__ == "__main__":
    main()
