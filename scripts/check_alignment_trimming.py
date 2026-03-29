#!/usr/bin/env python3
"""Diagnostic: count sequences shortened (or otherwise changed) by alignment.

For each record in an aligned FASTA with NNNNx.sha256hex IDs, compare:
  sha256(depadded_aligned_seq.upper())  vs  sha256 stored in the ID

A mismatch means the alignment trimmed or modified the sequence so that
its depadded form no longer matches the original pre-alignment sequence.

Usage:
    check_alignment_trimming.py --infilename=spikenuc1207.native2ascii.no_junk.counts.clean.fasta
"""
import hashlib
import sys
import argparse

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--infilename", required=True, help="Aligned FASTA with NNNNx.sha256 IDs")
parser.add_argument("--debug", action="store_true")
args = parser.parse_args()


def _extract_sha256(name):
    """Extract 64-char hex sha256 from NNNNx.sha256 ID, or None."""
    dot = name.find(".")
    if dot > 0 and len(name) >= dot + 65:
        candidate = name[dot + 1: dot + 65]
        if all(c in "0123456789abcdefABCDEF" for c in candidate):
            return candidate.lower()
    return None


n_total = n_with_sha = n_mismatch = 0
name = None
seq_parts = []


def _check(name, seq_parts):
    global n_with_sha, n_mismatch
    sha_in_id = _extract_sha256(name)
    if sha_in_id is None:
        return  # legacy ID without sha256 — skip
    n_with_sha += 1
    seq = "".join(seq_parts).replace("-", "").upper()
    sha_computed = hashlib.sha256(seq.encode()).hexdigest()
    if sha_computed != sha_in_id:
        n_mismatch += 1
        if args.debug:
            orig_len = len("".join(seq_parts).replace("-", ""))
            print(f"MISMATCH id={name}  id_sha={sha_in_id}  seq_sha={sha_computed}  depadded_len={orig_len}")


with open(args.infilename, "r", encoding="utf-8") as fh:
    for line in fh:
        line = line.rstrip("\r\n")
        if not line:
            continue
        if line.startswith(">"):
            if name is not None:
                _check(name, seq_parts)
                n_total += 1
                if n_total % 500_000 == 0:
                    print(f"Info: {n_total:,} records processed, {n_mismatch:,} mismatches so far", file=sys.stderr)
            name = line[1:].split()[0]
            seq_parts = []
        else:
            seq_parts.append(line)

if name is not None:
    _check(name, seq_parts)
    n_total += 1

print("\nSummary")
print(f"  Total records     : {n_total:,}")
print(f"  Records with sha256 in ID : {n_with_sha:,}")
print(f"  sha256 mismatches : {n_mismatch:,}  (sequences shortened/changed by alignment)")
print(f"  Match rate        : {(n_with_sha - n_mismatch) / n_with_sha * 100:.4f}%" if n_with_sha else "  N/A")
