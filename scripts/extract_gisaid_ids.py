#!/usr/bin/env python3
"""
Standalone utility designed to cleanly extract `|EPI_ISL_*|` tags from Pipeline `.sha256_to_descr_lines.tsv` mappings.

It supports two main operations:
1) Converting a generic description mapping TSV into a clean `.sha256_to_GISAID_ids.tsv` tree.
   If sequences lack an EPI_ISL_ marker, it falls back to their first-word ID.
2) Expanding a hash subset into a flat sequence array of IDs utilizing the newly formed GISAID TSV.

Usage Examples:
  # 2. Extract into final subset .txt targets
  python extract_gisaid_ids.py --gisaid-tsv-mapping spikenuc1207.sha256_to_GISAID_ids.tsv \
                               --subset-txt spikenuc1207.no_junk.counts.effectively_used_sha256_hashes.txt \
                               --out-subset-txt \\
                                 spikenuc1207.no_junk.counts.effectively_used_original_entry_GISAID_ids.txt

  # 3. Batch backfill across missing Phase 4 extraction pipelines
  # for tsv in *.effectively_used_original_entries.fasta.extraction_counts.tsv; do
  #     python3 extract_gisaid_ids.py --gisaid-tsv-mapping spikenuc1207.no_junk.sha256_to_GISAID_ids.tsv \
  #                                   --subset-txt "$tsv" \
  #                                   --out-subset-txt "${tsv%.fasta.extraction_counts.tsv}_entry_GISAID_ids.txt"
  # done
"""

import argparse
import os
import re
import sys
import time


def build_gisaid_tsv(descr_tsv: str, out_tsv: str):
    """Reads a .sha256_to_descr_lines.tsv and parses out EPI_ISL_ tags to build a GISAID TSV mapping."""
    if not os.path.exists(descr_tsv):
        print(f"Error: map file {descr_tsv} not found.", file=sys.stderr)
        return

    gisaid_re = re.compile(r'\|(EPI_ISL_[^|]+)\|')
    found_any = False

    start_t = time.time()
    n_lines = 0

    with open(descr_tsv, 'r', encoding='utf-8') as f_in, open(out_tsv, 'w', encoding='utf-8') as f_out:
        for line in f_in:
            n_lines += 1
            parts = line.rstrip('\r\n').split('\t')
            if len(parts) < 3:
                continue
            sha = parts[0]
            count = parts[1]
            out_parts = [sha, count]
            for hdr in parts[2:]:
                # Check for the GISAID EPI_ISL identifier
                match = gisaid_re.search(hdr)
                if match:
                    found_any = True
                    out_parts.append(match.group(1))
                else:
                    # Fallback cleanly to the very first string word of the sequence header natively
                    toks = hdr.replace('\\t', '\t').split()
                    first_word = toks[0] if toks else ""
                    out_parts.append(first_word)

            f_out.write('\t'.join(out_parts) + '\n')

    dur = time.time() - start_t
    print(f"Wrote {n_lines:,} mapped hashes to {out_tsv} in {dur:.2f}s.")

    if not found_any:
        print(f"Warning: No |EPI_ISL_*| identifiers were physically located in {descr_tsv}. "
              "Successfully fell back to primary native string sequence IDs for all entries.", file=sys.stderr)


def extract_subset_txt(gisaid_tsv: str, subset_txt: str, out_txt: str):
    """Reads a .discarded_sha256_hashes.txt and writes out all mapped GISAID flat list."""
    if not os.path.exists(gisaid_tsv):
        print(f"Error: GISAID mapping {gisaid_tsv} not found.", file=sys.stderr)
        return
    if not os.path.exists(subset_txt):
        print(f"Error: Subset file {subset_txt} not found.", file=sys.stderr)
        return

    print(f"Loading GISAID TSV mapping into memory: {gisaid_tsv} ...")
    gisaid_map = {}
    with open(gisaid_tsv, 'r', encoding='utf-8') as f_in:
        for line in f_in:
            parts = line.rstrip('\r\n').split('\t')
            if len(parts) >= 3:
                sha = parts[0]
                gisaid_map[sha] = parts[2:]

    start_t = time.time()
    n_found = 0
    n_missing = 0

    print(f"Resolving child array subset {subset_txt} ...")
    with open(subset_txt, 'r', encoding='utf-8') as f_in, open(out_txt, 'w', encoding='utf-8') as f_out:
        for line in f_in:
            if not line.strip():
                continue
            parts = line.strip().split('\t')

            # Explicitly catch TSV text headers
            if parts[0] == 'output_fasta' or (len(parts) > 1 and parts[1] == 'sha256'):
                continue

            # Support Phase 4 extraction TSVs correctly
            if len(parts) >= 3 and '.fasta' in parts[0]:
                sha = parts[1]
            else:
                tok = parts[0]
                # Strip NNNNx. prefix safely if it visually exists
                if tok.find('x.') >= 0:
                    sha = tok.split('x.')[1]
                elif tok.find('x') >= 0:
                    sha = tok.split('x')[1]
                else:
                    sha = tok

            targets = gisaid_map.get(sha, None)
            if targets is not None:
                for tgt in targets:
                    f_out.write(f"{tgt}\n")
                    n_found += 1
            else:
                n_missing += 1

    dur = time.time() - start_t
    print(f"Extracted {n_found:,} original GISAID mapped identifiers cleanly to {out_txt} in {dur:.2f}s.")
    if n_missing > 0:
        print(f"Warning: {n_missing:,} sequences within the subset array did NOT match any identical SHA256 "
              "string entries generically in the core mapping.", file=sys.stderr)


def crosscheck_fasta_vs_txt(fasta_path: str, generated_txt: str):
    """Mode 3: Parses biological .fasta dynamically cross-checking strictly against standard txt."""
    if not os.path.exists(fasta_path) or not os.path.exists(generated_txt):
        print("Error: Missing explicit crosscheck subsets.", file=sys.stderr)
        return

    gisaid_re = re.compile(r'\|(EPI_ISL_[^|]+)\|')
    fasta_ids = set()
    print(f"Scanning physical biological FASTA {os.path.basename(fasta_path)}...")
    with open(fasta_path, 'r', encoding='utf-8') as fh:
        for line in fh:
            if line.startswith('>'):
                match = gisaid_re.search(line)
                if match:
                    fasta_ids.add(match.group(1))

    txt_ids = set()
    print(f"Scanning generated text map {os.path.basename(generated_txt)}...")
    with open(generated_txt, 'r', encoding='utf-8') as fh:
        for line in fh:
            if line.strip():
                txt_ids.add(line.strip())

    d1 = fasta_ids - txt_ids
    d2 = txt_ids - fasta_ids
    if not d1 and not d2:
        print("100% Sequence Cross-Check Clean! Sets are perfectly mathematically identical computationally.")
    else:
        print("Validation Warning! Mismatched bounds dynamically mapped:")
        if d1:
            print(f"  {len(d1)} elements uniquely in FASTA (e.g. {list(d1)[0]}).")
        if d2:
            print(f"  {len(d2)} elements uniquely in Text (e.g. {list(d2)[0]}).")


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Standalone Extractor: Map and stream GISAID IDs natively from pipeline TSVs.")

    # Mode 1: Initial mapping generator
    parser.add_argument('--descr-tsv', help="Source TSV to parse (e.g. spikenuc1207.sha256_to_descr_lines.tsv)")
    parser.add_argument('--out-gisaid-tsv', help="Output mapping TSV (e.g. spikenuc1207.sha256_to_GISAID_ids.tsv)")

    # Mode 2: Sequential subset textual extractor
    parser.add_argument('--gisaid-tsv-mapping', help="Existing GISAID TSV generated from Mode 1.")
    parser.add_argument(
        '--subset-txt', help="Source subset hash text (e.g. spikenuc1207.no_junk.counts...extraction_counts.tsv)")
    parser.add_argument('--out-subset-txt', help="Output flat generic identifier list text file")

    # Mode 3: Validation Crosschecker
    parser.add_argument('--crosscheck', action='store_true', help="Enable orthornormal validation.")
    parser.add_argument('--fasta', help="Physical extracted FASTA array.")

    args = parser.parse_args()

    if args.crosscheck and args.fasta and args.subset_txt:
        crosscheck_fasta_vs_txt(args.fasta, args.subset_txt)
    elif args.descr_tsv and args.out_gisaid_tsv:
        build_gisaid_tsv(args.descr_tsv, args.out_gisaid_tsv)
    elif args.gisaid_tsv_mapping and args.subset_txt and args.out_subset_txt:
        extract_subset_txt(args.gisaid_tsv_mapping, args.subset_txt, args.out_subset_txt)

    else:
        parser.print_help()
        sys.exit(1)
