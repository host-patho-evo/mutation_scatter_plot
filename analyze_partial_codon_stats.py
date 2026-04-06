#!/usr/bin/env python3
import sys
import collections
import argparse
import os

def extract_stats(fasta_file, aa_start):
    if not os.path.exists(fasta_file):
        print(f"Error: Could not locate {fasta_file}")
        sys.exit(1)

    frame_start = aa_start - 1
    partial_codons = collections.Counter()
    total_internal_anomalies = 0

    with open(fasta_file, 'r') as f_in:
        seq_parts = []
        for line in f_in:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if seq_parts:
                    seq = "".join(seq_parts)
                    total_internal_anomalies += process_sequence(seq, frame_start, partial_codons)
                    seq_parts = []
            else:
                seq_parts.append(line)
        if seq_parts:
            total_internal_anomalies += process_sequence("".join(seq_parts), frame_start, partial_codons)

    print(f"Scan complete. Found {total_internal_anomalies} total internal mixed codons.")
    print("----- Frequency Table -----")
    for codon, count in partial_codons.most_common():
        print(f"{codon} : {count} times")

def process_sequence(seq, frame_start, partial_codons):
    codons = []
    for pos in range(frame_start, len(seq), 3):
        codon = seq[pos:pos+3]
        if len(codon) == 3:
            codons.append(codon)

    first_base_idx = -1
    last_base_idx = -1
    for i, codon in enumerate(codons):
        if any(c not in '-Nn?' for c in codon):
            if first_base_idx == -1:
                first_base_idx = i
            last_base_idx = i

    anomalies = 0
    if first_base_idx != -1 and first_base_idx + 1 < last_base_idx:
        for i in range(first_base_idx + 1, last_base_idx):
            codon = codons[i]
            has_gap = any(c in '-Nn?' for c in codon)
            has_base = any(c not in '-Nn?' for c in codon)

            if has_gap and has_base:
                partial_codons[codon] += 1
                anomalies += 1
    return anomalies

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Tally statistics of partial codons from the isolated fasta file")
    parser.add_argument("--fasta", required=True, help="Input failing_ragged_sequences.fasta file")
    parser.add_argument("--aa-start", type=int, default=1, help="1-based nucleotide index mapping frame")

    args = parser.parse_args()
    extract_stats(args.fasta, args.aa_start)
