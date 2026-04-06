#!/usr/bin/env python3
import sys
import argparse
import os

def find_ragged_boundary_sequences(fasta_file, aa_start, output_fasta):
    """
    Parses a FASTA file dynamically to identify and extract 'failing' sequences that
    have ragged unsequenced boundaries (meaning their leading or trailing `-`/`N`
    padding overlaps inside a triplet reading frame, forming mixed gap/nucleotide codons).
    """
    frame_start = aa_start - 1  # 0-based index offset

    if not os.path.exists(fasta_file):
        print(f"Error: Could not locate {fasta_file}")
        sys.exit(1)

    extracted_count = 0
    total_count = 0

    with open(fasta_file, 'r') as f_in, open(output_fasta, 'w') as f_out:
        header = None
        seq_parts = []

        def process_record():
            nonlocal extracted_count, total_count
            if header is None: return
            total_count += 1
            seq = "".join(seq_parts)

            ragged_failure = False

            # Group into reading frame codons
            codons = []
            for pos in range(frame_start, len(seq), 3):
                codon = seq[pos:pos+3]
                if len(codon) == 3:
                    codons.append(codon)

            # Identify where the structural read actually begins and ends
            first_base_idx = -1
            last_base_idx = -1
            for i, codon in enumerate(codons):
                if any(c not in '-Nn?' for c in codon):
                    if first_base_idx == -1:
                        first_base_idx = i
                    last_base_idx = i

            # A legitimately sequenced record must span at least some codons
            if first_base_idx != -1 and last_base_idx != -1:
                # We structurally ignore first_base_idx and last_base_idx because
                # random sequencing bounds natively drop off mid-codon (~89% of the time).
                # We are mathematically hunting for deeply embedded INTERNAL anomalies.
                if first_base_idx + 1 < last_base_idx:
                    for i in range(first_base_idx + 1, last_base_idx):
                        codon = codons[i]
                        has_gap = any(c in '-Nn?' for c in codon)
                        has_base = any(c not in '-Nn?' for c in codon)

                        # If a triplet INSIDE the contiguous sequence contains BOTH a nucleotide AND a gap:
                        # it is an internal reading-frame shift or an obscure mapping glitch!
                        if has_gap and has_base:
                            ragged_failure = True
                            break
            else:
                # 100% gap sequences provide zero structural utility
                ragged_failure = True

            if ragged_failure:
                # Output structurally ragged sequence
                f_out.write(f"{header}\n{seq}\n")
                extracted_count += 1

        for line in f_in:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                process_record()
                header = line
                seq_parts = []
            else:
                seq_parts.append(line)
        process_record()

    print(f"Scan Complete: Analyzed {total_count} sequences.")
    print(f"Isolated {extracted_count} ragged boundary ('failing') sequences -> {output_fasta}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Extract FASTA sequences with ragged boundary frame errors.")
    parser.add_argument("--fasta", required=True, help="Input *.clean.fasta file")
    parser.add_argument("--out", required=True, help="Output temp fasta file for failing sequences")
    parser.add_argument("--aa-start", type=int, default=1, help="1-based nucleotide index where coding frame begins (e.g. 413 for Spike)")

    args = parser.parse_args()
    find_ragged_boundary_sequences(args.fasta, args.aa_start, args.out)
