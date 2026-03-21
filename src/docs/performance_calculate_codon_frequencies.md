# Performance Optimizations for calculate_codon_frequencies

This document details the performance optimizations implemented to handle large FASTA alignments efficiently.

## 1. Streaming & Global Grouping
- **Streaming Input**: Replaced `AlignIO.read` (which loads the entire file into memory) with `SeqIO.parse`.
- **Memory Efficiency**: Memory usage is now proportional to the number of *unique* sequences in the alignment, rather than the total number of entries. This allows the tool to process files with hundreds of thousands of sequences.
- **Global Grouping**: The tool builds a dictionary of unique sequences and their global occurrence counts during the initial streaming pass.

## 2. Local Codon-Grouping (Site-wise)
- **Algorithmic Reduction**: At each alignment position, sequences are grouped by their specific 3-nucleotide codon and action state (e.g., gaps or insertions).
- **Complexity**: Reduced from $O(N_{sequences})$ to $O(N_{unique\_codons})$ per site.
- **Logic Hoisting**: Complex mutation, insertion, and deletion logic is executed once per unique codon (typically 1–64) instead of once per sequence.

## 3. Implementation Details
- **Sequence Caching**: Sequences are pre-processed and stored as efficient string objects.
- **Parsing Loop**: The inner loop now iterates over unique codon groups, significantly reducing Python overhead.
- **Verification**: All optimizations were verified for byte-for-byte identity against original baseline outputs using the integration test suite.
