# Profiling and Latency Analysis Report

This report summarizes the profiling results and specific speedups achieved during the optimization phase.

## 1. `calculate_codon_frequencies` Profiling
- **Baseline**: 23.89 seconds for `test2.fasta`.
- **Primary Bottleneck**: $O(N^2)$ complexity in `parse_alignment` due to redundant string operations (regex, cleaning) inside the tightest inner loop.
- **Optimization 1 (Hoisting)**: Static sequence properties (depadding, regex searches) were moved outside the codon evaluation loop.
- **Optimization 2 (Caching)**: Applied `functools.lru_cache` to `alt_translate` to eliminate recursive BioPython translation overhead.
- **Final Result**: **0.60 seconds** (approx. 40x total speedup, 119x for `parse_alignment` specifically).

## 2. `mutation_scatter_plot` Profiling
- **Baseline**: 17.83 seconds for typical visual rendering.
- **Primary Bottleneck**: Redundant `df.loc` queries in Pandas for empty matrix sites.
- **Optimization**: Hoisted data extraction logic behind a frequency threshold guard. Dictionary lookups and Pandas slicing now only occur for non-zero data points.
- **Final Result**: **12.26 seconds** (31% total script speedup).
- **Remaining Latency**: Analysis shows remaining time is dominated by library imports (Pandas, Matplotlib) and PNG zlib encoding, with no remaining algorithmic bottlenecks.

## 3. Summary of Algorithmic Complexity
- **Alignment Parsing**: Reduced from $O(N_{seq} \cdot L)$ to $O(N_{unique\_seq\_parts} \cdot L)$.
- **Scatter Plotting**: Reduced from $O(N_{sites} \cdot 64)$ to $O(N_{visible\_dots})$.
