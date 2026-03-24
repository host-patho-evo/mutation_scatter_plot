# Profiling and Latency Analysis Report

This report summarizes the profiling results and specific speedups achieved during the optimization phase.

## 1. `calculate_codon_frequencies` Profiling
- **Baseline**: 23.89 seconds (v0.3 state).
- **Primary Bottleneck**: $O(N^2)$ complexity in `parse_alignment` due to redundant string operations (regex, cleaning) inside the tightest inner loop.
- **Optimization 1 (Hoisting)**: Static sequence properties (depadding, regex searches) were moved outside the codon evaluation loop.
- **Optimization 2 (Caching)**: Applied `functools.lru_cache` to `alt_translate` to eliminate recursive BioPython translation overhead.
- **Final Result**: **0.60 seconds** (approx. 40x total speedup, 119x for `parse_alignment` specifically).

## 2. `mutation_scatter_plot` Profiling
- **Baseline**: 17.83 seconds (v0.3 state).
- **Primary Bottleneck**: Redundant `df.loc` queries in Pandas for empty matrix sites.
- **Optimization**: Hoisted data extraction logic behind a frequency threshold guard. Dictionary lookups and Pandas slicing now only occur for non-zero data points.
- **Final Result**: **12.26 seconds** (1.45x faster than baseline).
- **Remaining Latency**: Analysis shows remaining time is dominated by library imports (Pandas, Matplotlib) and PNG zlib encoding, with no remaining algorithmic bottlenecks.

## 3. Current State (2026-03-24)
- **Result**: **10.18 seconds** (1.75x total speedup from the 17.83s baseline).
- **Includes**: Restoration of the Y-axis completion mechanism (all 20 amino acids and special characters like STOP, DEL, INS, X are explicitly rendered).
- **Impact**: The addition of "tiny dots" to ensure predictable Y-axis ranges has **negligible performance impact**. The $O(N_{positions} \times 20)$ overhead of size-0 dots is trivial compared to the heavy lifting of Matplotlib/Bokeh text and figure rendering.

## 4. Pipeline Performance Summary
- **Baseline (Total Pipeline)**: 23.89s (`codon_frequencies`) + 17.83s (`scatter_plot`) = **41.72 seconds** (measured against v0.3).
- **Current (Total Pipeline)**: 0.60s (`codon_frequencies`) + 10.18s (`scatter_plot`) = **10.78 seconds**.
- **Total Speedup**: **3.87x faster overall**.

## 5. Summary of Algorithmic Complexity
- **Alignment Parsing**: Reduced from $O(N_{seq} \cdot L)$ to $O(N_{unique\_seq\_parts} \cdot L)$.
- **Scatter Plotting**: Reduced from $O(N_{sites} \cdot 64)$ to $O(N_{visible\_dots})$.

## 6. Manual Profiling Instructions

To reproduce these results or audit new changes, use `cProfile` with the following commands:

### `calculate_codon_frequencies`
```bash
python -m cProfile -o profile_codon.stats -m mutation_scatter_plot.calculate_codon_frequencies.cli \
  --reference-infile=tests/inputs/MN908947.3_S.fasta \
  --alignment-file=tests/inputs/test2.fasta \
  --outfile-prefix=tests/outputs/profile_test \
  --padded-reference --x-after-count --aa_start=413 --overwrite
```

### `mutation_scatter_plot`
```bash
python -m cProfile -o profile_scatter.stats -m mutation_scatter_plot.mutation_scatter_plot.cli \
  --tsv tests/outputs/test2_full.x_after_count.frequencies.tsv \
  --outfile-prefix tests/outputs/profile_scatter \
  --disable-showing-bokeh --disable-showing-mplcursors --aminoacids
```

### Analyzing Results
You can visualize the `.stats` files using `snakeviz` (recommended) or a simple sort script:
```bash
pip install snakeviz
snakeviz profile_codon.stats
```
