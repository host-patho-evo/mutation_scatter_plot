# Integration Testing Framework

The repository features a robust integration testing framework designed to ensure functional correctness and prevent regression across both calculation and visualization components.

## 1. Calculation Tests (`test_calculate_codon_frequencies.py`)
- **Automated Regression**: Uses `unittest` to execute the `calculate_codon_frequencies` binary across multiple parameter sets.
- **Deterministic Validation**: Compares generated TSV outputs against "golden baselines" stored in `tests/outputs/` using byte-for-byte file comparison (`filecmp.cmp`).
- **Temporary Execution**: Utilizes `tempfile.TemporaryDirectory` to ensure a clean environment and safe cleanup of artifacts between test cases.

## 2. Visualization Tests (`test_mutation_scatter_plot.py`)
- **rendering Robustness**: Evaluates Matplotlib (PNG/PDF) and Bokeh (HTML) generation.
- **Deterministic Bokeh Validation**: Since Bokeh HTML contains non-deterministic UUIDs and timestamps, the test suite includes a structural JSON extractor. It pulls the raw `ColumnDataSource` payload from the HTML to verify the plotted data points match the expected values exactly.
- **Headless Execution**: Forcing `MPLBACKEND=Agg` ensures tests run reliably in CI environments without attempting to launch graphical windows.

## 3. Test Coverage
- **Edge Cases**: Includes tests for synaptic mutations, stop codons, insertions, and deletions.
- **CLI Options**: Covers various combinations of `--x-after-count`, `--min_start`, and `--aminoacids`.
