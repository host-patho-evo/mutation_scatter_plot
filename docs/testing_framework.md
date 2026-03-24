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
- **Non-Interactive Flags**: The `--disable-showing-bokeh` and `--disable-showing-mplcursors` flags are utilized during testing to prevent the browser or interactive cursors from obstructing the automated flow.

## 3. Geometric Scaling Tests (`test_scaling_logic.py`)
- **Area-to-Frequency Verification**: Programmatically audits the Matplotlib and Bokeh internal attributes (circle sizes) to ensure they follow the 1:10:100:300 area-scaling law.
- **Precision Auditing**: Compares computed properties (e.g., `_offsets` and `_sizes` in Matplotlib) against theoretical ratios to maintain 95%+ geometric fidelity.

## 4. Regression Suite & Golden Refresh
- **Dataset Coverage**: The suite manages **16 distinct regression datasets** (38+ artifacts including PNG, PDF, HTML, and TSV).
- **`REGENERATE_GOLDENS=1`**: Setting this environment variable before running `pytest` will overwrite the existing golden baseline files with the current output. This is the standard procedure for synchronizing artifacts after intentional rendering logic changes.
    ```bash
    REGENERATE_GOLDENS=1 pytest tests/test_mutation_scatter_plot.py
    ```

## 5. Test Coverage
- **CLI Options**: Covers various combinations of `--x-after-count`, `--min_start`, and `--aminoacids`.

## 6. Execution Guide: Copy & Paste Commands

### Running All Tests
```bash
PYTHONPATH=$(pwd)/src MPLBACKEND=Agg pytest tests/
```

### Running Specific Test Modules
- **Mutation Scatter Plot (Regression)**:
  ```bash
  PYTHONPATH=$(pwd)/src MPLBACKEND=Agg pytest tests/test_mutation_scatter_plot.py
  ```
- **Codon Frequencies**:
  ```bash
  PYTHONPATH=$(pwd)/src pytest tests/test_calculate_codon_frequencies.py
  ```
- **Geometric Scaling**:
  ```bash
  PYTHONPATH=$(pwd)/src MPLBACKEND=Agg pytest tests/test_scaling_logic.py
  ```
- **Render Utility Functions**:
  ```bash
  PYTHONPATH=$(pwd)/src pytest tests/test_render_functions.py
  ```

### Code Coverage
```bash
PYTHONPATH=$(pwd)/src pytest --cov=mutation_scatter_plot tests/
```

### Static Analysis (Pylint)
To verify the 10.0/10.0 score across the codebase:
```bash
pylint $(git ls-files '*.py')
```

### Synchronizing Golden Artifacts
If you intentionally changed the rendering logic, regenerate the baselines:
```bash
REGENERATE_GOLDENS=1 PYTHONPATH=$(pwd)/src MPLBACKEND=Agg pytest tests/test_mutation_scatter_plot.py
```
