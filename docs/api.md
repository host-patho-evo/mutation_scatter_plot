# Programmatic API Reference

`mutation_scatter_plot` exposes a high-level Python API for rendering mutation
scatter plots, timeline plots, and computing codon frequencies — without using
the command-line interface.

All functions are importable directly from the top-level package:

```python
from mutation_scatter_plot import (
    render_scatter,
    render_timeline,
    calculate_frequencies,
    scatter_options,
    timeline_options,
    frequency_options,
    alt_translate,
)
```

> [!NOTE]
> These imports are **lazy**: importing `mutation_scatter_plot` loads only
> `Bio.Seq` (~300 ms).  Heavy dependencies (matplotlib, pandas, numpy) are
> deferred until you actually call a function.

---

## Table of Contents

- [render\_scatter()](#render_scatter)
- [render\_timeline()](#render_timeline)
- [calculate\_frequencies()](#calculate_frequencies)
- [Options Factories](#options-factories)
- [Low-Level Access](#low-level-access)
- [Usage Examples](#usage-examples)
  - [Example 1: Jupyter Notebook — Interactive Mutation Explorer](#example-1-jupyter-notebook--interactive-mutation-explorer)
  - [Example 2: Flask Web App — On-Demand Plot Generation](#example-2-flask-web-app--on-demand-plot-generation)
  - [Example 3: Automated Surveillance Pipeline](#example-3-automated-surveillance-pipeline)
  - [Example 4: BLOSUM Scoring Matrix Comparison](#example-4-blosum-scoring-matrix-comparison)
  - [Example 5: Codon Usage Bias Analyser](#example-5-codon-usage-bias-analyser)

---

## `render_scatter()`

Render a mutation scatter plot from a `.frequencies.tsv` file.

```python
render_scatter(
    tsv_path: str,
    outfile_prefix: str = '',
    *,
    aminoacids: bool = False,
    xmin: int = 0,
    xmax: int = 0,
    colormap: str = 'amino_acid_changes',
    matrix: str = 'BLOSUM80',
    linear_scaling: bool = False,
    threshold: float = 0.001,
    include_synonymous: bool = False,
    show_bokeh: bool = False,
    show_mplcursors: bool = False,
    dpi: int = 600,
    title: str = '',
    **extra_options,
) -> dict
```

### Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `tsv_path` | `str` | *required* | Path to a `.frequencies.tsv` file |
| `outfile_prefix` | `str` | `''` | Output path prefix; when empty, no files are saved |
| `aminoacids` | `bool` | `False` | Amino acid mode (`True`) or codon mode (`False`) |
| `xmin`, `xmax` | `int` | `0` | X-axis range (padded positions); 0 = auto-detect |
| `colormap` | `str` | `'amino_acid_changes'` | Colormap name |
| `matrix` | `str` | `'BLOSUM80'` | BLOSUM substitution matrix |
| `linear_scaling` | `bool` | `False` | Linear circle-size scaling instead of area |
| `threshold` | `float` | `0.001` | Minimum frequency to display |
| `include_synonymous` | `bool` | `False` | Include synonymous mutations in amino acid mode |
| `show_bokeh` | `bool` | `False` | Open Bokeh HTML in browser |
| `show_mplcursors` | `bool` | `False` | Open interactive matplotlib window |
| `dpi` | `int` | `600` | Raster output resolution |
| `title` | `str` | `''` | Figure title (auto-derived if empty) |

### Returns

A `dict` with keys:

| Key | Type | Description |
|-----|------|-------------|
| `'figure'` | `matplotlib.figure.Figure` | The rendered figure object |
| `'axes'` | `tuple` | `(ax1, ax2, ax3, ax4)` — scatter, bar, colorbar, empty |
| `'dataframe'` | `pandas.DataFrame` | The cleaned frequency data |
| `'files_written'` | `list[str]` | Paths of saved files (empty if no `outfile_prefix`) |

### Basic Usage

```python
from mutation_scatter_plot import render_scatter

# Get just the Figure — no files written
result = render_scatter(
    tsv_path='spike.frequencies.tsv',
    aminoacids=True,
    xmin=430, xmax=528,
    threshold=0.01,
)
fig = result['figure']
fig.savefig('custom_output.png', dpi=300)

# Save all output files (PNG, PDF, HTML, colors TSV)
result = render_scatter(
    tsv_path='spike.frequencies.tsv',
    outfile_prefix='output/rbd_mutations',
    aminoacids=True,
    xmin=430, xmax=528,
)
print(result['files_written'])
# ['output/rbd_mutations.BLOSUM80.area_scaling.amino_acid_changes.pdf',
#  'output/rbd_mutations.BLOSUM80.area_scaling.amino_acid_changes.png', ...]
```

---

## `render_timeline()`

Render a timeline scatter plot from per-month frequency TSV files.

```python
render_timeline(
    directory: str,
    positions: list[str],
    outfile_prefix: str = '',
    *,
    aminoacids: bool = False,
    threshold: float = 0.0,
    colormap: str = 'coolwarm_r',
    matrix: str = 'BLOSUM80',
    show_bokeh: bool = False,
    dpi: int = 600,
    **extra_options,
) -> dict
```

### Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `directory` | `str` | *required* | Directory with per-month `.frequencies.tsv` files |
| `positions` | `list[str]` | *required* | Position specs, e.g. `['N501Y', '498[RHQ]', '484']` |
| `outfile_prefix` | `str` | `''` | Output prefix; auto-inferred from filenames if empty |
| `aminoacids` | `bool` | `False` | Amino acid mode (`True`) or codon mode |
| `threshold` | `float` | `0.0` | Minimum frequency to include |
| `colormap` | `str` | `'coolwarm_r'` | Colormap name |
| `matrix` | `str` | `'BLOSUM80'` | BLOSUM matrix  |
| `show_bokeh` | `bool` | `False` | Open Bokeh HTML in browser |
| `dpi` | `int` | `600` | Raster output resolution |

### Position Specification Syntax

Positions support several formats:

| Format | Meaning | Example |
|--------|---------|---------|
| Numeric | All mutations at that position | `'614'` |
| Range | All mutations in a range | `'480-530'` |
| Specific mutation | Single mutation | `'D614G'`, `'N501Y'` |
| Bracket notation | Multiple mutants at one position | `'498[RHQ]'` |

### Returns

A `dict` with keys:

| Key | Type | Description |
|-----|------|-------------|
| `'data'` | `TimelineData` | Dataclass with `.points`, `.months`, `.positions` |
| `'files_written'` | `list[str]` | Paths of saved PNG/PDF files |

### Basic Usage

```python
from mutation_scatter_plot import render_timeline

result = render_timeline(
    directory='per_month_results/',
    positions=['N501Y', '498[RHQ]', 'D614G'],
    outfile_prefix='output/timeline',
    aminoacids=True,
)

# Inspect collected data
for point in result['data'].points[:5]:
    print(f"{point.month}: {point.label} freq={point.frequency:.4f}")

print(f"Saved: {result['files_written']}")
```

---

## `calculate_frequencies()`

Calculate codon/amino acid frequencies from a FASTA alignment
and return the results as a pandas DataFrame.

```python
calculate_frequencies(
    alignment_file: str,
    reference_file: str,
    outfile_prefix: str,
    *,
    padded_reference: bool = True,
    x_after_count: bool = False,
    print_unchanged_sites: bool = True,
    threads: int = 0,
    translation_table: int = 1,
    overwrite: bool = True,
    **extra_options,
) -> pandas.DataFrame
```

### Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `alignment_file` | `str` | *required* | Path to padded multi-FASTA alignment |
| `reference_file` | `str` | *required* | Path to reference FASTA (single sequence) |
| `outfile_prefix` | `str` | *required* | Output prefix for `.frequencies.tsv` etc. |
| `padded_reference` | `bool` | `True` | Whether the reference is padded with dashes |
| `x_after_count` | `bool` | `False` | FASTA IDs contain `NNNNx.` count prefixes |
| `print_unchanged_sites` | `bool` | `True` | Also write unchanged codons TSV |
| `threads` | `int` | `0` | Worker processes (0 = auto-detect) |
| `translation_table` | `int` | `1` | NCBI genetic code table |
| `overwrite` | `bool` | `True` | Overwrite existing output files |

### Returns

A `pandas.DataFrame` with columns:

| Column | Type | Description |
|--------|------|-------------|
| `padded_position` | `int` | Alignment-padded codon position |
| `position` | `int` | Original (unpadded) codon position |
| `original_aa` | `str` | Reference amino acid |
| `mutant_aa` | `str` | Observed amino acid |
| `frequency` | `float` | Mutation frequency (0.0–1.0) |
| `original_codon` | `str` | Reference codon triplet |
| `mutant_codon` | `str` | Observed codon triplet |
| `observed_codon_count` | `int` | Number of sequences with this codon |
| `total_codons_per_site` | `int` | Total sequences aligned at this site |

### Basic Usage

```python
from mutation_scatter_plot import calculate_frequencies

df = calculate_frequencies(
    alignment_file='spike_aligned.fasta',
    reference_file='MN908947.3_S.3873.fasta',
    outfile_prefix='output/spike',
    padded_reference=True,
    threads=4,
)

# Filter for high-frequency mutations
high_freq = df[df['frequency'] > 0.05]
print(f"Found {len(high_freq)} high-frequency variants")
print(high_freq[['position', 'original_aa', 'mutant_aa', 'frequency']])
```

---

## Options Factories

For advanced use cases that require direct access to the internal pipeline
functions, factory functions create properly-defaulted `argparse.Namespace`
objects:

```python
from mutation_scatter_plot import scatter_options, timeline_options, frequency_options
```

### `scatter_options(**overrides)`

```python
opts = scatter_options(
    tsv_file_path='data.tsv',
    outfile_prefix='my_plot',
    aminoacids=True,
    xmin=430, xmax=528,
    disable_showing_bokeh=True,
    disable_showing_mplcursors=True,
)
```

### `timeline_options(**overrides)`

```python
opts = timeline_options(
    directory='per_month/',
    positions=['N501Y', '484'],
    outfile_prefix='output/timeline',
)
```

### `frequency_options(**overrides)`

```python
opts = frequency_options(
    alignment_infilename='aligned.fasta',
    reference_infilename='reference.fasta',
    outfileprefix='output/freq',
    padded_reference=True,
    threads=4,
)
```

> [!WARNING]
> Unknown keyword arguments raise a `TypeError` with the list of valid keys.
> This catches typos early.

```python
>>> scatter_options(aminacids=True)   # typo!
TypeError: Unknown option 'aminacids'. Valid options are: ['aminoacids', ...]
```

---

## Low-Level Access

For maximum control, the internal pipeline components are accessible via
fully-qualified imports:

```python
# Substitution matrix loading and scoring
from mutation_scatter_plot.mutation_scatter_plot import load_matrix
from mutation_scatter_plot.mutation_scatter_plot.core import get_score, get_colormap

# DataFrame manipulation
from mutation_scatter_plot.mutation_scatter_plot import (
    load_and_clean_dataframe,
    build_frequency_tables,
    setup_matplotlib_figure,
    collect_scatter_data,
    render_bokeh,
    render_matplotlib,
)

# Timeline components
from mutation_scatter_plot.mutation_scatter_plot.timeline import (
    scan_directory,
    parse_positions,
    collect_timeline_data,
    render_timeline_matplotlib,
    render_timeline_bokeh,
    TimelineData,
    PositionSpec,
)

# Translation utility
from mutation_scatter_plot import alt_translate
alt_translate('ATG')  # → 'M'
alt_translate('TC-')  # → 'S'  (gap treated as N → TCN → Ser)
alt_translate('---')  # → '-'  (full gap → gap amino acid)
```

---

## Usage Examples

### Example 1: Jupyter Notebook — Interactive Mutation Explorer

```python
"""Interactive exploration of SARS-CoV-2 RBD mutations in a Jupyter notebook."""

from mutation_scatter_plot import render_scatter
import matplotlib.pyplot as plt

# Render with inline display, no file/GUI popups
result = render_scatter(
    tsv_path='spikenuc0719.no_junk.frequencies.tsv',
    outfile_prefix='notebook_output/rbd_mutations',
    aminoacids=True,
    xmin=430, xmax=528,
    colormap='amino_acid_changes',
    threshold=0.01,
    show_mplcursors=False,
    show_bokeh=False,
    dpi=150,  # lower DPI for notebook display
    title='SARS-CoV-2 Spike RBD Mutations',
)

fig = result['figure']
fig.set_size_inches(14, 8)
plt.show()

# Access the underlying data for custom analysis
df = result['dataframe']
high_freq = df[df['frequency'] > 0.1]
print(f"Positions with >10% mutation frequency:\n{high_freq}")
```

---

### Example 2: Flask Web App — On-Demand Plot Generation

```python
"""Flask endpoint that generates mutation scatter plots on demand."""

from flask import Flask, request, send_file
from mutation_scatter_plot import render_scatter
import tempfile
import os

app = Flask(__name__)

@app.route('/api/scatter', methods=['POST'])
def generate_scatter():
    """POST a TSV file, get back a PNG scatter plot."""
    tsv_file = request.files['tsv']
    xmin = int(request.form.get('xmin', 0))
    xmax = int(request.form.get('xmax', 0))
    mode = request.form.get('mode', 'aminoacids')  # or 'codons'

    with tempfile.TemporaryDirectory() as tmpdir:
        tsv_path = os.path.join(tmpdir, 'upload.frequencies.tsv')
        tsv_file.save(tsv_path)

        result = render_scatter(
            tsv_path=tsv_path,
            outfile_prefix=os.path.join(tmpdir, 'plot'),
            aminoacids=(mode == 'aminoacids'),
            xmin=xmin, xmax=xmax,
            show_bokeh=False,
            show_mplcursors=False,
            dpi=300,
        )

        png_path = [f for f in result['files_written'] if f.endswith('.png')][0]
        return send_file(png_path, mimetype='image/png')
```

---

### Example 3: Automated Surveillance Pipeline

```python
"""Monthly GISAID surveillance: calculate frequencies, render timeline,
detect emerging mutations, and print alerts."""

from mutation_scatter_plot import calculate_frequencies, render_timeline

# Step 1: Calculate codon frequencies from new month's alignment
df = calculate_frequencies(
    alignment_file='2026-03/spike_aligned.fasta',
    reference_file='reference/MN908947.3_S.3873.fasta',
    outfile_prefix='2026-03/spike',
    padded_reference=True,
    x_after_count=True,
    threads=8,
)

# Step 2: Detect emerging mutations (>5% frequency, non-synonymous)
emerging = df[
    (df['frequency'] > 0.05)
    & (df['original_aa'] != df['mutant_aa'])
]
for _, row in emerging.iterrows():
    codon_change = f"{row['original_codon']}→{row['mutant_codon']}"
    aa_change = f"{row['original_aa']}{row['position']}{row['mutant_aa']}"
    print(f"ALERT: {aa_change} ({codon_change}) at {row['frequency']:.1%}")

# Step 3: Render the full timeline across all months
result = render_timeline(
    directory='.',
    positions=['498[RHQ]', 'N501Y', '484'],
    outfile_prefix='timeline/spike_rbd',
    aminoacids=True,
    threshold=0.001,
)
print(f"Timeline saved to: {result['files_written']}")
```

---

### Example 4: BLOSUM Scoring Matrix Comparison

```python
"""Compare how different substitution matrices score the same mutations."""

from mutation_scatter_plot import scatter_options
from mutation_scatter_plot.mutation_scatter_plot import load_matrix
from mutation_scatter_plot.mutation_scatter_plot.core import get_score

matrices = ['BLOSUM45', 'BLOSUM62', 'BLOSUM80', 'BLOSUM90']
mutation = ('Q', 'R')  # Q498R — key Omicron RBD mutation

for matrix_name in matrices:
    opts = scatter_options(matrix=matrix_name, outfile_prefix='/tmp/dummy')
    matrix, _name, _min_s, _max_s, _prefix = load_matrix(opts)
    score = get_score(opts, matrix, False, mutation[0], mutation[1])
    print(f"{matrix_name}: Q→R score = {score}")
```

Output:
```
BLOSUM45: Q→R score = 1
BLOSUM62: Q→R score = 1
BLOSUM80: Q→R score = 1
BLOSUM90: Q→R score = 1
```

---

### Example 5: Codon Usage Bias Analyser

```python
"""Analyse codon usage bias per genomic region using the frequency pipeline."""

from mutation_scatter_plot import calculate_frequencies

regions = {
    'NTD': (1, 305),
    'RBD': (319, 541),
    'S2':  (686, 1273),
}

for region_name, (start, end) in regions.items():
    df = calculate_frequencies(
        alignment_file='spike_aligned.fasta',
        reference_file='MN908947.3_S.3873.fasta',
        outfile_prefix=f'codon_bias/{region_name}',
        padded_reference=True,
    )

    # Filter to region
    region_df = df[
        (df['position'] >= start) & (df['position'] <= end)
    ]

    # Calculate codon diversity per site
    site_diversity = (
        region_df.groupby('position')
        .agg(n_codons=('mutant_codon', 'nunique'),
             max_freq=('frequency', 'max'))
        .reset_index()
    )

    conserved = (site_diversity['max_freq'] > 0.99).sum()
    variable = (site_diversity['n_codons'] > 3).sum()
    print(f"{region_name} ({start}–{end}): "
          f"{conserved} conserved sites, {variable} highly variable sites")
```

---

## Complete Pipeline: From FASTA to Figure

A typical end-to-end workflow combining `calculate_frequencies`
and `render_scatter`:

```python
from mutation_scatter_plot import calculate_frequencies, render_scatter

# Step 1: Compute frequencies from a raw alignment
df = calculate_frequencies(
    alignment_file='aligned_sequences.fasta',
    reference_file='reference.fasta',
    outfile_prefix='output/my_analysis',
    padded_reference=True,
    threads=4,
)
print(f"Computed {len(df)} frequency entries")

# Step 2: Render the scatter plot from the TSV that was just written
result = render_scatter(
    tsv_path='output/my_analysis.frequencies.tsv',
    outfile_prefix='output/my_analysis',
    aminoacids=True,
    xmin=319, xmax=541,     # RBD domain
    threshold=0.005,
    dpi=300,
    title='Spike RBD Mutation Frequencies',
)

print(f"Figure saved to: {result['files_written']}")

# Step 3: Custom analysis on the DataFrame
rbd_hotspots = df[
    (df['position'].between(319, 541))
    & (df['frequency'] > 0.1)
    & (df['original_aa'] != df['mutant_aa'])
].sort_values('frequency', ascending=False)

print(f"\nTop RBD mutations (>10% frequency):")
for _, row in rbd_hotspots.head(10).iterrows():
    print(f"  {row['original_aa']}{row['position']}{row['mutant_aa']}: "
          f"{row['frequency']:.1%}")
```
