## Calculate in each codon position frequency of codons and respective amino acids from a multiple sequence alignment and draw an interactive scatter plot in GUI and as HTML with Javascript annotation and static images in JPG, PNG and PDF

The software code and data contained in this folder were used during the work published in **In Vitro and Viral Evolution Convergence Reveal the Selective Pressures Driving Omicron Emergence** publication by Shoshany et al. (submitted, see [bioRxiv](https://www.biorxiv.org/content/10.1101/2025.04.23.650148v1)). Original input data, calculated frequencies in TSV files and also the ready-made figures in JPG/PNG/PDF/HTML+Javascript can be found at [DOI:10.5281/zenodo.18601695](https://zenodo.org/records/18601695) meanwhile. One does not need to install these two utilities to study the results. However, we provide our code to facilitate similar studies of other datasets.

We developed two standalone programs:
`calculate_codon_frequencies` takes a multiple sequence alignment file in a FASTA format (possibly padded with `-` as gaps) and calculates frequencies of the codons. It parses a reference nucleotide sequence of the respective protein to stay in the reading frame (only reading frame +1 is supported). Therefore it is best to map sequencing reads (notably including amplification primer sequences to anchor the alignment ends perfectly) to a complete open reading frame (ORF) sequence encoding the protein (incl. START and STOP codons) although only part of it may have been studied. When this is followed the nucleotide or amino acid positions can be easily calculated from the padding with dashes (`-`) in a multi-FASTA 2-line file (input multiple-sequence alignment file). Otherwise the program allows to specify an arbitrary offset (to be added to the output codon position values) to output native coordinates despite short sequence provided. One can also specify some narrow regions (multiples of three) so that not all columns are to be inspected and codon frequencies calculated.

`mutation_scatter_plot` when run on a graphical terminal opens an interactive matplotlib window with the figures. One can use computer mouse to `hover()` the bubbles to see detailed annotation. After closing the interactive window it also outputs files in JPG, PNG and PDF format and finally also JSON and an interactive HTML+javascript output files are created. One can disable the interactive matplotlib window being raised by adjusting the `MPLBACKEND` environment variable (read further below).

Interactive figures can be visualized and individual frequencies of mutations inspected using python with [matplotlib](https://matplotlib.org/) and [mplcursors](https://pypi.org/project/mplcursors/) installed and concurrently, interactive HTML files are rendered thanks to [Bokeh](https://bokeh.org/).


## System requirements

Any `python-3.10` or newer will work fine. We use python-3.11 and 3.12.

We developed the software using the following versions:
```
>=python-3.10
>=biopython-1.81
>=matplotlib-3.8.4
>=mplcursors-0.5.3
>=pandas-2.2.2
>=numpy-2.2.5
>=blosum-2.0.3
>=bokeh-3.7.2
```

**Installation**

The package is installed with a single command that handles all dependencies automatically. Follow either A.1 (pip) or A.2 (conda) route.

After installation the following commands are available on your `PATH`: `mutation_scatter_plot`, `calculate_codon_frequencies`, `count_motifs_in_sequences`, `alignment2dots`, and `split_fasta_entries_by_lengths`.

The top-level `scripts/` directory contains shell helper scripts (`render-figures.sh`, `per_residue_frequencies.sh`, `create_per_residue_tables.sh`, `split_fasta_entries_by_lengths.sh`) that are **not** installed by `pip install .` — they are provided as reference workflows and must be run directly from the `scripts/` directory or copied manually to a location on your `PATH`.

**A.1 Install using pip into a virtual environment**

```bash
git clone https://github.com/host-patho-evo/mutation_scatter_plot.git
cd mutation_scatter_plot
python3 -m venv .venv
source .venv/bin/activate
pip install .
```

Later you can re-enter this virtual environment by `source .venv/bin/activate`. To leave it call `deactivate`.

To yield interactive matplotlib figures one of the following backends must also be installed:

```bash
pip install wxpython
pip install pyqt5
pip install pyqt6
pip install pycairo
pip install cairocffi
```

**A.2 Install using conda**

```bash
git clone https://github.com/host-patho-evo/mutation_scatter_plot.git
cd mutation_scatter_plot
conda create -n mutation_scatter_plot python=3.13
conda activate mutation_scatter_plot
pip install .
```

To yield interactive figures, one of the matplotlib backends listed above must also be available. For the wxagg backend `wxpython` is needed:

```bash
conda install wxpython
```

You can leave the conda environment by `conda deactivate` and re-enter it later with `conda activate mutation_scatter_plot`.


## Usage

After installation via `pip install .` both commands are available directly on your `PATH`:

The mplcursors `hover()` works under the following backends. Pick just a single command from the following:

```
export MPLBACKEND=gtk3cairo
export MPLBACKEND=wxagg
export MPLBACKEND=qtagg
export MPLBACKEND=gtk4agg
export MPLBACKEND=tkagg # is rather slow
```

The following starts a tornado webserver on a localhost but you do not need it because the script renders interactive HTML with with javascript inside:
```
export MPLBACKEND=webagg
```


To disable interactive matplotlib figures raised on your X11/Wayland display (interactive HTML with with javascript inside are still rendered) use:
```
export MPLBACKEND=agg
```

Below we explain what happens further upon the execution. Not only there is helpful information reported on STDOUT but also a web browser is opened (or a TAB in an existing browser is opened with the HTML+Javascript interactivity) but also an interactive Matplotlib window is opened. Once the Matplotliob window is closed the `mutation_scatter_plot` can finish. In the next paragraphs we describe some of these steps in more detail.

**Simple Testcase to calculate codon frequencies and to render the figures**

One can calculate the codon frequencies from a provided FASTA input file. The sequence must be in frame +1. In version 0.1 the padding dashes `-` are ignored and the software keeps fetches next() nucleotides as long until there are three nucleotides representing the codon to be translated. In the `main` branch there is now a version which requires the inut sequence to be fully aligned to codons. 

```bash
calculate_codon_frequencies --reference-infile=tests/inputs/MN908947.3_S_full.fasta --alignment-file=tests/inputs/test2_full.fasta \
    --outfile-prefix=tests/outputs/test1.default --padded-reference --print-unchanged-sites

prefix='tests/outputs/test1.default'
mutation_scatter_plot --xmin 340 --xmax 516 --tsv "${prefix}.tsv" --outfile-prefix "${prefix}.aa" --aminoacids
mutation_scatter_plot --xmin 340 --xmax 516 --tsv "${prefix}.tsv" --outfile-prefix "${prefix}.codon"
```

**More complex testing**

We have several native alignment scenarios tracked locally in `tests/inputs/`. The following commands illustrate how testing edge case variations are parsed, specifically involving missing reference frames or padded alignments triggering `--aa_start` and `--min_start` coordinate shifts:

```bash
# test2.fasta evaluation:
calculate_codon_frequencies --reference-infile=tests/inputs/MN908947.3_S.fasta --alignment-file=tests/inputs/test2.fasta --outfile-prefix=tests/outputs/test2.x_after_count --padded-reference --x-after-count --aa_start=413

# test2.fasta (with an internal alignment shift via min_start):
calculate_codon_frequencies --reference-infile=tests/inputs/MN908947.3_S.fasta --alignment-file=tests/inputs/test2.fasta --outfile-prefix=tests/outputs/test2.x_after_count_and_min_start --padded-reference --x-after-count --min_start=7 --aa_start=413

# test3.fasta (overcoming a fake early codon sequence):
calculate_codon_frequencies --reference-infile=tests/inputs/MN908947.3_S.fasta --alignment-file=tests/inputs/test3.fasta --outfile-prefix=tests/outputs/test3.default --padded-reference --x-after-count --min_start=4 --aa_start=413
```

You can view the resulting output coordinate lists populated cleanly inside the `tests/outputs/` directory.


**More realistic usage example**

Although we provide already the input, intermediate and resulting files in their respective ZIP bundles for download, to repeat the work or process other data one can take the following procedure to re-create our results. Download real data from [Zenodo https://doi.org/10.5281/zenodo.18601695](https://zenodo.org/records/18601695/files/per_sample_observed_codon_frequencies.zip?download=1). Unpack the ZIP file and pick any from the TSV files, for example `data/intermediates/BA2-4th-round-of-sort__G6.BA2.WTref.frequencies.tsv`.

```
```bash
curl -o per_sample_observed_codon_frequencies.zip "https://zenodo.org/records/18601695/files/per_sample_observed_codon_frequencies.zip?download=1"
unzip per_sample_observed_codon_frequencies.zip
prefix='BA2-4th-round-of-sort__G6.BA2.WTref'
mkdir -p data/outputs/aa/
mkdir -p data/outputs/codon/
mutation_scatter_plot --xmin 430 --xmax 528 --tsv data/intermediates/"$prefix".frequencies.tsv --outfile-prefix data/outputs/aa/"$prefix".aa --aminoacids --disable-showing-bokeh --disable-showing-mplcursors
mutation_scatter_plot --xmin 430 --xmax 528 --tsv data/intermediates/"$prefix".frequencies.tsv --outfile-prefix data/outputs/codon/"$prefix".codon --disable-showing-bokeh --disable-showing-mplcursors
```

**Automated Regression Testing**

The test suite includes a comprehensive regression gallery. When intentional rendering changes are made, the "golden" baselines can be synchronized using:

```bash
REGENERATE_GOLDENS=1 pytest tests/test_mutation_scatter_plot.py
```

For more details on the testing infrastructure, see [src/docs/testing_framework.md](src/docs/testing_framework.md).
```

We also provide a utility to count motifs in [per_sample_unique_sequences_in_FASTA.zip](https://zenodo.org/records/18601695/files/per_sample_unique_sequences_in_FASTA.zip?download=1)
```
curl -o per_sample_unique_sequences_in_FASTA.zip "https://zenodo.org/records/18601695/files/per_sample_unique_sequences_in_FASTA.zip?download=1"
unzip per_sample_unique_sequences_in_FASTA.zip
count_motifs_in_sequences --infilename=data/intermediates/"$prefix".scores_above_84.fastp.amplicons.clean.prot.counts.fasta --motif=RPTY
```

## Sequence Deduplication & Traceability Pipeline

The `scripts/` directory contains a set of standalone helper scripts for
processing large FASTA datasets (e.g. GISAID downloads with millions of
records) through a multi-stage deduplication and traceability pipeline.

### Pipeline scripts

| Script | Purpose |
|---|---|
| `count_same_sequences.py` | Deduplicate FASTA; write `NNNNx.sha256` IDs + sha256→ID mapping TSV |
| `create_list_of_discarded_sequences.py` | Expand a counts FASTA back to original FASTA IDs; identify discarded records |
| `kick.py` | Split a counts FASTA into three files by sequence length (exactly / shorter / longer than N nt) |
| `summarize_fasta_pipeline.py` | Audit the entire pipeline: count records, compute NNNNx sums, show per-step deltas |
| `check_alignment_trimming.py` | Diagnostic: count sequences shortened during alignment (causing sha256 mismatch) |

### File naming convention

Each pipeline stage appends a dot-separated suffix to the output name,
making the parent→child relationship machine-readable:

```
filename_prefix.fasta                          ← raw input
filename_prefix.counts.fasta                   ← after count_same_sequences.py
filename_prefix.counts.clean.fasta             ← after alignment filter
filename_prefix.counts.clean.exactly_N.fasta   ← after kick.py
filename_prefix.counts.clean.exactly_N.discarded_original_ids.txt
```

### Quick start

```bash
# Deduplicate and build sha256 mapping
count_same_sequences.py \
    --infilename=filename_prefix.fasta \
    --outfile-prefix=filename_prefix

# (align externally, producing filename_prefix.counts.clean.fasta)

# Split by target length
kick.py \
    --infile=filename_prefix.counts.clean.fasta \
    --outfile-prefix=filename_prefix.counts.clean \
    --full-length=3822

# List what was discarded at the length-filter step
create_list_of_discarded_sequences.py \
    --infilename=filename_prefix.counts.clean.exactly_3822.fasta \
    --original-infilename=filename_prefix.fasta \
    --inverted

# Audit the whole pipeline at once (uses cached TSV/TXT files when available)
summarize_fasta_pipeline.py . filename_prefix
```

All scripts implement a **make-style timestamp guard**: if all output files
are already newer than all input files the script prints `Info: up-to-date,
skipping` and exits 0.  Add `--overwrite` to force regeneration.

For full documentation see [docs/sequence_deduplication_pipeline.md](docs/sequence_deduplication_pipeline.md).

## Performance & Precision Optimizations

The `calculate_codon_frequencies` and `mutation_scatter_plot` pipelines have been significantly enhanced with the following optimizations:

### Key Improvements
- **Speedup 4: NumPy Vectorized Slicing**: Refactored core alignment parsing to use NumPy arrays and `np.bincount` aggregation, eliminating O(N_sequences) Python loops.
- **Speedup 5: Multi-processing Parallelization**: Decoupled site-specific calculations into a picklable unit, distributed across available CPU cores via `multiprocessing.Pool.starmap`.
- **Speedup 6: High-Precision Decimal Model**: Standardized the entire numerical pipeline on `decimal.Decimal` from the point of ingestion to ensure bit-identical stability across platforms and avoid floating-point drift.

### Optimization Summary

| Speedup | Component | Optimization Type | Result |
| :--- | :--- | :--- | :--- |
| **1** | `build_frequency_tables` | Vectorized `groupby().sum()` | O(1) row aggregation |
| **2** | Matplotlib Hover | O(1) Dictionary Lookup | Instant interactive hover |
| **3** | Bokeh Hover | O(1) Pre-formatted Payload | Eliminated per-dot rendering lag |
| **4** | Codon Slicing | NumPy Vectorization | Rapid alignment column extraction |
| **5** | Codon Processing | Multi-core Parallelization | Scalable performance for large FASTA |
| **6** | Numeric Model | High-Precision `Decimal(str)` | Bit-identical, platform-independent output |
| **7** | **Total Pipeline** | **End-to-End Optimization** | **3.87x overall speedup** (vs v0.3) |

Detailed performance audit results are available in [src/docs/profiling_report.md](src/docs/profiling_report.md).

### Graphical Scaling & Rendering Architecture

The tool supports three circle-size scaling modes that control how mutation
frequencies map to circle sizes across both the Matplotlib and Bokeh backends.
The active mode is encoded in every output filename between the matrix name and
the colormap name (e.g. `.BLOSUM80.area_scaling.amino_acid_changes.pdf`):

| Mode | Flag | Matplotlib `s=` | Bokeh `size=` | Result |
|:---|:---|:---|:---|:---|
| **area_scaling** (default) | *(none)* | `freq × 5000` | `√freq × 100` | area ∝ freq — recommended for bubble charts |
| **linear_scaling** | `--linear-circle-size` | `freq² × 5000` | `freq × 100` | radius ∝ freq — matches v0.2 Bokeh |
| Legacy | `--disable-bokeh-sqrt-size` | `freq × 5000` | `freq × 100` | inconsistent between backends |

The default **area_scaling** mode implements a perceptually linear area-to-frequency
scaling law (1 : 10 : 100 : 300 for frequencies 0.1 % : 1 % : 10 % : 30 %).  The
**linear_scaling** mode makes low-frequency mutations appear relatively larger and
matches the classic v0.2 Bokeh rendering.  Detailed mathematical specifications
and verification results can be found in [docs/rendering_scaling.md](docs/rendering_scaling.md).

## Detailed Technical Documentation

For in-depth analysis of specific modules and subsystems, refer to the following documentation:

- **[Graphical Scaling Architecture](docs/rendering_scaling.md)**: Three scaling modes (area_scaling, linear_scaling, legacy), mathematical derivations, and pixel-level rasterization verification.
- **[Performance Profiling Report](docs/profiling_report.md)**: End-to-end latency analysis and comparison against the v0.3 baseline.
- **[Codon Frequency Optimizations](docs/performance_calculate_codon_frequencies.md)**: Technical breakdown of streaming, grouping, and caching strategies.
- **[Testing Framework Guide](docs/testing_framework.md)**: Overview of the regression suite, golden file management, and automated verification protocols.
- **[Sequence Deduplication Pipeline](docs/sequence_deduplication_pipeline.md)**: Scripts for deduplication, traceability, length filtering, and pipeline auditing of large FASTA datasets.

## Project History

The git commit history has been enriched with the full technical context for major milestones, ensuring that the rationale and implementation details of every optimization (including the 3.87x pipeline speedup and area-proportional scaling laws) are preserved alongside the codebase.


### Parallelization & Environment Variables

The `calculate_codon_frequencies` tool is designed for high-performance multi-core execution. By default, it automatically detects the available CPU allocation from common HPC environment variables, falling back to all available cores if none are set.

**Supported Environment Variables (checked in order):**
1. `PBS_NUM_PPN`
2. `PBS_NCPUS`
3. `OMP_NUM_THREADS`

You can explicitly control the thread count using the `--threads` option:
```bash
# Force use of 8 threads
calculate_codon_frequencies --threads 8 [other options]
```

**Benchmarking Results (4-core system):**
With a high-diversity synthetic dataset (50,000 unique sequences), we measured the following efficiency:| Threads | Time (s) | Speedup | Efficiency |
| :--- | :--- | :--- | :--- |
| 1 | 27.31 | 1.00x | 100% |
| 2 | 28.08 | 0.97x | 49% |
| 4 | 28.32 | 0.96x | 24% |

**Note on Parallelization Efficiency:** The tool is now so heavily optimized that for most datasets, the **sequential grouping phase** (FASTA parsing and record-to-unique grouping) and **NumPy vectorization** dominate the runtime. The parallelized positional processing typically takes less than 1% of the total execution time. Multi-threading is most beneficial for extremely high-diversity alignments with hundreds of thousands of *unique* sequences.
e grouping optimization (Speedup 4) is so effective that it handles the majority of the computational load by collapsing redundant sequences. For most datasets, the processing is now dominated by the sequential FASTA parsing and grouping phase (~26s for 50k sequences), while the parallelized positional analysis takes a fraction of a second (~0.26s). Multi-core processing provides additional scaling headroom for extremely large alignments with extremely high unique sequence diversity.

The `mutation_scatter_plot` tool provides instant interactive hover performance and O(1) data aggregation.


## Example static output images (without interactive features)

The following were rendered by `mutation_scatter_plot`:
![BA2-4th-round-of-sort__G6.BA2.WTref.aa.frequencies.jpg](data/outputs/aa/jpg/BA2-4th-round-of-sort__G6.BA2.WTref.aa.frequencies.jpg)
![BA2-4th-round-of-sort__G6.BA2.WTref.codon.frequencies.jpg](data/outputs/codon/jpg/BA2-4th-round-of-sort__G6.BA2.WTref.codon.frequencies.jpg)

All figures we prepared for our new publication are at [https://host-patho-evo.github.io/mutation_scatter_plot](https://host-patho-evo.github.io/mutation_scatter_plot/).

## Example of a dynamic output of matplotlib window (with interactive features)

Once the interactive window is raised by matplotlib user can point the computer mouse to any circular object and a bubble with more detailed annotation is raised. Once the window is closed `mutation_scatter_plot` can continue and finish.

![BA2-4th-round-of-sort__G6.BA2.WTref.gofasta.aa.matplotlib.screenshot.aa.F490S.jpg](images/BA2-4th-round-of-sort__G6.BA2.WTref.gofasta.aa.matplotlib.screenshot.aa.F490S.jpg)
![BA2-4th-round-of-sort__G6.BA2.WTref.gofasta.aa.matplotlib.screenshot.aa.N501Y.jpg](images/BA2-4th-round-of-sort__G6.BA2.WTref.gofasta.aa.matplotlib.screenshot.aa.N501Y.jpg)
![BA2-4th-round-of-sort__G6.BA2.WTref.gofasta.codon.matplotlib.screenshot.codon.Q498R.jpg](images/BA2-4th-round-of-sort__G6.BA2.WTref.gofasta.codon.matplotlib.screenshot.codon.Q498R.jpg)
![BA2-4th-round-of-sort__G6.BA2.WTref.gofasta.codon.matplotlib.screenshot.codon.N501Y.jpg](images/BA2-4th-round-of-sort__G6.BA2.WTref.gofasta.codon.matplotlib.screenshot.codon.N501Y.jpg)
![spikenuc1207.native2ascii.no_junk.clean.mafft.frequencies.aa15INS.jpg](images/spikenuc1207.native2ascii.no_junk.clean.mafft.frequencies.aa15INS.jpg)
![spikenuc1207.native2ascii.no_junk.clean.mafft.frequencies.codon.15INS.jpg](images/spikenuc1207.native2ascii.no_junk.clean.mafft.frequencies.codon.15INS.jpg)
![spikenuc1207.native2ascii.no_junk.clean.mafft.frequencies.codon.N501Y.jpg](images/spikenuc1207.native2ascii.no_junk.clean.mafft.frequencies.codon.N501Y.jpg)
![spikenuc1207.native2ascii.no_junk.clean.mafft.frequencies.codon.Q498R.zoom.jpg](images/spikenuc1207.native2ascii.no_junk.clean.mafft.frequencies.codon.Q498R.zoom.jpg)
![spikenuc1207.native2ascii.no_junk.clean.mafft.frequencies.codon.N501Y.zoom.jpg](images/spikenuc1207.native2ascii.no_junk.clean.mafft.frequencies.codon.N501Y.zoom.jpg)
![spikenuc1207.native2ascii.no_junk.clean.mafft.frequencies.codon.INS15P.jpg](images/spikenuc1207.native2ascii.no_junk.clean.mafft.frequencies.codon.INS15P.jpg)
![spikenuc1207.native2ascii.no_junk.clean.mafft.frequencies.codon.P26DEL.jpg](images/spikenuc1207.native2ascii.no_junk.clean.mafft.frequencies.codon.P26DEL.jpg)
![spikenuc1207.native2ascii.no_junk.clean.mafft.frequencies.codon.V143DEL.jpg](images/spikenuc1207.native2ascii.no_junk.clean.mafft.frequencies.codon.V143DEL.jpg)
![spikenuc1207.native2ascii.no_junk.clean.mafft.frequencies.codon.V483DEL.jpg](images/spikenuc1207.native2ascii.no_junk.clean.mafft.frequencies.codon.V483DEL.jpg)


## Example of a dynamic output HTML pages with images (with interactive features)

The following were rendered by `mutation_scatter_plot`. They are not so visually appealing like the figures from matplotlib but they work in any www browser. The following figure is from Firefox.

![BA2-4th-round-of-sort__G6.BA2.WTref.gofasta.aa.firefox.screenshot.aa.F490S.jpg](images/BA2-4th-round-of-sort__G6.BA2.WTref.gofasta.aa.firefox.screenshot.aa.F490S.jpg)
![spikenuc1207.native2ascii.no_junk.clean.mafft.frequencies.firefox.screenshot.codon.png](images/spikenuc1207.native2ascii.no_junk.clean.mafft.frequencies.firefox.screenshot.codon.png)

Browse the live HTML with Javascript file contents at:

[data/outputs/aa/html/BA2-4th-round-of-sort__G6.BA2.WTref.aa.frequencies.html](https://host-patho-evo.github.io/mutation_scatter_plot/aa/html/BA2-4th-round-of-sort__G6.BA2.WTref.aa.frequencies.html)

[data/outputs/codon/html/BA2-4th-round-of-sort__G6.BA2.WTref.codon.frequencies.html](https://host-patho-evo.github.io/mutation_scatter_plot/codon/html/BA2-4th-round-of-sort__G6.BA2.WTref.codon.frequencies.html)

#### Command line output

Both `calculate_codon_frequencies` and `mutation_scatter_plot` output some helpful _Info:_ text on their STDOUT. One can also enable `--debug` option with debugging level.

When matplotlib raises its interactive image window and user points the mouse pointer some circular object in the chart (triggering the mouse `hover()` event) the `mutation_scatter_plot` writes on the STDOUT the values parsed for the codon or amino acid, for example:

```
$ calculate_codon_frequencies --alignment-file tests/inputs/test2_full.fasta --outfile-prefix tests/outputs/test2_full.x_after_count --padded-reference --reference-infile tests/inputs/MN908947.3_S_full.fasta --x-after-count --overwrite
Info: consensus = ACAGGCTGCGTTATAGCTTGGAATTCTAACAAGCTTGATTCTAAGGTTGGTGGTAATTATAATTACCTGTATAGATTGTTTAGGAAGTCTAATCTCAAACCTTTTGAGAGAGATATTTCAACTGAAATCTATCAGGCCGGTAACAAACCTTGTAATGGTGTTGCAGGTTTTAATTGTTACTTTCCTTTACAATCATATGGTTTCCGACCCACTTATGGTGTTGGTCACCAACCATACAGAGTAGTAGTACTTTCTTTTGAACTTCTACATGCACCAGCAACTGTTTGTGGACCTAAA
Info: Sample consensus sequence should roughly match substring inside the reference ATG---...
```

```
$ mutation_scatter_plot --tsv tests/outputs/test2_full.x_after_count.frequencies.tsv --outfile-prefix tests/outputs/test2_full.x_after_count.scatter_codons --show-STOP --show-X --show-DEL --show-INS --threshold=0.001
Info: Using BLOSUM80 matrix now. Theoretical minimum score is -6.0, theoretical maximum score is 11.0, values are {0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 11.0, -1.0, -6.0, -5.0, -4.0, -3.0, -2.0}
Info: Parsing input file tests/outputs/test2_full.x_after_count.frequencies.tsv
Info: Autodetected new TSV file format with a header in tests/outputs/test2_full.x_after_count.frequencies.tsv
Info: The file tests/outputs/test2_full.x_after_count.frequencies.tsv contains now these columns: Index(['padded_position', 'position', 'original_aa', 'mutant_aa', 'frequency',
       'original_codon', 'mutant_codon', 'observed_codon_count',
       'total_codons_per_site'],
      dtype='object')
Info: Originally there were 23 rows but after discarding codons with [N n] there are only 21 left
Info: Writing into tests/outputs/test2_full.x_after_count.scatter_codons.actually_rendered.tsv
Info: Title will be tests/outputs/test2_full.x_after_count
Info: Writing into tests/outputs/test2_full.x_after_count.scatter_codons.BLOSUM80.area_scaling.codon.frequencies.colors.tsv
Info: The following values were collected from matrix BLOSUM80 based on the actual data (some values from matrix might not be needed for your data, hence are not listed here): [-6, -4, -3, -2, -1, 0, 1, 2] . Range spans 9 values (before symmetrization).
Info: Writing into tests/outputs/test2_full.x_after_count.scatter_codons.BLOSUM80.area_scaling.amino_acid_changes.html
Info: Writing into tests/outputs/test2_full.x_after_count.scatter_codons.BLOSUM80.area_scaling.amino_acid_changes.png, figure size is [16.  9.] inches and [1600.  900.] dpi
Info: Writing into tests/outputs/test2_full.x_after_count.scatter_codons.BLOSUM80.area_scaling.amino_acid_changes.pdf, figure size is [16.  9.] inches and [1600.  900.] dpi

$ # Example hover output for position 498
Info: _padded_position=515, ypos=21
Info: _new_codon=CGA, _padded_position=515, _position_in_protein=498, _frequency=0.700161
Info: 1 aa residues observed in position 498:
       padded_position  position original_aa mutant_aa  frequency original_codon mutant_codon  observed_codon_count  total_codons_per_site
21                 515       498           Q         R   0.700161            CAA          CGA               1103909                1576651
```
`

## Command line arguments

```
$ calculate_codon_frequencies --help
usage: calculate_codon_frequencies [-h] [--infilename INFILENAME] [--infile-format INFILE_FORMAT] [--out-prefix OUT_PREFIX] [--region-start REGION_START] [--region-end REGION_END]
                                  [--x-after-count] [--print-unchanged-sites] [--disable-print-unchanged-sites] [--discard-this-many-leading-nucs DISCARD_THIS_MANY_LEADING_NUCS]
                                  [--discard-this-many-trailing-nucs DISCARD_THIS_MANY_TRAILING_NUCS] [--minimum-alignments-length MINIMUM_ALN_LENGTH] [--debug DEBUG] [--overwrite]
                                  [--threads THREADS]

Identify codons at each alignment position of input sequences and output frequency of every mutant codon. Parent codon/amino acid is identified as the most frequent one.

options:
  -h, --help            show this help message and exit
  --infilename INFILENAME
                        Input FASTA/Q file path. (default: )
  --infile-format INFILE_FORMAT
                        Input FASTA/Q file format [default: fasta]. Outfile is always fasta. (default: fasta)
  --out-prefix OUT_PREFIX
                        Output file prefix, eventually also with path. (default: )
  --region-start REGION_START
                        Start parsing of the sequence [default: none] eventually after discarding leading nucleotides (see '--discard-this-many-leading-nucs'). Default: 0 (default: 0)
  --region-end REGION_END
                        Stop parsing of the sequence [default: none]. Default: 0 (parse until the very end) (default: 0)
  --x-after-count       The FASTA file ID contains the count value followed by lowercase 'x' (default: False)
  --print-unchanged-sites
                        Print out also sites with unchanged codons in to unchanged_codons.tsv file [default] (default: True)
  --disable-print-unchanged-sites
                        Do NOT print out sites with unchanged codons to unchanged_codons.tsv file (default: False)
  --discard-this-many-leading-nucs DISCARD_THIS_MANY_LEADING_NUCS
                        Specify how many offending nucleotides are at the front of the FASTA sequences shifting the reading frame of the input FASTA file from frame +1 to either of the two
                        remaining. Count the leading dashes and eventual nucleotides of incomplete codons too and check if it can be divided by 3.0 without slack. By default reading
                        frame +1 is expected and hence no leading nucleotides are discarded. Default: 0 (default: 0)
  --discard-this-many-trailing-nucs DISCARD_THIS_MANY_TRAILING_NUCS
                        Specify how many offending nucleotides are at the end of each sequence. Default: 0 (default: 0)
  --minimum-alignments-length MINIMUM_ALN_LENGTH
                        Minimum length of aligned NGS read to be used for calculations (default: 50)
  --debug DEBUG         Set debug level to some real number (default: 0)
  --overwrite           Overwrite existing output files instead of raising RuntimeError (default: False)
  --threads THREADS     Number of CPU threads/processes to use for parallelization. If 0 [default], the tool automatically detects the allocation from environment variables (tried in
                        order: PBS_NUM_PPN, PBS_NCPUS, OMP_NUM_THREADS) or falls back to all available cores via multiprocessing.cpu_count(). (default: 0)
```

```
$ mutation_scatter_plot --help
usage: mutation_scatter_plot [options] [-h] [--tsv TSV_FILE_PATH] [--column COLUMN_WITH_FREQUENCIES] [--outfile-prefix OUTFILE_PREFIX] [--offset OFFSET] [--xmin XMIN] [--xmax XMAX]
                            [--x-axis-bins XAXIS_BINS] [--x-axis-major-ticks-spacing XAXIS_MAJOR_TICKS_SPACING] [--x-axis-minor-ticks-spacing XAXIS_MINOR_TICKS_SPACING]
                            [--x-axis-label-start XAXIS_LABEL_START] [--aminoacids] [--show-STOP] [--show-INS] [--show-DEL] [--show-X] [--enable-colorbar]
                            [--disable-short-legend] [--include-synonymous] [--threshold THRESHOLD] [--title TITLE] [--disable-2nd-Y-axis] [--disable-showing-bokeh]
                            [--disable-showing-mplcursors] [--legend] [--matrix MATRIX] [--matrix-file MATRIX_FILE] [--colormap COLORMAP] [--dpi DPI] [--backend BACKEND]
                             [--debug DEBUG] [--disable-bokeh-sqrt-size] [--linear-circle-size] [--show-invisible-placeholder-dots]

Render scatter figures of amino acid or codon mutation frequencies from TSV files using Matplotlib and Bokeh

options:
  -h, --help            show this help message and exit
  --tsv TSV_FILE_PATH   Path to a TAB separated file with 5-columns: ['position', 'original_aa', 'mutant_aa', 'frequency', 'original_codon', 'mutant_codon', 'coverage_per_codon'] or more
                        up to 11-columns ['padded_position', 'position', 'original_aa', 'mutant_aa', 'frequency', 'original_codon', 'mutant_codon', 'observed_codon_count',
                        'total_codons_per_site', 'frequency_parent', 'frequency_selected'] as we kept extending the file format (default: )
  --column COLUMN_WITH_FREQUENCIES
                        Name of a column in input TSV to be used for rendering frequencies [frequency] (default: frequency)
  --outfile-prefix OUTFILE_PREFIX
                        Output file prefix, eventually also with path. Output files will be PNG, JPG or PDF, HTML (default: )
  --offset OFFSET       Define offset value to be added to every amino acid position in the first column of the TSV input file at runtime. This is to obtain real amino acid position
                        within the protein. Normally protein starts at position 1. Check the first aa in TSV file and provide whatever number to be added to it to get the desired
                        amino acid position in the full-length protein. (default: 0)
  --xmin XMIN           Define minimum X-axis value. This should be the position in the padded alignment (using '-'). (default: 0)
  --xmax XMAX           Define maximum X-axis value. This should be the position in the padded alignment (using '-'). (default: 0)
  --x-axis-bins XAXIS_BINS
                        Set number of bins (labels) on the X-axis. Could be used to override the amount of major ticks in a different way. Use 20 for example. [default is 0] (default:
                        0)
  --x-axis-major-ticks-spacing XAXIS_MAJOR_TICKS_SPACING
                        Set distance between the major ticks on X-axis (default: 10)
  --x-axis-minor-ticks-spacing XAXIS_MINOR_TICKS_SPACING
                        Set distance between the minor ticks on X-axis (default: 1)
  --x-axis-label-start XAXIS_LABEL_START
                        Set first label value on the X-axis (default: 0)
  --aminoacids          Draw chart with amino acid residues on Y-axis instead of codons. [default is False] (default: False)
  --show-STOP           Include STOP codons or '*' in charts on Y-axis. [default is False] (default: False)
  --show-INS            Include INS in charts on Y-axis. [default is False] (default: False)
  --show-DEL            Include DEL in charts on Y-axis. [default is False] (default: False)
  --show-X              Include X in charts on Y-axis in --aminoacids mode. [default is False] (default: False)
  --enable-colorbar     Enable colorbar [is Disabled by default] (default: False)
  --disable-short-legend
                        Disable short legend in charts on X-axis. [is Enabled by default] (default: False)
  --include-synonymous  Include synonymous changes in --aminoacids output as green diamonds. In codon output they are always shown. [default is False] (default: False)
  --threshold THRESHOLD
                        Define minimum frequency threshold to display a pictogram in the output. For codon mode use 0.001 and for aa mode use 0.01. [default: 0.001] (default: 0.001)
  --title TITLE         Set title for the figures, by default trailing '.frequencies.tsv' is stripped from the end of the input filename (default: )
  --disable-2nd-Y-axis  Disable 2nd Y-axis on the right side of the plot. [default is False] (default: False)
  --disable-showing-bokeh
                        Disable calling a GUI or text-based browser TO show the rendered Bokeh HTML+Javascript. [default is False] (default: False)
  --disable-showing-mplcursors
                        Disable rendering an interactive window WITH hover() events recognized by mplcursors. [default is False] (default: False)
  --legend              Draw legend chart. [default is False] (default: False)
  --matrix MATRIX       BLOSUM matrix: BLOSUM45,BLOSUM50,BLOSUM62,BLOSUM80,BLOSUM90 [default is BLOSUM80] (default: BLOSUM80)
  --matrix-file MATRIX_FILE
                        Matrix file compatible with BLOSUM matrices, e.g. prot26050-sup-0004-Supinfo03.sm from https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8641535/bin/NIHMS1664401-
                        supplement-supinfo.rar if you do not like default BLOSUM62 (default: None)
  --colormap COLORMAP   Pick a colormap recognized by matplotlib. See https://i.sstatic.net/cmk1J.png [default is coolwarm_r but seismic_r and coolwarm_r are great too] (default:
                        amino_acid_changes)
  --dpi DPI             DPI resolution for images (default: 600)
  --backend BACKEND     Matplotlib backend to render resulting figures: agg, wxpython, pyqt5, pyqt6, pycairo, cairocffi [default: unset] To disable Matplotlib interactive window being
                        raised up you can set MPLBACKEND=agg env variable. (default: )
  --debug DEBUG         Set debug to some real number (default: 0)
  --disable-bokeh-sqrt-size
                        Disable sqrt scaling for Bokeh circle sizes; size (diameter) will be proportional to frequency, area proportional to frequency². By default sqrt scaling is on,
                        matching the perceptual appearance of the matplotlib figure. (default: True)
  --linear-circle-size  Render circles with radius proportional to frequency in BOTH Matplotlib and Bokeh. Matplotlib: s = freq² × 5000 so that radius ∝ frequency. Bokeh: diameter =
                        freq × 100. Implies --disable-bokeh-sqrt-size so that both figure types use the same linear-radius formula and remain visually in sync. Matches the v0.2 Bokeh
                        rendering. The scaling mode is encoded in output filenames as 'linear_scaling' (default: False)
  --show-invisible-placeholder-dots
                        Include below-threshold dots in the plot. [default: False] (default: False)
```

```
$ count_motifs_in_sequences --help
Usage: count_motifs_in_sequences [options]

Options:
  --version             show program's version number and exit
  -h, --help            show this help message and exit
  --infilename=INFILENAME
                        Input FASTA/Q file path.
  --infile-format=INFILE_FORMAT
                        Input FASTA/Q file format [default: fasta]. Outfile is
                        always fasta.
  --motif=MOTIF         Protein or nucleotide motif to count. No REGEXPs
                        allowed (yet)
  --start-position=STARTPOS
                        Exact start position of the query
  --end-position=ENDPOS
                        Exact end position of the query
  --debug=DEBUG         Set debug to some value
```

```
$ alignment2dots --help
Usage: alignment2dots [options]

Options:
  --version             show program's version number and exit
  -h, --help            show this help message and exit
  --reference-infile=FILE
                        FASTA formatted input file with reference padded
                        sequence or not
  --alignment-file=FILE
                        FASTA formatted input file with padded sequence or not
  --aminoacids          FASTA formatted input files are protein sequences
                        instead of DNA [default: DNA]
  --outfilename=FILE    Output filename. If the filename ends with .tsv the
                        output lines will be split using TABs and original
                        FASTA ID will be on the same line with the split
                        sequence
  --aln_start=ALN_START
                        First nucleotide of the ORF region of interest to be
                        sliced out from the input sequences
  --aln_stop=ALN_STOP   Last nucleotide of the last codon of interest to be
                        sliced out from the input sequences
  --print-fasta-ids     The FASTA file ID should be printed for each aligned
                        entry
  --threshold=THRESHOLD
                        Set the minimum absolute count of different characters
                        in sequence to be output. Set this to 1 or higher if
                        you want to see at least 2 aa residues being changed
                        and in turn, do not want to see just dashes in some
                        cases for synonymous changes [default: 0].
  --relative_threshold=RELATIVE_THRESHOLD
                        Set the minimum relative incidence threshold of the
                        whole sequence hash. Maybe you want something like
                        0.001 [default: 0].
  --top_n=TOP_N         Write only first N entries containing some difference
                        [default: all]
  --debug=DEBUG         Set debug level to some real number
```

```
$ split_fasta_entries_by_lengths --help
Usage: split_fasta_entries_by_lengths [options]

Options:
  --version             show program's version number and exit
  -h, --help            show this help message and exit
  --infile=INFILE       Input FASTA/Q file path.
  --outfile-prefix=OUTFILE_PREFIX
                        Output file path prefix.
  --full-length=FULL_LENGTH
                        Full length required for perfect alignment [0]
  --format=FORMAT       Input file format.
  --debug=DEBUG         Set debug to some value
```

## Project Evolution & Key Improvements

This project has undergone several major architectural and performance transformations through a series of focused optimization sessions. Below is a summary of the key milestones and improvements:

| Stage | Focus Area | Key Architectural & Performance Milestones |
| :--- | :--- | :--- |
| **1. Quality** | **Refactoring** | Migrated to modern **f-strings**, updated **Matplotlib 3.8+ colormap** syntax, and achieved a **10.0/10 Pylint score**. |
| **2. UX** | **CLI & Docs** | Comprehensive **`--help`** refresh, synchronized realistic **`test2_full`** examples in README, and automated help-output verification. |
| **3. Reliability** | **Testing** | Established a **"Golden File" regression suite** with 120+ artifacts and standardized naming. Refactored **`alt_translate`** to the core module. |
| **4. Core Speed** | **Algorithms** | **Streaming FASTA** input, **Global/Local grouping**, and **Render Pruning**, achieving a **119x core speedup**. |
| **5. Scale** | **Parallelism** | **NumPy Vectorization**, **Multi-core multiprocessing**, and **High-speed FASTA parser** with **Raw Grouping**. |

### Major Technical Milestones
- **Streaming & Memory Efficiency**: Transitioned from loading entire alignments to streaming unique sequence hashes, allowing the tool to process files with millions of reads in constant memory.
- **Precision Hardening**: Implemented a project-wide `decimal.Decimal` model to eliminate floating-point drift and ensure identical results across different OS platforms.
- **Static Analysis**: Integrated native **Pyright/Pylance** support in `pyproject.toml` for robust IDE development and import resolution.
- **Render Optimization**: Implemented $O(1)$ hover metadata lookups for both Matplotlib and Bokeh, eliminating rendering lag in interactive plots with thousands of mutations.

## Version 0.2 supports DELetions and INSertions relative to the (padded) reference sequence

We further improved the software to be able to report DELetions and INSertions appearing in sample data (multiple-sequence alignment). The `calculate_codon_frequencies` now also reports total counts of reads covering each codon (per-site coverage) in additional columns 8 and 9 of the TSV output file. The more detailed file can be parsed by `mutation_scatter_plot`. Another significant change was the requirement for padded alignment at input which must have exactly same length as the reference (which might need to be padded as well). To observe DELetions in the sample sequence one does not need to adjust the reference sequence, because in the aligned sample sequence will be just `---` for a DELeted codon. Obviously, if an INSertion is to be reflected in the sample sequence, the reference sequence must be inflated by paddings. The reason for that is that we use a pairwise NCBI blastn to create the alignment and do not create a multiple-sequence alignment at all (the hints mentioning `gofasta` are now removed).

![GISAID SARS-CoV-2 spikenuc1207 frequencies of mutations in S protein](images/spikenuc1207.native2ascii.no_junk.clean.mafft.frequencies.aa.jpg)
![GISAID SARS-CoV-2 spikenuc1207 frequencies of mutations in S protein at aa position 501](images/spikenuc1207.native2ascii.no_junk.clean.mafft.frequencies.aa.jpg)

## Citation

Please cite the following article if you use our data or software in your research:

Shoshany A., Tian R., Padilla-Blanco M., Hruška A., Baxová K., Zoler E., Mokrejš M., Schreiber G., Zahradník J. (submitted) In Vitro and Viral Evolution Convergence Reveal the Selective Pressures Driving Omicron Emergence. [bioRxiv](https://www.biorxiv.org/content/10.1101/2025.04.23.650148v1)


## Website

https://github.com/host-patho-evo/mutation_scatter_plot


## License

This work © 2025-2026 by Jiří Zahradník and Martin Mokrejš (First Medical Faculty - Charles University in Prague) is licensed under Creative Commons Attribution 4.0 International (CC BY 4.0). To view a copy of this license, visit https://creativecommons.org/licenses/by/4.0/


## Acknowledgements
[This project was supported by the National Institute of Virology and Bacteriology (Programme EXCELES, LX22NPO5103) - funded by the European Union - NextGenerationEU](https://nivb.cz/en/)

![logos/loga_hlavicka_colour_ENG.png](logos/loga_hlavicka_colour_ENG.png)Funded by the European Union NextGenerationEU -- Czech Recovery Plan -- Ministry of Education, Youth and Sports
