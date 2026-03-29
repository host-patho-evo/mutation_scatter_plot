# Sequence Deduplication Pipeline

This document describes the set of standalone Python helper scripts in `scripts/`
that implement a multi-stage sequence deduplication and traceability pipeline
for large FASTA datasets (e.g. GISAID SARS-CoV-2 spike nucleotide sequences).

---

## Pipeline Overview

```
filename_prefix.fasta                        (raw GISAID download, millions of records)
        │
        │  count_same_sequences.py
        ▼
filename_prefix.counts.fasta                 (deduplicated; NNNNx.sha256 IDs)
filename_prefix.sha256_to_ids.tsv            (sha256 → original FASTA ID mapping)
        │
        │  external alignment (e.g. MAFFT/BLASTN), then
        │  filter for successful alignments
        ▼
filename_prefix.counts.clean.fasta           (alignment-filtered subset)
        │
        │  kick.py --full-length=N
        ▼
filename_prefix.counts.clean.exactly_N.fasta  (exactly N nt, after padding removal)
filename_prefix.counts.clean.shorter_N.fasta  (shorter than N nt)
filename_prefix.counts.clean.longer_N.fasta   (longer than N nt)
        │
        │  create_list_of_discarded_sequences.py --inverted
        ▼
filename_prefix.counts.clean.exactly_N.discarded_original_ids.txt
```

The naming convention is purposeful: each pipeline stage **appends a
dot-separated suffix** to the output name, making the parent→child
relationship machine-readable (see `summarize_fasta_pipeline.py`).

---

## File Naming Convention

All scripts follow a consistent naming rule:

| File pattern | Produced by | Content |
|---|---|---|
| `filename_prefix.fasta` | upstream source | raw FASTA (GISAID format) |
| `filename_prefix.counts.fasta` | `count_same_sequences.py` | deduplicated; `NNNNx.sha256hex` IDs |
| `filename_prefix.sha256_to_ids.tsv` | `count_same_sequences.py` | sha256 → original FASTA ID map |
| `filename_prefix.counts.clean.fasta` | alignment filter | sequences that aligned successfully |
| `filename_prefix.counts.clean.exactly_N.fasta` | `kick.py` | sequences of exactly N nt |
| `filename_prefix.counts.clean.shorter_N.fasta` | `kick.py` | sequences shorter than N nt |
| `filename_prefix.counts.clean.longer_N.fasta` | `kick.py` | sequences longer than N nt |
| `filename_prefix.counts.clean.exactly_N.discarded_original_ids.txt` | `create_list_of_discarded_sequences.py` | original IDs of discarded sequences |

> **Note:** Replace `filename_prefix` with your actual file base name. These
> scripts never embed literal dataset names — the prefix is always provided via
> command-line arguments.

---

## Scripts

### `count_same_sequences.py`

**Purpose:** Deduplicate a FASTA file and annotate each unique sequence with
its occurrence count and SHA-256 checksum.

**Mechanism:** Pipes the input through `reformat.sh` (BBTools) which uppercases
sequences and strips alignment dashes, then through `sort | uniq -c | sort -nr`
to count duplicates. Each unique sequence is written as:

```
>NNNNx.sha256hex
SEQUENCE
```

where `NNNN` is the count and `sha256hex` is the 64-character hexadecimal
SHA-256 of the uppercase, dash-stripped sequence.

**Simultaneously** it performs a second pass over the original FASTA to build
a `sha256 → [original_ID_1, original_ID_2, …]` mapping TSV
(`filename_prefix.sha256_to_ids.tsv`). This TSV enables tracing any
deduplicated record back to all original sequences that share its content.

**SHA-256 normalization:** Both passes normalize sequences identically —
uppercase, alignment dashes stripped, no trailing newline — to ensure the
hash embedded in the `NNNNx.sha256hex` ID matches the TSV lookup key.

**Key options:**

| Option | Default | Description |
|---|---|---|
| `--infilename` | required | Input FASTA/Q |
| `--outfile-prefix` | derived from `--infilename` | Prefix for `.counts.fasta` |
| `--mapping-outfile` | `{prefix}.sha256_to_ids.tsv` | TSV output path |
| `--overwrite` | off | Clobber existing outputs |
| `--min-count N` | 0 | Only write sequences seen ≥ N times |
| `--top-n N` | 0 | Only write the top N most-frequent sequences |
| `--sort-bucket-size` | `30%` | `sort -S` argument (RAM or %, e.g. `10G`) |

**Input validation:**
- Rejects `--outfile-prefix` or `--mapping-outfile` values that contain `=`,
  which catches the common typo
  `--outfile-prefix=infilename=filename_prefix` (should be
  `--outfile-prefix=filename_prefix`).

**Timestamp-aware guard (make-style):**

| State | Behaviour |
|---|---|
| Both outputs exist and are *newer* than input | `Info: up-to-date, skipping` → exit 0 |
| Any output exists and is *older* than input | `RuntimeError: stale … use --overwrite` |
| `--overwrite` is set | Always proceeds and clobbers |

---

### `create_list_of_discarded_sequences.py`

**Purpose:** Given a deduplicated counts FASTA and (optionally) the original
FASTA, produce a text file listing the original FASTA IDs corresponding to
sequences that were either *kept* or *discarded* at a pipeline stage.

**Two modes:**

| Mode | `--infilename` is | Output is |
|---|---|---|
| **Default** | the *discarded* file | original IDs whose sha256 **is** in `--infilename` |
| **Inverted** (`--inverted`) | the *kept* file | original IDs whose sha256 is **absent** from `--infilename` |

**Source for original IDs — three paths (fastest first):**

1. **TSV fast path** (`--mapping-outfile`): reads the pre-built
   `sha256_to_ids.tsv` table produced by `count_same_sequences.py`. If not
   supplied explicitly, the script auto-detects it next to `--infilename`.
2. **FASTA scan** (`--original-infilename`): scans the original FASTA,
   computes sha256 per record, and compares. Required for `--inverted` when
   the TSV is absent.
3. **`--inverted` + TSV**: `--mapping-outfile` can also be used in inverted
   mode — the TSV is scanned to find sha256s *not* present in `--infilename`.

**Output:** A plain-text file with one original FASTA ID per line, defaulting
to `{stem}.discarded_original_ids.txt`.

**Key options:**

| Option | Description |
|---|---|
| `--infilename` | Deduplicated counts FASTA (kept or discarded file) |
| `--original-infilename` | Original (pre-deduplication) FASTA |
| `--mapping-outfile` | TSV mapping file (`sha256_to_ids.tsv`) — fast path |
| `--inverted` | Emit IDs *absent* from `--infilename` (discarded) |
| `--outfile` | Output path; defaults to `{stem}.discarded_original_ids.txt` |
| `--overwrite` | Clobber existing output |

**Timestamp-aware guard:** Same three-state make-style logic as
`count_same_sequences.py`. The guard considers the mtime of all input files
(`--infilename`, `--original-infilename`, `--mapping-outfile` if given).
Passing `--outfile=/dev/null` bypasses the guard entirely (used internally by
`summarize_fasta_pipeline.py` for stats-only calls — the output message in
this mode reads *"IDs not written anywhere"* rather than *"wrote … to /dev/null"*).

---

### `kick.py`

**Purpose:** Split a deduplicated counts FASTA by sequence length into three
output files: `exactly_N`, `shorter_N`, and `longer_N`.

**Input:** A FASTA file in `NNNNx.sha256` or any standard format, where
sequences are padded with dashes (`-`). Padding is **preserved** in the output
so that downstream alignment coordinates remain valid.

**Key options:**

| Option | Description |
|---|---|
| `--infile` | Input FASTA/Q |
| `--outfile-prefix` | Prefix for the three output files |
| `--full-length N` | Target length in nucleotides (required) |
| `--format` | Input format (`fasta` by default) |
| `--overwrite` | Clobber existing outputs |

**Timestamp-aware guard:** All three outputs must be newer than the input for
the skip to trigger.

---

### `summarize_fasta_pipeline.py`

**Purpose:** Audit the entire pipeline in one run — scan every matching FASTA
file, compute record counts and NNNNx sum, infer the parent→child pipeline
topology from filename suffixes, and (optionally) invoke
`create_list_of_discarded_sequences.py` for each step to show discarded-ID
statistics inline.

**Usage:**

```bash
summarize_fasta_pipeline.py <search_path> <filename_prefix> [options]
```

**Table columns** (printed to stdout after the scanning phase):

| Column | Description |
|---|---|
| `File` | Relative path from `<search_path>` |
| `Modified` | File modification time (`YYYY-MM-DD HH:MM`) |
| `Records` | Number of `>` header lines |
| `ΔRecords` | Record-count change vs. direct parent stage |
| `Sum of NNNNx` | Sum of the integer count prefixes across all IDs |
| `ΔSum` | NNNNx-sum change vs. direct parent stage |

**Parent→child inference:**

A file `B` is considered a direct child of `A` when
`strip_fasta_suffix(basename(B)).startswith(strip_fasta_suffix(basename(A)) + ".")`.
The *longest* matching ancestor is chosen as the direct parent (so
`filename_prefix.counts.clean` is a child of `filename_prefix.counts`, not of
`filename_prefix`).

**Discard statistics — four-tier cache (fastest first):**

1. **Existing `.discarded_original_ids.txt`** — newer than both FASTAs:
   read line count and NNNNx sum directly. Near-instant.
2. **Existing `.sha256_to_ids.tsv`** — newer than parent FASTA:
   invoke `create_list_of_discarded_sequences.py --mapping-outfile`. Fast.
3. **Auto-generated `.sha256_to_ids.tsv`** — if no fresh TSV exists, the
   parent FASTA is scanned in-process (no external tools required) and a new
   `{parent_base}.sha256_to_ids.tsv` is written alongside the FASTA. The TSV
   is then used immediately and serves as a cache for all subsequent calls.
4. **Fallback FASTA scan**: `create_list_of_discarded_sequences.py
   --original-infilename`. Only reached if TSV auto-generation fails.

All cache lookups compare file modification times (make-style): a file that
exists but is older than its source is treated as stale and the next tier is
used.

**Options:**

| Option | Description |
|---|---|
| `--no-discard-stats` | Print only the table; skip the discard-script calls |
| `--add-missing-checksums-to-fasta-files` | When a pipeline FASTA has legacy `NNNNx` IDs (no sha256 embedded), rewrite it in-place with `NNNNx.sha256hex` IDs and rename the original to `.fasta.orig`. Skipped silently if a `.fasta.orig` already exists. Makes all future analyses faster: sha256 can be read from the ID directly instead of being recomputed from the sequence. |

---

### `check_alignment_trimming.py`

**Purpose:** Diagnostic tool to quantify how many sequences were **shortened**
during the alignment step. After alignment the padded-then-stripped sequence
may differ from the pre-alignment sequence, causing the SHA-256 embedded in
the `NNNNx.sha256` ID to no longer match the sequence content. This explains
any residual discrepancy between pre- and post-alignment counts when using the
TSV fast path.

Run it on the aligned FASTA to count mismatches between the ID-embedded sha256
and the sha256 recomputed from the aligned sequence:

```bash
check_alignment_trimming.py --infilename=filename_prefix.counts.clean.fasta
```

---

## Mixed-Encoding Handling

GISAID FASTA files contain non-ASCII characters in sample descriptions (country
names, institute names, researcher names). These appear as Latin-1 or Latin-2
encoded bytes (e.g. `0xED = í`, `0xE9 = é`) which are not valid UTF-8.

All FASTA reading in this pipeline is done in **binary mode** with a
**per-line UTF-8 → Latin-1 fallback**:

```python
def _decode_fasta_line(raw: bytes) -> str:
    try:
        return raw.decode("utf-8")
    except UnicodeDecodeError:
        return raw.decode("latin-1")  # never fails — every byte is valid
```

Since only the **sequence** (pure ASCII: ACGT + N + IUPAC) and the **first
word of the header** (pure ASCII accession) are used for SHA-256 computation
and output, the encoding of the description fields does not affect correctness.

---

## Common Workflow

```bash
# Step 1 – deduplicate
count_same_sequences.py \
    --infilename=filename_prefix.fasta \
    --outfile-prefix=filename_prefix \
    --mapping-outfile=filename_prefix.sha256_to_ids.tsv

# Step 2 – align (external tool, e.g. blastn or mafft)
# produces: filename_prefix.counts.clean.fasta

# Step 3 – split by length
kick.py \
    --infile=filename_prefix.counts.clean.fasta \
    --outfile-prefix=filename_prefix.counts.clean \
    --full-length=3822

# Step 4 – list discarded IDs for the length-filtered subset
create_list_of_discarded_sequences.py \
    --infilename=filename_prefix.counts.clean.exactly_3822.fasta \
    --original-infilename=filename_prefix.counts.clean.fasta \
    --inverted

# Step 5 – audit the whole pipeline at once
summarize_fasta_pipeline.py . filename_prefix
```

Re-running any step is safe: scripts skip silently when all outputs are
up-to-date and fail with a clear error when outputs are stale, prompting you
to add `--overwrite` only when you genuinely intend to regenerate.
