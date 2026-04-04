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
        │  split_fasta_entries_by_lengths --full-length=N
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
| `filename_prefix.counts.clean.exactly_N.fasta` | `split_fasta_entries_by_lengths` | sequences of exactly N nt |
| `filename_prefix.counts.clean.shorter_N.fasta` | `split_fasta_entries_by_lengths` | sequences shorter than N nt |
| `filename_prefix.counts.clean.longer_N.fasta` | `split_fasta_entries_by_lengths` | sequences longer than N nt |
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

### Post-Alignment Header Structure (`.clean.fasta`)
If the sequences are subsequently processed through external alignment pipelines like `blastn` (e.g., via `processing4.sh`), the FASTA header is expanded into multiple tab-delimited columns. 

Crucially, the **very last item** appended to the end of the header line is the **padded reference sequence (`sseq`)** extracted directly from the true BLAST High-Scoring Segment Pair (HSP) alignment block.

*Example:*
```text
>NNNNx.sha256hex	<other_blast_columns>...	ATG--TTT-GTT...
A-G--T-T-GTT...
```

> [!WARNING]
> **Important Sequence Distinction:** 
> Beware that the sequence appended to the end of the FASTA header line is **NOT** the raw sequence from the root (e.g., `spikenuc1207.fasta`), nor is it the sequence from the initial deduplication step (`spikenuc1207.no_junk.counts.fasta`), nor is it the plain, non-padded reference sequence of the Spike protein. It was excised from the HSP pairwise-alignment block.
> 
> **It is a padded segment cut out from the aligned region** with its own unique length and dash placement. Furthermore, the DNA sequence data directly below it (the `SEQUENCE` body) is *also* padded, as it is extracted from that same HSP pairwise-alignment block.

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
this mode reads *"IDs not written anywhere"*; to persist them for a specific
step, run `create_list_of_discarded_sequences.py` directly for that pair).

---

### `split_fasta_entries_by_lengths`

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
| `Unique DNA entries` | Number of `>` header lines |
| `ΔDNA entries` | Record-count change vs. direct parent stage |
| `Sum of NNNNx` | Sum of the integer count prefixes across all IDs |
| `ΔSumToParent` | NNNNx-sum change vs. direct parent stage |
| `Discarded original unique Entries` | Unique sha256 groups discarded at this step |
| `Sum of discarded NNNNx` | Total individual sequences discarded (NNNNx expansion) |
| `Novel sha256s` | sha256 hashes in child but absent from parent's TSV |
| `Unique protein entries` | Unique sequences in companion `.prot.fasta` (if present) |
| `Sum of protein NNNNx` | NNNNx sum of the companion protein file |

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
| `--extract-discarded-fasta` | For each child file with a `.discarded_sha256_hashes.txt`, extract the corresponding original FASTA records from the root ancestor FASTA. The lookup uses the first GISAID-level ancestor's `sha256_to_ids.tsv` (mapping sha256 → original accession IDs). Output: `{child_stem}.discarded_original_entries.fasta`. Uses `filterbyname.sh` (BBTools) when available, otherwise falls back to a native Python byte-streaming extractor. |
| `--add-missing-checksums-to-fasta-files` | When a pipeline FASTA has legacy `NNNNx` IDs (no sha256 embedded), rewrite it in-place with `NNNNx.sha256hex` IDs and rename the original to `.fasta.orig`. Skipped silently if a `.fasta.orig` already exists. Makes all future analyses faster: sha256 can be read from the ID directly instead of being recomputed from the sequence. |
| `--cpu-bind {local,spread,none}` | NUMA CPU binding policy (default: `local`). `local` binds all threads to the NUMA node with the most free RAM before any I/O begins. `spread` leaves OS placement unchanged. `none` disables autobind. Silently skipped on single-node hosts and when a job scheduler has already set placement via `GOMP_CPU_AFFINITY`, `KMP_AFFINITY`, `SLURM_CPU_BIND`, `SGE_BINDING`, or `OMP_PROC_BIND`. |
| `--jobs N` | Number of parallel workers for Phase 0 (sha256 verification) and Phase 1 (FASTA scanning). **Default: 1 (sequential).** No auto-detection from environment variables: must be set explicitly. Recommended value: 4–8 for a pipeline with 5–20 FASTA files over NFS. |

> **Why `--jobs` and not `--threads`?**  
> `summarize_fasta_pipeline.py` parallelises *across files* using
> `concurrent.futures.ThreadPoolExecutor` (one thread per FASTA file).
> In Phase 0 and Phase 1, FASTA reads utilize C-level operations (e.g., `bytes.split()`)
> which natively release Python’s GIL, allowing true multi-core speedup for I/O bounds.
> The flag name follows the GNU `make -j N` / `parallel -j N` convention.
>
> In Phase 2 (`_verify_integrity_of_pairs_and_classify`), alignment verification relies
> heavily on Python string manipulation bytecode. The GIL forces sequential bytecode
> execution, causing concurrent Phase 2 validations to context-switch and effectively
> behave as if single-threaded (maxing out at 100% CPU). However, deploying
> `ThreadPoolExecutor` instead of `ProcessPoolExecutor` guarantees the 11.8+ GB 
> `sha256_sets` dictionary is safely shared among all parallel verifications.
> If processes were used, RAM cloning would trigger Out of Memory (OOM) kernel kills.

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

## Parallelization, NUMA & NFS Strategy

The pipeline is heavily optimized for multi-threaded setups via the `--jobs N`
flag in `summarize_fasta_pipeline.py`.

1. **Across-File Parallelism (ThreadPoolExecutor)**: Operations heavily
   rely on standard I/O (often fetching files that take gigabytes of disk space).
   Phase 1 drops to C-level byte processing which safely releases the Global Interpreter Lock (GIL),
   allowing threads to saturate total RAID or network bandwidth parallelly across multiple cores.
   Phase 2 runs sequential Python bytecode, forcing the GIL to fall back to effectively
   single-threaded throughput per Node. However, Threads strictly bypass memory duplication
   and securely share the global 12GB+ tracking RAM map.

2. **Within-File Design Avoided for NFS Performance**: The code explicitly
   **avoids** reading the *same* file via concurrent threads.  Performing
   concurrent distinct read-passes on NFSv4 wreaks havoc on the kernel’s
   read-ahead cache and causes severe network thrashing.  By completing one
   file entirely in a single thread pass, the pipeline allows NFS to stream
   giant files sequentially at maximum throughput.

3. **NUMA Auto-Bind** (`--cpu-bind local` default): On multi-socket NUMA
   hosts, the library automatically detects the NUMA topology via
   `/sys/devices/system/node/` and binds all threads + child processes to the
   NUMA node with the most free RAM before any phase begins.  On single-socket
   machines this is a silent no-op with zero overhead.  Binding stack (tried
   in priority order):

   | Back-end | Requires | Covers |
   |---|---|---|
   | `ctypes → libnuma.so.1` | `libnuma1` (Debian) / `numactl-libs` (RHEL) | CPU set + memory policy, in-process |
   | `numactl --cpunodebind=N --localalloc` | `numactl` in PATH | CPU set + memory policy, re-exec |
   | `os.sched_setaffinity()` | Python stdlib | CPU set only, fork-inherited |

   HPC job schedulers that already set placement (`GOMP_CPU_AFFINITY`,
   `KMP_AFFINITY`, `SLURM_CPU_BIND`, `SGE_BINDING`, `OMP_PROC_BIND`) cause
   autobind to be skipped silently.


---

## Encoding Pre-processing: `fix_fasta_encoding.py`

GISAID FASTA exports aggregate submissions from thousands of labs with
different software, operating systems, and locale settings.  Raw files can
contain a mixture of encoding artefacts that break downstream tools.  Run
`fix_fasta_encoding.py` on the raw FASTA **before** any other pipeline step:

```bash
fix_fasta_encoding.py spikenuc1207.fasta
# writes clean UTF-8 to spikenuc1207.fasta
# renames original to spikenuc1207.fasta.orig
```

It normalises three encoding forms in a single streaming pass:

| Artefact | Example bytes | Example glyph | Fix |
|---|---|---|---|
| Literal `\uXXXX` escapes | `\u00e9` (6 ASCII chars) | `é` | `_unescape_unicode()` |
| Raw Latin-1 bytes | `0xe9` (1 byte) | `é` | surrogate-escape decode |
| C0 control characters | `\u0003` (6 ASCII chars) | *ETX* | unescape → strip |

### The `\u0003` issue and `filterbyname.sh` v39.62

GISAID record **EPI_ISL_2016759** (Belgium, 2021-04-12) has a literal
`\u0003` (backslash-u-0-0-0-3 as six ASCII bytes) embedded in its header:

```
>Spike|...|EPI_ISL_2016759|...|'[Run20210508\u0003template bulk upload 20210508.xlsx]...|Belgium
ATGTTTGTTTTTCTT...
```

Verified with `od -c` on the raw file:

```
0000160   2   1   0   5   0   8   \   u   0   0   0   3   t   e   m   p
```

The display shows `\ u 0 0 0 3` — **six separate ASCII characters** — not
`\003` (the single ETX byte). The header and the sequence are on separate
lines; the FASTA format is well-formed.

**Why `filterbyname.sh` (BBTools ≥ v39.06) mishandles this:**

`filterbyname.sh` is a Java/Kotlin tool.  Java resolves `\uXXXX` escape
sequences **natively during string parsing**: when Java reads the six
characters `\u0003` from the file, its internal string layer silently
converts them to the actual Unicode codepoint U+0003 (ETX = End of Text).
filterbyname.sh then uses ETX as a field or record delimiter, producing
three failure modes:

1. **Header truncation** — filterbyname.sh registers a truncated record
   name (everything before the ETX: `>...|'[Run20210508`).  Depending on
   how filterbyname.sh matches names, the accession `EPI_ISL_2016759` may
   never be found in the names list.

2. **`ignorejunk=t` bypass** — if filterbyname.sh classifies the record
   as "junk" due to the ETX codepoint, `ignorejunk=t` forwards it to the
   output stream **without consulting the names list** — meaning the record
   passes through a junk-removal step it was supposed to be filtered by.

3. **Malformed output FASTA** — when filterbyname.sh writes the extracted
   record, it emits a real LF (`\n`) where the ETX was, splitting the
   header.  The tail of the original header (`template bulk upload...|Belgium`)
   is concatenated directly onto the nucleotide sequence with no newline
   separator.  The SHA-256 of this mangled "sequence" does not match the
   expected sha256.  The `summarize_fasta_pipeline.py` audit then reports:
   - **1 spurious *Novel sha256*** (the mangled sequence is a "new" hash)
   - **63 discarded entries instead of 64** (the real entry is miscounted)

**The fix applied by `fix_fasta_encoding.py`:**

```
Step 3 (_unescape_unicode):  \u0003  →  \x03   (real ETX byte, 1 byte)
Step 4 (_strip_c0_controls): \x03   →  ""      (stripped)
Result:  >...|'[Run20210508template bulk upload 20210508.xlsx]...|Belgium
```

The clean header has no control characters, filterbyname.sh reads it as a
plain ASCII string, and all three failure modes are eliminated.

### Backup file exclusion in `summarize_fasta_pipeline.py`

`fix_fasta_encoding.py` renames the original to `.fasta.orig` before writing
the cleaned file.  `summarize_fasta_pipeline.py` scans for all FASTA variants
including `.fasta.orig`.  Without a guard, **both** the clean `spikenuc1207.fasta`
and the original `spikenuc1207.fasta.orig` would appear in the file list with
the same stem, causing:

- Spurious duplicate entries in the pipeline table
- Ambiguous parent→child resolution (two roots with identical stems)
- Potential use of the *unfixed* `.fasta.orig` as the `root_path` for
  `_extract_discarded_to_fasta`, re-triggering the filterbyname.sh
  `\u0003` failure

Since commit `29c48fb`, the scanner silently skips any `.fasta.orig` /
`.fasta.ori` / `.fasta.old` file when a corresponding `.fasta` is present,
printing a diagnostic to stderr:

```
  Skipping backup 'spikenuc1207.fasta.orig' (clean 'spikenuc1207.fasta' present).
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
split_fasta_entries_by_lengths \
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

---

## Performance Benchmarks (GISAID 17M Dataset)

Following native Python byte-parsing refactors (dropping external `grep` pipelines), the `summarize_fasta_pipeline.py` script yields extreme high-throughput metrics on the full GISAID tracking sets.

**Dataset Footprint (April 2026):**
- 17,039,211 raw sequences (`spikenuc1207.fasta` / 63 GB)
- 4,723,081 distinct sequence combinations deduplicated (`spikenuc1207.no_junk.counts.fasta` / 17 GB)
- Cumulative `NNNNx` sequence expansion verified: 17,039,148 
- Mathematical exact-match tracking confirmed lossless (63 explicitly junked sequences filtered).

**Hardware Metrics (NVMe SSD with `--jobs 5`)**
Based on profiling of the full pipeline execution natively avoiding legacy shell forks:
* `Phase 0: Identity Check`: **2m 05s**, 0.1 GB Peak RAM (CPU Avg: 452%)
* `Phase 1: Validation & Accounting`: **1h 10m 00s**, 4.4 GB Peak RAM (CPU Avg: 320%)
* `Phase 2: SHA256 Verification`: **24m 33s**, 11.9 GB Peak RAM (CPU Avg: 304%)
* `Phase 3: Discard Statistics`: **17m 08s**, 16.8 GB Peak RAM (CPU Avg: 687%)
* `Phase 4: Extract Discarded Records`: **8m 30s**, 18.4 GB Peak RAM (CPU Avg: 324%)

**Total Time**: ~2 hours. 
**(Note: Following the integration of `--use-nnnx-counts` in the upcoming `processing5.sh` runs, Phase 1 times are projected to further drop radically as `O(N)` heavy sequence regex counts are dynamically bypassed).**

The native integration of `concurrent.futures.ThreadPoolExecutor` directly passing I/O buffers natively into C-space allows spanning multiple cores effectively through Python's GIL. Peak RAM reached 18.4GB across Phase 4 mapping operations, while iterating dynamically over **210+ GB of raw sequence data** efficiently in streaming chunks.
