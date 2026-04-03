# Translation Logic in `mutation_scatter_plot`

This document explains how nucleotide sequences are translated to amino acids
throughout the pipeline, why standard Biopython cannot be used directly on
alignment-padded FASTA files, and how the `scripts/translate.py` command-line
tool should be used.

## Table of Contents

1. [Background — the Biopython limitation](#background)
2. [The per-codon splitting workaround](#workaround)
3. [Codon resolution semantics](#semantics)
4. [Translation table summary (all 4096 IUPAC+gap codons)](#table)
5. [Old vs new: the NNN regression and its fix](#old-vs-new)
6. [Using `scripts/translate.py`](#cli)
7. [Selecting an NCBI genetic code table](#tables)
8. [Implementation notes for `calculate_codon_frequencies`](#codon-freq)

---

## 1. Background — the Biopython limitation {#background}

Standard Biopython `Bio.Seq.translate()` is designed for *clean* nucleotide
sequences.  Multi-FASTA alignment files produced by tools such as `mafft`,
`muscle`, or `clustalw` contain **gap characters** (`-`) used to preserve
reading-frame alignment across insertions and deletions.  When Biopython
is given such a sequence it fails in two ways:

### 1a. Whole-sequence gap handling

Calling `translate()` on a padded sequence requires the explicit `gap='-'`
keyword argument so that Biopython recognises full-gap codons (`---`):

```python
from Bio.Seq import translate

translate("TCA---TCG")        # raises TranslationError — unknown codon '---'
translate("TCA---TCG", gap='-')  # OK → "S-S"
```

### 1b. Partial-gap codons (the hard case)

Even with `gap='-'` Biopython only handles the *full-gap* codon `---`.  If a
gap character falls within a codon (a "partial-gap codon") — which is common
in gapped NGS reads aligned to a reference — Biopython silently returns `'X'`
rather than exploiting the *known* nucleotides:

```python
translate("TC-", gap='-')  # → 'X'  (wrong — should be 'S', see below)
```

This is the bug described in
[Biopython issue #5036](https://github.com/biopython/biopython/issues/5036)
and the related
[PR #4992](https://github.com/biopython/biopython/pull/4992), which the
Biopython maintainers declined to merge.

### 1c. Whole-sequence IUPAC translation is unreliable for alignments

Calling `translate()` on a sequence longer than 3 nucleotides works for clean
sequences, but fails for aligned sequences because the gap handling applies
only at call-level, not per-codon.  The ambiguity resolution (`N`, `R`, `Y`,
…) is per-codon — a mixture of standard and IUPAC nucleotides in the same
sequence can only be handled correctly if each codon is evaluated
independently.

### 1d. Incomplete trailing codons

NGS reads that have been trimmed or partially aligned may end with 1 or 2
nucleotides that do not form a complete codon.  Calling `translate()` on such
a fragment causes:

```
BiopythonWarning: Partial codon, len(sequence) not a multiple of three.
Explicitly trim the sequence or add trailing N before translation.
This may become an error in future.
```

This warning becomes a hard error in recent Biopython development branches
and must be avoided in production pipelines.

---

## 2. The per-codon splitting workaround {#workaround}

The solution is to split the sequence into 3-nt chunks **before** calling
Biopython, handle each codon individually, and join the results.  This is
implemented in `mutation_scatter_plot.alt_translate()`:

```python
@functools.lru_cache(maxsize=128)
def alt_translate(seq, table=1):
    result = []
    for i in range(0, len(seq), 3):
        codon = seq[i:i + 3]
        if not codon or len(codon) < 3:   # skip incomplete trailing codon
            continue
        if '-' in codon and codon != '---':
            codon = codon.replace('-', 'N')  # TC- → TCN
        result.append(translate(codon, table=table, gap='-'))
    return ''.join(result)
```

Key design decisions:

| Decision | Reason |
|---|---|
| **Per-codon loop** | Lets Biopython resolve IUPAC ambiguity independently for each codon |
| **`'-'` → `'N'` substitution** | An alignment gap means *unknown nucleotide*, semantically identical to `N` |
| **Skip if `len(codon) < 3`** | Avoids the `BiopythonWarning` for incomplete trailing codons |
| **`lru_cache`** | Most alignment files contain thousands of identical codons; caching gives ~10× speedup |
| **`table` parameter in cache key** | Different NCBI genetic codes are cached independently |

The same logic is used in `scripts/translate.py` via `_translate_seq()`.
The difference is that `_translate_seq()` uses a pre-built 64-entry dict for
the *fast path* (standard ATGC codons) and only calls Biopython for the
*fallback path* (IUPAC or partial-gap codons):

```python
def _translate_seq(seq, codon_table, ignore_gaps, table_id=1):
    for i in range(0, length - length % 3, 3):  # skip trailing nts cleanly
        codon = seq[i:i + 3]
        aa = codon_table.get(codon)    # O(1) dict lookup — no Biopython call
        if aa is None:                 # IUPAC or partial-gap: fall through
            codon_n = codon.replace('-', 'N')
            aa = _bio_translate(codon_n, table=table_id, gap='-')
        result.append(aa)
```

This two-path design means that for typical virus spike protein datasets
(≥95% standard ATGC codons) Biopython is called for fewer than 5% of codons.

---

## 3. Codon resolution semantics {#semantics}

### 3a. Standard codons (fast dict path, O(1))

All 64 ATGC codons (upper and lower case) are pre-built into a dict at tool
startup via `Bio.Data.CodonTable.unambiguous_dna_by_id[table_id]`.  No
Biopython call is made at translation time.

### 3b. Full-gap codon

| Input | Output | Meaning |
|-------|--------|---------|
| `---` | `-`    | Entire codon is a deletion relative to the reference |

### 3c. Partial-gap codons — the interesting case

An alignment gap in a non-gap-only codon position means the nucleotide at that
position is *unknown* — exactly the same semantics as the IUPAC code `N`.  The
substitution `'-' → 'N'` lets Biopython's per-codon ambiguity resolution
determine whether the known nucleotides uniquely constrain the amino acid.

**Complete list of partial-gap codons that resolve to a specific amino acid
(standard genetic code, table 1):**

| Codon | Intermediate | Amino acid | Why |
|-------|--------------|------------|-----|
| `TC-` | `TCN`        | `S` (Serine)   | All four TCA/TCC/TCG/TCT encode Serine — third position is four-fold degenerate |
| `AC-` | `ACN`        | `T` (Threonine)| All four ACA/ACC/ACG/ACT encode Threonine |
| `CC-` | `CCN`        | `P` (Proline)  | All four CCA/CCC/CCG/CCT encode Proline |
| `CG-` | `CGN`        | `R` (Arginine) | All four CGA/CGC/CGG/CGT encode Arginine |
| `CT-` | `CTN`        | `L` (Leucine)  | All four CTA/CTC/CTG/CTT encode Leucine |
| `GC-` | `GCN`        | `A` (Alanine)  | All four GCA/GCC/GCG/GCT encode Alanine |
| `GG-` | `GGN`        | `G` (Glycine)  | All four GGA/GGC/GGG/GGT encode Glycine |
| `GT-` | `GTN`        | `V` (Valine)   | All four GTA/GTC/GTG/GTT encode Valine |

**All other partial-gap codons resolve to `X` (ambiguous amino acid):**

| Codon example | Intermediate | Why |
|---------------|--------------|-----|
| `AT-`         | `ATN`        | ATA/ATC/ATT=Ile, ATG=Met — ambiguous |
| `A--`         | `ANN`        | Too many unknowns |
| `-TC`         | `NTC`        | First position unknown — covers multiple codons |
| `---`         | —            | Full-gap → `-` (handled separately) |

### 3d. Pure IUPAC ambiguity codons (no gap)

Same logic applies — Biopython resolves each IUPAC code to its set of possible
nucleotides and returns `X` if the amino acid is ambiguous across those
combinations.

| Codon | Amino acid | Explanation |
|-------|------------|-------------|
| `TCN` | `S`  | All possible third bases give Serine |
| `CTN` | `L`  | All possible third bases give Leucine |
| `GGN` | `G`  | Four-fold degenerate → Glycine |
| `GCN` | `A`  | Four-fold degenerate → Alanine |
| `GTN` | `V`  | Four-fold degenerate → Valine |
| `CCN` | `P`  | Four-fold degenerate → Proline |
| `ACN` | `T`  | Four-fold degenerate → Threonine |
| `CGN` | `R`  | Four-fold degenerate → Arginine |
| `ATN` | `X`  | ATA=I, ATC=I, ATT=I, ATG=M — Ile/Met ambiguous |
| `AGN` | `X`  | AGA=R, AGC=S, AGG=R, AGT=S — Arg/Ser ambiguous |
| `NNN` | `X`  | Completely ambiguous |

### 3e. Incomplete trailing codons (respect-alignment mode)

In alignment-respecting mode, a 1-nt or 2-nt fragment at the end of a
sequence represents a real codon position in the alignment where the read was
trimmed.  The amino acid is unknown but the position is not absent — so `'X'`
is returned, *not* an empty string.  Calling Biopython on a sub-3-nt fragment
causes a `BiopythonWarning` and may raise in future versions; this was one of
the original motivations for `alt_translate()` to exist.

| Input    | Output | Note |
|----------|--------|------|
| `TCA`    | `S`    | Complete codon |
| `TCATC`  | `SX`   | TCA=S + TC (incomplete trailing codon) → X |
| `TC`     | `X`    | Incomplete — position exists, AA unknown |
| `T`      | `X`    | Incomplete — position exists, AA unknown |
| `T-`     | `X`    | Incomplete (2 nt total) — position exists, AA unknown |
| `""`     | `""`   | Empty input — no positions at all |

---

## 4. Translation table summary (all 4096 IUPAC+gap codons) {#table}

All 4096 codons from `{A,C,G,T,R,Y,S,W,K,M,B,D,H,V,N,-}³` were tested in
`tests/test_translation.py::TestExhaustiveIUPACGapCodons`.

Distribution of results (NCBI table 1, standard genetic code):

| Result | Count | Meaning |
|--------|------:|---------|
| `X`    | 3896  | Ambiguous — multiple amino acids possible |
| `-`    |     1 | Full-gap codon `---` |
| `A`    |    16 | Alanine (four-fold degenerate family + IUPAC variants) |
| `C`    |     3 | Cysteine |
| `D`    |     3 | Aspartate |
| `E`    |     3 | Glutamate |
| `F`    |     3 | Phenylalanine |
| `G`    |    16 | Glycine |
| `H`    |     3 | Histidine |
| `I`    |     7 | Isoleucine |
| `J`    |     9 | Ile or Leu (IUPAC ambiguity) |
| `K`    |     3 | Lysine |
| `L`    |    22 | Leucine |
| `M`    |     1 | Methionine (ATG only) |
| `N`    |     3 | Asparagine |
| `P`    |    16 | Proline |
| `Q`    |     3 | Glutamine |
| `R`    |    22 | Arginine |
| `S`    |    19 | Serine |
| `T`    |    16 | Threonine |
| `V`    |    16 | Valine |
| `W`    |     1 | Tryptophan (TGG only) |
| `Y`    |     3 | Tyrosine |
| `Z`    |     3 | Glu or Gln (IUPAC ambiguity) |
| `B`    |     3 | Asp or Asn (IUPAC ambiguity) |
| `*`    |     5 | Stop codon |
| **Total** | **4096** | |

> **Note**: `J`, `B`, `Z` are IUPAC protein ambiguity codes returned by
> Biopython when the nucleotide ambiguity code does not uniquely determine
> the amino acid but constrains it to a known pair.

---

## 5. Original PR intent vs the current `'-'→'N'` extension {#old-vs-new}

> **Reference**: [Biopython PR #4992](https://github.com/biopython/biopython/pull/4992) —
> *Add more ways to translate nucleotide sequence with dashes* (authored by
> Martin Mokrejš, the author of this pipeline).

### 5a. What the PR proposed

The PR proposed two new optional arguments for `Bio.Seq.translate()`:

| Mode | Behaviour |
|------|-----------|
| `respect_alignment=True` | `---` → `-`; **any partial-gap codon** (e.g. `AT-`, `TC-`) → `X` |
| `ignore_gaps=True` | Internally does `.replace('-', '')` then translates (re-frames) |

The PR was **not accepted** by the Biopython maintainers, who suggested
workarounds such as `.replace("-", "").translate()` or custom codon tables.
These workarounds are all wrong for alignment-aware translation (see §1).

The PR explicitly noted a desirable **future extension** (not implemented in
the patch):
> *"If you want to make biopython clever, there are a few cases when one could
> safely translate CTN into L, CGN into R, TCN into S, GGN into G, GCN into A,
> GTN into V, ACN into T and CCN into P. I do not add this functionality in
> my patch but it would be handy extension…"*

### 5b. What our implementation does (the extension)

Our `'-'→'N'` per-gap substitution is precisely that future extension —
applied per codon before calling Biopython's IUPAC-aware translation.

| Codon | PR intent | Our result | Difference |
|-------|-----------|-----------|------------|
| `TC-` | `X`       | **`S`**   | TCN → all Serine (four-fold degenerate) |
| `AC-` | `X`       | **`T`**   | ACN → all Threonine |
| `CC-` | `X`       | **`P`**   | CCN → all Proline |
| `CG-` | `X`       | **`R`**   | CGN → all Arginine |
| `CT-` | `X`       | **`L`**   | CTN → all Leucine |
| `GC-` | `X`       | **`A`**   | GCN → all Alanine |
| `GG-` | `X`       | **`G`**   | GGN → all Glycine |
| `GT-` | `X`       | **`V`**   | GTN → all Valine |
| `AT-` | `X`       | **`X`**   | ATN ambiguous (Ile/Met) — same as PR |
| `AA-` | `X`       | **`X`**   | AAN ambiguous (Lys/Asn) — same as PR |
| `---` | `-`       | **`-`**   | Full-gap codon — same as PR |

The 8 cases where we differ from the raw PR design are exactly the
four-fold degenerate third-position families in the standard genetic code.
All other partial-gap codons still yield `X`, consistent with the PR.

### 5c. The old NNN workaround (commit dd812e8)

Before the `-→N` approach, the code used a **whole-codon NNN substitution**:

```python
codons = ("NNN" if "-" in codon and codon != "---" else codon for codon in codons)
```

This always gives `X` for any partial-gap codon — equivalent to the original
PR's `respect_alignment=True` intent, but worse than the `-→N` extension for
the 8 four-fold degenerate cases above.

### 5d. The Biopython string/Seq inconsistency

An important subtlety discovered in the PR thread: `translate("---")` (string
argument) raises `TranslationError: Codon '---' is invalid` on some Biopython
versions, while `Seq("---").translate()` returns `Seq('-')`.  This is why
`alt_translate` calls `translate(codon, gap='-')` per-codon rather than on
the whole sequence — the per-codon call with `gap='-'` handles `---` → `-`
reliably across versions.

### 5e. Why global substitution strategies are phase-wrong

The PR thread shows why `seq.replace("---", "???").replace("-", "N").replace("???", "---")` fails: `---` inside an alignment may not be in-phase with the codon frame.  Consider `AA---A`:

```
Codon frame: AA- | --A
Per-codon:   AA- → AAN → X
             --A → NNA → X
Result: XX  ✓
```

A global `---`→placeholder approach treats `---` as a unit, but the actual
codon boundaries may split it differently.  Per-codon processing is the only
correct approach.

---

## 6. Using `scripts/translate.py` {#cli}

`translate.py` is a **high-performance standalone CLI tool** for translating
nucleotide FASTA files to protein FASTA files.  It is designed to process
millions of records without loading the entire file into memory.

### Why not use `seqkit translate` or Biopython directly?

| Tool | Problem |
|------|---------|
| `seqkit translate` | Does not handle alignment gap characters in codons |
| `biopython translate(seq)` | Crashes on `---` without `gap='-'`; returns wrong `X` for partial-gap codons like `TC-`; emits warnings for incomplete trailing codons |
| `biopython translate(seq, gap='-')` | Handles `---`→`-` only; still wrong for `TC-`→`X` |
| `translate.py` | Correct per-codon semantics, streaming I/O, no memory limit |

### Basic usage

```bash
# Translate a FASTA file (reads from stdin if --infile not given):
python scripts/translate.py \
    --infile  sequences.nucleotide.fasta \
    --outfile sequences.protein.fasta

# From a pipeline (stdin → stdout):
zcat sequences.fasta.gz | python scripts/translate.py --infile - --outfile - | gzip > out.fasta.gz

# Ignore alignment gaps (re-frame gaps away before slicing):
python scripts/translate.py --infile aligned.fasta --outfile out.fasta --ignore-gaps

# Use vertebrate mitochondrial code (NCBI table 2):
python scripts/translate.py --infile mito.fasta --outfile mito_prot.fasta --translation-table 2
```

### All options

| Option | Default | Description |
|--------|---------|-------------|
| `--infile FILE` | stdin | Input FASTA file (`-` for stdin) |
| `--outfile FILE` | stdout | Output FASTA file (`-` for stdout) |
| `--infileformat` | `fasta` | Input format (passed to format detection) |
| `--outfileformat` | `fasta-2line` | Output format |
| `--ignore-gaps` | off | Strip all `-` before translation (destroys alignment) |
| `--respect-alignment` | on | Keep `-` in alignment; partial-gap codons use `'-'→'N'` |
| `--translation-table N` | `1` | NCBI genetic code table number |

### Performance

`translate.py` uses a streaming raw-bytes FASTA parser that avoids all
`SeqRecord` allocation.  The core hot path is a single `dict.get()` lookup
per codon — no Biopython overhead for the ~95% of standard ATGC codons.
Biopython is only called for IUPAC/partial-gap fallback codons.

Benchmark on a 17 M-record FASTA file (GISAID spike protein dataset):
- Throughput: ≈ 400 MB/s on a single core
- Memory: O(1) — constant regardless of file size

### Why the tool is not importable as a Python module

`translate.py` calls `myparser.parse_args()` at **module level** (a legacy of
the `optparse`-based design).  This means `import translate` triggers argument
parsing, which fails under pytest and other tools that pass their own
`sys.argv`.

The test suite works around this by temporarily patching `sys.argv` before
import:

```python
import sys
sys.argv = ['translate']
try:
    from translate import _build_codon_table, _translate_seq, parse_input
finally:
    sys.argv = _saved_argv
```

A future refactor should move `parse_args()` inside `_main()`.

---

## 7. Selecting an NCBI genetic code table {#tables}

The `--translation-table N` option (default: `1`) selects the NCBI genetic
code.  Notable alternatives:

| Table | Name | Key differences from table 1 |
|-------|------|-------------------------------|
| 1 | Standard | (baseline) |
| 2 | Vertebrate Mitochondrial | `TGA`=Trp (not stop); `ATA`=Met (not Ile); `AGA`/`AGG`=stop (not Arg) |
| 4 | Mold/Protozoan Mitochondrial | `TGA`=Trp |
| 11 | Bacterial/Archaeal | Same codons as table 1 but different initiator set |

Full list: <https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi>

The table selection propagates through the entire pipeline:

```
translate.py --translation-table 2
  └─ _build_codon_table(table_id=2)      # fast-path dict for standard codons
  └─ _translate_seq(..., table_id=2)
       └─ _bio_translate(codon_n, table=2, gap='-')   # fallback: IUPAC/gaps

calculate_codon_frequencies --translation-table 2
  └─ alt_translate(reference_seq, table=2)             # reference protein
  └─ parse_alignment(..., translation_table=2)
       └─ write_tsv_line(..., translation_table=2)
            └─ alt_translate(codon, table=2)           # each observed codon
```

---

## 8. Implementation notes for `calculate_codon_frequencies` {#codon-freq}

`calculate_codon_frequencies` uses `alt_translate()` via `write_tsv_line()`.
One important difference from the standalone `translate.py` tool:

| Context | Incomplete codon (len < 3) handling |
|---------|-------------------------------------|
| `alt_translate()` (respect-alignment) | Returns `"X"` — position exists, AA unknown |
| `write_tsv_line()` | Returns `"X"` — consistent with `alt_translate` |
| `_translate_seq(..., ignore_gaps=False)` | Silently drops trailing fragment (no `X`) |

In the codon-frequency context, a codon shorter than 3 nucleotides comes from
an NGS read that was trimmed mid-codon.  Mapping it to `X` (unknown) rather
than dropping it is the correct biological treatment — the read *does* provide
evidence of sequence at that position, just insufficient to determine the amino
acid.

Exhaustive testing of the `write_tsv_line` translation path is in
`tests/test_translation.py::TestCalculateCodonFrequenciesTranslation`.
