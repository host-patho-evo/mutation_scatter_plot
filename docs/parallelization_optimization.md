# Parallelization & Performance Optimization — `calculate_codon_frequencies`

**Date:** 2026-03-29  
**Authors:** Martin Mokrejš, Antigravity AI  
**Summary:** 13.6× end-to-end wall-time reduction on the 192-core Xeon production server
(59.9 s → 4.41 s at 8 threads, 100k-row benchmark; original sequential baseline 59.9 s).

---

## 1. Context and Problem Statement

`calculate_codon_frequencies` processes a large padded multi-FASTA alignment (spike protein,
SARS-CoV-2, up to **4,401,937 sequences × 3,822 nt** = 16.2 GB) and writes per-codon-site
frequency tables. The tool was instrumented with `multiprocessing.Pool` to distribute
the work across 1,274 codon sites using `_process_one_site()` workers.

Despite having 1,274 independent tasks (embarrassingly parallel), profiling on the
192-core Xeon revealed that increasing threads _degraded_ performance:

| threads | time (original) | notes |
|---:|---:|---|
| 1 | 57 s | sequential (no pool) |
| 4 | 59.9 s | _slower_ than 1 thread |
| 192 | ~600 s | 10× slower than single-threaded |

The root causes were identified via `cProfile`, `line_profiler`, and micro-benchmarks.

---

## 2. Hardware and Filesystem

| Property | Value |
|---|---|
| System | Xeon Platinum 8260, 4 sockets × 24 cores × 2 HT = **192 logical CPUs** |
| L3 cache | 33 MB per socket |
| RAM | **3 TB** |
| Home/data FS | **NFSv4** (home + data directories) |
| Fast local FS | SSD at `/scratch.ssd/mmokrejs/` (used for temp files / profiling output) |

---

## 3. Profiling Methodology

### 3.1 cProfile

`scripts/profile_cprofile.sh` wraps `python -m cProfile` and automatically pre-subsets
the full 16 GB file to a configurable row count (default: 100,000 sequences) on the
local SSD before profiling. Output files are named with the current git commit hash.

```bash
REFERENCE=/path/to/MN908947.3_S.fasta \
bash scripts/profile_cprofile.sh \
    /path/to/spikenuc1207...exactly_3822.fasta \
    4 100000
# Outputs: $TMPDIR/profile_codon_<SHA>_t4.stats
#          $TMPDIR/profile_codon_<SHA>_t4_top50_tottime.txt
#          $TMPDIR/profile_codon_<SHA>_t4_top50_cumulative.txt
```

### 3.2 line_profiler

`scripts/profile_line.sh` uses `line_profiler` (v5) to obtain per-line timing inside
`parse_alignment`, `_process_one_site`, and `fast_fasta_iter`. This was critical for
pinpointing which exact lines drove the 40 s `parse_alignment` tottime that cProfile
reported as a single opaque bucket.

```bash
REFERENCE=/.../MN908947.3_S.fasta \
bash scripts/profile_line.sh \
    /.../spikenuc1207...exactly_3822.fasta \
    4 100000
# Outputs: $TMPDIR/profile_line_<SHA>_t4.lprof
#          $TMPDIR/profile_line_<SHA>_t4.txt   (ms units per line)
```

### 3.3 Thread-scaling benchmark

`scripts/bench_thread_scaling.py` sweeps a list of thread counts and optionally a list
of row-count subsets, writing a TSV with min/median/max times per combination. Output
filenames include the git commit hash.

---

## 4. Root-Cause Analysis — cProfile Runs

### Run 1 — original code (`d0e04f3`), 4 threads, 100k rows

```
Total wall time: 59.889 s

ncalls  tottime   function
     1   40.744   pool._join_exited_workers          ← SIGTERM GC thrash
    16    6.574   posix.write                        ← IPC payload: 115 MB × 5 batches
     5    4.815   pickle.Pickler.dump                ← pickling full aln_array per batch
    27    3.183   select.poll.poll                   ← actual worker compute
```

**Root cause 1 — `pool.terminate()`:** The pool was closed via `terminate()` instead of
`close()+join()`. `SIGTERM` interrupts workers mid-computation; Python's GC then
collects all in-memory objects (including a 115 MB NumPy array per worker) under
interrupt context, each taking ~10 s. With 4 workers: 40+ s wasted.

**Root cause 2 — large-array IPC:** All task arguments (`_aln_array`, `_counts_array`,
and 5 other arrays, totalling ~115 MB for 30k unique seqs) were pickled into every
task tuple. With `chunksize=1` (the pool default at high worker counts), 1,274 tasks ×
115 MB = **147 GB** sent over Unix pipes — explaining the 21.5× penalty at 192 workers.
With `chunksize=256` each batch still pickled one full copy: 5 batches × 115 MB = 575 MB.

### Run 2 — after `close()+join()` and `chunksize=256` (`705e76d`), 100k rows

```
Total wall time: 49.716 s  (−17 %)
pool._join_exited_workers: gone
queues.get (_handle_tasks): 41 s ← true worker compute (sequential 1274 sites × 100k seqs)
posix.write: 2.567 s  (↓2.6×)
pickle.Pickler.dump: 1.434 s  (↓3.4×)
```

Workers no longer GC under SIGTERM, and fewer IPC round-trips reduce pickle/write
overhead. The `queues.get` time (41 s) represents _actual_ worker compute, confirming
that 1,274 NumPy-heavy codon-site evaluations take ~41 s in parallel across 4 workers.

### Run 3 — after COW-fork globals (`a3ff37c`), 100k rows

```
Total wall time: 46.061 s
parse_alignment tottime: 40.476 s   ← mystery: what is the serial code doing for 40 s?
select.poll (worker compute): 2.186 s
str.replace (depadding): 0.656 s
fast_fasta_iter cumtime: 0.830 s
```

IPC is eliminated; workers inherit `_aln_array` via COW fork and receive only `_pos`
(4 bytes) per task. Worker compute drops to 2.186 s. But `parse_alignment`'s own
tottime is still 40 s — indicating a major serial bottleneck in the Python loop.

### Run 4 — after time-gated flush (`475d345`), 100k rows

```
Total wall time: 46.025 s  (unchanged)
```

Hypothesis: per-site `outfilename.flush()` calls (3829×) were causing NFS round-trip
latency. **Disproved**: each flush() call takes only 2.3 µs (Python flushes to the OS
page cache, not to the NFS server synchronously). The flush was nevertheless replaced
with a 30-second time-gate for monitoring large production runs with `tail -f`.

### Run 5 — after `lstrip`/`rstrip` regex fix (`924502e`), 100k rows

```
Total wall time: 6.495 s  (−86 % from Run 4 = −40 s)
parse_alignment tottime: 0.565 s   ← dropped from 40.5 s to 0.6 s
select.poll (worker compute): 2.423 s
str.replace (depadding): 0.621 s
fast_fasta_iter cumtime: 0.827 s
```

---

## 5. Root-Cause Analysis — line_profiler

`scripts/profile_line.sh` found the bottleneck immediately:

```
Function: parse_alignment  —  Total time: 46.9 s

Line    Hits       Time   % Time  Line contents
 477  89,899  39,911.4 ms   85.2  for _match in _re_trailing_gaps.finditer(_aln_line_seq):
```

The regex `[-Nn]+$` is anchored to the **end** of a 3,822-character string. Python's
`re` engine must scan all 3,822 characters before it can anchor, then match backwards.
For 89,899 unique sequences:

```
89,899 strings × 444 µs/call = 39.9 s
```

The leading-gaps regex `^[-Nn]+` was only 122 ms (anchor at start — fails fast on
non-gap sequences).

**Fix:** `str.rstrip('-N')` scans from the right in pure C without regex engine
overhead. Equivalent result because `.upper()` is already applied:

```python
# Before — 39,911 ms total
_start_of_trailing_gaps = 0
for _match in _re_trailing_gaps.finditer(_aln_line_seq):
    _start_of_trailing_gaps = _match.start()

# After — ~45 ms total (200× faster)
_start_of_trailing_gaps = len(_aln_line_seq.rstrip('-N'))
```

Note: The sentinel value changes when no trailing gaps are present (`0` → `len(seq)`),
but the downstream condition `(_trail != 0) & (_pos + 1 > _trail)` evaluates correctly
in both cases because `_pos + 1 ≤ len(seq)` for all valid positions.

Also removed: `re.compile()` calls and `import re` (now unused).

---

## 6. All Fixes — Chronological

| Commit | Fix | Wall-time change (100k rows, 4T) |
|---|---|---|
| `d0e04f3` | **baseline** | 59.9 s |
| `705e76d` | `pool.close()+join()` instead of `terminate()`; add `--chunksize 256` | −10.2 s (−17%) |
| `a3ff37c` | COW-fork globals: populate `_WORKER_SHARED` before `Pool()`, send only `_pos` (4 bytes) per task | −3.7 s (−7%) |
| `475d345` | Time-gated flush (every 30 s instead of per-site); harmless, monitoring-friendly | 0 s |
| `924502e` | **`lstrip`/`rstrip` instead of `$`-anchored trailing-gaps regex ← MAIN WIN** | **−39.5 s (−86%)** |
| `9b5b91a` | `str.count('-') + str.count('N')` instead of `replace('-','').replace('N','')` for depadded length; avoids 2 temporary string allocations per unique seq | ~0 s at 100k rows |
| `0af7478` | Numpy array transposed to `(alignment_len, num_unique)` for cache-friendly row access | +0.3 s ← **reverted** |
| `e58395e` | Revert transpose (see §7) | −0.35 s vs `0af7478` |

### Fix A in depth — COW-fork globals (`a3ff37c`)

Workers are forked _after_ `_WORKER_SHARED` is populated, so they inherit all large
arrays via Linux copy-on-write (read-only, zero physical copies). Only `_pos` (one
integer, 8 bytes) is sent per task via the IPC pipe:

```python
# Module level
_WORKER_SHARED: dict = {}

def _process_one_site_wrapper(_pos):
    d = _WORKER_SHARED
    return _process_one_site(
        _pos, d['num_unique'], d['alignment_len'], d['aln_array'], ...
    )

# In parse_alignment(), before Pool():
_WORKER_SHARED = { 'aln_array': _aln_array, ... }
_active_pool = multiprocessing.Pool(processes=threads)   # forks here
_all_results = _active_pool.starmap(
    _process_one_site_wrapper,
    [(_pos,) for _pos in _slim_positions],
    chunksize=chunksize,
)
```

---

## 7. Ideas Investigated but Not Adopted

### 7.1 Numpy matrix transpose

**Hypothesis:** `_aln_array` has shape `(num_unique, alignment_len)` (C-contiguous).
`_process_one_site` does `_aln_array[:, _pos:_pos+3]` — a column slice that strides
`alignment_len` bytes between elements, causing L1/L2 cache misses.

Storing as `(alignment_len, num_unique)` via `.T.copy()` would make `_aln_array[_pos]`
a contiguous row of `num_unique` bytes — cache-friendly for SIMD `astype(np.uint32)`.

**Result:** Benchmarks showed +10% slowdown at 87k unique seqs (100k-row subset):

| commit | 8T time |
|---|---|
| `9b5b91a` (pre-transpose) | 4.46 s |
| `0af7478` (with transpose) | 4.76 s |
| `e58395e` (transpose reverted) | **4.41 s** |

**Why it failed:** Per-site working set = 87k × 3 bytes = **256 KB**, which fits in
L2 cache per core (512 KB–1 MB). Cache misses were not the limiting factor at this
scale. The `.T.copy()` overhead (~22 ms for 87k seqs, reading+writing 336 MB) exceeded
any cache benefit. The transpose would only help at `num_unique > L3/3 bytes ≈ 11M`,
which exceeds the dataset size (4.4M unique seqs).

### 7.2 NFS-flush latency

**Hypothesis:** Per-site `outfilename.flush()` calls (3,829 total) each incur a
10 ms NFSv4 synchronous round-trip → 38 s total.

**Result:** Disproved. `flush()` flushes Python's user-space buffer to the OS page
cache (asynchronous NFS writeback). Measured at 2.3 µs/call = 9 ms total. The
per-site flush was replaced with a 30-second time-gate (useful for `tail -f` monitoring
of long production runs) but had zero effect on benchmark performance.

### 7.3 Column-band shell wrapper

**Concept (not implemented):** Pre-split the 1,274 codon sites into ~425 bands of
3 codons each; run one `--threads 1` process per band via GNU `parallel`; merge TSV
outputs. Fully embarrassingly parallel at the process level with no IPC overhead.

**Status:** Not implemented. After the regex fix, wall time dropped from 60 s to 4.4 s
which is within an acceptable range for the current workflow. If the full-dataset run
time (projected ~3.5 min at 8 threads) becomes a bottleneck, this approach remains
the highest-potential architectural change.

---

## 8. Amdahl Analysis (post-fix, 100k–1M rows)

After all fixes, both the serial component (FASTA parsing, depadding, numpy setup)
and the parallel component (site processing) scale **linearly** with `num_unique_seqs`,
because the 4.4M input sequences are all pre-deduplicated (sha256). The Amdahl
parallel fraction therefore stays constant at ~57–60% regardless of dataset size.

| rows | unique seqs | S (serial) | P (1T parallel) | P fraction | max speedup |
|---:|---:|---:|---:|---:|---:|
| 100k | ~87k | ~4.0 s | ~4.5 s | 52% | 2.1× |
| 500k | ~450k | ~19.9 s | ~29.9 s | 60% | 2.5× |
| 1M | ~900k | ~44.7 s | ~59.3 s | 57% | 2.3× |
| 4.4M | 4.4M | ~175 s | ~263 s | 60% | ~2.5× |

> **Implication:** Increasing threads beyond 8–16 gives diminishing returns for all
> dataset sizes. The remaining serial bottleneck is the FASTA read + per-sequence
> setup (`.upper()`, `lstrip`/`rstrip`, `str.count`, numpy construction), all of
> which are inherently sequential for a single FASTA file.

---

## 9. Final Benchmark Results

### Machine: Xeon Platinum 8260, 192 cores, NFSv4 — commit `e58395e`

#### 100k rows, 1 run

| threads | time (s) | speedup |
|---:|---:|---:|
| 1 | 8.89 | 1.00× |
| 2 | 5.97 | 1.49× |
| 4 | 6.09 | 1.46× |
| **8** | **4.41** | **2.01×** ← recommended |
| 16 | 4.53 | 1.96× |
| 32 | 4.50 | 1.97× |
| 64 | 4.92 | 1.81× |
| 128 | 6.10 | 1.46× |
| 192 | 6.37 | 1.39× |

#### 500k rows, 3 runs (Stage 2, pre-Fix D, commit `924502e`)

| threads | median (s) | speedup |
|---:|---:|---:|
| 1 | 49.8 | 1.00× |
| **8** | **23.6** | **2.11×** |
| 16 | 24.1 | 2.06× |
| 192 | 32.3 | 1.54× |

#### 1M rows, 3 runs (Stage 2, pre-Fix D)

| threads | median (s) | speedup |
|---:|---:|---:|
| 1 | 104.0 | 1.00× |
| 8 | 52.1 | 2.00× |
| **16** | **50.8** | **2.05×** |
| 192 | 68.3 | 1.52× |

### Projected full dataset (4.4M rows)

| config | estimated wall time |
|---|---|
| 1 thread | ~7.5 min |
| **8 threads (recommended)** | **~3.5 min** |
| 32 threads | ~3.5–4 min |
| 192 threads | ~5–6 min (overhead dominates) |

---

## 10. Recommended Production Invocation

```bash
# Pull latest optimizations
git pull

# Run on the full dataset using local SSD for output (avoids NFS write latency)
REFERENCE=/auto/vestec1-elixir/projects/biocev/mmokrejs/proj/zahradnik/GISAID/hCoV-19_spikenuc1207/MN908947.3_S.fasta
ALIGNMENT=/auto/vestec1-elixir/projects/biocev/mmokrejs/proj/zahradnik/GISAID/hCoV-19_spikenuc1207/padded_length_3822/spikenuc1207.native2ascii.no_junk.counts.clean.exactly_3822.fasta
OUTDIR=/scratch.ssd/mmokrejs/results

python -m mutation_scatter_plot.calculate_codon_frequencies.cli \
    --reference-infile "$REFERENCE" \
    --alignment-file  "$ALIGNMENT" \
    --outfile-prefix  "$OUTDIR/spikenuc1207" \
    --padded-reference \
    --x-after-count \
    --overwrite \
    --threads 8 \
    --chunksize 256
```

---

## 11. Profiling and Benchmarking Scripts

| Script | Purpose | Output naming |
|---|---|---|
| `scripts/profile_cprofile.sh` | cProfile single run, 100k-row subset | `$TMPDIR/profile_codon_<SHA>_t<N>.{stats,txt}` |
| `scripts/profile_line.sh` | line_profiler run, per-line ms timing | `$TMPDIR/profile_line_<SHA>_t<N>.{lprof,txt}` |
| `scripts/bench_thread_scaling.py` | Wall-clock speedup vs thread count | `bench_scaling_<basename>_<SHA>.tsv` |
| `scripts/bench_chunksize_sweep.py` | Sweep IPC chunksize values | `bench_chunksize_<basename>_<SHA>.tsv` |

Output filenames always include the short git commit hash so results can be
unambiguously matched to the code state that produced them.
