#!/usr/bin/env python3
"""Scan monthly frequency data to find positions with significant mutation dynamics.

**Strategy** (two-phase):

1. **Candidate detection** — parse ``*.unchanged_codons.tsv`` files and collect
   every AA position whose unchanged-codon frequency dropped below a threshold
   (default 90 %) in *at least one month* since 2020-01.  These are sites
   where the reference codon was significantly displaced at some point.

2. **Dynamics characterisation** — for each candidate, walk through the
   ``*.frequencies.tsv`` files month-by-month and build a rich profile
   tracking total mutant frequency, per-AA trajectories, dominant-mutation
   switches, sweep counts, and volatility.

   Classifications:
   - **replacement**    — dominant mutant AA changed over time (variant replacement,
     e.g. E484K → E484A → E484K)
   - **recurrence**     — frequency swept up through 20 %% at least twice
   - **fixation**       — rose to ≥90 %% and stayed there
   - **transient**      — peaked ≥20 %% then fell back below 5 %%
   - **growing**        — trending upward in recent months
   - **declining**      — was high, now dropping
   - **persistent_low** — never exceeded 20 %%
   - **variable**       — significant fluctuation not fitting other categories

The output is a sorted list of AA positions ready for
``mutation_timeline_plot --positions …``.

Usage
-----
::

    python find_dynamic_positions.py \\
        --dir /path/to/monthly_tsv_files \\
        [--unchanged-pattern '*.unchanged_codons.tsv'] \\
        [--freq-pattern '*.frequencies.tsv'] \\
        [--threshold 0.10] \\
        [--start-month 2020-01] \\
        [--end-month 2025-07] \\
        [--min-months 3] \\
        [--output-format positions|tsv|both]
"""

from __future__ import annotations

import argparse
import csv
import glob
import os
import re
import sys
from collections import defaultdict
from dataclasses import dataclass, field as dc_field


# ──────────────────────────────────────────────────────────────────────────────
# Helpers
# ──────────────────────────────────────────────────────────────────────────────

_MONTH_RE = re.compile(r'\.(\d{4}-\d{2})\.')


def _extract_month(path: str) -> str | None:
    """Return YYYY-MM from the filename, or None."""
    m = _MONTH_RE.search(os.path.basename(path))
    return m.group(1) if m else None


def _month_in_range(month: str, start: str, end: str) -> bool:
    """Check if YYYY-MM is within [start, end] inclusive (string comparison)."""
    return start <= month <= end


# ──────────────────────────────────────────────────────────────────────────────
# Phase 1: scan unchanged_codons.tsv for candidate positions
# ──────────────────────────────────────────────────────────────────────────────

def find_candidate_positions(
    directory: str,
    pattern: str,
    threshold: float,
    start_month: str,
    end_month: str,
) -> dict[int, dict]:
    """Return {position: info_dict} for positions that dropped below *threshold*.

    info_dict keys:
        - min_freq:   minimum unchanged frequency observed
        - min_month:  the month in which min_freq occurred
        - months_below: number of months below the threshold
    """
    files = sorted(glob.glob(os.path.join(directory, pattern)))
    if not files:
        print(f"Warning: no files matched '{pattern}' in {directory}",
              file=sys.stderr)
        return {}

    # position -> {min_freq, min_month, months_below}
    candidates: dict[int, dict] = {}

    for fpath in files:
        month = _extract_month(fpath)
        if month is None:
            continue
        # Skip YYYY-00 aggregate files
        if month.endswith('-00'):
            continue
        if not _month_in_range(month, start_month, end_month):
            continue

        with open(fpath, newline='') as fh:
            reader = csv.DictReader(fh, delimiter='\t')
            for row in reader:
                try:
                    pos = int(row['position'])
                    freq = float(row['frequency'])
                except (ValueError, KeyError):
                    continue

                if freq < (1.0 - threshold):
                    # This position has ≥ threshold mutant frequency
                    if pos not in candidates:
                        candidates[pos] = {
                            'min_freq': freq,
                            'min_month': month,
                            'months_below': 1,
                        }
                    else:
                        c = candidates[pos]
                        c['months_below'] += 1
                        if freq < c['min_freq']:
                            c['min_freq'] = freq
                            c['min_month'] = month

    return candidates


# ──────────────────────────────────────────────────────────────────────────────
# Phase 2: build rich profiles from frequencies.tsv
# ──────────────────────────────────────────────────────────────────────────────

@dataclass
class PositionProfile:
    """Rich trajectory data for one AA position."""
    position: int
    total_series: list[tuple[str, float]] = dc_field(default_factory=list)
    per_aa_series: dict[str, list[tuple[str, float]]] = dc_field(default_factory=dict)
    dominant_history: list[tuple[str, str | None]] = dc_field(default_factory=list)

    # ── Derived metrics (filled by compute_metrics) ──
    max_mutant_freq: float = 0.0
    last_mutant_freq: float = 0.0
    n_sweeps: int = 0              # times freq crossed 0.20 going up
    n_dominant_mutations: int = 0  # distinct AAs ever dominant (>10%)
    dominant_switches: int = 0     # times the dominant AA changed
    volatility: float = 0.0       # mean |Δfreq| per month
    classification: str = 'unknown'

    def compute_metrics(self) -> None:
        """Compute derived metrics from the raw series."""
        freqs = [f for _, f in self.total_series]
        if not freqs:
            return
        self.max_mutant_freq = max(freqs)
        self.last_mutant_freq = freqs[-1]

        # ── Sweep count: transitions across 0.20 going upward ──
        for i in range(1, len(freqs)):
            if freqs[i - 1] < 0.20 <= freqs[i]:
                self.n_sweeps += 1

        # ── Volatility: mean |Δfreq| per month ──
        if len(freqs) >= 2:
            deltas = [abs(freqs[i] - freqs[i - 1]) for i in range(1, len(freqs))]
            self.volatility = sum(deltas) / len(deltas)

        # ── Dominant-mutation analysis ──
        dominants_seen: set[str] = set()
        prev_dom: str | None = None
        for _month, dom_aa in self.dominant_history:
            if dom_aa is not None:
                dominants_seen.add(dom_aa)
                if prev_dom is not None and dom_aa != prev_dom:
                    self.dominant_switches += 1
                prev_dom = dom_aa

        self.n_dominant_mutations = len(dominants_seen)

        # ── Classification ──
        self.classification = _classify(freqs, self)


def _classify(freqs: list[float], prof: PositionProfile) -> str:
    """Classify a mutation trajectory using rich metrics.

    Categories (from most to least specific):
        'replacement'    - dominant AA switched ≥1 time (variant replacement)
        'recurrence'     - ≥2 sweeps (frequency rose, fell, rose again)
        'fixation'       - rose to ≥90% and stayed; single dominant mutation
        'transient'      - peaked ≥20% then fell back below 5%
        'growing'        - trending upward in recent months
        'declining'      - was high, now dropping
        'persistent_low' - never exceeded 20%
        'variable'       - everything else
    """
    if not freqs:
        return 'unknown'

    last_n = freqs[-3:] if len(freqs) >= 3 else freqs
    first_n = freqs[:3] if len(freqs) >= 3 else freqs
    avg_last = sum(last_n) / len(last_n)
    avg_first = sum(first_n) / len(first_n)

    # Variant replacement: the dominant mutant AA changed at this position
    if prof.dominant_switches >= 1 and prof.max_mutant_freq >= 0.50:
        return 'replacement'

    # Recurrence: multiple sweeps through 20%
    if prof.n_sweeps >= 2:
        return 'recurrence'

    # Fixation: reached high frequency and stayed
    if prof.max_mutant_freq >= 0.90 and avg_last >= 0.85:
        return 'fixation'

    # Transient: peaked then vanished
    if prof.max_mutant_freq >= 0.20 and avg_last < 0.05:
        return 'transient'

    # Growing: trending upward
    if avg_last > avg_first + 0.10 and avg_last >= 0.10:
        return 'growing'

    # Declining: was high, now lower
    if avg_first >= 0.20 and avg_last < avg_first * 0.5:
        return 'declining'

    # Persistent low
    if prof.max_mutant_freq < 0.20:
        return 'persistent_low'

    return 'variable'


def build_profiles(
    directory: str,
    pattern: str,
    positions: set[int],
    start_month: str,
    end_month: str,
) -> dict[int, PositionProfile]:
    """Build rich PositionProfile objects from monthly frequency files."""
    files = sorted(glob.glob(os.path.join(directory, pattern)))

    # position -> month -> {mutant_aa: sum_freq}
    pos_month_aa: dict[int, dict[str, dict[str, float]]] = defaultdict(
        lambda: defaultdict(lambda: defaultdict(float))
    )
    all_months: set[str] = set()

    for fpath in files:
        month = _extract_month(fpath)
        if month is None or month.endswith('-00'):
            continue
        if not _month_in_range(month, start_month, end_month):
            continue
        all_months.add(month)

        with open(fpath, newline='') as fh:
            reader = csv.DictReader(fh, delimiter='\t')
            for row in reader:
                try:
                    pos = int(row['position'])
                    freq = float(row['frequency'])
                    mut_aa = row['mutant_aa']
                except (ValueError, KeyError):
                    continue
                if pos in positions:
                    pos_month_aa[pos][month][mut_aa] += freq

    # Build profiles
    sorted_months = sorted(all_months)
    profiles: dict[int, PositionProfile] = {}

    for pos in sorted(positions):
        prof = PositionProfile(position=pos)

        # Collect all mutant AAs seen at this position
        all_aas: set[str] = set()
        for m in sorted_months:
            for aa in pos_month_aa[pos].get(m, {}):
                all_aas.add(aa)

        for m in sorted_months:
            aa_freqs = pos_month_aa[pos].get(m, {})
            total = sum(aa_freqs.values())
            prof.total_series.append((m, total))

            # Dominant mutant AA this month (highest freq, must exceed 10%)
            if aa_freqs:
                best_aa = max(aa_freqs, key=aa_freqs.get)
                if aa_freqs[best_aa] >= 0.10:
                    prof.dominant_history.append((m, best_aa))
                else:
                    prof.dominant_history.append((m, None))
            else:
                prof.dominant_history.append((m, None))

        # Per-AA series
        for aa in sorted(all_aas):
            prof.per_aa_series[aa] = [
                (m, pos_month_aa[pos].get(m, {}).get(aa, 0.0))
                for m in sorted_months
            ]

        prof.compute_metrics()
        profiles[pos] = prof

    return profiles


# ──────────────────────────────────────────────────────────────────────────────
# Main
# ──────────────────────────────────────────────────────────────────────────────

def main() -> None:
    parser = argparse.ArgumentParser(
        description='Find AA positions with significant mutation dynamics.',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    parser.add_argument(
        '--dir', required=True,
        help='Directory containing monthly TSV files.',
    )
    parser.add_argument(
        '--unchanged-pattern',
        default='*.frequencies.unchanged_codons.tsv',
        help='Glob pattern for unchanged-codons files '
             '[default: %(default)s].',
    )
    parser.add_argument(
        '--freq-pattern',
        default='*.frequencies.tsv',
        help='Glob pattern for mutation-frequencies files '
             '[default: %(default)s].',
    )
    parser.add_argument(
        '--threshold', type=float, default=0.10,
        help='Minimum total mutant frequency to flag a position '
             '(i.e. unchanged < 1-threshold) [default: %(default)s].',
    )
    parser.add_argument(
        '--start-month', default='2020-01',
        help='Earliest month to consider [default: %(default)s].',
    )
    parser.add_argument(
        '--end-month', default='2025-07',
        help='Latest month to consider [default: %(default)s].',
    )
    parser.add_argument(
        '--min-months', type=int, default=3,
        help='Minimum number of months a position must be below threshold '
             'to be reported [default: %(default)s].',
    )
    parser.add_argument(
        '--output-format', choices=['positions', 'tsv', 'both'],
        default='both',
        help='Output format: "positions" for a space-separated list suitable '
             'for --positions, "tsv" for detailed table, "both" for both '
             '[default: %(default)s].',
    )

    args = parser.parse_args()

    # Phase 1: find candidates from unchanged_codons files
    print(f"Phase 1: Scanning unchanged codons (threshold={args.threshold})...",
          file=sys.stderr)
    candidates = find_candidate_positions(
        args.dir, args.unchanged_pattern, args.threshold,
        args.start_month, args.end_month,
    )

    # Filter by minimum number of months
    candidates = {
        pos: info for pos, info in candidates.items()
        if info['months_below'] >= args.min_months
    }

    if not candidates:
        print("No positions found matching criteria.", file=sys.stderr)
        sys.exit(0)

    print(f"  Found {len(candidates)} candidate positions.", file=sys.stderr)

    # Phase 2: build rich profiles from frequencies files
    print("Phase 2: Building mutation profiles...", file=sys.stderr)
    profiles = build_profiles(
        args.dir, args.freq_pattern, set(candidates.keys()),
        args.start_month, args.end_month,
    )

    # Merge candidate info with profiles and build results
    results = []
    for pos in sorted(candidates.keys()):
        info = candidates[pos]
        prof = profiles.get(pos)
        if prof is None:
            continue

        # Summarise dominant mutations seen
        dom_aas = sorted({
            aa for _, aa in prof.dominant_history if aa is not None
        })
        dom_summary = ','.join(dom_aas) if dom_aas else '-'

        results.append({
            'position': pos,
            'classification': prof.classification,
            'max_mutant_freq': prof.max_mutant_freq,
            'last_mutant_freq': prof.last_mutant_freq,
            'min_unchanged': info['min_freq'],
            'min_month': info['min_month'],
            'months_below': info['months_below'],
            'n_sweeps': prof.n_sweeps,
            'dom_switches': prof.dominant_switches,
            'n_dom_aas': prof.n_dominant_mutations,
            'dom_aas': dom_summary,
            'volatility': prof.volatility,
        })

    # Sort by max_mutant_freq descending (most dynamic first)
    results.sort(key=lambda r: -r['max_mutant_freq'])

    # Output
    # TSV table: sorted by max_mutant_freq (most dynamic first).
    # Positions list: numerically sorted (matches timeline renderer order).
    positions_list = sorted([r['position'] for r in results])

    if args.output_format in ('tsv', 'both'):
        header = [
            'position', 'classification', 'max_mutant_freq',
            'last_mutant_freq', 'min_unchanged_freq', 'worst_month',
            'months_below_threshold', 'n_sweeps', 'dom_switches',
            'n_dom_AAs', 'dominant_AAs', 'volatility',
        ]
        print('\t'.join(header))
        for r in results:
            print('\t'.join([
                str(r['position']),
                r['classification'],
                f"{r['max_mutant_freq']:.4f}",
                f"{r['last_mutant_freq']:.4f}",
                f"{r['min_unchanged']:.4f}",
                r['min_month'],
                str(r['months_below']),
                str(r['n_sweeps']),
                str(r['dom_switches']),
                str(r['n_dom_aas']),
                r['dom_aas'],
                f"{r['volatility']:.4f}",
            ]))

    if args.output_format in ('positions', 'both'):
        if args.output_format == 'both':
            print(file=sys.stderr)
        print(f"# mutation_timeline_plot --positions for {len(results)} sites:",
              file=sys.stderr)
        print(' '.join(str(p) for p in positions_list), file=sys.stderr)

        # Also print a copy-pasteable command
        print(f"\n# Ready-to-use command:", file=sys.stderr)
        print(
            f"mutation_timeline_plot --dir {args.dir} "
            f"--positions {' '.join(str(p) for p in positions_list)} "
            f"--pattern '*.frequencies.tsv' "
            f"--threshold=0.0001",
            file=sys.stderr,
        )


if __name__ == '__main__':
    main()
