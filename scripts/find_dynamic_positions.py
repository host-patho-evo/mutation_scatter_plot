#!/usr/bin/env python3
"""Scan monthly frequency data to find positions with significant mutation dynamics.

**Strategy** (two-phase):

1. **Candidate detection** — parse ``*.unchanged_codons.tsv`` files and collect
   every AA position whose unchanged-codon frequency dropped below a threshold
   (default 90 %) in *at least one month* since 2020-01.  These are sites
   where the reference codon was significantly displaced at some point.

2. **Dynamics characterisation** — for each candidate, walk through the
   ``*.frequencies.tsv`` files month-by-month and compute the *total mutant
   frequency* (1 − unchanged_freq) to build a trajectory.  Classify the
   trajectory as "rose to fixation", "transient sweep", or "persistent low".

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

import argparse
import csv
import glob
import os
import re
import sys
from collections import defaultdict


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
# Phase 2: build trajectories from frequencies.tsv
# ──────────────────────────────────────────────────────────────────────────────

def build_trajectories(
    directory: str,
    pattern: str,
    positions: set[int],
    start_month: str,
    end_month: str,
) -> dict[int, list[tuple[str, float]]]:
    """For each position, collect (month, total_mutant_freq) time-series.

    total_mutant_freq is the sum of all non-reference mutation frequencies
    at that position in a given month.
    """
    files = sorted(glob.glob(os.path.join(directory, pattern)))

    # position -> month -> sum of mutant frequencies
    pos_month_freq: dict[int, dict[str, float]] = defaultdict(
        lambda: defaultdict(float)
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
                except (ValueError, KeyError):
                    continue
                if pos in positions:
                    pos_month_freq[pos][month] += freq

    # Convert to sorted time-series
    sorted_months = sorted(all_months)
    trajectories: dict[int, list[tuple[str, float]]] = {}
    for pos in sorted(positions):
        series = []
        for m in sorted_months:
            series.append((m, pos_month_freq[pos].get(m, 0.0)))
        trajectories[pos] = series

    return trajectories


def classify_trajectory(series: list[tuple[str, float]]) -> str:
    """Classify a mutation frequency trajectory.

    Returns one of:
        'fixation'       - rose to ≥90% and stayed
        'transient'      - peaked ≥20% then fell back below 5%
        'persistent_low' - stayed mostly below 20% but above threshold
        'growing'        - still increasing in latest months
        'declining'      - was high, now declining
        'variable'       - significant fluctuation
    """
    if not series:
        return 'unknown'

    freqs = [f for _, f in series]
    max_freq = max(freqs)
    last_n = freqs[-3:] if len(freqs) >= 3 else freqs
    first_n = freqs[:3] if len(freqs) >= 3 else freqs
    avg_last = sum(last_n) / len(last_n)
    avg_first = sum(first_n) / len(first_n)

    if max_freq >= 0.90 and avg_last >= 0.85:
        return 'fixation'
    if max_freq >= 0.20 and avg_last < 0.05:
        return 'transient'
    if avg_last > avg_first + 0.10 and avg_last >= 0.10:
        return 'growing'
    if avg_first >= 0.20 and avg_last < avg_first * 0.5:
        return 'declining'
    if max_freq < 0.20:
        return 'persistent_low'
    return 'variable'


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

    # Phase 2: build trajectories from frequencies files
    print("Phase 2: Building mutation trajectories...", file=sys.stderr)
    trajectories = build_trajectories(
        args.dir, args.freq_pattern, set(candidates.keys()),
        args.start_month, args.end_month,
    )

    # Classify and sort
    results = []
    for pos in sorted(candidates.keys()):
        info = candidates[pos]
        traj = trajectories.get(pos, [])
        classification = classify_trajectory(traj)
        max_mutant_freq = max((f for _, f in traj), default=0.0)
        last_freq = traj[-1][1] if traj else 0.0
        results.append({
            'position': pos,
            'min_unchanged': info['min_freq'],
            'min_month': info['min_month'],
            'months_below': info['months_below'],
            'max_mutant_freq': max_mutant_freq,
            'last_mutant_freq': last_freq,
            'classification': classification,
        })

    # Sort by max_mutant_freq descending (most dynamic first)
    results.sort(key=lambda r: -r['max_mutant_freq'])

    # Output
    positions_list = [str(r['position']) for r in results]

    if args.output_format in ('tsv', 'both'):
        header = [
            'position', 'classification', 'max_mutant_freq',
            'last_mutant_freq', 'min_unchanged_freq', 'worst_month',
            'months_below_threshold',
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
            ]))

    if args.output_format in ('positions', 'both'):
        if args.output_format == 'both':
            print(file=sys.stderr)
        print(f"# mutation_timeline_plot --positions for {len(results)} sites:",
              file=sys.stderr)
        print(' '.join(positions_list), file=sys.stderr)

        # Also print a copy-pasteable command
        print(f"\n# Ready-to-use command:", file=sys.stderr)
        print(
            f"mutation_timeline_plot --dir {args.dir} "
            f"--positions {' '.join(positions_list)} "
            f"--pattern '*.frequencies.tsv' "
            f"--threshold=0.0001",
            file=sys.stderr,
        )


if __name__ == '__main__':
    main()
