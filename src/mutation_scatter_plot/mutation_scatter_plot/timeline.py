# -*- coding: utf-8 -*-
# vim: set fileencoding=utf-8 ts=4 sw=4 expandtab :
# This work © 2025-2026 by Jiří Zahradník and Martin Mokrejš
# (First Medical Faculty - Charles University in Prague) is licensed under
# Creative Commons Attribution 4.0 International. To view a copy of this
# license, visit https://creativecommons.org/licenses/by/4.0/

"""Timeline scatter plot for mutation frequencies across monthly datasets.

This module scans a directory of per-month ``.frequencies.tsv`` files produced
by ``calculate_codon_frequencies`` and renders a timeline scatter plot where:

- **X-axis** = time (YYYY-MM parsed from filenames)
- **Y-axis** = selected amino acid / codon positions
- **Circle colour** = BLOSUM substitution score (same scheme as mutation_scatter_plot)
- **Circle size** = mutation frequency (scaled down vs scatter plots)

The coloring, scoring, and BLOSUM matrix logic are reused from ``core.py``.
"""

import glob
import os
import re
import typing
from collections import defaultdict
from dataclasses import dataclass, field
from decimal import Decimal

import matplotlib.pyplot as plt
import pandas as pd

from .core import (
    resolve_codon_or_aa,
    adjust_size_and_color,
)
from .colorbar_helpers import setup_matplotlib_colorbar, add_bokeh_colorbar


# ── Data structures ──────────────────────────────────────────────────────

@dataclass
class PositionSpec:
    """A parsed position specification from the CLI.

    Examples
    --------
    ``PositionSpec(position=614, ref_aa=None, mutant_aas=None)``
        Track all mutations at position 614.
    ``PositionSpec(position=614, ref_aa='D', mutant_aas=['G'])``
        Track only the D614G mutation.
    ``PositionSpec(position=498, ref_aa=None, mutant_aas=['R', 'H', 'Q'])``
        Track mutations to R, H, or Q at position 498.
    """
    position: int
    ref_aa: typing.Optional[str] = None
    mutant_aas: typing.Optional[list[str]] = None

    def __str__(self):
        """Return human-readable position label, e.g. 'D614G' or '498[RHQ]'."""
        if self.ref_aa and self.mutant_aas and len(self.mutant_aas) == 1:
            return f"{self.ref_aa}{self.position}{self.mutant_aas[0]}"
        if self.mutant_aas:
            return f"{self.position}[{''.join(self.mutant_aas)}]"
        return str(self.position)


@dataclass
class TimelinePoint:
    """A single data point in the timeline scatter plot.

    Each point represents one observed mutation (old→new amino acid) at a
    specific genomic position in a specific month, with its frequency and
    BLOSUM substitution score.  The ``color`` and ``label`` fields are
    computed during data collection.
    """
    month: str           # 'YYYY-MM'
    position: int        # aa position
    ref_aa: str          # reference amino acid
    mutant_aa: str       # mutant amino acid
    ref_codon: str       # reference codon
    mutant_codon: str    # mutant codon
    frequency: Decimal   # observed frequency
    color: str = ''      # hex colour (computed)
    score: int = 0       # BLOSUM score (computed)
    label: str = ''      # e.g. 'D614G'


@dataclass
class TimelineData:
    """Container for parsed per-month, per-position frequency data.

    Aggregates all :class:`TimelinePoint` instances together with sorted
    month labels, unique positions, and the original position specifications
    from the CLI.  Passed to both ``render_timeline_matplotlib`` and
    ``render_timeline_bokeh``.
    """
    points: list[TimelinePoint] = field(default_factory=list)
    months: list[str] = field(default_factory=list)    # sorted YYYY-MM
    positions: list[int] = field(default_factory=list)  # sorted unique positions seen
    position_specs: list[PositionSpec] = field(default_factory=list)


def recolor_timeline_data(
    data: TimelineData,
    myoptions: typing.Any,
    norm: typing.Any,
    colors: list[str],
) -> None:
    """Recompute the ``.color`` attribute on every point in *data*.

    This allows rendering the same collected data with a different colormap
    without re-reading files or re-computing BLOSUM scores.  The function
    replicates the colour-assignment logic of
    :func:`~core.adjust_size_and_color` (synonymous → dark green,
    DEL/INS → red, X/NNN → gray, otherwise index into *colors*), but
    reads the already-stored ``.score`` instead of calling ``get_score()``.

    Parameters
    ----------
    data : TimelineData
        Collected data whose ``.color`` fields will be mutated in-place.
    myoptions : argparse.Namespace
        CLI options (used for ``cmap_vmin``).
    norm : matplotlib.colors.BoundaryNorm or None
        Discrete colour normaliser, or None for continuous colormaps.
    colors : list of str
        Hex colour palette.
    """

    for pt in data.points:
        # Synonymous sentinel (score == 12) → dark green
        if pt.score == 12:
            pt.color = '#219f11'
            continue

        # DEL / INS / stop → red
        if pt.mutant_codon.upper() in ('---', 'DEL', 'INS'):
            pt.color = '#ff0000'
            continue

        # Ambiguous / unknown → gray
        if pt.ref_codon.upper() in ('NNN',) or pt.mutant_codon.upper() in ('NNN',):
            pt.color = '#808080'
            continue
        if pt.ref_aa.upper() in ('X',) or pt.mutant_aa.upper() in ('X',):
            pt.color = '#808080'
            continue

        # Normal scoring path — map score to colour index
        if norm is not None:
            _colorindex = norm(pt.score)
        else:
            _vmin = getattr(myoptions, 'cmap_vmin', -11)
            _colorindex = max(0, min(len(colors) - 1, pt.score - _vmin))

        pt.color = colors[_colorindex]


def filter_timeline_data(
    data: TimelineData,
    positions: list[int],
) -> TimelineData:
    """Return a new :class:`TimelineData` containing only *positions*.

    All original months are preserved so the X-axis remains consistent
    across paginated outputs.

    Parameters
    ----------
    data : TimelineData
        The full collected data.
    positions : list of int
        Subset of positions to keep.

    Returns
    -------
    TimelineData
        Filtered copy (points list is new; month list is shared).
    """
    pos_set = set(positions)
    return TimelineData(
        points=[pt for pt in data.points if pt.position in pos_set],
        months=list(data.months),
        positions=sorted(p for p in data.positions if p in pos_set),
        position_specs=[s for s in data.position_specs if s.position in pos_set],
    )


# ── File scanning ────────────────────────────────────────────────────────

_MONTH_RE = re.compile(r'\.(\d{4})-(\d{2})\.')


def scan_directory(
    dirpath: str,
    pattern: str = '*.frequencies.tsv',
    include_unknown_month: bool = False,
) -> list[tuple[str, str]]:
    """Scan *dirpath* for files matching *pattern* and extract YYYY-MM dates.

    Parameters
    ----------
    dirpath : str
        Directory to scan.
    pattern : str
        Glob pattern for input files.
    include_unknown_month : bool
        If False (default), skip files where month is '00' (unknown).

    Returns
    -------
    list of (month_str, filepath) tuples, sorted chronologically.
    """
    results: list[tuple[str, str]] = []
    for fpath in sorted(glob.glob(os.path.join(dirpath, pattern))):
        m = _MONTH_RE.search(os.path.basename(fpath))
        if not m:
            print(f"Info: Skipping {os.path.basename(fpath)} — no YYYY-MM found in filename")
            continue
        year, month = m.group(1), m.group(2)
        if month == '00' and not include_unknown_month:
            print(f"Info: Skipping {os.path.basename(fpath)} — month '00' (use --include-unknown-month)")
            continue
        month_str = f"{year}-{month}"
        results.append((month_str, fpath))
    results.sort(key=lambda x: x[0])
    return results


def infer_common_prefix(
    files: list[tuple[str, str]],
    dirpath: str,
) -> str:
    """Infer a common output prefix from scanned filenames.

    Extracts the portion of each filename that precedes the ``YYYY-MM``
    date stamp, then finds the longest common prefix among those stems.
    The result includes the input directory so that output files land
    in the same directory as the input unless overridden.

    Examples
    --------
    >>> files = [
    ...   ('2020-01', '/data/spikenuc0719.no_junk.2020-01.counts.tsv'),
    ...   ('2020-02', '/data/spikenuc0719.no_junk.2020-02.counts.tsv'),
    ... ]
    >>> infer_common_prefix(files, '/data')
    '/data/spikenuc0719.no_junk'
    """
    if not files:
        return os.path.join(dirpath, 'timeline')

    stems: list[str] = []
    for _, fpath in files:
        basename = os.path.basename(fpath)
        m = _MONTH_RE.search(basename)
        if m:
            # Take everything before the .YYYY-MM. match, stripping trailing dot
            prefix = basename[:m.start()]
            if prefix.endswith('.'):
                prefix = prefix[:-1]
            stems.append(prefix)
        else:
            stems.append(basename)

    if not stems:
        return os.path.join(dirpath, 'timeline')

    # Find longest common prefix
    common = stems[0]
    for s in stems[1:]:
        while not s.startswith(common):
            common = common[:-1]
            if not common:
                break
    # Strip trailing dots/separators from the common prefix
    common = common.rstrip('.')

    if not common:
        common = 'timeline'

    return os.path.join(dirpath, common)


# ── Position parsing ─────────────────────────────────────────────────────

_MUTATION_RE = re.compile(r'^([A-Z])(\d+)([A-Z])$')
_BRACKET_RE = re.compile(r'^(\d+)\[([A-Z]+)\]$')
_REF_BRACKET_RE = re.compile(r'^([A-Z])(\d+)\[([A-Z]+)\]$')


def parse_positions(position_args: list[str]) -> list[PositionSpec]:
    """Parse position specifications from CLI arguments.

    Accepts:
    - Numeric positions: ``'614'``, ``'501'``
    - Ranges: ``'480-530'``
    - Mutation labels: ``'D614G'``, ``'N501Y'``
    - Bracket notation: ``'498[RHQ]'`` (mutations to R, H, or Q at pos 498)
    - Ref+bracket: ``'Q498[RH]'`` (Q→R or Q→H at pos 498)

    Returns a list of :class:`PositionSpec` objects.
    """
    specs: list[PositionSpec] = []
    for arg in position_args:
        arg_upper = arg.upper()

        # Try ref + bracket notation (e.g. Q498[RH])
        m = _REF_BRACKET_RE.match(arg_upper)
        if m:
            specs.append(PositionSpec(
                position=int(m.group(2)),
                ref_aa=m.group(1),
                mutant_aas=list(m.group(3)),
            ))
            continue

        # Try bracket notation without ref (e.g. 498[RHQ])
        m = _BRACKET_RE.match(arg_upper)
        if m:
            specs.append(PositionSpec(
                position=int(m.group(1)),
                ref_aa=None,
                mutant_aas=list(m.group(2)),
            ))
            continue

        # Try mutation label (e.g. D614G)
        m = _MUTATION_RE.match(arg_upper)
        if m:
            specs.append(PositionSpec(
                position=int(m.group(2)),
                ref_aa=m.group(1),
                mutant_aas=[m.group(3)],
            ))
            continue
        # Try range (e.g. 480-530)
        if '-' in arg and not arg.startswith('-'):
            parts = arg.split('-', 1)
            try:
                start, end = int(parts[0]), int(parts[1])
                for pos in range(start, end + 1):
                    specs.append(PositionSpec(position=pos))
                continue
            except ValueError:
                pass
        # Try single numeric position
        try:
            specs.append(PositionSpec(position=int(arg)))
        except ValueError:
            print(f"Warning: Cannot parse position specification '{arg}' — skipping")

    # ── Deduplicate: merge specs for the same position into a union ──
    # Rules:
    #   - Unfiltered spec (mutant_aas=None) subsumes all filtered specs for that position.
    #   - Multiple filtered specs: union of mutant_aas, keep ref_aa only if all agree.
    merged: dict[int, PositionSpec] = {}
    for s in specs:
        if s.position not in merged:
            merged[s.position] = PositionSpec(
                position=s.position,
                ref_aa=s.ref_aa,
                mutant_aas=list(s.mutant_aas) if s.mutant_aas else None,
            )
        else:
            existing = merged[s.position]
            if existing.mutant_aas is None:
                # Already unfiltered — subsumes everything
                continue
            if s.mutant_aas is None:
                # New spec is unfiltered — subsumes existing
                existing.mutant_aas = None
                existing.ref_aa = None
                continue
            # Both filtered: merge mutant_aas union
            combined = set(existing.mutant_aas) | set(s.mutant_aas)
            existing.mutant_aas = sorted(combined)
            # Keep ref_aa only if both agree
            if existing.ref_aa != s.ref_aa:
                existing.ref_aa = None

    result = list(merged.values())
    # Preserve position order (sorted numerically)
    result.sort(key=lambda s: s.position)

    if len(result) < len(specs):
        print(f"Info: Merged {len(specs)} position specs into {len(result)} unique positions")

    return result


# ── Data collection ──────────────────────────────────────────────────────

def collect_timeline_data(
    files: list[tuple[str, str]],
    specs: list[PositionSpec],
    myoptions: typing.Any,
    matrix: typing.Any,
    norm: typing.Any,
    colors: list[str],
) -> TimelineData:
    """Load data from monthly TSV files for the requested positions.

    For each file in *files*, reads the ``.frequencies.tsv`` table, filters
    rows to the amino acid positions listed in *specs*, and computes the
    BLOSUM substitution score and colour using :func:`adjust_size_and_color`
    from ``core.py``.

    Parameters
    ----------
    files : list of (month_str, filepath)
        Output of :func:`scan_directory`.
    specs : list of PositionSpec
        Parsed position filters from :func:`parse_positions`.
    myoptions : argparse.Namespace
        CLI options (threshold, debug, aminoacids, etc.).
    matrix : blosum.BLOSUM
        Active substitution-score matrix.
    norm : matplotlib.colors.BoundaryNorm or None
        Discrete colour normaliser (for ListedColormap), or None.
    colors : list of str
        Hex colour palette, indexed by normalised score.

    Returns
    -------
    TimelineData
        Populated container with all matching data points, sorted month
        list, and sorted position list.
    """
    all_positions = set()
    for s in specs:
        all_positions.add(s.position)

    # Build a lookup: position → list of PositionSpec (for mutation-label filtering)
    pos_to_specs: dict[int, list[PositionSpec]] = defaultdict(list)
    for s in specs:
        pos_to_specs[s.position].append(s)

    data = TimelineData(position_specs=specs)
    seen_months: set[str] = set()
    seen_positions: set[int] = set()

    for month_str, fpath in files:
        try:
            df = pd.read_csv(fpath, sep='\t', dtype={'position': int})
        except Exception as exc:
            print(f"Warning: Cannot read {fpath}: {exc}")
            continue

        # Normalise column names (handle legacy/new formats)
        if 'position' not in df.columns and 'padded_position' in df.columns:
            df['position'] = df['padded_position']

        # Filter to requested positions
        df_filtered = df[df['position'].isin(all_positions)]
        if df_filtered.empty:
            continue

        seen_months.add(month_str)

        for _, row in df_filtered.iterrows():
            pos = int(row['position'])
            ref_codon = str(row.get('original_codon', ''))
            mut_codon = str(row.get('mutant_codon', ''))
            ref_aa = str(row.get('original_aa', ''))
            mut_aa = str(row.get('mutant_aa', ''))
            freq_val = Decimal(str(row.get('frequency', 0)))

            # Apply mutation-label or bracket filter if specified
            matched_specs = pos_to_specs.get(pos, [])
            if matched_specs:
                # Check if any spec has a mutant filter
                has_filter = any(s.mutant_aas is not None for s in matched_specs)
                if has_filter:
                    matches = False
                    for s in matched_specs:
                        if s.mutant_aas is None:
                            # Unfiltered spec for this position — accept all
                            matches = True
                            break
                        # Check ref_aa if specified
                        if s.ref_aa and s.ref_aa != ref_aa.upper():
                            continue
                        # Check mutant_aa against the list
                        if mut_aa.upper() in s.mutant_aas:
                            matches = True
                            break
                    if not matches:
                        continue

            # Apply frequency threshold
            if abs(freq_val) < myoptions.threshold:
                continue

            # Skip synonymous changes in aminoacid mode (unless --include-synonymous)
            if (getattr(myoptions, 'aminoacids', False)
                    and not getattr(myoptions, 'include_synonymous', False)
                    and ref_aa and mut_aa and ref_aa == mut_aa):
                continue
            # Compute colour/score using core functions
            codon_on_input, _old, _new = resolve_codon_or_aa(
                myoptions, ref_codon, mut_codon
            )
            _score, _freq, _color = adjust_size_and_color(
                myoptions, freq_val, codon_on_input,
                ref_codon, mut_codon, _old, _new,
                matrix, norm, colors,
            )

            label = f"{ref_aa}{pos}{mut_aa}"

            pt = TimelinePoint(
                month=month_str,
                position=pos,
                ref_aa=ref_aa,
                mutant_aa=mut_aa,
                ref_codon=ref_codon,
                mutant_codon=mut_codon,
                frequency=freq_val,
                color=_color,
                score=_score,
                label=label,
            )
            data.points.append(pt)
            seen_positions.add(pos)

    data.months = sorted(seen_months)
    data.positions = sorted(seen_positions)
    return data


# ── Matplotlib rendering ─────────────────────────────────────────────────

# Timeline-specific scaling (much smaller than scatter plot's 3000-5000)
TIMELINE_CIRCLE_SCALE = 800
TIMELINE_MIN_SIZE = 10
TIMELINE_MAX_SIZE = 800


def _month_to_float(month_str: str, all_months: list[str]) -> float:
    """Convert a YYYY-MM string to a float index for plotting."""
    return all_months.index(month_str) if month_str in all_months else 0.0


def _compute_band_spacing(data: TimelineData) -> float:
    """Dynamically compute vertical band spacing based on data density.

    Examines the maximum number of mutations sharing a single
    (month, position) slot and scales spacing so that vertically
    offset circles do not bleed into neighbouring bands.

    Returns
    -------
    float
        Band spacing in data units (minimum 1.2).
    """
    if not data.points:
        return 1.5

    # Count max overlap per (month, position)
    overlap_counts: dict[tuple[str, int], int] = defaultdict(int)
    for pt in data.points:
        overlap_counts[(pt.month, pt.position)] += 1
    max_overlap = max(overlap_counts.values()) if overlap_counts else 1

    # Scale: 1 mutation → 2.0, 2 → 2.6, 3 → 3.2, 5 → 4.4, ...
    spacing = max(2.0, max_overlap * 0.6 + 1.4)
    return spacing


def _compute_intra_band_spread(band_spacing: float) -> float:
    """Compute the ± vertical offset range within a band.

    Uses 60% of the half-band width so circles stay well inside
    the band borders.
    """
    return band_spacing * 0.5 * 0.6


def _slot_key(pt: TimelinePoint) -> str:
    """Return a unique key for vertical slot assignment within a band.

    Uses codon-level identity so that different codon mutations producing
    the same AA change (e.g. GAT→GGT vs GAT→GGC, both D614G) get
    distinct vertical positions.
    """
    return f"{pt.ref_codon}{pt.position}{pt.mutant_codon}"


def _prepare_layout(
    data: TimelineData,
    myoptions: typing.Any = None,
) -> tuple[dict[int, float], dict[int, float],
              dict[int, float],
              dict[tuple[str, int], list['TimelinePoint']],
              dict[int, dict[str, float]]]:
    """Compute shared layout parameters for both renderers.

    Parameters
    ----------
    data : TimelineData
        Collected timeline data.
    myoptions : argparse.Namespace, optional
        CLI options.  If provided, ``band_spacing_factor`` is used to
        scale the auto-computed band spacing.

    Returns
    -------
    tuple of (pos_heights, pos_spread, pos_to_y, grouped, label_offsets)
        - pos_heights: ``{position: height}`` — per-position band height.
        - pos_spread: ``{position: spread}`` — ± vertical offset range
          within each band.
        - pos_to_y: mapping from position int to y-coordinate (band centre).
        - grouped: data points grouped by (month, position).
        - label_offsets: ``{position: {label: y_offset}}`` — fixed
          vertical slot for each unique mutation label within a band,
          consistent across all months.
    """
    # Group data points
    grouped: dict[tuple[str, int], list[TimelinePoint]] = defaultdict(list)
    for pt in data.points:
        grouped[(pt.month, pt.position)].append(pt)

    # Collect unique codon mutations per position (needed for height + offset)
    # Uses codon-level identity so different codons producing the same AA
    # change get distinct vertical slots instead of overlapping.
    slots_per_pos: dict[int, list[str]] = {}
    for pos in data.positions:
        _unique: set[str] = set()
        for (_m, _p), pts in grouped.items():
            if _p == pos:
                for pt in pts:
                    _unique.add(_slot_key(pt))
        slots_per_pos[pos] = sorted(_unique)

    # Apply user-specified scaling factor
    _factor = getattr(myoptions, 'band_spacing_factor', 1.0) if myoptions else 1.0

    # Compute per-position band heights.
    # Height scales with the number of unique labels at each position,
    # so dense positions get taller bands while simple positions stay compact.
    # The constants below must be large enough for the circle + 3-line
    # annotation text at TIMELINE_CIRCLE_SCALE=800 and 600 DPI.
    _SLOT_HEIGHT = 1.4    # data units per slot (circle + annotation)
    _MIN_BAND_HEIGHT = 2.0  # minimum for positions with 1 mutation

    pos_heights: dict[int, float] = {}
    for pos in data.positions:
        n_slots = len(slots_per_pos[pos])
        # Height = padding + n_slots * per-slot height
        pos_heights[pos] = max(_MIN_BAND_HEIGHT, n_slots * _SLOT_HEIGHT + 0.6) * _factor

    # Per-position spread: ± vertical offset range within each band
    pos_spread: dict[int, float] = {}
    for pos in data.positions:
        pos_spread[pos] = _compute_intra_band_spread(pos_heights[pos])

    # Compute cumulative y positions (variable spacing)
    pos_to_y: dict[int, float] = {}
    y_cursor = 0.0
    for pos in data.positions:
        half_h = pos_heights[pos] * 0.5
        pos_to_y[pos] = y_cursor + half_h
        y_cursor += pos_heights[pos]

    # Pre-compute consistent vertical slots for each mutation label within
    # each position band.  Every label gets the same y-offset regardless of
    # how many other labels appear in a given month, preventing circles
    # from jumping vertically.
    label_offsets: dict[int, dict[str, float]] = {}
    for pos in data.positions:
        _sorted = slots_per_pos[pos]
        _n = len(_sorted)
        _spr = pos_spread[pos]
        offsets: dict[str, float] = {}
        if _n <= 1:
            for key in _sorted:
                offsets[key] = 0.0
        else:
            for j, key in enumerate(_sorted):
                offsets[key] = -_spr + (2 * _spr * j / (_n - 1))
        label_offsets[pos] = offsets

    return pos_heights, pos_spread, pos_to_y, grouped, label_offsets


def _compute_ytick_labels(
    positions: list[int],
    pos_to_y: dict[int, float],
    grouped: dict[tuple[str, int], list['TimelinePoint']],
    pos_spread: dict[int, float],
    merge_threshold: float = 0.15,
) -> tuple[list[float], list[str]]:
    """Compute merged Y-axis tick positions and labels.

    For each position band, collects all unique mutation labels, merges
    labels that would be too close, then **evenly distributes** the
    resulting labels across the band height.  This prevents overlap
    regardless of the actual circle positions within the band.

    Returns
    -------
    tuple of (tick_positions, tick_labels)
        Parallel lists ready for ``ax.set_yticks`` / ``ax.set_yticklabels``
        or Bokeh ``yaxis.ticker`` / ``major_label_overrides``.
    """
    tick_y: list[float] = []
    tick_labels: list[str] = []

    for pos in positions:
        y_base = pos_to_y[pos]
        _spread = pos_spread.get(pos, 0.6)
        # Collect unique codon mutations (slot keys) with their display labels.
        # Uses _slot_key so that ytick positions match circle vertical slots.
        seen_slots: dict[str, str] = {}  # slot_key → display label
        for (_month, p_key), pts in grouped.items():
            if p_key != pos:
                continue
            for pt in pts:
                key = _slot_key(pt)
                if key not in seen_slots:
                    seen_slots[key] = pt.label

        slot_keys = sorted(seen_slots.keys())

        if not slot_keys:
            tick_y.append(y_base)
            tick_labels.append(str(pos))
            continue

        # Evenly distribute labels within the band [y_base - _spread, y_base + _spread]
        n_labels = len(slot_keys)
        if n_labels == 1:
            tick_y.append(y_base)
            tick_labels.append(seen_slots[slot_keys[0]])
        else:
            for j, key in enumerate(slot_keys):
                y = y_base - _spread + (2 * _spread * j / (n_labels - 1))
                tick_y.append(y)
                tick_labels.append(seen_slots[key])

    return tick_y, tick_labels


def _format_pct(freq: float) -> str:
    """Format a frequency as a percentage string.

    Examples: 0.14 → '14%', 0.015 → '1.5%', 0.002 → '0.20%'.
    """
    pct = freq * 100
    if pct >= 10:
        return f"{pct:.0f}%"
    if pct >= 1:
        return f"{pct:.1f}%"
    return f"{pct:.2f}%"


def render_timeline_matplotlib(
    data: TimelineData,
    myoptions: typing.Any,
    norm: typing.Any,
    cmap: typing.Any,
    colors: typing.Any,
    outfile_prefix: str,
) -> None:
    """Render the timeline scatter plot using matplotlib.

    Produces PNG and PDF files at ``{outfile_prefix}.png`` and ``.pdf``.

    Visual elements
    ---------------
    * Circles coloured by BLOSUM score (via *cmap* / *norm*).
    * Circle size proportional to mutation frequency (area or linear scaling).
    * Percentage labels annotated next to each circle.
    * Y-axis tick labels show mutation names (e.g. ``N501Y``, ``E484A, E484K``).
    * Colourbar trimmed to the data-driven score range.
    * Interactive hover via ``mplcursors`` (if installed).

    Parameters
    ----------
    data : TimelineData
        Collected timeline data points.
    myoptions : argparse.Namespace
        CLI options.
    norm : matplotlib.colors.BoundaryNorm or None
        Discrete colour normaliser, or None for continuous colormaps.
    cmap : matplotlib.colors.Colormap
        Colormap instance.
    colors : list or None
        Resolved palette list from ``get_colormap``.
    outfile_prefix : str
        Output path prefix (without extension).
    """
    if not data.points:
        print("Warning: No data points to render in timeline plot")
        return

    months = data.months
    positions = data.positions

    # Shared layout: per-position heights, spread, y-mapping, grouping, label offsets
    pos_heights, pos_spread, pos_to_y, grouped, label_offsets = _prepare_layout(data, myoptions)

    # Prepare scatter data
    x_vals: list[float] = []
    y_vals: list[float] = []
    sizes: list[float] = []
    colors_hex: list[str] = []
    freqs: list[float] = []
    labels: list[str] = []
    mut_labels: list[str] = []
    codon_changes: list[str] = []

    for (month, pos), pts in grouped.items():
        x = _month_to_float(month, months)
        y_base = pos_to_y[pos]
        _pos_offsets = label_offsets.get(pos, {})

        for pt in pts:
            # Fixed vertical slot for this mutation label
            y_offset = _pos_offsets.get(_slot_key(pt), 0.0)

            x_vals.append(x)
            y_vals.append(y_base + y_offset)

            # Scale circle size by frequency
            freq_f = float(pt.frequency)
            if getattr(myoptions, 'linear_circle_size', False):
                # Linear scaling: size proportional to frequency
                raw_size = freq_f * TIMELINE_CIRCLE_SCALE
            else:
                # Area scaling (default): size proportional to sqrt(frequency)
                # so that circle *area* is proportional to frequency
                raw_size = (freq_f ** 0.5) * TIMELINE_CIRCLE_SCALE
            sizes.append(max(TIMELINE_MIN_SIZE, min(TIMELINE_MAX_SIZE, raw_size)))
            colors_hex.append(pt.color)
            freqs.append(freq_f)
            mut_labels.append(pt.label)
            codon_changes.append(f"{pt.ref_codon}{pt.position}{pt.mutant_codon}")
            # Build hover text matching mutation_scatter_plot format
            _matrix_name = getattr(myoptions, 'matrix', 'BLOSUM80')
            _hover = (
                f"Position: {pt.position}\n"
                f"Original Codon: {pt.ref_codon} ({pt.ref_aa})\n"
                f"New Codon: {pt.mutant_codon} ({pt.mutant_aa})\n"
                f"{_matrix_name} score: {pt.score}\n"
                f"Frequency: {float(pt.frequency):.6f}"
            )
            labels.append(_hover)

    # Create figure with dedicated colorbar and legend columns
    n_pos = len(positions)
    # Total vertical extent = sum of all band heights
    _total_height = sum(pos_heights.get(p, 2.0) for p in positions)
    fig_height = max(5, _total_height * 0.5 + 2)
    fig_width = max(12, len(months) * 0.8 + 5)
    fig, (ax, ax_cb, ax_leg) = plt.subplots(
        1, 3, figsize=(fig_width, fig_height),
        width_ratios=[50, 1, 8],
    )
    fig.subplots_adjust(left=0.06, right=0.86, top=0.94, bottom=0.15)

    # Dynamic heading font: scales with the smaller figure dimension so that
    # title, axis labels, and 2nd-y-axis position numbers form a visually
    # consistent group.  Ranges from 11pt (compact plots with few bands) to
    # 20pt (spacious plots with many bands).
    _min_dim = min(fig_width, fig_height)
    _heading_fontsize = max(11, min(20, _min_dim * 1.8))

    title = getattr(myoptions, 'title', '') or f"Mutation Timeline ({len(months)} months, {n_pos} positions)"
    ax.set_title(title, fontsize=_heading_fontsize, fontweight='bold', pad=15)

    # Scatter plot
    scatter = ax.scatter(
        x_vals, y_vals,
        s=sizes,
        c=colors_hex,
        alpha=0.8,
        edgecolors='#333333',
        linewidths=0.5,
        zorder=5,
    )

    # Percentage + mutation labels next to circles
    for xi, yi, fi, ml, cc in zip(x_vals, y_vals, freqs, mut_labels, codon_changes):
        _annot_text = f"{_format_pct(fi)}\n{ml}\n{cc}"
        ax.annotate(_annot_text, (xi, yi), textcoords='offset points',
                    xytext=(4, 3), fontsize=5,
                    color='black', zorder=6)

    # X-axis: months
    ax.set_xticks(range(len(months)))
    ax.set_xticklabels(months, rotation=45, ha='right', fontsize=8)
    ax.set_xlabel('Month', fontsize=_heading_fontsize)
    # Y-axis: use merged mutation labels as tick labels
    _all_tick_y, _all_tick_labels = _compute_ytick_labels(
        positions, pos_to_y, grouped, pos_spread,
    )
    ax.set_yticks(_all_tick_y)
    ax.set_yticklabels(_all_tick_labels, fontsize=7)
    ax.set_ylabel('AA Position', fontsize=_heading_fontsize)

    # Grid and styling
    ax.set_xlim(-0.5, len(months) - 0.5)
    # y-limits: half-band padding around first and last position
    _first_half = pos_heights.get(positions[0], 2.0) * 0.5 if positions else 1.0
    _last_half = pos_heights.get(positions[-1], 2.0) * 0.5 if positions else 1.0
    _y_min = pos_to_y.get(positions[0], 0.0) - _first_half - 0.2 if positions else -0.2
    _y_max = pos_to_y.get(positions[-1], 0.0) + _last_half + 0.2 if positions else 0.2
    ax.set_ylim(_y_min, _y_max)
    ax.grid(axis='x', alpha=0.3, linestyle='--')
    # Draw horizontal band borders above and below each position
    # (not through the centre where data points are)
    for pos in positions:
        y_center = pos_to_y[pos]
        _half = pos_heights[pos] * 0.5
        ax.axhline(y=y_center - _half, color='#cccccc', linewidth=0.5, alpha=0.4, zorder=1)
        ax.axhline(y=y_center + _half, color='#cccccc', linewidth=0.5, alpha=0.4, zorder=1)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    # ── Second Y-axis (right): position numbers ──
    ax2 = ax.twinx()
    ax2.set_ylim(ax.get_ylim())
    _right_ticks = [pos_to_y[p] for p in positions]
    _right_labels = [str(p) for p in positions]
    ax2.set_yticks(_right_ticks)
    ax2.set_yticklabels(_right_labels, fontsize=_heading_fontsize,
                        fontweight='bold', color='black', alpha=0.7,
                        ha='left')
    # Small positive pad pushes labels just outside the right spine, into the
    # gap between the plot and the colorbar column.
    ax2.tick_params(axis='y', length=0, pad=20)
    ax2.spines['top'].set_visible(False)
    ax2.spines['right'].set_visible(False)

    # Colourbar — rendered in dedicated ax_cb column via shared helper
    _vmin = getattr(myoptions, 'cmap_vmin', -11)
    _vmax = getattr(myoptions, 'cmap_vmax', 11)
    _matrix_name = getattr(myoptions, 'matrix', 'BLOSUM80')
    # Same conditional label as mutation_scatter_plot: mention synonymous
    # codon changes only when they are actually present in the data.
    if getattr(myoptions, 'include_synonymous', False) or not getattr(myoptions, 'aminoacids', False):
        _cb_label = f'{_matrix_name} score values (synonymous codon changes shown in dark green)'
    else:
        _cb_label = f'{_matrix_name} score values'
    setup_matplotlib_colorbar(
        fig, ax_cb, norm, cmap, colors, _vmin, _vmax,
        label=_cb_label, alpha=0.5, label_fontsize=_heading_fontsize,
    )

    # Size legend — rendered in dedicated ax_leg column
    ax_leg.set_axis_off()
    _legend_freqs = [0.01, 0.1, 0.5, 1.0]
    _legend_sizes = []
    for f in _legend_freqs:
        if getattr(myoptions, 'linear_circle_size', False):
            raw = f * TIMELINE_CIRCLE_SCALE
        else:
            raw = (f ** 0.5) * TIMELINE_CIRCLE_SCALE
        _legend_sizes.append(max(TIMELINE_MIN_SIZE, min(TIMELINE_MAX_SIZE, raw)))
    for f, s in zip(_legend_freqs, _legend_sizes):
        ax_leg.scatter([], [], s=s, c='gray', alpha=0.5, edgecolors='#333',
                       linewidths=0.5, label=f'{f:.1%}')
    ax_leg.legend(loc='center', fontsize=7, frameon=False, title='Circle size',
                  title_fontsize=8,
                  labelspacing=3.5, handletextpad=1.5, borderpad=1.2,
                  scatterpoints=1)

    # NOTE: do NOT call tight_layout() here — it collapses the colorbar axes.
    # subplots_adjust is set at figure creation to match mutation_scatter_plot.

    # ── Auto-fix Y-axis label overlaps ──
    # Draw the canvas to compute actual text bounding boxes, then check
    # for overlapping Y-axis tick labels and increase figure height if needed.
    for _attempt in range(4):
        fig.canvas.draw()
        renderer = fig.canvas.get_renderer()
        ytick_bboxes = []
        for label in ax.get_yticklabels():
            bb = label.get_window_extent(renderer=renderer)
            ytick_bboxes.append(bb)
        # Check for pairwise overlaps among consecutive labels (sorted by y)
        sorted_bbs = sorted(ytick_bboxes, key=lambda b: b.y0)
        has_overlap = False
        max_overlap_px = 0.0
        for i in range(len(sorted_bbs) - 1):
            gap = sorted_bbs[i + 1].y0 - sorted_bbs[i].y1
            if gap < 0:
                has_overlap = True
                max_overlap_px = max(max_overlap_px, -gap)
        if not has_overlap:
            break
        # Increase figure height proportionally to fix the overlap
        old_h = fig.get_size_inches()[1]
        # Scale up by the ratio needed to eliminate the worst overlap
        # Each label pair needs (overlap + 2px padding) more space
        total_label_span = sorted_bbs[-1].y1 - sorted_bbs[0].y0
        if total_label_span > 0:
            scale = 1.0 + (max_overlap_px + 2.0) * len(sorted_bbs) / total_label_span
        else:
            scale = 1.3
        scale = min(scale, 1.5)  # cap per-iteration growth
        fig.set_size_inches(fig.get_size_inches()[0], old_h * scale)
        fig.subplots_adjust(left=0.06, right=0.86, top=0.94, bottom=0.15)

    # ── mplcursors hover support ──
    if not getattr(myoptions, 'disable_showing_mplcursors', False):
        try:
            import mplcursors
            _cursor = mplcursors.cursor(scatter, hover=True)

            @_cursor.connect("add")
            def _on_add(sel):
                """Display mutation details in the hover annotation."""
                idx = sel.index
                if 0 <= idx < len(labels):
                    sel.annotation.set_text(labels[idx])
                else:
                    sel.annotation.set_text("?")
        except ImportError:
            pass  # mplcursors not installed
        except Exception:  # pylint: disable=broad-exception-caught
            pass  # graceful fallback for non-interactive backends

    # Embed the git version string as a small footnote so that PNG and PDF
    # files can always be traced back to the exact code version that produced
    # them.  Matches the pattern in mutation_scatter_plot.render_matplotlib.
    from . import VERSION, _GIT_VERSION
    _version_label = (
        f"mutation_timeline_plot v{VERSION}  git:{_GIT_VERSION}  "
        f"{os.path.basename(outfile_prefix)}"
    )
    fig.text(
        0.99, 0.002, _version_label,
        ha='right', va='bottom',
        fontsize=6, color='#808080',
        transform=fig.transFigure,
    )

    # Save outputs
    for ext in ('png', 'pdf'):
        outpath = f"{outfile_prefix}.{ext}"
        _dpi = getattr(myoptions, 'dpi', 600)
        fig.savefig(outpath, dpi=_dpi, bbox_inches='tight', facecolor='white')
        print(f"Info: Saved {outpath}")

    plt.close(fig)


# ── Bokeh rendering ──────────────────────────────────────────────────────

def render_timeline_bokeh(
    data: TimelineData,
    myoptions: typing.Any,
    norm: typing.Any,
    cmap: typing.Any,
    colors: typing.Any,
    outfile_prefix: str,
) -> None:
    """Render interactive Bokeh HTML timeline.

    Produces an HTML file at ``{outfile_prefix}.html`` with:

    * Full-page responsive layout (``sizing_mode='stretch_width'``).
    * Hover tooltips showing: Mutation, Position, Month, Original Codon (AA),
      New Codon (AA), BLOSUM score, and Frequency.
    * Percentage labels next to each circle.
    * Y-axis tick labels showing mutation names.
    * Band borders separating position groups.

    Parameters
    ----------
    data : TimelineData
        Collected timeline data points.
    myoptions : argparse.Namespace
        CLI options.
    outfile_prefix : str
        Output path prefix (without extension).
    """
    if not data.points:
        print("Warning: No data points to render in Bokeh timeline")
        return

    try:
        from bokeh.plotting import figure, output_file, save
        from bokeh.models import ColumnDataSource, HoverTool
    except ImportError:
        print("Warning: bokeh not available, skipping Bokeh timeline output")
        return

    months = data.months
    positions = data.positions

    # Shared layout: per-position heights, spread, y-mapping, grouping, label offsets
    pos_heights, pos_spread, pos_to_y, grouped, label_offsets = _prepare_layout(data, myoptions)

    x_vals: list[float] = []
    y_vals: list[float] = []
    sizes_px: list[float] = []
    colors_hex: list[str] = []
    hover_labels: list[str] = []
    hover_freqs: list[str] = []
    hover_months: list[str] = []
    hover_codons: list[str] = []
    hover_positions: list[str] = []
    hover_ref_codons: list[str] = []
    hover_new_codons: list[str] = []
    hover_scores: list[str] = []

    for (month, pos), pts in grouped.items():
        x = _month_to_float(month, months)
        y_base = pos_to_y[pos]
        _pos_offsets = label_offsets.get(pos, {})

        for pt in pts:
            y_offset = _pos_offsets.get(_slot_key(pt), 0.0)

            x_vals.append(x)
            y_vals.append(y_base + y_offset)

            raw_size = float(pt.frequency) * TIMELINE_CIRCLE_SCALE
            if getattr(myoptions, 'linear_circle_size', False):
                sizes_px.append(max(3, min(25, raw_size ** 0.5 * 2)))
            else:
                sizes_px.append(max(3, min(25, (raw_size ** 0.5) ** 0.5 * 4)))
            colors_hex.append(pt.color)
            hover_labels.append(pt.label)
            hover_freqs.append(f"{float(pt.frequency):.6f}")
            hover_months.append(month)
            hover_codons.append(f"{pt.ref_codon}{pt.position}{pt.mutant_codon}")
            hover_positions.append(str(pt.position))
            hover_ref_codons.append(f"{pt.ref_codon} ({pt.ref_aa})")
            hover_new_codons.append(f"{pt.mutant_codon} ({pt.mutant_aa})")
            hover_scores.append(str(pt.score))

    source = ColumnDataSource(data={
        "x": x_vals,
        "y": y_vals,
        "size": sizes_px,
        "color": colors_hex,
        "label": hover_labels,
        "freq": hover_freqs,
        "month": hover_months,
        "codon": hover_codons,
        "position": hover_positions,
        "ref_codon": hover_ref_codons,
        "new_codon": hover_new_codons,
        "score": hover_scores,
    })

    title = getattr(myoptions, 'title', '') or f"Mutation Timeline ({len(months)} months, {len(positions)} positions)"

    n_pos = len(positions)
    _total_height = sum(pos_heights.get(p, 2.0) for p in positions)
    # y-limits: half-band padding around first and last position
    _first_half = pos_heights.get(positions[0], 2.0) * 0.5 if positions else 1.0
    _last_half = pos_heights.get(positions[-1], 2.0) * 0.5 if positions else 1.0
    _y_min = pos_to_y.get(positions[0], 0.0) - _first_half - 0.2 if positions else -0.2
    _y_max = pos_to_y.get(positions[-1], 0.0) + _last_half + 0.2 if positions else 0.2
    bokeh_fig = figure(
        title=title,
        width=2000,
        height=max(600, int(_total_height * 27 + 200)),
        x_range=(-0.5, len(months) - 0.5),
        y_range=(_y_min, _y_max),
        tools="pan,wheel_zoom,box_zoom,reset,save",
        sizing_mode='stretch_width',
    )

    bokeh_fig.scatter(
        'x', 'y',
        source=source,
        size='size',
        color='color',
        alpha=0.8,
        line_color='#333333',
        line_width=0.5,
    )

    _matrix_name = getattr(myoptions, 'matrix', 'BLOSUM80')
    hover = HoverTool(tooltips=[
        ("Mutation", "@label"),
        ("Position", "@position"),
        ("Month", "@month"),
        ("Original Codon", "@ref_codon"),
        ("New Codon", "@new_codon"),
        (f"{_matrix_name} score", "@score"),
        ("Frequency", "@freq"),
    ])
    bokeh_fig.add_tools(hover)

    # Axis labels
    bokeh_fig.xaxis.ticker = list(range(len(months)))
    bokeh_fig.xaxis.major_label_overrides = dict(enumerate(months))
    bokeh_fig.xaxis.major_label_orientation = 0.785  # 45 degrees
    bokeh_fig.xaxis.axis_label = "Month"

    # Y-axis: use merged mutation labels as tick labels
    _bokeh_tick_y, _bokeh_tick_labels = _compute_ytick_labels(
        positions, pos_to_y, grouped, pos_spread,
    )
    bokeh_fig.yaxis.ticker = _bokeh_tick_y
    bokeh_fig.yaxis.major_label_overrides = dict(zip(_bokeh_tick_y, _bokeh_tick_labels))
    bokeh_fig.yaxis.axis_label = "AA Position"

    # Grid
    bokeh_fig.xgrid.grid_line_alpha = 0.3
    bokeh_fig.ygrid.grid_line_alpha = 0.0

    # ── Bokeh Colorbar (shared helper, same approach as mutation_scatter_plot) ──
    _bk_vmin = getattr(myoptions, 'cmap_vmin', -11)
    _bk_vmax = getattr(myoptions, 'cmap_vmax', 11)
    if getattr(myoptions, 'include_synonymous', False) or not getattr(myoptions, 'aminoacids', False):
        _matrix_label = f"{_matrix_name} score values (synonymous codon changes shown in dark green)"
    else:
        _matrix_label = f"{_matrix_name} score values"
    try:
        add_bokeh_colorbar(
            bokeh_fig, norm, cmap, colors, _bk_vmin, _bk_vmax,
            alpha=0.8, label=_matrix_label,
        )
    except Exception:  # pylint: disable=broad-exception-caught
        pass

    # Percentage + mutation labels next to circles
    try:
        from bokeh.models import LabelSet
        annot_texts = [
            f"{_format_pct(float(f))}\n{lbl}\n{cod}"
            for f, lbl, cod in zip(hover_freqs, hover_labels, hover_codons)
        ]

        pct_source = ColumnDataSource(data={
            "x": x_vals, "y": y_vals, "text": annot_texts,
        })
        pct_labels = LabelSet(
            x='x', y='y', text='text', source=pct_source,
            text_font_size='7pt', text_color='black',
            x_offset=4, y_offset=3,
        )
        bokeh_fig.add_layout(pct_labels)
    except Exception:  # pylint: disable=broad-exception-caught
        pass

    # Embed the git version string below the x-axis, mirroring render_bokeh
    # in mutation_scatter_plot.
    from . import VERSION, _GIT_VERSION
    import bokeh.models
    _version_label_text = f"mutation_timeline_plot v{VERSION}  git:{_GIT_VERSION}"
    _version_label = bokeh.models.Title(
        text=_version_label_text,
        align='right',
        text_font_size='9pt', text_color='#808080',
        text_font_style='normal',
    )
    bokeh_fig.add_layout(_version_label, 'below')

    # Save HTML — embed version in browser tab title
    _html_title = (
        f"{os.path.basename(outfile_prefix)} "
        f"| mutation_timeline_plot v{VERSION} git:{_GIT_VERSION}"
    )
    html_path = f"{outfile_prefix}.html"
    output_file(html_path, title=_html_title)
    save(bokeh_fig)
    print(f"Info: Saved {html_path}")

    # Post-process: set/replace <title> tag unconditionally so the version
    # string is visible in the browser tab regardless of Bokeh version.
    if os.path.exists(html_path):
        with open(html_path, 'r', encoding='utf-8') as _fh:
            _html_content = _fh.read()
        _title_tag = f'<title>{_html_title}</title>'
        if '<title>' in _html_content:
            _html_content = re.sub(
                r'<title>[^<]*</title>', _title_tag, _html_content, count=1,
            )
        else:
            _html_content = _html_content.replace(
                '<head>', f'<head>\n  {_title_tag}', 1,
            )
        with open(html_path, 'w', encoding='utf-8') as _fh:
            _fh.write(_html_content)
