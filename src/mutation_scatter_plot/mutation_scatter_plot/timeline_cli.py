# -*- coding: utf-8 -*-
# vim: set fileencoding=utf-8 ts=4 sw=4 expandtab :
"""Command-line interface for mutation_timeline_plot.

This module provides the ``mutation_timeline_plot`` console script entry
point.  It parses command-line arguments, scans a directory of per-month
``.frequencies.tsv`` files, collects mutation data at user-specified amino
acid positions, and renders both static (matplotlib PNG/PDF) and
interactive (Bokeh HTML) timeline scatter plots.

Usage
-----
::

    mutation_timeline_plot --dir <input_dir> --positions 498[RHQ] N501Y 484 \
                           --threshold 0.001
"""

import argparse
import sys

from ..profiler import PROFILER
from .core import get_colormap, load_matrix
from .timeline import (
    collect_timeline_data,
    filter_timeline_data,
    infer_common_prefix,
    parse_positions,
    recolor_timeline_data,
    render_timeline_bokeh,
    render_timeline_matplotlib,
    scan_directory,
)


class NoWrapFormatter(argparse.RawDescriptionHelpFormatter, argparse.ArgumentDefaultsHelpFormatter):
    """Help formatter that does not wrap long lines, preserving URLs."""

    def _split_lines(self, text, width):
        """Preserve explicit newlines in help text (e.g. for URL formatting)."""
        return text.splitlines()


def build_option_parser():
    """Build and return the command-line argument parser."""
    myparser = argparse.ArgumentParser(
        description="mutation_timeline_plot — timeline scatter plot of mutation frequencies across monthly datasets",
        formatter_class=NoWrapFormatter,
    )
    myparser.add_argument(
        "--dir", action="store", type=str, dest="input_dir",
        required=True,
        help="Directory containing monthly .frequencies.tsv files to scan",
    )
    myparser.add_argument(
        "--positions", nargs='+', type=str, dest="positions",
        required=True,
        help="Positions to track. Accepts numeric (614), mutation labels (D614G), or ranges (480-530)",
    )
    myparser.add_argument(
        "--pattern", action="store", type=str, dest="pattern",
        default='*.frequencies.tsv',
        help="Glob pattern for input TSV files",
    )
    myparser.add_argument(
        "--aminoacids", action="store_true", dest="aminoacids",
        default=False,
        help="Use amino acid mode instead of codon mode",
    )
    myparser.add_argument(
        "--colormap", nargs='+', type=str, dest="colormaps",
        default=['coolwarm_r'],
        help="Colormap name(s) for BLOSUM scoring."
             " Multiple names can be given to render all in one pass"
             " (e.g. --colormap coolwarm_r amino_acid_changes)",
    )
    myparser.add_argument(
        "--matrix", action="store", type=str, dest="matrix",
        default='BLOSUM80',
        help="Substitution scoring matrix",
    )
    myparser.add_argument(
        "--matrix-file", action="store", type=str, dest="matrix_file",
        default='',
        help="Path to a custom substitution matrix file",
    )
    myparser.add_argument(
        "--outfile-prefix", action="store", type=str, dest="outfile_prefix",
        default='',
        help="Output file prefix. If omitted, inferred from common root of input filenames",
    )
    myparser.add_argument(
        "--title", action="store", type=str, dest="title",
        default='',
        help="Figure title (auto-generated if empty)",
    )
    myparser.add_argument(
        "--offset", action="store", type=int, dest="offset",
        default=0,
        help="Position offset to add to amino acid positions from TSV",
    )
    myparser.add_argument(
        "--threshold", action="store", type=float, dest="threshold",
        default=0.0,
        help="Minimum frequency to display a data point",
    )
    myparser.add_argument(
        "--include-unknown-month", action="store_true",
        dest="include_unknown_month", default=False,
        help="Include files with month '00' (unknown month)",
    )
    myparser.add_argument(
        "--debug", action="store_true", dest="debug",
        default=False,
        help="Enable debug output",
    )
    myparser.add_argument(
        "--linear-circle-size", action="store_true", dest="linear_circle_size",
        default=False,
        help="Use linear (not area-based) circle scaling",
    )
    myparser.add_argument(
        "--both-scaling", action="store_true", dest="both_scaling",
        default=False,
        help="Render both area-scaling and linear-scaling variants in one pass."
             " Overrides --linear-circle-size.",
    )
    myparser.add_argument(
        "--spread-colormap-over-virtual-matrix", action="store_true",
        dest="spread_colormap_virtual_matrix", default=False,
        help="Spread colormap over full theoretical matrix range instead of actual data range",
    )
    myparser.add_argument(
        "--dpi", action="store", type=int, dest="dpi", default=600,
        help="DPI resolution for PNG/PDF output images",
    )
    myparser.add_argument(
        "--include-synonymous", action="store_true",
        dest="include_synonymous", default=False,
        help="Include synonymous changes in --aminoacids output."
             " In codon mode they are always shown [default: False]",
    )
    myparser.add_argument(
        "--disable-showing-bokeh", action="store_true",
        dest="disable_showing_bokeh", default=False,
        help="Skip generating the interactive Bokeh HTML output [default: False]",
    )
    myparser.add_argument(
        "--disable-showing-mplcursors", action="store_true",
        dest="disable_showing_mplcursors", default=False,
        help="Disable interactive mplcursors hover on matplotlib output [default: False]",
    )
    myparser.add_argument(
        "--matplotlib-backend", action="store", type=str, dest="matplotlib_backend", default='',
        help="Matplotlib backend (agg, wxpython, pyqt5, etc.). "
             "Set MPLBACKEND=agg to prevent interactive windows [default: unset]",
    )
    myparser.add_argument(
        "--band-spacing-factor", action="store", type=float,
        dest="band_spacing_factor", default=1.0,
        help="Multiplier for vertical spacing between Y-axis position bands."
             " Values >1 increase spacing (useful for dense datasets)."
             " [default: 1.0]",
    )
    myparser.add_argument(
        "--positions-per-page", action="store", type=int,
        dest="positions_per_page", default=10,
        help="Maximum number of Y-axis position bands per output page."
             " When the total number of positions exceeds this value the"
             " output is split into numbered pages (e.g. prefix.page1.png)."
             " Set to 0 to disable splitting."
             " [default: 10]",
    )
    return myparser


def main():
    """Entry point for the mutation_timeline_plot CLI.

    Workflow:

    1. Parse command-line arguments via :func:`build_option_parser`.
    2. Scan the input directory for monthly ``.frequencies.tsv`` files.
    3. Infer a common filename prefix for output naming.
    4. Load the BLOSUM substitution matrix.
    5. Collect timeline data once (scores are matrix-dependent, not
       colormap-dependent).
    6. For each colormap × scaling combination:
       a. Derive the colormap and recolour the data points.
       b. Render static matplotlib output (PNG + PDF).
       c. Render interactive Bokeh output (HTML).
    """
    myparser = build_option_parser()
    myoptions = myparser.parse_args()

    # Apply matplotlib backend before any plotting imports occur
    if myoptions.matplotlib_backend:
        import matplotlib
        matplotlib.use(myoptions.matplotlib_backend)

    # Start profiler
    PROFILER.start()
    PROFILER.mark_phase_start("timeline_init")

    # Build the list of colormaps and scaling modes to render
    colormaps = myoptions.colormaps
    if myoptions.both_scaling:
        scaling_modes = [False, True]  # area_scaling, then linear_scaling
    else:
        scaling_modes = [myoptions.linear_circle_size]

    n_combos = len(colormaps) * len(scaling_modes)

    print(f"mutation_timeline_plot: scanning {myoptions.input_dir}")
    print(f"  pattern:   {myoptions.pattern}")
    print(f"  positions: {' '.join(myoptions.positions)}")
    print(f"  matrix:    {myoptions.matrix}")
    print(f"  colormaps: {' '.join(colormaps)}")
    print(f"  scaling:   {'both' if myoptions.both_scaling else ('linear' if myoptions.linear_circle_size else 'area')}")
    print(f"  variants:  {n_combos} (= {len(colormaps)} colormap(s) × {len(scaling_modes)} scaling(s))")
    print(f"  dpi:       {myoptions.dpi}")

    # Parse position specifications
    specs = parse_positions(myoptions.positions)
    if not specs:
        print("Error: No valid positions specified", file=sys.stderr)
        sys.exit(1)
    print(f"  parsed:    {', '.join(str(s) for s in specs)}")

    # Scan directory for TSV files
    files = scan_directory(
        myoptions.input_dir,
        myoptions.pattern,
        myoptions.include_unknown_month,
    )
    if not files:
        print(f"Error: No matching files found in {myoptions.input_dir}", file=sys.stderr)
        sys.exit(1)
    print(f"  files:     {len(files)} monthly datasets found")

    # Auto-infer output prefix from filenames if not specified
    if not myoptions.outfile_prefix:
        myoptions.outfile_prefix = infer_common_prefix(files, myoptions.input_dir)
        print(f"  prefix:    {myoptions.outfile_prefix} (auto-inferred)")
    else:
        print(f"  prefix:    {myoptions.outfile_prefix}")

    summary = PROFILER.pop_phase_summary()
    if summary:
        print(summary)

    # Append .timeline.aa/.codon to the base prefix (before matrix/scaling/colormap)
    _mode_label = 'aa' if myoptions.aminoacids else 'codon'
    _base_prefix = myoptions.outfile_prefix + '.timeline.' + _mode_label

    # Load matrix (once — scoring is matrix-dependent, not colormap-dependent)
    # Temporarily set colormap and scaling so load_matrix() constructs its
    # outfile_prefix — we won't use that prefix; we build our own per-combo.
    myoptions.outfile_prefix = _base_prefix
    myoptions.colormap = colormaps[0]
    PROFILER.mark_phase_start("matrix_load")
    _matrix, _matrix_name, _min_score, _max_score, _ = load_matrix(myoptions)

    # Initial colormap for data collection (first colormap in the list)
    myoptions.cmap_actual_vmin = -10
    myoptions.cmap_actual_vmax = 10
    _norm, _cmap, _colors = get_colormap(myoptions, colormaps[0])

    summary = PROFILER.pop_phase_summary()
    if summary:
        print(summary)

    # Collect timeline data ONCE (scores are matrix-dependent, colours will be
    # recomputed per colormap via recolor_timeline_data)
    PROFILER.mark_phase_start("data_collection")
    data = collect_timeline_data(files, specs, myoptions, _matrix, _norm, _colors)
    print(f"  data:      {len(data.points)} points across {len(data.months)} months, {len(data.positions)} positions")

    summary = PROFILER.pop_phase_summary()
    if summary:
        print(summary)

    if not data.points:
        print("Warning: No data points found for the specified positions.")
        sys.exit(0)

    # Compute actual score range from data (shared across all combos)
    _actual_scores = [pt.score for pt in data.points if pt.score != 12]
    _actual_vmin = min(_actual_scores) if _actual_scores else -11
    _actual_vmax = max(_actual_scores) if _actual_scores else 11
    _actual_bound = max(abs(_actual_vmin), abs(_actual_vmax))
    myoptions.cmap_actual_vmin = -_actual_bound
    myoptions.cmap_actual_vmax = _actual_bound
    print(f"  scores:    actual range [{_actual_vmin}, {_actual_vmax}], "
          f"colorbar range [{myoptions.cmap_actual_vmin}, {myoptions.cmap_actual_vmax}]")

    # ── Render each colormap × scaling combination ──
    _combo_idx = 0

    # Pagination: split positions into pages of N bands each
    _per_page = myoptions.positions_per_page
    if 0 < _per_page < len(data.positions):
        _chunks = [
            data.positions[i:i + _per_page]
            for i in range(0, len(data.positions), _per_page)
        ]
    else:
        _chunks = [data.positions]
    _n_pages = len(_chunks)
    if _n_pages > 1:
        print(f"  pages:     {_n_pages} (= {len(data.positions)} positions "
              f"/ {_per_page} per page)")

    for cmap_name in colormaps:
        # Derive colormap with data-driven bounds
        # Clear stale cmap_vmin/vmax so get_colormap picks up fresh values
        if hasattr(myoptions, 'cmap_vmin'):
            delattr(myoptions, 'cmap_vmin')
        if hasattr(myoptions, 'cmap_vmax'):
            delattr(myoptions, 'cmap_vmax')
        myoptions.colormap = cmap_name
        _norm, _cmap, _colors = get_colormap(myoptions, cmap_name)

        # Recolour all data points for this colormap
        recolor_timeline_data(data, myoptions, _norm, _colors)

        for linear in scaling_modes:
            _combo_idx += 1
            myoptions.linear_circle_size = linear
            _scaling_suffix = 'linear_scaling' if linear else 'area_scaling'
            _outfile_prefix = f"{_base_prefix}.{_matrix_name}.{_scaling_suffix}.{cmap_name}"
            print(f"Info: _outfile_prefix={_outfile_prefix}")

            if n_combos > 1:
                print(f"  ── variant {_combo_idx}/{n_combos}: "
                      f"{cmap_name} + {_scaling_suffix} ──")

            for page_idx, page_positions in enumerate(_chunks, 1):
                if _n_pages > 1:
                    page_data = filter_timeline_data(data, page_positions)
                    _page_prefix = f"{_outfile_prefix}.page{page_idx}"
                    print(f"  ── page {page_idx}/{_n_pages}: "
                          f"positions {page_positions[0]}..{page_positions[-1]} "
                          f"({len(page_data.points)} points) ──")
                else:
                    page_data = data
                    _page_prefix = _outfile_prefix

                # Render matplotlib (PNG + PDF + JPG)
                PROFILER.mark_phase_start("render_matplotlib")
                render_timeline_matplotlib(
                    page_data, myoptions, _norm, _cmap, _colors, _page_prefix,
                )

                summary = PROFILER.pop_phase_summary()
                if summary:
                    print(summary)

                # Render Bokeh (HTML)
                if not myoptions.disable_showing_bokeh:
                    PROFILER.mark_phase_start("render_bokeh")
                    render_timeline_bokeh(
                        page_data, myoptions, _norm, _cmap, _colors, _page_prefix,
                    )

                summary = PROFILER.pop_phase_summary()
                if summary:
                    print(summary)

    print("mutation_timeline_plot: done")


if __name__ == '__main__':
    main()
