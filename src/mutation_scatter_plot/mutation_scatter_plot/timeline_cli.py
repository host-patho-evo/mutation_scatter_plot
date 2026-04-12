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
    TimelineData,
    collect_timeline_data,
    infer_common_prefix,
    parse_positions,
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
        "--colormap", action="store", type=str, dest="colormap",
        default='coolwarm_r',
        help="Colormap name for BLOSUM scoring",
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
        "--backend", action="store", type=str, dest="backend", default='',
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
    return myparser


def main():
    """Entry point for the mutation_timeline_plot CLI.

    Workflow:

    1. Parse command-line arguments via :func:`build_option_parser`.
    2. Scan the input directory for monthly ``.frequencies.tsv`` files.
    3. Infer a common filename prefix for output naming.
    4. Load the BLOSUM substitution matrix and colormap.
    5. Collect timeline data for the specified positions.
    6. Compute the actual BLOSUM score range from the data and re-derive
       the colormap with data-driven bounds (unless
       ``--spread-colormap-virtual-matrix`` forces the full theoretical range).
    7. Render static matplotlib output (PNG + PDF).
    8. Render interactive Bokeh output (HTML).
    """
    myparser = build_option_parser()
    myoptions = myparser.parse_args()

    # Apply matplotlib backend before any plotting imports occur
    if myoptions.backend:
        import matplotlib
        matplotlib.use(myoptions.backend)

    # Start profiler
    PROFILER.start()
    PROFILER.mark_phase_start("timeline_init")

    print(f"mutation_timeline_plot: scanning {myoptions.input_dir}")
    print(f"  pattern:   {myoptions.pattern}")
    print(f"  positions: {' '.join(myoptions.positions)}")
    print(f"  matrix:    {myoptions.matrix}")
    print(f"  colormap:  {myoptions.colormap}")
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

    # Append .timeline and .aa/.codon before load_matrix() adds .BLOSUM80.area_scaling.coolwarm_r
    # so the final order is: prefix.timeline.aa.BLOSUM80.area_scaling.coolwarm_r.ext
    _mode_label = 'aa' if myoptions.aminoacids else 'codon'
    myoptions.outfile_prefix = myoptions.outfile_prefix + '.timeline.' + _mode_label

    # Load matrix
    PROFILER.mark_phase_start("matrix_load")
    _matrix, _matrix_name, _min_score, _max_score, _outfile_prefix = load_matrix(myoptions)

    # Get colormap
    # Set cmap bounds from actual data (will be refined by collect_timeline_data)
    myoptions.cmap_actual_vmin = -10
    myoptions.cmap_actual_vmax = 10
    _norm, _cmap, _colors = get_colormap(myoptions, myoptions.colormap)

    summary = PROFILER.pop_phase_summary()
    if summary:
        print(summary)

    # Collect timeline data
    PROFILER.mark_phase_start("data_collection")
    data = collect_timeline_data(files, specs, myoptions, _matrix, _norm, _colors)
    print(f"  data:      {len(data.points)} points across {len(data.months)} months, {len(data.positions)} positions")

    summary = PROFILER.pop_phase_summary()
    if summary:
        print(summary)

    if not data.points:
        print("Warning: No data points found for the specified positions.")
        sys.exit(0)

    # Compute actual score range from data and update colormap bounds
    # Exclude synonymous sentinel score (+12) from range computation
    _actual_scores = [pt.score for pt in data.points if pt.score != 12]
    _actual_vmin = min(_actual_scores) if _actual_scores else -11
    _actual_vmax = max(_actual_scores) if _actual_scores else 11
    # Make symmetric around zero
    _actual_bound = max(abs(_actual_vmin), abs(_actual_vmax))
    myoptions.cmap_actual_vmin = -_actual_bound
    myoptions.cmap_actual_vmax = _actual_bound
    print(f"  scores:    actual range [{_actual_vmin}, {_actual_vmax}], "
          f"colorbar range [{myoptions.cmap_actual_vmin}, {myoptions.cmap_actual_vmax}]")

    # Re-derive colormap with data-driven bounds (unless --spread-colormap-virtual-matrix)
    # Clear stale cmap_vmin/vmax so get_colormap picks up the updated cmap_actual_vmin/vmax
    if hasattr(myoptions, 'cmap_vmin'):
        delattr(myoptions, 'cmap_vmin')
    if hasattr(myoptions, 'cmap_vmax'):
        delattr(myoptions, 'cmap_vmax')
    _norm, _cmap, _colors = get_colormap(myoptions, myoptions.colormap)

    # Render matplotlib (PNG + PDF)
    PROFILER.mark_phase_start("render_matplotlib")
    render_timeline_matplotlib(data, myoptions, _norm, _cmap, _outfile_prefix)

    summary = PROFILER.pop_phase_summary()
    if summary:
        print(summary)

    # Render Bokeh (HTML)
    if not myoptions.disable_showing_bokeh:
        PROFILER.mark_phase_start("render_bokeh")
        render_timeline_bokeh(data, myoptions, _outfile_prefix)

    summary = PROFILER.pop_phase_summary()
    if summary:
        print(summary)

    print("mutation_timeline_plot: done")


if __name__ == '__main__':
    main()
