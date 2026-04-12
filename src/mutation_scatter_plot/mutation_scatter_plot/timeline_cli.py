# -*- coding: utf-8 -*-
# vim: set fileencoding=utf-8 ts=4 sw=4 expandtab :
"""Command-line interface for mutation_timeline_plot."""

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
    return myparser


def main():
    """Entry point for mutation_timeline_plot CLI."""
    myparser = build_option_parser()
    myoptions = myparser.parse_args()

    # Start profiler
    PROFILER.start()
    PROFILER.mark_phase_start("timeline_init")

    print(f"mutation_timeline_plot: scanning {myoptions.input_dir}")
    print(f"  pattern:   {myoptions.pattern}")
    print(f"  positions: {' '.join(myoptions.positions)}")
    print(f"  matrix:    {myoptions.matrix}")
    print(f"  colormap:  {myoptions.colormap}")

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

    # Append .timeline before load_matrix() adds .BLOSUM80.area_scaling.coolwarm_r
    # so the final order is: prefix.timeline.BLOSUM80.area_scaling.coolwarm_r.ext
    myoptions.outfile_prefix = myoptions.outfile_prefix + '.timeline'

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

    # Render matplotlib (PNG + PDF)
    PROFILER.mark_phase_start("render_matplotlib")
    render_timeline_matplotlib(data, myoptions, _norm, _cmap, _outfile_prefix)

    summary = PROFILER.pop_phase_summary()
    if summary:
        print(summary)

    # Render Bokeh (HTML)
    PROFILER.mark_phase_start("render_bokeh")
    render_timeline_bokeh(data, myoptions, _outfile_prefix)

    summary = PROFILER.pop_phase_summary()
    if summary:
        print(summary)

    print("mutation_timeline_plot: done")


if __name__ == '__main__':
    main()
