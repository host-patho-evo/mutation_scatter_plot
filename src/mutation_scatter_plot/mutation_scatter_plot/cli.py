# -*- coding: utf-8 -*-
# vim: set fileencoding=utf-8 ts=4 sw=4 expandtab :
"""Command-line interface for mutation_scatter_plot."""

import argparse
import os
import sys

from ..profiler import PROFILER
from . import (VERSION, build_frequency_tables, collect_scatter_data,
               load_and_clean_dataframe, load_matrix, render_bokeh,
               render_matplotlib, setup_matplotlib_figure)


class NoWrapFormatter(argparse.RawDescriptionHelpFormatter, argparse.ArgumentDefaultsHelpFormatter):
    """Help formatter that does not wrap long lines, preserving URLs."""

    def _split_lines(self, text, width):
        return text.splitlines()


def build_option_parser():
    """Build and return the command-line argument parser."""
    myparser = argparse.ArgumentParser(
        description=f"mutation_scatter_plot version {VERSION}",
        formatter_class=NoWrapFormatter,
    )
    myparser.add_argument(
        "--tsv", action="store", type=str, dest="tsv_file_path",
        default='',
        help="Path to a TAB separated file with 5-columns: ['position',"
             " 'original_aa', 'mutant_aa', 'frequency', 'original_codon',"
             " 'mutant_codon', 'coverage_per_codon'] or more up to 11-columns"
             " ['padded_position', 'position', 'original_aa', 'mutant_aa',"
             " 'frequency', 'original_codon', 'mutant_codon',"
             " 'observed_codon_count', 'total_codons_per_site',"
             " 'frequency_parent', 'frequency_selected'] as we kept extending"
             " the file format",
    )
    myparser.add_argument(
        "--column", action="store", type=str,
        dest="column_with_frequencies", default='frequency',
        help="Name of a column in input TSV to be used for rendering"
             " frequencies [frequency]",
    )
    myparser.add_argument(
        "--outfile-prefix", action="store", type=str,
        dest="outfile_prefix", default='',
        help="Output file prefix, eventually also with path. Output files"
             " will be PNG, JPG or PDF, HTML",
    )
    myparser.add_argument(
        "--offset", action="store", type=int, dest="offset", default=0,
        help="Define offset value to be added to every amino acid position"
             " in the first column of the TSV input file at runtime. This is"
             " to obtain real amino acid position within the protein. Normally"
             " protein starts at position 1. Check the first aa in TSV file"
             " and provide whatever number to be added to it to get the"
             " desired amino acid position in the full-length protein.",
    )
    myparser.add_argument(
        "--xmin", action="store", type=int, dest="xmin", default=0,
        help="Define minimum X-axis value. This should be the position in the padded alignment (using '-').",
    )
    myparser.add_argument(
        "--xmax", action="store", type=int, dest="xmax", default=0,
        help="Define maximum X-axis value. This should be the position in the padded alignment (using '-').",
    )
    myparser.add_argument(
        "--x-axis-bins", action="store", type=int, dest="xaxis_bins",
        default=0,
        help="Set number of bins (labels) on the X-axis. Could be used to"
             " override the amount of major ticks in a different way. Use 20"
             " for example. [default is 0]",
    )
    myparser.add_argument(
        "--x-axis-major-ticks-spacing", action="store", type=int,
        dest="xaxis_major_ticks_spacing", default=10,
        help="Set distance between the major ticks on X-axis",
    )
    myparser.add_argument(
        "--x-axis-minor-ticks-spacing", action="store", type=int,
        dest="xaxis_minor_ticks_spacing", default=5,
        help="Set distance between the minor ticks on X-axis",
    )
    myparser.add_argument(
        "--x-axis-label-start", action="store", type=int,
        dest="xaxis_label_start", default=0,
        help="Set first label value on the X-axis",
    )
    myparser.add_argument(
        "--aminoacids", action="store_true", dest="aminoacids", default=False,
        help="Draw chart with amino acid residues on Y-axis instead of"
             " codons. [default is False]",
    )
    myparser.add_argument(
        "--show-STOP", action="store_true", dest="showstop", default=False,
        help="Include STOP codons or '*' in charts on Y-axis."
             " [default is False]",
    )
    myparser.add_argument(
        "--show-INS", action="store_true", dest="showins", default=False,
        help="Include INS in charts on Y-axis. [default is False]",
    )
    myparser.add_argument(
        "--show-DEL", action="store_true", dest="showdel", default=False,
        help="Include DEL in charts on Y-axis. [default is False]",
    )
    myparser.add_argument(
        "--show-X", action="store_true", dest="showx", default=False,
        help="Include X in charts on Y-axis in --aminoacids mode."
             " [default is False]",
    )
    myparser.add_argument(
        "--enable-colorbar", action="store_true", dest="colorbar",
        default=False,
        help="Enable colorbar [is Disabled by default]",
    )
    myparser.add_argument(
        "--disable-short-legend", action="store_false", dest="shortlegend",
        default=True,
        help="Disable short legend in charts on X-axis."
             " [is Enabled by default]",
    )
    myparser.add_argument(
        "--include-synonymous", action="store_true",
        dest="include_synonymous", default=False,
        help="Include synonymous changes in --aminoacids output as green"
             " diamonds. In codon output they are always shown. Supplying this"
             " flag causes the total-frequency top bar heights in --aminoacids"
             " mode to mathematically mirror the taller bar heights natively"
             " reached in codon mode."
             " [default is False]",
    )
    myparser.add_argument(
        "--threshold", action="store", type=float, dest="threshold",
        default=0.001,
        help="Define minimum frequency threshold to display a pictogram in"
             " the output. For codon mode use 0.001 and for aa mode use 0.01."
             " [default: 0.001]",
    )
    myparser.add_argument(
        "--title", action="store", type=str, dest="title", default='',
        help="Set title for the figures, by default trailing"
             " '.frequencies.tsv' is stripped from the end of the input"
             " filename",
    )
    myparser.add_argument(
        "--disable-2nd-Y-axis", action="store_true",
        dest="disable_2nd_Y_axis", default=False,
        help="Disable 2nd Y-axis on the right side of the plot."
             " [default is False]",
    )
    myparser.add_argument(
        "--disable-showing-bokeh", action="store_true",
        dest="disable_showing_bokeh", default=False,
        help="Disable calling a GUI or text-based browser TO show the rendered Bokeh HTML+Javascript."
             " [default is False]",
    )
    myparser.add_argument(
        "--disable-showing-mplcursors", action="store_true",
        dest="disable_showing_mplcursors", default=False,
        help="Disable rendering an interactive window WITH hover() events recognized by mplcursors."
             " [default is False]",
    )
    myparser.add_argument(
        "--disable-padded-x-axis", action="store_true",
        dest="disable_padded_x_axis", default=False,
        help="Plot exclusively on natural 'position' instead of 'padded_position'."
             " Makes the X-axis narrower if insertions are present. Automatically toggles"
             " ON if no INS events exist within the [xmin, xmax] viewing window."
             " [default is False]",
    )
    myparser.add_argument(
        "--legend", action="store_true", dest="legend", default=False,
        help="Draw legend chart. [default is False]",
    )
    myparser.add_argument(
        "--matrix", action="store", type=str, dest="matrix",
        default="BLOSUM80",
        help="BLOSUM matrix: BLOSUM45,BLOSUM50,BLOSUM62,BLOSUM80,BLOSUM90"
             " [default is BLOSUM80]",
    )
    myparser.add_argument(
        "--matrix-file", action="store", type=str, dest="matrix_file",
        default=None,
        help="Matrix file compatible with BLOSUM matrices, e.g."
             " prot26050-sup-0004-Supinfo03.sm from"
             " https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8641535/bin/"
             "NIHMS1664401-supplement-supinfo.rar"
             " if you do not like default BLOSUM62",
    )
    myparser.add_argument(
        "--colormap", action="store", type=str, dest="colormap",
        default='amino_acid_changes',
        help="Pick a colormap recognized by matplotlib."
             " See https://i.sstatic.net/cmk1J.png"
             " [default is amino_acid_changes; coolwarm_r and seismic_r are"
             " great alternatives]",
    )
    myparser.add_argument(
        "--dpi", action="store", type=int, dest="dpi", default=600,
        help="DPI resolution for images",
    )
    myparser.add_argument(
        "--backend", action="store", type=str, dest="backend", default='',
        help="Matplotlib backend to render resulting figures: agg, wxpython,"
             " pyqt5, pyqt6, pycairo, cairocffi [default: unset]\n"
             "To disable Matplotlib interactive window being raised up you"
             " can set MPLBACKEND=agg env variable.",
    )
    myparser.add_argument(
        "--debug", action="store", type=int, dest="debug", default=0,
        help="Set debug to some real number",
    )
    myparser.add_argument(
        "--disable-bokeh-sqrt-size", action="store_false",
        dest="bokeh_sqrt_size", default=True,
        help="Disable sqrt scaling for Bokeh circle sizes; size (diameter)"
             " will be proportional to frequency, area proportional to"
             " frequency². By default sqrt scaling is on, matching the"
             " perceptual appearance of the matplotlib figure.",
    )
    myparser.add_argument(
        "--linear-circle-size", action="store_true",
        dest="linear_circle_size", default=False,
        help="Render circles with radius proportional to frequency in BOTH"
             " Matplotlib and Bokeh.  Matplotlib: s = freq\u00b2 \u00d7 5000 so that"
             " radius \u221d frequency.  Bokeh: diameter = freq \u00d7 100."
             " Implies --disable-bokeh-sqrt-size so that both figure types"
             " use the same linear-radius formula and remain visually in sync."
             " Matches the v0.2 Bokeh rendering. [default: False]",
    )
    myparser.add_argument(
        "--show-invisible-placeholder-dots", action="store_true",
        dest="show_invisible_placeholder_dots", default=False,
        help="Include below-threshold dots in the plot. [default: False]",
    )
    return myparser


def main():  # pylint: disable=too-many-locals
    """Entry point: parse options and run the full scatter-plot pipeline."""
    PROFILER.start()
    PROFILER.mark_phase_start("Phase 1: Input parsing and processing")
    myparser = build_option_parser()
    myoptions = myparser.parse_args()
    # --linear-circle-size implies --disable-bokeh-sqrt-size so that both the
    # Matplotlib and Bokeh figures use the same linear-radius formula and
    # remain visually in sync with each other.
    if myoptions.linear_circle_size:
        myoptions.bokeh_sqrt_size = False

    # Exception removed: We allow combining --disable-padded-x-axis and --show-INS
    # Programmatically, INS points will be masked from rendering below if disable-padded-x-axis=True.

    _matrix, _matrix_name, _min_theoretical_score, _max_theoretical_score, \
        _outfile_prefix = load_matrix(myoptions)

    # create a conversion dictionary
    _padded_position2position = {}

    # parse the .frequencies.tsv contents and fill-in the conversion dictionary
    _df, _padded_position2position = load_and_clean_dataframe(
        myoptions, myoptions.tsv_file_path, _padded_position2position,
    )

    if not myoptions.disable_padded_x_axis:
        _effective_xmin = myoptions.xmin if myoptions.xmin else min(_padded_position2position.keys())
        _effective_xmax = myoptions.xmax if myoptions.xmax else max(_padded_position2position.keys())
        _window = {k: v for k, v in _padded_position2position.items() if _effective_xmin <= k <= _effective_xmax}
        if _window:
            _differences = set(k - v for k, v in _window.items())
            # Algebraically, if _differences has only 1 unique continuous mapping
            # (e.g. k - v = offset 0 exclusively), no structural INS strings disrupt the map.
            # If multiple gaps overlap (e.g. padded 335 -> 333, padded 336 -> 333),
            # the subset produces varying offsets (2, 3), scaling len > 1 natively.
            if len(_differences) <= 1:
                myoptions.disable_padded_x_axis = True
                print(
                    "Info: Auto-detected no INS events altering coordinate width in the current region. "
                    "Switching to natural position X-axis.",
                    file=sys.stderr
                )

    print(f"Info: Writing into {_outfile_prefix}.actually_rendered.tsv")
    # Speedup 6: Export Decimal to string for reporting while keeping precision
    _df_to_save = _df.copy()
    _df_to_save[myoptions.column_with_frequencies] = _df_to_save[
        myoptions.column_with_frequencies
    ].apply(lambda x: f"{x:.6f}")
    _df_to_save.to_csv(
        f"{_outfile_prefix}.actually_rendered.tsv",
        sep='\t', header=True, index=False, float_format='{:.6f}'.format,
    )

    if '.frequencies.tsv' in myoptions.tsv_file_path:
        _count_filename = f"{_outfile_prefix}.count"
        if os.path.exists(_count_filename):
            try:
                with open(_count_filename, encoding="utf-8") as _aln_handle:
                    _aln_rows = _aln_handle.readline()
            except OSError:
                _aln_rows = '0'
        else:
            _aln_rows = '0'
    else:
        _aln_rows = '0'

    if not myoptions.title:
        _title_data = myoptions.tsv_file_path.replace('.frequencies.tsv', '')
    else:
        _title_data = myoptions.title

    print(f"Info: Title will be {_title_data}")

    # supplement the conversion dictionary with values from .frequencies.unchanged_codons.tsv
    _unchanged_tsv = myoptions.tsv_file_path.replace(
        '.frequencies.tsv', '.frequencies.unchanged_codons.tsv'
    )
    if os.path.exists(_unchanged_tsv):
        _df_frequencies_unchanged_codons, _padded_position2position = load_and_clean_dataframe(
            myoptions, _unchanged_tsv, _padded_position2position,
        )
        del _df_frequencies_unchanged_codons
    else:
        sys.stderr.write(
            f"Warning: File {_unchanged_tsv} not found, "
            f"skipping supplemental conversion data.{os.linesep}"
        )

    (
        _amino_acids, _codons_whitelist, _codons_whitelist2,
        _final_sorted_whitelist,
        _unique_aa_padded_positions, _unique_codon_padded_positions,
        _old_aa_table, _new_aa_table, _old_codon_table, _new_codon_table,
        _calculated_aa_offset, _padded_position2position,
    ) = build_frequency_tables(myoptions, _df, _padded_position2position)

    _figure, _ax1, _ax2, _ax3, _ax4, _xmin, _xmax = \
        setup_matplotlib_figure(
            myoptions,
            _title_data, _aln_rows, _matrix_name, _amino_acids,
            _codons_whitelist, _final_sorted_whitelist,
            _unique_aa_padded_positions, _unique_codon_padded_positions,
            _new_aa_table, _new_codon_table, _padded_position2position,
        )

    _table = _new_aa_table if myoptions.aminoacids else _new_codon_table

    (
        _norm, _cmap, _colors, _used_colors, _matrix_values,
        _mutations,
        _circles_bokeh, _circles_matplotlib, _markers, _dots,
        _hover_text_bokeh,
    ) = collect_scatter_data(
        myoptions,
        _df, _table, _outfile_prefix, _matrix, _amino_acids,
        _codons_whitelist2, _padded_position2position,
        _xmin, _xmax,
    )

    _xlabel = _ax1.get_xlabel()

    _prof_sum = PROFILER.pop_phase_summary()
    if _prof_sum:
        print()
        print(_prof_sum)

    if _circles_bokeh:
        render_bokeh(
            myoptions,
            _outfile_prefix, _xmin, _xmax, _amino_acids,
            _final_sorted_whitelist,
            _circles_bokeh, _mutations, _hover_text_bokeh,
            _title_data, _xlabel,
            _matrix_name, _colors, _norm, _cmap,
            _padded_position2position,
            show=not myoptions.disable_showing_bokeh,
        )

    render_matplotlib(
        myoptions,
        _figure, _ax1, _ax2, _ax3, _ax4, _outfile_prefix,
        _circles_matplotlib, _markers, _dots, _cmap, _norm, _colors,
        _matrix, _matrix_name,
        show=not myoptions.disable_showing_mplcursors,
    )


if __name__ == "__main__":
    main()

# vim:ts=4:sw=4:expandtab:smartindent
