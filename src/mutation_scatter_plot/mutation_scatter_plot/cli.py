"""Command-line interface for mutation_scatter_plot."""

import os
from optparse import OptionParser, IndentedHelpFormatter
from . import (
    VERSION,
    load_matrix,
    load_and_clean_dataframe,
    build_frequency_tables,
    setup_matplotlib_figure,
    collect_scatter_data,
    render_bokeh,
    render_matplotlib,
)


class NoWrapFormatter(IndentedHelpFormatter):
    """Help formatter that does not wrap long lines, preserving URLs."""

    def format_description(self, description):
        """Return description with a trailing newline, unwrapped."""
        return f"{description}\n" if description else ""

    def format_option(self, option):
        """Format a single option without wrapping the help text."""
        result = []
        opts = self.option_strings[option]
        opt_width = self.help_position - self.current_indent - 2
        if len(opts) > opt_width:
            opts = f"{'':>{self.current_indent}}{opts}\n"
            indent_first = self.help_position
        else:
            opts = f"{'':>{self.current_indent}}{opts:<{opt_width}}  "
            indent_first = 0
        result.append(opts)
        if option.help:
            help_text = self.expand_default(option)
            # Do not wrap — emit the full help string as a single line
            result.append(f"{'':>{indent_first}}{help_text}\n")
        elif opts[-1] != "\n":
            result.append("\n")
        return "".join(result)


def build_option_parser():
    """Build and return the command-line option parser."""
    myparser = OptionParser(
        version=f"%prog version {VERSION}",
        formatter=NoWrapFormatter(),
    )
    myparser.add_option(
        "--tsv", action="store", type="string", dest="tsv_file_path",
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
    myparser.add_option(
        "--column", action="store", type="string",
        dest="column_with_frequencies", default='frequency',
        help="Name of a column in input TSV to be used for rendering"
             " frequencies [frequency]",
    )
    myparser.add_option(
        "--outfile-prefix", action="store", type="string",
        dest="outfile_prefix", default='',
        help="Output file prefix, eventually also with path. Output files"
             " will be PNG, JPG or PDF, HTML",
    )
    myparser.add_option(
        "--offset", action="store", type="int", dest="offset", default=0,
        help="Define offset value to be added to every amino acid position"
             " in the first column of the TSV input file at runtime. This is"
             " to obtain real amino acid position within the protein. Normally"
             " protein starts at position 1. Check the first aa in TSV file"
             " and provide whatever number to be added to it to get the"
             " desired amino acid position in the full-length protein.",
    )
    myparser.add_option(
        "--xmin", action="store", type="int", dest="xmin", default=0,
        help="Define minimum X-axis value. This should be the position in the padded alignment (using '-').",
    )
    myparser.add_option(
        "--xmax", action="store", type="int", dest="xmax", default=0,
        help="Define maximum X-axis value. This should be the position in the padded alignment (using '-').",
    )
    myparser.add_option(
        "--x-axis-bins", action="store", type="int", dest="xaxis_bins",
        default=0,
        help="Set number of bins (labels) on the X-axis. Could be used to"
             " override the amount of major ticks in a different way. Use 20"
             " for example. [default is 0]",
    )
    myparser.add_option(
        "--x-axis-major-ticks-spacing", action="store", type="int",
        dest="xaxis_major_ticks_spacing", default=10,
        help="Set distance between the major ticks on X-axis",
    )
    myparser.add_option(
        "--x-axis-minor-ticks-spacing", action="store", type="int",
        dest="xaxis_minor_ticks_spacing", default=5,
        help="Set distance between the minor ticks on X-axis",
    )
    myparser.add_option(
        "--x-axis-label-start", action="store", type="int",
        dest="xaxis_label_start", default=0,
        help="Set first label value on the X-axis",
    )
    myparser.add_option(
        "--aminoacids", action="store_true", dest="aminoacids", default=False,
        help="Draw chart with amino acid residues on Y-axis instead of"
             " codons. [default is False]",
    )
    myparser.add_option(
        "--show-STOP", action="store_true", dest="showstop", default=False,
        help="Include STOP codons or '*' in charts on Y-axis."
             " [default is False]",
    )
    myparser.add_option(
        "--show-INS", action="store_true", dest="showins", default=False,
        help="Include INS in charts on Y-axis. [default is False]",
    )
    myparser.add_option(
        "--show-DEL", action="store_true", dest="showdel", default=False,
        help="Include DEL in charts on Y-axis. [default is False]",
    )
    myparser.add_option(
        "--show-X", action="store_true", dest="showx", default=False,
        help="Include X in charts on Y-axis in --aminoacids mode."
             " [default is False]",
    )
    myparser.add_option(
        "--enable-colorbar", action="store_true", dest="colorbar",
        default=False,
        help="Enable colorbar [is Disabled by default]",
    )
    myparser.add_option(
        "--disable-short-legend", action="store_false", dest="shortlegend",
        default=True,
        help="Disable short legend in charts on X-axis."
             " [is Enabled by default]",
    )
    myparser.add_option(
        "--include-synonymous", action="store_true",
        dest="include_synonymous", default=False,
        help="Include synonymous changes in --aminoacids output as green"
             " diamonds. In codon output they are always shown."
             " [default is False]",
    )
    myparser.add_option(
        "--threshold", action="store", type="float", dest="threshold",
        default=0.001,
        help="Define minimum frequency threshold to display a pictogram in"
             " the output. For codon mode use 0.001 and for aa mode use 0.01."
             " [default: 0.001]",
    )
    myparser.add_option(
        "--title", action="store", type="string", dest="title", default='',
        help="Set title for the figures, by default trailing"
             " '.frequencies.tsv' is stripped from the end of the input"
             " filename",
    )
    myparser.add_option(
        "--disable-2nd-Y-axis", action="store_true",
        dest="disable_2nd_Y_axis", default=False,
        help="Disable rendering of the 2nd Y-axis showing sequencing coverage",
    )
    myparser.add_option(
        "--legend", action="store_true", dest="legend", default=False,
        help="Draw legend chart. [default is False]",
    )
    myparser.add_option(
        "--matrix", action="store", type="str", dest="matrix",
        default="BLOSUM80",
        help="BLOSUM matrix: BLOSUM45,BLOSUM50,BLOSUM62,BLOSUM80,BLOSUM90"
             " [default is BLOSUM80]",
    )
    myparser.add_option(
        "--matrix-file", action="store", type="string", dest="matrix_file",
        default=None,
        help="Matrix file compatible with BLOSUM matrices, e.g."
             " prot26050-sup-0004-Supinfo03.sm from"
             " https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8641535/bin/"
             "NIHMS1664401-supplement-supinfo.rar"
             " if you do not like default BLOSUM62",
    )
    myparser.add_option(
        "--colormap", action="store", type="string", dest="colormap",
        default='amino_acid_changes',
        help="Pick a colormap recognized by matplotlib."
             " See https://i.sstatic.net/cmk1J.png"
             " [default is coolwarm_r but seismic_r and coolwarm_r are"
             " great too]",
    )
    myparser.add_option(
        "--dpi", action="store", type="int", dest="dpi", default=600,
        help="DPI resolution for images",
    )
    myparser.add_option(
        "--backend", action="store", type="string", dest="backend", default='',
        help="Matplotlib backend to render resulting figures: agg, wxpython,"
             " pyqt5, pyqt6, pycairo, cairocffi [default: unset]\n"
             "To disable Matplotlib interactive window being raised up you"
             " can set MPLBACKEND=agg env variable.",
    )
    myparser.add_option(
        "--debug", action="store", type="int", dest="debug", default=0,
        help="Set debug to some real number",
    )
    myparser.add_option(
        "--disable-bokeh-sqrt-size", action="store_false",
        dest="bokeh_sqrt_size", default=True,
        help="Disable sqrt scaling for Bokeh circle sizes; size (diameter)"
             " will be proportional to frequency, area proportional to"
             " frequency². By default sqrt scaling is on, matching the"
             " perceptual appearance of the matplotlib figure.",
    )
    myparser.add_option(
        "--show-invisible-placeholder-dots", action="store_true",
        dest="show_invisible_placeholder_dots", default=False,
        help="Include below-threshold dots in the plot. [default: False]",
    )
    return myparser


def main():  # pylint: disable=too-many-locals
    """Entry point: parse options and run the full scatter-plot pipeline."""
    myparser = build_option_parser()
    myoptions, _ = myparser.parse_args()

    _matrix, _matrix_name, _min_theoretical_score, _max_theoretical_score, \
        _outfile_prefix = load_matrix(myoptions)

    # create a conversion dictionary
    _padded_position2position = {}

    # parse the .frequencies.tsv contents and fill-in the conversion dictionary
    _df, _padded_position2position = load_and_clean_dataframe(
        myoptions, myoptions.tsv_file_path, _padded_position2position,
    )

    print(f"Info: Writing into {_outfile_prefix}.actually_rendered.tsv")
    # Speedup 6: Export Decimal to string for reporting while keeping precision
    _df_to_save = _df.copy()
    _df_to_save[myoptions.column_with_frequencies] = _df_to_save[myoptions.column_with_frequencies].apply(lambda x: f"{x:.6f}")
    _df_to_save.to_csv(
        f"{_outfile_prefix}.actually_rendered.tsv",
        sep='\t', header=None, index=False,
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
    _df_frequencies_unchanged_codons, _padded_position2position = load_and_clean_dataframe(
        myoptions, _unchanged_tsv, _padded_position2position,
    )
    del _df_frequencies_unchanged_codons

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
            _new_aa_table, _new_codon_table,
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

    if _circles_bokeh:
        render_bokeh(
            myoptions,
            _outfile_prefix, _xmin, _xmax, _amino_acids,
            _final_sorted_whitelist,
            _circles_bokeh, _mutations, _hover_text_bokeh,
            _title_data, _xlabel,
            _matrix_name, _colors, _norm, _cmap,
        )

    render_matplotlib(
        myoptions,
        _figure, _ax1, _ax2, _ax3, _ax4, _outfile_prefix,
        _circles_matplotlib, _markers, _dots, _cmap, _norm, _colors,
        _matrix, _matrix_name,
        _new_aa_table, _new_codon_table, _df, _codons_whitelist2,
        _final_sorted_whitelist,
        _padded_position2position,
    )


if __name__ == "__main__":
    main()

# vim:ts=4:sw=4:expandtab:smartindent
