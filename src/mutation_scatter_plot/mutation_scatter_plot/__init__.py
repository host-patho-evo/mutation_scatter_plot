# This work © 2025 by Jiří Zahradník and Martin Mokrejš
# (First Medical Faculty - Charles University in Prague) is licensed under
# Creative Commons Attribution 4.0 International. To view a copy of this
# license, visit https://creativecommons.org/licenses/by/4.0/

"""This program can be used to render figures from input TSV files using
Matplotlib and Bokeh. The figures show all amino acid residues or their
underlying codons in a protein-coding region per each position. It counts
naturally from 1, not from zero like python does internally.

Matplotlib provides PNG, JPG, PDF and interactive outputs whereas Bokeh
provides JSON and HTML+javascript outputs providing similar interactivity
upon mouse hover().

If the running python instance and its Matplotlib has access to an X11 or
Wayland screen or other renderer it runs an interactive window with the
figure and one can use computer mouse to hover() above the circles in the
scatter plot. It shows a bubble-like pop-up legend of the particular data
item but also, the script prints out on STDOUT of the command line terminal
the underlying Pandas dataframe lines for the amino acid residue and its
codon.

It can be run either in codon mode or in amino acid mode, default display
thresholds differ for these.

To respect reading frame calculate_codon_frequencies.py follows reference
sequence and parses always 3 nucleotides at once from the ALN input file.
If there was a gap (padding) character '-' if picks an extra character from
the input until it has 3 nucleotides. Then it calculates frequencies for
all possible codons. The results are stored in a TAB-separated TSV file for
easy post-processing by mutation_scatter_plot.py which on-the-fly discards
codons containing unknown (N) nucleotides using python Pandas library and
finally draws interactive figures using Matplotlib and Bokeh graphical
libraries.

To color the figures mutation_scatter_plot.py uses the widely used BLOSUM62
matrix from
https://www.ncbi.nlm.nih.gov/IEB/ToolBox/C_DOC/lxr/source/data/BLOSUM62 to
discern evolutionarily conservative amino acid changes from those more
drastic. For each amino acid replacement pair there is a value in the matrix.
Notably, the values are not centered around zero. Here is an overview:

BLOSUM45:     [-5,     -3, -2, -1, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10,    12]     Range spans 18 values but we need 25 for centering.
BLOSUM50:     [-5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5, 6, 7, 8,    10,        13] Range spans 19 values but we need 27 for centering.
BLOSUM62:         [-4, -3, -2, -1, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9]                Range spans 14 values but we need 19 for centering.
BLOSUM80: [-6,     -4, -3, -2, -1, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9]                Range spans 16 values but we need 19 for centering.
BLOSUM90: [-6, -5, -4, -3, -2, -1, 0, 1, 2, 3,    5, 6, 7, 8, 9]                Range spans 16 values but we need 19 for centering.


Colormaps

We prefer either discerning colormaps of two colors separated by some narrow
band of a color in the middle  (for example represented by gray color in the
"coolwarm_r" color palette).
The color palette "coolwarm_r" visible at https://i.sstatic.net/cmk1J.png is
used for drawing ranges from red to blue (in dark red should be rather
pronounced change while in dark blue should be functionally similar amino
acid change).
Users can discriminate about 3 colors only. That is not much but the figures
are easy to interpret.

Please cite our article if you use our data or software in your work:

Shoshany A., Tian R., Padilla-Blanco M., Hruška A., Baxova K., Zoler E., Mokrejš M., Schreiber G., Zahradník J. (submitted) In Vitro and Viral Evolution Convergence Reveal the Selective Pressures Driving Omicron Emergence. [bioRxiv](https://www.biorxiv.org/content/10.1101/2025.04.23.650148v1)
"""

import os
import sys
import re
import json
import typing

from decimal import Decimal, ExtendedContext, setcontext, getcontext
import pandas as pd

import matplotlib
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

import numpy as np
import mplcursors

# https://pypi.org/project/blosum/
import blosum

# https://docs.bokeh.org/en/latest/docs/user_guide/interaction/tools.html#ug-interaction-tools-hover-tool
import bokeh.plotting
import bokeh.models

from .. import alt_translate

setcontext(ExtendedContext)
c = getcontext()
c.prec = 99

VERSION = "0.3"



__all__ = [
    "VERSION",
    "alt_translate",
    "get_colormap",
    "resolve_codon_or_aa",
    "get_score",
    "adjust_size_and_color",
    "adjust_size_and_color_neutralized_escape",
    "adjust_size_and_color_weighted",
    "load_matrix",
    "load_and_clean_dataframe",
    "build_frequency_tables",
    "build_conversion_table",
    "setup_matplotlib_figure",
    "collect_scatter_data",
    "render_bokeh",
    "render_matplotlib",
]


def get_colormap(myoptions, colormapname):
    _norm = None
    _colors = None
    try:
        _cmap = matplotlib.colormaps.get_cmap(colormapname)
    except Exception:
        try:
            _cmap = plt.get_cmap(colormapname)
        except Exception:
            try:
                import nclcmaps
                _cmap = nclcmaps.cm.get_cmap(colormapname)
            except (ImportError, Exception):
                try:
                    import palettable
                    _wished_cmapname_prefix, _wished_cmapname_num = '_'.join(myoptions.colormap.split('_')[:-1]), int(myoptions.colormap.split('_')[-1])
                    print(f"Debug: _wished_cmapname_prefix={_wished_cmapname_prefix}, _wished_cmapname_num={_wished_cmapname_num}")
                    _cmap_from_palettable = palettable.colorbrewer.get_map(_wished_cmapname_prefix, 'diverging', _wished_cmapname_num)
                    _cmap = matplotlib.colors.ListedColormap(_cmap_from_palettable.mpl_colors)
                except (ImportError, Exception):
                    try:
                        import palettable
                        _cmap_from_palettable = palettable.colorbrewer.get_map(_wished_cmapname_prefix, 'sequential', _wished_cmapname_num)
                        _cmap = matplotlib.colors.ListedColormap(_cmap_from_palettable.mpl_colors)
                    except (ImportError, Exception):
                        try:
                            import palettable
                            _cmap_from_palettable = palettable.scientific.get_map(_wished_cmapname_prefix, 'diverging', _wished_cmapname_num)
                            _cmap = matplotlib.colors.ListedColormap(_cmap_from_palettable.mpl_colors)
                        except (ImportError, Exception):
                            try:
                                import palettable
                                _cmap_from_palettable = palettable.scientific.get_map(_wished_cmapname_prefix, 'sequential', _wished_cmapname_num)
                                _cmap = matplotlib.colors.ListedColormap(_cmap_from_palettable.mpl_colors)
                            except (ImportError, Exception):
                                try:
                                    import palettable
                                    _cmap_from_palettable = palettable.tableau.get_map(myoptions.colormap)
                                    _cmap = matplotlib.colors.ListedColormap(_cmap_from_palettable.mpl_colors)
                                except (ImportError, Exception):
                                    _micro_cvd_gray = ["#F5F5F5", "#D6D6D6", "#B7B7B7", "#8B8B8B","#616161"]
                                    _micro_cvd_green = ["#DDFFA0",  "#BDEC6F",  "#97CE2F", "#6D9F06","#4E7705"]
                                    _micro_cvd_orange = ["#FFD5AF",  "#FCB076","#F09163", "#C17754", "#9D654C"]
                                    _micro_cvd_blue = ["#E7F4FF", "#BCE1FF", "#7DCCFF", "#56B4E9","#098BD9"]
                                    _micro_cvd_turquoise = ["#A3E4D7", "#48C9B0",  "#43BA8F",  "#009E73", "#148F77"]
                                    _micro_cvd_purple = ["#EFB6D6", "#E794C1", "#CC79A7", "#A1527F", "#7D3560"]
                                    _micro_cvd_gray_r = _micro_cvd_gray[::-1]
                                    _micro_cvd_green_r = _micro_cvd_green[::-1]
                                    _micro_cvd_orange_r = _micro_cvd_orange[::-1]
                                    _micro_cvd_blue_r = _micro_cvd_blue[::-1]
                                    _micro_cvd_turquoise_r = _micro_cvd_turquoise[::-1]
                                    _micro_cvd_purple_r = _micro_cvd_purple[::-1]

                                    if myoptions.colormap == 'microshades_cvd_palettes':
                                        _colors = list(_micro_cvd_orange) + list(_micro_cvd_turquoise) + list(_micro_cvd_blue) + list(_micro_cvd_purple) + list(_micro_cvd_green) + list(_micro_cvd_gray)
                                        _cmap = matplotlib.colors.ListedColormap(_colors, 'microshades_cvd_palettes')
                                    elif myoptions.colormap == 'microshades_cvd_palettes_r':
                                        _colors = list(_micro_cvd_orange_r) + list(_micro_cvd_turquoise_r) + list(_micro_cvd_blue_r) + list(_micro_cvd_purple_r) + list(_micro_cvd_green_r) + list(_micro_cvd_gray_r)
                                        _cmap = matplotlib.colors.ListedColormap(_colors, 'microshades_cvd_palettes_r')
                                    elif myoptions.colormap == 'merged':
                                        _colors1 = plt.cm.gist_heat(np.linspace(0.2, 0.8, 9))
                                        _colors2 = plt.cm.cividis_r(np.linspace(0.2, 0.8, 9))
                                        _colors3 = plt.cm.GnBu_r(np.linspace(0.2, 0.8, 9))
                                        _colors = np.vstack((_colors1, _colors2, _colors3))
                                        _cmap = matplotlib.colors.LinearSegmentedColormap.from_list('merged', _colors)
                                    elif myoptions.colormap == 'merged1':
                                        _colors1 = plt.cm.viridis(np.linspace(0.2, 0.8, 9))
                                        _colors2 = plt.cm.coolwarm_r(np.linspace(0., 1, 9))
                                        _colors3 = plt.cm.GnBu_r(np.linspace(0.2, 0.8, 9))
                                        _colors = np.vstack((_colors1, _colors2, _colors3))
                                        _cmap = matplotlib.colors.LinearSegmentedColormap.from_list('merged', _colors)
                                    elif myoptions.colormap == 'dkeenan_26cols':
                                        _colors = ["#00B7FF", "#004DFF", "#00FFFF", "#826400", "#580041", "#FF00FF", "#00FF00", "#C500FF", "#B4FFD7", "#FFCA00", "#969600", "#B4A2FF", "#C20078", "#000000", "#0000C1", "#FF8B00", "#FFC8FF", "#666666", "#FF0000", "#CCCCCC", "#009E8F", "#D7A870", "#8200FF", "#960000", "#BBFF00", "#FFFF00", "#006F00"]
                                        _cmap = matplotlib.colors.ListedColormap(_colors, "dkeenan_26cols")
                                        myoptions.colormap = 'dkeenan_26cols'
                                        _norm = matplotlib.colors.BoundaryNorm(np.arange(-13, 13, 1), _cmap.N)
                                    else:
                                        _colors = ["#930000", "#930000", "#930000", "#930000", "#930000", "#930000", "#960000", "#580041", "#8200ff", "#c500ff", "#ff00fd", "#CC79A7", "#eea1d0", "#cc0000", "#ff0000", "#ff4f00", "#ff7c7c", "#ff9999", "#c58a24", "#9c644b", "#ffff00", "#ffcc00", "#ffa200", "#7DCCFF", "#0042ff", "#0000ff", "#D6D6D6", "#B7B7B7", "#8B8B8B", "palegreen", "#bbff00", "#97CE2F", "#219f11", "#930000", "#930000", "#930000", "#930000", "#930000", "#930000"]
                                        _cmap = matplotlib.colors.ListedColormap(_colors, "amino_acid_changes", len(_colors))
                                        myoptions.colormap = 'amino_acid_changes'
                                        _norm = matplotlib.colors.BoundaryNorm(np.arange(-19, 19, 1), _cmap.N)
                                        return _norm, _cmap, _colors

    return _norm, _cmap, _colors


def resolve_codon_or_aa(myoptions, old_codon_or_aa, new_codon_or_aa):
    _len_new_codon_or_aa = len(new_codon_or_aa)

    if myoptions.aminoacids:
        _old_codon_or_aa = alt_translate(old_codon_or_aa)
        if new_codon_or_aa == 'NNN':
            _new_codon_or_aa = 'X'
            _len_new_codon_or_aa = len(_new_codon_or_aa)
        elif _len_new_codon_or_aa > 1 and new_codon_or_aa not in ('---', 'DEL', 'INS'):
            print(f"Info: Weird, the new_codon_or_aa={new_codon_or_aa}")
            _new_codon_or_aa = new_codon_or_aa.upper()
        else:
            _new_codon_or_aa = new_codon_or_aa.upper()
    else:
        _old_codon_or_aa = old_codon_or_aa.upper()
        _new_codon_or_aa = new_codon_or_aa.upper()

    if len(_old_codon_or_aa) > 1 and old_codon_or_aa not in ('---', 'DEL', 'INS'):
        _codon_on_input = True
    elif len(_new_codon_or_aa) > 1 and new_codon_or_aa not in ('---', 'DEL', 'INS'):
        _codon_on_input = True
    else:
        _codon_on_input = False
    return _codon_on_input, _old_codon_or_aa, _new_codon_or_aa


# Module-level cache for substitution scores — cleared at the start of each plot run.
# Avoids repeated BLOSUM matrix lookups for the same (old_aa, new_aa) pairs.
_score_cache: dict = {}


def get_score(myoptions, matrix, codon_on_input, old_codon_or_aa, new_codon_or_aa):
    if codon_on_input:
        _old_aa = alt_translate(old_codon_or_aa)
        _new_aa = alt_translate(new_codon_or_aa)
    else:
        _old_aa = old_codon_or_aa
        _new_aa = new_codon_or_aa
    _cache_key = (_old_aa, _new_aa)
    if _cache_key in _score_cache:
        return _score_cache[_cache_key]
    try:
        _score = int(matrix[_old_aa][_new_aa])
    except OverflowError:
        _score = -11
    except KeyError as exc:
        raise ValueError(f"Cannot get a score for myoptions.matrix='{myoptions.matrix}', old_codon_or_aa='{old_codon_or_aa}', new_codon_or_aa='{new_codon_or_aa}'") from exc
    _score_cache[_cache_key] = _score
    return _score


def adjust_size_and_color(myoptions, frequency, codon_on_input, old_codon_or_aa, new_codon_or_aa, _old_codon_or_aa, _new_codon_or_aa, matrix, norm, colors):
    """Define colors for basic type of scatterplot figures."""

    if not codon_on_input:
        if not _old_codon_or_aa:
            raise ValueError(f"Aieee, old_codon_or_aa='{str(old_codon_or_aa)}' is empty")
        elif not _new_codon_or_aa:
            raise ValueError(f"Aieee, new_codon_or_aa='{str(new_codon_or_aa)}' is empty")
        else:
            if myoptions.debug:
                print(f"Info: some single-letter but neither asterisk nor dash amino acid residue: _old_codon_or_aa={str(_old_codon_or_aa)}, new_codon_or_aa='{str(new_codon_or_aa)}'")

    if _new_codon_or_aa in ('---', 'DEL', 'INS'):
        _score = -11
    else:
        _score = get_score(myoptions, matrix, codon_on_input, _old_codon_or_aa, _new_codon_or_aa)

    _colorindex = norm(_score)

    if old_codon_or_aa.upper() == new_codon_or_aa.upper():
        _color = '#00ff04'
        _score = 12
    elif _old_codon_or_aa.upper() == _new_codon_or_aa.upper():
        _color = '#00ff04'
        _score = 12
    elif old_codon_or_aa in ('---', 'DEL', 'INS', '*') or new_codon_or_aa in ('---', 'DEL', 'INS', '*', 'TGA', 'TAA', 'TAG'):
        _color = 'red'
        _score = -6
    elif _old_codon_or_aa in ('X', 'NNN') or new_codon_or_aa in ('X', 'NNN'):
        _color = 'gray'
    elif codon_on_input:
        if alt_translate(_old_codon_or_aa) == alt_translate(_new_codon_or_aa):
            _color = '#00ff04'
            _score = 12
        elif alt_translate(_new_codon_or_aa) == 'X':
            _color = 'gray'
        else:
            if myoptions.debug:
                sys.stdout.write(f"Info: Translating {_old_codon_or_aa} to {_new_codon_or_aa} to fetch color from _colors,")
            _color = colors[_colorindex]
            if myoptions.debug:
                sys.stdout.write(f" which has yielded {_colorindex} and {str(_color)}{os.linesep}")
    else:
        if myoptions.debug:
            sys.stdout.write(f"Info: Fetching color for {_old_codon_or_aa} to {_new_codon_or_aa} directly from _colors,")
        _color = colors[_colorindex]
        if myoptions.debug:
            sys.stdout.write(f" which has yielded {_colorindex} and {str(_color)}{os.linesep}")

    return _score, frequency, _color


def adjust_size_and_color_neutralized_escape(neutralized_parent_difference, codon_on_input):
    "Used only for neutralized_parent_difference and escape_parent_difference figure types."

    if neutralized_parent_difference < -0.001:
        _size = abs(neutralized_parent_difference)
        _color = 'red'
    elif neutralized_parent_difference > 0.001:
        if not codon_on_input:
            _size, _color = 0, '#00ff04'
        else:
            _size = neutralized_parent_difference
            _color = '#00ff04'
    else:
        _size = 0
        _color = 'black'
    return _size, _color


def adjust_size_and_color_weighted(weighted_diff_escape_neutralized, generic_circle_size=5000, weighted_diff_escape_neutralized_size1=3000, weighted_diff_escape_neutralized_size2=2000):
    "Used only for weighted_diff_escape_neutralized figure type."

    if weighted_diff_escape_neutralized > 1:
        _size = (weighted_diff_escape_neutralized) * generic_circle_size
        _color = 'black'
    elif weighted_diff_escape_neutralized > 0.1:
        _size = (weighted_diff_escape_neutralized) * weighted_diff_escape_neutralized_size1
        _color = 'darkblue'
    elif weighted_diff_escape_neutralized > 0.01:
        _size = (weighted_diff_escape_neutralized) * weighted_diff_escape_neutralized_size2
        _color = 'skyblue'
    else:
        _size = abs(weighted_diff_escape_neutralized) * 0
        _color = 'black'
    return _size, _color


def load_matrix(myoptions):
    """Load BLOSUM or custom substitution matrix and set up the output file prefix."""
    if myoptions.matrix_file and os.path.exists(myoptions.matrix_file):
        print(myoptions.matrix)
        if not myoptions.matrix:
            myoptions.matrix = re.sub(r'.*/', '', myoptions.matrix_file)
        print(myoptions.matrix)
        _matrix_type, _matrix_num = re.sub(r'\d+', '', myoptions.matrix), int(re.sub(r'[a-zA-Z]+', '', myoptions.matrix))
        _matrix = blosum.BLOSUM(myoptions.matrix_file)
        _matrix_name = myoptions.matrix_file.split(os.path.sep)[-1]
        myoptions.matrix = _matrix_name
    else:
        _matrix_type, _matrix_num = re.sub(r'\d+', '', myoptions.matrix), int(re.sub(r'[a-zA-Z]+', '', myoptions.matrix))
        if _matrix_type == 'BLOSUM':
            _matrix = blosum.BLOSUM(_matrix_num)
            _matrix_name = "BLOSUM%d" % (_matrix_num)
        else:
            sys.stderr.write(f"Warning: Unexpected matrix type {str(_matrix_type)}, falling back to BLOSUM\n")
            _matrix = blosum.BLOSUM(_matrix_num)
            _matrix_name = "BLOSUM%d" % (_matrix_num)
            myoptions.matrix = _matrix_name

    if not myoptions.outfile_prefix:
        raise RuntimeError("Please provide output filename prefix via --outfile-prefix")
    else:
        _outfile_prefix = myoptions.outfile_prefix + '.' + _matrix_name + '.' + myoptions.colormap
        print(f"Info: _outfile_prefix={_outfile_prefix}")

    _theoretical_scores = set()
    for _aa in _matrix.keys():
        for _score in _matrix[_aa].values():
            _theoretical_scores.add(_score)
    print(f"Info: Using {_matrix_name} matrix now. Theoretical minimum score is {min(_theoretical_scores)}, theoretical maximum score is {max(_theoretical_scores)}, values are {str(_theoretical_scores)}")
    _min_theoretical_score, _max_theoretical_score = int(min(_theoretical_scores)), int(max(_theoretical_scores))

    return _matrix, _matrix_name, _min_theoretical_score, _max_theoretical_score, _outfile_prefix


def load_and_clean_dataframe(myoptions, infilename, outfile_prefix, padded_position2position):
    """Load the input TSV, normalise column names for legacy formats, and filter noisy rows.
    Parse only rows with codons [ATGCatgc-], so not those with more exotic IUPAC codes.
    Provided we discard some data on-the-fly, we have to create the mapping dictionary
    ideally before that.

    By design the real _aa_position or _codon_position will always be reused for adjacent (rightwards)
    codons. Only the number in the _padded_position can be trusted to lead back to the actual
    codon and its frequency:

    243     243     T       T       0.526316        ACC     ACG     20      38
    243     243     T       V       0.105263        ACC     GTC     4       38
    243     243     T       DEL     0.342105        ACC     ---     13      38
    244     244     Q       H       0.052632        CAG     CAC     2       38
    244     244     Q       DEL     0.342105        CAG     ---     13      38
    245     244     INS     V       0.105263        ---     GTC     4       38
    246     245     R       DEL     0.342105        AGG     ---     13      38

    Attempts to get back to CAG -> CAC change at real codon position 244 will return more than one result.
    It would return much much results if there were more INSertions after each other.
    """

    print(f"Info: Parsing input file {myoptions.tsv_file_path}")
    if not myoptions.tsv_file_path:
        raise RuntimeError("Please provide an input TSV file via --tsv")
    if not os.path.exists(myoptions.tsv_file_path):
        raise RuntimeError(f"Input file not found: {myoptions.tsv_file_path}")
    if os.path.getsize(myoptions.tsv_file_path) == 0:
        raise RuntimeError(f"Input file is empty: {myoptions.tsv_file_path}")

    df = pd.read_csv(myoptions.tsv_file_path, sep='\t', header='infer', na_filter=False, na_values=[None])

    if 'position' in df.columns:
        print(f"Info: Autodetected new TSV file format with a header in {myoptions.tsv_file_path}")
    else:
        print(f"Info: Autodetected old TSV file format without a header in {myoptions.tsv_file_path}, assigning default column names")
        df = pd.read_csv(myoptions.tsv_file_path, sep='\t', header=None, na_filter=False, na_values=[None])
        _count_columns = len(df.columns.values)
        if _count_columns == 9:
            df.columns = ['padded_position', 'position', 'original_aa', 'mutant_aa', 'frequency', 'original_codon', 'mutant_codon', 'observed_codon_count', 'total_codons_per_site']
        elif _count_columns == 11:
            df.columns = ['padded_position', 'position', 'original_aa', 'mutant_aa', 'frequency', 'original_codon', 'mutant_codon', 'observed_codon_count', 'total_codons_per_site', 'frequency_parent', 'frequency_selected']
        elif _count_columns == 10:
            df.columns = ['position', 'original_aa', 'mutant_aa', 'frequency', 'original_codon', 'mutant_codon', 'observed_codon_count', 'total_codons_per_site', 'frequency_parent', 'frequency_selected']
        elif _count_columns == 8:
            df.columns = ['position', 'original_aa', 'mutant_aa', 'frequency', 'original_codon', 'mutant_codon', 'observed_codon_count', 'total_codons_per_site']
        elif _count_columns == 6:
            df.columns = ['position', 'original_aa', 'mutant_aa', 'frequency', 'original_codon', 'mutant_codon']
        else:
            raise RuntimeError(f"Unexpected number of columns in the {myoptions.tsv_file_path} file")

    print(f"Info: The file {myoptions.tsv_file_path} contains now these columns: {str(df.columns)}")

    if 'padded_position' not in df.columns:
        df['padded_position'] = df['position']

    if myoptions.offset:
        df['position'] = df['position'] + int(myoptions.offset)
        df['padded_position'] = df['padded_position'] + int(myoptions.offset)

    build_conversion_table(df, padded_position2position) # parse the data

    _mutant_codon = df['mutant_codon']
    _before = len(_mutant_codon)
    try:
        df = df.loc[_mutant_codon.str.match('[ATGCatgc-][ATGCatgc-][ATGCatgc-]')]
    except Exception as _exc:
        raise ValueError("Cannot parse column df['mutant_codon']=%s containing %d values" % (str(_mutant_codon), len(_mutant_codon))) from _exc

    _mutant_aa = df['mutant_aa']
    _aas_to_filter = ['X']
    if not myoptions.showstop:
        _aas_to_filter.append('*')
    if not myoptions.showdel:
        _aas_to_filter.append('DEL')
    if not myoptions.showins:
        _aas_to_filter.append('INS')
    df = df.loc[~_mutant_aa.isin(_aas_to_filter)]
    _after = len(df['mutant_codon'])
    if myoptions.showstop:
        if myoptions.showdel:
            print("Info: Originally there were %d rows but after discarding codons with [N n] there are only %d left" % (_before, _after))
        else:
            print("Info: Originally there were %d rows but after discarding codons with [N n DEL] there are only %d left" % (_before, _after))
    elif myoptions.showdel:
        print("Info: Originally there were %d rows but after discarding codons with [N n *] there are only %d left" % (_before, _after))
    else:
        print("Info: Originally there were %d rows but after discarding codons with [N n DEL] there are only %d left" % (_before, _after))
    return df, padded_position2position


def build_frequency_tables(myoptions, df, padded_position2position):
    """Build amino-acid and codon frequency tables from the cleaned dataframe.
    It is at the beginning initialized with Decimal(0).

    Use 'padded_position' to label the columns instead of 'position', so that thing
    do not break if we trip over INSertions. Due to that transition we need
    padded_position2position to convert the positions.
    """

    _amino_acids = ['C', 'R', 'K', 'E', 'Q', 'D', 'N', 'T', 'S', 'H', 'M', 'P', 'W', 'Y', 'F', 'V', 'L', 'I', 'A', 'G']
    if myoptions.showins:
        _amino_acids.append('INS')
    if myoptions.showx:
        _amino_acids.append('X')
    if myoptions.showstop:
        _amino_acids.append('*')
    if myoptions.showdel:
        _amino_acids.append('DEL')

    _unique_padded_aa_positions = sorted(padded_position2position.keys())
    if myoptions.debug:
        print(f"Debug: len(_unique_padded_aa_positions)={len(_unique_padded_aa_positions)}, _unique_padded_aa_positions: {str(_unique_padded_aa_positions)}")
    _unique_padded_codon_positions = list(set(_unique_padded_aa_positions))
    if myoptions.debug:
        print(f"Debug: len(_unique_padded_codon_positions)={len(_unique_padded_codon_positions)}, _unique_padded_codon_positions: {str(_unique_padded_codon_positions)}")

    _codons_whitelist = ['TTT', 'TTC', 'TTA', 'TTG', 'CTT', 'CTC', 'CTA', 'CTG', 'ATT', 'ATC', 'ATA', 'ATG', 'GTT', 'GTC', 'GTA', 'GTG', 'TCT', 'TCC', 'TCA', 'TCG', 'CCT', 'CCC', 'CCA', 'CCG', 'ACT', 'ACC', 'ACA', 'ACG', 'GCT', 'GCC', 'GCA', 'GCG', 'TAT', 'TAC', 'TAA', 'TAG', 'CAT', 'CAC', 'CAA', 'CAG', 'AAT', 'AAC', 'AAA', 'AAG', 'GAT', 'GAC', 'GAA', 'GAG', 'TGT', 'TGC', 'TGA', 'TGG', 'CGT', 'CGC', 'CGA', 'CGG', 'AGT', 'AGC', 'AGA', 'AGG', 'GGT', 'GGC', 'GGA', 'GGG']
    if not myoptions.showstop:
        for _stopcodon in ('TGA', 'TAA', 'TAG'):
            _codons_whitelist.remove(_stopcodon)

    if myoptions.debug:
        print("Debug: Translating pre-made all possible codons")
    _codons_whitelist_aa = [alt_translate(_codon) for _codon in _codons_whitelist]
    if 'DEL' in _amino_acids:
        _codons_whitelist.append('---')
        _codons_whitelist_aa.append('DEL')

    if len(_codons_whitelist) != len(_codons_whitelist_aa):
        raise ValueError("Length of _codons_whitelist is %d which is not equal to _codons_whitelist_aa with %d" % (len(_codons_whitelist), len(_codons_whitelist_aa)))
    _sorted_whitelist = sorted(zip(_codons_whitelist, _codons_whitelist_aa), key=lambda x: x[1])

    if myoptions.debug:
        print(f"Debug: _sorted_whitelist={str(_sorted_whitelist)}")

    _final_sorted_whitelist = [_tuple for x in _amino_acids for _tuple in _sorted_whitelist if _tuple[1] == x]
    if myoptions.debug:
        print(f"Debug: _final_sorted_whitelist={str(_final_sorted_whitelist)}")

    _codons_whitelist2 = [x[0] for x in _final_sorted_whitelist]

    if myoptions.debug:
        print(f"Debug: _codons_whitelist2={str(_codons_whitelist2)}")

    _old_aa_table = pd.DataFrame(Decimal(0), index=_amino_acids, columns=_unique_padded_aa_positions)
    _new_aa_table = pd.DataFrame(Decimal(0), index=_amino_acids, columns=_unique_padded_aa_positions) # the table contains 64+1 rows and len(_unique_padded_codon_positions) columns
    _old_codon_table = pd.DataFrame(Decimal(0), index=_codons_whitelist2, columns=_unique_padded_codon_positions)
    _new_codon_table = pd.DataFrame(Decimal(0), index=_codons_whitelist2, columns=_unique_padded_codon_positions) # the table contains 64+1 rows and len(_unique_padded_codon_positions) columns

    _very_leftmost_aa_pos = None
    _calculated_aa_offset = 0
    # make tables with yet another number of rows summing up eventually the frequencies
    for _df_index, _row in df.iterrows():
        _padded_position = _row['padded_position']
        # It is not necessary to skip N-containing codons as we anyway draw just those in codons list. Skipping some rows would make new_aa_table and new_codon_table have a different amount of rows, breaking slicing
        if _very_leftmost_aa_pos is None:
            _very_leftmost_aa_pos = int(_padded_position)
            # # if AA positions in the input file do NOT start from the first-one the numbering of sites gets shifted, so calculate the offset
            _calculated_aa_offset = _very_leftmost_aa_pos - myoptions.offset + 1
            if myoptions.debug:
                print(f"Debug: calculated offset is {_calculated_aa_offset}")
        _old_amino_acid = _row['original_aa']
        _new_amino_acid = _row['mutant_aa']
        _frequency = Decimal(_row[myoptions.column_with_frequencies])
        _old_codon = _row['original_codon'].upper()
        _new_codon = _row['mutant_codon'].upper()

        if (not myoptions.aminoacids or _old_amino_acid != _new_amino_acid or myoptions.include_synonymous) and not np.abs(_frequency) < myoptions.threshold:
            _old_value = Decimal(0)
            try:
                _old_value = _new_aa_table.at[_new_amino_acid, _padded_position]
            except KeyError:
                _new_aa_table.at[_new_amino_acid, _padded_position] = Decimal(_frequency)
            except TypeError as exc:
                raise TypeError(f"Weird value {_new_aa_table.at[_new_amino_acid, _padded_position]}") from exc
            else:
                _new_aa_table.at[_new_amino_acid, _padded_position] = Decimal(_old_value) + Decimal(_frequency)

            _old_value = Decimal(0)
            try:
                _old_value = _old_aa_table.at[_old_amino_acid, _padded_position]
            except KeyError:
                _old_aa_table.at[_old_amino_acid, _padded_position] = Decimal(_frequency)
            except TypeError as exc:
                raise TypeError(f"Weird value {_old_aa_table.at[_old_amino_acid, _padded_position]}") from exc
            else:
                _old_aa_table.at[_old_amino_acid, _padded_position] = Decimal(_old_value) + Decimal(_frequency)

            _old_value = Decimal(0)
            try:
                _old_value = Decimal(_old_codon_table.at[_old_codon, _padded_position])
            except KeyError:
                _old_codon_table.at[_old_codon, _padded_position] = Decimal(_frequency)
            except TypeError as exc:
                raise TypeError(f"Weird value {_old_codon_table.at[_old_codon, _padded_position]}") from exc
            else:
                _old_codon_table.at[_old_codon, _padded_position] = Decimal(_old_codon_table.at[_old_codon, _padded_position]) + Decimal(_frequency)

            _old_value = Decimal(0)
            _new_codon_table.at[_new_codon, _padded_position] = Decimal(_frequency)
            if myoptions.debug:
                print(f"Debug: OriginalDataFrameRowNumber: {_df_index}, Old: {_old_amino_acid}, New: {_new_amino_acid}, Frequency: {_frequency}")

    if myoptions.debug:
        for t in (_old_aa_table, _new_aa_table, _old_codon_table, _new_codon_table):
            print(f"Debug: len({t})={len(t)}")

    if myoptions.debug:
        if myoptions.aminoacids:
            print(_old_aa_table)
            print(_new_aa_table)
        else:
            print(_old_codon_table)
            print(_new_codon_table)

    return (
        _amino_acids, _codons_whitelist, _codons_whitelist2, _final_sorted_whitelist,
        _unique_padded_aa_positions, _unique_padded_codon_positions,
        _old_aa_table, _new_aa_table, _old_codon_table, _new_codon_table,
        _calculated_aa_offset, padded_position2position
    )


def build_conversion_table(df, padded_position2position):
    """Create a conversion table from padded_position to position from the data parsed from the input TSV file.
    The integers need not be contiguous, because the outfile_prefix.frequencies.tsv are filtered listings of
    mutations above some frequency threshold. One can supplement them with values from
    outfile_prefix.frequencies.unchanged_codons.tsv file.
    """

    # df.columns = ['padded_position', 'position', 'original_aa', 'mutant_aa', 'frequency', 'original_codon', 'mutant_codon', 'observed_codon_count', 'total_codons_per_site', 'frequency_parent', 'frequency_selected']
    for _df_index, _row in df.iterrows():
        _padded_position = int(_row['padded_position'])
        _position = int(_row['position'])
        if _padded_position not in padded_position2position:
            padded_position2position[_padded_position] = _position

    return padded_position2position


def setup_matplotlib_figure(
    myoptions,
    title_data, aln_rows, matrix_name, amino_acids, codons_whitelist,
    final_sorted_whitelist, unique_aa_padded_positions, unique_padded_codon_positions,
    new_aa_table, new_codon_table, padded_position2position,
):
    """Configure matplotlib figure, axes, labels, ticks, and the frequency bar chart."""
    matplotlib.rcParams['font.family'] = 'monospace'
    plt.rcParams["pdf.use14corefonts"] = True
    plt.rcParams["ps.useafm"] = True

    print(f"Info: matplotlib.get_backend={matplotlib.get_backend()}")

    _figure, (_ax1, _ax3, _ax4) = plt.subplots(1, 3, figsize=(16, 9), width_ratios=[55, 1, 6])

    if myoptions.aminoacids:
        if myoptions.shortlegend:
            _xlabel = 'Padded amino acid position'
        else:
            _xlabel = f'Padded amino acid position{os.linesep}based on {aln_rows.strip(os.linesep)} ALN rows, matrix {matrix_name}, colormap {myoptions.colormap}, mutation_scatter_plot {VERSION}'
    else:
        if myoptions.shortlegend:
            _xlabel = 'Padded codon position'
        else:
            _xlabel = f'Padded codon position{os.linesep}based on {aln_rows.strip(os.linesep)} ALN rows, matrix {matrix_name}, colormap {myoptions.colormap}, mutation_scatter_plot {VERSION}'
    _ax1.set_xlabel(_xlabel, fontsize=14)
    if myoptions.aminoacids:
        _ax1.set_ylabel('Introduced amino acid changes', fontsize=14)
        _ax1.set_title(title_data, fontsize=18)
        if myoptions.xmin:
            _xmin = myoptions.xmin
        else:
            _xmin = min(unique_aa_padded_positions) - 1

        if myoptions.xmax:
            _xmax = myoptions.xmax
        else:
            _xmax = max(unique_aa_padded_positions) + 1 # this should be the position in the padded alignment
    else:
        _ax1.set_ylabel('Introduced codon changes', fontsize=14)
        _ax1.set_title(title_data, fontsize=14)
        if myoptions.xmin:
            _xmin = myoptions.xmin
        else:
            _xmin = min(unique_padded_codon_positions) - 1

        if myoptions.xmax:
            _xmax = myoptions.xmax
        else:
            _xmax = max(unique_padded_codon_positions) + 1

    _ax1.set_xlim(myoptions.xaxis_label_start or _xmin, _xmax)
    _ax1.xaxis.set_major_locator(ticker.MultipleLocator(myoptions.xaxis_major_ticks_spacing))
    _ax1.xaxis.set_minor_locator(ticker.MultipleLocator(myoptions.xaxis_minor_ticks_spacing))

    if myoptions.debug:
        print(f"Debug: X-axis1: {_xmin}-{_xmax}")

    if myoptions.aminoacids:
        _total_frequencies = np.sum(np.abs(new_aa_table), axis=0)
    else:
        _total_frequencies = np.sum(np.abs(new_codon_table), axis=0)

    if myoptions.debug:
        print(f"Debug: len(_total_frequencies)={len(_total_frequencies)}, _total_frequencies={str(_total_frequencies.to_list())}")

    if myoptions.aminoacids:
        _ax1.xaxis.set_tick_params(labelsize=14)
        _ax1.tick_params(axis='x', which='both', labelsize=14)
        _y_ticks = np.arange(len(amino_acids))
        _ax1.set_yticks(_y_ticks)
        _ax1.set_yticklabels(amino_acids, fontsize=14)
    else:
        _ax1.xaxis.set_tick_params(labelsize=14)
        _ax1.tick_params(axis='x', which='both', labelsize=14)
        _y_ticks = np.arange(len(codons_whitelist))
        _ax1.set_yticks(_y_ticks)
        _ax1.set_yticklabels([_pairs[0] + ' (' + _pairs[1] + ')  ' if _pairs[1] not in ('INS', 'DEL') else _pairs[0] + ' (' + _pairs[1] + ')' for _pairs in final_sorted_whitelist], fontsize=8)
        _ax1.tick_params(axis='y', which='major', labelsize=8)

    if myoptions.xaxis_bins:
        plt.locator_params(axis='x', nbins=myoptions.xaxis_bins)

    _ax1.grid(True, linestyle='--', alpha=0.3, color='gray')

    _ax2 = None
    if not myoptions.disable_2nd_Y_axis:
        _ax2 = _ax1.twinx()
        _ax2.set_xlim(_xmin, _xmax)
        _ax2.set_ylim(0, 1)

        x1, x2 = _ax2.get_xlim()
        _ax2.set_ylabel(f'Cumulative frequency of mutations above threshold {myoptions.threshold:.2%} per codon', fontsize=12)
        _ax1.figure.canvas.draw()
        _ax2.figure.canvas.draw()

        if myoptions.aminoacids:
            _ax2.bar(unique_aa_padded_positions, _total_frequencies, color='black', alpha=0.5, width=0.8, align='center')
        else:
            _ax2.bar(unique_padded_codon_positions, _total_frequencies, color='black', alpha=0.5, width=0.8, align='center')

        x1, x2 = _ax2.get_xlim()

    if myoptions.aminoacids:
        _table = new_aa_table
    else:
        _table = new_codon_table
    if myoptions.debug:
        print(f"Debug: Printing the _table: {str(_table)}")

    if myoptions.debug:
        print("Debug: Printing the _table:")
        print(_table)

    return _figure, _ax1, _ax2, _ax3, _ax4, _xmin, _xmax


def collect_scatter_data(
    myoptions: typing.Any,
    df: typing.Any, table: typing.Any, outfile_prefix: str, matrix: typing.Any,
    amino_acids: list[str], codons_whitelist2: list[str], padded_position2position: dict[int, int],
    xmin: int, xmax: int,
):
    "Iterate over frequency tables and collect scatter plot data, labels, and colors."

    _real_aa_positions = sorted(padded_position2position.values()) # better extract it dynamically from the table then from the dictionary just in case some values would be discarded while parsing on-the-fly
    if myoptions.debug:
        print(f"Debug: _real_aa_positions={str(_real_aa_positions)}")

    if not myoptions.aminoacids:
        if list(table.index) != codons_whitelist2:
            raise ValueError("Both lists should be equal: {}={}, {}={}".format('table.index', str(table.index), 'codons_whitelist2', codons_whitelist2))
    else:
        if list(table.index) != amino_acids:
            raise ValueError("Both lists should be equal: {}={}, {}={}".format('table.index', str(table.index), 'amino_acids', amino_acids))

    _used_colors = set()
    _norm, _cmap, _colors = get_colormap(myoptions, myoptions.colormap)
    _labels: list[str] = []
    _label_padded_positions: list[str] = []
    _label_codon_positions: list[str] = []
    _label_original_amino_acids: list[str] = []
    _label_new_amino_acids: list[str] = []
    _label_cumulative_frequencies: list[str] = []
    _label_codon_frequencies: list[str] = []
    _label_observed_codon_counts: list[typing.Any] = []
    _label_observed_codon_count_sum: list[typing.Any] = []
    _label_total_codons_per_site: list[typing.Any] = []
    _label_scores: list[typing.Any] = []

    _html_labels: list[str] = []
    _mutations: list[str] = []
    _circles_bokeh: list[tuple[typing.Any, ...]] = []
    _circles_matplotlib: list[tuple[typing.Any, ...]] = []
    _markers: list[tuple[typing.Any, ...]] = []
    _dots: list[tuple[typing.Any, ...]] = []
    _warn_once: list[int] = []
    _matrix_values: set[int] = set()

    if myoptions.aminoacids:
        _outfilename = outfile_prefix + '.aa.frequencies.colors.tsv'
    else:
        _outfilename = outfile_prefix + '.codon.frequencies.colors.tsv'
    with open(_outfilename, 'w', encoding="utf-8") as _color_file:
        print(f"Info: Writing into {_outfilename}")

        # the tables were constructed with the following in build_frequency_tables()
        #     _new_aa_table = pd.DataFrame(Decimal(0), index=_amino_acids, columns=_unique_padded_aa_positions)
        #     _new_codon_table = pd.DataFrame(Decimal(0), index=_codons_whitelist2, columns=_unique_padded_codon_positions)
        _pos_to_old_codon = df.groupby('padded_position')['original_codon'].first().to_dict()
        _pos_to_old_aa = df.groupby('padded_position')['original_aa'].first().to_dict()
        _mut_col = 'mutant_aa' if myoptions.aminoacids else 'mutant_codon'
        _df_indexed = df.set_index(['padded_position', _mut_col])
        # Pre-build a dict for O(1) per-cell lookup instead of repeated .loc[] calls
        _df_groups: dict = {}
        for _key, _sub in _df_indexed.groupby(level=[0, 1]):
            _df_groups[_key] = _sub
        # Clear the score cache so each plot run starts fresh (matrix may differ)
        _score_cache.clear()
        for i, _some_codon_or_aa in enumerate(table.index): # so _some_codon_or_aa contains the index specified when the table was constructed
            for j, _padded_position in enumerate(table.columns): # so _aa_position contains the real aa_position
                if not (xmin <= _padded_position <= xmax):
                    continue
                if myoptions.debug:
                    print(f"Debug: i: {str(i)}, j: {str(j)}, _padded_position column: {str(_padded_position)}")
                try:
                    _aa_position = padded_position2position[_padded_position]
                except KeyError:
                    continue
                if myoptions.debug:
                    print(f"Debug0: _padded_positions (typically will not be contiguous and will contain multiplicates): {str(sorted(list(table.columns)))}{os.linesep}")
                    print(f"Debug0:     _aa_positions (typically will not be contiguous and will contain multiplicates): {sorted(padded_position2position.values())}{os.linesep}")
                _frequency = table.loc[_some_codon_or_aa, _padded_position]
                if myoptions.debug and _frequency:
                    print(f"Debug0: _padded_position={_padded_position}, _aa_position={_aa_position}, _some_codon_or_aa={_some_codon_or_aa}, _frequency={str(_frequency)}")
                try:
                    _old_codon = _pos_to_old_codon[_padded_position]
                except KeyError:
                    if _frequency and myoptions.debug:
                        print(f"Debug0b: _padded_position={_padded_position}, _some_codon_or_aa={_some_codon_or_aa}, _frequency={str(_frequency)}")
                    if _padded_position not in _warn_once:
                        sys.stderr.write(f"Warning: Cannot determine original codon for position {_padded_position}, seems missing from input TSV{os.linesep}")
                        _warn_once.append(_padded_position)
                    continue
                _codon_on_input, _old_codon_or_aa, _new_codon_or_aa = resolve_codon_or_aa(myoptions, _old_codon, _some_codon_or_aa)
                _score = None
                if not np.abs(Decimal(_frequency)) < myoptions.threshold:
                    if myoptions.column_with_frequencies in ['neutralized_parent_difference', 'escape_parent_difference']:
                        _size, _color = adjust_size_and_color_neutralized_escape(Decimal(_frequency), _codon_on_input)
                        _score = -12
                    elif myoptions.column_with_frequencies in ['weighted_diff_escape_neutralized']:
                        _size, _color = adjust_size_and_color_weighted(Decimal(_frequency))
                        _score = -12
                    else:
                        _score, _size, _color = adjust_size_and_color(myoptions, Decimal(_frequency), _codon_on_input, _old_codon, _some_codon_or_aa, _old_codon_or_aa, _new_codon_or_aa, matrix, _norm, _colors)
                        _matrix_values.add(_score)
                    if myoptions.debug:
                        print(f"Debug: Padded AA position: {_padded_position}, Real AA position: {_padded_position}, observed codon: {_some_codon_or_aa}, _frequency: {_frequency}, _size: {_size}, color: {_color}")
                    _bokeh_size = float(np.sqrt(np.abs(_size)) * 100) if myoptions.bokeh_sqrt_size else float(np.abs(_size) * 100)
                    if myoptions.aminoacids:
                        if _score < 0:
                            _circles_bokeh.append((_aa_position, _some_codon_or_aa, _bokeh_size, 'circle', _color, 0.5, _score, _aa_position, _padded_position))
                        else:
                            _circles_bokeh.append((_aa_position, _some_codon_or_aa, _bokeh_size, 'hex', _color, 0.5, _score, _aa_position, _padded_position))
                    else:
                        if _score < 0:
                            _circles_bokeh.append((_aa_position, _some_codon_or_aa + ' (' + alt_translate(_some_codon_or_aa) + ')', _bokeh_size, 'circle_x', _color, 0.5, _score, _aa_position, _padded_position))
                        else:
                            _circles_bokeh.append((_aa_position, _some_codon_or_aa + ' (' + alt_translate(_some_codon_or_aa) + ')', _bokeh_size, 'hex', _color, 0.5, _score, _aa_position, _padded_position))
                    if _score < 0:
                        _circles_matplotlib.append((_aa_position, i, float(np.abs(_size) * 5000), 'circle_x', _color, 0.5, _score, _aa_position, _padded_position)) # 'o' circle shape
                        _markers.append((_aa_position, i, 1, 'dot', 'black', 0.5))
                    else:
                        _circles_matplotlib.append((_aa_position, i, float(np.abs(_size) * 5000), 'circle', _color, 0.5, _score, _aa_position, _padded_position)) # 'h' hex shape
                        _markers.append((_aa_position, i, 1, 'circle', 'black', 0.5))
                else:
                    _size, _color = 0, 'black'
                    _score = get_score(myoptions, matrix, _codon_on_input, _old_codon_or_aa, _new_codon_or_aa)
                    if getattr(myoptions, 'show_invisible_placeholder_dots', False):
                        _dots.append((_padded_position, i, _size, 'dot', _color, 0.5, _score))
                    if myoptions.debug:
                        print(f"Debug: Invisible dot. Real AA position: {_padded_position}, observed codon: {_some_codon_or_aa}, _frequency: {_frequency}, _size: {_size}, color: {_color}")

                if not np.abs(Decimal(_frequency)) < myoptions.threshold:
                    if _padded_position not in _warn_once:
                        try:
                            _base_df = _df_groups.get((_padded_position, _some_codon_or_aa), pd.DataFrame())
                        except KeyError:
                            _base_df = pd.DataFrame()

                        if not _base_df.empty:
                            _sub_df = _base_df[_base_df[myoptions.column_with_frequencies] >= myoptions.threshold]
                        else:
                            _sub_df = pd.DataFrame()

                        _frequencies = [Decimal(x) for x in _sub_df[myoptions.column_with_frequencies].to_list()] if not _sub_df.empty else []

                        try:
                            _old_amino_acid = _pos_to_old_aa[_padded_position]
                        except KeyError:
                            print(f"Error: Cannot slice {str(df.loc[df['padded_position'] == _padded_position]['original_aa'].to_list())}")
                            _old_amino_acid = _pos_to_old_aa[_padded_position]

                        _new_codons = _sub_df['mutant_codon'].to_list() if not _sub_df.empty and 'mutant_codon' in _sub_df.columns else []
                        if myoptions.aminoacids:
                            if _frequency != table.at[_some_codon_or_aa, _padded_position]:
                                raise ValueError(f"Values _frequency={_frequency} and table.at[_some_codon_or_aa, _padded_position]={table.at[_some_codon_or_aa, _padded_position]} should be equal")
                            if len(_new_codons) != len(_frequencies):
                                raise ValueError(f"len(_new_codons) != len(_frequencies), specifically: {len(_new_codons)} != {len(_frequencies)}")
                            if 'observed_codon_count' in df.columns.values:
                                _observed_codon_counts = _sub_df['observed_codon_count'].to_list() if not _sub_df.empty else []
                                _total_codons_per_site_list = _sub_df['total_codons_per_site'].to_list() if not _sub_df.empty else []
                                if _total_codons_per_site_list:
                                    _total_codons_per_site = _total_codons_per_site_list[0]
                            else:
                                _observed_codon_counts = []
                                _total_codons_per_site = 0
                            _observed_codon_count_sum = sum(_observed_codon_counts)

                            if not _frequency < myoptions.threshold:
                                _labels.append(f"Padded position: {_padded_position}\nPosition: {_aa_position}\nOriginal Amino Acid: {_old_amino_acid} ({_old_codon})\nNew Amino Acid: {_some_codon_or_aa} {_new_codons}\n{myoptions.matrix} score: {_score}\nCumulative Frequency: {sum(_frequencies):.6f}\nCodon Frequencies: {[f'{x:.6f}' for x in _frequencies]}\nObserved codon counts: {_observed_codon_counts}\nObserved codon counts sum: {_observed_codon_count_sum}\nTotal codons per site: {_total_codons_per_site}")
                                _label_padded_positions.append(f"{_padded_position}")
                                _label_codon_positions.append(f"{_aa_position}")
                                _label_original_amino_acids.append(f"{_old_amino_acid} ({_old_codon})")
                                _label_new_amino_acids.append(f"{_some_codon_or_aa} {_new_codons}")
                                _label_cumulative_frequencies.append(f"{sum(_frequencies):.6f}")
                                _label_codon_frequencies.append(f"{[f'{x:.6f}' for x in _frequencies]}")
                                if 'observed_codon_count' in df.columns.values:
                                    _label_observed_codon_counts.append(f"{_observed_codon_counts}")
                                    _label_observed_codon_count_sum.append(f"{sum(_observed_codon_counts)}")
                                    _label_total_codons_per_site.append(f"{_total_codons_per_site}")
                                else:
                                    _label_observed_codon_counts.append(_observed_codon_counts)
                                    _label_observed_codon_count_sum.append(sum(_observed_codon_counts))
                                    _label_total_codons_per_site.append(_total_codons_per_site)
                                _label_scores.append(_score)
                                _html_labels.append(f"Padded position: {_padded_position}<br>Position: {_aa_position}<br>Original Amino Acid: {_old_amino_acid} ({_old_codon})<br>New Amino Acid: {_some_codon_or_aa} {_new_codons}<br>{myoptions.matrix} score: {_score}<br>Cumulative Frequency: {sum(_frequencies):.6f}<br>Codon Frequencies: {[f'{x:.6f}' for x in _frequencies]}<br>Observed codon counts: {_observed_codon_counts}<br>Observed codon count sum: {_observed_codon_count_sum}<br>Total codons per site: {_total_codons_per_site}")
                                _mutations.append(f"{_old_amino_acid}{_aa_position}{_some_codon_or_aa}")
                        else:
                            if _frequency != table.at[_some_codon_or_aa, _padded_position]:
                                raise ValueError(f"Values _frequency={_frequency} and table.at[_some_codon_or_aa, _padded_position]={table.at[_some_codon_or_aa, _padded_position]} should be equal")
                            if 'observed_codon_count' in df.columns.values:
                                _observed_codon_counts = _sub_df['observed_codon_count'].to_list() if not _sub_df.empty else []
                                if _observed_codon_counts and len(_observed_codon_counts) < 2:
                                    _observed_codon_counts = _observed_codon_counts[0]
                                    _observed_codon_count_sum = _observed_codon_counts
                                else:
                                    _observed_codon_count_sum = sum(_observed_codon_counts)
                                _total_codons_per_site_list = _sub_df['total_codons_per_site'].to_list() if not _sub_df.empty else []
                                if _total_codons_per_site_list:
                                    _total_codons_per_site = _total_codons_per_site_list[0]
                            else:
                                _observed_codon_counts = []
                                _total_codons_per_site = 0

                            _observed_aminoacids = df.loc[df['padded_position'] == _padded_position]['mutant_aa'].to_list()
                            if myoptions.debug:
                                print(f"Info: {len(_observed_aminoacids)} aa residues observed in position {_aa_position}:{os.linesep} {str(df.loc[df['padded_position'] == _padded_position][0:])}{os.linesep}")

                            try:
                                # make sure we dot not fetch also INSertions, which could be even multiple rows in addition to the row with changed_codon, especially if the reference protein and is padded on the right with dashes
                                # try to switch to 'padded_position' instead of 'position' to fetch a row from Pandas
                                _new_amino_acid_list = _base_df['mutant_aa'].to_list() if not _base_df.empty and 'mutant_aa' in _base_df.columns else []
                                if _new_amino_acid_list:
                                    _new_amino_acid = _new_amino_acid_list[0]
                                else:
                                    _new_amino_acid = None
                                if myoptions.debug and _new_amino_acid:
                                    print(f"Info: Parsed existing _new_amino_acid '{_new_amino_acid}' for codon position {_aa_position}")
                            except Exception as exc:
                                _new_amino_acid = None

                            if _new_amino_acid:
                                try:
                                    # make sure we dot not fetch also INSertions, which could be even multiple rows in addition to the row with changed_codon
                                    _some_freq_list = _base_df[myoptions.column_with_frequencies].to_list() if not _base_df.empty and myoptions.column_with_frequencies in _base_df.columns else []
                                    _some_frequency = Decimal(_some_freq_list[0]) if _some_freq_list else Decimal('0.00000000009')
                                except (IndexError, ValueError, TypeError):
                                    _some_frequency = Decimal('0.00000000009')
                                if _some_frequency != _frequency and (_some_frequency != Decimal('0.00000000009') and _frequency != 0):
                                    # 334	333	R	R	0.432432	AGG	AGA	16	37 # .frequencies.tsv
                                    # 338	333	INS	R	0.027027	---	AGA	1	37 # .frequencies.tsv
                                    # 344	333	INS	R	0.648649	---	AGA	24	37 # .frequencies.tsv
                                    #
                                    # 344	333	INS	R	0.648649	---	AGA	24	37 # .frequencies.unchanged_codons.tsv
                                    raise ValueError(f"Frequency new_codon_table.at[_some_codon_or_aa, _padded_position]={_frequency} _some_codon_or_aa={_some_codon_or_aa}, _padded_position={_padded_position} not same as df.loc[(df['padded_position'] == _padded_position) & (df['mutant_codon'] == _some_codon_or_aa)][myoptions.column_with_frequencies].to_list()[0]={_some_frequency}. Are multiple rows matching? We picked just the first one: df.loc[(df['padded_position'] == _padded_position) & (df['mutant_codon'] == _some_codon_or_aa)][myoptions.column_with_frequencies].to_list()={_some_freq_list} for _old_codon={_old_codon}, _old_amino_acid={_old_amino_acid}")

                                if _old_amino_acid and not _frequency < myoptions.threshold and _size:
                                    _label_scores.append(_score)

                                    if myoptions.column_with_frequencies == 'neutralized_parent_difference':
                                        _labels.append(f"Padded position: {_padded_position}\nPosition: {_aa_position}\nOriginal Codon: {_old_codon} ({_old_amino_acid})\nNew Codon: {_some_codon_or_aa} ({_new_amino_acid})\n{myoptions.matrix} score: {_score}\nDifference neutralized2parent: {_frequency:.6f}")
                                        _html_labels.append(f"Padded position: {_padded_position}<br>Position: {_aa_position}<br>Original Codon: {_old_codon} ({_old_amino_acid})<br>New Codon: {_some_codon_or_aa} ({_new_amino_acid})<br>{myoptions.matrix} score: {_score}<br>Difference neutralized2parent: {_frequency:.6f}")
                                    elif myoptions.column_with_frequencies == 'escape_parent_difference':
                                        _labels.append(f"Padded position: {_padded_position}\nPosition: {_aa_position}\nOriginal Codon: {_old_codon} ({_old_amino_acid})\nNew Codon: {_some_codon_or_aa} ({_new_amino_acid})\n{myoptions.matrix} score: {_score}\nDifference escape2parent: {_frequency:.6f}")
                                        _html_labels.append(f"Padded position: {_padded_position}<br>Position: {_aa_position}<br>Original Codon: {_old_codon} ({_old_amino_acid})<br>New Codon: {_some_codon_or_aa} ({_new_amino_acid})<br>{myoptions.matrix} score: {_score}<br>Difference escape2parent: {_frequency:.6f}")
                                    elif myoptions.column_with_frequencies == 'weighted_diff_escape_neutralized':
                                        _labels.append(f"Padded position: {_padded_position}\nPosition: {_aa_position}\nOriginal Codon: {_old_codon} ({_old_amino_acid})\nNew Codon: {_some_codon_or_aa} ({_new_amino_acid})\n{myoptions.matrix} score: {_score}\nWeighted difference escape2neutralized: {_frequency:.6f}")
                                        _html_labels.append(f"Padded position: {_padded_position}<br>Position: {_aa_position}<br>Original Codon: {_old_codon} ({_old_amino_acid})<br>New Codon: {_some_codon_or_aa} ({_new_amino_acid})<br>{myoptions.matrix} score: {_score}<br>Weighted difference escape2neutralized: {_frequency:.6f}")
                                    elif myoptions.column_with_frequencies == 'frequency':
                                        _labels.append(f"Padded position: {_padded_position}\nPosition: {_aa_position}\nOriginal Codon: {_old_codon} ({_old_amino_acid})\nNew Codon: {_some_codon_or_aa} ({_new_amino_acid})\n{myoptions.matrix} score: {_score}\nFrequency: {_frequency:.6f}\nObserved codon count: {_observed_codon_counts}\nTotal codons per site: {_total_codons_per_site}")
                                        _html_labels.append(f"Padded position: {_padded_position}<br>Position: {_aa_position}<br>Original Codon: {_old_codon} ({_old_amino_acid})<br>New Codon: {_some_codon_or_aa} ({_new_amino_acid})<br>{myoptions.matrix} score: {_score}<br>Frequency: {_frequency:.6f}<br>Observed codon counts: {_observed_codon_counts}<br>Total codons per site: {_total_codons_per_site}")
                                        _label_padded_positions.append(f"{_padded_position}")
                                        _label_codon_positions.append(f"{_aa_position}")
                                        _label_original_amino_acids.append(f"{_old_amino_acid} ({_old_codon})")
                                        _label_new_amino_acids.append(f"{_new_amino_acid} ({_some_codon_or_aa})")
                                        _label_cumulative_frequencies.append(f"{_frequency:.6f}")
                                        _label_codon_frequencies.append(f"{_frequency:.6f}")
                                    else:
                                        _labels.append(f"Padded position: {_padded_position}\nPosition: {_aa_position}\nOriginal Codon: {_old_codon} ({_old_amino_acid})\nNew Codon: {_some_codon_or_aa} ({_new_amino_acid})\n{myoptions.matrix} score: {_score}\nFrequency: {_frequency:.6f}\nObserved codon count: {_observed_codon_counts}\nTotal codons per site: {_total_codons_per_site}")
                                        _html_labels.append(f"Padded position: {_padded_position}<br>Position: {_aa_position}<br>Original Codon: {_old_codon} ({_old_amino_acid})<br>New Codon: {_some_codon_or_aa} ({_new_amino_acid})<br>{myoptions.matrix} score: {_score}<br>Frequency: {_frequency:.6f}<br>Observed codon counts: {_observed_codon_counts}<br>Total codons per site: {_total_codons_per_site}")
                                        _label_padded_positions.append(f"{_padded_position}")
                                        _label_codon_positions.append(f"{_aa_position}")
                                        _label_original_amino_acids.append(f"{_old_amino_acid} ({_old_codon})")
                                        _label_new_amino_acids.append(f"{_new_amino_acid} ({_some_codon_or_aa})")
                                        _label_cumulative_frequencies.append(f"{_frequency:.6f}")
                                        _label_codon_frequencies.append(f"{_frequency:.6f}")
                                    if not _frequency < myoptions.threshold:
                                        _mutations.append(f"{_old_codon}{_aa_position}{_some_codon_or_aa}")
                                    else:
                                        _mutations.append('')
                                    if 'observed_codon_count' in df.columns.values:
                                        _label_observed_codon_counts.append(f"{_observed_codon_counts}")
                                        _label_observed_codon_count_sum.append(f"{_observed_codon_count_sum}")
                                        _label_total_codons_per_site.append(f"{_total_codons_per_site}")
                                    else:
                                        _label_observed_codon_counts.append('')
                                        _label_observed_codon_count_sum.append('')
                                        _label_total_codons_per_site.append('')

                        _original_aas = df.loc[(df['padded_position'] == _padded_position)]['original_aa'].to_list()
                        if _original_aas:
                            _original_aa = _original_aas[0]
                        else:
                            _original_aa = None
                        _new_amino_acid = _some_codon_or_aa
                        _mutant_aa = None
                        if myoptions.aminoacids:
                            _mutant_codons = df.loc[(df['padded_position'] == _padded_position) & (df['mutant_aa'] == _new_amino_acid)]['mutant_codon'].to_list()
                        else:
                            _mutant_codon = _some_codon_or_aa
                            _mutant_codons = [_mutant_codon]
                            if _some_codon_or_aa == '---':
                                _mutant_aa = 'DEL'
                            elif _some_codon_or_aa == 'INS':
                                _mutant_aa = 'INS'
                            elif _some_codon_or_aa == 'DEL':
                                _mutant_aa = 'DEL'
                            else:
                                _mutant_aa = alt_translate(_some_codon_or_aa)
                            _mutant_codons = df.loc[(df['padded_position'] == _padded_position) & (df['mutant_codon'] == _some_codon_or_aa)]['mutant_codon'].to_list()
                        if _original_aa and not _frequency < myoptions.threshold and _size:
                            if len(_mutant_codons) > 1:
                                _color_file.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}{}".format(_aa_position, _old_codon, str(_mutant_codons), _original_aa, _some_codon_or_aa, f'{_frequency:.6f}', _color, _score, os.linesep))
                            elif not myoptions.aminoacids:
                                _color_file.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}{}".format(_aa_position, _old_codon, _mutant_codon, _original_aa, _mutant_aa, f'{_frequency:.6f}', _color, _score, os.linesep))
                            elif _mutant_codons:
                                _color_file.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}{}".format(_aa_position, _old_codon, _mutant_codons[0], _original_aa, _some_codon_or_aa, f'{_frequency:.6f}', _color, _score, os.linesep))
                        elif myoptions.debug:
                            if myoptions.aminoacids:
                                sys.stderr.write("Debug: Skipped line for _aa_position={} _original_aa={} _old_codon={} _mutant_codons={} _some_codon_or_aa={} _frequency={} score={} myoptions.threshold={}\n".format(_aa_position, _original_aa, _old_codon, str(_mutant_codons), _some_codon_or_aa, f'{_frequency:.6f}', _score, myoptions.threshold))
                            else:
                                sys.stderr.write("Debug: Skipped line for _aa_position={} _original_aa={} _old_codon={} _mutant_codons={} _some_codon_or_aa={} _frequency={} score={} myoptions.threshold={}\n".format(_aa_position, _original_aa, _old_codon, _mutant_codon, _mutant_aa, f'{_frequency:.6f}', _score, myoptions.threshold))
                    _used_colors.add(_color)

    print("Info: The following values were collected from matrix %s based on the actual data (some values from matrix might not be needed for your data, hence are not listed here): %s . Range spans %d values (before symmetrization)." % (myoptions.matrix, str(sorted(_matrix_values)), abs(min(_matrix_values)) + 1 + max(_matrix_values)))
    if myoptions.debug:
        print(f"Debug: {len(_used_colors)} _used_colors used: {str(_used_colors)}")

    if myoptions.debug:
        print(f"Debug: {len(_labels)} _labels={str(_labels)}")
        print(f"Debug: {len(_html_labels)} _html_labels={str(_html_labels)}")
        print(f"Debug: {len(_circles_bokeh)} _circles_bokeh={str(_circles_bokeh)}")
        print(f"Debug: {len(_markers)} _markers={str(_markers)}")
        print(f"Debug: {len(_dots)} _dots={str(_dots)}")

    return (
        _norm, _cmap, _colors, _used_colors, _matrix_values,
        _labels, _html_labels, _mutations,
        _circles_bokeh, _circles_matplotlib, _markers, _dots, _label_padded_positions,
        _label_codon_positions, _label_original_amino_acids, _label_new_amino_acids,
        _label_cumulative_frequencies, _label_codon_frequencies,
        _label_observed_codon_counts, _label_observed_codon_count_sum,
        _label_total_codons_per_site, _label_scores,
    )


def pretty_print_bokeh_html(filename):
    """Pretty-print the JSON block inside a Bokeh-generated HTML file."""
    if not os.path.exists(filename):
        return
    with open(filename, 'r', encoding='utf-8') as f:
        content = f.read()

    # Bokeh embeds JSON in a script tag. Use a regex to find it.
    # Pattern: <script type="application/json" id="...">JSON_HERE</script>
    pattern = r'(<script type="application/json"[^>]*>)(.*?)(</script>)'

    def replacer(match):
        prefix = match.group(1)
        json_str = match.group(2)
        suffix = match.group(3)
        try:
            data = json.loads(json_str)
            pretty_json = json.dumps(data, indent=4)
            return f"{prefix}\n{pretty_json}\n{suffix}"
        except Exception:
            return match.group(0)

    new_content = re.sub(pattern, replacer, content, flags=re.DOTALL)

    if new_content != content:
        with open(filename, 'w', encoding='utf-8') as f:
            f.write(new_content)


def render_bokeh(
    myoptions,
    outfile_prefix, xmin, xmax, amino_acids, final_sorted_whitelist,
    circles_bokeh, labels, html_labels, mutations, label_padded_positions,
    label_codon_positions, label_original_amino_acids, label_new_amino_acids,
    label_cumulative_frequencies, label_codon_frequencies,
    label_observed_codon_counts, label_observed_codon_count_sum,
    label_total_codons_per_site, label_scores,
    title_data, xlabel,
    matrix_name, colors, norm, cmap, padded_position2position,
):
    """Build and save the interactive Bokeh HTML scatter plot.

    Parameters
    ----------
    myoptions : argparse.Namespace
        Parsed command-line options.  Relevant flags: ``aminoacids``,
        ``matrix``, ``bokeh_sqrt_size``.
    outfile_prefix : str
        Path prefix for the output file.  ``'.html'`` is appended automatically.
    xmin, xmax : int
        Genomic / protein position range for the x-axis.
    amino_acids : list[str]
        Y-axis categories when ``myoptions.aminoacids`` is True.
    final_sorted_whitelist : list[tuple[str, str]]
        Y-axis codon-change pairs when running in codon mode.
    circles_bokeh : list[tuple]
        One tuple per scatter point with layout::

            (x, y, size, marker, colour_hex, alpha, score, aa_position, padded_position)
              0  1    2      3           4       5      6

        ``colour_hex`` is a hex string already produced by the matplotlib
        ``cmap(norm(score))`` pipeline in ``collect_scatter_data``.
        ``score`` is the integer substitution-matrix score for that mutation.
    labels, html_labels, mutations : list[str]
        Per-point annotation strings.
    label_codon_positions … label_scores : list
        Per-point tooltip field values (see TOOLTIPS definition in the body).
    title_data, xlabel : str
        Figure title and x-axis label.
    matrix_name : str
        Human-readable name of the substitution matrix (e.g. ``'BLOSUM62'``,
        ``'amino_acid_changes'``).  Used as the colorbar title.
    colors : list or None
        The raw colour list returned by ``get_colormap()``.  Its **length**
        encodes the designed score range: a list of N entries covers scores
        ``-N//2 .. +N//2`` (e.g. 39 entries → scores -19 … +19).  Used only
        as a fallback when ``norm`` is ``None`` (continuous colormaps).
    norm : matplotlib.colors.BoundaryNorm or None
        Normalisation object returned by ``get_colormap()``.  Present only for
        the discrete palettes ``amino_acid_changes`` (BoundaryNorm -19…+18)
        and ``dkeenan`` (BoundaryNorm -13…+12).  ``None`` for all other
        colormaps.
    cmap : matplotlib.colors.Colormap or None
        Colormap object returned by ``get_colormap()``.

    Colorbar implementation notes
    ------------------------------
    Bokeh's ``LinearColorMapper`` maps the continuous interval ``[low, high]``
    linearly across the palette list.  With a palette of N discrete colours we
    need each colour band to be exactly 1 score-unit wide and centred on its
    integer score value so that the tick label (an integer) falls at the visual
    centre of the band, and so that score 0 always maps to the correct colour
    (yellow in the ``amino_acid_changes`` palette).

    The key relationships are:

    * ``low  = -half - 0.5``  (half = N // 2)
    * ``high = +half + 0.5``
    * Total range = N units → band width = 1.0 exactly.
    * Score *s* occupies the band ``[s - 0.5, s + 0.5]``, centred on *s*.
    * Ticks are placed at plain integers via ``FixedTicker``.

    Why ``cmap(norm(s))`` rather than indexing ``colors`` directly:
        The ``amino_acid_changes`` palette has 39 colour entries but
        ``BoundaryNorm(np.arange(-19, 19, 1), 39)`` defines only 38 boundary
        edges (37 bins between them) and then rescales those into 39 colour
        slots.  As a result, colour index 20 (yellow, ``#ffff00``) maps to
        score 0 via the norm, while a naive ``colors[19]`` (index 19) would
        return the wrong brown colour (``#9c644b``).  Evaluating
        ``cmap(norm(s))`` replicates the exact same rescaling that matplotlib
        uses for the scatter points, guaranteeing that the colorbar colours
        match the plotted circles for every score value.

    Fallback for continuous colormaps (``norm`` is ``None``):
        The palette is sampled linearly across the cmap at ``N`` equally spaced
        positions.  If ``cmap`` is also ``None`` a flat grey palette is used.

    Circle size scaling
    -------------------
    Bokeh's ``scatter(size=...)`` interprets the value as a **diameter in
    screen pixels**, so area ∝ size².  A naive linear mapping
    (``size = frequency × 100``) therefore makes large frequencies appear
    disproportionately large.

    Two modes are available via ``--disable-bokeh-sqrt-size``:

    sqrt mode (default, ``myoptions.bokeh_sqrt_size = True``)::

        _bokeh_size = float(np.sqrt(np.abs(_size)) * 100)

    Diameter ∝ sqrt(frequency), so area ∝ frequency — matching human
    perception and matching the visual appearance of the matplotlib figure.

    Linear mode (``--disable-bokeh-sqrt-size``)::

        _bokeh_size = float(np.abs(_size) * 100)

    Diameter ∝ frequency, area ∝ frequency².

    Relationship to matplotlib scatter:
        ``matplotlib.scatter(s=...)`` interprets *s* as **area in points²**, so
        ``s = frequency × 5000`` already gives perceived radius ∝
        sqrt(frequency) without an explicit sqrt.  The Bokeh sqrt mode mirrors
        this behaviour.  Use ``--disable-bokeh-sqrt-size`` to revert to
        diameter-proportional scaling in Bokeh only.

    X-axis tick layout
    ------------------
    Major ticks (with labels) and minor ticks (unlabelled) are placed using
    ``FixedTicker``, mirroring the ``ticker.MultipleLocator`` behaviour used
    by ``render_matplotlib``.

    The first major tick is the smallest multiple of
    ``myoptions.xaxis_major_ticks_spacing`` that is **≥ xmin** (or ≥
    ``myoptions.xaxis_label_start`` when that option is set).  This means
    the first label always aligns with the actual data range rather than
    being back-projected to zero.

    Example: if the first data point is at position 331 and
    ``xaxis_major_ticks_spacing=10``, the first label is at **340**::

        _first_major = int(_x_start) + (-int(_x_start) % _major_spacing)
        # 331 + (-331 % 10) = 331 + 9 = 340

    Minor ticks are placed at every ``myoptions.xaxis_minor_ticks_spacing``
    position starting from the first multiple of that spacing that is ≥ xmin,
    excluding positions already occupied by a major tick.  They are rendered
    via a second ``LinearAxis`` overlay with invisible labels and a smaller
    tick length.

    Relevant options (shared with ``render_matplotlib``):

    * ``--x-axis-major-ticks-spacing``  (default: 10)
    * ``--x-axis-minor-ticks-spacing``  (default: 5)
    * ``--x-axis-label-start``          (default: 0, meaning use xmin)
    """
    if circles_bokeh:
        _circles_x, _circles_y, _circles_size, _circles_marker, _circles_color, _circles_alpha, _circles_score, _circles_aa_pos, _circles_padded_pos = zip(*circles_bokeh)
    else:
        _circles_x, _circles_y, _circles_size, _circles_marker, _circles_color, _circles_alpha, _circles_score, _circles_aa_pos, _circles_padded_pos = [], [], [], [], [], [], [], [], []

    _mysource = bokeh.models.ColumnDataSource(data=dict(
        x=_circles_x,
        y=_circles_y,
        s=_circles_size,
        m=_circles_marker,
        c=_circles_color,
        a=_circles_alpha,
        score=_circles_score,
        aaposition=_circles_aa_pos,
        label=labels,
        label1=label_padded_positions,
        label2=label_original_amino_acids,
        label3=label_new_amino_acids,
        label4=label_cumulative_frequencies,
        label5=label_codon_frequencies,
        label6=label_observed_codon_counts,
        label7=label_observed_codon_count_sum,
        label8=label_total_codons_per_site,
        label9=label_scores,
        label10=label_codon_positions,
        padded_pos=_circles_padded_pos,
        mutation=mutations,
    ))

    if myoptions.aminoacids:
        TOOLTIPS = [
            ("Padded Codon Position", "@label1"),
            ("Codon Position", "@label10"),
            ("Original Amino Acid", "@label2"),
            ("New Amino Acid", "@label3"),
            ("Cumulative Frequency", "@label4"),
            ("Codon Frequencies", "@label5"),
            ("Observed codon counts", "@label6"),
            ("Observed codon count sum", "@label7"),
            ("Total codons per site", "@label8"),
            (f"{myoptions.matrix} score", "@label9"),
        ]
    else:
        TOOLTIPS = [
            ("Padded Codon Position", "@label1"),
            ("Codon Position", "@label10"),
            ("Original Amino Acid", "@label2"),
            ("New Amino Acid", "@label3"),
            ("Cumulative Frequency", "@label4"),
            ("Codon Frequencies", "@label5"),
            ("Observed codon count", "@label6"),
            ("Total codons per site", "@label8"),
            (f"{myoptions.matrix} score", "@label9"),
        ]
    if myoptions.aminoacids:
        _p = bokeh.plotting.figure(x_range=(xmin, xmax), y_range=amino_acids, tooltips=TOOLTIPS, title=title_data, x_axis_label=xlabel, y_axis_label='Introduced amino acid changes', width=2000, height=1200, sizing_mode='stretch_width')
    else:
        _p = bokeh.plotting.figure(x_range=(xmin, xmax), y_range=[_pairs[0] + ' (' + _pairs[1] + ')' for _pairs in final_sorted_whitelist], tooltips=TOOLTIPS, title=title_data, x_axis_label=xlabel, y_axis_label='Introduced codon changes', height=1200, width=2000, sizing_mode='stretch_width')

    # Mirror matplotlib's xaxis_major_ticks_spacing / xaxis_minor_ticks_spacing.
    # Major ticks carry labels; minor ticks are unlabelled.
    _x_start = myoptions.xaxis_label_start or xmin
    _major_spacing = myoptions.xaxis_major_ticks_spacing  # default 10
    _minor_spacing = myoptions.xaxis_minor_ticks_spacing  # default 5
    # First major tick: the smallest multiple of _major_spacing that is >= _x_start.
    # This ensures labels are aligned to the data, not to an arbitrary origin.
    _first_major = int(_x_start) + (-int(_x_start) % _major_spacing)
    _major_ticks = list(range(_first_major, int(xmax) + _major_spacing, _major_spacing))
    # First minor tick: smallest multiple of _minor_spacing that is >= xmin,
    # excluding positions already covered by a major tick.
    _first_minor = int(xmin) + (-int(xmin) % _minor_spacing)
    _minor_ticks = [
        x for x in range(_first_minor, int(xmax) + _minor_spacing, _minor_spacing)
        if x not in _major_ticks
    ]
    _p.xaxis.ticker = bokeh.models.FixedTicker(ticks=_major_ticks)
    _p.xaxis.minor_tick_line_color = "black"
    _p.xgrid.minor_grid_line_color = None
    # Attach minor ticks via a second overlay axis sharing the same range
    _p.add_layout(bokeh.models.LinearAxis(
        x_range_name="default",
        ticker=bokeh.models.FixedTicker(ticks=_minor_ticks),
        major_tick_in=4,
        major_tick_out=0,
        major_label_text_font_size="0px",
        axis_line_color=None,
    ), "below")

    _p.xaxis.major_label_orientation = 'vertical'
    _p.axis.axis_label_text_font_size = "12px"
    _p.axis.major_label_text_font_size = "12px"
    _p.scatter(x='x', y='y', size='s', marker='m', color='c', alpha='a', source=_mysource)

    # Build the palette by evaluating cmap(norm(s)) for every integer score in
    # the full theoretical range. This is the same computation used to colour
    # the scatter points, so every colour band is guaranteed to be correct.
    # norm and cmap are matplotlib objects but are only used here at Python
    # build time to produce plain hex strings — nothing matplotlib is passed
    # to Bokeh itself.
    # Guard: norm is only set for amino_acid_changes / dkeenan colormaps.
    # For other colormaps colors may be None or a continuous cmap without a
    # discrete norm, so fall back to sampling cmap linearly.
    _n = len(colors) if colors is not None else 256
    _half = _n // 2
    _score_range = range(-_half, _half + 1)
    if norm is not None and cmap is not None:
        _score_palette = [
            matplotlib.colors.to_hex(cmap(norm(s)))
            for s in _score_range
        ]
    elif cmap is not None:
        _score_palette = [
            matplotlib.colors.to_hex(cmap(i / max(1, _n - 1)))
            for i in range(_n)
        ]
    else:
        _score_palette = ['#aaaaaa'] * _n
    # low/high extended by ±0.5 so each band is exactly 1 unit wide and the
    # integer score tick sits at the centre of its band.
    _color_mapper = bokeh.models.LinearColorMapper(
        palette=_score_palette,
        low=-_half - 0.5,
        high= _half + 0.5,
    )
    _tick_positions = list(_score_range)
    _colorbar = bokeh.models.ColorBar(
        color_mapper=_color_mapper,
        label_standoff=8,
        title=f"{matrix_name} score values",
        title_standoff=10,
        location=(0, 0),
        ticker=bokeh.models.FixedTicker(ticks=_tick_positions),
    )
    _p.add_layout(_colorbar, 'right')

    _p.text(x='x', y='y', text='mutation', text_color='c', text_font_size='10px', text_font_style='bold', source=_mysource)
    _p.title.align = 'center'
    _p.title.text_font_size = '14pt'
    print(f"Info: Writing into {outfile_prefix} + '.html'")
    bokeh.plotting.output_file(outfile_prefix + '.html')
    bokeh.plotting.show(_p)
    pretty_print_bokeh_html(outfile_prefix + '.html')


def render_matplotlib(
    myoptions,
    figure, ax1, ax2, ax3, ax4, outfile_prefix,
    circles_matplotlib, markers, dots, cmap, norm, colors,
    matrix, matrix_name,
    new_aa_table, new_codon_table, df, codons_whitelist2, final_sorted_whitelist,
    calculated_aa_offset, padded_position2position,
):
    """Render the matplotlib scatter figure with hover callbacks and save to PNG/PDF.

    Circle size scaling (matplotlib):
        matplotlib's scatter(s=...) interprets s as **area in points²**.
        Sizes are stored as float(np.abs(_size) * 5000) where _size is the raw
        frequency value returned by adjust_size_and_color().
        Because s is area, perceived radius ∝ sqrt(s) ∝ sqrt(frequency * 5000),
        so the visual radius is proportional to sqrt(frequency).
        This gives perceptually linear scaling — a frequency twice as large appears
        ~1.41× wider — without any explicit sqrt transformation.

    Relationship to Bokeh scatter:
        Bokeh's scatter(size=...) interprets size as **diameter in screen pixels**,
        so the same raw frequency * 100 mapping gives diameter ∝ frequency and
        area ∝ frequency², making large frequencies appear disproportionately large.
        By default --disable-bokeh-sqrt-size is enabled, so Bokeh applies sqrt scaling and
        matches the perceptual appearance of this matplotlib figure. Use
        --disable-bokeh-sqrt-size to revert to diameter-proportional scaling in Bokeh.

    X-axis tick layout
    ------------------
    Major ticks (with labels) are placed by ``ticker.MultipleLocator(
    myoptions.xaxis_major_ticks_spacing)`` and minor ticks (unlabelled) by
    ``ticker.MultipleLocator(myoptions.xaxis_minor_ticks_spacing)``.
    ``MultipleLocator`` automatically places the first tick at the smallest
    multiple of the spacing that falls within the view range set by
    ``set_xlim()``, so labels always align with the data.

    Example: if the first data point is at position 331 and
    ``xaxis_major_ticks_spacing=10``, the view starts at 331 and
    ``MultipleLocator(10)`` places the first label at **340** — the first
    multiple of 10 that is ≥ 331.

    ``render_bokeh`` replicates this behaviour with an explicit ceiling
    calculation::

        _first_major = int(_x_start) + (-int(_x_start) % _major_spacing)
        # 331 + (-331 % 10) = 331 + 9 = 340

    Relevant options:

    * ``--x-axis-major-ticks-spacing``  (default: 10)
    * ``--x-axis-minor-ticks-spacing``  (default: 5)
    * ``--x-axis-label-start``          (default: 0, meaning use xmin)
    * ``--x-axis-bins``                 overrides spacing via ``locator_params``

    Hover callback (mplcursors)
    ---------------------------
    The cursor is attached to ``_mpl_scatterplot`` (circles only), NOT to the
    full ``ax1`` axes object.  Attaching to ``ax1`` would also cover the
    ``markers`` and ``dots`` scatter collections; dots intentionally have
    ``size=0`` and represent below-threshold positions — hovering over them
    would produce an IndexError because no above-threshold row exists in ``df``
    for that ``(padded_position, codon/aa)`` combination.

    Inside ``on_add``, point identity is resolved via ``sel.index`` (the integer
    index of the nearest point in the PathCollection), NOT via ``sel.target``
    (which returns the raw mouse cursor position in data coordinates and is
    therefore subject to floating-point imprecision and off-by-one errors when
    converted back to a table column index).  ``circles_matplotlib[sel.index]``
    gives the exact ``(_padded_position, row_index, ...)`` tuple that was used
    to plot the point, guaranteeing that the hover annotation always refers to
    the correct position and codon/amino-acid.
    """

    if circles_matplotlib:
        cm_x, cm_y, cm_s, _, _, _, cm_c, _, _ = zip(*circles_matplotlib)
        _mpl_scatterplot = ax1.scatter(cm_x, cm_y, marker='o', s=cm_s, alpha=0.5, c=cm_c, cmap=cmap, norm=norm)
    else:
        _mpl_scatterplot = ax1.scatter([], [], marker='o', s=[], alpha=0.5, c=[], cmap=cmap, norm=norm)

    _colorbar = figure.colorbar(_mpl_scatterplot, cax=ax3, label=f"{matrix_name} score values", location='right', pad=-0.1, alpha=0.5)
    _colorbar.ax.set_yticks(np.arange(-18.5, 18.5, 1), np.arange(-19, 18, 1))
    _colorbar.ax.tick_params(axis='y', which='minor', length=0)
    
    if markers:
        mk_x, mk_y, mk_s, _, _, _ = zip(*markers)
        ax1.scatter(mk_x, mk_y, s=mk_s, marker='x', color='black', alpha=0.5)
    
    if dots:
        dt_x, dt_y, dt_s, _, _, _, _ = zip(*dots)
        ax1.scatter(dt_x, dt_y, s=dt_s, marker='.', color='black', alpha=0.5)

    for _label in ax1.get_xticklabels():
        _label.set_rotation(90)
        _label.set_ha("center")

    _cursor = mplcursors.cursor(_mpl_scatterplot, hover=True)
    if myoptions.aminoacids:
        @_cursor.connect("add")
        def on_add(sel):
            # sel.index is the index into _mpl_scatterplot's data array, i.e.
            # into circles_matplotlib.  Do NOT use sel.target: it returns the
            # raw mouse position in data coordinates, which produces wrong
            # padded_position values when converted back via column arithmetic.
            _pt = circles_matplotlib[sel.index]
            _padded_position = _pt[0]  # stored as _padded_position when appended
            ypos = _pt[1]              # stored as i (row index in new_aa_table)
            _new_amino_acid = new_aa_table.index[ypos]
            _position_in_protein = padded_position2position[_padded_position]
            _frequency = new_aa_table.at[_new_amino_acid, _padded_position]
            _frequencies = [Decimal(x) for x in df.loc[(df['padded_position'] == _padded_position) & (df['mutant_aa'] == _new_amino_acid) & (df[myoptions.column_with_frequencies] >= myoptions.threshold)][myoptions.column_with_frequencies].to_list()]
            _old_amino_acid = df.loc[df['padded_position'] == _padded_position]['original_aa'].to_list()[0]
            _old_codon = df.loc[df['padded_position'] == _padded_position]['original_codon'].to_list()[0]
            _new_codons = df.loc[(df['padded_position'] == _padded_position) & (df['mutant_aa'] == _new_amino_acid) & (df[myoptions.column_with_frequencies] >= myoptions.threshold)]['mutant_codon'].to_list()
            _new_codon = _new_codons[0]
            _observed_codon_counts = df.loc[(df['padded_position'] == _padded_position) & (df['mutant_aa'] == _new_amino_acid) & (df[myoptions.column_with_frequencies] >= myoptions.threshold)]['observed_codon_count'].to_list()
            _observed_codon_count_sum = sum(_observed_codon_counts)
            _total_codons_per_site = df.loc[(df['padded_position'] == _padded_position) & (df['mutant_codon'] == _new_codon)]['total_codons_per_site'].to_list()
            if len(_total_codons_per_site):
                _total_codons_per_site = _total_codons_per_site[0]
            _score = matrix[_old_amino_acid][_new_amino_acid]
            if _new_codon not in _new_codons:
                raise ValueError(f"The new codon {_new_codon} is not in the list of all codons {str(_new_codons)} encoding this aa {_new_amino_acid}")

            if myoptions.column_with_frequencies == 'neutralized_parent_difference':
                sel.annotation.set_text(f"Padded position: {_padded_position}\nPosition: {_position_in_protein}\nOriginal Amino Acid: {_old_amino_acid} ({_old_codon})\nNew Amino Acid: {_new_amino_acid} ({_new_codons})\n{matrix_name} score: {_score}\nCumulative Frequency: {_frequency:.6f}\nCodon Frequencies: {[f'{x:.6f}' for x in _frequencies]}")
            elif myoptions.column_with_frequencies == 'escape_parent_difference':
                sel.annotation.set_text(f"Padded position: {_padded_position}\nPosition: {_position_in_protein}\nOriginal Amino Acid: {_old_amino_acid} ({_old_codon})\nNew Amino Acid: {_new_amino_acid} ({_new_codons})\n{matrix_name} score: {_score}\nCumulative Frequency: {_frequency:.6f}\nCodon Frequencies: {[f'{x:.6f}' for x in _frequencies]}")
            elif myoptions.column_with_frequencies == 'weighted_diff_escape_neutralized':
                sel.annotation.set_text(f"Padded position: {_padded_position}\nPosition: {_position_in_protein}\nOriginal Amino Acid: {_old_amino_acid} ({_old_codon})\nNew Amino Acid: {_new_amino_acid} ({_new_codons})\n{matrix_name} score: {_score}\nCumulative Frequency: {_frequency:.6f}\nCodon Frequencies: {[f'{x:.6f}' for x in _frequencies]}")
            elif myoptions.column_with_frequencies == 'frequency':
                sel.annotation.set_text(f"Padded position: {_padded_position}\nPosition: {_position_in_protein}\nOriginal Amino Acid: {_old_amino_acid} ({_old_codon})\nNew Amino Acid: {_new_amino_acid} ({_new_codons})\n{matrix_name} score: {_score}\nCumulative Frequency: {_frequency:.6f}\nCodon Frequencies: {[f'{x:.6f}' for x in _frequencies]}\nObserved codon counts: {_observed_codon_counts}\nObserved codon count sum: {_observed_codon_count_sum}\nTotal codons per site: {_total_codons_per_site}")
            else:
                sel.annotation.set_text(f"Padded position: {_padded_position}\nPosition: {_position_in_protein}\nOriginal Amino Acid: {_old_amino_acid} ({_old_codon})\nNew Amino Acid: {_new_amino_acid} ({_new_codons})\n{matrix_name} score: {_score}\nCumulative Frequency: {_frequency:.6f}\nCodon Frequencies: {[f'{x:.6f}' for x in _frequencies]}\nObserved codon counts: {_observed_codon_counts}\nObserved codon count sum: {_observed_codon_count_sum}\nTotal codons per site: {_total_codons_per_site}")
    else:
        @_cursor.connect("add")
        def on_add(sel):
            # sel.index is the index into _mpl_scatterplot's data array, i.e.
            # into circles_matplotlib.  Do NOT use sel.target: it returns the
            # raw mouse position in data coordinates, which produces wrong
            # padded_position values when converted back via column arithmetic.
            _pt = circles_matplotlib[sel.index]
            _padded_position = _pt[0]  # stored as _padded_position when appended
            ypos = _pt[1]              # stored as i (row index in new_codon_table)
            print(f"Info: _padded_position={_padded_position}, ypos={ypos}")
            _new_codon = new_codon_table.index[ypos]
            _position_in_protein = padded_position2position[_padded_position]
            _frequency = new_codon_table.at[_new_codon, _padded_position]
            print(f"Info: _new_codon={_new_codon}, _padded_position={_padded_position}, _position_in_protein={_position_in_protein}, _frequency={_frequency}")
            _old_codon = df.loc[df['padded_position'] == _padded_position]['original_codon'].to_list()[0]
            _old_amino_acid = df.loc[df['padded_position'] == _padded_position]['original_aa'].to_list()[0]
            _observed_codon_count = df.loc[(df['padded_position'] == _padded_position) & (df['mutant_codon'] == _new_codon) & (df[myoptions.column_with_frequencies] >= myoptions.threshold)]['observed_codon_count'].to_list()[0]
            _total_codons_per_site = df.loc[(df['padded_position'] == _padded_position) & (df['mutant_codon'] == _new_codon) & (df[myoptions.column_with_frequencies] >= myoptions.threshold)]['total_codons_per_site'].to_list()[0]
            print(f"Info: {len(df.loc[df['padded_position'] == _padded_position]['position'])} aa residues observed in position {_position_in_protein}:{os.linesep} {str(df.loc[df['padded_position'] == _padded_position][0:])}{os.linesep}")
            try:
                _new_amino_acid = df.loc[(df['padded_position'] == _padded_position) & (df['mutant_codon'] == _new_codon)]['mutant_aa'].to_list()[0]
            except IndexError:
                _new_amino_acid = f"Failed to find, the {_new_codon} is wrong and too far from original codon"
            _score = matrix[_old_amino_acid][_new_amino_acid]
            try:
                _some_frequency = Decimal(df.loc[(df['padded_position'] == _padded_position) & (df['mutant_codon'] == _new_codon)][myoptions.column_with_frequencies].to_list()[0])
            except (IndexError, ValueError, TypeError):
                _some_frequency = 0.00000000009
            if _some_frequency != _frequency and (_some_frequency != 0.00000000009 and _frequency != 0) and not np.abs(_frequency) < myoptions.threshold:
                raise ValueError(f"Frequency new_codon_table.at[_new_codon, _padded_position]={_frequency} not same as df.loc[(df['padded_position'] == _padded_position) & (df['mutant_codon'] == _new_codon)][myoptions.column_with_frequencies].to_list()[0]={_frequency}")

            if myoptions.debug:
                print(f"Debug: Padded position: {_padded_position} Position: {_position_in_protein} Original Codon: {_old_codon} ({_old_amino_acid}) New Codon: {_new_codon} ({_new_amino_acid}) {matrix_name} score: {_score} Frequency: {_frequency:.6f}")

            if myoptions.column_with_frequencies == 'neutralized_parent_difference':
                sel.annotation.set_text(f"Padded position: {_padded_position}\nPosition: {_position_in_protein}\nOriginal Codon: {_old_codon} ({_old_amino_acid})\nNew Codon: {_new_codon} ({_new_amino_acid})\n{matrix_name} score: {_score}\nDifference neutralized2parent: {_frequency:.6f}")
            elif myoptions.column_with_frequencies == 'escape_parent_difference':
                sel.annotation.set_text(f"Padded position: {_padded_position}\nPosition: {_position_in_protein}\nOriginal Codon: {_old_codon} ({_old_amino_acid})\nNew Codon: {_new_codon} ({_new_amino_acid})\n{matrix_name} score: {_score}\nDifference escape2parent: {_frequency:.6f}")
            elif myoptions.column_with_frequencies == 'weighted_diff_escape_neutralized':
                sel.annotation.set_text(f"Padded position: {_padded_position}\nPosition: {_position_in_protein}\nOriginal Codon: {_old_codon} ({_old_amino_acid})\nNew Codon: {_new_codon} ({_new_amino_acid})\n{matrix_name} score: {_score}\nWeighted difference escape2neutralized: {_frequency:.6f}")
            elif myoptions.column_with_frequencies == 'frequency':
                sel.annotation.set_text(f"Padded position: {_padded_position}\nPosition: {_position_in_protein}\nOriginal Codon: {_old_codon} ({_old_amino_acid})\nNew Codon: {_new_codon} ({_new_amino_acid})\n{matrix_name} score: {_score}\nFrequency: {_frequency:.6f}\nObserved codon count: {_observed_codon_count}\nTotal codons per site: {_total_codons_per_site}")
            else:
                sel.annotation.set_text(f"Padded position: {_padded_position}\nPosition: {_position_in_protein}\nOriginal Codon: {_old_codon} ({_old_amino_acid})\nNew Codon: {_new_codon} ({_new_amino_acid})\n{matrix_name} score: {_score}\nFrequency: {_frequency:.6f}\nObserved codon count: {_observed_codon_count}\nTotal codons per site: {_total_codons_per_site}")
            if myoptions.debug:
                print(f"Debug: final_sorted_whitelist={str(final_sorted_whitelist)}")
                print(f"Debug: codons_whitelist2={str(codons_whitelist2)}")

    _handles, _labels = [], []
    if myoptions.aminoacids:
        _junk = 'NNN'
        _codon_on_input = False
    else:
        _junk = 'NNN'
        _codon_on_input = True
    for _freq in [0.001, 0.01, 0.1, 0.3]:
        if myoptions.column_with_frequencies in ['neutralized_parent_difference', 'escape_parent_difference']:
            _size, _color = adjust_size_and_color_neutralized_escape(Decimal(_freq), _codon_on_input)
        elif myoptions.column_with_frequencies in ['weighted_diff_escape_neutralized']:
            _size, _color = adjust_size_and_color_weighted(Decimal(_freq))
        else:
            _score, _size, _color = adjust_size_and_color(myoptions, Decimal(_freq), _codon_on_input, _junk, _junk, _junk, _junk, matrix, norm, colors)
        _handle = ax2.scatter(_size, - 400 + _freq, s=float(_freq * 5000), color='gray', alpha=0.5, label=f'Frequency {_freq:.1%}')
        _label = str(_freq)
        _handles.append(_handle)
        _labels.append(_label)
    ax4.set_axis_off()
    ax2.legend(loc='upper center', bbox_to_anchor=(1.25, 1.00), labelspacing=3, frameon=False, handletextpad=1.5)

    for _ext in ('.png', '.pdf'):
        _wholefig = plt.gcf()
        _figsize = _wholefig.get_size_inches()*_wholefig.dpi
        print(f"Info: Writing into {outfile_prefix + _ext}, figure size is {_wholefig.get_size_inches()} inches and {_figsize} dpi")
        figure.savefig(outfile_prefix + _ext, dpi=myoptions.dpi)
    plt.show()
    plt.clf()



# vim:ts=4:sw=4:expandtab:smartindent
