# pylint: disable=too-many-lines
# This work © 2025-2026 by Jiří Zahradník and Martin Mokrejš
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
from matplotlib import ticker

import numpy as np
import mplcursors

# https://pypi.org/project/blosum/
import blosum

# bokeh is imported lazily inside render_bokeh()

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
    """Retrieve or create a matplotlib colormap and normalization object."""

    _norm = None
    _colors = None

    _micro_cvd_gray = ["#F5F5F5", "#D6D6D6", "#B7B7B7", "#8B8B8B", "#616161"]
    _micro_cvd_green = ["#DDFFA0", "#BDEC6F", "#97CE2F", "#6D9F06", "#4E7705"]
    _micro_cvd_orange = ["#FFD5AF", "#FCB076", "#F09163", "#C17754", "#9D654C"]
    _micro_cvd_blue = ["#E7F4FF", "#BCE1FF", "#7DCCFF", "#56B4E9", "#098BD9"]
    _micro_cvd_turquoise = ["#A3E4D7", "#48C9B0", "#43BA8F", "#009E73", "#148F77"]
    _micro_cvd_purple = ["#EFB6D6", "#E794C1", "#CC79A7", "#A1527F", "#7D3560"]

    _custom_found = True
    # Handle custom hardcoded colormaps first to avoid "not found" warnings
    if colormapname == 'amino_acid_changes':
        _colors = ["#930000", "#930000", "#930000", "#930000", "#930000", "#930000", "#960000", "#580041", "#8200ff", "#c500ff", "#ff00fd", "#CC79A7", "#eea1d0", "#cc0000", "#ff0000", "#ff4f00", "#ff7c7c", "#ff9999", "#c58a24", "#9c644b", "#ffff00", "#ffcc00", "#ffa200", "#7DCCFF", "#0042ff", "#0000ff", "#D6D6D6", "#B7B7B7", "#8B8B8B", "palegreen", "#bbff00", "#97CE2F", "#219f11", "#930000", "#930000", "#930000", "#930000", "#930000", "#930000"]
        _cmap = matplotlib.colors.ListedColormap(_colors, "amino_acid_changes", len(_colors))
        myoptions.colormap = 'amino_acid_changes'
        _norm = matplotlib.colors.BoundaryNorm(np.arange(-19, 19, 1), _cmap.N)

    elif colormapname == 'dkeenan_26cols':
        _colors = ["#00B7FF", "#004DFF", "#00FFFF", "#826400", "#580041", "#FF00FF", "#00FF00", "#C500FF", "#B4FFD7", "#FFCA00", "#969600", "#B4A2FF", "#C20078", "#000000", "#0000C1", "#FF8B00", "#FFC8FF", "#666666", "#FF0000", "#CCCCCC", "#009E8F", "#D7A870", "#8200FF", "#960000", "#BBFF00", "#FFFF00", "#006F00"]
        _cmap = matplotlib.colors.ListedColormap(_colors, "dkeenan_26cols")
        myoptions.colormap = 'dkeenan_26cols'
        _norm = matplotlib.colors.BoundaryNorm(np.arange(-13, 13, 1), _cmap.N)

    elif colormapname == 'microshades_cvd_palettes':
        _colors = list(_micro_cvd_orange) + list(_micro_cvd_turquoise) + list(_micro_cvd_blue) + list(_micro_cvd_purple) + list(_micro_cvd_green) + list(_micro_cvd_gray)
        _cmap = matplotlib.colors.ListedColormap(_colors, 'microshades_cvd_palettes')

    elif colormapname == 'microshades_cvd_palettes_r':
        _colors = list(_micro_cvd_orange[::-1]) + list(_micro_cvd_turquoise[::-1]) + list(_micro_cvd_blue[::-1]) + list(_micro_cvd_purple[::-1]) + list(_micro_cvd_green[::-1]) + list(_micro_cvd_gray[::-1])
        _cmap = matplotlib.colors.ListedColormap(_colors, 'microshades_cvd_palettes_r')

    elif colormapname == 'merged':
        _colors1 = matplotlib.colormaps['gist_heat'](np.linspace(0.2, 0.8, 9))
        _colors2 = matplotlib.colormaps['cividis_r'](np.linspace(0.2, 0.8, 9))
        _colors3 = matplotlib.colormaps['GnBu_r'](np.linspace(0.2, 0.8, 9))
        _colors = np.vstack((_colors1, _colors2, _colors3))
        _cmap = matplotlib.colors.LinearSegmentedColormap.from_list('merged', _colors)

    elif colormapname == 'merged1':
        _colors1 = matplotlib.colormaps['viridis'](np.linspace(0.2, 0.8, 9))
        _colors2 = matplotlib.colormaps['coolwarm_r'](np.linspace(0., 1, 9))
        _colors3 = matplotlib.colormaps['GnBu_r'](np.linspace(0.2, 0.8, 9))
        _colors = np.vstack((_colors1, _colors2, _colors3))
        _cmap = matplotlib.colors.LinearSegmentedColormap.from_list('merged', _colors)
    else:
        _custom_found = False

    if _custom_found:
        return _norm, _cmap, _colors

    _wished_cmapname_prefix = ''
    _wished_cmapname_num = 0
    if '_' in colormapname:
        try:
            _wished_cmapname_prefix = '_'.join(colormapname.split('_')[:-1])
            _wished_cmapname_num = int(colormapname.split('_')[-1])
        except (ValueError, IndexError):
            _wished_cmapname_prefix = colormapname
            _wished_cmapname_num = 0

    try:
        _cmap = matplotlib.colormaps.get_cmap(colormapname)
    except Exception:
        try:
            _cmap = plt.get_cmap(colormapname)
        except Exception:
            try:
                import palettable
                _cmap_from_palettable = palettable.colorbrewer.get_map(_wished_cmapname_prefix, 'diverging', _wished_cmapname_num)
                _cmap = matplotlib.colors.ListedColormap(_cmap_from_palettable.mpl_colors)
            except (ImportError, Exception):
                try:
                    import palettable
                    _cmap_from_palettable = palettable.colorbrewer.get_map(_wished_cmapname_prefix, 'sequential', _wished_cmapname_num)
                    _cmap = matplotlib.colors.ListedColormap(_cmap_from_palettable.mpl_colors)
                except (ImportError, Exception):
                    try:
                        import palettable.scientific.diverging
                        _cmap_from_palettable = palettable.scientific.diverging.get_map(_wished_cmapname_prefix)
                        _cmap = matplotlib.colors.ListedColormap(_cmap_from_palettable.mpl_colors)
                    except (ImportError, Exception):
                        try:
                            import palettable.scientific.sequential
                            _cmap_from_palettable = palettable.scientific.sequential.get_map(_wished_cmapname_prefix)
                            _cmap = matplotlib.colors.ListedColormap(_cmap_from_palettable.mpl_colors)
                        except (ImportError, Exception):
                            try:
                                import palettable.mygbm
                                _cmap_from_palettable = palettable.mygbm.get_map(_wished_cmapname_prefix) # pylint: disable=no-member
                                _cmap = matplotlib.colors.ListedColormap(_cmap_from_palettable.mpl_colors)
                            except (ImportError, Exception):
                                print(f"Warning: Colormap {colormapname} not found, falling back to coolwarm_r")
                                _cmap = matplotlib.colormaps.get_cmap('coolwarm_r')

    return _norm, _cmap, _colors


def resolve_codon_or_aa(myoptions, old_codon_or_aa, new_codon_or_aa):
    """Determine the codon or amino acid labels to use based on mode."""

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


# Module-level cache for BLOSUM substitution scores.
#
# Performance rationale
# ---------------------
# blosum.BLOSUM is a subclass of collections.defaultdict and stores scores in a
# nested dict of dicts. Individual lookups (matrix[aa1][aa2]) are O(1) at the
# Python dict level, so the structure itself is not the bottleneck.
#
# The problem is that our caller loop is O(N_codons_or_aa x N_positions), and
# with standard amino acid alphabets there are at most 20x20 = 400 distinct
# (old_aa, new_aa) pairs regardless of how many positions the alignment has.
# Before this cache was added, every cell in that loop independently called
# matrix[old_aa][new_aa], meaning the same pairs were looked up thousands of
# times (e.g. 2331 calls observed during profiling of a real dataset).
#
# Additionally, when codon_on_input=True, each call also invoked alt_translate()
# twice to convert codons to amino acids before doing the matrix lookup. While
# alt_translate() is itself cached via @functools.lru_cache, two extra function
# calls and cache lookups per score request still accumulate when multiplied
# across all positions.
#
# The dictionary is keyed by (old_aa, new_aa) string tuples — always single
# uppercase amino acid letters after translation. functools.lru_cache is not
# used here because the matrix object (blosum.BLOSUM) is not hashable, so a
# module-level plain dict is the simplest correct solution.
#
# The cache is cleared at the start of each collect_scatter_data() call to
# ensure correctness when different matrices (e.g. BLOSUM62 vs BLOSUM80) are
# used in successive calls within the same Python process.
_score_cache: dict = {}


def get_score(myoptions, matrix, codon_on_input, old_codon_or_aa, new_codon_or_aa):
    """Return the integer substitution score for an amino acid (or codon) pair.

    Parameters
    ----------
    myoptions : argparse.Namespace
        CLI options; only used for error messages (myoptions.matrix name).
    matrix : blosum.BLOSUM
        The substitution scoring matrix (a nested defaultdict subclass).
        blosum.BLOSUM stores precomputed scores in plain Python dicts, so
        individual lookups are O(1) and the library has no caching of its own
        — it simply does not know whether the caller will repeat the same pair.
    codon_on_input : bool
        If True, old_codon_or_aa and new_codon_or_aa are 3-nt codons that must
        be translated to amino acids before the matrix lookup.
    old_codon_or_aa : str
        The reference codon or amino acid (before mutation).
    new_codon_or_aa : str
        The observed codon or amino acid (after mutation).

    Returns
    -------
    int
        The BLOSUM (or custom matrix) substitution score, clipped to -11 on
        arithmetic overflow.

    Performance notes
    -----------------
    Results are cached in the module-level ``_score_cache`` dict, keyed by the
    (old_aa, new_aa) pair after translation. With a standard 20-amino-acid
    alphabet there are at most 400 unique pairs, so the cache saturates quickly
    and subsequent calls for the same pair cost only a single dict lookup.
    The cache must be cleared (via ``_score_cache.clear()``) before each
    independent plot run when the matrix may change.
    """
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
            raise ValueError(f"Aieee, old_codon_or_aa='{old_codon_or_aa}' is empty")
        if not _new_codon_or_aa:
            raise ValueError(f"Aieee, new_codon_or_aa='{new_codon_or_aa}' is empty")

        if myoptions.debug:
            print(f"Info: some single-letter but neither asterisk nor dash amino acid residue: _old_codon_or_aa='{_old_codon_or_aa}', new_codon_or_aa='{new_codon_or_aa}'")

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
                sys.stdout.write(f"Info: Translating {_old_codon_or_aa} to {_new_codon_or_aa} to fetch color from _colors," )
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
        _size = weighted_diff_escape_neutralized * generic_circle_size
        _color = 'black'
    elif weighted_diff_escape_neutralized > 0.1:
        _size = weighted_diff_escape_neutralized * weighted_diff_escape_neutralized_size1
        _color = 'darkblue'
    elif weighted_diff_escape_neutralized > 0.01:
        _size = weighted_diff_escape_neutralized * weighted_diff_escape_neutralized_size2
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
            _matrix_name = f"BLOSUM{_matrix_num}"
        else:
            sys.stderr.write(f"Warning: Unexpected matrix type {str(_matrix_type)}, falling back to BLOSUM\n")
            _matrix = blosum.BLOSUM(_matrix_num)
            _matrix_name = f"BLOSUM{_matrix_num}"
            myoptions.matrix = _matrix_name

    if not myoptions.outfile_prefix:
        raise RuntimeError("Please provide output filename prefix via --outfile-prefix")

    # Clean up the output prefix by stripping common extensions if user provided a full filename
    _clean_prefix = myoptions.outfile_prefix
    while True:
        _changed = False
        for _ext in ['.frequencies.tsv', '.frequencies', '.tsv', '.png', '.jpg', '.pdf', '.html']:
            if _clean_prefix.endswith(_ext):
                _clean_prefix = _clean_prefix[:-len(_ext)]
                _changed = True
                break
        if not _changed:
            break

    _outfile_prefix = _clean_prefix + '.' + _matrix_name + '.' + myoptions.colormap
    print(f"Info: _outfile_prefix={_outfile_prefix}")

    _theoretical_scores = set()
    for _scores_dict in _matrix.values():
        for _score in _scores_dict.values():
            _theoretical_scores.add(_score)
    print(f"Info: Using {_matrix_name} matrix now. Theoretical minimum score is {min(_theoretical_scores)}, theoretical maximum score is {max(_theoretical_scores)}, values are {str(_theoretical_scores)}")
    _min_theoretical_score, _max_theoretical_score = int(min(_theoretical_scores)), int(max(_theoretical_scores))

    return _matrix, _matrix_name, _min_theoretical_score, _max_theoretical_score, _outfile_prefix


def load_and_clean_dataframe(myoptions, infilename, padded_position2position):
    """
    Load the input TSV, normalise column names for legacy formats, and filter noisy rows.

    Optimization:
    - Speedup 6: Uses Decimal(str) for all frequency parsing to maintain
      high precision and bit-identical output across platforms.

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

    print(f"Info: Parsing input file {infilename}")
    if not infilename:
        raise RuntimeError("Please provide an input TSV file")
    if not os.path.exists(infilename):
        raise RuntimeError(f"Input file not found: {infilename}")
    if os.path.getsize(infilename) == 0:
        raise RuntimeError(f"Input file is empty: {infilename}")

    df = pd.read_csv(infilename, sep='\t', header='infer', na_filter=False, na_values=[None])

    _freq_col = myoptions.column_with_frequencies
    if _freq_col in df.columns:
        # High precision: Convert to Decimal from string to avoid float64 drift
        df[_freq_col] = df[_freq_col].astype(str).apply(Decimal)

    if 'position' in df.columns:
        print(f"Info: Autodetected new TSV file format with a header in {infilename}")
    else:
        print(f"Info: Autodetected old TSV file format without a header in {infilename}, assigning default column names")
        df = pd.read_csv(infilename, sep='\t', header=None, na_filter=False, na_values=[None])
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
            raise RuntimeError(f"Unexpected number of columns in the {infilename} file")

    print(f"Info: The file {infilename} contains now these columns: {str(df.columns)}")

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
        raise ValueError(f"Cannot parse column df['mutant_codon']={_mutant_codon} containing {len(_mutant_codon)} values") from _exc

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
            print(f"Info: Originally there were {_before} rows but after discarding codons with [N n] there are only {_after} left")
        else:
            print(f"Info: Originally there were {_before} rows but after discarding codons with [N n DEL] there are only {_after} left")
    elif myoptions.showdel:
        print(f"Info: Originally there were {_before} rows but after discarding codons with [N n *] there are only {_after} left")
    else:
        print(f"Info: Originally there were {_before} rows but after discarding codons with [N n DEL] there are only {_after} left")
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
        raise ValueError(f"Length of _codons_whitelist is {len(_codons_whitelist)} which is not equal to _codons_whitelist_aa with {len(_codons_whitelist_aa)}")
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
    # Pass 1: Filter and group by mutations
    # We use float casting for the threshold check to avoid Decimal vs float comparison overhead
    _freq_col = myoptions.column_with_frequencies
    # Condition 1: Keep row if (not aminoacid mode) OR (AA changed) OR (including synonymous)
    if myoptions.include_synonymous or not myoptions.aminoacids:
        _cond_mutation = True
    else:
        _cond_mutation = df['original_aa'] != df['mutant_aa']

    # Condition 2: Frequencies above threshold
    _cond_threshold = df[_freq_col].astype(float).abs() >= myoptions.threshold

    _mask = _cond_mutation & _cond_threshold
    _filtered_df = df.loc[_mask]

    if myoptions.debug:
        print(f"Debug: build_frequency_tables processing {len(_filtered_df)} valid rows out of {len(df)}")

    # Pass 2: Vectorized aggregation
    # Grouping by (mutant/original, position) and summing frequencies.
    # We use these intermediate series to update our pre-allocated Decimal tables.
    _new_aa_sums = _filtered_df.groupby(['mutant_aa', 'padded_position'])[_freq_col].sum()
    _old_aa_sums = _filtered_df.groupby(['original_aa', 'padded_position'])[_freq_col].sum()
    _new_codon_sums = _filtered_df.groupby(['mutant_codon', 'padded_position'])[_freq_col].sum()
    _old_codon_sums = _filtered_df.groupby(['original_codon', 'padded_position'])[_freq_col].sum()

    # Pass 3: Update pre-allocated Decimal tables
    # While we still iterate over unique (aa/codon, pos) pairs, this is far fewer visits
    # than the original O(N_rows) iterrows() loop.
    for (_aa, _pos), _val in _new_aa_sums.items():
        if _aa in _new_aa_table.index and _pos in _new_aa_table.columns:
            _new_aa_table.at[_aa, _pos] = Decimal(_val)

    for (_aa, _pos), _val in _old_aa_sums.items():
        if _aa in _old_aa_table.index and _pos in _old_aa_table.columns:
            _old_aa_table.at[_aa, _pos] = Decimal(_val)

    for (_codon, _pos), _val in _new_codon_sums.items():
        _codon = _codon.upper()
        if _codon in _new_codon_table.index and _pos in _new_codon_table.columns:
            _new_codon_table.at[_codon, _pos] = Decimal(_val)

    for (_codon, _pos), _val in _old_codon_sums.items():
        _codon = _codon.upper()
        if _codon in _old_codon_table.index and _pos in _old_codon_table.columns:
            _old_codon_table.at[_codon, _pos] = Decimal(_val)

    # Finalize offset and debug info
    if not _filtered_df.empty:
        _first_row = _filtered_df.iloc[0]
        _very_leftmost_aa_pos = int(_first_row['padded_position'])
        _calculated_aa_offset = _very_leftmost_aa_pos - myoptions.offset + 1
        if myoptions.debug:
            print(f"Debug: calculated offset is {_calculated_aa_offset}")

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
    new_aa_table, new_codon_table,
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

        _ax2.set_ylabel(f'Cumulative frequency of mutations above threshold {myoptions.threshold:.2%} per codon', fontsize=12)
        _ax2.yaxis.set_major_formatter(ticker.PercentFormatter(1.0))
        _ax2.yaxis.set_major_locator(ticker.MultipleLocator(0.1))
        _ax1.figure.canvas.draw()
        _ax2.figure.canvas.draw()

        if myoptions.aminoacids:
            _ax2.bar(unique_aa_padded_positions, _total_frequencies, color='black', alpha=0.5, width=0.8, align='center')
        else:
            _ax2.bar(unique_padded_codon_positions, _total_frequencies, color='black', alpha=0.5, width=0.8, align='center')


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
    """
    Iterate over frequency tables and collect scatter plot data, labels, and colors.
    Synonymous mutations will be represented like 'D1146D' or 'I68I'.

    Optimization:
    - Speedup 3: Implements an O(1) hover text pipeline for Bokeh by pre-formatting
      metadata in ColumnDataSource, avoiding heavy per-dot calculation during rendering.
    """

    _real_aa_positions = sorted(padded_position2position.values()) # better extract it dynamically from the table then from the dictionary just in case some values would be discarded while parsing on-the-fly
    if myoptions.debug:
        print(f"Debug: _real_aa_positions={str(_real_aa_positions)}")

    if not myoptions.aminoacids:
        if list(table.index) != codons_whitelist2:
            raise ValueError(f"Both lists should be equal: table.index={str(table.index)}, codons_whitelist2={codons_whitelist2}")
    else:
        if list(table.index) != amino_acids:
            raise ValueError(f"Both lists should be equal: table.index={str(table.index)}, amino_acids={amino_acids}")

    _used_colors = set()
    _norm, _cmap, _colors = get_colormap(myoptions, myoptions.colormap)
    _mutations: list[str] = []
    _circles_bokeh: list[tuple[typing.Any, ...]] = []
    _circles_matplotlib: list[tuple[typing.Any, ...]] = []
    _markers: list[tuple[typing.Any, ...]] = []
    _dots: list[tuple[typing.Any, ...]] = []
    _warn_once: list[int] = []
    _matrix_values: set[int] = set()
    _hover_text_bokeh: list[str] = []

    if myoptions.aminoacids:
        _colors_tsv_filename = outfile_prefix + '.aa.frequencies.colors.tsv'
    else:
        _colors_tsv_filename = outfile_prefix + '.codon.frequencies.colors.tsv'

    _colors_tsv_rows = []

    # the tables were constructed with the following in build_frequency_tables()
    #     _new_aa_table = pd.DataFrame(Decimal(0), index=_amino_acids, columns=_unique_padded_aa_positions)
    #     _new_codon_table = pd.DataFrame(Decimal(0), index=_codons_whitelist2, columns=_unique_padded_codon_positions)
    _pos_to_old_codon = df.groupby('padded_position')['original_codon'].first().to_dict()
    _pos_to_old_aa = df.groupby('padded_position')['original_aa'].first().to_dict()
    _mut_col = 'mutant_aa' if myoptions.aminoacids else 'mutant_codon'
    _df_indexed = df.set_index(['padded_position', _mut_col])
    # Pre-build a plain dict mapping (padded_position, codon_or_aa) -> DataFrame
    # for O(1) per-cell lookup instead of repeated pandas .loc[] calls.
    #
    # Performance rationale
    # ---------------------
    # The main loop below is O(N_rows x N_positions) where N_rows is the
    # number of distinct codon or amino acid rows in the frequency table.
    # Profiling showed 2676 calls to pandas DataFrame.loc/__getitem__ during
    # a single plot run (test3.default, aminoacids mode), each of which
    # traverses the pandas index machinery even though the result is always
    # a single row or an empty frame.
    #
    # pandas.DataFrame.loc[] is designed for flexible, label-based access
    # with type coercion and index alignment, all of which add overhead on
    # every call. Using a pre-built plain Python dict eliminates this:
    # _df_groups.get(key, pd.DataFrame()) is a single hash lookup with no
    # pandas indexing overhead, and construction of the dict (one pass
    # through groupby) is O(N_rows_in_df), done once before the main loop.
    _df_groups: dict = {}
    for _key, _sub in _df_indexed.groupby(level=[0, 1]):
        _df_groups[_key] = _sub

    # Clear the score cache so each plot run starts fresh (matrix may differ)
    _score_cache.clear()
    _matrix_values = set()
    matrix_name = myoptions.matrix

    for i, _some_codon_or_aa in enumerate(table.index): # so _some_codon_or_aa contains the index specified when the table was constructed
        for j, _padded_position in enumerate(table.columns): # so _aa_position contains the real aa_position
            if not xmin <= _padded_position <= xmax:
                continue
            if myoptions.debug:
                print(f"Debug: i: {i}, j: {j}, _padded_position column: {_padded_position}")
            try:
                _aa_position = padded_position2position[_padded_position]
            except KeyError:
                continue
            if myoptions.debug:
                print(f"Debug0: _padded_positions (typically will not be contiguous and will contain multiplicates): {sorted(list(table.columns))}{os.linesep}")
                print(f"Debug0:     _aa_positions (typically will not be contiguous and will contain multiplicates): {sorted(padded_position2position.values())}{os.linesep}")
            _frequency = table.loc[_some_codon_or_aa, _padded_position]
            if myoptions.debug and _frequency:
                print(f"Debug0: _padded_position={_padded_position}, _aa_position={_aa_position}, _some_codon_or_aa={_some_codon_or_aa}, _frequency={_frequency}")
            try:
                _old_codon = _pos_to_old_codon[_padded_position]
            except KeyError:
                if _frequency and myoptions.debug:
                    print(f"Debug0b: _padded_position={_padded_position}, _some_codon_or_aa={_some_codon_or_aa}, _frequency={_frequency}")
                if _padded_position not in _warn_once:
                    sys.stderr.write(f"Warning: Cannot determine original codon for position {_padded_position}, seems missing from input TSV{os.linesep}")
                    _warn_once.append(_padded_position)
                continue
            _codon_on_input, _old_codon_or_aa, _new_codon_or_aa = resolve_codon_or_aa(myoptions, _old_codon, _some_codon_or_aa)
            if myoptions.aminoacids and not myoptions.include_synonymous and _old_codon_or_aa == _new_codon_or_aa:
                continue

            # Data lookup and variable initialization for metadata/reporting
            _old_amino_acid = _pos_to_old_aa.get(_padded_position, "Unknown")
            _base_df = _df_groups.get((_padded_position, _some_codon_or_aa), pd.DataFrame())
            _sub_df = _base_df[_base_df[myoptions.column_with_frequencies].astype(float).abs() >= myoptions.threshold] if not _base_df.empty else pd.DataFrame()

            if myoptions.aminoacids:
                _new_amino_acid = _some_codon_or_aa
                _new_codons = _sub_df['mutant_codon'].to_list() if not _sub_df.empty and 'mutant_codon' in _sub_df.columns else []
            else:
                _new_amino_acid = _sub_df['mutant_aa'].iloc[0] if not _sub_df.empty and 'mutant_aa' in _sub_df.columns else "N/A"
                if _new_amino_acid == "N/A":
                    _new_amino_acid = alt_translate(_some_codon_or_aa)
                _new_codons = [_some_codon_or_aa]

            _score = None
            if not abs(Decimal(_frequency)) < myoptions.threshold:
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
                _bokeh_size = float(np.sqrt(abs(_size)) * 100) if myoptions.bokeh_sqrt_size else float(abs(_size) * 100)

                # --- O(1) Hover Metadata Reconstruction (Matplotlib) ---
                if myoptions.aminoacids:
                    _frequencies = [Decimal(x) for x in _sub_df[myoptions.column_with_frequencies].to_list()] if not _sub_df.empty else []
                    _observed_codon_counts = _sub_df['observed_codon_count'].to_list() if not _sub_df.empty and 'observed_codon_count' in _sub_df.columns else []
                    _total_codons_per_site = _sub_df['total_codons_per_site'].iloc[0] if not _sub_df.empty and 'total_codons_per_site' in _sub_df.columns else 0
                    _observed_codon_count_sum = sum(_observed_codon_counts)

                    try:
                        _hover_score = matrix[_old_amino_acid][_new_codon_or_aa]
                        _score_str = f"{float(_hover_score):.1f}" if _hover_score != float('-inf') else "-inf"
                    except (KeyError, TypeError):
                        _score_str = f"{float(_score):.1f}"
                    _hover_text = (
                        f"Padded position: {_padded_position}\n"
                        f"Position: {_aa_position}\n"
                        f"Original Amino Acid: {_old_amino_acid} ({_old_codon})\n"
                        f"New Amino Acid: {_some_codon_or_aa} ({_new_codons})\n"
                        f"{matrix_name} score: {_score_str}\n"
                    )
                    if myoptions.column_with_frequencies == 'neutralized_parent_difference':
                        _hover_text += f"Cumulative difference neutralized2parent: {sum(_frequencies):.6f}\n"
                    elif myoptions.column_with_frequencies == 'escape_parent_difference':
                        _hover_text += f"Cumulative difference escape2parent: {sum(_frequencies):.6f}\n"
                    elif myoptions.column_with_frequencies == 'weighted_diff_escape_neutralized':
                        _hover_text += f"Cumulative weighted difference escape2neutralized: {sum(_frequencies):.6f}\n"
                    else:
                        _hover_text += f"Cumulative Frequency: {sum(_frequencies):.6f}\n"

                    if 'observed_codon_count' in df.columns.values:
                        _hover_text += (
                            f"Codon Frequencies: {[f'{x:.6f}' for x in _frequencies]}\n"
                            f"Observed codon counts: {_observed_codon_counts}\n"
                            f"Observed codon count sum: {_observed_codon_count_sum}\n"
                            f"Total codons per site: {_total_codons_per_site}"
                        )
                else:
                    _observed_codon_count = _sub_df['observed_codon_count'].iloc[0] if not _sub_df.empty and 'observed_codon_count' in _sub_df.columns else 0
                    _total_codons_per_site = _sub_df['total_codons_per_site'].iloc[0] if not _sub_df.empty and 'total_codons_per_site' in _sub_df.columns else 0
                    if _new_amino_acid == "N/A":
                        _new_amino_acid = alt_translate(_some_codon_or_aa)

                    try:
                        _hover_score = matrix[_old_amino_acid][_new_amino_acid]
                        _score_str = f"{float(_hover_score):.1f}" if _hover_score != float('-inf') else "-inf"
                    except (KeyError, TypeError):
                        _score_str = f"{float(_score):.1f}"

                    _hover_text = (
                        f"Padded position: {_padded_position}\n"
                        f"Position: {_aa_position}\n"
                        f"Original Codon: {_old_codon} ({_old_amino_acid})\n"
                        f"New Codon: {_some_codon_or_aa} ({_new_amino_acid})\n"
                        f"{matrix_name} score: {_score_str}\n"
                    )
                    if myoptions.column_with_frequencies == 'neutralized_parent_difference':
                        _hover_text += f"Difference neutralized2parent: {_frequency:.6f}"
                    elif myoptions.column_with_frequencies == 'escape_parent_difference':
                        _hover_text += f"Difference escape2parent: {_frequency:.6f}"
                    elif myoptions.column_with_frequencies == 'weighted_diff_escape_neutralized':
                        _hover_text += f"Weighted difference escape2neutralized: {_frequency:.6f}"
                    else:
                        _hover_text += f"Frequency: {_frequency:.6f}\n"
                        _hover_text += (
                            f"Observed codon count: {_observed_codon_count}\n"
                            f"Total codons per site: {_total_codons_per_site}"
                        )

                # --- End Hover Reconstruction ---

                if myoptions.aminoacids:
                    _mutation_str = f"{_old_amino_acid}{_aa_position}{_some_codon_or_aa}"
                    if _score < 0:
                        _circles_bokeh.append((_padded_position, _some_codon_or_aa, _bokeh_size, 'circle', _color, 0.5, _score, _aa_position, _padded_position))
                    else:
                        _circles_bokeh.append((_padded_position, _some_codon_or_aa, _bokeh_size, 'hex', _color, 0.5, _score, _aa_position, _padded_position))
                else:
                    _mutation_str = f"{_old_codon}{_aa_position}{_some_codon_or_aa}"
                    _y_label_bokeh = _some_codon_or_aa + ' (' + alt_translate(_some_codon_or_aa) + ')'
                    if _score < 0:
                        _circles_bokeh.append((_padded_position, _y_label_bokeh, _bokeh_size, 'circle_x', _color, 0.5, _score, _aa_position, _padded_position))
                    else:
                        _circles_bokeh.append((_padded_position, _y_label_bokeh, _bokeh_size, 'hex', _color, 0.5, _score, _aa_position, _padded_position))

                _mutations.append(_mutation_str)
                _hover_text_bokeh.append(_hover_text)

                if _score < 0:
                    _circles_matplotlib.append((_padded_position, i, float(np.abs(_size) * 5000), 'circle_x', _color, 0.5, _score, _aa_position, _padded_position, _hover_text))
                    _markers.append((_padded_position, i, 1, 'dot', 'black', 0.5))
                else:
                    _circles_matplotlib.append((_padded_position, i, float(np.abs(_size) * 5000), 'circle', _color, 0.5, _score, _aa_position, _padded_position, _hover_text))
                    _markers.append((_padded_position, i, 1, 'circle', 'black', 0.5))
                _used_colors.add(_color)

                # Record for .colors.tsv
                if _old_amino_acid and _size:
                    _mutant_codons_list = _new_codons if myoptions.aminoacids else [_some_codon_or_aa]
                    if len(_mutant_codons_list) > 1:
                        _colors_tsv_rows.append((_padded_position, _aa_position, _old_codon, str(_mutant_codons_list), _old_amino_acid, _new_amino_acid, _frequency, _color, _score))
                    elif _mutant_codons_list:
                        _colors_tsv_rows.append((_padded_position, _aa_position, _old_codon, _mutant_codons_list[0], _old_amino_acid, _new_amino_acid, _frequency, _color, _score))
            else:
                # just draw some tiny dot otherwise pandas will drop empty Y-rows for unused amino acids or codons,
                # which sucks and it btw does happen for charts with codons too although they have more data and
                # supposedly are less likely to run into this issue but it does happen too
                _size, _color = 0.00000000009, 'black'
                _score = get_score(myoptions, matrix, _codon_on_input, _old_codon_or_aa, _new_codon_or_aa)
                _dots.append((_padded_position, i, _size, 'dot', _color, 0.5, _score))

                if myoptions.debug:
                    print(f"Debug: Invisible dot. Real AA position: {_padded_position}, observed codon: {_some_codon_or_aa}, _frequency: {_frequency}, _size: {_size}, color: {_color}")

                _used_colors.add(_color)

    with open(_colors_tsv_filename, 'w', encoding="utf-8") as _color_file:
        print(f"Info: Writing into {_colors_tsv_filename}")
        _color_file.write(f"padded_position\tposition\toriginal_codon\tmutant_codon\toriginal_aa\tmutant_aa\tfrequency\tcolor\tscore{os.linesep}")
        # Sort by padded_position (x[0]), then position (x[1]), then mutant_codon (x[3])
        for _row in sorted(_colors_tsv_rows, key=lambda x: (x[0], x[1], x[3])):
            _color_file.write(f"{_row[0]}\t{_row[1]}\t{_row[2]}\t{_row[3]}\t{_row[4]}\t{_row[5]}\t{_row[6]:.6f}\t{_row[7]}\t{_row[8]}{os.linesep}")

    if _matrix_values:
        print(f"Info: The following values were collected from matrix {myoptions.matrix} based on the actual data (some values from matrix might not be needed for your data, hence are not listed here): {str(sorted(_matrix_values))} . Range spans {abs(min(_matrix_values)) + 1 + max(_matrix_values)} values (before symmetrization).")
    else:
        print(f"Info: No mutations were found in the specified range (xmin={myoptions.xmin}, xmax={myoptions.xmax}) or with threshold {myoptions.threshold}.")
        _norm = matplotlib.colors.Normalize(vmin=-5, vmax=5)
        _cmap = matplotlib.colormaps.get_cmap('coolwarm_r')
        _used_colors = []
        _mutations = []
        return _norm, _cmap, _colors, _used_colors, _matrix_values, _mutations, _circles_bokeh, _circles_matplotlib, _markers, _dots, _hover_text_bokeh
    if myoptions.debug:
        print(f"Debug: {len(_used_colors)} _used_colors used: {str(_used_colors)}")

    if myoptions.debug:
        print(f"Debug: {len(_circles_bokeh)} _circles_bokeh={str(_circles_bokeh)}")
        print(f"Debug: {len(_markers)} _markers={str(_markers)}")
        print(f"Debug: {len(_dots)} _dots={str(_dots)}")

    return (
        _norm, _cmap, _colors, _used_colors, _matrix_values,
        _mutations,
        _circles_bokeh, _circles_matplotlib, _markers, _dots,
        _hover_text_bokeh,
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
    circles_bokeh, mutations, hover_texts,
    title_data, xlabel,
    matrix_name, colors, norm, cmap,
    show=True,
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
    mutations : list[str]
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
    import bokeh.plotting
    import bokeh.models

    if circles_bokeh:
        _circles_x, _circles_y, _circles_size, _circles_marker, _circles_color, _circles_alpha, _circles_score, _circles_aa_pos, _circles_padded_pos = zip(*circles_bokeh)
    else:
        _circles_x, _circles_y, _circles_size, _circles_marker, _circles_color, _circles_alpha, _circles_score, _circles_aa_pos, _circles_padded_pos = [], [], [], [], [], [], [], [], []

    _mysource = bokeh.models.ColumnDataSource(data={
        "x": _circles_x,
        "y": _circles_y,
        "s": _circles_size,
        "m": _circles_marker,
        "c": _circles_color,
        "a": _circles_alpha,
        "score": _circles_score,
        "aaposition": _circles_aa_pos,
        "hover_text": [x.replace('\n', '<br>') for x in hover_texts],
        "padded_pos": _circles_padded_pos,
        "mutation": mutations,
    })


    if myoptions.aminoacids:
        _tooltips = "@hover_text{safe}"
    else:
        _tooltips = "@hover_text{safe}"
    if myoptions.aminoacids:
        _p = bokeh.plotting.figure(x_range=(xmin, xmax), y_range=amino_acids, tooltips=_tooltips, title=title_data, x_axis_label=xlabel, y_axis_label='Introduced amino acid changes', width=2000, height=1200, sizing_mode='stretch_width')
    else:
        _p = bokeh.plotting.figure(x_range=(xmin, xmax), y_range=[_pairs[0] + ' (' + _pairs[1] + ')' for _pairs in final_sorted_whitelist], tooltips=_tooltips, title=title_data, x_axis_label=xlabel, y_axis_label='Introduced codon changes', height=1200, width=2000, sizing_mode='stretch_width')

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
        title=f"{matrix_name} score values (for synonymous changes forcibly set to +12 (dark green))",
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
    if show:
        bokeh.plotting.show(_p)
    else:
        bokeh.plotting.save(_p)
    pretty_print_bokeh_html(outfile_prefix + '.html')


def render_matplotlib(
    myoptions,
    figure, ax1, ax2, ax3, ax4, outfile_prefix,
    circles_matplotlib, markers, dots, cmap, norm, colors,
    matrix, matrix_name,
    show=True,
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
        cm_x, cm_y, cm_s, _, _, _, cm_c, _, _, _ = zip(*circles_matplotlib)
        _mpl_scatterplot = ax1.scatter(cm_x, cm_y, marker='o', s=cm_s, alpha=0.5, c=cm_c, cmap=cmap, norm=norm)
    else:
        _mpl_scatterplot = ax1.scatter([], [], marker='o', s=[], alpha=0.5, c=[], cmap=cmap, norm=norm)

    _colorbar = figure.colorbar(_mpl_scatterplot, cax=ax3, label=f"{matrix_name} score values (for synonymous changes forcibly set to +12 (dark green))", location='right', pad=-0.1, alpha=0.5)
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

    def on_add(sel):
        # sel.index is the index into _mpl_scatterplot's data array, i.e.
        # into circles_matplotlib.  Do NOT use sel.target: it returns the
        # raw mouse position in data coordinates, which produces wrong
        # padded_position values when converted back via column arithmetic.
        # Speedup: circles_matplotlib now contains the pre-formatted hover text
        # in the last element (index 9), providing O(1) lookup speed.
        try:
            _pt = circles_matplotlib[sel.index]
            _hover_text = _pt[9]
            sel.annotation.set_text(_hover_text)
        except (IndexError, TypeError):
            sel.annotation.set_text("Error: Hover metadata not found")

    try:
        _cursor = mplcursors.cursor(_mpl_scatterplot, hover=True)
        _cursor.connect("add", on_add)
    except (ImportError, Exception): # pylint: disable=broad-exception-caught
        pass

    _handles, _labels = [], []
    if myoptions.aminoacids:
        _junk = 'NNN'
        _codon_on_input = False
    else:
        _junk = 'NNN'
        _codon_on_input = True
    for _freq in [1.0, 0.5, 0.3, 0.1, 0.01, 0.001]:
        if myoptions.column_with_frequencies in ['neutralized_parent_difference', 'escape_parent_difference']:
            _size, _color = adjust_size_and_color_neutralized_escape(Decimal(_freq), _codon_on_input)
        elif myoptions.column_with_frequencies in ['weighted_diff_escape_neutralized']:
            _size, _color = adjust_size_and_color_weighted(Decimal(_freq))
        else:
            _score, _size, _color = adjust_size_and_color(myoptions, Decimal(_freq), _codon_on_input, _junk, _junk, _junk, _junk, matrix, norm, colors)
        _handle = ax2.scatter(0, - 400 + _freq, s=float(_freq * 5000), color='gray', alpha=0.5, label=f'{_freq:.1%}') # Frequency
        _label = str(_freq)
        _handles.append(_handle)
        _labels.append(_label)
    ax4.set_axis_off()
    ax2.legend(loc='upper center', bbox_to_anchor=(1.25, 1.00), labelspacing=5.4, frameon=False, handletextpad=2.8)

    if myoptions.debug or os.environ.get('PYTEST_CURRENT_TEST'):
        _mpl_hovers = []
        class MockAnnotation: # pylint: disable=too-few-public-methods
            """Internal mock of an annotation object for hover generation."""
            def __init__(self):
                """Initialize with empty text."""
                self.text = ""
            def set_text(self, t):
                """Set the annotation text."""
                self.text = t
        class MockSel: # pylint: disable=too-few-public-methods
            """Internal mock of a selection object for hover generation."""
            def __init__(self, idx):
                """Initialize with a specific index and mock annotation."""
                self.index = idx
                self.annotation = MockAnnotation()

        for i in range(len(circles_matplotlib)):
            mock_sel = MockSel(i)
            try:
                on_add(mock_sel)
                # Only include entries that are clearly described mutations (not "N/A", "Valid Unchanged Site", or "ProgrammingError")
                if all(x not in mock_sel.annotation.text for x in ["(N/A)", "N/A", "Valid Unchanged Site", "ProgrammingError"]):
                    _mpl_hovers.append({'index': i, 'text': mock_sel.annotation.text})
            except Exception as e:
                # Skip errors as well, as they indicate broken/incomplete hover data
                if myoptions.debug:
                    print(f"Debug: Skipping hover index {i} due to error: {str(e)}")

        with open(outfile_prefix + ".matplotlib_hovers.json", "w", encoding="utf-8") as f:
            json.dump(_mpl_hovers, f, indent=2)

    for _ext in ('.png', '.pdf'):
        _wholefig = plt.gcf()
        _figsize = _wholefig.get_size_inches()*_wholefig.dpi
        print(f"Info: Writing into {outfile_prefix + _ext}, figure size is {_wholefig.get_size_inches()} inches and {_figsize} dpi")
        figure.savefig(outfile_prefix + _ext, dpi=myoptions.dpi)
    if show:
        plt.show()
    plt.close(figure)



# vim:ts=4:sw=4:expandtab:smartindent
