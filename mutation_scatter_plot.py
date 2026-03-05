#! /usr/bin/env python3

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

Red to Blue color palettes (practically 3 to 2 colors discernible):
cool_r (practically almost 3 colors discernible) > cet_CET_C4 > coolwarm_r (practically 2 to 3 colors discernible) > cet_diverging_bwr_20_95_c54_r > seismic_r > cet_CET_D1A_r

Red to Green but diverging around zero:
Carbone42_r = SCook18_r (practically 5 colors discernible) > Carbone17_r > cet_rainbow_bgyr_10_90_c83_r ~= cet_rainbow4_r ~= cet_CET_R4_r > cet_rainbow_r

Red via narrow Green to Blue (practically 3 colors discernible):
LangRainbow12_r >> cmc.roma > cet_isoluminant_cgo_70_c39_r > cet_isoluminant_cgo_80_c38_r ~ cet_isolum_r

Red via Orange to Green:
Wild25_r (practically 4 to 5 colors discernible) > turbo_r > cet_rainbow_bgyr_35_85_c72_r > cet_rainbow_bgyr_35_85_c73_r


Blue to Orange (practically 3 colors discernible):
cet_rainbow_bgyrm_35_85_c71 ~ cet_rainbow_bgyr_35_85_c72 ~= cet_rainbow_bgyr_35_85_c73 ~= cet_rainbow

Violet to Green (practically 3 colors discernible):
cet_gouldian ~ cet_CET_L20 > cet_linear_kbgoy_20_95_c57 > cet_linear_kbgyw_10_98_c63 > cet_linear_bgyw_15_100_c68 > cet_linear_bgyw_15_100_c67 > cet_kbgyw > cet_linear_bgy_10_95_c74 > cet_CET_L9 ~ cet_CET_L16

Blue to Green (practically 3 to 2 colors discernible):
cet_linear_bgyw_20_98_c66 > cet_bgy > cmc.bukavu (is a bit yellow or brown) > cet_isoluminant_cgo_70_c39 > cet_isoluminant_cgo_80_c38 > cmc.bamako (practically 2 colors discernible)


Color Vision Deficiency:

Somewhat warm Violet to Blue (practically 3 to 4 colors discernible):
cet_linear_tritanopic_krjcw_5_95_c24

Somewhat colder Violet to Blue BUT values around zero have too intense color, almost like when at the right end:
cet_linear_tritanopic_kcw_5_95_c22 > cet_linear_protanopic_deuteranopic_kbw_5_98_c40 > cet_linear_protanopic_deuteranopic_kbw_5_95_c34


Somewhat warm Blue to Orange, better than cet_linear_protanopic_deuteranopic_kbjyw_5_95_c25 and cet_linear_protanopic_deuteranopic_kyw_5_95_c49 but suffers similar problem around zero value:
cmc.batlow

Somewhat warm Violet to Brown BUT values around zero have too intense color, almost like when at the right end:
cet_linear_protanopic_deuteranopic_kbjyw_5_95_c25 > cet_linear_protanopic_deuteranopic_kyw_5_95_c49

Blue to Brown (practically 2 colors discernible):
cet_diverging_linear_bjy_30_90_c45 > cet_diverging_linear_protanopic_deuteranopic_bjy_57_89_c34

Blue to Yellow (practically 3 colors discernible):
cet_CET_R4 > HomeyerRainbow

Single-color (practically 2 colors discernible):
cet_linear_blue_5_95_c73

Green (practically 2 colors discernible):
cet_linear_kgy_5_95_c69 ~ cet_linear_green_5_95_c69




Getting out the most from the figures

If one is after discrimination of different colorings in the figure, then several
colormaps performed well and allowed 5 or even 7 colors to be identified from the
datapoints in our suboptimal testcases (having probably just 7 different data values
with narrow spacing while rest of the space unoccupied).

cet_glasbey_Hv (practically 7 colors discernible) > cet_glasbey_dark_r (practically 6 colors discernible, somewhat light) ~ cet_glasbey_bw_minc_20_maxl_70_r (practically 6 colors discernible) > cet_glasbey_bw_minc_20_r (practically 5 colors discernible)


In default (codon) mode pale green is used to show synonymous changes.

Unfortunately, colors for values around zero are somewhat light blue and
they might appear to some readers as almost gray. Moreover, for some amino
acid changes defined in the BLOSUM62 matrix, the worst (drastic change)
value is e.g. -2 for T change to either of G, H, F, W or Y while for some
other it is -4 (W change to N, D, F or V).
Therefore, one cannot directly compare intensity of the (red) color of
multiple changes because their lowest possible minimum is just different
for each original amino acid. Still we find the color-coding of changes
helpful though not ideal.
The program mutation_scatter_plot.py developed to render the figures can
display either amino acid changes or codon changes in a single run, no both.

To omit less abundant changes from rendering there is a minimal threshold
option available.


Please cite our article if you use our data or software in your work:

Shoshany A., Tian R., Padilla-Blanco M., Hruška A., Baxova K., Zoler E., Mokrejš M., Schreiber G., Zahradník J. (submitted) In Vitro and Viral Evolution Convergence Reveal the Selective Pressures Driving Omicron Emergence. [bioRxiv](https://www.biorxiv.org/content/10.1101/2025.04.23.650148v1)


Other approaches

A completely different approach to display a few mutations is taken by:

Mutplot: An easy-to-use online tool for plotting complex mutation data with flexibility
https://doi.org/10.1371%2Fjournal.pone.0215838

Lollipops in the Clinic: Information Dense Mutation Plots for Precision Medicine
https://doi.org/10.1371%2Fjournal.pone.0160519
https://github.com/pbnjay/lollipops

Mutations Needle Plot
https://github.com/bbglab/muts-needle-plot
DOI:10.5281/zenodo.14561

Plot protein: visualization of mutations
https://doi.org/10.1186%2F2043-9113-3-14
https://sites.google.com/site/plotprotein/home
"""

import os
import sys
import re

from decimal import *
from optparse import OptionParser
import pandas as pd

# https://jiffyclub.github.io/palettable/
import palettable

# https://github.com/akvas/nclcmaps
# is not in conda-forge nor PyPi, you have to install it manually from the source directory using 'python -m pip install .'
# import nclcmaps

# https://colorcet.com/index.html
# https://github.com/holoviz/colorcet
import colorcet

# https://arm-doe.github.io/pyart/examples/plotting/plot_choose_a_colormap.html
# JJ Helmus and SM Collis, JORS 2016, doi: 10.5334/jors.119
# import pyart

# https://www.fabiocrameri.ch/colourmaps/
# https://zenodo.org/records/8409685
# import cmcrameri

import matplotlib
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

import numpy as np
from Bio.Seq import translate
import mplcursors

# https://pypi.org/project/blosum/
import blosum

# https://docs.bokeh.org/en/latest/docs/user_guide/interaction/tools.html#ug-interaction-tools-hover-tool
import bokeh.plotting
import bokeh.models

setcontext(ExtendedContext)
c = getcontext()
c.prec = 99

version = 202509041145

myparser = OptionParser(version="%s version %s" % ('%prog', version))
myparser.add_option("--tsv", action="store", type="string", dest="tsv_file_path", default='',
    help="Path to a TAB separated file with 5-columns: ['position', 'original_aa', 'mutant_aa', 'frequency', 'original_codon', 'mutant_codon', 'coverage_per_codon'] or more up to 11-columns ['padded_position', 'position', 'original_aa', 'mutant_aa', 'frequency', 'original_codon', 'mutant_codon', 'observed_codon_count', 'total_codons_per_site', 'frequency_parent', 'frequency_selected'] as we kept extending the file format")
myparser.add_option("--column", action="store", type="string", dest="column_with_frequencies", default='frequency',
    help="Name of a column in input TSV to be used for rendering frequencies [frequency]")
myparser.add_option("--outfile-prefix", action="store", type="string", dest="outfile_prefix", default='',
    help="Output file prefix, eventually also with path. Output files will be PNG, JPG or PDF, HTML")
myparser.add_option("--offset", action="store", type="int", dest="offset", default=0,
    help="Define offset value to be added to every amino acid position in the first column of the TSV input file at runtime. This is to obtain real amino acid position within the protein. Normally protein starts at position 1. Check the first aa in TSV file and provide whatever number to be added to it to get the desired amino acid position in the full-length protein.")
myparser.add_option("--xmin", action="store", type="int", dest="xmin", default=0,
    help="Define minimum X-axis value")
myparser.add_option("--xmax", action="store", type="int", dest="xmax", default=0,
    help="Define maximum X-axis value")
myparser.add_option("--x-axis-bins", action="store", type="int", dest="xaxis_bins", default=0,
    help="Set number of bins (labels) on the X-axis. Could be used to override the amount of major ticks in a different way. Use 20 for example. [default is 0]")
myparser.add_option("--x-axis-major-ticks-spacing", action="store", type="int", dest="xaxis_major_ticks_spacing", default=10,
    help="Set distance between the major ticks on X-axis")
myparser.add_option("--x-axis-minor-ticks-spacing", action="store", type="int", dest="xaxis_minor_ticks_spacing", default=5,
    help="Set distance between the minor ticks on X-axis")
myparser.add_option("--x-axis-label-start", action="store", type="int", dest="xaxis_label_start", default=0,
    help="Set first label value on the X-axis")
myparser.add_option("--aminoacids", action="store_true", dest="aminoacids", default=False,
    help="Draw chart with amino acid residues on Y-axis instead of codons. [default is False]")
myparser.add_option("--show-STOP", action="store_true", dest="showstop", default=False,
    help="Include STOP codons or '*' in charts on Y-axis. [default is False]")
myparser.add_option("--show-INS", action="store_true", dest="showins", default=False,
    help="Include INS in charts on Y-axis. [default is False]")
myparser.add_option("--show-DEL", action="store_true", dest="showdel", default=False,
    help="Include DEL in charts on Y-axis. [default is False]")
myparser.add_option("--show-X", action="store_true", dest="showx", default=False,
    help="Include X in charts on Y-axis in --aminoacids mode. [default is False]")
myparser.add_option("--enable-colorbar", action="store_true", dest="colorbar", default=False,
    help="Enable colorbar [is Disabled by default]")
myparser.add_option("--disable-short-legend", action="store_false", dest="shortlegend", default=True,
    help="Disable short legend in charts on X-axis. [is Enabled by default]")
myparser.add_option("--include-synonymous", action="store_true", dest="include_synonymous", default=False,
    help="Include synonymous changes in --aminoacids output as green diamonds. In codon output they are always shown. [default is False]")
myparser.add_option("--threshold", action="store", type="float", dest="threshold", default=0.001,
    help="Define minimum frequency threshold to display a pictogram in the output. For codon mode use 0.001 and for aa mode use 0.01. [default: 0.001]")
myparser.add_option("--title", action="store", type="string", dest="title", default='',
    help="Set title for the figures, by default trailing '.frequencies.tsv' is stripped from the end of the input filename")
myparser.add_option("--disable-2nd-Y-axis", action="store_true", dest="disable_2nd_Y_axis", default=False,
    help="Disable rendering of the 2nd Y-axis showing sequencing coverage")
myparser.add_option("--legend", action="store_true", dest="legend", default=False,
    help="Draw legend chart. [default is False]")
myparser.add_option("--matrix", action="store", type="str", dest="matrix", default="BLOSUM80",
    help="BLOSUM matrix: BLOSUM45,BLOSUM50,BLOSUM62,BLOSUM80,BLOSUM90 [default is BLOSUM80]")
myparser.add_option("--matrix-file", action="store", type="string", dest="matrix_file", default=None,
    help="Matrix file compatible with BLOSUM matrices, e.g. prot26050-sup-0004-Supinfo03.sm from https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8641535/bin/NIHMS1664401-supplement-supinfo.rar if you do not like default BLOSUM62")
myparser.add_option("--colormap", action="store", type="string", dest="colormap", default='amino_acid_changes_cvd',
    help="Pick a colormap recognized by matplotlib. See https://i.sstatic.net/cmk1J.png [default is coolwarm_r but seismic_r and coolwarm_r are great too]")
myparser.add_option("--dpi", action="store", type="int", dest="dpi", default=600,
    help="DPI resolution for images")
myparser.add_option("--backend", action="store", type="string", dest="backend", default='',
    help="Matplotlibg backend to render resulting figures: agg, wxpython, pyqt5, pyqt6, pycairo, cairocffi [default: unset]\nTo disable Matplolib interactive window being raised up you can set MPLBACKEND=agg env variable.")
myparser.add_option("--debug", action="store", type="int", dest="debug", default=0,
    help="Set debug to some real number")
(myoptions, myargs) = myparser.parse_args()

def alt_translate(seq):
    """Biopython cannot sometimes translate a sequence but one can get around by
    splitting it into codons and merging later back
    https://github.com/biopython/biopython/pull/4992
    https://github.com/biopython/biopython/pull/4992#issuecomment-2865429105
    """

    codons = (seq[i:i+3] for i in range(0, len(seq), 3))
    codons = ("NNN" if "-" in codon and codon != "---" else codon for codon in codons)
    return "".join(translate(codon, gap='-') for codon in codons)


def get_colormap(colormapname):
    try:
        _cmap = matplotlib.colormaps.get_cmap(colormapname) # LinearSegmentedColormap object
        # tab20b = mpl.colormaps['tab20b'].resampled(27)
        # _colors = [matplotlib.colors.to_hex(mpl.colormaps['tab20b'](i), keep_alpha=True) for i in np.linspace(0, 1, 27)]
        # https://matplotlib.org/stable/api/_as_gen/matplotlib.colors.to_hex.html
        # [matplotlib.colors.to_hex(cmap(i), keep_alpha=True) for i in np.linspace(0, 1, _number_of_colors_needed)] # use #rrggbbaa instead of #rrggbb
        # ['#393b79ff', '#393b79ff', '#5254a3ff', '#6b6ecfff', '#9c9edeff', '#9c9edeff', '#637939ff', '#8ca252ff', '#b5cf6bff', '#b5cf6bff', '#cedb9cff', '#8c6d31ff', '#bd9e39ff', '#e7ba52ff', '#e7ba52ff', '#e7cb94ff', '#843c39ff', '#ad494aff', '#ad494aff', '#d6616bff', '#e7969cff', '#7b4173ff', '#7b4173ff', '#a55194ff', '#ce6dbdff', '#de9ed6ff', '#de9ed6ff']
    except:
        try:
            _cmap = plt.get_cmap(colormapname) # matplotlib.colors.ListedColormap object
        except:
            try:
                # Use nclcmaps.cm.colormaps() to list all available colormaps
                _cmap = nclcmaps.cm.get_cmap(colormapname) # matplotlib.colors.ListedColormap object
            except:
                try:
                    _wished_cmapname_prefix, _wished_cmapname_num = '_'.join(myoptions.colormap.split('_')[:-1]), int(myoptions.colormap.split('_')[-1])
                    print(f"Debug: _wished_cmapname_prefix={_wished_cmapname_prefix}, _wished_cmapname_num={_wished_cmapname_num}")
                    _cmap_from_palettable = palettable.colorbrewer.get_map(_wished_cmapname_prefix, 'diverging', _wished_cmapname_num)
                    _cmap = matplotlib.colors.ListedColormap(_cmap_from_palettable.mpl_colors)
                except:
                    try:
                        _cmap_from_palettable = palettable.colorbrewer.get_map(_wished_cmapname_prefix, 'sequential', _wished_cmapname_num)
                        _cmap = matplotlib.colors.ListedColormap(_cmap_from_palettable.mpl_colors)
                    except:
                        try:
                            _cmap_from_palettable = palettable.scientific.get_map(_wished_cmapname_prefix, 'diverging', _wished_cmapname_num)
                            _cmap = matplotlib.colors.ListedColormap(_cmap_from_palettable.mpl_colors)
                        except:
                            try:
                                _cmap_from_palettable = palettable.scientific.get_map(_wished_cmapname_prefix, 'sequential', _wished_cmapname_num)
                                _cmap = matplotlib.colors.ListedColormap(_cmap_from_palettable.mpl_colors)
                            except:
                                try:
                                    _cmap_from_palettable = palettable.tableau.get_map(myoptions.colormap)
                                    _cmap = matplotlib.colors.ListedColormap(_cmap_from_palettable.mpl_colors)
                                except:
                                    # #EFB6D6 pink
                                    # #CC79A7 darker pink
                                    # #ff9900 orange
                                    # #F09163 salmon
                                    # #9D654C brown
                                    # #ffcc00 yellow
                                    # #cccc00 yellow to green
                                    # #A3E4D7 light blue
                                    # #7DCCFF light blue
                                    # #43BA8F green to blue
                                    # #DDFFA0 light green
                                    # #97CE2F green with yellow
                                    # 
                                    # https://github.com/KarstensLab/microshades
                                    _micro_cvd_gray = ("#F5F5F5", "#D6D6D6", "#B7B7B7", "#8B8B8B","#616161")
                                    _micro_cvd_green = ("#DDFFA0",  "#BDEC6F",  "#97CE2F", "#6D9F06","#4E7705")
                                    _micro_cvd_orange = ("#FFD5AF",  "#FCB076","#F09163", "#C17754", "#9D654C")
                                    _micro_cvd_blue = ("#E7F4FF", "#BCE1FF", "#7DCCFF", "#56B4E9","#098BD9")
                                    _micro_cvd_turquoise = ("#A3E4D7", "#48C9B0",  "#43BA8F",  "#009E73", "#148F77")
                                    _micro_cvd_purple = ("#EFB6D6", "#E794C1", "#CC79A7", "#A1527F", "#7D3560")
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
                                        # https://stackoverflow.com/questions/31051488/combining-two-matplotlib-colormaps
                                        # sample the colormaps that you want to use. Use 128 from each so we get 256
                                        # colors in total
                                        _colors1 = plt.cm.gist_heat(np.linspace(0.2, 0.8, 9))
                                        _colors2 = plt.cm.cividis_r(np.linspace(0.2, 0.8, 9))
                                        _colors3 = plt.cm.GnBu_r(np.linspace(0.2, 0.8, 9))

                                        # combine them and build a new colormap
                                        _colors = np.vstack((_colors1, _colors2, _colors3))
                                        _cmap = matplotlib.colors.LinearSegmentedColormap.from_list('merged', _colors)
                                    elif myoptions.colormap == 'merged1':
                                        # https://stackoverflow.com/questions/31051488/combining-two-matplotlib-colormaps
                                        # sample the colormaps that you want to use. Use 128 from each so we get 256
                                        # colors in total
                                        _colors1 = plt.cm.viridis(np.linspace(0.2, 0.8, 9))
                                        _colors2 = plt.cm.coolwarm_r(np.linspace(0., 1, 9))
                                        _colors3 = plt.cm.GnBu_r(np.linspace(0.2, 0.8, 9))

                                        # combine them and build a new colormap
                                        _colors = np.vstack((_colors1, _colors2, _colors3))
                                        _cmap = matplotlib.colors.LinearSegmentedColormap.from_list('merged', _colors)
                                    elif myoptions.colormap == 'dkeenan_26cols':
                                        # https://graphicdesign.stackexchange.com/questions/3682/where-can-i-find-a-large-palette-set-of-contrasting-colors-for-coloring-many-d
                                        # https://dkeenan.com/Colour/
                                        _colors = ["#00B7FF", "#004DFF", "#00FFFF", "#826400", "#580041", "#FF00FF", "#00FF00", "#C500FF", "#B4FFD7", "#FFCA00", "#969600", "#B4A2FF", "#C20078", "#0000C1", "#FF8B00", "#FFC8FF", "#666666", "#FF0000", "#CCCCCC", "#009E8F", "#D7A870", "#8200FF", "#960000", "#BBFF00", "#FFFF00", "#006F00"]
                                        _cmap = matplotlib.colors.ListedColormap(_colors, "dkeenan_26cols") # 26 colors
                                        myoptions.colormap = 'dkeenan_26cols'
                                    else:
                                        # _amino_acid_changes_cvd_cmap = ["#F5F5F5", "#D6D6D6", "#B7B7B7", "#8B8B8B", "#ff0000", "#f16264", "#ffc8ff", "#EFB6D6", "#CC79A7", "#ff9900", "#F09163", "#9c644b", "#ffcc00", "#cccc00", "#A3E4D7", "#7DCCFF", "#0042ff", "#DDFFA0", "#97CE2F", "#219f11", "#00ff04", "#c20078", "#ff00fd", "#c500ff", "#8200ff", "#960000", "#580041"]
                                        _colors = ["#930000", "#930000", "#930000", "#930000", "#930000", "#930000", "#960000", "#580041", "#8200ff", "#c500ff", "#ff00fd", "#eea1d0", "#CC79A7", "#cc0000", "#ff0000", "#ff4f00", "#ff7c7c", "#ff9999", "#c58a24", "#9c644b", "#ffff00", "#ffcc00", "#cccc00", "#7DCCFF", "#0042ff", "#0000ff", "#D6D6D6", "#B7B7B7", "#8B8B8B", "#97CE2F", "#bbff00", "#219f11", "#00ff04", "#930000", "#930000", "#930000", "#930000", "#930000", "#930000"]
                                        _cmap = matplotlib.colors.ListedColormap(_colors, "amino_acid_changes_cvd", len(_colors)) # 39 colors
                                        #plt.cm.register_cmap(name='amino_acid_changes_cvd', cmap=_cmap)
                                        myoptions.colormap = 'amino_acid_changes_cvd'
                                        return _cmap, _colors

    return _cmap, _colors


def adjust_size_and_color(frequency, old_codon_or_aa, new_codon_or_aa, matrix, min_score, max_score, cmap):
    """Define colors for basic type of scatterplot figures. The circles have variable diameter
    in a range from 0 to 100%, each. The total sum of amino acid residues in a particular position
    in a protein is 100%. Provided there can be multiple aa residues encoded by even more codons,
    the diameter is not always a 100% (not a full-size).

    The checks need to anticipate these scenarios:
    old_codon_or_aa = 'NNN'
    new_codon_or_aa = 'N' # incomplete new codon

    old_codon_or_aa = 'NNN'
    new_codon_or_aa = 'DEL'

    old_codon_or_aa = '--'
    new_codon_or_aa = 'NNN'

    https://stackoverflow.com/questions/25408393/getting-individual-colors-from-a-color-map-in-matplotlib
    https://stackoverflow.com/questions/34314356/how-to-view-all-colormaps-available-in-matplotlib/68317686#68317686
    https://i.sstatic.net/cmk1J.png

    https://matplotlib.org/stable/gallery/color/named_colors.html

    https://www.ncbi.nlm.nih.gov/IEB/ToolBox/C_DOC/lxr/source/data/BLOSUM45
    https://www.ncbi.nlm.nih.gov/IEB/ToolBox/C_DOC/lxr/source/data/BLOSUM50
    https://www.ncbi.nlm.nih.gov/IEB/ToolBox/C_DOC/lxr/source/data/BLOSUM62
    https://www.ncbi.nlm.nih.gov/IEB/ToolBox/C_DOC/lxr/source/data/BLOSUM80
    https://www.ncbi.nlm.nih.gov/IEB/ToolBox/C_DOC/lxr/source/data/BLOSUM90
    """

    _len_old_codon_or_aa = len(old_codon_or_aa)
    _len_new_codon_or_aa = len(new_codon_or_aa)

    if myoptions.aminoacids:
        # in aa mode we always pass to this routine a codon instead of the original amino acid, dunno why but at least when rendering a static figure legend we pass down 'NNN'
        _old_codon_or_aa = alt_translate(old_codon_or_aa)
        _len_old_codon_or_aa = len(_old_codon_or_aa)
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

    # define color palette, use 'coolwarm nr. 90' from https://i.sstatic.net/cmk1J.png or 'rainbow nr. 138' or 'seismic_r 141' or 'RdYIBu_r 51'
    # https://datascientyst.com/full-list-named-colors-pandas-python-matplotlib/

    # do not get confused by the 3-letter DEL and INS values
    if len(_old_codon_or_aa) > 1 and old_codon_or_aa not in ('---', 'DEL', 'INS'):
        _codon_on_input = True
    elif len(_new_codon_or_aa) > 1 and new_codon_or_aa not in ('---', 'DEL', 'INS'):
        _codon_on_input = True
    else:
        _codon_on_input = False

    if _codon_on_input:
        if str(_old_codon_or_aa)[0].upper() != str(_new_codon_or_aa)[0].upper(): # the 1st nucleotide in a codon is different
            _significant_change = True
        elif str(_old_codon_or_aa)[1].upper() != str(_new_codon_or_aa)[1].upper(): # the 2nd nucleotide in a codon is different
            _significant_change = True
        else:
            _significant_change = False
    elif not _old_codon_or_aa:
        raise ValueError("Aieee, old_codon_or_aa='%s' is empty" % str(old_codon_or_aa))
    elif not _new_codon_or_aa:
        raise ValueError("Aieee, new_codon_or_aa='%s' is empty" % str(new_codon_or_aa))
    else:
        if new_codon_or_aa in ('---', 'DEL', 'INS', '*'):
            _significant_change = True
        elif _old_codon_or_aa in ('---', 'DEL', 'INS', '*'):
            _significant_change = True
        else:
            # some single-letter amino acid residue or '*' or 'DEL' or 'INS'
            if myoptions.debug:
                print("Info: some single-letter but neither asterisk nor dash amino acid residue: _old_codon_or_aa=%s, new_codon_or_aa='%s'" % (str(_old_codon_or_aa), str(new_codon_or_aa)))
            _significant_change = False

    if frequency < Decimal(0):
        _size, _color = frequency, 'red'
    elif frequency < myoptions.threshold:
        _size, _color = 0, 'blue'  # ignore values below given threshold, e.g. below 0.001, use blue color but it is ignored later on anyway
    elif frequency >= Decimal(0.05) and frequency <= Decimal(0.3):
        if _significant_change:
            _color = 'magenta'
        else:
            _color = 'gray'
        _size = frequency
    elif frequency > Decimal(0.3):
        if _significant_change:
            _color = 'paleturquoise'
        else:
            _color = 'black'
        _size = frequency
    else:
        if _significant_change:
            _color = 'yellow'  # other frequencies
        else:
            _color = 'skyblue'  # other frequencies, skyblue color
        _size = frequency

    # zap the color with a color from gradient, make sure the middle of the gradient should be symmetrical around zero (white) [-11 ... 0 ... +11] which is actually [min_score ... 0 ... max_score]
    _half_size = int(max(abs(min_score), max_score))
    _number_of_colors_needed = int(max(abs(min_score), max_score) * 2 + 1) # calculate how many colors are needed to cover all values in the matrix while 0 would be in the middle and convert to int from float

    # make sure the BLOSUM matrix contains the AA pair otherwise blosum module returns -inf (infinity) triggering OverflowError
    if _new_codon_or_aa in ('---', 'DEL', 'INS'):
        _score = -11
    elif _codon_on_input:
        _score = int(matrix[alt_translate(_old_codon_or_aa)][alt_translate(_new_codon_or_aa)])
    else:
        try:
            _score = int(matrix[_old_codon_or_aa][_new_codon_or_aa])
        except:
            raise ValueError("Cannot get a score for myoptions.matrix='%s', _old_codon_or_aa='%s', _new_codon_or_aa='%s'" % (myoptions.matrix, _old_codon_or_aa, _new_codon_or_aa))
        
    if not _score:
        _pos = max(abs(min_score), max_score) # get the position of the 0 in the middle of the list of colors
    elif _score < 0:
        _pos = max(abs(min_score), max_score) - abs(int(_score)) # do not add 1 for the 0 in the middle of the list of colors because python counts items anyway from 0
    else:
        _pos = max(abs(min_score), max_score) + int(_score) # do not add 1 for the 0 in the middle of the list of colors because python counts items anyway from 0

    if _size:
        # https://matplotlib.org/stable/api/_as_gen/matplotlib.colors.to_hex.html
        _colors = [matplotlib.colors.to_hex(cmap(i), keep_alpha=True) for i in np.linspace(0, 1, _number_of_colors_needed)] # use #rrggbbaa instead of #rrggbb
        if myoptions.debug: print([x for x in _colors])

        if old_codon_or_aa.upper() == new_codon_or_aa.upper():
            _size, _color = 0, 'palegreen' # do not draw circles for unchanged amino acids
        elif _old_codon_or_aa.upper() == _new_codon_or_aa.upper():
            _size, _color = 0, 'palegreen' # do not draw circles for unchanged amino acids
        elif old_codon_or_aa in ('---', 'DEL', 'INS', '*') or new_codon_or_aa in ('---', 'DEL', 'INS', '*', 'TGA', 'TAA', 'TAG'):
            _color = 'red'
        elif _old_codon_or_aa in ('X', 'NNN') or new_codon_or_aa in ('X', 'NNN'):
            _color = 'gray'
        elif _codon_on_input:
            if alt_translate(_old_codon_or_aa) == alt_translate(_new_codon_or_aa):
                _color = 'palegreen'
            elif alt_translate(_new_codon_or_aa) == 'X':
                _color = 'gray'
            else:
                if myoptions.debug:
                    sys.stdout.write("Info: Translating %s to %s to fetch color from cmap," % (_old_codon_or_aa, _new_codon_or_aa))
                _color = _colors[_pos] # get the color by its index position in the ListedColormap
                if myoptions.debug:
                    sys.stdout.write(" which has yielded %s and %s%s" % (_pos, str(_color), os.linesep))
        else:
            if myoptions.debug:
                sys.stdout.write("Info: Fetching color for %s to %s directly from cmap," % (_old_codon_or_aa, _new_codon_or_aa))
            _color = _colors[_pos]
            if myoptions.debug:
                sys.stdout.write(" which has yielded %s and %s%s" % (_pos, str(_color), os.linesep))
        if _pos < 0 or _pos > _number_of_colors_needed:
            raise ValueError("_pos=%s, Only %d colors are expected, widen the range int(max(abs(%d), %d) * 2 + 1) if you need more." % (str(_pos), _number_of_colors_needed, min_score, max_score))

    if not _size:
        return _score, 0.00005, _color
    elif not frequency < myoptions.threshold:
        if old_codon_or_aa.upper() == new_codon_or_aa.upper(): # do not draw circles for unchanged amino acids
            if myoptions.debug: print("Debug: For _old_codon_or_aa = '%s', _new_codon_or_aa = '%s' returning _size = '%s', _color = '%s'" % (_old_codon_or_aa, _new_codon_or_aa, 0, 'palegreen'))
            return _score, 0.00005, 'palegreen'
        elif _old_codon_or_aa.upper() == _new_codon_or_aa.upper(): # do not draw circles for unchanged amino acids
            if myoptions.debug: print("Debug: For _old_codon_or_aa = '%s', _new_codon_or_aa = '%s' returning _size = '%s', _color = '%s'" % (_old_codon_or_aa, _new_codon_or_aa, 0, 'palegreen'))
            return _score, 0.00005, 'palegreen'
        else:
            if myoptions.debug: print("Debug: For _old_codon_or_aa = '%s', _new_codon_or_aa = '%s' returning _size = '%s', _color = '%s'" % (_old_codon_or_aa, _new_codon_or_aa, frequency, _color))
            return _score, frequency, _color # do not draw circles for frequencies below threshold
    else:
        if myoptions.debug: print("Debug: For _old_codon_or_aa = '%s', _new_codon_or_aa = '%s' returning _size = '%s', _color = '%s'" % (_old_codon_or_aa, _new_codon_or_aa, 0, _color))
        return _score, 0.00005, _color


def adjust_size_and_color_neutralized_escape(neutralized_parent_difference, old_codon_or_aa, new_codon_or_aa, matrix):
    "Used only for neutralized_parent_difference and escape_parent_difference figure types."

    if neutralized_parent_difference < -0.001:
        _size = abs(neutralized_parent_difference)
        _color = 'red'
    elif neutralized_parent_difference > 0.001:
        if not _codon_on_input:
            _size, _color = 0, 'palegreen' # do not draw circles for unchanged amino acids
        else:
            _size = neutralized_parent_difference
            _color = 'palegreen'
    else:        # Values close to zero (between -0.001 and 0.001) non-visible; use size 5 to retain a controlling handle
        _size = 0
        _color = 'black'
    return _size, _color


def adjust_size_and_color_weighted(weighted_diff_escape_neutralized, old_codon_or_aa, new_codon_or_aa, matrix, generic_circle_size=5000, weighted_diff_escape_neutralized_size1=3000, weighted_diff_escape_neutralized_size2=2000): # use for 'weighted_diff_escape_neutralized'
    "Used only for weighted_diff_escape_neutralized figure type."

    if weighted_diff_escape_neutralized > 1:
        # Values above 1 will be scaled more prominently and colored dark blue
        _size = (weighted_diff_escape_neutralized) * generic_circle_size
        _color = 'black'
    elif weighted_diff_escape_neutralized > 0.1:
        # Values between 0.001 and 1 will be scaled and colored skyblue
        _size = (weighted_diff_escape_neutralized) * weighted_diff_escape_neutralized_size1
        _color = 'darkblue'
    elif weighted_diff_escape_neutralized > 0.01:
        # Values between 0.001 and 1 will be scaled and colored skyblue
        _size = (weighted_diff_escape_neutralized) * weighted_diff_escape_neutralized_size2
        _color = 'skyblue'
    else:
        _size = abs(weighted_diff_escape_neutralized) * 0
        _color = 'black'
    return _size, _color


def main():
    if myoptions.matrix_file and os.path.exists(myoptions.matrix_file):
        print(myoptions.matrix)
        if not myoptions.matrix:
            myoptions.matrix = re.sub(r'.*/', '', myoptions.matrix_file)
        print(myoptions.matrix)
        _matrix_type, _matrix_num = re.sub(r'\d+', '', myoptions.matrix), int(re.sub(r'[a-zA-Z]+', '', myoptions.matrix))
        _matrix = blosum.BLOSUM(myoptions.matrix_file) # if compatible with BLOSUM, import the custom file
        _matrix_name = myoptions.matrix_file.split(os.path.sep)[-1]
        _commonly_used_matrix = False
        myoptions.matrix = _matrix_name
    else:
        _matrix_type, _matrix_num = re.sub(r'\d+', '', myoptions.matrix), int(re.sub(r'[a-zA-Z]+', '', myoptions.matrix))
        if _matrix_type == 'BLOSUM':
            _matrix = blosum.BLOSUM(_matrix_num) # 45,50,62,80,90
            _matrix_name = "BLOSUM%d" % (_matrix_num)
        else:
            sys.stderr.write("Warning: Unexpected matrix type %s, falling back to BLOSUM\n" % str(_matrix_type))
            _matrix = blosum.BLOSUM(_matrix_num) # 45,50,62,80,90
            _matrix_name = "BLOSUM%d" % (_matrix_num)
            myoptions.matrix = _matrix_name
        _commonly_used_matrix = True

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
    _min_theoretical_score, _max_theoretical_score = int(min(_theoretical_scores)), int(max(_theoretical_scores)) # also convert from float to int

    print(f"Info: Parsing input file {myoptions.tsv_file_path}")
    if not myoptions.tsv_file_path:
        raise RuntimeError("Please provide an input TSV file via --tsv")
    df = pd.read_csv(myoptions.tsv_file_path, sep='\t', header='infer', na_filter=False, na_values=[None])#, nrows=500)

    if 'position' not in df.columns:
        del(df)
        print("Info: Autodetected old TSV file format without a header in %s, assigning default column names" % myoptions.tsv_file_path)
        df = pd.read_csv(myoptions.tsv_file_path, sep='\t', header=0, na_filter=False, na_values=[None])#, nrows=500)
        # Assign column names
        # released results on Zenodo: 430	T	K	0.000386	ACA	AAA
        print(f"Info: The file {myoptions.tsv_file_path} contained initially these columns: {str(df.columns)}")
        _count_columns = len(df.columns.values)
        if _count_columns == 9:
            # 340     340     E       E       0.012149        GAA     GAG     97849   8054365
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
            raise RuntimeError("Unexpected number of columns in the %s file" % myoptions.tsv_file_path)
    else:
        print(f"Info: Autodetected new TSV file format with a header in {myoptions.tsv_file_path}")
    print(f"Info: The file {myoptions.tsv_file_path} contains now these columns: {str(df.columns)}")

    #print(f"Info: df={str(df)}")
    # sanitize input aa positions and add myoptions.offset to each value
    df['position'] = df['position'] + int(myoptions.offset)

    _before = len(df['mutant_codon'])
    try:
        df = df.loc[df['mutant_codon'].str.match('[ATGCatgc-][ATGCatgc-][ATGCatgc-]')] # discard [nN] but keep '-', beware GISAID also contains other IUPAC codes
    except Exception as _exc:
        raise ValueError("Cannot parse column df['mutant_codon']=%s containing %d values" % (str(df['mutant_codon']), len(df['mutant_codon'])))
    if not myoptions.showstop:
        # discard sequencing errors appearing as STOPs
        df = df.loc[~df['mutant_aa'].isin(['*'])]
    if not myoptions.showdel:
        # discard sequencing errors appearing as DEL
        df = df.loc[~df['mutant_aa'].isin(['DEL'])]
    if not myoptions.showins:
        df = df.loc[~df['mutant_aa'].isin(['INS'])]
    df = df.loc[~df['mutant_aa'].isin(['X'])] # discard frame-breaking changes (codons contain one or two dashes), TODO: maybe we should check if any dash was in the reference_sequence ?
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
    # print(df)
    print(f"Info: Writing into {_outfile_prefix}.actually_rendered.tsv")
    df.to_csv(f"{_outfile_prefix}.actually_rendered.tsv", sep='\t', header=None, index=False, float_format='{:7.6f}'.format) # dump back the edited dataframe for eventual checks and disable scientific notation of float numbers
    if '.frequencies.tsv' in myoptions.tsv_file_path:
        _count_filename = _outfile_prefix + '.count'
        if os.path.exists(_count_filename):
            try:
                _aln_handle = open(_count_filename)
            except:
                _aln_rows = '0'
                _aln_handle = None
        else:
            _aln_rows = '0'
            _aln_handle = None
    else:
        _aln_rows = '0'
        _aln_handle = None

    if not myoptions.title:
        title_data = myoptions.tsv_file_path.replace('.frequencies.tsv', '')
    else:
        title_data = myoptions.title

    if _aln_handle:
        _aln_rows = _aln_handle.readline()
        _aln_handle.close()

    print(f"Info: Title will be {title_data}")

    # create an empty table pre-filled with zeroes, without ['B', 'Z']
    #amino_acids = ['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V', '*']
    amino_acids = ['C', 'R', 'K', 'E', 'Q', 'D', 'N', 'T', 'S', 'H', 'M', 'P', 'W', 'Y', 'F', 'V', 'L', 'I', 'A', 'G'] # 'INS', 'X', 'DEL'] #, None] # sorted by physicochemical properties
    if myoptions.showins:
        amino_acids.append('INS')
    if myoptions.showx:
        amino_acids.append('X')
    if myoptions.showstop:
        amino_acids.append('*')
    if myoptions.showdel:
        amino_acids.append('DEL')

    unique_aa_positions = [x for x in df['position'].unique()]
    _min_aa_pos = min(unique_aa_positions)
    _max_aa_pos = max(unique_aa_positions)
    #_number_of_insertions = _padded_reference_sequence.count('-') # count number of insertions based on the number of dashes in a padded reference sequence
    _number_of_insertions = 0
    unique_aa_positions = [x for x in range(_min_aa_pos, _max_aa_pos + 1 + _number_of_insertions)] # ensure we have a continuous range, eventually inject INS positions, add +1 to accommodate for the pandas dataframe index column
    if myoptions.debug: print(f"Debug: len(unique_aa_positions)={len(unique_aa_positions)}, unique_aa_positions: {str(unique_aa_positions)}")
    unique_codon_positions = list(unique_aa_positions)
    if myoptions.debug: print(f"Debug: len(unique_codon_positions)={len(unique_codon_positions)}, unique_codon_positions: {str(unique_codon_positions)}")

    # these are all theoretically appearing codons in a sample
    codons_whitelist = ['TTT', 'TTC', 'TTA', 'TTG', 'CTT', 'CTC', 'CTA', 'CTG', 'ATT', 'ATC', 'ATA', 'ATG', 'GTT', 'GTC', 'GTA', 'GTG', 'TCT', 'TCC', 'TCA', 'TCG', 'CCT', 'CCC', 'CCA', 'CCG', 'ACT', 'ACC', 'ACA', 'ACG', 'GCT', 'GCC', 'GCA', 'GCG', 'TAT', 'TAC', 'TAA', 'TAG', 'CAT', 'CAC', 'CAA', 'CAG', 'AAT', 'AAC', 'AAA', 'AAG', 'GAT', 'GAC', 'GAA', 'GAG', 'TGT', 'TGC', 'TGA', 'TGG', 'CGT', 'CGC', 'CGA', 'CGG', 'AGT', 'AGC', 'AGA', 'AGG', 'GGT', 'GGC', 'GGA', 'GGG']
    if not myoptions.showstop:
        for _stopcodon in ('TGA', 'TAA', 'TAG'):
            codons_whitelist.remove(_stopcodon)
    
    if myoptions.debug: print("Debug: Translating pre-made all possible codons")
    codons_whitelist_aa = [alt_translate(_codon) for _codon in codons_whitelist]
    # ensure the resulting codons_whitelist2 will have same length as the codons later parsed from TSV file into _table.index
    if 'DEL' in amino_acids:
        codons_whitelist.append('---')
        codons_whitelist_aa.append('DEL')

    if len(codons_whitelist) != len(codons_whitelist_aa):
        raise ValueError("Length of codons_whitelist is %d which is not equal to codons_whitelist_aa with %d" % (len(codons_whitelist), len(codons_whitelist_aa)))
    sorted_whitelist = sorted(zip(codons_whitelist, codons_whitelist_aa), key=lambda x: x[1])

    if myoptions.debug: print(f"Debug: sorted_whitelist={str(sorted_whitelist)}")
    # Info: sorted_whitelist=[('TAA', '*'), ('TAG', '*'), ('TGA', '*'), ('GCT', 'A'), ('GCC', 'A'), ('GCA', 'A'), ('GCG', 'A'), ('TGT', 'C'), ('TGC', 'C'), ('GAT', 'D'), ('GAC', 'D'), ('GAA', 'E'), ('GAG', 'E'), ('TTT', 'F'), ('TTC', 'F'), ('GGT', 'G'), ('GGC', 'G'), ('GGA', 'G'), ('GGG', 'G'), ('CAT', 'H'), ('CAC', 'H'), ('ATT', 'I'), ('ATC', 'I'), ('ATA', 'I'), ('AAA', 'K'), ('AAG', 'K'), ('TTA', 'L'), ('TTG', 'L'), ('CTT', 'L'), ('CTC', 'L'), ('CTA', 'L'), ('CTG', 'L'), ('ATG', 'M'), ('AAT', 'N'), ('AAC', 'N'), ('CCT', 'P'), ('CCC', 'P'), ('CCA', 'P'), ('CCG', 'P'), ('CAA', 'Q'), ('CAG', 'Q'), ('CGT', 'R'), ('CGC', 'R'), ('CGA', 'R'), ('CGG', 'R'), ('AGA', 'R'), ('AGG', 'R'), ('TCT', 'S'), ('TCC', 'S'), ('TCA', 'S'), ('TCG', 'S'), ('AGT', 'S'), ('AGC', 'S'), ('ACT', 'T'), ('ACC', 'T'), ('ACA', 'T'), ('ACG', 'T'), ('GTT', 'V'), ('GTC', 'V'), ('GTA', 'V'), ('GTG', 'V'), ('TGG', 'W'), ('TAT', 'Y'), ('TAC', 'Y')]
    
    # >>> [pairs[0] + ' (' + pairs[1] + ')' for pairs in sorted_whitelist]
    # ['TAA (*)', 'TAG (*)', 'TGA (*)', 'GCT (A)', 'GCC (A)', 'GCA (A)', 'GCG (A)', 'TGT (C)', 'TGC (C)', 'GAT (D)', 'GAC (D)', 'GAA (E)', 'GAG (E)', 'TTT (F)', 'TTC (F)', 'GGT (G)', 'GGC (G)', 'GGA (G)', 'GGG (G)', 'CAT (H)', 'CAC (H)', 'ATT (I)', 'ATC (I)', 'ATA (I)', 'AAA (K)', 'AAG (K)', 'TTA (L)', 'TTG (L)', 'CTT (L)', 'CTC (L)', 'CTA (L)', 'CTG (L)', 'ATG (M)', 'AAT (N)', 'AAC (N)', 'CCT (P)', 'CCC (P)', 'CCA (P)', 'CCG (P)', 'CAA (Q)', 'CAG (Q)', 'CGT (R)', 'CGC (R)', 'CGA (R)', 'CGG (R)', 'AGA (R)', 'AGG (R)', 'TCT (S)', 'TCC (S)', 'TCA (S)', 'TCG (S)', 'AGT (S)', 'AGC (S)', 'ACT (T)', 'ACC (T)', 'ACA (T)', 'ACG (T)', 'GTT (V)', 'GTC (V)', 'GTA (V)', 'GTG (V)', 'TGG (W)', 'TAT (Y)', 'TAC (Y)']
    
    # sort the pairs according to physicochemical properties in amino_acids
    final_sorted_whitelist = [tuple for x in amino_acids for tuple in sorted_whitelist if tuple[1] == x]
    if myoptions.debug: print(f"Debug: final_sorted_whitelist={str(final_sorted_whitelist)}")

    # >>> [tuple for x in amino_acids for tuple in sorted_whitelist if tuple[1] == x]
    # [('TGT', 'C'), ('TGC', 'C'), ('CGT', 'R'), ('CGC', 'R'), ('CGA', 'R'), ('CGG', 'R'), ('AGA', 'R'), ('AGG', 'R'), ('AAA', 'K'), ('AAG', 'K'), ('GAA', 'E'), ('GAG', 'E'), ('CAA', 'Q'), ('CAG', 'Q'), ('GAT', 'D'), ('GAC', 'D'), ('AAT', 'N'), ('AAC', 'N'), ('ACT', 'T'), ('ACC', 'T'), ('ACA', 'T'), ('ACG', 'T'), ('TCT', 'S'), ('TCC', 'S'), ('TCA', 'S'), ('TCG', 'S'), ('AGT', 'S'), ('AGC', 'S'), ('CAT', 'H'), ('CAC', 'H'), ('ATG', 'M'), ('CCT', 'P'), ('CCC', 'P'), ('CCA', 'P'), ('CCG', 'P'), ('TGG', 'W'), ('TAT', 'Y'), ('TAC', 'Y'), ('TTT', 'F'), ('TTC', 'F'), ('GTT', 'V'), ('GTC', 'V'), ('GTA', 'V'), ('GTG', 'V'), ('TTA', 'L'), ('TTG', 'L'), ('CTT', 'L'), ('CTC', 'L'), ('CTA', 'L'), ('CTG', 'L'), ('ATT', 'I'), ('ATC', 'I'), ('ATA', 'I'), ('GCT', 'A'), ('GCC', 'A'), ('GCA', 'A'), ('GCG', 'A'), ('GGT', 'G'), ('GGC', 'G'), ('GGA', 'G'), ('GGG', 'G'), ('TAA', '*'), ('TAG', '*'), ('TGA', '*')]
    # >>> 

    codons_whitelist2 = [x[0] for x in final_sorted_whitelist] # get back the list of codons in the order used in final_sorted_whitelist but also in sync with the order on Y-axis further down
    
    if myoptions.debug: print(f"Debug: codons_whitelist2={str(codons_whitelist2)}")
    
    # Define the genetic code and codons with an additional parameter (1 to 65)
    genetic_code = {
        "TTT": ("F", 1), "TTC": ("F", 2), "TTA": ("L", 3), "TTG": ("L", 4),
        "CTT": ("L", 5), "CTC": ("L", 6), "CTA": ("L", 7), "CTG": ("L", 8),
        "ATT": ("I", 9), "ATC": ("I", 10), "ATA": ("I", 11), "ATG": ("M", 12),
        "GTT": ("V", 13), "GTC": ("V", 14), "GTA": ("V", 15), "GTG": ("V", 16),
        "TCT": ("S", 17), "TCC": ("S", 18), "TCA": ("S", 19), "TCG": ("S", 20),
        "CCT": ("P", 21), "CCC": ("P", 22), "CCA": ("P", 23), "CCG": ("P", 24),
        "ACT": ("T", 25), "ACC": ("T", 26), "ACA": ("T", 27), "ACG": ("T", 28),
        "GCT": ("A", 29), "GCC": ("A", 30), "GCA": ("A", 31), "GCG": ("A", 32),
        "TAT": ("Y", 33), "TAC": ("Y", 34), "TAA": ("Stop", 35), "TAG": ("Stop", 36),
        "CAT": ("H", 37), "CAC": ("H", 38), "CAA": ("Q", 39), "CAG": ("Q", 40),
        "AAT": ("N", 41), "AAC": ("N", 42), "AAA": ("K", 43), "AAG": ("K", 44),
        "GAT": ("D", 45), "GAC": ("D", 46), "GAA": ("E", 47), "GAG": ("E", 48),
        "TGT": ("C", 49), "TGC": ("C", 50), "TGA": ("Stop", 51), "TGG": ("W", 52),
        "CGT": ("R", 53), "CGC": ("R", 54), "CGA": ("R", 55), "CGG": ("R", 56),
        "AGT": ("S", 57), "AGC": ("S", 58), "AGA": ("R", 59), "AGG": ("R", 60),
        "GGT": ("G", 61), "GGC": ("G", 62), "GGA": ("G", 63), "GGG": ("G", 64),
        "NNN": ("Other", 65),  # Placeholder for anything else
    }
    old_aa_table = pd.DataFrame(Decimal(0), index=amino_acids, columns=unique_aa_positions)
    new_aa_table = pd.DataFrame(Decimal(0), index=amino_acids, columns=unique_aa_positions)
    old_codon_table = pd.DataFrame(Decimal(0), index=codons_whitelist2, columns=unique_codon_positions)
    new_codon_table = pd.DataFrame(Decimal(0), index=codons_whitelist2, columns=unique_codon_positions)

    _very_leftmost_aa_pos = None
    _calculated_aa_offset = 0
    # make tables with yet another number of rows summing up eventually the frequencies
    for df_index, row in df.iterrows():
        # 1       G       DEL      0.012625       GGA     ---
        #print(f"Row: {str(row)}")

        # It is not necessary to skip N-containing codons as we anyway draw just those
        # in codons list. Skipping some rows would make new_aa_table and new_codon_table
        # have a different amount of rows, breaking slicing
        position = row['position'] # real position_in_protein
        if _very_leftmost_aa_pos is None:
            _very_leftmost_aa_pos = int(position)
            _calculated_aa_offset = _very_leftmost_aa_pos - myoptions.offset + 1 # if AA positions in the input file do NOT start from the first-one the numbering of sites gets shifted, so calculate the offset
            if myoptions.debug: print(f"Debug: calculated offset is {_calculated_aa_offset}")
        old_amino_acid = row['original_aa']
        new_amino_acid = row['mutant_aa']
        frequency = Decimal(row[myoptions.column_with_frequencies])
        old_codon = row['original_codon'].upper()
        new_codon = row['mutant_codon'].upper()

        if (not myoptions.aminoacids or old_amino_acid != new_amino_acid or myoptions.include_synonymous) and not np.abs(frequency) < myoptions.threshold: # discard sequencing noise by zapping the values instead of deleting altogether
            _old_value = Decimal(0)
            try:
                _old_value = new_aa_table.at[new_amino_acid, position]
            except KeyError:
                new_aa_table.at[new_amino_acid, position] = Decimal(frequency)
            except TypeError:
                raise TypeError("Weird value %s" % new_aa_table.at[new_amino_acid, position])
            else:
                new_aa_table.at[new_amino_acid, position] = Decimal(_old_value) + Decimal(frequency) # the position with the new amino acid may appear multiple times so sum up the values

            _old_value =Decimal( 0)
            try:
                _old_value = old_aa_table.at[old_amino_acid, position]
            except KeyError:
                old_aa_table.at[old_amino_acid, position] = Decimal(frequency)
            except TypeError:
                raise TypeError("Weird value %s" % old_aa_table.at[old_amino_acid, position])
            else:
                old_aa_table.at[old_amino_acid, position] = Decimal(_old_value) + Decimal(frequency) # the position with the original amino acid will appear multiple times so sum up the values

            _old_value = Decimal(0)
            try:
                _old_value = Decimal(old_codon_table.at[old_codon, position])
            except KeyError:
                old_codon_table.at[old_codon, position] = Decimal(frequency) # the position with the original codon will appear multiple times so sum up the values
            except TypeError:
                raise TypeError("Weird value %s" % old_codon_table.at[old_codon, position])
            else:
                old_codon_table.at[old_codon, position] = Decimal(old_codon_table.at[old_codon, position]) + Decimal(frequency) # rewrite np.float64 into Decimal

            _old_value = Decimal(0)
    #        try:
    #            new_codon_table.at[new_codon, position] = Decimal(new_codon_table.at[new_codon, position]) + Decimal(frequency)
    #        except KeyError:
            new_codon_table.at[new_codon, position] = Decimal(frequency) # the position and a codon should appear just once so no reason to sum up new value with previous
    #        except TypeError:
    #            # sanitize nan values in the table
    #            raise TypeError("Weird value %s" % new_codon_table.at[new_codon, position])
            if myoptions.debug: print(f"Debug: OriginalDataFrameRowNumber: {df_index}, Old: {old_amino_acid}, New: {new_amino_acid}, Frequency: {frequency}")

    if myoptions.debug:
        # print tables with summed up frequencies, notably each has different amount of rows which differ from the original TSV row number as well
        for t in (old_aa_table, new_aa_table, old_codon_table, new_codon_table):
            print(f"Debug: len({t})={len(t)}")

    if myoptions.debug:
        if myoptions.aminoacids:
            print(old_aa_table)
            print(new_aa_table)
        else:
            print(old_codon_table)
            print(new_codon_table)


    #print(plt.rcParams["font.sans-serif"][0])
    #print(plt.rcParams["font.monospace"][0])
    matplotlib.rcParams['font.family'] = 'monospace'
    #matplotlib.rcParams["font.monospace"] = ["DejavuSans"]
    # trigger core fonts for PDF backend
    plt.rcParams["pdf.use14corefonts"] = True
    # trigger core fonts for PS backend
    plt.rcParams["ps.useafm"] = True
    #plt.text(2, 0.65, 'Version %s' % version)

    print(f"Info: matplotlib.get_backend={matplotlib.get_backend()}")

    _figure, (_ax1, _ax3, _ax4) = plt.subplots(1, 3, figsize=(16, 9), width_ratios=[55, 1, 6]) # default was 6.4 x 4.8 in and 100dpi

    # _aln_rows contain number of rows in the ALN file typically each read is counted so roughly divide it by two to get at about number of amplicons sequenced
    if myoptions.aminoacids:
        if myoptions.shortlegend:
            _xlabel = 'Amino acid position'
        else:
            _xlabel = 'Amino acid position%sbased on %s ALN rows, matrix %s, colormap %s, mutation_scatter_plot.py %s' % (os.linesep, _aln_rows.strip(os.linesep), _matrix_name, myoptions.colormap, version)
    else:
        if myoptions.shortlegend:
            _xlabel = 'Codon position'
        else:
            _xlabel = 'Codon position%sbased on %s ALN rows, matrix %s, colormap %s, mutation_scatter_plot.py %s' % (os.linesep, _aln_rows.strip(os.linesep), _matrix_name, myoptions.colormap, version)
    _ax1.set_xlabel(_xlabel, fontsize=14)
    if myoptions.aminoacids:
        _ax1.set_ylabel('Introduced amino acid changes', fontsize=14)
        _ax1.set_title(title_data, fontsize=18)
        if myoptions.xmin:
            # --offset 299 --xmin 305 --xmax 485 ## respect xmin and do not make the figure wider spanning into primer region
            _xmin = myoptions.xmin
        else:
            _xmin = min(unique_aa_positions) - 1

        if myoptions.xmax:
            _xmax = myoptions.xmax
        else:
            _xmax = max(unique_aa_positions) + 1
    else:
        _ax1.set_ylabel('Introduced codon changes', fontsize=14)
        _ax1.set_title(title_data, fontsize=14)
        if myoptions.xmin:
            # --offset 299 --xmin 305 --xmax 485 ## respect xmin and do not make the figure wider spanning into primer region
            _xmin = myoptions.xmin
        else:
            _xmin = min(unique_codon_positions) - 1

        if myoptions.xmax:
            _xmax = myoptions.xmax
        else:
            _xmax = max(unique_codon_positions) + 1

    _ax1.set_xlim(myoptions.xaxis_label_start or _xmin, _xmax) # start X-axis from 1, not zero
    _ax1.xaxis.set_major_locator(ticker.MultipleLocator(myoptions.xaxis_major_ticks_spacing))
    _ax1.xaxis.set_minor_locator(ticker.MultipleLocator(myoptions.xaxis_minor_ticks_spacing))

    if myoptions.debug: print(f"Debug: X-axis1: {_xmin}-{_xmax}")

    # add a barchart with relative changes
    #print(str(unique_aa_positions)) # these are one-based positions
    #print(f"Length of unique_aa_positions is {len(unique_aa_positions)}")
    if myoptions.aminoacids:
        total_frequencies = np.sum(np.abs(new_aa_table), axis=0)
    else:
        total_frequencies = np.sum(np.abs(new_codon_table), axis=0)

    if myoptions.debug: print(f"Debug: len(total_frequencies)={len(total_frequencies)}, total_frequencies={str(total_frequencies.to_list())}")

    # BUG: we cannot change X-axis xticks on one X-axis without resetting xmin and xmax for the twinx plot
    #     It makes no differenceif we use "ax2 = ax.twinx()" or if we use "fig, (ax, ax2) = plt.subplots(1, 2, sharex=True)"
    # https://github.com/matplotlib/matplotlib/issues/6860
    # https://github.com/matplotlib/matplotlib/pull/7904
    # https://github.com/matplotlib/matplotlib/pull/7528
    # https://github.com/chartjs/Chart.js/issues/3484
    if myoptions.aminoacids:
    #    ax.set_xticks(np.arange(0, max(unique_aa_positions), round(round(max(unique_aa_positions)/20)/10.0)*10)) # this breaks syncing of ax and ax2 axes to be out-of-sync, worked around via plt.locator_params()
        #ax.set_xticklabels(np.arange(0, max(unique_aa_positions)), rotation=90, ha='right')  # rotate X-axis legend
        # https://stackoverflow.com/questions/53747298/how-to-format-axis-tick-labels-from-number-to-thousands-or-millions-125-436-to#53747693
        #ax.xaxis.set_major_formatter(ticker.FuncFormatter(lambda x, pos: '{:,.0f}'.format(x/1))) # get rid of non-real numbers on the X-axis,like 2.5
        _ax1.xaxis.set_tick_params(labelsize=14)
        _ax1.tick_params(axis='x', which='both', labelsize=14)
        _y_ticks = np.arange(len(amino_acids))
        _ax1.set_yticks(_y_ticks)
        _ax1.set_yticklabels(amino_acids, fontsize=14)
    else:
    #    ax.set_xticks(np.arange(0, max(unique_codon_positions), round(round(max(unique_aa_positions)/20)/10.0)*10)) # this breaks syncing of ax and ax2 axes to be out-of-sync, worked around via plt.locator_params()
        #ax.set_xticklabels(np.arange(0, max(unique_codon_positions)), rotation=90, ha='right')  # rotate X-axis legend
        # https://stackoverflow.com/questions/53747298/how-to-format-axis-tick-labels-from-number-to-thousands-or-millions-125-436-to#53747693
        #ax.xaxis.set_major_formatter(ticker.FuncFormatter(lambda x, pos: '{:,.0f}'.format(x/1))) # get rid of non-real numbers on the X-axis,like 2.5
        _ax1.xaxis.set_tick_params(labelsize=14)
        _ax1.tick_params(axis='x', which='both', labelsize=14)
        _y_ticks = np.arange(len(codons_whitelist))
        _ax1.set_yticks(_y_ticks)
        _ax1.set_yticklabels([pairs[0] + ' (' + pairs[1] + ')  ' if pairs[1] not in ('INS', 'DEL') else pairs[0] + ' (' + pairs[1] + ')' for pairs in final_sorted_whitelist], fontsize=8) # sorted by amino acids
        _ax1.tick_params(axis='y', which='major', labelsize=8)

    #plt.xticks(rotation=90)
    # work around the bug with xmin being reset to zero when ax.set_xticks(() is used and adjust the spacing of ticks in a different way
    # https://www.geeksforgeeks.org/how-to-change-the-number-of-ticks-in-matplotlib/
    if myoptions.xaxis_bins:
        plt.locator_params(axis='x', nbins=myoptions.xaxis_bins )

    # add the grid in gray
    _ax1.grid(True, linestyle='--', alpha=0.3, color='gray')

    if not myoptions.disable_2nd_Y_axis:
        _ax2 = _ax1.twinx()
        _ax2.set_xlim(_xmin, _xmax) # start X-axis from 1, not zero
        ## Set the y-axis limit for the second axis (_ax2)
        _ax2.set_ylim(0, 1)

        x1, x2 = _ax2.get_xlim()
        _ax2.set_ylabel(f'Cumulative frequency of mutations above threshold {myoptions.threshold:.2%} per codon', fontsize=12)
        _ax1.figure.canvas.draw()
        _ax2.figure.canvas.draw()

        if myoptions.aminoacids:
            _ax2.bar(unique_aa_positions, total_frequencies, color='black', alpha=0.5, width=0.8, align='center') # , label='AA Frequencies')
        else:
            _ax2.bar(unique_codon_positions, total_frequencies, color='black', alpha=0.5, width=0.8, align='center') # , label='Codon Frequencies')

        x1, x2 = _ax2.get_xlim()
    
    if myoptions.aminoacids:
        _table = new_aa_table
    else:
        _table = new_codon_table
    if myoptions.debug:
        print(f"Debug: Printing the _table: {str(_table)}")
        # the _observed_codon_or_aa contains e.g. ANA, ANN, etc.

    if myoptions.debug:
        print("Debug: Printing the _table:")
        print(_table)

    #for row in _table.itertuples():
    #>>> for row in df.itertuples():
    #...    print(row)
    #Pandas(Index=4, A=0, B=2, C=3)

    #for row in _table.itertuples(index=False):
    #>>> for row in df.itertuples():
    #...    print(row)
    #Pandas(A=0, B=2, C=3)


    #for row in _table.items():
    #>>> for row in df.items():
    #...     print(row)
    #...
    #('A', 4     0
    #5     0
    #6    10
    #Name: A, dtype: int64)
    #('B', 4     2
    #5     4
    #6    20
    #Name: B, dtype: int64)
    #('C', 4     3
    #5     1
    #6    30
    #Name: C, dtype: int64)
    #>>>

    #for i, _some_codon_or_aa in enumerate(_table.index): # causes off-by-one error
    #    for j, _aa_position in enumerate(_table.columns): # causes off-by-one error

    _real_aa_positions = list(_table.columns)
    if myoptions.debug: print(f"Debug: _real_aa_positions={str(_real_aa_positions)}")
    #for row in _table.itertuples(index=True):
    #    print(f"Info: len(row)={len(row)}; row={str(row)}")
    #    _some_codon_or_aa = row['position']
    #    for _index, frequency in enumerate(row[1:]):

    if not myoptions.aminoacids:
        if list(_table.index) != codons_whitelist2: # ignore codons with 'N'
            # check the number of codons on Y-xis to match parsed data
            raise ValueError("Both lists should be equal: %s=%s, %s=%s" % ('_table.index', str(_table.index), 'codons_whitelist2', codons_whitelist2))
    else:
        if list(_table.index) != amino_acids: # check if the input were parsed correctly and codons with Ns dropped resulting in X amino acids
            # check the number of amino acid residues on Y-xis to match parsed data
            raise ValueError("Both lists should be equal: %s=%s, %s=%s" % ('_table.index', str(_table.index), 'amino_acids', amino_acids))

    # record the colors used for rendering in a TSV file
    if myoptions.aminoacids:
        _outfilename = _outfile_prefix + '.aa.frequencies.colors.tsv'
        _color_file = open(_outfilename, 'w')
        print(f"Info: Writing into {_outfilename}")
    else:
        _outfilename = _outfile_prefix + '.codon.frequencies.colors.tsv'
        _color_file = open(_outfilename, 'w')
        print(f"Info: Writing into {_outfilename}")

    _used_colors = set()
    _cmap, _colors = get_colormap(myoptions.colormap)
    _labels = [] # list of all simple label strings as they are added to a scatter plot axis
    _label_codon_positions = [] # position
    _label_original_amino_acids = [] # original_amino_acid
    _label_new_amino_acids = [] # new amino acid
    _label_cumulative_frequencies = [] # cumulative frequency
    _label_codon_frequencies = [] # codon frequencies
    _label_observed_codon_counts = [] # observed_codon_count
    _label_observed_codon_count_sum = [] # sum of the synonymous codon counts
    _label_total_codons_per_site = [] # total_codons_per_site
    _label_scores = [] # keep score for each position

    _html_labels = [] # they should be strings of unescaped HTML
    _mutations = []
    _points = [] # any of _circles, _markers, _dots
    _circles = [] # scaled by 100 for bokeh figures
    _circles5000 = [] # scaled by 5000 for matplotlib figures
    _markers = []
    _dots = []
    _warn_once = []
    _matrix_values = set() # unique matrix values collected for the data points
    for i, _some_codon_or_aa in enumerate(_table.index): # introduces off-by-one error, iterate over amino_acids which were used as the index using 'pd.DataFrame(Decimal(0), index=amino_acids, columns=unique_aa_positions)'
        for j, _aa_position in enumerate(_table.columns): # introduces off-by-one error
            if myoptions.debug: print(f"Debug: i: {str(i)}, j: {str(j)}, _aa_position column: {str(_aa_position)}")
            _real_aa_position = _real_aa_positions[j]
            if _real_aa_position != _aa_position:
                raise ValueError("Values _aa_position=%s and _real_aa_position=%s should be equal" % (_aa_position, _real_aa_position)) # both approaches should work
            frequency = _table.loc[_some_codon_or_aa, _aa_position]
            if myoptions.debug and frequency:
                # Debug0: _aa_position=351, _some_codon_or_aa=C, frequency=0.00109599999999999996903865540076594697893597185611724853515625
                # 7-SA-AA.SAPA@SA-Lib_F+SA-Lib_R.None+None.Wuhan.gofasta.frequencies.tsv
                # 351	Y	C	0.001096	TAT	TGT
                print(f"Debug0: _aa_position={_aa_position}, _some_codon_or_aa={_some_codon_or_aa}, frequency={str(frequency)}")
            try:
                _old_codon = df.loc[df['position'] == _aa_position]['original_codon'].to_list()[0] # pick any/first row to infer original_codon
            except IndexError:
                # do not die in regions where the are no data, for example in gaps between amplicons
                # H4/7-OM-EE.IOME@OM-Lib_F+OM-Lib_R.None+None.Wuhan.gofasta.frequencies.tsv at position 341
                # raise IndexError("Cannot determine original codon for position %s" % _aa_position)
                if frequency and myoptions.debug:
                    print(f"Debug0b: _aa_position={_aa_position}, _some_codon_or_aa={_some_codon_or_aa}, frequency={str(frequency)}")
                if _aa_position not in _warn_once:
                    sys.stderr.write("Warning: Cannot determine original codon for position %s, seems missing from input TSV or cannot split list %s%s" % (_aa_position, str(df.loc[df['position'] == _aa_position]['original_codon'].to_list()), os.linesep))
                    _warn_once.append(_aa_position)
                next
            if not np.abs(Decimal(frequency)) < myoptions.threshold:
                if myoptions.column_with_frequencies in ['neutralized_parent_difference', 'escape_parent_difference']:
                    size, _color = adjust_size_and_color_neutralized_escape(Decimal(frequency), _old_codon, _some_codon_or_aa, _matrix)
                    _score = -22 # TODO: junk value
                elif myoptions.column_with_frequencies in ['weighted_diff_escape_neutralized']:
                    size, _color = adjust_size_and_color_weighted(Decimal(frequency), _old_codon, _some_codon_or_aa, _matrix)
                    _score = -22 # TODO: junk value
                else:
                    #  we always pass-down the original codon, even when running in aa mode which confuses the code in adjust_size_and_color()
                    _score, size, _color = adjust_size_and_color(Decimal(frequency), _old_codon, _some_codon_or_aa, _matrix, _min_theoretical_score, _max_theoretical_score, _cmap) # the _color was derived via np.linspace mapping to range(0,1)
                    _matrix_values.add(_score) # update the set() with unique score values actually present in the dataset
                if myoptions.debug: print(f"Debug: Real AA position: {_aa_position}, observed codon: {_some_codon_or_aa}, frequency: {frequency}, size: {size}, color: {_color}")
                # for bokeh plots use _some_codon_or_aa index
                # TODO: try to emphasize differences by interpreting negative and positive matrix values, but these are later ignored anyway
                if myoptions.aminoacids:
                    if _score < 0:
                        _circles.append((_aa_position, _some_codon_or_aa, float(np.abs(size) * 100), 'circle_x', _color, 0.5, _score))
                    else:
                        _circles.append((_aa_position, _some_codon_or_aa, float(np.abs(size) * 100), 'circle_x', _color, 0.5, _score))
                else:
                    # print later codons followed by the amino acid they encode on Y-axis
                    _circles.append((_aa_position, _some_codon_or_aa + ' (' + alt_translate(_some_codon_or_aa) + ')', float(np.abs(size) * 100), 'circle_x', _color, 0.5))
                # for matplotlib figures use i index instead
                # TODO: try to emphasize differences by interpreting negative and positive matrix values, but these are later ignored anyway
                if _score < 0:
                    _circles5000.append((_aa_position, i, float(np.abs(size) * 5000), 'circle_x', _color, 0.5, _score)) 
                    _markers.append((_aa_position, i, 1, 'dot', 'black', 0.5)) # for matplotlib figures provide j, i pointers instead of _aa_position, _some_codon_or_aa for some reason
                else:
                    _circles5000.append((_aa_position, i, float(np.abs(size) * 5000), 'circle_x', _color, 0.5, _score))
                    _markers.append((_aa_position, i, 1, 'circle', 'black', 0.5)) # for matplotlib figures provide j, i pointers instead of _aa_position, _some_codon_or_aa for some reason
            else:
                # just draw some tiny dot otherwise pandas will drop empty Y-rows for unused amino acids or codons, which sucks and it btw does happen for charts with codons too although they have more data and supposedly are less likely to run into this issue but it does happen too
                size, _color = 0, 'black' # if we force the size to 0.00000000009 the invisible dots are drawn in matplotlib figs
                # for matplotlib figures use i index
                _dots.append((_aa_position, i, size, 'dot', _color, 0.5, _score))
                if myoptions.debug: print(f"Debug: Invisible dot. Real AA position: {_aa_position}, observed codon: {_some_codon_or_aa}, frequency: {frequency}, size: {size}, color: {_color}")
            # df.columns = ['position', 'original_aa', 'mutant_aa', 'frequency', 'original_codon', 'mutant_codon']


            # prepare labels for all datapoints
            # discard sequencing noise
            if _aa_position not in _warn_once:
                frequencies = [Decimal(x) for x in df.loc[(df['position'] == _aa_position) & (df['mutant_aa'] == _some_codon_or_aa) & (df[myoptions.column_with_frequencies] >= myoptions.threshold)][myoptions.column_with_frequencies].to_list()] # convert strings to Decimal values
                try:
                    old_amino_acid = df.loc[df['position'] == _aa_position]['original_aa'].to_list()[0] # this sometimes fails to slice 0 element
                except IndexError:
                    print(f"Error: Cannot slice {str(df.loc[df['position'] == _aa_position]['original_aa'].to_list())}")
                    old_amino_acid = df.loc[df['position'] == _aa_position]['original_aa'].to_list()[0]

                new_codons = df.loc[(df['position'] == _aa_position) & (df['mutant_aa'] == _some_codon_or_aa) & (df[myoptions.column_with_frequencies] >= myoptions.threshold)]['mutant_codon'].to_list()
                if myoptions.aminoacids:
                    if frequency != _table.at[_some_codon_or_aa, _aa_position]:
                        raise ValueError("Values frequency=%s and _table.at[_some_codon_or_aa, _aa_position]=%s should be equal" % (frequency, _table.at[_some_codon_or_aa, _aa_position])) # both approaches should work
                    if len(new_codons) != len(frequencies):
                        raise ValueError("len(new_codons) != len(frequencies), specifically: %s != %s" % (len(new_codons), len(frequencies)))
                    if 'observed_codon_count' in df.columns.values:
                        _observed_codon_counts = df.loc[(df['position'] == _aa_position) & (df['mutant_aa'] == _some_codon_or_aa) & (df[myoptions.column_with_frequencies] >= myoptions.threshold)]['observed_codon_count'].to_list()
                        _total_codons_per_site = df.loc[(df['position'] == _aa_position) & (df['mutant_aa'] == _some_codon_or_aa) & (df[myoptions.column_with_frequencies] >= myoptions.threshold)]['total_codons_per_site'].to_list() # sometimes this yields multiple exactly same numbers like [6863138, 6863138]
                        if _total_codons_per_site:
                            _total_codons_per_site = _total_codons_per_site[0] # anyway pick just the very first value, they are all same
                    else:
                        _observed_codon_counts = []
                        _total_codons_per_site = 0
                    _observed_codon_count_sum = sum(_observed_codon_counts)

                    # filtered_codons, filtered_frequencies, _observed_codon_count_sum = filter_codons(_some_codon_or_aa, new_codons, frequencies, _observed_codon_counts, myoptions.threshold)
                    if not frequency < myoptions.threshold: # append to the lists only if the frequency is above threshold
                        _labels.append(f"Position: {_aa_position}\nOriginal Amino Acid: {old_amino_acid} ({_old_codon})\nNew Amino Acid: {_some_codon_or_aa} {new_codons}\n{myoptions.matrix} score: {_score}\nCumulative Frequency: {sum(frequencies):.6f}\nCodon Frequencies: {['{:.6f}'.format(x) for x in frequencies]}\nObserved codon counts: {_observed_codon_counts}\nObserved codon counts sum: {_observed_codon_count_sum}\nTotal codons per site: {_total_codons_per_site}")
                        _label_codon_positions.append(f"{_aa_position}")
                        _label_original_amino_acids.append(f"{old_amino_acid} ({_old_codon})")
                        _label_new_amino_acids.append(f"{_some_codon_or_aa} {new_codons}")
                        _label_cumulative_frequencies.append(f"{sum(frequencies):.6f}")
                        _label_codon_frequencies.append(f"{['{:.6f}'.format(x) for x in frequencies]}")
                        if 'observed_codon_count' in df.columns.values:
                            _label_observed_codon_counts.append(f"{_observed_codon_counts}")
                            _label_observed_codon_count_sum.append(f"{sum(_observed_codon_counts)}")
                            _label_total_codons_per_site.append(f"{_total_codons_per_site}")
                        else:
                            _label_observed_codon_counts.append(_observed_codon_counts)
                            _label_observed_codon_count_sum.append(sum(_observed_codon_counts))
                            _label_total_codons_per_site.append(_total_codons_per_site)
                        _label_scores.append(_score)
                        _html_labels.append(f"Position: {_aa_position}<br>Original Amino Acid: {old_amino_acid} ({_old_codon})<br>New Amino Acid: {_some_codon_or_aa} {new_codons}<br>{myoptions.matrix} score: {_score}<br>Cumulative Frequency: {sum(frequencies):.6f}<br>Codon Frequencies: {['{:.6f}'.format(x) for x in frequencies]}<br>Observed codon counts: {_observed_codon_counts}<br>Observed codon count sum: {_observed_codon_count_sum}<br>Total codons per site: {_total_codons_per_site}")
                        _mutations.append(f"{old_amino_acid}{_aa_position}{_some_codon_or_aa}")
                else:
                    if frequency != _table.at[_some_codon_or_aa, _aa_position]:
                        raise ValueError("Values frequency=%s and _table.at[_some_codon_or_aa, _aa_position]=%s should be equal" % (frequency, _table.at[_some_codon_or_aa, _aa_position])) # both approaches should work

                    if 'observed_codon_count' in df.columns.values:
                        _observed_codon_counts = df.loc[(df['position'] == _aa_position) & (df['mutant_codon'] == _some_codon_or_aa) & (df[myoptions.column_with_frequencies] >= myoptions.threshold)]['observed_codon_count'].to_list()
                        if _observed_codon_counts and len(_observed_codon_counts) < 2:
                            _observed_codon_counts = _observed_codon_counts[0]
                            _observed_codon_count_sum = _observed_codon_counts
                        else:
                            _observed_codon_count_sum = sum(_observed_codon_counts)
                        _total_codons_per_site = df.loc[(df['position'] == _aa_position) & (df['mutant_codon'] == _some_codon_or_aa) & (df[myoptions.column_with_frequencies] >= myoptions.threshold)]['total_codons_per_site'].to_list()
                        if _total_codons_per_site and len(_total_codons_per_site) < 2:
                            _total_codons_per_site = _total_codons_per_site[0]
                    else:
                        _observed_codon_counts = []
                        _total_codons_per_site = 0

                    #old_amino_acid = df.loc[df['position'] == _aa_position]['original_aa'].to_list()[0]

                    if myoptions.debug:
                        # print relevant lines from df matching a particular codon column
                        print(f"Info: {len(df.loc[df['position'] == _aa_position]['position'])} aa residues observed in position {_aa_position}:{os.linesep} {str(df.loc[df['position'] == _aa_position][0:])}")
                    #print(f"Info: Looking for new_codon={_some_codon_or_aa}")
                    try:
                        new_amino_acid = df.loc[(df['position'] == _aa_position) & (df['mutant_codon'] == _some_codon_or_aa)]['mutant_aa'].to_list()[0]
                    except:
                        # new_amino_acid = "Failed to find, the %s is wrong and too far from original codon" % _some_codon_or_aa
                        new_amino_acid = None

                    if new_amino_acid:
                        try:
                            _frequency = Decimal(df.loc[(df['position'] == _aa_position) & (df['mutant_codon'] == _some_codon_or_aa)][myoptions.column_with_frequencies].to_list()[0]) # convert strings to Decimal values
                        except:
                            _frequency = 0.00000000009 # enter junk value but prevent matplotlib from crashing
                        if _frequency != frequency and (_frequency != 0.00000000009 and frequency != 0):
                            raise ValueError("Frequency new_codon_table.at[_some_codon_or_aa, _aa_position]=%s _some_codon_or_aa=%s, _aa_position=%s not same as df.loc[(df['position'] == _aa_position) & (df['mutant_codon'] == _some_codon_or_aa)][myoptions.column_with_frequencies].to_list()[0]=%s" % (frequency, _some_codon_or_aa, _aa_position, _frequency))

                        if old_amino_acid and not frequency < myoptions.threshold and size:
                            # print colors used only for data points above the threshold, because those below are not rendered
                            # also do not draw data points if the limits were applied manually with size=0 and palegreen color enforced

                            _label_scores.append(_score)

                            if myoptions.column_with_frequencies == 'neutralized_parent_difference':
                                _labels.append(f"Position: {_aa_position}\nOriginal Codon: {_old_codon} ({old_amino_acid})\nNew Codon: {_some_codon_or_aa} ({new_amino_acid})\n{myoptions.matrix} score: {_score}\nDifference neutralized2parent: {_frequency:.6f}")
                                _html_labels.append(f"Position: {_aa_position}<br>Original Codon: {_old_codon} ({old_amino_acid})<br>New Codon: {_some_codon_or_aa} ({new_amino_acid})<br>{myoptions.matrix} score: {_score}<br>Difference neutralized2parent: {_frequency:.6f}")
                            elif myoptions.column_with_frequencies == 'escape_parent_difference':
                                _labels.append(f"Position: {_aa_position}\nOriginal Codon: {_old_codon} ({old_amino_acid})\nNew Codon: {_some_codon_or_aa} ({new_amino_acid})\n{myoptions.matrix} score: {_score}\nDifference escape2parent: {_frequency:.6f}")
                                _html_labels.append(f"Position: {_aa_position}<br>Original Codon: {_old_codon} ({old_amino_acid})<br>New Codon: {_some_codon_or_aa} ({new_amino_acid})<br>{myoptions.matrix} score: {_score}<br>Difference escape2parent: {_frequency:.6f}")
                            elif myoptions.column_with_frequencies == 'weighted_diff_escape_neutralized':
                                _labels.append(f"Position: {_aa_position}\nOriginal Codon: {_old_codon} ({old_amino_acid})\nNew Codon: {_some_codon_or_aa} ({new_amino_acid})\n{myoptions.matrix} score: {_score}\nWeighted difference escape2neutralized: {_frequency:.6f}")
                                _html_labels.append(f"Position: {_aa_position}<br>Original Codon: {_old_codon} ({old_amino_acid})<br>New Codon: {_some_codon_or_aa} ({new_amino_acid})<br>{myoptions.matrix} score: {_score}<br>Weighted difference escape2neutralized: {_frequency:.6f}")
                            elif myoptions.column_with_frequencies == 'frequency':
                                _labels.append(f"Position: {_aa_position}\nOriginal Codon: {_old_codon} ({old_amino_acid})\nNew Codon: {_some_codon_or_aa} ({new_amino_acid})\n{myoptions.matrix} score: {_score}\nFrequency: {_frequency:.6f}\nObserved codon count: {_observed_codon_counts}\nTotal codons per site: {_total_codons_per_site}")
                                _html_labels.append(f"Position: {_aa_position}<br>Original Codon: {_old_codon} ({old_amino_acid})<br>New Codon: {_some_codon_or_aa} ({new_amino_acid})<br>{myoptions.matrix} score: {_score}<br>Frequency: {_frequency:.6f}<br>Observed codon counts: {_observed_codon_counts}<br>Total codons per site: {_total_codons_per_site}")
                                _label_codon_positions.append(f"{_aa_position}")
                                _label_original_amino_acids.append(f"{old_amino_acid} ({_old_codon})")
                                _label_new_amino_acids.append(f"{new_amino_acid} ({_some_codon_or_aa})")
                                _label_cumulative_frequencies.append(f"{_frequency:.6f}")
                                _label_codon_frequencies.append(f"{_frequency:.6f}")
                            else:
                                _labels.append(f"Position: {_aa_position}\nOriginal Codon: {_old_codon} ({old_amino_acid})\nNew Codon: {_some_codon_or_aa} ({new_amino_acid})\n{myoptions.matrix} score: {_score}\nFrequency: {_frequency:.6f}\nObserved codon count: {_observed_codon_counts}\nTotal codons per site: {_total_codons_per_site}")
                                _html_labels.append(f"Position: {_aa_position}<br>Original Codon: {_old_codon} ({old_amino_acid})<br>New Codon: {_some_codon_or_aa} ({new_amino_acid})<br>{myoptions.matrix} score: {_score}<br>Frequency: {_frequency:.6f}<br>Observed codon counts: {_observed_codon_counts}<br>Total codons per site: {_total_codons_per_site}")
                                _label_codon_positions.append(f"{_aa_position}")
                                _label_original_amino_acids.append(f"{old_amino_acid} ({_old_codon})")
                                _label_new_amino_acids.append(f"{new_amino_acid} ({_some_codon_or_aa})")
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
#                        else:
#                            # fill in the labels to get same amount of items in the list
#                            _frequency = Decimal(df.loc[(df['position'] == _aa_position) & (df['mutant_codon'] == _some_codon_or_aa)][myoptions.column_with_frequencies].to_list()[0])
#                            _labels.append(f"Position: {_aa_position}\nOriginal Codon: {_old_codon} ({old_amino_acid})\nNew Codon: {_some_codon_or_aa} ({new_amino_acid})\nFrequency: {_frequency:.6f}")
#                            _html_labels.append(f"Position: {_aa_position}<br>Original Codon: {_old_codon} ({old_amino_acid})<br>New Codon: {_some_codon_or_aa} ({new_amino_acid})<br>Frequency: {_frequency:.6f}")
#                    else:
#                        # fill in the labels to get same amount of items in the list
#                        try:
#                            _frequency = Decimal(df.loc[(df['position'] == _aa_position) & (df['mutant_codon'] == _some_codon_or_aa)][myoptions.column_with_frequencies].to_list()[0])
#                        except:
#                            _frequency = 0
#                        _labels.append(f"Position: {_aa_position}\nOriginal Codon: {_old_codon} ({old_amino_acid})\nNew Codon: {_some_codon_or_aa} ({new_amino_acid})\nFrequency: {_frequency}")
#                        _html_labels.append(f"Position: {_aa_position}<br>Original Codon: {_old_codon} ({old_amino_acid})<br>New Codon: {_some_codon_or_aa} ({new_amino_acid})<br>Frequency: {_frequency}")


                # write colors used while in amino acid mode
                position_in_protein = _aa_position
                _original_aas = df.loc[(df['position'] == position_in_protein)]['original_aa'].to_list()
                if _original_aas:
                    _original_aa = _original_aas[0]
                else:
                    _original_aa = None
                new_amino_acid = _some_codon_or_aa
                if myoptions.aminoacids:
                    _mutant_codons = df.loc[(df['position'] == position_in_protein) & (df['mutant_aa'] == new_amino_acid)]['mutant_codon'].to_list() # this returns empty list if not running in amino acid mode
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
                    # sanitize empty Data Frame for codon mode
                    _mutant_codons = df.loc[(df['position'] == position_in_protein) & (df['mutant_codon'] == _some_codon_or_aa)]['mutant_codon'].to_list()
                if _original_aa and not frequency < myoptions.threshold and size:
                    # print colors used only for data points above the threshold, because those below are not rendered
                    # also do not draw data points if the limits were applied manually with size=0 and palegreen color enforced
                    if len(_mutant_codons) > 1:
                        #  multiple codons encoding same mutated_aa do appear: TCT -> ['TCA', 'TCC', 'TCG']
                        _color_file.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s%s" % (_aa_position, _old_codon, str(_mutant_codons), _original_aa, _some_codon_or_aa, '{0:.6f}'.format(frequency), _color, _score, os.linesep))
                    elif not myoptions.aminoacids:
                        _color_file.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s%s" % (_aa_position, _old_codon, _mutant_codon, _original_aa, _mutant_aa, '{0:.6f}'.format(frequency), _color, _score, os.linesep))
                        # Note: beware there will be multiple lines per each mutated codon appearing. To sort the changes from worst to synonymous by the BLOSUM score use:
                        # sort -t$'\t' -k 1,1n -k 8,8n -k 6,6n prefix.gofasta.aa.frequencies.colors.tsv
                    elif _mutant_codons:
                        _color_file.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s%s" % (_aa_position, _old_codon, _mutant_codons[0], _original_aa, _some_codon_or_aa, '{0:.6f}'.format(frequency), _color, _score, os.linesep))
                    # else:
                    #     # do not draw an empty [] if the mutation did not happen and hence does not exist in the input TSV file
                elif myoptions.debug:
                    if myoptions.aminoacids:
                        sys.stderr.write("Debug: Skipped line for _aa_position=%s _original_aa=%s _old_codon=%s _mutant_codons=%s _some_codon_or_aa=%s frequency=%s score=%s myoptions.threshold=%s\n" % (_aa_position, _original_aa, _old_codon, str(_mutant_codons), _some_codon_or_aa, '{0:.6f}'.format(frequency), _score, myoptions.threshold))
                    else:
                        sys.stderr.write("Debug: Skipped line for _aa_position=%s _original_aa=%s _old_codon=%s _mutant_codons=%s _some_codon_or_aa=%s frequency=%s score=%s myoptions.threshold=%s\n" % (_aa_position, _original_aa, _old_codon, _mutant_codon, _mutant_aa, '{0:.6f}'.format(frequency), _score, myoptions.threshold))
                _used_colors.add(_color)
    _color_file.close()

    _half_size = int(max(abs(_min_theoretical_score), _max_theoretical_score))
    # cax=_ax1.inset_axes([- _half_size, - _half_size / 2, _half_size / 2, _half_size], transform=_ax1.transData)

    # https://www.bomberbot.com/python/exploring-the-power-and-versatility-of-matplotlibs-listedcolormap/
    #_figure.subplots_adjust(right=0.8)
    #_cbar_ax = _figure.add_axes([0.92, 0.1, 0.02, 0.8]) # [left, bottom, width, height]
    #_ax3.xaxis.set_major_locator(ticker.MultipleLocator(2))
    # _figure.tight_layout(h_pad=0)

    for label in _ax1.get_xticklabels():
        label.set_rotation(90)
        label.set_ha("center")


    print("Info: The following values were collected from matrix %s based on the actual data (some values from matrix might not be needed for your data, hence are not listed here): %s . Range spans %d values (before symmetrization)." % (myoptions.matrix, str(sorted(_matrix_values)), abs(min(_matrix_values)) + 1 + max(_matrix_values)))
    _merged_lists = _circles # _markers + _dots # do not add _markers as they are just the central dots inside circles
    if myoptions.debug:
        print(f"Debug: {len(_used_colors)} _used_colors used: {str(_used_colors)}")

    if myoptions.debug:
        print(f"Debug: {len(_labels)} _labels={str(_labels)}")
        print(f"Debug: {len(_html_labels)} _html_labels={str(_html_labels)}")
        print(f"Debug: {len(_circles)} _circles={str(_circles)}")
        print(f"Debug: {len(_markers)} _markers={str(_markers)}")
        print(f"Debug: {len(_dots)} _dots={str(_dots)}")


    # https://docs.bokeh.org/en/latest/docs/user_guide/interaction/tools.html#ug-interaction-tools-hover-tool
    mysource = bokeh.models.ColumnDataSource(data=dict(
        x=[x[0] for x in _circles],
        y=[y[1] for y in _circles],
        s=[s[2] for s in _circles], # size
        m=[m[3] for m in _circles], # marker
        c=[c[4] for c in _circles], # _color
        a=[a[5] for a in _circles], # alpha
        label=_labels,
        label1=_label_codon_positions,
        label2=_label_original_amino_acids,
        label3=_label_new_amino_acids,
        label4=_label_cumulative_frequencies,
        label5=_label_codon_frequencies,
        label6=_label_observed_codon_counts,
        label7=_label_observed_codon_count_sum,
        label8=_label_total_codons_per_site,
        label9=_label_scores,
        mutation=_mutations,
    ))

    if myoptions.aminoacids:
        TOOLTIPS = [
            #("index", "$index"),
            #("(x,y)", "($x, $y)"),
            # ("(snap_x,snap_y)", "($snap_x, $snap_y)"),
            #("size", "@s"),
            # ("marker", "@m"),
            # ("color", "@c"),
            # ("alpha", "@a"),
            # ("label", "@label"),
            ("Codon Position", "@label1"),
            ("Original Amino Acid", "@label2"),
            ("New Amino Acid", "@label3"),
            ("Cumulative Frequency", "@label4"),
            ("Codon Frequencies", "@label5"),
            ("Observed codon counts", "@label6"),
            ("Observed codon count sum", "@label7"),
            ("Total codons per site", "@label8"),
            ("Score", "@label9"),
        ]
    else:
        TOOLTIPS = [
            #("index", "$index"),
            #("(x,y)", "($x, $y)"),
            # ("(snap_x,snap_y)", "($snap_x, $snap_y)"),
            #("size", "@s"),
            # ("marker", "@m"),
            # ("color", "@c"),
            # ("alpha", "@a"),
            # ("label", "@label"),
            ("Codon Position", "@label1"),
            ("Original Amino Acid", "@label2"),
            ("New Amino Acid", "@label3"),
            ("Cumulative Frequency", "@label4"),
            ("Codon Frequencies", "@label5"),
            ("Observed codon count", "@label6"),
            ("Total codons per site", "@label8"),
            ("Score", "@label9"),
        ]
    if myoptions.aminoacids:
        p = bokeh.plotting.figure(x_range=(_xmin, _xmax), y_range=amino_acids, tooltips=TOOLTIPS, title=title_data, x_axis_label=_xlabel, y_axis_label='Introduced amino acid changes', x_minor_ticks=10, width=2000, height=1200, sizing_mode='stretch_width') # width=5000, height=5000)
    else:
        # codons sorted by amino acids
        p = bokeh.plotting.figure(x_range=(_xmin, _xmax), y_range=[pairs[0] + ' (' + pairs[1] + ')' for pairs in final_sorted_whitelist], tooltips=TOOLTIPS, title=title_data, x_axis_label=_xlabel, y_axis_label='Introduced codon changes', x_minor_ticks=10, height=1200, width=2000, sizing_mode='stretch_width') # width=5000, height=5000)

    # print(f"ticks={str(_ticks)}")
    # p.xaxis.ticker = _x_ticks
    #
    # axis_label, axis_label_align, axis_label_orientation, axis_label_standoff, axis_label_text_align, axis_label_text_alpha, axis_label_text_baseline, axis_label_text_color, axis_label_text_font, axis_label_text_font_size, axis_label_text_font_style, axis_label_text_line_height, axis_label_text_outline_color, axis_line_alpha, axis_line_cap, axis_line_color, axis_line_dash, axis_line_dash_offset, axis_line_join, axis_line_width, background_fill_alpha, background_fill_color, bounds, context_menu, coordinates, css_classes, css_variables, dimension, face, fixed_location, formatter, group, js_event_callbacks, js_property_callbacks, level, major_label_orientation, major_label_overrides, major_label_policy, major_label_standoff, major_label_text_align, major_label_text_alpha, major_label_text_baseline, major_label_text_color, major_label_text_font, major_label_text_font_size, major_label_text_font_style, major_label_text_line_height, major_label_text_outline_color, major_tick_in, major_tick_line_alpha, major_tick_line_cap, major_tick_line_color, major_tick_line_dash, major_tick_line_dash_offset, major_tick_line_join, major_tick_line_width, major_tick_out, minor_tick_in, minor_tick_line_alpha, minor_tick_line_cap, minor_tick_line_color, minor_tick_line_dash, minor_tick_line_dash_offset, minor_tick_line_join, minor_tick_line_width, minor_tick_out, name, propagate_hover, styles, stylesheets, subscribed_events, syncable, tags, ticker, visible, x_range_name or y_range_name
    p.xaxis.major_label_orientation = 'vertical'
    p.axis.axis_label_text_font_size = "12px"
    p.axis.major_label_text_font_size = "12px" # major_label_text_font_size, axis_label_text_font_size or major_label_text_font_style
    _scatter = p.scatter(x='x', y='y', size='s', marker='m', color='c', alpha='a', source=mysource)
    
    # https://github.com/bokeh/bokeh/issues/12660
    #TICKERS = [
    #bokeh.models.AdaptiveTicker(
    #    min_interval = 1,
    #    max_interval = 2
    #    )
    #]
    # p.xaxis.ticker = TICKERS # BUG: does not work

    # p.text() accepts:  anchor, angle, angle_units, background_fill_alpha, background_fill_color, background_hatch_alpha, background_hatch_color, background_hatch_extra, background_hatch_pattern, background_hatch_scale, background_hatch_weight, border_line_alpha, border_line_cap, border_line_color, border_line_dash, border_line_dash_offset, border_line_join, border_line_width, border_radius, decorations, js_event_callbacks, js_property_callbacks, name, outline_shape, padding, subscribed_events, syncable, tags, text, text_align, text_alpha, text_baseline, text_color, text_font, text_font_size, text_font_style, text_line_height, text_outline_color, x, x_offset, y or y_offset
    p.text(x='x', y='y', text='mutation', text_color='c', text_font_size='10px', text_font_style='bold', source=mysource)
    p.title.align = 'center'
    p.title.text_font_size = '14pt'
    print(f"Info: Writing into {_outfile_prefix + '.html'}")
    bokeh.plotting.output_file(_outfile_prefix + '.html')
    bokeh.plotting.show(p)

    
    #Some known testing position
    #index_position_in_TSV = df.loc[df['position'] == 145][0].to_list()[0]
    #old_amino_acid =  df.loc[df['position'] == 145][1].to_list()[0]
    #print(f"Original AA at position {index_position_in_TSV} was {old_amino_acid}")

    _norm = matplotlib.colors.BoundaryNorm(np.arange(-19, 19, 1), _cmap.N) # to account for BoundaryNorm choosing based on the lower bound
    _mpl_scatterplot = _ax1.scatter([x[0] for x in _circles5000], [x[1] for x in _circles5000], s=[x[2] for x in _circles5000], alpha=0.5, c=[x[6] for x in _circles5000], cmap=_cmap, norm=_norm)
    _colorbar = _figure.colorbar(_mpl_scatterplot, cax=_ax3, label="%s values" % myoptions.matrix, location='right', pad=-0.1, alpha=0.5)
    _colorbar.ax.set_yticks(np.arange(-18.5, 18.5, 1), np.arange(-18, 19, 1)) # colorbar is only -11 to +10, instead of -11 to +11
    _colorbar.ax.tick_params(axis='y', which='minor', length=0) # to clean up old ticks
    _ax1.scatter([x[0] for x in _markers], [x[1] for x in _markers], s=[x[2] for x in _markers], marker='x', color='black', alpha=0.5) # TODO: this applies same marker='x' to negative and positive values, somehow cannot find a way to utilize [x[3] for x in _markers] containing a mix of prepared dots and circles
    _ax1.scatter([x[0] for x in _dots], [x[1] for x in _dots], s=[x[2] for x in _dots], marker='.', color='black', alpha=0.5, cmap=_cmap)


    # display info on mouse hover()
    cursor = mplcursors.cursor(_ax1, hover=True)
    if myoptions.aminoacids:
        @cursor.connect("add")
        def on_add(sel):
            ypos, xpos = int(sel.target[1]), int(sel.target[0])
            #old_amino_acid = old_aa_table.index[ypos]
            new_amino_acid = new_aa_table.index[ypos]
            position_in_protein = new_aa_table.columns[xpos+1-myoptions.offset-_calculated_aa_offset] # undo the off-by-one error removal attempt on the ax.scatter line
            #print(f"Position: {position_in_protein}")
            frequency = new_aa_table.at[new_amino_acid, position_in_protein]
            frequencies = [Decimal(x) for x in df.loc[(df['position'] == position_in_protein) & (df['mutant_aa'] == new_amino_acid) & (df[myoptions.column_with_frequencies] >= myoptions.threshold)][myoptions.column_with_frequencies].to_list()] # convert strings to Decimal values
            old_amino_acid = df.loc[df['position'] == position_in_protein]['original_aa'].to_list()[0]
            old_codon = df.loc[df['position'] == position_in_protein]['original_codon'].to_list()[0]
            new_codons = df.loc[(df['position'] == position_in_protein) & (df['mutant_aa'] == new_amino_acid) & (df[myoptions.column_with_frequencies] >= myoptions.threshold)]['mutant_codon'].to_list()
            new_codon = new_codons[0]
            _observed_codon_counts = df.loc[(df['position'] == position_in_protein) & (df['mutant_aa'] == new_amino_acid) & (df[myoptions.column_with_frequencies] >= myoptions.threshold)]['observed_codon_count'].to_list() # fetch coverage counts for all synonymous codons
            _observed_codon_count_sum = sum(_observed_codon_counts)
            _total_codons_per_site = df.loc[(df['position'] == position_in_protein) & (df['mutant_codon'] == new_codon)]['total_codons_per_site'].to_list()
            if len(_total_codons_per_site):
                # multiple rows with same values must have been returned by pandas so pick just the first number
                _total_codons_per_site = _total_codons_per_site[0]
            _score = _matrix[old_amino_acid][new_amino_acid]
            # discard sequencing noise
            # filtered_codons, filtered_frequencies, _observed_codon_count_sum = filter_codons(_some_codon_or_aa, new_codons, frequencies, _observed_codon_counts, myoptions.threshold)
            if new_codon not in new_codons:
                raise ValueError("The new codon %s is not in the list of all codons %s encoding this aa %s" % (new_codon, str(new_codons), new_amino_acid))

            if myoptions.column_with_frequencies == 'neutralized_parent_difference':
                sel.annotation.set_text(f"Position: {position_in_protein}\nOriginal Amino Acid: {old_amino_acid} ({old_codon})\nNew Amino Acid: {new_amino_acid} ({new_codons})\n{myoptions.matrix} score: {_score}\nCumulative Frequency: {frequency:.6f}\nCodon Frequencies: {['{:.6f}'.format(x) for x in frequencies]}")
            elif myoptions.column_with_frequencies == 'escape_parent_difference':
                sel.annotation.set_text(f"Position: {position_in_protein}\nOriginal Amino Acid: {old_amino_acid} ({old_codon})\nNew Amino Acid: {new_amino_acid} ({new_codons})\n{myoptions.matrix} score: {_score}\nCumulative Frequency: {frequency:.6f}\nCodon Frequencies: {['{:.6f}'.format(x) for x in frequencies]}")
            elif myoptions.column_with_frequencies == 'weighted_diff_escape_neutralized':
                sel.annotation.set_text(f"Position: {position_in_protein}\nOriginal Amino Acid: {old_amino_acid} ({old_codon})\nNew Amino Acid: {new_amino_acid} ({new_codons})\n{myoptions.matrix} score: {_score}\nCumulative Frequency: {frequency:.6f}\nCodon Frequencies: {['{:.6f}'.format(x) for x in frequencies]}")
            elif myoptions.column_with_frequencies == 'frequency':
                sel.annotation.set_text(f"Position: {position_in_protein}\nOriginal Amino Acid: {old_amino_acid} ({old_codon})\nNew Amino Acid: {new_amino_acid} ({new_codons})\n{myoptions.matrix} score: {_score}\nCumulative Frequency: {frequency:.6f}\nCodon Frequencies: {['{:.6f}'.format(x) for x in frequencies]}\nObserved codon counts: {_observed_codon_counts}\nObserved codon count sum: {_observed_codon_count_sum}\nTotal codons per site: {_total_codons_per_site}")
            else:
                sel.annotation.set_text(f"Position: {position_in_protein}\nOriginal Amino Acid: {old_amino_acid} ({old_codon})\nNew Amino Acid: {new_amino_acid} ({new_codons})\n{myoptions.matrix} score: {_score}\nCumulative Frequency: {frequency:.6f}\nCodon Frequencies: {['{:.6f}'.format(x) for x in frequencies]}\nObserved codon counts: {_observed_codon_counts}\nObserved codon count sum: {_observed_codon_count_sum}\nTotal codons per site: {_total_codons_per_site}")
    else:
        @cursor.connect("add")
        def on_add(sel):
            ypos, xpos = int(sel.target[1]), int(sel.target[0])
            print(f"Info: xpos={xpos}, ypos={ypos}, _calculated_aa_offset={_calculated_aa_offset}")
            new_codon = new_codon_table.index[ypos]
            #_old_codon = old_codon_table.index[ypos]
            position_in_protein = new_codon_table.columns[xpos+1-myoptions.offset-_calculated_aa_offset]
            frequency = new_codon_table.at[new_codon, position_in_protein]
            print(f"Info: position_in_protein={position_in_protein}, frequency={frequency}")
            old_codon = df.loc[df['position'] == position_in_protein]['original_codon'].to_list()[0]
            old_amino_acid = df.loc[df['position'] == position_in_protein]['original_aa'].to_list()[0]
            #new_codon = df.loc[(df['position'] == position_in_protein)][xpos-myoptions.offset].to_list()['mutant_codon']
            _observed_codon_count = df.loc[(df['position'] == position_in_protein) & (df['mutant_codon'] == new_codon) & (df[myoptions.column_with_frequencies] >= myoptions.threshold)]['observed_codon_count'].to_list()[0]
            _observed_codon_count_sum = sum(_observed_codon_counts)
            _total_codons_per_site = df.loc[(df['position'] == position_in_protein) & (df['mutant_codon'] == new_codon) & (df[myoptions.column_with_frequencies] >= myoptions.threshold)]['total_codons_per_site'].to_list()[0]
            # print relevant lines from df matching a particular codon column
            print(f"Info: {len(df.loc[df['position'] == position_in_protein]['position'])} aa residues observed in position {position_in_protein}:{os.linesep} {str(df.loc[df['position'] == position_in_protein][0:])}")
            #print(f"Info: Looking for new_codon={new_codon}")
            try:
                new_amino_acid = df.loc[(df['position'] == position_in_protein) & (df['mutant_codon'] == new_codon)]['mutant_aa'].to_list()[0]
            except:
                new_amino_acid = "Failed to find, the %s is wrong and too far from original codon" % new_codon
            _score = _matrix[old_amino_acid][new_amino_acid]
            try:
                _frequency = Decimal(df.loc[(df['position'] == position_in_protein) & (df['mutant_codon'] == new_codon)][myoptions.column_with_frequencies].to_list()[0]) # convert strings to Decimal values
            except:
                _frequency = 0.00000000009 # enter junk value but prevent matplotlib from crashing
            if _frequency != frequency and (_frequency != 0.00000000009 and _frequency != 0) and not np.abs(frequency) < myoptions.threshold:
                raise ValueError("Frequency new_codon_table.at[new_codon, position_in_protein]=%s not same as df.loc[(df['position'] == position_in_protein) & (df['mutant_codon'] == new_codon)][myoptions.column_with_frequencies].to_list()[0]=%s" % (frequency, _frequency))

            if myoptions.debug: print(f"Debug: Position: {position_in_protein }\nOriginal Codon: {old_codon} ({old_amino_acid})\nNew Codon: {new_codon} ({new_amino_acid})\n{myoptions.matrix} score: {_score}\nFrequency: {frequency:.6f}")

            if myoptions.column_with_frequencies == 'neutralized_parent_difference':
                sel.annotation.set_text(f"Position: {position_in_protein}\nOriginal Codon: {old_codon} ({old_amino_acid})\nNew Codon: {new_codon} ({new_amino_acid})\n{myoptions.matrix} score: {_score}\nDifference neutralized2parent: {frequency:.6f}")
            elif myoptions.column_with_frequencies == 'escape_parent_difference':
                sel.annotation.set_text(f"Position: {position_in_protein}\nOriginal Codon: {old_codon} ({old_amino_acid})\nNew Codon: {new_codon} ({new_amino_acid})\n{myoptions.matrix} score: {_score}\nDifference escape2parent: {frequency:.6f}")
            elif myoptions.column_with_frequencies == 'weighted_diff_escape_neutralized':
                sel.annotation.set_text(f"Position: {position_in_protein}\nOriginal Codon: {old_codon} ({old_amino_acid})\nNew Codon: {new_codon} ({new_amino_acid})\n{myoptions.matrix} score: {_score}\nWeighted difference escape2neutralized: {frequency:.6f}")
            elif myoptions.column_with_frequencies == 'frequency':
                sel.annotation.set_text(f"Position: {position_in_protein}\nOriginal Codon: {old_codon} ({old_amino_acid})\nNew Codon: {new_codon} ({new_amino_acid})\n{myoptions.matrix} score: {_score}\nFrequency: {frequency:.6f}\nObserved codon count: {_observed_codon_count}\nTotal codons per site: {_total_codons_per_site}")
            else:
                sel.annotation.set_text(f"Position: {position_in_protein}\nOriginal Codon: {old_codon} ({old_amino_acid})\nNew Codon: {new_codon} ({new_amino_acid})\n{myoptions.matrix} score: {_score}\nFrequency: {frequency:.6f}\nObserved codon count: {_observed_codon_count}\nTotal codons per site: {_total_codons_per_site}")
            if myoptions.debug:
                print(f"Debug: final_sorted_whitelist={str(final_sorted_whitelist)}")
                print(f"Debug: codons_whitelist2={str(codons_whitelist2)}")

    # draw legend
    handles, labels = [], []
    if myoptions.aminoacids:
        _junk = 'NNN'
    else:
        _junk = 'NNN'
    for _freq in [0.001, 0.01, 0.1, 0.3]:
        if myoptions.column_with_frequencies in ['neutralized_parent_difference', 'escape_parent_difference']:
            _size, _color = adjust_size_and_color_neutralized_escape(Decimal(_freq), _junk, _junk, _matrix)
            # _size = _size * 5000
        elif myoptions.column_with_frequencies in ['weighted_diff_escape_neutralized']:
            _size, _color = adjust_size_and_color_weighted(Decimal(_freq), _junk, _junk, _matrix) # size is already multiplied by either 2000 or 3000 or 5000
        else:
            _score, _size, _color = adjust_size_and_color(Decimal(_freq), _junk, _junk, _matrix, _min_theoretical_score, _max_theoretical_score, _cmap)
            # _size = _size * 5000
        # print(f"Info: Freq is {_freq.__round__(3)}")
        handle = _ax2.scatter(_size, - 400 + _freq, s=float(_freq * 5000), color='gray', alpha=0.5, label=f'Frequency {_freq:.1%}')
        label = str(_freq)
        handles.append(handle)
        labels.append(label)
    # https://stackoverflow.com/questions/4700614/how-to-put-the-legend-outside-the-plot#4701285
    # Shrink current x-axis and y-axis by 50%
    #_box = _ax4.get_position()
    #_ax4.set_position([_box.x0, _box.y0, _box.width * 1.00, _box.height * 1.00])
    _ax4.set_axis_off()
    #plt.xticks(rotation=90)
    # 'center top' is not a valid value for loc; supported values are 'best', 'upper right', 'upper left', 'lower left', 'lower right', 'right', 'center left', 'center right', 'lower center', 'upper center', 'center'
    _ax2.legend(loc='upper center', bbox_to_anchor=(1.25, 1.00), labelspacing=3, frameon=False, handletextpad=1.5)

    for _ext in ('.png', '.pdf'): # '.jpg'
        _wholefig = plt.gcf()
        _figsize= _wholefig.get_size_inches()*_wholefig.dpi # size in pixels
        print("Info: Writing into %s, figure size is %s inches and %s dpi" % (_outfile_prefix + _ext, _wholefig.get_size_inches(), _figsize))
        _figure.savefig(_outfile_prefix + _ext, dpi=myoptions.dpi)
    plt.show()
    plt.clf()


#    # draw legends in a separate figure, just in case
#    _figure_legend = plt.figure()
#    _ax5 = plt.subplot(111)
#    handles, labels = [], []
#    if myoptions.aminoacids:
#        _junk = 'NNN'
#    else:
#        _junk = 'NNN'
#    for _freq in [0.001, 0.01, 0.1, 0.3]:
#        if myoptions.column_with_frequencies in ['neutralized_parent_difference', 'escape_parent_difference']:
#            _size, _color = adjust_size_and_color_neutralized_escape(Decimal(_freq), _junk, _junk, _matrix)
#            # _size = _size * 5000
#        elif myoptions.column_with_frequencies in ['weighted_diff_escape_neutralized']:
#            _size, _color = adjust_size_and_color_weighted(Decimal(_freq), _junk, _junk, _matrix) # size is already multiplied by either 2000 or 3000 or 5000
#        else:
#            _score, _size, _color = adjust_size_and_color(Decimal(_freq), _junk, _junk, _matrix, _min_theoretical_score, _max_theoretical_score, _cmap)
#            # _size = _size * 5000
#        # print(f"Info: Freq is {_freq.__round__(3)}")
#        handle = _ax5.scatter(400 + _freq, _size, s=float(_freq * 5000), color='magenta', alpha=0.6, label=f'Frequency {_freq:.3f}')
#        label = str(_freq)
#        handles.append(handle)
#        labels.append(label)
#    # https://stackoverflow.com/questions/4700614/how-to-put-the-legend-outside-the-plot#4701285
#    # Shrink current x-axis by 70%
#    _box = _ax5.get_position()
#    _ax5.set_position([_box.x0, _box.y0, _box.width * 0.7, _box.height])
#    plt.xticks(rotation=90)
#
#    # https://matplotlib.org/stable/api/_as_gen/matplotlib.figure.Figure.legend.html
#    # https://matplotlib.org/stable/users/explain/axes/legend_guide.html#legend-guide
#    # print(str(dir(_figure_legend)))
#    #_figure_legend.legend(bbox_to_anchor=(0.15, 0.15), loc='lower left') #, bbox_transform=_figure_legend.transFigure)
#    # supported values are 'best', 'upper right', 'upper left', 'lower left', 'lower right', 'right', 'center left', 'center right', 'lower center', 'upper center', 'center'
#    _ax5.legend(loc='center left', bbox_to_anchor=(1.05, 0.5), labelspacing=3, frameon=False, handletextpad=1.5)
#    print(f"Info: Writing into {_outfile_prefix + '.legend.png'}")
#    _figure_legend.savefig(_outfile_prefix + '.legend.png', dpi=myoptions.dpi)



#    if not myoptions.colorbar:
#        plt.clf()
#        _figure_colorbar = plt.figure()
#        _ax6 = plt.subplot(111)
#        handles, labels = [], []
#        if myoptions.aminoacids:
#            _junk = 'NNN'
#        else:
#            _junk = 'NNN'
#        for _freq in [0.001, 0.01, 0.1, 0.3]:
#            if myoptions.column_with_frequencies in ['neutralized_parent_difference', 'escape_parent_difference']:
#                _size, _color = adjust_size_and_color_neutralized_escape(Decimal(_freq), _junk, _junk, _matrix)
#                # _size = _size * 5000
#            elif myoptions.column_with_frequencies in ['weighted_diff_escape_neutralized']:
#                _size, _color = adjust_size_and_color_weighted(Decimal(_freq), _junk, _junk, _matrix) # size is already multiplied by either 2000 or 3000 or 5000
#            else:
#                _score, _size, _color = adjust_size_and_color(Decimal(_freq), _junk, _junk, _matrix, _min_theoretical_score, _max_theoretical_score, _cmap)
#                # _size = _size * 5000
#            # print(f"Info: Freq is {_freq.__round__(3)}")
#            handle = _ax6.scatter(400 + _freq, _size, s=float(_freq * 5000), color='magenta', alpha=0.6, label=f'Frequency {_freq:.3f}')
#            label = str(_freq)
#            handles.append(handle)
#            labels.append(label)
#        # https://stackoverflow.com/questions/4700614/how-to-put-the-legend-outside-the-plot#4701285
#        # Shrink current x-axis by 70%
#        _box = _ax6.get_position()
#        _ax6.set_position([_box.x0, _box.y0, _box.width * 0.7, _box.height])
#        plt.xticks(rotation=90)
#
#        _colorbar = _figure_colorbar.colorbar(plt.cm.ScalarMappable(norm=matplotlib.colors.Normalize(- _half_size, _half_size, _half_size * 2 + 1), cmap=_cmap), ax=_ax6, label="%s values" % myoptions.matrix, location='left', pad=0.15, alpha=0.5)
#        print(f"Info: Writing into {_outfile_prefix + '.colorbar.png'}")
#        _figure_colorbar.savefig(_outfile_prefix + '.colorbar.png', dpi=myoptions.dpi)
#        print(f"Info: Writing into {_outfile_prefix + '.colorbar.pdf'}")
#        _figure_colorbar.savefig(_outfile_prefix + '.colorbar.pdf', dpi=myoptions.dpi)
#        plt.clf()


if __name__ == "__main__":
    main()

# vim:ts=4:sw=4:expandtab:smartindent
