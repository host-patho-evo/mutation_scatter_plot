import os
import sys
sys.path.insert(0, os.path.abspath('src'))
from mutation_scatter_plot.mutation_scatter_plot import get_colormap, adjust_size_and_color

class Options:
    def __init__(self):
        self.colormap = 'amino_acid_changes'
        self.aminoacids = True
        self.spread_colormap_virtual_matrix = False
        self.matrix_min_theoretical = -6
        self.matrix_max_theoretical = 11
        self.cmap_actual_vmin = -4
        self.cmap_actual_vmax = 4
        self.debug = False
        self.synonymous_size = 1000
        self.synonymous_size2 = 1000
        self.linear_circle_size = False
        self.matrix = 'BLOSUM80'

myoptions = Options()

_norm, _cmap, _colors = get_colormap(myoptions, 'amino_acid_changes')

print(f"Colors length: {len(_colors)}")
print(f"Colors array: {_colors}")

matrix = {'V': {'H': -4}}
print("Testing adjust_size_and_color for score -4:")
_score, freq, _hex = adjust_size_and_color(myoptions, 0.5, False, "V", "H", "V", "H", matrix, _norm, _colors)
print(f"Returned Hex: {_hex}")

