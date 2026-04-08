# -*- coding: utf-8 -*-
# vim: set fileencoding=utf-8 ts=4 sw=4 expandtab :
"""
Tests verifying the dynamic and static-virtual colormap bounding constraint logic.
"""
import os
import unittest

import matplotlib
import matplotlib.pyplot as plt

# Import core functions
# pylint: disable=import-error
from mutation_scatter_plot.mutation_scatter_plot import (
    build_frequency_tables, collect_scatter_data, load_and_clean_dataframe,
    load_matrix, setup_matplotlib_figure)


class TestDynamicColormap(unittest.TestCase):
    """Verify that continuous colormaps are correctly dynamically mapped or forcibly unbound."""

    def setUp(self):
        matplotlib.use('Agg')
        self.tests_dir = os.path.dirname(os.path.abspath(__file__))
        self.outputs_dir = os.path.join(self.tests_dir, "outputs")
        self.tsv_path = os.path.join(self.outputs_dir, "test2_full.x_after_count.frequencies.tsv")

        if not os.path.exists(self.tsv_path):
            self.skipTest(f"Input file {self.tsv_path} not found.")

    def _get_mock_options(self, dynamic=True):
        class Options:
            # pylint: disable=too-few-public-methods
            """Mock options class for CLI parameters."""

            def __init__(self):
                self.cmap_vmin = None
                self.cmap_vmax = None
                self.tsv_file_path = None
                self.matrix = 'BLOSUM80'
                self.matrix_file = None
                self.outfile_prefix = 'scaling_test'
                self.offset = 0
                self.column_with_frequencies = 'frequency'
                self.threshold = 0.001
                self.aminoacids = True
                self.showstop = False
                self.showins = False
                self.showdel = False
                self.showx = False
                self.debug = 0
                self.xaxis_major_ticks_spacing = 10
                self.xaxis_minor_ticks_spacing = 5
                self.xaxis_label_start = 0
                self.xaxis_bins = 0
                self.shortlegend = True
                self.title = 'Colormap Test'
                self.disable_2nd_Y_axis = False
                self.include_synonymous = False
                self.colormap = 'coolwarm_r'
                self.dpi = 100
                self.xmin = 0
                self.xmax = 0
                self.disable_showing_bokeh = True
                self.disable_showing_mplcursors = True
                self.bokeh_sqrt_size = True
                self.linear_circle_size = False
                self.disable_padded_x_axis = False

                # Activate flag if dynamic scaling is disabled
                self.spread_colormap_virtual_matrix = not dynamic

        opts = Options()
        opts.tsv_file_path = self.tsv_path
        return opts

    def test_colormap_bounds_divergence(self):
        """Compare the resulting palettes from dynamic dataset bounding vs static theoretical bounding."""

        # 1. Dynamic Mode
        options_dyn = self._get_mock_options(dynamic=True)
        matrix, matrix_name, min_theoretical, max_theoretical, outfile_prefix = load_matrix(options_dyn)
        df, p2p = load_and_clean_dataframe(options_dyn, options_dyn.tsv_file_path, {})
        (aa, cw, cw2, fsw, uap, ucp, _, nat, _, nct, _, p2p) = \
            build_frequency_tables(options_dyn, df, p2p)

        fig, ax1, ax2, ax3, ax4, xmin, xmax = setup_matplotlib_figure(
            options_dyn, 'Test', '0', matrix_name, aa, cw, fsw, uap, ucp, nat, nct, p2p)

        _, _, _, _, _, _, _, circles_matplotlib_dyn, _, _, _ = \
            collect_scatter_data(options_dyn, df, nat, outfile_prefix, matrix, aa, cw2, p2p, xmin, xmax)

        dyn_vmin = options_dyn.cmap_vmin
        dyn_vmax = options_dyn.cmap_vmax
        plt.close(fig)

        # 2. Virtual Static Mode
        options_stat = self._get_mock_options(dynamic=False)
        _, _, _, _, _ = load_matrix(options_stat)
        fig, ax1, ax2, ax3, ax4, xmin, xmax = setup_matplotlib_figure(
            options_stat, 'Test', '0', matrix_name, aa, cw, fsw, uap, ucp, nat, nct, p2p)

        _, _, _, _, _, _, _, circles_matplotlib_stat, _, _, _ = \
            collect_scatter_data(options_stat, df, nat, outfile_prefix, matrix, aa, cw2, p2p, xmin, xmax)

        stat_vmin = options_stat.cmap_vmin
        stat_vmax = options_stat.cmap_vmax
        plt.close(fig)

        # Assert physical boundary constraints are perfectly symmetric extremes
        bound_abs = max(abs(min_theoretical), abs(max_theoretical))
        self.assertEqual(stat_vmin, -bound_abs)
        self.assertEqual(stat_vmax, bound_abs)

        self.assertGreaterEqual(dyn_vmin, -bound_abs)
        self.assertLessEqual(dyn_vmax, bound_abs)

        # Compare color scaling for a specific observed score map
        # If the dataset doesn't strictly span the exact theoretical boundaries, the hex palettes
        # will wildly diverge for the same score.
        score_to_color_dyn = {c[6]: c[4] for c in circles_matplotlib_dyn}
        score_to_color_stat = {c[6]: c[4] for c in circles_matplotlib_stat}

        # Save explicitly for reviewing artifact properties
        import json
        with open(os.path.join(self.outputs_dir, "test_colormap_bounds.dyn.colors.json"), "w") as f:
            json.dump(score_to_color_dyn, f, indent=4)
        with open(os.path.join(self.outputs_dir, "test_colormap_bounds.stat.colors.json"), "w") as f:
            json.dump(score_to_color_stat, f, indent=4)

        diverged_scores = 0
        for score in score_to_color_dyn:
            if score == 12:  # Skip the sentinel dark-green
                continue
            color_dyn = score_to_color_dyn[score]
            color_stat = score_to_color_stat[score]
            if color_dyn != color_stat:
                diverged_scores += 1

        # We expect color values to diverge if dyn_vmin != stat_vmin or dyn_vmax != stat_vmax
        if dyn_vmin != stat_vmin or dyn_vmax != stat_vmax:
            self.assertGreater(diverged_scores, 0,
                               "Color mappings identically collided despite distinct spread boundaries!")


if __name__ == "__main__":
    unittest.main()
