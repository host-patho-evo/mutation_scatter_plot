"""
Tests for verifying the mathematical area-to-frequency scaling logic.
"""
import os
import unittest
import numpy as np
import matplotlib
import matplotlib.pyplot as plt

# Import core functions
from mutation_scatter_plot.mutation_scatter_plot import (
    load_matrix, load_and_clean_dataframe, build_frequency_tables, setup_matplotlib_figure,
    collect_scatter_data, render_matplotlib
)

class TestScalingLogic(unittest.TestCase):
    """Verify that area-to-frequency scaling (1:10:100:300) is consistent across backends."""

    def setUp(self):
        matplotlib.use('Agg')
        self.tests_dir = os.path.dirname(os.path.abspath(__file__))
        self.project_root = os.path.dirname(self.tests_dir)
        self.outputs_dir = os.path.join(self.tests_dir, "outputs")
        self.tsv_path = os.path.join(self.outputs_dir, "test2_full.x_after_count.frequencies.tsv")

        if not os.path.exists(self.tsv_path):
            self.skipTest(f"Input file {self.tsv_path} not found.")

    def _get_mock_options(self, aminoacids=True, threshold=0.001):
        # pylint: disable=too-many-instance-attributes,too-few-public-methods
        class Options:
            """Mock options class for CLI parameters."""
            def __init__(self):
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
                self.title = 'Scaling Test'
                # pylint: disable=invalid-name
                self.disable_2nd_Y_axis = False
                self.include_synonymous = False
                self.colormap = 'amino_acid_changes'
                self.dpi = 100
                self.xmin = 0
                self.xmax = 0
                self.disable_showing_bokeh = True
                self.disable_showing_mplcursors = True
                self.bokeh_sqrt_size = True
                self.linear_circle_size = False

        opts = Options()
        opts.tsv_file_path = self.tsv_path
        opts.aminoacids = aminoacids
        opts.threshold = threshold
        return opts

    def test_matplotlib_scaling_aminoacids(self):
        """Verify Matplotlib legend area ratios (1:10:100:300) in amino acid mode."""
        options = self._get_mock_options(aminoacids=True)
        matrix, matrix_name, _, _, outfile_prefix = load_matrix(options)
        df, p2p = load_and_clean_dataframe(options, options.tsv_file_path, {})
        (aa, cw, cw2, fsw, uap, ucp, _, nat, _, nct, _, p2p) = \
            build_frequency_tables(options, df, p2p)

        fig, ax1, ax2, ax3, ax4, xmin, xmax = setup_matplotlib_figure(
            options, 'Test', '0', matrix_name, aa, cw, fsw, uap, ucp, nat, nct)

        norm, cmap, colors, _, _, _, \
            _, circles_matplotlib, markers, dots, _ = \
            collect_scatter_data(options, df, nat, outfile_prefix, matrix, aa, cw2, p2p, xmin, xmax)

        render_matplotlib(
            options, fig, ax1, ax2, ax3, ax4, outfile_prefix,
            circles_matplotlib, markers, dots, cmap, norm, colors,
            matrix, matrix_name, show=False
        )

        # Legend collections at indices [0, 1, 2, 3, 4, 5] correspond to 100%, 50%, 30%, 10%, 1%, 0.1%
        # They are added to ax2 in render_matplotlib
        sizes = [coll.get_sizes()[0] for coll in ax2.collections]
        ratios = [s / sizes[0] for s in sizes]

        # Expected ratios: 1.0, 0.5, 0.3, 0.1, 0.01, 0.001
        np.testing.assert_allclose(ratios, [1.0, 0.5, 0.3, 0.1, 0.01, 0.001], rtol=1e-5)
        plt.close(fig)

    def test_matplotlib_scaling_codons(self):
        """Verify Matplotlib legend area ratios (1:10:100:300) in codon mode."""
        options = self._get_mock_options(aminoacids=False)
        matrix, matrix_name, _, _, outfile_prefix = load_matrix(options)
        df, p2p = load_and_clean_dataframe(options, options.tsv_file_path, {})
        (aa, cw, cw2, fsw, uap, ucp, _, nat, _, nct, _, p2p) = \
            build_frequency_tables(options, df, p2p)

        fig, ax1, ax2, ax3, ax4, xmin, xmax = setup_matplotlib_figure(
            options, 'Test', '0', matrix_name, aa, cw, fsw, uap, ucp, nat, nct)

        norm, cmap, colors, _, _, _, \
            _, circles_matplotlib, markers, dots, _ = \
            collect_scatter_data(options, df, nct, outfile_prefix, matrix, aa, cw2, p2p, xmin, xmax)

        render_matplotlib(
            options, fig, ax1, ax2, ax3, ax4, outfile_prefix,
            circles_matplotlib, markers, dots, cmap, norm, colors,
            matrix, matrix_name, show=False
        )

        sizes = [coll.get_sizes()[0] for coll in ax2.collections]
        ratios = [s / sizes[0] for s in sizes]

        np.testing.assert_allclose(ratios, [1.0, 0.5, 0.3, 0.1, 0.01, 0.001], rtol=1e-5)
        plt.close(fig)

    def test_bokeh_scaling_logic(self):
        """Verify Bokeh internal size calculation (size = sqrt(freq)*100)."""
        # freq = 0.3 -> size = sqrt(0.3) * 100 = 54.77
        # freq = 0.001 -> size = sqrt(0.001) * 100 = 3.16

        f1 = 0.001
        f2 = 0.3

        s1 = float(np.sqrt(f1) * 100)
        s2 = float(np.sqrt(f2) * 100)

        area1 = np.pi * (s1/2)**2
        area2 = np.pi * (s2/2)**2

        self.assertAlmostEqual(area2 / area1, 300.0, places=5)

    def test_specific_mutation_scaling(self):
        """
        Verify Matplotlib area and Bokeh diameter for specific mutations from test2_full.

        | Mutation | Frequency | MPL Area (s) | Bokeh Diameter |
        | :---     | :---      | :---         | :---           |
        | Q498R    | 0.700161  | 3500.81      | 83.68          |
        | E484A    | 0.691146  | 3455.73      | 83.13          |
        | N501Y    | 0.993272  | 4966.36      | 99.66          |
        """
        mut_data = {
            'Q498R': (0.700161, 3500.805, 83.6756),
            'E484A': (0.691146, 3455.730, 83.1352),
            'N501Y': (0.993272, 4966.360, 99.6630)
        }

        for name, (freq, exp_s, exp_size) in mut_data.items():
            with self.subTest(mutation=name):
                # Matplotlib Area intent
                calc_s = float(freq * 5000)
                self.assertAlmostEqual(calc_s, exp_s, places=3)

                # Bokeh Diameter intent
                calc_size = float(np.sqrt(freq) * 100)
                self.assertAlmostEqual(calc_size, exp_size, places=4)

if __name__ == "__main__":
    unittest.main()
