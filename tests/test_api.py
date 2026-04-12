# -*- coding: utf-8 -*-
# vim: set fileencoding=utf-8 ts=4 sw=4 expandtab :
"""Tests for the high-level programmatic API (api.py and options.py)."""

import os
import tempfile
import unittest

import matplotlib
matplotlib.use('Agg')


class TestScatterOptions(unittest.TestCase):
    """Test the scatter_options factory function."""

    def test_defaults(self):
        """scatter_options() returns sensible defaults matching CLI."""
        from mutation_scatter_plot.mutation_scatter_plot.options import scatter_options
        opts = scatter_options(tsv_file_path='test.tsv', outfile_prefix='out')
        self.assertEqual(opts.matrix, 'BLOSUM80')
        self.assertEqual(opts.colormap, 'amino_acid_changes')
        self.assertEqual(opts.threshold, 0.001)
        self.assertEqual(opts.dpi, 600)
        self.assertFalse(opts.aminoacids)
        self.assertFalse(opts.include_synonymous)
        self.assertFalse(opts.linear_circle_size)

    def test_overrides(self):
        """scatter_options() applies keyword overrides correctly."""
        from mutation_scatter_plot.mutation_scatter_plot.options import scatter_options
        opts = scatter_options(
            tsv_file_path='test.tsv', outfile_prefix='out',
            aminoacids=True, xmin=430, xmax=528,
            matrix='BLOSUM62', threshold=0.01,
        )
        self.assertTrue(opts.aminoacids)
        self.assertEqual(opts.xmin, 430)
        self.assertEqual(opts.xmax, 528)
        self.assertEqual(opts.matrix, 'BLOSUM62')
        self.assertEqual(opts.threshold, 0.01)

    def test_unknown_key_raises(self):
        """scatter_options() raises TypeError on unknown keys."""
        from mutation_scatter_plot.mutation_scatter_plot.options import scatter_options
        with self.assertRaises(TypeError) as ctx:
            scatter_options(tsv_file_path='', outfile_prefix='', nonexistent=True)
        self.assertIn('nonexistent', str(ctx.exception))


class TestTimelineOptions(unittest.TestCase):
    """Test the timeline_options factory function."""

    def test_defaults(self):
        """timeline_options() returns sensible defaults."""
        from mutation_scatter_plot.mutation_scatter_plot.options import timeline_options
        opts = timeline_options(directory='/tmp', positions=['501'])
        self.assertEqual(opts.matrix, 'BLOSUM80')
        self.assertEqual(opts.threshold, 0.0)


class TestFrequencyOptions(unittest.TestCase):
    """Test the frequency_options factory function."""

    def test_defaults(self):
        """frequency_options() returns sensible defaults."""
        from mutation_scatter_plot.mutation_scatter_plot.options import frequency_options
        opts = frequency_options(
            alignment_infilename='aln.fasta',
            reference_infilename='ref.fasta',
            outfileprefix='out',
        )
        self.assertEqual(opts.translation_table, 1)
        self.assertFalse(opts.x_after_count)
        self.assertTrue(opts.print_unchanged_sites)


class TestLazyReExports(unittest.TestCase):
    """Test that lazy re-exports from package __init__.py work."""

    def test_scatter_options_importable(self):
        """scatter_options is importable from mutation_scatter_plot."""
        from mutation_scatter_plot import scatter_options
        opts = scatter_options(tsv_file_path='t.tsv', outfile_prefix='o')
        self.assertEqual(opts.matrix, 'BLOSUM80')

    def test_render_scatter_importable(self):
        """render_scatter is importable from mutation_scatter_plot."""
        from mutation_scatter_plot import render_scatter
        self.assertTrue(callable(render_scatter))

    def test_calculate_frequencies_importable(self):
        """calculate_frequencies is importable from mutation_scatter_plot."""
        from mutation_scatter_plot import calculate_frequencies
        self.assertTrue(callable(calculate_frequencies))

    def test_unknown_attr_raises(self):
        """Accessing a non-existent attribute from the package raises AttributeError."""
        import mutation_scatter_plot
        with self.assertRaises(AttributeError):
            _ = mutation_scatter_plot.nonexistent_function


class TestRenderScatterAPI(unittest.TestCase):
    """Integration test: render_scatter() produces correct outputs."""

    def setUp(self):
        """Set up paths to test fixtures."""
        self.tests_dir = os.path.dirname(os.path.abspath(__file__))
        self.outputs_dir = os.path.join(self.tests_dir, "outputs")
        self.tsv_input = os.path.join(
            self.outputs_dir, "test1.default.frequencies.tsv"
        )

    def test_render_scatter_aminoacids(self):
        """render_scatter() in aminoacids mode produces Figure and files."""
        from mutation_scatter_plot.api import render_scatter

        with tempfile.TemporaryDirectory() as tmpdir:
            outfile_prefix = os.path.join(tmpdir, "api_test")
            result = render_scatter(
                tsv_path=self.tsv_input,
                outfile_prefix=outfile_prefix,
                aminoacids=True,
                threshold=0.001,
                show_bokeh=False,
                show_mplcursors=False,
            )

            # Check return structure
            self.assertIn('figure', result)
            self.assertIn('axes', result)
            self.assertIn('dataframe', result)
            self.assertIn('files_written', result)

            # Figure is a matplotlib Figure
            self.assertIsNotNone(result['figure'])
            self.assertEqual(len(result['axes']), 4)

            # Files were written
            self.assertGreater(len(result['files_written']), 0)

            # PNG and HTML exist
            extensions = {os.path.splitext(f)[1] for f in result['files_written']}
            self.assertIn('.png', extensions)
            self.assertIn('.html', extensions)

            # DataFrame has expected columns
            df = result['dataframe']
            self.assertIn('padded_position', df.columns)
            self.assertIn('frequency', df.columns)

    def test_render_scatter_codons_coolwarm(self):
        """render_scatter() in codon mode with coolwarm_r colormap."""
        from mutation_scatter_plot.api import render_scatter

        with tempfile.TemporaryDirectory() as tmpdir:
            outfile_prefix = os.path.join(tmpdir, "api_test_codons")
            result = render_scatter(
                tsv_path=self.tsv_input,
                outfile_prefix=outfile_prefix,
                aminoacids=False,
                colormap='coolwarm_r',
                threshold=0.001,
                show_bokeh=False,
                show_mplcursors=False,
            )
            self.assertIsNotNone(result['figure'])
            self.assertGreater(len(result['files_written']), 0)

    def test_render_scatter_no_outfile(self):
        """render_scatter() without outfile_prefix returns Figure without writing files."""
        from mutation_scatter_plot.api import render_scatter

        result = render_scatter(
            tsv_path=self.tsv_input,
            aminoacids=True,
            threshold=0.001,
        )
        self.assertIsNotNone(result['figure'])
        self.assertEqual(result['files_written'], [])

    def test_render_scatter_linear_scaling(self):
        """render_scatter() with linear_scaling=True works."""
        from mutation_scatter_plot.api import render_scatter

        with tempfile.TemporaryDirectory() as tmpdir:
            outfile_prefix = os.path.join(tmpdir, "api_test_linear")
            result = render_scatter(
                tsv_path=self.tsv_input,
                outfile_prefix=outfile_prefix,
                aminoacids=True,
                linear_scaling=True,
                threshold=0.001,
                show_bokeh=False,
                show_mplcursors=False,
            )
            self.assertIsNotNone(result['figure'])
            self.assertGreater(len(result['files_written']), 0)


class TestCalculateFrequenciesAPI(unittest.TestCase):
    """Integration test: calculate_frequencies() with real test data."""

    def setUp(self):
        """Set up paths to test fixtures."""
        self.tests_dir = os.path.dirname(os.path.abspath(__file__))
        self.inputs_dir = os.path.join(self.tests_dir, "inputs")

    def test_calculate_frequencies_basic(self):
        """calculate_frequencies() produces a DataFrame with correct columns."""
        from mutation_scatter_plot.api import calculate_frequencies

        alignment = os.path.join(self.inputs_dir, "test.fasta")
        reference = os.path.join(self.inputs_dir, "MN908947.3_S.fasta")

        if not (os.path.exists(alignment) and os.path.exists(reference)):
            self.skipTest("Test fixtures not available")

        with tempfile.TemporaryDirectory() as tmpdir:
            outfile_prefix = os.path.join(tmpdir, "api_freq_test")
            df = calculate_frequencies(
                alignment_file=alignment,
                reference_file=reference,
                outfile_prefix=outfile_prefix,
                padded_reference=True,
                threads=1,
            )

            # DataFrame should have the expected columns
            expected_cols = [
                'padded_position', 'position', 'original_aa', 'mutant_aa',
                'frequency', 'original_codon', 'mutant_codon',
                'observed_codon_count', 'total_codons_per_site',
            ]
            for col in expected_cols:
                self.assertIn(col, df.columns, f"Missing column: {col}")

            # Should have data
            self.assertGreater(len(df), 0)

            # TSV file should exist on disk
            self.assertTrue(os.path.exists(
                os.path.join(tmpdir, "api_freq_test.tsv")
            ))


class TestRenderTimelineAPI(unittest.TestCase):
    """Integration test: render_timeline() produces correct outputs."""

    def setUp(self):
        """Set up paths to test fixtures."""
        self.tests_dir = os.path.dirname(os.path.abspath(__file__))
        self.timeline_dir = os.path.join(
            self.tests_dir, "inputs", "timeline"
        )

    def test_render_timeline_basic(self):
        """render_timeline() produces files from per-month TSVs."""
        from mutation_scatter_plot.api import render_timeline

        if not os.path.isdir(self.timeline_dir):
            self.skipTest("Timeline fixtures not available")

        with tempfile.TemporaryDirectory() as tmpdir:
            outfile_prefix = os.path.join(tmpdir, "api_timeline")
            result = render_timeline(
                directory=self.timeline_dir,
                positions=['N501Y', '484', 'D614G'],
                outfile_prefix=outfile_prefix,
                aminoacids=True,
                show_bokeh=False,
            )

            self.assertIn('data', result)
            self.assertIn('files_written', result)
            self.assertGreater(len(result['data'].points), 0)

            # Files should have been written
            self.assertGreater(len(result['files_written']), 0)

            # PNG should exist
            extensions = {os.path.splitext(f)[1] for f in result['files_written']}
            self.assertIn('.png', extensions)

    def test_render_timeline_no_data(self):
        """render_timeline() with position that has no data returns empty."""
        from mutation_scatter_plot.api import render_timeline

        if not os.path.isdir(self.timeline_dir):
            self.skipTest("Timeline fixtures not available")

        with tempfile.TemporaryDirectory() as tmpdir:
            outfile_prefix = os.path.join(tmpdir, "api_timeline_empty")
            result = render_timeline(
                directory=self.timeline_dir,
                positions=['999'],  # position not in data
                outfile_prefix=outfile_prefix,
                aminoacids=True,
            )

            self.assertEqual(len(result['data'].points), 0)
            self.assertEqual(result['files_written'], [])

    def test_render_timeline_missing_dir(self):
        """render_timeline() with nonexistent directory raises FileNotFoundError."""
        from mutation_scatter_plot.api import render_timeline

        with self.assertRaises(FileNotFoundError):
            render_timeline(
                directory='/nonexistent/path',
                positions=['501'],
            )


class TestExample4MatrixComparison(unittest.TestCase):
    """Test that Example 4 from the implementation plan works."""

    def test_blosum_matrix_comparison(self):
        """load_matrix + get_score works for multiple BLOSUM matrices."""
        from mutation_scatter_plot import scatter_options
        from mutation_scatter_plot.mutation_scatter_plot import load_matrix
        from mutation_scatter_plot.mutation_scatter_plot.core import get_score

        matrices = ['BLOSUM45', 'BLOSUM62', 'BLOSUM80', 'BLOSUM90']
        scores = {}
        for matrix_name in matrices:
            opts = scatter_options(
                matrix=matrix_name, outfile_prefix='/tmp/dummy',
            )
            matrix, _name, _min_s, _max_s, _prefix = load_matrix(opts)
            score = get_score(opts, matrix, False, 'Q', 'R')
            scores[matrix_name] = score

        # All should return integer scores
        for name, score in scores.items():
            self.assertIsInstance(score, (int, float), f"{name} score is not numeric")

        # BLOSUM matrices produce consistent Q->R scores
        self.assertEqual(scores['BLOSUM62'], scores['BLOSUM80'],
                         "BLOSUM62 and BLOSUM80 should give same Q→R score")


if __name__ == "__main__":
    unittest.main()

