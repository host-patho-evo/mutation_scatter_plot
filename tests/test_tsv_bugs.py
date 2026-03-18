
import unittest
import os
import subprocess
import pandas as pd

class TestTsvBugs(unittest.TestCase):
    def tearDown(self):
        # The outfile_prefix is modified by the matrix and colormap in the real app
        expected_rendered = "old.BLOSUM80.amino_acid_changes.actually_rendered.tsv"
        files_to_remove = [
            "old.frequencies.tsv", "old.frequencies.unchanged_codons.tsv",
            expected_rendered,
            "old.BLOSUM80.amino_acid_changes.aa.frequencies.colors.tsv",
            "old.BLOSUM80.amino_acid_changes.codon.frequencies.colors.tsv",
            "old.BLOSUM80.amino_acid_changes.png",
            "old.BLOSUM80.amino_acid_changes.pdf",
            "old.BLOSUM80.amino_acid_changes.html"
        ]
        for f in files_to_remove:
            if os.path.exists(f):
                os.remove(f)

    def test_legacy_tsv_loading(self):
        # 6-column headerless TSV
        tsv_content = "1\tM\tV\t0.5\tATG\tGTG\n"
        with open("old.frequencies.tsv", "w") as f:
            f.write(tsv_content)

        with open("old.frequencies.unchanged_codons.tsv", "w") as f:
            f.write("1\tM\tM\t0.5\tATG\tATG\n")

        env = os.environ.copy()
        env["MPLBACKEND"] = "agg"

        cmd = [
            "mutation_scatter_plot",
            "--tsv=old.frequencies.tsv",
            "--outfile-prefix=old",
            "--xmin=1",
            "--xmax=10"
        ]

        subprocess.run(cmd, capture_output=True, text=True, env=env)

        expected_rendered = "old.BLOSUM80.amino_acid_changes.actually_rendered.tsv"
        self.assertTrue(os.path.exists(expected_rendered))

        rendered_df = pd.read_csv(expected_rendered, sep='\t', header=None)
        # Bug 1: header=0 caused loss of the only row
        self.assertEqual(len(rendered_df), 1)
        # Bug 2: 'padded_position' KeyError
        self.assertEqual(rendered_df.iloc[0, 0], 1)
        self.assertEqual(rendered_df.iloc[0, 1], 'M')

    def test_legacy_tsv_hover_logic(self):
        # 6-column headerless TSV (no observed_codon_count)
        tsv_content = "1\tM\tV\t0.5\tATG\tGTG\n"
        with open("old_hover.frequencies.tsv", "w") as f:
            f.write(tsv_content)
        with open("old_hover.frequencies.unchanged_codons.tsv", "w") as f:
            f.write("1\tM\tM\t0.5\tATG\tATG\n")

        import mutation_scatter_plot.mutation_scatter_plot as msp
        from argparse import Namespace

        myoptions = Namespace(
            tsv_file_path="old_hover.frequencies.tsv",
            column_with_frequencies="frequency",
            offset=0,
            threshold=0.001,
            showstop=False,
            showdel=False,
            showins=False,
            showx=False,
            matrix="BLOSUM80",
            matrix_file=None,
            outfile_prefix="old_hover",
            colormap="amino_acid_changes",
            aminoacids=False,
            debug=0,
            xaxis_bins=0,
            xaxis_major_ticks_spacing=10,
            xaxis_minor_ticks_spacing=5,
            xaxis_label_start=0,
            xmin=0,
            xmax=0,
            shortlegend=True,
            disable_2nd_Y_axis=False,
            bokeh_sqrt_size=True,
            dpi=600
        )

        # Minimal pipeline to get to render_matplotlib
        _matrix, _matrix_name, _, _, _outfile_prefix = msp.load_matrix(myoptions)
        _padded_position2position = {}
        _df, _padded_position2position = msp.load_and_clean_dataframe(myoptions, myoptions.tsv_file_path, _outfile_prefix, _padded_position2position)

        (
            _amino_acids, _codons_whitelist, _codons_whitelist2,
            _final_sorted_whitelist,
            _unique_aa_padded_positions, _unique_codon_padded_positions,
            _old_aa_table, _new_aa_table, _old_codon_table, _new_codon_table,
            _calculated_aa_offset, _padded_position2position,
        ) = msp.build_frequency_tables(myoptions, _df, _padded_position2position)

        import matplotlib.pyplot as plt
        import matplotlib
        matplotlib.use('agg')
        _figure, _ax1, _ax2, _ax3, _ax4, _xmin, _xmax = msp.setup_matplotlib_figure(
            myoptions, "title", "0", _matrix_name, _amino_acids,
            _codons_whitelist, _final_sorted_whitelist,
            _unique_aa_padded_positions, _unique_codon_padded_positions,
            _new_aa_table, _new_codon_table, _padded_position2position
        )

        _table = _new_codon_table
        data = msp.collect_scatter_data(myoptions, _df, _table, _outfile_prefix, _matrix, _amino_acids, _codons_whitelist2, _padded_position2position)

        # This should NOT crash with KeyError: 'observed_codon_count'
        try:
            msp.render_matplotlib(
                myoptions, _figure, _ax1, _ax2, _ax3, _ax4, _outfile_prefix,
                data[8], data[9], data[10], data[1], data[0], data[2],
                _matrix, _matrix_name,
                _new_aa_table, _new_codon_table, _df, _codons_whitelist2,
                _final_sorted_whitelist,
                _calculated_aa_offset, _padded_position2position
            )
        except Exception as e:
            self.fail(f"render_matplotlib crashed with legacy TSV: {e}")
        finally:
            plt.close('all')
            for f in ["old_hover.frequencies.tsv", "old_hover.frequencies.unchanged_codons.tsv", "old_hover.BLOSUM80.amino_acid_changes.actually_rendered.tsv", "old_hover.BLOSUM80.amino_acid_changes.codon.frequencies.colors.tsv", "old_hover.BLOSUM80.amino_acid_changes.png", "old_hover.BLOSUM80.amino_acid_changes.pdf"]:
                if os.path.exists(f):
                    os.remove(f)

if __name__ == "__main__":
    unittest.main()
