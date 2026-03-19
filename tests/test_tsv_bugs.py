import unittest
import os
import subprocess
import pandas as pd
import numpy as np

class TestTsvBugs(unittest.TestCase):
    def test_legacy_tsv_loading(self):
        # 6-column headerless TSV
        # position, original_aa, mutant_aa, frequency, original_codon, mutant_codon
        tsv_content = "1\tM\tV\t0.5\tATG\tGTG\n"
        with open("old.frequencies.tsv", "w") as f:
            f.write(tsv_content)

        with open("old.frequencies.unchanged_codons.tsv", "w") as f:
            f.write("1\tM\tM\t0.5\tATG\tATG\n")

        env = os.environ.copy()
        env["MPLBACKEND"] = "agg"
        env["PYTHONPATH"] = "src"

        cmd = [
            "python3", "-m", "mutation_scatter_plot.mutation_scatter_plot.cli",
            "--tsv=old.frequencies.tsv",
            "--outfile-prefix=old",
            "--xmin=1",
            "--xmax=10"
        ]

        subprocess.run(cmd, capture_output=True, text=True, env=env)

        expected_rendered = "old.BLOSUM80.amino_acid_changes.actually_rendered.tsv"
        self.assertTrue(os.path.exists(expected_rendered))

    def test_legacy_tsv_hover_and_colorbar_logic(self):
        # 6-column headerless TSV (no observed_codon_count)
        # position, original_aa, mutant_aa, frequency, original_codon, mutant_codon
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
        self.assertGreater(len(data), 0)

if __name__ == "__main__":
    unittest.main()
