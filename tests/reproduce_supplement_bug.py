
import os
import subprocess
import pandas as pd

def main():
    # Main frequencies file
    with open("bug.frequencies.tsv", "w") as f:
        f.write("padded_position\tposition\toriginal_aa\tmutant_aa\tfrequency\toriginal_codon\tmutant_codon\n")
        f.write("1\t1\tM\tV\t0.5\tATG\tGTG\n")

    # Unchanged codons file with a position NOT in the main file
    with open("bug.frequencies.unchanged_codons.tsv", "w") as f:
        f.write("padded_position\tposition\toriginal_aa\tmutant_aa\tfrequency\toriginal_codon\tmutant_codon\n")
        f.write("2\t2\tA\tA\t1.0\tGCC\tGCC\n")

    print("Running mutation_scatter_plot to check if it supplements from unchanged_codons.tsv")
    env = os.environ.copy()
    env["MPLBACKEND"] = "agg"

    cmd = [
        "mutation_scatter_plot",
        "--tsv=bug.frequencies.tsv",
        "--outfile-prefix=bug",
        "--xmin=1",
        "--xmax=10",
        "--debug=1"
    ]

    result = subprocess.run(cmd, capture_output=True, text=True, env=env)
    print(f"STDOUT:\n{result.stdout}")

    # If it correctly supplemented, _padded_position2position should have entry for position 2.
    # We can check the actually_rendered.tsv if it was written or just debug output.
    # The debug output of collect_scatter_data should show _real_aa_positions
    if "Debug: _real_aa_positions=[1, 2]" in result.stdout:
        print("SUCCESS: Supplemented position 2 from unchanged_codons.tsv")
    elif "Debug: _real_aa_positions=[1]" in result.stdout:
        print("BUG CONFIRMED: Did NOT supplement from unchanged_codons.tsv")
    else:
        print("Could not determine result from output")

    # Cleanup
    for f in [
        "bug.frequencies.tsv", "bug.frequencies.unchanged_codons.tsv",
        "bug.BLOSUM80.amino_acid_changes.actually_rendered.tsv",
        "bug.BLOSUM80.amino_acid_changes.aa.frequencies.colors.tsv",
        "bug.BLOSUM80.amino_acid_changes.codon.frequencies.colors.tsv",
        "bug.BLOSUM80.amino_acid_changes.png",
        "bug.BLOSUM80.amino_acid_changes.pdf",
        "bug.BLOSUM80.amino_acid_changes.html"
    ]:
        if os.path.exists(f):
            os.remove(f)

if __name__ == "__main__":
    main()
