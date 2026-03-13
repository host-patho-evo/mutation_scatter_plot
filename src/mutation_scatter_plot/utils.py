# This work © 2025 by Jiří Zahradník and Martin Mokrejš
# (First Medical Faculty - Charles University in Prague) is licensed under
# Creative Commons Attribution 4.0 International. To view a copy of this
# license, visit https://creativecommons.org/licenses/by/4.0/

"""Shared utilities used by both mutation_scatter_plot and calculate_codon_frequencies."""

from Bio.Seq import translate


def alt_translate(seq):
    """Biopython cannot sometimes translate a sequence but one can get around by
    splitting it into codons and merging later back.
    https://github.com/biopython/biopython/pull/4992
    https://github.com/biopython/biopython/pull/4992#issuecomment-2865429105
    """
    codons = (seq[i:i+3] for i in range(0, len(seq), 3))
    codons = ("NNN" if "-" in codon and codon != "---" else codon for codon in codons)
    return "".join(translate(codon, gap='-') for codon in codons)
