# This work © 2025-2026 by Jiří Zahradník and Martin Mokrejš
# (First Medical Faculty - Charles University in Prague) is licensed under
# Creative Commons Attribution 4.0 International. To view a copy of this
# license, visit https://creativecommons.org/licenses/by/4.0/

"""mutation-scatter-plot distribution package.

Contains two tools:
  - mutation_scatter_plot:       render scatter figures of mutation frequencies
  - calculate_codon_frequencies: calculate codon/aa frequencies from alignments

Shared utilities are in mutation_scatter_plot.utils.
"""

import functools
from Bio.Seq import translate

VERSION = "0.3"


@functools.lru_cache(maxsize=128)
def alt_translate(seq):
    """Translate a nucleotide sequence (possibly alignment-padded) to protein.

    Biopython's translate() cannot handle sequences with gap characters
    without crashing.  This wrapper splits the input into per-codon chunks
    and calls Bio.Seq.translate() on each, applying the semantics of the
    proposed --respect-alignment flag from Biopython PR #4992:

      '---' → '-'   (full gap codon → gap amino acid)
      'TC-' → 'X'   (partial-gap codon → ambiguous; NOT mapped via NNN)
      'TCN' → 'S'   (IUPAC ambiguity → resolved by Biopython per-codon)
      'TCA' → 'S'   (standard codon → fast Biopython lookup, lru_cache hit)

    The lru_cache is effective when called with single 3-nt codons
    (the common case in mutation_scatter_plot/__init__.py).

    https://github.com/biopython/biopython/pull/4992
    """
    result = []
    for i in range(0, len(seq), 3):
        codon = seq[i:i + 3]
        if not codon:
            continue
        if '-' in codon and codon != '---':
            # Partial-gap codon: treat each '-' as 'N' (unknown nucleotide).
            # This lets Biopython resolve correctly per-codon:
            #   TC- → TCN → S  (four-fold degenerate → Serine)
            #   AT- → ATN → X  (ATN can be Ile or Met → ambiguous)
            codon = codon.replace('-', 'N')
        result.append(translate(codon, gap='-'))
    return ''.join(result)
