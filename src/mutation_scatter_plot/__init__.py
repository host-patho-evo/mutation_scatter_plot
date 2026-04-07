# -*- coding: utf-8 -*-
# vim: set fileencoding=utf-8 ts=4 sw=4 expandtab :
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
def alt_translate(seq, table=1):
    """Translate a nucleotide sequence (possibly alignment-padded) to protein.

    Biopython's translate() cannot handle sequences with gap characters
    without crashing.  This wrapper splits the input into per-codon chunks
    and calls Bio.Seq.translate() on each, applying the semantics of the
    proposed --respect-alignment flag from Biopython PR #4992:

      '---' → '-'   (full gap codon → gap amino acid)
      'TC-' → 'TCN' → 'S'  ('-' treated as 'N'; Biopython resolves per-codon)
      'AT-' → 'ATN' → 'X'  (ATN ambiguous: Ile or Met)
      'TCN' → 'S'   (IUPAC ambiguity → resolved by Biopython per-codon)
      'TCA' → 'S'   (standard codon → Biopython lookup, lru_cache hit)

    *table* selects the NCBI genetic code table (default: 1 = standard).
    The lru_cache key includes *table* so different codes are cached
    independently.  The cache is most effective when called with single
    3-nt codon strings (the common case in mutation_scatter_plot/__init__.py).

    https://github.com/biopython/biopython/pull/4992
    """
    result = []
    for i in range(0, len(seq), 3):
        codon = seq[i:i + 3]
        if not codon:
            continue
        if len(codon) < 3:
            # Incomplete trailing codon (1 or 2 nt): we know a codon position
            # exists but cannot determine the amino acid — return 'X'.
            # Calling Biopython on a sub-3-nt fragment causes a
            # BiopythonWarning and may raise in future versions; this was one
            # of the original motivations for alt_translate() to exist.
            result.append('X')
            continue
        if '-' in codon and codon != '---':
            # Partial-gap codon: treat each '-' as 'N' (unknown nucleotide).
            # A gap in an alignment column means 'unknown', same as N.
            # This lets Biopython resolve correctly per-codon:
            #   TC- → TCN → S  (four-fold degenerate → Serine)
            #   AT- → ATN → X  (ATN can be Ile or Met → ambiguous)
            codon = codon.replace('-', 'N')
        result.append(translate(codon, table=table, gap='-'))
    return ''.join(result)
