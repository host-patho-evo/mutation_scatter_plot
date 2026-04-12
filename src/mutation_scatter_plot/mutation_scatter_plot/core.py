# -*- coding: utf-8 -*-
# vim: set fileencoding=utf-8 ts=4 sw=4 expandtab :
# pylint: disable=too-many-lines
# This work © 2025-2026 by Jiří Zahradník and Martin Mokrejš
# (First Medical Faculty - Charles University in Prague) is licensed under
# Creative Commons Attribution 4.0 International. To view a copy of this
# license, visit https://creativecommons.org/licenses/by/4.0/

"""Shared scoring, coloring, and matrix-loading logic.

This module is the single source for functions that are shared between
``mutation_scatter_plot`` (the existing scatter-plot tool) and
``mutation_timeline_plot`` (the new timeline tool).

By extracting these functions into a separate module, we avoid duplicating
the BLOSUM scoring, colormap construction, and size/colour assignment logic.

The parent ``__init__.py`` re-exports everything defined here so that all
existing ``from mutation_scatter_plot.mutation_scatter_plot import …``
imports continue to work unmodified.
"""

import os
import re
import sys
import typing

# https://pypi.org/project/blosum/
import blosum
import matplotlib
import matplotlib.pyplot as plt
import numpy as np

from .. import alt_translate

def get_colormap(myoptions, colormapname):
    """Retrieve or create a matplotlib colormap and normalization object."""

    _norm = None
    _colors = None
    _cmap = None

    _micro_cvd_gray = ["#F5F5F5", "#D6D6D6", "#B7B7B7", "#8B8B8B", "#616161"]
    _micro_cvd_green = ["#DDFFA0", "#BDEC6F", "#97CE2F", "#6D9F06", "#4E7705"]
    _micro_cvd_orange = ["#FFD5AF", "#FCB076", "#F09163", "#C17754", "#9D654C"]
    _micro_cvd_blue = ["#E7F4FF", "#BCE1FF", "#7DCCFF", "#56B4E9", "#098BD9"]
    _micro_cvd_turquoise = ["#A3E4D7", "#48C9B0", "#43BA8F", "#009E73", "#148F77"]
    _micro_cvd_purple = ["#EFB6D6", "#E794C1", "#CC79A7", "#A1527F", "#7D3560"]

    _custom_found = True
    # Handle custom hardcoded colormaps first to avoid "not found" warnings
    if colormapname == 'amino_acid_changes':
        # All entries are explicit hex strings — no CSS named colours — so
        # that every downstream consumer (_blend_with_white, Bokeh palette
        # builder, matplotlib ListedColormap) receives a uniform format.
        # Previously index 29 was the CSS name "palegreen"; it is now its
        # exact hex equivalent #98fb98 (R=152 G=251 B=152).
        _colors = ["#930000", "#930000", "#930000", "#930000", "#930000", "#930000", "#960000", "#580041", "#8200ff", "#c500ff", "#ff00fd", "#CC79A7", "#eea1d0", "#cc0000", "#ff0000", "#ff4f00", "#ff7c7c", "#ff9999", "#c58a24",
                   "#9c644b", "#ffff00", "#ffcc00", "#ffa200", "#7DCCFF", "#0042ff", "#0000ff", "#D6D6D6", "#B7B7B7", "#8B8B8B", "#98fb98", "#bbff00", "#97CE2F", "#219f11", "#930000", "#930000", "#930000", "#930000", "#930000", "#930000"]

        if getattr(myoptions, 'spread_colormap_virtual_matrix', False):
            _local_vmin = getattr(myoptions, 'matrix_min_theoretical', -19)
            _local_vmax = getattr(myoptions, 'matrix_max_theoretical', 19)
        else:
            _local_vmin = getattr(myoptions, 'cmap_actual_vmin', -19)
            _local_vmax = getattr(myoptions, 'cmap_actual_vmax', 19)
        if _local_vmin >= _local_vmax:
            _local_vmin = _local_vmax - 1
        myoptions.cmap_vmin = _local_vmin
        myoptions.cmap_vmax = _local_vmax

        _cmap = matplotlib.colors.ListedColormap(_colors, "amino_acid_changes", len(_colors))
        myoptions.colormap = 'amino_acid_changes'
        _norm = matplotlib.colors.BoundaryNorm(np.arange(-19, 19, 1), len(_colors))

    elif colormapname == 'dkeenan_26cols':
        _colors = ["#00B7FF", "#004DFF", "#00FFFF", "#826400", "#580041", "#FF00FF", "#00FF00", "#C500FF", "#B4FFD7", "#FFCA00", "#969600", "#B4A2FF", "#C20078",
                   "#000000", "#0000C1", "#FF8B00", "#FFC8FF", "#666666", "#FF0000", "#CCCCCC", "#009E8F", "#D7A870", "#8200FF", "#960000", "#BBFF00", "#FFFF00", "#006F00"]

        if getattr(myoptions, 'spread_colormap_virtual_matrix', False):
            _local_vmin = getattr(myoptions, 'matrix_min_theoretical', -13)
            _local_vmax = getattr(myoptions, 'matrix_max_theoretical', 13)
        else:
            _local_vmin = getattr(myoptions, 'cmap_actual_vmin', -13)
            _local_vmax = getattr(myoptions, 'cmap_actual_vmax', 13)
        if _local_vmin >= _local_vmax:
            _local_vmin = _local_vmax - 1
        myoptions.cmap_vmin = _local_vmin
        myoptions.cmap_vmax = _local_vmax

        _cmap = matplotlib.colors.ListedColormap(_colors, "dkeenan_26cols")
        myoptions.colormap = 'dkeenan_26cols'
        _norm = matplotlib.colors.BoundaryNorm(np.arange(-13, 13, 1), len(_colors))

    elif colormapname == 'microshades_cvd_palettes':
        _colors = list(_micro_cvd_orange) + list(_micro_cvd_turquoise) + list(_micro_cvd_blue) + \
            list(_micro_cvd_purple) + list(_micro_cvd_green) + list(_micro_cvd_gray)
        _cmap = matplotlib.colors.ListedColormap(_colors, 'microshades_cvd_palettes')

    elif colormapname == 'microshades_cvd_palettes_r':
        _colors = list(_micro_cvd_orange[::-1]) + list(_micro_cvd_turquoise[::-1]) + list(_micro_cvd_blue[::-1]) + \
            list(_micro_cvd_purple[::-1]) + list(_micro_cvd_green[::-1]) + list(_micro_cvd_gray[::-1])
        _cmap = matplotlib.colors.ListedColormap(_colors, 'microshades_cvd_palettes_r')

    elif colormapname == 'merged':
        _colors1 = matplotlib.colormaps['gist_heat'](np.linspace(0.2, 0.8, 9))
        _colors2 = matplotlib.colormaps['cividis_r'](np.linspace(0.2, 0.8, 9))
        _colors3 = matplotlib.colormaps['GnBu_r'](np.linspace(0.2, 0.8, 9))
        _colors = np.vstack((_colors1, _colors2, _colors3))
        _cmap = matplotlib.colors.LinearSegmentedColormap.from_list('merged', _colors)

    elif colormapname == 'merged1':
        _colors1 = matplotlib.colormaps['viridis'](np.linspace(0.2, 0.8, 9))
        _colors2 = matplotlib.colormaps['coolwarm_r'](np.linspace(0., 1, 9))
        _colors3 = matplotlib.colormaps['GnBu_r'](np.linspace(0.2, 0.8, 9))
        _colors = np.vstack((_colors1, _colors2, _colors3))
        _cmap = matplotlib.colors.LinearSegmentedColormap.from_list('merged', _colors)
    else:
        _custom_found = False

    if _custom_found:
        return _norm, _cmap, _colors

    _wished_cmapname_prefix = ''
    _wished_cmapname_num = 0
    if '_' in colormapname:
        try:
            _wished_cmapname_prefix = '_'.join(colormapname.split('_')[:-1])
            _wished_cmapname_num = int(colormapname.split('_')[-1])
        except (ValueError, IndexError):
            _wished_cmapname_prefix = colormapname
            _wished_cmapname_num = 0

    try:
        _cmap = matplotlib.colormaps.get_cmap(colormapname)
    except Exception:  # pylint: disable=broad-except
        try:
            _cmap = plt.get_cmap(colormapname)
        except Exception:  # pylint: disable=broad-except
            try:
                import palettable  # pylint: disable=import-error
                _cmap_from_palettable = palettable.colorbrewer.get_map(_wished_cmapname_prefix, 'diverging', _wished_cmapname_num)
                _cmap = matplotlib.colors.ListedColormap(_cmap_from_palettable.mpl_colors)
            except Exception as _e3:  # pylint: disable=broad-except
                print(f"Info: colorbrewer diverging lookup failed for '{_wished_cmapname_prefix}': {_e3}")
                try:
                    import palettable  # pylint: disable=import-error
                    _cmap_from_palettable = palettable.colorbrewer.get_map(_wished_cmapname_prefix, 'sequential', _wished_cmapname_num)
                    _cmap = matplotlib.colors.ListedColormap(_cmap_from_palettable.mpl_colors)
                except Exception as _e4:  # pylint: disable=broad-except
                    print(f"Info: colorbrewer sequential lookup failed for '{_wished_cmapname_prefix}': {_e4}")
                    try:
                        import palettable.scientific.diverging  # pylint: disable=import-error
                        _cmap_from_palettable = palettable.scientific.diverging.get_map(_wished_cmapname_prefix)
                        _cmap = matplotlib.colors.ListedColormap(_cmap_from_palettable.mpl_colors)
                    except Exception as _e5:  # pylint: disable=broad-except
                        print(f"Info: palettable scientific diverging lookup failed for '{_wished_cmapname_prefix}': {_e5}")
                        try:
                            import palettable.scientific.sequential  # pylint: disable=import-error
                            _cmap_from_palettable = palettable.scientific.sequential.get_map(_wished_cmapname_prefix)
                            _cmap = matplotlib.colors.ListedColormap(_cmap_from_palettable.mpl_colors)
                        except Exception as _e6:  # pylint: disable=broad-except
                            print(f"Info: palettable scientific sequential lookup failed for '{_wished_cmapname_prefix}': {_e6}")
                            try:
                                import palettable.mygbm  # pylint: disable=import-error,no-name-in-module
                                _cmap_from_palettable = palettable.mygbm.get_map(_wished_cmapname_prefix)  # pylint: disable=no-member
                                _cmap = matplotlib.colors.ListedColormap(_cmap_from_palettable.mpl_colors)
                            except Exception as _e7:  # pylint: disable=broad-except
                                print(f"Info: palettable mygbm lookup failed for '{_wished_cmapname_prefix}': {_e7}")
                                print(f"Warning: Colormap {colormapname} not found, falling back to coolwarm_r")
                                _cmap = matplotlib.colormaps.get_cmap('coolwarm_r')
    # Standard matplotlib colormaps (non-custom / no ListedColormap): sample
    # the cmap uniformly over [0, 1] to build a discrete colour list — the
    # original v0.3 approach.  _norm is intentionally kept None to activate the
    # direct-index code path in adjust_size_and_color and a separate
    # continuous-Normalize colorbar in render_matplotlib.
    if _norm is None and _cmap is not None:
        _n_std = 23  # covers scores −11 … +11 (BLOSUM80 non-stop range; DEL/INS use _min_theoretical_score)
        _colors = [
            matplotlib.colors.to_hex(_cmap(i / (_n_std - 1)))
            for i in range(_n_std)
        ]
        # _norm intentionally stays None — signals the separate code path

    if getattr(myoptions, 'spread_colormap_virtual_matrix', False):
        if getattr(myoptions, 'cmap_vmin', None) is None:
            myoptions.cmap_vmin = getattr(myoptions, 'matrix_min_theoretical', -11)
            myoptions.cmap_vmax = getattr(myoptions, 'matrix_max_theoretical', 11)
    else:
        if getattr(myoptions, 'cmap_vmin', None) is None:
            myoptions.cmap_vmin = getattr(myoptions, 'cmap_actual_vmin', -11)
            myoptions.cmap_vmax = getattr(myoptions, 'cmap_actual_vmax', 11)

    return _norm, _cmap, _colors


def resolve_codon_or_aa(myoptions, old_codon_or_aa, new_codon_or_aa):
    """Determine the codon or amino acid labels to use based on mode."""

    if myoptions.aminoacids:
        _old_codon_or_aa = alt_translate(old_codon_or_aa)
        if new_codon_or_aa == 'NNN':
            _new_codon_or_aa = 'X'
        elif len(new_codon_or_aa) > 1 and new_codon_or_aa not in ('---', 'DEL', 'INS'):
            print(f"Info: Weird, the new_codon_or_aa={new_codon_or_aa}")
            _new_codon_or_aa = new_codon_or_aa.upper()
        else:
            _new_codon_or_aa = new_codon_or_aa.upper()
    else:
        _old_codon_or_aa = old_codon_or_aa.upper()
        _new_codon_or_aa = new_codon_or_aa.upper()

    if len(_old_codon_or_aa) > 1 and old_codon_or_aa not in ('---', 'DEL', 'INS'):
        _codon_on_input = True
    elif len(_new_codon_or_aa) > 1 and new_codon_or_aa not in ('---', 'DEL', 'INS'):
        _codon_on_input = True
    else:
        _codon_on_input = False
    return _codon_on_input, _old_codon_or_aa, _new_codon_or_aa


# Module-level cache for BLOSUM substitution scores.
#
# Performance rationale
# ---------------------
# blosum.BLOSUM is a subclass of collections.defaultdict and stores scores in a
# nested dict of dicts. Individual lookups (matrix[aa1][aa2]) are O(1) at the
# Python dict level, so the structure itself is not the bottleneck.
#
# The problem is that our caller loop is O(N_codons_or_aa x N_positions), and
# with standard amino acid alphabets there are at most 20x20 = 400 distinct
# (old_aa, new_aa) pairs regardless of how many positions the alignment has.
# Before this cache was added, every cell in that loop independently called
# matrix[old_aa][new_aa], meaning the same pairs were looked up thousands of
# times (e.g. 2331 calls observed during profiling of a real dataset).
#
# Additionally, when codon_on_input=True, each call also invoked alt_translate()
# twice to convert codons to amino acids before doing the matrix lookup. While
# alt_translate() is itself cached via @functools.lru_cache, two extra function
# calls and cache lookups per score request still accumulate when multiplied
# across all positions.
#
# The dictionary is keyed by (old_aa, new_aa) string tuples — always single
# uppercase amino acid letters after translation. functools.lru_cache is not
# used here because the matrix object (blosum.BLOSUM) is not hashable, so a
# module-level plain dict is the simplest correct solution.
#
# The cache is cleared at the start of each collect_scatter_data() call to
# ensure correctness when different matrices (e.g. BLOSUM62 vs BLOSUM80) are
# used in successive calls within the same Python process.
_score_cache: dict = {}

# Module-level worst-case score, set by load_matrix() from the active matrix.
#
# Used as the sentinel score for DEL/INS entries and OverflowError fallbacks
# in get_score() / adjust_size_and_color() without requiring those functions to
# iterate the matrix themselves.  Left as None until load_matrix() is called so
# that any code path that scores mutations before loading a matrix fails loudly
# rather than silently using a stale -11 default.
_min_theoretical_score: typing.Optional[int] = None  # pylint: disable=invalid-name


def get_score(myoptions, matrix, codon_on_input, old_codon_or_aa, new_codon_or_aa):
    """Return the integer substitution score for an amino acid (or codon) pair.

    Parameters
    ----------
    myoptions : argparse.Namespace
        CLI options; only used for error messages (myoptions.matrix name).
    matrix : blosum.BLOSUM
        The substitution scoring matrix (a nested defaultdict subclass).
        blosum.BLOSUM stores precomputed scores in plain Python dicts, so
        individual lookups are O(1) and the library has no caching of its own
        — it simply does not know whether the caller will repeat the same pair.
    codon_on_input : bool
        If True, old_codon_or_aa and new_codon_or_aa are 3-nt codons that must
        be translated to amino acids before the matrix lookup.
    old_codon_or_aa : str
        The reference codon or amino acid (before mutation).
    new_codon_or_aa : str
        The observed codon or amino acid (after mutation).

    Returns
    -------
    int
        The BLOSUM (or custom matrix) substitution score.  On arithmetic
        overflow the score falls back to ``_min_theoretical_score`` — the
        stop-codon penalty from the active matrix (e.g. −6 for BLOSUM80),
        set by ``load_matrix()``.

    Performance notes
    -----------------
    Results are cached in the module-level ``_score_cache`` dict, keyed by the
    (old_aa, new_aa) pair after translation. With a standard 20-amino-acid
    alphabet there are at most 400 unique pairs, so the cache saturates quickly
    and subsequent calls for the same pair cost only a single dict lookup.
    The cache must be cleared (via ``_score_cache.clear()``) before each
    independent plot run when the matrix may change.
    """
    if codon_on_input:
        _old_aa = alt_translate(old_codon_or_aa)
        _new_aa = alt_translate(new_codon_or_aa)
    else:
        _old_aa = old_codon_or_aa
        _new_aa = new_codon_or_aa
    _cache_key = (_old_aa, _new_aa)
    if _cache_key in _score_cache:
        return _score_cache[_cache_key]
    try:
        _score = int(matrix[_old_aa][_new_aa])
    except OverflowError:
        # Fall back to the pre-computed module-level worst score rather than
        # a hardcoded constant.  _min_theoretical_score is set by load_matrix()
        # from the active matrix (e.g. -6 for BLOSUM80).
        if _min_theoretical_score is None:
            raise RuntimeError(
                "_min_theoretical_score is not initialised — call load_matrix() before scoring"
            ) from None
        _score = _min_theoretical_score
    except KeyError as exc:
        raise ValueError(
            f"Cannot get a score for myoptions.matrix='{myoptions.matrix}', old_codon_or_aa='{old_codon_or_aa}', new_codon_or_aa='{new_codon_or_aa}'") from exc
    _score_cache[_cache_key] = _score
    return _score


def adjust_size_and_color(myoptions, frequency, codon_on_input, old_codon_or_aa, new_codon_or_aa, _old_codon_or_aa, _new_codon_or_aa, matrix, norm, colors):
    """Assign BLOSUM score, relative size, and colour to one scatter data-point.

    Parameters
    ----------
    myoptions : argparse.Namespace
        Parsed CLI options (used for ``debug`` and ``colormap`` context).
    frequency : Decimal
        Observed mutation frequency; controls the rendered circle area.
    codon_on_input : bool
        True when old_codon_or_aa / _old_codon_or_aa are 3-nt codon strings
        rather than single-letter amino acid codes.
    old_codon_or_aa : str
        Raw reference string from the TSV: a codon (e.g. ``'GCC'``) in codon
        mode or an amino acid letter (e.g. ``'A'``) in aminoacids mode.
    new_codon_or_aa : str
        Raw observed string from the TSV (codon or amino acid).
    _old_codon_or_aa : str
        Resolved reference: codon string if ``codon_on_input`` else translated
        amino acid (output of ``resolve_codon_or_aa``).
    _new_codon_or_aa : str
        Resolved observed codon or amino acid.
    matrix : blosum.BLOSUM
        Substitution-score matrix used for non-synonymous mutations.
    norm : matplotlib.colors.BoundaryNorm or None
        Discrete colour normaliser for ``amino_acid_changes`` colormap;
        ``None`` for continuous colormaps (``coolwarm_r``, etc.).
    colors : list[str]
        Palette as hex strings, indexed by ``norm(score)`` or by the
        centre-offset formula.

    Returns
    -------
    tuple[int, Decimal, str]
        ``(_score, frequency, _color)`` where ``_score`` is the BLOSUM
        substitution score (or a sentinel value, see below) and ``_color``
        is a hex colour string.

    Dark-green (``#219f11``) override — when and why
    -------------------------------------------------
    Synonymous substitutions carry no amino acid change but are scientifically
    meaningful (e.g. codon-usage bias, RNA secondary structure).  To distinguish
    them visually from missense mutations they are painted in dark green
    (``#219f11``) and assigned the sentinel score **+12**, which places them in
    a dedicated band at the top of the colour scale.

    The override is applied under **three mutually exclusive conditions**,
    checked in priority order:

    1. **Identical raw strings** (``old_codon_or_aa == new_codon_or_aa``)
       The most common case in codon mode: the unchanged-codon TSV contains
       entries where original and mutant codon are the same letter-for-letter
       (e.g. ``GCC → GCC``).  In aminoacids mode this catches ``A → A``.

    2. **Identical resolved strings** (``_old_codon_or_aa == _new_codon_or_aa``)
       After ``resolve_codon_or_aa()``, if the resolved labels are equal but the
       raw strings differed (e.g. raw ``'GCC'`` vs resolved ``'A'`` in aminoacids
       mode), this branch triggers.

    3. **Different codons encoding the same amino acid** (codon mode only)
       ``alt_translate(old_codon) == alt_translate(new_codon)`` catches true
       synonymous *substitutions* (e.g. ``GCC → GCT``, both encoding Ala).
       This is the case the user most likely means by "synonymous mutation".

    **Sentinel score +12**
        For the ``amino_acid_changes`` ListedColormap, ``BoundaryNorm`` maps
        +12 to index 32, where ``colors[32] == '#219f11'`` (dark green).  For
        continuous colormaps (``coolwarm_r`` etc.) the pre-resolved
        ``_color = '#219f11'`` hex string is passed directly to the scatter
        call (``c=cm_hex`` path), so the same colour appears regardless of
        what score +12 would map to in that colormap's gradient.

    **When the override is NOT applied**
        - In ``--aminoacids`` mode without ``--include-synonymous``: the entry
          is silently **skipped** in ``collect_scatter_data()`` before this
          function is ever called (line 1101 in that function).  No circle is
          drawn at all — not even in the wrong colour.
        - In ``--aminoacids`` mode **with** ``--include-synonymous``: the entry
          reaches this function, branch 1 or 2 fires, and the dark-green circle
          is rendered as a diamond marker (``'D'``) instead of the default hex.
        - In codon mode (no ``--aminoacids``): synonymous codon entries from
          both the main frequencies TSV (branch 3, e.g. GCC→GCT) and the
          unchanged-codons TSV (branch 1, e.g. GCC→GCC) are **always** drawn
          in dark green, regardless of ``--include-synonymous``.
    """

    if not codon_on_input:
        if not _old_codon_or_aa:
            raise ValueError(f"Aieee, old_codon_or_aa='{old_codon_or_aa}' is empty")
        if not _new_codon_or_aa:
            raise ValueError(f"Aieee, new_codon_or_aa='{new_codon_or_aa}' is empty")

        if myoptions.debug:
            print(
                f"Info: some single-letter but neither asterisk nor dash amino acid residue: _old_codon_or_aa='{_old_codon_or_aa}', new_codon_or_aa='{new_codon_or_aa}'")

    if _new_codon_or_aa in ('---', 'DEL', 'INS'):
        # Use the pre-computed module-level worst score (set by load_matrix())
        # as the sentinel for deletions/insertions instead of a hardcoded -11.
        # For BLOSUM80 this equals matrix['A']['*'] = -6; other matrices differ.
        if _min_theoretical_score is None:
            raise RuntimeError(
                "_min_theoretical_score is not initialised — call load_matrix() before scoring"
            )
        _score = _min_theoretical_score
    else:
        _score = get_score(myoptions, matrix, codon_on_input, _old_codon_or_aa, _new_codon_or_aa)

    if norm is not None:
        _colorindex = norm(_score)
    else:
        # Standard matplotlib cmap path (original v0.3 approach modified for dynamic symmetry):
        # Because bounds are perfectly symmetric, centre slot (index len(colors)//2) exactly represents score 0.
        _vmin = getattr(myoptions, 'cmap_vmin', -11)
        _colorindex = max(0, min(len(colors) - 1, int(_score) - _vmin))

    if old_codon_or_aa.upper() == new_codon_or_aa.upper():
        # Dark green override — Branch 1: identical raw strings.
        # Covers unchanged-codon TSV entries (e.g. GCC→GCC) and aminoacids
        # entries where original_aa == mutant_aa (A→A).
        # Sentinel score +12 maps to colors[32]='#219f11' in amino_acid_changes;
        # for continuous colormaps the pre-resolved #219f11 hex is used directly.
        _color = '#219f11'  # dark green — synonymous sentinel colour (scores[32] in amino_acid_changes)
        _score = 12
    elif _old_codon_or_aa.upper() == _new_codon_or_aa.upper():
        # Dark green override — Branch 2: identical *resolved* labels.
        # Fires when the raw strings differ but resolve to the same amino acid
        # label after resolve_codon_or_aa() (rare edge case).
        _color = '#219f11'  # dark green — synonymous sentinel colour
        _score = 12
    elif old_codon_or_aa in ('---', 'DEL', 'INS', '*') or new_codon_or_aa in ('---', 'DEL', 'INS', '*', 'TGA', 'TAA', 'TAG'):
        # DEL/INS/STOP: score is already set to _min_theoretical_score (above). No override here.
        _color = '#ff0000'  # pure red — deletion, insertion, or stop codon
    elif _old_codon_or_aa in ('X', 'NNN') or new_codon_or_aa in ('X', 'NNN'):
        _color = '#808080'  # medium gray — ambiguous/unknown codon or amino acid
    elif codon_on_input:
        if alt_translate(_old_codon_or_aa) == alt_translate(_new_codon_or_aa):
            # Dark green override — Branch 3: different codons, same amino acid.
            # True synonymous codon substitution, e.g. GCC→GCT (both Ala).
            # This is the primary case for codon-mode scatter plots.
            _color = '#219f11'  # dark green — synonymous sentinel colour
            _score = 12
        elif alt_translate(_new_codon_or_aa) == 'X':
            _color = '#808080'  # medium gray — codon translates to an ambiguous/stop residue
        else:
            if myoptions.debug:
                sys.stdout.write(f"Info: Translating {_old_codon_or_aa} to {_new_codon_or_aa} to fetch color from _colors,")
            _color = colors[_colorindex]
            if myoptions.debug:
                sys.stdout.write(f" which has yielded {_colorindex} and {str(_color)}{os.linesep}")
    else:
        if myoptions.debug:
            sys.stdout.write(f"Info: Fetching color for {_old_codon_or_aa} to {_new_codon_or_aa} directly from _colors,")
        _color = colors[_colorindex]
        if myoptions.debug:
            sys.stdout.write(f" which has yielded {_colorindex} and {str(_color)}{os.linesep}")

    return _score, frequency, _color


def adjust_size_and_color_neutralized_escape(neutralized_parent_difference, codon_on_input):
    """Assign size and colour for neutralized_parent_difference / escape_parent_difference figures.

    Used only for the ``neutralized_parent_difference`` and
    ``escape_parent_difference`` figure types, not for the standard
    mutation-frequency scatter plots.

    Colour legend
    -------------
    * ``#ff0000``  (pure red)       — negative parent difference (< −0.001);
                                      indicates reduced neutralisation / escape.
    * ``#00ff04``  (bright lime green) — positive parent difference (> +0.001);
                                        indicates enhanced neutralisation / escape.
    * ``#000000``  (black)          — near-zero difference (−0.001 … +0.001);
                                      drawn as a zero-size placeholder so the
                                      position is not dropped from the figure.
    """

    if neutralized_parent_difference < -0.001:
        _size = abs(neutralized_parent_difference)
        _color = '#ff0000'  # pure red — reduced neutralisation / escape
    elif neutralized_parent_difference > 0.001:
        if not codon_on_input:
            _size, _color = 0, '#00ff04'  # bright lime green — enhanced in aa mode (no size)
        else:
            _size = neutralized_parent_difference
            _color = '#00ff04'  # bright lime green — enhanced in codon mode
    else:
        _size = 0
        _color = '#000000'  # black — near-zero; zero-size placeholder
    return _size, _color


def adjust_size_and_color_weighted(weighted_diff_escape_neutralized, generic_circle_size=5000, weighted_diff_escape_neutralized_size1=3000, weighted_diff_escape_neutralized_size2=2000):
    """Assign size and colour for weighted_diff_escape_neutralized figures.

    Used only for the ``weighted_diff_escape_neutralized`` figure type, not
    for the standard mutation-frequency scatter plots.

    Colour legend (three-tier, largest signal first)
    -------------------------------------------------
    * ``#000000``  (black)     — very high signal (> 1); largest circles.
    * ``#00008b``  (dark blue) — medium signal (0.1 … 1).
    * ``#87ceeb``  (sky blue)  — low signal (0.01 … 0.1).
    * ``#000000``  (black)     — negligible signal (≤ 0.01); zero-size placeholder.
    """

    if weighted_diff_escape_neutralized > 1:
        _size = weighted_diff_escape_neutralized * generic_circle_size
        _color = '#000000'  # black — very high signal
    elif weighted_diff_escape_neutralized > 0.1:
        _size = weighted_diff_escape_neutralized * weighted_diff_escape_neutralized_size1
        _color = '#00008b'  # dark blue — medium signal
    elif weighted_diff_escape_neutralized > 0.01:
        _size = weighted_diff_escape_neutralized * weighted_diff_escape_neutralized_size2
        _color = '#87ceeb'  # sky blue — low signal
    else:
        _size = abs(weighted_diff_escape_neutralized) * 0
        _color = '#000000'  # black — negligible signal; zero-size placeholder
    return _size, _color


def load_matrix(myoptions):
    """Load BLOSUM or custom substitution matrix and set up the output file prefix."""
    if myoptions.matrix_file and os.path.exists(myoptions.matrix_file):
        print(myoptions.matrix)
        if not myoptions.matrix:
            myoptions.matrix = re.sub(r'.*/', '', myoptions.matrix_file)
        print(myoptions.matrix)
        _matrix_type, _matrix_num = re.sub(r'\d+', '', myoptions.matrix), int(re.sub(r'[a-zA-Z]+', '', myoptions.matrix))
        _matrix = blosum.BLOSUM(myoptions.matrix_file)
        _matrix_name = myoptions.matrix_file.split(os.path.sep)[-1]
        myoptions.matrix = _matrix_name
    else:
        try:
            _matrix_type = re.sub(r'\d+', '', myoptions.matrix)
            _matrix_num = int(re.sub(r'[a-zA-Z_]+', '', myoptions.matrix))
        except ValueError:
            sys.stderr.write(
                f"Warning: Cannot parse '{myoptions.matrix}' as a BLOSUM-style matrix name "
                f"(expected e.g. BLOSUM80); falling back to BLOSUM80 for scoring "
                f"but keeping '{myoptions.matrix}' in the output filename.{os.linesep}"
            )
            _matrix = blosum.BLOSUM(80)
            _matrix_name = myoptions.matrix
        else:
            if _matrix_type == 'BLOSUM':
                _matrix = blosum.BLOSUM(_matrix_num)
                _matrix_name = f"BLOSUM{_matrix_num}"
            else:
                sys.stderr.write(f"Warning: Unexpected matrix type {str(_matrix_type)}, falling back to BLOSUM{os.linesep}")
                _matrix = blosum.BLOSUM(_matrix_num)
                _matrix_name = f"BLOSUM{_matrix_num}"
                myoptions.matrix = _matrix_name

    if not myoptions.outfile_prefix:
        raise RuntimeError("Please provide output filename prefix via --outfile-prefix")

    # Clean up the output prefix by stripping common extensions if user provided a full filename
    _clean_prefix = myoptions.outfile_prefix
    while True:
        _changed = False
        for _ext in ['.frequencies.tsv', '.frequencies', '.tsv', '.png', '.jpg', '.pdf', '.html']:
            if _clean_prefix.endswith(_ext):
                _clean_prefix = _clean_prefix[:-len(_ext)]
                _changed = True
                break
        if not _changed:
            break

    _scaling_suffix = 'linear_scaling' if myoptions.linear_circle_size else 'area_scaling'
    _outfile_prefix = _clean_prefix + '.' + _matrix_name + '.' + _scaling_suffix + '.' + myoptions.colormap
    print(f"Info: _outfile_prefix={_outfile_prefix}")

    _theoretical_scores = set()
    for _scores_dict in _matrix.values():
        for _score in _scores_dict.values():
            _theoretical_scores.add(_score)
    print(f"Info: Using {_matrix_name} matrix now. Theoretical minimum score is {min(_theoretical_scores)}, theoretical maximum score is {max(_theoretical_scores)}, values are {str(_theoretical_scores)}")
    # Publish to the module-level sentinel so get_score() / adjust_size_and_color()
    # can use it without a matrix dict lookup.
    global _min_theoretical_score  # pylint: disable=global-statement
    _min_theoretical_score = int(min(_theoretical_scores))
    _max_theoretical_score = int(max(_theoretical_scores))

    _bound_abs = max(abs(_min_theoretical_score), abs(_max_theoretical_score))
    myoptions.matrix_min_theoretical = -_bound_abs
    myoptions.matrix_max_theoretical = _bound_abs

    return _matrix, _matrix_name, _min_theoretical_score, _max_theoretical_score, _outfile_prefix
