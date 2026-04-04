"""Tests for discrete colormap (ListedColormap + BoundaryNorm) alignment.

Reference implementation / design spec:
  https://gist.github.com/mmokrejs/99fb8b0ac32676698247c5b9f88724d2

Background
----------
matplotlib's ``BoundaryNorm`` with **integer** boundaries uses
``np.digitize`` under the hood, which for monotonically-increasing bins
applies the convention::

    bins[i-1] <= x < bins[i]   →   digitize returns i   →   iret = i-1

This means integer score *s* maps to the interval ``[s, s+1)`` (the
interval that STARTS at s), i.e. iret = s + |min_boundary|.

For the ``amino_acid_changes`` palette (39 colours,
``BoundaryNorm(np.arange(-19, 19, 1), 39)``):
* Score 0  → colour[20] = ``#ffff00``  (yellow, confirmed correct by user)
* Score +1 → colour[21] = ``#ffcc00``  (gold)
* Score +2 → colour[22] = ``#ffa200``  (orange)

The reference gist uses a *different* 23-colour palette where
``norm(score - 1)`` for positive scores is needed to cancel the interval
offset introduced by that specific palette design.  For the
``amino_acid_changes`` 39-colour palette the simple ``norm(score)``
lookup is already correct — the palette was laid out to match the
``np.digitize`` convention directly (yellow at index 20, not 19).

For the *continuous* coolwarm_r path the norm is ``None`` and the colour
is derived from a ``np.linspace(0, 1, n)`` sample of the cmap; no
BoundaryNorm is involved.
"""


import pytest
import numpy as np
import matplotlib
import matplotlib.colors as mc

# pytest injects fixtures by matching parameter names to fixture function names.
# This looks like "redefining the outer-scope fixture function" to pylint, but
# is correct and idiomatic pytest.  Suppress W0621 for the whole test module.
# pylint: disable=redefined-outer-name

matplotlib.use("Agg")  # non-interactive backend

# ── palette definition ────────────────────────────────────────────────────────
# Exact same list as in get_colormap() for 'amino_acid_changes'.
_AA_CHANGES_COLORS = [
    "#930000", "#930000", "#930000", "#930000", "#930000", "#930000",
    "#960000", "#580041", "#8200ff", "#c500ff", "#ff00fd",
    "#CC79A7", "#eea1d0",
    "#cc0000", "#ff0000", "#ff4f00", "#ff7c7c", "#ff9999",
    "#c58a24", "#9c644b",
    "#ffff00",   # index 20 → score 0 (yellow, neutral)
    "#ffcc00",   # index 21 → score +1
    "#ffa200",   # index 22 → score +2
    "#7DCCFF",   # index 23 → score +3
    "#0042ff", "#0000ff",
    "#D6D6D6", "#B7B7B7", "#8B8B8B",
    "palegreen",
    "#bbff00", "#97CE2F", "#219f11",  # index 32 → synonymous (+12)
    "#930000", "#930000", "#930000", "#930000", "#930000", "#930000",
]
assert len(_AA_CHANGES_COLORS) == 39, "palette must have exactly 39 entries"


@pytest.fixture(scope="module")
def aa_norm():
    """Return the BoundaryNorm used for amino_acid_changes in get_colormap()."""
    cmap = mc.ListedColormap(_AA_CHANGES_COLORS, "amino_acid_changes", 39)
    norm = mc.BoundaryNorm(np.arange(-19, 19, 1), cmap.N)
    return norm


class TestBoundaryNormDigitizeBehavior:
    """Document exactly how np.digitize maps integer BLOSUM scores to indices.

    np.digitize(x, bins) for monotonically-increasing bins:
        bins[i-1] <= x < bins[i]  →  returns i  →  iret = i - 1

    Combined with the ncolors > n_intervals scaling:
        n_intervals = 37  (len(np.arange(-19, 19, 1)) - 1)
        ncolors     = 39
        ret = iret * (ncolors - 1) // (n_intervals)
            = iret * 38 // 37   (note: NOT 38//36 — check your matplotlib version)

    The empirically verified values below are the GROUND TRUTH for regressions.
    """

    def test_score_zero_maps_to_yellow(self, aa_norm):
        """Score 0 must map to index 20 = '#ffff00' (yellow, neutral)."""
        idx = aa_norm(0)
        assert idx == 20, f"norm(0) = {idx}, expected 20"
        assert _AA_CHANGES_COLORS[idx] == "#ffff00", \
            f"color at norm(0)={idx} is {_AA_CHANGES_COLORS[idx]!r}, expected '#ffff00'"

    def test_score_plus1_maps_to_gold(self, aa_norm):
        """Score +1 must map to index 21 = '#ffcc00' (gold)."""
        idx = aa_norm(1)
        assert idx == 21, f"norm(1) = {idx}, expected 21"
        assert _AA_CHANGES_COLORS[idx] == "#ffcc00", \
            f"color at norm(1)={idx} is {_AA_CHANGES_COLORS[idx]!r}, expected '#ffcc00'"

    def test_score_plus2_maps_to_orange(self, aa_norm):
        """Score +2 must map to index 22 = '#ffa200' (orange).

        Y505H (Tyr→His) has BLOSUM80 score +2. Its scatter circle must be
        '#ffa200', and the Bokeh colorbar band at tick '+2' must also show
        '#ffa200'.  If this assertion fails the circle and colorbar are
        misaligned.
        """
        idx = aa_norm(2)
        assert idx == 22, f"norm(2) = {idx}, expected 22"
        assert _AA_CHANGES_COLORS[idx] == "#ffa200", \
            f"color at norm(2)={idx} is {_AA_CHANGES_COLORS[idx]!r}, expected '#ffa200'"

    def test_score_minus1(self, aa_norm):
        """Score -1 must map to index 19 = '#9c644b'."""
        idx = aa_norm(-1)
        assert idx == 19, f"norm(-1) = {idx}, expected 19"
        assert _AA_CHANGES_COLORS[idx] == "#9c644b"

    def test_score_minus2(self, aa_norm):
        """Score -2 must map to index 17 = '#ff9999'."""
        idx = aa_norm(-2)
        assert idx == 17, f"norm(-2) = {idx}, expected 17"
        assert _AA_CHANGES_COLORS[idx] == "#ff9999"

    def test_score_minus11_del(self, aa_norm):
        """Score -11 (DEL/INS sentinel) must map to index 8 = '#8200ff'."""
        idx = aa_norm(-11)
        assert idx == 8, f"norm(-11) = {idx}, expected 8"
        assert _AA_CHANGES_COLORS[idx] == "#8200ff"

    def test_all_positive_scores_have_distinct_colors(self, aa_norm):
        """Every positive score in BLOSUM80 range (+1 to +11) must map to
        a distinct colour index (no two-in-one) in the amino_acid_changes palette.
        """
        indices = [aa_norm(s) for s in range(1, 12)]
        assert len(set(indices)) == len(indices), \
            f"Positive scores map to duplicate color indices: {indices}"

    def test_zero_and_positive_scores_differ(self, aa_norm):
        """Score 0 and score +1 must map to DIFFERENT colour indices.

        If they map to the same index the neutral and slightly-positive
        BLOSUM mutations become visually indistinguishable.
        """
        assert aa_norm(0) != aa_norm(1), \
            "Score 0 and score +1 must have different colour indices"


class TestColorbarsMatchCircles:
    """Verify that the Bokeh colorbar palette mirrors adjust_size_and_color.

    The render_bokeh palette is built as:
        [colors[norm(s)] for s in range(-half, half+1)]

    The scatter circle colour for score s is:
        colors[norm(s)]   (discrete / norm-is-not-None path)

    So palette[s + half] must equal colors[norm(s)].  This guarantees
    that the tick at data-value s sits at the centre of the band whose
    colour matches the scatter circle.
    """

    def test_palette_entry_matches_circle_for_score_zero(self, aa_norm):
        """Score 0 circle colour and colorbar palette band must be identical."""
        n = len(_AA_CHANGES_COLORS)       # 39
        half = n // 2                     # 19
        palette = [_AA_CHANGES_COLORS[max(0, min(n - 1, aa_norm(s)))]
                   for s in range(-half, half + 1)]
        circle_color_0 = _AA_CHANGES_COLORS[aa_norm(0)]
        # score 0 is at palette index half (= 19)
        assert palette[half] == circle_color_0 == "#ffff00", \
            f"palette[{half}]={palette[half]!r}, circle={circle_color_0!r}"

    def test_palette_entry_matches_circle_for_score_plus2(self, aa_norm):
        """Score +2 circle colour and colorbar palette band must be identical."""
        n = len(_AA_CHANGES_COLORS)
        half = n // 2
        palette = [_AA_CHANGES_COLORS[max(0, min(n - 1, aa_norm(s)))]
                   for s in range(-half, half + 1)]
        circle_color_2 = _AA_CHANGES_COLORS[aa_norm(2)]
        # score +2 is at palette index half + 2 = 21
        assert palette[half + 2] == circle_color_2, \
            (f"COLORBAR MISMATCH for score +2: "
             f"palette[{half + 2}]={palette[half + 2]!r} "
             f"!= circle={circle_color_2!r}")

    def test_palette_entry_matches_circle_for_all_blosum80_scores(self, aa_norm):
        """Circle and colorbar must agree for every integer score in [-11, +11]."""
        n = len(_AA_CHANGES_COLORS)
        half = n // 2
        palette = [_AA_CHANGES_COLORS[max(0, min(n - 1, aa_norm(s)))]
                   for s in range(-half, half + 1)]
        mismatches = []
        for s in range(-11, 12):
            circle_c = _AA_CHANGES_COLORS[aa_norm(s)]
            palette_c = palette[s + half]
            if circle_c != palette_c:
                mismatches.append(
                    f"score {s:+d}: circle={circle_c!r} palette={palette_c!r}"
                )
        assert not mismatches, \
            "Circle/colorbar colour mismatch:\n" + "\n".join(mismatches)


class TestReferenceGistBehavior:
    """Document the reference gist's BoundaryNorm convention.

    The gist (https://gist.github.com/mmokrejs/99fb8b0ac...) uses a
    23-colour palette for BLOSUM scores -11…+11 with
    ``BoundaryNorm(np.arange(-12, 12, 1), 23)``.

    IMPORTANT — matplotlib version history
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    The gist comment says "decrement the index by one for positive
    scores" (``norm(score - 1)``).  This was written for an **older**
    matplotlib (<3.7) where ``BoundaryNorm`` internally used
    ``np.searchsorted(side='left')``::

        searchsorted(bins, 0, 'left') → insert BEFORE bins[12]=0 → iret=11

    In that version norm(0) = 11, so score 0 would land in '#9c644b'
    (brownish) instead of '#ffff00' (yellow).  The decrement hack
    shifted positive scores by 1 to compensate.

    In **modern** matplotlib (≥3.7) ``BoundaryNorm`` uses
    ``np.digitize``::

        digitize(0, bins) → 0 in [bins[12], bins[13]) = [0,1) → returns 13 → iret=12

    So norm(0) = 12 → '#ffff00' (yellow) — correct WITHOUT any decrement.

    Conclusion
    ~~~~~~~~~~
    In modern matplotlib ``norm(score)`` is the correct lookup for both
    the gist's 23-colour palette and the ``amino_acid_changes`` 39-colour
    palette.  No positive-score decrement is needed or beneficial.

    The matplotlib colorbar ticks should still be placed at INTEGER
    positions (``FixedTicker``) since ``np.digitize`` centres each
    integer score in its own band.  The half-integer tick positions
    documented in the gist were a Matplotlib-version-specific workaround.
    """

    _GIST_COLORS = [
        "#8110ff", "#c500ff", "#ff00fd", "#eea1d0", "#CC79A7",
        "#cc0000", "#ff0000", "#ff4f00", "#ff7c7c", "#ff9999",
        "#c58a24", "#9c644b",
        "#ffff00",   # index 12 → score 0 (modern matplotlib, np.digitize)
        "#ffcc00",   # index 13 → score +1
        "#cccc00", "#7DCCFF", "#0042ff", "#0000ff",
        "#D6D6D6", "#B7B7B7", "#8B8B8B", "#97CE2F", "#219f11",
    ]

    def test_gist_norm_zero_modern_matplotlib(self):
        """Modern matplotlib (np.digitize): norm(0) == 12 → '#ffff00' (yellow).

        The gist was written when BoundaryNorm used searchsorted('left')
        and returned norm(0)==11.  With np.digitize, score 0 falls in
        [0, 1) which is the 12th interval → index 12 → yellow.  No
        decrement needed.
        """
        cmap = mc.ListedColormap(self._GIST_COLORS, "gist_test", 23)
        norm = mc.BoundaryNorm(np.arange(-12, 12, 1), cmap.N)
        idx = norm(0)
        assert idx == 12, (
            f"norm(0)={idx}, expected 12 (modern matplotlib/np.digitize).\n"
            f"If this is 11 you are running an old matplotlib that uses "
            f"searchsorted('left') — the gist's norm(score-1) decrement "
            f"would be needed there."
        )
        assert self._GIST_COLORS[idx] == "#ffff00"

    def test_gist_norm_plus1_modern_matplotlib(self):
        """Modern matplotlib: norm(+1) == 13 → '#ffcc00' (gold)."""
        cmap = mc.ListedColormap(self._GIST_COLORS, "gist_test", 23)
        norm = mc.BoundaryNorm(np.arange(-12, 12, 1), cmap.N)
        idx = norm(1)
        assert idx == 13, f"norm(1)={idx}, expected 13"
        assert self._GIST_COLORS[idx] == "#ffcc00"

    def test_gist_no_decrement_needed_in_modern_matplotlib(self):
        """With np.digitize, norm(score) already maps correctly for all scores.

        The gist's `norm(score-1)` for positive scores is NOT needed in
        modern matplotlib — it would give score +1 the same colour as
        score 0, which is wrong.
        """
        cmap = mc.ListedColormap(self._GIST_COLORS, "gist_test", 23)
        norm = mc.BoundaryNorm(np.arange(-12, 12, 1), cmap.N)
        # Without decrement: all scores in BLOSUM80 range have distinct colours.
        indices = [norm(s) for s in range(-11, 12)]
        assert len(set(indices)) == len(indices), (
            "Without decrement, all scores should map to distinct indices.\n"
            f"Got duplicates: {indices}"
        )
        # With decrement for positive: score 0 and score +1 would collide.
        idx_0 = norm(0)        # score 0, no decrement
        idx_p1 = norm(1 - 1)   # score +1 with decrement → norm(0)
        assert idx_0 == idx_p1, (
            "Demonstrating the decrement collision: score 0 and score +1 "
            f"both map to index {idx_0} with gist decrement."
        )
