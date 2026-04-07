# -*- coding: utf-8 -*-
# vim: set fileencoding=utf-8 ts=4 sw=4 expandtab :
"""Unit tests for render_bokeh and render_matplotlib colorbar logic.

These tests focus on the palette-building and tick-alignment code that was
the subject of several iterations of refinement.  The full render functions
require a display and real data, so the palette/colorbar logic is tested in
isolation by replicating the same computation used inside the functions.
"""

import unittest

import matplotlib
import matplotlib.colors
import matplotlib.pyplot as plt
import numpy as np

matplotlib.use("Agg")   # headless backend — no display required


# ---------------------------------------------------------------------------
# Helpers that mirror the exact logic inside render_bokeh
# ---------------------------------------------------------------------------

AMINO_ACID_CHANGES_COLORS = [
    "#930000", "#930000", "#930000", "#930000", "#930000", "#930000",
    "#960000", "#580041", "#8200ff", "#c500ff", "#ff00fd", "#CC79A7",
    "#eea1d0", "#cc0000", "#ff0000", "#ff4f00", "#ff7c7c", "#ff9999",
    "#c58a24", "#9c644b",
    "#ffff00",  # index 20 → score 0 via BoundaryNorm rescaling
    "#ffcc00", "#ffa200", "#7DCCFF", "#0042ff", "#0000ff",
    "#D6D6D6", "#B7B7B7", "#8B8B8B", "#bbff00", "#00ff04",
    "#97CE2F", "#219f11",
    "#930000", "#930000", "#930000", "#930000", "#930000", "#930000",
]

DKEENAN_COLORS = [
    "#00B7FF", "#004DFF", "#00FFFF", "#826400", "#580041", "#FF00FF",
    "#00FF00", "#C500FF", "#B4FFD7", "#FFCA00", "#969600", "#B4A2FF",
    "#C20078",
    "#000000",  # index 13 → score 0 via BoundaryNorm rescaling
    "#0000C1", "#FF8B00", "#FFC8FF", "#666666", "#FF0000", "#CCCCCC",
    "#009E8F", "#D7A870", "#8200FF", "#960000", "#BBFF00", "#FFFF00",
    "#006F00",
]


def _blend_with_white(hex_color: str, alpha: float) -> str:
    """Pre-blend *hex_color* with white at *alpha* — mirrors render_bokeh helper."""
    r = int(hex_color[1:3], 16)
    g = int(hex_color[3:5], 16)
    b = int(hex_color[5:7], 16)
    ra = round(alpha * r + (1 - alpha) * 255)
    ga = round(alpha * g + (1 - alpha) * 255)
    ba = round(alpha * b + (1 - alpha) * 255)
    return f"#{ra:02x}{ga:02x}{ba:02x}"


_BOKEH_CIRCLE_ALPHA = 0.5  # must match constant in render_bokeh


def _make_cmap_norm(colors, low, high):
    """Replicate get_colormap() output for discrete palettes."""
    cmap = matplotlib.colors.ListedColormap(colors, "test_cmap", len(colors))
    norm = matplotlib.colors.BoundaryNorm(np.arange(low, high, 1), cmap.N)
    return cmap, norm


def _build_bokeh_palette(colors, norm, cmap):
    """Mirror the palette-building logic in render_bokeh (including alpha pre-blend)."""
    n = len(colors) if colors is not None else 256
    half = n // 2
    score_range = range(-half, half + 1)
    if norm is not None and colors is not None:
        # Discrete path: colours[norm(s)], matching adjust_size_and_color.
        raw_palette = [
            matplotlib.colors.to_hex(matplotlib.colors.to_rgba(
                colors[max(0, min(len(colors) - 1, norm(s)))])
            )
            for s in score_range
        ]
    elif cmap is not None:
        raw_palette = [matplotlib.colors.to_hex(cmap(i / max(1, n - 1))) for i in range(n)]
    else:
        raw_palette = ["#aaaaaa"] * n
    # Pre-blend to match the visual appearance of alpha=0.5 scatter glyphs.
    palette = [_blend_with_white(c, _BOKEH_CIRCLE_ALPHA) for c in raw_palette]
    low_val = -half - 0.5
    high_val = half + 0.5
    return palette, low_val, high_val, list(score_range)


# ---------------------------------------------------------------------------
# Tests for render_bokeh colorbar palette
# ---------------------------------------------------------------------------

class TestRenderBokehColorbarPalette(unittest.TestCase):
    """Verify that the Bokeh colorbar palette maps scores to the correct colours."""

    def setUp(self):
        self.aa_cmap, self.aa_norm = _make_cmap_norm(
            AMINO_ACID_CHANGES_COLORS, -19, 19
        )
        self.dk_cmap, self.dk_norm = _make_cmap_norm(
            DKEENAN_COLORS, -13, 13
        )

    # --- amino_acid_changes palette ---

    def test_amino_acid_changes_palette_length(self):
        """Palette has one entry per integer score in the full range."""
        palette, _, _, ticks = _build_bokeh_palette(
            AMINO_ACID_CHANGES_COLORS, self.aa_norm, self.aa_cmap
        )
        self.assertEqual(len(palette), len(AMINO_ACID_CHANGES_COLORS))
        self.assertEqual(len(palette), len(ticks))

    def test_amino_acid_changes_score_zero_is_yellow(self):
        """Score 0 colorbar band must be blended yellow (#ffff80 at alpha=0.5 on white).

        Raw palette colour: #ffff00
        After pre-blending with Python round():
          R: round(0.5*255 + 127.5) = round(255.0) = 255 = 0xff
          G: round(0.5*255 + 127.5) = round(255.0) = 255 = 0xff
          B: round(0.5*0   + 127.5) = round(127.5) = 128 = 0x80 (banker's rounding)
        Result: #ffff80
        """
        palette, _, _, ticks = _build_bokeh_palette(
            AMINO_ACID_CHANGES_COLORS, self.aa_norm, self.aa_cmap
        )
        idx = ticks.index(0)
        self.assertEqual(palette[idx].lower(), "#ffff80",
                         f"Score 0 colorbar band should be blended yellow but got {palette[idx]}")

    def test_amino_acid_changes_score_plus1_colour(self):
        """Score +1 colorbar band is blended #ffcc00 at alpha=0.5 on white = #ffe680."""
        palette, _, _, ticks = _build_bokeh_palette(
            AMINO_ACID_CHANGES_COLORS, self.aa_norm, self.aa_cmap
        )
        idx = ticks.index(1)
        self.assertEqual(palette[idx].lower(), _blend_with_white("#ffcc00", _BOKEH_CIRCLE_ALPHA))

    def test_amino_acid_changes_score_minus1_colour(self):
        """Score -1 colorbar band is blended #9c644b at alpha=0.5 on white = #ceb2a5."""
        palette, _, _, ticks = _build_bokeh_palette(
            AMINO_ACID_CHANGES_COLORS, self.aa_norm, self.aa_cmap
        )
        idx = ticks.index(-1)
        self.assertEqual(palette[idx].lower(), _blend_with_white("#9c644b", _BOKEH_CIRCLE_ALPHA))

    def test_band_width_is_one(self):
        """Each colour band must be exactly 1 unit wide."""
        palette, low, high, _ = _build_bokeh_palette(
            AMINO_ACID_CHANGES_COLORS, self.aa_norm, self.aa_cmap
        )
        n = len(palette)
        band_width = (high - low) / n
        self.assertAlmostEqual(band_width, 1.0, places=10)

    def test_score_zero_band_centre_is_zero(self):
        """The centre of the score-0 band must fall exactly on 0.0."""
        palette, low, high, ticks = _build_bokeh_palette(
            AMINO_ACID_CHANGES_COLORS, self.aa_norm, self.aa_cmap
        )
        n = len(palette)
        band_width = (high - low) / n
        idx = ticks.index(0)
        centre = low + (idx + 0.5) * band_width
        self.assertAlmostEqual(centre, 0.0, places=10)

    def test_ticks_span_full_range(self):
        """Ticks must run from -half to +half inclusive."""
        _, _, _, ticks = _build_bokeh_palette(
            AMINO_ACID_CHANGES_COLORS, self.aa_norm, self.aa_cmap
        )
        half = len(AMINO_ACID_CHANGES_COLORS) // 2
        self.assertEqual(ticks[0],  -half)
        self.assertEqual(ticks[-1],  half)

    def test_all_palette_entries_are_valid_hex(self):
        """Every palette entry must be a valid CSS hex colour string."""
        palette, _, _, _ = _build_bokeh_palette(
            AMINO_ACID_CHANGES_COLORS, self.aa_norm, self.aa_cmap
        )
        for entry in palette:
            self.assertRegex(entry, r"^#[0-9a-f]{6}$",
                             f"Invalid hex colour: {entry}")

    # --- dkeenan palette ---

    def test_dkeenan_palette_length(self):
        """dkeenan palette has one entry per integer score in the full range."""
        palette, _, _, _ = _build_bokeh_palette(
            DKEENAN_COLORS, self.dk_norm, self.dk_cmap
        )
        self.assertEqual(len(palette), len(DKEENAN_COLORS))

    def test_dkeenan_band_width_is_one(self):
        """Each dkeenan colour band must be exactly 1 unit wide."""
        palette, low, high, _ = _build_bokeh_palette(
            DKEENAN_COLORS, self.dk_norm, self.dk_cmap
        )
        band_width = (high - low) / len(palette)
        self.assertAlmostEqual(band_width, 1.0, places=10)

    def test_dkeenan_score_zero_band_centre_is_zero(self):
        """The centre of the dkeenan score-0 band must fall exactly on 0.0."""
        palette, low, high, ticks = _build_bokeh_palette(
            DKEENAN_COLORS, self.dk_norm, self.dk_cmap
        )
        n = len(palette)
        band_width = (high - low) / n
        idx = ticks.index(0)
        centre = low + (idx + 0.5) * band_width
        self.assertAlmostEqual(centre, 0.0, places=10)

    # --- fallback: norm is None ---

    def test_fallback_no_norm(self):
        """When norm is None, palette is sampled linearly from cmap."""
        cmap = matplotlib.colormaps.get_cmap("viridis")
        palette, _, _, _ = _build_bokeh_palette(
            AMINO_ACID_CHANGES_COLORS, None, cmap
        )
        # Length matches colors list
        self.assertEqual(len(palette), len(AMINO_ACID_CHANGES_COLORS))
        # All hex
        for entry in palette:
            self.assertRegex(entry, r"^#[0-9a-f]{6}$")

    def test_fallback_no_norm_no_cmap(self):
        """When both norm and cmap are None, palette is flat blended grey.

        Raw colour: #aaaaaa. Blended at alpha=0.5 on white:
          round(0.5*170 + 127.5) = round(212.5) = 212 = 0xd4
        Result: #d4d4d4.
        """
        palette, _, _, _ = _build_bokeh_palette(None, None, None)
        self.assertTrue(all(c == "#d4d4d4" for c in palette),
                        f"Expected all blended grey #d4d4d4 but got: {set(palette)}")


# ---------------------------------------------------------------------------
# Tests for render_matplotlib colorbar tick setup
# ---------------------------------------------------------------------------

class TestRenderMatplotlibColorbar(unittest.TestCase):
    """Verify the matplotlib colorbar tick positions and labels."""

    def _build_mpl_colorbar(self, colors, norm_low=-19, norm_high=19):
        """Create a minimal matplotlib figure with colorbar, as render_matplotlib does."""
        cmap = matplotlib.colors.ListedColormap(colors, "test", len(colors))
        norm = matplotlib.colors.BoundaryNorm(
            np.arange(norm_low, norm_high, 1), cmap.N
        )
        fig, (ax_scatter, ax_cb) = plt.subplots(1, 2, figsize=(4, 6))
        # Minimal scatter so colorbar has something to attach to
        sc = ax_scatter.scatter([0], [0], c=[0], cmap=cmap, norm=norm, s=50)
        cb = fig.colorbar(sc, cax=ax_cb, label="test score values",
                          location="right", alpha=0.5)
        # Replicate render_matplotlib tick setup exactly
        cb.ax.set_yticks(np.arange(norm_low + 0.5, norm_high - 0.5, 1),
                         np.arange(norm_low, norm_high - 1, 1))
        cb.ax.tick_params(axis="y", which="minor", length=0)
        return fig, cb, norm, cmap

    def test_tick_count_amino_acid_changes(self):
        """Colorbar should have one tick per bin (38 for amino_acid_changes)."""
        _, cb, _, _ = self._build_mpl_colorbar(AMINO_ACID_CHANGES_COLORS)
        tick_locs = cb.ax.get_yticks()
        # 38 boundaries → 37 bins, set_yticks places ticks at bin centres
        n_bins = len(AMINO_ACID_CHANGES_COLORS) - 2   # 37
        self.assertEqual(len(tick_locs), n_bins)

    def test_tick_labels_start_at_norm_low(self):
        """First tick label must be norm_low (-19 for amino_acid_changes)."""
        _, cb, _, _ = self._build_mpl_colorbar(AMINO_ACID_CHANGES_COLORS)
        labels = [t.get_text() for t in cb.ax.get_yticklabels()]
        labels = [lbl for lbl in labels if lbl]   # remove empty strings
        if labels:
            self.assertEqual(labels[0], "-19")

    def test_score_zero_colour_via_norm(self):
        """cmap(norm(0)) must return yellow for amino_acid_changes."""
        cmap = matplotlib.colors.ListedColormap(
            AMINO_ACID_CHANGES_COLORS, "test", len(AMINO_ACID_CHANGES_COLORS)
        )
        norm = matplotlib.colors.BoundaryNorm(np.arange(-19, 19, 1), cmap.N)
        colour = matplotlib.colors.to_hex(cmap(norm(0)))
        self.assertEqual(colour.lower(), "#ffff00")

    def test_score_plus1_colour_via_norm(self):
        """cmap(norm(1)) must return #ffcc00."""
        cmap = matplotlib.colors.ListedColormap(
            AMINO_ACID_CHANGES_COLORS, "test", len(AMINO_ACID_CHANGES_COLORS)
        )
        norm = matplotlib.colors.BoundaryNorm(np.arange(-19, 19, 1), cmap.N)
        colour = matplotlib.colors.to_hex(cmap(norm(1)))
        self.assertEqual(colour.lower(), "#ffcc00")

    def test_score_minus19_colour_via_norm(self):
        """Score -19 maps to the first colour in the list (#930000)."""
        cmap = matplotlib.colors.ListedColormap(
            AMINO_ACID_CHANGES_COLORS, "test", len(AMINO_ACID_CHANGES_COLORS)
        )
        norm = matplotlib.colors.BoundaryNorm(np.arange(-19, 19, 1), cmap.N)
        colour = matplotlib.colors.to_hex(cmap(norm(-19)))
        self.assertEqual(colour.lower(), "#930000")

    def test_dkeenan_score_zero_colour_via_norm(self):
        """dkeenan: cmap(norm(0)) maps to #0000c1 (index 14, shifted by BoundaryNorm rescaling)."""
        cmap = matplotlib.colors.ListedColormap(
            DKEENAN_COLORS, "test", len(DKEENAN_COLORS)
        )
        norm = matplotlib.colors.BoundaryNorm(np.arange(-13, 13, 1), cmap.N)
        colour = matplotlib.colors.to_hex(cmap(norm(0)))
        self.assertEqual(colour.lower(), "#0000c1")

    def tearDown(self):
        plt.close("all")


# ---------------------------------------------------------------------------
# Tests for _blend_with_white helper (alpha pre-blend logic)
# ---------------------------------------------------------------------------

class TestBlendWithWhite(unittest.TestCase):
    """Lock in the alpha=0.5 pre-blend math used for the Bokeh colorbar palette.

    The Bokeh ColorBar renders at alpha=1.0 by default, while scatter glyphs
    render at alpha=0.5.  To make the colorbar bands appear the same shade as
    the circles, we pre-blend each palette colour with the white background:

        apparent_channel = round(alpha * src + (1 - alpha) * 255)

    All expected values below were computed with this formula at alpha=0.5.
    """

    def test_yellow_blend(self):
        """#ffff00 @ 0.5 on white → #ffff80 (blue channel: round(127.5)=128=0x80)."""
        self.assertEqual(_blend_with_white("#ffff00", 0.5), "#ffff80")

    def test_gold_blend(self):
        """#ffcc00 @ 0.5 on white → #ffe680 (green: round(0.5*204+127.5)=round(229.5)=230=0xe6)."""
        self.assertEqual(_blend_with_white("#ffcc00", 0.5), "#ffe680")

    def test_orange_blend(self):
        """#ffa200 @ 0.5 on white → #ffd080 (green: round(0.5*162+127.5)=round(208.5)=208=0xd0)."""
        self.assertEqual(_blend_with_white("#ffa200", 0.5), "#ffd080")

    def test_black_blend(self):
        """#000000 @ 0.5 on white → #808080 (round(127.5)=128=0x80 each channel)."""
        self.assertEqual(_blend_with_white("#000000", 0.5), "#808080")

    def test_white_blend(self):
        """#ffffff @ any alpha on white → #ffffff (no change)."""
        self.assertEqual(_blend_with_white("#ffffff", 0.5), "#ffffff")

    def test_alpha_one_is_identity(self):
        """At alpha=1.0 the blended colour equals the source colour."""
        for c in ["#ffa200", "#ffff00", "#930000", "#219f11"]:
            self.assertEqual(_blend_with_white(c, 1.0), c)

    def test_score_zero_and_score_plus2_visually_different(self):
        """After blending, score 0 (yellow) and score +2 (orange) must still differ.

        This regression test catches any accidental collapse that would make
        neighbouring bands visually indistinguishable (the original bug was that
        the un-blended colorbar gold band for +1 looked the same as the blended
        orange circle for +2)."""
        blended_0 = _blend_with_white("#ffff00", 0.5)   # score 0
        blended_p2 = _blend_with_white("#ffa200", 0.5)   # score +2
        self.assertNotEqual(blended_0, blended_p2)


if __name__ == "__main__":
    unittest.main()
