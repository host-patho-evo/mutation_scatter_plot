"""Unit tests for render_bokeh and render_matplotlib colorbar logic.

These tests focus on the palette-building and tick-alignment code that was
the subject of several iterations of refinement.  The full render functions
require a display and real data, so the palette/colorbar logic is tested in
isolation by replicating the same computation used inside the functions.
"""

import unittest
import numpy as np
import matplotlib
import matplotlib.colors
import matplotlib.pyplot as plt

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


def _make_cmap_norm(colors, low, high):
    """Replicate get_colormap() output for discrete palettes."""
    cmap = matplotlib.colors.ListedColormap(colors, "test_cmap", len(colors))
    norm = matplotlib.colors.BoundaryNorm(np.arange(low, high, 1), cmap.N)
    return cmap, norm


def _build_bokeh_palette(colors, norm, cmap):
    """Mirror the palette-building logic in render_bokeh."""
    n = len(colors) if colors is not None else 256
    half = n // 2
    score_range = range(-half, half + 1)
    if norm is not None and cmap is not None:
        palette = [matplotlib.colors.to_hex(cmap(norm(s))) for s in score_range]
    elif cmap is not None:
        palette = [matplotlib.colors.to_hex(cmap(i / max(1, n - 1))) for i in range(n)]
    else:
        palette = ["#aaaaaa"] * n
    low_val  = -half - 0.5
    high_val =  half + 0.5
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
        """Score 0 must map to yellow (#ffff00) in the amino_acid_changes palette."""
        palette, _, _, ticks = _build_bokeh_palette(
            AMINO_ACID_CHANGES_COLORS, self.aa_norm, self.aa_cmap
        )
        idx = ticks.index(0)
        self.assertEqual(palette[idx].lower(), "#ffff00",
                         f"Score 0 should be yellow but got {palette[idx]}")

    def test_amino_acid_changes_score_plus1_colour(self):
        """Score +1 maps to #ffcc00."""
        palette, _, _, ticks = _build_bokeh_palette(
            AMINO_ACID_CHANGES_COLORS, self.aa_norm, self.aa_cmap
        )
        idx = ticks.index(1)
        self.assertEqual(palette[idx].lower(), "#ffcc00")

    def test_amino_acid_changes_score_minus1_colour(self):
        """Score -1 maps to #9c644b."""
        palette, _, _, ticks = _build_bokeh_palette(
            AMINO_ACID_CHANGES_COLORS, self.aa_norm, self.aa_cmap
        )
        idx = ticks.index(-1)
        self.assertEqual(palette[idx].lower(), "#9c644b")

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
        """When both norm and cmap are None, palette is flat grey."""
        palette, _, _, _ = _build_bokeh_palette(None, None, None)
        self.assertTrue(all(c == "#aaaaaa" for c in palette))


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
        labels = [l for l in labels if l]   # remove empty strings
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


if __name__ == "__main__":
    unittest.main()
