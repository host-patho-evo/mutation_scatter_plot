"""Shared colorbar helpers for matplotlib and Bokeh renderers.

Both matplotlib and Bokeh have design limitations that prevent correct colorbar
rendering out of the box for our use case (discrete integer score bands with
semi-transparent colours).  The workarounds below were discovered empirically
in ``mutation_scatter_plot`` and extracted here so that
``mutation_timeline_plot`` produces identical colorbars without duplicating
the tricks.

Matplotlib workarounds
----------------------
1. **BoundaryNorm tick centering.**  When ``BoundaryNorm`` is used to map
   integer scores to discrete colour bands, matplotlib places tick marks at the
   *boundary* values (e.g. -3, -2, -1, 0, 1, 2, 3, 4).  But the user expects
   the integer *labels* to sit at the **visual centre** of each band, not at
   a boundary edge.  We solve this by constructing the boundaries at every
   integer from ``vmin`` to ``vmax + 1`` (N + 1 boundaries → N bands, each 1
   unit wide), then explicitly replacing the auto-ticks with positions shifted
   by **+0.5**::

       ticks at  vmin + 0.5,  vmin + 1.5,  …,  vmax + 0.5
       labels    vmin,        vmin + 1,     …,  vmax

   This places each label at the midpoint of its band.

2. **ListedColormap slicing.**  The full palette returned by ``get_colormap``
   may span a wider range than the displayed score range (e.g. -19…+19 for
   amino_acid_changes while the data only covers -3…+3).  We slice the palette
   to ``[vmin, vmax]`` using the same ``colors[norm(score)]`` lookup that the
   scatter circles use, ensuring colour consistency between glyphs and
   colorbar.

3. **Dedicated cax axes.**  ``tight_layout()`` compresses auxiliary axes and
   makes colorbar bands too small.  The caller must use ``subplots_adjust()``
   with fixed margins instead, and pass a dedicated ``cax`` axes created via
   ``gridspec`` ``width_ratios``.

Bokeh workarounds
-----------------
4. **Alpha pre-blending.**  Bokeh's ``ColorBar`` does not support an ``alpha``
   property on its colour bands.  To match the semi-transparent appearance of
   scatter-plot glyphs (drawn with e.g. ``alpha=0.5``), we analytically
   pre-blend each palette colour with white::

       apparent_channel = alpha × source + (1 − alpha) × 255

   The resulting opaque hex colours visually match the translucent circles when
   rendered over a white background.

5. **LinearColorMapper ±0.5 range offset.**  Bokeh's ``LinearColorMapper``
   maps a continuous range ``[low, high]`` to the palette by dividing it into
   ``len(palette)`` equal bins.  To make an integer tick (e.g. ``2``) land at
   the *centre* of the corresponding colour band, we extend ``low`` and
   ``high`` by 0.5 beyond the integer range::

       low  = vmin − 0.5   (e.g. -3.5)
       high = vmax + 0.5   (e.g.  3.5)

   Each band then spans exactly one integer unit (e.g. [1.5, 2.5)), and the
   integer tick at 2.0 sits at its geometric midpoint.

6. **FixedTicker.**  Bokeh's automatic tickers produce non-integer or
   misaligned ticks.  We use ``FixedTicker`` with the explicit integer list
   ``[vmin, …, vmax]`` to guarantee one tick per score band.
"""

import typing

import matplotlib
import matplotlib.cm
import matplotlib.colors
import numpy as np


def blend_with_white(hex_color: str, alpha: float) -> str:
    """Pre-blend *hex_color* with white at *alpha* opacity (Bokeh workaround).

    Bokeh's ``ColorBar`` widget renders palette bands as fully opaque
    rectangles — there is no way to set ``alpha`` on individual bands.  To
    reproduce the washed-out appearance of semi-transparent scatter glyphs
    (which are drawn over a white figure background), we compute the colour
    that *would* result from alpha-compositing the source colour onto white::

        apparent = alpha × source + (1 − alpha) × 255

    This is the standard Porter-Duff "source over destination" formula with a
    white destination (RGB = 255, 255, 255).

    Parameters
    ----------
    hex_color : str
        Six-digit CSS hex colour string, e.g. ``'#ffa200'``.  Must include the
        leading ``#`` and contain exactly 6 hex digits (no alpha digit).
    alpha : float
        Opacity in ``[0, 1]``.  ``0`` → fully white; ``1`` → original colour.

    Returns
    -------
    str
        Pre-blended hex colour string, e.g. ``'#ffd17f'``.
    """
    r = int(hex_color[1:3], 16)
    g = int(hex_color[3:5], 16)
    b = int(hex_color[5:7], 16)
    ra = round(alpha * r + (1 - alpha) * 255)
    ga = round(alpha * g + (1 - alpha) * 255)
    ba = round(alpha * b + (1 - alpha) * 255)
    return f"#{ra:02x}{ga:02x}{ba:02x}"


def setup_matplotlib_colorbar(
    fig: typing.Any,
    cax: typing.Any,
    norm: typing.Any,
    cmap: typing.Any,
    colors: typing.Any,
    vmin: int,
    vmax: int,
    label: str = 'BLOSUM score',
    alpha: float = 0.5,
) -> None:
    """Create a matplotlib colorbar in the dedicated *cax* axes.

    **Always builds a discrete ``ListedColormap`` + ``BoundaryNorm``**, even
    for continuous colormaps like ``coolwarm_r``.  Since BLOSUM/PAM scores are
    always integers, discrete colour bands are more informative than a smooth
    gradient — the reader can instantly see the distinct colour assigned to
    each integer score value.

    **Why a dedicated cax?**
    Using ``fig.colorbar(mappable, ax=some_axes)`` steals space from the data
    axes and makes layout unpredictable.  Instead, a separate ``cax`` is
    created via ``gridspec`` / ``width_ratios``.  The caller must also use
    ``fig.subplots_adjust(…)`` instead of ``plt.tight_layout()``, because
    ``tight_layout`` compresses auxiliary axes and makes the colour bands too
    narrow to visually verify tick centering (this was the root cause of the
    "ticks not centred" issue in the timeline renderer).

    Tick centering trick
    ~~~~~~~~~~~~~~~~~~~~
    ``BoundaryNorm`` with boundaries ``[vmin, vmin+1, …, vmax, vmax+1]``
    creates N = (vmax − vmin + 1) colour bands, each one unit wide.  The band
    for score ``s`` spans ``[s, s+1)`` in data coordinates.  matplotlib's
    default tick placement would put labels at the *boundary* values, so label
    ``s`` would sit at the left/bottom *edge* of its band.

    To centre each label, we replace the auto-ticks with positions at
    ``s + 0.5`` (the midpoint of each 1-unit band)::

        tick positions:  vmin + 0.5,  vmin + 1.5,  …,  vmax + 0.5
        tick labels:     vmin,        vmin + 1,     …,  vmax

    For pre-built palettes (``amino_acid_changes``, ``dkeenan``), the sliced
    palette is built by looking up ``colors[norm(s)]`` for each integer score,
    which is the same lookup used to colour scatter glyphs — guaranteeing
    colour consistency between the data points and the colorbar.

    For continuous colormaps (``coolwarm_r``, etc.), the cmap is sampled at
    evenly-spaced fractions to produce one colour per integer score, then
    wrapped in a ``ListedColormap``.  The same ``BoundaryNorm`` + tick
    centering trick is applied.

    Parameters
    ----------
    fig : matplotlib.figure.Figure
        The parent figure (needed for ``fig.colorbar(…)``).
    cax : matplotlib.axes.Axes
        Dedicated axes for the colorbar, created via ``gridspec`` or
        ``width_ratios``.  Must **not** be a data-plot axes.
    norm : matplotlib.colors.BoundaryNorm or None
        ``BoundaryNorm`` from ``get_colormap`` for pre-built discrete
        palettes, or ``None`` for continuous colormaps (which will be
        discretised here).
    cmap : matplotlib.colors.Colormap
        Colormap instance (e.g. ``coolwarm_r``, or a ``ListedColormap``).
    colors : list or None
        Resolved palette list from ``get_colormap``.  Required for the
        pre-built discrete palette path (each element is an RGBA tuple or hex
        string that matches the scatter-circle colouring).  May be ``None``
        for continuous cmaps, which will be sampled automatically.
    vmin, vmax : int
        Integer score range to display on the colorbar (inclusive on both
        ends).  Typically derived from ``cmap_vmin`` / ``cmap_vmax`` on
        ``myoptions``, which is clamped to the data-driven score range.
    label : str
        Text label for the colorbar axis.
    alpha : float
        Alpha transparency applied to the colorbar bands via
        ``fig.colorbar(…, alpha=…)``.
    """
    n_bands = vmax - vmin + 1

    if norm is not None and colors is not None:
        # ── Pre-built discrete palette (amino_acid_changes, dkeenan) ──
        #
        # Slice the full palette to only [vmin, vmax].
        # Use the same colors[norm(s)] lookup as the scatter circles so that
        # the colorbar colours are pixel-identical to the glyph colours.
        # Clamp the index to valid bounds to handle edge cases where norm(s)
        # falls outside the palette range.
        _cb_sliced = [
            colors[max(0, min(len(colors) - 1, norm(s)))]
            for s in range(vmin, vmax + 1)
        ]
    elif cmap is not None:
        # ── Continuous cmap (coolwarm_r etc.) — discretise into bands ──
        #
        # Sample the continuous cmap at evenly-spaced fractions to produce
        # one flat-colour band per integer score.  This makes it easy for
        # the reader to distinguish individual score values instead of
        # seeing a smooth gradient where neighbouring integers blend.
        _cb_sliced = [
            cmap(i / max(1, n_bands - 1))
            for i in range(n_bands)
        ]
    else:
        return  # nothing to render

    # Build a fresh ListedColormap + BoundaryNorm for the (possibly sliced)
    # palette.  The BoundaryNorm boundaries are [vmin, vmin+1, …, vmax+1],
    # i.e. N+1 boundaries for N colour bands, where N = vmax − vmin + 1.
    _cb_cmap = matplotlib.colors.ListedColormap(_cb_sliced, "sliced")
    _cb_norm = matplotlib.colors.BoundaryNorm(
        np.arange(vmin, vmax + 2, 1), len(_cb_sliced),
    )
    _cb_sm = matplotlib.cm.ScalarMappable(cmap=_cb_cmap, norm=_cb_norm)
    _cb_sm.set_array([])  # required by matplotlib even though we don't map data

    # Render the colorbar.
    # location='right' is ignored when cax is given, but kept for clarity.
    # pad=-0.1 reduces whitespace between the plot and colorbar.
    # NOTE: We do NOT pass label= to fig.colorbar() here.  The default
    # placement of the label is a rotated text along the colorbar axis,
    # which gets obscured by the large position labels on the timeline's
    # secondary Y-axis (ax2.twinx).  Instead, we set the label as a title
    # above the colorbar axes, where it is always visible.
    _colorbar = fig.colorbar(
        _cb_sm, cax=cax, location='right', pad=-0.1,
        alpha=alpha,
    )
    # Place the label as a title above the colorbar so it is never obscured
    # by neighbouring elements (e.g. the large position labels on the
    # timeline's secondary Y-axis).
    _colorbar.ax.set_title(label, fontsize=9, pad=4)

    # TICK CENTERING TRICK.
    # Each band for score s spans [s, s+1) in data coordinates.
    # The midpoint is at s + 0.5.  We set custom tick positions at these
    # midpoints and label them with the integer scores.
    # Without this, matplotlib would put tick labels at boundary positions
    # (i.e. at band *edges*), making -3 appear at the bottom edge and +3
    # at the top edge of their bands rather than centred.
    _colorbar.ax.set_yticks(
        np.arange(vmin + 0.5, vmax + 1.5, 1),
        np.arange(vmin, vmax + 1, 1),
    )
    # Suppress minor ticks that matplotlib may auto-add.
    _colorbar.ax.tick_params(axis='y', which='minor', length=0)


def build_bokeh_colorbar_palette(
    norm: typing.Any,
    cmap: typing.Any,
    colors: typing.Any,
    vmin: int,
    vmax: int,
    alpha: float = 0.5,
) -> list[str]:
    """Build a Bokeh-ready pre-blended hex palette for the colorbar.

    Returns a list of hex colour strings, one per integer score value in
    ``[vmin, vmax]``, pre-blended with white to simulate the given *alpha*
    on a white background.

    **Why pre-blend?**  Bokeh's ``ColorBar`` widget renders each palette entry
    as a fully opaque rectangle.  Unlike matplotlib's ``fig.colorbar(…,
    alpha=…)``, there is no ``alpha`` parameter on Bokeh's ``ColorBar`` or
    ``LinearColorMapper``.  To make the colorbar visually match the
    semi-transparent scatter circles (which are drawn with e.g.
    ``circle(…, alpha=0.5)`` on a white ``Plot.background_fill_color``), we
    analytically pre-compute the composited colour using
    :func:`blend_with_white`.

    The palette is built using the same ``colors[norm(s)]`` lookup (discrete
    path) or ``cmap(fraction)`` sampling (continuous path) that the scatter
    glyphs use, so the colorbar colours are consistent with the data points.

    Parameters
    ----------
    norm : BoundaryNorm or None
        Discrete normaliser from ``get_colormap``, or ``None`` for continuous.
    cmap : Colormap
        Colormap instance.
    colors : list or None
        Resolved palette list for the discrete path (from ``get_colormap``).
    vmin, vmax : int
        Integer score range (inclusive).
    alpha : float
        Circle opacity to simulate in the colorbar bands.

    Returns
    -------
    list[str]
        One hex colour string per score band, suitable as ``palette`` arg
        for ``bokeh.models.LinearColorMapper``.
    """
    n_bands = vmax - vmin + 1
    palette: list[str] = []
    if norm is not None and colors is not None:
        # Pre-built discrete palette: lookup and pre-blend each score's colour.
        for s in range(vmin, vmax + 1):
            idx = max(0, min(len(colors) - 1, norm(s)))
            hex_c = matplotlib.colors.to_hex(
                matplotlib.colors.to_rgba(colors[idx]),
            )
            palette.append(blend_with_white(hex_c, alpha))
    elif cmap is not None:
        # Continuous cmap — discretise into per-integer bands, same as the
        # matplotlib helper.  Sample at evenly-spaced fractions.
        for i in range(n_bands):
            hex_c = matplotlib.colors.to_hex(cmap(i / max(1, n_bands - 1)))
            palette.append(blend_with_white(hex_c, alpha))
    else:
        # Fallback: neutral grey bands.
        palette = ['#aaaaaa'] * n_bands
    return palette


def add_bokeh_colorbar(
    bokeh_fig: typing.Any,
    norm: typing.Any,
    cmap: typing.Any,
    colors: typing.Any,
    vmin: int,
    vmax: int,
    alpha: float = 0.5,
    label: str = 'BLOSUM score',
) -> None:
    """Add a ``ColorBar`` to a Bokeh figure with centred ticks and alpha pre-blend.

    This is the Bokeh counterpart of :func:`setup_matplotlib_colorbar`.  It
    applies two workarounds for Bokeh's ``ColorBar`` design limitations:

    1. **Alpha pre-blending** — Bokeh ``ColorBar`` has no ``alpha`` property.
       Palette colours are pre-blended with white via
       :func:`build_bokeh_colorbar_palette` so the colorbar visually matches
       the semi-transparent scatter-circle glyphs.

    2. **±0.5 range offset for tick centering** — ``LinearColorMapper`` divides
       the continuous range ``[low, high]`` into ``len(palette)`` equal-width
       bins.  If we used ``low=vmin, high=vmax+1`` (the natural N-boundary
       range), the bins would be ``[vmin, vmin+1), [vmin+1, vmin+2), …`` and
       an integer tick at ``vmin`` would sit at the *left/bottom edge* of the
       first bin rather than its centre.

       By setting ``low = vmin − 0.5`` and ``high = vmax + 0.5``, each bin
       spans ``[s − 0.5, s + 0.5)`` and the integer tick at ``s`` lands
       exactly at the bin's geometric midpoint::

           bin for score 0:  [-0.5, +0.5)  →  tick at 0.0 = centre  ✓
           bin for score 3:  [ 2.5,  3.5)  →  tick at 3.0 = centre  ✓

    3. **FixedTicker** — Bokeh's automatic tick algorithms produce non-integer
       or unevenly-spaced ticks for small ranges (e.g. [-3, 3]).
       ``FixedTicker`` with the explicit integer list guarantees one tick per
       score band.

    Parameters
    ----------
    bokeh_fig : bokeh.plotting.Figure
        Target Bokeh figure to which the ColorBar will be added (``'right'``
        layout position).
    norm : BoundaryNorm or None
        Discrete normaliser from ``get_colormap``, or ``None`` for continuous
        cmaps.
    cmap : Colormap
        Colormap instance.
    colors : list or None
        Resolved palette list for the discrete path.
    vmin, vmax : int
        Integer score range (inclusive).
    alpha : float
        Circle opacity to simulate in the colorbar.
    label : str
        Colorbar title text (displayed above the colorbar).
    """
    import bokeh.models  # pylint: disable=import-outside-toplevel

    palette = build_bokeh_colorbar_palette(norm, cmap, colors, vmin, vmax, alpha)

    # ±0.5 offset trick: extend the mapper range by half a unit on each side
    # so that the integer ticks (placed by FixedTicker below) fall at the
    # geometric centre of each colour band, not at a band edge.
    mapper = bokeh.models.LinearColorMapper(
        palette=palette,
        low=vmin - 0.5,
        high=vmax + 0.5,
    )
    colorbar = bokeh.models.ColorBar(
        color_mapper=mapper,
        label_standoff=8,
        title=label,
        title_standoff=10,
        location=(0, 0),
        # Explicit integer ticks — Bokeh's auto-ticker may choose non-integer
        # or unevenly-spaced values for small ranges like [-3, 3].
        ticker=bokeh.models.FixedTicker(ticks=list(range(vmin, vmax + 1))),
    )
    bokeh_fig.add_layout(colorbar, 'right')
