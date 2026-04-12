"""Shared colorbar helpers for matplotlib and Bokeh renderers.

These functions encapsulate the workarounds needed for correct colorbar
rendering (tick centering for BoundaryNorm, alpha pre-blending for Bokeh)
so that both ``mutation_scatter_plot`` and ``mutation_timeline_plot`` produce
identical, correctly centred colorbars.
"""

import typing

import matplotlib
import matplotlib.cm
import matplotlib.colors
import numpy as np


def blend_with_white(hex_color: str, alpha: float) -> str:
    """Return *hex_color* pre-blended with white at *alpha* opacity.

    Simulates how a semi-transparent glyph appears when composited over a
    white background::

        apparent_channel = round(alpha * source + (1 - alpha) * 255)

    Used for Bokeh ``ColorBar`` palette bands which cannot render with
    alpha natively.

    Parameters
    ----------
    hex_color : str
        Six-digit CSS hex colour string, e.g. ``'#ffa200'``.
    alpha : float
        Opacity in [0, 1].
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

    Handles both the discrete ``BoundaryNorm`` path (e.g. amino_acid_changes)
    and the continuous cmap path (e.g. coolwarm_r), with proper tick centering.

    This is the single source of truth for colorbar setup — extracted from
    ``render_matplotlib`` in ``mutation_scatter_plot`` to be shared by the
    timeline renderer.

    Parameters
    ----------
    fig : matplotlib.figure.Figure
        The parent figure.
    cax : matplotlib.axes.Axes
        Dedicated axes for the colorbar (from gridspec).
    norm : matplotlib.colors.BoundaryNorm or None
        Discrete normaliser, or ``None`` for continuous colormaps.
    cmap : matplotlib.colors.Colormap
        Colormap instance.
    colors : list or None
        Resolved palette list (from ``get_colormap``).  Required for the
        discrete ``BoundaryNorm`` path; may be ``None`` for continuous cmaps.
    vmin, vmax : int
        Score range to display on the colorbar.
    label : str
        Colorbar label text.
    alpha : float
        Alpha transparency applied to the colorbar bands.
    """
    if norm is not None and colors is not None:
        # Discrete BoundaryNorm path (amino_acid_changes, dkeenan).
        # Build a sliced Mappable covering only [vmin, vmax], using the
        # same colors[norm(s)] lookup as the scatter circles.
        _cb_sliced = [
            colors[max(0, min(len(colors) - 1, norm(s)))]
            for s in range(vmin, vmax + 1)
        ]
        _cb_cmap = matplotlib.colors.ListedColormap(_cb_sliced, "sliced")
        _cb_norm = matplotlib.colors.BoundaryNorm(
            np.arange(vmin, vmax + 2, 1), len(_cb_sliced),
        )
        _cb_sm = matplotlib.cm.ScalarMappable(cmap=_cb_cmap, norm=_cb_norm)
        _cb_sm.set_array([])

        _colorbar = fig.colorbar(
            _cb_sm, cax=cax, label=label, location='right', pad=-0.1,
            alpha=alpha,
        )

        # Centre each integer label inside its colour band with +0.5 offset.
        _colorbar.ax.set_yticks(
            np.arange(vmin + 0.5, vmax + 1.5, 1),
            np.arange(vmin, vmax + 1, 1),
        )
        _colorbar.ax.tick_params(axis='y', which='minor', length=0)
    elif cmap is not None:
        # Continuous cmap path (coolwarm_r etc.).
        _cb_norm = matplotlib.colors.Normalize(vmin=vmin, vmax=vmax)
        _sm = matplotlib.cm.ScalarMappable(cmap=cmap, norm=_cb_norm)
        _sm.set_array([])
        _colorbar = fig.colorbar(
            _sm, cax=cax, label=label, location='right', pad=-0.1,
            alpha=alpha,
        )
        _colorbar.ax.set_yticks(np.arange(vmin, vmax + 1, 1))
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

    Parameters
    ----------
    norm : BoundaryNorm or None
        Discrete normaliser, or ``None`` for continuous colormaps.
    cmap : Colormap
        Colormap instance.
    colors : list or None
        Resolved palette list for the discrete path.
    vmin, vmax : int
        Score range.
    alpha : float
        Circle opacity to match in the colorbar.
    """
    palette: list[str] = []
    if norm is not None and colors is not None:
        for s in range(vmin, vmax + 1):
            idx = max(0, min(len(colors) - 1, norm(s)))
            hex_c = matplotlib.colors.to_hex(
                matplotlib.colors.to_rgba(colors[idx]),
            )
            palette.append(blend_with_white(hex_c, alpha))
    elif cmap is not None:
        n = vmax - vmin + 1
        for i in range(n):
            hex_c = matplotlib.colors.to_hex(cmap(i / max(1, n - 1)))
            palette.append(blend_with_white(hex_c, alpha))
    else:
        palette = ['#aaaaaa'] * (vmax - vmin + 1)
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
    """Add a ColorBar to a Bokeh figure with centred ticks and alpha pre-blend.

    Parameters
    ----------
    bokeh_fig : bokeh.plotting.Figure
        Target Bokeh figure.
    norm, cmap, colors : same as ``build_bokeh_colorbar_palette``.
    vmin, vmax : int
        Score range.
    alpha : float
        Circle opacity to simulate.
    label : str
        Colorbar title.
    """
    import bokeh.models  # pylint: disable=import-outside-toplevel

    palette = build_bokeh_colorbar_palette(norm, cmap, colors, vmin, vmax, alpha)

    # low/high extended by ±0.5 so each band is 1 score-unit wide and
    # the integer tick coordinate falls at the geometric centre.
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
        ticker=bokeh.models.FixedTicker(ticks=list(range(vmin, vmax + 1))),
    )
    bokeh_fig.add_layout(colorbar, 'right')
