# -*- coding: utf-8 -*-
# vim: set fileencoding=utf-8 ts=4 sw=4 expandtab :
# This work © 2025-2026 by Jiří Zahradník and Martin Mokrejš
# (First Medical Faculty - Charles University in Prague) is licensed under
# Creative Commons Attribution 4.0 International. To view a copy of this
# license, visit https://creativecommons.org/licenses/by/4.0/

"""Factory functions for creating option namespaces programmatically.

These functions return pre-populated ``argparse.Namespace`` objects with
the same defaults as the CLI parsers.  They are the recommended way to
configure the library when calling functions from Python code rather
than from the command line.

Example
-------
>>> from mutation_scatter_plot.mutation_scatter_plot.options import scatter_options
>>> opts = scatter_options(
...     tsv_file_path='data.tsv',
...     outfile_prefix='my_plot',
...     aminoacids=True,
...     xmin=430, xmax=528,
... )
"""

import argparse


def scatter_options(**overrides) -> argparse.Namespace:
    """Create a scatter plot options namespace with sensible defaults.

    All CLI defaults are applied first, then any keyword arguments
    override them.  Unknown keys raise ``TypeError`` so that typos
    are caught early.

    Parameters
    ----------
    **overrides
        Any attribute of the scatter plot CLI namespace.  Common ones:

        - ``tsv_file_path`` (str): path to the input ``.frequencies.tsv``
        - ``outfile_prefix`` (str): output path prefix (PNG, PDF, HTML, …)
        - ``aminoacids`` (bool): amino-acid mode if ``True``, codon if ``False``
        - ``xmin``, ``xmax`` (int): X-axis range
        - ``colormap`` (str): matplotlib colormap name
        - ``matrix`` (str): BLOSUM matrix name
        - ``linear_circle_size`` (bool): linear scaling mode
        - ``threshold`` (float): minimum frequency to display
        - ``disable_showing_bokeh`` (bool): suppress browser popup
        - ``disable_showing_mplcursors`` (bool): suppress interactive window
        - ``dpi`` (int): output image resolution

    Returns
    -------
    argparse.Namespace
        A fully populated options namespace.

    Example
    -------
    >>> opts = scatter_options(
    ...     tsv_file_path='spike.frequencies.tsv',
    ...     outfile_prefix='output/spike',
    ...     aminoacids=True,
    ...     xmin=430, xmax=528,
    ...     disable_showing_bokeh=True,
    ...     disable_showing_mplcursors=True,
    ... )
    """
    from .cli import build_option_parser
    parser = build_option_parser()
    defaults = parser.parse_args([
        '--tsv', overrides.pop('tsv_file_path', ''),
        '--outfile-prefix', overrides.pop('outfile_prefix', ''),
    ])
    _apply_overrides(defaults, overrides)
    return defaults


def timeline_options(**overrides) -> argparse.Namespace:
    """Create a timeline plot options namespace with sensible defaults.

    Parameters
    ----------
    **overrides
        Any attribute of the timeline CLI namespace.  Common ones:

        - ``directory`` (str): path to directory with per-month TSV files
        - ``positions`` (list[str]): position specs like ``['N501Y', '484']``
        - ``outfile_prefix`` (str): output path prefix
        - ``threshold`` (float): minimum frequency to display
        - ``colormap`` (str): matplotlib colormap name (convenience alias;
          stored as ``colormaps`` list internally)
        - ``disable_showing_bokeh`` (bool): suppress browser popup

    Returns
    -------
    argparse.Namespace
    """
    from .timeline_cli import build_option_parser
    parser = build_option_parser()
    defaults = parser.parse_args([
        '--dir', overrides.pop('directory', ''),
        '--positions', *(overrides.pop('positions', ['501'])),
    ])
    # Accept colormap= (str) as convenience alias for colormaps= (list)
    if 'colormap' in overrides:
        _cmap = overrides.pop('colormap')
        defaults.colormaps = [_cmap] if isinstance(_cmap, str) else list(_cmap)
        defaults.colormap = defaults.colormaps[0]
    _apply_overrides(defaults, overrides)
    return defaults


def frequency_options(**overrides) -> argparse.Namespace:
    """Create a codon frequency calculator options namespace.

    Parameters
    ----------
    **overrides
        Any attribute of the frequency CLI namespace.  Common ones:

        - ``alignment_infilename`` (str): path to the alignment FASTA
        - ``reference_infilename`` (str): path to the reference FASTA
        - ``outfileprefix`` (str): output path prefix
        - ``padded_reference`` (bool): reference is padded
        - ``x_after_count`` (bool): IDs contain ``NNNNx.`` count prefixes
        - ``threads`` (int): number of worker processes (0 = auto)
        - ``translation_table`` (int): genetic code table number

    Returns
    -------
    argparse.Namespace
    """
    from ..calculate_codon_frequencies.cli import build_option_parser
    parser = build_option_parser()
    defaults = parser.parse_args([
        '--alignment-file', overrides.pop('alignment_infilename', ''),
        '--reference-infile', overrides.pop('reference_infilename', ''),
        '--outfile-prefix', overrides.pop('outfileprefix', ''),
    ])
    _apply_overrides(defaults, overrides)
    return defaults


def _apply_overrides(namespace: argparse.Namespace, overrides: dict) -> None:
    """Apply keyword overrides to a namespace, raising on unknown keys."""
    for key, value in overrides.items():
        if not hasattr(namespace, key):
            valid = sorted(vars(namespace).keys())
            raise TypeError(
                f"Unknown option {key!r}. Valid options are: {valid}"
            )
        setattr(namespace, key, value)
