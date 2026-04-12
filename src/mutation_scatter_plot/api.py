# -*- coding: utf-8 -*-
# vim: set fileencoding=utf-8 ts=4 sw=4 expandtab :
# This work © 2025-2026 by Jiří Zahradník and Martin Mokrejš
# (First Medical Faculty - Charles University in Prague) is licensed under
# Creative Commons Attribution 4.0 International. To view a copy of this
# license, visit https://creativecommons.org/licenses/by/4.0/

"""High-level programmatic API for mutation_scatter_plot.

These functions wrap the low-level rendering pipeline into single-call
entry points suitable for library consumers, Jupyter notebooks, and
3rd-party applications.

Quick start
-----------
>>> from mutation_scatter_plot.api import render_scatter
>>> result = render_scatter(
...     tsv_path='data.frequencies.tsv',
...     outfile_prefix='output/my_plot',
...     aminoacids=True,
...     xmin=430, xmax=528,
... )
>>> result['figure'].savefig('custom.png')   # matplotlib Figure
>>> result['files_written']                  # list of saved paths
"""

import os
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    import pandas as pd



def render_scatter(
    tsv_path: str,
    outfile_prefix: str = '',
    *,
    aminoacids: bool = False,
    xmin: int = 0,
    xmax: int = 0,
    colormap: str = 'amino_acid_changes',
    matrix: str = 'BLOSUM80',
    linear_scaling: bool = False,
    threshold: float = 0.001,
    include_synonymous: bool = False,
    show_bokeh: bool = False,
    show_mplcursors: bool = False,
    dpi: int = 600,
    title: str = '',
    **extra_options,
) -> dict:
    """Render a mutation scatter plot from a frequency TSV file.

    If *outfile_prefix* is provided, PNG/PDF/HTML files are written and
    their paths are returned.  If omitted, only the matplotlib ``Figure``
    is created (no files are written).

    Parameters
    ----------
    tsv_path : str
        Path to a ``.frequencies.tsv`` file produced by
        ``calculate_codon_frequencies``.
    outfile_prefix : str, optional
        Output path prefix.  When empty, no files are saved.
    aminoacids : bool
        Amino-acid mode (``True``) or codon mode (``False``).
    xmin, xmax : int
        X-axis range (padded alignment positions).
    colormap : str
        Matplotlib colormap name (default ``'amino_acid_changes'``).
    matrix : str
        BLOSUM matrix name (default ``'BLOSUM80'``).
    linear_scaling : bool
        Use linear circle-size scaling instead of area scaling.
    threshold : float
        Minimum frequency to display a data point (default 0.001).
    include_synonymous : bool
        Include synonymous mutations in amino-acid mode.
    show_bokeh : bool
        Open the Bokeh HTML in a browser after rendering.
    show_mplcursors : bool
        Open the interactive matplotlib window with hover support.
    dpi : int
        Resolution for raster outputs (default 600).
    title : str
        Figure title; auto-derived from filename if empty.
    **extra_options
        Any additional options forwarded to the scatter options namespace.

    Returns
    -------
    dict
        - ``'figure'``: ``matplotlib.figure.Figure``
        - ``'axes'``: tuple ``(ax1, ax2, ax3, ax4)``
        - ``'dataframe'``: ``pandas.DataFrame`` — the cleaned input data
        - ``'files_written'``: ``list[str]`` — paths of saved files (empty
          if *outfile_prefix* was not given)

    Example
    -------
    >>> result = render_scatter(
    ...     tsv_path='spike.frequencies.tsv',
    ...     aminoacids=True, xmin=430, xmax=528,
    ... )
    >>> fig = result['figure']
    >>> fig.savefig('my_custom_output.png', dpi=300)
    """
    # pylint: disable=too-many-locals
    from .mutation_scatter_plot.options import scatter_options
    from .mutation_scatter_plot import (
        load_and_clean_dataframe,
        build_frequency_tables,
        setup_matplotlib_figure,
        collect_scatter_data,
        render_bokeh,
        render_matplotlib,
        load_matrix,
    )

    # When no outfile_prefix is given, we still need one for load_matrix()
    # which requires it.  Use a tempdir that we'll clean up at the end.
    import tempfile as _tempfile
    _tmpdir = None
    _save_files = bool(outfile_prefix)
    if not outfile_prefix:
        _tmpdir = _tempfile.mkdtemp()
        outfile_prefix = os.path.join(_tmpdir, 'api_tmp')

    opts = scatter_options(
        tsv_file_path=tsv_path,
        outfile_prefix=outfile_prefix,
        aminoacids=aminoacids,
        xmin=xmin, xmax=xmax,
        colormap=colormap,
        matrix=matrix,
        linear_circle_size=linear_scaling,
        threshold=threshold,
        include_synonymous=include_synonymous,
        disable_showing_bokeh=not show_bokeh,
        disable_showing_mplcursors=not show_mplcursors,
        dpi=dpi,
        title=title,
        **extra_options,
    )

    # --linear-circle-size implies --disable-bokeh-sqrt-size (same as CLI)
    if opts.linear_circle_size:
        opts.bokeh_sqrt_size = False

    # ── 1. Load matrix ───────────────────────────────────────────────────
    (_matrix, _matrix_name, _min_theoretical_score,
     _max_theoretical_score, _outfile_prefix) = load_matrix(opts)

    # ── 2. Load and clean the TSV ────────────────────────────────────────
    _padded_position2position = {}
    _df, _padded_position2position = load_and_clean_dataframe(
        opts, tsv_path, _padded_position2position,
    )

    # ── 3. Auto-detection logic (mirrors cli.py:main()) ──────────────────
    _effective_xmin = opts.xmin if opts.xmin else min(_padded_position2position.keys())
    _effective_xmax = opts.xmax if opts.xmax else max(_padded_position2position.keys())

    if not opts.disable_padded_x_axis:
        _window = {k: v for k, v in _padded_position2position.items()
                   if _effective_xmin <= k <= _effective_xmax}
        if _window:
            _differences = set(k - v for k, v in _window.items())
            if len(_differences) <= 1:
                opts.disable_padded_x_axis = True

    if opts.showins is None:
        _ins_in_range = _df[
            (_df['original_aa'] == 'INS') &
            (_df['padded_position'] >= _effective_xmin) &
            (_df['padded_position'] <= _effective_xmax)
        ]
        opts.showins = len(_ins_in_range) > 0

    if opts.showdel is None:
        _del_in_range = _df[
            (_df['mutant_aa'] == 'DEL') &
            (_df['padded_position'] >= _effective_xmin) &
            (_df['padded_position'] <= _effective_xmax)
        ]
        opts.showdel = len(_del_in_range) > 0

    # ── 4. Supplement the conversion dictionary ──────────────────────────
    _unchanged_tsv = tsv_path.replace(
        '.frequencies.tsv', '.frequencies.unchanged_codons.tsv'
    )
    if os.path.exists(_unchanged_tsv):
        _df_unch, _padded_position2position = load_and_clean_dataframe(
            opts, _unchanged_tsv, _padded_position2position,
        )
        del _df_unch

    # ── 5. Build frequency tables ────────────────────────────────────────
    (
        _amino_acids, _codons_whitelist, _codons_whitelist2,
        _final_sorted_whitelist,
        _unique_aa_padded_positions, _unique_codon_padded_positions,
        _old_aa_table, _new_aa_table, _old_codon_table, _new_codon_table,
        _calculated_aa_offset, _padded_position2position,
    ) = build_frequency_tables(opts, _df, _padded_position2position)

    # ── 6. Title ─────────────────────────────────────────────────────────
    if not title:
        _title_data = os.path.basename(tsv_path).replace('.frequencies.tsv', '')
    else:
        _title_data = title

    # Determine aln_rows from .count file
    _aln_rows = '0'
    if '.frequencies.tsv' in tsv_path:
        _count_filename = f"{_outfile_prefix}.count"
        if os.path.exists(_count_filename):
            try:
                with open(_count_filename, encoding="utf-8") as _fh:
                    _aln_rows = _fh.readline().strip()
            except OSError:
                pass

    # ── 7. Set up matplotlib figure ──────────────────────────────────────
    _figure, _ax1, _ax2, _ax3, _ax4, _xmin, _xmax = \
        setup_matplotlib_figure(
            opts,
            _title_data, _aln_rows, _matrix_name, _amino_acids,
            _codons_whitelist, _final_sorted_whitelist,
            _unique_aa_padded_positions, _unique_codon_padded_positions,
            _new_aa_table, _new_codon_table, _padded_position2position,
        )

    # ── 8. Collect scatter data ──────────────────────────────────────────
    _table = _new_aa_table if aminoacids else _new_codon_table
    (
        _norm, _cmap, _colors, _used_colors, _matrix_values,
        _mutations,
        _circles_bokeh, _circles_matplotlib, _markers, _dots,
        _hover_text_bokeh,
    ) = collect_scatter_data(
        opts,
        _df, _table, _outfile_prefix, _matrix, _amino_acids,
        _codons_whitelist2, _padded_position2position,
        _xmin, _xmax,
    )

    files_written = []

    # ── 9. Render ────────────────────────────────────────────────────────
    _xlabel = _ax1.get_xlabel()

    if _circles_bokeh and _save_files:
        render_bokeh(
            opts,
            _outfile_prefix, _xmin, _xmax, _amino_acids,
            _final_sorted_whitelist,
            _circles_bokeh, _mutations, _hover_text_bokeh,
            _title_data, _xlabel,
            _matrix_name, _colors, _norm, _cmap,
            _padded_position2position,
            show=show_bokeh,
        )

    if _save_files:
        render_matplotlib(
            opts,
            _figure, _ax1, _ax2, _ax3, _ax4, _outfile_prefix,
            _circles_matplotlib, _markers, _dots, _cmap, _norm, _colors,
            _matrix, _matrix_name,
            show=show_mplcursors,
        )

        # Discover all files that were written
        out_dir = os.path.dirname(_outfile_prefix) or '.'
        prefix_base = os.path.basename(_outfile_prefix)
        for fname in os.listdir(out_dir):
            if fname.startswith(prefix_base):
                files_written.append(os.path.join(out_dir, fname))

    # Clean up temp directory if we created one
    if _tmpdir is not None:
        import shutil
        shutil.rmtree(_tmpdir, ignore_errors=True)

    return {
        'figure': _figure,
        'axes': (_ax1, _ax2, _ax3, _ax4),
        'dataframe': _df,
        'files_written': sorted(files_written),
    }


def render_timeline(
    directory: str,
    positions: list[str],
    outfile_prefix: str = '',
    *,
    aminoacids: bool = False,
    threshold: float = 0.0,
    colormap: str = 'coolwarm_r',
    matrix: str = 'BLOSUM80',
    show_bokeh: bool = False,
    dpi: int = 600,
    **extra_options,
) -> dict:
    """Render a timeline scatter plot from per-month frequency TSV files.

    Scans *directory* for ``.frequencies.tsv`` files named with a
    ``YYYY-MM`` month tag and renders a timeline where the X-axis
    represents time and circle size/colour encode mutation frequency
    and BLOSUM score.

    Parameters
    ----------
    directory : str
        Directory containing per-month ``.frequencies.tsv`` files.
    positions : list[str]
        Position specifications, e.g. ``['N501Y', '498[RHQ]', '484']``.
    outfile_prefix : str, optional
        Output path prefix.  When empty, auto-inferred from filenames.
    aminoacids : bool
        Amino-acid mode (``True``) or codon mode (``False``).
    threshold : float
        Minimum frequency to include in the timeline (default 0.0).
    colormap : str
        Matplotlib colormap name (default ``'coolwarm_r'``).
    matrix : str
        BLOSUM matrix name (default ``'BLOSUM80'``).
    show_bokeh : bool
        Open the Bokeh HTML in a browser.
    dpi : int
        Resolution for raster outputs (default 600).
    **extra_options
        Additional options forwarded to the timeline options namespace.

    Returns
    -------
    dict
        - ``'data'``: ``TimelineData`` object with all collected points
        - ``'files_written'``: ``list[str]``
    """
    # pylint: disable=too-many-locals
    from .mutation_scatter_plot.options import timeline_options
    from .mutation_scatter_plot import load_matrix
    from .mutation_scatter_plot.core import get_colormap
    from .mutation_scatter_plot.timeline import (
        scan_directory,
        infer_common_prefix,
        parse_positions,
        collect_timeline_data,
        render_timeline_matplotlib,
        render_timeline_bokeh,
    )

    opts = timeline_options(
        directory=directory,
        positions=positions,
        outfile_prefix=outfile_prefix,
        aminoacids=aminoacids,
        threshold=threshold,
        colormap=colormap,
        matrix=matrix,
        dpi=dpi,
        disable_showing_bokeh=not show_bokeh,
        **extra_options,
    )

    # ── 1. Discover files ────────────────────────────────────────────────
    files = scan_directory(
        opts.input_dir, opts.pattern, opts.include_unknown_month,
    )
    if not files:
        raise FileNotFoundError(
            f"No .frequencies.tsv files found in {directory!r}"
        )

    # ── 2. Parse positions ───────────────────────────────────────────────
    specs = parse_positions(positions)

    # ── 3. Auto-infer output prefix ──────────────────────────────────────
    if not opts.outfile_prefix:
        opts.outfile_prefix = infer_common_prefix(files, opts.input_dir)

    # Append .timeline.aa/.codon before load_matrix adds .BLOSUM80.area_scaling.coolwarm_r
    _mode_label = 'aa' if aminoacids else 'codon'
    opts.outfile_prefix = opts.outfile_prefix + '.timeline.' + _mode_label

    # ── 4. Load matrix ───────────────────────────────────────────────────
    (_matrix, _matrix_name, _min_score,
     _max_score, _outfile_prefix) = load_matrix(opts)

    # ── 5. Initial colormap (will be refined after data collection) ─────
    opts.cmap_actual_vmin = -10
    opts.cmap_actual_vmax = 10
    _norm, _cmap, _colors = get_colormap(opts, colormap)

    # ── 6. Collect timeline data ─────────────────────────────────────────
    data = collect_timeline_data(files, specs, opts, _matrix, _norm, _colors)

    if not data.points:
        return {
            'data': data,
            'files_written': [],
        }

    # ── 7. Re-derive colormap with actual score range ────────────────────
    _actual_scores = [pt.score for pt in data.points if pt.score != 12]
    _actual_vmin = min(_actual_scores) if _actual_scores else -11
    _actual_vmax = max(_actual_scores) if _actual_scores else 11
    _actual_bound = max(abs(_actual_vmin), abs(_actual_vmax))
    opts.cmap_actual_vmin = -_actual_bound
    opts.cmap_actual_vmax = _actual_bound

    if hasattr(opts, 'cmap_vmin'):
        delattr(opts, 'cmap_vmin')
    if hasattr(opts, 'cmap_vmax'):
        delattr(opts, 'cmap_vmax')
    _norm, _cmap, _colors = get_colormap(opts, colormap)

    # ── 8. Render ────────────────────────────────────────────────────────
    files_written = []

    render_timeline_matplotlib(
        data, opts, _norm, _cmap, _colors, _outfile_prefix,
    )
    for ext in ('.png', '.pdf'):
        candidate = f"{_outfile_prefix}{ext}"
        if os.path.exists(candidate):
            files_written.append(candidate)

    if show_bokeh or outfile_prefix:
        render_timeline_bokeh(
            data, opts, _norm, _cmap, _colors, _outfile_prefix,
        )
        html_candidate = f"{_outfile_prefix}.html"
        if os.path.exists(html_candidate):
            files_written.append(html_candidate)

    return {
        'data': data,
        'files_written': sorted(files_written),
    }


def calculate_frequencies(
    alignment_file: str,
    reference_file: str,
    outfile_prefix: str,
    *,
    padded_reference: bool = True,
    x_after_count: bool = False,
    print_unchanged_sites: bool = True,
    threads: int = 0,
    translation_table: int = 1,
    overwrite: bool = True,
    **extra_options,
) -> 'pd.DataFrame':
    """Calculate codon frequencies from a FASTA alignment.

    Runs the full ``calculate_codon_frequencies`` pipeline: parses the
    alignment against the reference, writes the standard TSV output
    files, and returns the result as a pandas DataFrame.

    Parameters
    ----------
    alignment_file : str
        Path to the padded multi-FASTA alignment file.
    reference_file : str
        Path to the reference FASTA file (single sequence).
    outfile_prefix : str
        Output path prefix.  Files ``{prefix}.frequencies.tsv`` and
        ``{prefix}.unchanged_codons.tsv`` are created.
    padded_reference : bool
        Whether the reference is padded with dashes (default ``True``).
    x_after_count : bool
        Whether FASTA IDs contain ``NNNNx.`` count prefixes.
    print_unchanged_sites : bool
        Also write unchanged codons to a separate TSV.
    threads : int
        Number of worker processes (0 = auto-detect).
    translation_table : int
        NCBI genetic code table number (default 1 = standard).
    overwrite : bool
        Overwrite existing output files (default ``True``).
    **extra_options
        Additional options forwarded to the frequency options namespace.

    Returns
    -------
    pandas.DataFrame
        The frequency table with columns: ``padded_position``,
        ``position``, ``original_aa``, ``mutant_aa``, ``frequency``,
        ``original_codon``, ``mutant_codon``, ``observed_codon_count``,
        ``total_codons_per_site``.

    Example
    -------
    >>> from mutation_scatter_plot.api import calculate_frequencies
    >>> df = calculate_frequencies(
    ...     alignment_file='spike_aligned.fasta',
    ...     reference_file='MN908947.3_S.3873.fasta',
    ...     outfile_prefix='output/spike',
    ... )
    >>> high_freq = df[df['frequency'] > 0.1]
    """
    from .mutation_scatter_plot.options import frequency_options
    from .calculate_codon_frequencies import (
        get_codons,
        parse_alignment,
        open_file,
    )
    from Bio import SeqIO

    opts = frequency_options(
        alignment_infilename=alignment_file,
        reference_infilename=reference_file,
        outfileprefix=outfile_prefix,
        padded_reference=padded_reference,
        x_after_count=x_after_count,
        print_unchanged_sites=print_unchanged_sites,
        threads=threads,
        translation_table=translation_table,
        overwrite=overwrite,
        **extra_options,
    )

    # Read reference — mirrors cli.py: SeqIO.parse, take first record
    for _record in SeqIO.parse(reference_file, 'fasta'):
        padded_ref_dna = str(_record.seq)
        break

    # Translate reference to protein (whole padded sequence, same as CLI)
    from . import alt_translate
    ref_protein = alt_translate(padded_ref_dna, table=translation_table)
    ref_codons = get_codons(padded_ref_dna)

    # Build output file paths — matches CLI's naming convention
    if outfile_prefix.endswith('.tsv'):
        freq_tsv_path = outfile_prefix
        unchanged_tsv_path = f"{outfile_prefix[:-4]}.unchanged_codons.tsv"
    else:
        freq_tsv_path = f"{outfile_prefix}.tsv"
        unchanged_tsv_path = f"{outfile_prefix}.unchanged_codons.tsv"

    count_path = (
        f"{'.'.join(alignment_file.split('.')[:-1])}.count"
    )

    _tsv_header = (
        "padded_position\tposition\toriginal_aa\tmutant_aa\tfrequency\t"
        "original_codon\tmutant_codon\tobserved_codon_count\t"
        "total_codons_per_site\n"
    )

    outfile = open_file(freq_tsv_path, overwrite)
    outfile.write(_tsv_header)

    outfile_unchanged = None
    if print_unchanged_sites:
        outfile_unchanged = open_file(unchanged_tsv_path, overwrite)
        outfile_unchanged.write(_tsv_header)

    outfile_count = open(count_path, 'w', encoding='utf-8')

    _aa_start = (opts.aa_start - 1) if opts.aa_start else 0
    _min_start = (opts.min_start - 1) if opts.min_start else 0
    _max_stop = opts.max_stop if opts.max_stop else 0
    _threads = threads if threads > 0 else None

    try:
        parse_alignment(
            opts, alignment_file, padded_ref_dna, ref_protein, ref_codons,
            outfile, outfile_unchanged, outfile_count,
            _aa_start, _min_start, _max_stop,
            threads=_threads,
            translation_table=translation_table,
        )
    finally:
        outfile.close()
        if outfile_unchanged:
            outfile_unchanged.close()
        outfile_count.close()

    # Read back and return as DataFrame
    import pandas as pd  # noqa: PLC0415  (deferred for import-time savings)
    df = pd.read_csv(freq_tsv_path, sep='\t')
    return df

