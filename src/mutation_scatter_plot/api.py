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
    from .mutation_scatter_plot.options import scatter_options
    from .mutation_scatter_plot import (
        load_and_clean_dataframe,
        build_frequency_tables,
        build_conversion_table,
        setup_matplotlib_figure,
        collect_scatter_data,
        render_bokeh,
        render_matplotlib,
    )
    from .mutation_scatter_plot.core import load_matrix, get_colormap

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

    # 1. Load data
    padded_position2position: dict = {}
    df = load_and_clean_dataframe(opts, tsv_path, padded_position2position)
    new_aa_table, new_codon_table = build_frequency_tables(opts, df, padded_position2position)
    amino_acids, codons_whitelist = build_conversion_table(df, padded_position2position)

    # 2. Set up BLOSUM matrix and colormap
    matrix_obj = load_matrix(opts)
    matrix_name = opts.matrix
    colors, norm, cmap = get_colormap(opts, opts.colormap)

    # 3. Set up matplotlib figure
    (figure, ax1, ax2, ax3, ax4,
     final_sorted_whitelist,
     unique_aa_padded_positions,
     unique_padded_codon_positions,
     title_data, aln_rows, xlabel) = setup_matplotlib_figure(
        opts, title or os.path.basename(tsv_path).replace('.frequencies.tsv', ''),
        0,  # aln_rows placeholder
        matrix_name,
        amino_acids, codons_whitelist, [],
        [], [], new_aa_table, new_codon_table,
        padded_position2position,
    )

    # 4. Collect scatter data
    (circles_bokeh, circles_matplotlib, mutations,
     hover_texts, markers, dots) = collect_scatter_data(
        opts, df, new_aa_table if aminoacids else new_codon_table,
        outfile_prefix, matrix_obj, amino_acids,
        codons_whitelist, padded_position2position, xmin, xmax,
    )

    files_written = []

    # 5. Render matplotlib
    if outfile_prefix:
        render_matplotlib(
            opts, figure, ax1, ax2, ax3, ax4,
            outfile_prefix, circles_matplotlib, markers, dots,
            cmap, norm, colors, matrix_obj, matrix_name,
            show_mplcursors,
        )
        # Collect written files
        for ext in ('.png', '.pdf'):
            candidate = f"{outfile_prefix}.{matrix_name}"
            scaling_tag = 'linear_scaling' if linear_scaling else 'area_scaling'
            cmap_tag = colormap
            for suffix in (f'.{scaling_tag}.{cmap_tag}{ext}',):
                path = candidate + suffix
                if os.path.exists(path):
                    files_written.append(path)

        # 6. Render Bokeh
        if not opts.disable_showing_bokeh or outfile_prefix:
            render_bokeh(
                opts, outfile_prefix, xmin, xmax,
                amino_acids, final_sorted_whitelist,
                circles_bokeh, mutations, hover_texts,
                title_data, xlabel, matrix_name,
                colors, norm, cmap, padded_position2position,
                show_bokeh,
            )
            html_path = f"{outfile_prefix}.{matrix_name}"
            scaling_tag = 'linear_scaling' if linear_scaling else 'area_scaling'
            cmap_tag = colormap
            html_candidate = f"{html_path}.{scaling_tag}.{cmap_tag}.html"
            if os.path.exists(html_candidate):
                files_written.append(html_candidate)

    return {
        'figure': figure,
        'axes': (ax1, ax2, ax3, ax4),
        'dataframe': df,
        'files_written': files_written,
    }


def render_timeline(
    directory: str,
    positions: list[str],
    outfile_prefix: str = '',
    *,
    threshold: float = 0.001,
    colormap: str = 'coolwarm_r',
    matrix: str = 'BLOSUM80',
    show_bokeh: bool = False,
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
        Output path prefix.  When empty, only the ``Figure`` is returned.
    threshold : float
        Minimum frequency to include in the timeline.
    colormap : str
        Matplotlib colormap name (default ``'coolwarm_r'``).
    matrix : str
        BLOSUM matrix name (default ``'BLOSUM80'``).
    show_bokeh : bool
        Open the Bokeh HTML in a browser.
    **extra_options
        Additional options forwarded to the timeline options namespace.

    Returns
    -------
    dict
        - ``'figure'``: ``matplotlib.figure.Figure``
        - ``'data'``: ``TimelineData`` object with all collected points
        - ``'files_written'``: ``list[str]``
    """
    from .mutation_scatter_plot.options import timeline_options
    from .mutation_scatter_plot.core import load_matrix, get_colormap
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
        threshold=threshold,
        colormap=colormap,
        matrix=matrix,
        disable_showing_bokeh=not show_bokeh,
        **extra_options,
    )

    # 1. Discover files
    files = scan_directory(directory, opts.pattern, opts.include_unknown_month)
    if not files:
        raise FileNotFoundError(
            f"No .frequencies.tsv files found in {directory!r}"
        )

    # 2. Parse positions and load scoring
    specs = parse_positions(positions)
    matrix_obj = load_matrix(opts)
    colors, norm, cmap = get_colormap(opts, colormap)

    # 3. Collect data
    data = collect_timeline_data(files, specs, opts, matrix_obj, norm, colors)
    prefix = outfile_prefix or infer_common_prefix(files, directory)

    files_written = []

    # 4. Render
    fig = render_timeline_matplotlib(
        data, opts, norm, cmap, colors, prefix,
    )
    if prefix:
        for ext in ('.png', '.pdf'):
            for candidate in [f"{prefix}.timeline{ext}"]:
                if os.path.exists(candidate):
                    files_written.append(candidate)

    render_timeline_bokeh(data, opts, norm, cmap, colors, prefix)
    html_candidate = f"{prefix}.timeline.html"
    if os.path.exists(html_candidate):
        files_written.append(html_candidate)

    return {
        'figure': fig,
        'data': data,
        'files_written': files_written,
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
) -> pd.DataFrame:
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
    ...     reference_file='MN908947.3_S_full.fasta',
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

    # Read reference
    ref_record = SeqIO.read(reference_file, 'fasta')
    ref_seq = str(ref_record.seq).upper()

    if padded_reference:
        padded_ref_dna = ref_seq
    else:
        padded_ref_dna = ref_seq

    # Translate reference to protein
    from . import alt_translate
    ref_depadded = ref_seq.replace('-', '')
    ref_protein = ''
    for i in range(0, len(ref_depadded) - 2, 3):
        codon = ref_depadded[i:i + 3]
        ref_protein += alt_translate(codon, table=translation_table)

    ref_codons = get_codons(ref_seq)

    # Open output files
    freq_tsv_path = f"{outfile_prefix}.frequencies.tsv"
    unchanged_tsv_path = f"{outfile_prefix}.unchanged_codons.tsv"
    count_path = f"{outfile_prefix}.count"

    outfile = open_file(freq_tsv_path, overwrite)
    outfile_unchanged = open_file(unchanged_tsv_path, overwrite) if print_unchanged_sites else None
    outfile_count = open_file(count_path, overwrite)

    # Write headers
    header = ("padded_position\tposition\toriginal_aa\tmutant_aa\t"
              "frequency\toriginal_codon\tmutant_codon\t"
              "observed_codon_count\ttotal_codons_per_site\n")
    outfile.write(header)
    if outfile_unchanged:
        outfile_unchanged.write(header)

    try:
        parse_alignment(
            opts, alignment_file, padded_ref_dna, ref_protein, ref_codons,
            outfile, outfile_unchanged, outfile_count,
            aa_start=opts.aa_start, min_start=opts.min_start,
            max_stop=opts.max_stop, threads=threads,
            translation_table=translation_table,
        )
    finally:
        outfile.close()
        if outfile_unchanged:
            outfile_unchanged.close()

    # Read back and return as DataFrame
    df = pd.read_csv(freq_tsv_path, sep='\t')
    return df
