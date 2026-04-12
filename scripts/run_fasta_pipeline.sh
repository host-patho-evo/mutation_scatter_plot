#!/bin/sh

# run_fasta_pipeline.sh — FASTA processing pipeline: raw download → mutation scatter plots.
#
# Usage:
#   run_fasta_pipeline.sh --infile=spikenuc1207.fasta --reference=MN908947.3_S.fasta [OPTIONS]
#
# Logged run:
#   nohup bash -c "run_fasta_pipeline.sh --infile=spikenuc1207.fasta --reference=MN908947.3_S.fasta 2>&1 | ts '[%Y-%m-%d %H:%M:%S]'" > pipeline_run.log &
#
# Real-time monitoring:
#   while true; do COLUMNS=512 top -c -b -w 512 -u "$USER" -n1 | grep -E 'grep|awk|python|blast|sed' | head -35; sleep 10; done

VERSION="202604112100"
GIT_VERSION=$(git -C "$(dirname "$0")" rev-parse --short=12 HEAD 2>/dev/null || echo "unknown")

set -eo pipefail

# ──────────────────────────────────────────────────────────────────────────────
# Defaults
# ──────────────────────────────────────────────────────────────────────────────
infile=""
reference=""
discard_junk=""
xmin=""
xmax=""
threshold=0.001
jobs=5
sort_bucket_size="40%"
codon_freq_threads=8
compare_frequencies=""
old_alignment_file=""
realign_covid19_spike=false
full_length=""
extract_discarded_fasta=false
full_fasta_header=false
write_original_descr_lines=false
use_nnnx_counts=false
split_gisaid_by_month=false
xranges=""
title=""

# ──────────────────────────────────────────────────────────────────────────────
# Parse command-line arguments
# ──────────────────────────────────────────────────────────────────────────────
usage() {
    cat <<'USAGE'
Usage: run_fasta_pipeline.sh --infile=FILE --reference=FILE [OPTIONS]

Required:
  --infile=FILE             Input FASTA file (e.g. spikenuc1207.fasta).
                            The prefix (basename minus .fasta) is used for all
                            intermediate and output file names.
  --reference=FILE          Reference FASTA file (single entry, e.g. MN908947.3_S.fasta).
                            Its sequence length determines the full-length filter.

Optional:
  --discard-junk=FILE      Plaintext file listing FASTA headers to filter out.
                            If provided, a .no_junk filtering step is performed
                            and all downstream files include a '.no_junk' segment.
                            If omitted, filtering is skipped and filenames are
                            shorter (no '.no_junk' segment) [default: skip].
  --full-length=INT         Override the full-length filter value.  By default,
                            this is the reference sequence length.  Use when the
                            reference is padded but you want to filter for the
                            unpadded length (e.g. ref=3870, filter=3822).
  --xmin=INT                Start codon for scatter plots [default: full range].
  --xmax=INT                End codon for scatter plots [default: full range].
  --xranges=RANGES          Comma-separated min-max pairs for scatter plots,
                            e.g. '1-1274,331-528'.  Overrides --xmin/--xmax.
                            Each range produces a separate set of figures.
                            [default: use --xmin/--xmax or full range].
  --title=STRING            Figure title for scatter plots.  In per-month
                            mode the month tag is appended automatically.
                            [default: auto-derived by mutation_scatter_plot].
  --threshold=FLOAT         Minimum frequency for filtered TSVs [default: 0.001].
  --jobs=INT                Parallel jobs for summarize_fasta_pipeline.py [default: 5].
  --sort-bucket-size=SIZE   Memory fraction for sort buckets [default: 40%].
  --codon-freq-threads=INT  Threads for calculate_codon_frequencies [default: 8].
  --compare-frequencies=FILE
                            Full path to an old frequencies.tsv for comparison.
                            The corresponding .unchanged_codons.tsv is derived
                            automatically. If set, produces .comparison.txt files.
                            [default: skip comparison].
  --old-alignment-file=FILE
                            Full path to a pre-existing cleaned alignment FASTA
                            (e.g. from MAFFT or a previous pipeline run).  When
                            set, the file is symlinked as the
                            .counts.${REFLEN}.clean.fasta and the BLAST pipeline
                            (Stage 3) is skipped entirely.  The downstream
                            summarize_fasta_pipeline.py will later update headers
                            if they lack sha256 checksums [default: run BLAST].
  --version                 Show version and exit.
  --realign-covid19-spike   Enable the SARS-CoV-2 spike-specific InDel
                            realignment step (manual_SARS-CoV-2_InDel_realignment.py)
                            in the BLAST pipeline.  If omitted, the realignment
                            step is skipped, making the pipeline generic for
                            any organism [default: disabled].
  --extract-discarded-fasta Pass --extract-discarded-fasta to summarize_fasta_pipeline.py.
  --full-fasta-header       Pass --full-fasta-header to summarize_fasta_pipeline.py.
  --write-original-descr-lines
                            Pass --write-original-descr-lines to summarize_fasta_pipeline.py.
  --use-nnnx-counts         Pass --use-nnnx-counts to summarize_fasta_pipeline.py.
  --split-gisaid-by-month   After the no_junk filtering step, split the FASTA
                            into per-month files (YYYY-MM) using
                            split_GISAID_sequences_by_month.py.  Each monthly
                            file then runs through the full pipeline (Stages
                            2-9).  The combined file is NOT processed.
                            [default: disabled].
  --help                    Show this help message and exit.
USAGE
    exit "${1:-0}"
}

for arg in "$@"; do
    case "$arg" in
        --infile=*)              infile="${arg#*=}" ;;
        --reference=*)           reference="${arg#*=}" ;;
        --discard-junk=*)        discard_junk="${arg#*=}" ;;
        --xmin=*)                xmin="${arg#*=}" ;;
        --xmax=*)                xmax="${arg#*=}" ;;
        --threshold=*)           threshold="${arg#*=}" ;;
        --jobs=*)                jobs="${arg#*=}" ;;
        --sort-bucket-size=*)    sort_bucket_size="${arg#*=}" ;;
        --codon-freq-threads=*)  codon_freq_threads="${arg#*=}" ;;
        --compare-frequencies=*)    compare_frequencies="${arg#*=}" ;;
        --old-alignment-file=*)     old_alignment_file="${arg#*=}" ;;
        --version)               echo "run_fasta_pipeline.sh  version $VERSION  git:$GIT_VERSION"; exit 0 ;;
        --realign-covid19-spike) realign_covid19_spike=true ;;
        --full-length=*)         full_length="${arg#*=}" ;;
        --extract-discarded-fasta) extract_discarded_fasta=true ;;
        --full-fasta-header)     full_fasta_header=true ;;
        --write-original-descr-lines) write_original_descr_lines=true ;;
        --use-nnnx-counts)       use_nnnx_counts=true ;;
        --split-gisaid-by-month) split_gisaid_by_month=true ;;
        --xranges=*)             xranges="${arg#*=}" ;;
        --title=*)               title="${arg#*=}" ;;
        --help)                  usage 0 ;;
        *)                       echo "Error: unknown argument: $arg" >&2; usage 1 ;;
    esac
done

# Validate required arguments
if [ -z "$infile" ]; then
    echo "Error: --infile is required." >&2
    usage 1
fi
if [ -z "$reference" ]; then
    echo "Error: --reference is required." >&2
    usage 1
fi
if [ ! -f "$infile" ] && [ ! -L "$infile" ]; then
    echo "Error: input file not found: $infile" >&2
    exit 1
fi
if [ ! -f "$reference" ]; then
    echo "Error: reference file not found: $reference" >&2
    exit 1
fi

# ──────────────────────────────────────────────────────────────────────────────
# Derived variables
# ──────────────────────────────────────────────────────────────────────────────
# Prefix: strip .fasta (or .fasta.gz, .fa, etc.) from the input filename.
prefix="$(basename "$infile" .fasta)"
prefix="${prefix%.fa}"   # also handle .fa extension

export reference
export PATH=/auto/vestec1-elixir/projects/biocev/mmokrejs/proj/mutation_scatter_plot/scripts:$PATH

# Compute reference sequence length from the single-entry FASTA.
reference_length=$(grep -v '^>' "$reference" | tr -d '\n' | wc -c)

# Determine full-length filter.
#   --old-alignment-file has highest authority (first sequence length).
#   --full-length must agree if both are given; otherwise it is a hard error.
#   --full-length must also match reference_length if no old alignment is given.
if [ -n "$old_alignment_file" ] && [ -f "$old_alignment_file" -o -L "$old_alignment_file" ]; then
    somelen=$(sed -n '2p' "$old_alignment_file" | tr -d '\n' | wc -c)
    echo "  somelen=$somelen (derived from first sequence in $old_alignment_file)"
    if [ -n "$full_length" ] && [ "$full_length" -ne "$somelen" ]; then
        echo "Error: --full-length=$full_length contradicts --old-alignment-file" >&2
        echo "  The first sequence in $old_alignment_file has length $somelen." >&2
        echo "  Remove --full-length or fix the old alignment file." >&2
        exit 1
    fi
elif [ -n "$full_length" ]; then
    if [ "$full_length" -ne "$reference_length" ]; then
        echo "Error: --full-length=$full_length does not match reference length $reference_length." >&2
        echo "  Fix the reference file or remove --full-length." >&2
        exit 1
    fi
    somelen=$full_length
else
    somelen=$reference_length
fi

# Filtered prefix: includes '.no_junk' segment when junk filtering is active.
if [ -n "$discard_junk" ]; then
    fp="${prefix}.no_junk"
else
    fp="${prefix}"
fi
export fp

echo "═══════════════════════════════════════════════════════════════════════"
echo "  run_fasta_pipeline.sh  version $VERSION  git:$GIT_VERSION"
echo "  Input:     $infile  (prefix: $prefix)"
echo "  Reference: $reference  (length: $reference_length)"
if [ -n "$discard_junk" ]; then
    echo "  Discard junk: $discard_junk  (filtered prefix: $fp)"
else
    echo "  Discard junk: (disabled — no .no_junk filtering)"
fi
echo "  Full-length filter: $somelen"
if [ -n "$xranges" ]; then
    echo "  Plot ranges: $xranges, threshold: $threshold"
elif [ -n "$xmin" ] && [ -n "$xmax" ]; then
    echo "  Plot range: codons ${xmin}–${xmax}, threshold: $threshold"
else
    echo "  Plot range: full range, threshold: $threshold"
fi
echo "  Jobs: $jobs, codon-freq threads: $codon_freq_threads"
if [ -n "$compare_frequencies" ]; then
    echo "  Compare with:  $compare_frequencies"
fi
if [ -n "$old_alignment_file" ]; then
    echo "  Old alignment: $old_alignment_file (skip BLAST)"
fi
if [ "$realign_covid19_spike" = true ]; then
    echo "  SARS-CoV-2 spike realignment: enabled"
fi
if [ "$split_gisaid_by_month" = true ]; then
    echo "  Split by month: enabled (per-month pipeline, combined file skipped)"
fi
if [ -n "$title" ]; then
    echo "  Title: $title"
fi
echo "═══════════════════════════════════════════════════════════════════════"

# ──────────────────────────────────────────────────────────────────────────────
# Stage 1: Filter junk headers and normalize encoding (optional)
# ──────────────────────────────────────────────────────────────────────────────
if [ -n "$discard_junk" ]; then
    # BBTools filterbyname.sh internally crashes on UTF-8 bytes and corrupts
    # standard characters (e.g. í to ýý).  We now normalize the encodings and
    # filter the junk lists simultaneously in a single pure Python native byte-stream.
    if [ ! -e "${fp}.fasta" ]; then
        /usr/bin/time -v -o time_filterbyname.log \
            fix_fasta_encoding.py \
                --infile="${prefix}.fasta" \
                --filter-away-fasta-headers="${discard_junk}" \
                --outfile="${fp}.fasta"
    fi
fi

# ──────────────────────────────────────────────────────────────────────────────
# Stage 1b: Split FASTA by month (optional, GISAID-specific)
# ──────────────────────────────────────────────────────────────────────────────
if $split_gisaid_by_month; then
    # Guard: skip if per-month files already exist
    _existing_months=$(ls -1 "${fp}".????-??.fasta 2>/dev/null | head -1)
    if [ -z "$_existing_months" ]; then
        echo "Info: Splitting ${fp}.fasta by YYYY-MM into per-month files..."
        split_GISAID_sequences_by_month.py \
            --infilename="${fp}.fasta" \
            --outfile-prefix="${fp}"
        echo "Info: Per-month files created."
    else
        echo "Info: Per-month files already exist, skipping split step."
    fi
fi

# ══════════════════════════════════════════════════════════════════════════════
# Helper: run Stages 2–9 for a single FASTA prefix.
# Usage: _run_stages "$fp_current" ["$title_override"]
# ══════════════════════════════════════════════════════════════════════════════
_run_stages() {
    local _fp="$1"
    local _title_override="${2:-}"

    echo "───────────────────────────────────────────────────────────────────────"
    echo "  Processing: ${_fp}.fasta"
    echo "───────────────────────────────────────────────────────────────────────"

    # ──────────────────────────────────────────────────────────────────────
    # Stage 2: Deduplicate identical sequences (NNNNx. counting)
    # ──────────────────────────────────────────────────────────────────────
    if [ ! -e "${_fp}.counts.fasta" ]; then
        /usr/bin/time -v -o "time_count_same_sequences_$(basename "${_fp}").log" \
            count_same_sequences.py \
                --infilename="${_fp}.fasta" \
                --outfile-prefix="${_fp}" \
                --sort-bucket-size="${sort_bucket_size}"
    fi

    # ──────────────────────────────────────────────────────────────────────
    # Stage 3: BLAST + realignment pipeline → cleaned FASTA
    # ──────────────────────────────────────────────────────────────────────
    local _clean_fasta="${_fp}.counts.${somelen}.clean.fasta"
    if [ -n "$old_alignment_file" ] && [ ! -e "$_clean_fasta" ]; then
        if [ ! -f "$old_alignment_file" ] && [ ! -L "$old_alignment_file" ]; then
            echo "Error: old alignment file not found: $old_alignment_file" >&2
            exit 1
        fi
        echo "  Injecting old alignment: $old_alignment_file → $_clean_fasta"
        ln -s "$(realpath "$old_alignment_file")" "$_clean_fasta"
    fi

    if [ ! -e "$_clean_fasta" ]; then
        export somelen

        # Gracefully auto-detect NUMA topology capabilities
        if command -v numactl >/dev/null 2>&1 && numactl --hardware >/dev/null 2>&1 && [ "$(numactl --hardware | grep -c '^node [0-9]')" -gt 1 ]; then
            export BEST_NODE=$(numactl --hardware | awk '/free:/ {print $2, $4}' | sort -k2,2rn | head -n1 | awk '{print $1}')
            export NODE_CORES=$(numactl --hardware | grep -E "^node ${BEST_NODE} cpus:" | awk '{print NF-3}')
            NUMACTL_CMD="numactl --cpunodebind=${BEST_NODE} --localalloc "
        else
            export BEST_NODE="N/A"
            export NODE_CORES=$(nproc 2>/dev/null || echo 8)
            NUMACTL_CMD=""
        fi

        # Reserve cores for downstream pipeline operators.
        if [ "$NODE_CORES" -gt 24 ]; then
            export REALIGN_CORES=$(( (NODE_CORES * 3) / 10 ))
            export BLASTN_CORES=$((NODE_CORES - REALIGN_CORES - 4))
        elif [ "$NODE_CORES" -gt 8 ]; then
            export REALIGN_CORES=$((NODE_CORES / 4))
            export BLASTN_CORES=$((NODE_CORES - REALIGN_CORES - 2))
        elif [ "$NODE_CORES" -gt 4 ]; then
            export REALIGN_CORES=1
            export BLASTN_CORES=$((NODE_CORES - 2))
        else
            export REALIGN_CORES=1
            export BLASTN_CORES=1
        fi

        # Common pipeline prefix (written to temp script).
        # The heredoc uses unquoted PIPEEOF so shell variables are expanded
        # at write time (they are all known at this point).
        local _pipe_script="tmp_blastn_pipe_${somelen}_$(basename "${_fp}").sh"
        cat << PIPEEOF > "$_pipe_script"
#!/bin/bash
blastn -task blastn -reward 2 -max_hsps 1 -num_alignments 1 -word_size 4 -num_threads "${BLASTN_CORES}" -dust no -evalue 1e-5 -outfmt "6 qacc sstart send sstrand qstart qend slen qlen evalue bitscore score length pident nident mismatch positive gapopen gaps ppos qseq sseq" -db "$reference" -query "${_fp}.counts.fasta" 2>/dev/null | awk '{print \$1,\$2,\$3,\$4,\$5,\$6,\$7,\$8,\$9,\$10,\$11,\$12,\$13,\$14,\$15,\$16,\$17,\$18,\$19,\$20,\$21}' | drop_erroneous_insertions.py --infile=- --outfile=- --drop-misinsertions | reversecomplement_reads_on_minus.py --infile=- --reference="$reference" --min_start=1 --max_stop="${somelen}" \\
PIPEEOF

        # Append the tail: realignment step or direct output.
        if [ "$realign_covid19_spike" = true ]; then
            echo "| manual_SARS-CoV-2_InDel_realignment.py --cores \"${REALIGN_CORES}\" > \"${_fp}.counts.${somelen}.clean.fasta\"" >> "$_pipe_script"
        else
            echo "> \"${_fp}.counts.${somelen}.clean.fasta\"" >> "$_pipe_script"
        fi

        /usr/bin/time -v -o "time_blastn_pipeline_${somelen}_$(basename "${_fp}").log" ${NUMACTL_CMD} bash "$_pipe_script"
    fi

    # ──────────────────────────────────────────────────────────────────────
    # Stage 4: Split by length + smart symlink deduplication
    # ──────────────────────────────────────────────────────────────────────
    local _parent="${_fp}.counts.${somelen}.clean.fasta"
    local _exactly="${_fp}.counts.${somelen}.clean.exactly_${somelen}.fasta"
    local _shorter="${_fp}.counts.${somelen}.clean.shorter_${somelen}.fasta"
    local _longer="${_fp}.counts.${somelen}.clean.longer_${somelen}.fasta"
    local _combined="${_fp}.counts.${somelen}.clean.exactly_or_shorter_${somelen}.fasta"

    if [ ! -e "$_combined" ]; then
        /usr/bin/time -v -o "time_split_fasta_entries_by_lengths_${somelen}_$(basename "${_fp}").log" \
            split_fasta_entries_by_lengths \
                --infile="$_parent" \
                --outfile-prefix="${_fp}.counts.${somelen}.clean" \
                --full-length="${somelen}" \
                --convert-to-upper

        if [ ! -s "$_longer" ] && [ ! -s "$_shorter" ]; then
            echo "  split: exactly_${somelen} ≡ parent (shorter and longer both empty) → symlink"
            rm -f "$_exactly"
            ln -s "$(basename "$_parent")" "$_exactly"
        fi

        if [ ! -s "$_shorter" ]; then
            echo "  split: exactly_or_shorter_${somelen} = exactly_${somelen} (shorter empty) → symlink"
            ln -s "$(basename "$_exactly")" "$_combined"
        elif [ ! -s "$_exactly" ]; then
            echo "  split: exactly_or_shorter_${somelen} = shorter_${somelen} (exactly empty) → symlink"
            ln -s "$(basename "$_shorter")" "$_combined"
        else
            echo "  split: cat exactly_${somelen} + shorter_${somelen} → exactly_or_shorter_${somelen}"
            cat "$_exactly" "$_shorter" > "$_combined"
        fi

        if [ ! -s "$_longer" ]; then
            echo "  split: exactly_or_shorter_${somelen} ≡ parent (longer empty)"
        fi
    fi

    # ──────────────────────────────────────────────────────────────────────
    # Stage 5: Protein translation
    # ──────────────────────────────────────────────────────────────────────
    if [ ! -e "${_fp}.counts.${somelen}.clean.exactly_or_shorter_${somelen}.prot.fasta" ]; then
        echo "translate.py  invoked" >&2
        /usr/bin/time -v -o "time_translate_$(basename "${_fp}").log" \
            translate.py \
                --infile="${_fp}.counts.${somelen}.clean.exactly_or_shorter_${somelen}.fasta" \
                --respect-alignment \
                --outfile="${_fp}.counts.${somelen}.clean.exactly_or_shorter_${somelen}.prot.fasta"
    fi

    # ──────────────────────────────────────────────────────────────────────
    # Stage 6: Codon frequency calculation
    # ──────────────────────────────────────────────────────────────────────
    local _freq_prefix="${_fp}.counts.${somelen}.clean.exactly_or_shorter_${somelen}"

    if [ ! -e "${_freq_prefix}.frequencies.tsv" ]; then
        echo "calculate_codon_frequencies  invoked" >&2
        /usr/bin/time -v -o "time_calculate_codon_frequencies_${somelen}_$(basename "${_fp}").log" \
            calculate_codon_frequencies \
                --reference-infile="$reference" \
                --alignment-file="${_freq_prefix}.fasta" \
                --print-unchanged-sites \
                --outfile-prefix="${_freq_prefix}.frequencies" \
                --x-after-count \
                --padded-reference \
                --threads "$codon_freq_threads"
    fi

    # ──────────────────────────────────────────────────────────────────────
    # Stage 7: Frequency filtering
    # ──────────────────────────────────────────────────────────────────────
    local _threshold_tag
    _threshold_tag=$(echo "$threshold" | tr -d '.')

    awk -v t="$threshold" '{if ($5 >= t) print}' "${_freq_prefix}.frequencies.tsv" \
        > "${_freq_prefix}.frequencies.${_threshold_tag}.tsv"
    awk -v t="$threshold" '{if ($5 >= t) print}' "${_freq_prefix}.frequencies.unchanged_codons.tsv" \
        > "${_freq_prefix}.frequencies.unchanged_codons.${_threshold_tag}.tsv"

    # ──────────────────────────────────────────────────────────────────────
    # Stage 9: Mutation scatter plots
    # ──────────────────────────────────────────────────────────────────────
    local _interactive="--disable-showing-bokeh --disable-showing-mplcursors"

    # Build the list of (xmin, xmax) ranges to iterate over.
    local _ranges=""
    if [ -n "$xranges" ]; then
        _ranges="$xranges"
    elif [ -n "$xmin" ] && [ -n "$xmax" ]; then
        _ranges="${xmin}-${xmax}"
    fi

    if [ -n "$_ranges" ]; then
        # Loop over each comma-separated range
        echo "$_ranges" | tr ',' '\n' | while IFS='-' read -r _rmin _rmax; do
            local _aa_suffix=".aa${_rmin}-${_rmax}"
            local _codon_suffix=".codon${_rmin}-${_rmax}"
            local _xrange_args="--xmin $_rmin --xmax $_rmax"

            # Use wider tick spacing for full-protein ranges (>500 codons)
            local _tick_args=""
            if [ $(( _rmax - _rmin )) -gt 500 ]; then
                _tick_args="--x-axis-major-ticks-spacing=40 --x-axis-minor-ticks-spacing=20"
            fi

            local _title_arg=""
            if [ -n "$_title_override" ]; then
                _title_arg="--title=$_title_override"
            fi

            for _scaling in '--linear-circle-size' ''; do
                mutation_scatter_plot $_scaling \
                    --tsv="${_freq_prefix}.frequencies.tsv" \
                    --outfile-prefix="${_freq_prefix}${_aa_suffix}" \
                    --aminoacids --show-DEL --show-INS \
                    $_xrange_args --threshold "$threshold" \
                    $_tick_args ${_title_arg:+"$_title_arg"} \
                    --colormap=coolwarm_r $_interactive

                mutation_scatter_plot $_scaling \
                    --tsv="${_freq_prefix}.frequencies.tsv" \
                    --outfile-prefix="${_freq_prefix}${_codon_suffix}" \
                    --show-DEL --show-INS \
                    $_xrange_args --threshold "$threshold" \
                    $_tick_args ${_title_arg:+"$_title_arg"} \
                    --include-synonymous --colormap=coolwarm_r $_interactive

                mutation_scatter_plot $_scaling \
                    --tsv="${_freq_prefix}.frequencies.tsv" \
                    --outfile-prefix="${_freq_prefix}${_codon_suffix}" \
                    --show-DEL --show-INS \
                    $_xrange_args --threshold "$threshold" \
                    $_tick_args ${_title_arg:+"$_title_arg"} \
                    --include-synonymous $_interactive

                mutation_scatter_plot $_scaling \
                    --tsv="${_freq_prefix}.frequencies.tsv" \
                    --outfile-prefix="${_freq_prefix}${_aa_suffix}" \
                    --aminoacids --show-DEL --show-INS \
                    $_xrange_args --threshold "$threshold" \
                    $_tick_args ${_title_arg:+"$_title_arg"} \
                    $_interactive
            done
    done
    else
        # No ranges specified — render full range.
        # Extract natural min/max from TSV column 2 (natural position).
        local _nat_range=""
        if [ -f "${_freq_prefix}.frequencies.unchanged_codons.tsv" ]; then
            _nat_range=$(awk 'NR>1 {if(min=="" || $2+0<min+0) min=$2; if($2+0>max+0) max=$2} END {print min, max}' \
                         "${_freq_prefix}.frequencies.unchanged_codons.tsv")
        elif [ -f "${_freq_prefix}.frequencies.tsv" ]; then
            _nat_range=$(awk 'NR>1 {if(min=="" || $2+0<min+0) min=$2; if($2+0>max+0) max=$2} END {print min, max}' \
                         "${_freq_prefix}.frequencies.tsv")
        fi

        local _aa_suffix=".aa"
        local _codon_suffix=".codon"
        if [ -n "$_nat_range" ] && [ "$_nat_range" != " " ]; then
            local _nat_min=${_nat_range% *}
            local _nat_max=${_nat_range#* }
            _aa_suffix=".aa${_nat_min}-${_nat_max}"
            _codon_suffix=".codon${_nat_min}-${_nat_max}"
        fi

        for _scaling in '--linear-circle-size' ''; do
            mutation_scatter_plot $_scaling \
                --tsv="${_freq_prefix}.frequencies.tsv" \
                --outfile-prefix="${_freq_prefix}${_aa_suffix}" \
                --aminoacids --show-DEL --show-INS \
                --threshold "$threshold" \
                --colormap=coolwarm_r $_interactive

            mutation_scatter_plot $_scaling \
                --tsv="${_freq_prefix}.frequencies.tsv" \
                --outfile-prefix="${_freq_prefix}${_codon_suffix}" \
                --show-DEL --show-INS \
                --threshold "$threshold" \
                --include-synonymous --colormap=coolwarm_r $_interactive

            mutation_scatter_plot $_scaling \
                --tsv="${_freq_prefix}.frequencies.tsv" \
                --outfile-prefix="${_freq_prefix}${_codon_suffix}" \
                --show-DEL --show-INS \
                --threshold "$threshold" \
                --include-synonymous $_interactive

            mutation_scatter_plot $_scaling \
                --tsv="${_freq_prefix}.frequencies.tsv" \
                --outfile-prefix="${_freq_prefix}${_aa_suffix}" \
                --aminoacids --show-DEL --show-INS \
                --threshold "$threshold" \
                $_interactive
        done
    fi
}

# ══════════════════════════════════════════════════════════════════════════════
# Main dispatch: per-month loop vs. single combined processing
# ══════════════════════════════════════════════════════════════════════════════
if $split_gisaid_by_month; then
    # ──────────────────────────────────────────────────────────────────────
    # Per-month processing: run Stages 2–9 for each monthly file
    # ──────────────────────────────────────────────────────────────────────
    for month_fasta in "${fp}".????-??.fasta; do
        [ -f "$month_fasta" ] || continue
        fp_month="${month_fasta%.fasta}"
        month_tag="${fp_month##*.}"   # e.g. "2025-06"
        if [ -n "$title" ]; then
            _month_title="${title} (${month_tag})"
        else
            _month_title="Mutations in SARS-CoV-2 per month in GISAID ${prefix} (${month_tag})"
        fi
        _run_stages "$fp_month" "$_month_title"
    done

    # ──────────────────────────────────────────────────────────────────────
    # Stage 12: Animated GIF assembly (excluding YYYY-00 months)
    # ──────────────────────────────────────────────────────────────────────
    if command -v magick >/dev/null 2>&1; then
        echo "Info: Assembling animated GIFs from per-month PNGs..."
        # Build the list of ranges for GIF assembly
        _gif_ranges=""
        if [ -n "$xranges" ]; then
            _gif_ranges="$xranges"
        elif [ -n "$xmin" ] && [ -n "$xmax" ]; then
            _gif_ranges="${xmin}-${xmax}"
        fi

        if [ -n "$_gif_ranges" ]; then
            echo "$_gif_ranges" | tr ',' '\n' | while IFS='-' read -r _rmin _rmax; do
                for _mode in "aa${_rmin}-${_rmax}" "codon${_rmin}-${_rmax}"; do
                    for _stype in area_scaling linear_scaling; do
                        for _cmap in coolwarm_r amino_acid_changes; do
                            # Collect PNGs for this specific (mode, scaling, colormap)
                            # combination, sorted by YYYY-MM, excluding month-00.
                            _files=$(ls -1 ${fp}.????-??.counts.*.${_mode}.*.${_stype}.${_cmap}.png 2>/dev/null \
                                     | grep -v '\.[0-9]\{4\}-00\.' | sort)
                            if [ -n "$_files" ]; then
                                _outgif="${fp}.${_mode}.BLOSUM80.${_stype}.${_cmap}.animated.gif"
                                echo "Info: Creating $_outgif"
                                magick -dispose previous -delay 300 $_files gif:"$_outgif"
                            fi
                        done
                    done
                done
            done
        fi
    else
        echo "Warning: 'magick' (ImageMagick) not found, skipping GIF assembly." >&2
    fi
else
    # ──────────────────────────────────────────────────────────────────────
    # Standard combined-file processing
    # ──────────────────────────────────────────────────────────────────────
    _run_stages "$fp" "$title"

    # ──────────────────────────────────────────────────────────────────────
    # Stage 10: Pipeline summary + traceability
    # ──────────────────────────────────────────────────────────────────────
    summarize_args=""
    [ "$extract_discarded_fasta" = true ] && summarize_args="$summarize_args --extract-discarded-fasta"
    [ "$full_fasta_header" = true ] && summarize_args="$summarize_args --full-fasta-header"
    [ "$write_original_descr_lines" = true ] && summarize_args="$summarize_args --write-original-descr-lines"
    [ "$use_nnnx_counts" = true ] && summarize_args="$summarize_args --use-nnnx-counts"

    summarize_fasta_pipeline.py . "$prefix" \
        $summarize_args \
        --jobs "$jobs"
fi

# ──────────────────────────────────────────────────────────────────────────────
# Stage 11: Post-pipeline deduplication
# ──────────────────────────────────────────────────────────────────────────────
# Replace verified-identical files with symlinks.  Each call is guarded:
# only runs if the target exists and the file is not already a symlink
# (idempotent on re-runs).

_symlink_if_dup() {
    local target="$1" link="$2"
    if [ -f "$target" ] && [ ! -L "$link" ]; then
        rm -f "$link"
        ln -s "$target" "$link"
        echo "  dedup: $link → $target"
    fi
}

# Group A: effectively_used_original_entries.fasta (.no_junk vs .counts)
_symlink_if_dup "${fp}.effectively_used_original_entries.fasta" \
                "${fp}.counts.effectively_used_original_entries.fasta"

# Group C: traceability TSVs (shorter/longer discarded + prot survived)
# All identical when shorter and longer split files are empty.
_symlink_if_dup "${fp}.counts.${somelen}.clean.exactly_or_shorter_${somelen}.prot.sha256_to_original_descr_lines_of_survived.tsv" \
                "${fp}.counts.${somelen}.clean.longer_${somelen}.sha256_to_original_descr_lines_of_discarded.tsv"
_symlink_if_dup "${fp}.counts.${somelen}.clean.exactly_or_shorter_${somelen}.prot.sha256_to_original_descr_lines_of_survived.tsv" \
                "${fp}.counts.${somelen}.clean.shorter_${somelen}.sha256_to_original_descr_lines_of_discarded.tsv"

# Group D: survived descr_lines TSVs (.no_junk vs .counts)
# Identical when .counts deduplication discards 0 entries.
_symlink_if_dup "${fp}.sha256_to_original_descr_lines_of_survived.tsv" \
                "${fp}.counts.sha256_to_original_descr_lines_of_survived.tsv"

# Group F: effectively_used_original_ids.txt (exactly vs exactly_or_shorter)
# Identical when shorter is empty.
_symlink_if_dup "${fp}.counts.${somelen}.clean.exactly_${somelen}.effectively_used_original_ids.txt" \
                "${fp}.counts.${somelen}.clean.exactly_or_shorter_${somelen}.effectively_used_original_ids.txt"

# Group G: effectively_used_sha256_hashes.txt
# Three identical when clean ≡ exactly ≡ exactly_or_shorter.
_symlink_if_dup "${fp}.counts.${somelen}.clean.effectively_used_sha256_hashes.txt" \
                "${fp}.counts.${somelen}.clean.exactly_or_shorter_${somelen}.effectively_used_sha256_hashes.txt"
_symlink_if_dup "${fp}.counts.${somelen}.clean.effectively_used_sha256_hashes.txt" \
                "${fp}.counts.${somelen}.clean.exactly_${somelen}.effectively_used_sha256_hashes.txt"

echo "  dedup: done"

# ──────────────────────────────────────────────────────────────────────────────
# Notes
# ──────────────────────────────────────────────────────────────────────────────
# The table columns with --classify-mismatches now read:
#   'Altered due to end-clipping'  'Altered inside the sequence'
# Both names are self-explanatory: the first covers expected alignment trimming,
# and the second flags anything where bases inside the sequence were actually
# changed — the category you'd want to investigate further.
#
# Get original header lines for a sha256 checksum:
#   awk -F'\t' '$1 == "SHA256HEX" { for (i=1; i<=NF; i++) print $i }' \
#       ${fp}.counts.${somelen}.clean.exactly_or_shorter_${somelen}.sha256_to_original_descr_lines_of_survived.tsv | wc -l
