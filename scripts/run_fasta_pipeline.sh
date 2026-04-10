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
#   while true; do COLUMNS=512 top -c -b -w 512 -u "$USER" -n1 | grep -E 'grep|awk|python|blast|sed' | head -15; sleep 10; done

VERSION="202604101800"
GIT_VERSION=$(git -C "$(dirname "$0")" rev-parse --short=12 HEAD 2>/dev/null || echo "unknown")

set -eo pipefail

# ──────────────────────────────────────────────────────────────────────────────
# Defaults
# ──────────────────────────────────────────────────────────────────────────────
infile=""
reference=""
discard_junk=""
xmin=430
xmax=528
threshold=0.001
jobs=5
sort_bucket_size="40%"
codon_freq_threads=8
compare_frequencies=""
old_alignment_file=""
realign_covid19_spike=false

# ──────────────────────────────────────────────────────────────────────────────
# Parse command-line arguments
# ──────────────────────────────────────────────────────────────────────────────
usage() {
    cat <<'USAGE'
Usage: processing5.sh --infile=FILE --reference=FILE [OPTIONS]

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
  --xmin=INT                Start codon for scatter plots [default: 430].
  --xmax=INT                End codon for scatter plots [default: 528].
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
somelen=$reference_length

# Filtered prefix: includes '.no_junk' segment when junk filtering is active.
if [ -n "$discard_junk" ]; then
    fp="${prefix}.no_junk"
else
    fp="${prefix}"
fi

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
echo "  Plot range: codons ${xmin}–${xmax}, threshold: $threshold"
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
# Stage 2: Deduplicate identical sequences (NNNNx. counting)
# ──────────────────────────────────────────────────────────────────────────────
if [ ! -e "${fp}.counts.fasta" ]; then
    /usr/bin/time -v -o time_count_same_sequences.log \
        count_same_sequences.py \
            --infilename="${fp}.fasta" \
            --outfile-prefix="${fp}" \
            --sort-bucket-size="${sort_bucket_size}"
fi

# ──────────────────────────────────────────────────────────────────────────────
# Stage 3: BLAST + realignment pipeline → cleaned FASTA
# ──────────────────────────────────────────────────────────────────────────────
# If --old-alignment-file was provided, symlink it instead of running BLAST.
clean_fasta="${fp}.counts.${somelen}.clean.fasta"
if [ -n "$old_alignment_file" ] && [ ! -e "$clean_fasta" ]; then
    if [ ! -f "$old_alignment_file" ] && [ ! -L "$old_alignment_file" ]; then
        echo "Error: old alignment file not found: $old_alignment_file" >&2
        exit 1
    fi
    echo "  Injecting old alignment: $old_alignment_file → $clean_fasta"
    ln -s "$(realpath "$old_alignment_file")" "$clean_fasta"
fi

if [ ! -e "$clean_fasta" ]; then
    export somelen

    # ────────────────────────────────────────────────────────────────────────
    # PIPELINE SCALING & NUMA PERFORMANCE NOTES:
    # ────────────────────────────────────────────────────────────────────────
    # Originally, the `blastn` command was explicitly hardcoded to request
    # `-num_threads 64` natively on a 192 CPU NUMA host, completely without
    # assessing its actual hardware partitioning.  Because each isolated NUMA
    # node on this platform physically caps out strictly at 48 hardware cores,
    # spinning 64 threads inherently forced execution to brutally cross NUMA
    # socket boundaries via sheer oversubscription.
    #
    # Performance Solution:
    # We now dynamically parse hardware footprints on-the-fly to isolate the
    # single NUMA node holding the largest chunk of free memory and cap blastn
    # threads to fit inside that node.
    # ────────────────────────────────────────────────────────────────────────

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
    cat << 'PIPEEOF' > tmp_blastn_pipe_"${somelen}".sh
#!/bin/bash
blastn -task blastn -reward 2 -max_hsps 1 -num_alignments 1 -word_size 4 -num_threads "${BLASTN_CORES}" -dust no -evalue 1e-5 -outfmt "6 qacc sstart send sstrand qstart qend slen qlen evalue bitscore score length pident nident mismatch positive gapopen gaps ppos qseq sseq" -db "$reference" -query "${fp}.counts.fasta" 2>/dev/null | awk '{print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18,$19,$20,$21}' | drop_erroneous_insertions.py --infile=- --outfile=- --drop-misinsertions | reversecomplement_reads_on_minus.py --infile=- --reference="$reference" --min_start=1 --max_stop="${somelen}" \
PIPEEOF

    # Append the tail: realignment step or direct output.
    if [ "$realign_covid19_spike" = true ]; then
        echo '| manual_SARS-CoV-2_InDel_realignment.py --cores "${REALIGN_CORES}" > "${fp}.counts.${somelen}.clean.fasta"' >> tmp_blastn_pipe_"${somelen}".sh
    else
        echo '> "${fp}.counts.${somelen}.clean.fasta"' >> tmp_blastn_pipe_"${somelen}".sh
    fi

    /usr/bin/time -v -o time_blastn_pipeline_"${somelen}".log ${NUMACTL_CMD} bash tmp_blastn_pipe_"${somelen}".sh
fi

# ──────────────────────────────────────────────────────────────────────────────
# Stage 4: Split by length + smart symlink deduplication
# ──────────────────────────────────────────────────────────────────────────────
parent="${fp}.counts.${somelen}.clean.fasta"
exactly="${fp}.counts.${somelen}.clean.exactly_${somelen}.fasta"
shorter="${fp}.counts.${somelen}.clean.shorter_${somelen}.fasta"
longer="${fp}.counts.${somelen}.clean.longer_${somelen}.fasta"
combined="${fp}.counts.${somelen}.clean.exactly_or_shorter_${somelen}.fasta"

if [ ! -e "$combined" ]; then
    /usr/bin/time -v -o time_split_fasta_entries_by_lengths_"${somelen}".log \
        split_fasta_entries_by_lengths \
            --infile="$parent" \
            --outfile-prefix="${fp}.counts.${somelen}.clean" \
            --full-length="${somelen}" \
            --convert-to-upper

    # Smart deduplication after splitting:
    #   exactly + shorter + longer = parent  (partition by length)
    #
    # If longer is empty AND shorter is empty → exactly ≡ parent → symlink.
    # For combining exactly + shorter:
    #   - If shorter is empty → exactly_or_shorter = exactly → symlink.
    #   - If exactly is empty → exactly_or_shorter = shorter → symlink.
    #   - Otherwise → cat (genuine concatenation needed).

    # Step 1: If a split file equals the parent, symlink it.
    if [ ! -s "$longer" ] && [ ! -s "$shorter" ]; then
        echo "  split: exactly_${somelen} ≡ parent (shorter and longer both empty) → symlink"
        rm -f "$exactly"
        ln -s "$(basename "$parent")" "$exactly"
    fi

    # Step 2: Build exactly_or_shorter — avoid cat when possible.
    if [ ! -s "$shorter" ]; then
        echo "  split: exactly_or_shorter_${somelen} = exactly_${somelen} (shorter empty) → symlink"
        ln -s "$(basename "$exactly")" "$combined"
    elif [ ! -s "$exactly" ]; then
        echo "  split: exactly_or_shorter_${somelen} = shorter_${somelen} (exactly empty) → symlink"
        ln -s "$(basename "$shorter")" "$combined"
    else
        echo "  split: cat exactly_${somelen} + shorter_${somelen} → exactly_or_shorter_${somelen}"
        cat "$exactly" "$shorter" > "$combined"
    fi

    if [ ! -s "$longer" ]; then
        echo "  split: exactly_or_shorter_${somelen} ≡ parent (longer empty)"
    fi
fi

# ──────────────────────────────────────────────────────────────────────────────
# Stage 5: Protein translation
# ──────────────────────────────────────────────────────────────────────────────
find . -name "*.clean.exactly_or_shorter_${somelen}.fasta" | while read f; do
    p=$(echo "$f" | sed -e 's/.fasta//')
    if [ ! -e "$p".prot.fasta ]; then
        echo "translate.py  invoked: translate.py --infile=\"$f\" --respect-alignment --outfile=\"$p\".prot.fasta" >&2
        nohup /usr/bin/time -v -o time_translate_"$(basename "$p")".log \
            translate.py --infile="$f" --respect-alignment --outfile="$p".prot.fasta &
    fi
done
wait  # wait for background translate.py jobs

# ──────────────────────────────────────────────────────────────────────────────
# Stage 6: Codon frequency calculation
# ──────────────────────────────────────────────────────────────────────────────
freq_prefix="${fp}.counts.${somelen}.clean.exactly_or_shorter_${somelen}"

if [ ! -e "${freq_prefix}.frequencies.tsv" ]; then
    echo "calculate_codon_frequencies  invoked" >&2
    /usr/bin/time -v -o time_calculate_codon_frequencies_"${somelen}".log \
        calculate_codon_frequencies \
            --reference-infile="$reference" \
            --alignment-file="${freq_prefix}.fasta" \
            --print-unchanged-sites \
            --outfile-prefix="${freq_prefix}.frequencies" \
            --x-after-count \
            --padded-reference \
            --threads "$codon_freq_threads"
fi

# ──────────────────────────────────────────────────────────────────────────────
# Stage 7: Frequency filtering
# ──────────────────────────────────────────────────────────────────────────────
# Convert threshold to a filename-safe string (e.g. 0.001 → 0001)
threshold_tag=$(echo "$threshold" | tr -d '.')

awk -v t="$threshold" '{if ($5 >= t) print}' "${freq_prefix}.frequencies.tsv" \
    > "${freq_prefix}.frequencies.${threshold_tag}.tsv"
awk -v t="$threshold" '{if ($5 >= t) print}' "${freq_prefix}.frequencies.unchanged_codons.tsv" \
    > "${freq_prefix}.frequencies.unchanged_codons.${threshold_tag}.tsv"

# ──────────────────────────────────────────────────────────────────────────────
# Stage 8: Comparison with old frequencies (optional)
# ──────────────────────────────────────────────────────────────────────────────
if [ -n "$compare_frequencies" ]; then
    # Derive the unchanged_codons variant from the given frequencies file.
    old_freqs="$compare_frequencies"
    old_unchanged="$(echo "$compare_frequencies" | sed 's/\.tsv$/.unchanged_codons.tsv/')"

    if [ -f "$old_freqs" ]; then
        for n in $(seq "$xmin" "$xmax"); do
            echo ""
            echo "$n OLD $old_freqs"
            awk -F'\t' '$2=='"$n" "$old_freqs"
            echo "$n NEW ${freq_prefix}.frequencies.${threshold_tag}.tsv"
            awk -F'\t' '$2=='"$n" "${freq_prefix}.frequencies.${threshold_tag}.tsv"
        done > "${freq_prefix}.frequencies.comparison.txt"
    fi

    if [ -f "$old_unchanged" ]; then
        for n in $(seq "$xmin" "$xmax"); do
            echo ""
            echo "$n OLD $old_unchanged"
            awk -F'\t' '$2=='"$n" "$old_unchanged"
            echo "$n NEW ${freq_prefix}.frequencies.unchanged_codons.${threshold_tag}.tsv"
            awk -F'\t' '$2=='"$n" "${freq_prefix}.frequencies.unchanged_codons.${threshold_tag}.tsv"
        done > "${freq_prefix}.frequencies.unchanged_codons.comparison.txt"
    fi
fi

# ──────────────────────────────────────────────────────────────────────────────
# Stage 9: Mutation scatter plots
# ──────────────────────────────────────────────────────────────────────────────
interactive="--disable-showing-bokeh --disable-showing-mplcursors"
for scaling in '--linear-circle-size' ''; do
    # Amino acid plots with coolwarm_r
    mutation_scatter_plot $scaling \
        --tsv="${freq_prefix}.frequencies.tsv" \
        --outfile-prefix="${freq_prefix}.aa${xmin}-aa${xmax}" \
        --aminoacids --show-DEL --show-INS \
        --xmin "$xmin" --xmax "$xmax" --threshold "$threshold" \
        --colormap=coolwarm_r $interactive

    # Codon plots with coolwarm_r
    mutation_scatter_plot $scaling \
        --tsv="${freq_prefix}.frequencies.tsv" \
        --outfile-prefix="${freq_prefix}.codon${xmin}-codon${xmax}" \
        --show-DEL --show-INS \
        --xmin "$xmin" --xmax "$xmax" --threshold "$threshold" \
        --include-synonymous --colormap=coolwarm_r $interactive

    # Codon plots with default colormap (amino_acid_changes)
    mutation_scatter_plot $scaling \
        --tsv="${freq_prefix}.frequencies.tsv" \
        --outfile-prefix="${freq_prefix}.codon${xmin}-codon${xmax}" \
        --show-DEL --show-INS \
        --xmin "$xmin" --xmax "$xmax" --threshold "$threshold" \
        --include-synonymous $interactive

    # Amino acid plots with default colormap
    mutation_scatter_plot $scaling \
        --tsv="${freq_prefix}.frequencies.tsv" \
        --outfile-prefix="${freq_prefix}.aa${xmin}-aa${xmax}" \
        --aminoacids --show-DEL --show-INS \
        --xmin "$xmin" --xmax "$xmax" --threshold "$threshold" \
        $interactive
done

# ──────────────────────────────────────────────────────────────────────────────
# Stage 10: Pipeline summary + traceability
# ──────────────────────────────────────────────────────────────────────────────
summarize_fasta_pipeline.py . "$prefix" \
    --extract-discarded-fasta \
    --full-fasta-header \
    --write-original-descr-lines \
    --use-nnnx-counts \
    --jobs "$jobs"

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
