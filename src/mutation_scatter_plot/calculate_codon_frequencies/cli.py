"""Command-line interface for calculate_codon_frequencies."""
# This work © 2025-2026 by Jiří Zahradník and Martin Mokrejš
# (First Medical Faculty - Charles University in Prague) is licensed under
# Creative Commons Attribution 4.0 International. To view a copy of this
# license, visit https://creativecommons.org/licenses/by/4.0/

import os
from contextlib import ExitStack
import argparse

from Bio import SeqIO

from . import (
    VERSION,
    get_codons,
    parse_alignment,
    open_file,
)
from .. import alt_translate
from ..profiler import PROFILER


class NoWrapFormatter(argparse.RawDescriptionHelpFormatter, argparse.ArgumentDefaultsHelpFormatter):
    """Help formatter that does not wrap long lines, preserving URLs."""

    def _split_lines(self, text, width):
        return text.splitlines()


def build_option_parser():
    """Build the argument parser for calculate_codon_frequencies."""
    myparser = argparse.ArgumentParser(
        description=__import__('mutation_scatter_plot.calculate_codon_frequencies',
                                fromlist=['']).__doc__,
        formatter_class=NoWrapFormatter,
    )
    myparser.add_argument("--version", action="version", version=f"%(prog)s {VERSION}")
    myparser.add_argument("--reference-infile", action="store", type=str,
        dest="reference_infilename", default=None, metavar="FILE",
        help="FASTA formatted input file with reference padded sequence or not")
    myparser.add_argument("--padded-reference", action="store_true",
        dest="padded_reference", default=False,
        help="By default we do NOT require the reference sequence to be padded with '-' characters to match the alignment delineating INSertions. If it is not padded [default case] then INSertion will not be reported but gaps parsed in the alignment will be skipped as long until 3 nucleotides are available for codon translation. Regardless of this --padded-reference setting, length of the reference sequence must match length of each alignment line.")
    myparser.add_argument("--alignment-file", action="store", type=str,
        dest="alignment_infilename", default=None, metavar="FILE",
        help="ALIGNMENT file in FASTA format with - (minus) chars to adjust the alignment to the --reference-infile")
    myparser.add_argument("--outfile-prefix", action="store", type=str,
        dest="outfileprefix", default=None, metavar="FILE",
        help="It assumes *.frequencies.fasta files. The prefix specified should end with .frequencies . The .tsv and .unchanged_codons.tsv will be appended to the prefix.")
    myparser.add_argument("--left-reference-offset", action="store", type=int,
        dest="left_reference_offset", default=0,
        help="First nucleotide of the ORF region of the REFERENCE of interest to be sliced out from the input sequences. This requires 0-based numbering.")
    myparser.add_argument("--right-reference-offset", action="store", type=int,
        dest="right_reference_offset", default=0,
        help="Last nucleotide of the last codon of the REFERENCE of interest to be sliced out from the input sequences. This requires 0-based numbering.")
    myparser.add_argument("--aa_start", action="store", type=int,
        dest="aa_start", default=0,
        help="Adjust (padded) real position of the very first codon unless (1 for an initiator ATG). This value is added to the codon position reported in the output TSV file (the ATG position minus one). Use this if you cannot use --left-reference-offset nor --right-reference-offset which would have been used for slicing the input reference. The value provided is decremented by one to match pythonic 0-based numbering.")
    myparser.add_argument("--min_start", action="store", type=int,
        dest="min_start", default=0,
        help="Start parsing the alignment since this position of the ALIGNMENT file. This requires 1-based numbering. This is to speedup parsing of input sequences and of the reference by skipping typical leading and trailing padding dashes. Default: 0 (parse since the beginning)")
    myparser.add_argument("--max_stop", action="store", type=int,
        dest="max_stop", default=0,
        help="Stop parsing the alignment at this position of the ALIGNMENT file. This requires 1-based numbering. This is to speedup parsing of input sequences and of the reference by skipping typical leading and trailing padding dashes. Default: 0 (parse until the very end)")
    myparser.add_argument("--x-after-count", action="store_true",
        dest="x_after_count", default=False,
        help="The FASTA file ID contains the count value followed by lowercase 'x'")
    myparser.add_argument("--print-unchanged-sites", action="store_true",
        dest="print_unchanged_sites", default=True,
        help="Print out also sites with unchanged codons in to unchanged_codons.tsv file [default]")
    myparser.add_argument("--disable-print-unchanged-sites", action="store_false",
        dest="print_unchanged_sites",
        help="Do NOT print out sites with unchanged codons to unchanged_codons.tsv file")
    myparser.add_argument("--discard-this-many-leading-nucs", action="store", type=int,
        dest="discard_this_many_leading_nucs", default=0,
        help="Specify how many offending nucleotides are at the front of the FASTA sequences shifting the reading frame of the input FASTA file from frame +1 to either of the two remaining. Count the leading dashes and eventual nucleotides of incomplete codons too and check if it can be divided by 3.0 without slack. By default reading frame +1 is expected and hence no leading nucleotides are discarded. Default: 0")
    myparser.add_argument("--discard-this-many-trailing-nucs", action="store", type=int,
        dest="discard_this_many_trailing_nucs", default=0,
        help="Specify how many offending nucleotides are at the end of each sequence. Default: 0")
    myparser.add_argument("--minimum-alignments-length", action="store", type=int,
        dest="minimum_aln_length", default=50,
        help="Minimum length of aligned NGS read to be used for calculations")
    myparser.add_argument("--debug", action="store", type=int,
        dest="debug", default=0,
        help="Set debug level to some real number")
    myparser.add_argument("--overwrite", action="store_true",
        dest="overwrite", default=False,
        help="Overwrite existing output files instead of raising RuntimeError")

    # Auto-detect default threads from environment variables
    _default_threads = 0
    for _env_var in ["PBS_NUM_PPN", "PBS_NCPUS", "OMP_NUM_THREADS"]:
        if _env_var in os.environ:
            try:
                _default_threads = int(os.environ[_env_var])
                break
            except ValueError:
                continue

    myparser.add_argument("--threads", action="store", type=int,
        dest="threads", default=_default_threads,
        help="Number of CPU worker *processes* to use (via multiprocessing.Pool). "
             "If 0 [default], the tool auto-detects the allocation from environment "
             "variables (tried in order: PBS_NUM_PPN, PBS_NCPUS, OMP_NUM_THREADS) "
             "or falls back to all available cores via multiprocessing.cpu_count(). "
             "Recommended: always set --threads 8 explicitly; benchmarks show 8 "
             "threads is the sweet spot (2x speedup); beyond 16 there are "
             "diminishing returns due to Amdahl's serial fraction. "
             "Note: named '--threads' (not '--jobs') because these are CPU-bound "
             "worker *processes* forked by multiprocessing.Pool to parallelise "
             "1,274 codon-site computations -- process isolation avoids GIL "
             "contention on CPU-heavy workloads. "
             "Contrast with summarize_fasta_pipeline.py --jobs N, which dispatches "
             "I/O-bound file-level tasks via ThreadPoolExecutor (OS threads, not "
             "processes, because FASTA reads release the GIL).")
    myparser.add_argument("--chunksize", action="store", type=int,
        dest="chunksize", default=256,
        help="Number of codon sites dispatched per IPC round-trip to worker processes. "
             "Larger values reduce socket overhead; smaller values improve load-balance. "
             "Default 256 was empirically optimal on a 192-core Xeon (4 workers, 100k rows). "
             "Use 0 to let Python choose automatically.")
    myparser.add_argument("--cpu-bind", choices=['local', 'spread', 'none'],
        dest="cpu_bind", default="local", metavar="POLICY",
        help="NUMA CPU binding policy.  'local' (default) binds all worker "
             "processes to the single NUMA node with the most free RAM — best "
             "for CPU-bound multiprocessing workloads.  'spread' leaves the OS "
             "to use all nodes.  'none' disables autobind.  Autobind is "
             "silently skipped on single-node hosts and when a job scheduler "
             "has already set the CPU affinity (GOMP_CPU_AFFINITY, "
             "KMP_AFFINITY, SLURM_CPU_BIND, SGE_BINDING, OMP_PROC_BIND).")
    myparser.add_argument("--translation-table", action="store", type=int,
        dest="translation_table", default=1, metavar="N",
        help="NCBI genetic code table number to use for translation "
             "(default: 1 = standard genetic code). "
             "See https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi for all tables.")
    return myparser


def main():
    """Main function for calculate_codon_frequencies CLI."""
    myparser = build_option_parser()
    myoptions = myparser.parse_args()

    # parse the reference DNA
    if not myoptions.reference_infilename:
        raise ValueError("Error: Please specify --reference-infile with FASTA sequence")
    if not os.path.exists(myoptions.reference_infilename):
        raise ValueError(f"Error: File {myoptions.reference_infilename} does not exist")
    if os.path.getsize(myoptions.reference_infilename) == 0:
        raise ValueError(f"Error: File {myoptions.reference_infilename} is empty")

    for _record in SeqIO.parse(myoptions.reference_infilename, "fasta"):
        _padded_reference_dna_seq = str(_record.seq)
        break  # parse only the first and supposedly the only entry
    _reference_protein_seq = alt_translate(_padded_reference_dna_seq,
                                           table=myoptions.translation_table)
    _reference_as_codons = get_codons(_padded_reference_dna_seq,
                                      debug=myoptions.debug)

    if not myoptions.alignment_infilename:
        raise RuntimeError("Please specify --alignment-file")
    if not myoptions.outfileprefix:
        raise RuntimeError("Please specify output filename prefix via --outfile-prefix")

    with ExitStack() as stack:
        _count_filename = f"{'.'.join(myoptions.alignment_infilename.split('.')[:-1])}.count"
        _alnfilename_count_handle = stack.enter_context(open(_count_filename, 'w', encoding="utf-8"))

        _tsv_header = "padded_position\tposition\toriginal_aa\tmutant_aa\tfrequency\toriginal_codon\tmutant_codon\tobserved_codon_count\ttotal_codons_per_site\n"

        if myoptions.outfileprefix.endswith('.tsv'):
            _out_name = myoptions.outfileprefix
            _unchanged_name = f"{myoptions.outfileprefix[:-4]}.unchanged_codons.tsv"
        else:
            _out_name = f"{myoptions.outfileprefix}.tsv"
            _unchanged_name = f"{myoptions.outfileprefix}.unchanged_codons.tsv"

        _outfilename_handle = stack.enter_context(
            open_file(_out_name, overwrite=myoptions.overwrite, encoding="utf-8")
        )
        _outfilename_handle.write(_tsv_header)

        if myoptions.print_unchanged_sites:
            _outfilename_unchanged_codons_handle = stack.enter_context(
                open_file(_unchanged_name, overwrite=myoptions.overwrite, encoding="utf-8")
            )
            _outfilename_unchanged_codons_handle.write(_tsv_header)
        else:
            _outfilename_unchanged_codons_handle = None

        _aa_start = (myoptions.aa_start - 1) if myoptions.aa_start else 0
        _min_start = (myoptions.min_start - 1) if myoptions.min_start else 0
        _max_stop  = (myoptions.max_stop  + 1) if myoptions.max_stop  else 0
        _threads   = myoptions.threads if myoptions.threads > 0 else None
        _chunksize = myoptions.chunksize if myoptions.chunksize > 0 else None

        # ── NUMA auto-bind (before multiprocessing Pool creation) ────────────
        # Silently skipped on single-node / non-NUMA hosts and when a
        # job scheduler has already configured CPU placement.
        from .. import numa_bind  # noqa: PLC0415  (import inside function is intentional)
        numa_bind.autobind(mode=myoptions.cpu_bind)

        # parse_alignment manages its own Pool internally so that _WORKER_SHARED
        # (the COW-fork global) is populated *before* the Pool is created.
        # This means each task sends only _pos (4 bytes) over IPC instead of the
        # full alignment arrays (~115 MB), reducing IPC from ~575 MB to ~5 KB.
        if os.path.exists(myoptions.alignment_infilename):
            if os.path.getsize(myoptions.alignment_infilename) == 0:
                raise RuntimeError(f"Input file {myoptions.alignment_infilename} is empty")

            PROFILER.start()
            PROFILER.mark_phase_start("Phase 1: Parse sequences")
            parse_alignment(
                myoptions,
                myoptions.alignment_infilename,
                _padded_reference_dna_seq,
                _reference_protein_seq,
                _reference_as_codons,
                _outfilename_handle,
                _outfilename_unchanged_codons_handle,
                _alnfilename_count_handle,
                _aa_start,
                _min_start,
                _max_stop,
                threads=_threads,
                chunksize=_chunksize,
                translation_table=myoptions.translation_table,
            )

            _prof_sum = PROFILER.pop_phase_summary()
            if _prof_sum:
                print()
                print(_prof_sum)
        else:
            raise RuntimeError(
                f"Input file {str(myoptions.alignment_infilename)} does not exist or is not defined"
            )


if __name__ == "__main__":
    main()
