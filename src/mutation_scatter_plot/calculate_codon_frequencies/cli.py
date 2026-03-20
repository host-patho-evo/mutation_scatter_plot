# This work © 2025 by Jiří Zahradník and Martin Mokrejš
# (First Medical Faculty - Charles University in Prague) is licensed under
# Creative Commons Attribution 4.0 International. To view a copy of this
# license, visit https://creativecommons.org/licenses/by/4.0/

import os

from optparse import OptionParser, IndentedHelpFormatter

from Bio import SeqIO

from . import (
    VERSION,
    get_codons,
    parse_alignment,
    open_file,
)
from ..utils import alt_translate


class NoWrapFormatter(IndentedHelpFormatter):
    """Help formatter that does not wrap long lines, preserving URLs."""

    def format_description(self, description):
        return f"{description}\n" if description else ""

    def format_option(self, option):
        result: list[str] = []
        opts = self.option_strings[option]
        opt_width = self.help_position - self.current_indent - 2
        if len(opts) > opt_width:
            opts = "%*s%s\n" % (self.current_indent, "", opts)
            indent_first = self.help_position
        else:
            opts = "%*s%-*s  " % (self.current_indent, "", opt_width, opts)
            indent_first = 0
        result.append(opts)
        if option.help:
            help_text = self.expand_default(option)
            # Do not wrap — emit the full help string as a single line
            result.append("%*s%s\n" % (indent_first, "", help_text))
        elif opts[-1] != "\n":
            result.append("\n")
        return "".join(result)


def build_option_parser():
    myparser = OptionParser(
        version="{} version {}".format('%prog', VERSION),
        formatter=NoWrapFormatter(),
        description=__import__('mutation_scatter_plot.calculate_codon_frequencies',
                               fromlist=['']).__doc__,
    )
    myparser.add_option("--reference-infile", action="store", type="string",
        dest="reference_infilename", default=None, metavar="FILE",
        help="FASTA formatted input file with reference padded sequence or not")
    myparser.add_option("--padded-reference", action="store_true",
        dest="padded_reference", default=False,
        help="By default we do NOT require the reference sequence to be padded with '-' characters to match the alignment delineating INSertions. If it is not padded [default case] then INSertion will not be reported but gaps parsed in the alignment will be skipped as long until 3 nucleotides are available for codon translation. Regardless of this --padded-reference setting, length of the reference sequence must match length of each alignment line.")
    myparser.add_option("--alignment-file", action="store", type="string",
        dest="alignment_infilename", default=None, metavar="FILE",
        help="ALIGNMENT file in FASTA format with - (minus) chars to adjust the alignment to the --reference-infile")
    myparser.add_option("--outfile-prefix", action="store", type="string",
        dest="outfileprefix", default=None, metavar="FILE",
        help="It assumes *.frequencies.fasta files. The prefix specified should end with .frequencies . The .tsv and .unchanged_codons.tsv will be appended to the prefix.")
    myparser.add_option("--left-reference-offset", action="store", type="int",
        dest="left_reference_offset", default=0,
        help="First nucleotide of the ORF region of the REFERENCE of interest to be sliced out from the input sequences. This requires 0-based numbering.")
    myparser.add_option("--right-reference-offset", action="store", type="int",
        dest="right_reference_offset", default=0,
        help="Last nucleotide of the last codon of the REFERENCE of interest to be sliced out from the input sequences. This requires 0-based numbering.")
    myparser.add_option("--aa_start", action="store", type="int",
        dest="aa_start", default=0,
        help="Adjust (padded) real position of the very first codon unless (1 for an initiator ATG). This value is added to the codon position reported in the output TSV file (the ATG position minus one). Use this if you cannot use --left-reference-offset nor --right-reference-offset which would have been used for slicing the input reference. The value provided is decremented by one to match pythonic 0-based numbering.")
    myparser.add_option("--min_start", action="store", type="int",
        dest="min_start", default=0,
        help="Start parsing the alignment since this position of the ALIGNMENT file. This requires 1-based numbering. This is to speedup parsing of input sequences and of the reference by skipping typical leading and trailing padding dashes. Default: 0 (parse since the beginning)")
    myparser.add_option("--max_stop", action="store", type="int",
        dest="max_stop", default=0,
        help="Stop parsing the alignment at this position of the ALIGNMENT file. This requires 1-based numbering. This is to speedup parsing of input sequences and of the reference by skipping typical leading and trailing padding dashes. Default: 0 (parse until the very end)")
    myparser.add_option("--x-after-count", action="store_true",
        dest="x_after_count", default=False,
        help="The FASTA file ID contains the count value followed by lowercase 'x'")
    myparser.add_option("--print-unchanged-sites", action="store_true",
        dest="print_unchanged_sites", default=False,
        help="Print out also sites with unchanged codons in to unchanged_codons.tsv file")
    myparser.add_option("--discard-this-many-leading-nucs", action="store", type="int",
        dest="discard_this_many_leading_nucs", default=0,
        help="Specify how many offending nucleotides are at the front of the FASTA sequences shifting the reading frame of the input FASTA file from frame +1 to either of the two remaining. Count the leading dashes and eventual nucleotides of incomplete codons too and check if it can be divided by 3.0 without slack. By default reading frame +1 is expected and hence no leading nucleotides are discarded. Default: 0")
    myparser.add_option("--discard-this-many-trailing-nucs", action="store", type="int",
        dest="discard_this_many_trailing_nucs", default=0,
        help="Specify how many offending nucleotides are at the end of each sequence. Default: 0")
    myparser.add_option("--minimum-alignments-length", action="store", type="int",
        dest="minimum_aln_length", default=50,
        help="Minimum length of aligned NGS read to be used for calculations")
    myparser.add_option("--debug", action="store", type="int",
        dest="debug", default=0,
        help="Set debug level to some real number")
    myparser.add_option("--overwrite", action="store_true",
        dest="overwrite", default=False,
        help="Overwrite existing output files instead of raising RuntimeError")
    return myparser


def main():
    myparser = build_option_parser()
    myoptions, myargs = myparser.parse_args()

    # parse the reference DNA
    if not myoptions.reference_infilename:
        raise ValueError("Error: Please specify --reference-infile with FASTA sequence")
    elif not os.path.exists(myoptions.reference_infilename):
        raise ValueError(f"Error: File {myoptions.reference_infilename} does not exist")
    elif os.path.getsize(myoptions.reference_infilename) == 0:
        raise ValueError(f"Error: File {myoptions.reference_infilename} is empty")
    else:
        for _record in SeqIO.parse(myoptions.reference_infilename, "fasta"):
            _padded_reference_dna_seq = str(_record.seq)
            break  # parse only the first and supposedly the only entry
        _reference_protein_seq = alt_translate(_padded_reference_dna_seq)
        _reference_as_codons = get_codons(_padded_reference_dna_seq,
                                          debug=myoptions.debug)

    if myoptions.alignment_infilename:
        _alnfilename_count_handle = open(
            f"{'.'.join(myoptions.alignment_infilename.split('.')[:-1])}.count", 'w'
        )
    else:
        raise RuntimeError("Please specify --alignment-file")

    if myoptions.outfileprefix:
        if myoptions.outfileprefix.endswith('.tsv'):
            _outfilename_handle = open_file(myoptions.outfileprefix, overwrite=myoptions.overwrite)
            _outfilename_unchanged_codons_handle = open_file(
                f"{myoptions.outfileprefix[:-4]}.unchanged_codons.tsv",
                overwrite=myoptions.overwrite
            )
        else:
            _outfilename_handle = open_file(f"{myoptions.outfileprefix}.tsv", overwrite=myoptions.overwrite)
            _outfilename_unchanged_codons_handle = open_file(
                f"{myoptions.outfileprefix}.unchanged_codons.tsv",
                overwrite=myoptions.overwrite
            )
    else:
        raise RuntimeError("Please specify output filename prefix via --outfile-prefix")

    _aa_start = (myoptions.aa_start - 1) if myoptions.aa_start else 0
    _min_start = (myoptions.min_start - 1) if myoptions.min_start else 0
    _max_stop  = (myoptions.max_stop  + 1) if myoptions.max_stop  else 0

    if myoptions.alignment_infilename and os.path.exists(myoptions.alignment_infilename):
        if os.path.getsize(myoptions.alignment_infilename) == 0:
            raise RuntimeError(f"Input file {myoptions.alignment_infilename} is empty")
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
        )
    else:
        raise RuntimeError(
            f"Input file {str(myoptions.alignment_infilename)} does not exist or is not defined"
        )


if __name__ == "__main__":
    main()
