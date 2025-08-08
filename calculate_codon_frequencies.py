#! /usr/bin/env python3

# This work © 2025 by Jiří Zahradník and Martin Mokrejš
# (First Medical Faculty - Charles University in Prague) is licensed under
# Creative Commons Attribution 4.0 International. To view a copy of this
# license, visit https://creativecommons.org/licenses/by/4.0/

"""This program calculate_codon_frequencies.py can parse multi-FASTA
files (FASTA) with padding dashes '-' and calculates frequencies of all up
to 64 codons (or up to 20 amino acids) in the respective column of the FASTA
file. It also counts number of INS-, DEL-, X- and STOP-events in the
sample multi-FASTA file (compared to pre-aligned reference).

It uses decimal module to keep numeric precision at maximum and only converts
to floating point numbers during output if not outputting just a string
representation of the number.

To prepare the input data we typically map NGS short sequences from Illumina
platform using 'bwa mem' to a reference protein-coding ORF sequence and then
convert the SAM file format to FASTA using `gofasta sam toMultiAlign` command.
The tool gofasta can be obtained from https://github.com/virus-evolution/gofasta .

We aligned the sequences using BLASTN to the de-padded reference sequence
and extracted padded alignment of the best match. This multi-FASTA file however
has varying length of the entries, caused by INSertion event in the sample
sequences compared tot he reference. To get a well aligned alignment with all
of them with same length one needs to extend the original reference sequence
with extra dashes to yield same length. Then this program can report INSertion
in the sample properly.

The utility discards internally codons containing unknown ('N') nucleotides.
It uses Biopython for parsing and translation. Similarly it discards incomplete
codon containing just one or two dashes.

The output is a TAB-separated TSV file with frequencies, e.g.:

146	N	D	0.137285	GGT	GAC	6750538
146	N	D	0.638309	GGT	GAT	6750538
146	N	X	0.001593	GGT	GNT	6750538
146	N	X	0.007141	GGT	GRT	6750538
146	N	V	0.000231	GGT	GTT	6750538
147	D	X	0.000809	GTT	GT-	5826744
147	D	X	0.000437	GTT	NN-	5826744
147	D	F	0.000144	GTT	TTT	5826744
147	D	DEL	0.137225	GTT	---	5826744
147	INS	X	0.446303	---	--N	9926
147	INS	X	0.297502	---	--T	9926
147	INS	T	0.256196	---	ACT	9926

If the reference sequence contains no padding symbols '-', it will not report
any INSertions in the sample but should work otherwise. To respect reading frame
it follows reference sequence and parses always 3 nucleotides at once from the FASTA
input file.
If there was a gap (padding) character '-' if picks an extra character from
the input until it has 3 nucleotides. Then it calculates frequencies for
all possible codons. The results are stored in a TAB-separated TSV file for
easy post-processing by mutation_scatter_plot.py which on-the-fly discards
codons containing unknown (N) nucleotides using python Pandas library and
finally draws interactive figures using Matplotlib and Bokeh graphical
libraries.


Please cite the following article if you use our data or software in your research:

Shoshany A., Tian R., Padilla-Blanco M., Hruška A., Baxova K., Zoler E., Mokrejš M., Schreiber G., Zahradník J. (submitted) In Vitro and Viral Evolution Convergence Reveal the Selective Pressures Driving Omicron Emergence. [bioRxiv](https://www.biorxiv.org/content/10.1101/2025.04.23.650148v1)
"""

import os
import re
import sys
from subprocess import Popen, PIPE, call

from collections import Counter
import hashlib

from optparse import OptionParser
from decimal import Decimal
from numpy import median, average

from Bio import SeqIO
from Bio.Seq import Seq, reverse_complement, translate
from Bio import AlignIO

version = 202508082050

myparser = OptionParser()
myparser.add_option("--reference-infile", action="store", type="string", dest="reference_infilename", default=None, metavar="FILE",
    help="FASTA formatted input file with reference padded sequence or not")
myparser.add_option("--padded-reference", action="store_true", dest="padded_reference", default=False,
    help="By default we do NOT require the reference sequence to be padded with '-' characters to match the alignment delineating INSertions. If it is not padded [default case] then INSertion will not be reported but gaps parsed in the alignment will be skipped as long until 3 nucleotides are available for codon translation. Regardless of this --padded-reference setting, length of the reference sequence must match length of each alignment line.")
myparser.add_option("--alignment-file", action="store", type="string", dest="alignment_infilename", default=None, metavar="FILE",
    help="Alignment file in FASTA format with - (minus) chars to adjust the alignment to the --reference-infile")
myparser.add_option("--outfile-prefix", action="store", type="string", dest="outfileprefix", default=None, metavar="FILE",
    help="It assumes *.frequencies.fasta files. The prefix specified should end with .frequencies . The .tsv and .unchanged_codons.tsv will be appended to the prefix.")
myparser.add_option("--left-reference-offset", action="store", type="int", dest="left_reference_offset", default=0,
    help="First nucleotide of the ORF region of the REFERENCE of interest to be sliced out from the input sequences. This requires 0-based numbering.")
myparser.add_option("--right-reference-offset", action="store", type="int", dest="right_reference_offset", default=0,
    help="Last nucleotide of the last codon of the REFERENCE of interest to be sliced out from the input sequences. This requires 0-based numbering.")
myparser.add_option("--aa_start", action="store", type="int", dest="aa_start", default=0,
    help="Real position of the very first codon unless (1 for an initiator ATG). This value is added to the codon position reported in the output TSV file (the ATG position minus one). Use this if you cannot use --left-reference-offset nor --right-reference-offset which would have been used for slicing the input reference. The value provided is decremented by one to match pythonic 0-based numbering.")
myparser.add_option("--min_start", action="store", type="int", dest="min_start", default=0,
    help="Start parsing the alignment since this position of the amplicon region. This requires 1-based numbering. This is to speedup parsing of input sequences and of the reference by skipping typical leading and trailing padding dashes. Default: 0 (parse since the beginning)")
myparser.add_option("--max_stop", action="store", type="int", dest="max_stop", default=0,
    help="Stop parsing the alignment at this position of the amplicon region. This requires 1-based numbering. This is to speedup parsing of input sequences and of the reference by skipping typical leading and trailing padding dashes. Default: 0 (parse until the very end)")
myparser.add_option("--x-after-count", action="store_true", dest="x_after_count", default=False,
    help="The FASTA file ID contains the count value followed by lowercase 'x'")
myparser.add_option("--print-unchanged-sites", action="store_true", dest="print_unchanged_sites", default=False,
    help="Print out also sites with unchanged codons in to unchanged_codons.tsv file")
myparser.add_option("--discard-this-many-leading-nucs", action="store", type="int", dest="discard_this_many_leading_nucs", default=0,
    help="Specify how many offending nucleotides are at the front of the FASTA sequences shifting the reading frame of the input FASTA file from frame +1 to either of the two remaining. Count the leading dashes and eventual nucleotides of incomplete codons too and check if it can be divided by 3.0 without slack. By default reading frame +1 is expected and hence no leading nucleotides are discarded. Default: 0")
myparser.add_option("--discard-this-many-trailing-nucs", action="store", type="int", dest="discard_this_many_trailing_nucs", default=0,
    help="Specify how many offending nucleotides are at the end of each sequence. Default: 0")
myparser.add_option("--minimum-alignments-length", action="store", type="int", dest="minimum_aln_length", default=50,
    help="Minimum length of aligned NGS read to be used for calculations")
myparser.add_option("--debug", action="store", type="int", dest="debug", default=0,
    help="Set debug level to some real number")

(myoptions, myargs) = myparser.parse_args()


def get_codons(seq):
    """This function is only used to parse reference sequence and if it was not
    already in frame +1 and divisible by 3 without remainder, any padding dashes
    are removed and sequence is split into triplets (codons)."""

    if len(seq) % 3 == 0:
        #split list into codons of three
        codons = [seq[i:i+3] for i in range(0, len(seq), 3)]
    elif len(seq.replace('-','')) % 3 == 0:
        #split list into codons of three
        codons = [seq[i:i+3] for i in range(0, len(seq.replace('-','')), 3)]
        if myoptions.debug: print("Debug: Detected %s minus signs in the sequence but after all the nucleotide sequence can be divided by three when they are omitted, good." % seq.count('-'))
    else:
        raise ValueError("Error: Sequence %s cannot be divided by 3 and removing minus signs does not help either" % seq)
    return codons

def alt_translate(seq):
    """Biopython cannot sometimes translate a sequence but one can get around by
    splitting it into codons and merging later back
    https://github.com/biopython/biopython/pull/4992
    https://github.com/biopython/biopython/pull/4992#issuecomment-2865429105
    """

    codons = (seq[i:i+3] for i in range(0, len(seq), 3))
    codons = ("NNN" if "-" in codon and codon != "---" else codon for codon in codons)
    return "".join(translate(codon, gap='-') for codon in codons)

def write_tsv_line(outfilename, codons, natural_codon_position_padded, natural_codon_position_depadded, reference_aa, total_codons_per_site_sum, reference_codon, debug=False):
    if not total_codons_per_site_sum:
        _total_codons_per_site_sum = 0
    else:
        _total_codons_per_site_sum = total_codons_per_site_sum
    for _some_codon in codons:
        if len(_some_codon) < 3:
            _some_aa = 'X'
        else:
            _some_aa = alt_translate(_some_codon)
        _observed_codon_count = Decimal(codons[_some_codon])
        if not _observed_codon_count:
            _observed_codon_count2 = 0
        else:
            _observed_codon_count2 = _observed_codon_count
        outfilename.write("{}\t{}\t{}\t{}\t{:8.6f}\t{}\t{}\t{}\t{}\n".format(natural_codon_position_padded, natural_codon_position_depadded, reference_aa, _some_aa, _observed_codon_count2 / Decimal(total_codons_per_site_sum), reference_codon, _some_codon, _observed_codon_count2, _total_codons_per_site_sum))
        if debug: print("TESTING1:\t{}\t{}\t{}\t{}\t{:8.6f}\t{}\t{}\t{}\t{}".format(natural_codon_position_padded, natural_codon_position_depadded, reference_aa, _some_aa, _observed_codon_count2 / Decimal(total_codons_per_site_sum), reference_codon, _some_codon, _observed_codon_count2, _total_codons_per_site_sum))
    outfilename.flush()


def parse_alignment(alignment_file, padded_reference_dna_seq, reference_protein_seq, reference_as_codons, outfilename, outfilename_unchanged_codons, alnfilename_count, aa_start, min_start, max_stop, calculate_checksums=False):
    """left_reference_offset and right_reference_offset are used to slice the reference
    discard_this_many_leading_nucs and discard_this_many_trailing_nucs are used to discard some leading/trailing nucs if the aligned region was a bit wider, aspecially if not starting at frame +1
    """

    if myoptions.debug: print("Debug0: Depadded reference sequence has length %d, padded reference sequence has length %d, each entry from %s must also have same padded length %d" % (len(padded_reference_dna_seq.replace('-','')), len(padded_reference_dna_seq), alignment_file, len(padded_reference_dna_seq)))
    if myoptions.debug: print("Debug1: reference_protein_seq=%s with length %d" % (str(reference_protein_seq), len(reference_protein_seq)))
    try:
        _align = AlignIO.read(alignment_file, "fasta") # read whole alignment into memory
    except ValueError:
        raise ValueError("Error: one of the entries in the %s file has different length" % alignment_file)
    if myoptions.left_reference_offset or myoptions.right_reference_offset:
        _padded_reference_dna_seq = padded_reference_dna_seq[max(myoptions.left_reference_offset - 1, 0) : min(len(padded_reference_dna_seq), myoptions.right_reference_offset)] # discard leading sequence of the reference not contained in FASTA file, if the left_reference_offset is correct
        _reference_protein_seq = reference_protein_seq[int(max(myoptions.left_reference_offset - 1, 0) / 3): int(min(len(reference_protein_seq), myoptions.right_reference_offset / 3))]
        _reference_as_codons = reference_as_codons[int(max(myoptions.left_reference_offset - 1, 0) / 3): int(min(len(reference_as_codons), myoptions.right_reference_offset / 3))]
        if myoptions.debug:
            print("Info: len(padded_reference_dna_seq)=%s, myoptions.left_reference_offset=%s, myoptions.right_reference_offset=%s" % (len(padded_reference_dna_seq), myoptions.left_reference_offset, myoptions.right_reference_offset))
            print("Info: After cutting input using offset position: %s" % _padded_reference_dna_seq)
            print("Info: After cutting input using offset position: %s" % _reference_protein_seq)
            print("Info: After cutting input using offset position: %s" % _reference_as_codons)
        if not _reference_protein_seq:
            raise ValueError("Error: No _reference_protein_seq provided or left. Was the slicing using myoptions.left_reference_offset=%s, myoptions.right_reference_offset=%s wrong?" % (myoptions.left_reference_offset, myoptions.right_reference_offset))
        if not _reference_as_codons:
            raise ValueError("Error: No _reference_as_codons provided or left. Was the slicing using myoptions.left_reference_offset=%s, myoptions.right_reference_offset=%s wrong?" % (myoptions.left_reference_offset, myoptions.right_reference_offset))

        if myoptions.debug: print("Debug2: Depadded reference sequence has length %d, padded reference sequence has length %d, each padded entry from %s must also have same padded length" % (len(_padded_reference_dna_seq.replace('-','')), len(_padded_reference_dna_seq), alignment_file))
    else:
        _padded_reference_dna_seq = padded_reference_dna_seq
        _reference_protein_seq = reference_protein_seq
        _reference_as_codons = reference_as_codons

    _zero_based_padded_reference_codon_index = 0 # used to fetch next codon or a dash
    _zero_based_padded_reference_aa_index = 0 # used to fetch next aa residue or a dash from already sliced-out reference sequence, will be increased by one per every 3 nt
    _reference_aa = _reference_protein_seq[_zero_based_padded_reference_aa_index] # fetch the very first aa residue, supposedly 'M' but above we already chopped the sequence and cut only a region of interest
    _previous_gaps = 0
    _new_gaps_in_reference = 0
    _re_leading_gaps = re.compile("^[-Nn]+")
    _re_trailing_gaps = re.compile("[-Nn]+$")
    _already_checked_starts = [] # keep list of the pythonic codon _starts of the for loop below so that some check are performed only once per iteration, not for every input alignment line
    _top_most_codons = []
    _total_aln_entries_used = 0
    _start_from = myoptions.discard_this_many_leading_nucs # cut away some eventually offending leading nucleotides to get the alignment into frame +1
    _stop_to = myoptions.discard_this_many_trailing_nucs # cut away some eventually offending trailing nucleotides
    for _zero_based_codon_startpos in range(min_start, max_stop or len(_padded_reference_dna_seq), 3): # generate pythonic numbers so the pythonic slicing below works, the input sequence must be in frame +1 (if not, use discard_this_many_leading_nucs and discard_this_many_trailing_nucs)
        _unchanged_codons = Counter()
        _changed_codons = Counter()
        _deleted_reference_codons = Counter()
        _inserted_codons = Counter()
        _unchanged_aa_residues = Counter()
        _changed_aa_residues = Counter()
        _deleted_reference_aa_residues = Counter()
        _new_gaps_in_reference = 0
        _is_deletion = False
        _is_insertion = False
        _new_aa_residue = None
        _reference_aa = _reference_protein_seq[_zero_based_padded_reference_aa_index].upper() # fetch the very first aa residue, supposedly 'M' but above we already chopped the sequence and cut only a region of interest
        _current_codon_position = _zero_based_codon_startpos / 3 + 1
        if myoptions.debug: print("Debug3: _start=%s, _padded_reference_dna_seq=%s" % (_zero_based_codon_startpos, _padded_reference_dna_seq))
        for _aln_line in _align:
            _record_id = _aln_line.id
            if myoptions.x_after_count:
                if 'x.' in _record_id:
                    _record_count, _checksum = _record_id.split('x.')
                    _record_count = int(_record_count)
                else:
                    _record_count = int(_record_id.replace('x',''))
                    if calculate_checksums:
                        _hash_func = hashlib.sha256() # always create a new object other .update() just appends to it new string
                        _hash_func.update(_record_id.seq.encode()) # calculate the hash without a trailing newline
                        _checksum = _hash_func.hexdigest() # TODO: _checksum is unusued
            else:
                _record_count = 1
                if calculate_checksums:
                    _hash_func = hashlib.sha256() # always create a new object other .update() just appends to it new string
                    _hash_func.update(_record_id.seq.encode()) # calculate the hash without a trailing newline
                    _checksum = _hash_func.hexdigest() # TODO: _checksum is unusued
            if not _zero_based_codon_startpos or not _total_aln_entries_used:
                # count lines and actually the counts mentioned in the FASTA ID only when iterating over the first codon
                _total_aln_entries_used += _record_count
            if _start_from:
                if _stop_to:
                    _aln_line_seq = str(_aln_line.seq)[_start_from:-_stop_to] # discard leading one or two nucleotides to get into frame +1 while reference was supposedly in +1 frame already
                else:
                    _aln_line_seq = str(_aln_line.seq)[_start_from:] # discard leading one or two nucleotides to get into frame +1 while reference was supposedly in +1 frame already
            elif _stop_to:
                _aln_line_seq = str(_aln_line.seq)[:-_stop_to] # discard some trailing nucleotides
            else:
                _aln_line_seq = str(_aln_line.seq)
            _amplicon_length = (max_stop or len(_padded_reference_dna_seq)) - min_start
            _padded_aln_line_length = len(_aln_line_seq) # padded length
            _depadded_aln_line_length = len(_aln_line_seq.replace('-','').lstrip('Nn').rstrip('Nn'))
            _new_gaps_in_reference = 0 # zap value from previous codon search if it contained a minus char
            _codon = None
            _deleted_reference_codon = None
            _deleted_reference_codon_depadded = None
            _reference_codon = _padded_reference_dna_seq[3*_zero_based_padded_reference_aa_index:3*_zero_based_padded_reference_aa_index + 3].upper()
            _reference_codon_depadded = _reference_codon.replace('-','')
            _rough_sample_codon = _aln_line_seq[_zero_based_codon_startpos:_zero_based_codon_startpos + 3].upper()
            _sample_codon_depadded = _rough_sample_codon.replace('-','')
            if len(_rough_sample_codon) != len(_sample_codon_depadded):
                _sample_codon_contained_pad = True
            else:
                _sample_codon_contained_pad = False
            if len(_reference_codon) != len(_reference_codon_depadded):
                _reference_codon_contained_pad = True
            else:
                _reference_codon_contained_pad = False
            _new_aa_residue = None

            _start_of_trailing_gaps = 0
            _end_of_leading_gaps = 0
            for _match in _re_leading_gaps.finditer(_aln_line_seq):
                _end_of_leading_gaps = _match.end()
            for _match in _re_trailing_gaps.finditer(_aln_line_seq):
                _start_of_trailing_gaps = _match.start()
            if myoptions.debug: print(f"Debug4: End of leading gaps is at {_end_of_leading_gaps:6d}, start of trailing gaps is at {_start_of_trailing_gaps:6d}, {_aln_line_seq}")

            if myoptions.minimum_aln_length and _depadded_aln_line_length < myoptions.minimum_aln_length:
                pass # skip alignments shorter than 50nt
                if myoptions.debug: print("Debug5:Here1, _end_of_leading_gaps=%s, _start_of_trailing_gaps=%s, _zero_based_codon_startpos=%s" % (_end_of_leading_gaps, _start_of_trailing_gaps, _zero_based_codon_startpos))
            elif not _depadded_aln_line_length:
                pass # skip empty alignments
                if myoptions.debug: print("Debug6:Here2, _end_of_leading_gaps=%s, _start_of_trailing_gaps=%s, _zero_based_codon_startpos=%s" % (_end_of_leading_gaps, _start_of_trailing_gaps, _zero_based_codon_startpos))
            elif _end_of_leading_gaps and _zero_based_codon_startpos < _end_of_leading_gaps - 1:
                pass # skip leading gaps so they do not get counted as DELetions
                if myoptions.debug: print("Debug7:Here3, _end_of_leading_gaps=%s, _start_of_trailing_gaps=%s, _zero_based_codon_startpos=%s" % (_end_of_leading_gaps, _start_of_trailing_gaps, _zero_based_codon_startpos))
            elif _start_of_trailing_gaps and _zero_based_codon_startpos + 1 > _start_of_trailing_gaps:
                pass # skip trailing gaps so they do not get counted as DELetions
                if myoptions.debug: print("Debug8:Here4, _end_of_leading_gaps=%s, _start_of_trailing_gaps=%s, _zero_based_codon_startpos=%s" % (_end_of_leading_gaps, _start_of_trailing_gaps, _zero_based_codon_startpos))
            elif not _end_of_leading_gaps and _zero_based_codon_startpos + 1 > _start_of_trailing_gaps and _start_of_trailing_gaps:
                pass # skip trailing gaps so they do not get counted as DELetions but do not enter here is there are no leading or trailing dashes at all causing _start_of_trailing_gaps=0 _end_of_leading_gaps=0
                # Q409del
                if myoptions.debug: print("Debug9:Here5, _end_of_leading_gaps=%s, _start_of_trailing_gaps=%s, _zero_based_codon_startpos=%s" % (_end_of_leading_gaps, _start_of_trailing_gaps, _zero_based_codon_startpos))
            elif not _start_of_trailing_gaps and _zero_based_codon_startpos < _end_of_leading_gaps + 1 and _end_of_leading_gaps:
                if myoptions.debug: print("Debug10:Here6, _end_of_leading_gaps=%s, _start_of_trailing_gaps=%s, _zero_based_codon_startpos=%s" % (_end_of_leading_gaps, _start_of_trailing_gaps, _zero_based_codon_startpos))
                pass # skip leading gaps so they do not get counted as DELetions but do not enter here when there are no leading or trailing dashes at all causing _start_of_trailing_gaps=0 _end_of_leading_gaps=0
            elif _padded_aln_line_length + 1 > _zero_based_codon_startpos + 3: # avoid infinite loop when increasing slice end position behind the existing input sequence
                if myoptions.debug: print("Debug11: Rough sample codon at %d is %s, reference codon is %s, amplicon region length %s (nt), sliced [%d:%d] (pythonic slice numbering), _sample_codon_depadded is %s, _reference_codon_depadded is %s, _sample_codon_contained_pad=%s, _reference_codon_contained_pad=%s" % (_current_codon_position, _rough_sample_codon, _reference_codon, _amplicon_length, _zero_based_codon_startpos+_previous_gaps, _zero_based_codon_startpos+_previous_gaps+3+_new_gaps_in_reference, _sample_codon_depadded, _reference_codon_depadded, str(_sample_codon_contained_pad), str(_reference_codon_contained_pad)))
                if _rough_sample_codon == 'NNN' or _reference_codon == 'NNN':
                    if myoptions.debug: print("Debug11: Skipping a codon containing NNN")
                if _rough_sample_codon == '---' and _reference_codon != '---' and _reference_codon != 'NNN':
                    _is_deletion = True
                    _deleted_reference_codon = str(_reference_codon) # make an extra variable as we re-write it in the below loop
                    _deleted_reference_codon_depadded = _deleted_reference_codon.replace('-','')
                    if myoptions.debug > 1: print("Debug12: Padded deleted codon at %d is %s, original reference codon is %s, amplicon region length %s (nt), end of slice is %d" % (_current_codon_position, _deleted_reference_codon, _reference_codon, _amplicon_length, _zero_based_codon_startpos + _previous_gaps + 3 + _new_gaps_in_reference))
                    _deleted_reference_codons[_deleted_reference_codon] += _record_count
                    _deleted_aa_residue = alt_translate(_deleted_reference_codon)
                    _deleted_reference_aa_residues[_deleted_aa_residue] += _record_count
                    _codon = '---'
                    _sample_codon_depadded = ''
                    if myoptions.debug: print("Debug13: Added Deleted codon at %d is %s, reference codon is %s, amplicon region length %s (nt), sliced [%d:%d] (pythonic slice numbering), aa_residue is %s" % (_current_codon_position, '---', _reference_codon, _amplicon_length, _zero_based_codon_startpos+_previous_gaps, _zero_based_codon_startpos + _previous_gaps + 3 + _new_gaps_in_reference, _deleted_aa_residue))
                elif _reference_codon == '---' and _rough_sample_codon != '---':
                    _is_insertion = True
                    _sample_codon_depadded = _rough_sample_codon.replace('-','')
                    _reference_aa = 'INS'
                    _inserted_codons[_sample_codon_depadded] += _record_count # use depadded representation to compress --N and -N- and N-- together
                    _new_aa_residue = alt_translate(_rough_sample_codon)
                    if _new_aa_residue:
                        _changed_aa_residues[_new_aa_residue] += _record_count
                    else:
                        # _new_aa_residue = '-'
                        pass # do not count junk if reference was anyway just gaps
                    if myoptions.debug > 1: print("Debug14: Added Deleted codon at %d is %s, reference codon is %s, amplicon region length %s (nt), sliced [%d:%d] (pythonic slice numbering), aa_residue is %s" % (_current_codon_position, '---', _reference_codon, _amplicon_length, _zero_based_codon_startpos+_previous_gaps, _zero_based_codon_startpos + _previous_gaps + 3 + _new_gaps_in_reference, _new_aa_residue))
                else:
                    _sample_codon_depadded = _rough_sample_codon.replace('-','')
                    if myoptions.debug > 1: print("Debug15: Rough sample codon at %d is %s, reference codon is %s, amplicon region length %s (nt), sliced [%d:%d] (pythonic slice numbering), _sample_codon_depadded is %s, _reference_codon_depadded is %s, _sample_codon_contained_pad=%s, _reference_codon_contained_pad=%s" % (_current_codon_position, _rough_sample_codon, _reference_codon, _amplicon_length, _zero_based_codon_startpos+_previous_gaps, _zero_based_codon_startpos+_previous_gaps+3+_new_gaps_in_reference, _sample_codon_depadded, _reference_codon_depadded, str(_sample_codon_contained_pad), str(_reference_codon_contained_pad)))

                    if (_sample_codon_contained_pad and _reference_codon_contained_pad) and not _sample_codon_depadded and not _reference_codon:
                        # rough codon is ---, reference codon is ---
                        _unchanged_codons['---'] += _record_count
                        _unchanged_aa_residues['-'] += _record_count
                        if myoptions.debug > 1: print("Debug16: Skipping unchanged codon at %d due to a padding dash in both %s and %s" % (_current_codon_position, _rough_sample_codon, _reference_codon))
                    elif _rough_sample_codon == '---' and _reference_codon == '---':
                        _unchanged_codons['---'] += _record_count
                        _unchanged_aa_residues['-'] += _record_count
                        if myoptions.debug > 1: print("Debug17: Skipping unchanged codon at %d due to a padding dash in both sample and reference codon %s" % (_current_codon_position, _rough_sample_codon))
                    elif (_sample_codon_contained_pad or _reference_codon_contained_pad) and _sample_codon_depadded == _reference_codon_depadded:
                        # skip empty codon parsed just after the end of the sequence
                        # skip 'AT-C' versus 'ATC' but increase the incidence counter
                        _unchanged_codons[_sample_codon_depadded] += _record_count
                        _unchanged_aa_residues[alt_translate(_sample_codon_depadded)] += _record_count
                        if myoptions.debug > 1: print("Debug18: Skipping unchanged codon at %d due to a padding dash in both %s and %s" % (_current_codon_position, _rough_sample_codon, _reference_codon))
                    elif _sample_codon_depadded == _reference_codon_depadded:
                        _unchanged_codons[_sample_codon_depadded] += _record_count
                        _unchanged_aa_residues[alt_translate(_sample_codon_depadded)] += _record_count
                        if myoptions.debug > 1: print("Debug19: Unchanged codon at %d is %s" % (_current_codon_position, _sample_codon_depadded))
                    elif _sample_codon_depadded != _reference_codon_depadded:
                        _changed_codons[_sample_codon_depadded] += _record_count
                        _new_aa_residue = alt_translate(_rough_sample_codon)
                        if 'N' in _sample_codon_depadded and _new_aa_residue != 'X':
                            if _sample_codon_depadded not in ('TCN', 'CTN', 'GTN', 'CCN', 'ACN', 'GCN', 'CGN', 'GGN'):
                                raise ValueError("Error: biopython translated codon %s at %d into %s" % (_sample_codon_depadded, _current_codon_position, _new_aa_residue))
                            else:
                                if myoptions.debug > 1: print("Debug20: biopython translated codon %s at %d into %s" % (_sample_codon_depadded, _current_codon_position, _new_aa_residue))
                        if _new_aa_residue:
                            if _new_aa_residue != _reference_aa:
                                _changed_aa_residues[_new_aa_residue] += _record_count
                            else:
                                _unchanged_aa_residues[_new_aa_residue] += _record_count
                        else:
                            _changed_aa_residues['None'] += _record_count
                            if myoptions.debug > 1: print("Debug21: Added None to _changed_aa_residues for codon at %s" % (_current_codon_position))
                        if myoptions.debug: print("Debug22: Final codon at %d is %s, reference codon is %s, amplicon region length %s (nt), sliced [%d:%d] (pythonic slice numbering), aa_residue is %s" % (_current_codon_position, _rough_sample_codon, _reference_codon, _amplicon_length, _zero_based_codon_startpos+_previous_gaps, _zero_based_codon_startpos+_previous_gaps+3+_new_gaps_in_reference, _new_aa_residue))
                    else:
                        _unchanged_codons[_rough_sample_codon] += _record_count
                        _new_aa_residue = alt_translate(_rough_sample_codon)
                        if _new_aa_residue:
                            if _new_aa_residue != _reference_aa:
                                _changed_aa_residues[_new_aa_residue] += _record_count
                            else:
                                _unchanged_aa_residues[_new_aa_residue] += _record_count
                        else:
                            _changed_aa_residues['None'] += _record_count
                        if myoptions.debug > 1: print("Debug23: Including frame-breaking codon %s at %d into per-codon coverage, reference codon is %s, amplicon region length %s (nt), sliced [%d:%d] (pythonic slice numbering), aa_residue is %s" % (_rough_sample_codon, _current_codon_position, _reference_codon, _amplicon_length, _zero_based_codon_startpos+_previous_gaps, _zero_based_codon_startpos+_previous_gaps+3+_new_gaps_in_reference, _new_aa_residue))
            else:
                if myoptions.debug: print("Debug10b:Here7. Leaking sample codon %s at position %d, _amplicon_length=%d" % (_rough_sample_codon, _current_codon_position, _amplicon_length))

        # all entries in padded FASTA file were processed for the current codon column
        if myoptions.debug: print("Debug24: Lists at %d of changed codons and aminoacid residues: %s-%s %s %s" % (_current_codon_position, _zero_based_padded_reference_aa_index + _zero_based_codon_startpos + 1, _zero_based_padded_reference_aa_index + _zero_based_codon_startpos + 3, _changed_codons, _changed_aa_residues)) # undo off-by-one error due to pythonic counting
        if myoptions.debug: print("Debug25: Lists at %d of unchanged codons and aminoacid residues: %s-%s %s %s" % (_current_codon_position, _zero_based_padded_reference_aa_index + _zero_based_codon_startpos + 1, _zero_based_padded_reference_aa_index + _zero_based_codon_startpos + 3, _unchanged_codons, _unchanged_aa_residues)) # undo off-by-one error due to pythonic counting

        if myoptions.debug: print("Debug26a: Deleted codons at %d: %s" % (_current_codon_position, _deleted_reference_codons))
        if myoptions.debug: print("Debug26b: Inserted codons at %d: %s" % (_current_codon_position, _inserted_codons))

        _total_codons_per_site_counts = _inserted_codons + _deleted_reference_codons + _changed_codons + _unchanged_codons
        _total_codons_per_site_sum = sum(_total_codons_per_site_counts.values())
        # _coverage_per_codon = sum(Counter(_inserted_codons + _deleted_reference_codons + _changed_codons + _unchanged_codons).values())

        if myoptions.debug: print("Debug26c: _reference_codon was %s" % str(_reference_codon))
        #_reference_codon = _padded_reference_dna_seq[_zero_based_codon_startpos+int(_previous_gaps/3.0):_zero_based_codon_startpos + int(_previous_gaps/3.0) + 3 + _new_gaps_in_reference]
        #if myoptions.debug: print("Debug26d: _reference_codon changed to %s" % str(_reference_codon))

        if len(_reference_protein_seq) > _zero_based_padded_reference_aa_index:
            _reference_aa = _reference_protein_seq[_zero_based_padded_reference_aa_index]
            # otherwise keep previous values (actually keep the C-term residue)
        else:
            sys.stderr.write("Error: Probably _zero_based_padded_reference_aa_index=%s at the end of the padded reference sequence which was %s  long\n" % (_zero_based_padded_reference_aa_index, len(_reference_protein_seq)))
            sys.stderr.flush()

        # this would performed for every input alignment line so it slows down the execution, hence perform the check only once per reference codon
        if _zero_based_codon_startpos not in _already_checked_starts and len(_reference_codon) == 3:
            if _reference_codon == '---':
                _reference_aa = '-'
                _current_aa = '-'
            else:
                try:
                    _current_aa = alt_translate(_reference_codon)
                except:
                    # Bio.Data.CodonTable.TranslationError: Codon '---' is invalid
                    _current_aa = 'Biopython failure'
            if _current_aa != _reference_aa:
                if not min_start:
                    # ValueError: The reference aa should be 'INS' but parsed reference codon '---' encodes 'None'. The _zero_based_codon_startpos = 45
                    raise ValueError("Error: The reference aa at %d should be '%s' but parsed reference codon '%s' encodes '%s'. The _zero_based_codon_startpos = %s" % (_current_codon_position, _reference_aa, _reference_codon, _current_aa, _zero_based_codon_startpos))
                    # but do not die if the alignment does not start at initiator methionine
            else:
                _already_checked_starts.append(_zero_based_codon_startpos)

        if myoptions.debug:
            if myoptions.debug > 1:
                for _key in _changed_codons:
                    print("Debug27: {}: {} = {:8.6f} {} {}".format(_key, _changed_codons[_key], _changed_codons[_key] / _total_codons_per_site_sum, set([x for x in _changed_codons]), set([x for x in _changed_aa_residues])))
        #if sum([ x for x in _changed_codons]) != sum([x for x in _changed_aa_residues]):
        #    raise ValueError("Error: len(_changed_codons)=%d != len(_changed_aa_residues)=%d, _changed_codons=%s, _changed_aa_residues=%s" % (len(_changed_codons), len(_changed_aa_residues), str(_changed_codons), str(_changed_aa_residues)))

        if myoptions.debug: print("Debug27a: _new_gaps_in_reference=%s, _reference_codon_contained_pad=%s, _sample_codon_contained_pad=%s, _reference_codon=%s, _rough_sample_codon=%s" %(_new_gaps_in_reference, _reference_codon_contained_pad, _sample_codon_contained_pad, _reference_codon, _rough_sample_codon))
        # count the pads now
        if _new_gaps_in_reference:
            raise ValueError("Error: _new_gaps_in_reference should be zero but is %s instead" % _new_gaps_in_reference)
        _new_gaps_in_reference = _reference_codon.count('-')

        _natural_codon_position_padded = _zero_based_padded_reference_aa_index + 1 + int( myoptions.left_reference_offset / 3.0) + aa_start
        _natural_codon_position_depadded = _natural_codon_position_padded - int((_previous_gaps + _new_gaps_in_reference) / 3.0)

        write_tsv_line(outfilename_unchanged_codons, _unchanged_codons, _natural_codon_position_padded, _natural_codon_position_depadded, _reference_aa, _total_codons_per_site_sum, _reference_codon, debug=myoptions.debug)

        # TODO: print contents of _unchanged_aa_residues eventually too
        if len(_changed_codons):
            write_tsv_line(outfilename, _changed_codons, _natural_codon_position_padded, _natural_codon_position_depadded, _reference_aa, _total_codons_per_site_sum, _reference_codon, debug=myoptions.debug)

        if len(_inserted_codons):
            write_tsv_line(outfilename, _inserted_codons, _natural_codon_position_padded, _natural_codon_position_depadded, 'INS', _total_codons_per_site_sum, _reference_codon, debug=myoptions.debug)

        if _is_deletion:
            for _some_deleted_codon in _deleted_reference_codons:
                _observed_codon_count = Decimal(_deleted_reference_codons[_some_deleted_codon])
                if not _observed_codon_count:
                    _observed_codon_count2 = 0
                else:
                    _observed_codon_count2 = _observed_codon_count
                outfilename.write("{}\t{}\t{}\t{}\t{:8.6f}\t{}\t{}\t{}\t{}\n".format(_natural_codon_position_padded, _natural_codon_position_depadded, _reference_aa, 'DEL', Decimal(_deleted_reference_codons[_some_deleted_codon]) / Decimal(_total_codons_per_site_sum), _some_deleted_codon, '---', _observed_codon_count2, _total_codons_per_site_sum))
                if myoptions.debug: print("TESTING2:\t{}\t{}\t{}\t{}\t{:8.6f}\t{}\t{}\t{}\t{}".format(_natural_codon_position_padded, _natural_codon_position_depadded, _reference_aa, 'DEL', Decimal(_deleted_reference_codons[_some_deleted_codon]) / Decimal(_total_codons_per_site_sum), _some_deleted_codon, '---', _observed_codon_count2, _total_codons_per_site_sum))
        outfilename.flush()

        _previous_gaps += _new_gaps_in_reference
        _new_gaps_in_reference = 0

        # get the top-most codon in this FASTA column
        try:
            _top_most_codon, _top_most_count = _total_codons_per_site_counts.most_common()[0]
        except IndexError:
            pass
        else:
            _top_most_codons.append(_top_most_codon)

#        if len(_reference_protein_seq) > _zero_based_padded_reference_aa_index + 1: # undo off-by-one error while stop increasing after the C-terminus of the reference protein
#            if _is_insertion:
#                _zero_based_padded_reference_aa_index += 1 # increase the value of the counter for every INS traversed too
#                if myoptions.debug: print("Debug33: Increased _zero_based_padded_reference_aa_index to %d, moving to next codon" % _zero_based_padded_reference_aa_index)
#            elif _is_deletion:
#                _zero_based_padded_reference_aa_index += 1 # increase the counter in extra for every DEL traversed
#                if myoptions.debug: print("Debug31: Increased _zero_based_padded_reference_aa_index to %d, moving to next codon" % _zero_based_padded_reference_aa_index)
#            else:
#                _zero_based_padded_reference_aa_index += 1 # increase the counter anyway since we support INS in the reference
#                if myoptions.debug: print("Debug32: Increased _zero_based_padded_reference_aa_index to %d, moving to next codon" % _zero_based_padded_reference_aa_index)
        if myoptions.debug: print("Debug33: After increments: _natural_codon_position_depadded=%d _previous_gaps=%d _new_gaps_in_reference=%d" % (_natural_codon_position_depadded, _previous_gaps, _new_gaps_in_reference))
        _zero_based_padded_reference_aa_index += 1 # increase the padded aa counter
        _zero_based_padded_reference_codon_index += 3

    del _align
    alnfilename_count.write("%s\n" % _total_aln_entries_used)
    alnfilename_count.close()
    _consensus = ''.join(_top_most_codons)
    print("Info: consensus = %s" % str(_consensus))
    if _consensus in _padded_reference_dna_seq:
        print("Info: Sample consensus sequence should roughly match substring inside the reference and IT DOES: %s" % str(_top_most_codons))
    else:
        print("Info: Sample consensus sequence should roughly match substring inside the reference BUT IT DOES NOT, maybe due to some true major mutations in some codons: %s" % str(_top_most_codons))


def open_file(outfilename):
    if os.path.exists(outfilename):
        raise RuntimeError("The file %s already exists, will not overwrite it." % outfilename)
    else:
        _outfilename_handle = open(outfilename, 'w')
        return _outfilename_handle

def main():
    # parse the reference DNA
    if not myoptions.reference_infilename:
        raise ValueError("Error: Please specify --reference-infile with FASTA sequence")
    elif os.path.exists(myoptions.reference_infilename):
        for _record in SeqIO.parse(myoptions.reference_infilename, "fasta"):
            _padded_reference_dna_seq = str(_record.seq)
            break # parse only the first and supposedly the only entry
        _reference_dna_seq = _padded_reference_dna_seq.replace('-', '')
        _reference_protein_seq = alt_translate(_padded_reference_dna_seq)
        _reference_as_codons = get_codons(_padded_reference_dna_seq)
    else:
        raise ValueError("Error: File %s does not exist" % myoptions.reference_infilename)

    if myoptions.alignment_infilename and os.path.exists(myoptions.alignment_infilename):
        _alnfilename_count_handle = open('.'.join(myoptions.alignment_infilename.split('.')[:-1]) + '.count', 'w')
    else:
        _alnfilename_count_handle = open('.'.join(myoptions.alignment_infilename.split('.')[:-1]) + '.count', 'w')

    if myoptions.outfileprefix:
        if myoptions.outfileprefix.endswith('.tsv'):
            _outfilename_handle = open_file(myoptions.outfileprefix)
            _outfilename_unchanged_codons_handle = open_file(myoptions.outfileprefix[:-4] + '.unchanged_codons.tsv')
        else:
            _outfilename_handle = open_file(myoptions.outfileprefix + '.tsv')
            _outfilename_unchanged_codons_handle = open_file(myoptions.outfileprefix + '.unchanged_codons.tsv')
    else:
        raise RuntimeError("Please specify output filename prefix via --outfile-prefix")
    if myoptions.aa_start:
        _aa_start = myoptions.aa_start - 1 # sanitize the input to be one value less than the real position
    else:
        _aa_start = 0

    if myoptions.min_start:
        _min_start = myoptions.min_start - 1 # alter it to pythonic 0-based numbering
    else:
        _min_start = 0
    if myoptions.max_stop:
        _max_stop = myoptions.max_stop # incidentally the 1-based value is matching pythonic slicing expectations
    else:
        _max_stop = 0

    if myoptions.alignment_infilename and os.path.exists(myoptions.alignment_infilename):
        parse_alignment(myoptions.alignment_infilename, _padded_reference_dna_seq, _reference_protein_seq, _reference_as_codons, _outfilename_handle, _outfilename_unchanged_codons_handle, _alnfilename_count_handle, _aa_start, _min_start, _max_stop)
    else:
        raise RuntimeError("Input file %s does not exist or is not defined" % str(myoptions.alignment_infilename))


if __name__ == "__main__":
    main()

# vim:ts=4:sw=4:expandtab:smartindent
