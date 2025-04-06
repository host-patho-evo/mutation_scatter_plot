#! /usr/bin/env python3

# This work © 2025 by Jiří Zahradník and Martin Mokrejš
# (First Medical Faculty - Charles University in Prague) is licensed under
# Creative Commons Attribution 4.0 International. To view a copy of this
# license, visit https://creativecommons.org/licenses/by/4.0/

"""This program calculate_codon_frequencies.py can parse multi-FASTA *.aln
files (ALN) of protein-coding regions along with GFF3 annotation file to
read a codon (or amino acid residue) from the reference sequence and calculate
frequencies of all up to 64 codons (or up to 20 amino acids) in the respective
column of the ALN file. It also counts number of DELeletion events in the
sample (compared to reference) and also of STOP codons ('*' or 'X' is displayed
in figures depending on the runtime mode).

It uses decimal module to keep numeric precision at maximum and only converts
to floating point numbers during output if not outputting just a string
representation of the number.

To prepare the input data we typically map NGS short sequences from Illumina
platform using 'bwa mem' to a reference protein-coding ORF sequence and then
convert the SAM file format to ALN using `gofasta sam toMultiAlign` command.
The tool gofasta can be obtained from https://github.com/virus-evolution/gofasta .

It does not keep record of INSertions appearing relative to the reference
sequence from the sample sequences but provided in our assays these are just
sequencing errors and not real mutations (as they would break the reading frame
and the proteins would not be synthesized in expression host) we can safely
ignore gofasta discarding these erroneous INSertions.
In principle, 3nt, 6nt, ... long INSertions would be lost but they cannot
appear in our experiments with single nucleotide changes and very very rarely
with double-nucleotide changes.

The utility discards internally codons containing unknown ('N') nucleotides.
It uses Biopython for parsing and translation. The output is a TAB-separated
TSV file with frequencies, e.g.:

439	N	K	0.000310	AAC	AAA

To respect reading frame calculate_codon_frequencies.py follows reference
sequence and parses always 3 nucleotides at once from the ALN input file.
If there was a gap (padding) character '-' if picks an extra character from
the input until it has 3 nucleotides. Then it calculates frequencies for
all possible codons. The results are stored in a TAB-separated TSV file for
easy post-processing by mutation_scatter_plot.py which on-the-fly discards
codons containing unknown (N) nucleotides using python Pandas library and
finally draws interactive figures using Matplotlib and Bokeh graphical
libraries.


Please cite the following article if you use our data or software in your research:

Shoshany A., Tian R., Padilla-Blanco M., Hruška A., Baxova K., Zoler E., Mokrejš M., Schreiber G., Zahradník J. (submitted) Convergence of Directed and Viral Evolution Reveals the Selective Pressures Driving Omicron Emergence.
"""

import os, re
from subprocess import Popen, PIPE, call

from collections import Counter

from optparse import OptionParser
from decimal import Decimal
from numpy import median, average

from Bio import SeqIO
from Bio.Seq import Seq, reverse_complement, translate
from Bio import AlignIO

from BCBio import GFF

version = 2025004050820

myparser = OptionParser()
myparser.add_option("--reference-infile", action="store", type="string", dest="reference_infilename", default=None, metavar="FILE",
    help="FASTA formatted input file with reference padded sequence matching GFF3 contents and the FASTA aln file")
myparser.add_option("--reference-gff3-infile", action="store", type="string", dest="gff3_infilename", default=None, metavar="FILE",
    help="GFF3 formatted annotation of the reference")
myparser.add_option("--alignment-file", action="store", type="string", dest="alignment_infilename", default=None, metavar="FILE",
    help="Alignment file in FASTA format with - (minus) chars to adjust the alignment to the --reference-infile")
myparser.add_option("--outfile-prefix", action="store", type="string", dest="outfileprefix", default=None, metavar="FILE",
    help="Prefix to derive output TSV filename to be created, .tsv and .unchanged_codons.tsv will be appeded to it")
myparser.add_option("--left-offset", action="store", type="int", dest="loffset", default=0,
    help="First nucleotide of the ORF region of interest to be sliced out from the input sequences")
myparser.add_option("--right-offset", action="store", type="int", dest="roffset", default=0,
    help="Last nucleotide of the last codon of interest to be sliced out from the input sequences")
myparser.add_option("--aa_start", action="store", type="int", dest="aa_start", default=0,
    help="Real position of the very first codon unless (1 for an initiator ATG). Add this value to shift the codon position in the output TSV file (the ATG position minus one). Use this if you cannot use --left-offset nor --right-offset which would have been used for slicing the input reference. The value provided is decremented by one to match pythonic 0-based numbering.")
myparser.add_option("--print-unchanged-sites", action="store_true", dest="print_unchanged_sites", default=False,
    help="Print out also sites with unchanged codons in to unchanged_codons.tsv file")
myparser.add_option("--discard-this-many-leading-nucs", action="store", type="int", dest="discard_this_many_leading_nucs", default=0,
    help="Specify how many offending nucleotides are at the front of the ALN sequences shifting the reading frame of the input ALN file from frame +1 so either of two remaining. Count the leading dashes and eventual nucleotides of incomplete codons too and check if it can be divided by 3.0 without slack. By default reading frame +1 is expected and hence no leading nucleotides are discarded.")
myparser.add_option("--minimum-alignments-length", action="store", type="int", dest="minimum_aln_length", default=50,
    help="Minimum length of aligned NGS read to be used for calculations") 
myparser.add_option("--debug", action="store", type="int", dest="debug", default=0,
    help="Set debug level to some integer value")

(myoptions, myargs) = myparser.parse_args()


def get_codons(seq):
    if len(seq) % 3 == 0:
        #split list into codons of three
        codons = [seq[i:i+3] for i in range(0, len(seq), 3)]
    elif len(seq.replace('-','')) % 3 == 0:
        #split list into codons of three
        codons = [seq[i:i+3] for i in range(0, len(seq.replace('-','')), 3)]
        if myoptions.debug: print("Debug: Detected %s minus signs in the sequence but after all the nucleotide sequence can be divided by three when they are omitted, good." % seq.count('-'))
    else:
        raise ValueError("Error: sequence %s cannot be divided by 3 and removing minus signs does not help either" % seq)
    return codons


def parse_alignment(alignment_file, padded_reference_dna_seq, reference_protein_seq, reference_as_codons, outfilename, outfilename_unchanged_codons, alnfilename_count, aa_start):
    if myoptions.debug: print("Debug: Depadded reference sequence has length %d, padded reference sequence has length %d, each entry from %s must also have same padded length" % (len(padded_reference_dna_seq.replace('-','')), len(padded_reference_dna_seq), alignment_file))
    if myoptions.debug: print("Debug: reference_protein_seq=%s with length %d" % (str(reference_protein_seq), len(reference_protein_seq)))
    try:
        align = AlignIO.read(alignment_file, "fasta")
    except ValueError:
        raise ValueError("Error: one of the entries in the %s file has different length" % alignment_file)
    if myoptions.loffset or myoptions.roffset:
        padded_reference_dna_seq = padded_reference_dna_seq[max(myoptions.loffset - 1, 0) : min(len(padded_reference_dna_seq), myoptions.roffset)] # discard leading sequence of the reference not contained in ALN file, if the loffset is correct
        reference_protein_seq = reference_protein_seq[int(max(myoptions.loffset - 1, 0) / 3): int(min(len(reference_protein_seq), myoptions.roffset / 3))]
        reference_as_codons = reference_as_codons[int(max(myoptions.loffset - 1, 0) / 3): int(min(len(reference_as_codons), myoptions.roffset / 3))]
        if True or myoptions.debug:
            print("Info: len(padded_reference_dna_seq)=%s, myoptions.loffset=%s, myoptions.roffset=%s" % (len(padded_reference_dna_seq), myoptions.loffset, myoptions.roffset))
            print("Info: After cutting input using offset position: %s" % padded_reference_dna_seq)
            print("Info: After cutting input using offset position: %s" % reference_protein_seq)
            print("Info: After cutting input using offset position: %s" % reference_as_codons)
        if not reference_protein_seq:
            raise ValueError("Error: No reference_protein_seq provided or left. Was the slicing using myoptions.loffset=%s, myoptions.roffset=%s wrong?" % (myoptions.loffset, myoptions.roffset))
        if not reference_as_codons:
            raise ValueError("Error: No reference_as_codons provided or left. Was the slicing using myoptions.loffset=%s, myoptions.roffset=%s wrong?" % (myoptions.loffset, myoptions.roffset))

        if True or myoptions.debug: print("Debug: Depadded reference sequence has length %d, padded reference sequence has length %d, each entry from %s must also have same padded length" % (len(padded_reference_dna_seq.replace('-','')), len(padded_reference_dna_seq), alignment_file))
    if len(padded_reference_dna_seq) != len(align[0].seq):
        raise ValueError("Error: Length %s of the padded reference sequence %s does not match length %s of the first alignment item %s read from thr ALN file" % (len(padded_reference_dna_seq), padded_reference_dna_seq, len(align[0].seq), align[0].seq))
    _reference_aa_index = 0
    _reference_aa = reference_protein_seq[_reference_aa_index] # fetch the very first aa residue, supposedly 'M' but above we already chopped the sequence and cut only a region of interest
    _previous_gaps = 0
    _new_gaps = 0
    _new_deleted_gaps = 0
    _re_leading_gaps = re.compile("^[-Nn]+")
    _re_trailing_gaps = re.compile("[-Nn]+$")
    _already_checked_starts = [] # keep list of the pythonic codon _starts of the for loop below so that some check are performed only once per iteration, not for every input alignment line
    _top_most_codons = []
    _total_aln_entries_used = 0
    for _start in range(0, len(padded_reference_dna_seq), 3): # generate pythonic numbers so the pythonic slicing below works, the input sequence must be in frame +1
        _unchanged_codons = []
        _changed_codons = []
        _deleted_codons = []
        _unchanged_aa_residues = []
        _changed_aa_residues = []
        _deleted_aa_residues = []
        _new_gaps = 0
        _new_deleted_gaps = 0
        _is_deletion = False
        _is_insertion = False
        _new_aa_residue = None
        if myoptions.debug: print("Debug0: _start=%s, padded_reference_dna_seq=%s" % (_start, padded_reference_dna_seq))
        for _aln_line in align:
            if not _start:
                # count lines only when iterating over the first codon
                _total_aln_entries_used += 1
            if myoptions.discard_this_many_leading_nucs:
                _aln_line_seq = str(_aln_line.seq)[myoptions.discard_this_many_leading_nucs:]
            else:
                _aln_line_seq = str(_aln_line.seq)
            _padded_aln_line_length = len(_aln_line_seq) # padded length
            _depadded_aln_line_length = len(_aln_line_seq.replace('-','').lstrip('Nn').rstrip('Nn'))
            _new_gaps = 0 # zap value from previous codon search if it contained a minus char
            _new_deleted_gaps = 0 # zap value from previous codon search if it contained a minus char
            _codon = None
            _deleted_codon = None
            _deleted_codon_depadded = None
            _rough_sample_codon = None
            _reference_codon = None
            _new_aa_residue = None

            _end_of_leading_gaps = 0
            _start_of_trailing_gaps = 0
            for _match in _re_leading_gaps.finditer(_aln_line_seq):
                _end_of_leading_gaps = _match.end()
            for _match in _re_trailing_gaps.finditer(_aln_line_seq):
                _start_of_trailing_gaps = _match.start()
            if myoptions.debug: print("Debug0: End of leading gaps is at %s, start of trailing gaps is at %s, %s" % (_end_of_leading_gaps, _start_of_trailing_gaps, _aln_line_seq))

            if myoptions.minimum_aln_length and _depadded_aln_line_length < myoptions.minimum_aln_length:
                pass # skip alignments shorter than 50nt
            elif not _depadded_aln_line_length:
                pass # skip empty alignments
            elif _end_of_leading_gaps and _start < _end_of_leading_gaps + 1:
                pass # skip leading gaps so they do not get counted as DELetions
            elif _start_of_trailing_gaps and _start + 1 > _start_of_trailing_gaps:
                pass # skip trailing gaps so they do not get counted as DELetions
            elif not _end_of_leading_gaps and _start + 1 > _start_of_trailing_gaps and _start_of_trailing_gaps:
                pass # skip trailing gaps so they do not get counted as DELetions but do not enter here is there are no leading or trailing dashes at all causing _start_of_trailing_gaps=0 _end_of_leading_gaps=0
                # Q409del
            elif not _start_of_trailing_gaps and _start < _end_of_leading_gaps + 1 and _end_of_leading_gaps:
                pass # skip leading gaps so they do not get counted as DELetions but do not enter here is there are no leading or trailing dashes at all causing _start_of_trailing_gaps=0 _end_of_leading_gaps=0
            elif _padded_aln_line_length + 1 > _start + _previous_gaps + 3 + _new_gaps: # avoid infinite loop when increasing slice end position behind the existing input sequence
                _rough_sample_codon = _aln_line_seq[_start+_previous_gaps:_start + _previous_gaps + 3 + _new_gaps]
                _reference_codon = padded_reference_dna_seq[_start+_previous_gaps:_start + _previous_gaps + 3 + _new_gaps]
                if myoptions.debug: print("Debug1: Rough codon is %s, reference codon is %s, refseq length %s, sliced [%d:%d] (pythonic slice numbering)" % (_rough_sample_codon, _reference_codon, _padded_aln_line_length, _start+_previous_gaps, _start + _previous_gaps + 3 + _new_gaps))
                if _rough_sample_codon == 'NNN' or _reference_codon == 'NNN':
                    pass
                elif _rough_sample_codon == '---' and _reference_codon != '---' and _reference_codon != 'NNN':
                    _new_deleted_gaps = 0 # make no sense to increase it
                    _deleted_codon = padded_reference_dna_seq[_start+_previous_gaps:_start + _previous_gaps + 3 + _new_deleted_gaps]
                    _deleted_codon_depadded = _deleted_codon.replace('-','')
                    if myoptions.debug: print("Debug2: Deleted codon is %s, reference codon is %s, refseq length %s, end of slice is %d" % (_deleted_codon, _reference_codon, _padded_aln_line_length, _start + _previous_gaps + 3 + _new_gaps))
                    while len(_deleted_codon_depadded) > 0 and len(_deleted_codon_depadded) < 3: # slice wider if we hit gaps and have less than 3 nucleotides
                        _new_deleted_gaps += 1
                        if myoptions.debug: print("Debug2: Increased _new_deleted_gaps by one to %d" % _new_deleted_gaps)
                        _deleted_codon = padded_reference_dna_seq[_start+_previous_gaps:_start + _previous_gaps + 3 + _new_deleted_gaps]
                        _deleted_codon_depadded = _deleted_codon.replace('-','')
                    if _deleted_codon_depadded:
                        if _deleted_codon_depadded != 'NNN' and len(_deleted_codon_depadded) == 3:
                            _deleted_codons += [ _deleted_codon_depadded ]
                            _new_aa_residue = translate(_deleted_codon_depadded)
                            _deleted_aa_residues += [ _new_aa_residue ]
                            #_deleted_codons += [ '---' ]
                            #_new_aa_residue = 'None'
                            #_deleted_aa_residues += [ 'None' ]
                            _is_deletion = True
                            _codon = '---'
                            _codon_depadded = 'None'
                            _new_aa_residue = 'None'
                            if myoptions.debug: print("Debug3: Added Deleted codon is %s, reference codon is %s, refseq length %s, sliced [%d:%d] (pythonic slice numbering), aa_residue is %s" % ('---', _reference_codon, _padded_aln_line_length, _start+_previous_gaps, _start + _previous_gaps + 3 + _new_gaps, 'None'))
                elif _reference_codon == '---' and _rough_sample_codon != '---':
                    _is_insertion = True
                    _rough_sample_codon = _aln_line_seq[_start+_previous_gaps:_start+_previous_gaps+3+_new_gaps]
                    _codon_depadded = _rough_sample_codon.replace('-','')
                    _reference_codon = padded_reference_dna_seq[_start+_previous_gaps:_start + _previous_gaps + 3 + _new_gaps]
                    _reference_codon_depadded = _reference_codon.replace('-','')
                    _changed_codons += [ _codon_depadded ]
                    _new_aa_residue = translate(_codon_depadded)
                    if _new_aa_residue:
                        _changed_aa_residues += [ _new_aa_residue ]
                    else:
                        pass # do not count junk if reference was anyway just gaps
                    _is_insertion = True
                else:
                    _rough_sample_codon = _aln_line_seq[_start+_previous_gaps:_start+_previous_gaps+3+_new_gaps]
                    _codon_depadded = _rough_sample_codon.replace('-','')
                    _reference_codon = padded_reference_dna_seq[_start+_previous_gaps:_start + _previous_gaps + 3 + _new_gaps]
                    _reference_codon_depadded = _reference_codon.replace('-','')
                    if myoptions.debug: print("Debug: _padded_aln_line_length=%d" % _padded_aln_line_length)
                    if myoptions.debug: print("Debug: Rough codon is %s, reference codon is %s, refseq length %s, sliced [%d:%d] (pythonic slice numbering), _codon_depadded is %s, _reference_codon_depadded is %s" % (_rough_sample_codon, _reference_codon, _padded_aln_line_length, _start+_previous_gaps, _start+_previous_gaps+3+_new_gaps, _codon_depadded, _reference_codon_depadded))
                    while len(_reference_codon_depadded) < 3 and len(_codon_depadded) > 0 and len(_codon_depadded) < 3 and _padded_aln_line_length > _start+_previous_gaps+3+_new_gaps: # slice wider if we hit gaps and have less than 3 nucleotides
                        _new_gaps += 1
                        _rough_sample_codon = _aln_line_seq[_start+_previous_gaps:_start+_previous_gaps+3+_new_gaps]
                        _codon_depadded = _rough_sample_codon.replace('-','')
                        _reference_codon = padded_reference_dna_seq[_start+_previous_gaps:_start + _previous_gaps + 3 + _new_gaps]
                        _reference_codon_depadded = _reference_codon.replace('-','')
                        if myoptions.debug: print("Debug4: Rough codon is %s, reference codon is %s, refseq length %s, sliced [%d:%d] (pythonic slice numbering)" % (_rough_sample_codon, _reference_codon, _padded_aln_line_length, _start+_previous_gaps, _start+_previous_gaps+3+_new_gaps))
                    if myoptions.debug:
                        print("Debug5: Rough codon is %s, reference codon is %s, refseq length %s, sliced [%d:%d] (pythonic slice numbering), _codon_depadded is %s, _reference_codon_depadded is %s" % (_rough_sample_codon, _reference_codon, _padded_aln_line_length, _start+_previous_gaps, _start+_previous_gaps+3+_new_gaps, _codon_depadded, _reference_codon_depadded))

                    # skip empty codon parsed just after the end of the sequence
                    if _codon_depadded.upper() == _reference_codon_depadded.upper():
                        # skip 'AT-C' versus 'ATC' but increase the incidence counter
                        _unchanged_codons += [ _codon_depadded ]
                        _unchanged_aa_residues += translate(_codon_depadded)
                    elif _codon_depadded and len(_codon_depadded) == 3:
                        _changed_codons += [ _codon_depadded ]
                        _new_aa_residue = translate(_codon_depadded)
                        if 'N' in _codon_depadded and _new_aa_residue != 'X':
                            if _codon_depadded not in ('TCN', 'CTN', 'GTN', 'CCN', 'ACN', 'GCN', 'CGN', 'GGN'):
                                raise ValueError("Error: biopython translated %s into %s" % (_codon_depadded, _new_aa_residue))
                            else:
                                if myoptions.debug: print("Debug6: biopython translated %s into %s" % (_codon_depadded, _new_aa_residue))
                        if _new_aa_residue:
                            if _new_aa_residue != _reference_aa:
                                _changed_aa_residues += [ _new_aa_residue ]
                            else:
                                _unchanged_aa_residues += [ _reference_aa ]
                        else:
                            _changed_aa_residues += [ 'None' ]
                        _rough_sample_codon = _aln_line_seq[_start+_previous_gaps:_start + _previous_gaps + 3 + _new_gaps]
                        if myoptions.debug: print("Debug7: Final codon is %s, reference codon is %s, refseq length %s, sliced [%d:%d] (pythonic slice numbering), aa_residue is %s" % (_rough_sample_codon, _reference_codon, _padded_aln_line_length, _start+_previous_gaps, _start+_previous_gaps+3+_new_gaps, _new_aa_residue))
                    else:
                        _new_aa_residue = translate(_codon_depadded)
                        if myoptions.debug: print("Debug8: Discarded frame-breaking codon %s, reference codon is %s, refseq length %s, sliced [%d:%d] (pythonic slice numbering), aa_residue is %s" % (_rough_sample_codon, _reference_codon, _padded_aln_line_length, _start+_previous_gaps, _start+_previous_gaps+3+_new_gaps, _new_aa_residue))
        # all entries in ALN file were processed for the current codon column
        if myoptions.debug: print("Debug9: Lists of codons and aminoacid residues: %s-%s %s %s" % (_start+1, _start+3, _changed_codons, _changed_aa_residues)) # undo off-by-one error due to pythonic counting

        if myoptions.debug: print("Debug10: Deleted codons: %s" % _deleted_codons)
        if myoptions.print_unchanged_sites:
            _unchanged_codon_freq = Counter(_unchanged_codons)
            _unchanged_aa_freq = Counter(_unchanged_aa_residues) # TODO: unused

        _deleted_codon_freq = Counter(_deleted_codons)
        _deleted_aa_freq = Counter(_deleted_aa_residues)

        _changed_codon_freq = Counter(_changed_codons)
        _changed_aa_freq = Counter(_changed_aa_residues) # the values also include 'None' eventually

        _total_codon_freq = Counter(_deleted_codons + _changed_codons + _unchanged_codons)
        _total_codons_counted = sum(_total_codon_freq.values())

        _reference_codon = padded_reference_dna_seq[_start+_previous_gaps:_start + _previous_gaps + 3 + _new_gaps]
        if len(reference_protein_seq) > _reference_aa_index:
            _reference_aa = reference_protein_seq[_reference_aa_index]
            # otherwise keep previous values (actually keep the C-term residue)

        # this is performed for every input ALN line so it slows down the execution, perform the check only once per reference codon
        if _start not in _already_checked_starts and len(_reference_codon) == 3:
            if translate(_reference_codon) != _reference_aa:
                raise ValueError("The reference aa should be '%s' but parsed reference codon '%s' encodes '%s'. The _start = %s" % (_reference_aa, _reference_codon, translate(_reference_codon), _start))
            else:
                _already_checked_starts.append(_start)

        if myoptions.debug:
            if myoptions.debug > 1:
                for _key in _changed_codon_freq:
                    print("Debug: {}: {} = {:8.6f} {} {}".format(_key, _changed_codon_freq[_key], _changed_codon_freq[_key] / _total_codons_counted, set([x for x in Counter(_changed_codons)]), set([x for x in Counter(_changed_aa_residues)])))
        #if sum([ x for x in _changed_codon_freq]) != sum([x for x in _changed_aa_freq]):
        #    raise ValueError("Error: len(_changed_codon_freq)=%d != len(_changed_aa_freq)=%d, _changed_codon_freq=%s, _changed_aa_freq=%s" % (len(_changed_codon_freq), len(_changed_aa_freq), str(_changed_codon_freq), str(_changed_aa_freq)))

        if myoptions.print_unchanged_sites:
            for _some_codon in _unchanged_codon_freq:
                _some_aa = translate(_some_codon)
                outfilename_unchanged_codons.write("{}\t{}\t{}\t{:8.6f}\t{}\t{}\n".format(_reference_aa_index + 1 + int( myoptions.loffset / 3.0) + aa_start, _reference_aa, _some_aa, Decimal(_unchanged_codon_freq[_some_codon]) / Decimal(_total_codons_counted), _reference_codon, _some_codon))
        #
        # print contents of _unchanged_aa_freq
        if len(_changed_codon_freq):
            if not _is_insertion: # print only sites where there is > 1 codon used
                for _some_codon in sorted(_changed_codon_freq):
                    _some_aa = translate(_some_codon)
                    if myoptions.debug: print("Debug11: _reference_aa_index = %s, len(reference_protein_seq)=%d" % (_reference_aa_index, len(reference_protein_seq)))
                    if (_reference_codon.upper() != _some_codon.upper() and translate(padded_reference_dna_seq[_reference_aa_index * 3]) != _some_aa):
                        # BUG: _reference_aa_index is 482 in aa numbering of depadded protein sequence 1446/3.0=482
                        #      what we are after is a position in the protein but with leading Ns omitted
                        #      until that is fixed, reference_protein_seq[_reference_aa_index] does reference_protein_seq[482] which is outside
                        outfilename.write("{}\t{}\t{}\t{:8.6f}\t{}\t{}\n".format(_reference_aa_index + 1 + int( myoptions.loffset / 3.0) + aa_start, _reference_aa, _some_aa, Decimal(_changed_codon_freq[_some_codon]) / Decimal(_total_codons_counted), _reference_codon, _some_codon))
                    elif translate(padded_reference_dna_seq[_reference_aa_index * 3]) == _some_aa:
                        if myoptions.debug: print("Debug12: translate(padded_reference_dna_seq[_reference_aa_index * 3])=%s, _some_aa=%s" % (translate(padded_reference_dna_seq[_reference_aa_index * 3]), _some_aa))
                        outfilename.write("{}\t{}\t{}\t{:8.6f}\t{}\t{}\n".format(_reference_aa_index + 1 + int( myoptions.loffset / 3.0) + aa_start, _reference_aa, _some_aa, Decimal(_changed_codon_freq[_some_codon]) / Decimal(_total_codons_counted), _reference_codon, _some_codon))
                    else:
                        outfilename.write("{}\t{}\t{}\t{:8.6f}\t{}\t{}\t{}\n".format(_reference_aa_index + 1 + int( myoptions.loffset / 3.0) + aa_start, _reference_aa, _some_aa, Decimal(_changed_codon_freq[_some_codon]) / Decimal(_total_codons_counted), _reference_codon, _some_codon, 'LOST'))
                        raise ValueError()
            # handle insertion in sample if reference contained a gap otherwise reference_protein_seq[_reference_aa_index] fails with IndexError
            if _is_insertion:
                for _some_codon in sorted(_changed_codon_freq):
                    _some_aa = translate(_some_codon)
                    if myoptions.debug: print("Debug13: Reference contains a full gap spanning a whole aminoacid position in padded range %d-%d, NOT increasing _reference_aa_index now. _previous_gaps=%d, _new_gaps=%s" % (_start+_previous_gaps + 1, _start + _previous_gaps + 3 + _new_gaps, _previous_gaps, _new_gaps)) # undo off-by-one error
                    outfilename.write("{}\t{}\t{}\t{:8.6f}\t{}\t{}\n".format(_reference_aa_index + int( myoptions.loffset / 3.0) + aa_start, 'INS', _some_aa, Decimal(_changed_codon_freq[_some_codon]) / Decimal(_total_codons_counted), _reference_codon, _some_codon))

        if _is_deletion:
            for _some_deleted_codon, _some_deleted_aa in zip(_deleted_codon_freq, _deleted_aa_residues):
                outfilename.write("{}\t{}\t{}\t{:8.6f}\t{}\t{}\n".format(_reference_aa_index + 1 + int( myoptions.loffset / 3.0) + aa_start, _reference_aa, 'DEL', Decimal(_deleted_codon_freq[_some_deleted_codon]) / Decimal(_total_codons_counted), _some_deleted_codon, '---'))

        # get the top-most codon in this ALN column
        try:
            _top_most_codon, _top_most_count = _total_codon_freq.most_common()[0]
        except IndexError:
            pass
        else:
            _top_most_codons.append(_top_most_codon)

        if not _is_deletion and len(reference_protein_seq) > _reference_aa_index + 1: # undo off-by-one error while stop increasing after the C-terminus of the reference protein
            _reference_aa_index += 1
            if myoptions.debug: print("Debug14: Increased _reference_aa_index to %d, moving to next codon" % _reference_aa_index)
        if _is_deletion:
            _reference_aa_index += 1 # increase the counter in extra for every DEL traversed
            if myoptions.debug: print("Debug15: Increased _reference_aa_index to %d, moving to next codon" % _reference_aa_index)
        if _is_insertion:
            _reference_aa_index -= 1 # decrement the counter for every INS traversed
            if myoptions.debug: print("Debug16: Decreased _reference_aa_index to %d, moving to next codon" % _reference_aa_index)

        _previous_gaps += _new_gaps
        _new_gaps = 0
        _new_deleted_gaps = 0
    alnfilename_count.write("%s\n" % _total_aln_entries_used)
    _consensus = ''.join(_top_most_codons)
    print("Info: consensus = %s" % str(_consensus))
    if _consensus in padded_reference_dna_seq:
        print("Info: Sample consensus sequence should roughly match substring inside the reference and IT DOES: %s" % str(_top_most_codons))
    else:
         print("Info: Sample consensus sequence should roughly match substring inside the reference BUT IT DOES NOT, maybe due to some true major mutations in some codons: %s" % str(_top_most_codons))


def main():
    # parse the reference DNA
    if not myoptions.reference_infilename:
        raise ValueError("Error: Please specify --reference-infile with FASTA sequence matching --reference-gff3-infile")
    elif os.path.exists(myoptions.reference_infilename):
        _padded_reference_dna_seq_fh = open(myoptions.reference_infilename)
        _padded_reference_dna_seqs = SeqIO.to_dict(SeqIO.parse(_padded_reference_dna_seq_fh, "fasta"))
        _padded_reference_dna_seq_fh.close()
    else:
        raise ValueError("Error: File %s does not exist" % myoptions.reference_infilename)

    # parse GFF3 annotation and decorate further the already existing Seq objects
    if myoptions.gff3_infilename and os.path.exists(myoptions.gff3_infilename):
        in_handle = open(myoptions.gff3_infilename)
        for rec in GFF.parse(in_handle, base_dict=_padded_reference_dna_seqs):
            #print(rec)
            #print(rec.features[0])
            pass
        in_handle.close()
        if myoptions.debug: print("Debug15: Parsed %d entries from %s and %s" % (len(_padded_reference_dna_seqs), myoptions.reference_infilename, myoptions.gff3_infilename))

        #print(dir(rec))
        # ['_AnnotationsDict', '_AnnotationsDictValue', '__add__', '__annotations__', '__bool__', '__bytes__', '__class__', '__contains__', '__delattr__', '__dict__', '__dir__', '__doc__', '__eq__', '__format__', '__ge__', '__getattribute__', '__getitem__', '__getstate__', '__gt__', '__hash__', '__init__', '__init_subclass__', '__iter__', '__le__', '__len__', '__lt__', '__module__', '__ne__', '__new__', '__radd__', '__reduce__', '__reduce_ex__', '__repr__', '__setattr__', '__sizeof__', '__slotnames__', '__str__', '__subclasshook__', '__weakref__', '_per_letter_annotations', '_seq', '_set_per_letter_annotations', '_set_seq', 'annotations', 'count', 'dbxrefs', 'description', 'features', 'format', 'id', 'islower', 'isupper', 'letter_annotations', 'lower', 'name', 'reverse_complement', 'seq', 'translate', 'upper']
        # /sequence-region=[('7-WU-FF', 978, 1440)]
    
        #print(str(rec.features))
        # [SeqFeature(SimpleLocation(ExactPosition(978), ExactPosition(1440), strand=1), type='gene', id='gene:S', qualifiers=...)]

        #print(str(rec.features[0]))

        #print(dir(rec.features[0]))
        # ['__bool__', '__class__', '__contains__', '__delattr__', '__dict__', '__dir__', '__doc__', '__eq__', '__format__', '__ge__', '__getattribute__', '__getstate__', '__gt__', '__hash__', '__init__', '__init_subclass__', '__iter__', '__le__', '__len__', '__lt__', '__module__', '__ne__', '__new__', '__reduce__', '__reduce_ex__', '__repr__', '__setattr__', '__sizeof__', '__str__', '__subclasshook__', '__weakref__', '_flip', '_get_ref', '_get_ref_db', '_get_strand', '_set_ref', '_set_ref_db', '_set_strand', '_shift', 'extract', 'id', 'location', 'qualifiers', 'ref', 'ref_db', 'strand', 'sub_features', 'translate', 'type']

        #print(dir(rec.features[0].location))

        _start = rec.features[0].location.start # 1st column if gff3 does not match reference FASTA ID
        _stop = rec.features[0].location.end
        _strand = rec.features[0].location.strand
        #print(rec.seq[_start:_stop])

        _padded_reference_dna_seq = str(rec.seq)
        _reference_dna_seq = _padded_reference_dna_seq.replace('-', '')
        _reference_protein_seq = translate(_reference_dna_seq)
        _reference_as_codons = get_codons(_reference_dna_seq)

        print("Info: Before cutting input using offset position _padded_reference_dna_seq: %s" % _padded_reference_dna_seq)
        print("Info: Before cutting input using offset position _reference_protein_seq: %s" % _reference_protein_seq)
        print("Info: Before cutting input using offset position _reference_as_codons: %s" % _reference_as_codons)
    else:
        raise ValueError("Error: Please specify existing file instead of --reference-gff3-infile=%s" % str(myoptions.gff3_infilename))

    # this skips many unmapped read from the input SAM stream
    # samtools view -h /tmp/testcases.sam | ./gofasta sam toMultiAlign --start 978 --end 1440 --pad -o /tmp/testcases.aln

    if myoptions.alignment_infilename and os.path.exists(myoptions.alignment_infilename):
        _alnfilename_count_handle = open(myoptions.alignment_infilename + '.count', 'w')
    else:
        _alnfilename_count_handle = open(myoptions.alignment_infilename + '.count', 'w')

    if myoptions.outfileprefix:
        if myoptions.outfileprefix.endswith('.tsv'):
            _outfilename = open(myoptions.outfileprefix, 'w')
            _outfilename_unchanged_codons = open(myoptions.outfileprefix[:-4] + '.unchanged_codons.tsv', 'w')
        else:
            _outfilename = open(myoptions.outfileprefix + '.tsv', 'w')
            _outfilename_unchanged_codons = open(myoptions.outfileprefix + '.unchanged_codons.tsv', 'w')
    else:
        raise RuntimeError("Please specify output filename prefix via --outfile-prefix")
    if myoptions.aa_start:
        _aa_start = myoptions.aa_start - 1 # sanitize the input to be one value less than the real position
    else:
        _aa_start = 0
    if myoptions.alignment_infilename and os.path.exists(myoptions.alignment_infilename):
        parse_alignment(myoptions.alignment_infilename, _padded_reference_dna_seq, _reference_protein_seq, _reference_as_codons, _outfilename, _outfilename_unchanged_codons, _alnfilename_count_handle, _aa_start)
    else:
        raise RuntimeError("Input file %s does not exist or is not defined" % str(myoptions.alignment_infilename))


if __name__ == "__main__":
    main()

# vim:ts=4:sw=4:expandtab:smartindent
