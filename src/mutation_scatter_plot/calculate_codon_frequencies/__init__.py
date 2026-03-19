# This work © 2025 by Jiří Zahradník and Martin Mokrejš
# (First Medical Faculty - Charles University in Prague) is licensed under
# Creative Commons Attribution 4.0 International. To view a copy of this
# license, visit https://creativecommons.org/licenses/by/4.0/

"""calculate_codon_frequencies — parse multi-FASTA alignments and calculate
codon/amino-acid frequencies per site relative to a reference sequence.

Please cite our article if you use our data or software in your work:

Shoshany A., Tian R., Padilla-Blanco M., Hruška A., Baxova K., Zoler E.,
Mokrejš M., Schreiber G., Zahradník J. (submitted) In Vitro and Viral
Evolution Convergence Reveal the Selective Pressures Driving Omicron
Emergence. https://www.biorxiv.org/content/10.1101/2025.04.23.650148v1
"""

import os
import re
import sys

from collections import Counter
from decimal import Decimal

from Bio import AlignIO

from ..utils import alt_translate

VERSION = "0.3"

__all__ = [
    "VERSION",
    "get_codons",
    "write_tsv_line",
    "parse_alignment",
    "open_file",
]


def get_codons(seq, debug=False):
    """Parse reference sequence into a list of codon triplets.

    If the sequence is not already divisible by 3, padding dashes are removed
    before splitting.
    """
    _seq_len = len(seq)
    if _seq_len % 3 == 0:
        _codons = [seq[_i:_i+3] for _i in range(0, _seq_len, 3)]
    else:
        _seq_depadded = seq.replace('-', '')
        _seq_depadded_len = len(_seq_depadded)
        if _seq_depadded_len % 3 == 0:
            _codons = [seq[_i:_i+3] for _i in range(0, _seq_depadded_len, 3)]
            if debug:
                print("Debug: Detected %s minus signs in the sequence but after all "
                      "the nucleotide sequence can be divided by three when they are "
                      "omitted, good." % seq.count('-'))
        else:
            raise ValueError(
                "Error: Sequence %s cannot be divided by 3 and removing minus "
                "signs does not help either" % seq
            )
    return _codons


def write_tsv_line(outfilename, codons, natural_codon_position_padded,
                   natural_codon_position_depadded, reference_aa,
                   total_codons_per_site_sum, reference_codon, debug=False):
    """Write one or more TSV lines for all codons observed at a site."""
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
        outfilename.write(
            "{}\t{}\t{}\t{}\t{:8.6f}\t{}\t{}\t{}\t{}\n".format(
                natural_codon_position_padded,
                natural_codon_position_depadded,
                reference_aa,
                _some_aa,
                _observed_codon_count2 / Decimal(total_codons_per_site_sum),
                reference_codon,
                _some_codon,
                _observed_codon_count2,
                _total_codons_per_site_sum,
            )
        )
        if debug:
            print("TESTING1:\t{}\t{}\t{}\t{}\t{:8.6f}\t{}\t{}\t{}\t{}".format(
                natural_codon_position_padded,
                natural_codon_position_depadded,
                reference_aa,
                _some_aa,
                _observed_codon_count2 / Decimal(total_codons_per_site_sum),
                reference_codon,
                _some_codon,
                _observed_codon_count2,
                _total_codons_per_site_sum,
            ))
    outfilename.flush()


def parse_alignment(myoptions, alignment_file, padded_reference_dna_seq,
                    reference_protein_seq, reference_as_codons,
                    outfilename, outfilename_unchanged_codons,
                    alnfilename_count, aa_start, min_start, max_stop):
    """Parse a padded multi-FASTA alignment and write codon frequency TSV files.

    left_reference_offset and right_reference_offset are used to slice the
    reference. discard_this_many_leading_nucs and
    discard_this_many_trailing_nucs are used to discard offending
    leading/trailing nucleotides.

    If the reference protein is shorter than the sample entries, the trailing codons
    after reference protein terminated are treated as INSertions with aa position
    of the aminoacid residue corresponding to the last aminoacid codon of the
    reference protein. But provided there are still codon columns to be parsed
    there are multiple lines in the resulting TSV files, which later on breaks
    mutation_scatter_plot because multiple lines are retrieved via Pandas instead
    of just a single row. A workaround so far is to ignore lines with INSertion
    event from the input .frequency.tsv file.

    333	332	R	M	0.108108	AGG	ATG	4	37
    333	332	R	T	0.027027	AGG	ACC	1	37
    333	332	R	T	0.027027	AGG	ACG	1	37
    334	333	R	R	0.432432	AGG	AGA	16	37
    334	333	R	S	0.027027	AGG	AGC	1	37
    335	333	X	R	0.405405	GA-	AGG	15	37
    335	333	X	X	0.108108	GA-	AG	4	37
    335	333	X	X	0.027027	GA-	AC	1	37
    335	333	X	X	0.027027	GA-	GT	1	37
    336	333	INS	E	0.405405	---	GAG	15	37
    336	333	INS	X	0.243243	---	G	9	37
    336	333	INS	X	0.351351	---	T	13	37
    337	333	INS	X	0.891892	---	C	33	37
    337	333	INS	X	0.081081	---	T	3	37
    337	333	INS	P	0.027027	---	CCC	1	37
    338	333	INS	R	0.027027	---	AGA	1	37
    339	333	INS	W	0.027027	---	TGG	1	37
    340	333	INS	X	0.675676	---	AT	25	37
    340	333	INS	X	0.243243	---	AA	9	37
    340	333	INS	X	0.054054	---	AC	2	37
    340	333	INS	G	0.027027	---	GGT	1	37
    341	333	INS	P	0.864865	---	CCT	32	37
    341	333	INS	H	0.108108	---	CAT	4	37
    341	333	INS	L	0.027027	---	CTT	1	37
    342	333	INS	L	0.972973	---	CTG	36	37
    342	333	INS	M	0.027027	---	ATG	1	37
    343	333	INS	L	0.513514	---	CTA	19	37
    343	333	INS	L	0.027027	---	TTA	1	37
    343	333	INS	A	0.081081	---	GCA	3	37
    343	333	INS	S	0.027027	---	TCA	1	37
    343	333	INS	D	0.297297	---	GAC	11	37
    343	333	INS	Y	0.054054	---	TAC	2	37
    344	333	INS	R	0.648649	---	AGA	24	37
    344	333	INS	T	0.270270	---	ACA	10	37
    344	333	INS	A	0.081081	---	GCA	3	37
    345	333	INS	L	0.540541	---	CTG	20	37
    345	333	INS	P	0.108108	---	CCG	4	37
    345	333	INS	F	0.351351	---	TTC	13	37
    346	333	INS	R	0.513514	---	AGG	19	37
    346	333	INS	K	0.027027	---	AAG	1	37
    346	333	INS	M	0.108108	---	ATG	4	37
    346	333	INS	X	0.351351	---	A	13	37
    347	333	INS	P	0.459459	---	CCT	17	37
    347	333	INS	L	0.054054	---	CTT	2	37
    347	333	INS	A	0.027027	---	GCT	1	37
    347	333	INS	X	0.108108	---	GC	4	37
    348	333	INS	P	0.351351	---	CCT	13	37
    348	333	INS	H	0.189189	---	CAT	7	37
    349	333	INS	L	0.540541	---	CTG	20	37
    350	333	INS	X	0.432432	---	CA	16	37
    350	333	INS	X	0.081081	---	CG	3	37
    350	333	INS	H	0.027027	---	CAT	1	37
    351	333	INS	H	0.027027	---	CAT	1	37
    352	333	INS	X	0.621622	---	T	23	37
    352	333	INS	P	0.027027	---	CCT	1	37
    353	333	INS	X	0.135135	---	C	5	37
    353	333	INS	P	0.378378	---	CCT	14	37
    353	333	INS	L	0.027027	---	CTT	1	37
    353	333	INS	H	0.027027	---	CAT	1	37
    353	333	INS	X	0.054054	---	CYT	2	37
    353	333	INS	P	0.027027	---	CCA	1	37
    354	333	INS	P	0.513514	---	CCA	19	37
    354	333	INS	X	0.324324	---	CA	12	37
    354	333	INS	X	0.027027	---	TA	1	37
    355	333	INS	X	0.135135	---	CT	5	37
    355	333	INS	A	0.837838	---	GCT	31	37
    355	333	INS	A	0.027027	---	GCW	1	37
    356	333	INS	L	1.000000	---	CTG	37	37
    357	333	INS	C	0.810811	---	TGT	30	37
    357	333	INS	W	0.189189	---	TGG	7	37
    358	333	INS	L	0.972973	---	CTC	36	37
    358	333	INS	F	0.027027	---	TTC	1	37
    359	333	INS	L	1.000000	---	CTG	37	37
    360	333	INS	H	1.000000	---	CAT	37	37
    """
    if myoptions.debug:
        print("Debug0: Depadded reference sequence has length %d, padded "
              "reference sequence has length %d, each entry from %s must also "
              "have same padded length %d" % (
                  len(padded_reference_dna_seq.replace('-', '')),
                  len(padded_reference_dna_seq),
                  alignment_file,
                  len(padded_reference_dna_seq)))
        print("Debug1: reference_protein_seq=%s with length %d" % (
            str(reference_protein_seq), len(reference_protein_seq)))
    if not os.path.exists(alignment_file):
        raise RuntimeError("Alignment file not found: %s" % alignment_file)
    if os.path.getsize(alignment_file) == 0:
        raise RuntimeError("Alignment file is empty: %s" % alignment_file)
    try:
        _align = AlignIO.read(alignment_file, "fasta")
    except ValueError as exc:
        raise ValueError(
            'Error: one of the entries in the %s file has different length'
            % alignment_file
        ) from exc

    if myoptions.left_reference_offset or myoptions.right_reference_offset:
        _padded_reference_dna_seq = padded_reference_dna_seq[
            max(myoptions.left_reference_offset - 1, 0):
            min(len(padded_reference_dna_seq), myoptions.right_reference_offset)
        ]
        _reference_protein_seq = reference_protein_seq[
            int(max(myoptions.left_reference_offset - 1, 0) / 3):
            int(min(len(reference_protein_seq), myoptions.right_reference_offset / 3))
        ]
        _reference_as_codons = reference_as_codons[
            int(max(myoptions.left_reference_offset - 1, 0) / 3):
            int(min(len(reference_as_codons), myoptions.right_reference_offset / 3))
        ]
        if myoptions.debug:
            print("Info: len(padded_reference_dna_seq)=%s, "
                  "myoptions.left_reference_offset=%s, "
                  "myoptions.right_reference_offset=%s" % (
                      len(padded_reference_dna_seq),
                      myoptions.left_reference_offset - 1,
                      myoptions.right_reference_offset))
            print("Info: After cutting input using offset position: %s" % _padded_reference_dna_seq)
            print("Info: After cutting input using offset position: %s" % _reference_protein_seq)
            print("Info: After cutting input using offset position: %s" % _reference_as_codons)
        if not _reference_protein_seq:
            raise ValueError(
                "Error: No _reference_protein_seq provided or left. Was the "
                "slicing using myoptions.left_reference_offset-1=%s, "
                "myoptions.right_reference_offset=%s wrong?" % (
                    myoptions.left_reference_offset - 1,
                    myoptions.right_reference_offset))
        if not _reference_as_codons:
            raise ValueError(
                "Error: No _reference_as_codons provided or left. Was the "
                "slicing using myoptions.left_reference_offset-1=%s, "
                "myoptions.right_reference_offset=%s wrong?" % (
                    myoptions.left_reference_offset - 1,
                    myoptions.right_reference_offset))
        if myoptions.debug:
            print("Debug2: Depadded reference sequence has length %d, padded "
                  "reference sequence has length %d, each padded entry from %s "
                  "must also have same padded length" % (
                      len(_padded_reference_dna_seq.replace('-', '')),
                      len(_padded_reference_dna_seq),
                      alignment_file))
    else:
        _padded_reference_dna_seq = padded_reference_dna_seq
        _reference_protein_seq = reference_protein_seq
        _reference_as_codons = reference_as_codons

    _zero_based_padded_reference_aa_index = 0
    _reference_aa = _reference_protein_seq[_zero_based_padded_reference_aa_index]
    _reference_codon = _padded_reference_dna_seq[3*_zero_based_padded_reference_aa_index:3*_zero_based_padded_reference_aa_index + 3].upper()
    _reference_codon_depadded = _reference_codon.replace('-', '')
    _previous_gaps = 0
    _new_gaps_in_reference = 0
    _re_leading_gaps = re.compile("^[-Nn]+")
    _re_trailing_gaps = re.compile("[-Nn]+$")
    _already_checked_starts = []
    _top_most_codons = []
    _total_aln_entries_used = 0
    _start_from = myoptions.discard_this_many_leading_nucs
    _stop_to = myoptions.discard_this_many_trailing_nucs

    for _zero_based_codon_startpos in range(
        min_start, max_stop or len(_padded_reference_dna_seq), 3
    ):
        _unchanged_codons = Counter()
        _changed_codons = Counter()
        _deleted_reference_codons = Counter()
        _inserted_codons = Counter()
        _unchanged_aa_residues = Counter()
        _changed_aa_residues = Counter()
        _inserted_aa_residues = Counter() # practically unused
        _deleted_reference_aa_residues = Counter()
        _new_gaps_in_reference = 0
        _is_deletion = False
        _is_insertion = False
        _new_aa_residue = None
        try:
            _reference_aa = _reference_protein_seq[
                _zero_based_padded_reference_aa_index
            ].upper()
        except IndexError as exc:
            raise IndexError(
                'Error: Cannot slice reference protein sequence %s at '
                'position %s' % (_reference_protein_seq,
                                 _zero_based_padded_reference_aa_index)
            ) from exc
        _current_codon_position = _zero_based_codon_startpos / 3 + 1
        _reference_codon = _padded_reference_dna_seq[3*_zero_based_padded_reference_aa_index:3*_zero_based_padded_reference_aa_index + 3].upper()
        _reference_codon_depadded = _reference_codon.replace('-', '')
        if myoptions.debug:
            print("Debug3: _start=%s, _padded_reference_dna_seq=%s" % (
                _zero_based_codon_startpos, _padded_reference_dna_seq))

        for _aln_line in _align:
            _record_id = _aln_line.id
            if myoptions.x_after_count:
                if 'x.' in _record_id:
                    _record_count, _ = _record_id.split('x.')
                    _record_count = int(_record_count)
                else:
                    try:
                        _record_count = int(_record_id.replace('x', ''))
                    except ValueError:
                        # we cannot make an integer from supposedly a string, probably user just enabled x-after-id option but there are just no counts
                        _record_count = 1
            else:
                _record_count = 1

            if not _zero_based_codon_startpos or not _total_aln_entries_used:
                _total_aln_entries_used += _record_count

            if _start_from:
                if _stop_to:
                    _aln_line_seq = str(_aln_line.seq)[_start_from:-_stop_to]
                else:
                    _aln_line_seq = str(_aln_line.seq)[_start_from:]
            elif _stop_to:
                _aln_line_seq = str(_aln_line.seq)[:-_stop_to]
            else:
                _aln_line_seq = str(_aln_line.seq)

            _amplicon_length = (max_stop or len(_padded_reference_dna_seq)) - min_start
            _padded_aln_line_length = len(_aln_line_seq)
            _depadded_aln_line_length = len(
                _aln_line_seq.replace('-', '').lstrip('Nn').rstrip('Nn')
            )
            _new_gaps_in_reference = 0
            _deleted_reference_codon = None
            _rough_sample_codon = _aln_line_seq[_zero_based_codon_startpos:_zero_based_codon_startpos + 3].upper()
            _sample_codon_depadded = _rough_sample_codon.replace('-', '')
            _sample_codon_contained_pad = len(_rough_sample_codon) != len(_sample_codon_depadded)
            _reference_codon_contained_pad = len(_reference_codon) != len(_reference_codon_depadded)
            _new_aa_residue = None

            _start_of_trailing_gaps = 0
            _end_of_leading_gaps = 0
            for _match in _re_leading_gaps.finditer(_aln_line_seq):
                _end_of_leading_gaps = _match.end()
            for _match in _re_trailing_gaps.finditer(_aln_line_seq):
                _start_of_trailing_gaps = _match.start()
            if myoptions.debug:
                print(f"Debug4: End of leading gaps is at {_end_of_leading_gaps:6d}, "
                      f"start of trailing gaps is at {_start_of_trailing_gaps:6d}, "
                      f"{_aln_line_seq}")

            if myoptions.minimum_aln_length and _depadded_aln_line_length < myoptions.minimum_aln_length:
                if myoptions.debug:
                    print("Debug5:Here1, _end_of_leading_gaps=%s, "
                          "_start_of_trailing_gaps=%s, "
                          "_zero_based_codon_startpos=%s" % (
                              _end_of_leading_gaps, _start_of_trailing_gaps,
                              _zero_based_codon_startpos))
            elif not _depadded_aln_line_length:
                if myoptions.debug:
                    print("Debug6:Here2, _end_of_leading_gaps=%s, "
                          "_start_of_trailing_gaps=%s, "
                          "_zero_based_codon_startpos=%s" % (
                              _end_of_leading_gaps, _start_of_trailing_gaps,
                              _zero_based_codon_startpos))
            elif _end_of_leading_gaps and _zero_based_codon_startpos < _end_of_leading_gaps - 1:
                if myoptions.debug:
                    print("Debug7:Here3, _end_of_leading_gaps=%s, "
                          "_start_of_trailing_gaps=%s, "
                          "_zero_based_codon_startpos=%s" % (
                              _end_of_leading_gaps, _start_of_trailing_gaps,
                              _zero_based_codon_startpos))
            elif _start_of_trailing_gaps and _zero_based_codon_startpos + 1 > _start_of_trailing_gaps:
                if myoptions.debug:
                    print("Debug8:Here4, _end_of_leading_gaps=%s, "
                          "_start_of_trailing_gaps=%s, "
                          "_zero_based_codon_startpos=%s" % (
                              _end_of_leading_gaps, _start_of_trailing_gaps,
                              _zero_based_codon_startpos))
            elif not _end_of_leading_gaps and _zero_based_codon_startpos + 1 > _start_of_trailing_gaps and _start_of_trailing_gaps:
                if myoptions.debug:
                    print("Debug9:Here5, _end_of_leading_gaps=%s, "
                          "_start_of_trailing_gaps=%s, "
                          "_zero_based_codon_startpos=%s" % (
                              _end_of_leading_gaps, _start_of_trailing_gaps,
                              _zero_based_codon_startpos))
            elif not _start_of_trailing_gaps and _zero_based_codon_startpos < _end_of_leading_gaps + 1 and _end_of_leading_gaps:
                if myoptions.debug:
                    print("Debug10:Here6, _end_of_leading_gaps=%s, "
                          "_start_of_trailing_gaps=%s, "
                          "_zero_based_codon_startpos=%s" % (
                              _end_of_leading_gaps, _start_of_trailing_gaps,
                              _zero_based_codon_startpos))
            elif _padded_aln_line_length + 1 > _zero_based_codon_startpos + 3:
                if myoptions.debug:
                    print("Debug11: Rough sample codon at %d is %s, reference "
                          "codon is %s, amplicon region length %s (nt), sliced "
                          "[%d:%d] (pythonic slice numbering), "
                          "_sample_codon_depadded is %s, "
                          "_reference_codon_depadded is %s, "
                          "_sample_codon_contained_pad=%s, "
                          "_reference_codon_contained_pad=%s" % (
                              _current_codon_position, _rough_sample_codon,
                              _reference_codon, _amplicon_length,
                              _zero_based_codon_startpos + _previous_gaps,
                              _zero_based_codon_startpos + _previous_gaps + 3 + _new_gaps_in_reference,
                              _sample_codon_depadded, _reference_codon_depadded,
                              str(_sample_codon_contained_pad),
                              str(_reference_codon_contained_pad)))
                if _rough_sample_codon == 'NNN' or _reference_codon == 'NNN':
                    if myoptions.debug:
                        print("Debug11: Skipping a codon containing NNN")
                if _rough_sample_codon == '---' and _reference_codon != '---' and _reference_codon != 'NNN':
                    _is_deletion = True
                    _deleted_reference_codon = str(_reference_codon)
                    if myoptions.debug > 1:
                        print("Debug12: Padded DELeted codon at %d is %s, "
                              "original reference codon is %s, amplicon region "
                              "length %s (nt), end of slice is %d" % (
                                  _current_codon_position,
                                  _deleted_reference_codon,
                                  _reference_codon, _amplicon_length,
                                  _zero_based_codon_startpos + _previous_gaps + 3 + _new_gaps_in_reference))
                    _deleted_reference_codons[_deleted_reference_codon] += _record_count
                    _deleted_aa_residue = alt_translate(_deleted_reference_codon)
                    _deleted_reference_aa_residues[_deleted_aa_residue] += _record_count
                    _sample_codon_depadded = ''
                    if myoptions.debug:
                        print("Debug13: Added DELeted codon at %d is %s, "
                              "reference codon is %s, amplicon region length %s "
                              "(nt), sliced [%d:%d] (pythonic slice numbering), "
                              "aa_residue is %s" % (
                                  _current_codon_position, '---',
                                  _reference_codon, _amplicon_length,
                                  _zero_based_codon_startpos + _previous_gaps,
                                  _zero_based_codon_startpos + _previous_gaps + 3 + _new_gaps_in_reference,
                                  _deleted_aa_residue))
                elif _reference_codon == '---' and _rough_sample_codon != '---':
                    _is_insertion = True
                    _sample_codon_depadded = _rough_sample_codon.replace('-', '')
                    _reference_aa = 'INS'
                    _inserted_codons[_sample_codon_depadded] += _record_count
                    _new_aa_residue = alt_translate(_rough_sample_codon)
                    if _new_aa_residue:
                        _inserted_aa_residues[_new_aa_residue] += _record_count # practically unused so far
                    if myoptions.debug > 1:
                        print("Debug14: Added INSerted codon at %d is %s, "
                              "reference codon is %s, amplicon region length %s "
                              "(nt), sliced [%d:%d] (pythonic slice numbering), "
                              "aa_residue is %s" % (
                                  _current_codon_position, '---',
                                  _reference_codon, _amplicon_length,
                                  _zero_based_codon_startpos + _previous_gaps,
                                  _zero_based_codon_startpos + _previous_gaps + 3 + _new_gaps_in_reference,
                                  _new_aa_residue))
                else:
                    _sample_codon_depadded = _rough_sample_codon.replace('-', '')
                    if myoptions.debug > 1:
                        print("Debug15: Rough sample codon at %d is %s, "
                              "reference codon is %s, amplicon region length %s "
                              "(nt), sliced [%d:%d] (pythonic slice numbering), "
                              "_sample_codon_depadded is %s, "
                              "_reference_codon_depadded is %s, "
                              "_sample_codon_contained_pad=%s, "
                              "_reference_codon_contained_pad=%s" % (
                                  _current_codon_position, _rough_sample_codon,
                                  _reference_codon, _amplicon_length,
                                  _zero_based_codon_startpos + _previous_gaps,
                                  _zero_based_codon_startpos + _previous_gaps + 3 + _new_gaps_in_reference,
                                  _sample_codon_depadded, _reference_codon_depadded,
                                  str(_sample_codon_contained_pad),
                                  str(_reference_codon_contained_pad)))

                    if (_sample_codon_contained_pad and _reference_codon_contained_pad) and not _sample_codon_depadded and not _reference_codon:
                        _unchanged_codons['---'] += _record_count
                        _unchanged_aa_residues['-'] += _record_count
                        if myoptions.debug > 1:
                            print("Debug16: Skipping unchanged codon at %d due "
                                  "to a padding dash in both %s and %s" % (
                                      _current_codon_position,
                                      _rough_sample_codon, _reference_codon))
                    elif _rough_sample_codon == '---' and _reference_codon == '---':
                        _unchanged_codons['---'] += _record_count
                        _unchanged_aa_residues['-'] += _record_count
                        if myoptions.debug > 1:
                            print("Debug17: Skipping unchanged codon at %d due "
                                  "to a padding dash in both sample and "
                                  "reference codon %s" % (
                                      _current_codon_position, _rough_sample_codon))
                    elif (_sample_codon_contained_pad or _reference_codon_contained_pad) and _sample_codon_depadded == _reference_codon_depadded:
                        _unchanged_codons[_sample_codon_depadded] += _record_count
                        _unchanged_aa_residues[alt_translate(_sample_codon_depadded)] += _record_count
                        if myoptions.debug > 1:
                            print("Debug18: Skipping unchanged codon at %d due "
                                  "to a padding dash in both %s and %s" % (
                                      _current_codon_position,
                                      _rough_sample_codon, _reference_codon))
                    elif _sample_codon_depadded == _reference_codon_depadded:
                        _unchanged_codons[_sample_codon_depadded] += _record_count
                        _unchanged_aa_residues[alt_translate(_sample_codon_depadded)] += _record_count
                        if myoptions.debug > 1:
                            print("Debug19: Unchanged codon at %d is %s" % (
                                _current_codon_position, _sample_codon_depadded))
                    elif _sample_codon_depadded != _reference_codon_depadded:
                        _changed_codons[_sample_codon_depadded] += _record_count
                        _new_aa_residue = alt_translate(_rough_sample_codon)
                        if 'N' in _sample_codon_depadded and _new_aa_residue != 'X':
                            if _sample_codon_depadded not in ('TCN', 'CTN', 'GTN', 'CCN', 'ACN', 'GCN', 'CGN', 'GGN'):
                                raise ValueError(
                                    "Error: biopython translated codon %s at "
                                    "%d into %s" % (
                                        _sample_codon_depadded,
                                        _current_codon_position,
                                        _new_aa_residue))
                            else:
                                if myoptions.debug > 1:
                                    print("Debug20: biopython translated codon "
                                          "%s at %d into %s" % (
                                              _sample_codon_depadded,
                                              _current_codon_position,
                                              _new_aa_residue))
                        if _new_aa_residue:
                            if _new_aa_residue != _reference_aa:
                                _changed_aa_residues[_new_aa_residue] += _record_count
                            else:
                                _unchanged_aa_residues[_new_aa_residue] += _record_count
                        else:
                            _changed_aa_residues['None'] += _record_count
                            if myoptions.debug > 1:
                                print("Debug21: Added None to "
                                      "_changed_aa_residues for codon at %s" % (
                                          _current_codon_position))
                        if myoptions.debug:
                            print("Debug22: Final codon at %d is %s, reference "
                                  "codon is %s, amplicon region length %s (nt), "
                                  "sliced [%d:%d] (pythonic slice numbering), "
                                  "aa_residue is %s" % (
                                      _current_codon_position, _rough_sample_codon,
                                      _reference_codon, _amplicon_length,
                                      _zero_based_codon_startpos + _previous_gaps,
                                      _zero_based_codon_startpos + _previous_gaps + 3 + _new_gaps_in_reference,
                                      _new_aa_residue))
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
                        if myoptions.debug > 1:
                            print("Debug23: Including frame-breaking codon %s "
                                  "at %d into per-codon coverage, reference "
                                  "codon is %s, amplicon region length %s (nt), "
                                  "sliced [%d:%d] (pythonic slice numbering), "
                                  "aa_residue is %s" % (
                                      _rough_sample_codon,
                                      _current_codon_position, _reference_codon,
                                      _amplicon_length,
                                      _zero_based_codon_startpos + _previous_gaps,
                                      _zero_based_codon_startpos + _previous_gaps + 3 + _new_gaps_in_reference,
                                      _new_aa_residue))
            else:
                if myoptions.debug:
                    print("Debug10b:Here7. Leaking sample codon %s at position "
                          "%d, _amplicon_length=%d" % (
                              _rough_sample_codon, _current_codon_position,
                              _amplicon_length))

        # all entries processed for the current codon column
        if myoptions.debug:
            print("Debug24: Lists at %d of changed codons and aminoacid "
                  "residues: %s-%s %s %s" % (
                      _current_codon_position,
                      _zero_based_padded_reference_aa_index + _zero_based_codon_startpos + 1,
                      _zero_based_padded_reference_aa_index + _zero_based_codon_startpos + 3,
                      _changed_codons, _changed_aa_residues))
            print("Debug25: Lists at %d of unchanged codons and aminoacid "
                  "residues: %s-%s %s %s" % (
                      _current_codon_position,
                      _zero_based_padded_reference_aa_index + _zero_based_codon_startpos + 1,
                      _zero_based_padded_reference_aa_index + _zero_based_codon_startpos + 3,
                      _unchanged_codons, _unchanged_aa_residues))
            print("Debug26a: DELeted codons at %d: %s" % (
                _current_codon_position, _deleted_reference_codons))
            print("Debug26b: INSerted codons at %d: %s" % (
                _current_codon_position, _inserted_codons))

        _total_codons_per_site_counts = (
            _inserted_codons + _deleted_reference_codons
            + _changed_codons + _unchanged_codons
        )
        _total_codons_per_site_sum = sum(_total_codons_per_site_counts.values())

        if myoptions.debug:
            print("Debug26c: _reference_codon was %s" % str(_reference_codon))

        if len(_reference_protein_seq) > _zero_based_padded_reference_aa_index:
            _reference_aa = _reference_protein_seq[_zero_based_padded_reference_aa_index]
        else:
            sys.stderr.write(
                "Error: Probably _zero_based_padded_reference_aa_index=%s at "
                "the end of the padded reference sequence which was %s long\n" % (
                    _zero_based_padded_reference_aa_index,
                    len(_reference_protein_seq)))
            sys.stderr.flush()

        if _zero_based_codon_startpos not in _already_checked_starts and len(_reference_codon) == 3:
            if _reference_codon == '---':
                _reference_aa = '-'
                _current_aa = '-'
            else:
                try:
                    _current_aa = alt_translate(_reference_codon)
                except Exception:
                    _current_aa = 'Biopython failure'
            if _current_aa != _reference_aa:
                if not min_start:
                    raise ValueError(
                        "Error: The reference aa at %d should be '%s' but "
                        "parsed reference codon '%s' encodes '%s'. The "
                        "_zero_based_codon_startpos = %s" % (
                            _current_codon_position, _reference_aa,
                            _reference_codon, _current_aa,
                            _zero_based_codon_startpos))
            else:
                _already_checked_starts.append(_zero_based_codon_startpos)

        if myoptions.debug and myoptions.debug > 1:
            for _key in _changed_codons:
                print("Debug27: {}: {} = {:8.6f} {} {}".format(
                    _key, _changed_codons[_key],
                    _changed_codons[_key] / _total_codons_per_site_sum,
                    set([x for x in _changed_codons]),
                    set([x for x in _changed_aa_residues])))

        if myoptions.debug:
            print("Debug27a: _new_gaps_in_reference=%s, "
                  "_reference_codon_contained_pad=%s, "
                  "_sample_codon_contained_pad=%s, _reference_codon=%s, "
                  "_rough_sample_codon=%s" % (
                      _new_gaps_in_reference, _reference_codon_contained_pad,
                      _sample_codon_contained_pad, _reference_codon,
                      _rough_sample_codon))

        if _new_gaps_in_reference:
            raise ValueError(
                "Error: _new_gaps_in_reference should be zero but is %s "
                "instead" % _new_gaps_in_reference)
        _new_gaps_in_reference = _reference_codon.count('-')

        _natural_codon_position_padded = (
            _zero_based_padded_reference_aa_index + 1
            + int(myoptions.left_reference_offset / 3.0)
            + aa_start
        )
        _natural_codon_position_depadded = (
            _natural_codon_position_padded
            - int((_previous_gaps + _new_gaps_in_reference) / 3.0)
        )

        write_tsv_line(
            outfilename_unchanged_codons, _unchanged_codons,
            _natural_codon_position_padded, _natural_codon_position_depadded,
            _reference_aa, _total_codons_per_site_sum, _reference_codon,
            debug=myoptions.debug,
        )

        if len(_changed_codons):
            write_tsv_line(
                outfilename, _changed_codons,
                _natural_codon_position_padded, _natural_codon_position_depadded,
                _reference_aa, _total_codons_per_site_sum, _reference_codon,
                debug=myoptions.debug,
            )

        if len(_inserted_codons): # possibly use a different output file for INSertion so that they do not collide with changed codons and add multiple lines
            write_tsv_line(
                outfilename, _inserted_codons,
                _natural_codon_position_padded, _natural_codon_position_depadded,
                'INS', _total_codons_per_site_sum, _reference_codon,
                debug=myoptions.debug,
            )

        # we add the DEL lines into the main frequencies.tsv file here even without calling write_tsv_line()
        # the write_tsv_line() accepts list of codons in the sample but here we have a 
        if _is_deletion:
            for _some_deleted_codon in _deleted_reference_codons:
                _observed_codon_count = Decimal(
                    _deleted_reference_codons[_some_deleted_codon]
                )
                _observed_codon_count2 = _observed_codon_count if _observed_codon_count else 0
                outfilename.write(
                    "{}\t{}\t{}\t{}\t{:8.6f}\t{}\t{}\t{}\t{}\n".format(
                        _natural_codon_position_padded,
                        _natural_codon_position_depadded,
                        _reference_aa, 'DEL',
                        Decimal(_deleted_reference_codons[_some_deleted_codon])
                        / Decimal(_total_codons_per_site_sum),
                        _some_deleted_codon, '---',
                        _observed_codon_count2, _total_codons_per_site_sum,
                    )
                )
                if myoptions.debug:
                    print("TESTING2:\t{}\t{}\t{}\t{}\t{:8.6f}\t{}\t{}\t{}\t{}".format(
                        _natural_codon_position_padded,
                        _natural_codon_position_depadded,
                        _reference_aa, 'DEL',
                        Decimal(_deleted_reference_codons[_some_deleted_codon])
                        / Decimal(_total_codons_per_site_sum),
                        _some_deleted_codon, '---',
                        _observed_codon_count2, _total_codons_per_site_sum,
                    ))
        outfilename.flush()

        _previous_gaps += _new_gaps_in_reference
        _new_gaps_in_reference = 0

        try:
            _top_most_codon, _top_most_count = _total_codons_per_site_counts.most_common()[0]
        except IndexError:
            pass
        else:
            _top_most_codons.append(_top_most_codon)

        if len(_reference_protein_seq) > _zero_based_padded_reference_aa_index + 1:
            if _is_insertion:
                _zero_based_padded_reference_aa_index += 1
                if myoptions.debug:
                    print("Debug33: Increased _zero_based_padded_reference_aa_index "
                          "to %d, moving to next codon" % _zero_based_padded_reference_aa_index)
            elif _is_deletion:
                _zero_based_padded_reference_aa_index += 1
                if myoptions.debug:
                    print("Debug31: Increased _zero_based_padded_reference_aa_index "
                          "to %d, moving to next codon" % _zero_based_padded_reference_aa_index)
            else:
                _zero_based_padded_reference_aa_index += 1
                if myoptions.debug:
                    print("Debug32: Increased _zero_based_padded_reference_aa_index "
                          "to %d, moving to next codon" % _zero_based_padded_reference_aa_index)
            if myoptions.debug:
                print("Debug33: After increments: "
                      "_natural_codon_position_depadded=%d "
                      "_previous_gaps=%d _new_gaps_in_reference=%d" % (
                          _natural_codon_position_depadded,
                          _previous_gaps, _new_gaps_in_reference))
        else:
            break

    del _align
    alnfilename_count.write("%s\n" % _total_aln_entries_used)
    alnfilename_count.close()
    _consensus = ''.join(_top_most_codons).upper()
    print("Info: consensus = %s" % str(_consensus))
    if _consensus in _padded_reference_dna_seq.upper():
        print("Info: Sample consensus sequence should roughly match substring "
              "inside the reference %s and IT DOES: %s" % (
                  str(_padded_reference_dna_seq), str(_top_most_codons)))
    else:
        print("Info: Sample consensus sequence should roughly match substring "
              "inside the reference %s BUT IT DOES NOT, maybe due to some true "
              "major mutations in some codons: %s" % (
                  str(_padded_reference_dna_seq), str(_top_most_codons)))


def open_file(outfilename):
    """Open a new file for writing, raising an error if it already exists."""
    if os.path.exists(outfilename):
        raise RuntimeError(
            "The file %s already exists, will not overwrite it." % outfilename
        )
    return open(outfilename, 'x')
