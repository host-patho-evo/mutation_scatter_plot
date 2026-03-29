# This work © 2025-2026 by Jiří Zahradník and Martin Mokrejš
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
import sys
import time
import typing

import multiprocessing
from collections import Counter
from decimal import Decimal

import numpy as np

from .. import alt_translate

VERSION = "0.3"

# Module-level dict populated by parse_alignment *before* forking the worker
# pool.  Forked workers inherit all entries via copy-on-write so tasks only
# need to send a single integer (_pos) over the IPC pipe.
_WORKER_SHARED: dict = {}

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
                print(f"Debug: Detected {seq.count('-')} minus signs in the sequence but after all "
                      "the nucleotide sequence can be divided by three when they are "
                      "omitted, good.")
        else:
            raise ValueError(
                f"Error: Sequence {seq} cannot be divided by 3 and removing minus "
                "signs does not help either"
            )
    return _codons


def write_tsv_line(outfilename, codons, natural_codon_position_padded,
                   natural_codon_position_depadded, reference_aa,
                   total_codons_per_site_sum, reference_codon, debug=False):
    """Write one or more TSV lines for all codons observed at a site."""
    if not outfilename:
        return
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
            f"{natural_codon_position_padded}\t{natural_codon_position_depadded}\t{reference_aa}\t{_some_aa}\t{Decimal(_observed_codon_count2) / Decimal(total_codons_per_site_sum):.6f}\t{reference_codon}\t{_some_codon}\t{_observed_codon_count2}\t{_total_codons_per_site_sum}\n"
        )
        if debug:
            print(f"TESTING1:\t{natural_codon_position_padded}\t{natural_codon_position_depadded}\t{reference_aa}\t{_some_aa}\t{Decimal(_observed_codon_count2) / Decimal(total_codons_per_site_sum):.6f}\t{reference_codon}\t{_some_codon}\t{_observed_codon_count2}\t{_total_codons_per_site_sum}")
    # Flush is done by the caller on a time-gate to avoid expensive NFS round-trips.


def _process_one_site(
    _pos, _num_unique, _alignment_len, _aln_array, _counts_array,
    _padded_len_array, _depadded_len_array, _leading_gap_array, _trailing_gap_array,
    _padded_reference_dna_seq, _reference_protein_seq, _ref_gaps_cumulative,
    aa_start, myoptions_minimum_aln_length, myoptions_left_reference_offset
):
    """
    Process a single codon position and return the aggregated results.

    This function implements Speedup 4 (NumPy Vectorization) by using vectorized
    column slicing and np.bincount to group unique sequences by their (codon, action)
    tuple. It is designed to be picklable for Speedup 5 (Multi-processing).
    """
    _zero_based_padded_reference_aa_index = int(_pos / 3)

    _reference_codon = _padded_reference_dna_seq[3*_zero_based_padded_reference_aa_index:3*_zero_based_padded_reference_aa_index + 3]
    _reference_codon_depadded = _reference_codon.replace('-', '')

    _reference_aa = _reference_protein_seq[_zero_based_padded_reference_aa_index] if len(_reference_protein_seq) > _zero_based_padded_reference_aa_index else '?'

    # Pass 1: Local Grouping (Vectorized with NumPy)
    _local_groups = {}
    if _num_unique > 0:
        _actions = np.zeros(_num_unique, dtype=np.uint8)
        if myoptions_minimum_aln_length:
            _actions[_depadded_len_array < myoptions_minimum_aln_length] = 1
        _actions[_depadded_len_array == 0] = 2

        _is_masked = _actions == 0
        if np.any(_is_masked):
            _lead = _leading_gap_array
            _trail = _trailing_gap_array
            _mask = _is_masked & (_lead != 0) & (_pos < _lead - 1)
            _actions[_mask] = 3
            _is_masked &= ~_mask
            _mask = _is_masked & (_trail != 0) & (_pos + 1 > _trail)
            _actions[_mask] = 4
            _is_masked &= ~_mask
            _mask = _is_masked & (_padded_len_array + 1 <= _pos + 3)
            _actions[_mask] = 5

        _chunks = _aln_array[:, _pos:_pos + 3]
        _packed = (_chunks[:, 0].astype(np.uint32) << 16 |
                   _chunks[:, 1].astype(np.uint32) << 8 |
                   _chunks[:, 2].astype(np.uint32))
        _combined = _packed | (_actions.astype(np.uint32) << 24)
        _unique_combined, _first_indices, _inverse = np.unique(_combined, return_index=True, return_inverse=True)
        _agg_counts = np.bincount(_inverse, weights=_counts_array)
        _appearance_order = np.argsort(_first_indices)
        _unique_combined = _unique_combined[_appearance_order]
        _agg_counts = _agg_counts[_appearance_order]

        _action_map = {0: "regular", 1: "min_len_fail", 2: "empty_seq",
                       3: "leading_gap", 4: "trailing_gap", 5: "too_short"}

        for _val, _count in zip(_unique_combined, _agg_counts):
            if _count == 0:
                continue
            _act_code = int(_val >> 24)
            _pck = _val & 0xFFFFFF
            _codon = bytes([(_pck >> 16) & 0xFF, (_pck >> 8) & 0xFF, _pck & 0xFF]).decode()
            _action_str = _action_map[_act_code]
            _local_groups[(_codon, _action_str)] = int(_count)

    _unchanged_codons = Counter()
    _changed_codons = Counter()
    _deleted_reference_codons = Counter()
    _inserted_codons = Counter()
    _is_deletion = False
    _is_insertion = False

    for (_rough_sample_codon, _action), _record_count in _local_groups.items():
        _sample_codon_contained_pad = '-' in _rough_sample_codon
        _sample_codon_depadded = _rough_sample_codon.replace('-', '') if _sample_codon_contained_pad else _rough_sample_codon

        if _action == "regular":
            if _rough_sample_codon == '---' and _reference_codon != '---' and _reference_codon != 'NNN':
                _is_deletion = True
                _deleted_reference_codons[_reference_codon] += _record_count
            elif _reference_codon == '---' and _rough_sample_codon != '---':
                _is_insertion = True
                _inserted_codons[_sample_codon_depadded] += _record_count
            elif _sample_codon_depadded == _reference_codon_depadded:
                _unchanged_codons[_sample_codon_depadded] += _record_count
            else:
                _changed_codons[_sample_codon_depadded] += _record_count

    _total_counts = _inserted_codons + _deleted_reference_codons + _changed_codons + _unchanged_codons
    _total_sum = sum(_total_counts.values())
    _new_gaps = _reference_codon.count('-')
    _prev_gaps = int(_ref_gaps_cumulative[_pos - 1]) if _pos > 0 else 0
    _nat_pos_padded = _zero_based_padded_reference_aa_index + 1 + int(myoptions_left_reference_offset / 3.0) + aa_start
    _nat_pos_depadded = _nat_pos_padded - int((_prev_gaps + _new_gaps) / 3.0)

    return {
        'pos': _pos,
        'nat_padded': _nat_pos_padded,
        'nat_depadded': _nat_pos_depadded,
        'ref_aa': _reference_aa,
        'ref_codon': _reference_codon,
        'total_sum': int(_total_sum),
        'unchanged': _unchanged_codons,
        'changed': _changed_codons,
        'inserted': _inserted_codons,
        'deleted': _deleted_reference_codons,
        'is_deletion': _is_deletion,
        'is_insertion': _is_insertion,
        'counts': _total_counts,
    }


def fast_fasta_iter(handle):
    """
    A lightweight FASTA iterator that yields (name, sequence, lineno) tuples.
    Significantly faster than Bio.SeqIO.parse for large files.

    *lineno* is the 1-based line number of the ``>`` header line so that
    callers can report the exact location of malformed records.

    Embedded carriage-returns (``\\r``, ``^M``) inside sequence lines are
    stripped so that files with mixed or Windows/old-Mac line endings do not
    produce sequences with unexpected lengths.
    """
    name, seq = None, []
    lineno = 0
    start_lineno = 0
    for line in handle:
        lineno += 1
        if not line.strip():
            continue
        if line[0] == ">":
            if name is not None:
                yield name, "".join(seq), start_lineno
            # Mimic Bio.SeqIO.parse: id is the first word after the '>'
            _parts = line[1:].split()
            name = _parts[0] if _parts else ""
            start_lineno = lineno
            seq = []
        else:
            # Strip trailing whitespace AND any embedded \r (^M) that appears
            # when files have Windows/old-Mac line endings within sequence data.
            seq.append(line.strip().replace('\r', ''))
    if name is not None:
        yield name, "".join(seq), start_lineno


def _process_one_site_wrapper(_pos: int):
    """Thin IPC wrapper used with the COW-fork optimisation.

    All heavy data (alignment array, reference sequences, etc.) lives in the
    module-level ``_WORKER_SHARED`` dict which forked workers inherit via
    copy-on-write.  Only the codon position index ``_pos`` travels over the
    IPC pipe — reducing per-task payload from ~115 MB to 4 bytes.
    """
    _d = _WORKER_SHARED
    return _process_one_site(
        _pos,
        _d['num_unique'],
        _d['alignment_len'],
        _d['aln_array'],
        _d['counts_array'],
        _d['padded_len_array'],
        _d['depadded_len_array'],
        _d['leading_gap_array'],
        _d['trailing_gap_array'],
        _d['padded_reference_dna_seq'],
        _d['reference_protein_seq'],
        _d['ref_gaps_cumulative'],
        _d['aa_start'],
        _d['minimum_aln_length'],
        _d['left_reference_offset'],
    )

def parse_alignment(myoptions: typing.Any, alignment_file: str, padded_reference_dna_seq: str,
                    reference_protein_seq: str, reference_as_codons: list[str],
                    outfilename: typing.Any, outfilename_unchanged_codons: typing.Any,
                    alnfilename_count: typing.Any, aa_start: int, min_start: int, max_stop: int,
                    threads: int = None,
                    pool: typing.Optional["multiprocessing.pool.Pool"] = None,
                    chunksize: typing.Optional[int] = None):
    """
    Parse a padded multi-FASTA alignment and write codon frequency TSV files.

    Optimizations:
    - Speedup 4: NumPy-based vectorized column slicing for O(1) site extraction.
    - Speedup 5: Multi-processing parallelization using multiprocessing.Pool.
    - Speedup 6: High-precision Decimal(str) formatting for bit-identical output.

    left_reference_offset and right_reference_offset are used to slice the
    reference. discard_this_many_leading_nucs and
    discard_this_many_trailing_nucs are used to discard offending
    leading/trailing nucleotides.

    If the reference protein is shorter than the sample entries, the trailing codons
    after reference protein terminated are treated as INSertions with aa position
    of the aminoacid residue corresponding to the last aminoacid codon of the
    reference protein. We had to switch to using padded-reference instead.

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
        print(f"Debug0: Depadded reference sequence has length {len(padded_reference_dna_seq.replace('-', ''))}, padded "
              f"reference sequence has length {len(padded_reference_dna_seq)}, each entry from {alignment_file} must also "
              f"have same padded length {len(padded_reference_dna_seq)}")
        print(f"Debug1: reference_protein_seq={str(reference_protein_seq)} with length {len(reference_protein_seq)}")
    if not os.path.exists(alignment_file):
        raise RuntimeError(f"Alignment file not found: {alignment_file}")
    if os.path.getsize(alignment_file) == 0:
        raise RuntimeError(f"Alignment file is empty: {alignment_file}")
    # We use SeqIO.parse for streaming grouping to handle huge files
    # while only keeping unique sequences in memory.

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
            print(f"Info: len(padded_reference_dna_seq)={len(padded_reference_dna_seq)}, "
                  f"myoptions.left_reference_offset={myoptions.left_reference_offset - 1}, "
                  f"myoptions.right_reference_offset={myoptions.right_reference_offset}")
            print(f"Info: After cutting input using offset position: {_padded_reference_dna_seq}")
            print(f"Info: After cutting input using offset position: {_reference_protein_seq}")
            print(f"Info: After cutting input using offset position: {_reference_as_codons}")
        if not _reference_protein_seq:
            raise ValueError(
                "Error: No _reference_protein_seq provided or left. Was the "
                f"slicing using myoptions.left_reference_offset-1={myoptions.left_reference_offset - 1}, "
                f"myoptions.right_reference_offset={myoptions.right_reference_offset} wrong?")
        if not _reference_as_codons:
            raise ValueError(
                "Error: No _reference_as_codons provided or left. Was the "
                f"slicing using myoptions.left_reference_offset-1={myoptions.left_reference_offset - 1}, "
                f"myoptions.right_reference_offset={myoptions.right_reference_offset} wrong?")
        if myoptions.debug:
            print(f"Debug2: Depadded reference sequence has length {len(_padded_reference_dna_seq.replace('-', ''))}, padded "
                  f"reference sequence has length {len(_padded_reference_dna_seq)}, each padded entry from {alignment_file} "
                  f"must also have same padded length")
    else:
        _padded_reference_dna_seq = padded_reference_dna_seq
        _reference_protein_seq = reference_protein_seq
        _reference_as_codons = reference_as_codons

    _padded_reference_dna_seq = _padded_reference_dna_seq.upper()
    _reference_protein_seq = _reference_protein_seq.upper()
    _top_most_codons = []
    _total_aln_entries_used: int = 0
    _start_from: int = int(myoptions.discard_this_many_leading_nucs) if myoptions.discard_this_many_leading_nucs else 0
    _stop_to: int = int(myoptions.discard_this_many_trailing_nucs) if myoptions.discard_this_many_trailing_nucs else 0

    _expected_aln_len = len(_padded_reference_dna_seq)
    _raw_groups = {}
    with open(alignment_file, "r", encoding="utf-8") as _handle:
        if myoptions.x_after_count:
            for _record_id, _seq_str, _lineno in fast_fasta_iter(_handle):
                # Parse count from ID prefix: format is NNNNx.SHA256HASH or NNNNx
                # e.g. '576521x.7cbee25f...' → count=576521
                # We extract digits before the first 'x'; fall back to 1 if absent.
                _x_pos = _record_id.find('x')
                if _x_pos > 0 and _record_id[:_x_pos].isdigit():
                    _record_count = int(_record_id[:_x_pos])
                else:
                    _record_count = 1

                if _seq_str in _raw_groups:
                    _raw_groups[_seq_str][0] += _record_count
                else:
                    _raw_groups[_seq_str] = [_record_count, _record_id, _lineno]
        else:
            for _, _seq_str, _lineno in fast_fasta_iter(_handle):
                if _seq_str in _raw_groups:
                    _raw_groups[_seq_str][0] += 1
                else:
                    _raw_groups[_seq_str] = [1, "", _lineno]

    _parsed_alignments = {}
    for _raw_seq, (_record_count, _record_id, _lineno) in _raw_groups.items():
        _total_aln_entries_used += _record_count

        # Validate sequence length early, while we still have the record ID and
        # line number available for a useful diagnostic message.
        _raw_len = len(_raw_seq)
        if _raw_len != _expected_aln_len and not (_start_from or _stop_to):
            print(
                f"Warning: skipping record at line {_lineno} "
                f"(ID: '{_record_id}', sequence length {_raw_len} "
                f"!= reference length {_expected_aln_len}). "
                f"Use: grep -n -A 1 '^>{_record_id}' '{alignment_file}'",
                file=sys.stderr,
            )
            _total_aln_entries_used -= _record_count
            continue

        # Slicing and upper-casing logic
        if _start_from or _stop_to:
            if _stop_to:
                _aln_line_seq = _raw_seq[_start_from:-_stop_to].upper()
            else:
                _aln_line_seq = _raw_seq[_start_from:].upper()
        else:
            _aln_line_seq = _raw_seq.upper()

        if _aln_line_seq in _parsed_alignments:
            _parsed_alignments[_aln_line_seq]['count'] += _record_count
        else:
            _padded_aln_line_length = len(_aln_line_seq)
            # str.count scans in C without allocating any temporary strings;
            # ~2× faster than .replace('-','').replace('N','') which creates
            # two 3822-char heap allocations per sequence (87k × 2 = 174k allocs).
            _depadded_aln_line_length = (
                _padded_aln_line_length
                - _aln_line_seq.count('-')
                - _aln_line_seq.count('N')
            )

            # lstrip/rstrip in pure C: ~200× faster than $-anchored regex on 3822-char strings.
            # _aln_line_seq is already .upper() so 'n' → 'N' — strip set covers both gap chars.
            _end_of_leading_gaps = len(_aln_line_seq) - len(_aln_line_seq.lstrip('-N'))
            _start_of_trailing_gaps = len(_aln_line_seq.rstrip('-N'))

            _parsed_alignments[_aln_line_seq] = {
                'count': _record_count,
                'seq': _aln_line_seq,
                'padded_len': _padded_aln_line_length,
                'depadded_len': _depadded_aln_line_length,
                'end_of_leading_gaps': _end_of_leading_gaps,
                'start_of_trailing_gaps': _start_of_trailing_gaps
            }

    _parsed_alignments_list = list(_parsed_alignments.values())
    _num_unique = len(_parsed_alignments_list)
    _alignment_len = len(_padded_reference_dna_seq)

    # Pre-allocate NumPy arrays for vectorized operations
    if _num_unique > 0:
        # Optimization: use np.frombuffer on a single string join instead of list comprehension
        _all_seqs_str = "".join(item['seq'] for item in _parsed_alignments_list)
        # Safety check: wrong-length sequences should have been caught and skipped
        # in the grouping loop above, but guard here in case slicing changes things.
        if len(_all_seqs_str) != _num_unique * _alignment_len:
            raise ValueError(
                f"Internal error: joined sequence buffer length {len(_all_seqs_str)} "
                f"!= {_num_unique} * {_alignment_len}. "
                "This indicates sequences with unexpected lengths slipped through filtering."
            )
        _aln_array = np.frombuffer(_all_seqs_str.encode('ascii'), dtype=np.uint8).reshape(_num_unique, _alignment_len)


        _counts_array = np.array([item['count'] for item in _parsed_alignments_list], dtype=np.int64)
        _padded_len_array = np.array([item['padded_len'] for item in _parsed_alignments_list], dtype=np.int32)
        _depadded_len_array = np.array([item['depadded_len'] for item in _parsed_alignments_list], dtype=np.int32)
        _leading_gap_array = np.array([item['end_of_leading_gaps'] for item in _parsed_alignments_list], dtype=np.int32)
        _trailing_gap_array = np.array([item['start_of_trailing_gaps'] for item in _parsed_alignments_list], dtype=np.int32)
    else:
        _aln_array = np.array([], dtype=np.uint8).reshape(0, _alignment_len)
        _counts_array = np.array([], dtype=np.int64)
        _padded_len_array = np.array([], dtype=np.int32)
        _depadded_len_array = np.array([], dtype=np.int32)
        _leading_gap_array = np.array([], dtype=np.int32)
        _trailing_gap_array = np.array([], dtype=np.int32)

    # Pre-calculate previous gaps in reference to make each loop iteration independent
    _ref_aln_array = np.frombuffer(_padded_reference_dna_seq.encode('ascii'), dtype=np.uint8)
    _ref_gaps_cumulative = np.cumsum(_ref_aln_array == ord('-'))

    # Pre-allocate counters once and clear() each iteration to avoid 8 * N_sites constructor calls
    _top_most_codons = []

    # ── COW-fork dispatch ───────────────────────────────────────────────────
    # Populate the module global *before* creating (forking) the Pool so that
    # worker processes inherit all large arrays via copy-on-write.  Each task
    # then only sends _pos (one int, 4 bytes) over the IPC pipe instead of
    # the full alignment array (~115 MB with 30k unique seqs × 3822 cols).
    _slim_positions = list(range(min_start, max_stop or _alignment_len, 3))

    if _slim_positions:
        global _WORKER_SHARED  # pylint: disable=global-statement
        _WORKER_SHARED = {
            'num_unique':             _num_unique,
            'alignment_len':          _alignment_len,
            'aln_array':              _aln_array,
            'counts_array':           _counts_array,
            'padded_len_array':       _padded_len_array,
            'depadded_len_array':     _depadded_len_array,
            'leading_gap_array':      _leading_gap_array,
            'trailing_gap_array':     _trailing_gap_array,
            'padded_reference_dna_seq': _padded_reference_dna_seq,
            'reference_protein_seq':  _reference_protein_seq,
            'ref_gaps_cumulative':    _ref_gaps_cumulative,
            'aa_start':               aa_start,
            'minimum_aln_length':     myoptions.minimum_aln_length,
            'left_reference_offset':  myoptions.left_reference_offset,
        }
        # Pool is created AFTER the global is set so forked workers see it.
        _external_pool = pool is not None
        _all_results = []
        if _external_pool:
            _active_pool = pool
            try:
                _all_results = _active_pool.starmap(
                    _process_one_site_wrapper,
                    [(_pos,) for _pos in _slim_positions],
                    chunksize=chunksize,
                )
            finally:
                _WORKER_SHARED = {}
        else:
            with multiprocessing.Pool(processes=threads) as _active_pool:
                try:
                    _all_results = _active_pool.starmap(
                        _process_one_site_wrapper,
                        [(_pos,) for _pos in _slim_positions],
                        chunksize=chunksize,
                    )
                finally:
                    _WORKER_SHARED = {}

        _last_flush_t = time.monotonic()
        for _res in _all_results:
            # Write results
            write_tsv_line(outfilename_unchanged_codons, _res['unchanged'], _res['nat_padded'], _res['nat_depadded'], _res['ref_aa'], _res['total_sum'], _res['ref_codon'], debug=myoptions.debug)
            if _res['changed']:
                write_tsv_line(outfilename, _res['changed'], _res['nat_padded'], _res['nat_depadded'], _res['ref_aa'], _res['total_sum'], _res['ref_codon'], debug=myoptions.debug)
            if _res['inserted']:
                write_tsv_line(outfilename, _res['inserted'], _res['nat_padded'], _res['nat_depadded'], 'INS', _res['total_sum'], _res['ref_codon'], debug=myoptions.debug)

            if _res['is_deletion']:
                for _some_deleted_codon in _res['deleted']:
                    _count = _res['deleted'][_some_deleted_codon]
                    outfilename.write(f"{_res['nat_padded']}\t{_res['nat_depadded']}\t{_res['ref_aa']}\tDEL\t{Decimal(_count) / Decimal(_res['total_sum']):.6f}\t{_some_deleted_codon}\t---\t{_count}\t{_res['total_sum']}\n")

            # Time-gated flush: avoids one NFS round-trip (~10 ms) per site.
            # On large datasets (hours runtime) this still flushes every 30 s so
            # output can be monitored with tail -f; on short benchmark runs it
            # flushes 0-1 times instead of once per site (was 3829x = ~38 s).
            _now = time.monotonic()
            if _now - _last_flush_t >= 30.0:
                if outfilename:
                    outfilename.flush()
                if outfilename_unchanged_codons:
                    outfilename_unchanged_codons.flush()
                _last_flush_t = _now

            # Update top most codons
            if _res['counts']:
                _top_most_codon, _count = _res['counts'].most_common(1)[0]
                _top_most_codons.append(_top_most_codon)

        # Final flush to ensure the last partial buffer is written.
        if outfilename:
            outfilename.flush()
        if outfilename_unchanged_codons:
            outfilename_unchanged_codons.flush()

    alnfilename_count.write(f"{_total_aln_entries_used}\n")
    alnfilename_count.close()
    _consensus = ''.join(_top_most_codons).upper()
    print(f"Info: consensus = {_consensus!s}")
    if _consensus in _padded_reference_dna_seq.upper():
        print("Info: Sample consensus sequence should roughly match substring "
              f"inside the reference {str(_padded_reference_dna_seq)} and IT DOES: {str(_top_most_codons)}")
    else:
        print("Info: Sample consensus sequence should roughly match substring "
              f"inside the reference {str(_padded_reference_dna_seq)} BUT IT DOES NOT, maybe due to some true "
              f"major mutations in some codons: {str(_top_most_codons)}")


def open_file(outfilename, overwrite=False, encoding=None):
    """Open a new file for writing, raising an error if it already exists and overwrite is False."""
    if os.path.exists(outfilename) and not overwrite:
        raise RuntimeError(
            f"The file {outfilename} already exists, will not overwrite it."
        )
    if overwrite and os.path.exists(outfilename):
        return open(outfilename, 'w', encoding=encoding)
    return open(outfilename, 'x', encoding=encoding)
