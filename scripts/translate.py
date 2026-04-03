#! /usr/bin/env python3

# >A00808:1538:HWL7GDRX3:2:2101:7491:1125 minus
# AGAAAGTACTACTAATCTGTATGGTAGGTGACCAACACCATAAGTGGGTCGGAAACCATATGATTGTAAAGGAGAGTAACAATTAGGACCTGCAACACAAATACAAGGTTTGATACCGGCCAGATAGATATCAGTTGAAATATATCTCTCAAAAGGTTTCAGTTTAGACTTCCTAAACAATATATACCAGTAATCATAAATACCACTATGCTTAGAATCAAGCATGTTAGAATTCCAAGCTATATCGCAGCCTTTTAAATCATCTGGTAATTTTTAATTATAATCAGCAATATTTCCAGTTTGCCCTGGAGCGATTTGGCTGACTTCATTACCTTTAATTACAAATGAATCTGCATAGACATGAGTAAAGCAGTGATCATTTATTTTAGTAGGAGACACTCC
# >A00808:1538:HWL7GDRX3:2:2101:7419:1157 minus
# AGAAAGTACTACTACTCTGTATGGTTGGTGACCAACACCATAAGTGGGTCGGAAACCATATGATTGTAAAGGAGAGTAACAATTAGGACCTGCAACACCATTACAAGGTTTGTTACCGGCCTGATAGATTTCAGTTGAAATATCTCTCTCAAAAGGGTTGAGTTTAGACTTCCTAAACAATCTATACCAGTAATCATAATTACCACTATGCTTAGAATCAAGCTTGTTAGAATTCCAAGCTATAACGCAGCCTGTAAAATCATCTGGTAATTTATAATTATAATCAGCAATATTTCCAGTTTGCCCTGGAGCGATTTGGCTGACTTCATTACCTTTAATTACAAATGAATCTGCATAGACATTAGTAAAGCAGAGATCATTTAATTTAGTAGGAGACACTCC
# >A00808:1538:HWL7GDRX3:2:2101:21296:1204 plus
# GGAGTGTCTCCTACTAAATTAAATGATCTCTGCTTTACTAATGTCTATGCAGATTCATTTGTAATTAAAGGTAATGAAGTCAGCCAAATCGCTCCAGGGCAAACTGGAAATATTGCTGATTATAATTATAAATTACCAGATGATTTTACAGGCTGCGTTATAGCTTGGAATTCTAACAAGCTTGATTCTAAGCATAGTGGTAATTATGATTACTGGTATAGATTGTTTAGGAAGTCTAAACTCAAACCTTTTGAGAGAGATATTTCAACTGAAATCTATCAGGCCGGTAACAAACCTTGTAAAGG---TAAAGGTCCTAATTGTTACTTTCCTTTACAATCATATGGTTTCCGACCCACTTATGGTGTTGGTCACCAACCATACAGAGTAGTAGTACTTTCT
# >A00808:1538:HWL7GDRX3:2:2101:24831:1219 minus
# AGAAAGTACTACTACTCTGTATGGTTGGTGACCAACACCATAAGTGGGTCGGAAACCATATGATTGTAAAGGAAAGTAACAATTAGGACCTTTA---CCTTTACAAGGTTTGTTACCGGCCTGATAGATTTCAGTTGAAATATCTCTCTCAAAAGGTTTGAGTTTAGACTTCCTAAACAATCTATACCAGTAATCATAATTACCACTATGCTTAGAATCAAGCTTGTTAGAATTCCAAGCTATAACGCAGCCTGTAAAATCATCTGGTAATTTATAATTATAATCAGCAATATTTCCAGTTTGCCCTGGAGCGATTTGGCTGACTTCATTACCTTTAATTACAAATGAATCTGCATAGACATTAGTAAAGCAGAGATCATTTAATTTAGTAGGAGACACTCC
# >A00808:1538:HWL7GDRX3:2:2101:25021:1235 minus
# AGAAAGTACTACTACTCTGTATGGTTGGTGACCAACACCATAAGTGGGTCGGAAACCATATGATTGTAAAGGAAAGTAACAATTAGGACCTTTA---CCTTTACAAGGTTTGTTACCGGCCTGATAGATTTCAGTTGAAATATCTCTCTCAAAAGGTTTGAGTTTAGACTTCCTAAACAATCTATACCAGTAATCATAATTACCACTATGCTTAGAATCAAGCTTGTTAGAATTCCAAGCTATAACGCAGCCTGTAAAATCATCTGGTAATTTATAATTATAATCAGCAATATTTCCAGTTTGCCCTGGAGCGATTTGGCTGACTTCATTACCTTTAATTACAAATGAATCTGCATAGACATTAGTAAAGCAGAGATCATTTAATTTAGTAGGAGACACTCC

# testcases:
# /auto/vestec1-elixir/projects/biocev/mmokrejs/proj/zahradnik/sekvenace_wastewater_202403/demultiplex/processing_in_4_steps_by_barcode_and_primer/merged_resulting_files/split_by_flowcell_and_lane/WW_waste-water_samples__Sample17.VSCHT@JZ268_WW_F.GGCCTGCTACGTCGA.HWL7GDRX3.2.fastp.amplicons.clean.faa
# /auto/vestec1-elixir/projects/biocev/mmokrejs/proj/zahradnik/sekvenace_wastewater_202403/demultiplex/processing_in_4_steps_by_barcode_and_primer/merged_resulting_files/split_by_flowcell_and_lane/WW_waste-water_samples__RBDopure.IFQ491_1C1A.3@JZ271_WW_F.TTACCCCAGCCAGTC.HYKFKDRX3.1.fastp.amplicons.clean.faa
# /auto/vestec1-elixir/projects/biocev/mmokrejs/proj/zahradnik/sekvenace_wastewater_202403/demultiplex/processing_in_4_steps_by_barcode_and_primer/merged_resulting_files/split_by_flowcell_and_lane/WW_waste-water_samples__RBDopure.IFQ491_1C1A.3@JZ271_WW_F.TTACCCCAGCCAGTC.HWL7GDRX3.2.fastp.amplicons.clean.faa
# /auto/vestec1-elixir/projects/biocev/mmokrejs/proj/zahradnik/sekvenace_wastewater_202403/demultiplex/processing_in_4_steps_by_barcode_and_primer/merged_resulting_files/split_by_flowcell_and_lane/WW_waste-water_samples__WW-P1.parental@JZ240_50_F.GTGGGTTGGAGAT.HWLFTDRX3.2.fastp.amplicons.clean.faa
# /auto/vestec1-elixir/projects/biocev/mmokrejs/proj/zahradnik/sekvenace_wastewater_202403/demultiplex/processing_in_4_steps_by_barcode_and_primer/merged_resulting_files/split_by_flowcell_and_lane/WW_waste-water_samples__RBDopure.IFQ491_1C1A.3@JZ271_WW_F.TTACCCCAGCCAGTC.HWLFTDRX3.2.fastp.amplicons.clean.faa
# /auto/vestec1-elixir/projects/biocev/mmokrejs/proj/zahradnik/sekvenace_wastewater_202403/demultiplex/processing_in_4_steps_by_barcode_and_primer/merged_resulting_files/split_by_flowcell_and_lane/WW_waste-water_samples__Sample19.VSCHT@JZ270_WW_F.GATTGGTGGTTGCTA.HWL7GDRX3.2.fastp.amplicons.clean.faa
# /auto/vestec1-elixir/projects/biocev/mmokrejs/proj/zahradnik/sekvenace_wastewater_202403/demultiplex/processing_in_4_steps_by_barcode_and_primer/merged_resulting_files/split_by_flowcell_and_lane/WW_waste-water_samples__WW-A1.affinity@JZ190_0_F.ATCGATGTAATAA.HWLFTDRX3.2.fastp.amplicons.clean.faa

"""Example
# Note: the unique deltaV483 deletion aka CKG-KG
# Note: re-align the deletion around BA.2.86 region using sed -e 's/AAGG---TAAA/AAGGT---AAA/' included in /auto/vestec1-elixir/projects/biocev/mmokrejs/proj/zahradnik/Izrael/bin/fix_SARS-CoV2_S-protein_indel_misalignments.sh

gzip -dc /auto/vestec1-elixir/projects/biocev/mmokrejs/proj/zahradnik/sekvenace_wastewater_202403/demultiplex/processing_in_4_steps_by_barcode_and_primer/merged_resulting_files/split_by_flowcell_and_lane/WW_waste-water_samples__WW-P1.parental@JZ244_54_F.CATGCTCAATAGC.HWL7GDRX3.2.fastp.fastq.gz | fastq_to_fasta | blastn -task blastn -reward 2 -max_hsps 1 -num_alignments 1 -word_size 4 -num_threads 64 -dust no -evalue 1e-5 -outfmt "6 qacc sstart send sstrand qseq sseq" -db MN908947.3_S.fasta -query - 2>/dev/null | awk '{print $1,$2,$3,$4,$5,$6}' | python /auto/vestec1-elixir/projects/biocev/mmokrejs/proj/zahradnik/Izrael/bin/drop_erroneous_insertions.py --infile=- outfile=- | python /auto/vestec1-elixir/projects/biocev/mmokrejs/proj/zahradnik/Izrael/bin/reversecomplement_reads_on_minus.py --infile=- | /auto/vestec1-elixir/projects/biocev/mmokrejs/proj/zahradnik/Izrael/bin/fix_SARS-CoV2_S-protein_indel_misalignments.sh | python  /auto/vestec1-elixir/projects/biocev/mmokrejs/proj/zahradnik/Izrael/bin/translate.py --infile=- --respect-alignment | grep -v '^>' | sort | uniq -c | sort -nr
"""

import sys
import os
import io
from optparse import OptionParser

from Bio.Seq import translate as _bio_translate

VERSION = "202504031900"

myparser = OptionParser(version=f"%prog version {VERSION}")
myparser.add_option("--infile", action="store", type="string", dest="infile", default='',
    help="Input FASTA file with word plus or minus as the second word of the header, use minus (dash) for stdin")
myparser.add_option("--outfile", action="store", type="string", dest="outfile", default='',
    help="Output filename with words plus and minus removed from FASTA header lines, use minus (dash) for stdout")
myparser.add_option("--infileformat", action="store", type="string", dest="infileformat", default='fasta',
    help="Input file format (Default: fasta)")
myparser.add_option("--outfileformat", action="store", type="string", dest="outfileformat", default='fasta-2line',
    help="Output file format (Default: fasta-2line)")
myparser.add_option("--ignore-gaps", action="store_true", dest="ignore_gaps", default=False,
    help="Ignore eventual gaps in sequences and drop them before translation. If you padded the sequence alignment to keep it in-frame do not enable this.")
myparser.add_option("--respect-alignment", action="store_true", dest="respect_alignment", default=False,
    help="Respect padded alignment as input and do not die with an incomplete codon 'AA-' but return 'X'")
(myoptions, myargs) = myparser.parse_args()


# ── codon lookup table ────────────────────────────────────────────────────────

def _build_codon_table() -> dict:
    """Build a flat codon→amino-acid dict from the NCBI standard genetic code.

    Special cases handled up front so the hot translation loop needs only a
    single dict.get(codon, 'X') call per codon — no branching, no Biopython
    function-call overhead, no generator chaining.

    Codon resolution priority (applied once at build time):
      1. '---' → '-'           (alignment gap codon)
      2. codon with '-' but not '---':
           respect_alignment → 'X'  (partial gap in aligned input)
           ignore_gaps       → handled by stripping '-' before slicing
      3. Codon with 'N' → 'X'  (ambiguous nucleotide)
      4. Standard 64 codons    → looked up from Biopython CodonTable once
    """
    from Bio.Data import CodonTable  # pylint: disable=import-outside-toplevel
    std = CodonTable.unambiguous_dna_by_id[1]

    table: dict[str, str] = {}

    # Standard 64 sense codons (upper-case ACGTacgt)
    for codon, aa in std.forward_table.items():
        table[codon] = aa
        table[codon.lower()] = aa
    # Stop codons → '*'
    for codon in std.stop_codons:
        table[codon] = '*'
        table[codon.lower()] = '*'

    # Gap codon
    table['---'] = '-'

    # Any codon containing '-' (but not '---') → 'X' when respecting alignment
    # (These are built lazily in _translate_seq via the dict.get default.)

    return table


# ── raw-bytes streaming FASTA parser ─────────────────────────────────────────

def _iter_fasta_raw(fh):
    """Yield (header_str, seq_str) pairs from a binary file object.

    Avoids all SeqRecord allocation.  Sequences are concatenated from raw
    lines decoded as ASCII (non-ASCII bytes silently dropped — FASTA sequence
    data is always ASCII).
    """
    hdr: str | None = None
    parts: list[str] = []
    for raw in fh:
        if raw.startswith(b'>'):
            if hdr is not None:
                yield hdr, ''.join(parts)
            hdr = raw[1:].rstrip(b'\r\n').decode('utf-8', errors='replace')
            parts = []
        elif hdr is not None:
            parts.append(raw.rstrip(b'\r\n').decode('ascii', errors='ignore'))
    if hdr is not None:
        yield hdr, ''.join(parts)


# ── translation ───────────────────────────────────────────────────────────────

def _translate_seq(seq: str, table: dict, ignore_gaps: bool) -> str:
    """Translate *seq* to amino acids using the pre-built *table*.

    Three codon resolution paths, fastest first:

    1. Fast dict lookup (64 standard codons + '---'):
       Standard ATGC codons and the full-gap codon are resolved via a
       single dict.get() with no Biopython overhead.

    2. Partial-gap codon ('TC-', 'AT-', 'A-C', …) or pure IUPAC ('TCN'):
       Each '-' is treated as 'N' (unknown nucleotide — same semantics as
       an alignment gap), then Bio.Seq.translate() resolves per-codon.
       This correctly gives TC- → TCN → S and AT- → ATN → X (ambiguous).

    When *ignore_gaps* is True, all '-' are stripped before slicing so the
    reading frame is preserved across gapped regions.

    Any trailing incomplete codon is silently dropped.
    """
    if ignore_gaps:
        seq = seq.replace('-', '')
    length = len(seq)
    result: list[str] = []
    for i in range(0, length - length % 3, 3):
        codon = seq[i:i + 3]
        aa = table.get(codon)
        if aa is None:
            # Not in the fast dict: partial-gap (TC-) or IUPAC (TCN).
            # Treat '-' as 'N' — a gap in an alignment column means
            # "unknown nucleotide", the same semantics as N.  This lets
            # Biopython resolve the codon correctly per-codon:
            #   TC- → TCN → S  (all four nucleotides give Serine)
            #   AT- → ATN → X  (ATN can be Ile or Met → ambiguous)
            codon_n = codon.replace('-', 'N')
            try:
                aa = _bio_translate(codon_n, gap='-')
            except Exception:  # pylint: disable=broad-except
                aa = 'X'
        result.append(aa)
    return ''.join(result)


# ── main I/O loop ─────────────────────────────────────────────────────────────

def parse_input(
        infile,
        source_name: str,
        outfileh,
        ignore_gaps: bool,
        respect_alignment: bool,  # pylint: disable=unused-argument
) -> None:
    """Stream *infile* (binary), translate each record, write to *outfileh* (binary).

    *respect_alignment* is accepted for API compatibility.  The lookup-table
    implementation already handles partial-gap codons (any codon containing
    '-' that is not '---') by returning 'X' via the dict default, so the
    flag has no effect on the output.
    """
    table = _build_codon_table()
    for header, seq in _iter_fasta_raw(infile):
        try:
            aa_seq = _translate_seq(seq, table, ignore_gaps)
        except Exception as exc:
            raise ValueError(
                f"Cannot translate record '{header.split()[0]}' in {source_name}"
            ) from exc
        outfileh.write(b'>')
        outfileh.write(header.encode('utf-8', errors='replace'))
        outfileh.write(b'\n')
        outfileh.write(aa_seq.encode('ascii'))
        outfileh.write(b'\n')


# ── entry point ───────────────────────────────────────────────────────────────

def _main() -> None:
    """Parse CLI options and run the translation pipeline."""
    # ── open input ────────────────────────────────────────────────────────────
    if not myoptions.infile or myoptions.infile == '-':
        infileh = sys.stdin.buffer
        source_name = 'STDIN'
    elif os.path.exists(myoptions.infile):
        infileh = open(myoptions.infile, 'rb')  # pylint: disable=consider-using-with
        source_name = myoptions.infile
    else:
        raise RuntimeError(f"File {myoptions.infile} does not exist")

    # ── open output ───────────────────────────────────────────────────────────
    if not myoptions.outfile or myoptions.outfile == '-':
        outfileh: io.RawIOBase = sys.stdout.buffer
    elif os.path.exists(myoptions.outfile):
        raise RuntimeError(f"Error: File {myoptions.outfile} already exists")
    else:
        outfileh = open(myoptions.outfile, 'xb')  # pylint: disable=consider-using-with

    # Wrap in a large BufferedWriter to amortise write() syscall overhead.
    buf_out = io.BufferedWriter(outfileh, buffer_size=4 << 20)  # 4 MB

    try:
        parse_input(
            infileh,
            source_name,
            buf_out,
            ignore_gaps=myoptions.ignore_gaps,
            respect_alignment=myoptions.respect_alignment,
        )
    finally:
        buf_out.flush()
        if infileh is not sys.stdin.buffer:
            infileh.close()
        if outfileh is not sys.stdout.buffer:
            outfileh.close()


if __name__ == "__main__":
    _main()

