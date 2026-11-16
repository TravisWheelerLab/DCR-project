#!/usr/bin/env python3
"""esl_translate_replica.py

Pure-Python replica of `esl-translate -l 1 --watson` (Easel / HMMER miniapp,
Eddy/Rivas Lab), used in the DCR conservation pipeline to enumerate top-strand
ORFs. Validated to reproduce the exact ORF set of the pipeline's shipped
esl-translate output for all 1599 mouse-ortholog DCR sequences (100% match).

Semantics (verified):
  - Top strand only (--watson): no reverse complement.
  - Standard genetic code (NCBI transl_table 1).
  - ORFs are maximal stop-to-stop runs in each of the 3 forward frames; NO
    ATG-start requirement (esl-translate default). A leading partial ORF (frame
    start -> first stop) and a trailing partial ORF (last stop -> seq end) are
    both emitted.
  - -l 1: minimum ORF length 1 amino acid.
  - Coords are 1-based inclusive nt positions of the CODING region (stop codon
    excluded), on the input (top) strand.

Purpose     stand-in for the esl-translate binary, which could not be installed/built on the
            analysis machine; used by the revision controls (GC-matched controls, self-shuffle,
            latent potential, definition reconciliation, X/X' filter -- see ../archive_prior_session/)
            and re-implemented inline in selection/selection_v3.py (esl_orfs).
Inputs      FASTA (CLI) or a sequence string (esl_translate_watson()).
Outputs     CLI: TSV on stdout -- name, start1, end1, frame (1-3), peptide. This is a flat table,
            not esl-translate's FASTA output format.
Usage       python esl_translate_replica.py seqs.fa [--min-aa 1] > orfs.tsv
            from esl_translate_replica import esl_translate_watson
Manuscript  Methods (ORF enumeration, Figure 2B(i)-(ii) controls).
"""

_BASES = "TCAG"
_CODONS = [a+b+c for a in _BASES for b in _BASES for c in _BASES]
_AAS = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
CODON_TABLE = dict(zip(_CODONS, _AAS))
STOPS = {"TAA", "TAG", "TGA"}


def translate_codon(c):
    return CODON_TABLE.get(c.upper(), "X")


def esl_translate_watson(seq, min_aa=1):
    """Return list of (start1, end1, frame, pep) for top-strand ORFs.
    start1/end1 are 1-based inclusive nt coords of the coding region."""
    s = seq.upper()
    n = len(s)
    orfs = []
    for frame in (1, 2, 3):
        off = frame - 1
        orf_start = off
        pep = []
        i = off
        while i + 3 <= n:
            codon = s[i:i+3]
            if codon in STOPS:
                if pep:
                    orfs.append((orf_start + 1, i, frame, "".join(pep)))
                pep = []
                orf_start = i + 3
            else:
                pep.append(translate_codon(codon))
            i += 3
        if pep:
            orfs.append((orf_start + 1, orf_start + 3*len(pep), frame, "".join(pep)))
    return [o for o in orfs if len(o[3]) >= min_aa]


def read_fasta(path):
    d = {}; name = None; buf = []
    with open(path) as f:
        for line in f:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if name is not None:
                    d[name] = "".join(buf)
                name = line[1:].split()[0]; buf = []
            else:
                buf.append(line)
        if name is not None:
            d[name] = "".join(buf)
    return d


if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser(description="pure-python replica of esl-translate -l 1 --watson")
    ap.add_argument("fasta")
    ap.add_argument("--min-aa", type=int, default=1)
    a = ap.parse_args()
    for name, seq in read_fasta(a.fasta).items():
        for st, en, fr, pep in esl_translate_watson(seq, a.min_aa):
            print(f"{name}\t{st}\t{en}\t{fr}\t{pep}")
