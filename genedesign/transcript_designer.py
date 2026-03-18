"""
TranscriptDesigner (improved)
==============================
Reverse-translates a protein sequence into a CDS that:
  1. Maximises Codon Adaptation Index (CAI) for E. coli K-12
  2. Minimises RNA hairpin count
  3. Contains no forbidden restriction-enzyme sites
  4. Contains no internal σ70 promoter-like elements

Strategy
--------
We use a sliding-window Monte-Carlo / greedy repair loop:

  1. Start from the highest-CAI codon at every position (greedy).
  2. Run all checkers on the full CDS.
  3. If any checker fails, attempt local codon swaps at the offending
     region (trying synonymous codons in descending CAI order) until all
     checkers pass or the iteration budget is exhausted.
  4. Return the best design found (fewest checker violations, then
     highest CAI).
"""

from __future__ import annotations

import random
from typing import Dict, List, Tuple

# These imports assume the standard genedesign package layout.
# Adjust if your repo uses a different namespace.
from genedesign.rbs_chooser import RBSChooser
from genedesign.models.transcript import Transcript

from genedesign.checkers.hairpin_checker import HairpinChecker
from genedesign.checkers.forbidden_sequence_checker import ForbiddenSequenceChecker
from genedesign.checkers.codon_checker import CodonChecker
from genedesign.checkers.internal_promoter_checker import PromoterChecker
from genedesign.seq_utils.reverse_complement import reverse_complement


# ---------------------------------------------------------------------------
# E. coli K-12 codon table
# Each amino acid maps to a list of (codon, relative_adaptiveness) tuples
# sorted by relative adaptiveness descending (w = f/f_max, 0 < w ≤ 1).
# Values derived from Carbone et al. / Codon Usage Database (kazusa.or.jp).
# ---------------------------------------------------------------------------

CODON_TABLE: Dict[str, List[Tuple[str, float]]] = {
    'A': [("GCG", 1.00), ("GCC", 0.76), ("GCA", 0.59), ("GCT", 0.58)],
    'C': [("TGC", 1.00), ("TGT", 0.70)],
    'D': [("GAT", 1.00), ("GAC", 0.62)],
    'E': [("GAA", 1.00), ("GAG", 0.52)],
    'F': [("TTC", 1.00), ("TTT", 0.74)],
    'G': [("GGT", 1.00), ("GGC", 0.72), ("GGA", 0.34), ("GGG", 0.30)],
    'H': [("CAC", 1.00), ("CAT", 0.78)],
    'I': [("ATC", 1.00), ("ATT", 0.86), ("ATA", 0.20)],
    'K': [("AAA", 1.00), ("AAG", 0.51)],
    'L': [("CTG", 1.00), ("TTA", 0.46), ("TTG", 0.43), ("CTT", 0.39),
          ("CTC", 0.37), ("CTA", 0.20)],
    'M': [("ATG", 1.00)],
    'N': [("AAC", 1.00), ("AAT", 0.88)],
    'P': [("CCG", 1.00), ("CCA", 0.54), ("CCT", 0.47), ("CCC", 0.38)],
    'Q': [("CAG", 1.00), ("CAA", 0.58)],
    'R': [("CGT", 1.00), ("CGC", 0.72), ("CGG", 0.42), ("CGA", 0.33),
          ("AGA", 0.32), ("AGG", 0.25)],
    'S': [("TCT", 1.00), ("TCC", 0.89), ("AGC", 0.87), ("TCA", 0.70),
          ("AGT", 0.65), ("TCG", 0.56)],
    'T': [("ACC", 1.00), ("ACA", 0.74), ("ACT", 0.73), ("ACG", 0.56)],
    'V': [("GTT", 1.00), ("GTC", 0.80), ("GTA", 0.48), ("GTG", 0.47)],
    'W': [("TGG", 1.00)],
    'Y': [("TAC", 1.00), ("TAT", 0.84)],
    '*': [("TAA", 1.00), ("TGA", 0.70), ("TAG", 0.35)],   # stop codons
}

# Pre-build reverse map: codon -> amino acid
_CODON_TO_AA: Dict[str, str] = {
    codon: aa
    for aa, pairs in CODON_TABLE.items()
    for codon, _ in pairs
}

MAX_ITER = 500   # repair iterations before giving up


class TranscriptDesigner:
    """
    Reverse translates a protein sequence into an optimised DNA CDS and
    selects a ribosome binding site.
    """

    def __init__(self):
        self.rbsChooser: RBSChooser = None
        self._hairpin   = HairpinChecker()
        self._codon = CodonChecker()
        self._forbidden = ForbiddenSequenceChecker()
        self._promoter  = PromoterChecker()

    # ------------------------------------------------------------------
    # Initialisation
    # ------------------------------------------------------------------

    def initiate(self) -> None:
        self.rbsChooser = RBSChooser()
        self.rbsChooser.initiate()
        self._hairpin.initiate()
        self._codon.initiate()
        self._forbidden.initiate()
        self._promoter.initiate()

    # ------------------------------------------------------------------
    # Main entry point
    # ------------------------------------------------------------------

    def run(self, peptide: str, ignores: set) -> Transcript:
        """
        Design a CDS for *peptide* that passes all quality checks.

        Parameters
        ----------
        peptide : str
            Amino acid sequence (single-letter codes, NO stop symbol).
        ignores : set
            RBS candidates to skip (passed through to RBSChooser).

        Returns
        -------
        Transcript
            The best design found within the iteration budget.
        """
        # Start with greedy highest-CAI assignment
        codons = self._greedy_codons(peptide)
        codons.append("TAA")           # stop codon

        best_codons = list(codons)
        best_violations = self._count_violations(codons)

        for _ in range(MAX_ITER):
            if best_violations == 0:
                break
            codons = self._targeted_repair(list(best_codons), peptide)
            v = self._count_violations(codons)
            if v <= best_violations:
                best_violations = v
                best_codons = list(codons)

        cds = "".join(best_codons)
        selected_rbs = self.rbsChooser.run(cds, ignores)
        return Transcript(selected_rbs, peptide, best_codons)

    # ------------------------------------------------------------------
    # Internal helpers
    # ------------------------------------------------------------------

    def _greedy_codons(self, peptide: str) -> List[str]:
        """Pick the highest-CAI codon for each amino acid."""
        return [CODON_TABLE[aa][0][0] for aa in peptide]

    def _cds(self, codons: List[str]) -> str:
        return "".join(codons)

    def _count_violations(self, codons: List[str]) -> int:
        cds = self._cds(codons)
        total = 0
        for checker in (self._hairpin, self._forbidden, self._promoter):
            passed, _ = checker.run(cds)
            if not passed:
                total += 1
        passed_codon, _, _, _ = self._codon.run(codons) 
        if not passed_codon: 
          total += 1 
        return total

    def _swap_codon(self, codons: List[str], idx: int, peptide: str, prefer_first: bool = False) -> List[str]:
        """Replace codon at idx with a different synonymous codon (not the current one)."""
        aa = peptide[idx]
        synonyms = CODON_TABLE[aa]
        if len(synonyms) == 1:
          return codons
        current = codons[idx]
        alternatives = [(c, w) for c, w in synonyms if c != current]
        if not alternatives:
          return codons
        if prefer_first:
          codons[idx] = alternatives[0][0]  # highest CAI non-current
        else:
          weights = [w for _, w in alternatives]
          codons[idx] = random.choices([c for c, _ in alternatives], weights=weights)[0]
        return codons

    def _targeted_repair(self, codons: List[str], peptide: str) -> List[str]:
        """
        Check each issue and attempt a targeted fix:
        - For forbidden/promoter: find the offending region in the CDS, swap
          a codon overlapping that region.
        - For hairpin: swap a random codon in the first half of the sequence
          (hairpins tend to form in GC-rich regions).
        - For codon diversity: swap the most-repeated codon to a synonym.
        - Falls back to random swap if nothing targeted applies.
        """
        cds = self._cds(codons)
        n = len(peptide)

        # --- Forbidden sequence ---
        passed_f, site_f = self._forbidden.run(cds)
        if not passed_f and site_f:
            clean_site = site_f.replace(" (reverse strand)", "")
            pos = cds.find(clean_site)
            if pos < 0:
                pos = cds.find(reverse_complement(clean_site))
            if pos >= 0:
                codon_idx = min(pos // 3, n - 1)
                for offset in range(-2, 5):
                    idx = max(0, min(codon_idx + offset, n - 1))
                    if len(CODON_TABLE[peptide[idx]]) > 1:
                      return self._swap_codon(codons, idx, peptide, prefer_first=True)
        # --- Internal promoter ---
        passed_p, site_p = self._promoter.run(cds)
        if not passed_p and site_p:
          anchor = site_p[:6]
          pos = cds.find(anchor)
          if pos < 0:
            pos = cds.find(reverse_complement(anchor))
          if pos >= 0:
            codon_idx = min(pos // 3, n - 1)
            for offset in range(-2, 5):
                idx = max(0, min(codon_idx + offset, n - 1))
                if len(CODON_TABLE[peptide[idx]]) > 1:
                  return self._swap_codon(codons, idx, peptide, prefer_first=True)


        # --- Hairpin: swap in a window around a random position ---
        passed_h, _ = self._hairpin.run(cds)
        if not passed_h:
            # Try swapping codons in a random window of 10
            start = random.randint(0, max(0, n - 10))
            for idx in range(start, min(start + 10, n)):
                if len(CODON_TABLE[peptide[idx]]) > 1:
                    return self._swap_codon(codons, idx, peptide)

        # --- Codon diversity: swap the most duplicated codon ---
        passed_c, _, _, _ = self._codon.run(codons)
        if not passed_c:
            from collections import Counter
            counts = Counter(codons[:n])  # exclude stop
            most_common_codon, _ = counts.most_common(1)[0]
            # Find positions where this codon is used and swap one
            positions = [i for i, c in enumerate(codons[:n]) if c == most_common_codon]
            if positions:
                idx = random.choice(positions)
                aa = peptide[idx]
                if len(CODON_TABLE[aa]) > 1:
                    return self._swap_codon(codons, idx, peptide)

        # --- Fallback: random swap ---
        idx = random.randrange(n)
        return self._swap_codon(codons, idx, peptide)

    # ------------------------------------------------------------------
    # Utility: compute CAI for a list of codons
    # ------------------------------------------------------------------

    @staticmethod
    def compute_cai(codons: List[str]) -> float:
        """
        Return the geometric-mean CAI for *codons* (stop codon excluded).
        """
        import math
        weights = []
        for codon in codons:
            aa = _CODON_TO_AA.get(codon.upper())
            if aa is None or aa == "*":
                continue
            for c, w in CODON_TABLE[aa]:
                if c == codon.upper():
                    weights.append(w)
                    break
        if not weights:
            return 0.0
        return math.exp(sum(math.log(w) for w in weights) / len(weights))


# ---------------------------------------------------------------------------
# Quick smoke-test
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    peptide  = "MYPFIRTARMTV"
    designer = TranscriptDesigner()
    designer.initiate()
    transcript = designer.run(peptide, set())
    print(transcript)
    cai = TranscriptDesigner.compute_cai(transcript.codons[:-1])  # exclude stop
    print(f"CAI: {cai:.4f}")
