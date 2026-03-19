import random
import os
from genedesign.rbs_chooser import RBSChooser
from genedesign.models.transcript import Transcript
from genedesign.checkers.forbidden_sequence_checker import ForbiddenSequenceChecker
from genedesign.checkers.hairpin_checker import hairpin_checker
from genedesign.checkers.internal_promoter_checker import PromoterChecker
from genedesign.checkers.rna_interference_checker import RNAInterferenceChecker
from genedesign.checkers.codon_checker import CodonChecker
from genedesign.seq_utils.translate import Translate


class TranscriptDesigner:
    """
    Reverse translates a protein sequence into a DNA sequence using
    CAI-weighted random codon selection, validating against forbidden sequences.
    """

    def __init__(self):
        self.aminoAcidToCodon = {}
        self.aa_to_codons = {}  # {amino_acid: [(codon, weight), ...]}
        self.rbsChooser = None
        self.forbidden_checker = None
        self.translator = None

    def initiate(self) -> None:
        self.rbsChooser = RBSChooser()
        self.rbsChooser.initiate()

        # Checks for restriction enzyme sites (BsaI, BbsI, etc.) that would be cut during cloning
        self.forbidden_checker = ForbiddenSequenceChecker()
        self.forbidden_checker.initiate()

        # Detects internal -10/-35 promoter-like sequences that could cause spurious transcription
        self.promoter_checker = PromoterChecker()
        self.promoter_checker.initiate()

        # Checks for reverse-complement duplexes with native E. coli mRNAs that trigger RNase III degradation
        self.rna_interference_checker = RNAInterferenceChecker()
        self.rna_interference_checker.initiate({})  # starts empty, builds up as genes are designed

        # Validates CAI, codon diversity, and absence of rare codons that stall translation due to low tRNA availability
        self.codon_checker = CodonChecker()
        self.codon_checker.initiate()

        # Used to verify the CDS back-translates to the original protein sequence
        self.translator = Translate()
        self.translator.initiate()

        # Codon table with highest CAI codon for each amino acid (for E. coli)
        self.aminoAcidToCodon = {
            'A': "GCG", 'C': "TGC", 'D': "GAT", 'E': "GAA", 'F': "TTC",
            'G': "GGT", 'H': "CAC", 'I': "ATC", 'K': "AAA", 'L': "CTG",
            'M': "ATG", 'N': "AAC", 'P': "CCG", 'Q': "CAG", 'R': "CGT",
            'S': "TCT", 'T': "ACC", 'V': "GTT", 'W': "TGG", 'Y': "TAC"
        }

        # Build weighted codon table from codon_usage.txt for guided random (Monte Carlo) selection
        self.aa_to_codons = {}
        data_path = os.path.join(os.path.dirname(__file__), 'data', 'codon_usage.txt')
        with open(data_path, 'r') as f:
            for line in f:
                parts = line.split()
                if len(parts) < 3:
                    continue
                codon = parts[0]
                aa = parts[1]
                freq = float(parts[2])
                if aa == '*':
                    continue
                if aa not in self.aa_to_codons:
                    self.aa_to_codons[aa] = []
                self.aa_to_codons[aa].append((codon, freq))

    def guided_random_codon(self, amino_acid: str) -> str:
        """Monte Carlo codon selection weighted by CAI frequencies from the host organism's codon usage table."""
        options = self.aa_to_codons[amino_acid]
        codons = [c for c, _ in options]
        weights = [w for _, w in options]
        return random.choices(codons, weights=weights, k=1)[0]

    def reverse_translate(self, peptide: str) -> str:
        """Reverse translates a protein sequence into a CDS using guided random codons."""
        codons = [self.guided_random_codon(aa) for aa in peptide]
        codons.append("TAA")
        return ''.join(codons)

    def check_window(self, cds_fragment):
        """Checks a CDS fragment for forbidden sequences, hairpins, and promoters."""
        if not self.forbidden_checker.run(cds_fragment)[0]:
            return False
        if not hairpin_checker(cds_fragment)[0]:
            return False
        if not self.promoter_checker.run(cds_fragment)[0]:
            return False
        return True

    def run(self, peptide: str, ignores: set) -> Transcript:
        # Every ORF must begin with ATG; M (methionine) is the only amino acid that encodes ATG
        if peptide[0] != 'M':
            raise ValueError(f"Peptide must start with M (methionine), got '{peptide[0]}'")

        max_global = 100   # Number of full re-rolls before giving up
        max_window = 50    # Max attempts to fix each failing window via synonymous re-randomization
        window_size = 17   # Codons per window (~51bp, covers hairpin and promoter check windows)

        for _ in range(max_global):
            # Start with a guided random encoding of the full protein
            codons = [self.guided_random_codon(aa) for aa in peptide]
            codons.append("TAA")  # stop codon

            # Sliding window: scan the CDS and re-randomize failing regions with synonymous codons
            fixed = True
            for i in range(len(codons) - window_size + 1):
                window_dna = ''.join(codons[i:i + window_size])
                if self.check_window(window_dna):
                    continue

                # Window failed — try re-randomizing just these codons (same amino acids, different codons)
                window_fixed = False
                for _ in range(max_window):
                    for j in range(i, min(i + window_size, len(peptide))):
                        codons[j] = self.guided_random_codon(peptide[j])
                    window_dna = ''.join(codons[i:i + window_size])
                    if self.check_window(window_dna):
                        window_fixed = True
                        break

                if not window_fixed:
                    fixed = False
                    break

            if not fixed:
                continue

            cds = ''.join(codons)

            # Check RNA interference against previously designed genes
            if not self.rna_interference_checker.run(cds)[0]:
                continue

            # Sanity check: confirm the generated CDS translates back to the original protein
            retranslated = self.translator.run(cds)
            if retranslated != peptide:
                continue

            # All passed — register this CDS so future genes are checked for cross-talk with it
            self.rna_interference_checker.native_genes[peptide[:10]] = cds

            selectedRBS = self.rbsChooser.run(cds, ignores)
            return Transcript(selectedRBS, peptide, codons)

        raise ValueError(f"Failed to generate valid CDS for peptide after {max_global} global attempts")


if __name__ == "__main__":
    # Longer test peptide (47 aa) to stress-test checkers on a realistically-sized protein
    peptide = "MYPFIRTARMTVKDELSQGHNWCPVAILFYEKDGSTRMQNHWCPVAL"

    designer = TranscriptDesigner()
    designer.initiate()
    ignores = set()

    # Show guided random codon function producing different codons
    print("Guided random codons for 'A':")
    for _ in range(5):
        print(f"  {designer.guided_random_codon('A')}")

    # Show reverse_translate producing different CDS each time
    print(f"\nReverse translate '{peptide}':")
    for _ in range(3):
        print(f"  {designer.reverse_translate(peptide)}")

    # Run full design with validation
    transcript = designer.run(peptide, ignores)
    cds = ''.join(transcript.codons)

    # Verify translation
    translator = Translate()
    translator.initiate()
    retranslated = translator.run(cds)
    print(f"\nCDS: {cds}")
    print(f"Peptide:      {peptide}")
    print(f"Retranslated: {retranslated}")
    print(f"Match: {peptide == retranslated}")
