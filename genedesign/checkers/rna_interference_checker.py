from genedesign.seq_utils.reverse_complement import reverse_complement


class RNAInterferenceChecker:
    """
    Checks for potential duplex formation between the designed RNA and native
    host RNA sequences. Reverse-complementary regions can cause expression
    interference by forming stable duplexes that trigger degradation or block
    translation.

    Returns True if no problematic duplexes are found.
    Returns False, the problematic sequence, and the gene name if a duplex is detected.
    """

    def __init__(self):
        self.native_genes = {}  # {gene_name: cds_sequence}
        self.min_duplex_length = 20  # minimum bp for stable duplex formation

    def initiate(self, native_genes: dict) -> None:
        """
        Initializes the checker with native gene CDS sequences.

        Parameters:
            native_genes (dict): {gene_name: cds_sequence} for the host organism.
        """
        self.native_genes = {name: seq.upper() for name, seq in native_genes.items()}

    def run(self, designed_seq: str) -> tuple:
        """
        Checks if the designed sequence could form duplexes with any native gene.

        Parameters:
            designed_seq (str): The designed DNA/CDS sequence.

        Returns:
            tuple: (bool, str or None, str or None)
                - (True, None, None) if no duplex is predicted.
                - (False, problematic_sequence, gene_name) if a duplex is detected.
        """
        designed_seq = designed_seq.upper()
        rc_designed = reverse_complement(designed_seq)

        # Slide a window of min_duplex_length across the reverse complement
        # and check if it appears in any native gene
        for i in range(len(rc_designed) - self.min_duplex_length + 1):
            window = rc_designed[i:i + self.min_duplex_length]

            for gene_name, native_seq in self.native_genes.items():
                if window in native_seq:
                    return False, window, gene_name

        return True, None, None


if __name__ == "__main__":
    checker = RNAInterferenceChecker()

    # Example native genes
    native = {
        "geneA": "ATGAAAGCGATCGATCGATCGATCGATCGATCGATCGTAA",
        "geneB": "ATGCCCGGGAAATTTCCCGGGAAATTTCCCGGGAAATAA",
    }
    checker.initiate(native)

    # Sequence that is reverse complement of part of geneA
    rc_of_geneA = reverse_complement("GCGATCGATCGATCGATCGATCG")
    print(f"Testing RC of geneA region: {rc_of_geneA}")
    result = checker.run(rc_of_geneA)
    print(f"Result: passed={result[0]}, seq={result[1]}, gene={result[2]}")

    # Sequence with no match
    print("\nTesting unrelated sequence:")
    result = checker.run("ATGTACTACTACTACTACTACTACTACTACTAC")
    print(f"Result: passed={result[0]}, seq={result[1]}, gene={result[2]}")
