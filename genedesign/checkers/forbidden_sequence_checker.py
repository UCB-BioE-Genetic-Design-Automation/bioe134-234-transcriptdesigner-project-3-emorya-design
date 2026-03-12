from genedesign.seq_utils.reverse_complement import reverse_complement


class ForbiddenSequenceChecker:
    def __init__(self):
        self.forbidden = []

    def initiate(self):
        self.forbidden = [
            "AAAAAAAA",  # poly(A)
            "TTTTTTTT",  # poly(T)
            "CCCCCCCC",  # poly(C)
            "GGGGGGGG",  # poly(G)
            "ATATATAT",  # poly(AT)
            "CAATTG",    # MfeI
            "GAATTC",    # EcoRI
            "GGATCC",    # BamHI
            "AGATCT",    # BglII
            "ACTAGT",    # SpeI
            "TCTAGA",    # XbaI
            "GGTCTC",    # BsaI
            "CGTCTC",    # BsmBI
            "CACCTGC",   # AarI
            "CTGCAG",    # PstI
            "CTCGAG",    # XhoI
            "GCGGCCGC",  # NotI
            "AAGCTT",    # HindIII
        ]

    def run(self, dnaseq):
        rc = reverse_complement(dnaseq)
        combined = (dnaseq + "x" + rc).upper()
        for site in self.forbidden:
            if site in combined:
                return False, site
        return True, None
