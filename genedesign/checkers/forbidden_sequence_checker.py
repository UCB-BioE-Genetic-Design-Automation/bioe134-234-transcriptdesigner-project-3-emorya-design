from genedesign.seq_utils.reverse_complement import reverse_complement


class ForbiddenSequenceChecker:

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
            # Asymmetric sites + their explicit reverse complements
            "GAGACC",    # BsaI RC
            "GAGACG",    # BsmBI RC
            "GCAGGTG",   # AarI RC
        ]

    def run(self, dnaseq):
        dnaseq = dnaseq.upper()
        rc = reverse_complement(dnaseq)
        for site in self.forbidden:
            if site in dnaseq:
                return False, site
            if site in rc:
                return False, f"{site} (reverse strand)" 
        return True, None
