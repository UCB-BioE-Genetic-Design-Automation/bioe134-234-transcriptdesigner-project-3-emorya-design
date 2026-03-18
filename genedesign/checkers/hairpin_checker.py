from genedesign.seq_utils.hairpin_counter import hairpin_counter


class HairpinChecker:
    CHUNK_SIZE = 50
    OVERLAP    = 25
    MIN_STEM   = 3
    MIN_LOOP   = 4
    MAX_LOOP   = 9

    def initiate(self):
        pass

    def run(self, dna):
        starts = list(range(0, max(1, n - self.CHUNK_SIZE + 1), self.OVERLAP))
    if n > self.CHUNK_SIZE and starts[-1] + self.CHUNK_SIZE < n:
        starts.append(n - self.CHUNK_SIZE)
    
    for i in starts:
        chunk = dna[i : i + self.CHUNK_SIZE]
        count, hairpin_str = hairpin_counter(
            chunk, self.MIN_STEM, self.MIN_LOOP, self.MAX_LOOP
        )
        if count > 0:
            return False, hairpin_str
    return True, None
def hairpin_checker(dna):
    checker = HairpinChecker()
    checker.initiate()
    return checker.run(dna)
