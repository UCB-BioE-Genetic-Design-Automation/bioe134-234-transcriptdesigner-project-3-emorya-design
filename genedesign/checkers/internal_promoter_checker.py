import math
from genedesign.seq_utils.reverse_complement import reverse_complement


class PromoterChecker:

    def initiate(self):
        pfm = [
            [0, 0, 0, 12, 0, 12, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 0, 12, 0, 12, 12, 0],
            [0, 0, 0,  0,12,  0, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 0,  0, 0,  0,  0, 0],
            [0, 0,12,  0, 0,  0, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 0,  0, 0,  0,  0, 0],
            [12,12, 0,  0, 0,  0, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3,12,  0,12,  0,  0,12],
        ]
        ncols = len(pfm[0])
        self.pwm = [[0] * ncols for _ in range(4)]
        for x in range(ncols):
            total = sum(pfm[y][x] for y in range(4))
            for y in range(4):
                freq = pfm[y][x]
                prob_base = 0.25
                w = (math.log(
                    (freq + math.sqrt(total) * prob_base) /
                    (total + math.sqrt(total)) / prob_base
                )) / math.log(2)
                self.pwm[y][x] = w

    def run(self, seq):
        seq = seq.upper()
        rc = reverse_complement(seq)
        combined = seq + "x" + rc
        sliding_frame = 29
        threshold = 9.134
        for i in range(len(combined) - sliding_frame + 1):
            score = 0.0
            partseq = combined[i:i + sliding_frame]
            for x, base in enumerate(partseq):
                y = {'A': 0, 'C': 1, 'G': 2, 'T': 3}.get(base, -1)
                if y != -1:
                    score += self.pwm[y][x]
            if score >= threshold:
                return False, partseq
        return True, None
