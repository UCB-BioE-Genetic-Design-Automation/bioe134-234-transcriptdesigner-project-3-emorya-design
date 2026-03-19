import pytest
from genedesign.checkers.rna_interference_checker import RNAInterferenceChecker
from genedesign.seq_utils.reverse_complement import reverse_complement


@pytest.fixture
def checker():
    c = RNAInterferenceChecker()
    native_genes = {
        "dnaA": "ATGTCACTTTCGCTTTGGCAGCAGTGTCTTGCCCGATTGCAGGATGAGTTACCAGCCGGATATT",
        "rpoB": "ATGGTTTACTCCTATACCGAGAAAGCGACATTTGACAAACGTCTGATTGCTCAGGTCATGGAAGA",
        "ftsZ": "ATGTTTGAACCAATGGAACTTACCAATGACGCGATCGTCAATATTATTGGTGTCGGTATCGCAC",
    }
    c.initiate(native_genes)
    return c


def test_no_duplex_passes(checker):
    """A sequence with no RC match to native genes should pass."""
    seq = "ATGTACTACTACTACTACTACTACTACTACTACTACTACTAC"
    passed, prob_seq, gene = checker.run(seq)
    assert passed is True
    assert prob_seq is None
    assert gene is None


def test_duplex_with_native_gene_fails(checker):
    """A sequence that is reverse complement of a native gene region should fail."""
    # Take a 25bp chunk from dnaA and reverse complement it
    dnaa_region = "ATGTCACTTTCGCTTTGGCAGCAGTG"
    rc_region = reverse_complement(dnaa_region)
    passed, prob_seq, gene = checker.run(rc_region)
    assert passed is False
    assert gene == "dnaA"
    assert prob_seq is not None


def test_short_match_passes(checker):
    """A match shorter than min_duplex_length should pass."""
    # Take only 10bp from dnaA — below 20bp threshold
    short_region = "ATGTCACTTT"
    rc_region = reverse_complement(short_region)
    passed, prob_seq, gene = checker.run(rc_region)
    assert passed is True


def test_exact_threshold_length_fails(checker):
    """A match exactly at min_duplex_length (20bp) should fail."""
    # Take exactly 20bp from rpoB
    rpob_20bp = "ATGGTTTACTCCTATACCGA"
    rc_region = reverse_complement(rpob_20bp)
    passed, prob_seq, gene = checker.run(rc_region)
    assert passed is False
    assert gene == "rpoB"


def test_one_below_threshold_passes(checker):
    """A match at 19bp (one below threshold) should pass."""
    rpob_19bp = "ATGGTTTACTCCTATACCG"
    rc_region = reverse_complement(rpob_19bp)
    passed, prob_seq, gene = checker.run(rc_region)
    assert passed is True


def test_empty_sequence_passes(checker):
    """An empty designed sequence should pass (no duplex possible)."""
    passed, prob_seq, gene = checker.run("")
    assert passed is True


def test_no_native_genes_passes():
    """With no native genes loaded, any sequence should pass."""
    c = RNAInterferenceChecker()
    c.initiate({})
    passed, prob_seq, gene = c.run("ATGAAAGCGATCGATCGATCGATCGATCGATCG")
    assert passed is True


def test_identifies_correct_gene(checker):
    """Should report the specific gene that the duplex matches."""
    # RC of a region from ftsZ
    ftsz_region = "ATGTTTGAACCAATGGAACTTACC"
    rc_region = reverse_complement(ftsz_region)
    passed, prob_seq, gene = checker.run(rc_region)
    assert passed is False
    assert gene == "ftsZ"
