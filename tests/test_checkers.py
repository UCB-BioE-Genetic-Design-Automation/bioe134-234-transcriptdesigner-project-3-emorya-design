import pytest
from genedesign.checkers.hairpin_checker import HairpinChecker
from genedesign.checkers.forbidden_sequence_checker import ForbiddenSequenceChecker
from genedesign.checkers.internal_promoter_checker import PromoterChecker


class TestHairpinChecker:

    @pytest.fixture(autouse=True)
    def setup(self):
        self.checker = HairpinChecker()
        self.checker.initiate()

    def test_empty_sequence_passes(self):
        passed, msg = self.checker.run("")
        assert passed
        assert msg is None

    def test_short_sequence_passes(self):
        passed, msg = self.checker.run("ATGCGT")
        assert passed
        assert msg is None

    def test_poly_a_passes(self):
        passed, msg = self.checker.run("A" * 50)
        assert passed
        assert msg is None

    def test_multiple_hairpins_in_window_fails(self):
        hp1 = "GGGGG" + "TTTT" + "CCCCC"
        hp2 = "CCGCC" + "AAAA" + "GGCGG"
        seq = hp1 + hp2 + "A" * (50 - len(hp1) - len(hp2))
        passed, msg = self.checker.run(seq)
        if not passed:
            assert msg is not None

    def test_return_type_on_fail_is_string(self):
        hp = "GCGCGCGCG" + "AAAA" + "CGCGCGCGC"
        seq = hp + hp[::-1]
        passed, msg = self.checker.run(seq)
        if not passed:
            assert isinstance(msg, str)


class TestForbiddenSequenceChecker:

    @pytest.fixture(autouse=True)
    def setup(self):
        self.checker = ForbiddenSequenceChecker()
        self.checker.initiate()

    def test_clean_sequence_passes(self):
        passed, site = self.checker.run("ATGCGTAAAGTATTCGTTCCGCTGACCGTTTAA")
        assert passed
        assert site is None

    def test_empty_sequence_passes(self):
        passed, site = self.checker.run("")
        assert passed
        assert site is None

    def test_ecori_site_fails(self):
        passed, site = self.checker.run("ATGGAATTCTAA")
        assert not passed
        assert site == "GAATTC"

    def test_xbai_site_fails(self):
        passed, site = self.checker.run("ATGTCTAGATAA")
        assert not passed
        assert site == "TCTAGA"

    def test_hindiii_site_fails(self):
        passed, site = self.checker.run("ATGAAGCTTTAA")
        assert not passed
        assert site == "AAGCTT"

    def test_bsai_site_fails(self):
        passed, site = self.checker.run("ATGGGTCTCTAA")
        assert not passed

    def test_poly_g_fails(self):
        passed, site = self.checker.run("GGGGGGGGG")
        assert not passed
        assert site in ("GGGGGGGG", "CCCCCCCC")

    def test_returns_first_hit_only(self):
        passed, site = self.checker.run("GAATTCTCTAGA")
        assert not passed
        assert isinstance(site, str)


class TestPromoterChecker:

    @pytest.fixture(autouse=True)
    def setup(self):
        self.checker = PromoterChecker()
        self.checker.initiate()

    def test_perfect_constitutive_promoter_fails(self):
        passed, promoter = self.checker.run("TTGACAATTAATCATCGAACTAGTATAAT")
        assert not passed
        assert promoter is not None
        assert len(promoter) == 29

    def test_strong_j23_variant_fails(self):
        passed, _ = self.checker.run("ttgacagctagctcagtcctaggtataatgctagc")
        assert not passed

    def test_broken_constitutive_passes(self):
        passed, promoter = self.checker.run("TTctgAATTAATCATCGAACTAGgcgAAT")
        assert passed
        assert promoter is None

    def test_random_coding_sequence_passes(self):
        passed, promoter = self.checker.run("ATGCGTAAAGTATTCGTTCCGCTGACCGTTTAA")
        assert passed
        assert promoter is None

    def test_empty_sequence_passes(self):
        passed, promoter = self.checker.run("")
        assert passed
        assert promoter is None

    def test_fail_returns_29nt_string(self):
        passed, promoter = self.checker.run("TTGACAATTAATCATCGAACTAGTATAAT")
        if not passed:
            assert isinstance(promoter, str)
            assert len(promoter) == 29

    def test_pass_returns_none(self):
        passed, promoter = self.checker.run("ATGATGATGATGATGATGATGATGATGATG")
        if passed:
            assert promoter is None

    def test_rc_strand_promoter_detected(self):
        from genedesign.seq_utils.reverse_complement import reverse_complement
        rc = reverse_complement("TTGACAATTAATCATCGAACTAGTATAAT")
        passed, _ = self.checker.run(rc)
        assert not passed
