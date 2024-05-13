# -*- coding: utf-8 -*-
import os
import unittest

import pytest
from support import CACHE

import hgvs
from hgvs.pretty_print import PrettyPrint


@pytest.mark.quick
@pytest.mark.models
class Test_SimplePosition(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.hp = hgvs.parser.Parser()
        cls.hdp = hgvs.dataproviders.uta.connect(
            mode=os.environ.get("HGVS_CACHE_MODE", "run"), cache=CACHE
        )
        cls.pp = PrettyPrint(
            cls.hdp,
        )
        cls.pp.useColor = False

        cls.atta_expected_results = (
            "          :   123,346,500         123,346,520         123,346,540\n"
            + "chrom pos :   |    .    |    .    |    .    |    .    |    .\n"
            + "seq    -> : ATAAAGCTTTTCCAAATGTTATTAATTACTGGCATTGCTTTTTGCCAA\n"
            + "region    :                     |------|                    \n"
            + "aa seq <- :                      TerAsnSerAlaAsnSerLysAlaLeu\n"
            + "tx seq <- : TATTTCGAAAAGGTTTACAATAATTAATGACCGTAACGAAAAACGGTT\n"
            + "tx pos    :  |    .    |    .   |   |    .    |    .    |   \n"
            + "          :  *20       *10      *1  2880      2870      2860\n"
            + "ref>alt   : ATTAATTA>ATTAATTAATTA\n"
        )

    def test_var_c1_forward(self):
        """test c1 on -> strand"""

        hgvs_c = "NM_198689.2:c.1="
        var_c = self.hp.parse(hgvs_c)
        result = self.pp.display(var_c)
        print(result)
        result = result.split("\n")
        expected_str = (
            "hgvs      : NC_000021.8:g.46020522=\n"
            + "hgvs      : NM_198689.2:c.1=\n"
            + "          :         46,020,510          46,020,530\n"
            + "chrom pos :    .    |    .    |    .    |    .    |  \n"
            + "seq    -> : CCTCCAGTTCAATCCCCAGCATGGCCGCGTCCACTATGTCT\n"
            + "tx ref dif:                          X               \n"
            + "region    :                     =                    \n"
            + "aa seq -> :                     MetAlaAlaSerThrMetSer\n"
            + "tx seq -> : cctccagttcaatccccagcATGGCTGCGTCCACTATGTCT\n"
            + "tx pos    : |    .    |    .    |   .    |    .    | \n"
            + "          : -20       -10       1        10        20\n"
        ).split("\n")
        for r, e in zip(result, expected_str):
            self.assertEqual(e, r)

    def test_var_c1_reverse(self):
        """test c1 on <- strand"""

        hgvs_c = "NM_001177507.2:c.1="
        var_c = self.hp.parse(hgvs_c)
        result = self.pp.display(var_c)
        print(result)
        result = result.split("\n")
        expected_str = (
            "hgvs      : NC_000007.13:g.36763753=\n"
            + "hgvs      : NM_001177507.2:c.1=\n"
            + "          :        36,763,740          36,763,760\n"
            + "chrom pos :   .    |    .    |    .    |    .    |   \n"
            + "seq    -> : GATTTTCCAGGGGGACTGCATCTCCGAGCTATGCACCCCAA\n"
            + "region    :                     =                    \n"
            + "aa seq <- : IleLysTrpProSerGlnMet                    \n"
            + "tx seq <- : CTAAAAGGTCCCCCTGACGTAgaggctcgatacgtggggtt\n"
            + "tx pos    :  |    .    |    .   |    .    |    .    |\n"
            + "          :  20        10       1         -10       -20\n"
        ).split("\n")
        for r, e in zip(result, expected_str):
            self.assertEqual(e, r)

    def test_var_g_substitution(self):
        hgvs_g = "NC_000007.13:g.36561662C>T"
        var_g = self.hp.parse(hgvs_g)

        result = self.pp.display(var_g)
        print(result)
        result = result.split("\n")
        expected_str = (
            "hgvs      : NC_000007.13:g.36561662C>T\n"
            + "hgvs      : NM_001177507.2:c.1486G>A\n"
            + "          :         36,561,650          36,561,670\n"
            + "chrom pos :    .    |    .    |    .    |    .    |  \n"
            + "seq    -> : TACCTCGTTGGGGTGGAATCCATCCACGGGCTCGATGAGCT\n"
            + "region    :                     T                    \n"
            + "aa seq <- :    GluAsnProHisPheGlyAspValProGluIleLeuGl\n"
            + "tx seq <- :    GAGCAACCCCACCTTAGGTAGGTGCCCGAGCTACTCGA\n"
            + "tx pos    :       |    .    |    .    |    .    |    \n"
            + "          :       1500      1490      1480      1470\n"
        ).split("\n")
        for r, e in zip(result, expected_str):
            self.assertEqual(e, r)

    def test_var_g_ins(self):
        """[ATTA]x2 -> x3"""
        hgvs_g = "NC_000005.10:g.123346517_123346518insATTA"
        var_g = self.hp.parse(hgvs_g)

        result = self.pp.display(var_g)
        print(result)
        result = result.split("\n")
        expected_str = (
            "hgvs      : NC_000005.10:g.123346517_123346518insATTA\n"
            + "hgvs      : NM_001166226.1:c.*1_*2insTAAT\n"
            + self.atta_expected_results
        ).split("\n")
        for r, e in zip(result, expected_str):
            self.assertEqual(e, r)

    def test_var_g_dup(self):
        hgvs_g = "NC_000005.10:g.123346522_123346525dup"
        var_g = self.hp.parse(hgvs_g)

        result = self.pp.display(var_g)
        print(result)
        result = result.split("\n")
        expected_str = (
            "hgvs      : NC_000005.10:g.123346522_123346525dup\n"
            + "hgvs      : NM_001166226.1:c.2880_2883dup\n"
            + self.atta_expected_results
        ).split("\n")
        for r, e in zip(result, expected_str):
            self.assertEqual(e, r)

    def test_insertion(self):
        "A shuffleable insertion, shuffleable unit: TCGTCATC additional residues: G"
        hgvs_g = "NC_000004.11:g.1643284_1643285insTCGTCATCG"
        var_g = self.hp.parse(hgvs_g)

        result = self.pp.display(var_g)
        print(result)
        result = result.split("\n")
        expected_str = (
            "hgvs      : NC_000004.11:g.1643284_1643285insTCGTCATCG\n"
            + "hgvs      : NM_001174070.2:c.932_933insCGATGACGA\n"
            + "          :      1,643,270 1,643,280 1,643,290 1,643,300 1,643,310\n"
            + "chrom pos : .    |    .    |    .    |    .    |    .    |  \n"
            + "seq    -> : TCACTGGGGTGTCATCCTCATCGTCATCTTCGTAATTGAGGGAGCAAA\n"
            + "region    :                     |------|                    \n"            
            + "aa seq <- : sValProThrAspAspGluAspAspAspGluTyrAsnLeuSerCysLe\n"
            + "tx seq <- : AGTGACCCCACAGTAGGAGTAGCAGTAGAAGCATTAACTCCCTCGTTT\n"
            + "tx pos    :   |    .    |    .    |    .    |    .    |    .\n"
            + "          :   950       940       930       920       910\n"
            + "ref>alt   : TCGTCATC>TCGTCATCGTCGTCATC\n"
        ).split("\n")
        for r, e in zip(result, expected_str):
            self.assertEqual(e, r)

    def test_insertion_size_1(self):
        hgvs_g = "NC_000007.13:g.36561662_36561663insT"
        var_g = self.hp.parse(hgvs_g)

        result = self.pp.display(var_g)
        print(result)

        result = result.split("\n")
        expected_str = (
            "hgvs      : NC_000007.13:g.36561662_36561663insT\n"
            + "hgvs      : NM_001177507.2:c.1485_1486insA\n"
            + "          :        36,561,650          36,561,670\n"
            + "chrom pos :   .    |    .    |    .    |    .    |  \n"
            + "seq    -> : ACCTCGTTGGGGTGGAATCCATCCACGGGCTCGATGAGCT\n"
            + "region    :                    ^^                   \n"
            + "aa seq <- :   GluAsnProHisPheGlyAspValProGluIleLeuGl\n"
            + "tx seq <- :   GAGCAACCCCACCTTAGGTAGGTGCCCGAGCTACTCGA\n"  
            + "tx pos    :      |    .    |    .    |    .    |    \n"
            + "          :      1500      1490      1480      1470\n"
        ).split("\n")
        for r, e in zip(result, expected_str):
            self.assertEqual(e, r)

    def test_del_2bp(self):
        hgvs_g = "NC_000007.13:g.36561662_36561663del"
        var_g = self.hp.parse(hgvs_g)

        result = self.pp.display(var_g)
        print(result)

        result = result.split("\n")
        expected_str = (
            "hgvs      : NC_000007.13:g.36561662_36561663del\n"
            + "hgvs      : NM_001177507.2:c.1485_1486del\n"
            + "          :         36,561,650          36,561,670\n"
            + "chrom pos :    .    |    .    |    .    |    .    |   \n"
            + "seq    -> : TACCTCGTTGGGGTGGAATCCATCCACGGGCTCGATGAGCTG\n"
            + "region    :                     xx                    \n"
            + "aa seq <- :    GluAsnProHisPheGlyAspValProGluIleLeuGln\n"
            + "tx seq <- :    GAGCAACCCCACCTTAGGTAGGTGCCCGAGCTACTCGAC\n"
            + "tx pos    :       |    .    |    .    |    .    |    .\n"
            + "          :       1500      1490      1480      1470\n"
            
        ).split("\n")
        for r, e in zip(result, expected_str):
            self.assertEqual(e, r)

    def test_del_1bp_shuffleable(self):
        hgvs_g = "NC_000007.13:g.36561662del"
        var_g = self.hp.parse(hgvs_g)

        result = self.pp.display(var_g)
        print(result)

        result = result.split("\n")
        expected_str = (
            "hgvs      : NC_000007.13:g.36561662del\n"
            + "hgvs      : NM_001177507.2:c.1487del\n"
            + "          :          36,561,650          36,561,670\n"
            + "chrom pos :     .    |    .    |    .    |    .    |  \n"
            + "seq    -> : TTACCTCGTTGGGGTGGAATCCATCCACGGGCTCGATGAGCT\n"
            + "region    :                     xx                    \n"
            + "aa seq <- :     GluAsnProHisPheGlyAspValProGluIleLeuGl\n"
            + "tx seq <- :     GAGCAACCCCACCTTAGGTAGGTGCCCGAGCTACTCGA\n"
            + "tx pos    :        |    .    |    .    |    .    |    \n"
            + "          :        1500      1490      1480      1470\n"
            
            + "ref>alt   : CC>C\n"
        ).split("\n")
        for r, e in zip(result, expected_str):
            self.assertEqual(e, r)

    def test_del_1bp(self):
        hgvs_g = "NC_000007.13:g.36561663del"
        var_g = self.hp.parse(hgvs_g)

        result = self.pp.display(var_g)
        print(result)

        result = result.split("\n")
        expected_str = (
            "hgvs      : NC_000007.13:g.36561663del\n"
            + "hgvs      : NM_001177507.2:c.1485del\n"
            + "          :        36,561,650          36,561,670\n"
            + "chrom pos :   .    |    .    |    .    |    .    |   \n"
            + "seq    -> : ACCTCGTTGGGGTGGAATCCATCCACGGGCTCGATGAGCTG\n"
            + "region    :                     x                    \n"
            + "aa seq <- :   GluAsnProHisPheGlyAspValProGluIleLeuGln\n"
            + "tx seq <- :   GAGCAACCCCACCTTAGGTAGGTGCCCGAGCTACTCGAC\n"
            + "tx pos    :      |    .    |    .    |    .    |    .\n"
            + "          :      1500      1490      1480      1470\n"
                        
        ).split("\n")
        for r, e in zip(result, expected_str):
            self.assertEqual(e, r)

    def test_dup_1bp_shuffleable(self):
        hgvs_g = "NC_000007.13:g.36561662dup"
        var_g = self.hp.parse(hgvs_g)

        result = self.pp.display(var_g)
        print(result)

        result = result.split("\n")
        expected_str = (
            "hgvs      : NC_000007.13:g.36561662dup\n"
            + "hgvs      : NM_001177507.2:c.1487dup\n"
            + "          :          36,561,650          36,561,670\n"
            + "chrom pos :     .    |    .    |    .    |    .    |  \n"
            + "seq    -> : TTACCTCGTTGGGGTGGAATCCATCCACGGGCTCGATGAGCT\n"
            + "region    :                     ||                    \n"
            + "aa seq <- :     GluAsnProHisPheGlyAspValProGluIleLeuGl\n"
            + "tx seq <- :     GAGCAACCCCACCTTAGGTAGGTGCCCGAGCTACTCGA\n"
            + "tx pos    :        |    .    |    .    |    .    |    \n"
            + "          :        1500      1490      1480      1470\n"
            + "ref>alt   : CC>CCC\n"
        ).split("\n")
        for r, e in zip(result, expected_str):
            self.assertEqual(e, r)

    def test_dup_1bp(self):
        hgvs_g = "NC_000007.13:g.36561663dup"
        var_g = self.hp.parse(hgvs_g)

        result = self.pp.display(var_g)
        print(result)

        result = result.split("\n")
        expected_str = (
            "hgvs      : NC_000007.13:g.36561663dup\n"
            + "hgvs      : NM_001177507.2:c.1485dup\n"
            + "          :        36,561,650          36,561,670\n"
            + "chrom pos :   .    |    .    |    .    |    .    |   \n"
            + "seq    -> : ACCTCGTTGGGGTGGAATCCATCCACGGGCTCGATGAGCTG\n"
            + "region    :                     |                    \n"
            + "aa seq <- :   GluAsnProHisPheGlyAspValProGluIleLeuGln\n"
            + "tx seq <- :   GAGCAACCCCACCTTAGGTAGGTGCCCGAGCTACTCGAC\n"
            + "tx pos    :      |    .    |    .    |    .    |    .\n"
            + "          :      1500      1490      1480      1470\n"
            + "ref>alt   : A>AA\n"
        ).split("\n")
        for r, e in zip(result, expected_str):
            self.assertEqual(e, r)

    def test_identity(self):

        hgvs_g = "NC_000007.13:g.36561663="
        var_g = self.hp.parse(hgvs_g)

        result = self.pp.display(var_g)
        print(result)

        result = result.split("\n")
        expected_str = (
            "hgvs      : NC_000007.13:g.36561663=\n"
            + "hgvs      : NM_001177507.2:c.1485=\n"
            + "          :        36,561,650          36,561,670\n"
            + "chrom pos :   .    |    .    |    .    |    .    |   \n"
            + "seq    -> : ACCTCGTTGGGGTGGAATCCATCCACGGGCTCGATGAGCTG\n"
            + "region    :                     =                    \n"
            + "aa seq <- :   GluAsnProHisPheGlyAspValProGluIleLeuGln\n"
            + "tx seq <- :   GAGCAACCCCACCTTAGGTAGGTGCCCGAGCTACTCGAC\n"
            + "tx pos    :      |    .    |    .    |    .    |    .\n"
            + "          :      1500      1490      1480      1470\n"
            
        ).split("\n")
        for r, e in zip(result, expected_str):
            self.assertEqual(e, r)

    def test_tiny(self):
        """Test a variant with bad input."""
        hgvs_g = "NC_000005.10:g.123346517_123346518insATTA"
        var_g = self.hp.parse(hgvs_g)

        tiny_pp = PrettyPrint(self.hdp, padding_left=0, padding_right=0)

        result = tiny_pp.display(var_g)
        print(result)

        result = result.split("\n")

        expected_str = (
            "hgvs      : NC_000005.10:g.123346517_123346518insATTA\n"
            + "hgvs      : NM_001166226.1:c.*1_*2insTAAT\n"
            + "          :   123,346,520\n"
            + "chrom pos :   |    .\n"
            + "seq    -> : ATTAATTA\n"
            + "region    : |------|\n"
            + "aa seq <- :  TerAsnS\n"
            + "tx seq <- : TAATTAAT\n"
            + "tx pos    : |   |   \n"
            + "          : *1  2880\n"
            + "ref>alt   : ATTAATTA>ATTAATTAATTA"
        ).split("\n")

        for r, e in zip(result, expected_str):
            self.assertEqual(e, r)

    @pytest.mark.skip(reason="CNVs not implemented yet")
    def test_cnv(self):
        """Test a CNV variant. TODO: make display compact"""
        hgvs_g = "NC_000005.10:g.123345517_123346518del"
        var_g = self.hp.parse(hgvs_g)

        result = self.pp.display(var_g)

        print(result)

    def test_hgvs_c(self):
        """Test a hgvs_c variant overlapping start codon on reverse strand."""
        hgvs_c = "NM_004572.3:c.-9_12dup"
        var_c = self.hp.parse(hgvs_c)
        pp = PrettyPrint(
            self.hdp, padding_left=10, padding_right=110, useColor=False, showLegend=False
        )
        result = pp.display(var_c)

        print(result)

        result = result.split("\n")
        expected_str = (
            "NC_000012.11:g.33049660_33049680dup\n"
            + "NM_004572.3:c.-9_12dup\n"
            + "      33,049,650          33,049,670          33,049,690          33,049,710          33,049,730          33,049,750          33,049,770          33,049,790\n"
            + " .    |    .    |    .    |    .    |    .    |    .    |    .    |    .    |    .    |    .    |    .    |    .    |    .    |    .    |    .    |\n"
            + "CTGGGGCGCCGGGGGCTGCCATGGGGCCGGTGGGGGCGACCGAGCTGCTCGCCTGCCTCTGGACTCGCGGGCGAAGCCGCCACGGAGCTGGGGGCGCTGGCGCGAGCCCCGCCCCGCTCGAGTCCGGCCCCGCCCCTGGCCCGCCCC\n"
            + "          |-------------------------|                                                                                                              \n"
            + "aProAlaGlyProAlaAlaMet                                                                                                                             \n"
            + "GACCCCGCGGCCCCCGACGGTAccccggccacccccgctggctcgacgagcggacggagacctgagcgcccgcttcggcggtgcctcgacccccgcgaccgcgctcggggcggggcgagctcaggccggggcgggga          \n"
            + "  |    .    |    .   |    .    |    .    |    .    |    .    |    .    |    .    |    .    |    .    |    .    |    .    |    .    |    .          \n"
            + "  20        10       1         -10       -20       -30       -40       -50       -60       -70       -80       -90       -100      -110\n"
            + "GGGGGCTGCCATGGGGCCGGTGGGGGC>GGGGGCTGCCATGGGGCCGGTGGGGGCTGCCATGGGGCCGGTGGGGGC"
        ).split("\n")
        for r, e in zip(result, expected_str):
            self.assertEqual(e, r)

    def test_ref_disagree(self):
        """Test a tx ref disagree variant."""
        hgvs_g = "NM_001111.4:c.298G>A"
        var_g = self.hp.parse(hgvs_g)
        result = self.pp.display(var_g)

        print(result)

        result = result.split("\n")

        # note the X in the transscript sequence
        expected_str = (
            "hgvs      : NC_000001.10:g.154574820=\n"
            + "hgvs      : NM_001111.4:c.298G>A\n"
            + "          : 154,574,800         154,574,820         154,574,840\n"
            + "chrom pos : |    .    |    .    |    .    |    .    |\n"
            + "seq    -> : TCTCTGGAGCCCCTGACTTCTGAGATGCACGCCCCTGGGGA\n"
            + "tx ref dif:                     X                    \n"
            + "region    :                     T                    \n"
            + "aa seq <- : ArgGlnLeuGlyGlnSerGlyLeuHisValGlyArgProVa\n"
            + "tx seq <- : AGAGACCTCGGGGACTGAAGGCTCTACGTGCGGGGACCCCT\n"
            + "tx pos    :    .    |    .    |    .    |    .    |  \n"
            + "          :         310       300       290       280\n"
                        
        ).split("\n")
        for r, e in zip(result, expected_str):
            self.assertEqual(e, r)

    def test_ref_disagree_ref_ins(self):
        """A ref disagree with a region inserted inref, that is missing in transcript"""
        # hgvs_g = "NC_000001.10:g.154574820_154574821delinsCA"
        # var_g = self.hp.parse(hgvs_g)
        # hgvs_c = "NM_020469.2:c.188_189="
        # hgvs_c = "NM_003777.3:c.5475dup" # an I variant
        # hgvs_c =

        hgvs_c = "NM_198689.2:c.124_135="
        # this would match chromosome: "NM_198689.2:c.124_135insCTGCTGCGCCCCCAG"
        var_c = self.hp.parse(hgvs_c)
        pp = PrettyPrint(self.hdp, infer_hgvs_c=True, padding_left=30, padding_right=40)
        result = pp.display(var_c)
        print(result)
        result = result.split("\n")
        expected_str = (
            "hgvs      : NC_000021.8:g.46020668_46020682del\n"
            "hgvs      : NM_198689.2:c.124_135=\n"
            "          :     46,020,630          46,020,650          46,020,670          46,020,690          46,020,710\n"
            "chrom pos :     |    .    |    .    |    .    |    .    |    .    |    .    |    .    |    .    |    .    |  \n"
            "seq    -> : CGACTGCCCAGAGAGCTGCTGCGAGCCCCCCTGCTGCGCCCCCAGCTGCTGCGCCCCGGCCCCCTGCCTGAGCCTGGTCTGCACCCCAGTGAGCCGT\n"
            "tx ref dif:                               IIIIIIIIIIIIIII                                                 XX \n"
            "region    :                               =-------------------------=                                        \n"
            "aa seq -> : pAspCysProGluSerCysCysGluProPr---------------oCysCysAlaProAlaProCysLeuSerLeuValCysThrProValSerTyr\n"
            "tx seq -> : CGACTGCCCAGAGAGCTGCTGCGAGCCCCC---------------CTGCTGCGCCCCGGCCCCCTGCCTGAGCCTGGTCTGCACCCCAGTGAGCTAT\n"
            "tx pos    : .    |    .    |    .    |                   .    |    .    |    .    |    .    |    .    |    . \n"
            "          :      110       120       130                      140       150       160       170       180\n"
            "ref>alt   : CTGCTGCGCCCCCAGCTGCTGCGCCCC>CTGCTGCGCCCC\n"
        ).split("\n")

        for r, e in zip(result, expected_str):
            self.assertEqual(e, r)

    def test_ref_disagree_del(self):
        # hgvs_g = "NC_000001.10:g.154574820_154574821delinsCA" - one base is a svn relative to the tx and part of the variant -> NM_001025107.2:c.-589C>T 
        # var_g = self.hp.parse(hgvs_g)
        # hgvs_c = "NM_020469.2:c.188_189=" is  NC_000009.11:g.136135237_136135238delinsGC in ref


        hgvs_c = "NM_000682.6:c.901_911del"  # a del variant

        var_c = self.hp.parse(hgvs_c)
        pp = PrettyPrint(self.hdp, infer_hgvs_c=True, padding_left=30, padding_right=40)
        result = pp.display(var_c)
        print(result)
        result = result.split("\n")

        expected_str = (
            "hgvs      : NC_000002.11:g.96780987_96780997del\n"
            + "hgvs      : NM_000682.6:c.901_911del\n"
            + "          :     96,780,960          96,780,980                   96,781,000          96,781,020\n"
            + "chrom pos :     |    .    |    .    |    .    |    .  _________  |    .    |    .    |    .    |    .  \n"
            + "seq    -> : TGCCTGGGGTTCACACTCTTCCTCCTCCTCCTCCTCCTCTTC.........GGCTTCATCCTCTGGAGATGCCCCACAAACACCCTCCTTC\n"
            + "tx ref dif:                                           DDDDDDDDD                               \n"
            + "region    :                               x----------x                                                 \n"
            + "aa seq <- : AlaGlnProGluCysGluGluGluGluGluGluGluGluGluGluGluGluAlaGluAspGluProSerAlaGlyCysValGlyGluLysG\n"
            + "tx seq <- : ACGGACCCCAAGTGTGAGAAGGAGGAGGAGGAGGAGGAGAAGGAGGAGAAGTCGAAGTAGGAGACCTCTACGGGGTGTTTGTGGGAGGAAG\n"
            + "tx pos    :   |    .    |    .    |    .    |    .    |    .    |    .    |    .    |    .    |    .   \n"
            + "          :   940       930       920       910       900       890       880       870       860\n"
            + "ref>alt   : CTCCTCCTCTTC>C\n"
        ).split("\n")
        for r, e in zip(result, expected_str):
            self.assertEqual(e, r)

    def test_ref_disagree_tx_ins(self):
        """ A variant that is a tx ref disagree with an insertion in the ref genome"""
        aa!


    @pytest.mark.skip(reason="actually not that special, but still a nice variant.")
    def test_exon_boundary_overlap_forward_strand(self):
        hgvs_c = "NM_001283009.2:c.1228_1266+39del"
        var_c = self.hp.parse(hgvs_c)
        pp = PrettyPrint(self.hdp, showLegend=True, useColor=False)

        result = pp.display(var_c)

        print(result)

        assert False

    def test_ruler(self):
        """Test the ruler display option turned on."""
        hgvs_c = "NM_001111.4:c.298G>A"
        var_c = self.hp.parse(hgvs_c)
        pp = PrettyPrint(self.hdp, showLegend=False)

        result = pp.display(var_c)

        print(result)

        result = result.split("\n")

        # note the X in the transscript ref-disagree row
        expected_str = (
            "NC_000001.10:g.154574820=\n"
            + "NM_001111.4:c.298G>A\n"
            + "154,574,800         154,574,820         154,574,840\n"
            + "|    .    |    .    |    .    |    .    |\n"
            + "TCTCTGGAGCCCCTGACTTCTGAGATGCACGCCCCTGGGGA\n"
            + "                    X                    \n"
            + "                    T                    \n"
            + "ArgGlnLeuGlyGlnSerGlyLeuHisValGlyArgProVa\n"
            + "AGAGACCTCGGGGACTGAAGGCTCTACGTGCGGGGACCCCT\n"
            + "   .    |    .    |    .    |    .    |  \n"
            + "        310       300       290       280\n"
        ).split("\n")
        for r, e in zip(result, expected_str):
            self.assertEqual(e, r)
