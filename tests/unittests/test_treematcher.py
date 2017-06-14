from __future__ import absolute_import
from ete3 import Tree, PhyloTree, NCBITaxa
from treematcher import TreePattern
import unittest


class Test_TreePattern(unittest.TestCase):
    def test_ete_params(self):

        """
        tests exact match is working
        tests basic ete parameters like node.name, node.dist, node.support are working

        """
        pattern0 = """
             ( bye , kk );
            """
        pattern1 = """
            ( hello , kk );
            """
        pattern2 = """
            (
            'len(@.children) > 2 and @.name in ("bye","kk") '
            )
            '(len(@.name) > 3) and @.dist >= 0.5';
            """

        pattern0 = TreePattern(pattern0)
        pattern1 = TreePattern(pattern1)
        pattern2 = TreePattern(pattern2)

        tree = Tree("(hello,(1,2,3)kk)pasa:1;", format=1)
        match0 = pattern0.find_match(tree)
        match1 = pattern1.find_match(tree)
        match2 = pattern2.find_match(tree)


        self.assertEqual(list(match0), [])
        self.assertEqual(list(match1)[0].name, 'pasa')
        self.assertEqual(list(match2)[0].name, 'pasa')

    def test_maxhits(self):

        """
        tests syntax we've created like node.leaves, node.children are working
        """
        pattern1 = """' "Pan_troglodytes_1" in leaves(@)';"""
        pattern1 = TreePattern(pattern1)

        tree = PhyloTree(
            "((((Anolis_carolinensis_1:1, Gallus_gallus_1:1), (Felis_catus_1:1, (Homo_sapiens_1:1, Pan_troglodytes_1:1)p4)p3)p2, ((Danio_rerio_1:1, (Xenopus_laevis_1:1, Anolis_carolinensis_1:1)), Saccharomyces_cerevisiae_2:1))p1, Saccharomyces_cerevisiae_1:1)root;",
            format=1)

        # note that maxhits<1 not allowed, uses default which is one match
        one_match = list(pattern1.find_match(tree))
        self.assertEqual(len(one_match), 1)
        self.assertEqual(one_match[0].name, 'root')

        three_matches= list(pattern1.find_match(tree, maxhits=3))
        self.assertEqual(len(three_matches), 3)

        pan_ancestor = list(pattern1.find_match(tree, maxhits=None))  #note that leaf constraint also returns leaf


        self.assertEqual(len(pan_ancestor), 6)
        self.assertEqual(pan_ancestor[0].name, 'root')
        self.assertEqual(pan_ancestor[1].name, 'p1')
        self.assertEqual(pan_ancestor[2].name, 'p2')
        self.assertEqual(pan_ancestor[3].name, 'p3')
        self.assertEqual(pan_ancestor[4].name, 'p4')
        self.assertEqual(pan_ancestor[5].name, 'Pan_troglodytes_1')

        tree2 = PhyloTree(
            "((((Anolis_carolinensis_1:1, Gallus_gallus_1:1), (Felis_catus_1:1, (Homo_sapiens_1:1, Pan_troglodytes_2:1))), ((Danio_rerio_1:1, (Xenopus_laevis_1:1, Anolis_carolinensis_1:1)), Saccharomyces_cerevisiae_2:1)), Saccharomyces_cerevisiae_1:1);",
            format=1)

        self.assertEqual(len(list(pattern1.find_match(tree2))), 0)


    def test_species(self):
        """
        tests if node.species and ncbi_query are working
        """

        # test node.species

        species_tree = PhyloTree(
            """(Felis_catus_1:1,
                (Homo_sapiens_1:1, Pan_troglodytes_1:1),
                Saccharomyces_cerevisiae_1:1);""",
            format=1)
        species_tree.set_species_naming_function(lambda n: n.name.split("_")[1] if "_" in n.name else '')

        pattern0 = """('',
                       (' len(set(["sapiens","pygmaeus"]) & species(@))>0',
                       Pan_troglodytes_1)
                       );"""

        pattern0 = TreePattern(pattern0)


        root = species_tree.get_tree_root()
        self.assertEqual(list(pattern0.find_match(species_tree)), [root])

        # test ncbi taxonomy

        ncbi = NCBITaxa()
        taxonomy_tree = PhyloTree("((9598, 9606), 10090);", sp_naming_function=lambda name: name)
        taxonomy_tree.annotate_ncbi_taxa()
        root = taxonomy_tree.get_tree_root()

        pattern1 = """ '  @.sci_name == "Euarchontoglires" ';"""
        pattern2 = """
          (( '@.sci_name=="Homo sapiens"' , '9526 in @.lineage ' )' @.rank=="subfamily" and @.taxid == 207598 ')
          '  @.sci_name == "Euarchontoglires" and "cellular organisms" in @.named_lineage';
          """

        pattern1 = TreePattern(pattern1)
        pattern2 = TreePattern(pattern2)

        match1 = pattern1.find_match(taxonomy_tree)
        match2 = pattern2.find_match(taxonomy_tree)

        self.assertEqual(list(match1), [root])
        self.assertEqual(list(match2), [root])


    def test_lineages(self):
        """
        Search trees (naming format: NumericTaxid.SequenceName)
        for nodes containing branches that separate two groups of primate genes where,
        in one side, the human gene has been lost,
        and the branch support value of the matching node is higher than 0.9.

                                  /-Any primate taxid (9443 in lineage)
        support >= 0.9--|
                                  \-Any primate taxid except human

        """
        t1 = PhyloTree("(9601.ENSPPYP00000022176:1,9593.ENSGGOP00000009720:1);")
        t2 = PhyloTree("(9361.ENSDNOP00000016844:1,9258.ENSOANP00000032529:1);")
        t3 = PhyloTree(
            "(((((37347.ENSTBEP00000010698:0.120098,(9361.ENSDNOP00000000113:0.0697238,(9785.ENSLAFP00000009564:0.0297499,(9371.ENSETEP00000002412:0.0588324,9813.ENSPCAP00000006440:0.026638)0.985184:0.0242194)0.99985:0.0211882)0.99706:0.0161759)0.756:0.00666819,((132908.ENSPVAP00000002358:0.0439546,59463.ENSMLUP00000004598:0.0635161)0.994843:0.00885432,(9796.ENSECAP00000009809:0.0292517,((9685.ENSFCAP00000004938:0.056779,(9615.ENSCAFP00000008559:0.039179,(9823.ENSSSCP00000024070:0.126803,(9669.ENSMPUP00000010096:0.0341928,9646.ENSAMEP00000005906:0.0189746)0.995231:0.00951966)0.915476:0.0046099)0.949664:0.00417374)0.99985:0.0133593,(9739.ENSTTRP00000009464:0.0664336,9913.ENSBTAP00000001687:0.036632)0.99985:0.0236174)0.939309:0.00508062)0.991475:0.00823937)0.99985:0.0107263)0.99985:0.0100107,((9986.ENSOCUP00000014919:0.0830612,10141.ENSCPOP00000005291:0.12195)0.99985:0.0202639,((9483.ENSCJAP00000047968:0.0446865,(9544.ENSMMUP00000007168:0.0201746,((9593.ENSGGOP00000005929:0.00916494,(9606.ENSP00000294053:1.3e-07,9598.ENSPTRP00000006940:0.0068176)0.955193:0.00220905)0.99985:0.00778854,(9601.ENSPPYP00000004174:0.00495163,61853.ENSNLEP00000020892:0.179569)0.290072:0.00153447)0.998732:0.00889714)0.99985:0.0144864)0.99985:0.0344562,(9478.ENSTSYP00000006073:0.129349,(30608.ENSMICP00000010690:0.0852248,30611.ENSOGAP00000013738:0.0467206)0.99985:0.0188861)0.232709:0.00179852)0.99985:0.00929928)0.51042:0.00516905)0.367617:0.00813494,(43179.ENSSTOP00000004287:0.0599707,(10020.ENSDORP00000000618:0.138502,(10116.ENSRNOP00000026665:0.0528487,10090.ENSMUSP00000001884:0.0307781)0.99985:0.089983)0.99985:0.018366)0.698647:0.00414256)0.995833:0.06629,(9258.ENSOANP00000012946:0.33344,(13616.ENSMODP00000032549:0.0348012,(9315.ENSMEUP00000011030:0.0138664,9305.ENSSHAP00000003293:0.0185119)0.570293:0.0137766)0.99985:0.143897)0.995833:0.06629);")
        t4 = PhyloTree("(9593.ENSGGOP00000025542:1,9601.ENSPPYP00000004907:1);")
        t5 = PhyloTree(
            "(9371.ENSETEP00000005103:0.0955875,(9785.ENSLAFP00000014743:0.0214619,(9813.ENSPCAP00000005573:0.0376639,(9796.ENSECAP00000019319:0.0196571,(37347.ENSTBEP00000012329:0.0242927,((9361.ENSDNOP00000011716:0.0676669,(9606.ENSP00000374323:9e-07,(9593.ENSGGOP00000028731:0.00246332,(61853.ENSNLEP00000002377:0.0030064,(9601.ENSPPYP00000015233:0.0112606,(9598.ENSPTRP00000026129:0.00246268,9483.ENSCJAP00000015834:0.0290829)0:1.2e-07)0:6.5e-07)0.146278:0.00614181)0.146329:0.00485474)0.991187:0.014264)0.763764:0.00352544,((10020.ENSDORP00000008692:0.0259566,(30608.ENSMICP00000002718:0.0380742,9478.ENSTSYP00000009200:0.0174548)0.197348:0.00155005)0.99985:0.0110622,((((132908.ENSPVAP00000013183:0.0099908,59463.ENSMLUP00000014424:0.0115111)0.99985:0.00655941,(10141.ENSCPOP00000003417:0.0535498,((9669.ENSMPUP00000002651:0.0156675,(9646.ENSAMEP00000014393:0.0142536,9615.ENSCAFP00000013394:0.00243184)0.930921:0.00345947)0.99985:0.015828,(9913.ENSBTAP00000053531:0.0545233,9739.ENSTTRP00000001508:0.0344514)0.985783:0.00536759)0:1.1e-07)0:1.1e-07)0.99985:0.00795592,(10090.ENSMUSP00000066734:0.0572278,(43179.ENSSTOP00000020881:0.021661,30611.ENSOGAP00000000479:0.00876016)0.955042:0.00724791)0.992776:0.0044053)0:3.4e-07,(9258.ENSOANP00000012014:0.10692,(9315.ENSMEUP00000001901:0.0451997,13616.ENSMODP00000021214:0.00830289)0.994926:0.0229072)0.99985:0.0500253)0.981032:0.00621499)0:9e-08)0.723103:0.00185076)0.580248:0.00162611)0.99985:0.0167207)0.863552:0.00574499)1:0.0955875);")
        t6 = PhyloTree(
            "((9305.ENSSHAP00000010229:0.0607855,13616.ENSMODP00000009656:0.0615237)0.99985:0.0877765,(9785.ENSLAFP00000028174:0.0885004,(((9823.ENSSSCP00000002806:0.0860827,9823.ENSSSCP00000002780:0.0111508)0.99985:0.122086,((9913.ENSBTAP00000038896:0.050358,(9685.ENSFCAP00000017257:0.0778567,(9986.ENSOCUP00000017975:0.161424,(9615.ENSCAFP00000020783:0.056902,(9646.ENSAMEP00000019763:0.0857189,9669.ENSMPUP00000019474:0.0325693)0.99985:0.0314116)0.875671:0.00690881)0.942895:0.0136375)0.798192:0.00741364)0.967573:0.0100004,(59463.ENSMLUP00000020576:0.0755216,9796.ENSECAP00000004613:0.0777605)0.799782:0.00471384)0.911021:0.00832673)0.659845:0.00664335,((43179.ENSSTOP00000021465:0.123042,9593.ENSGGOP00000020601:0.0781752)0.987812:0.0311266,(30611.ENSOGAP00000021055:0.090792,(10116.ENSRNOP00000016702:0.0112116,10090.ENSMUSP00000050705:0.0330259)0.99985:0.134681)0.972881:0.0174783)0.998643:0.0179346)0.901179:0.017737)0.99985:0.0877765);")
        t7 = PhyloTree(
            "(9258.ENSOANP00000017269:0.144169,(((10090.ENSMUSP00000089169:0.0424834,10116.ENSRNOP00000026070:0.0151696)0.99985:0.0742333,(((((132908.ENSPVAP00000008558:0.0138473,(30608.ENSMICP00000004293:1.5e-07,((9986.ENSOCUP00000020707:0.0691049,37347.ENSTBEP00000002617:0.0138881)0:1.2e-07,(9371.ENSETEP00000012957:0.0515389,(9785.ENSLAFP00000009919:0.0260641,9813.ENSPCAP00000013834:0.0329521)0.741149:0.0041225)0.998768:0.00855745)0.99985:0.0111961)0.867255:0.00524663)0:4.3e-07,(9361.ENSDNOP00000010929:0.0359312,(9739.ENSTTRP00000015818:0.0267351,9796.ENSECAP00000009501:0.0168218)0.868862:0.00355516)0:8e-08)0.99985:0.0056594,(9913.ENSBTAP00000012912:0.0231165,(9669.ENSMPUP00000002012:0.00320767,9823.ENSSSCP00000023102:0.0629927)0.99134:0.00309237)0.988361:0.00284581)0:1.5e-07,((59463.ENSMLUP00000015155:0.0360776,9615.ENSCAFP00000002053:0.00579656)0.961397:0.00553059,(9685.ENSFCAP00000023114:0.0115974,9646.ENSAMEP00000004090:0.00575272)0.959045:0.00279601)0.988458:0.00279093)0.998008:0.00284847,(30611.ENSOGAP00000001383:0.00849776,((9483.ENSCJAP00000006698:0.0114709,(9544.ENSMMUP00000006654:0.00568623,(61853.ENSNLEP00000004122:0.00566385,(9601.ENSPPYP00000021653:0.00853215,(9593.ENSGGOP00000020462:1.8e-07,(9598.ENSPTRP00000035990:1e-08,9606.ENSP00000365550:1e-08)0.99985:0.00282071)0.996162:0.00281965)0:1.7e-07)0:8e-08)0.954037:0.0027827)0.99985:0.00818313,(43179.ENSSTOP00000012068:0.0109022,(9478.ENSTSYP00000008441:0.0132658,10141.ENSCPOP00000000986:0.0564111)0.314526:0.00294575)0:7e-08)0.980721:0.00309462)0.991529:0.00280168)0:1.6e-07)0.99985:0.0483405,(9315.ENSMEUP00000015273:0.00839008,(9305.ENSSHAP00000020642:0.00542335,13616.ENSMODP00000010568:0.101485)0:2.1e-07)0.99985:0.0336521)1:0.144169);")
        t8 = PhyloTree(
            "(((9371.ENSETEP00000003671:0.0131637,(9258.ENSOANP00000006745:0.117598,(132908.ENSPVAP00000001122:0.0159907,(30611.ENSOGAP00000013217:0.0071702,(((9823.ENSSSCP00000000042:0.0144457,(9646.ENSAMEP00000009872:0.0154876,9361.ENSDNOP00000012437:0.0817179)0:1e-06)0.998538:0.00765581,(9544.ENSMMUP00000001765:1e-08,(10116.ENSRNOP00000010491:0.0292686,(9669.ENSMPUP00000016236:0.340739,9615.ENSCAFP00000001415:4e-07)0.989009:0.00985882)0:8.7e-07)0:8.7e-07)0.99736:0.00973955,(((9606.ENSP00000379704:1e-08,(9601.ENSPPYP00000013264:0.00772278,9598.ENSPTRP00000024873:1e-08)0:2.3e-07)0.996569:0.00720502,(9913.ENSBTAP00000017531:0.0145949,9739.ENSTTRP00000016448:0.00723237)0.996503:0.00710774)0:4.2e-07,((9593.ENSGGOP00000008768:0.270021,(9785.ENSLAFP00000013194:0.00881524,9478.ENSTSYP00000011482:6.1e-07)0.482225:0.00675219)0.500314:0.00675139,(((59463.ENSMLUP00000002337:0.0319341,30608.ENSMICP00000003266:6.2e-07)0.987498:0.010619,(9796.ENSECAP00000021110:0.0073991,(9986.ENSOCUP00000007142:0.0196352,37347.ENSTBEP00000000333:0.0989537)0:9.5e-07)0:1.09e-06)0.873107:0.00951386,((9685.ENSFCAP00000000826:3e-07,(43179.ENSSTOP00000011619:0.00863897,10090.ENSMUSP00000023095:1e-08)0:1e-08)0.99985:0.132958,(10020.ENSDORP00000013215:0.0339132,10141.ENSCPOP00000011894:4.1e-07)0:4.1e-07)0.524756:0.00714334)0:8.1e-07)0.99985:0.00971634)0:7e-08)0:7e-08)0.772739:0.0177399)0.992096:0.0404786)0.817723:0.0310407)0.522416:0.072068,(9305.ENSSHAP00000014579:0.246289,9315.ENSMEUP00000008760:0.0666798)0.977479:0.195421)0.99985:1.2587,((((37347.ENSTBEP00000000946:0.0956163,(9483.ENSCJAP00000024301:0.0743892,(9593.ENSGGOP00000012469:0.00721405,(9606.ENSP00000391249:1e-08,9606.ENSP00000461549:1e-08)0:1.3e-07)0.993649:0.00856538)0.99985:0.0230549)0.975176:0.0143781,(30611.ENSOGAP00000003324:0.104251,30608.ENSMICP00000007369:0.0381575)0.990656:0.0183563)0.916137:0.00581305,(9823.ENSSSCP00000018191:0.0558998,((10020.ENSDORP00000010153:0.197695,((9796.ENSECAP00000018039:0.0363101,132908.ENSPVAP00000013461:0.0941126)0.892367:0.013635,((9739.ENSTTRP00000004783:0.0138565,9913.ENSBTAP00000003415:0.0166473)0.99985:0.0326524,((9371.ENSETEP00000006140:0.107709,(9785.ENSLAFP00000006435:0.170692,9813.ENSPCAP00000005503:0.0655274)0:2.68e-06)0.99985:0.0526328,(9258.ENSOANP00000002804:0.150016,(9315.ENSMEUP00000001056:0.0197146,(13616.ENSMODP00000002021:0.0382813,9305.ENSSHAP00000007534:0.0357616)0.99985:0.0843541)0.99985:0.115238)0.99985:0.133971)0.964252:0.0135998)0.99559:0.0163904)0.732303:0.00993157)0.99985:0.0470037,(9685.ENSFCAP00000008713:0.124988,(9615.ENSCAFP00000007771:0.0225216,(9646.ENSAMEP00000014479:0.0718956,9669.ENSMPUP00000013273:0.0487162)0.99985:0.0148769)0:9.2e-07)0.99985:0.0433867)0.99277:0.027679)0.99985:0.0134312)0:4.7e-07,(43179.ENSSTOP00000019919:0.152642,((10116.ENSRNOP00000003891:0.158016,10090.ENSMUSP00000091435:0.0102936)0.99985:0.0704992,(10141.ENSCPOP00000011436:0.130601,9986.ENSOCUP00000015843:0.529405)0:5.42e-06)0.909203:0.011833)0.428577:0.0186403)0.99985:1.2587);")
        t9 = PhyloTree("(9305.ENSSHAP00000009662:1,9305.ENSSHAP00000009620:1);")
        t10 = PhyloTree("((9315.ENSMEUP00000008285:0.899711,9258.ENSOANP00000027752:0.559777)0.99985:0.11989,((9739.ENSTTRP00000010720:0.164873,9913.ENSBTAP00000003500:0.298158)0.99985:0.109903,((9685.ENSFCAP00000006440:0.239731,(9615.ENSCAFP00000042310:0.122399,(9646.ENSAMEP00000002314:0.18278,9669.ENSMPUP00000005544:0.270727)0.6117:0.0396991)0.99985:0.0702148)0.99985:0.082488,(132908.ENSPVAP00000014833:0.488081,(9796.ENSECAP00000022144:0.310699,(((9785.ENSLAFP00000009512:0.187095,9813.ENSPCAP00000004417:0.493329)0.99985:0.359095,(30611.ENSOGAP00000016876:0.334272,(9483.ENSCJAP00000021314:0.178043,(9601.ENSPPYP00000003401:0.0415077,((61853.ENSNLEP00000003253:0.196659,9544.ENSMMUP00000037769:0.326984)0.835225:0.0989423,(9593.ENSGGOP00000004740:0.101826,9606.ENSP00000182290:0.0204981)0.997196:0.020731)0.307827:0.0046059)0.99985:0.0991112)0.99985:0.162323)0.972253:0.0380139)0.70642:0.0193389,((10141.ENSCPOP00000016274:0.272126,43179.ENSSTOP00000015376:0.458416)0.996119:0.0901785,(37347.ENSTBEP00000013312:0.328061,(10020.ENSDORP00000010739:0.398341,(10116.ENSRNOP00000051746:0.0455948,10090.ENSMUSP00000009396:0.0811741)0.99985:0.269525)0.791467:0.0577236)0.536676:0.0461933)0.99985:0.0620583)0.99985:0.0788824)0.969465:0.0395994)0.635969:0.0171601)0.702925:0.0283261)0.99985:0.11989);")

        trees = [(t1, "t1", True), (t2, "t2", False), (t3, "t3", True),
                 (t4, "t4", True), (t5, "t5", True), (t6, "t6", False),
                 (t7, "t7", True), (t8, "t8", True), (t9, "t9", False),
                 (t10, "t10", True)]
        for tree, tree_name, has_matches in trees:
            tree.set_species_naming_function(lambda n: n.name.split(".")[0] if "." in n.name else '')
            tree.annotate_ncbi_taxa()
            # Has support for two primates where at least one is not Homo sapiens
            pattern = """
                ( ' 9443 in @.lineage ' , ' 9443 in @.lineage and @.name!=9606 ' )' @.support >= 0.9 ';
                """
            pattern = TreePattern(pattern)
            if not has_matches:
                self.assertEqual(list(pattern.find_match(tree)), [])
            else:
                match = pattern.find_match(tree).next()
                self.assertEqual(match.support >= 0.9, True)
                test_status = (9443 in match.children[0].lineage and \
                               9443 in match.children[1].lineage and \
                               match.children[1].name != '9606')
                # permute children and check again
                test_status2 = (9443 in match.children[1].lineage and \
                               9443 in match.children[0].lineage and \
                               match.children[0].name != '9606')
                self.assertEqual(test_status, True)
                self.assertEqual(test_status2, True)



    def test_cached_attributes(self):
        pattern0 = """  '"Gallus_gallus_1" in leaves(@)' ;"""
        pattern1 = """( '"Hom" in species(@) and n_leaves(@) > 2')'"Pan_troglodytes_1" in leaves(@)';"""

        pattern0 = TreePattern(pattern0)
        pattern1 = TreePattern(pattern1)

        tree = PhyloTree(
            "((((Anolis_carolinensis_1:1, Gallus_gallus_1:1), (Felis_catus_1:1, (Homo_sapiens_1:1, Pan_troglodytes_1:1)primates)primates), ((Danio_rerio_1:1, (Xenopus_laevis_1:1, Anolis_carolinensis_1:1)), Saccharomyces_cerevisiae_2:1)), Saccharomyces_cerevisiae_1:1);",
            format=1)
        root = tree.get_tree_root()

        pattern0_match = list(pattern0.find_match(tree, maxhits=None))
        self.assertEqual(len(pattern0_match), 5)  # returns leaf itself
        self.assertEqual(pattern0_match[0], root)
        self.assertEqual(pattern0_match[4].name, "Gallus_gallus_1")

        pattern1_match = list(pattern1.find_match(tree, maxhits=None))
        self.assertEqual(len(pattern1_match), 3)
        self.assertEqual(pattern1_match[0], root)
        self.assertEqual(pattern1_match[2].children[1].children[1].children[0].name, "Homo_sapiens_1")

    def test_shortcut_functions(self):
        t = PhyloTree(
            """((((Human_1, Chimp_1), (Human_2, (Chimp_2, Chimp_3))),
            ((Fish_1, (Human_3, Fish_3)), Yeast_2)), Yeast_1);""")
        t.set_species_naming_function(lambda node: node.name.split("_")[0])
        t.get_descendant_evol_events()  # DDDSSSDDS

        root = t.get_tree_root()
        # Detects two consecutive nodes with duplications
        pattern0 = """('n_duplications(@) > 0')'n_duplications(@) > 0 '; """
        pattern1 = """( 'contains_leaves(@, ["Chimp_2", "Chimp_3"])'); """
        pattern2 = """'n_speciations(@) > 3 '; """

        pattern0 = TreePattern(pattern0)
        pattern1 = TreePattern(pattern1)
        pattern2 = TreePattern(pattern2)

        pattern0_match = list(pattern0.find_match(t, maxhits=None))
        pattern1_match = list(pattern1.find_match(t, maxhits=None))
        pattern2_match = list(pattern2.find_match(t, maxhits=None))

        self.assertEqual(len(pattern0_match), 5)

        self.assertEqual(len(pattern1_match), 4)
        self.assertEqual(pattern1_match[0], root)

        self.assertEqual(len(pattern2_match), 2)
        self.assertEqual(pattern2_match[0], root)
        self.assertEqual(pattern2_match[1], root.children[0])

    def regular_expression_syntax_test(self):
            '''

            The + represents one or more nodes.

            '''

            # should not match any pattern
            t1 = PhyloTree(""" ((c,g)a) ; """, format=8, quoted_node_names=False)

            # should not match any pattern
            # does not match pattern 1 because there must be at least one node between a and (c,d)
            t2 = PhyloTree(""" ((c,d)a) ; """, format=8, quoted_node_names=False)

            # should match patterns 1,2
            t3 = PhyloTree(""" ((d,c)b)a ; """, format=8, quoted_node_names=False)

            # should match patterns 1,2
            t4 = PhyloTree(""" ((c,d),(e,f)b)a ; """, format=8, quoted_node_names=False)

            # should match pattern 1,2
            t5 = PhyloTree(""" (((e,f)dum,(c,d)dee)b)a ; """, format=8, quoted_node_names=False)

            # should match only 1
            # does not match pattern 2 since (c,g) does not match (c,d)
            t6 = PhyloTree(""" (((e,f),(c,g)b)b)a ; """, format=8, quoted_node_names=False)

            # should match 1,2,3
            t7 = PhyloTree(""" (((e,f,g)d,(e,f,i)c)b)a ; """, format=8, quoted_node_names=False)

            # should match 1,2,3
            t8 = PhyloTree(""" (((e,f,i)d,(e,f,g)c)b)a ; """, format=8, quoted_node_names=False)

            # should match 1,2 not 3
            t9 = PhyloTree(""" (((e,f,i)d,(e,f,j)c)b)a ; """, format=8, quoted_node_names=False)

            # Should match 1,2,3
            # does not match pattern4 because ('e','f','g') should come from sibling of b
            t10 = PhyloTree(""" (b,((g,h,i)b,(e,f,g)c)d)a ; """, format=8, quoted_node_names=False)

            # should match 1,3,4
            # does not match pattern 2 because (c,c) does not match (c,d)
            t11 = PhyloTree("""  ( ((e, f, g) c) b, ((g, h, i)c) d) a ; """, format=8, quoted_node_names=False)

            pattern1 = TreePattern(""" ((c)+)a ;""", quoted_node_names=False)
            pattern2 = TreePattern(""" (('c','d')'+') 'a' ;""", quoted_node_names=True)
            pattern3 = TreePattern(""" (('e','f','g')'+') 'a' ;""", quoted_node_names=True)
            pattern4 = TreePattern(""" ((('g','h','i')+)'d',('e','f','g')'+') 'a' ;""", quoted_node_names=True)

            pattern1_match = [3, 4, 5, 6, 7, 8, 9, 10, 11]
            pattern2_match = [3, 4, 5, 7, 8, 9, 10]
            pattern3_match = [7, 8, 10, 11]
            pattern4_match = [11]
            true_match = [pattern1_match, pattern2_match, pattern3_match, pattern4_match]

            for p_num,pattern in enumerate([pattern1, pattern2, pattern3, pattern4]):
                pattern_match = []
                for tree_num, tree in enumerate([t1, t2, t3, t4, t5, t6, t7, t8, t9, t10, t11]):
                    if list(pattern.find_match(tree, maxhits=None)):
                        pattern_match += [tree_num+1]
                self.assertEqual(pattern_match, true_match[p_num])







def run():
    unittest.main()

if __name__ == '__main__':
    run()
