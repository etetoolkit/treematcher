import unittest
from ete3 import  Tree
from treematcher.treematcher import TreePattern
from copy import deepcopy
#class Test_strict_match():
class Test_strict_match(unittest.TestCase):

    def setUp(self):
        t1 = Tree(""" ((c,g)a) ; """, format=8, quoted_node_names=False)
        t2 = Tree(""" ((c,d)a) ; """, format=8, quoted_node_names=False)
        t3 = Tree(""" ((d,c)b)a ; """, format=8, quoted_node_names=False)
        t4 = Tree(""" ((c,d),(e,f)b)a ; """, format=8, quoted_node_names=False)
        t5 = Tree(""" (((e,f)dum,(c,d)dee)b)a ; """, format=8, quoted_node_names=False)
        t6 = Tree(""" (((e,f),(c,g)b)b)a ; """, format=8, quoted_node_names=False)
        t7 = Tree(""" (((e,f,g)d,(e,f,i)c)b)a ; """, format=8, quoted_node_names=False)
        t8 = Tree(""" (((e,f,i)d,(e,f,g)c)b)a ; """, format=8, quoted_node_names=False)
        t9 = Tree(""" (((e,f,i)d,(e,f,j)c)b)a ; """, format=8, quoted_node_names=False)
        t10 = Tree(""" (b,((g,h,i)b,(e,f,g)c)d)a ; """, format=8, quoted_node_names=False)
        t11 = Tree("""  ( ((e, f, g) c) b, ((g, (w)h, i)c) d) a ; """, format=8, quoted_node_names=False)

        self.trees = [t1, t2, t3, t4, t5, t6, t7, t8, t9, t10, t11]

    def test_one_terminal_node(self):
        # The c node exists as leaf (not internal)
        pattern = TreePattern(" (c)^; ")
        true_match = [1, 2, 3, 4, 5, 6]
        matches = []
        for num, tree in enumerate(self.trees):

            one_use = deepcopy(pattern)

            result = one_use.find_match(tree)
            try:
                res = next(result)
            except:
                res = None

            if res:
                matches += [num + 1]

        self.assertTrue(matches == true_match)

    def test_two_terminal_nodes(self):
        # The presense of leaves e, f as sister nodes.
        pattern = TreePattern(" (e, f)^; ")
        true_match = [4, 5, 6, 7, 8, 9, 10, 11]
        matches = []
        for num, tree in enumerate(self.trees):
            one_use = deepcopy(pattern)

            result = one_use.find_match(tree)
            try:
                res = next(result)
            except:
                res = None

            if res:
                matches += [num + 1]

        self.assertTrue(matches == true_match)

    def test_simple_paren_two_children(self):
        tree = Tree(" (((b, c)a, (b, c)a), (e, f)d) ;", format=1)
        pattern = TreePattern("(b,c)a ;")
        result = pattern.find_match(tree)
        self.assertTrue(len(list(result)) == 2)

    def test_simple_parent_two_children_false(self):
        tree = Tree(" (((b, c)a, (b, c)a), (e, f)d) ;", format=1)
        pattern = TreePattern("(b,c)qq ;")
        result = pattern.find_match(tree)
        self.assertTrue(len(list(result)) == 0)

    def test_simple_complete_topology(self):
        pattern = TreePattern("((e, i, f)d)^ ; ")
        true_match = [8, 9]
        match = []

        for num, tree in enumerate(self.trees):
            result = pattern.find_match(tree)
            if (len(list(result)) > 0):
                match += [num+1]

        self.assertTrue(true_match == match)

    # Waits for update
    # def test_incomplete_topology(self):
    #    pattern = TreePattern(" (e, i, f)d ;")
    #
    #    for i in pattern.find_match(self.trees[8]):
    #        print i

    #    self.assertTrue(True)

class Test_metacharacters_at_terminal_nodes(unittest.TestCase):
    def test_simple_plus(self):
        tree = Tree(" (((a, a, b, qq), (a, b, c, ww)), (b, b, a, ee));", format=8)
        pattern = TreePattern(" (qq, a+)^ ;")

        result = pattern.find_match(tree)
        expected = (tree&'qq').up
        self.assertTrue(len(list(result)) > 0 )

    def test_double_match(self):
        tree = Tree(" (((a, a, b), (c, c, d) ), (e, e, f), (g, h, i)) ; ")
        pattern = TreePattern( " ((a+, b)^, (e+, f)^);")

        result = (pattern.find_match(tree))

        #self.assertTrue(len(list(result)) > 0 )
        self.assertEqual(next(result), tree)

    def test_simple_zero_or_more(self):
        tree = Tree(" ((((a, a, b)), (c, d), (e, f)), (g, h, i)) ; ")
        pattern = TreePattern(" (a, a, b, d*) ;")

        result = (pattern.find_match(tree))
        self.assertTrue(len(list(result)) > 0)

    def test_skipped_zero_or_more(self):
        tree = Tree(" ((((a, a, b)), (c, d), (e, f)), (g, h, i)) ; ")
        pattern = TreePattern(" ( a, b, d*) ;")

        result = (pattern.find_match(tree))
        # false test
        self.assertTrue(len(list(result)) == 0)

    def test_constraintes_pattern(self):
        tree = Tree(" ((((a, b)), (c, d), (e, f)), (g, h, i)) ; ")

        (tree&'a').dist = 0.2
        (tree&'b').dist = 0.4
        (tree&'c').dist = 0.5
        (tree&'d').dist = 0.6
        (tree&'e').dist = 0.7
        (tree&'f').dist = 0.8
        (tree&'g').dist = 0.9


        pattern = TreePattern(" ('@.dist > 0.5+'); ", quoted_node_names=True)

        result = pattern.find_match(tree)
        expected = [(tree&'e').up, (tree&'g').up]

        found = True
        count = 0

        for node in result:
            count += 1
            found &= node in expected

        found &= count == 2

        self.assertTrue(found)

    def test_constraints_and_loose(self):
        tree = Tree(" (((((a, b), (c, d), (e, f)), (g, h, i)))) ; ")

        (tree&'a').dist = 0.2
        (tree&'b').dist = 0.4
        (tree&'c').dist = 0.5
        (tree&'d').dist = 0.6
        (tree&'e').dist = 0.7
        (tree&'f').dist = 0.7
        (tree&'g').dist = 0.9

        pattern = TreePattern(""" ('@.dist == 0.2', 'b')'^', ('@.dist > 0.5', '@.dist == 0.7+')'^' ; """, quoted_node_names=True)
        result = pattern.find_match(tree)
        res = next(result)
        #self.assertEqual(res, ((tree&'f').up).up)
        self.assertTrue( len(list(result)) > 0 )

    def test_star_and_logical_constraints(self):
        tree = Tree(" (((a, b), (c, d), (e, f)), (g, h, i)) ; ")

        (tree&'a').dist = 0.2
        (tree&'b').dist = 0.4
        (tree&'c').dist = 0.5
        (tree&'d').dist = 0.6
        (tree&'e').dist = 0.7
        (tree&'f').dist = 0.8
        (tree&'g').dist = 0.9

        pattern = TreePattern(""" ('g', '@.dist == 1+', 'fls_node*'); """, quoted_node_names=True)
        result = pattern.find_match(tree)

        self.assertEqual(next(result), (tree&'g').up)

    def test_exact_number_of_repeat(self):
        tree = Tree("((a, a, a, b, c), (d, d, qq), (e, e, e, ww, e, e, e, e, e)); ")
        p1 = TreePattern(" (b, c, 'a{1,3}') ;")
        p2 = TreePattern(" (b, c, 'a{2,3}') ;")
        p3 = TreePattern(" (b, c, 'a{3,3}') ;")
        p4 = TreePattern(" (b, c, 'a{4,5}') ;")
        p5 = TreePattern(" (ww, 'e{1,8}') ;")
        p6 = TreePattern(" (ww, 'e{7,9}') ;")
        p7 = TreePattern(" (ww, 'e{1,3}') ;")

        patterns = [p1, p2, p3, p4, p5, p6, p7]
        true_match = [True, True, True, False, True, True, False]
        match = True

        for num, pattern in enumerate(patterns):
            result = pattern.find_match(tree)
            found = len(list(result)) > 0
            match &= found == true_match[num]

        self.assertTrue(match)

    def test_exact_number_and_topology(self):
        tree = Tree(" ((a, a, b)p1, ((c, c, c, d)p2, (e, f, g)p3)p4)p5 ;", format=1)
        p1 = TreePattern(" ('a{2,2}', 'b')'p1' ;", quoted_node_names=True)
        p2 = TreePattern(" ('c{1,5}', 'd')'p2' ;",quoted_node_names=True)
        p3 = TreePattern(" ('c{2,3}', d, 'ww{0,3}')p2 ;")
        p4 = TreePattern(" ('c{3,3}', 'd{0,5}', 'ww{0,3}')p2;")
        p5 = TreePattern(" ('c{1,2}', 'd{0,1}', 'ww*')p2;")

        patterns = [p1, p2, p3, p4, p5]
        true_match = [True, True, True, True, False]
        match = True

        for num, pattern in enumerate(patterns):
            result = pattern.find_match(tree)
            found = len(list(result)) > 0
            match &= found == true_match[num]

        self.assertTrue(match)


class Test_basic_tests(unittest.TestCase):
    def test_all(self):

        test = True

        t1 = Tree(" ((a, a, b)p1, ((c, c, c, d)p2, (e, f, g)p3)p4)p5 ;", format=1)
        p1 = TreePattern(" ('c+', 'd')'p2' ;",quoted_node_names=True)
        test &= len(list(p1.find_match(t1))) > 0


        # Should  match
        t1 = Tree(" (((F, G)E, (C, D)B), A);", format=8)
        p1  = TreePattern("('@.support > 0', '@.support > 0')'B' ;")
        test &= len(list(p1.find_match(t1))) > 0



        # Should NOT match
        t1 = Tree(" (((F, G)E, (C, D)B), A);", format=8)
        p1  = TreePattern("('@.support > 0', '@.support > 0{2,3}')'B' ;")
        test &= len(list(p1.find_match(t1))) == 0



        # Should  match
        t1 = Tree(" (((F, G)E, (C, D)B), A);", format=8)
        p1  = TreePattern("('C', '@.support > 0')'B' ;")
        test &= len(list(p1.find_match(t1))) > 0


        # Should not match
        t1 = Tree("(((A, A, A), (B,C)), K);")
        p1 = TreePattern("(((A, A+, A, A), (B,C)), K);")
        test &= len(list(p1.find_match(t1))) == 0


        # Should match
        t1 = Tree("(((A, A, A), (B,C)), K);")
        p1 = TreePattern("(((A, A+, A), (B,C)), K);")
        test &= len(list(p1.find_match(t1))) > 0


        # Should match
        t1 = Tree("(((A, A, A), (B,C)), K);")
        p1 = TreePattern("(((A, A+), (B,C)), K);")
        test &= len(list(p1.find_match(t1))) > 0


        # ^ after a ) means that the two children of that node can be connected by
        # any number of internal up/down nodes
        t1 = Tree("(  ((B,Z), (D,F)), G);")
        p1 = TreePattern("( (B,Z), G)^;")
        test &= len(list(p1.find_match(t1))) > 0


        t1 = Tree("(  ((G, ((B,Z),A)), (D,G)), C);")
        p1 = TreePattern("(((B,Z)^,C), G)^;")
        test &= len(list(p1.find_match(t1))) == 0


        t1 = Tree("(  ((G, ((B,Z),A)), (D,G)), C);")
        p1 = TreePattern("(((B,Z)^,G), C)^;")
        test &= len(list(p1.find_match(t1))) > 0


        t1 = Tree("(((A, (B,C,D)), ((B,C), A)), F);")
        p1 = TreePattern("((C,B,D*), A);")
        test &= len(list(p1.find_match(t1))) > 0


        t1 = Tree("(((A, (B,C,D, D, D)), ((B,C), A)), F);")
        p1 = TreePattern("((C,B,'D{2,3}'), A);")
        test &= len(list(p1.find_match(t1))) > 0


        t1 = Tree("(a, b, b, a);")
        p1 = TreePattern("(a+, b+);")
        test &= len(list(p1.find_match(t1))) > 0


        t1 = Tree("((a, b), c);")
        p1 = TreePattern("((a, b, d*), c);")
        test &= len(list(p1.find_match(t1))) > 0


        t1 = Tree("(  (((B,H), (B,B,H), C), A), (K, J));")
        p1 = TreePattern("((C, (B+,H)+), A);")
        test &= len(list(p1.find_match(t1))) > 0


        self.assertTrue(test)


if __name__ == '__main__':
    unittest.main()
