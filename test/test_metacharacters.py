import unittest

#assuming ete3 is not installed. they will be imported from local directory on developement modules.
from ete3 import PhyloTree, Tree, NCBITaxa
from treematcher import TreePattern


# plus symbol acts "one_or_more" matches
# trees from treematcher.py test function, except t6.

class Test_one_or_more_functionality_test(unittest.TestCase):
    def setUp(self):
        self.t1 = PhyloTree(""" ((c,g)a) ; """, format=8, quoted_node_names=False)
        self.t2 = PhyloTree(""" ((c,d)a) ; """, format=8, quoted_node_names=False)
        self.t3 = PhyloTree(""" ((d,c)b)a ; """, format=8, quoted_node_names=False)
        self.t4 = PhyloTree(""" ((c,d),(e,f)b)a ; """, format=8, quoted_node_names=False)
        self.t5 = PhyloTree(""" (((e,f)dum,(c,d)dee)b)a ; """, format=8, quoted_node_names=False)
        self.t6 = PhyloTree(""" ((d,c)a)a ; """, format=8, quoted_node_names=False)
        self.t11 = PhyloTree("""  ( ((e, f, g) c) b, ((g, h, i)c) d) a ; """, format=8, quoted_node_names=False)
        self.t12 = PhyloTree("""  ((( ((e, f, g) c) b, (((g, h, i)c)n) d)k)m) a ; """, format=8, quoted_node_names=False)

        self.pattern = TreePattern(""" ((c)+)a ;""", quoted_node_names=False)

    # there is no intermediate node between 'c' and 'a' nodes.
    def test_no_match(self):
        result = len(list(self.pattern.find_match(self.t1)))
        self.assertTrue(result == 0)

    # there is one intermediate node
    def test_one_intermediate_node(self):
        result = len(list(self.pattern.find_match(self.t3)))
        self.assertTrue(result > 0)

    #one intermediate node. pattern's 'c' (lowest node) has siter
    def test_one_intermediate_node_with_sister(self):
        result = len(list(self.pattern.find_match(self.t4)))
        self.assertTrue(result > 0)

    # two intermediate nodes. check the 'more' functionality
    def test_two_intermediate_nodes(self):
        result = len(list(self.pattern.find_match(self.t5)))
        self.assertTrue(result > 0)

    # one intermediate node. same name as the match.
    def test_one_intermediate_node_with_the_same_name(self):
        result = len(list(self.pattern.find_match(self.t6)))
        self.assertTrue(result == 1)

    # one intermediate node. intermediate node has siter node
    def test_one_intermediate_node_more_complex_tree(self):
        result = len(list(self.pattern.find_match(self.t11)))
        self.assertTrue(result > 0)

    # four intermediate nodes.
    def test_four_intermediate_nodes(self):
        result = len(list(self.pattern.find_match(self.t12)))
        self.assertTrue(result > 0)

class Test_zero_or_more_functionality(unittest.TestCase):
    def setUp(self):
        self.t1 = PhyloTree(""" ((c,g)a) ; """, format=8, quoted_node_names=False)
        self.t2 = PhyloTree(""" ((c,d)a) ; """, format=8, quoted_node_names=False)
        self.t3 = PhyloTree(""" ((d,c)b)a ; """, format=8, quoted_node_names=False)
        self.t4 = PhyloTree(""" ((c,d),(e,f)b)a ; """, format=8, quoted_node_names=False)
        self.t5 = PhyloTree(""" (((e,f)dum,(c,d)dee)b)a ; """, format=8, quoted_node_names=False)
        self.t6 = PhyloTree(""" ((d,c)a)a ; """, format=8, quoted_node_names=False)
        self.t11 = PhyloTree("""  ( ((e, f, g) c) b, ((g, h, i)c) d) a ; """, format=8, quoted_node_names=False)
        self.t12 = PhyloTree("""  ((( ((e, f, g) c) b, (((g, h, i)c)n) d)k)m) a ; """, format=8, quoted_node_names=False)
        self.t13 = PhyloTree("""  ( ((e, f, g) c) b, ((g, (w)h, i)c) d) a ; """, format=8, quoted_node_names=False)

        self.pattern  = TreePattern(""" ((c)*)a ;""", quoted_node_names=False)
        self.pattern2  = TreePattern(""" (((d)c)*)a ;""", quoted_node_names=False)
        self.pattern3 = TreePattern(""" ((((w)*)c)*)a ;""", quoted_node_names=False)
        self.pattern4 = TreePattern(""" ((d)*)a ;""")

    # there is no intermediate node between 'c' and 'a' nodes.
    def test_no_match(self):
        result = len(list(self.pattern.find_match(self.t1)))
        self.assertTrue(result > 0)

    # there is one intermediate node
    def test_one_intermediate_node(self):
        result = len(list(self.pattern.find_match(self.t3)))
        self.assertTrue(result > 0)

    #one intermediate node. pattern's 'c' (lowest node) has siter
    def test_one_intermediate_node_with_sister(self):
        result = len(list(self.pattern.find_match(self.t4)))
        self.assertTrue(result > 0)

    # two intermediate nodes. check the 'more' functionality
    def test_two_intermediate_nodes(self):
        result = len(list(self.pattern.find_match(self.t5)))
        self.assertTrue(result > 0)

    # one intermediate node. same name as the match.
    def test_one_intermediate_node_with_the_same_name(self):
        result = len(list(self.pattern.find_match(self.t6, maxhits=None)))
        self.assertTrue(result == 2)

    # one intermediate node. intermediate node has siter node
    def test_one_intermediate_node_more_complex_tree(self):
        result = len(list(self.pattern.find_match(self.t11)))
        self.assertTrue(result > 0)

    # four intermediate nodes.
    def test_four_intermediate_nodes(self):
        result = len(list(self.pattern.find_match(self.t12)))
        self.assertTrue(result > 0)

    def test_no_result(self):
        result = len(list(self.pattern2.find_match(self.t1)))
        self.assertTrue(result  == 0)

    def test_double_symbols(self):
        result = len(list(self.pattern3.find_match(self.t13)))
        self.assertTrue(result  > 0)

    def test_direct_connection(self):
        result = len(list(self.pattern3.find_match(self.t13)))
        self.assertTrue(result > 0)

class Test_leaf_and_root_metacharacters(unittest.TestCase):
    def setUp(self):
        self.t1 = PhyloTree(""" ((c,g)a) ; """, format=8, quoted_node_names=False)
        self.t2 = PhyloTree(""" ((c,d)a) ; """, format=8, quoted_node_names=False)
        self.t3 = PhyloTree(""" ((d,c)b)a ; """, format=8, quoted_node_names=False)
        self.t4 = PhyloTree(""" ((c,d),(e,f)b)a ; """, format=8, quoted_node_names=False)
        self.t5 = PhyloTree(""" (((e,f)dum,(c,d)dee)b)a ; """, format=8, quoted_node_names=False)
        self.t6 = PhyloTree(""" ((d,c)a)a ; """, format=8, quoted_node_names=False)
        self.t11 = PhyloTree("""  ( ((e, f, g) c) b, ((g, h, i)c) d) a ; """, format=8, quoted_node_names=False)
        self.t12 = PhyloTree("""  ((( ((e, f, g) c) b, (((g, h, i)c)n) d)k)m) a ; """, format=8, quoted_node_names=False)
        self.trees = [self.t1, self.t2, self.t3, self.t4, self.t5, self.t6, self.t11, self.t12]

    def test_simple_only_root(self):
        # should match all except t1 and t2 a is not root
        pattern = TreePattern(""" ((c)*)a^ ;""", quoted_node_names=False)
        count = 0
        for tree in self.trees:
            if any(True for result in pattern.find_match(tree) ):
                count += 1
        self.assertEqual(count, len(self.trees) - 2)

    def test_simple_only_leaf(self):
        # should match all except t11 and t12 'c' is not leaf
        pattern = TreePattern(""" ((c$)*)a ;""", quoted_node_names=False)
        count = 0
        for tree in self.trees:
            if any(True for result in pattern.find_match(tree) ):
                count += 1
        self.assertEqual(count, len(self.trees) - 2)

    def test_both_root_and_leaf(self):
        # should match all except t1, t2 'a' is not root, t11 and t12 'c' is not leaf
        pattern = TreePattern(""" ((c$)*)a^ ;""", quoted_node_names=False)
        count = 0
        for tree in self.trees:
            if any(True for result in pattern.find_match(tree) ):
                count += 1
        self.assertEqual(count, len(self.trees) - 4)

    def test_simple_false_root(self):
        # should not match any tree
        pattern = TreePattern(""" ((c^)*)a ;""", quoted_node_names=False)
        count = 0
        for tree in self.trees:
            if any(True for result in pattern.find_match(tree) ):
                count += 1
        self.assertEqual(count, 0)

    def test_simple_false_leaf(self):
        # should not match any tree
        pattern = TreePattern(""" ((c)*)a$ ;""", quoted_node_names=False)
        count = 0
        for tree in self.trees:
            if any(True for result in pattern.find_match(tree) ):
                count += 1
        self.assertEqual(count, 0)

    def test_false_root_true_leaf(self):
        # should not match any tree
        pattern = TreePattern(""" (((e$)c^)*)a ;""", quoted_node_names=False)
        count = 0
        for tree in self.trees:
            if any(True for result in pattern.find_match(tree) ):
                count += 1
        self.assertEqual(count, 0)

    def test_true_root_false_leaf(self):
        # should not match any tree
        pattern = TreePattern(""" (((e)c$)*)a^ ;""", quoted_node_names=False)
        count = 0
        for tree in self.trees:
            if any(True for result in pattern.find_match(tree) ):
                count += 1
        self.assertEqual(count, 0)

    def test_multiple_leaves(self):
         # should not match any tree
        pattern = TreePattern("""((( ((e$, f$, g$) c) b, (((g$, h$, i$)c)n) d)k)m) a^ ;;""", quoted_node_names=False)
        count = 0
        for tree in self.trees:
            if any(True for result in pattern.find_match(tree) ):
                count += 1
        self.assertEqual(count, 1)

class Test_controlled_skipping(unittest.TestCase):
    def setUp(self):
        self.t1 = PhyloTree(""" ((c,g)a) ; """, format=8, quoted_node_names=False)
        self.t2 = PhyloTree(""" ((c,d)a) ; """, format=8, quoted_node_names=False)
        self.t3 = PhyloTree(""" ((d,c)b)a ; """, format=8, quoted_node_names=False)
        self.t4 = PhyloTree(""" ((c,d),(e,f)b)a ; """, format=8, quoted_node_names=False)
        self.t5 = PhyloTree(""" (((e,f)dum,(c,d)dee)b)a ; """, format=8, quoted_node_names=False)
        self.t6 = PhyloTree(""" (((e,f),(c,g)b)b)a ; """, format=8, quoted_node_names=False)
        self.t7 = PhyloTree(""" (((e,f,g)d,(e,f,i)c)b)a ; """, format=8, quoted_node_names=False)
        self.t8 = PhyloTree(""" (((e,f,i)d,(e,f,g)c)b)a ; """, format=8, quoted_node_names=False)
        self.t9 = PhyloTree(""" (((e,f,i)d,(e,f,j)c)b)a ; """, format=8, quoted_node_names=False)
        self.t10 = PhyloTree(""" (b,((g,h,i)b,(e,f,g)c)d)a ; """, format=8, quoted_node_names=False)
        self.t11 = PhyloTree("""  ( ((e, f, g) c) b, ((g, h, i)c) d) a ; """, format=8, quoted_node_names=False)
        self.t12 = PhyloTree("""  ((( ((e, f, g) c) b, (((g, h, i)c)n) d)k)m) a ; """, format=8, quoted_node_names=False)
        self.t13 = PhyloTree(""" ((d,c)a)a ; """, format=8, quoted_node_names=False)

        self.trees = [self.t1, self.t2, self.t3, self.t4, self.t5, self.t6, self.t7, self.t8, self.t9, self.t10, self.t11, self.t12, self.t13]

    def compare(self, patterns, true_matches):
        matches = []
        for num, pattern in enumerate(patterns):
            current_matches = []
            for t_num, tree in enumerate(self.trees):
                if list(pattern.find_match(tree, maxhits=None)):
                    current_matches += [t_num + 1]
            matches += [current_matches]
        return  true_matches == matches

    def test_exact_number(self):
        pt1 = TreePattern(" ((c){0})a ; ", quoted_node_names=False)
        pt1_match = [1, 2, 13]
        pt2 = TreePattern(" ((c){1})a ; ", quoted_node_names=False)
        pt2_match = [3, 4, 7, 8, 9, 10, 11, 13]
        pt3 = TreePattern(" ((c){2})a ; ", quoted_node_names=False)
        pt3_match = [5, 6]
        pt4 = TreePattern(" ((c){4})a ; ", quoted_node_names=False)
        pt4_match = [12]

        patterns = [pt1, pt2, pt3, pt4]
        true_matches = [pt1_match, pt2_match, pt3_match, pt4_match]

        self.assertTrue(self.compare(patterns, true_matches))

    def test_only_low_bound(self):
        pt1 = TreePattern(" ((c){0-})a ;", quoted_node_names=False)
        pt1_match = range(1, 14) #all trees

        pt2 = TreePattern(" ((c){1-})a; ", quoted_node_names=False)
        pt2_match = range(3, 14) # all except 1, 2

        pt3 = TreePattern(" ((c){2-})a ;", quoted_node_names=False)
        pt3_match = [5, 6, 12]

        patterns = [pt1, pt2, pt3]
        true_matches = [pt1_match, pt2_match, pt3_match]

        self.assertTrue(self.compare(patterns, true_matches))

    def test_only_high_bound(self):
        pt1 = TreePattern(" ((c){-0})a; ", quoted_node_names=False)
        pt1_match = [1, 2, 13]

        pt2 = TreePattern(" ((c){-1})a; ", quoted_node_names=False)
        pt2_match = [1, 2, 3, 4, 7, 8, 9, 10, 11, 13]

        pt3 = TreePattern(" ((c){-2})a; ", quoted_node_names=False)
        pt3_match = range(1, 12) + [13] # all except 12

        pt4 = TreePattern(" ((c){-3})a; ", quoted_node_names=False)
        pt4_match = range(1,14)

        patterns = [pt1, pt2, pt3, pt4]
        true_matches = [pt1_match, pt2_match, pt3_match, pt4_match]

        self.assertTrue(self.compare(patterns, true_matches))

    def test_both_bounds(self):
        pt1 = TreePattern(" ((c){0-1})a; ", quoted_node_names=False)
        pt1_match = [1, 2, 3, 4, 7, 8, 9, 10, 11, 13]
        pt2 = TreePattern(" ((c){1-2})a; ", quoted_node_names=False)
        pt2_match = [3, 4, 5, 6, 7, 8, 9, 10, 11, 13]
        pt3 = TreePattern(" ((c){2-3})a; ", quoted_node_names=False)
        pt3_match = [5, 6, 12]
        pt4 = TreePattern(" ((c){3-4})a; ", quoted_node_names=False)
        pt4_match = [12]
        pt5 = TreePattern(" ((c){4-5})a; ", quoted_node_names=False)
        pt5_match = [12]

        patterns = [pt1, pt2, pt3, pt4, pt5]
        true_matches = [pt1_match, pt2_match, pt3_match, pt4_match, pt5_match]

        self.assertTrue(self.compare(patterns, true_matches))

class Test_zero_or_one_symbol(unittest.TestCase):
    def setUp(self):
        self.t1 = PhyloTree(""" ((c,g)a) ; """, format=8, quoted_node_names=False)
        self.t2 = PhyloTree(""" ((c,d)a) ; """, format=8, quoted_node_names=False)
        self.t3 = PhyloTree(""" ((d,c)b)a ; """, format=8, quoted_node_names=False)
        self.t4 = PhyloTree(""" ((c,d),(e,f)b)a ; """, format=8, quoted_node_names=False)
        self.t5 = PhyloTree(""" (((e,f)dum,(c,d)dee)b)a ; """, format=8, quoted_node_names=False)
        self.t6 = PhyloTree(""" (((e,f),(c,g)b)b)a ; """, format=8, quoted_node_names=False)
        self.t7 = PhyloTree(""" (((e,f,g)d,(e,f,i)c)b)a ; """, format=8, quoted_node_names=False)
        self.t8 = PhyloTree(""" (((e,f,i)d,(e,f,g)c)b)a ; """, format=8, quoted_node_names=False)
        self.t9 = PhyloTree(""" (((e,f,i)d,(e,f,j)c)b)a ; """, format=8, quoted_node_names=False)
        self.t10 = PhyloTree(""" (b,((g,h,i)b,(e,f,g)c)d)a ; """, format=8, quoted_node_names=False)
        self.t11 = PhyloTree("""  ( ((e, f, g) c) b, ((g, h, i)c) d) a ; """, format=8, quoted_node_names=False)
        self.t12 = PhyloTree("""  ((( ((e, f, g) c) b, (((g, h, i)c)n) d)k)m) a ; """, format=8, quoted_node_names=False)
        self.t13 = PhyloTree(""" ((d,c)a)a ; """, format=8, quoted_node_names=False)

        self.trees = [self.t1, self.t2, self.t3, self.t4, self.t5, self.t6, self.t7, self.t8, self.t9, self.t10, self.t11, self.t12, self.t13]

    def test_symbol(self):
        pt1 = TreePattern(" ((c)?)a; ", quoted_node_names=False)
        pt1_match = [1, 2, 3, 4, 7, 8, 9, 10, 11, 13]

        patterns = [pt1]
        true_matches = [pt1_match]

        matches = []
        for num, pattern in enumerate(patterns):
            current_matches = []
            for t_num, tree in enumerate(self.trees):
                if list(pattern.find_match(tree, maxhits=None)):
                    current_matches += [t_num + 1]
            matches += [current_matches]
        self.assertEqual(true_matches,  matches)

if __name__ == '__main__':
    unittest.main()
