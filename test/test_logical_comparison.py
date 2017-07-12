import unittest
import random
import sys

#assuming ete3 is not installed. they will be imported from local directory on developement modules.
from ete3 import PhyloTree, Tree, NCBITaxa
from treematcher import TreePattern

class test_sets_logical_comparison(unittest.TestCase):
    def setUp(self):
        if len(sys.argv) > 1 and sys.argv[1]:
          self.number_of_itteration = int(sys.argv[1])
        else:
          self.number_of_itteration = 1000

        if len(sys.argv) > 2 and sys.argv[2]:
          self.number_of_leaves = int(sys.argv[2])
        else:
          self.number_of_leaves = 30

    def logical_comparison_distance(self, treePattern, maxhits=None):
        # define a tree
        a = Tree()
        a.populate(self.number_of_leaves)

        # find the number of nodes
        length = 0
        for node in a.traverse():
            length += 1

        # create a list of defined random distances, so the max can be verified.
        maximum = 1000.0 # allow a bigger possibility of mistake, by adding a lot of values.
        dists = []
        for i in range(0, length):
            dists += [round(random.uniform(0.1, maximum), 1)]

        # set nodes distances to the nodes
        cursor = 0
        for node in a.traverse():
            node.dist = dists[cursor]
            cursor += 1

        matched = []
        not_matched = []

        result = list(treePattern.find_match(a, maxhits=None))

        for node in a.traverse():
            if node in result:
                matched += [node]
            else:
                not_matched += [node]

        return [matched, not_matched]


    def test_less_than_children(self):
        pattern = TreePattern(""" '@.dist < [:children:].dist' ; """, quoted_node_names=True)
        test_matches = True
        for i in range(1,self.number_of_itteration+1):
            result = self.logical_comparison_distance(pattern)
            matched = result[0]
            not_matched = result[1]

            for node in matched:
                test_matches &= all(node.dist < child.dist for child in node.children)
            for node in not_matched:
                test_child = True
                for child in node.children:
                    test_child &= node.dist < child.dist
                test_matches &= not test_child

        self.assertTrue(test_matches)

    def test_more_than_children(self):
        pattern = TreePattern(""" '@.dist > [:children:].dist' ; """, quoted_node_names=True)
        test_matches = True
        for i in range(1,self.number_of_itteration+1):
            result = self.logical_comparison_distance(pattern)
            matched = result[0]
            not_matched = result[1]

            for node in matched:
                test_matches &= all(node.dist > child.dist for child in node.children)
            for node in not_matched:
                test_child = True
                for child in node.children:
                    test_child &= node.dist > child.dist
                test_matches &= not test_child

        self.assertTrue(test_matches)

    def test_more_than_any_child(self):
        pattern = TreePattern(""" '@.dist > [:any_child:].dist' ; """, quoted_node_names=True)
        test_matches = True
        for i in range(1,self.number_of_itteration+1):
            result = self.logical_comparison_distance(pattern)
            matched = result[0]
            not_matched = result[1]

            for node in matched:
                test_matches &= any(node.dist > child.dist for child in node.children)
            for node in not_matched:
                test_matches &= all(node.dist <= child.dist for child in node.children)

        self.assertTrue(test_matches)

    def test_less_than_any_child(self):
        pattern =  TreePattern(""" '@.dist < [:any_child:].dist' ; """, quoted_node_names=True)
        test_matches = True
        for i in range(1,self.number_of_itteration+1):
            result = self.logical_comparison_distance(pattern)
            matched = result[0]
            not_matched = result[1]

            for node in matched:
                test_matches &= any(node.dist < child.dist for child in node.children)
            for node in not_matched:
                test_matches &= all(node.dist >= child.dist for child in node.children)

        self.assertTrue(test_matches)


class test_exrteme_cases(unittest.TestCase):
    def setUp(self):
        if len(sys.argv) > 1 and sys.argv[1]:
          self.number_of_itteration = int(sys.argv[1])
        else:
          self.number_of_itteration = 1000

        if len(sys.argv) > 2 and sys.argv[2]:
          self.number_of_leaves = int(sys.argv[2])
        else:
          self.number_of_leaves = 30

    def compare_extreme_values(self):
        # define a tree
        a = Tree()
        a.populate(self.number_of_leaves)

        # find the number of nodes
        length = 0
        for node in a.traverse():
            length += 1

        # create a list of defined random distances, so the max can be verified.
        maximum = 1000.0 # allow a bigger possibility of mistake, by adding a lot of values.
        dists = []
        for i in range(0, length):
            dists += [round(random.uniform(0.1, maximum), 1)]

        # set nodes distances to the nodes
        cursor = 0
        for node in a.traverse():
            node.dist = dists[cursor]
            cursor += 1

        var_max = max(dists)
        var_min = min(dists) # The correct mimimum value should be this using python function.

        # Define the patterns and run the match.
        b = TreePattern(""" '@.dist > [:all_nodes:].dist' ; """, quoted_node_names=True)
        c = TreePattern(""" '@.dist < [:all_nodes:].dist' ; """, quoted_node_names=True)

        for i in b.find_match(a): tree_max = i.dist
        for i in c.find_match(a): tree_min = i.dist

        # Retrieve again all the values and find the max again for verification.
        dists_sec = []
        for node in a.traverse():
            dists_sec += [node.dist]
            again_max = max(dists_sec)
            again_min = min(dists_sec)

        if not sum(dists) == sum(dists_sec):
            return 2

        # get results
        if var_max == tree_max == again_max and var_min == tree_min == again_min:
            return 0
        else:
            return 1

    def test_exrteme_cases(self):
        num = self.compare_extreme_values()
        correct = 0
        fails = 0

        for i in range(1, self.number_of_itteration+1):
            num = self.compare_extreme_values()
            if num == 0: correct += 1
            else: fails += 1

        self.assertEqual(correct, self.number_of_itteration)


if __name__ == '__main__':
    unittest.main()
