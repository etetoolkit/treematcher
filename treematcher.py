#system modules
import re
from itertools import permutations

# third party modules
import ast
import six
import copy
from ete3 import PhyloTree, Tree, NCBITaxa
from symbols import SYMBOL, SET

# internal modules
#...


class TreePatternCache(object):
    def __init__(self, tree):
        """ Creates a cache for attributes that require multiple tree
        traversal when using complex TreePattern queries.

        :param tree: a regular ETE tree instance
         """
        # Initialize cache (add more stuff as needed)
        self.leaves_cache = tree.get_cached_content()
        self.all_node_cache = tree.get_cached_content(leaves_only=False)

    def get_cached_attr(self, attr_name, node, leaves_only=False):
        """
        easy access to cached attributes on trees.

        :param attr_name: any attribute cached in tree nodes (e.g., species,
         name, dist, support, etc.)

        :param node: The pattern tree containing the cache

        :leaves_only: If True, return cached values only from leaves

        :return: cached values for the requested attribute (e.g., Homo sapiens,
         Human, 1.0, etc.)

        """
        #print("USING CACHE")
        cache = self.leaves_cache if leaves_only else self.all_node_cache
        values = [getattr(n, attr_name, None) for n in cache[node]]
        return values

    def get_leaves(self, node):
        return self.leaves_cache[node]

    def get_descendants(self, node):
        return self.all_node_cache[node]


class _FakeCache(object):
    """TreePattern cache emulator."""
    def __init__(self):
        pass

    def get_cached_attr(self, attr_name, node, leaves_only=False):
        """ Helper function to mimic the behaviour of a cache, so functions can
        refer to a cache even when one has not been created, thus simplifying code
        writing. """
        if leaves_only:
            iter_nodes = node.iter_leaves
        else:
            iter_nodes = node.traverse

        values = [getattr(n, attr_name, None) for n in iter_nodes()]
        return values

    def get_leaves(self, node):
        return node.get_leaves()

    def get_descendants(self, node):
        return node.get_descendants()


class PatternSyntax(object):
    def __init__(self):
        # Creates a fake cache to ensure all functions below are functioning
        # event if no real cache is provided
        self.__fake_cache = _FakeCache()
        self.__cache = None

    def __get_cache(self):
        if self.__cache:
            return self.__cache
        else:
            return self.__fake_cache

    def __set_cache(self, value):
        self.__cache = value

    cache = property(__get_cache, __set_cache)

    def leaves(self, target_node):
        return sorted([name for name in self.cache.get_cached_attr(
            'name', target_node, leaves_only=True)])

    def descendants(self, target_node):
        return sorted([name for name in self.cache.get_cached_attr(
            'name', target_node)])

    def species(self, target_node):
        return set([name for name in self.cache.get_cached_attr(
            'species', target_node, leaves_only=True)])

    def contains_species(self, target_node, species_names):
        """
        Shortcut function to find the species at a node and any of it's descendants.
        """
        if isinstance(species_names, six.string_types):
            species_names = set([species_names])
        else:
            species_names = set(species_names)

        found = 0
        for sp in self.cache.get_cached_attr('species', target_node, leaves_only=True):
            if sp in species_names:
                found += 1
        return found == len(species_names)

    def contains_leaves(self, target_node, node_names):
        """ Shortcut function to find if a node contains at least one of the
        node names provided. """

        if isinstance(node_names, six.string_types):
            node_names = set([node_names])
        else:
            node_names = set(node_names)

        found = 0
        for name in self.cache.get_cached_attr('name', target_node, leaves_only=True):
            if name in node_names:
                found += 1
        return found == len(node_names)

    def n_species(self, target_node):
        """ Shortcut function to find the number of species within a node and
        any of it's descendants. """

        species = self.cache.get_cached_attr('species', target_node, leaves_only=True)
        return len(set(species))

    def n_leaves(self, target_node):
        """ Shortcut function to find the number of leaves within a node and any
                of it's descendants. """
        return len(self.cache.get_leaves(target_node))

    def n_duplications(self, target_node):
        """
            Shortcut function to find the number of duplication events at or below a node.
            :param target_node: Node to be evaluated, given as @.
            :return: True if node is a duplication, otherwise False.
        """
        events = self.cache.get_cached_attr('evoltype', target_node)
        return(events.count('D'))

    def n_speciations(self, target_node):
        """
            Shortcut function to find the number of speciation events at or below a node.
        """
        events = self.cache.get_cached_attr('evoltype', target_node)
        return(events.count('S'))


    # TODO: review next function

    def smart_lineage(self, constraint):
        """ Get names instead of tax ids if a string is given before the "in
        @.linage" in a query. Otherwise, returns Taxonomy ids. Function also
        works for constraint that contains something besides the given target
        node (e.g., @.children[0].lineage).

        :param constraint: Internal use.

        :return:  Returns list of lineage tax ids if taxid is searched,
        otherwise returns names in lineage. """

        parsedPattern = ast.parse(constraint, mode='eval')

        lineage_node = [n for n in ast.walk(parsedPattern)
                        if hasattr(n, 'comparators') and type(n.comparators[0]) == ast.Attribute
                        and n.comparators[0].attr == "lineage"]

        index = 0
        for lineage_search in lineage_node:
            if hasattr(lineage_node[index].left,'s'):
                # retrieve what is between __target and .lineage
                found_target = (re.search(r'__target[^ ]*\.lineage', constraint).span())
                extracted_target = constraint[found_target[0]: found_target[1]]

                syntax = "(ncbi.get_taxid_translator(" + \
                         str(extracted_target) + ")).values()"
                if index == 0:
                    constraint = constraint.replace(str(extracted_target), syntax, 1)
                else:

                    constraint = re.sub(r'^((.*?' + extracted_target + r'.*?){' + str(index) + r'})' + extracted_target,
                             r'\1' + syntax, constraint)

            index += 1

        return constraint


class TreePattern(Tree):
    def __str__(self):
        return self.get_ascii(show_internal=True, attributes=["name"])

    def __init__(self, newick=None, format=1, dist=None, support=None,
                 name=None, quoted_node_names=True, syntax=None):
        """ Creates a tree pattern instance that can be used to search within
        other trees.

        :param newick: Path to the file containing the tree or, alternatively,
            the text string containing the same information.
        :param format: subnewick format that is a number between 0 (flexible)
            and 9 (leaf name only), or 100 (topology only).
        :param quoted_node_names: Set to True if node names are quoted with
            stings, otherwise False.
        :syntax: Syntax controller class containing functions to use as
            constraints within the pattern.
        """
        # Load the pattern string as a normal ETE tree, where node names are
        # python expressions
        super(TreePattern, self).__init__(newick, format, dist, support, name, quoted_node_names)

        # Set a default syntax controller if a custom one is not provided
        self.syntax = syntax if syntax else PatternSyntax()

        # check the tree syntax.

    def is_in_bounds(self, bound):
        """
        Checks in case of indirect connection the current node is inside the bounds.
        :param bound: specifies what bound to test
        """

        if bound == 'low':
            return self.controller["skipped"] >= self.controller["low"]
        elif bound == 'high':
            if self.controller["high"] >= 0:
                # high bound has to be less but not equal, because it performs a 'can skip one more' test.
                return self.controller["skipped"] < self.controller["high"]
        return True

    def find_extreme_case(self, match_list):
        """
        Given a single match pattern and a list of (matched) nodes
        returns the one that fulfills the constraint. That constraint should
        describe an axtreme condition such as minimum or maxhimum.

        :param match_list: a list of (matched) nodes
        """

        # find the node in self that has the comparison expression.
        # now is useless. In case this comparison can be made inside other
        # patterns, the every pattern should compare their corresponding node of
        # the other patterns.
        for node in self.traverse():
            if node.controller["single_match"]:
                pattern = node
                constraint = pattern.controller["single_match_contstraint"]

        # for all nodes in the matched pattern find the one that fits the best
        # the constraint expression.
        # assumes that the root has to be tested.
        correct_node = match_list[0]

        st = False # truth flag for eval()

        for node in match_list:
            # for every node in the list compare the best match so far with
            # the current node in the list priority.
            constraint = pattern.controller["single_match_contstraint"]
            constraint_scope = {attr_name: getattr(self.syntax, attr_name)
                                for attr_name in dir(self.syntax)}
            constraint_scope.update({"__target_node": node})
            constraint_scope.update({"__correct_node": correct_node})
            constraint = constraint.replace("@", "__target_node")
            constraint = constraint.replace(SET["all_nodes"], "__correct_node")

            try:
                st = eval(constraint, constraint_scope)

            except ValueError:
                raise ValueError("not a boolean result: . Check quoted_node_names.")

            except (AttributeError, IndexError) as err:
                raise ValueError('Constraint evaluation failed at %s: %s' %
                                 (target_node, err))
            except NameError:
                try:
                    # temporary fix. Can not access custom syntax on all nodes. Get it from the root node.
                    root_syntax = self.get_tree_root().syntax
                    constraint_scope = {attr_name: getattr(root_syntax, attr_name)
                                        for attr_name in dir(root_syntax)}
                    constraint_scope.update({"__target_node": node})
                    constraint_scope.update({"__correct_node": correct_node})

                    st = eval(constraint, constraint_scope)
                    if st: correct_node = node

                except NameError as err:
                    raise NameError('Constraint evaluation failed at %s: %s' %
                             (target_node, err))
            else:
                if st: correct_node = node
        return correct_node

    def decode_repeat_symbol(self, bounds):
        """
        Extracts valuable information from the controlled skipping case.

        :param bounds: a string that contains the {x-y} pattern.

        returns a list with the lower and the highest bounds in respectively order.
        """
        bounds = bounds[0:len(bounds)-1]
        if '-' in bounds:
            split = bounds.split("-")
            low = int(split[0]) if split[0] else 0
            high = int(split[1] ) if split[1] else -1
        else:
            low = high = int(bounds)
        return [low, high]

    def set_controller(self):
        """
        Creates a dictionary that contains information about a node.
        That information is about how a node interacts with the tree topology.
        It describes how the metacharacter connects with the rest of nodes and
        if it is leaf or root.
        """
        controller = {}

        # bounds and (already) skipped nodes values
        controller["low"] = 0
        controller["high"] = -1
        controller["skipped"] = 0
        controller["single_match"] = False

        # update controller in case of root or leaf metacharacters and set the node name
        if SYMBOL["is_root"] in self.name:
            controller["root"] = True
            self.name = self.name.split(SYMBOL["is_root"])[0]
        if SYMBOL["is_leaf"] in self.name:
            controller["leaf"] = True
            self.name = self.name.split(SYMBOL["is_leaf"])[0]

        # transform sets to the corresponding code
        if SET["any_child"] in self.name:
            self.name = " any( " + self.name.split("[")[0] + " " + ("[" + self.name.split("[")[1]).replace(SET["any_child"], "x") + " for x in __target_node.children)"
        if SET["children"] in self.name:
            self.name = " all( " + self.name.split("[")[0] + " " + ("[" + self.name.split("[")[1]).replace(SET["children"], "x") + " for x in __target_node.children)"
        if SET["all_nodes"] in self.name:
            controller["single_match"] = True
            controller["single_match_contstraint"] = self.name
            self.name = '@'

        # update controller according to metacharacter connection properties.
        if self.name == SYMBOL["one_or_more"]:
            controller["allow_indirect_connection"] = True
            controller["direct_connection_first"] = False
        elif self.name == SYMBOL["zero_or_more"]:
            controller["allow_indirect_connection"] = True
            controller["direct_connection_first"] = True
        elif self.name == SYMBOL["zero_or_one"]:
            controller["direct_connection_first"] = True
            controller["allow_indirect_connection"] = True
            controller["high"] = 1
        elif SYMBOL["defined_number_set"] in self.name:
            split = self.name.split(SYMBOL["defined_number_set"])
            self.name = split[0]
            bounds = self.decode_repeat_symbol(split[1])
            controller["low"] = bounds[0]
            controller["high"] = bounds[1]
            controller["allow_indirect_connection"] = True
            if controller["low"]  == 0: controller["direct_connection_first"] = True
            else: controller["direct_connection_first"] = False
        else:
            controller["allow_indirect_connection"] = False
            controller["direct_connection_first"] = False

        self.controller = controller
        return self.controller["single_match"]

    # FUNCTIONS EXPOSED TO USERS START HERE
    def match(self, node, cache=None):
        """
        Check all constraints interactively on the target node.

        :param node: A tree (node) to be searched for a given pattern.

        :param local_vars:  Dictionary of TreePattern class variables and
        functions for constraint evaluation.

        :return: True if a match has been found, otherwise False.
        """
        self.syntax.cache = cache
        # does the target node match the root node of the pattern?

        #check the zero intermediate node case.
        #assumes that SYMBOL["zero_or_more"] has only one child.
        if self.controller["direct_connection_first"]:
            self = self.children[0]

        status = self.is_local_match(node, cache)

        if not status:
            if self.up is not None and self.up.controller["allow_indirect_connection"] and self.up.is_in_bounds("high"):  # skip node by resetting pattern
                status = True
                self = self.up
                self.controller["skipped"] += 1
        elif self.controller["allow_indirect_connection"] and self.controller["skipped"] == 0:
            self.controller["skipped"] += 1


        # if so, continues evaluating children pattern nodes against target node
        # children
        if status and self.children:

                #if the number of children do not match, find where they do and check that
                nodes = []

                if len(node.children) < len(self.children):
                    if self.controller["allow_indirect_connection"]:
                        count = 0
                        for skip_to_node in node.traverse(strategy="levelorder"):
                            # skip to node with correct number of children
                            if len(skip_to_node.children) >= len(self.children):
                                count += 1
                                nodes += [skip_to_node]
                                sisters = skip_to_node.get_sisters()
                                if len(sisters) > 0:
                                   for sister in sisters:
                                       nodes += [sister]

                                break
                        if count < 1:
                            status = False

                    else:
                        #print("setting status to false")
                        status = False

                # If pattern node expects children nodes, tries to find a
                # combination of target node children that match the pattern

                if len(nodes) == 0:
                    nodes = [node]

                for node in nodes:
                    sub_status_count = 0
                    if len(node.children) >= len(self.children):
                        for candidate in permutations(node.children):
                            sub_status = True

                            for i in range(len(self.children)):
                                st = self.children[i].match(candidate[i], cache)
                                if st and not self.is_in_bounds("low"):
                                    # in case it matches, but has exited the lower bound (in case it exist), force False the match
                                    st = False
                                if st == False and self.controller["allow_indirect_connection"] and len(candidate[i].children) > 0:
                                    pass

                                else:
                                    sub_status &= st
                                    if sub_status == True: sub_status_count += 1

                            if sub_status and sub_status_count > 0:
                                status = True
                                break
                            else:
                                status = False
                    if status and sub_status_count > 0:
                        break

        # 'skipped' tracks the maximum skipped node. So only in case of not match, it decreases
        if not status and self.controller["allow_indirect_connection"]: self.controller["skipped"] -= 1
        return status


    def is_local_match(self, target_node, cache):
        """ Evaluate if this pattern constraint node matches a target tree node.

        TODO: args doc here...
        """

        # Creates a local scope containing function names, variables and other
        # stuff referred within the pattern expressions. We use Syntax() as a
        # container of those custom functions and shortcuts.

        constraint_scope = {attr_name: getattr(self.syntax, attr_name)
                            for attr_name in dir(self.syntax)}
        constraint_scope.update({"__target_node": target_node})
        constraints = []

        if "root" in self.controller and self.controller["root"]: constraints.append("__target_node.is_root()")
        if "leaf" in self.controller and self.controller["leaf"]: constraints.append("__target_node.is_leaf()")


        if not self.name:
            # empty pattern node should match any target node
            constraints.append('')
        elif '@' in self.name:
            # use as any node if used alone or
            # converts references to node itself if it's in an expression.
            if len(self.name) == 1: constraints.append('')
            else: constraints.append(self.name.replace('@', '__target_node'))
        elif self.controller["allow_indirect_connection"]:
            # pattern nodes that allow indirect connection should match any target node
            constraints.append('')

        else:
            # if no references to itself, let's assume we search an exact name
            # match (allows using regular newick string as patterns)
            constraints.append('__target_node.name == "%s"' % self.name)

        try:
            st = True
            for constraint in constraints:
                if constraint:
                    st &= eval(constraint, constraint_scope)
                else: st &= True

        except ValueError:
            raise ValueError("not a boolean result: . Check quoted_node_names.")

        except (AttributeError, IndexError) as err:
            raise ValueError('Constraint evaluation failed at %s: %s' %
                             (target_node, err))
        except NameError:
            try:
                # temporary fix. Can not access custom syntax on all nodes. Get it from the root node.
                root_syntax = self.get_tree_root().syntax
                constraint_scope = {attr_name: getattr(root_syntax, attr_name)
                                    for attr_name in dir(root_syntax)}
                constraint_scope.update({"__target_node": target_node})

                st = True
                for constraint in constraints:
                    if constraint:
                        st &= eval(constraint, constraint_scope)
                    else: st &= True
                return st
            except NameError as err:
                raise NameError('Constraint evaluation failed at %s: %s' %
                         (target_node, err))
        else:
            return st


    def find_match(self, tree, maxhits=1, cache=None, target_traversal="preorder"):
        """ A pattern search continues until the number of specified matches are
        found.
        :param tree: tree to be searched for a matching pattern.
        :param maxhits: Number of matches to be searched for.
        :param cache: When a cache is provided, preloads target tree information
                               (recommended only for large trees)
        :param None maxhits: Pattern search will continue until all matches are found.
        :param target_traversal: choose the traversal order (e.g. preorder, levelorder, etc)
        :param nested: for each match returned,

        """

        # the controller is for only one use, because it changes the node names.
        # for that reason a deepcopy of the pattern and not the pattern
        # is traversed. Very slow operation.
        one_use_pattern = copy.deepcopy(self)

        # indicaes whether the pattern is a pattern of single match
        # such as maximum or minimum value of an attribute
        single_match_pattern = False

        # sets the controller for every node.
        # should chagne if only is a metacharacter acording to benchmarking tests.
        # if one node indicates single_match, the whole pattern is single_match.
        # no functionality to support other than pattern's root single_match yet
        for node in one_use_pattern.traverse():
            single_match_pattern |= node.set_controller()

        # in case of single_match pattern save the nodes and filter them
        # in other case find the requested match
        if single_match_pattern:
            matched_nodes = []
            for node in tree.traverse(target_traversal):
                if one_use_pattern.match(node, cache):
                    matched_nodes += [node]
            yield one_use_pattern.find_extreme_case(matched_nodes)
        else:
            num_hits = 0
            for node in tree.traverse(target_traversal):
                if one_use_pattern.match(node, cache):
                    num_hits += 1
                    yield node
                if maxhits is not None and num_hits >= maxhits:
                    break


def test():
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
    t11 = PhyloTree("""  ( ((e, f, g) c) b, ((g, (w)h, i)c) d) a ; """, format=8, quoted_node_names=False)

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
        print(pattern_match == true_match[p_num])


if __name__ == '__main__':
    test()
