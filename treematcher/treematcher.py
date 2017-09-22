import re
import itertools
from collections import defaultdict, OrderedDict

import six
from copy import deepcopy
from ete3 import PhyloTree, Tree, NCBITaxa

from pprint import pprint

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

class TreePattern(Tree):
    def __str__(self):
        return self.get_ascii(show_internal=True, attributes=["name"])

    def __repr__(self):
        return "Pattern node '%s' (%s)" %(self.name, hex(self.__hash__()))

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


    def parse_metacharacters(self, raw_constraint):
        """Takes a string as node name, extracts metacharacters and interpret them as
        min and max occurrences. Assumes that all metacharacters are defined at
        the end of the string.

        """
        raw_constraint = raw_constraint.strip()
        if raw_constraint.endswith('+'):
            self.min_occur = 1
            self.max_occur = 9999999
            raw_constraint = raw_constraint[:-1]
        elif raw_constraint.endswith('*'):
            self.min_occur = 0
            self.max_occur = 9999999
            raw_constraint = raw_constraint[:-1]
        elif raw_constraint.endswith('}'):
            exp = '\{\s*(\d+)\s*,\s*(\d+)\s*\}'
            s = re.search(exp, raw_constraint)
            if s:
                minv, maxv = map(int, s.groups())
                self.min_occur = minv
                self.max_occur = maxv
                raw_constraint = re.sub(exp, '', raw_constraint)
            else:
                self.min_occur = 1
                self.max_occur = 1

        else:
            self.min_occur = 1
            self.max_occur = 1

        if raw_constraint.startswith('~') and self.children:
            self.loose_children = True
            raw_constraint = raw_constraint[1:]
        else:
            self.loose_children = False

        return raw_constraint.strip()

    def parse_node_name(self):
        """transforms node.name into an evaluable python expression. Metachars are also
        extracted and parsed as min and max occurrence information.
        """
        clean_name = self.parse_metacharacters(self.name)

        # Translate alias and shortcut expressions in clean names
        if '@' not in clean_name:
            constraint = '__target_node.name == "%s"' %clean_name
        elif clean_name:
            constraint = clean_name.replace('@', '__target_node')
        else:
            constraint = 'True'

        if self.children:
            constraint = '(%s) and __target_node.children' %constraint
        else:
            constraint = '(%s) and not __target_node.children' %constraint

        return constraint

    def init_controller(self):
        """
        Creates a dictionary that contains information about a node.
        That information is about how a node interacts with the tree topology.
        It describes how the metacharacter connects with the rest of nodes and
        if it is leaf or root.
        """
        # Interpret node name to python expression
        self.constraint = self.parse_node_name()


    def is_local_match(self, target_node, cache):
        """ Evaluate if a tree nodes matches the constraints in this pattern node.  """

        # Creates a local scope containing function names, variables and other
        # stuff referred within the pattern expressions. We use Syntax() as a
        # container of those custom functions and shortcuts.
        constraint_scope = {attr_name: getattr(self.syntax, attr_name)
                            for attr_name in dir(self.syntax)}
        constraint_scope.update({"__target_node": target_node})

        try:
            if self.constraint:
                st = eval(self.constraint, constraint_scope)
            else:
                st = True

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

    def find_match(self, t):
        return find_matches(t, self)



# NEW APPROACH
def compute_match_matrix(pattern, tree):
    '''Computes a dictionary where keys are all the constraints observed in a
    pattern and values all nodes matching those patterns.'''

    c2nodes = defaultdict(set)
    for n in tree.traverse():
        for cn in pattern.traverse():
            if cn.is_local_match(n, None):
                c2nodes[cn.constraint].add(n)
    return c2nodes

def children_match(tnode, pnode, c2nodes, loose_constraint=None):
    '''returns True if a subtree (tnode) matches recursively a given pattern
    (pnode), handling min and max number of occurrences. pnode should not
    contain loose connections
    '''

    # If no children expected in pattern node, return True, as local
    # conditions have already been checked
    if not pnode.children:
        return True

    t_children = set(tnode.children)

    matches = []
    matched_children = set()
    constraint2max_occur = defaultdict(lambda: [set(), 0, 0])
    for pnode_ch in pnode.children:
        match_nodes = c2nodes[pnode_ch.constraint] & t_children
        constraint2max_occur[pnode_ch.constraint][1] += pnode_ch.min_occur
        constraint2max_occur[pnode_ch.constraint][2] += pnode_ch.max_occur
        constraint2max_occur[pnode_ch.constraint][0].update(match_nodes)

        # check min nodes each pattern constraint
        if pnode_ch.min_occur > 0 and len(match_nodes) < pnode_ch.min_occur:
            #print 1
            return False

        # Keep track all children nodes with matches
        matched_children.update(match_nodes)

        # And prepare all permutations of matches in node with minimum occurrences
        if not match_nodes and pnode_ch.min_occur == 0:
            matches.append([[None]])
        else:
            matches.append(list(itertools.permutations(match_nodes, max(pnode_ch.min_occur, 1))))

    # there should be nodes without a match
    if len(matched_children) < len(t_children):
        #print 2
        return False

    # Let's check if there is a non-overlapping combination of nodes matches
    # satisfying patterns. For instance, avoid cases where one node matches the
    # two required patterns
    for comb in itertools.product(*matches):
        valid = set()
        potential_match = comb
        for x in comb:
            x = set(x)
            inter = valid & set(x)
            if not inter:
                valid.update(x)
            else:
                potential_match = None
                # let's check next comb
                break

        # Validate max number of occurrences assuming current valid combination
        if potential_match:
            for k, (matches, mino, maxo) in constraint2max_occur.items():
                remain_match_nodes = matches - valid
                if len(remain_match_nodes) > (maxo - mino):
                    potential_match = None
                    break

        # Validate descendants
        if potential_match:
            match = True
            # Let's check inside the node
            for i, pnode_ch in enumerate(pnode.children):
                for tnode_ch in potential_match[i]:
                    if tnode_ch is None:
                        continue
                    if not children_match(tnode_ch, pnode_ch, c2nodes):
                        match = False
                        break
                if not match:
                    break
            if match:
                return True

    #print 3
    return False

def split_by_loose_nodes(pattern):
    '''split a pattern tree into all subpatterns connected through loose connections
    (allowing multiple intermediate between them). '''

    pnode2content = pattern.get_cached_content(leaves_only=False)

    # split pattern by loose nodes and store a list of partitions that can be
    # used for strict matches
    to_visit = set()
    for n in list(pattern.traverse()):
        if n.loose_children:
            for ch in n.get_children():
                if not ch.loose_children:
                    to_visit.add(ch.detach())
                else:
                    ch.detach()
        elif n is pattern:
            to_visit.add(pattern)

    # Calculate expected groupings of the split partitions
    expected_groups = set()
    for p, content in six.iteritems(pnode2content):
        c = frozenset(content & to_visit)
        if len(c) > 1:
            expected_groups.add(c)

    return to_visit, sorted(expected_groups, key=lambda x: len(x))


def find_matches(tree, pattern):
    '''Iterate over all possible matches of pattern in tree'''
    pattern = deepcopy(pattern)
    for n in pattern.traverse():
        n.init_controller()

    c2nodes = compute_match_matrix(pattern, tree)
    root2matches = OrderedDict()
    to_visit, expected_groups = split_by_loose_nodes(pattern)

    for proot in to_visit:
        matches = []
        for match_node in c2nodes[proot.constraint]:
            if children_match(match_node, proot, c2nodes):
                matches.append(match_node)
        if not matches:
            raise StopIteration

        root2matches[proot]=matches

    if len(root2matches) == 1:
        for match in root2matches[proot]:
            yield match
        raise StopIteration

    p2index = {p:i for i,p in enumerate(root2matches.keys())}
    to_visit = set(to_visit)
    for nodes in itertools.product(*root2matches.values()):
        ancestors = list()
        if len(nodes) != len(set(nodes)):
            continue
        is_match = True
        for group in expected_groups:
            observed_group = [nodes[p2index[v]] for v in group]
            anc = tree.get_common_ancestor(observed_group)
            if anc not in ancestors:
                ancestors.append(anc)
            else:
                is_match = False
                break
        if is_match:
            yield ancestors[-1]

def expand_loose_connection_aliases(nw):
    def find_first_unmatched_closing_par(string):
        open_par = 0
        for i, ch in enumerate(string):
            if ch == '(':
                open_par += 1
            elif ch == ')':
                if open_par == 0:
                    return i
                open_par -= 1
        return -1
    chunks = []
    i = 0
    while True:
        next_split = nw[i:].find('~')
        if next_split == -1:
            chunks.append(nw[i:])
            break

        chunks.append(nw[i:next_split])
        chunks.append(',')

        i += next_split + 1
        ipoint = find_first_unmatched_closing_par(nw[i:])
        chunks.append(nw[i:i+ipoint+1])
        chunks.append('~')
        i += ipoint + 1

    return ''.join(chunks)



# Facts:
#   1) A pattern is a tree structure expressed in newick format.
#
#   2) Every pattern node is associated to an attributes constraint, which is
#   expressed and pattern node name and can be a python expression.
#
#   3) The pattern topology is also used as a constraint. If no special symbols
#   are used, a strict match of topology (parent-node relationships) is
#   necessary to have a match
#
#   4) The symbol '+' at the end of a pattern node constraint expression
#   indicates that ONE OR MORE nodes like the one specified can exist as SISTER
#   nodes
#
#   5) The syntax '{min,max}' at the end of a pattern node constraint
#   expression indicates that between MIN and MAX number of nodes matching the
#   constraint can exist as SISTER nodes
#
#   6) The symbol '*' at the end of a pattern node constraint expression
#   indicates that the between 0 and INFINITY number of nodes matching the
#   constraint can exist as SISTER nodes
#
#   7) The modifier symbols cab be used both in terminal and internal nodes.
#   When used in internal nodes, they refer to the whole partition defined by
#   the node.
#
#   8) To express that two terminal or internal nodes can be linked by an
#   arbitrary number of connections (so called a LOOSE CONNECTION), the symbol
#   '~' should be added to their parent nodes
#
# How the matching algorithm works:
#
# The problem of finding matches is subdivided into three sub-problems:
#
#  A1) finding matches of the strict pattern expressions (no looses connections contained).
#
#  A1.1) The original pattern is split into pieces (sub-patterns) that contain
#  no loose connections. This is mostly done by the `split_by_loose_nodes()` function
#
#  A1.2) Each sub-pattern is treated individually and used as input for A2.
#
#  A2) find exact matches of a sub-pattern in a tree
#
#   A2.1) Recursive search starting from the root node of the sub-pattern and
#   descending through. Each pattern node is evaluated to check if their
#   children combination exist in the target tree.
#
#   Which TREE nodes match each children PATTERN node.
#   Fails if:
#     if any child pattern nodes cannot find min and max number of nodes in the target target.
#
#
#
# to be continued...
