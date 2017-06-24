# -*- coding: utf-8 -*-

from __future__ import print_function, division
import timeit

# Explanation:
# Anything in a "setup_..." variable is excluded from the total time
# Anything in a "test_..." variable is timed

# 5 tests are done:
#  * With no warmup nor cache
#  * With pattern warmup but no cache
#  * With pattern warmup and cache
#  * With cache and cache warmup
#  * With pattern warmup, cache and cache warmup

# Warmup here means that the code is executed once to ensure all caches are populated
# i.e. ignore the time it takes to create any cache

# If the time it takes to create a cache is significant the "warmed up" variants
# should be slightly faster. In the current example the difference is not
# measurable which is puzzling.

# It should take about a couple of seconds to run this.
# If you want to have more precise numbers at the end, increase the number of repeats
repeat = 1000


tree = " '((((Human_1, Chimp_1), (Human_2, (Chimp_2, Yeast_3))), ((FuguRubripes_1, (Pig, Fish_3)), Whale)), Frog);' "

setup_base = """
from ete3 import PhyloTree
from treematcher import TreePattern, TreePatternCache

t = PhyloTree(""" + tree + """)

t.set_species_naming_function(lambda node: node.name.split("_")[0])
t.get_descendant_evol_events()

"""

setup_pattern_warmup = """
# This 'warms up' the tree
pattern.match(t)
"""

setup_cache = """
cache = TreePatternCache(t)
"""

setup_cache_warmup = """
# This 'warms up' the tree and cache object
pattern.match(t, cache)
"""

test_nocache = '''
pattern.match(t)
'''

test_cache = '''
pattern.match(t, cache)
'''

# Matching patterns
match_pattern1 = '''
pattern = TreePattern("""('"Chimp" in species(@)', ''); """)
'''

match_pattern2 = '''
pattern = TreePattern("""('"Macaque" not in species(@)', ''); """)
'''

match_pattern3 = '''
pattern = TreePattern("""('"Yeast" in species(@)', ''); """)
'''

match_pattern4 = '''
pattern = TreePattern("""('"Frog" in species(@)', ''); """)
'''

match_pattern5 = '''
pattern = TreePattern("""('"Human" in species(@)', ''); """)
'''

match_pattern6 = '''
pattern = TreePattern("""('"Fugu" in species(@)', ''); """)
'''

match_pattern7 = '''
pattern = TreePattern("""('"Pig" in species(@)', ''); """)
'''

# Non-matching patterns
nomatch_pattern1 = '''
pattern = TreePattern("""('"James" in species(@)', ''); """)
'''

nomatch_pattern2 = '''
pattern = TreePattern("""('"Worm" in species(@)', ''); """)
'''

nomatch_pattern3 = '''
pattern = TreePattern("""('"Human" not in species(@)', ''); """)
'''

patterns = (
    # Pattern, Expected to match
    (match_pattern1, True),
    (match_pattern2, True),
    (match_pattern3, True),
    (match_pattern4, True),
    (match_pattern5, True),
    (match_pattern6, True),
    (match_pattern7, True),
    (nomatch_pattern1, False),
    (nomatch_pattern2, False),
    (nomatch_pattern3, False),
)

# This just makes sure that the pattern actually produces the expected result
# Prints a warning otherwise
pattern_test = """
if pattern.match(t) != {0}:
    if {0} == True:
        print '''\033[91m>>>>>>> WARNING Pattern does NOT match and should: {1}tree: {2}\033[0m'''
    elif {0} == False:
        print '''\033[91m>>>>>>> WARNING Pattern should NOT match but it does: {1}tree: {2}\033[0m'''
else:
    print '''Using pattern: {1}on tree: {2}'''
"""

for pattern, expected in patterns:
    setup = setup_base + pattern

    # Run the sanity check here
    timeit.timeit(pattern_test.format(expected, pattern, tree), setup=setup, number=1)

    setup_nocache_nowarm = setup
    setup_nocache_warm = setup + setup_pattern_warmup
    setup_cache_nopatwarm = setup + setup_cache
    setup_cache_nocachewarm = setup + setup_pattern_warmup + setup_cache
    setup_cache_warm = setup + setup_pattern_warmup + setup_cache + setup_cache_warmup

    time_nocache_nowarm = timeit.timeit(test_nocache, setup=setup_nocache_nowarm,
                                        number=repeat) / repeat
    time_nocache_warm = timeit.timeit(test_nocache, setup=setup_nocache_warm,
                                      number=repeat) / repeat
    time_cache_nopatwarm = timeit.timeit(test_cache, setup=setup_cache_nopatwarm,
                                         number=repeat) / repeat
    time_cache_nocachewarm = timeit.timeit(test_cache, setup=setup_cache_nocachewarm,
                                           number=repeat) / repeat
    time_cache_warm = timeit.timeit(test_cache, setup=setup_cache_warm,
                                    number=repeat) / repeat

    print("Took {:.3f}".format(time_nocache_nowarm * 1e6), "µs per pattern to search",
          repeat, "patterns without cache nor pattern warmup")
    print("Took {:.3f}".format(time_nocache_warm * 1e6), "µs per pattern to search",
          repeat, "patterns with pattern warmup but no cache")
    print("Took {:.3f}".format(time_cache_nopatwarm * 1e6), "µs per pattern to search",
          repeat, "patterns with cache but no pattern nor cache warmup")
    print("Took {:.3f}".format(time_cache_nocachewarm * 1e6), "µs per pattern to search",
          repeat, "patterns with cache and pattern warmup but no cache warmup")
    print("Took {:.3f}".format(time_cache_warm * 1e6), "µs per pattern to search",
          repeat, "patterns with cache and full warmup")

# vim: ai sts=4 et sw=4
