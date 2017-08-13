#!/usr/bin/env python

import sys
import logging
import os.path
from argparse import ArgumentParser
from ete3.tools.common import src_tree_iterator
from ete3.phylo import PhyloTree
from treematcher import TreePattern, TreePatternCache

class match_stats(object):
    def __init__(self, name=""):
        self.name = name
        self.total = 0
        self.num_of_patterns = 0 # > 0 is for Summarized results
        self.num_of_trees = 0    # > 0 is for Summarized results
        self.matched = 0
        self.not_matched = 0
        self.errors = 0

    def __str__(self):
        printable = "{}\n".format(self.name)
        printable += "Total comparisons: {}\n".format(self.total)
        printable +="Matches found: {}\n".format(self.matched)
        printable +="Not matches: {}\n".format(self.not_matched)
        if self.num_of_patterns > 0:
            printable += "Number of patterns: {}\n".format(self.num_of_patterns)
        if self.num_of_trees > 0:
            printable += "Number of trees: {}\n".format(self.num_of_trees)

        printable +="Errors: {}\n".format(self.errors)
        return printable

DESC=''

#ete3 treematcher --pattern "(hello, kk);" --pattern-format 8 --tree-format 8 --trees "(hello,(1,2,3)kk);" --quoted-node-names

def populate_args(treematcher_args_p):

    treematcher_args = treematcher_args_p.add_argument_group('TREEMATCHER OPTIONS')


    treematcher_args.add_argument("--quoted_node_names", dest="quoted_node_names",
                              action="store_true",
                              help="True if using quotes to designate node names in the pattern. Otherwise False.")
    treematcher_args.add_argument("--tree_format", dest="tree_format",
                              type=int,
                              default=0,
                              help="A number 0-8 designating Newick format.")
    treematcher_args.add_argument("--tab", dest="taboutput",
                              action="store_true",
                              help="output results in tab delimited format")
    treematcher_args.add_argument("--ascii", dest="asciioutput",
                              action="store_true",
                              help="output results in ascii format")
    treematcher_args.add_argument("-t", "--tree", dest="src_trees", type=str,
                                nargs="*", help=("a list of trees in newick format (filenames or"
                                "quoted strings) to be used as target tree(s)"))
    treematcher_args.add_argument("--target_tree_list", dest="src_tree_list",
                              type=str,
                              help=("path to a file containing many pattern trees, one per line"))
    treematcher_args.add_argument("-p", dest='pattern_trees',
                              type=str, nargs="*",
                              help=("a list of trees in newick format (filenames or"
                              "quoted strings) to be used as pattern tree(s)"))
    treematcher_args.add_argument("--pattern_tree_list", dest="pattern_tree_list",
                              type=str,
                              help=("path to a file containing many pattern trees, one per line"))
    treematcher_args.add_argument("-o", "--output", dest="output", type=str,
                                help=("specify an output file"))
    treematcher_args.add_argument("--render", dest="render",
                               type=str,
                               help="filename (.SVG, .PDF, or .PNG), to render the tree")
    treematcher_args.add_argument("-r", "--root", dest="whole_tree", action="store_true",
                                    help=("Returns the tree from root if match found. Is used as\
                                    flag to indicate match presence rather than match it self."))
    treematcher_args.add_argument("-v", "--verbosity", dest="verbosity",
                                    type=int, nargs=1,
                                    help=("A number between 1-4. The verbosity level.\
                                    1: print matches (default), 2: print statistsics, \
                                    3: print the pattern, 4: print statistsics for \
                                    each pattern."))

def run(args):
    # a list of stats objects. one for every pattern
    all_stats = []

    if vars(args)["src_trees"] is None and vars(args)["src_tree_list"] is None:
        logging.error('Please specify a tree to search (i.e. -t) ')
        sys.exit(-1)
    if not vars(args)["pattern_trees"] and not vars(args)["pattern_tree_list"]:
        logging.error('Please specify a pattern to search for. (i.e. -p)')
        sys.exit(-1)


    pattern_length = len(list(pattern_tree_iterator(args)))

    for pattern_num, p in enumerate(pattern_tree_iterator(args)):
        try :
            pattern = TreePattern(p, quoted_node_names=vars(args)["quoted_node_names"])
        except:
            logging.error("Could not create pattern from newick.")
            continue

        stats = match_stats("pattern_" + str(pattern_num))

        # handle file creation
        if vars(args)["output"]:
            filename = vars(args)["output"]
            if pattern_length > 1:
                if '.' in vars(args)["output"]:
                    filename = filename.replace('.', str(pattern_num) + '.')
                else:
                    filename += str(pattern_num)

            outputfile = open(filename, 'w')

        if vars(args)["verbosity"] and int(vars(args)["verbosity"][0]) > 2:
            print("pattern_{} is: ".format(pattern_num))
            print(pattern)

        # for every tree
        if vars(args)["verbosity"] and vars(args)["verbosity"][0] > 2 and not vars(args)["output"]:
            print("match(es) for pattern_{}:".format(pattern_num))

        for n, nw in enumerate(src_tree_iterator(args)):
            stats.total += 1
            try:
                t = PhyloTree(nw, format=args.tree_format)
            except:
                logging.error("Could not creat tree from newick format.")
                stats.errors += 1
                continue

            matches = list(pattern.find_match(t))
            match_length=len(matches)
            if match_length > 0:
                stats.matched += 1
            else:
                stats.not_matched += 1

            if args.render:
                image = args.render
                if pattern_length > 1:  # multiple patterns
                    if match_length > 1:  # one file per match on each pattern
                        for m, match in enumerate(matches):
                            if '.' in image:
                                image = image.replace('.', str(pattern_num) + '_' + str(m) + '.')
                            else:
                                image += str(pattern_num) + str(m)
                            match.render(image)
                    elif match_length == 1:  # One match on multiple patterns
                        if '.' in image:
                            image = image.replace('.', str(pattern_num) + '.')
                        else:
                            image += str(pattern_num)
                        matches[0].render(image)
                    else:
                        if vars(args)["verbosity"] and vars(args)["verbosity"][0] > 1:
                            print("No matches for pattern {} tree {}".format(pattern_num, n))
                else:  # one pattern
                    if match_length > 1:  # one file per match on one pattern
                        for m, match in enumerate(matches):
                            if '.' in image:
                                image = image.replace('.', '_' + str(m) + '.')
                            else:
                                image += str(m)
                            match.render(image)
                    elif match_length == 1:  # one file for one match
                        matches[0].render(image)
                    else:
                        if vars(args)["verbosity"][0] > 1:
                            print("No matches for tree {}".format(n))

            if vars(args)["output"]:
                if vars(args)["asciioutput"]:
                    if vars(args)["whole_tree"] and match_length > 0:
                        outputfile.write(str(t))
                    else:
                        for match in matches:
                            outputfile.write(str(match) + '\n')
                else:  #args.taboutput
                    if vars(args)["whole_tree"]:
                        outputfile.whrite(t.write(features=[]))
                    else:
                        outputfile.write('\t'.join([match.write(features=[]) for match in matches]))

            if not vars(args)["output"] and not args.render:
                if vars(args)["asciioutput"]:
                    if vars(args)["whole_tree"] and match_length > 0:
                        print(t)
                    else:
                        for match in matches:
                            print(match)
                else:
                    if vars(args)["whole_tree"] and match_length > 0:
                        print(t.write(features=[]))
                    else:
                        for match in matches:
                            print(match.write(features=[]))

        all_stats += [stats]
        if vars(args)["verbosity"] and vars(args)["verbosity"][0] > 3:
            print("{}".format(stats))

        if vars(args)["output"]:
            outputfile.close()

    concentrated = match_stats("\nSummarize")
    concentrated.total = sum([ stat.total for stat in all_stats])
    concentrated.num_of_patterns = len(all_stats)
    concentrated.num_of_trees = concentrated.total / concentrated.num_of_patterns
    concentrated.matched = sum([stat.matched for stat in all_stats])
    concentrated.not_matched = sum([stat.not_matched for stat in all_stats])
    concentrated.errors = sum([stat.errors for stat in all_stats])

    if vars(args)["verbosity"] and vars(args)["verbosity"][0] > 1:
        print("{}".format(concentrated))

def pattern_tree_iterator(args):
    if not vars(args)["pattern_trees"] and not sys.stdin.isatty():
        vars(args)["pattern_trees"] = sys.stdin
    if vars(args)["pattern_trees"]:
        for p_tree in vars(args)["pattern_trees"]:
            yield p_tree.strip()
    elif vars(args)["pattern_tree_list"]:
        for line in open(vars(args)["pattern_tree_list"]):
            line = line.strip()
            if line:
                yield line




if __name__ == "__main__":
    parser = ArgumentParser()
    populate_args(parser)
    args = parser.parse_args(sys.argv[1:])
    run(args)







    '''
    #add in ete.py

     from . import (ete_split, ete_expand, ete_annotate, ete_ncbiquery, ete_view,
               ete_generate, ete_mod, ete_extract, ete_compare, ete_evol,
               ete_maptrees, ete_treematcher)


    # -treematcher-
    treematcher_args_p = subparser.add_parser("treematcher", parents=[source_args_p, main_args_p],
                                           description=ete_treematcher.DESC)
    treematcher_args_p.set_defaults(func=ete_treematcher.run)
    ete_treematcher.populate_args(treematcher_args_p)


    '''
