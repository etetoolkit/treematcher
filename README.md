# TreeMatcher: A new tool for creating Python-based queries on trees

TreeMatcher is a subproject of ete3.
TreeMatcher will be merged with ete3. Until then treematcher have to be in PYTHONPATH environmental variable.

Insert that line to `~/.bashrc` file so it can be loaded automatically on every login
in case you use bash shell.

```
path_to_treematcher='path_to_treematcher_repository'
export PYTHONPATH=$PYTHONPATH":${path_to_treematcher}"
```

### Program Description

In mathematics, a standard way of representing graphical trees with edge lengths is the Newick format. The TreeMatcher module extends the Newick format to define a tree pattern and includes rules and filters with a Python-based vocabulary. These patterns are then searched for using a tree traversal algorithm. A pattern can be written by accessing the attributes and functions available to an ETE Tree (see Tree and PhyloTree classes), using Python code directly, or through custom functions and syntax.

### How to use treematcher

The simplest way to begin using treematcher is to create a pattern on a single node. In the following example, a string defines the pattern and a TreePattern instance is created. If an attribute is not specified, the node name is assumed by default.

```
# Example 1: Find a node named "sample_1"
pattern1 = ' sample_1 ; '	 # begin with a string
pattern1 = TreePattern(pattern1)  # create a TreePattern Instance

```

Now that you know how to search for the name of a single node, you may be tempted to access other nodes through constraints like @.children[0].name=="sample_1" and @.children[1].name=="sample_2" but calling a node's descendants in this way restricts the order in which they are considered a match. For example, the permutation @.children[0].name=="sample_2" and @.children[1].name=="sample_1" would not be returned as a match. Using the Newick format ensures that both permutations of children are matched.


Note that the format type is set to 1 as the default which does not allow internal node names. Access other Newick format types using the format argument.

```
# Example 2: Find a tree where sample_1 and sample_2 are siblings under ancestor_a
tree = Tree("((sample_1,sample_2)ancestor_a,(sample_1,sample_2)ancestor_b)root;", format = 8)
pattern2 = TreePattern(' (sample_1, sample_2)ancestor_a ; ', format=8)
```

### Quoted node names and the node symbol @
[More about quoted_node_names](sdoc/tutorial/tutorial_treematcher.md#quoted-node-names-and-the-node-symbol-)

In order to differentiate the parentheses of a function call from the parentheses defining Newick structure, quoted node names are used.
That means simply that you enclose each node name in quotes.
In order to access a method on a node, use the @ symbol to represent the node.

Be sure that these quotes are different from those of the overall pattern definition.
For example:
`TreePattern(""" ('the_quoted_pattern'); """, quoted_node_names=True ) `
and not
`TreePattern(" ("the_quoted_pattern") ;", quoted_node_names=True)`


```

# Example 3: Find a tree where sample_2 and another leaf are siblings where a leaf is determined by number of children.
pattern2 = TreePattern(""" ('len(@.children)==0', 'sample_2')ancestor_a ; """, quoted_node_names=True)

# Example 4: Find a tree where sample_1 and another leaf are siblings by accessing the the is_leaf method.
pattern3 = TreePattern(""" ('sample_1', '@.is_leaf()')ancestor_a ; """, quoted_node_names=True)
```

### Relax matches
[More about relax matches](sdoc/tutorial/tutorial_treematcher.md#relax_patterns)

Treematcher allows to test against relax matched patterns.

By setting ` ^ ` as ancestor you enable the loose connection ability.
Loose connection means that the ` ^ ` children may connect
loosely via any number of intermediate nodes.

```
# example 5: test if tips A, B, C exists in the same tree
pattern5 = TreePattern(" (A, B, C)^ ;")

#example 6: test if (A, B) and (C, D) are conected via any number of nodes
pattern6 = TreePattern(" ((A,B)^), ((C,D)^) ;")

```

You can test for a relax number of children too.
You can use:
` * `, ` + ` and ` {min, max} `
They borrow their meaning from regular expressions.

```
#example 7: zero or more ocuurances of a node
TreePattern( " ('A', '@.dist > 0.5*')'^' ;", quoted_node_names=True)

#example 8: exact number
TreePattern(" ('A', '@.dist > 0,5{2, 5}') ;", quoted_node_names=True)

```


### To Run
To run, use the find_match() function. By default, find_match will look for one match.
To find the number of matches returned, use len().

```

# Example 9: Find the parent node of the siblings sample_1 and sample_2
tree = Tree("((sample_1,sample_2)ancestor_a,(sample_1,sample_2)ancestor_b)root;", format = 8)
pattern = TreePattern(' (sample_1) ; ', quoted_node_names=False)
solution = list(pattern.find_match(tree))
print("The number of solutions are: ", len(solution))

```

For more details on how to use treematcher read the tutorial.
[Treematcher tutorial](sdoc/tutorial/tutorial_treematcher.md)

### Advanced Topics

To make treematcher perform faster, break complex patterns into smaller searches. If conditional statements are used, try putting the part of the search that you think will be faster first.

####  Custom Functions
You can use your own custom functions and syntax in treematcher.  In the following example, a custom function is created in a custom class called MySyntax.

```
# Example 8: Expanding vocabulary
class MySyntax(PatternSyntax):
	def my_nice_function(self, node):
		return node.species == 'Chimp'

my_syntax = MySyntax()

pattern = """ 'my_nice_function(@)'; """
t_pattern = TreePattern(pattern, syntax=my_syntax)
for match in t_pattern.find_match(t):
	print(list(match))

```

### Command line tool

ete_search is the command line interface to treematcher. Using ete_search you can run multiple
pattern comparisons to multiple trees using files or text and retrieve some basic statistics and
save results for later use.

[Tutorial on command line tool](sdoc/tutorial/tutorial_treematcher.rst#ete_search-command-line-tool)

examples:
Read patterns from a file called MyPatterns.txt and apply to each tree in MyTargetTrees.txt, output the results of each pattern in separate files called treematches0.txt, treematches1.txt, etc
If there is only one pattern, the result file will not be numbered.

`python -m treematcher.tools.ete_search --pattern_tree_list "MyPattern.txt" --tree_format 8 --src_tree_list "MyTargetTrees.txt" -o treematches.txt `

Provide the pattern and tree as strings and print the result to the terminal.
`python -m treematcher.tools.ete_search -p "(e,d);" --tree_format 8 -t "(c,(d,e)b)a;" `


Count how many trees matches a pattern from a list of trees.
` python -m treematcher.tools.ete_search -p "(the, pattern)" --src_tree_list trees.file --root | wc -l`


The render option will save each match as an image. If there are multiple patterns, numbers will be used to designate each pattern starting from 0.
If there are multiple matches, and underscore is used with a number for each match starting with 0. If I had two

`python -m treematcher.tools.ete_search --pattern_tree_list "MyPatterns.txt" --tree_format 8 --src_tree_list "MyTargetTrees.txt" --render treematches.png `
