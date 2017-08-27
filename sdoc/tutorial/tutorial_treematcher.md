# TreeMatcher: A new tool for creating Python-based queries on trees

### Program Description

In mathematics, a standard way of representing graphical trees with edge lengths is the Newick format. The TreeMatcher module extends the Newick format to define a tree pattern and includes rules and filters with a Python-based vocabulary.
These patterns are then searched for. A pattern can be written by accessing the attributes and functions available to an ETE Tree (see Tree and PhyloTree classes), using Python code directly, or through custom functions and syntax.
Treematcher's result is a generator object containing sub-trees with nodes that matched the pattern.


### How to use treematcher


The simplest way to begin using treematcher is to create a pattern on a single node.
In the following example, a string defines the pattern and a TreePattern instance is created.

If an attribute is not specified, the node name is assumed by default.

` Example with a simple match `

```
# Example 1: Find a node named "sample_1"

# create a demo tree.
tree = Tree(" ((sample_1, sample_2), sample_3); ")

# create the patternt to search for.
pattern1 = ' sample_1 ; '	 # begin with a string
pattern1 = TreePattern(pattern1)  # create a TreePattern Instance

# search for the pattrn
result = pattern1.find_match(tree) # returns a generator object

```


#### Types of nodes

There are two types of nodes that can be written in a TreePattern.

- Node name
- Node reference (access attributes)

Treematcher uses ` @ ` to access a node reference and reach the node attributes.
Any other string is considered a target node name (meaning name as attribute).

```
# Example 2: Find a tree where sample_1 and sample_2 are siblings.

# create a demo tree.
tree2 = Tree(" ((sample_1, sample_2), sample_3); ")

pattern2 = TreePattern(' (sample_1, sample_2) ; ')

# search for the pattern
result = pattern2.find_match(tree2) # returns a generator object
```

###### Quoted node names and the node symbol @
In order to differentiate the parentheses of a function call from the parentheses defining Newick structure, quoted node names are used. The quotes surrounding each node will be removed and the contents inside will be processed as python code. If no quotes are present, as in the previous examples, all parenthesis are assumed to be part of the Newick structure. In order to ensure that the pattern is being processed correctly, set the quoted_node_names to False when not quoting node names. You can set quoted_node_names to True when they are used but this is assumed by default. In order to access a method on a node, use the @ symbol to represent the node.

Be sure that these quotes are different from those of the overall pattern definition.
For example:
`TreePattern(""" ('the_quoted_pattern'); """, quoted_node_names=True ) `
and not
`TreePattern(" ("the_quoted_pattern") ;", quoted_node_names=True)`


```
# Example 3: search a node using an attrubute
tree3 = Tree(" (((D, E), (B, C)), A);")
(tree3&'E').dist = 0.3 # Tree()&node_name returns the node with that name
pattern3 = TreePattern("('@.dist == 0.3', 'D');",  quoted_node_names=True)
result = pattern3.find_match(tree3)
```

#### Types of patterns

- Strict patterns
- Relax patters

Both types of nodes (Node names and Node reference) can be used at both types of
patterns and can be mixed in a pattern.

The actual node name is a truth condition. So the target is to have a valid truth
condition for each node. When a node name is written it's being translated
to test the @.name attribute.

##### Strict patterns

A strict pattern is written normally as a (sub)tree and will match sub-trees with
the same topology and order as itself.

```
#Example 4: Search for a strict match
tree4 = Tree(" (((F, G)E, (C, D)B), A);", format=8)
pattern4 = TreePattern(" (D,C)B ;", format=8)
result = pattern4.find_match(tree4)
```

```
# Example 5: Search for strict match, but using  node reference
tree5 = Tree(" (((F, G)E, (C, D)B), A);", format=8)
pattern5 = TreePattern(""" ('C', '@.support > 0+')'B' ;""", format=8, quoted_node_names=True)
result = pattern5.find_match(tree5)
```

##### Relax patterns

A relax pattern is written as a newick tree, but in addition it has some symbols
mimicking the behavior of regular expressions' metacharacters  (characters with
special meaning).

These metacharacters have impact on the tree topology on how nodes are connected
or the number and/or the attributes of the children a node may have.

##### Relax topology

###### Connect via any number of nodes.

Relax topology uses ` ^ ` and it's written as an ancestor node.
That symbol indicates that each of it's child/children may connect to the rest of pattern
via any number of intermediate nodes. That symbol overrides the pattern's topology
meaning that it's child/children may exists at a higher or lower level.

```
# Example 6: Search for relax match between two nodes
  tree6 = Tree("((A,B), (C,D),((E,F), (G,H)));")
  pattern6 = TreePattern("(G,F)^;")

  result = pattern6.find_match(tree6)
```

```
# Example 7: Search for relax match between subtress
  tree7 = Tree("((A,B), (C,D),((E,F), (G,H)));")
  pattern7 = TreePattern("((C,D)^), (E,F)^;")

  result = pattern7.find_match(tree7)
```

###### Relax number of children

Relax number of children has symbols to indicate continuous repeat of a node
and uses:
- ` * ` zero or more
- ` + ` one or more
- ` {min, max} ` defined number

```
# Examples 8, 9, 10 : showing the use of metacharacters.

tree8 = Tree(" (a, a, b), ((c, d), (e, f, g, g, g)) ;")
pattern8 = TreePattern("(b, a+)^;")
pattern9 = TreePattern(" (a+, b*)^;")
pattern10 = TreePattern( """ ('e', 'f', 'g{1,3}'); """, quoted_node_names=True)


result = pattern8.find_match(tree8)
result = pattern9.find_match(tree8)
result = pattern10.find_match(tree8)
```
These symbols borrows it's behavior from regular expressions' meaning.


A short list of commonly used constraints is given in the following table.

Table 1: Examples of common constraints.

|  type                     |custom|  syntax example       						            | example meaning       				        |  Comments																        |
| --------------------------|:-:|:---------------------------------------------------------:|:---------------------------------------------:|:-----------------------------------------------------------------------------:|
| node                      |   | @	            						                    |a  node, default for nodes left blank	        | Use @.attribute to access attribute, function(@) to access function           |
| node name                 |   | node_name, "node_name", or @.name=="node_name"	        | when attribute not specified, name is assumed | Looking for multiple names, use list: @.name in ("sample1","sample2")         |
| distance                  |   | @.dist >= 0.5     					                    | branch length no less than 0.5		        | Use any of the following: <, <=, ==, >=, !=								    |
| support                   |   | @.support > 0.9	            		                    | Has a support value greater than 0.90	        | 																		        |
| species                   |   | @.species=="Homo sapiens"	    		                    | Homo sapiens is species of node		        | See set_species_naming_function()	for details							        |
| scientific name           |   | @.sci_name == Euarchontoglires 		                    | scientific name is Euarchontoglires	        | See annotate_ncbi_taxa() function for details 						        |
| rank                      |   | @.rank == subfamily 					                    | node is ranked at the subfamily level	        | See annotate_ncbi_taxa() function for details							        |
| taxonomic id              |   | @.taxid == 207598						                    | 20758 is the taxid of the node			    | See annotate_ncbi_taxa() function for details							        |
| number of children        |   | len(@.children)						                    | binary tree internal node has 2, leaf hss 0   | use quoted node names to differentiate parentheses from Newick Structure	    |
| size of subtree           |   | @.get_descendants() > 5							        | size of tree  is greater than 5               | number of descendants												            |
| number of leaves          |   | len(@) > 2                                                | number of leaves is greater than 2            | number of leaves descending from a node                                       |
| lineage                   | * | 9606 in @.lineage or "Homo sapiens" in @.named_lineage    | Homo sapiens in @.lineage				        | Find NCBI taxonomy ID or the full scientific name	in a node's lineage	        |
| species in descendant node| * | contains_species(@, ["Pan troglodytes", "Homo sapiens"])	| Find species in the last at or below node     | species at a node and any of it's descendants							        |
| leaf name                 | * | contains_leaves(@, ["Chimp_2", "Chimp_3"])		        | Pan_troglodytes_1 is descendant leaf name	    | Find the leaf name within a list of leaf names                                |
| number of duplications    | * |  		n_duplications(@) > 0                               | Number of duplications beyond and including this node is greater than zero.	    | number of duplication events at or below a node  |
| number of speciations     | * |  		n_speciations(@) > 0                                | Number of speciations beyond and including this node is greater than zero.	    | number of speciation events at or below a node  |

* functions do not exist outside of treematcher classes.



### Advanced Topics

#### optimization

Virtually any attribute available in ETE can be searched for on a tree, however, the larger the structure the more complex the pattern is, the more computationally intensive the search will be. Large Newick trees with complex conditional statements calling functions that require several tree traversals is not recommended.
Instead, break complex patterns into smaller searches. If conditional statements are used, try putting the part of the search that you think will be faster first.



## ete_search command line tool.

|  argument       						| meaning       						                                                  |
| --------------------------------------|:---------------------------------------------------------------------------------------:|
| -p								    | a list of patterns in newick format (filenames with one per file or quoted strings)     |
| -t								    | a list of trees in newick format (filenames or quoted strings)                          |
|-v                                     | prints the current pattern, prints which trees (by number) do not match the pattern     |
| --tree_format							| format for trees, default = 1	                            		                      |
| --quoted_node_names 					| default = True					                            	                      |
| -o, --output                  | output file for search results
| --src_tree_list                       | path to a file containing many target trees, one per line                               |
| --pattern_tree_list                   | path to a file containing many pattern trees, one per line                              |
|-r, --root                             | flag to return the root of the tree if at least a match was found
| --render                              | filename (.SVG, .PDF, or .PNG), to render the tree image                                |
| --tab                                 | output results in tab delimited format, default if -o used and ascii not specified      |
| --ascii                               | output results in ascii format                                                          |



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
