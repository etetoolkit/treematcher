#!/usr/bin/env python

from ete3 import PhyloTree, Tree, NCBITaxa
from treematcher import TreePattern


a = Tree("""  (((a,b)c)d, ((g,h)k)t)l ;""", format=1)


print a

b = TreePattern( """  ((a,b)+)l ;""")

print b

lala  = b.find_match(a)

for i in lala:
		print i