import pdb
from ete3 import Tree
from treematcher import TreePattern

t = Tree('((aaaaaaaaad:1,(aaaaaaaaae:1,aaaaaaaaaf:1)1:1)1:1,((aaaaaaaaag:1,aaaaaaaaah:1)1:1,((aaaaaaaaai:1,aaaaaaaaaj:1)1:1,(aaaaaaaaaa:1,(aaaaaaaaab:1,aaaaaaaaac:1)1:1)1:1)1:1)1:1);')

#print t
"""
      /-aaaaaaaaad
   /-|
  |  |   /-aaaaaaaaae
  |   \-|
  |      \-aaaaaaaaaf
--|
  |      /-aaaaaaaaag
  |   /-|
  |  |   \-aaaaaaaaah
  |  |
   \-|      /-aaaaaaaaai
     |   /-|
     |  |   \-aaaaaaaaaj
      \-|
        |   /-aaaaaaaaaa
         \-|
           |   /-aaaaaaaaab
            \-|
               \-aaaaaaaaac
"""
#pt = '((aaaaaaaaaa)*,(aaaaaaaaaj)*)@;'
#tpt = TreePattern(pt)
#print "As expected: " + str(len(list(tpt.find_match(t))) > 0)

"""
      /-aaaaaaaaad
   /-|
  |  |   /-aaaaaaaaae
  |   \-|
  |      \-aaaaaaaaaf
--|
  |      /-aaaaaaaaag
  |   /-|
  |  |   \-aaaaaaaaah
  |  |
   \-|      /-aaaaaaaaai
     |   /-|
     |  |   \-aaaaaaaaaj
      \-|
        |   /-aaaaaaaaaa
         \-|
           |   /-aaaaaaaaab
            \-|
               \-aaaaaaaaac
"""

#pdb.set_trace()
pt = '((aaaaNOTaaaaaa)+,(aaaaaaNOTaaaj)+)@;'
tpt = TreePattern(pt)
test = tpt.find_match(t)
print "As expected: " + str(len(list(test)) == 0)



"""
      /-aaaaaaaaad
   /-|
  |  |   /-aaaaaaaaae
  |   \-|
  |      \-aaaaaaaaaf
--|
  |      /-aaaaaaaaag
  |   /-|
  |  |   \-aaaaaaaaah
  |  |
   \-|      /-aaaaaaaaai
     |   /-|
     |  |   \-aaaaaaaaaj
      \-|
        |   /-aaaaaaaaaa
         \-|
           |   /-aaaaaaaaab
            \-|
               \-aaaaaaaaac
"""


#for node in t.traverse():
#	print node.name