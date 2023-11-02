import sys
sys.setrecursionlimit(10000)
import ete3
import pdb 
pdb.set_trace()
tree = ete3.Tree("closest_with_outliers.nwk", format = 1)
isols = []
for n in tree:
    if n.name != "AZPAE14941_AE004091":
        isols.append(n.name)
    else:
        tree.set_outgroup(n)
tree.prune(isols)
tree.write(outfile = "closest_with_outliers_rooted.nwk", format = 1)

