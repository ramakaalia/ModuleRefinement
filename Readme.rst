**REFMOD**

**About**

This module refines network modules from a protein protein interaction network. Large modules obtained by a modularity optimization algorithm like 
Louvain [1] or Clauset-Newman-Moore Greedy maximization algorithm [2] can be further remodularized and refined using our algorithm. 
The algorithm works by relaxing the maximal modularity values and mining submodules of larger modules in suboptimal zone of modularity.
For details, see our publication "Refining modules to determine functionally significant clusters in molecular networks. Rama Kaalia and Jagath C. Rajapakse. BMC Supplements"

**Usage**

```

#Modularity optimization can be done using two algotihms:
    #method: greedy for Clauset-Newman-Moore greedy modularity maximization
    #method: louvain for Louvian Modularity optimization
#Input: Following files are needed
    #edgelist: Tab spearated file of interactions between proteins or genes
#Output: original and refined node membership files, community size files, modularity values and refined modules for different iterations
import import_ipynb
import refmod_app
import networkx as nx
import timeit
import numpy as np
from networkx.algorithms.community import greedy_modularity_communities


start= timeit.default_timer()
#Read edgelist as a graph using networkx
G = nx.read_edgelist("edgelist.tsv", delimiter='\t', nodetype=str)
#method="greedy"
method="louvain"

iterations=np.array(list(range(1,6)))
resolution=np.arange(1,3,1)
#note: Iterations should be outer loop and resolution should be inner loop
for i in iterations:
    for r in resolution:
        ref_main(G,method,it=i,g=r,writeorig=True)

stop= timeit.default_timer()
print("runtime : "+str(stop-start))

```


**References:**

[1] Blondel, V. D., Guillaume, J.-L., Lambiotte, R. & Lefebvre, E. Fast unfolding of communities in large networks. J. Stat. Mech. theory Exp. 2008, P10008 (2008).
[2] Clauset, A., Newman, M. E. J. & Moore, C. Finding community structure in very large networks. Phys. Rev. E 70, 66111 (2004).
