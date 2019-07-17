**REFMOD: Module Refinement Algorithm**

**About:**

    This module refines network modules from a protein protein interaction network. Large modules obtained by a modularity optimization algorithm like Louvain [1] or Clauset-Newman-Moore Greedy maximization algorithm [2] can be further remodularized and refined using our algorithm. 
    The algorithm works by relaxing the maximal modularity values and mining submodules of larger modules in suboptimal zone of modularity. The observed modules are shown to be biologically more relevant.
    For details, see our publication *"Refining modules to determine functionally significant clusters in molecular networks. Rama Kaalia and Jagath C. Rajapakse. BMC Supplements"*

**Description:**


    Modularity optimization can be done using two algorithms:
        1: greedy for Clauset-Newman-Moore greedy modularity maximization
        2: louvain for Louvian Modularity optimization
    Input: Following files are needed
        edgelist: Tab spearated file of interactions between proteins or genes
    Output: 
        Original and refined node membership files, community size files, modularity values and refined modules for different iterations
    Source code: refmod.py

**Code Usage: Example one**

    import networkx as nx

    import numpy as np

    from networkx.algorithms.community import greedy_modularity_communities

    #Read edgelist as a graph using networkx

    G = nx.read_edgelist("edgelist.tsv", delimiter='\t', nodetype=str)

    method="louvain"

    iterations=np.array(list(range(1,6)))

    resolution=np.arange(1,3,1)

    #note: Iterations should be outer loop and resolution should be inner loop

    for i in iterations:

            for r in resolution:
    
                  ref_main(G,method,it=i,g=r,writeorig=True)
                  
                  

**Code Usage: Example two**

    import networkx as nx

    import numpy as np

    from networkx.algorithms.community import greedy_modularity_communities

    #Read edgelist as a graph using networkx

    G = nx.read_edgelist("edgelist.tsv", delimiter='\t', nodetype=str)

    method="greedy"

    ref_main(G,method,writeorig=True)
    
 
 
 **Data:**
 
 human_ppi :
    The interactions between proteins or genes are provided in a tab separated file. These interactions can be derived on the basis of gene expression, protein similarities, co-regulation or physical interactions between proteins. The edgelist for protein-protein interactions in human used for analysis published in publication is given in human_ppi. This file contains physical and functional interactions of human proteome collected from HPRD, STRING, BioGRID and IMEx consortium databases. It contains 78705 interactions in 12022 unique human proteins.
    

**Requirements:**

The following packages softwares are needed to run this algorithm successfully
    [1] Louvain Community Detection: The python implementation of Louvain algorithm [1] can be installed from pip using pip install python-louvain. For more details see https://github.com/taynaud/python-louvain .
    
    [2] Networkx
    
    [3] Python3

**References:**

    [1] Blondel, V. D., Guillaume, J.-L., Lambiotte, R. & Lefebvre, E. Fast unfolding of communities in large networks. J. Stat. Mech. theory Exp. 2008, P10008 (2008).
    
    [2] Clauset, A., Newman, M. E. J. & Moore, C. Finding community structure in very large networks. Phys. Rev. E 70, 66111 (2004).
