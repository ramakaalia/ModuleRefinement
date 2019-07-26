**RefMod: An Algorithm for Refining Modules in biological networks**

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
        Download and save the code as refmod.py.

**Code Usage: Example one**


    import networkx as nx
    
    import numpy as np
    
    import refmod

    from networkx.algorithms.community import greedy_modularity_communities

    #Read edgelist as a graph using networkx

    G = nx.read_edgelist("samplefile", delimiter='\t', nodetype=str)

    method="louvain"

    iterations=np.array(list(range(1,6)))

    resolution=np.arange(1,3,1)

    #note: Iterations should be outer loop and resolution should be inner loop

    for i in iterations:

            for r in resolution:
    
                  refmod.ref_main(G,method,it=i,g=r,writeorig=True)
                  
                  

**Code Usage: Example two**

    import networkx as nx

    import numpy as np
    
    import refmod

    from networkx.algorithms.community import greedy_modularity_communities

    #Read edgelist as a graph using networkx

    G = nx.read_edgelist("samplefile", delimiter='\t', nodetype=str)

    method="greedy"

    refmod.ref_main(G,method,writeorig=True)
    
 
 **Note:**
 
    The input for the code is a tab separated edgelist. The output of the code is written in a refmod_output directory (Warning: Any existing files in this directory are rewritten at each run). If the option writeorig is set to True, the original partition are also written in output files orig_membership, orig_community and orig_modularity. Default is set to None. 
    
    Output files for original partitioning from Louvain or Greedy modularity maximization are:
    orig_membership: tab separated file for nodes and their communities
    orig_community: tab separated file for community and their size
    orig_modularity: tab separated file for every partition.
    
    Output files for refined modules are:
    ref_membership: tab separated file for nodes and their final refined communities
    ref_community: tab separated file for refined community and their sizes 
    ref_modularity: tab separated file for modularities for refined communities at every remodularization step
    ref_stats: refined and nonrefined modules at every remodularization step of the refinement algorithm
 
    

**Requirements:**

The following packages softwares are needed to run this algorithm successfully
    [1] Louvain Community Detection: The python implementation of Louvain algorithm [1] can be installed from pip using pip install python-louvain. For more details see https://github.com/taynaud/python-louvain .
    
    [2] Networkx
    
    [3] Python3


**Data:**
 
 human_ppi :
    The interactions between proteins or genes are provided in a tab separated file. These interactions can be derived on the basis of gene expression, protein similarities, co-regulation or physical interactions between proteins. The edgelist for protein-protein interactions in human used for analysis published in publication is given in human_ppi. This file contains physical and functional interactions of human proteome collected from HPRD, STRING, BioGRID and IMEx consortium databases. It contains 78705 interactions in 12022 unique human proteins.
 
 Partition files used for comparing different algorithms in paper: *"Refining modules to determine functionally significant clusters in molecular networks. Rama Kaalia and Jagath C. Rajapakse. BMC Supplements"*
 
 louvain_membership: The best iteration partition (out of 25 iterations) from louvain community detection at resolution parameter 2 (out of 1 to 10).
 
 louvref_membership: The partition after refinement of the best iteration partition from louvain community detection at resolution parameter 2 (i.e. from louvain_membership).
  
 greedy_membership: The best partition from greedy community detection.
 
 greedyref_membership: The partition after refinement of the best iteration partition from greedy community detection (i.e. from greedy_membership).
    
 asy_membership: The best partition from Asymptotic Surprise community detection.
      
 mcode_membership: The best partition from MCODE community detection.
 
 dpclus_membership: The best partition from DPCLUS community detection.
 
 labelprop_membership: The best partition from Label Propagation community detection.
 
 
    
**References:**

    [1] Blondel, V. D., Guillaume, J.-L., Lambiotte, R. & Lefebvre, E. Fast unfolding of communities in large networks. J. Stat. Mech. theory Exp. 2008, P10008 (2008).
    
    [2] Clauset, A., Newman, M. E. J. & Moore, C. Finding community structure in very large networks. Phys. Rev. E 70, 66111 (2004).
