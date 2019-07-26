#refmod.py
"""
This module implements refmod, the module refinement algorithm.
Presently refmod offers refinement on modules observed through Lovain and Clauset-Newman-Moore greedy modularity maximization.
"""
import networkx as nx
import community
import math
import os
from itertools import chain
from networkx.algorithms.community import greedy_modularity_communities


__author__ = """Rama Kaalia (ramakaalia@gmail.com, rkaalia@ntu.edu.sg)"""
#    Copyright (C) 2019 by
#    Rama Kaalia <ramakaalia@gmail.com>
#    All rights reserved.


def partition_from_set(comset):
    #returns a dictionary for partition  from a list of frozensets  
    #get node_membership
    i=0
    partition={}
    for com in comset:
        i+=1
        for node in com:
            partition[node]=i      
    return partition


def get_modularity_set(graph,m):
    #Description:returns modularity value (float) for a partition 'm' of 'graph'
    #Parameters:
    #graph= supergraph or super module (networkx Graph object) and 
    #m= partition or set of modules(list of frozensets)
    graph_link = graph.size() #total edges in graph
    Q=0.
    for com in m:
        (din,dout) = module_degree(graph,com) #calculate indegree & outdegree of module
        dtot=din+dout
        m_subgraph = graph.subgraph(list(com)) #edges within the module
        m_link = m_subgraph.size()  #in degree of module
        Q += (m_link/graph_link) - (dtot/(2*graph_link))**2  #modularity of a module
    return Q

def greedy(G,gamma=1.0):
    #if method=greedy, modularization is done using this function
    c = list(greedy_modularity_communities(G))
    partition=partition_from_set(c)
    modularity=get_modularity_set(G,c)
    return partition,modularity

def louvain(G,gamma=1.0):
    #if method=louvain, modularization is done using this function
    partition= community.best_partition(G,randomize=True,resolution=gamma)
    best_mod = community.modularity(partition,G)
    return (partition,best_mod)

def module_degree(graph, sub):
    #Returns in_degree and out_degree of a module wrt a super-module
    #graph: supergraph or super module and 
    #sub: module (frozenset) whose in degree and out degree has to be calculated
    
    nodelist = list(sub)  #nodes of module
    module_graph = graph.subgraph(nodelist) #render the module as a graph
    total_degree = sum([val for (node, val) in graph.degree(nodelist)])  #total degree of module
    in_degree = sum([in_val for (in_node, in_val) in module_graph.degree(nodelist)])
    out_degree = total_degree - in_degree   
    return (in_degree,out_degree)

def get_modularity(graph,m):
    #graph= supergraph or super module and m= module (frozenset) whose modularity has to be calculated
    graph_link = len(nx.edges(graph))
    (din,dout) = module_degree(graph,m)
    dtot = din + dout
    m_subgraph = graph.subgraph(list(m))
    m_link = len(nx.edges(m_subgraph))
    m_Q = (m_link/graph_link) - (dtot/(2*graph_link))**2  #modularity of a module
    return m_Q

def compare(s, t): #compare two nonhashable lists of modules
    t = list(t)   # make a mutable copy
    try:
        for elem in s:
            t.remove(elem)
    except ValueError:
        return False
    return not t

def nodelist_from_multimodules(A): 
    #returns a unique list of nodes from list of communities (frozenset)
    B=[]
    for x in A:
        y=list(x)
        for z in y:
            B.append(z)
    return(list(set(B)))

def inc_modularize(giant,method,M,ref_stats):
    """
    Incremental modularization: The original partition is incrementally modularized untill there are no more refinable 
    modules left to be remodularized.
    
    """
    M_final=[]
    M_unref = []  #unrefined module set: modules that can be remodularized
    M_ref = []  #new refined module set : modules after remodularization
    M_non = [] #non-refinable module set : modules that can't be remodularized
    L = len(nx.edges(giant))
    
    i=1
    orig_mod = len(M)
    for m in M:
        #create networkx graph out of each community
        m_nodes = list(m)
        #Render first generation module as subgraph
        m_graph = giant.subgraph(m_nodes) 
        lm = len(nx.edges(m_graph)) 
        
        #Check for remodularization: A module is refinable if it clubs two or more sub-modules, 
        #thus having at least √2L intra-module edges
        
        if lm<math.sqrt(2*L):   
            M_non.append(m) 
            print("Module",i,"can't be remodularized!")
            print("...")
        else: 
            #Remodularize the modules which may include sub-modules 'M1'
            print("Submitting Module",i,"for remodularization...")
            print("...")
            mod_cmd=method+"(m_graph)" #if method=louvain, louvain(g) returns (partition,best modularity)
            modularization=eval(mod_cmd)
            partition1 = modularization[0]
            Q1 = modularization[1]
            #Next generation partition
            M1 = [frozenset([nodes for nodes in partition1.keys() if partition1[nodes] == module]) for module in set(partition1.values())] 
            print("Submodules in next generation:",len(M1))

            #Remodularization of next generation modules
            unref_submodule=[]
            ref_submodule=[]
            for mm in M1:
                (in_deg,out_deg)=module_degree(m_graph,mm)
                qmm = get_modularity(m_graph,mm)
                # Refine Modules if::
                # 1) Module has postive modularity and 2) intra module connections are more than intermodule connections
                # i.e. module is topologically true
                # Note: 
                if qmm >0 and out_deg - in_deg < out_deg*0.3:
                    ref_submodule.append(mm)
                else:
                    unref_submodule.append(mm)   
            print("Newly refined sub-modules:",len(ref_submodule))
            print("Modules not refined:",len(unref_submodule))
            print('...')
            print('...')
            #merge unrefined modules in a subgraph after the incremental step to submit again for remodularization
            if unref_submodule:
                unref_nodes= nodelist_from_multimodules(unref_submodule)
                unref_bigmodule=frozenset(unref_nodes)
                M_unref.append(unref_bigmodule)          
        
            if ref_submodule:           
                M_ref.append(ref_submodule)


        i+=1
    if M_non:
        M_final.append(M_non)
    if M_ref:
        M_ref = list(chain.from_iterable(M_ref))
        M_final.append(M_ref)
    if M_unref:
        #print("unref:",M_unref)
        M_final.append(M_unref)

    M_final = list(chain.from_iterable(M_final))
    print("Original_modules:",orig_mod,"non_refinable:", len(M_non),"unrefined:",len(M_unref),"new_refined:",len(M_ref))
    otp=str(orig_mod)+"\t"+str(len(M_non))+"\t"+str(len(M_unref))+"\t"+str(len(M_ref))+"\n"
    ref_stats.write(otp)
    #print(M_final)
    return M_final

def converge_inc_mod(giant,method,M1,M2,p,iterations,it_modfile,ref_stats):
    '''This function converges the algorithm i.e. finalizes partition M2 if loss in modularity of the partition q(M2)-q(M1) 
    is less than threshold modularity loss p. In case of no significant change in modularity of the partition, 
    if there are any unrefined modules in M2, remodularize them using inc_modularize function. Keep refining and 
    remodularizing unless the modularity loss is more than threshold p.
    
    '''
    
    if compare(M1, M2): #returns true if no change in two partitions
        #write change in modularity for every iteration/run
        mfin_sum=0
        for mfin in M2:
            mfin_Q = get_modularity(giant,mfin)
            mfin_sum += mfin_Q
        qm_final = mfin_sum
        it_output=str(iterations)+"\t"+str(qm_final)+"\n"
        it_modfile.write(it_output)
        print("Modularity after remodularization run",it_output)
        return M2
    else:      
        msum=0
        for m in M1:
            m_Q = get_modularity(giant,m)
            msum += m_Q
        qm = msum

        mfin_sum=0
        for mfin in M2:
            mfin_Q = get_modularity(giant,mfin)
            mfin_sum += mfin_Q
        qmfinal = mfin_sum
        print("change in modularity:",qm," to ",qmfinal)
        
        
        #if loss in modularity is less than threshold modularity loss
        if qm - qmfinal <= (qm*p*0.01):
            it_output=str(iterations)+"\t"+str(qmfinal)+"\n"
            it_modfile.write(it_output)
            print("Modularity after remodularization run ",it_output)
            
            iterations+=1
            print("Mining submodules within suboptimal region: run ",iterations)
            print("...")
       
            M2_buffer=[]
            M2_buffer = inc_modularize(giant,method,M2,ref_stats)
            M1=M2  
            M2=M2_buffer
            return(converge_inc_mod(giant,method,M1,M2,p,iterations,it_modfile,ref_stats))
        else:
            mfin_sum=0
            for mfin in M2:
                mfin_Q = get_modularity(giant,mfin)
                mfin_sum += mfin_Q
            qm_final = mfin_sum
            #write change in modularity for every iteration/run
            it_output=str(iterations)+"\t"+str(qm_final)+"\n"
            it_modfile.write(it_output)
            print("Modularity after remodularization run",it_output)
            return M2

def write_finalpartition(partition,i,r):
    """
    Write the final partition after refinement.
    Node membership in ref_membership file, community sizes in ref_community file 
    and modularity for different iterations of modularization of edgelist in ref_modularity file.
  
    Parameters
    ----------
    partition: dictionary
        Final partition after various runs of refinement on the original partition
        Partition is in the form of dictionary where keys are nodes and values are their communities.

    i: integer, iterations
        The iteration at which modules were detected using community detection method (louvain or greedy).
        Default is 1. The corresponding partition at this iteration is on 
    
    r: float, resolution
        The resolution for greedy, r=1. For louvain, resolution r can be set as any postive float value, default is set to 1.
    
    Returns
    -----
    ref_membership: tab separated file for nodes and their final refined communities
    ref_community: tab separated file for refined community and their sizes 
    
    """
    if i==1 and r==1:
        nodefile="./refmod_output/ref_membership"
        comfile="./refmod_output/ref_community"
    else:
        nodefile = "./refmod_output/ref_membership"+"_"+str(i)+"_"+str(r)
        comfile = "./refmod_output/ref_community"+"_"+str(i)+"_"+str(r)
    nodefh = open(nodefile,"w")
    comfh = open(comfile,"w")
    comfh.write("communityid\tsize\n")
    nodefh.write("nodeid\tcommunityid\n")
    
#get node_membership
    i=0
    for y in partition:
        i+=1
        for z in y:
            x1=str(z)+"\t"+str(i)+"\n"
            nodefh.write(x1)
#get community_size
    j=0
    for y in partition:
        j+=1
        x2=str(j)+"\t"+str(len(y))+"\n"
        comfh.write(x2)
    nodefh.close()
    comfh.close()

def write_originalpartition(modln,i,r):
    """
    Write the original partition. 
    Node membership in orig_membership file, community sizes in orig_community file 
    and modularity for different iterations of modularization of edgelist in orig_modularity file.
    
    Parameters
    ----------
    modln: tuple,(partition,modularity)
        Returned from community detection method louvain or greedy. 
        Partition is in the form of dictionary where keys are nodes and values are their communities.

    i: integer, iterations
        The iterations for detecting modules to be run using community detection method (louvain or greedy).
        Default is 1
    
    r: float, resolution
        The resolution for greedy, r=1. For louvain, resolution r can be set as any postive float value, default is set to 1.
    
    Returns
    -----
    orig_membership: tab separated file for nodes and their communities
    orig_community: tab separated file for community and their size
    orig_modularity: tab separated file for every partition.
        It reports iteration, resolution, modularity of partition, number of communities.
    """
    i_partition=modln[0]
    if i==1 and r==1:
        memfile="./refmod_output/orig_membership"
        comfile="./refmod_output/orig_community"
    else:
        memfile="./refmod_output/orig_membership"+"_"+str(i)+"_"+str(r)
        comfile="./refmod_output/orig_community"+"_"+str(i)+"_"+str(r)

    modfile=open("./refmod_output/orig_modularity","a")
    modfile.write("iteration"+"\t"+"resolution"+"\t"+"modularity"+"\t"+"communities"+"\n")


        #stats from original algorithm for different iterations
    
    file1=open(comfile,"w") #size of different communities for best partition
    file2=open(memfile,"w") #memberships of nodes in best partition
    #print number and size of communities

    size = float(len(set(i_partition.values())))
    
    file1.write("communityid\tsize\n")
    #print ("\nCommunity","Size")

    for com in set(i_partition.values()):
        list_nodes = [nodes for nodes in i_partition.keys()
                      if i_partition[nodes] == com]
        com_size = float(len(list_nodes))
        x1=str(com)+"\t"+str(com_size)+"\n"
        file1.write(x1)

    file2.write("nodeid\tcommunityid\n")
    #print memberships
    for nodes in set(i_partition.keys()):
        com= i_partition[nodes]
        #print (nodes,com)
        x2=str(nodes)+"\t"+str(com)+"\n"
        file2.write(x2)

    file1.close()
    file2.close()
    best_mod = modln[1]
    x3=str(i)+"\t"+str(r)+"\t"+str(best_mod)+"\t"+str(size)+"\n"
    modfile.write(x3)
    modfile.close()
        
##Main module
def ref_main(g1,method,it=1,g=1.0,p=30.0,writeorig=False):
    """
    Modularize the network using a modularity based community detection method 
    and relax the modularity levels, remodularize & refine the observed modules to obtain smaller modules.
    
    Parameters
    ----------
    g1: networkx.Graph
        The graph which has to be partitioned into communities/modules
    
    method: louvain | greedy
        Modularity optimization can be done using two algorithms:
        method=greedy for Clauset-Newman-Moore greedy modularity maximization
        method=louvain for Louvian Modularity optimization
        
    it: integer
        The index of iteration at which the modules were detected using community detection method (louvain or greedy).
        Default is 1. It's useful in case of stochasticity based modularization such as louvain 
        where each run gives different partition and multiple iterations of modularization are required.
        
    g: float, resolution
        The resolution for greedy, r=1. 
        For louvain, resolution r can be set as any postive float value, default is set to 1.
        
    p: float
    The threshold modularity loss in percentage, default value is 30 for 30% loss
        
    writeorig: True|False, optional
        If set to True, the original partition are also written in output files orig_membership, orig_community
        and orig_modularity. Default is set to None.


    """
    #returns original louvain partition and modularity if writeorig=True
    
    #Final partition
    M_f=[]  #list of communities(frozenset)

    #get First Generation Modules
 
    #output files
    #create output directory and delete existing files
    if os.path.exists("./refmod_output"):       
        flist = os.listdir("./refmod_output")
        # Remove old output files if any
        for filePath in flist:
            try:
                if os.path.exists(filePath):
                    os.remove(filePath)
            except:
                print("Error while deleting old file : ", filePath)    
    else:
        os.mkdir("./refmod_output")
                           
    ref_stats=open("./refmod_output/ref_stats","a") #refined and nonrefined modules at every refinement step
    hdr="iteration="+str(it)+" resolution="+str(g)+"\n"
    ref_stats.write(hdr)
    ref_stats.write("original\tnon_refinable\tunrefined\tnew_refined\n")
    #modularities for refined communities at every refinement run for each iteration at resolution g
    it_modfile=open("./refmod_output/ref_modularity","a")
    it_modfile.write(hdr)
    it_modfile.write("run\tmodularity\n")

    giant= max(nx.connected_component_subgraphs(g1), key=len) #main component of graph to be modularized

    print("giant component has ", giant.__len__(), " nodes out of total ",g1.__len__()," nodes")
    
    #First modularization
    mod_cmd=method+"(giant,gamma=g)" #if method=louvain, louvain(g) returns (partition,best modularity)
    modularization=eval(mod_cmd)
    partition = modularization[0]
    Q = modularization[1]

    #first generation partition
    M = [frozenset([nodes for nodes in partition.keys() if partition[nodes] == module]) for module in set(partition.values())] 
    print("Modules in first generation:",len(M))
    print('...')
    
    #Incremental modularization
    global iterations
    iterations=1
    print("Mining submodules within suboptimal region: run ",iterations)
    print("...")
    M_f=inc_modularize(giant,method,M,ref_stats)
    M_f = converge_inc_mod(giant,method,M,M_f,p,iterations,it_modfile,ref_stats)

    write_finalpartition(M_f,it,g)
    
    it_modfile.close()
    if writeorig==True:
        write_originalpartition(modularization,it,g)




