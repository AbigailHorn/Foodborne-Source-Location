#############################################################################################
# IMPORT TOOLS
#############################################################################################
import numpy as np
from numpy.linalg import inv

#############################################################################################
# TASK:    sparse flow matrix
#           create table with ID of 
#            outgoing node, incoming node, and flow probability
#
# INPUT:   network flow matrix <flows>
#           Total number of nodes <nrofnodes>
#
# OUTPUT:  table with size (nrofnodes,3) <sparsedFlows>
#############################################################################################
def sparseFlow(flows,nrofnodes):
    nrofLinks = np.count_nonzero(flows)
    sparsedFlows = np.zeros((nrofLinks,3))
    
    running = 0
    for i in range(nrofnodes):
        for j in range(nrofnodes):
            if flows[i,j]>0:
                sparsedFlows[running,0] = i
                sparsedFlows[running,1] = j
                sparsedFlows[running,2] = flows[i,j]
                running += 1
                
    return sparsedFlows

#############################################################################################
# TASK:    find the corresponding parent nodes to list of given nodes
#           
# INPUT:   network flow matrix <flows>
#           list of nodes from flow matrix <setofNodes>
#
# OUTPUT:  list of parent nodes <setofParentNodes>
#############################################################################################
def parentNodes(setofNodes,flows):
    setofParentNodes = []
    
    for i in range(len(setofNodes)):
        parentNodes = np.argwhere(flows[:,setofNodes[i]] > 0).tolist()
        setofParentNodes.extend(parentNodes)
        # maybe later use?
        #nrofParentLinks = np.count_nonzero(flows[:,setofNodes[i]])
    
    
    # to be sorted and filter only unique nodes, we need to prep the list
    setofParentNodes = np.array(setofParentNodes).reshape(len(setofParentNodes),)
    
    # ordering the list 
    setofParentNodes = [int(j) for j in sorted(set(setofParentNodes))]
    return setofParentNodes


#############################################################################################
# TASK:    calculate the probability distribution of feasible sources to be true source
#           Formula is A = (I_t - Q)^{-1}*R
#
# INPUT:   list of nodes, possible sources <feasibleSources>
#           list of absorbing nodes, outbreak cases <contamReports>
#            network flow matrix <flows>
#             ID of the first node of each layer <firstNode>
#              ID of the last node of each layer <lastNode>
#               total number of layers <nroflayers>
#
# OUTPUT:  normalized probability distribution of feasible sources <pmf>
#############################################################################################
def exactVolumeComponent(feasibleSources, contamReports, flows, firstNode, lastNode, nroflayers):
    
    pmf = np.zeros(len(feasibleSources))

    # size: (transient nodes, transient nodes) 
    I_t = np.identity(lastNode[nroflayers-2])
    
    # size: (transient nodes, transient nodes) 
    Q = flows[firstNode[0]:lastNode[nroflayers-2],firstNode[0]:lastNode[nroflayers-2]]
    
    # size: (transient nodes, absorbing nodes)
    R = flows[firstNode[0]:lastNode[nroflayers-2],firstNode[nroflayers-1]:lastNode[nroflayers-1]]
       
    # size: (transient nodes, absorbing nodes)
    A = np.dot(inv(I_t - Q),R)
    
    
    for i in range(len(feasibleSources)):
        pathLikelihood = 1.
       
        if contamReports.size>1:    #catching numpy doesnt like lists of length 1
             for j in range(contamReports.size):
                pathLikelihood *= A[i,contamReports[j]-lastNode[nroflayers-2]]
    
        else:
            pathLikelihood *= A[i,contamReports-lastNode[nroflayers-2]] 
        pmf[i] = pathLikelihood    
        
    #normalize
    pmf /= pmf.sum()
    
    return pmf