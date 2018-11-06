#############################################################################################
##############################################################################################
# Back Tracing Algorithm
# updated: 
#        2017/11/30: finished code, no documentation yet
#        2017/12/01: documentation
#        2017/12/05: documentation and explanation of the code to the untrained eye
#############################################################################################
# Author: Abigail Horn
# Transferred from Matlab to Python for BfR: Marcel Fuhrmann
#
# Functionality: 
#        Given a network model and a number of cases, find the best estimate of
#        one true source based on Horn&Friedrich(2017):
#        "Locating the Source of Large-Scale Outbreaks of Foodborne Disease"
# INPUT: 
#        Network model: 
#            file: data format CSV-table: 
#                        rows original nodes, 
#                        columns destination nodes
#        Observations: (e.g. outbreak case data)
#            file: data format CSV-table:
#                        1st row: node ids, order of rows from network model
#
# OUTPUT: 
#        Feasible nodes as text:
#            prints a row of node ids of possible sources, no ranking
#        best estimate of true source:
#            gives a best estimate of true source and probability, given 1 true source
#
#    TODO:    
#        output as file, 
#        distribution of source probabilities, 
#        ranking of sources,
#                need info first what is needed or wanted
#        GUI?
#                for user who doesn't what to look at code, maybe transfer to KNIME?
#        prior knowledge
#                now the codes ignores all expert knowledge and number of observations
#        assumption is that all observations can reach the true source through the network
#                thereby ignoring all expert knowledge that may point to other sources: problem?
#############################################################################################
#############################################################################################
# IMPORT TOOLS
#############################################################################################
import numpy as np  # numerical analysis package
import pandas as pd # data reading

from functions import parentNodes, exactVolumeComponent







#############################################################################################
# DEFINE
# PARAMETERS         AND
# DATA SOURCE  
# TODO: nrofnodesinlayer should be readable by data itself, have to come up with something here
#############################################################################################
nrofnodesinlayer = [10, 5, 5, 10];

# gives the network model
df1=pd.read_csv('test1/flows.dat', sep=',',header=None)
# gives the outbreak cases
df2=pd.read_csv('test1/contam_reps.dat', sep=',',header=None)






#############################################################################################
# PREPROCESSING
# PARAMETERS DERIVED FROM DATA SOURCE AND PARAMETERS 
#############################################################################################

# total number of layers
nroflayers = len(nrofnodesinlayer)

# network flows read in from data 
flows = df1.values   
flows = np.array(flows) # transform into matrix!

# contaminated nodes read in from data
obs = df2.values  
obs = np.array(obs)     # transform into matrix!
obs = obs.astype(int)   # make it integers!

# number of nodes from flows
nrofnodes = flows.shape[0]  

# dont forget that in python index always starts with zero
obsID = obs[0]-1            

# ordered and unique set of contaminated nodes, have to be integer
setofContamNodes = [int(i) for i in sorted(set(obsID))]
print("IDs of contaminated nodes:",obsID)
print("Total number of nodes:",nrofnodes)

# creating the list for the last node of each layer
adresses = [0]*nroflayers
for i in range(nroflayers):
    adresses[i]=sum(nrofnodesinlayer[0:i+1])    #TODO: remind me why that is here
print("ID of last node of each layer:",adresses)

# defining 1st and last node of each layer 
firstNode = [0]*nroflayers
lastNode = [0]*nroflayers
for i in range(nroflayers):
    firstNode[i] = adresses[i] - nrofnodesinlayer[i]
    lastNode[i] = adresses[i]

# Function: sparsing flow on my own    
# (Deactivated: small matrices dont need to be sparsed, needs to be imported!)
#sFlow = sparseFlow(flows,nrofnodes)

# ids of maximum number of feasible sources aka. nodes of the 1st layer
# (Deactivated: using only feasible source, not all sources)
#firstLayerSources = np.arange(firstNode[0],lastNode[0])



#############################################################################################
# CONSISTENCY CHECK:
# IS DATA COMPATIBLE WITH DEFINED PARAMETERS?
#############################################################################################

# CHECK 1: checking for interaction between nodes of the same layer, 
# if data is consistent with assumption this should not happen
for i in range(nroflayers):
    if flows[firstNode[i]:lastNode[i],firstNode[i]:lastNode[i]].sum()!=0:
        print("on layer nr.",i+1,": interaction between nodes of the same layer detected!")

# CHECK 2: Is total number of nodes from flows consistent with the sum of the nodes in each layer?
if sum(nrofnodesinlayer)!=flows.shape[0]:
    print("Total number of nodes in file is not equal sum of nodes in each layer!")
#############################################################################################
# CORE PROGRAM
#############################################################################################
#
# TODO: dont forget about prior knowledge
#
#############################################################################################


# Function(recursive): Finding Parent Nodes
# (may be deactivated if number of nodes is small or all 1st layer nodes are feasible)
# note that i only look for feasible 1st-layer-nodes!
candidateNodes = np.array(setofContamNodes)
for i in range(nroflayers-1):
    candidateNodes = parentNodes(candidateNodes,flows)
print("Found all feasible nodes:",[x+1 for x in candidateNodes])

# Function: Back Tracing
pmf = exactVolumeComponent(candidateNodes, obsID, flows, firstNode, lastNode, nroflayers)
print("Most probable source is Node Nr.",np.argmax(pmf)+1,"with probability of ",round(np.amax(pmf)*100.,2),"%!")
#print("probability of each node of the 1st layer:", pmf)


