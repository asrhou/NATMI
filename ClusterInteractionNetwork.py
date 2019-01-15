#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 25 13:57:26 2018

@author: rhou
"""

import argparse
import pandas as pd
import numpy as np
import os, sys, json
import multiprocessing

#===================================================================#
#----------------------------MCL  METHOD----------------------------#
#===================================================================#

def Normalize(A):
    columnSums = A.sum(axis=0)
    for columnIndex in xrange(len(columnSums)):
        if columnSums[columnIndex] != 0:
            for rowIndex in xrange(len(A)):
                A[rowIndex,columnIndex] = A[rowIndex,columnIndex] / columnSums[columnIndex]
    newMatrix = np.nan_to_num(A)
    return newMatrix

# get clusters from the MCL output
def GetMCLClusters(A, nodenames):
    clusters = []
    for i, r in enumerate((A > 0).tolist()):
        if r[i]:
            clusters.append(A[i,:]>0)
    rstDict = {}
    clusterIndex = 0
    for cn, c in enumerate(clusters):
        nodeList = []
        for x in  [ i for i, x in enumerate(c) if x ]:
            nodeList.append(nodenames[x].split('\n')[0])
        if nodeList not in rstDict.values():
            clusterIndex += 1
            rstDict['Cluster %02d' % clusterIndex] = nodeList
    return rstDict

#MCL method
def MCL(adjM, expandPara, inflatePara, maxLoop):
    # get weights
    adjArray = adjM.values[:, :].T
    # add self loops
    adjArray = adjArray + np.identity(len(adjArray))
    # normalize values
    adjArray = Normalize(adjArray)
    
    for i in xrange(maxLoop):
        # expansion step
        adjArray = np.linalg.matrix_power(adjArray, expandPara)
        # inflation step
        adjArray = np.power(adjArray, inflatePara)
        # normalize values
        adjArray = Normalize(adjArray)
        
        # check if the matrix is converged
        if i%5 == 4: 
            temp = np.linalg.matrix_power(adjArray, expandPara) - adjArray
            if np.max(temp) - np.min(temp) == 0:
                print "Stop at iteration %s" % i
                break
    cltDict = GetMCLClusters(adjArray, adjM.columns)
    return cltDict

#===================================================================#
#--------------------------Louvain  METHOD--------------------------#
#===================================================================#

#generate node-community vector from node-neighbour matrix
def GenLVNodeCommunityV(nodeNeighM):
    #each community is a node
    nodeNumber = len(nodeNeighM)
    nodeCommV = np.arange(nodeNumber)
    
    #change node names to community names
    nodeNeighM.index = nodeCommV
    nodeNeighM.columns = nodeCommV
    #change nodeCommV will affect the nodeNeighM, use the individual copy to replace the nodeCommV
    nodeCommV = nodeCommV.copy()
    return nodeCommV, nodeNeighM
    
#check which nodes in the new communities
def CheckLVNodesOfComm(communtyIndexList, nodeCommV):
    mappingList = []
    for communtyIndex in communtyIndexList:
        nodes = np.where(nodeCommV == communtyIndex)[0]
        mappingList.append(nodes)
    return mappingList
    
#prepare the neighbours for each community(node), values are neighbours' index names and itself
def CheckLVNeighbours(nodeNeighM, communtyIndexList):
    #serch neighbours for each node
    commNeighbourList = []
    for commIndex in communtyIndexList:
        #non-zero values in the matrix are neighbours
        hneighbours = nodeNeighM.ix[commIndex,:] > 0
        vneighbours = nodeNeighM.ix[:,commIndex] > 0
        
        #elements in neighbours are bool, elements with true are wanted
        hneighbourIndexList = nodeNeighM.index[hneighbours].values
        vneighbourIndexList = nodeNeighM.index[vneighbours].values
        
        #merge the indices of nodes whose tail or head are the given node and itself
        neighboursIndexList = np.append(hneighbourIndexList, vneighbourIndexList)
        neighboursIndexList = np.append(neighboursIndexList, commIndex)
        
        #remove duplicate indices
        neighboursIndexList = np.unique(neighboursIndexList)
        commNeighbourList.append(neighboursIndexList)
        
    return commNeighbourList
    
#calculate the sum of the weights of the links inside community
def CalSumInComm(nodesInComm, nodeNeighM):
    return nodeNeighM.values[nodesInComm, nodesInComm].sum()
    
#calculate the sum of the weights of the links incident to nodes in community
def CalSumTot(nodesInComm, nodeNeighM):
    sumTotPara = 0.0
    #calculate the sum of the weights of the links that given nodes are heads
    sumTotPara += nodeNeighM.values[nodesInComm, :].sum()
    #calculate the sum of the weights of the links that given nodes are tails
    sumTotPara += nodeNeighM.values[:,nodesInComm].sum()
    return sumTotPara
    
#calculate the modularity
def CalModularity(nodeCommV, nodeNeighM, mPara):
    #initially modularity is 0.0
    mod = 0.0
    
    #find all non-empty communities
    commIndexList = np.unique(nodeCommV)
    
    #calculate the modularity of each commmunity
    for commIndex in commIndexList:
        #find nodes in the community
        nodeIndexList = np.where(nodeCommV == commIndex)[0]
        
        #calculate the sum weight of edges incident of nodes 
        totPara = CalSumTot(nodeIndexList, nodeNeighM)
        
        #add to the total modularity
        if  totPara > 0:
            #calculate the sum weight of edges in the community
            weightInComm = CalSumInComm(nodeIndexList, nodeNeighM)
            #update the sum of modularity
            mod += np.divide(weightInComm, mPara) - np.power(np.divide(totPara, 2.0 * mPara), 2)
            
    return (mod, nodeCommV)
    
#calculate the modularity of node in each possible community
def CalModularities(commOrder, neighbourCommList, currentCommunityPartion, nodeNeighM, coreNumToUse, mPara):
    #generate new partition
    partitionList = []
    for neighbourComm in neighbourCommList:
        newPartion = currentCommunityPartion.copy()
        newPartion[commOrder] = neighbourComm
        partitionList.append(newPartion)
    
    #set how many cores are used to calculate modularity, one core can calculate one kind of partition individually
    from multiprocessing import Pool
    pool = Pool(coreNumToUse)
    
    #map method can only take one argument, wrap the CalModularity with one fixed argument
    from functools import partial
    #prod_x has only one argument newPartion, nodeNeighM and mpapa are fixed
    prod_x=partial(CalModularity, nodeNeighM=nodeNeighM, mPara=mPara)
    modularityPartionList = pool.map(prod_x, partitionList)
    pool.close()
    pool.join()
    
    return modularityPartionList

# a greedy assignment of nodes to communities, favoring local optimizations of modularity
def GreedyLVExtending(nodeCommV, nodeNeighM, mPara, coreNumToUse, isRandom):
    
    #initialise the arguments
    roundCounter = 0
    
    #find the communities
    communtyIndexList = nodeNeighM.index.values
    
    #build the mapping between coarsen graph and original graph
    mappingList = CheckLVNodesOfComm(communtyIndexList, nodeCommV)
    coarsenCommPartitionV = nodeNeighM.index.values
    
    #record the current partition for comparing
    currentCommunityPartion = coarsenCommPartitionV
    
    #prepare the neighbours for each community
    commNeighbourList = CheckLVNeighbours(nodeNeighM, communtyIndexList)
    
    #limit the maximum steps of improvement
    while roundCounter < communtyIndexList.size:
        
        ##initial parameters
        #number of round has been done
        roundCounter += 1
        
        #check if there is improvement
        improvement = False
        
        #try to move nodes in given order 
        communties = np.arange(communtyIndexList.size)
        if isRandom:
            import random
            random.shuffle(communties)
        
        #try to merge two communities in given order
        for commOrder in communties:
            
            #check which communities the node can move to based on its neighbours, including the original one
            neighbourCommList = commNeighbourList[commOrder]
            
            #calculate the modularity of node in each possible community
            modularityPartionList = CalModularities(commOrder, neighbourCommList, currentCommunityPartion, nodeNeighM, coreNumToUse, mPara)
            
            #find the partion with the biggest modularity
            from operator import itemgetter
            #bestPartionPair is a tuple (mod, nodeCommV)
            bestPartionPair = max(modularityPartionList, key=itemgetter(0))
            bestPartion = bestPartionPair[1]
            bestModularity = bestPartionPair[0]
            
            #check if best partition is not the origignal one, the modularity is improved
            if not improvement and not np.array_equal(bestPartion, currentCommunityPartion):
                improvement = True
                
            #use the best partition for now to next round
            currentCommunityPartion = bestPartion
            currentModularity = bestModularity
                
        #if there is no improvement in this round, then stop. 
        if not improvement:
            break
    
    #update the new partition to the original graph
    for commOrder in xrange(currentCommunityPartion.size):
        affectedNodes = mappingList[commOrder]
        commIndex = currentCommunityPartion[commOrder]
        nodeCommV[affectedNodes] = commIndex
    return nodeCommV, currentModularity

#building a new network whose nodes are now the communities found during the first phase
def CoarsenLVGraph(currentNodeCommV, currentNodeNeighM):
    #find communityIndex, use them as the new graph's node names
    communityIndexList = np.unique(currentNodeCommV)
    newNodeNum = communityIndexList.size
    iniExist = np.zeros((newNodeNum, newNodeNum),dtype=np.float)
    newNodeNeighM = pd.DataFrame(iniExist, index=communityIndexList, columns=communityIndexList)
    
    #find nodes in the community
    commNodeList = CheckLVNodesOfComm(communityIndexList, currentNodeCommV)
    
    #calculate the sum weight of edges from one community to another one
    for headIndex in xrange(newNodeNum):
        for tailIndex in xrange(newNodeNum):
            headSet = commNodeList[headIndex]
            tailSet = commNodeList[tailIndex]
            sumValue = currentNodeNeighM.iloc[tailSet, headSet].sum().sum()
            newNodeNeighM.iloc[tailIndex, headIndex] = sumValue
            
    return newNodeNeighM

#greedy search partition
def GreedyLVSearch(nodeNeighM, coreNumToUse, isRandom):    
    level = 0
    
    #initial partition is each node is a community, which is the best one for now
    bestNodeCommV, nodeNeighM = GenLVNodeCommunityV(nodeNeighM)
    originalNodeNeighM = nodeNeighM.copy()
    
    #m is the sum of the weights of all the links, which will not changed in coarsen graphs
    mPara = nodeNeighM.unstack().sum()
        
    #greedy improve the modularity
    while 1:
        level += 1
        print 'level %s' % level
        
        #reassign nodes to communities based on the modularity
        newNodeCommV, newMod = GreedyLVExtending(bestNodeCommV.copy(), nodeNeighM.copy(), mPara, coreNumToUse, isRandom)
        
        #if the partition is firstly extended
        if level == 1:
            #check if there are changes, which means new partition generated
            if np.array_equal(bestNodeCommV, newNodeCommV):
                #there are no changes, return the initial partition
                return bestNodeCommV
            #new partition generated
            else:
                #there are changes, save as the best partition
                bestNodeCommV = newNodeCommV
                bestMod = newMod
                #if the best partition has two or less communities, terminate the extending
                if np.unique(bestNodeCommV).size <= 2:
                    return bestNodeCommV
                    
                #coarsen graph based on the new partitions in the originalNodeCommM
                nodeNeighM = CoarsenLVGraph(bestNodeCommV, originalNodeNeighM)
        else:
            #if there is no new partition generated, no improvement or new partition is worse, finish the searching
            if np.array_equal(bestNodeCommV, newNodeCommV) or newMod < bestMod:
                return bestNodeCommV
            #if the new partition is better, replace the current one and coarsen the graph
            else:
                #save the best partition for now
                bestNodeCommV = newNodeCommV
                bestMod = newMod
                
                #if the best partition has two or less communities, terminate the extending
                if np.unique(bestNodeCommV).size <= 2:
                    return bestNodeCommV
                    
                #coarsen graph based on the new partitions in the originalNodeCommM
                nodeNeighM = CoarsenLVGraph(bestNodeCommV, nodeNeighM)

#find the cluster list from the partion vector 
def GetModularityClusters(nodeCommV, nodenames):
    clusterIndexList = np.unique(nodeCommV)
    rstDict = {}
    clusterIndex = 0
    for clusterIndex in clusterIndexList:
        nodeList = list(np.where(nodeCommV == clusterIndex)[0])
        nodeList = [nodenames[i].split('\n')[0] for i in nodeList]
        if nodeList not in rstDict.values():
            clusterIndex += 1
            rstDict['Cluster %02d' % clusterIndex] = nodeList
    return rstDict

# Louvain method
def Louvain(adjM, coreNum, isRandom):
    nodeNames = adjM.columns
    greedyNodeCommV = GreedyLVSearch(adjM, coreNum, isRandom)
    rstDict = GetModularityClusters(greedyNodeCommV, nodeNames)
    return rstDict
    
#===================================================================#
#--------------------------Infomap  METHOD--------------------------#
#===================================================================#

#generate node-community vector from node-neighbour matrix
def GenIMNodeCommunityV(nodeNeighM):
    #each community is a node
    nodeNumber = len(nodeNeighM)
    nodeCommV = np.arange(nodeNumber)
    
    #change node names to community names
    nodeNeighM.index = nodeCommV
    nodeNeighM.columns = nodeCommV
    #change nodeCommV will affect the nodeNeighM, use the individual copy to replace the nodeCommV
    nodeCommV = nodeCommV.copy()
    return nodeCommV, nodeNeighM
    
#calculate the ergodic node visit frequencies of nodes, using PageRank algorithm
def CalErgodicNodeVisitFrequencies(nodeNeighM, tauPara):
    nodeNum = nodeNeighM.index.size
    #iniPPara is nodeNum rows and one column
    iniPPara = np.ones((nodeNum, 1))
    
    #normalise the memberships so the sum is one
    nodeNeighM = nodeNeighM.div(nodeNeighM.sum(axis=1), axis=0).fillna(0.0)
    
    #transitison matrix
    nodeNeighArray = nodeNeighM.values[:, :].T
    
    #start p matrix
    startPM = np.full((nodeNum, nodeNum), 1.0 /nodeNum, dtype = np.float)
    
    #flow matrix
    partOne = np.multiply((1 - tauPara), nodeNeighArray)
    partTwo = np.multiply(tauPara, startPM)
    flowM = np.add(partOne, partTwo)
    
    pPara = np.dot(flowM, iniPPara)
    
    while not np.array_equal(pPara, iniPPara):
        iniPPara = pPara 
        pPara = np.dot(flowM, iniPPara)
    pPara = pPara.T
    
    #pPara is a matrix, reduce it to a vector
    return pPara[0]
    
#check which nodes in the new communities
def CheckIMNodesOfComm(communtyIndexList, nodeCommV):
    mappingList = []
    for communtyIndex in communtyIndexList:
        nodes = np.where(nodeCommV == communtyIndex)[0]
        mappingList.append(nodes)
    return mappingList
    
#prepare the neighbours for each community(node), values are neighbours' index names and itself
def CheckIMNeighbours(nodeNeighM, communtyIndexList):
    #serch neighbours for each node
    commNeighbourList = []
    for commIndex in communtyIndexList:
        #non-zero values in the matrix are neighbours
        hneighbours = nodeNeighM.ix[commIndex,:] > 0
        vneighbours = nodeNeighM.ix[:,commIndex] > 0
        
        #elements in neighbours are bool, elements with true are wanted
        hneighbourIndexList = nodeNeighM.index[hneighbours].values
        vneighbourIndexList = nodeNeighM.index[vneighbours].values
        
        #merge the indices of nodes whose tail or head are the given node and itself
        neighboursIndexList = np.append(hneighbourIndexList, vneighbourIndexList)
        neighboursIndexList = np.append(neighboursIndexList, commIndex)
        
        #remove duplicate indices
        neighboursIndexList = np.unique(neighboursIndexList)
        commNeighbourList.append(neighboursIndexList)
        
    return commNeighbourList
    
#calculate exit probabilities
def CalExitProbab(pPara, tauPara, cellListOfCommunities, nodeNeighM):
    #normalise the memberships so the sum is one
    nodeNeighM = nodeNeighM.div(nodeNeighM.sum(axis=1), axis=0).fillna(0.0)
    
    exitProbabilities = []
    nodeNum = nodeNeighM.index.size
    #for each community
    for nodeList in cellListOfCommunities:
        partOneOne = np.divide(np.multiply(tauPara, (nodeNum - nodeList.size)), (nodeNum - 1))
        partOneTwo = pPara[nodeList].sum()
        partOne = np.multiply(partOneOne, partOneTwo)
        
        partTwoOne = 1.0 - tauPara
        partTwoTwo = 0.0
        outsideNodes = np.setdiff1d(np.arange(nodeNum), nodeList)
        for inNodeIndex in nodeList:
            for outsideNodeIndex in outsideNodes:
                partTwoTwo += np.multiply(pPara[inNodeIndex], nodeNeighM.iloc[inNodeIndex, outsideNodeIndex])
        partTwo = np.multiply(partTwoOne, partTwoTwo)
        
        iExitProb = partOne + partTwo
        exitProbabilities.append(iExitProb)
        
    #wrap up as a ndarray
    return np.asarray(exitProbabilities)
    
#calculate the map equation
def CalMapEquation(nodeCommV, nodeNeighM, pPara, tauPara):
    #find all non-empty communities
    commIndexList = np.unique(nodeCommV)
    
    #find the nodes in each community
    cellListOfCommunities = CheckIMNodesOfComm(commIndexList, nodeCommV)
    
    #calculate exit probabilities
    exitProbabilities = CalExitProbab(pPara, tauPara, cellListOfCommunities, nodeNeighM)
    nonZeroExitProbabilities = exitProbabilities[exitProbabilities>0]
    
    #for red part
    sumExitProbabilities = exitProbabilities.sum()
    if sumExitProbabilities > 0:
        redPart = np.multiply(sumExitProbabilities, np.log2(sumExitProbabilities))
    else:
        redPart = 0
    
    #for black part
    sumExitProbLog = 0.0
    loggedNonZeroExitProbabilities = np.log2(nonZeroExitProbabilities)
    multiplyLogExitProb = np.multiply(loggedNonZeroExitProbabilities, nonZeroExitProbabilities)
    sumExitProbLog = multiplyLogExitProb.sum()
    blackPart = np.multiply(2.0, sumExitProbLog)
    
    #for blue part one
    loggedPPara = np.where(pPara > 1e-15, np.log2(pPara), 0.0)
    multiplyLoggedPPara = np.multiply(loggedPPara, pPara)
    sumPParaLog = multiplyLoggedPPara.sum()

    #for blue part two     
    sumExitProbPParaLog = 0.0
    for commOrder in xrange(len(cellListOfCommunities)):
        itemExitProbPParaLog = pPara[cellListOfCommunities[commOrder]].sum() + exitProbabilities[commOrder]
        loggedItemExitProbPParaLog = np.where(itemExitProbPParaLog > 1e-15, np.log2(itemExitProbPParaLog), 0.0)
        sumExitProbPParaLog += np.multiply(itemExitProbPParaLog, loggedItemExitProbPParaLog)
        
    MP = redPart - blackPart - sumPParaLog + sumExitProbPParaLog
    return (MP, nodeCommV)
    
#calculate the map equation values of node in each possible community
def CalMapEquationValues(commOrder, neighbourCommList, currentCommunityPartion, nodeNeighM, pPara, tauPara, coreNumToUse):
    #generate new partition
    partitionList = []
    for neighbourComm in neighbourCommList:
        newPartion = currentCommunityPartion.copy()
        newPartion[commOrder] = neighbourComm
        partitionList.append(newPartion)
    
    #set how many cores are used to calculate map equation, one core can calculate one kind of partition individually
    from multiprocessing import Pool
    pool = Pool(coreNumToUse)
    
    #map method can only take one argument, wrap the CalMapEquation with one fixed argument
    from functools import partial
    #prod_x has only one argument newPartion, others are fixed
    prod_x=partial(CalMapEquation, nodeNeighM=nodeNeighM, pPara=pPara, tauPara=tauPara)
    mpPartionList = pool.map(prod_x, partitionList)
    pool.close()
    pool.join()
    
    return mpPartionList

# a greedy assignment of nodes to communities, favoring local optimizations of modularity
def GreedyIMExtending(nodeCommV, nodeNeighM, pPara, tauPara, coreNumToUse, isRandom):
    
    #initialise the arguments
    roundCounter = 0
    
    #find the communities
    communtyIndexList = nodeNeighM.index.values
    
    #build the mapping between coarsen graph and original graph
    mappingList = CheckIMNodesOfComm(communtyIndexList, nodeCommV)
    coarsenCommPartitionV = nodeNeighM.index.values
    
    #record the current partition for comparing
    currentCommunityPartion = coarsenCommPartitionV
    
    #prepare the neighbours for each community
    commNeighbourList = CheckIMNeighbours(nodeNeighM, communtyIndexList)
    
    #limit the maximum steps of improvement
    while roundCounter < communtyIndexList.size:
        
        ##initial parameters
        #number of round has been done
        roundCounter += 1
        
        #check if there is improvement
        improvement = False
        
        #try to move nodes in given order 
        communties = np.arange(communtyIndexList.size)
        if isRandom:
            import random
            random.shuffle(communties)
        
        #try to merge two communities in given order
        for commOrder in communties:
            
            #check which communities the node can move to based on its neighbours, including the original one
            neighbourCommList = commNeighbourList[commOrder]
            
            #calculate the map equation values of node in each possible community
            mpPartionList = CalMapEquationValues(commOrder, neighbourCommList, currentCommunityPartion, nodeNeighM, pPara, tauPara, coreNumToUse)
            
            #find the partion with the smallest map equation value
            from operator import itemgetter
            #bestPartionPair is a tuple (map equation value, nodeCommV)
            bestPartionPair = min(mpPartionList, key=itemgetter(0))
            bestPartion = bestPartionPair[1]
            bestMP = bestPartionPair[0]
            
            #check if best partition is not the origignal one, the modularity is improved
            if not improvement and not np.array_equal(bestPartion, currentCommunityPartion):
                improvement = True
                
            #use the best partition for now to next round
            currentCommunityPartion = bestPartion
            currentMP = bestMP
                
        #if there is no improvement in this round, then stop. 
        if not improvement:
            break
    
    #update the new partition to the original graph
    for commOrder in xrange(currentCommunityPartion.size):
        affectedNodes = mappingList[commOrder]
        commIndex = currentCommunityPartion[commOrder]
        nodeCommV[affectedNodes] = commIndex
        
    return nodeCommV, currentMP

#building a new network whose nodes are now the communities found during the first phase
def CoarsenIMGraph(currentNodeCommV, currentNodeNeighM):
    #find communityIndex, use them as the new graph's node names
    communityIndexList = np.unique(currentNodeCommV)
    newNodeNum = communityIndexList.size
    iniExist = np.zeros((newNodeNum, newNodeNum),dtype=np.float)
    newNodeNeighM = pd.DataFrame(iniExist, index=communityIndexList, columns=communityIndexList)
    
    #find nodes in the community
    commNodeList = CheckIMNodesOfComm(communityIndexList, currentNodeCommV)
    
    #calculate the sum weight of edges from one community to another one
    for headIndex in xrange(newNodeNum):
        for tailIndex in xrange(newNodeNum):
            headSet = commNodeList[headIndex]
            tailSet = commNodeList[tailIndex]
            sumValue = currentNodeNeighM.iloc[tailSet, headSet].sum().sum()
            newNodeNeighM.iloc[tailIndex, headIndex] = sumValue
    
    return newNodeNeighM
    
#greedy search partition
def GreedyIMSearch(nodeNeighM, tauPara, coreNumToUse, isRandom):    
    level = 0
    
    #initial partition is each node is a community, which is the best one for now
    bestNodeCommV, nodeNeighM = GenIMNodeCommunityV(nodeNeighM)
    originalNodeNeighM = nodeNeighM.copy()
    
    #greedy improve the modularity
    while 1:
        level += 1
        print 'level %s' % level
        
        #calculate the ergodic node visit frequencies of nodes in the current graph
        pPara = CalErgodicNodeVisitFrequencies(nodeNeighM.copy(), tauPara)
        
        #reassign nodes to communities based on the modularity
        newNodeCommV, newMP = GreedyIMExtending(bestNodeCommV.copy(), nodeNeighM.copy(), pPara, tauPara, coreNumToUse, isRandom)
        
        #if the partition is firstly extended
        if level == 1:
            #check if there are changes, which means new partition generated
            if np.array_equal(bestNodeCommV, newNodeCommV):
                #there are no changes, return the initial partition
                return bestNodeCommV
            #new partition generated
            else:
                #there are changes, save as the best partition
                bestNodeCommV = newNodeCommV
                bestMP = newMP
                
                #if the best partition has two or less communities, terminate the extending
                if np.unique(bestNodeCommV).size <= 2:
                    return bestNodeCommV
                    
                #coarsen graph based on the new partitions in the originalNodeCommM
                nodeNeighM = CoarsenIMGraph(bestNodeCommV, originalNodeNeighM)
        else:
            #if there is no new partition generated, no improvement or new partition is worse, finish the searching
            if np.array_equal(bestNodeCommV, newNodeCommV) or newMP > bestMP:
                return bestNodeCommV
            #if the new partition is better, replace the current one and coarsen the graph
            else:
                #save the best partition for now
                bestNodeCommV = newNodeCommV
                bestMP = newMP
                
                #if the best partition has two or less communities, terminate the extending
                if np.unique(bestNodeCommV).size <= 2:
                    return bestNodeCommV
                    
                #coarsen graph based on the new partitions in the originalNodeCommM
                nodeNeighM = CoarsenIMGraph(bestNodeCommV, originalNodeNeighM)

def Infomap(adjM, coreNum, isRandom, tauPara):
    nodeNames = adjM.columns
    greedyNodeCommV = GreedyIMSearch(adjM, tauPara, coreNum, isRandom)
    rstDict = GetModularityClusters(greedyNodeCommV, nodeNames)
    return rstDict


#start to cluster the cell-to-cell communication network
def main(sourceFolder, signalType, weightType, edgeWeightKind, specificityThreshold, method, expandPara, inflatePara, maxLoop, coreNum, isRandom, tauPara):
    # load data
    adjMFileName = os.path.join(sourceFolder, 'Adjacency%sMatrix_speficificity_above_%s_signalType_%s_weightType_%s.xlsx' % (edgeWeightKind.title(), specificityThreshold,signalType,weightType))
    print adjMFileName
    if os.path.exists(adjMFileName):
        adjM = pd.read_excel(adjMFileName)
        if method == 'mcl':
            cltDict = MCL(adjM, expandPara, inflatePara, maxLoop)
        elif method == 'louvain':
            cltDict = Louvain(adjM, coreNum, isRandom)
        elif method == 'infomap':
            cltDict = Infomap(adjM, coreNum, isRandom, tauPara)
        print '%s identified %s cluster(s)' % (method, len(cltDict.keys()))
        cltFileName = os.path.join(sourceFolder, '%s_clusters_speficificity_above_%s_signalType_%s_weightType_%s_edgeWeightKind_%s.txt' % (method,specificityThreshold,signalType,weightType,edgeWeightKind))
        with open(cltFileName, 'w') as outfile:  
            json.dump(cltDict, outfile, sort_keys=True, indent=4)
        
    else:
        sys.exit("The adjacency matrix file does not exist.")
    
#===================================================================#
#-------------------------------START-------------------------------#
#===================================================================#

'''
python ClusterInteractionNetwork.py --sourceFolder M6C_data.frame --signalType paracrine --specificityThreshold 0.01 --method mcl
'''

if __name__ == '__main__':
    #process arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('--sourceFolder', required=True, help="the path to the dataset's folder")
    parser.add_argument('--specificityThreshold', type=float, default=0.01, help='do not draw the interactions whose speficificities are lower than the threshold (default 0.01).')
    parser.add_argument('--signalType', default='paracrine', help='paracrine (default) | juxtacrine')
    parser.add_argument('--weightType', default='mean', help="mean (default) | sum | es (enrichment score)")
    parser.add_argument('--edgeWeightKind', default='weight', help='weight (default) | count')
    parser.add_argument('--method', default='mcl', help='MCL (default) | Louvain | InfoMap')
    #MCL parameters
    parser.add_argument('--expandPara', type=int, default=2, help='MCL: expansion parameter, default is two')
    parser.add_argument('--inflatePara', type=int, default=2, help='MCL: power that the elements of the matrix are raised to in the inflation step, default is two')
    parser.add_argument('--maxLoop', type=int, default=20, help='MCL: maximum number of loops after which the process is terminated, default is 20')
    #Louvain or InfoMap parameters
    parser.add_argument('--coreNum', type=int, default=1, help='Modularity-based: number of the cores to use, default is one')
    parser.add_argument('--useRandomOrder', default='n', help='Modularity-based: n(o) (default) | y(es)')
    #InfoMap parameter
    parser.add_argument('--tauPara', type=float, default=0.15, help='InfoMap: probability of teleporting (default is 0.15)')
    
    opt = parser.parse_args()
    #print(opt)
    
    #check sourceFolder
    if not os.path.exists(opt.sourceFolder):
        sys.exit("The source dataset's folder does not exist.")
        
    #check signalType
    avaSignalTypeList = ['paracrine', 'juxtacrine']
    if opt.signalType.lower() not in avaSignalTypeList:
        sys.exit("The signalType can only be 'paracrine' or 'juxtacrine'.")
    else:
        signalType = opt.signalType.lower()
        
    #check weightType
    avaWeightTypeList = ['mean', 'sum', 'es']
    if opt.weightType.lower() not in avaWeightTypeList:
        sys.exit("The weightType can only be 'mean' or 'sum'.")
    else:
        weightType = opt.weightType.lower()
        
    #check edgeWeightKind
    avaEdgeWeightKindList = ['weight', 'count']
    if opt.edgeWeightKind.lower() not in avaEdgeWeightKindList:
        sys.exit("The edgeWeightKind can only be 'weight' or 'count'.")
    else:
        edgeWeightKind = opt.edgeWeightKind.lower()
        
    #check specificityThreshold
    if opt.specificityThreshold > 1.0 or opt.specificityThreshold < 0.0:
        specificityThreshold = 0.01
    else:
        specificityThreshold = opt.specificityThreshold
        
    #check clustering method
    avaMethodList = ['mcl', 'louvain', 'infomap']
    if opt.method not in avaMethodList:
        sys.exit("The given clustering method has not been implemented.")
        
    #check coreNum
    maxCoreNum = multiprocessing.cpu_count()
    if opt.coreNum > maxCoreNum:
        sys.exit("There are only %s cores availble, less than %s cores." % (maxCoreNum, opt.coreNum))
    
    #check if use random order
    useRandomOrderStr = opt.useRandomOrder.lower()
    if useRandomOrderStr[0] == 'y':
        isRandom = True
    elif useRandomOrderStr[0] == 'n':
        isRandom = False
    else:
        sys.exit("useRandomOrder can only be 'y' or 'n'.")
        
    #check tauPara
    if opt.tauPara > 1.0 or opt.tauPara < 0.0:
        tauPara = 0.15
    else:
        tauPara = opt.tauPara

    #pass argument check, show input data
    print '==================================================='
    print 'Input data:'
    print 'The source dataset folder: %s' % opt.sourceFolder
    print 'The cell-to-cell signaling type: %s' % signalType
    print 'The weight type of cell-to-cell signaling: %s' % weightType
    print 'The kind of cell-to-cell signaling edge weight: %s' % edgeWeightKind
    print 'The speficificity for interactions to draw is: %s' % specificityThreshold
    print 'The clustering method: %s' % opt.method
    if opt.method == 'mcl':
        print 'MCL: expansion parameter: %s' % opt.expandPara
        print 'MCL: inflation parameter: %s' % opt.inflatePara
        print 'MCL: maximum number of loops: %s' % opt.maxLoop
    elif opt.method == 'louvain':
        print 'Louvain: use random order: %s' % useRandomOrderStr[0]
        print 'Louvain: number of cores to use: %s' % opt.coreNum
    elif opt.method == 'infomap':
        print 'InfoMap: use random order: %s' % useRandomOrderStr[0]
        print 'InfoMap: number of cores to use: %s' % opt.coreNum
        print 'InfoMap: probability of teleporting: %s' % tauPara
    print '==================================================='
    
    #start to cluster the cell-to-cell communication network
    main(opt.sourceFolder, signalType, weightType, edgeWeightKind, specificityThreshold, opt.method, opt.expandPara, opt.inflatePara, opt.maxLoop, opt.coreNum, isRandom, tauPara)