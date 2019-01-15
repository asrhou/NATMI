#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 19 13:19:43 2018

@author: rhou
"""

import argparse
import pandas as pd
import numpy as np
import os, sys
import seaborn as sns 
import igraph as ig
import networkx as nx
from matplotlib import pyplot as plt
plt.switch_backend('agg')
from matplotlib.figure import SubplotParams as sppars
from matplotlib import colorbar as cb
    
# build adjacency matrix
def BuildAdjM(edgeDF, origlabels, labels, specificityThreshold, weightThreshold):
    # only keep edges of interest
    if specificityThreshold == 0:
        specificityThreshold = 1e-9
    if 'delta specificity' in edgeDF.columns:
        edgeDF = edgeDF.ix[(edgeDF['delta specificity']>=specificityThreshold)&(edgeDF['delta weight']>weightThreshold),]
    else:
        edgeDF = edgeDF.ix[(edgeDF['product of specified']>=specificityThreshold)&(edgeDF['product of original']>weightThreshold),]
    if specificityThreshold == 1e-9:
        specificityThreshold = 0
    
            
    adjM = pd.DataFrame(0.0, index=origlabels, columns=origlabels)
    adjCountM = pd.DataFrame(0, index=origlabels, columns=origlabels)
    for idx in edgeDF.index:
        if 'delta specificity' in edgeDF.columns:
            adjM.ix[edgeDF.ix[idx,'sending cluster name'], edgeDF.ix[idx,'target cluster name']] += edgeDF.ix[idx,'delta weight']
        else:
            adjM.ix[edgeDF.ix[idx,'sending cluster name'], edgeDF.ix[idx,'target cluster name']] += edgeDF.ix[idx,'product of original']
        adjCountM.ix[edgeDF.ix[idx,'sending cluster name'], edgeDF.ix[idx,'target cluster name']] += 1

    adjM.index = labels
    adjM.columns = labels
    adjCountM.index = labels
    adjCountM.columns = labels
    nxgW = nx.MultiDiGraph(adjM)
    nxgC = nx.MultiDiGraph(adjCountM)
    return edgeDF, adjM, adjCountM, nxgW, nxgC

#def BuildInterClusterNetwork(origlabels, labels, cltSizes, edgeDF, specificityThreshold, signalType, weightType, layout, plotFormat, plotWidth, plotHeight, minClusterSize, maxClusterSize, fontSize, edgeWidth, edgeColorMap, arrowSize, arrowHead, sourceFolder, dataType=''):
def BuildInterClusterNetwork(origlabels, labels, cltSizes, edgeDF, specificityThreshold, weightThreshold, signalType, weightType, layout, plotFormat, plotWidth, plotHeight, fontSize, edgeWidth, sourceFolder, dataType=''):
    
    #colors = sns.color_palette("Paired",len(labels)).as_hex()
    colors = sns.hls_palette(n_colors=len(labels)).as_hex()
    
    # build adjacency matrix
    edgeDF, adjM, adjCountM, nxgW, nxgC = BuildAdjM(edgeDF, origlabels, labels, specificityThreshold, weightThreshold)
    
    if dataType == '':
        edgeDFFileName = os.path.join(sourceFolder, 'Communications_speficificity_above_%s_weight_above_%s_signalType_%s_weightType_%s.xlsx' % (specificityThreshold,weightThreshold,signalType, weightType))
        adjMFileName = os.path.join(sourceFolder, 'AdjacencyWeightMatrix_speficificity_above_%s_weight_above_%s_signalType_%s_weightType_%s.xlsx' % (specificityThreshold,weightThreshold,signalType, weightType))
        adjCountMFileName = os.path.join(sourceFolder, 'AdjacencyCountMatrix_speficificity_above_%s_weight_above_%s_signalType_%s_weightType_%s.xlsx' % (specificityThreshold,weightThreshold,signalType, weightType))
    else:
        edgeDFFileName = os.path.join(sourceFolder, '%s_communications_speficificity_above_%s_weight_above_%s_signalType_%s_weightType_%s.xlsx' % (dataType, specificityThreshold,weightThreshold,signalType, weightType))
        adjMFileName = os.path.join(sourceFolder, '%s_adjacencyWeightMatrix_speficificity_above_%s_weight_above_%s_signalType_%s_weightType_%s.xlsx' % (dataType, specificityThreshold,weightThreshold,signalType, weightType))
        adjCountMFileName = os.path.join(sourceFolder, '%s_adjacencyCountMatrix_speficificity_above_%s_weight_above_%s_signalType_%s_weightType_%s.xlsx' % (dataType, specificityThreshold,weightThreshold,signalType, weightType))
    edgeDF.to_excel(edgeDFFileName, index=False)
    adjM.to_excel(adjMFileName)
    adjCountM.to_excel(adjCountMFileName)
    
    ## draw use graphviz
    #============total weight
    # convert to a graphviz graph
    nxgW = nx.nx_agraph.to_agraph(nxgW)
    # set graph size
    nxgW.graph_attr.update(size = "%s,%s" % (adjM.shape))
    # set graph font
    nxgW.graph_attr.update(fontname = "Arial")
    nxgW.graph_attr.update(fontsize = fontSize)
    
    # set node color
    nxgW.node_attr['style']='filled'
    idx = 0 
    for nd in nxgW.nodes():
        nd.attr['fillcolor'] = colors[idx]+'80'
        idx += 1
    maxVal = adjM.max().max()
    for ed in nxgW.edges():
        sn = ed[0]
        tn = ed[1]
        if adjM.ix[sn, tn]*edgeWidth/maxVal < 1:
            ed.attr['penwidth'] = 1
        else:
            ed.attr['penwidth'] = int(adjM.ix[sn, tn]*edgeWidth/maxVal)
    if dataType == '':
        plotFileName = 'c2ccnetwork_speficificity_above_%s_weight_above_%s_signalType_%s_weightType_%s_unit_total_weight_layout_%s.%s' % (specificityThreshold,weightThreshold,signalType,weightType,layout,plotFormat)
    else:
        plotFileName = '%s_network_speficificity_above_%s_weight_above_%s_signalType_%s_weightType_%s_unit_total_weight_layout_%s.%s' % (dataType, specificityThreshold,weightThreshold,signalType,weightType,layout,plotFormat)
    plotFileName = os.path.join(sourceFolder, plotFileName)
    nxgW.draw(plotFileName,prog=layout) 
    
    #============edge count
    nxgC = nx.nx_agraph.to_agraph(nxgC)
    nxgC.graph_attr.update(size = "%s,%s" % (adjM.shape))
    nxgC.graph_attr.update(fontname = "Arial")
    nxgC.graph_attr.update(fontsize = fontSize)
    nxgC.node_attr['style']='filled'
    idx = 0 
    for nd in nxgC.nodes():
        nd.attr['fillcolor']= colors[idx]+'80'
        idx += 1
        
    maxVal = adjCountM.max().max()
    for ed in nxgC.edges():
        sn = ed[0]
        tn = ed[1]
        if adjCountM.ix[sn, tn]*edgeWidth/maxVal < 1:
            ed.attr['penwidth'] = 1
        else:
            ed.attr['penwidth'] = int(adjCountM.ix[sn, tn]*edgeWidth/maxVal)
    if dataType == '':
        plotFileName = 'c2ccnetwork_speficificity_above_%s_weight_above_%s_signalType_%s_weightType_%s_unit_edge_count_layout_%s.%s' % (specificityThreshold,weightThreshold,signalType,weightType,layout,plotFormat)
    else:
        plotFileName = '%s_network_speficificity_above_%s_weight_above_%s_signalType_%s_weightType_%s_unit_edge_count_layout_%s.%s' % (dataType, specificityThreshold,weightThreshold,signalType,weightType,layout,plotFormat)
    plotFileName = os.path.join(sourceFolder, plotFileName)
    nxgC.draw(plotFileName,prog=layout) 
    
    # save edge list
    edgeDict = {'sending cluster name':[], 'target cluster name':[], 'weight':[], 'count':[]}
    for send in adjM.index:
        for target in adjM.columns:
            if adjM.ix[send, target] > 0:
                edgeDict['sending cluster name'].append(send.replace('\n',' '))
                edgeDict['target cluster name'].append(target.replace('\n',' '))
                edgeDict['weight'].append(adjM.ix[send, target])
                edgeDict['count'].append(adjCountM.ix[send, target])
    if dataType == '':
        edgeFileName = os.path.join(sourceFolder, 'c2ccnetwork_edges_speficificity_above_%s_weight_above_%s_signalType_%s_weightType_%s.csv' % (specificityThreshold,weightThreshold,signalType,weightType))
    else:
        edgeFileName = os.path.join(sourceFolder, '%s_network_edges_speficificity_above_%s_weight_above_%s_signalType_%s_weightType_%s.csv' % (dataType,specificityThreshold,weightThreshold,signalType,weightType))
    
    edgeDF = pd.DataFrame(edgeDict)
    edgeDF['weight per pair'] = edgeDF['weight']/edgeDF['count']
    edgeDF.sort_values(by="weight", ascending=False).to_csv(edgeFileName,columns=['sending cluster name', 'target cluster name', 'weight', 'count', 'weight per pair'],index=False)
    return edgeDF

# prepare color mapping for chord plots
def PrepareChordPlotColors(clt2CltDF):
    cltSet = set()
    cltSet = cltSet.union(set(clt2CltDF['sending cluster name']))
    cltSet = cltSet.union(set(clt2CltDF['target cluster name']))
    cltList = sorted(list(cltSet))
    cltNum = len(cltList)

    # construct color dict for R code
    colors = sns.hls_palette(n_colors=cltNum, h=.1, l=.3).as_hex()
    colorMapStr = 'grid.col = c('
    for i in range(len(cltList)):
        colorMapStr = colorMapStr + '"' + cltList[i] +'" = "' + colors[i].upper() +'", '
    colorMapStr = colorMapStr[:-2]+')\n'
    
    return colorMapStr

# compose R script to draw chordplots
def ComposeR(labels, clt2CltEdgeDF, colorMapStr, specificityThreshold, weightThreshold, signalType, weightType, sourceFolder, plotFormat, dataType):
    tempclusterList = sorted(list(set(clt2CltEdgeDF.ix[:, 'sending cluster name']).union(set(clt2CltEdgeDF.ix[:, 'target cluster name']))))
    
    if dataType == '':
        RCodeFileName = os.path.join(sourceFolder, 'chord_plot_signalType_%s_weightType_%s_specificityThreshold_%s_weightThreshold_%s.r' % (signalType, weightType, specificityThreshold, weightThreshold))
    else:
        RCodeFileName = os.path.join(sourceFolder, 'chord_plot_kind_%s_signalType_%s_weightType_%s_specificityThreshold_%s_weightThreshold_%s.r' % (dataType, signalType, weightType, specificityThreshold, weightThreshold))
    RCodeFileName = os.path.abspath(RCodeFileName)
    
    with open(RCodeFileName, "wb") as outfile:
        outfile.write('#！/usr/bin/R\n')
        outfile.write('\n')
        outfile.write('library(circlize)\n')
        outfile.write(colorMapStr)
        
        # draw different type of plots
        outfile.write('order = c(%s)\n' % (str(tempclusterList)[1:-1])) 
        
        outfile.write('#load adjacency list\n')
        if dataType == '':
            edgeFileName = os.path.join(sourceFolder, 'c2ccnetwork_edges_speficificity_above_%s_weight_above_%s_signalType_%s_weightType_%s.csv' % (specificityThreshold, weightThreshold,signalType,weightType))
        else:
            edgeFileName = os.path.join(sourceFolder, '%s_network_edges_speficificity_above_%s_weight_above_%s_signalType_%s_weightType_%s.csv' % (dataType,specificityThreshold, weightThreshold,signalType,weightType))
        edgeFileName = os.path.abspath(edgeFileName)
        outfile.write('ad <- read.csv("%s")\n' % (edgeFileName))
    
        outfile.write('#draw total weight graph\n')
        if dataType == '':
            plotfilename = os.path.join(sourceFolder,'chordplot_signalType_%s_weightType_%s_specificityThreshold_%s_weightThreshold_%s_total_weight.%s' % (signalType, weightType, specificityThreshold, weightThreshold, plotFormat))
        else:
            plotfilename = os.path.join(sourceFolder,'chordplot_kind_%s_signalType_%s_weightType_%s_specificityThreshold_%s_weightThreshold_%s_total_weight.%s' % (dataType, signalType, weightType, specificityThreshold, weightThreshold, plotFormat))
        plotfilename = os.path.abspath(plotfilename)
        outfile.write('%s("%s", width=%s, height=%s)\n' % (plotFormat, plotfilename, len(tempclusterList), len(tempclusterList)))
        outfile.write('chordDiagram(ad[,c(1,2,3)], order=order, directional = 1, direction.type = c("diffHeight", "arrows"), grid.col = grid.col)\n')
        outfile.write('dev.off()\n')
        outfile.write('circos.clear()\n')
        
        outfile.write('#draw edge count graph\n')
        if dataType == '':
            plotfilename = os.path.join(sourceFolder,'chordplot_signalType_%s_weightType_%s_specificityThreshold_%s_weightThreshold_%s_edge_count.%s' % (signalType, weightType, specificityThreshold, weightThreshold, plotFormat))
        else:
            plotfilename = os.path.join(sourceFolder,'chordplot_kind_%s_signalType_%s_weightType_%s_specificityThreshold_%s_weightThreshold_%s_edge_count.%s' % (dataType, signalType, weightType, specificityThreshold, weightThreshold, plotFormat))
        plotfilename = os.path.abspath(plotfilename)
        outfile.write('%s("%s", width=%s, height=%s)\n' % (plotFormat, plotfilename, len(tempclusterList), len(tempclusterList)))
        outfile.write('chordDiagram(ad[,c(1,2,4)], order=order, directional = 1, direction.type = c("diffHeight", "arrows"), grid.col = grid.col)\n')
        outfile.write('dev.off()\n')
        outfile.write('circos.clear()\n')
        
        outfile.write('#draw average edge weight graph\n')
        if dataType == '':
            plotfilename = os.path.join(sourceFolder,'chordplot_signalType_%s_weightType_%s_specificityThreshold_%s_weightThreshold_%s_mean_edge_weight.%s' % (signalType, weightType, specificityThreshold, weightThreshold, plotFormat))
        else:
            plotfilename = os.path.join(sourceFolder,'chordplot_kind_%s_signalType_%s_weightType_%s_specificityThreshold_%s_weightThreshold_%s_total_mean_edge_weight.%s' % (dataType, signalType, weightType, specificityThreshold, weightThreshold, plotFormat))
        plotfilename = os.path.abspath(plotfilename)
        outfile.write('%s("%s", width=%s, height=%s)\n' % (plotFormat, plotfilename, len(tempclusterList), len(tempclusterList)))
        outfile.write('chordDiagram(ad[,c(1,2,5)], order=order, directional = 1, direction.type = c("diffHeight", "arrows"), grid.col = grid.col)\n')
        outfile.write('dev.off()\n')
        outfile.write('circos.clear()\n')
        
    #draw chord graphs
    os.system('R CMD BATCH %s' % RCodeFileName)

# draw cluster-to-cluster communicaiton chord plot
def DrawChordPlots(labels, cltSizes, clt2CltEdgeDF, specificityThreshold, weightThreshold, signalType, weightType, sourceFolder, plotFormat, dataType=''):
    # prepare color mapping for chord plots
    colorMapStr = PrepareChordPlotColors(clt2CltEdgeDF)
    
    # compose R script to draw chordplots
    ComposeR(labels, clt2CltEdgeDF, colorMapStr, specificityThreshold, weightThreshold, signalType, weightType, sourceFolder, plotFormat, dataType)

#start to construct the cell-to-cell communication network
def MainClt(sourceFolder, signalType, weightType, specificityThreshold, weightThreshold, drawChordPlot, layout, plotFormat, plotWidth, plotHeight, fontSize, edgeWidth):
    # load data
    #examine the type of dataset
    clusterMapFilename = os.path.join(sourceFolder, 'clusterMapping.csv')
    if os.path.exists(clusterMapFilename):
        # process node properties for single dataset
        clusterMapDF = pd.read_csv(clusterMapFilename, index_col=None, header=0)
        clusterSizeDF = clusterMapDF.groupby('cluster').count()
        origlabels = [str(i) for i in clusterSizeDF.index]
        labels = [str(i)+'\n(%s cells)' % clusterSizeDF.ix[i,'cell'] for i in clusterSizeDF.index]
        cltSizes = [float(clusterSizeDF.ix[i,'cell']) for i in clusterSizeDF.index]
        # load edge properties
        edgeDF = pd.read_csv(os.path.join(sourceFolder, '%sClusteredEM_%s_signaling interactions.csv' % (weightType,signalType)), index_col=None, header=0)
        
        #cluster to cluster weighted directed network
        print 'plotting cluster to cluster weighted directed network'
        clt2CltEdgeDF = BuildInterClusterNetwork(origlabels, labels, cltSizes, edgeDF, specificityThreshold, weightThreshold, signalType, weightType, layout, plotFormat, plotWidth, plotHeight, fontSize, edgeWidth, sourceFolder)
    
        # draw cluster-to-cluster communicaiton chord plot
        if drawChordPlot:
            print 'plotting cluster to cluster chord plot'
            DrawChordPlots(labels, cltSizes, clt2CltEdgeDF, specificityThreshold, weightThreshold, weightThreshold, signalType, weightType, sourceFolder, plotFormat)
        
    else:
        # process node properties for delta dataset
        sumClusterDF = pd.read_csv(os.path.join(sourceFolder, 'cluster_size_comparison.csv'), index_col=None, header=0)
        sumClusterDF['cluster'] = sumClusterDF['cluster'].astype(str)
        origlabels = sorted([i for i in set(sumClusterDF['cluster'])])
        labels = []
        cltSizes = []
        for clt in origlabels:
            cltDF = sumClusterDF.ix[sumClusterDF['cluster'] == clt,].sort_values(by=['population (cells)'])
            if len(cltDF) == 2:
                deltaSize = int(cltDF.iloc[1,1]) - int(cltDF.iloc[0,1])
            else:
                deltaSize = int(cltDF.iloc[0,1])
            cltSizes.append(deltaSize)
            labels.append(clt+'\n(%s cells)' % deltaSize)
        
        # load edge properties
        kindList = ['appeared', 'disappeared', 'up_regulated', 'down_regulated']
        for kind in kindList:
            edgeDF = pd.read_csv(os.path.join(sourceFolder, '%s_signalType_%s_weightType_%s_signaling interactions.csv' % (kind,signalType,weightType)), index_col=None, header=0)

            #cluster to cluster weighted directed network
            print 'plotting cluster to cluster weighted directed network for %s interactions' % kind
            BuildInterClusterNetwork(origlabels, labels, cltSizes, edgeDF, specificityThreshold, weightThreshold, signalType, weightType, layout, plotFormat, plotWidth, plotHeight, fontSize, edgeWidth, sourceFolder, kind)
    
            # draw cluster-to-cluster communicaiton chord plot
            if drawChordPlot:
                print 'plotting cluster to cluster chord plot'
                DrawChordPlots(origlabels, labels, cltSizes, clt2CltEdgeDF, specificityThreshold, weightThreshold, signalType, weightType, sourceFolder, plotFormat, kind)
        
#single-LR-based cluster to cluster weighted directed network
def BuildSingleLRInterClusterNetwork(origlabels, labels, cltSizes, edgeDF, specificityThreshold, weightThreshold, ls, rs, signalType, weightType, layout, plotFormat, plotWidth, plotHeight, fontSize, edgeWidth, sourceFolder, dataType = ''):
    colors = sns.hls_palette(n_colors=len(labels)).as_hex()
    
    # build adjacency matrix
    edgeDF, adjM, adjCountM, nxgW, nxgC = BuildAdjM(edgeDF, origlabels, labels, specificityThreshold, weightThreshold)
    
    if dataType == '':
        edgeDFFileName = os.path.join(sourceFolder, 'Communications_speficificity_above_%s_weight_above_%s_signalType_%s_weightType_%s_%s-%s-based.xlsx' % (specificityThreshold, weightThreshold,signalType, weightType, ls, rs))
        adjMFileName = os.path.join(sourceFolder, 'AdjacencyWeightMatrix_speficificity_above_weight_above_%s_%s_signalType_%s_weightType_%s_%s-%s-based.xlsx' % (specificityThreshold, weightThreshold,signalType, weightType, ls, rs))
    else:
        edgeDFFileName = os.path.join(sourceFolder, '%s_communications_speficificity_above_%s_weight_above_%s_signalType_%s_weightType_%s_%s-%s-based.xlsx' % (dataType, specificityThreshold, weightThreshold,signalType, weightType, ls, rs))
        adjMFileName = os.path.join(sourceFolder, '%s_adjacencyWeightMatrix_speficificity_above_%s_weight_above_%s_signalType_%s_weightType_%s_%s-%s-based.xlsx' % (dataType, specificityThreshold, weightThreshold,signalType, weightType, ls, rs))
    edgeDF.to_excel(edgeDFFileName, index=False)
    adjM.to_excel(adjMFileName)
    
    ## draw use graphviz
    #============total weight
    # convert to a graphviz graph
    nxgW = nx.nx_agraph.to_agraph(nxgW)
    # set graph size
    nxgW.graph_attr.update(size = "%s,%s" % (adjM.shape))
    # set graph font
    nxgW.graph_attr.update(fontname = "Arial")
    nxgW.graph_attr.update(fontsize = fontSize)
    
    # set node color
    nxgW.node_attr['style']='filled'
    idx = 0 
    for nd in nxgW.nodes():
        nd.attr['fillcolor'] = colors[idx]+'80'
        idx += 1
    maxVal = adjM.max().max()
    for ed in nxgW.edges():
        sn = ed[0]
        tn = ed[1]
        if adjM.ix[sn, tn]*edgeWidth/maxVal < 1:
            ed.attr['penwidth'] = 1
        else:
            ed.attr['penwidth'] = int(adjM.ix[sn, tn]*edgeWidth/maxVal)
    
    if dataType == '':
        plotFileName = 'c2ccnetwork_speficificity_above_%s_weight_above_%s_signalType_%s_weightType_%s_unit_total_weight_layout_%s_%s-%s-based.%s' % (specificityThreshold,weightThreshold,signalType,weightType,layout, ls, rs, plotFormat)
    else:
        plotFileName = '%s_network_speficificity_above_%s_weight_above_%s_signalType_%s_weightType_%s_unit_total_weight_layout_%s_%s-%s-based.%s' % (dataType, specificityThreshold,weightThreshold,signalType,weightType,layout, ls, rs, plotFormat)
    plotFileName = os.path.join(sourceFolder, plotFileName)
    nxgW.draw(plotFileName,prog=layout) 
    
    # save edge list
    edgeDict = {'sending cluster name':[], 'target cluster name':[], 'weight':[], 'count':[]}
    for send in adjM.index:
        for target in adjM.columns:
            if adjM.ix[send, target] > 0:
                edgeDict['sending cluster name'].append(send.replace('\n',' '))
                edgeDict['target cluster name'].append(target.replace('\n',' '))
                edgeDict['weight'].append(adjM.ix[send, target])
                edgeDict['count'].append(adjCountM.ix[send, target])
    if dataType == '':
        edgeFileName = os.path.join(sourceFolder, 'c2ccnetwork_edges_speficificity_above_%s_weight_above_%s_signalType_%s_weightType_%s_%s-%s-based.csv' % (specificityThreshold, weightThreshold,signalType, weightType, ls, rs))
    else:
        edgeFileName = os.path.join(sourceFolder, '%s_network_edges_speficificity_above_%s_weight_above_%s_signalType_%s_weightType_%s_%s-%s-based.csv' % (dataType, specificityThreshold, weightThreshold,signalType, weightType, ls, rs))
    
    edgeDF = pd.DataFrame(edgeDict)
    edgeDF['weight per pair'] = edgeDF['weight']/edgeDF['count']
    edgeDF.sort_values(by="weight", ascending=False).to_csv(edgeFileName,columns=['sending cluster name', 'target cluster name', 'weight', 'count', 'weight per pair'],index=False)
    return edgeDF
    
#start to draw single-LR network
def SimpleLRNetwork(sourceFolder, signalType, weightType, specificityThreshold, weightThreshold, ls, rs, plotFormat, layout, plotWidth, plotHeight, fontSize, edgeWidth):
    # load edge info
    edgeDF = pd.read_csv(os.path.join(sourceFolder, '%sClusteredEM_%s_signaling interactions.csv' % (weightType,signalType)), index_col=None, header=0)
    
    # only keep edges of interest
    if specificityThreshold == 0:
            specificityThreshold = 1e-9
    if 'delta specificity' in edgeDF.columns:
        edgeDF = edgeDF.ix[(edgeDF['delta specificity']>=specificityThreshold)&(edgeDF['delta weight']>weightThreshold),]
    else:
        edgeDF = edgeDF.ix[(edgeDF['product of specified']>=specificityThreshold)&(edgeDF['product of original']>weightThreshold),]
    if specificityThreshold == 1e-9:
        specificityThreshold = 0
            
    edgeDF = edgeDF.ix[(edgeDF['ligand']==ls) & (edgeDF['receptor']==rs),]
    
    # process node properties for single dataset
    clusterMapFilename = os.path.join(sourceFolder, 'clusterMapping.csv')
    clusterMapDF = pd.read_csv(clusterMapFilename, index_col=None, header=0)
    clusterSizeDF = clusterMapDF.groupby('cluster').count()
    origlabels = [str(i) for i in clusterSizeDF.index]
    labels = [str(i)+'\n(%s cells)' % clusterSizeDF.ix[i,'cell'] for i in clusterSizeDF.index]
    cltSizes = [float(clusterSizeDF.ix[i,'cell']) for i in clusterSizeDF.index]
    labelDict = {}
    for i in clusterSizeDF.index:
        labelDict[str(i)] = str(i)+' (%s cells)' % clusterSizeDF.ix[i,'cell']
    
    #single-LR-based cluster to cluster weighted directed network
    BuildSingleLRInterClusterNetwork(origlabels, labels, cltSizes, edgeDF, specificityThreshold, weightThreshold, ls, rs, signalType, weightType, layout, plotFormat, plotWidth, plotHeight, fontSize, edgeWidth, sourceFolder)
    
# build data df for plotting
def BuildBubbleChartData(curDF, ligandList, receptorList, dataType):
    resultDict = {'ligand':[], 'original ligand':[], 'specified ligand':[],'receptor':[], 'original receptor':[], 'specified receptor':[], 'specificity':[], 'original weight':[]}
    for ly in range(1, len(ligandList)+1):
        for rx in range(1, len(receptorList)+1):
            temp = curDF.ix[(curDF['ligand']==ligandList[ly-1])&(curDF['receptor']==receptorList[rx-1]),]
            if len(temp)>0:
                resultDict['ligand'].append(ly)
                resultDict['receptor'].append(rx)
                if 'delta specificity' in curDF.columns:
                    resultDict['original ligand'].append(temp.ix[temp.index[0], 'delta ligand expression'])
                    resultDict['specified ligand'].append(temp.ix[temp.index[0], 'delta ligand specificity'])
                    resultDict['original receptor'].append(temp.ix[temp.index[0], 'delta receptor expression'])
                    resultDict['specified receptor'].append(temp.ix[temp.index[0], 'delta receptor specificity'])
                    resultDict['specificity'].append(temp.ix[temp.index[0], 'delta specificity'])
                    resultDict['original weight'].append(temp.ix[temp.index[0], 'delta weight'])
                else:
                    resultDict['original ligand'].append(temp.ix[temp.index[0], 'original ligand'])
                    resultDict['specified ligand'].append(temp.ix[temp.index[0], 'specified ligand'])
                    resultDict['original receptor'].append(temp.ix[temp.index[0], 'original receptor'])
                    resultDict['specified receptor'].append(temp.ix[temp.index[0], 'specified receptor'])
                    resultDict['specificity'].append(temp.ix[temp.index[0], 'product of specified'])
                    resultDict['original weight'].append(temp.ix[temp.index[0], 'product of original'])
    rsltDF = pd.DataFrame(resultDict)
    return rsltDF

#draw bubble chart for ligand-receptor edges from one cluster to another one
def DrawBubbleChart(curDF, sendClt, targetClt, plotFormat, edgeColorMap, sourceFolder, signalType, weightType, specificityThreshold, weightThreshold, dataType):
    # collect basic info
    ligandList = sorted(list(set(curDF['ligand'])))
    receptorList = sorted(list(set(curDF['receptor'])))
    
    # build data df for plotting
    bcdata = BuildBubbleChartData(curDF, ligandList, receptorList, dataType)
    
    # plot the bubble chart
    cmap = plt.get_cmap(edgeColorMap)
    sns.set_style("whitegrid")
    try:
        if max(bcdata.shape) < 25:
            fig = sns.relplot(x="receptor", y="ligand", hue="original weight", size="specificity", aspect=float(len(receptorList))/len(ligandList),palette=cmap, data=bcdata)
        else:
            fig = sns.relplot(x="receptor", y="ligand", hue="original weight", size="specificity", height=len(ligandList)/4.0, aspect=float(len(receptorList))/len(ligandList),palette=cmap, data=bcdata)
    except Exception as e:
        print e
        return
    fig.set_xlabels('Receptors in '+targetClt)
    fig.set_ylabels('Ligands in '+sendClt)
    fig.set(xticks=range(1, len(receptorList)+1), yticks=range(1,len(ligandList)+1))
    fig.set_xticklabels(receptorList,rotation=90)
    fig.set_yticklabels(ligandList)
    
    sns.despine(trim=True)
    if dataType == '':
        plotFileName = '%s_mediated_%s_to_%s_communication_specificityThreshold_%s_weightThreshold_%s_weightType_%s.%s' % (signalType, sendClt.split(' (')[0].replace('/','\\'), targetClt.split(' (')[0].replace('/','\\'), specificityThreshold, weightThreshold, weightType, plotFormat)
    else:
        plotFileName = '%s_%s_mediated_%s_to_%s_communication_specificityThreshold_%s_weightThreshold_%s_weightType_%s.%s' % (dataType, signalType, sendClt.split(' (')[0].replace('/','\\'), targetClt.split(' (')[0].replace('/','\\'), specificityThreshold, weightThreshold, weightType, plotFormat)
    fig.savefig(os.path.join(sourceFolder,plotFileName), bbox_inches='tight')
    plt.close()
                
    #save associated data
    for idx in bcdata.index:
        bcdata.ix[idx, 'ligand'] = ligandList[bcdata.ix[idx, 'ligand']-1]
        bcdata.ix[idx, 'receptor'] = receptorList[bcdata.ix[idx, 'receptor']-1]
    bcdata["combined weight"] = bcdata["original weight"] * bcdata["specificity"]
    bcdata = bcdata.sort_values(by="combined weight", ascending=False)
    if dataType == '':
        dataFilenameName = '%s_mediated_%s_to_%s_communication_specificityThreshold_%s_weightThreshold_%s_weightType_%s.xlsx' % (signalType, sendClt.split(' (')[0].replace('/','\\'), targetClt.split(' (')[0].replace('/','\\'), specificityThreshold, weightThreshold, weightType)
    else:
        dataFilenameName = '%s_%s_mediated_%s_to_%s_communication_specificityThreshold_%s_weightThreshold_%s_weightType_%s.xlsx' % (dataType, signalType, sendClt.split(' (')[0].replace('/','\\'), targetClt.split(' (')[0].replace('/','\\'), specificityThreshold, weightThreshold, weightType)
    bcdata.to_excel(os.path.join(sourceFolder,dataFilenameName), columns = ['ligand', 'original ligand', 'specified ligand', 'receptor', 'original receptor', 'specified receptor', 'specificity', 'original weight',"combined weight"], index=False, header=True)

def BuildInterClusterBubbleChart(origlabels, labelDict, edgeDF, weightType, specificityThreshold, weightThreshold, plotFormat, edgeColorMap, sourceFolder, signalType, dataType=''):
        
    # only keep edges of interest
    if specificityThreshold == 0:
            specificityThreshold = 1e-9
    if 'delta specificity' in edgeDF.columns:
        edgeDF = edgeDF.ix[(edgeDF['delta specificity']>=specificityThreshold)&(edgeDF['delta weight']>weightThreshold),]
    else:
        edgeDF = edgeDF.ix[(edgeDF['product of specified']>=specificityThreshold)&(edgeDF['product of original']>weightThreshold),]
    if specificityThreshold == 1e-9:
            specificityThreshold = 0
        
    for sendClt in origlabels:
        for targetClt in origlabels:
            curDF = edgeDF.ix[(edgeDF['sending cluster name']==sendClt)&(edgeDF['target cluster name']==targetClt),]
            if len(curDF) > 0:
                print 'plotting "%s" to "%s" ligand-receptor bubble chart' % (sendClt, targetClt)
                DrawBubbleChart(curDF, labelDict[sendClt], labelDict[targetClt], plotFormat, edgeColorMap, sourceFolder, signalType, weightType, specificityThreshold, specificityThreshold, dataType)
    
#start to construct bubble charts
def MainLR(sourceFolder, signalType, weightType, specificityThreshold, weightThreshold, edgeColorMap, plotFormat):
    # load data
    #examine the type of dataset
    clusterMapFilename = os.path.join(sourceFolder, 'clusterMapping.csv')
    if os.path.exists(clusterMapFilename):
        # process node properties for single dataset
        clusterMapDF = pd.read_csv(clusterMapFilename, index_col=None, header=0)
        clusterSizeDF = clusterMapDF.groupby('cluster').count()
        origlabels = [str(i) for i in clusterSizeDF.index]
        labelDict = {}
        for i in clusterSizeDF.index:
            labelDict[str(i)] = str(i)+' (%s cells)' % clusterSizeDF.ix[i,'cell']
        # load edge properties
        edgeDF = pd.read_csv(os.path.join(sourceFolder, '%sClusteredEM_%s_signaling interactions.csv' % (weightType,signalType)), index_col=None, header=0)
        
        #cluster to cluster bubble chats
        print 'plotting cluster to cluster bubble charts'
        BuildInterClusterBubbleChart(origlabels, labelDict, edgeDF, weightType, specificityThreshold, weightThreshold, plotFormat, edgeColorMap, sourceFolder, signalType)
    else:
        # process node properties for delta dataset
        sumClusterDF = pd.read_csv(os.path.join(sourceFolder, 'cluster_size_comparison.csv'), index_col=None, header=0)
        sumClusterDF['cluster'] = sumClusterDF['cluster'].astype(str)
        origlabels = sorted([i for i in set(sumClusterDF['cluster'])])
        labelDict = {}
        for clt in origlabels:
            cltDF = sumClusterDF.ix[sumClusterDF['cluster'] == clt,].sort_values(by=['population (cells)'])
            if len(cltDF) == 2:
                deltaSize = int(cltDF.iloc[1,1]) - int(cltDF.iloc[0,1])
            else:
                deltaSize = int(cltDF.iloc[0,1])
            labelDict[clt] = clt+' (%s cells)' % deltaSize
            
        kindList = ['appeared', 'disappeared', 'up_regulated', 'down_regulated']
        for kind in kindList:
            edgeDF = pd.read_csv(os.path.join(sourceFolder, '%s_signalType_%s_weightType_%s_signaling interactions.csv' % (kind,signalType,weightType)), index_col=None, header=0)
            
            #cluster to cluster bubble chats
            print 'plotting cluster to cluster bubble charts for %s interactions' % kind
            BuildInterClusterBubbleChart(origlabels, labelDict, edgeDF, weightType, specificityThreshold, weightThreshold, plotFormat, edgeColorMap, sourceFolder, signalType, kind)
    
#===================================================================#
#-------------------------------START-------------------------------#
#===================================================================#

if __name__ == '__main__':
    #process arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('--sourceFolder', required=True, help="the path to the dataset's folder")
    parser.add_argument('--weightType', default='mean', help="mean (default) | sum | es (enrichment score)")
    
    parser.add_argument('--layout', default='circo', help="circo (default) | neato | dot | twopi | fdp | nop")
    parser.add_argument('--plotWidth', type=int, default=12, help='plot width (default 12 inches).')
    parser.add_argument('--plotHeight', type=int, default=10, help='plot height (default 10 inches).')
    parser.add_argument('--plotFormat', default='pdf', help="pdf (default) | png | svg")
    
    parser.add_argument('--minClusterSize', type=int, default=10, help='min radius of the clusters (default 10).')
    parser.add_argument('--maxClusterSize', type=int, default=5000, help='max radius of the clusters (default 5000).')
    parser.add_argument('--fontSize', type=int, default=12, help='font size for node labels (default 8).')
    parser.add_argument('--edgeWidth', type=float, default=6, help='scale for edge width (default 1.0).')
    parser.add_argument('--arrowSize', type=float, default=25.0, help='scale for edge width (default 25.0).')
    parser.add_argument('--arrowHead', default='-|>', help='"-|>" (default), "-" ,"->" ,"-[" ,"-|>" ,"<-" ,"<->" ,"<|-" ,"<|-|>" ,"]-" ,"]-[" ,"|-|", "fancy" ,"simple" ,"wedge" ,full list: https://matplotlib.org/api/_as_gen/matplotlib.patches.ArrowStyle.html')
    parser.add_argument('--colorMap', default='brg', help="brg (default) | Greys | gray, full list: https://matplotlib.org/examples/color/colormaps_reference.html")
    parser.add_argument('--specificityThreshold', type=float, default=0.01, help='do not draw the interactions whose speficificities are lower than the threshold (default 0.01).')
    parser.add_argument('--weightThreshold', type=float, default=0.00, help='do not draw the interactions whose weight are lower than the threshold (default 0).')
    parser.add_argument('--signalType', default='paracrine', help='paracrine (default) | juxtacrine')
    parser.add_argument('--drawBubbleChart', default='n', help='n(o) (default) | y(es)')
    parser.add_argument('--drawChordPlot', default='n', help='n(o) (default) | y(es)')
    parser.add_argument('--drawSingleLRNetwork', nargs='*', help="ligand and receptor's symbols")
    
    opt = parser.parse_args()
    
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
        
    #check specificityThreshold
    if opt.specificityThreshold > 1.0 or opt.specificityThreshold < 0.0:
        specificityThreshold = 0.01
    else:
        specificityThreshold = opt.specificityThreshold
        
    #check weightThreshold
    if opt.weightThreshold < 0.0:
        weightThreshold = 0
    else:
        weightThreshold = opt.weightThreshold
        
    #check plotformat
    avaFormatList = ["pdf", "png", "svg"]
    if opt.plotFormat.lower() not in avaFormatList:
        sys.exit("The given plot format is not permitted.")
    else:
        plotFormat = opt.plotFormat.lower()
    
    #check edgeColorMap
    avaCmapList = ['viridis', 'plasma', 'inferno', 'magma', 'Greys', 'Purples', 'Blues', 'Greens', 'Oranges', 'Reds',
            'YlOrBr', 'YlOrRd', 'OrRd', 'PuRd', 'RdPu', 'BuPu',
            'GnBu', 'PuBu', 'YlGnBu', 'PuBuGn', 'BuGn', 'YlGn', 'binary', 'gist_yarg', 'gist_gray', 'gray', 'bone', 'pink',
            'spring', 'summer', 'autumn', 'winter', 'cool', 'Wistia',
            'hot', 'afmhot', 'gist_heat', 'copper', 'PiYG', 'PRGn', 'BrBG', 'PuOr', 'RdGy', 'RdBu',
            'RdYlBu', 'RdYlGn', 'Spectral', 'coolwarm', 'bwr', 'seismic','Pastel1', 'Pastel2', 'Paired', 'Accent',
            'Dark2', 'Set1', 'Set2', 'Set3',
            'tab10', 'tab20', 'tab20b', 'tab20c','flag', 'prism', 'ocean', 'gist_earth', 'terrain', 'gist_stern',
            'gnuplot', 'gnuplot2', 'CMRmap', 'cubehelix', 'brg', 'hsv',
            'gist_rainbow', 'rainbow', 'jet', 'nipy_spectral', 'gist_ncar']
    if opt.colorMap not in avaCmapList:
        sys.exit("The given color map name is not permitted.")
    
    #check if draw bubble chart
    drawBubbleChartStr = opt.drawBubbleChart.lower()
    if drawBubbleChartStr[0] == 'y':
        drawBubbleChart = True
    elif drawBubbleChartStr[0] == 'n':
        drawBubbleChart = False
    else:
        sys.exit("drawBubbleChart can only be 'y' or 'n'.")
        
    #check if draw chord plot
    drawChordPlotStr = opt.drawChordPlot.lower()
    if drawChordPlotStr[0] == 'y':
        drawChordPlot = True
    elif drawChordPlotStr[0] == 'n':
        drawChordPlot = False
    else:
        sys.exit("drawChordPlot can only be 'y' or 'n'.")
    
    if drawBubbleChart:
        #pass argument check, show input data
        print '==================================================='
        print 'Input data:'
        print 'The source dataset folder: %s' % opt.sourceFolder
        print 'The cell-to-cell signaling type: %s' % signalType
        print 'The weight type of cell-to-cell signaling: %s' % weightType
        print 'The speficificity threshold for interactions to draw is: %s' % specificityThreshold
        print 'The weight threshold for interactions to draw is: %s' % weightThreshold
        print 'The color map for edge weights: %s' % opt.colorMap
        print 'The plot format: %s' % plotFormat
        print '==================================================='
        
        #draw bubble charts for cluster to cluster ligand-receptor edges
        #start to construct bubble charts
        MainLR(opt.sourceFolder, signalType, weightType, specificityThreshold, weightThreshold, opt.colorMap, plotFormat)
        
    else:
        #draw cluster to cluster weighted directed network
        #check layout
        avaLayoutList = ['fr', 'drl', 'kk', 'lgl', 'random', 'rt', 'sphere', 'circle']
        avaLayoutList = ['neato', 'dot', 'twopi', 'circo', 'fdp', 'nop']
        if opt.layout.lower() not in avaLayoutList:
            sys.exit("The given layout is not permitted.")
        else:
            layout = opt.layout.lower()
        
        #check arrowHead
        avaArrowHeadList = ['-', '->', '-[', '-|>', '<-', '<->', '<|-', '<|-|>', ']-', ']-[', 'fancy', 'simple', 'wedge', '|-|']
        if opt.arrowHead not in avaArrowHeadList:
            sys.exit("The given arrow head is not permitted.")
        else:
            arrowHead = opt.arrowHead.replace('~','-')
                
        #check minClusterSize
        if opt.minClusterSize < 3:
            minClusterSize = 3
        else:
            minClusterSize = opt.minClusterSize
            
        #check maxClusterSize
        if opt.maxClusterSize < minClusterSize:
            maxClusterSize = minClusterSize
        else:
            maxClusterSize = opt.maxClusterSize
            
        #check fontSize
        if opt.fontSize > minClusterSize-2:
            fontSize = minClusterSize-2
        else:
            fontSize = opt.fontSize
            
        #check edgeWidth
        if opt.edgeWidth > minClusterSize/2:
            edgeWidth = minClusterSize/2
        else:
            edgeWidth = opt.edgeWidth
        #check  arrowSize
        if opt.arrowSize < 10:
            arrowSize = 10
        else:
            arrowSize = opt.arrowSize
            
        #check if draw single-LR network
        if opt.drawSingleLRNetwork is not None:#len(opt.drawSingleLRNetwork) > 0:
            
            if len(opt.drawSingleLRNetwork) != 2:
                sys.exit("drawSingleLRNetwork can only accept a ligand-receptor pair.")
                
            #check if symbols are legal
            edgeDF = pd.read_csv(os.path.join(opt.sourceFolder, '%sClusteredEM_%s_signaling interactions.csv' % (weightType,signalType)), index_col=None, header=0)
            edgeDF = edgeDF.ix[:,['ligand', 'receptor']]
            if opt.drawSingleLRNetwork[0] in list(edgeDF['ligand']) and opt.drawSingleLRNetwork[1] in list(edgeDF['receptor']):
                ls = opt.drawSingleLRNetwork[0]
                rs = opt.drawSingleLRNetwork[1]
            elif opt.drawSingleLRNetwork[1] in edgeDF['ligand'] and opt.drawSingleLRNetwork[0] in edgeDF['receptor']:
                ls = opt.drawSingleLRNetwork[1]
                rs = opt.drawSingleLRNetwork[0]
            else:
                sys.exit("Cannot find %s-%s pair or %s-%s pair  in the dataset." % (opt.drawSingleLRNetwork[0], opt.drawSingleLRNetwork[1], opt.drawSingleLRNetwork[1], opt.drawSingleLRNetwork[0]))
                
            #pass argument check, show input data
            print '==================================================='
            print 'Input data:'
            print 'The source dataset folder: %s' % opt.sourceFolder
            print 'The cell-to-cell signaling type: %s' % signalType
            print 'The weight type of cell-to-cell signaling: %s' % weightType
            print 'The speficificity threshold for interactions to draw is: %s' % specificityThreshold
            print 'The weight threshold for interactions to draw is: %s' % weightThreshold
            print 'The network layout: %s' % layout
            print 'The font size for cluster labels: %s' % fontSize
            print 'The edge width: %s' % edgeWidth
            print 'The plot width: %s inches' % opt.plotWidth
            print 'The plot height: %s inches' % opt.plotHeight
            print 'Ligand: %s' % ls
            print 'Receptor: %s' % rs
            print 'The plot format: %s' % plotFormat
            print '==================================================='
            
            #start to draw single-LR network
            SimpleLRNetwork(opt.sourceFolder, signalType, weightType, specificityThreshold, weightThreshold, ls, rs, plotFormat, layout, opt.plotWidth, opt.plotHeight, fontSize, edgeWidth)
            print 'single-LR network for %s-%s pair has been drawn' % (ls, rs)
                
        else:
            #pass argument check, show input data
            print '==================================================='
            print 'Input data:'
            print 'The source dataset folder: %s' % opt.sourceFolder
            print 'The cell-to-cell signaling type: %s' % signalType
            print 'The weight type of cell-to-cell signaling: %s' % weightType
            print 'The speficificity threshold for interactions to draw is: %s' % specificityThreshold
            print 'The weight threshold for interactions to draw is: %s' % weightThreshold
            print 'The network layout: %s' % layout
            print 'The font size for cluster labels: %s' % fontSize
            print 'The edge width: %s' % edgeWidth
            print 'The plot width: %s inches' % opt.plotWidth
            print 'The plot height: %s inches' % opt.plotHeight
            print 'The plot format: %s' % plotFormat
            print 'Draw chord plots: %s' % drawChordPlotStr[0]
            print '==================================================='
            
            #start to construct the cell-to-cell communication network
            MainClt(opt.sourceFolder, signalType, weightType, specificityThreshold, weightThreshold, drawChordPlot, layout, plotFormat, opt.plotWidth, opt.plotHeight, fontSize, edgeWidth)
            
    