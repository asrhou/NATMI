#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 21 12:00:51 2019

@author: rhou
"""

import warnings
warnings.filterwarnings("ignore")
import pandas as pd
import numpy as np
import argparse, os, sys
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns

#build the folder to save the analysis results
def BuildSaveFolder(refFolder, targetFolder):
    ref = os.path.basename(refFolder)
    target = os.path.basename(targetFolder)
    resultFolder = os.path.join(os.path.dirname(refFolder), "delta_ref_%s_tgt_%s" % (ref, target))
    if not os.path.exists(resultFolder):
        os.mkdir(resultFolder)
    rmf = os.path.join(resultFolder,'README.txt')
    with open(rmf, 'w') as file_object:
        file_object.write('README\n')
        file_object.write('\n')
        file_object.write('Reference dataset: %s\n' % ref)
        file_object.write('Target dataset: %s\n' % target)
        file_object.write('\n')
    return resultFolder, rmf

# find population changes in clusters
def IdentifyPopulationChanges(refClusterMapDF, testClusterMapDF, resultFolder, rmf):
    refClusterSizeDF = refClusterMapDF.groupby('cluster').count()
    refClusterSizeDF['fraction'] = 100.0 * refClusterSizeDF['cell']/refClusterSizeDF['cell'].sum()
    refClusterSizeDF['class'] = ['reference dataset'] * len(refClusterSizeDF)
    refClusterSizeDF = refClusterSizeDF.reset_index()
    testClusterSizeDF = testClusterMapDF.groupby('cluster').count()
    testClusterSizeDF['fraction'] = 100.0 * testClusterSizeDF['cell']/testClusterSizeDF['cell'].sum()
    testClusterSizeDF['class'] = ['target dataset'] * len(testClusterSizeDF)
    testClusterSizeDF = testClusterSizeDF.reset_index()
    
    sumClusterDF = pd.concat([refClusterSizeDF, testClusterSizeDF])
    sumClusterDF.columns = ['Cluster', 'Population (cells)', 'Fraction (%)', 'Source']
    sumClusterDF = sumClusterDF.sort_values(by=['Cluster'])
    sumClusterDF.set_index('Cluster').to_excel(os.path.join(resultFolder,'cluster_comparison.xlsx'), index=True, header=True)
    
    # find all clusters
    allClusters = set(sumClusterDF['Cluster'])
    plotHeight = len(allClusters) * 2 / 5
    plotWidth = 2 * plotHeight
    
    # sort clusters based on the fold changes
    cltOrder = []
    for clt in allClusters:
        cltDF = sumClusterDF.ix[sumClusterDF['Cluster'] == clt,].sort_values(by=['Source'])
        if len(cltDF) == 2:
            cltOrder.append((clt,float(cltDF.iloc[1,2])/cltDF.iloc[0,2]))
        else:
            cltOrder.append((clt,0))
    cltOrder = sorted(cltOrder, key=lambda tup: tup[1], reverse=True)
    cltOrder = [i[0] for i in cltOrder]
    
    # visualise population changes in clusters
    sns.set_style("whitegrid")
    fig = plt.figure(figsize=(plotWidth, plotHeight))
    ax = sns.barplot(y="Cluster", x="Population (cells)", hue="Source", data=sumClusterDF, order=cltOrder, palette={'reference dataset':'r','target dataset':'lime'})
    fig = ax.get_figure()
    fig.savefig(os.path.join(resultFolder,'cluster_size_comparison.pdf'), bbox_inches='tight')
    fig = plt.figure(figsize=(plotWidth, plotHeight))
    ax = sns.barplot(y="Cluster", x="Fraction (%)", hue="Source", data=sumClusterDF, order=cltOrder, palette={'reference dataset':'r','target dataset':'lime'})
    fig = ax.get_figure()
    fig.savefig(os.path.join(resultFolder,'cluster_fraction_comparison.pdf'), bbox_inches='tight')
    
    with open(rmf, 'a') as file_object:
        file_object.write('cluster_comparison.xlsx: cluster size and fraction of each cluster in two datasets\n')
        file_object.write('cluster_size_comparison.pdf: comparison of the cluster sizes of each cluster in two datasets\n')
        file_object.write('cluster_fraction_comparison.pdf: comparison of the cluster fractions of each cluster in two datasets\n')
        file_object.write('\n')
        
    return sumClusterDF

# calculate changes in the cell-to-cell communications 
def IdentifyLREdgeChanges(sizeClusterDF, refEdgeDF, testEdgeDF, signalType, weightType, resultFolder, rmf):
    resultEdgeF = os.path.join(resultFolder, 'Delta_edges_%s' % (signalType))
    if not os.path.exists(resultEdgeF):
        os.mkdir(resultEdgeF)
    with open(rmf, 'a') as file_object:
        file_object.write('Delta_edges_xxx folder: variations in the edges of two datasets using the certain ligand-receptor pair list\n')
        
    mergedDF = pd.merge(refEdgeDF, testEdgeDF, how='outer', on=['sending cluster name', 'target cluster name', 'ligand', 'receptor']).fillna(0)
    oldColNames = ['sending cluster name', 'ligand', 'receptor', 'target cluster name', 'count ligand_x', 'frequency ligand_x', 'original ligand_x', 'specified ligand_x', 'count receptor_x', 'frequency receptor_x', 'original receptor_x', 'specified receptor_x', 'product of original_x', 'product of specified_x', 'count ligand_y', 'frequency ligand_y', 'original ligand_y', 'specified ligand_y', 'count receptor_y', 'frequency receptor_y', 'original receptor_y', 'specified receptor_y', 'product of original_y', 'product of specified_y']
    deltaColNames = ['Sending cluster', 'Ligand symbol', 'Receptor symbol', 'Target cluster', 'Delta ligand expressing cells', 'Delta ligand detection rate', 'Delta ligand expression', 'Delta ligand specificity', 'Delta receptor expressing cells', 'Delta receptor detection rate', 'Delta receptor expression', 'Delta receptor specificity', 'Delta edge expression weight', 'Delta edge specificity weight']

    disaperEdgeDF = mergedDF.ix[(mergedDF['product of original_y']==0)&(mergedDF['product of original_x']>0),['sending cluster name', 'ligand', 'receptor', 'target cluster name', 'count ligand_x', 'frequency ligand_x', 'original ligand_x', 'specified ligand_x', 'count receptor_x', 'frequency receptor_x', 'original receptor_x', 'specified receptor_x', 'product of original_x', 'product of specified_x']]
    disaperEdgeDF.columns = deltaColNames
    disaperEdgeDF = disaperEdgeDF.sort_values(by='Delta edge expression weight', ascending=False)
    disaperEdgeDF.to_csv(os.path.join(resultEdgeF, 'Disappeared_%s.csv' % (weightType)), index=False, header=True, columns=deltaColNames)
    print('%s edges disappeared' % "{:,}".format(len(disaperEdgeDF)))
    with open(rmf, 'a') as file_object:
        file_object.write('\t|->Disappeared_xxx.csv: edges that only detected in the reference dataset\n')
        
    aperEdgeDF = mergedDF.ix[(mergedDF['product of original_x']==0)&(mergedDF['product of original_y']>0),['sending cluster name', 'ligand', 'receptor', 'target cluster name', 'count ligand_y', 'frequency ligand_y', 'original ligand_y', 'specified ligand_y', 'count receptor_y', 'frequency receptor_y', 'original receptor_y', 'specified receptor_y', 'product of original_y', 'product of specified_y']]
    aperEdgeDF.columns = deltaColNames
    aperEdgeDF = aperEdgeDF.sort_values(by='Delta edge expression weight', ascending=False)
    aperEdgeDF.to_csv(os.path.join(resultEdgeF, 'Appeared_%s.csv' % (weightType)), index=False, header=True, columns=deltaColNames)
    print('%s edges appeared' % "{:,}".format(len(aperEdgeDF)))
    with open(rmf, 'a') as file_object:
        file_object.write('\t|->Appeared_xxx.csv: edges that only detected in the target dataset\n')
    
    oldColNames = ['sending cluster name', 'ligand', 'receptor', 'target cluster name', 'count ligand_x', 'frequency ligand_x', 'original ligand_x', 'specified ligand_x', 'count receptor_x', 'frequency receptor_x', 'original receptor_x', 'specified receptor_x', 'product of original_x', 'product of specified_x', 'count ligand_y', 'frequency ligand_y', 'original ligand_y', 'specified ligand_y', 'count receptor_y', 'frequency receptor_y', 'original receptor_y', 'specified receptor_y', 'product of original_y', 'product of specified_y']
    changedColNames = ['sending cluster name', 'ligand', 'receptor', 'target cluster name', 'original ligand count', 'original ligand frequency', 'original ligand expression', 'original ligand specificity', 'original receptor count', 'original receptor frequency', 'original receptor expression', 'original receptor specificity', 'original weight', 'original specificity', 'changed ligand count', 'changed ligand frequency', 'changed ligand expression', 'changed ligand specificity', 'changed receptor count', 'changed receptor frequency', 'changed receptor expression', 'changed receptor specificity', 'changed weight', 'changed specificity', 'delta ligand expression', 'fold change of ligand expression', 'log2 fold change of ligand expression', 'delta ligand specificity', 'fold change of ligand specificity', 'log2 fold change of ligand specificity', 'delta receptor expression', 'fold change of receptor expression', 'log2 fold change of receptor expression', 'delta receptor specificity', 'fold change of receptor specificity', 'log2 fold change of receptor specificity', 'delta edge expression weight', 'delta edge specificity weight', 'fold change of weight', 'log2 fold change of weight', 'fold change of specificity', 'log2 fold change of specificity']
    newdeltaColNames = ['Sending cluster', 'Ligand symbol', 'Receptor symbol', 'Target cluster', 'Ligand expressing cells in condition 1', 'Ligand detection rate in condition 1', 'Ligand expression in condition 1', 'Ligand specificity in condition 1', 'Receptor expressing cells in condition 1', 'Receptor detection rate in condition 1', 'Receptor expression in condition 1', 'Receptor specificity in condition 1', 'Edge expression weight in condition 1', 'Edge specificity weight in condition 1', 'Ligand expressing cells in condition 2', 'Ligand detection rate in condition 2', 'Ligand expression in condition 2', 'Ligand specificity in condition 2', 'Receptor expressing cells in condition 2', 'Receptor detection rate in condition 2', 'Receptor expression in condition 2', 'Receptor specificity in condition 2', 'Edge expression weight in condition 2', 'Edge specificity weight in condition 2', 'Delta ligand expression', 'Fold change of ligand expression', 'Log2-transformed fold change of ligand expression', 'Delta ligand specificity', 'Fold change of ligand specificity', 'Log2-transformed fold change of ligand specificity', 'Delta receptor expression', 'Fold change of receptor expression', 'Log2-transformed fold change of receptor expression', 'Delta receptor specificity', 'Fold change of receptor specificity', 'Log2-transformed fold change of receptor specificity', 'Delta edge expression weight', 'Delta edge specificity weight', 'Fold change of edge expression weight', 'Log2-transformed fold change of edge expression weight', 'Fold change of edge specificity weight', 'Log2-transformed fold change of edge specificity weight']

    upEdgeDF = mergedDF.ix[(mergedDF['product of original_x'] < mergedDF['product of original_y'])&(mergedDF['product of original_x']>0),oldColNames]
    upEdgeDF.columns = changedColNames[:-18]
    upEdgeDF['delta ligand expression'] = upEdgeDF['changed ligand expression'] - upEdgeDF['original ligand expression']
    upEdgeDF['fold change of ligand expression'] = upEdgeDF['original ligand expression'] / upEdgeDF['changed ligand expression']
    upEdgeDF['log2 fold change of ligand expression'] = np.log2(upEdgeDF['fold change of ligand expression'])
    upEdgeDF['delta ligand specificity'] = upEdgeDF['changed ligand specificity'] - upEdgeDF['original ligand specificity']
    upEdgeDF['fold change of ligand specificity'] = upEdgeDF['original ligand specificity'] / upEdgeDF['changed ligand specificity']
    upEdgeDF['log2 fold change of ligand specificity'] = np.log2(upEdgeDF['fold change of ligand specificity'])
    upEdgeDF['delta receptor expression'] = upEdgeDF['changed receptor expression'] - upEdgeDF['original receptor expression']
    upEdgeDF['fold change of receptor expression'] = upEdgeDF['original receptor expression'] / upEdgeDF['changed receptor expression']
    upEdgeDF['log2 fold change of receptor expression'] = np.log2(upEdgeDF['fold change of receptor expression'])
    upEdgeDF['delta receptor specificity'] = upEdgeDF['changed receptor specificity'] - upEdgeDF['original receptor specificity']
    upEdgeDF['fold change of receptor specificity'] = upEdgeDF['original receptor specificity'] / upEdgeDF['changed receptor specificity']
    upEdgeDF['log2 fold change of receptor specificity'] = np.log2(upEdgeDF['fold change of receptor specificity'])
    upEdgeDF['delta edge expression weight'] = upEdgeDF['changed weight'] - upEdgeDF['original weight']
    upEdgeDF['delta edge specificity weight'] = upEdgeDF['changed specificity'] - upEdgeDF['original specificity']
    upEdgeDF['fold change of weight'] = upEdgeDF['changed weight'] / upEdgeDF['original weight']
    upEdgeDF['log2 fold change of weight'] = np.log2(upEdgeDF['fold change of weight'])
    upEdgeDF['fold change of specificity'] = upEdgeDF['changed specificity'] / upEdgeDF['original specificity']
    upEdgeDF['log2 fold change of specificity'] = np.log2(upEdgeDF['fold change of specificity'])
    upEdgeDF = upEdgeDF.sort_values(by="fold change of weight", ascending=False)
    upEdgeDF = upEdgeDF.ix[:,changedColNames]
    upEdgeDF.columns = newdeltaColNames
    upEdgeDF.to_csv(os.path.join(resultEdgeF, 'UP-regulated_%s.csv' % (weightType)), index=False, header=True, columns=newdeltaColNames)
    print('%s edges are up-regulated' % "{:,}".format(len(upEdgeDF)))
    with open(rmf, 'a') as file_object:
        file_object.write('\t|->UP-regulated_xxx.csv: edges that have higher expression weights in the target dataset\n')
    
    downEdgeDF = mergedDF.ix[(mergedDF['product of original_x'] > mergedDF['product of original_y'])&(mergedDF['product of original_y']>0),oldColNames]
    downEdgeDF.columns = changedColNames[:-18]
    downEdgeDF['delta ligand expression'] = downEdgeDF['original ligand expression'] - downEdgeDF['changed ligand expression']
    downEdgeDF['fold change of ligand expression'] = downEdgeDF['original ligand expression'] / downEdgeDF['changed ligand expression']
    downEdgeDF['log2 fold change of ligand expression'] = np.log2(downEdgeDF['fold change of ligand expression'])
    downEdgeDF['delta ligand specificity'] = downEdgeDF['original ligand specificity'] - downEdgeDF['changed ligand specificity']
    downEdgeDF['fold change of ligand specificity'] = downEdgeDF['original ligand specificity'] / downEdgeDF['changed ligand specificity']
    downEdgeDF['log2 fold change of ligand specificity'] = np.log2(downEdgeDF['fold change of ligand specificity'])
    downEdgeDF['delta receptor expression'] = downEdgeDF['original receptor expression'] - downEdgeDF['changed receptor expression']
    downEdgeDF['fold change of receptor expression'] = downEdgeDF['original receptor expression'] / downEdgeDF['changed receptor expression']
    downEdgeDF['log2 fold change of receptor expression'] = np.log2(downEdgeDF['fold change of receptor expression'])
    downEdgeDF['delta receptor specificity'] = downEdgeDF['original receptor specificity'] - downEdgeDF['changed receptor specificity'] 
    downEdgeDF['fold change of receptor specificity'] = downEdgeDF['original receptor specificity'] / downEdgeDF['changed receptor specificity']
    downEdgeDF['log2 fold change of receptor specificity'] = np.log2(downEdgeDF['fold change of receptor specificity'])
    downEdgeDF['delta edge expression weight'] = downEdgeDF['original weight'] - downEdgeDF['changed weight']
    downEdgeDF['delta edge specificity weight'] = downEdgeDF['original specificity'] - downEdgeDF['changed specificity']
    downEdgeDF['fold change of weight'] = downEdgeDF['original weight'] / downEdgeDF['changed weight']
    downEdgeDF['log2 fold change of weight'] = np.log2(downEdgeDF['fold change of weight'])
    downEdgeDF['fold change of specificity'] = downEdgeDF['original specificity'] / downEdgeDF['changed specificity']
    downEdgeDF['log2 fold change of specificity'] = np.log2(downEdgeDF['fold change of specificity'])
    downEdgeDF = downEdgeDF.sort_values(by="fold change of weight", ascending=False)
    downEdgeDF = downEdgeDF.ix[:,changedColNames]
    downEdgeDF.columns = newdeltaColNames
    downEdgeDF.to_csv(os.path.join(resultEdgeF, 'DOWN-regulated_%s.csv' % (weightType)), index=False, header=True, columns=newdeltaColNames)
    print('%s edges are down-regulated' % "{:,}".format(len(downEdgeDF)))
    with open(rmf, 'a') as file_object:
        file_object.write('\t|->DOWN-regulated_xxx.csv: edges that have lower expression weights in the target dataset\n')
    
    eqEdgeDF = mergedDF.ix[(mergedDF['product of original_x'] == mergedDF['product of original_y'])&(mergedDF['product of original_x']>0)&(mergedDF['product of original_y']>0),oldColNames]
    eqEdgeDF.columns = changedColNames[:-18]
    eqEdgeDF['delta ligand expression'] = eqEdgeDF['original ligand expression'] - eqEdgeDF['changed ligand expression']
    eqEdgeDF['fold change of ligand expression'] = eqEdgeDF['original ligand expression'] / eqEdgeDF['changed ligand expression']
    eqEdgeDF['log2 fold change of ligand expression'] = np.log2(eqEdgeDF['fold change of ligand expression'])
    eqEdgeDF['delta ligand specificity'] = eqEdgeDF['original ligand specificity'] - eqEdgeDF['changed ligand specificity']
    eqEdgeDF['fold change of ligand specificity'] = eqEdgeDF['original ligand specificity'] / eqEdgeDF['changed ligand specificity']
    eqEdgeDF['log2 fold change of ligand specificity'] = np.log2(eqEdgeDF['fold change of ligand specificity'])
    eqEdgeDF['delta receptor expression'] = eqEdgeDF['original receptor expression'] - eqEdgeDF['changed receptor expression']
    eqEdgeDF['fold change of receptor expression'] = eqEdgeDF['original receptor expression'] / eqEdgeDF['changed receptor expression']
    eqEdgeDF['log2 fold change of receptor expression'] = np.log2(eqEdgeDF['fold change of receptor expression'])
    eqEdgeDF['delta receptor specificity'] = eqEdgeDF['original receptor specificity'] - eqEdgeDF['changed receptor specificity'] 
    eqEdgeDF['fold change of receptor specificity'] = eqEdgeDF['original receptor specificity'] / eqEdgeDF['changed receptor specificity']
    eqEdgeDF['log2 fold change of receptor specificity'] = np.log2(eqEdgeDF['fold change of receptor specificity'])
    eqEdgeDF['delta edge expression weight'] = eqEdgeDF['original weight'] - eqEdgeDF['changed weight']
    eqEdgeDF['delta edge specificity weight'] = eqEdgeDF['original specificity'] - eqEdgeDF['changed specificity']
    eqEdgeDF['fold change of weight'] = eqEdgeDF['original weight'] / eqEdgeDF['changed weight']
    eqEdgeDF['log2 fold change of weight'] = np.log2(eqEdgeDF['fold change of weight'])
    eqEdgeDF['fold change of specificity'] = eqEdgeDF['original specificity'] / eqEdgeDF['changed specificity']
    eqEdgeDF['log2 fold change of specificity'] = np.log2( eqEdgeDF['fold change of specificity'])
    eqEdgeDF = eqEdgeDF.sort_values(by="fold change of specificity", ascending=False)
    eqEdgeDF = eqEdgeDF.ix[:,changedColNames]
    eqEdgeDF.columns = newdeltaColNames
    eqEdgeDF.to_csv(os.path.join(resultEdgeF, 'Stable_%s.csv' % (weightType)), index=False, header=True, columns=newdeltaColNames)
    print('%s edges are not changed' % "{:,}".format(len(eqEdgeDF)))
    with open(rmf, 'a') as file_object:
        file_object.write('\t|->Stable_xxx.csv: edges that have identical expression weights in two datasets\n')
        
    mergedDF.columns = changedColNames[:-18]
    mergedDF['delta ligand expression'] = mergedDF['original ligand expression'] - mergedDF['changed ligand expression']
    mergedDF['delta ligand specificity'] = mergedDF['original ligand specificity'] - mergedDF['changed ligand specificity']
    mergedDF['fold change of ligand expression'] = mergedDF['original ligand expression'] / mergedDF['changed ligand expression']
    mergedDF['log2 fold change of ligand expression'] = np.log2(mergedDF['fold change of ligand expression'])
    mergedDF['fold change of ligand specificity'] = mergedDF['original ligand specificity'] / mergedDF['changed ligand specificity']
    mergedDF['log2 fold change of ligand specificity'] = np.log2(mergedDF['fold change of ligand specificity'])
    mergedDF['delta receptor expression'] = mergedDF['original receptor expression'] - mergedDF['changed receptor expression']
    mergedDF['delta receptor specificity'] = mergedDF['original receptor specificity'] - mergedDF['changed receptor specificity'] 
    mergedDF['fold change of receptor expression'] = mergedDF['original receptor expression'] / mergedDF['changed receptor expression']
    mergedDF['log2 fold change of receptor expression'] = np.log2(mergedDF['fold change of receptor expression'])
    mergedDF['fold change of receptor specificity'] = mergedDF['original receptor specificity'] / mergedDF['changed receptor specificity']
    mergedDF['log2 fold change of receptor specificity'] = np.log2(mergedDF['fold change of receptor specificity'])
    mergedDF['delta edge expression weight'] = mergedDF['original weight'] - mergedDF['changed weight']
    mergedDF['delta edge specificity weight'] = mergedDF['original specificity'] - mergedDF['changed specificity']
    mergedDF['fold change of weight'] = mergedDF['original weight'] / mergedDF['changed weight']
    mergedDF['log2 fold change of weight'] = np.log2(mergedDF['fold change of weight'])
    mergedDF['fold change of specificity'] = mergedDF['original specificity'] / mergedDF['changed specificity']
    mergedDF['log2 fold change of specificity'] = np.log2(mergedDF['fold change of specificity'])
    mergedDF = mergedDF.sort_values(by="fold change of specificity", ascending=False)
    mergedDF = mergedDF.ix[:,changedColNames]
    mergedDF.columns = newdeltaColNames
    mergedDF.to_csv(os.path.join(resultEdgeF, 'All_edges_%s.csv' % (weightType)), index=False, header=True, columns=newdeltaColNames)
    print('%s edges are detected in at least one dataset' % "{:,}".format(len(mergedDF)))
    
    with open(rmf, 'a') as file_object:
        file_object.write('\t|->All_edges_xxx.csv: all edges that detected in both dataset\n')
        file_object.write('\n')
        
        file_object.write('Sending cluster: clusters that expresses the ligand\n')
        file_object.write('Target cluster: clusters that expresses the receptor\n')
        file_object.write('Ligand/receptor symbol: official gene symbol of the detected ligand/receptor\n')
        file_object.write('Delta ligand/receptor expression: change of ligand/receptor expression values in two datasets\n')
        file_object.write('Delta ligand/receptor specificity: change of ligand/receptor expression value derived specificities in two datasets\n')
        file_object.write('Delta edge expression weight: change of the edge expression weights in two datasets\n')
        file_object.write('Delta edge specificity weight: change of the edge specificity weights in two datasets\n')
        file_object.write('Log2-transformed fold change of edge expression weight: ratio of the edge expression weight in the target dataset to that in the reference dataset\n')
        file_object.write('Log2-transformed fold change of edge specificity weight: ratio of the edge specificity weight in the target dataset to that in the reference dataset\n')
        file_object.write('\n')
    
    return aperEdgeDF, disaperEdgeDF, upEdgeDF, downEdgeDF, eqEdgeDF

#start to construct the delta network
def main(refFolder, targetFolder, signalType, weightType):
    #load data
    try:
        testClusterMapDF = pd.read_csv(os.path.join(targetFolder, 'ClusterMapping.csv'), index_col=None, header=0)
        refClusterMapDF = pd.read_csv(os.path.join(refFolder, 'ClusterMapping.csv'), index_col=None, header=0)
        testEdgeDF = pd.read_csv(os.path.join(targetFolder, 'Edges_%s.csv' % (signalType)), index_col=None, header=0)
        refEdgeDF = pd.read_csv(os.path.join(refFolder, 'Edges_%s.csv' % (signalType)), index_col=None, header=0)
        
        # only keep edges of interest
        if weightType == 'mean':
            selectCols = ['Sending cluster', 'Ligand symbol', 'Receptor symbol', 'Target cluster', 'Ligand-expressing cells', 'Ligand detection rate', 'Ligand average expression value', 
            'Ligand derived specificity of average expression value', 'Receptor-expressing cells', 'Receptor detection rate', 'Receptor average expression value', 
            'Receptor derived specificity of average expression value', 'Edge average expression weight', 'Edge average expression derived specificity']
        
        elif weightType == 'sum':
            selectCols =  ['Sending cluster', 'Ligand symbol', 'Receptor symbol', 'Target cluster', 'Ligand-expressing cells', 'Ligand detection rate', 'Ligand total expression value', 
            'Ligand derived specificity of total expression value', 'Receptor-expressing cells', 'Receptor detection rate', 'Receptor total expression value', 
            'Receptor derived specificity of total expression value', 'Edge total expression weight', 'Edge total expression derived specificity']
        
        refEdgeDF = refEdgeDF.ix[:,selectCols]
        testEdgeDF = testEdgeDF.ix[:,selectCols]
        newColNames = ["sending cluster name", "ligand", "receptor", "target cluster name", "count ligand", "frequency ligand", "original ligand", "specified ligand", "count receptor", "frequency receptor", "original receptor", "specified receptor", "product of original", "product of specified"]
        refEdgeDF.columns = newColNames
        testEdgeDF.columns = newColNames
        print('all data are loaded, reference dataset has %s edges, target dataset has %s edges' % (len(refEdgeDF), len(testEdgeDF)))
    
    except Exception as e: 
        print(e)
        return
        
    #build the folder to save the analysis results
    resultFolder, rmf = BuildSaveFolder(refFolder, targetFolder)
    print('the folder "%s" is made to save the analysis results' % resultFolder)
    with open(rmf, 'a') as file_object:
        file_object.write('Reference dataset has %s edges\n' % len(refEdgeDF))
        file_object.write('Target dataset has %s edges\n' % len(testEdgeDF))
        file_object.write('\n')
            
    # find population changes in clusters
    sizeClusterDF = IdentifyPopulationChanges(refClusterMapDF, testClusterMapDF, resultFolder, rmf)
    print('all population changes have been identified')
    
    # calculate changes in the cell-to-cell communications 
    aperEdgeDF, disaperEdgeDF, upEdgeDF, downEdgeDF, eqEdgeDF = IdentifyLREdgeChanges(sizeClusterDF, refEdgeDF, testEdgeDF, signalType, weightType, resultFolder, rmf)
    print('all signaling changes have been identified')
    

#===================================================================#
#-------------------------------START-------------------------------#
#===================================================================#

'''
python DiffEdges.py --refFolder /Users/rhou/Public/c2cnetworkwebtool/CLI/tbs/3m.em --targetFolder /Users/rhou/Public/c2cnetworkwebtool/CLI/tbs/18m.em --signalType LRC2.0r --coreNum 1
python DiffEdges.py --refFolder /Users/rhou/Public/c2cnetworkwebtool/CLI/tbs/18m.em --targetFolder /Users/rhou/Public/c2cnetworkwebtool/CLI/tbs/21m.em --signalType LRC2.0r --coreNum 1
'''

if __name__ == '__main__':
    #process arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('--refFolder', required=True, help='the path to the folder of the reference dataset')
    parser.add_argument('--targetFolder', required=True, help='the path to the folder of the target dataset')
    parser.add_argument('--signalType', default='lrc2p', help='lrc2p (default) | lrc2a, folder name of the interaction database')
    parser.add_argument('--weightType', default='mean', help="mean (default) | sum")
    
    opt = parser.parse_args()
    
    #check signalType
    signalType = os.path.basename(opt.signalType)
        
    #check weightType
    avaWeightTypeList = ['mean', 'sum']
    if opt.weightType.lower() not in avaWeightTypeList:
        sys.exit("The weightType can only be 'mean' or 'sum'.")
    else:
        weightType = opt.weightType.lower()
        
    #check refFolder
    if not os.path.exists(opt.refFolder):
        sys.exit("The folder of the reference dataset does not exist.")
    if not os.path.exists('%s/Edges_%s.csv' % (opt.refFolder,opt.signalType)):
        sys.exit("Cannot find the 'Edges_%s.csv' file in the reference dataset folder." % opt.signalType)
    if not os.path.exists(os.path.join(opt.refFolder,'ClusterMapping.csv')):
        sys.exit("Cannot find the 'ClusterMapping.csv' file in the reference dataset folder.")
    
    #check targetFolder
    if not os.path.exists(opt.targetFolder):
        sys.exit("The folder of the target dataset does not exist.")
    if not os.path.exists('%s/Edges_%s.csv' % (opt.targetFolder,opt.signalType)):
        sys.exit("Cannot find the 'Edges_%s.csv' file in the target dataset folder." % opt.signalType)
    if not os.path.exists(os.path.join(opt.targetFolder,'ClusterMapping.csv')):
        sys.exit("Cannot find the 'ClusterMapping.csv' file in the target dataset folder.")
    
    #pass argument check, show input data
    print('===================================================')
    print('Input data:')
    print('The reference dataset: %s' % opt.refFolder)
    print('The target dataset: %s' % opt.targetFolder)
    print('The cell-to-cell signaling type: %s' % signalType)
    print('The weight type of cell-to-cell signaling: %s' % weightType)
    print('===================================================')
    
    #start to construct the delta network
    main(opt.refFolder, opt.targetFolder, signalType, weightType)