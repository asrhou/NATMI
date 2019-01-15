#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Aug  9 21:40:04 2018

@author: rhou
"""

import argparse
import pandas as pd
import os, sys
import multiprocessing
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
matplotlib.style.use('ggplot')
import seaborn as sns

#build the folder to save the analysis results
def BuildSaveFolder(refFolder, targetFolder):
    ref = os.path.basename(refFolder)
    target = os.path.basename(targetFolder)
    resultFolder = os.path.join(os.path.dirname(refFolder), "delta_interactions_ref_%s_target_%s" % (ref, target))
    if not os.path.exists(resultFolder):
        os.mkdir(resultFolder)
    return resultFolder

# find population changes in clusters
def IdentifyPopulationChanges(refClusterMapDF, testClusterMapDF, resultFolder):
    refClusterSizeDF = refClusterMapDF.groupby('cluster').count()
    refClusterSizeDF['class'] = ['reference dataset'] * len(refClusterSizeDF)
    refClusterSizeDF = refClusterSizeDF.reset_index()
    
    testClusterSizeDF = testClusterMapDF.groupby('cluster').count()
    testClusterSizeDF['class'] = ['target dataset'] * len(testClusterSizeDF)
    testClusterSizeDF = testClusterSizeDF.reset_index()
    
    sumClusterDF = pd.concat([refClusterSizeDF, testClusterSizeDF])
    sumClusterDF.columns = ['cluster', 'population (cells)', 'source']
    sumClusterDF = sumClusterDF.sort_values(by=['cluster'])
    sumClusterDF.to_csv(os.path.join(resultFolder,'cluster_size_comparison.csv'), index=False, header=True)
    
    # find all clusters
    allClusters = set(sumClusterDF['cluster'])
    plotHeight = len(allClusters) * 2 / 5
    plotWidth = 2 * plotHeight
    
    # sort clusters based on the fold changes
    cltOrder = []
    for clt in allClusters:
        cltDF = sumClusterDF.ix[sumClusterDF['cluster'] == clt,].sort_values(by=['source'])
        if len(cltDF) == 2:
            cltOrder.append((clt,float(cltDF.iloc[1,1])/cltDF.iloc[0,1]))
        else:
            cltOrder.append((clt,0))
    cltOrder = sorted(cltOrder, key=lambda tup: tup[1], reverse=True)
    cltOrder = [i[0] for i in cltOrder]
    
    # visualise population changes in clusters
    sns.set_style("whitegrid")
    fig = plt.figure(figsize=(plotWidth, plotHeight))
    ax = sns.barplot(y="cluster", x="population (cells)", hue="source", data=sumClusterDF, order=cltOrder);
    ax.legend(bbox_to_anchor=(1.02, 0.1))
    fig = ax.get_figure()
    fig.savefig(os.path.join(resultFolder,'cluster_size_comparison.pdf'), bbox_inches='tight')

# calculate changes in the cell-to-cell communications 
def IdentifyLREdgeChanges(refEdgeDF, testEdgeDF, signalType, weightType, resultFolder):
    mergedDF = pd.merge(refEdgeDF, testEdgeDF, how='outer', on=['sending cluster name', 'target cluster name', 'ligand', 'receptor']).fillna(0)
    aperEdgeDF = mergedDF.ix[(mergedDF['product of original_x']==0)&(mergedDF['product of original_y']>0),['sending cluster name', 'target cluster name', 'ligand', 'receptor', 'original ligand_y', 'specified ligand_y', 'original receptor_y', 'specified receptor_y', 'product of original_y', 'product of specified_y']]
    aperEdgeDF.columns = ['sending cluster name', 'target cluster name', 'ligand', 'receptor', 'delta ligand expression', 'delta ligand specificity', 'delta receptor expression', 'delta receptor specificity', 'delta weight', 'delta specificity']
    aperEdgeDF = aperEdgeDF.sort_values(by="delta weight", ascending=False)
    aperEdgeDF.to_csv(os.path.join(resultFolder,'appeared_signalType_%s_weightType_%s_signaling interactions.csv' % (signalType, weightType)), index=False, header=True, columns=['sending cluster name', 'target cluster name', 'ligand', 'receptor', 'delta ligand expression', 'delta ligand specificity', 'delta receptor expression', 'delta receptor specificity', 'delta weight', 'delta specificity'])
    print '%s interactions appeared' % "{:,}".format(len(aperEdgeDF))
    
    disaperEdgeDF = mergedDF.ix[(mergedDF['product of original_y']==0)&(mergedDF['product of original_x']>0),['sending cluster name', 'target cluster name', 'ligand', 'receptor', 'original ligand_x', 'specified ligand_x', 'original receptor_x', 'specified receptor_x', 'product of original_x', 'product of specified_x']]
    disaperEdgeDF.columns = ['sending cluster name', 'target cluster name', 'ligand', 'receptor', 'delta ligand expression', 'delta ligand specificity', 'delta receptor expression', 'delta receptor specificity', 'delta weight', 'delta specificity']
    disaperEdgeDF = disaperEdgeDF.sort_values(by="delta weight", ascending=False)
    disaperEdgeDF.to_csv(os.path.join(resultFolder,'disappeared_signalType_%s_weightType_%s_signaling interactions.csv' % (signalType, weightType)), index=False, header=True, columns=['sending cluster name', 'target cluster name', 'ligand', 'receptor', 'delta ligand expression', 'delta ligand specificity', 'delta receptor expression', 'delta receptor specificity', 'delta weight', 'delta specificity'])
    print '%s interactions disappeared' % "{:,}".format(len(disaperEdgeDF))
    
    upEdgeDF = mergedDF.ix[(mergedDF['product of original_x'] < mergedDF['product of original_y'])&(mergedDF['product of original_x']>0),['sending cluster name', 'target cluster name', 'ligand', 'receptor', 'original ligand_x', 'specified ligand_x', 'original receptor_x', 'specified receptor_x', 'product of original_x', 'product of specified_x', 'original ligand_y', 'specified ligand_y', 'original receptor_y', 'specified receptor_y', 'product of original_y', 'product of specified_y']]
    upEdgeDF.columns = ['sending cluster name', 'target cluster name', 'ligand', 'receptor', 'original ligand expression', 'original ligand specificity', 'original receptor expression', 'original receptor specificity', 'original weight', 'original specificity', 'changed ligand expression', 'changed ligand specificity', 'changed receptor expression', 'changed receptor specificity', 'changed weight', 'changed specificity']
    upEdgeDF['delta ligand expression'] = upEdgeDF['changed ligand expression'] - upEdgeDF['original ligand expression']
    upEdgeDF['delta receptor expression'] = upEdgeDF['changed receptor expression'] - upEdgeDF['original receptor expression']
    upEdgeDF['delta ligand specificity'] = upEdgeDF['changed ligand specificity'] - upEdgeDF['original ligand specificity']
    upEdgeDF['delta receptor specificity'] = upEdgeDF['changed receptor specificity'] - upEdgeDF['original receptor specificity']
    upEdgeDF['delta weight'] = upEdgeDF['changed weight'] - upEdgeDF['original weight']
    upEdgeDF['delta specificity'] = upEdgeDF['changed specificity'] - upEdgeDF['original specificity']
    upEdgeDF['fold change of weight'] = upEdgeDF['changed weight'] / upEdgeDF['original weight']
    upEdgeDF['fold change of specificity'] = upEdgeDF['changed specificity'] / upEdgeDF['original specificity']
    upEdgeDF = upEdgeDF.sort_values(by="fold change of weight", ascending=False)
    upEdgeDF.to_csv(os.path.join(resultFolder,'up_regulated_signalType_%s_weightType_%s_signaling interactions.csv' % (signalType, weightType)), index=False, header=True, columns=['sending cluster name', 'target cluster name', 'ligand', 'receptor', 'original ligand expression', 'original ligand specificity', 'original receptor expression', 'original receptor specificity', 'original weight', 'original specificity', 'changed ligand expression', 'changed ligand specificity', 'changed receptor expression', 'changed receptor specificity', 'changed weight', 'changed specificity', 'delta ligand expression', 'delta ligand specificity', 'delta receptor expression', 'delta receptor specificity', 'delta weight', 'delta specificity', 'fold change of weight', 'fold change of specificity'])
    print '%s interactions are up-regulated' % "{:,}".format(len(upEdgeDF))
    
    downEdgeDF = mergedDF.ix[(mergedDF['product of original_x'] > mergedDF['product of original_y'])&(mergedDF['product of original_y']>0),['sending cluster name', 'target cluster name', 'ligand', 'receptor', 'original ligand_x', 'specified ligand_x', 'original receptor_x', 'specified receptor_x', 'product of original_x', 'product of specified_x', 'original ligand_y', 'specified ligand_y', 'original receptor_y', 'specified receptor_y', 'product of original_y', 'product of specified_y']]
    downEdgeDF.columns = ['sending cluster name', 'target cluster name', 'ligand', 'receptor', 'original ligand expression', 'original ligand specificity', 'original receptor expression', 'original receptor specificity', 'original weight', 'original specificity', 'changed ligand expression', 'changed ligand specificity', 'changed receptor expression', 'changed receptor specificity', 'changed weight', 'changed specificity']
    downEdgeDF['delta ligand expression'] = downEdgeDF['original ligand expression'] - downEdgeDF['changed ligand expression']
    downEdgeDF['delta receptor expression'] = downEdgeDF['original receptor expression'] - downEdgeDF['changed receptor expression']
    downEdgeDF['delta ligand specificity'] = downEdgeDF['original ligand specificity'] - downEdgeDF['changed ligand specificity']
    downEdgeDF['delta receptor specificity'] = downEdgeDF['original receptor specificity'] - downEdgeDF['changed receptor specificity'] 
    downEdgeDF['delta weight'] = downEdgeDF['original weight'] - downEdgeDF['changed weight']
    downEdgeDF['delta specificity'] = downEdgeDF['original specificity'] - downEdgeDF['changed specificity']
    downEdgeDF['fold change of weight'] = downEdgeDF['original weight'] / downEdgeDF['changed weight']
    downEdgeDF['fold change of specificity'] = downEdgeDF['original specificity'] / downEdgeDF['changed specificity']
    downEdgeDF = downEdgeDF.sort_values(by="fold change of weight", ascending=False)
    downEdgeDF.to_csv(os.path.join(resultFolder,'down_regulated_signalType_%s_weightType_%s_signaling interactions.csv' % (signalType, weightType)), index=False, header=True, columns=['sending cluster name', 'target cluster name', 'ligand', 'receptor', 'original ligand expression', 'original ligand specificity', 'original receptor expression', 'original receptor specificity', 'original weight', 'original specificity', 'changed ligand expression', 'changed ligand specificity', 'changed receptor expression', 'changed receptor specificity', 'changed weight', 'changed specificity', 'delta ligand expression', 'delta ligand specificity', 'delta receptor expression', 'delta receptor specificity', 'delta weight', 'delta specificity', 'fold change of weight', 'fold change of specificity'])
    print '%s interactions are down-regulated' % "{:,}".format(len(downEdgeDF))
    
    return aperEdgeDF, disaperEdgeDF, upEdgeDF, downEdgeDF

def GenDynamicCell2CellEdge(subDF):
    rstDict = {'sending cluster name':[], 'target cluster name':[], 'edge count':[], 'delta weight':[], 'delta specificity':[], 'most variable pair':[], 'biggest delta weight of a pair':[]}
    if len(subDF)>0:
        rstDict['sending cluster name'].append(subDF.ix[:,'sending cluster name'].values[0])
        rstDict['target cluster name'].append(subDF.ix[:,'target cluster name'].values[0])
        rstDict['edge count'].append(len(subDF))
        rstDict['delta weight'].append(subDF['delta weight'].sum())
        rstDict['delta specificity'].append(subDF['delta specificity'].sum())
        maxValue = subDF.sort_values(by=['delta weight']).iloc[-1,]
        rstDict['most variable pair'].append(maxValue['ligand'] +' -> '+ maxValue['receptor'])
        rstDict['biggest delta weight of a pair'].append(maxValue['delta weight'])
        
    rstDF = pd.DataFrame(rstDict).ix[:,['sending cluster name', 'target cluster name', 'edge count', 'delta weight', 'delta specificity', 'most variable pair', 'biggest delta weight of a pair']]
    return rstDF

# collect cluster-to-cluster communicaiton
def AnaClt2CltComm(aperEdgeDF, disaperEdgeDF, upEdgeDF, downEdgeDF, signalType, weightType, coreNum, resultFolder):
    dfDict = {'appeared':aperEdgeDF, 'disappeared':disaperEdgeDF, 'up_regulated':upEdgeDF, 'down_regulated':downEdgeDF}
    rstDict = {}
    for kind in dfDict.keys():
        df = dfDict[kind]
        grouped = df.groupby(['sending cluster name', 'target cluster name'])
        groupedDFs = [i[1] for i in grouped.__iter__()]
    
        p = multiprocessing.Pool(coreNum)
        resultList = p.map(GenDynamicCell2CellEdge, groupedDFs)
        p.close()
        p.join()
        
        if len(resultList) > 0:
            allRstDF = pd.concat(resultList)
            allRstDF = allRstDF.round({'delta weight':2, 'delta specificity':2, 'biggest delta weight of a pair':2})
            allRstDF.to_csv(os.path.join(resultFolder,'%s_signalType_%s_weightType_%s_cluster interactions.csv' % (kind, signalType, weightType)), index=False, header=True, columns=['sending cluster name', 'target cluster name', 'edge count', 'delta weight', 'delta specificity', 'most variable pair', 'biggest delta weight of a pair'])
            rstDict[kind] = allRstDF
        else:
            rstDict[kind] = pd.DataFrame()
    return rstDict

# prepare color mapping for chord plots
def PrepareChordPlotColors(clt2CltCommDict):
    cltSet = set()
    for kind in clt2CltCommDict.keys():
        df = clt2CltCommDict[kind]
        cltSet = cltSet.union(set(df['sending cluster name']))
        cltSet = cltSet.union(set(df['target cluster name']))
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
def ComposeR(colorMapStr, clt2CltCommDict, signalType, weightType, resultFolder, plotFormat):
    RCodeFileName = os.path.join(resultFolder, 'delta_signalType_%s_weightType_%s_chord.r' % (signalType, weightType))
    RCodeFileName = os.path.abspath(RCodeFileName)
    
    with open(RCodeFileName, "wb") as outfile:
        outfile.write('#！/usr/bin/R\n')
        outfile.write('\n')
        outfile.write('library(circlize)\n')
        outfile.write(colorMapStr)
        
        # draw different type of plots
        for kind in clt2CltCommDict.keys():
            df = clt2CltCommDict[kind]
            tempclusterList = sorted(list(set(df.ix[:, 'sending cluster name']).union(set(df.ix[:, 'target cluster name']))))
            outfile.write('%sorder = c(%s)\n' % (kind, str(tempclusterList)[1:-1])) 
            
            outfile.write('#load adjacency list\n')
            filename = os.path.join(resultFolder,'%s_signalType_%s_weightType_%s_cluster interactions.csv' % (kind, signalType, weightType))
            filename = os.path.abspath(filename)
            outfile.write('ad <- read.csv("%s")\n' % (filename))
        
            outfile.write('#draw edge count graph\n')
            plotfilename = os.path.join(resultFolder,'%s_%s_chordplot_edge_count.%s' % (kind, signalType, plotFormat))
            plotfilename = os.path.abspath(plotfilename)
            outfile.write('%s("%s", width=%s, height=%s)\n' % (plotFormat, plotfilename, len(tempclusterList), len(tempclusterList)))
            outfile.write('chordDiagram(ad[,c(1,2,3)], order=%sorder, directional = 1, direction.type = c("diffHeight", "arrows"), grid.col = grid.col)\n' % (kind))
            outfile.write('dev.off()\n')
            outfile.write('circos.clear()\n')
            
            outfile.write('#draw delta weight graph\n')
            plotfilename = os.path.join(resultFolder,'%s_%s_chordplot_delta_weight.%s' % (kind, signalType, plotFormat))
            plotfilename = os.path.abspath(plotfilename)
            outfile.write('%s("%s", width=%s, height=%s)\n' % (plotFormat, plotfilename, len(tempclusterList), len(tempclusterList)))
            outfile.write('chordDiagram(ad[,c(1,2,4)], order=%sorder, directional = 1, direction.type = c("diffHeight", "arrows"), grid.col = grid.col)\n' % (kind))
            outfile.write('dev.off()\n')
            outfile.write('circos.clear()\n')
            
            outfile.write('#draw delta specificity graph\n')
            plotfilename = os.path.join(resultFolder,'%s_%s_chordplot_delta_specificity.%s' % (kind, signalType, plotFormat))
            plotfilename = os.path.abspath(plotfilename)
            outfile.write('%s("%s", width=%s, height=%s)\n' % (plotFormat, plotfilename, len(tempclusterList), len(tempclusterList)))
            outfile.write('chordDiagram(ad[,c(1,2,5)], order=%sorder, directional = 1, direction.type = c("diffHeight", "arrows"), grid.col = grid.col)\n' % (kind))
            outfile.write('dev.off()\n')
            outfile.write('circos.clear()\n')
        
    #draw chord graphs
    os.system('R CMD BATCH %s' % RCodeFileName)
           
# draw cluster-to-cluster communicaiton chord plot
def DrawChordPlots(clt2CltCommDict, signalType, weightType, resultFolder, plotFormat):
    # prepare color mapping for chord plots
    colorMapStr = PrepareChordPlotColors(clt2CltCommDict)
    
    # compose R script to draw chordplots
    ComposeR(colorMapStr, clt2CltCommDict, signalType, weightType, resultFolder, plotFormat)
    
#start to construct the delta network
def main(refFolder, targetFolder, signalType, weightType, coreNum, drawChordPlot, plotFormat):
    #load data
    try:
        testClusterMapDF = pd.read_csv(os.path.join(targetFolder, 'clusterMapping.csv'), index_col=None, header=0)
        refClusterMapDF = pd.read_csv(os.path.join(refFolder, 'clusterMapping.csv'), index_col=None, header=0)
        testEdgeDF = pd.read_csv(os.path.join(targetFolder, '%sClusteredEM_%s_signaling interactions.csv' % (weightType, signalType)), index_col=None, header=0)
        refEdgeDF = pd.read_csv(os.path.join(refFolder, '%sClusteredEM_%s_signaling interactions.csv' % (weightType, signalType)), index_col=None, header=0)
           
    except Exception as e: 
        print e
    print 'all data are loaded'
    
        
    #build the folder to save the analysis results
    resultFolder = BuildSaveFolder(refFolder, targetFolder)
    print 'the folder "%s" is made to save the analysis results' % resultFolder
    
    # find population changes in clusters
    IdentifyPopulationChanges(refClusterMapDF, testClusterMapDF, resultFolder)
    return
    print 'all population changes have been identified'
    
    # calculate changes in the cell-to-cell communications 
    aperEdgeDF, disaperEdgeDF, upEdgeDF, downEdgeDF = IdentifyLREdgeChanges(refEdgeDF, testEdgeDF, signalType, weightType, resultFolder)
    print 'all signaling changes have been identified'
    
    # collect cluster-to-cluster communicaiton
    clt2CltCommDict = AnaClt2CltComm(aperEdgeDF, disaperEdgeDF, upEdgeDF, downEdgeDF, signalType, weightType, coreNum, resultFolder)
    
    # draw cluster-to-cluster communicaiton chord plot
    if drawChordPlot:
        DrawChordPlots(clt2CltCommDict, signalType, weightType, resultFolder, plotFormat)
        print 'all chord plots have been drawn'

#===================================================================#
#-------------------------------START-------------------------------#
#===================================================================#

if __name__ == '__main__':
    #process arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('--refFolder', required=True, help='the path to the folder of the reference dataset')
    parser.add_argument('--targetFolder', required=True, help='the path to the folder of the target dataset')
    parser.add_argument('--signalType', default='paracrine', help='paracrine (default) | juxtacrine')
    parser.add_argument('--weightType', default='mean', help="mean (default) | sum | es (enrichment score)")
    parser.add_argument('--coreNum', type=int, default=1, help='number of the cores to use, default is one')
    parser.add_argument('--drawChordPlot', default='n', help='n(o) (default) | y(es)')
    parser.add_argument('--plotFormat', default='pdf', help="pdf (default) | png | svg")
    
    opt = parser.parse_args()
    
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
        
    #check refFolder
    if not os.path.exists(opt.refFolder):
        sys.exit("The folder of the reference dataset does not exist.")
        
    #check targetFolder
    if not os.path.exists(opt.targetFolder):
        sys.exit("The folder of the target dataset does not exist.")
        
    #check coreNum
    maxCoreNum = multiprocessing.cpu_count()
    if opt.coreNum > maxCoreNum:
        sys.exit("There are only %s cores availble, less than %s cores." % (maxCoreNum, opt.coreNum))
        
    #check if draw chord plot
    drawChordPlotStr = opt.drawChordPlot.lower()
    if drawChordPlotStr[0] == 'y':
        drawChordPlot = True
    elif drawChordPlotStr[0] == 'n':
        drawChordPlot = False
    else:
        sys.exit("drawChordPlot can only be 'y' or 'n'.")
        
    if drawChordPlot:
        #check plotformat
        avaFormatList = ["pdf", "png", "svg"]
        if opt.plotFormat.lower() not in avaFormatList:
            sys.exit("The given plot format is not permitted.")
        else:
            plotFormat = opt.plotFormat.lower()
    else:
        plotFormat = ''
        
    #pass argument check, show input data
    print '==================================================='
    print 'Input data:'
    print 'The reference dataset: %s' % opt.refFolder
    print 'The target dataset: %s' % opt.targetFolder
    print 'The cell-to-cell signaling type: %s' % signalType
    print 'The weight type of cell-to-cell signaling: %s' % weightType
    print 'The number of cores to use: %s' % opt.coreNum
    print 'Draw chord plots: %s' % drawChordPlotStr[0]
    if drawChordPlot:
        print 'The plot format: %s' % plotFormat
    print '==================================================='
    
    #start to construct the delta network
    main(opt.refFolder, opt.targetFolder, signalType, weightType, opt.coreNum, drawChordPlot, plotFormat)