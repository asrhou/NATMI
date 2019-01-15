#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Aug  8 16:02:45 2018

@author: rhou
"""

import argparse
import pandas as pd
import os, sys
import multiprocessing
from functools import partial

# transfer hid to gene symbols of the certain species
def TransferToGeneSymbol(homoMapDir, speciestype, taxidCol, geneSymbolCol, hidCol, lrM):
    homoDF = pd.read_csv(homoMapDir, sep='\t', index_col=None, header=None)
    # reduce the list to genes of the given species
    humanDF = homoDF.ix[homoDF[taxidCol] == 9606,]
    homoDF = homoDF.ix[homoDF[taxidCol] == int(speciestype),]

    ligandGIDList = list(lrM.index.values)
    receptorGIDList = list(lrM.columns.values)
    
    lhumanDF = humanDF.ix[humanDF[geneSymbolCol].isin(ligandGIDList), [hidCol, geneSymbolCol]]
    lhumanDF.columns = ['hid', 'gene']
    lhumanDF = lhumanDF.drop_duplicates(subset=['hid']).drop_duplicates(subset=['gene'])
    rhumanDF = humanDF.ix[humanDF[geneSymbolCol].isin(receptorGIDList), [hidCol, geneSymbolCol]]
    rhumanDF.columns = ['hid', 'gene']
    rhumanDF = rhumanDF.drop_duplicates(subset=['hid']).drop_duplicates(subset=['gene'])
    lrM = lrM.ix[lhumanDF.ix[:, 'gene'], rhumanDF.ix[:, 'gene']]
    lrM.index = lhumanDF.ix[:, 'hid']
    lrM.columns = rhumanDF.ix[:, 'hid']
    
    lhomoDF = homoDF.ix[homoDF[hidCol].isin(lhumanDF.ix[:, 'hid']), [hidCol, geneSymbolCol]]
    lhomoDF.columns = ['hid', 'gene']
    lhomoDF = lhomoDF.drop_duplicates(subset=['hid']).drop_duplicates(subset=['gene'])
    rhomoDF = homoDF.ix[homoDF[hidCol].isin(rhumanDF.ix[:, 'hid']), [hidCol, geneSymbolCol]]
    rhomoDF.columns = ['hid', 'gene']
    rhomoDF = rhomoDF.drop_duplicates(subset=['hid']).drop_duplicates(subset=['gene'])
    
    lrM = lrM.ix[lhomoDF.ix[:, 'hid'], rhomoDF.ix[:, 'hid']]
    lrM.index = lhomoDF.ix[:, 'gene']
    lrM.columns = rhomoDF.ix[:, 'gene']
    return lrM

#cluster expression matrix
def ClusterAnnotateEM(resultDir, emDF, ann, weightType):
    # get cluster list
    clusterIdList = sorted(list(set(ann.ix[:, 'cluster'].tolist())))
    
    # calculate the expressions for each cluster
    counttableDFList = []
    for clusterId in clusterIdList:
        # get the sub dataframe of the cluster
        cellsInClusterList = list(ann.index[ann['cluster'] == clusterId])
        clusterDF = emDF.ix[:, cellsInClusterList]

        # replace headers for the cluster
        if weightType == 'mean':
            DF = clusterDF.mean(axis=1).to_frame(name=clusterId)
        elif weightType == 'sum':
            DF = clusterDF.sum(axis=1).to_frame(name=clusterId)
#        elif weightType == 'es':
#            DF = ((clusterDF.mean(axis=1)+1)/(emDF.median(axis=1)+1)).to_frame(name=clusterId)
            
        # add to the final dataframes
        counttableDFList.append(DF)

    # merge results and save
    counttableDF = pd.concat(counttableDFList, axis=1)
    # save result
    counttableFileDir = os.path.join(resultDir, '%sClusteredEM.csv' % weightType)
    counttableDF = counttableDF.round(4)
    counttableDF.to_csv(counttableFileDir)

    #return sumCounttableDF, meanCounttableDF
    return counttableDF

#generate ligand-receptor pair list
def GenLigandReceptorList(pairsDF):
    ligandList = []
    receptorList = []
    ligandNameList = pairsDF.index.values
    for ligandName in ligandNameList:
        pairedRecptors = list(pairsDF.ix[ligandName, pairsDF.ix[ligandName,] > 0].index.values)
        subLigand = [ligandName] * len(pairedRecptors)
        ligandList = ligandList +subLigand
        receptorList = receptorList + pairedRecptors
    pairListDF = pd.DataFrame({'ligand':ligandList, 'receptor':receptorList})
    return pairListDF

#split file to ligand and receptor file
def SplitIntoSinalProteins(dataDF, typeString, dataString, countTableFolder, ligandApprovedSymbolList, receptorApprovedSymbolList):
    #generate file name
    ligandFile = typeString + '_' + dataString + '_ligand.csv'
    receptorFile = typeString + '_' + dataString + '_receptor.csv'

    #split ligand and receptor
    ligandDF = dataDF.ix[ligandApprovedSymbolList,:].fillna(0.0)
    receptorDF = dataDF.ix[receptorApprovedSymbolList,:].fillna(0.0)
    
    #save file
    ligandDF = ligandDF.ix[ligandDF.sum(axis=1)>0,]
    ligandDF.to_csv(os.path.join(countTableFolder, ligandFile))
    receptorDF = receptorDF.ix[receptorDF.sum(axis=1)>0,]
    receptorDF.to_csv(os.path.join(countTableFolder, receptorFile))

    return ligandDF, receptorDF

#find cells for ligand and receptor
def FindCellsOfProtein(protein, proteinDF, specifiedProteinDF, combinedProteinDF):
    subproteinDF = proteinDF.ix[proteinDF.index == protein, :].T
    subspecifiedProteinDF = specifiedProteinDF.ix[specifiedProteinDF.index == protein, :].T
    subcombinedProteinDF = combinedProteinDF.ix[combinedProteinDF.index == protein, :].T
    return subproteinDF, subspecifiedProteinDF, subcombinedProteinDF

#find cells for ligand and receptor
def BuildHalfEdge(protein, proteinDF, specifiedProteinDF, combinedProteinDF, cellRole, proteinType):
    proteinCellListDF, proteinCellListSpecifiedDF, proteinCellListCombinedDF = FindCellsOfProtein(protein, proteinDF, specifiedProteinDF, combinedProteinDF)
    mergedProteinDF = pd.concat([proteinCellListDF, proteinCellListSpecifiedDF, proteinCellListCombinedDF], axis=1)

    #update headers
    newHeaders = ['original '+proteinType, 'specified '+proteinType, 'combined '+proteinType]
    mergedProteinDF.columns = newHeaders
    mergedProteinDF[cellRole+' name'] = mergedProteinDF.index.values
    mergedProteinDF[proteinType] = protein
    mergedProteinDF.reset_index(drop=True, inplace=True)
    
    return mergedProteinDF

#construct cell-to-cell edges
def GenSingleCell2CellEdge(ligandApprovedSymbolDict, receptorApprovedSymbolDict, pair):
    edgeList = []
    ligand = pair[0]
    receptor = pair[1]
    if ligand not in ligandApprovedSymbolDict.keys() or receptor not in receptorApprovedSymbolDict.keys():
        return edgeList
    ligandCellListDF = ligandApprovedSymbolDict[ligand]
    receptorCellListDF = receptorApprovedSymbolDict[receptor]
    for ligandCellListIndex in ligandCellListDF.index:
        for receptorCellListIndex in receptorCellListDF.index:
            ligandCellDF = ligandCellListDF.ix[ligandCellListIndex, :].to_frame().T.reset_index(drop=True)
            receptorCellDF = receptorCellListDF.ix[receptorCellListIndex, :].to_frame().T.reset_index(drop=True)
            mergedPairDF = pd.concat([ligandCellDF, receptorCellDF], axis=1)
            edgeList.append(mergedPairDF)
    return edgeList

#construct cell-to-cell edges
def GenCell2CellEdges(pairListDF, ligandApprovedSymbolDict, receptorApprovedSymbolDict, coreNum):
    pairList = [tuple(x) for x in pairListDF.values]
    
    p = multiprocessing.Pool(coreNum)
    func = partial(GenSingleCell2CellEdge, ligandApprovedSymbolDict, receptorApprovedSymbolDict)
    edgeList = p.map(func, pairList)
    p.close()
    p.join()
    
    edgeList = [item for sublist in edgeList for item in sublist]
    edgeListDF = pd.concat(edgeList, axis=0)
    #add prodect columns
    edgeListDF['product of original'] = edgeListDF['original ligand'] * edgeListDF['original receptor']
    edgeListDF['product of specified'] = edgeListDF['specified ligand'] * edgeListDF['specified receptor']
    edgeListDF['product of combined'] = edgeListDF['combined ligand'] * edgeListDF['combined receptor']
    #remove values equal zero
    edgeListDF = edgeListDF[edgeListDF['product of original'] > 0]
    #round values
    return edgeListDF

def GenSendingCellPairs(edgeListDF, targetCellList, sendingCell):
    cellPairDFTempList = []
    for targetCell in targetCellList:
        #find cell pair in the full list
        cellPairDF = edgeListDF[(edgeListDF['sending cluster name']==sendingCell) & (edgeListDF['target cluster name']==targetCell)]
        if len(cellPairDF) > 0:
            cellPairDF = cellPairDF[['product of original', 'product of specified', 'product of combined']]
            sumCellPairDF = cellPairDF.sum(axis=0).to_frame().T
            sumCellPairDF['sending cluster name'] = sendingCell
            sumCellPairDF['target cluster name'] = targetCell
            cellPairDFTempList.append(sumCellPairDF)
            
    #calculate the protortion of the sending cell to target cells
    cellPairTempDF = pd.concat(cellPairDFTempList)
    
    cellPairTempDF['proportion of original product'] = cellPairTempDF['product of original'].div(cellPairTempDF['product of original'].sum(), axis=0).fillna(0.0)
    cellPairTempDF['proportion of specified product'] = cellPairTempDF['product of specified'].div(cellPairTempDF['product of specified'].sum(), axis=0).fillna(0.0)
    cellPairTempDF['proportion of combined product'] = cellPairTempDF['product of combined'].div(cellPairTempDF['product of combined'].sum(), axis=0).fillna(0.0)        
    
    return cellPairTempDF
    
#construct and save cell-to-cell relations
def GenCellCellPairs(edgeListDF, coreNum):
    cellPairDFList = []
    sendingCellList = list(set(edgeListDF.ix[:,'sending cluster name'].tolist()))
    targetCellList = list(set(edgeListDF.ix[:,'target cluster name'].tolist()))
       
    p = multiprocessing.Pool(coreNum)
    func = partial(GenSendingCellPairs, edgeListDF, targetCellList)
    cellPairDFList = p.map(func, sendingCellList)
    p.close()
    p.join()
    pairListDF = pd.concat(cellPairDFList)
    
    return pairListDF

# generate cell-ligand-receptor-cell table
def GenerateCell2CellTable(emDF, typeString, pairsDF, resultDir, coreNum):
    # generate ligand-receptor pair list
    pairListDF = GenLigandReceptorList(pairsDF)
    ligandApprovedSymbolList = list(pairsDF.index.values)
    receptorApprovedSymbolList = list(pairsDF.columns.values)
    
    print '#### extract signaling proteins from %s' % typeString.split('_')[0]
    # split original and speciefied countable to ligand and receptor file
    ligandDF, receptorDF = SplitIntoSinalProteins(emDF, typeString, 'original', resultDir, ligandApprovedSymbolList, receptorApprovedSymbolList)

    # calculate specificity
    specifiedCounttableDF = emDF.div(emDF.sum(axis=1), axis=0).fillna(0.0)
    specifiedLigandDF, specifiedReceptorDF = SplitIntoSinalProteins(specifiedCounttableDF, typeString, 'specified', resultDir, ligandApprovedSymbolList, receptorApprovedSymbolList)

    # calculate combined expression levels
    combinedCounttableDF = emDF.mul(specifiedCounttableDF, fill_value=1).fillna(0.0)
    combinedLigandDF, combinedReceptorDF = SplitIntoSinalProteins(combinedCounttableDF, typeString, 'combined', resultDir, ligandApprovedSymbolList, receptorApprovedSymbolList)

    # find cells for ligand and receptor
    ligandApprovedSymbolDict = {}
    for ligandApprovedSymbol in ligandDF.index:#ligandApprovedSymbolList:
        mergedLigandDF = BuildHalfEdge(ligandApprovedSymbol, ligandDF, specifiedLigandDF, combinedLigandDF, 'sending cluster', 'ligand')
        ligandApprovedSymbolDict[ligandApprovedSymbol] = mergedLigandDF
    receptorApprovedSymbolDict = {}
    for receptorApprovedSymbol in receptorDF.index:#receptorApprovedSymbolList:
        mergedReceptorDF = BuildHalfEdge(receptorApprovedSymbol, receptorDF, specifiedReceptorDF, combinedReceptorDF, 'target cluster', 'receptor')
        receptorApprovedSymbolDict[receptorApprovedSymbol] = mergedReceptorDF

    print '#### construct signaling interactions from %s' % typeString.split('_')[0]
    # construct and save cell-to-cell edge
    edgeListDF = GenCell2CellEdges(pairListDF, ligandApprovedSymbolDict, receptorApprovedSymbolDict, coreNum)
    if len(edgeListDF) > 0:
        edgeListFileDir = os.path.join(resultDir, typeString + '_signaling interactions.csv')
        edgeListDF.round(4).to_csv(edgeListFileDir, index=False,
                               columns=["sending cluster name", "target cluster name", "ligand", "receptor", 
                                        "original ligand", "specified ligand", "combined ligand", 
                                        "original receptor", "specified receptor", "combined receptor", 
                                        "product of original", "product of specified", "product of combined"])

        print '#### construct cluster interactions from %s' % typeString.split('_')[0]
        # construct and save cell-to-cell relations
        pairListDF = GenCellCellPairs(edgeListDF, coreNum)
        pairListFileDir = os.path.join(resultDir, typeString + '_cluster interactions.csv')
        pairListDF.round(4).to_csv(pairListFileDir, index=False,
                                   columns=['sending cluster name', 'target cluster name', 'product of original',
                                            'product of specified', 'product of combined', 'proportion of original product',
                                            'proportion of specified product', 'proportion of combined product'])
    else:
        print '#### No signaling interactions found from %s' % typeString.split('_')[0]
    
    return edgeListDF, pairListDF

# generate cell-to-cell interaction files
def GenerateDataFiles(emDF, typeString, lrM, ann, resultDir, coreNum):
    
    # generate cell-ligand-receptor-cell table
    edgeListDF, pairListDF = GenerateCell2CellTable(emDF, typeString, lrM, resultDir, coreNum)
    
def main(species, emFile, annFile, signalType, weightType, coreNum):
    #load data
    emFileExtention = emFile.split('.')[-1]
    if emFileExtention == 'csv':
        em = pd.read_csv(emFile, index_col=0, header=0)
    elif emFileExtention == 'tsv' or emFileExtention == 'txt':
        em = pd.read_table(emFile, index_col=0, header=0)
    else:
        sys.exit("Cannot process the expression matrix file, please check the format of the expression matrix file.")
    
    annFileExtention = annFile.split('.')[-1]
    if annFileExtention == 'csv':
        ann = pd.read_csv(annFile, index_col=0, header=0)
    elif annFileExtention == 'tsv' or emFileExtention == 'txt':
        ann = pd.read_table(annFile, index_col=0, header=0)
    else:
        sys.exit("Cannot process the metafile, please check the format of the metafile.")
    ann.columns = ['cluster']
    
    #load interaction list
    lrM = pd.read_excel('%s/pairsM.xlsx' % signalType, index_col=0, header=0)
    
    # change gene symbols if necessary
    if species != '9606':
        #find taxonomy file
        ## HID (HomoloGene group id)[0] - Taxonomy ID[1] - Gene ID[2] - Gene Symbol[3] - Protein gi[4] - Protein accession[5]
        homoMapDir = 'homology/homologene.data'
        hidCol = 0
        taxidCol = 1
        geneSymbolCol = 3
        # transfer hid to gene symbols of the certain species
        lrM = TransferToGeneSymbol(homoMapDir, species, taxidCol, geneSymbolCol, hidCol, lrM)
    
    #build the folder to save the analysis results
    resultDir = emFile.split('.'+emFileExtention)[0]
    if not os.path.exists(resultDir):
        os.mkdir(resultDir)
    saveann = ann.reset_index()
    saveann.columns = ['cell', 'cluster']
    saveann.to_csv(os.path.join(resultDir,'clusterMapping.csv'), index=False, header=True)
    
    #cluster expression matrix
    print '#### cluster expression matrix'
    EMDF = ClusterAnnotateEM(resultDir, em, ann, weightType)
    
    # generate cell-to-cell interaction files
    GenerateDataFiles(EMDF, '%sClusteredEM_%s' % (weightType, os.path.basename(signalType)), lrM, ann, resultDir, coreNum)
    
#===================================================================#
#-------------------------------START-------------------------------#
#===================================================================#

'''
python NATMI.py --species mouse --emFile M6C_data.frame.txt --annFile M6C_metadata.frame.txt --signalType paracrine --coreNum 1
python NATMI.py --species mouse --emFile M6C_data.frame.txt --annFile M6C_metadata.frame.txt --signalType juxtacrine --coreNum 1
'''

if __name__ == '__main__':
    #process arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('--species', default='human', help='human(default) | mouse') 
    parser.add_argument('--emFile', required=True, help='the path to the file of the expression matrix with row names and column names')
    parser.add_argument('--annFile', required=True, help='the path to the metafile in which row one are the column names, column one are cellbarcodes and column two are corresponding cluster IDs')
    parser.add_argument('--signalType', default='paracrine', help='paracrine (default) | juxtacrine | folder name of the interaction database')
    parser.add_argument('--weightType', default='mean', help="mean (default) | sum | es (enrichment score)")
    parser.add_argument('--coreNum', type=int, default=1, help='number of the cores to use, default is one')
    
    opt = parser.parse_args()
    
    #check species
    avaSpecDict = {'human':'9606', 'mouse':'10090'}
    if opt.species.lower() not in avaSpecDict.keys():
        sys.exit("The species can only be 'human' or 'mouse'.")
    else:
        species = avaSpecDict[opt.species.lower()]
        
    #check signalType
    avaSignalTypeList = ['paracrine', 'juxtacrine']
    if opt.signalType.lower() not in avaSignalTypeList:
        if not os.path.exists('%s/pairsM.xlsx' % opt.signalType):
            sys.exit("The signalType can only be 'paracrine', 'juxtacrine' or the folder of interaction database.")
        else:
            signalType = opt.signalType
    else:
        signalType = opt.signalType.lower()
        
    #check weightType
    avaWeightTypeList = ['mean', 'sum', 'es']
    if opt.weightType.lower() not in avaWeightTypeList:
        sys.exit("The weightType can only be 'mean' or 'sum'.")
    else:
        weightType = opt.weightType.lower()
        
    #check emFile
    if not os.path.exists(opt.emFile):
        sys.exit("The expression matrix file does not exist.")
        
    #check annFile
    if not os.path.exists(opt.annFile):
        sys.exit("The metafile does not exist.")
        
    #check coreNum
    maxCoreNum = multiprocessing.cpu_count()
    if opt.coreNum > maxCoreNum:
        sys.exit("There are only %s cores availble, less than %s cores." % (maxCoreNum, opt.coreNum))
        
    #pass argument check, show input data
    print '==================================================='
    print 'Input data:'
    print 'Species: %s' % opt.species.lower()
    print 'The expression matrix file: %s' % opt.emFile
    print 'The metafile: %s' % opt.annFile
    print 'The cell-to-cell signaling type: %s' % signalType
    print 'The weight type of cell-to-cell signaling: %s' % weightType
    print 'The number of cores to use: %s' % opt.coreNum
    print '==================================================='
    
    #start to construct cell-to-cell network
    main(species, opt.emFile, opt.annFile, signalType, weightType, opt.coreNum)