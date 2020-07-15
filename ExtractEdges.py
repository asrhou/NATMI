#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sat Sep  7 23:17:45 2019

@author: rhou
"""

import warnings
warnings.filterwarnings("ignore")
import pandas as pd
import argparse, os, sys, glob
import multiprocessing
from functools import partial
    
# transfer hid to gene symbols of the certain species
def TransferToGeneSymbol(homoMapDir, speciestype, interSpeciestype, taxidCol, geneSymbolCol, hidCol, lrM):
    # load data
    homoDF = pd.read_csv(homoMapDir, sep='\t', index_col=None, header=None)
    # reduce the list to genes of the given species
    humanDF = homoDF.loc[homoDF[taxidCol] == int(interSpeciestype),]
    homoDF = homoDF.loc[homoDF[taxidCol] == int(speciestype),]

    ligandGIDList = list(lrM.index.values)
    receptorGIDList = list(lrM.columns.values)
    
    lhumanDF = humanDF.loc[humanDF[geneSymbolCol].isin(ligandGIDList), [hidCol, geneSymbolCol]]
    lhumanDF.columns = ['hid', 'gene']
    lhumanDF = lhumanDF.drop_duplicates(subset=['hid']).drop_duplicates(subset=['gene'])
    rhumanDF = humanDF.loc[humanDF[geneSymbolCol].isin(receptorGIDList), [hidCol, geneSymbolCol]]
    rhumanDF.columns = ['hid', 'gene']
    rhumanDF = rhumanDF.drop_duplicates(subset=['hid']).drop_duplicates(subset=['gene'])
    lrM = lrM.loc[lhumanDF.loc[:, 'gene'], rhumanDF.loc[:, 'gene']]
    lrM.index = lhumanDF.loc[:, 'hid']
    lrM.columns = rhumanDF.loc[:, 'hid']
    
    lhomoDF = homoDF.loc[homoDF[hidCol].isin(lhumanDF.loc[:, 'hid']), [hidCol, geneSymbolCol]]
    lhomoDF.columns = ['hid', 'gene']
    lhomoDF = lhomoDF.drop_duplicates(subset=['hid']).drop_duplicates(subset=['gene'])
    rhomoDF = homoDF.loc[homoDF[hidCol].isin(rhumanDF.loc[:, 'hid']), [hidCol, geneSymbolCol]]
    rhomoDF.columns = ['hid', 'gene']
    rhomoDF = rhomoDF.drop_duplicates(subset=['hid']).drop_duplicates(subset=['gene'])
    
    lrM = lrM.loc[lhomoDF.loc[:, 'hid'], rhomoDF.loc[:, 'hid']]
    lrM.index = lhomoDF.loc[:, 'gene']
    lrM.columns = rhomoDF.loc[:, 'gene']

    return lrM

def ClusterAnnotateEM(resultDir, emDF, ann):
    # get cluster list
    clusterIdList = sorted(list(set(ann.loc[:, 'cluster'].tolist())))
    
    # calculate the expressions for each cluster
    sumCounttableDFList = []
    meanCounttableDFList = []
    celltableDFList = []
    counttableDFList = []
    for clusterId in clusterIdList:
        # get the sub dataframe of the cluster
        cellsInClusterList = list(ann.index[ann['cluster'] == clusterId])
        clusterDF = emDF.loc[:, cellsInClusterList]
        
        # replace headers for the cluster
        sumDF = clusterDF.sum(axis=1).to_frame(name=clusterId)
        meanDF = clusterDF.mean(axis=1).to_frame(name=clusterId)
        # calculate the number of expressed cells
        cellDF = clusterDF[clusterDF>0].count(axis=1).astype(float).to_frame(name=clusterId)
        countDF = cellDF/len(cellsInClusterList)
        
        # add to the final dataframes
        sumCounttableDFList.append(sumDF)
        meanCounttableDFList.append(meanDF)
        celltableDFList.append(cellDF)
        counttableDFList.append(countDF)

    # merge results and save
    sumCounttableDF = pd.concat(sumCounttableDFList, axis=1)
    meanCounttableDF = pd.concat(meanCounttableDFList, axis=1)
    counttableDF = pd.concat(counttableDFList, axis=1)
    celltableDF = pd.concat(celltableDFList, axis=1)
    
    return sumCounttableDF, meanCounttableDF, counttableDF, celltableDF

def GenLigandReceptorList(pairsDF):
    ligandList = []
    receptorList = []
    ligandNameList = pairsDF.index.values
    for ligandName in ligandNameList:
        pairedRecptors = list(pairsDF.loc[ligandName, pairsDF.loc[ligandName,] > 0].index.values)
        subLigand = [ligandName] * len(pairedRecptors)
        ligandList = ligandList +subLigand
        receptorList = receptorList + pairedRecptors
    pairListDF = pd.DataFrame({'ligand':ligandList, 'receptor':receptorList})
    return pairListDF

def SplitIntoSinalProteins(sumEMDF, meanEMDF, countEMDF, cellEMDF, ligandApprovedSymbolList, receptorApprovedSymbolList):
    #split ligand and receptor
    ligandApprovedSymbolList = list(set(ligandApprovedSymbolList).intersection(set(sumEMDF.index)))
    receptorApprovedSymbolList = list(set(receptorApprovedSymbolList).intersection(set(sumEMDF.index)))
    sumligandDF = sumEMDF.loc[ligandApprovedSymbolList,:].fillna(0.0)
    sumreceptorDF = sumEMDF.loc[receptorApprovedSymbolList,:].fillna(0.0)
    meanligandDF = meanEMDF.loc[ligandApprovedSymbolList,:].fillna(0.0)
    meanreceptorDF = meanEMDF.loc[receptorApprovedSymbolList,:].fillna(0.0)
    countligandDF = countEMDF.loc[ligandApprovedSymbolList,:].fillna(0.0)
    countreceptorDF = countEMDF.loc[receptorApprovedSymbolList,:].fillna(0.0)
    cellligandDF = cellEMDF.loc[ligandApprovedSymbolList,:].fillna(0.0)
    cellreceptorDF = cellEMDF.loc[receptorApprovedSymbolList,:].fillna(0.0)
    
    # calculate specificity
    sumSpecifiedLigandDF = sumligandDF.div(sumligandDF.sum(axis=1), axis=0).fillna(0.0)
    sumSpecifiedReceptorDF = sumreceptorDF.div(sumreceptorDF.sum(axis=1), axis=0).fillna(0.0)
    meanSpecifiedLigandDF = meanligandDF.div(meanligandDF.sum(axis=1), axis=0).fillna(0.0)
    meanSpecifiedReceptorDF = meanreceptorDF.div(meanreceptorDF.sum(axis=1), axis=0).fillna(0.0)
    
    return cellligandDF, cellreceptorDF, countligandDF, countreceptorDF, sumligandDF, sumreceptorDF, meanligandDF, meanreceptorDF, sumSpecifiedLigandDF, sumSpecifiedReceptorDF, meanSpecifiedLigandDF, meanSpecifiedReceptorDF

def LRExpressions(typeString, ann, cellligandDF, cellreceptorDF, countligandDF, countreceptorDF, sumligandDF, sumreceptorDF, meanligandDF, meanreceptorDF, sumSpecifiedLigandDF, sumSpecifiedReceptorDF, meanSpecifiedLigandDF, meanSpecifiedReceptorDF, sourceFolder):
    origlabels = list(set(ann.iloc[:,0]))
    
    LRCountDict = {'Cluster':[],'Ligand count':[],'Receptor count':[]}
    
    writer = pd.ExcelWriter(os.path.join(sourceFolder, 'Ligands_Receptors_%s.xlsx' % (typeString.split('Edges_')[1])), engine='xlsxwriter')
    readmeDict = {'colA':['README'],'colB':['']}
    readmeDict['colA'].append('Ligand/receptor symbol')
    readmeDict['colB'].append('The official gene symbol of the detected ligand/receptor.')
    readmeDict['colA'].append('Total number of cells')
    readmeDict['colB'].append('Number of cells ligand/receptor is detected in.')
    readmeDict['colA'].append('Ligand/receptor detection rate')
    readmeDict['colB'].append('The ratio of cells that expressed the ligand/receptor to total cells in the cluster.')
    readmeDict['colA'].append('Ligand/receptor average expression value')
    readmeDict['colB'].append('The average expression level of the ligand/receptor in the cluster.')
    readmeDict['colA'].append('Ligand/receptor derived specificity of average expression value')
    readmeDict['colB'].append('The ratio of the average expression level of the ligand/receptor in the cluster to the sum of the average expression levels of the ligand/receptor in every cluster.')
    readmeDict['colA'].append('Ligand/receptor total expression value')
    readmeDict['colB'].append('The total expression level of the ligand/receptor in the cluster.')
    readmeDict['colA'].append('Ligand/receptor derived specificity of total expression value')
    readmeDict['colB'].append('The ratio of the total expression level of the ligand/receptor in the cluster to the sum of the total expression levels of the ligand/receptor in every cluster.')
    readmeDict['colA'].append('')
    readmeDict['colB'].append('')
    readmeDict['colA'].append('Summary of input file')
    readmeDict['colB'].append('')
    readmeDict['colA'].append('Cluster')
    readmeDict['colB'].append('Total number of cells')
    idmapDict = {}
    ididx = 0
    for origlabel in sorted(origlabels):
        ididx += 1
        idmapDict[origlabel] = 'Cluster %s' % ididx
        readmeDict['colA'].append(origlabel+'/'+idmapDict[origlabel])
        readmeDict['colB'].append(str(len(ann.index[ann['cluster'] == origlabel])))
    readmeDF = pd.DataFrame(readmeDict)
    readmeDF.to_excel(writer, sheet_name='README', index=False, header=False)
    
    for origlabel in sorted(origlabels):
        LRCountDict['Cluster'].append(origlabel)
        tempcellligandDF = cellligandDF.loc[:,[origlabel]]
        tempcountligandDF = countligandDF.loc[:,[origlabel]]
        tempsumligandDF = sumligandDF.loc[:,[origlabel]]
        tempmeanligandDF = meanligandDF.loc[:,[origlabel]]
        tempsumSpecifiedLigandDF = sumSpecifiedLigandDF.loc[:,[origlabel]]
        tempmeanSpecifiedLigandDF = meanSpecifiedLigandDF.loc[:,[origlabel]]
        tempLigandDF = pd.concat([tempcellligandDF, tempcountligandDF, tempmeanligandDF, tempmeanSpecifiedLigandDF, tempsumligandDF, tempsumSpecifiedLigandDF], axis=1)
        tempLigandDF = tempLigandDF.reset_index()
        tempLigandDF.columns = ['Ligand symbol', 'Total number of cells', 'Ligand detection rate', 'Ligand average expression value', 'Ligand derived specificity of average expression value', 'Ligand total expression value', 'Ligand derived specificity of total expression value']
        tempLigandDF['Total number of cells'] = tempLigandDF['Total number of cells'].astype(int)
        tempLigandDF = tempLigandDF.sort_values(by='Ligand detection rate', ascending=False)
        tempLigandDF = tempLigandDF.loc[tempLigandDF['Ligand total expression value']>0]
        LRCountDict['Ligand count'].append(len(tempLigandDF))
        tempcellreceptorDF = cellreceptorDF.loc[:,[origlabel]]
        tempcountreceptorDF = countreceptorDF.loc[:,[origlabel]]
        tempsumreceptorDF = sumreceptorDF.loc[:,[origlabel]]
        tempmeanreceptorDF = meanreceptorDF.loc[:,[origlabel]]
        tempmeanSpecifiedReceptorDF = meanSpecifiedReceptorDF.loc[:,[origlabel]]
        tempsumSpecifiedReceptorDF = sumSpecifiedReceptorDF.loc[:,[origlabel]]
        tempReceptorDF = pd.concat([tempcellreceptorDF, tempcountreceptorDF,tempmeanreceptorDF,tempmeanSpecifiedReceptorDF, tempsumreceptorDF,tempsumSpecifiedReceptorDF], axis=1)
        tempReceptorDF = tempReceptorDF.reset_index()
        tempReceptorDF.columns = ['Receptor symbol', 'Total number of cells', 'Receptor detection rate', 'Receptor average expression value', 'Receptor derived specificity of average expression value', 'Receptor total expression value', 'Receptor derived specificity of total expression value']
        tempReceptorDF['Total number of cells'] = tempReceptorDF['Total number of cells'].astype(int)
        tempReceptorDF = tempReceptorDF.sort_values(by='Receptor detection rate', ascending=False)
        tempReceptorDF = tempReceptorDF.loc[tempReceptorDF['Receptor total expression value']>0]
        LRCountDict['Receptor count'].append(len(tempReceptorDF))
        if len(tempLigandDF) > 0:
            tempLigandDF.to_excel(writer, sheet_name='Ligands in %s' % idmapDict[origlabel], index=False, header=True, columns = ['Ligand symbol', 'Total number of cells', 'Ligand detection rate', 'Ligand average expression value', 'Ligand derived specificity of average expression value', 'Ligand total expression value', 'Ligand derived specificity of total expression value'])
            worksheet = writer.sheets['Ligands in %s' % idmapDict[origlabel]]
            worksheet.conditional_format('C2:C%s' % (len(tempLigandDF)+1), {'type':'2_color_scale','min_color': '#FFFFFF','max_color': '#FF0000'})
            worksheet.conditional_format('D2:D%s' % (len(tempLigandDF)+1), {'type':'2_color_scale','min_color': '#FFFFFF','max_color': '#FF0000'})
            worksheet.conditional_format('E2:E%s' % (len(tempLigandDF)+1), {'type':'2_color_scale','min_color': '#FFFFFF','max_color': '#FF0000'})
            worksheet.conditional_format('F2:F%s' % (len(tempLigandDF)+1), {'type':'2_color_scale','min_color': '#FFFFFF','max_color': '#FF0000'})
            worksheet.conditional_format('G2:G%s' % (len(tempLigandDF)+1), {'type':'2_color_scale','min_color': '#FFFFFF','max_color': '#FF0000'})
            
        if len(tempReceptorDF) > 0:
            tempReceptorDF.to_excel(writer, sheet_name='Receptors in %s' % idmapDict[origlabel], index=False, header=True, columns = ['Receptor symbol', 'Total number of cells', 'Receptor detection rate', 'Receptor average expression value', 'Receptor derived specificity of average expression value', 'Receptor total expression value', 'Receptor derived specificity of total expression value'])
            worksheet = writer.sheets['Receptors in %s' % idmapDict[origlabel]]
            worksheet.conditional_format('C2:C%s' % (len(tempLigandDF)+1), {'type':'2_color_scale','min_color': '#FFFFFF','max_color': '#FF0000'})
            worksheet.conditional_format('D2:D%s' % (len(tempLigandDF)+1), {'type':'2_color_scale','min_color': '#FFFFFF','max_color': '#FF0000'})
            worksheet.conditional_format('E2:E%s' % (len(tempLigandDF)+1), {'type':'2_color_scale','min_color': '#FFFFFF','max_color': '#FF0000'})
            worksheet.conditional_format('F2:F%s' % (len(tempLigandDF)+1), {'type':'2_color_scale','min_color': '#FFFFFF','max_color': '#FF0000'})
            worksheet.conditional_format('G2:G%s' % (len(tempLigandDF)+1), {'type':'2_color_scale','min_color': '#FFFFFF','max_color': '#FF0000'})
    writer.save()
    

def FindCellsOfProtein(protein, proteinType, cellDF, countDF, sumDF, meanDF, sumSpecifiedDF, meanSpecifiedDF):
    subCellDF = cellDF.loc[cellDF.index == protein, :].T
    subCountDF = countDF.loc[countDF.index == protein, :].T
    subSumDF = sumDF.loc[sumDF.index == protein, :].T
    subMeanDF = meanDF.loc[meanDF.index == protein, :].T
    subSumSpecifiedDF = sumSpecifiedDF.loc[sumSpecifiedDF.index == protein, :].T
    subMeanSpecifiedDF = meanSpecifiedDF.loc[meanSpecifiedDF.index == protein, :].T
    mergedProteinDF = pd.concat([subCellDF, subCountDF, subSumDF, subMeanDF, subSumSpecifiedDF, subMeanSpecifiedDF], axis=1)
    newHeaders = ['cell '+proteinType, 'frequency '+proteinType, 'sum '+proteinType, 'mean '+proteinType, 'specified sum '+proteinType, 'specified mean '+proteinType]
    mergedProteinDF.columns = newHeaders
    
    return mergedProteinDF

def BuildHalfEdge(protein, cellDF, countDF, sumDF, meanDF, sumSpecifiedDF, meanSpecifiedDF, cellRole, proteinType):
    mergedProteinDF = FindCellsOfProtein(protein, proteinType, cellDF, countDF, sumDF, meanDF, sumSpecifiedDF, meanSpecifiedDF)
    
    #update headers
    mergedProteinDF[cellRole+' name'] = mergedProteinDF.index.values
    mergedProteinDF[proteinType] = protein
    mergedProteinDF.reset_index(drop=True, inplace=True)
    
    return mergedProteinDF

def GenSingleCell2CellEdge(ligandApprovedSymbolDict, receptorApprovedSymbolDict, resultDir, pair):
    ligand = pair[0]
    receptor = pair[1]
    if ligand not in ligandApprovedSymbolDict.keys() or receptor not in receptorApprovedSymbolDict.keys():
        return
    
    fn = os.path.join(resultDir,'%s-%s.xlsx' % (pair[0],pair[1]))
    
    edgeList = []
    ligandCellListDF = ligandApprovedSymbolDict[ligand]
    receptorCellListDF = receptorApprovedSymbolDict[receptor]
    for ligandCellListIndex in ligandCellListDF.index:
        for receptorCellListIndex in receptorCellListDF.index:
            ligandCellDF = ligandCellListDF.loc[ligandCellListIndex, :].to_frame().T.reset_index(drop=True)
            receptorCellDF = receptorCellListDF.loc[receptorCellListIndex, :].to_frame().T.reset_index(drop=True)
            mergedPairDF = pd.concat([ligandCellDF, receptorCellDF], axis=1)
            edgeList.append(mergedPairDF)
    
    edgeListDF = pd.concat(edgeList, axis=0)
    #add prodect columns
    edgeListDF['product of sum'] = edgeListDF['sum ligand'] * edgeListDF['sum receptor']
    edgeListDF['product of mean'] = edgeListDF['mean ligand'] * edgeListDF['mean receptor']
    edgeListDF['product of specified sum'] = edgeListDF['specified sum ligand'] * edgeListDF['specified sum receptor']
    edgeListDF['product of specified mean'] = edgeListDF['specified mean ligand'] * edgeListDF['specified mean receptor']
    
    #remove values equal zero
    edgeListDF = edgeListDF[edgeListDF['product of sum'] > 0]
    if len(edgeListDF) > 0:
        edgeListDF = edgeListDF.loc[:,['sending cluster name', 'ligand', 'receptor', 'target cluster name', 'cell ligand', 'frequency ligand', 'mean ligand', 'sum ligand', 'specified mean ligand', 'specified sum ligand', 
           'cell receptor', 'frequency receptor', 'mean receptor', 'sum receptor', 'specified mean receptor', 'specified sum receptor', 'product of mean', 'product of sum', 'product of specified mean', 'product of specified sum']]
        newCols = ['Sending cluster', 'Ligand symbol', 'Receptor symbol', 'Target cluster', 'Ligand-expressing cells', 'Ligand detection rate', 'Ligand average expression value', 'Ligand total expression value', 
            'Ligand derived specificity of average expression value', 'Ligand derived specificity of total expression value',  'Receptor-expressing cells', 'Receptor detection rate', 'Receptor average expression value', 'Receptor total expression value', 
            'Receptor derived specificity of average expression value', 'Receptor derived specificity of total expression value', 'Edge average expression weight', 'Edge total expression weight', 
            'Edge average expression derived specificity', 'Edge total expression derived specificity']
        edgeListDF.columns = newCols
        edgeListDF.to_excel(fn, header=True,index=False,columns=newCols)
    return

def GenCell2CellEdges(interDB, pairListDF, ligandApprovedSymbolDict, receptorApprovedSymbolDict, coreNum, resultDir):
    pairList = [tuple(x) for x in pairListDF.values]
    fn = os.path.join(resultDir,'LR-pairs_'+interDB)
    if not os.path.exists(fn):
        os.mkdir(fn)
    
    p = multiprocessing.Pool(coreNum)
    func = partial(GenSingleCell2CellEdge, ligandApprovedSymbolDict, receptorApprovedSymbolDict, fn)
    p.map(func, pairList)
    p.close()
    p.join()
    
    fl = []
    for p in pairList:
        f = os.path.join(fn,'%s-%s.xlsx' % (p[0],p[1]))
        if not os.path.exists(f):
            continue
        d = pd.read_excel(f, index_col=False, header=0)
        fl.append(d)
    edgeListDF = pd.concat(fl, axis=0)
    print('#### collected %s LR-mediated edges' % len(edgeListDF))
    
    return edgeListDF

def GenerateCell2CellTable(interDB, sumEMDF, meanEMDF, countEMDF, cellEMDF, typeString, pairsDF, ann, resultDir, coreNum):
    # generate ligand-receptor pair list
    pairListDF = GenLigandReceptorList(pairsDF)
    ligandApprovedSymbolList = list(pairsDF.index.values)
    receptorApprovedSymbolList = list(pairsDF.columns.values)
    
    print('#### extract signaling factors')
    # split original and speciefied countable to ligand and receptor file
    cellligandDF, cellreceptorDF, countligandDF, countreceptorDF, sumligandDF, sumreceptorDF, meanligandDF, meanreceptorDF, sumSpecifiedLigandDF, sumSpecifiedReceptorDF, meanSpecifiedLigandDF, meanSpecifiedReceptorDF = SplitIntoSinalProteins(sumEMDF, meanEMDF, countEMDF, cellEMDF, ligandApprovedSymbolList, receptorApprovedSymbolList)
    LRExpressions(typeString, ann, cellligandDF, cellreceptorDF, countligandDF, countreceptorDF, sumligandDF, sumreceptorDF, meanligandDF, meanreceptorDF, sumSpecifiedLigandDF, sumSpecifiedReceptorDF, meanSpecifiedLigandDF, meanSpecifiedReceptorDF, resultDir)
    
    # find cells for ligand and receptor
    ligandApprovedSymbolDict = {}
    for ligandApprovedSymbol in sumligandDF.index:#ligandApprovedSymbolList:
        mergedLigandDF = BuildHalfEdge(ligandApprovedSymbol, cellligandDF, countligandDF, sumligandDF, meanligandDF, sumSpecifiedLigandDF, meanSpecifiedLigandDF, 'sending cluster', 'ligand')
        ligandApprovedSymbolDict[ligandApprovedSymbol] = mergedLigandDF
    receptorApprovedSymbolDict = {}
    for receptorApprovedSymbol in sumreceptorDF.index:#receptorApprovedSymbolList:
        mergedReceptorDF = BuildHalfEdge(receptorApprovedSymbol, cellreceptorDF, countreceptorDF, sumreceptorDF, meanreceptorDF, sumSpecifiedReceptorDF, meanSpecifiedReceptorDF, 'target cluster', 'receptor')
        receptorApprovedSymbolDict[receptorApprovedSymbol] = mergedReceptorDF
    print('#### construct signaling interactions')
    
    # construct and save cell-to-cell edge
    edgeListDF = GenCell2CellEdges(interDB, pairListDF, ligandApprovedSymbolDict, receptorApprovedSymbolDict, coreNum, resultDir)
    
    if len(edgeListDF) > 0:
        return edgeListDF
    else:
        return []
        print('#### No signaling interactions found')

def GenerateDataFiles(interDB, sumEMDF, meanEMDF, countEMDF, cellEMDF, typeString, lrM, ann, resultDir, coreNum):
    # generate cell-ligand-receptor-cell table
    edgeListDF = GenerateCell2CellTable(interDB, sumEMDF, meanEMDF, countEMDF, cellEMDF, typeString, lrM, ann, resultDir, coreNum)
    
    #save result
    if len(edgeListDF) > 0:
        saveann = ann.reset_index()
        saveann.columns = ['cell', 'cluster']
        saveann.to_csv(os.path.join(resultDir,'ClusterMapping.csv'), index=False, header=True)
    
        edgeListFileDir = os.path.join(resultDir, typeString + '.csv')
        edgeListDF.to_csv(edgeListFileDir, index=False, header=True,columns=['Sending cluster', 'Ligand symbol', 'Receptor symbol', 'Target cluster', 'Ligand-expressing cells', 'Ligand detection rate', 'Ligand average expression value', 'Ligand total expression value', 
                'Ligand derived specificity of average expression value', 'Ligand derived specificity of total expression value',  'Receptor-expressing cells', 'Receptor detection rate', 'Receptor average expression value', 'Receptor total expression value', 
                'Receptor derived specificity of average expression value', 'Receptor derived specificity of total expression value', 'Edge average expression weight', 
                'Edge average expression derived specificity', 'Edge total expression weight', 'Edge total expression derived specificity'])
        with open(os.path.join(resultDir,'README.txt'), 'w') as file_object:
            file_object.write('README\n')
            file_object.write('\n')
            file_object.write('ClusterMapping.csv: cell-to-cluster mapping.\n')
            file_object.write('Ligands_Receptors_xx.xlsx: information about ligands and receptors in each cell-type/single-cell cluster.\n')
            file_object.write('Edges_xx.csv: all ligand-receptor-mediated communications.\n')
            file_object.write('LR-pairs_xx folder: all ligand-receptor-mediated communications via a ligand-receptor pair.\n')
            file_object.write('\n')
            file_object.write('Sending cluster: the cluster that expresses the ligand.\n')
            file_object.write('Target cluster: the cluster that expresses the receptor.\n')
            file_object.write('Ligand/receptor symbol: the official gene symbol of the detected ligand/receptor.\n')
            file_object.write('Ligand/receptor-expressing cells: number of cells ligand/receptor is detected in.\n')
            file_object.write('Ligand/receptor detection rate: the ratio of cells that expressed the ligand/receptor to total cells in the cluster.\n')
            file_object.write('Ligand/receptor average expression value: the average expression level of the ligand/receptor in the cluster.\n')
            file_object.write('Ligand/receptor derived specificity of average expression value: the ratio of the average expression level of the ligand/receptor in the cluster to the sum of the average expression levels of the ligand/receptor in every cluster.\n')
            file_object.write('Ligand/receptor total expression value: the total expression level of the ligand/receptor in the cluster.\n')
            file_object.write('Ligand/receptor derived specificity of total expression value: the ratio of the total expression level of the ligand/receptor in the cluster to the sum of the total expression levels of the ligand/receptor in every cluster.\n')
            file_object.write('Edge average expression weight: the product of average expression levels of the ligand and the receptor in the corresponding cluster(s).\n')
            file_object.write('Edge average expression derived specificity: the product of the ligand and receptor derived specificity of average expression values.\n')
            file_object.write('Edge total expression weight: the product of total expression levels of the ligand and the receptor in the corresponding cluster(s).\n')
            file_object.write('Edge total expression derived specificity: the product of the ligand and receptor derived specificity of total expression values.\n')
    else:
        with open(os.path.join(resultDir,'README.txt'), 'w') as file_object:
            file_object.write('README\n')
            file_object.write('\n')
            file_object.write('No signaling interactions found')

def main(species, emFile, annFile, idType, interDB, interSpecies, coreNum, outFolder):
    #load data
    emFileExtention = emFile.split('.')[-1]
    if emFileExtention == 'csv':
        em = pd.read_csv(emFile, index_col=0, header=0)
    elif emFileExtention == 'tsv' or emFileExtention == 'txt':
        em = pd.read_csv(emFile, index_col=0, header=0, sep='\t')
    elif emFileExtention == 'xls' or emFileExtention == 'xlsx':
        em = pd.read_excel(emFile, index_col=0, header=0)
    else:
        sys.exit("Cannot process the expression matrix file, please check the format of the expression matrix file, which can only be csv, tsv, txt, xls or xlsx.")

    em = em.loc[em.max(axis=1)>0,]
    
    #id conversion
    if idType != 'symbol' and idType != 'custom':
        if idType == 'hgnc' or idType == 'mgi':
            idType = 'ID'
        idT = pd.read_csv('ids/' + species + '_' + idType + '.csv', dtype=object, index_col=0, header=0)
        commonIDs = list(set(idT.index).intersection(set(em.index)))
        idT = idT.loc[commonIDs,]
        em = em.loc[commonIDs,]
        em.index = idT.loc[commonIDs,'Symbol']
        em = em.groupby(level=0).sum()
        
    if os.path.exists(opt.annFile):
        annFileExtention = annFile.split('.')[-1]
        if annFileExtention == 'csv':
            ann = pd.read_csv(annFile, index_col=0, header=0)
        elif annFileExtention == 'tsv' or emFileExtention == 'txt':
            ann = pd.read_csv(annFile, index_col=0, header=0, sep='\t')
        elif annFileExtention == 'xls' or emFileExtention == 'xlsx':
            ann = pd.read_excel(annFile, index_col=0, header=0)
        else:
            sys.exit("Cannot process the metafile, please check the format of the metafile, which can only be csv, tsv, txt, xls or xlsx.")
        
        #check annotation
        ann.columns = ['cluster']
        if len(set(em.columns).intersection(set(ann.index))) == 0:
            if len(set(em.columns).intersection(set(ann['cluster']))) == 0:
                sys.exit("Cannot find matched annotations for the barcodes.")
            elif len(set(em.columns).intersection(set(ann['cluster']))) < len(set(em.columns)):
                sys.exit("Cannot find matched annotations for all barcodes.")
            else:
                ann.columns = ['new.barcode']
                ann.reset_index().set_index('new.barcode')
                ann.columns = ['cluster']
        elif len(set(em.columns).intersection(set(ann.index))) < len(set(em.columns)):
                sys.exit("Cannot find matched annotations for all barcodes.")
    else:
        ann = pd.DataFrame({'cluster':list(em.columns)},index=em.columns)
            
    #load interaction list
    interDBExtention = interDB.split('.')[-1]
    if interDBExtention == 'csv':
        lrL = pd.read_csv(os.path.join('lrdbs', interDB), dtype=object, index_col=None, header=0)
    elif interDBExtention == 'tsv' or interDBExtention == 'txt':
        lrL = pd.read_csv(os.path.join('lrdbs', interDB), dtype=object, index_col=None, header=0, sep='\t')
    elif interDBExtention == 'xls' or interDBExtention == 'xlsx':
        lrL = pd.read_excel(os.path.join('lrdbs', interDB), dtype=object, index_col=None, header=0)
    else:
        sys.exit("Cannot process the ligand-receptor interaction database file, please check the format of the file, which can only be csv, tsv, txt, xls or xlsx.")
    
    if interSpecies == 'p':
        lrE = pd.read_csv(os.path.join('lrdbs', 'lrc2p.csv'), dtype=object, index_col=None, header=0)
    elif interSpecies == 'a':
        lrE = pd.read_csv(os.path.join('lrdbs', 'lrc2a.txt'), dtype=object, index_col=None, header=0, sep='\t')
    if interSpecies == 'p' or interSpecies == 'a':
        lrL = pd.concat([lrL,lrE],ignore_index=True)
        
    #to adj matrix
    lset = sorted(list(set(lrL['Ligand gene symbol'])))
    rset = sorted(list(set(lrL['Receptor gene symbol'])))
    lrM = pd.DataFrame(0,index=lset,columns=rset)
    for idx in lrL.index:
        lrM.loc[lrL.loc[idx,'Ligand gene symbol'], lrL.loc[idx,'Receptor gene symbol']] = 1
    interDB = interDB.split('.')[0]
    
    # change gene symbols if necessary
    if interSpecies in ['a','p']:
        interSpeciesType = '9606'
    else:
        interSpeciesType = interSpecies
    if idType != 'custom' and species != interSpeciesType:
        #find taxonomy file
        ## HID (HomoloGene group id)[0] - Taxonomy ID[1] - Gene ID[2] - Gene Symbol[3] - Protein gi[4] - Protein accession[5]
        homoMapDir = 'homology/homologene.data'
        hidCol = 0
        taxidCol = 1
        geneSymbolCol = 3
        lrM = TransferToGeneSymbol(homoMapDir, species, interSpeciesType, taxidCol, geneSymbolCol, hidCol, lrM)
        
    #build the folder to save the analysis results
    if outFolder != '':
        resultDir = os.path.abspath(outFolder)
    else:
        resultDir = emFile.split('.'+emFileExtention)[0]
    if not os.path.exists(resultDir):
        os.mkdir(resultDir)
    print('#### the folder "%s" has been created to save the analysis results' % os.path.abspath(resultDir))
    
    #cluster expression matrix
    print('#### cluster the expression matrix')
    sumEMDF, meanEMDF, countEMDF, cellEMDF = ClusterAnnotateEM(resultDir, em, ann)
    
    # generate cell-to-cell interaction files
    print('#### generate cell-to-cell interactions')
    GenerateDataFiles(interDB, sumEMDF, meanEMDF, countEMDF, cellEMDF, 'Edges_%s' % (os.path.basename(interDB)), lrM, ann, resultDir, coreNum)
    
    print('#### DONE')
    
#===================================================================#
#-------------------------------START-------------------------------#
#===================================================================#

if __name__ == '__main__':
    #process arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('--interDB', default='lrc2p', help='lrc2p (default) | lrc2a | name of the user-supplied interaction database file without extension')
    parser.add_argument('--interSpecies', default='human', help='human (default) | mouse | expandp | expanda')
    parser.add_argument('--emFile', required=True, help='the path to the file of the expression matrix with row names (gene identifiers) and column names (single-cell/cell-type identifiers)')
    parser.add_argument('--annFile', default='', help='the path to the metafile in which column one has single-cell identifiers and column two has corresponding cluster IDs (see file "toy.sc.ann.txt" as an example). This file is NOT required for bulk data')
    parser.add_argument('--species', default='human', help='human (default) | mouse | rat | zebrafish | fruitfly | chimpanzee | dog | monkey | cattle | chicken | frog | mosquito | nematode | thalecress | rice | riceblastfungus | bakeryeast | neurosporacrassa | fissionyeast | eremotheciumgossypii | kluyveromyceslactis') 
    parser.add_argument('--idType', default='symbol', help='symbol (default) | entrez | ensembl | uniprot | hgnc | mgi | custom (gene identifier used in the expression matrix)')
    parser.add_argument('--coreNum', type=int, default=1, help='the number of CPU cores used, default is one')
    parser.add_argument('--out', default='', help='the path to save the analysis results')
    
    opt = parser.parse_args()
    
    #check species
    avaSpecDict = {'human':'9606', 'mouse':'10090', 'chimpanzee':'9598', 'dog':'9615', 'monkey':'9544', 'cattle':'9913', 'rat':'10116', 'chicken':'9031', 'frog':'8364', 'zebrafish':'7955', 'fruitfly':'7227', 'mosquito':'7165', 'nematode':'6239', 'thalecress':'3702', 'rice':'4530', 'riceblastfungus':'318829', 'bakeryeast':'4932', 'neurosporacrassa':'5141', 'fissionyeast':'4896', 'eremotheciumgossypii':'33169', 'kluyveromyceslactis':'28985'}
    if opt.species.lower() not in avaSpecDict.keys() and opt.idType.lower() != 'custom':
        sys.exit("The species can only be 'human' or other 20 supported species for now.")
    else:
        species = avaSpecDict[opt.species.lower()]
    
    #check gene ids
    avaIDList = ['symbol', 'entrez', 'ensembl', 'uniprot', 'hgnc', 'mgi', 'custom']
    if opt.idType.lower() not in avaIDList:
        sys.exit("The gene identifiers can only be 'symbol', 'entrez', 'ensembl', 'uniprot', 'hgnc', or 'mgi' for now.")
    elif species == '9606' and opt.idType.lower() == 'mgi':
        sys.exit("Human expression data cannot use MGI IDs.")
    elif species == '10090' and opt.idType.lower() == 'hgnc':
        sys.exit("Mouse expression data cannot use HGNC IDs.")
    else:
        idType = opt.idType.lower()
    
    #check interDB
    stflist = glob.glob('lrdbs/%s.*' % opt.interDB)
    if len(stflist) == 0:
        sys.exit("Cannot find the ligand-receptor interaction database named '%s'." % opt.interDB)
    elif len(stflist) > 1:
        sys.exit("There are multiple ligand-receptor interaction databases with the same name '%s'." % opt.interDB)
    else:
        interDB = os.path.basename(stflist[0])
    
    #check interSpecies
    avaInterSpecDict = {'human':'9606', 'mouse':'10090', 'expandp':'p', 'expanda':'a'}
    if opt.interSpecies.lower() not in avaInterSpecDict.keys():
        sys.exit("'interSpecies' can only be 'human', 'mouse', 'expandp' or 'expanda' for now.")
    else:
        interSpecies = avaInterSpecDict[opt.interSpecies.lower()]
        
    #check emFile
    if not os.path.exists(opt.emFile):
        sys.exit("The expression matrix file does not exist.")
        
    #check annFile
    if not os.path.exists(opt.annFile):
        print("The metafile does not exist, each sample will be considered as a cluster.")
        
    #check coreNum
    maxCoreNum = multiprocessing.cpu_count()
    if opt.coreNum > maxCoreNum:
        sys.exit("There are only %s cores availble, less than %s cores." % (maxCoreNum, opt.coreNum))
        
    #pass argument check, show input data
    print('===================================================')
    print('Input data:')
    if idType != 'custom':
        print('Species to obtain abundance data: %s' % opt.species.lower())
    print('The expression matrix file: %s' % opt.emFile)
    if os.path.exists(opt.annFile):
        print('The metafile: %s' % opt.annFile)
    else:
        print('The metafile is not provided')
    print('The ligand-receptor interaction database: %s' % interDB)
    if opt.interSpecies.lower() in ['human', 'mouse']:
        print('Species to obtain ligand-receptor interactions: %s' % opt.interSpecies.lower())
    else:
        print('Ligand-receptor interactions are used to expand: lrc2%s' % interSpecies)
    print('The number of cores to use: %s' % opt.coreNum)
    if opt.out != '':
        print('The folder to save all results: %s' % os.path.abspath(opt.out))
    print('===================================================')
    
    #start to construct cell-to-cell network
    main(species, opt.emFile, opt.annFile, idType, interDB, interSpecies, opt.coreNum, opt.out)
