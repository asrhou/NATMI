# NATMI: A network analysis toolkit for multicellular interactions in single-cell transcriptome data

Recent development of high throughput single-cell sequencing technologies has made it cost-effective to profile thousands of cells from a complex sample. Examining ligand and receptor expression patterns in the cell types identified from these datasets allows prediction of cell-to-cell communication at the level of niches, tissues and organism-wide. Here, we developed NATMI (Network Analysis Toolkit for Multicellular Interactions), a Python-based toolkit for multi-cellular communication network construction and network analysis of multispecies single-cell and bulk gene expression and proteomic data. NATMI uses connectomeDB2020 (the most up-to-date manually curated ligand-receptor interaction list) and user supplied gene expression tables with cell type labels to predict and visualize cell-to-cell communication networks. By interrogating the Tabula Muris cell atlas we demonstrate the utility of NATMI to identify cellular communities and the key ligands and receptors involved. Notably, we confirm our previous predictions from bulk data that autocrine signalling is a major feature of cell-to-cell communication networks and for the first time ever show a substantial potential for self-signalling of individual cells through hundreds of co-expressed ligand-receptor pairs. Lastly, we identify age related changes in intercellular communication between the mammary gland of 3 and 18-month-old mice in the Tabula Muris dataset. NATMI and our updated ligand-receptor lists are freely available to the research community.

NATMI is maintained by Rui Hou [rui.hou@research.uwa.edu.au]

- [Download and Installation](#download-and-installation)
- [Required Data and Formats](#required-data-and-formats)
  * [Supported Species](#supported-species)
  * [Expression Data](#expression-data)
  * [Ligand-Receptor Interactions (connectomeDB2020 or user-specified)](#ligand-receptor-interactions-connectomeDB2020-or-user-supplied)
  * [Metafile (single cell analysis only)](#metafile-single-cell-analysis-only)
- [Command Line Utilities](#command-line-utilities)
  * [ExtractEdges.py](#extractedges-extracting-ligand-receptor-mediated-interactions-between-cell-types-in-the-input-transcriptome-data)
  * [DiffEdges.py](#diffedges-identification-of-changes-in-ligand-receptor-edge-weights-between-a-cell-type-pair-in-two-conditions)
  * [VisInteractions.py](#visinteractionspy-visualisation-of-the-network-analysis-results-from-extractedgespy-and-diffedgespy)
- [Example Workflow Simple (single-cell toy dataset)](#example-workflow-simple-single-cell-toy-dataset)
  * [Extract ligand-receptor-mediated interactions](#extract-ligand-receptor-mediated-interactions-in-toyscemtxt-and-save-results-to-test-folder-using-extractedgespy)
  * [Visualise cell-to-cell communication networks](#visualise-ligand-receptor-mediated-interaction-network-of-in-toyscemtxt-in-three-different-ways)
- [Example Workflow Advanced (Tabula Muris Senis dataset)](#example-workflow-advanced-tabula-muris-senis-dataset)
  * [Extract ligand-receptor-mediated interactions at two time-points.](#we-firstly-extract-edges-between-cells-of-the-3-and-18-month-old-mammary-glands-in-mice-using-extractedgespy)
  * [Identify variations in cell-to-cell signaling networks](#the-variations-in-cell-to-cell-signaling-between-3-month-old-and-18-month-old-murine-mammary-gland-are-then-identified-by-diffedgespy)
  * [Visualize the cell-to-cell communication networks (Figure 6 of the manuscript)](#we-visualize-up--and-down-regulated-edges-between-3-months-and-18-months-using-visinteractionspy-as-in-figure-6-of-the-manuscript)

## Download and Installation

In order to use NATMI, please make sure you match the following pre-requisites:

- Python 2.X or Python3.X

- Python libraries [seaborn](https://seaborn.pydata.org/), [igraph](https://igraph.org/python/), [NetworkX](https://networkx.github.io/) and [PyGraphviz](https://pygraphviz.github.io/) are required to visualise the cell-to-cell communication network at three distinct levels. 

NATMI was tested using python 2.7 and 3.7 versions and seaborn 0.8.1, igraph 0.7.1, NetworkX 2.1 and PyGraphviz 1.5 versions.

From the directory in which you wish to install NATMI, run the following command:
```bat
   git clone https://github.com/asrhou/NATMI.git
```

This tool currently provides command line utilities only.

## Required Data and Formats

To explore cell-cell to cell communication NATMI uses (1) user-supplied gene/protein abundance data, (2) ligand-receptor interactions (precomplied connectomeDB2020 or user-supplied interactions) and for the single cell data analysis it requires (3) the metafile describing mapping between each cell in the dataset and a cell-type label. By deaful NATMI uses official gene symbols from human and mouse, but it can also automatically convert HGNC IDs, MGI IDs, Entrez gene IDs, Ensembl gene IDs, and UniProt IDs using biomartXXX [short technical details here]. 

[We need to say what happens for the converted IDs -> will the output be in original IDs or in gene symbols?]


### Supported Species

Currently, NATMI can process expression data from human and mouse. The species of input data is specified by the argument '--species' in ExtractEdges.py.


[We might need to expand the above by a sentence or two. This is what I see in the methods: 

"The curated ligand-receptor database provided in this manuscript is based on interactions between human ligands and receptors. To run NATMI on other species we extracted the homologs of interacting pairs from NCBI HomoloGene Database [22]. The homologous pairs for 20 commonly requested species such as mouse, rat, zebrafish (https://www.ncbi.nlm.nih.gov/homologene/statistics/) can be automatically inferred by NATMI based on the input gene expression data. Note there are reported ligand-receptor pairs that are species specific, thus when applying NATMI to other species, please always check the literature to verify if a reported edge is only seen in humans."

So what happens if one let's say trie that on a zebrafish data?]

### Expression Data 

User-specified gene/protein abundance matrix files are supported in the following formats: csv, tsv, txt, xls or xlsx and require the gene/protein IDs to be one of the following: human and mouse official gene symbols (deafult) or HGNC IDs, MGI IDs, Entrez gene IDs, Ensembl gene IDs, and UniProt IDs (see [required data and formats](#required-data-and-formats)). 
Additionally, [Tabula Muris](https://tabula-muris.ds.czbiohub.org/), [Tabula Muris Senis](https://tabula-muris-senis.ds.czbiohub.org/) and [FANTOM5 cell atlas](http://fantom.gsc.riken.jp/5/suppl/Ramilowski_et_al_2015/) can also be explored. 

### Ligand-Receptor Interactions (connectomeDB2020 or user-supplied interactions)

In 2015, we publised a first draft of human cell interactions and a database of human ligand-receptor pairs ([Ramilowski, J. A., et al.  Nat Commun 6, 7866 (2015)](https://www.nature.com/articles/ncomms8866)). This database compiled 708 ligands and 691 receptors into 2,422 human ligand-receptor interacting pairs (1,894 pairs with primary literature support, and an additional 528 putative pairs with high-throughput protein-protein interaction evidence). In 2020, we made an updated and expanded database of 2,187 human ligand-receptor pairs with primary litarture support and additional 1,791 putative pairs and named it **connectomeDB2020** [\reference to the paper]. By default, ExtractEdges.py of NATMI extracts edges from input expression data based on the literature-supported ligand-receptor pairs from connectomeDB2020.

Alternatively, since ExtractEdges.py can also work with **user-supplied ligand-receptor interactions** (argument '--signalType'), we briefly describe requirements and the format of this dataset. Similar to the pre-compiled in the repository connectomeDB2020 datasets, an interaction data file 'pairsM.xlsx' must be stored in a folder at the same location as ExtractEdges.py (the script will search for the file named 'pairsM.xlsx' based on the folder name). More, the data file should be in a binary matrix form with row names denoting ligands and column names receptors (all represented by the appropriate human gene symbol). For an interacting ligand-receptor pair, the corresponding matrix element is 1, otherwise, it is 0. Following table provides a toy example of the correct interaction matrix.

'pairsM.xlsx' --> I wonder if we can think how to make it readible in other formats too. Seems like a punishment to ask for excel :)

    No, ExtractEdges.py will only search for 'pairsM.xlsx' only.

connectomeDB2020 --> it starts feeling like the name is real close to ConnectomeDB (we knew that...)

    In this FANTOM5 paper [Update of the FANTOM web resource: high resolution transcriptome of diverse cell types in mammals](https://academic.oup.com/nar/article/45/D1/D737/2333885), connectomeDB2015 was referred as 'Ligand Receptor Connectome'.

                        
||Receptor1|Receptor2|Receptor3|...|
|-:|:-:|:-:|:-:|:-|
|**Ligand1**|1|0|0|...|
|**Ligand2**|0|0|1|...|
|**Ligand3**|1|1|0|...|
|**...**|...|...|...|...|

### Metafile (single cell analysis only)

For single-cell gene expression data, the user needs to provide a metafile with the mapping between each cell in the dataset and a cell-type label. Following table displays the format of a metafile.

|Cell|Annotation|
|-:|:-|
|Barcode1|Cell-type1|
|Barcode2|Cell-type1|
|Barcode3|Cell-type2|
|...|...|

## Command Line Utilities

NATMI is a python-based tool (see [software requirements](#software-requirements)) to construct cell-to-cell ligand-receptor communication networks from multiomics data. It works with user-specified gene/protein abundance matrix files or can be used to explore Tabula Muris, Tabula Muris Senis] and FANTOM5 cell atlas (see [required data](#expression-data)). 

Path to NATMI scripts has to be specified or the following commands can be run from the instalation directory:

### ExtractEdges: Extracting ligand-receptor-mediated interactions between cell types in the input transcriptome data.
[transcriptome data--> once we agree on how we will word this in the paper, we can modify this]

*Optional arguments are enclosed in square brackets […]*

```
python ExtractEdges.py [-h] [--species SPECIES] --emFile EMFILE [--idType IDTYPE] [--annFile ANNFILE] [--signalType SIGNALTYPE] [--coreNum CORENUM] [--out OUT]

Arguments:
  -h, --help            show this help message and exit
  --species SPECIES     only human and mouse expression data are currently supported, default is "human"
  --emFile EMFILE       the path to the expression matrix file with row names (gene symbols) and column names (cell names/single-cell identifiers)
  --annFile ANNFILE     the path to the metafile in which column one has single-cell identifiers and column two has corresponding cluster IDs (see file 'toy.sc.ann.txt' as an example)
  --signalType SIGNALTYPE
                        lrc2p (default) has literature supported ligand-receptor pairs | lrc2a has putative and literature supported ligand-receptor pairs, folder name of the interaction database
  --idType IDTYPE       symbol (default) | entrez(https://www.ncbi.nlm.nih.gov/gene) | ensembl(https://www.ensembl.org/) | uniprot(https://www.uniprot.org/) | hgnc(https://www.genenames.org/) | mgi(http://www.informatics.jax.org/mgihome/nomen/index.shtml), gene identifier used in the expression matrix
  --coreNum CORENUM     the number of CPU cores used, default is one
  --out OUT             the path to save the analysis results [is this optonal or required?]
```

**Note**: Expression matrix and metafile are supported in csv, tsv, txt, xls or xlsx format.

Predict ligand-receptor-mediated interactions in a mouse single-cell RNA-seq dataset using literature supported ligand-receptor pairs and four CPUs:
```bat
   python ExtractEdges.py --species mouse --emFile toy.sc.em.txt --annFile toy.sc.ann.txt --signalType lrc2p --coreNum 4
```

Predict ligand-receptor-mediated interactions in a human bulk RNA-seq dataset using putative and literature supported ligand-receptor pairs and one CPU:
```bat
   python ExtractEdges.py --species human --emFile toy.bulk.em.xls --signalType lrc2a
```

ExtractEdges.py creates a folder using the name of the expression matrix or the user specified name. README.txt (in the output folder) contains information about other files in the folder.

### DiffEdges: Identification of changes in ligand-receptor edge weights between a cell-type pair in two conditions. 
*Optional arguments are enclosed in square brackets […]*

**Note**: Only weight changes across two condition from the same, or similar, datasets and with the same 'signalType' of ligand-receptor pairs (literature-supported with literature-supported or all with all) should be compared.

```
python DiffEdges.py [-h] --refFolder REFFOLDER --targetFolder TARGETFOLDER [--signalType SIGNALTYPE] [--weightType WEIGHTTYPE] [--out OUT]

Arguments:
  -h, --help            show this help message and exit
  --refFolder REFFOLDER
                        the path to the folder of the reference dataset
  --targetFolder TARGETFOLDER
                        the path to the folder of the target dataset
  --signalType SIGNALTYPE
                        lrc2p (default) | lrc2a, folder name of the interaction database
  --weightType WEIGHTTYPE
                        mean (default) | sum
  --out OUT
                        the path to save the analysis results
```

Detect changes in edge weight in two output folders (generated by ExtractEdges.py) using literature supported ligand-receptor pairs:
```bat
   python DiffEdges.py --refFolder /path/to/ExtractEdges.py's/output/folder/of/reference/dataset --targetFolder /path/to/ExtractEdges.py's/output/folder/of/target/dataset --signalType lrc2p
```

DiffEdges.py creates a folder from the names of the two datasets or the user specified name. README.txt (in the output folder) contains information about other files in the folder.

### VisInteractions.py: Visualisation of the network analysis results from ExtractEdges.py and DiffEdges.py.
*Optional arguments are enclosed in square brackets […]*

```
python VisInteractions.py [-h] --sourceFolder SOURCEFOLDER [--signalType SIGNALTYPE] [--weightType WEIGHTTYPE] [--specificityThreshold SPECIFICITYTHRESHOLD]
                          [--expressionThreshold EXPRESSIONTHRESHOLD] [--detectionThreshold DETECTIONTHRESHOLD]
                          [--keepTopEdge KEEPTOPEDGE] [--plotWidth PLOTWIDTH] [--plotHeight PLOTHEIGHT] [--plotFormat PLOTFORMAT]
                          [--edgeWidth EDGEWIDTH] [--clusterDistance CLUSTERDISTANCE] [--drawClusterPair DRAWCLUSTERPAIR]
                          [--layout LAYOUT] [--fontSize FONTSIZE] [--maxClusterSize MAXCLUSTERSIZE]
                          [--drawNetwork DRAWNETWORK] [--drawLRNetwork [LIGAND [RECEPTOR ...]]]

Arguments:
  -h, --help            show this help message and exit
  --sourceFolder SOURCEFOLDER
                        the path to the folder of extracted edges from ExtractEdges.py or DiffEdges.py
  --signalType SIGNALTYPE
                        lrc2p (default) | lrc2a, folder name of the interaction database
  --weightType WEIGHTTYPE
                        mean (default) | sum, method to calculate the expression level of a ligand/receptor in a cell type
  --specificityThreshold SPECIFICITYTHRESHOLD
                        do not draw the edges whose specicificities are not greater than the threshold (default 0).
  --expressionThreshold EXPRESSIONTHRESHOLD
                        do not draw the edges in which expression levels of the ligand and the receptor are not greater than the threshold (default 0).
  --detectionThreshold DETECTIONTHRESHOLD
                        do not draw the interactions in which detection rates of the ligand and the receptor are lower than the threshold (default 0.2).
  --keepTopEdge KEEPTOPEDGE
                        only draw top n interactions that passed the thresholds (default 0 means all interactions that passed the thresholds will be drawn).
  --plotWidth PLOTWIDTH
                        resulting plot's width (default 12).
  --plotHeight PLOTHEIGHT
                        resulting plot's height (default 10).
  --plotFormat PLOTFORMAT
                        pdf (default) | png | svg, format of the resulting plot(s)
  --edgeWidth EDGEWIDTH
                        maximum thickness of edges in the plot (default 0: edge weight is shown as a label around the edge).
  --clusterDistance CLUSTERDISTANCE
                        relative distance between clusters (default value is 1; if clusterDistance >1, the distance will be increased, if clusterDistance >0 and clusterDistance <1, the distance will be decreased).
  --drawClusterPair DRAWCLUSTERPAIR
                        n(o) (default) | y(es)
  --layout LAYOUT       kk (default) | circle | random | sphere; 'kk' stands for Kamada-Kawai force-directed algorithm
  --fontSize FONTSIZE   font size for node labels (default 8).
  --maxClusterSize MAXCLUSTERSIZE
                        maximum radius of the clusters (default 0: all clusters have identical radius).
  --drawNetwork DRAWNETWORK
                        y(es) (default) | n(o)
  --drawLRNetwork [DRAWLRNETWORK [DRAWLRNETWORK ...]]
                        ligand and receptor's symbols
```


Visualise cell-connectivity-summary networks from the results of ExtractEdges.py and DiffEdges.py:
```bat
   python VisInteractions.py --sourceFolder /path/to/result/folder --signalType lrc2p --weightType mean --detectionThreshold 0.2 --plotFormat pdf --drawNetwork y --plotWidth 12 --plotHeight 10 --layout kk --fontSize 8 --edgeWidth 0 --maxClusterSize 0 --clusterDistance 1
```

ExtractEdges.py creates a folder (in the result folder) containing networks with three different weights. DiffEdges.py creates a folder (in the result folder), containing networks with three different weights in reference and target datasets. Additionally, delta networks are drawn, where yellow edges are (non-significant) edges with the fold change of their weights in two conditions of two or less. For other edges, a red color indicates the edges with a weight higher in the reference dataset, and a green color indicates the edges with a weight higher in the target dataset. The color intensity scales with the degree of change.

Visualise cell-to-cell communication networks between two cell types using results of ExtractEdges.py and DiffEdges.py:

```bat
   python VisInteractions.py --sourceFolder /path/to/result/folder --signalType lrc2p --drawClusterPair y
```

ExtractEdges.py creates a folder (in the result folder) containing bipartite graphs with three different weights. DiffEdges.py creates a folder (in the result folder) containing four kinds of interactions. From a cell type to another cell type, each kind of interactions form a separate bipartite graph.

Visualise cell-to-cell communication networks via a ligand-receptor pair from the results of ExtractEdges.py:

```bat
   python VisInteractions.py --sourceFolder /path/to/result/folder --signalType lrc2p --drawLRNetwork LIGAND.SYMBOL RECEPTOR.SYMBOL
```

DiffEdges.py creates a folder (in the result folder) containing the simple graph and hypergraph for the given ligand-receptor pair in the dataset. 

## Example Workflow Simple (single-cell toy dataset)
This workflow shows how to extract and visualize intercellular communication using mouse single-cell RNA-seq dataset ('toy.sc.em.txt') and the corresponding annotation file ('toy.sc.ann.txt') and literature supported ligand-receptor pairs from connectomeDB2020.

### Extract ligand-receptor-mediated interactions in 'toy.sc.em.txt' and save results to 'test' folder using ExtractEdges.py. 
The first step of NATMI-based analysis is always to predict the potential ligand-receptor-mediated interactions between cells using the user-specified ligand-receptor pairs. Here, we use literature supported ligand-receptor pairs in connectomeDB2020 and ExtractEdges.py to extract interactions in the toy single-cell dataset.

```bat
   python ExtractEdges.py --species mouse --emFile toy.sc.em.txt --annFile toy.sc.ann.txt --signalType lrc2p --coreNum 4 --out test
```

### Visualise ligand-receptor-mediated interaction network of in 'toy.sc.em.txt' in three different ways. 
The output of ExtractEdges.py in 'test' folder are the predicted edges between three cell types. Visualisation of the extracted edges is a good place to start interrogating biological processes through these predicted edges. To have a complete view of the cell-to-cell communication network, we first visualise the cell-connectivity-summary network in 'test' folder.

```bat
   python VisInteractions.py --sourceFolder test --signalType lrc2p --weightType mean --detectionThreshold 0.2 --drawNetwork y --plotWidth 4 --plotHeight 4 --layout circle --fontSize 15 --edgeWidth 6 --maxClusterSize 0 --clusterDistance 0.6
```

Network in *test/Network_exp_0_spe_0_det_0.2_top_0_signal_lrc2p_weight_mean/network_total-specificity-based_layout_circle.pdf* is the cell-connectivity-summary network weighted by the sum of specificity weights for each ligand-receptor pair. Apparently, there are more specific ligand-receptor pairs from endothelial cell to itself. 

We then visualise top ten ligand-receptor pairs between three cell types.

```bat
   python VisInteractions.py --sourceFolder test --drawClusterPair y --keepTopEdge 15
```

Bipartite graph in *test/CltPair_exp_0_spe_0_det_0.2_top_15_signal_lrc2p_weight_mean/From_Endothelial Cells_to_Endothelial Cells_spe.pdf* shows fifteen most specific ligand-receptor pairs endothelial cell to itself. Based on the edge weights (specificities), eleven of them are only detected in endothelial cells. Since the specificity of Efnb2-Pecam1 pair in the graph is less than one, it is interesting to know which cell-type pairs are also connected by Efnb2-Pecam1 pair.

We then visualise the cell-to-cell communication network via Efnb2-Pecam1 pair.

```bat
   python VisInteractions.py --sourceFolder test --drawLRNetwork Efnb2 Pecam1 --plotWidth 4 --plotHeight 4 --layout circle --fontSize 15 --edgeWidth 6 --maxClusterSize 0 --clusterDistance 0.6
```

Network in *test/LRNetwork_Efnb2-Pecam1_exp_0_spe_0_det_0.2_top_0_signal_lrc2p_weight_mean/network_Efnb2-Pecam1_layout_circle.pdf* only has one edge. This means although other cell-type pairs are connected by edges of Efnb2-Pecam1 pair, only for endothelial cell, Efnb2 and Pecam1 are detected in > 20 % cells. Therefore, Efnb2-Pecam1 pair is only reliably detected in endothelial cell.

## Example Workflow Advanced (Tabula Muris Senis dataset)

To demonstrate the usage of delta network analysis, we show the analysis on Tabula Muris Senis as in our manuscript. Processed Tabula Muris Senis data 'Mammary_Gland_droplet.h5ad' was downloaded from figshare (https://figshare.com/projects/Tabula_Muris_Senis/64982). We then extracted 3 and 18-month-old mammary gland cells and normalized each expression profile by total number of unique molecular identifiers and then rescaled by multiplying by 1,000,000. Normalized gene expression data and annotations are available in figshare: https://figshare.com/s/7f45bf6352da453b3266.

### We firstly extract edges between cells of the 3 and 18-month-old mammary glands in mice using ExtractEdges.py. 
```bat
   python ExtractEdges.py --species mouse --emFile /path/to/3m.upm.em.csv --annFile /path/to/3m.ann.csv --signalType lrc2p --coreNum 4 --out 3m.mg

   python ExtractEdges.py --species mouse --emFile /path/to/18m.upm.em.csv --annFile /path/to/18m.ann.csv --signalType lrc2p --coreNum 4 --out 18m.mg
```

### The variations in cell-to-cell signaling between 3-month-old and 18-month-old murine mammary gland are then identified by DiffEdges.py.
```bat
   python DiffEdges.py --refFolder 3m.mg --targetFolder 18m.mg --signalType lrc2p --out 3m-18m
```

### We visualize up- and down-regulated edges between 3 months and 18 months using VisInteractions.py as in Figure 6 of the manuscript.

```bat
   python VisInteractions.py --sourceFolder 3m-18m --signalType lrc2p --weightType mean --detectionThreshold 0.2 --drawNetwork y --plotWidth 10 --plotHeight 10 --layout circle --fontSize 15 --edgeWidth 6 --maxClusterSize 0 --clusterDistance 0.6
```

*Resulting networks are in the folder '/path/to/3m-18m/Delta_Network_exp_0_spe_0_det_0.2_top_0_signal_lrc2p_weight_mean'*
