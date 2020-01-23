# NATMI: A network analysis toolkit for multicellular interactions in single-cell transcriptome data

Recent development of high throughput single-cell sequencing technologies has made it cost-effective to profile thousands of cells from a complex sample. Examining ligand and receptor expression patterns in the cell types identified from these datasets allows prediction of cell-to-cell communication at the level of niches, tissues and organism-wide. Here, we developed NATMI (Network Analysis Toolkit for Multicellular Interactions), a Python-based toolkit for multi-cellular communication network construction and network analysis of multispecies single-cell and bulk gene expression and proteomic data. NATMI uses connectomeDB2020 (the most up-to-date manually curated ligand-receptor interaction list) and user supplied gene expression tables with cell type labels to predict and visualize cell-to-cell communication networks. By interrogating the Tabula Muris cell atlas we demonstrate the utility of NATMI to identify cellular communities and the key ligands and receptors involved. Notably, we confirm our previous predictions from bulk data that autocrine signalling is a major feature of cell-to-cell communication networks and for the first time ever show a substantial potential for self-signalling of individual cells through hundreds of co-expressed ligand-receptor pairs. Lastly, we identify age related changes in intercellular communication between the mammary gland of 3 and 18-month-old mice in the Tabula Muris dataset. NATMI and our updated ligand-receptor lists are freely available to the research community.

NATMI is maintained by Rui Hou [rui.hou@research.uwa.edu.au]

## Download and Installation
```bat
   git clone https://github.com/asrhou/xxx.git
```

This tool currently provides command line utilities only.

## Command Line Utilities

NATMI is a python-based (both python 2 and 3 are supported) tool to construct cell-to-cell ligand-receptor communication networks from multiomics data. It works with with user-specified gene/protein abundance matrix files (csv, tsv, txt, xls or xlsx format) or can be used to explore [Tabula Muris](https://tabula-muris.ds.czbiohub.org/), [Tabula Muris Senis](https://tabula-muris-senis.ds.czbiohub.org/) and [FANTOM5 cell atlas](http://fantom.gsc.riken.jp/5/suppl/Ramilowski_et_al_2015/). NATMI was tested using python 2.7 and 3.7 versions. 

### ExtractEdges: Extracting ligand-receptor-mediated interactions between cell types in the input transcriptome data.
[transcriptome data--> once we agree on how we will word this in the paper, we can modify this]

```
ExtractEdges.py [-h] [--species SPECIES] --emFile EMFILE [--annFile ANNFILE] [--signalType SIGNALTYPE] [--coreNum CORENUM] [--out OUT]

optional arguments: [optional and required arguments are mixed]
  -h, --help            show this help message and exit
  --species SPECIES     only human and mouse expression data are currently supported, default is "human"
  --emFile EMFILE       the path to the expression matrix file with row names (gene symbols) and column names (cell names/single-cell identifiers)
  --annFile ANNFILE     the path to the metafile in which column one has single-cell identifiers and column two has corresponding cluster IDs (see file 'toy.sc.ann.txt' as an example)
  --signalType SIGNALTYPE
                        lrc2p (default) has literature supported ligand-receptor pairs | lrc2a has putative and literature supported ligand-receptor pairs, folder name of the interaction database
  --coreNum CORENUM     the number of CPU cores used, default is one
  --out OUT             the path to save the analysis results [is this optonal or required?]
```

*Note: Expression matrix and metafile are supported in csv, tsv, txt, xls or xlsx format.

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

*Note: Only weight changes across two condition from the same, or similar, datasets and with the same 'signalType' of ligand-receptor pairs (literature-supported with literature-supported or all with all) should be compared. 

```
DiffEdges.py [-h] --refFolder REFFOLDER --targetFolder TARGETFOLDER [--signalType SIGNALTYPE] [--weightType WEIGHTTYPE] [--out OUT]

optional arguments:  [optional and required arguments are mixed]
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

```
VisInteractions.py [-h] --sourceFolder SOURCEFOLDER [--signalType SIGNALTYPE] [--weightType WEIGHTTYPE] [--specificityThreshold SPECIFICITYTHRESHOLD]
                          [--expressionThreshold EXPRESSIONTHRESHOLD] [--detectionThreshold DETECTIONTHRESHOLD]
                          [--keepTopEdge KEEPTOPEDGE] [--plotWidth PLOTWIDTH] [--plotHeight PLOTHEIGHT] [--plotFormat PLOTFORMAT]
                          [--edgeWidth EDGEWIDTH] [--clusterDistance CLUSTERDISTANCE] [--drawClusterPair DRAWCLUSTERPAIR]
                          [--layout LAYOUT] [--fontSize FONTSIZE] [--maxClusterSize MAXCLUSTERSIZE]
                          [--drawNetwork DRAWNETWORK] [--drawLRNetwork [LIGAND [RECEPTOR ...]]]

optional arguments:
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

*Note: Python libraries [seaborn](https://seaborn.pydata.org/), [igraph](https://igraph.org/python/), [NetworkX](https://networkx.github.io/) and [PyGraphviz](https://pygraphviz.github.io/) are required to visualise the cell-to-cell communication network at three distinct levels. NATMI was tested using seaborn 0.8.1, igraph 0.7.1, NetworkX 2.1 and PyGraphviz 1.5 versions.

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


### Example workflow (Extract edges in xxx and yyy, identify and visualize resulting differential edges.)
(make a toy example of something you are showing in the paper so the reviewer can run it all at once and put it all to the "test folder". For VisInteractions.py, consider if your script can check if seaborn, igraph, NetworkX and PyGraphviz are installed and if not, warn the user and proceed to install):

*NotePython libraries [seaborn](https://seaborn.pydata.org/), [igraph](https://igraph.org/python/), [NetworkX](https://networkx.github.io/) and [PyGraphviz](https://pygraphviz.github.io/) are required. NATMI will proceed to instalation if they are missing. 

```bat
   python ExtractEdges.py --species mouse --emFile toy.sc.em.txt --annFile toy.sc.ann.txt --signalType lrc2p --coreNum 4
   python ExtractEdges.py --species mouse --emFile toy.sc.em.txt --annFile toy.sc.ann.txt --signalType lrc2p --coreNum 4
   python DiffEdges.py --refFolder /path/to/ExtractEdges.py's/output/folder/of/reference/dataset --targetFolder     /path/to/ExtractEdges.py's/output/folder/of/target/dataset --signalType lrc2p
   python VisInteractions.py --sourceFolder /path/to/result/folder --signalType lrc2p --weightType mean --detectionThreshold 0.2 --plotFormat pdf --drawNetwork y --plotWidth 12 --plotHeight 10 --layout kk --fontSize 8 --edgeWidth 0 --maxClusterSize 0 --clusterDistance 1
```
*Results are in the test folder
