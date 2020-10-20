# NetCore

NetCore is a network propagation approach based on node coreness, for phenotype-
genotype associations and module identification. NetCore addresses the node degree bias in PPI
networks by using node coreness in the random walk with restart procedure, and achieves better
re-ranking of genes after propagation. Furthermore, NetCore implements a semi-supervised
approach to identify network modules, which include both well-known genes together with new
candidates.

NetCore's workflow cosists of three parts:
1. Data initialization - includes the extraction of a high quality PPI network, data collection and extraction of aseed gene list from a manually curated database. 
2. Network propagation using node coreness.
3. Module identification in a semi-supervised fashion combining both network propagation results and the seed gene list.

## Step Up

### Download
	git clone https://github.molgen.mpg.de/barel/NetCore

### Install

The following software and packages are required for running NetCore:

* Linux/Unix
* [Python (3.6 or 3.7)](http://python.org/)
* [NumPy (1.16.3)](http://www.numpy.org/)
* [SciPy (1.2.1)](http://www.scipy.org/)
* [Pandas (0.24.2)](https://pandas.pydata.org/)
* [matplotlib (3.1.1)](https://matplotlib.org/3.1.3/)
* [seaborn (0.9.0)](https://seaborn.pydata.org/)
* [NetworkX (2.3)](http://networkx.github.io/)

After the directory was cloned, please run the following to install NetCore:

	python setup.py install

It is also possible to run:

	pip3 install NetCore


## Tutorial
Explanation on how to run NetCore is avilable in the [tutorial notebook](https://github.molgen.mpg.de/barel/NetCore/blob/master/Tutorial.ipynb)

## Data
Example data for running NetCore is in the data subdirectory.

### Protein Protein Interaction network
The CPDB PPI high confidence network [1] is provided as an edge list.
The CPDB PPI network can be downloaded via [ConsensusPathDB](http://cpdb.molgen.mpg.de/)

### Type II diabetes
GWAS data for Type II diabetes is provided from:
* [The GWAS Catalog](https://www.ebi.ac.uk/gwas/efotraits/EFO_0001360) - all the gene associations were downloaded and the p-values were converted to weights using -log10.
* [GWAS list](https://github.com/idekerlab/Network_Evaluation_Tools/blob/master/Data/GWAS_Catalog_genesets.txt) - the genes in this list as known to be associated with Type II Diabetes and can be used as seed genes for NetCore's module identification.

1 [Barel, G., & Herwig, R. (2018). Network and Pathway Analysis of Toxicogenomics Data. Frontiers in genetics, 9, 484.](https://doi.org/10.3389/fgene.2018.00484)