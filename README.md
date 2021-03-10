#  Testing limits of active module identification
## Table of contents
* [Python dependencies](#python_dependencies)
* [Reproduction of the manuscript results](#reproduction)
* [Run the pipeline](#run)
* [Test your own NEM](#test)
* [Citing](#citing)
* [Contact](#contact)

## Python dependencies <a name="python_dependencies"></a>

You need a Python3 interpreter for running the tests. Moreover, you need to install all Python packages listed in ``testsuite/requirements.txt``, e.g., using Anaconda.

## Reproduction of the manuscript results <a name="reproduction"></a>

### ClustEx2

- J. Gu, Y. Chen, S. Li, Y. Li (2010): "Identification of responsive gene modules by network-based gene clustering and extending: application to inflammation and angiogenesis", *BMC Syst. Biol.*, https://doi.org/10.1186/1752-0509-4-47.
- Z. Ding, W. Guo, J. Gu (2018): "ClustEx2: Gene module identification using density-based network hierarchical clustering", *CAC 2018*, pp. 2407-2412, https://doi.org/10.1109/CAC.2018.8623442.

Compilation under Linux:

```shell script
cd algorithms/clustex2
./install_dependencies.sh
./install_linux.sh
```

Compilation under macOS:

```shell script
cd algorithms/clustex2
./install_dependencies.sh
./install_macos.sh
```

### GXNA

- Ş. Nacu, R. Critchley-Thorne, P. Lee, S. Holmes (2007): "Gene  expression network analysis and applications to immunology", *Bioinformatics*, https://doi.org/10.1093/bioinformatics/btm019.

Compilation under Linux and macOS:

```shell script
cd algorithms/gxna
./install.sh
```

### DIAMOnD

- S.D. Ghiassian, J. Menche J, A.L. Barabási (2015): "A DIseAse MOdule Detection (DIAMOnD) algorithm derived from a  systematic analysis of connectivity patterns of disease proteins in the  human interactome", *PLOS Comp. Biol.*, https://doi.org/10.1371/journal.pcbi.1004120

Distributed as Python script, no compilation necessary.

### GiGA

- R. Breitling, A. Amtmann, P. Herzyk (2004): "Graph-based iterative Group Analysis enhances microarray interpretation", *BMC Bioinform.*, https://doi.org/10.1186/1471-2105-5-100.

Distributed as Perl script, no compilation necessary.


### KeyPathwayMiner

- N. Alcaraz, H. Kücük, J. Weile, A. Wipat, J. Baumbach (2011): "KeyPathwayMiner: Detecting case-specific biological pathways using expression data", *Internet Mathematics*, https://doi.org/10.1080/15427951.2011.604548.
- N. Alcaraz, J. Pauling, R. Batra, et al. (2014): "KeyPathwayMiner 4.0: condition-specific pathway analysis by combining multiple omics studies and networks with Cytoscape". *BMC Syst. Biol.*, https://doi.org/10.1186/s12918-014-0099-x.
- M. List, N. Alcaraz, M. Dissing-Hansen, H.J. Ditzel, J. Mollenhauer, J. Baumbach (2016): "KeyPathwayMinerWeb: online multi-omics network enrichment", *Nucleic Acids Res.*, https://doi.org/10.1093/nar/gkw373.

Distributed as jar executable, no compilation necessary.


### GrandForest 

- S.J. Larsen, H.H.H.W. Schmidt, J. Baumbach (2020): "De novo and supervised endophenotyping using network-guided ensemble learning", *Systems Medicine*, https://doi.org/10.1089/sysm.2019.0008.

Installation instructions can be found in the next paragraph.

### COSINE

- H. Ma, E.E. Schadt, L.M. Kaplan, H. Zhao (2011): "COSINE: COndition-SpecIfic sub-NEtwork identification using a global optimization method", *Bioinformatics*, https://doi.org/10.1093/bioinformatics/btr136.

To install GrandForest and COSINE, open an R console with sudo rights and run the following commands:

```R
install.packages(c("data.table", "tidyverse", "devtools", "COSINE", "tidyr"))
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
BiocManager::install("simpIntLists")
devtools::install_github("SimonLarsen/grandforest")
```

Under macOS, the installation of ``data.table`` will likely fail. Follow the instructions given here: https://github.com/Rdatatabl`e/data.table/wiki/Installation.

### DOMINO
Follow instructions for installations here https://github.com/Shamir-Lab/DOMINO
Please pay attention if `pcst-fast` package was succesfully installed. Otherwise install it with `python setup.py install` from here https://github.com/fraenkel-lab/pcst_fast

## To run the pipeline <a name="run"></a>

Navigate to the directory , activate the environment and run the pipeline in the parallel mode:
```shell
python run_tests.py parallel
```
To see all available options:

```shell
python run_tests.py parallel -h
```
You can also execute only one combination at a time by running the script in a sequential mode. To see all parameters run:
```shell
python run_tests.py sequential -h

```
## Test your own AMIM <a name="test"></a>

To test your own method you need the following:
1. Basic knowledge of python to follow provided instructions.
2. A way to run your method from a command-line.
3. All dependencies needed for your method.

###  Step-by-step instruction
* Download the repository and open the provided template *testsuite/custom_wrapper.py*.
* Do not change any lines from *testsuite/custom_wrapper.py*. Just add what you need to make your method work.
* For your convince the wrapper script is divided in 4 parts that are described in the example bellow:
 
Please keep the "prefix" variable exactly as shown in the example. It allows to distinguish this run from all others, parse results and remove all temp files.


1: Network is passed in the variable *ggi_network* as networkx graph where "GeneID" attribute correspond to entrez gene ID.
The network should be then saved in a format that is appropriate for you tool. For instance, if your tool requiers a .tsv file without a header:
```python
# 1. Write GGI network in format required by your method
path_to_network = f'../temp/{prefix}_custom_ggi.txt'
with open(path_to_network, 'w') as edge_list_file:
    gene_ids = nx.get_node_attributes(ggi_network, 'GeneID')
    for u, v in ggi_network.edges():
        edge_list_file.write(f'{gene_ids[u]}\t{gene_ids[v]}\n')
```
NOTE: here we discuss only formatting and not the actual content of your files. If you wish to include your own network, you need to save it as `CUSTOM.tsv` file without a header in the `data/networks/` directory and then run the pipeline with `--network CUSTOM`

2: Save gene expression and phenotype in the nessesary format. In the example bellow we save phenotype as one of columns in gene expression and save both as .txt. "phenotypes" variable is a list of phenotypes given for this condition 

```python
# 2. Write expression data and phenotype in the format required by your method
path_to_expression = f'../temp/{prefix}_custom_expr.txt'
expression_data["phenotype"] = phenotypes
expression_data.to_csv(path_to_expression, index = False)
```

3: 
Write the command for your tool. 
Don't forget to navigate to the corresponding folder and pass prefix to your method such that it can save the output to *path_to_output*.

```python
# 3. Insert the command to run your method, direct the output to path_to_output
path_to_output = f'../temp/{prefix}_custom_output.txt'
command = f'cd ../algorithms/your_method_directory/; ./your_method.sh {prefix}"
subprocess.call(command, shell = True, stdout=subprocess.PIPE)
```
4: Add the resulting gene as list of strings (Entrez ids)

```python
# 4. Process results such that they are formatted as a list of strings (entez IDs)
result_genes = []
with open(path_to_output, 'r') as results:
    for line in results:
        result_genes.append(line.strip())
```
* When the wrapper is done, all you need to do is to run the pipline:
```shell
#for sequential execution
python run_tests.py -sequential --method CUSTOM  
# OR for parallel execution
python run_tests.py -parallel --methods CUSTOM  

```
The results will be stored in the results folder in the following format: NETWORK_GENERATOR_CUSTOM.csv

### Visualize your results
Please note, that minimal information that is needed for the visualization includes 1 run on the original network (any network) and 1 run on any generator for any condition.
For instance:
```shell
#for sequential execution
python run_tests.py sequential --network HPRD --generator ORIGINAL --method CUSTOM  --condition GSE3790 --verbose
python run_tests.py sequential --network HPRD --generator RDPN --method CUSTOM  --condition GSE3790 --verbose

#OR for parallel
python run_tests.py parallel --networks HPRD --generators ORIGINAL RDPN --methods CUSTOM  --conditions GSE3790 --verbose
```
To visualize the results please run:

```shell
python show_plots.py
```

## Citing <a name="citing"></a>

If you use this test-suite for a scientific publication, please cite the following paper:

- O. Lazareva, J. Baumbach, M. List, D. B. Blumenthal: On the limits of active module identification, **Brief. Bioinform.**, 2021, [DOI: 10.1093/bib/bbab066](https://doi.org/10.1093/bib/bbab066).  

## Contact <a name="contact"></a>
If you experience any difficulties or need more options for your own tool evaluation, please reach out to us:

* [Olga Lazareva](mailto:olga.lazareva@wzw.tum.de)
* [David Blumenthal](mailto:david.blumenthal@wzw.tum.de)











