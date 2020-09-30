# Testing the one-network-fits-all hypothesis 

## Python dependencies and algorithms

### Python dependencies

You need to install all Python packages listed in ``testsuite/requirements.txt``.

### ClustEx2

- Z. Ding, W. Guo, J. Gu (2018): "ClustEx2: Gene module identification using density-based network hierarchical clustering", *CAC 2018*, pp. 2407-2412, https://doi.org/10.1109/CAC.2018.8623442.
- Gu J, Chen Y, Li S, Li Y (2010): "Identification of responsive gene modules by network-based gene clustering and extending: application to inflammation and angiogenesis", *BMC Syst. Biol.*, https://doi.org/10.1186/1752-0509-4-47.

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

### PinnacleZ

-  H.‐Y. Chuang, E. Lee, Y.-T. Liu, D. Lee, T. Ideker (2007): "Network‐based classification of breast cancer metastasis", *Mol. Syst. Biol.*, https://doi.org/:10.1038/msb4100180.

Distributed as jar executable, no compilation necessary.

### KeyPathwayMiner

- N. Alcaraz, H. Kücük, J. Weile, A. Wipat, J. Baumbach (2011): "KeyPathwayMiner: Detecting case-specific biological pathways using expression data", *Internet Mathematics*, https://doi.org/10.1080/15427951.2011.604548.
- M. List, N. Alcaraz, M. Dissing-Hansen, H.J. Ditzel, J. Mollenhauer, J. Baumbach (2016): "KeyPathwayMinerWeb: online multi-omics network enrichment", *Nucleic Acids Res.*, https://doi.org/10.1093/nar/gkw373.

Distributed as jar executable, no compilation necessary.


### GrandForest and COSINE

Open an R console with sudo rights and run the following commands:

```R
install.packages(c("data.table", "tidyverse", "devtools", "COSINE", "tidyr"))
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
BiocManager::install("simpIntLists")
devtools::install_github("SimonLarsen/grandforest")
```

Under macOS, the installation of ``data.table`` will likely fail. Follow the instructions given here: https://github.com/Rdatatabl`e/data.table/wiki/Installation.

