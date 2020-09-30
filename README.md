# Testing the one-network-fits-all hypothesis 

## Algorithms

### ClustEx2

- Z. Ding, W. Guo, J. Gu (2018): "ClustEx2: Gene module identification using density-based network hierarchical clustering", *CAC 2018*, pp. 2407-2412, https://doi.org/10.1109/CAC.2018.8623442.
- Gu J, Chen Y, Li S, Li Y (2010): "Identification of responsive gene modules by network-based gene clustering and extending: application to inflammation and angiogenesis", *BMC Syst. Biol.*, https://doi.org/10.1186/1752-0509-4-47.

Compilation under Linux:

```sh
cd algorithms/clustex2
./install_dependencies.sh
./install_linux.sh
```

Compilation under macOS:

```sh
cd algorithms/clustex2
./install_dependencies.sh
./install_macos.sh
```

### GXNA

- Ş. Nacu, R. Critchley-Thorne, P. Lee, S. Holmes (2007): "Gene  expression network analysis and applications to immunology", *Bioinformatics*, https://doi.org/10.1093/bioinformatics/btm019

Compilation under Linux and macOS:

```sh
cd algorithms/gxna
./install.sh
```

### DIAMOnD

- S.D. Ghiassian, J. Menche J, A.L. Barabási (2015): "A DIseAse MOdule Detection (DIAMOnD) algorithm derived from a  systematic analysis of connectivity patterns of disease proteins in the  human interactome", *PLOS Comp. Biol.*, https://doi.org/10.1371/journal.pcbi.1004120

Distributed as Python script, no compilation necessary.

### GrandForest and COSINE

```shell script
sudo R
install.packages(c("data.table", "tidyverse", "devtools", "COSINE", "tidyr"))
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("simpIntLists")
devtools::install_github("SimonLarsen/grandforest")
```

Under macOS, the installation of `data.table` will likely fail. Follow the instructions given here: https://github.com/Rdatatable/data.table/wiki/Installation.

