# Compared Network Enrichment Algorithms

## ClustEx2

- Z. Ding, W. Guo and J. Gu, "ClustEx2: Gene Module Identification using Density-Based Network Hierarchical Clustering," *2018 Chinese Automation Congress (CAC)*, Xi'an, China, 2018, pp. 2407-2412, https://doi.org/10.1109/CAC.2018.8623442.

Compilation under Linux:

```sh
cd clustex2
./install_dependencies.sh
./install_linux.sh
```

Compilation under macOS:

```sh
cd clustex2
./install_dependencies.sh
./install_macos.sh
```

## GXNA

- Şerban Nacu, Rebecca Critchley-Thorne, Peter Lee, Susan Holmes,  Gene  expression network analysis and applications to immunology, *Bioinformatics*, Volume 23, Issue 7, 1 April 2007, Pages 850–858, https://doi.org/10.1093/bioinformatics/btm019

Compilation under Linux and macOS:

```sh
cd gxna
./install.sh
```

## DIAMOnD

- Ghiassian SD, Menche J, Barabási AL   (2015)     A DIseAse MOdule Detection (DIAMOnD) Algorithm Derived from a  Systematic Analysis of Connectivity Patterns of Disease Proteins in the  Human Interactome. PLOS Computational Biology  11(4): e1004120. https://doi.org/10.1371/journal.pcbi.1004120

Distributed as Python script, no compilation necessary.

## Grand Forest and COSINE

install.packages(c("data.table", "tidyverse", "devtools", "COSINE", "tidyr"))
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("simpIntLists")
devtools::install_github("SimonLarsen/grandforest")