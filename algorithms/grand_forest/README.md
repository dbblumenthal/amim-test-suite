# Installation Instructions

Run the following commands in a terminal.

```shell script
sudo R
install.packages(c("data.table", "tidyverse", "devtools"))
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("simpIntLists")
devtools::install_github("SimonLarsen/grandforest")
```

Under macOS, the installation of `data.table` will likely fail. Follow the instructions given here: https://github.com/Rdatatable/data.table/wiki/Installation.
 