#! /usr/bin/Rscript
library(data.table)
suppressPackageStartupMessages(library(tidyverse))
library(simpIntLists)
library(grandforest)

setwd('.')
edges = read.table('../../temp/gf_ggi.txt', sep = "\t", header = TRUE)
D <- read.table("../../temp/gf_expr.txt", sep = ",", header= TRUE)
D$phenotype<- as.factor(D$phenotype)
colnames(D) <- gsub("^X", "",  colnames(D))

model <- grandforest(
  data = D,
  graph_data = edges,
  dependent.variable.name = "phenotype",
  num.trees = 10000
)

print(model)


results <- importance(model) %>%
  sort(decreasing=TRUE) %>%
  head(100) %>%
  enframe 

fileConn<-file("../../temp/gf_output.txt")
writeLines(results$name, fileConn)
close(fileConn)

