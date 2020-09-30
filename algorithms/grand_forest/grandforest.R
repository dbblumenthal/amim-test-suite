#! /usr/bin/Rscript
library(data.table)
suppressPackageStartupMessages(library(tidyverse))
library(simpIntLists)
library(grandforest)

setwd('.')
args = commandArgs(trailingOnly=TRUE)
prefix = args[1]
edges = read.table(paste("../../temp/", prefix, "_gf_ggi.txt", sep=""), sep = "\t", header = TRUE)
D <- read.table(paste("../../temp/", prefix, "_gf_expr.txt", sep=""), sep = ",", header= TRUE)
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

fileConn<-file(paste("../../temp/", prefix, "_gf_output.txt", sep=""))
writeLines(results$name, fileConn)
close(fileConn)

