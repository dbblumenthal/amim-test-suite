# STEP 1: generate FG nodes

load("/path_to/scripts/generate_foreground.rda")
# required packages
library(igraph)
library(parallel)
library(gdata)

# OPTION 1a : get fg genes using AVDk
# net.file: file with the network. sample file graph_hprd.txt
# fg.dir: directory where results should be stored
# algo.loc: location of the gbu executable
# n.u: upper limit of number of fgs
# n.l: lower limit of number of fgs
# n.by: icrement of n.l upto n.u
# k.u: upper limit of average shortest distance between the fgs
# k.l: lower limit of average shortest distance between the fgs
# kd: plus minus set k
# m: number of fg sets
# k.by: increment of k.l upto k.u
algo.loc <-"/path_to/tools/gbu/" 
net.file <- "/path_to/sample_network/graph_hprd.txt"
fg.dir <- "/path_to_fg_dir/"
CallGBU (net.file, fg.dir, algo="gbu", algo.loc, n.u=20, n.l=20, k.u=3, k.l=2, kd=2, m=3, n.by=0, k.by=1)

# OPTION 1b : get fg genes using subgraphs
# net.file: file with the network. sample file graph_hprd.txt
# fg.dir: directory where results should be stored
# snet.fize: number of fgs
# num.snet: number of fg sets
# num.core: if you want to generate it in parallel
net.file <- "/path_to/sample_network/graph_hprd.txt"
fg.dir <- "/path_to_fg_dir/"
num.snet <- 5
num.cores <- 10
snet.size <- 20
CallSgraphs (net.file, fg.dir, snet.size, num.snet, num.cores)
