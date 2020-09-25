
# STEP 2: simulate expression data
library(truncnorm)
library(psych)
library(parallel)
load("/path_to_/scripts/simulate_expression_data.rda")

# OPTION 2a: varying variance model
# net.file: file with the network. sample file graph_hprd.txt
# fg.dir: directory where results should be stored
# v.fgs: variance of fg nodes in case samples
# v.fgs.ctrl: variance of fg nodes in control samples
# v.bgs: variance of bg nodes in all the samples
# fg.files: files with fg nodes as generates in STEP1
# m.fgs: mean of fg nodes in case samples
# m.fgs.ctrl: mean of fg nodes in control samples
# m.bgs: mean of bg nodes in all the samples
# sim.num: 1 -> varying variance; 2 -> varying mean
# num.cores: cores to use for the task
fg.dir <- "/path_to_fg_dir/"
net.file <- "/path_to_/sample_network/graph_hprd.txt"
res.dir <- "path_to_fg_dir"
fg.files <- list.files(path=fg.dir, pattern=".txt", full.names=T)
m.bgs <- c(0, 0)
m.fgs <- c(0, 0)
m.fgs.ctrl <- c(0, 0)
v.bgs <- c(5, 0.2)
v.fgs <- c(5, 5)
v.fgs.ctrl <- c(5, 0.2)

CallModel (net.file, res.dir, v.fgs, v.fgs.ctrl, v.bgs, fg.files, model="model.exp.ctrl", m.fgs, m.fgs.ctrl, m.bgs, sim.num=1, num.cores=10)

# OPTION 2b: varying mean model
fg.dir <- "/path_to_fg_dir/"
net.file <- "/path_to_/sample_network/graph_hprd.txt"
res.dir <- "/path_to_fg_dir/"
fg.files <- list.files(path=fg.dir, pattern=".txt", full.names=T)
m.bgs <- c(0, 0)
m.fgs <- c(0, 2)
m.fgs.ctrl <- c(0, 1)
v.bgs <- c(1, 1)
v.fgs <- c(1, 1)
v.fgs.ctrl <- c(1, 1)

CallModel (net.file, res.dir, v.fgs, v.fgs.ctrl, v.bgs, fg.files, model="model.exp.ctrl", m.fgs, m.fgs.ctrl, m.bgs, sim.num=2, num.cores=10)
