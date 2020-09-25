# STEP 3: compute pathways
load("/path_to_/scripts/compute_pathways.rda")
library(BioNet)
library(igraph)
library(COSINE)
library(plyr)
library(parallel)
# for each tool mcores is for cores to use in parallel computation
######################## PINNACLE ###############################

pinnacle.dir <- "/path_to_/tools/pinnacle/" # directory with pinnacle tools
net.file  <- "/path_to_/sample_network/graph_hprd.sif" # network in sif
eps.dir <- "/path_to_/input_pinnacle/" # directory with files to use for pathway computation
path.dir <- "/path_to_pathway_dir/pinnacle/" # output dir
# get pathways, upM/lowM/stepM limit for M parameter of pinnacle that determines the size of output
#getPinnaclePathways (net.file, eps.dir, path.dir, pinnacle.dir, mcores=10, upM=10, lowM=10, stepM=0)

######################## KPM ###############################
net.file  <- "/path_to_/sample_network/graph_hprd.sif" # network in sif
kpm.loc <- "/path_to_/tools/kpm/" # directory with KPM jar file
input.dir <- "/path_to_/input_kpm/" # directory with files to use for pathway computation
res.dir <- "/path_to_pathway_dir/kpm/"  # output dir
# get pathways, upL/lowL/stepL  upK/lowK/stepK range for l7k parameter of kpm that determine the exceptions
# percentage : 0 --> the l values provided are real values; 1 --> l values are in percentages
# rest are KeypathwayMiner parameters, please refer to KPM help page for more information
#getKPMPathways (input.dir, kpm.loc, res.dir, algorithms=c("GREEDY"), strategies=c("GLONE"), combine.formula = "(L1 || L2)", 
#                combine.operation = "OR", percentage=0, lowL=250, upL=250, stepL=0,
#                lowK=0, upK=0, stepK=0, net.file)
######################## GXNA ###############################
gxna.dir <- "/path_to_/tools/gxna/" # tool
net.file  <- "/path_to_/tools/gxna/hprd.gra" #network
input.dir <- "/path_to_/input_gxna/" # input files
res.dir <- "/path_to_pathway_dir/gxna/"  # output dir
# get pathways, upM/lowM/stepM limit for M parameter of gxna that determines the size of output
#getGxnaPathways (net.file, input.dir, res.dir, gxna.dir, mcores=10, upM=10, lowM=10, stepM=0)

######################## GiGA ###############################
giga.dir <- "/path_to_/tools/giga/" # tool
setwd(giga.dir)
net.file  <- "/path_to_/sample_network/graph_hprd.txt" # network
ths <- c(10, 20, 30, 40, 50) # threshold parameter to determine the size of output
input.dir <- "/path_to_/input_giga/" # input
res.dir <- "/path_to_pathway_dir/giga/"  # output dir
#getGigaPathways(input.dir, res.dir, net.file, ths, num.cores=10)

######################## DEGAS ###############################
degas.dir <- "/path_to_/tools/degas/Expander/" # tool
net.file  <- "/path_to_/sample_network/graph_hprd.sif" # # network file with two columns indicating an edge in the network
input.dir <- "/path_to_/input_degas/" # input
res.dir <- "/path_to_pathway_dir/degas/"  # output dir
#upK/lowK/stepK range for K parameter of degas that determines the exceptions
#getDegasPathways (net.file, input.dir, res.dir, degas.dir, mcores=10, upK=4, lowK=4, stepK=0)

######################## COSINE ###############################
net.file  <- "/path_to_/sample_network/graph_hprd.txt" # # network file with two columns indicating an edge in the network
input.dir <- "/path_to_/input_cosine/" # input directory
res.dir <- "/path_to_pathway_dir/cosine/"  # output dir
exp.profiles <- list.files(input.dir, pattern=".txt", full.names=T) # list of mol profiles
snet.sizes <- c(20) # size of subgrpahs
res.files <- getCosinePathways (exp.profiles, net.file, res.dir, snet.sizes, mcores=10)

######################## BioNet ###############################
net.file  <- "/path_to_/sample_network/graph_hprd.txt" # # network file with two columns indicating an edge in the network
input.dir <- "/path_to_/input_bionet/" # input files
res.dir <- "/path_to_pathway_dir/bionet/"  # output dir
fdrs <- c(0.5) # false discovery rate of pathways
mol.profiles <- list.files(input.dir, pattern=".txt", full.names=T)
res.files <- getBionetPathways(mol.profiles, net.file, res.dir, fdrs, mcores=10)
