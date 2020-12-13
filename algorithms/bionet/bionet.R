library(BioNet)
library(DLBCL)
library(igraph)

data(dataLym)
data(interactome)

pvals <- cbind(t = dataLym$t.pval, s = dataLym$s.pval)
rownames(pvals) <- dataLym$label
pval <- aggrPvals(pvals, order = 2, plot = FALSE)

subnet <- subNetwork(dataLym$label, interactome)
subnet <- rmSelfLoops(subnet)
subnet
fb <- fitBumModel(pval, plot = FALSE)
scores <- scoreNodes(subnet, fb, fdr = 0.001)
module <- runFastHeinz(subnet, scores)
module@nodes

edges = as.matrix(read.table("../../temp/HPRD_ORIGINAL_bionet_ggi.graphml", sep = "\t", header = TRUE))

net = graph_from_edgelist(edges, directed = FALSE)
pv =read.csv("../../temp/HPRD_ORIGINAL_bionet_pvals.txt", sep = ",", header = TRUE, row.names="X")

fb <- fitBumModel(pv, plot = FALSE)
scores <- scoreNodes(net, fb, fdr = 0.001)


