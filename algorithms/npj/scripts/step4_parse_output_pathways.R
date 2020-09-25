library(gdata)
load("/path_to_/scripts/parse_pathways.rda")

# for each tool 
# res.dir is where pathways are stored
# fg.dir is the directory where fg is stored. with one to one correspondence to pathway file

######################## PINNACLE #############path_to_##################
res.dir <- "/path_to_pathway_dir/pinnacle/"
fg.dir <- "/path_to_fg_dir/"
res.mat <- getPinnacleResultMatrix(path.dir=res.dir, fg.dir)
save(res.mat, file=paste(res.dir, "res_mat.RData",sep=""))

######################## KPM ###############################
res.dir <- "/path_to_pathway_dir/kpm/"
fg.dir <- "/path_to_fg_dir/"
res.mat <- getKPMResultMatrix(g.path.dir=res.dir, 
                              algorithm="GREEDY", strategy="GLONE", 
                              fg.dir)
save(res.mat, file=paste(res.dir, "res_mat.RData",sep=""))

######################## GXNA ###############################
res.dir <- "/path_to_pathway_dir/gxna/"
fg.dir <- "/path_to_fg_dir/"
res.mat <- getGxnaResultMatrix(path.dir=res.dir, fg.dir)
save(res.mat, file=paste(res.dir, "res_mat.RData",sep=""))

######################## GiGA ###############################
res.dir <- "/path_to_pathway_dir/giga/"
fg.dir <- "/path_to_fg_dir/"
res.mat <- getGigaResultMatrix(path.dir=res.dir, fg.dir)
save(res.mat, file=paste(res.dir, "res_mat.RData",sep=""))

######################## DEGAS ###############################
res.dir <- "/path_to_pathway_dir/degas/"
fg.dir <- "/path_to_fg_dir/"
res.mat <- getDegasResultMatrix(path.dir=res.dir, 
                                fg.dir)
save(res.mat, file=paste(res.dir, "res_mat.RData",sep=""))

######################## COSINE ###############################
res.dir <- "/path_to_pathway_dir/cosine/"
fg.dir <- "/path_to_fg_dir/"
path.rdatas <- list.files(path=res.dir, pattern="pathway", full.names=T)
res.mats <- lapply(path.rdatas, FUN=function(path.rdata) getCosineResultMatrix(path.rdata, fg.dir))
res.mat <- data.frame(do.call(rbind, res.mats))
save(res.mat, file=paste(res.dir, "res_mat.RData",sep=""))

######################## BioNet ###############################
res.dir <- "/path_to_pathway_dir/bionet/"
fg.dir <- "/path_to_fg_dir/"
path.rdatas <- list.files(path=res.dir, pattern="pathway", full.names=T)
res.mats <- lapply(path.rdatas, FUN=function(path.rdata) getBionetResultMatrix(path.rdata, fg.dir))
res.mat <- data.frame(do.call(rbind, res.mats))
save(res.mat, file=paste(res.dir, "res_mat.RData",sep=""))
