library(geophyttertools)
library(ape)
library(phangorn)
library(ggplot2)

BHVTree_plot <- function(proj_trees, BHVTrees, outName, frechet_mean){

  tree_list_ori <- strsplit(proj.lines, split= " ")
  tree_list <- lapply(tree_list_ori, function(x) x[6])
  
  trees_bhv <- list()
  for(i in 1:length(tree_list)){
    trees_bhv[[i]] <- as.list(ape::read.tree(text=tree_list[[i]]))
  }
  
  type <- seq(length(trees_bhv))
  for(i in seq(length(trees_bhv))){
    if(i==type[i]) {
      for(j in seq(i, length(trees_bhv))) {
        if(RF.dist(trees_bhv[[i]],trees_bhv[[j]])==0) {
          type[j]=i
        }
      } 
    }
  }
  summary <- table(type)
  proj_bhv_trees <- BHVTrees
  proj_bhv_trees[4] <- frechet_mean 
  for(i in 1:length(summary)) proj_bhv_trees[i+4] <- trees_bhv[as.numeric(names(summary)[i])]
  # Replace the leaves label with new labels
  n <- length(trees_bhv[[1]]$tip.label) 
  
  
  for(i in 1:length(proj_bhv_trees)) proj_bhv_trees[[i]]$edge.length <- NULL  # ceiling(length(trees)/3)*3
  
  layout_row_n <- ceiling(length(proj_bhv_trees)/3)*3
  
  jpeg(paste("./BHVTrees",outName,".jpg",sep=""),width=800, height=300)
  par.old <- par(mar = c(0.2,0.25,1.1,0.25), xpd = NA, font.main = 1) # Define the parameters of the figure
  on.exit({
    par(par.old)
    layout(1)
  }) 
  layout(matrix(1:layout_row_n, nc = 3, byrow = TRUE))
  
  
  tip.orig <- c("Sample1", "Sample2", "Sample3","Sample4","Sample5")
  tip.abbrv <- c("1", "2", "3", "4", "5")
  for (i in seq_along(proj_bhv_trees)){
    tord <- match(proj_bhv_trees[[i]]$tip.label, tip.orig)
    proj_bhv_trees[[i]]$tip.label <- tip.abbrv[tord]
  } 
  
  topo_title <- c()
  for(i in seq_along(summary)) topo_title[i] <- paste("Topology ",i," (",summary[i],")", sep="")
  
  
  treelabs <- c("Tree 1", "Tree 2", "Tree 3","Fréchet mean", topo_title)
  
  # col_value assigns colors to the points
  
  col_value <- c("#000000","#000000","#000000","#000000","#8dd3c7","#696969","#bebada","#fb8072","#b3de69","#fdb462","#80b1d3","#fccde5","#d9d9d9")
  for(i in seq_along(proj_bhv_trees)) ape::plot.phylo(proj_bhv_trees[[i]], main=treelabs[i], type="c", direction="downwards", srt = 90, adj = 0.5,
                                             label.offset = 0.2, cex = 2, cex.main = 2,  edge.width = 2, font = 2, font.main = 2,col.main = col_value[i])
  dev.off()
}

# BHV Triangle
# The ".trop" file is an output of the java program, geophytterplus.DecomposeLFMTriangle.
# The ".colc" file is an output of the java program, geophytterplus.FitLFMTriangle.
# "read.topologies" function is from the R package geophyttertools
# "read.projections" function is from from R package geophyttertools

BHV_Triangle <- function(background, points, outName){
  
  jpeg(paste("./BHVTriangle",outName,".jpg",sep=""),width=400, height=300)
  pp <- ggplot(background) + corner_labels(offset=0.03, size = 6) +
    guides(fill=guide_legend(title="Topology")) +
    #scale_fill_brewer(palette="Set3") +
    scale_fill_manual(values=c("#8dd3c7","#ffffb3","#bebada","#fb8072","#b3de69","#fdb462","#80b1d3","#fccde5","#d9d9d9")) +
    geom_point(aes(x=x,y=y), points, inherit.aes = FALSE, size=0.5) +
    theme(legend.position = c(.82,.76),legend.text=element_text(size=16, face = "bold"), legend.title = element_text(size = 16,face = "bold"))
  print(pp)
  dev.off()
}

year <- 1993:2013
frechet_mean <- read.tree("./Frechet_mean.txt")

for(i in 1:21){
  # BHV Trees
  proj.lines <- readLines(paste("./ori_colc/N_NYh3n2_HA_20000_5_",year[i],".colc",sep=""))
  proj.lines <- proj.lines[2:(length(proj.lines)-9)]
  BHVTree_plot(proj.lines,BHVTrees,year[i],frechet_mean[i])
  
  background <- read.topologies(paste("./BHV/N_NYh3n2_HA_20000_5_",year[i],".trop",sep=""))
  points <- read.projections(paste("./ori_colc/N_NYh3n2_HA_20000_5_",year[i],".colc",sep=""))
  BHV_Triangle(background, points, year[i])
}
