library(geophyttertools)
library(ape)
library(ggplot2)


# NYh3n2_HA_20000_4_1993.txt: the original data.
# NYh3n2_HA_20000_4_1993.col: the output of geophytterplus.FitLFMTriangle. Separated into two parts, .tree and .colc. 
    # NYh3n2_HA_20000_4_1993.tree: part of .col file, BHV PCA trees.
    # NYh3n2_HA_20000_4_1993.colc: part of .col file, the projected trees.
# NYh3n2_HA_20000_4_1993.trop: the output of geophytterplus.DecomposeLFMTriangle, used as BHV triangle background.

# Check http://www.mas.ncl.ac.uk/~ntmwn/geophytterplus/index.html for more details.
# projected trees should be Newick format
proj.lines <- readLines("E:/working/5/NYh3n2_HA_20000_5_2011.colc")

proj.lines <- proj.lines[2:20000]
outliers <- c(NA)

# Trees: Should be Newick format
BHVTrees <- read.tree("E:/working/5/NYh3n2_HA_20000_5_2011.tree")


BHVTree_plot <- function(proj_trees, BHVTrees, outliers){
  if(sum(is.na(outliers)) == 0){
    proj.lines <- proj_trees[-c(outliers)]
  }
  
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
  for(i in 1:length(summary)) proj_bhv_trees[i+3] <- trees_bhv[as.numeric(names(summary)[i])]
  # Replace the leaves label with new labels
  n <- length(trees_bhv[[1]]$tip.label) 
  
  
  for(i in 1:length(proj_bhv_trees)) proj_bhv_trees[[i]]$edge.length <- NULL  # ceiling(length(trees)/3)*3
  
  layout_row_n <- ceiling(length(proj_bhv_trees)/3)*3
  
  par.old <- par(mar = c(0.2,0.25,1.1,0.25), xpd = NA, font.main = 1) # Define the parameters of the figure
  on.exit({
    par(par.old)
    layout(1)
  }) 
  layout(matrix(1:layout_row_n, nc = 3, byrow = TRUE))
  
  
  if(n == 4){
    tip.orig <- c("Sample1", "Sample2", "Sample3","Sample4")
    tip.abbrv <- c("1", "2", "3", "4")
    for (i in seq_along(proj_bhv_trees)){
      tord <- match(proj_bhv_trees[[i]]$tip.label, tip.orig)
      proj_bhv_trees[[i]]$tip.label <- tip.abbrv[tord]
    } 
  }else if(n == 5){
    tip.orig <- c("Sample1", "Sample2", "Sample3","Sample4","Sample5")
    tip.abbrv <- c("1", "2", "3", "4", "5")
    for (i in seq_along(proj_bhv_trees)){
      tord <- match(proj_bhv_trees[[i]]$tip.label, tip.orig)
      proj_bhv_trees[[i]]$tip.label <- tip.abbrv[tord]
    } 
  }
  
  topo_title <- c()
  for(i in seq_along(summary)) topo_title[i] <- paste("Topology ",i," (",summary[i],")", sep="")
  
  
  treelabs <- c("Tree 1", "Tree 2", "Tree 3", topo_title)
  
  # Plot the trees.
  for(i in seq_along(proj_bhv_trees)) ape::plot.phylo(proj_bhv_trees[[i]], main=treelabs[i], type="c", direction="downwards", srt = 90, adj = 0.5,
                                             label.offset = 0.2, cex = 2, cex.main = 2,  edge.width = 2, font = 2, font.main = 2)
}



# BHV Triangle
# The ".trop" file is define by user. Basically, it is the output of the JAVA program, geophytterplus.DecomposeLFMTriangle.
# The ".colc" file is define by user. It is a PART of the output of the JAVA program, geophytterplus.FitLFMTriangle.
# "read.topologies" function is from r package geophyttertools.
# "read.projections" function is from r package geophyttertools.

background <- read.topologies("E:/working/5/NYh3n2_HA_20000_5_2011.trop")
points <- read.projections("E:/working/5/NYh3n2_HA_20000_5_2011.colc")

BHV_Triangle <- function(background, points){
  pp <- ggplot(background) + corner_labels(offset=0.03, size = 6) +
    guides(fill=guide_legend(title="Topology")) +
    #scale_fill_brewer(palette="Set3") +
    scale_fill_manual(values=c("#8dd3c7","#ffffb3","#bebada","#fb8072","#b3de69","#fdb462","#80b1d3","#fccde5","#d9d9d9")) +
    geom_point(aes(x=x,y=y), points, inherit.aes = FALSE, size=0.5) +
    theme(legend.position = c(.82,.76),legend.text=element_text(size=16, face = "bold"), legend.title = element_text(size = 16,face = "bold"))
  pp
}
