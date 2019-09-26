library(ape)
library(phangorn)
library(parallel)

# Call the "func_ssh.R" file which includes auxiliary functions needed to run tropical PCA 
source("./func_ssh.R")

# Define the number of clusters used in parallel computing 
cl <- makeCluster(6)

tropTree_plot <- function(trees_ori, comb_set, outName, pcs=3){
  for(i in 1:length(trees_ori)){
    if(sum(trees_ori[[i]]$edge.length) > 0){
      trees_ori[[i]]$edge.length[trees_ori[[i]]$edge.length < 0] = 0
    }
  }
  
  n <- length(trees_ori[[1]]$tip.label) 
  to <- trees_ori[[1]]$tip.label 
  distVec_all <- distMat(trees_ori, tipOrder = to) 
  N <- length(distVec_all)  
  D_all <- matrix(unlist(distVec_all), ncol = N) 
  
  trees <- trees_ori[comb_set]
  
  base <- distMat(trees, tipOrder = to)
  D_base <- matrix(unlist(base), ncol = pcs)
  adj_proj_points <- lapply(distVec_all, project_pi,D_s = D_base)
  proj_trees_matrix <- lapply(adj_proj_points, make.matrix, n = length(to), tips = to)
  proj_trees <- lapply(proj_trees_matrix, upgma)
  
  type <- seq(length(proj_trees))
  for(i in seq(length(proj_trees))){
    if(i==type[i]) {
      for(j in seq(i, length(proj_trees))) {
        if(RF.dist(proj_trees[[i]],proj_trees[[j]]) == 0) {
          type[j] = i
        }
      } 
    }
  } 
  
  # We now classify the tree topologies.  First, we check how many tree topologies are identified.
  # trees[[1]], trees[[2]], trees[[3]] are the vertices for the tropical triangle, giving the second principal component.
  
  for(i in 1:length(table(type))){
    trees[[3+i]] <- proj_trees[[i]]
  }
  
  # Keep the tree topology.
  for(i in 1:length(trees)) trees[[i]]$edge.length <- NULL  # ceiling(length(trees)/3)*3
  
  layout_row_n <- ceiling(length(trees)/3)*3
  
  jpeg(paste("./tropTree",outName,".jpg",sep=""),width=700, height=400)
  par.old <- par(mar = c(0.2,0.25,1.1,0.25), xpd = NA, font.main = 1) # Define the parameters of the figure
  on.exit({
    par(par.old)
    layout(1)
  }) 
  layout(matrix(1:layout_row_n, nc = 3, byrow = TRUE))
  
  # Replace leaf labels with new labels
  tip.orig <- c("Sample1", "Sample2", "Sample3","Sample4","Sample5")
  tip.abbrv <- c("1", "2", "3", "4", "5")
  for (i in seq_along(trees)){
    tord <- match(trees[[i]]$tip.label, tip.orig)
    trees[[i]]$tip.label <- tip.abbrv[tord]
  } 
  
  topo_title <- c()
  for(i in seq_along(table(type))) topo_title[i] <- paste("Topology ",i," (",table(type)[i],")", sep="")
  
  treelabs <- c("PC1", "PC2", "PC3", topo_title)
  
  # Plot the trees.
  col_value <- c("#000000","#000000","#000000","#8dd3c7","#696969","#bebada","#fb8072","#b3de69","#fdb462","#80b1d3","#fccde5","#d9d9d9")
  
  for(i in seq_along(trees)) ape::plot.phylo(trees[[i]], main=treelabs[i], type="c", direction="downwards", srt = 90, adj = 0.5,
                                             label.offset = 0.2, cex = 2, cex.main = 2,  edge.width = 2, font = 2, font.main = 2,
                                             col.main = col_value[i]) # Plot the trees
  dev.off()
  
  new_base <- D_all[,comb_set]
  DD_base <- t(new_base)
  D_base <- normalize.ultrametrices(DD_base)
  proj_points_plot <- lapply(adj_proj_points, polytope_iso, D = DD_base)
  
  proj_plot_norm <- lapply(proj_points_plot, normalize.proj)
  proj_2D_plot_m <- matrix(unlist(proj_plot_norm), nrow = N, ncol = 3, byrow = T)
  
  type_ind <- as.matrix(table(type))
  
  tree_num <- as.numeric(rownames(type_ind))
  freq <- type
  
  # Specify colors
  for(i in 1:length(tree_num)){
    freq[type == tree_num[i]] <- i
  }
  
  jpeg(paste("./tropTriangle",outName,".jpg",sep=""),width=500, height=400)
  par(mfrow = c(1,1))
  par(mar = c(2,2,0.2,0.2), xpd = NA, font.main = 1)
  k <- ncol(D_base)
  
  plot(D_base[1,],D_base[2,], xlab = "", ylab = "",cex.axis = 1.3,font = 2)
  for(i in 1:(k - 1)){
    for(j in (i + 1):k){
      tseg1 <- tropical.geodesic.dim.2(D_base[,i], D_base[,j])
      tseg2 <- tropical.geodesic.dim.2(D_base[,i], D_base[,j], flag=1)
      if(tseg1[[2]] < tseg2[[2]]) tseg <- tseg1
      else tseg <- tseg2
      segments(tseg[[1]][1,1], tseg[[1]][2,1], tseg[[1]][1,2], tseg[[1]][2,2], col = 'black')
      segments(tseg[[1]][1,2], tseg[[1]][2,2], tseg[[1]][1,3], tseg[[1]][2,3], col = 'black')
    }
  }
  points(x = proj_2D_plot_m[,2], y = proj_2D_plot_m[,3], pch = 16, cex = 0.75, col = col_value[freq+3])
  dev.off()
}

# comb_set is the set includes the index of tree combination detected by the tropical PCA method.
comb_set <- matrix(c(13167, 19321, 33,   15922, 19825, 9155,  3882,  9324,  1175,
                     5747,  15559, 4815, 19778, 405,   2588,  19927, 16297, 7182,
                     5868,  6428,  1296, 11240, 10231, 5653,  5892,  11225, 13277,
                     4819,  18173, 3093, 9003,  4515,  13271, 4436,  7701,  3558,
                     5092,  17025, 1529, 13321, 13956, 19531, 5050,  5020,  12878,
                     17233, 1701,  4199, 5921,  4394,  159,   15820, 9341,  2506,
                     10879, 13188, 3336, 182,   5680,  3106,  10665, 19722, 373), ncol = 3,byrow = T)


dataList <- list.files("./data/")
year <- c(1993:2013)
for(i in 1:21){
  trees_ori <- read.tree(paste("./data/",dataList[i],sep=""))
  tropTree_plot(trees_ori, comb_set[i,], year[i], pcs=3)
}

# Appendix A 
# Combination list
# 1: 1, 13167 19321 33
# 2: 2, 15922 19825 9155
# 3: 1, 3882 9324 1175
# 4: 2, 5747 15559 4815
# 5: 2, 19778 405 2588
# 6: 2, 19927 16297 7182
# 7: 1, 5868 6428 1296
# 8: 1, 11240 10231 5653
# 9: 5, 5892 11225 13277
# 10: 4, 4819 18173 3093
# 11: 1, 9003 4515 13271
# 12: 3, 4436 7701 3558
# 13: 5, 5092 17025 1529
# 14: 5, 13321 13956 19531
# 15: 1, 5050 5020 12878
# 16: 2, 17233 1701 4199
# 17: 2, 5921 4394 159
# 18: 1, 15820 9341 2506
# 19: 2, 10879 13188 3336
# 20: 1, 182 5680 3106
# 21: 3, 10665 19722 373
