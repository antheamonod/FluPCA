library(lpSolveAPI)
library(parallel)
library(ape)
library(phangorn)

source("./program/func_ssh.R")
cl <- makeCluster(8) 

# Function to calculate proportion of variance explained (R^2)
# trees_ori: the original trees data, it should have a multiphy format
# comb_set:  the index of combination
# pcs:       number of vertices in tropical polytope 

r_square <- function(trees_ori, comb_set, pcs = 3,  outliers = NA){
  
  # Set negative branch lengths equal to 0
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
  
  new_base <- D_all[, comb_set]
  proj_points <- parLapply(cl, distVec_all, project_pi , D_s = new_base) 
  tropical_dist_vec <- mapply(tropical_dist, distVec_all, proj_points) 
  sum_dist <- sum(tropical_dist_vec) 
 
  r_proj_data <- matrix(unlist(proj_points), nrow = length(proj_points), byrow = T)
  
  apicom_fermat <- fermatweberdistance(r_proj_data)
  r_apicom <- apicom_fermat/(sum_dist + apicom_fermat)
  r_apicom
}
