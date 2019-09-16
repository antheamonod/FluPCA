# This file contains supplemental functions required to run BHV and tropical PCA, and calculate their respective proportions of variance explained (R-squared values)
# Authors: Qiwen Kang, Xu Zhang

project_pi<-function(D_s,D){
  if(is.null(dim(D_s))){
    lambda <- min(D - D_s) 
    pi_D <- c(t(lambda + t(D_s))) 
  }else{
    lambda <- apply(D - D_s, 2, min)#D_s by row
    pi_D <- apply(t(lambda + t(D_s)),1,max)
  }
  return(pi_D)
}

vec_fun<-function(x){
  m<-dim(x)[1]
  vecTreesVec<-rep(NA,choose(m,2)) 
  for(row.num in 1:(m-1)){
    for(col.num in (row.num+1):m){
      vecTreesVec[col.num-row.num+(m-1+(m-1-row.num+2))*(row.num-1)/2]<-x[row.num,col.num]  
    }
  }
  vecTreesVec
}

# Calculate the tropical distance of two points
#` D_1:  distance vector of tree 1
#` D_2:  distance vector of tree 2

tropical_dist<-function(D_1,D_2){
  e <- length(D_1)
  t_dist <- 0
  for(i in 1:(e-1)){
    for(j in (i+1):e){
      if(abs(D_1[i]-D_2[i]-D_1[j]+D_2[j])>t_dist){
        t_dist<-abs(D_1[i]-D_2[i]-D_1[j]+D_2[j])
      }
    }
  }
  t_dist
}

# Calculate the the sum of tropical distances of a set of trees
#`     pc_base:  distance vector usded to build tropical space
#` distVec_all:  all distance vectors of the whole trees

tropDist <- function(pc_base, distVec_all){
  
  proj_points <- parLapply(cl, distVec_all, project_pi , D_s = pc_base)
  tropical_dist_vec <- mapply(tropical_dist, distVec_all, proj_points)
  sum_dist <- sum(tropical_dist_vec)
  
  
  return(sum_dist)
}

# Transform the projected point vectors to plot them in 2 dimensions
#` D: of size (s by e) whose rows are vertices of tropical polytope, point 
#` P: in rowspan(D)
polytope_iso<-function(D, P){
  e = length(P)
  s = dim(D)[[1]]
  Q = mat.or.vec(1, s)
  for (i in seq(s)){
    maxvalue = D[i,1] - P[[1]]
    for (j in seq(e)){
      maxvalue = max(maxvalue, D[i,j] - P[[j]])
    }
    Q[[i]]=maxvalue
  }
  return(Q)
}

# Normalize the transformed vectors
# D: transfered projected points vector
normalize.proj<-function(D){
  r<-length(D)
  D.new<-rep(NA,r)
  for(i in 1:r){
    D.new[i]<-D[i] - D[1]
  }
  return(D.new)
}

# Calculate the distance matrix of the input tree data set 
#`    trees:  the input tree data set.
#` tipOrder:  the order of leasves.

distMat <- function(trees, tipOrder){
  if(class(trees)=="multiPhylo"){
    trees_root <- root(trees, outgroup = tipOrder[1],resolve.root=TRUE)
    
    chronotrees <- parLapply(cl, trees_root, chronos)
    dist_chrono <- parLapply(cl, chronotrees,cophenetic)

    dist_ordered <- parLapply(cl, dist_chrono, function(x) x[tipOrder, tipOrder])
    distVec_all <- parLapply(cl, dist_ordered,vec_fun)
    
  }else {
    treeOne <- root(trees, outgroup = tipOrder[1],resolve.root=TRUE)
    chronoTree <- chronos(treeOne)
    dist_chrono_one <- cophenetic(chronoTree)
    
    dist_ordered_one <- dist_chrono_one[tipOrder, tipOrder]
    distVec_all <- vec_fun(dist_ordered_one)
  }
  
  return(distVec_all)
}

# Normalize the ultrametric vector of a tree
#` D: the input ultrametric vector of a tree
normalize.ultrametrices <- function(D){
  k <- ncol(D)
  new.D <- matrix(rep(0, 2*k), nrow=2, ncol=k)
  for(i in 2:3)
    new.D[i-1, ] <- D[i, ] - D[1, ]
  return(new.D)
}

tropical.geodesic.dim.2 <- function(D1, D2, flag = 0){
  k <- length(D1)
  if(k != 2) warning("dimension has to be 2!")
  for(i in 1:k)
    D1[i] <- round(D1[i], 4)
  for(i in 1:k)
    D2[i] <- round(D2[i], 4)
  if(length(D2) != k)
    warning("dimension is wrong!")
  addd <- 0
  if(flag == 1){
    tmp.D <- D2
    D2 <- D1
    D1 <- tmp.D
  }
  tmp.metric <- (D2 - D1)
  sorted.tmp.metric <- sort.int(tmp.metric, index.return=TRUE)
 
  D <- rep(0, k)
  
  D[sorted.tmp.metric$ix[2]] <- D2[sorted.tmp.metric$ix[2]]
  D[sorted.tmp.metric$ix[1]] <- min(D2[sorted.tmp.metric$ix[2]] - D1[sorted.tmp.metric$ix[2]] + D1[sorted.tmp.metric$ix[1]], D1[sorted.tmp.metric$ix[1]])
  
  distance <- max(abs(D1 - D))
  distance <- distance + max(abs(D2 - D))
  
  segment <- matrix(rep(0, 6), nrow=2, ncol=3)
  segment[,1] <- D1
  segment[,2] <- D
  segment[,3] <- D2
  
  return(list(segment, distance))
}

# Calculate R-squared on lprec objects 
#` datamatrix: projected porints vectors
fermatweberdistance <- function(datamatrix) {
  n = dim(datamatrix)[1]
  m = dim(datamatrix)[2]
  lprec <- make.lp(0, n+m)
  objective = mat.or.vec(n+m,1)
  for (i in seq(n)) {
    objective[i] = 1
  }
  set.objfn(lprec, objective)
  for (i in seq(n)) {
    for (j in seq(m)) {
      for (k in seq(m)) {
        v = mat.or.vec(n+m,1)
        v[i] = 1
        v[n+k] = 1
        v[n+j] = -1
        add.constraint(lprec, v, ">=", datamatrix[i,k] - datamatrix[i,j])
      }
    }
  }
  solve(lprec)
  return((get.objective(lprec)))
}

# Transform the list of distance vectors to distance matrix
#`    D: list of distance vectors
#`    n: number of leaves
#` tips: order of leaves

make.matrix <- function(D, n, tips){
  dd <- matrix(rep(0, n*n), nrow=n, ncol=n, byrow=TRUE)
  count <- 1
  for(i in 1:(n-1))
    for(j in (i + 1):n){
      dd[i, j] <- D[count]
      dd[j, i] <- dd[i, j]
      count <- count+1
    }
  mymatrix <- matrix(dd, nrow=n, ncol=n, byrow=TRUE, dimnames=list(tips,tips))
  return(mymatrix)
}
