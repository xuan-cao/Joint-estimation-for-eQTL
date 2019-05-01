library(MASS)


GetDiagSub <- function(ind,U){
  return(U[ind,ind])
}

SpecialDet <- function(x){
  if(length(x)==0) return(1)
  if(length(x)==1) return(as.numeric(x))
  return(det(x))
}

lognorm <- function(Adj, U, alpha) {
  m <- nrow(Adj)
  if(m!=ncol(Adj)) stop("Adjacent matrix is not square!")
  pa.set <- apply(Adj,2,function(x){which(x!=0)})
  pa.count <- unlist(lapply(pa.set, length))
  alpha <- pa.count + alpha
  if(length(pa.set)>0){
    pa.count <- unlist(lapply(pa.set,length))
    if(sum((alpha-2)<=pa.count)>0) stop("The alpha must satisfy alpha_i > pa_i+2!")
    pa.mat <- lapply(pa.set,GetDiagSub,U=U)
    pa.mat.det <- unlist(lapply(pa.mat,SpecialDet))
    fa.set <- lapply(1:m,function(i){return(c(i,pa.set[[i]]))})
    fa.mat <- lapply(fa.set,GetDiagSub,U=U)
    #print(fa.mat)
    fa.mat.det <- unlist(lapply(fa.mat,SpecialDet))
  } else{
    pa.count <- rep(0,m)
    if(sum((alpha-2)<=pa.count)>0) stop("The alpha must satisfy alpha_i > pa_i+2!")
    pa.mat.det <- rep(1,m)
    fa.set <- lapply(1:m, function(i)return(i))
    fa.mat <- lapply(fa.set,GetDiagSub,U=U)
    fa.mat.det <- unlist(lapply(fa.mat,SpecialDet))
    
  }
  logz <- sum(lgamma(alpha/2-pa.count/2-1) + log(2)*(alpha/2-1) + log(pi)*pa.count/2 + log(pa.mat.det)*(alpha-pa.count-3)/2 - log(fa.mat.det)*(alpha-pa.count-1)/2)
  return(logz)
}

evaluation.dag <- function(Adj1, Adj2){
  #Adj1 <- as.matrix(get.adjacency(graph1))
  #Adj2 <- as.matrix(get.adjacency(graph2))
  true.index <- which(Adj1==1)
  #false.index <- setdiff(which(upper.tri(Adj1)),true.index)
  false.index <- which(Adj1==0)
  positive.index <- which(Adj2==1)
  # negative.index <- setdiff(which(upper.tri(Adj2)),positive.index)
  negative.index <- which(Adj2==0)
  
  TP <- length(intersect(true.index,positive.index))
  FP <- length(intersect(false.index,positive.index))
  FN <- length(intersect(true.index,negative.index))
  TN <- length(intersect(false.index,negative.index))
  
  MCC.denom <- sqrt(TP+FP)*sqrt(TP+FN)*sqrt(TN+FP)*sqrt(TN+FN)
  if(MCC.denom==0) MCC.denom <- 1
  MCC <- (TP*TN-FP*FN)/MCC.denom
  if((TN+FP)==0) MCC <- 1
  
  Precision <- TP/(TP+FP)
  if((TP+FP)==0) Precision <- 1
  Recall <- TP/(TP+FN)
  if((TP+FN)==0) Recall <- 1
  Sensitivity <- Recall
  Specific <- TN/(TN+FP)
  if((TN+FP)==0) Specific <- 1
  
  #graph1.graph <- igraph.to.graphNEL(graph1)
  #graph2.graph <- igraph.to.graphNEL(graph2)
  #SHD <- shd(graph1.graph,graph2.graph)
  return(list(Precision=Precision,Recall=Recall,Sensitivity=Sensitivity,Specific=Specific,MCC=MCC,TP=TP,FP=FP,TN=TN,FN=FN))
}

evaluation.var <- function(v1,v2){
  true.id <- which(v1==TRUE)
  false.id <- which(v1==FALSE)
  positive.id <- which(v2!=0)
  negative.id <- which(v2==0)
  
  TP <- length(intersect(true.id,positive.id))
  FP <- length(intersect(false.id,positive.id))
  FN <- length(intersect(true.id,negative.id))
  TN <- length(intersect(false.id,negative.id))
  

  
  Precision <- TP/(TP+FP)
  if((TP+FP)==0) Precision <- 1
  Recall <- TP/(TP+FN)
  if((TP+FN)==0) Recall <- 1
  Sensitivity <- Recall
  Specific <- TN/(TN+FP)
  if((TN+FP)==0) Specific <- 1
  
  #graph1.graph <- igraph.to.graphNEL(graph1)
  #graph2.graph <- igraph.to.graphNEL(graph2)
  #SHD <- shd(graph1.graph,graph2.graph)
  return(list(Precision=Precision,Recall=Recall,Sensitivity=Sensitivity,Specific=Specific,TP=TP,FP=FP,TN=TN,FN=FN))
}
  
