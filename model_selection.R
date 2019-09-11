library(reshape2)
library(ggplot2)
library(plotROC)

# load evaluation functions
source("functions.R")


# generate data -----------------------------------------------------------
# set dimensions
n <- 150
p <- 300
q <- 300
cat("n:", n, "\np:", p, "\nq:", q, "\n")

# set sparsity level in L
sparsity <- 0.05
cat("sparsity:", sparsity, "\n")

# set numbers of iterations
niter <- 15000
nburn <- 7000

# set cutoff value for inclusion probabilities
thres <- 0.5

# generate true adjancency matrix
row_1 <- c()
col_1 <- c()
for (j in 1:(q-1)) {
  for (i in (j+1):q) {
    row_1 <- c(row_1, i)
    col_1 <- c(col_1, j)
  }
}
true.adj <- matrix(0, q, q)
nonzero_id <- sample(length(row_1), sparsity * length(row_1))
for(i in nonzero_id){
  true.adj[row_1[i], col_1[i]] <- 1
}

# generate Cholesky factor and covariance matrix
true.L <- true.adj
diag(true.L) <- 1
true.Sigma <- solve(true.L %*% t(true.L))

# generate true variable indicators
true.z <- c(FALSE, TRUE, FALSE, FALSE, TRUE, rep(FALSE, p-5))

# generate coefficient matrix
c <- 0.3^2
true.B <- matrix(0, p, q)
true.B[true.z, ] <- mvrnorm(sum(true.z), rep(0, q), c * true.Sigma)

# generate design matrix and response matrix
X <- matrix(rnorm(n * p, 10, 1), n, p)
Y <- t(apply(X %*% true.B, 1, function(x) mvrnorm(1, x, true.Sigma)))


# JDAG --------------------------------------------------------------------
# initialize
zout <- rep(0, p)
dagout <- matrix(0, q, q)
z <- rep(FALSE, p)
z[sample(p, 5)] <- TRUE
adj <- t(chol(solve(t(Y) %*% Y / n + 0.1 * diag(q))))
diag(adj) <- 0
adj <- 1 * (adj > 0.2)

# run MCMC
for (iter in 1:niter) {
  Xz <- X[, z]
  Hz <- solve(diag(n) + c * Xz %*% t(Xz))
  
  # get z candidate
  zz <- z
  if (runif(1) > 0.5) {
    zz[sample(1 * which(z), 1)] <- FALSE
    Xzz <- X[, zz]
    Hzz <- solve(diag(n) + c * Xzz %*% t(Xzz))
    pa <- min(1, exp(lognorm(adj, diag(q)+t(Y)%*%Hzz%*%Y, n+3)-lognorm(adj, diag(q)+t(Y)%*%Hz%*%Y, n+3)+(log(det(Hzz))-log(det(Hz)))*(q/2))*p)
  } else {
    zz[sample(1 * which(!z), 1)] <- TRUE
    Xzz <- X[, zz]
    Hzz <- solve(diag(n) + c*Xzz%*%t(Xzz))
    pa <- min(1, exp(lognorm(adj, diag(q)+t(Y)%*%Hzz%*%Y, n+3)-lognorm(adj, diag(q)+t(Y)%*%Hz%*%Y, n+3)+(log(det(Hzz))-log(det(Hz)))*(q/2))/p)
  }
  if (runif(1)<pa) {
    cat('update z\n')
    z <- zz
    Hz <- Hzz
  }
  if(iter > nburn) zout <- z + zout
  
  # get adj candidate
  adjj <- adj
  row1 <- c()
  col1 <- c()
  row0 <- c()
  col0 <- c()
  for (j in 1:(q-1)) {
    for (i in (j+1):q) {
      if (adj[i, j]>0) {
        row1 <- c(row1, i)
        col1 <- c(col1, j)
      } else {
        row0 <- c(row0, i)
        col0 <- c(col0, j)
      }
    }
  }
  if (runif(1) > 0.5) {
    idx <- sample(length(row1), 1)
    adjj[row1[idx], col1[idx]] <- 0
    pb <- min(1, exp(lognorm(adjj, diag(q)+t(Y)%*%Hz%*%Y, n+3)-lognorm(adj, diag(q)+t(Y)%*%Hz%*%Y, n+3))*q^2)
  } else {
    idx <- sample(length(row0), 1)
    adjj[row0[idx], col0[idx]] <- 1
    pb <- min(1, exp(lognorm(adjj, diag(q)+t(Y)%*%Hz%*%Y, n+3)-lognorm(adj, diag(q)+t(Y)%*%Hz%*%Y, n+3))/q^2)
  }
  
  # accept/reject adjj
  if (runif(1)<pb) {
    cat('update adj\n')
    adj <- adjj
  }
  if(iter > nburn) dagout <- dagout + adj
  
  if (iter %% 100 == 0) {
    cat(iter, 'iterations completed.')
  }
}

# get estimated DAG and variable indicators
est.dag <- 1 * as.matrix(dagout / (niter - nburn) > thres)
est.z <- 1 * as.matrix(zout / (niter - nburn) > thres)


# evaluate the results ----------------------------------------------------
evaluation.dag(true.adj, est.dag)
evaluation.var(true.z, est.z)


# plot heatmaps and ROC curves --------------------------------------------
# estimated DAG
pdf(paste("est.dag",n,p,q,sparsity,".pdf", sep = ""))
ggplot(data = melt(est.dag), aes(x=Var1, y=Var2, fill=value)) +
  geom_tile()
dev.off()

# true DAG
pdf(paste("true.dag",n,p,q,sparsity,".pdf", sep = ""))
ggplot(data = melt(true.adj), aes(x=Var1, y=Var2, fill=value)) +
  geom_tile()
dev.off()

# ROC for DAG selection
pdf(paste("roc.dag",n,p,q,sparsity,".pdf", sep = ""))
ggplot(data=data.frame(D=true.adj[lower.tri(true.adj)], M=est.dag[lower.tri(est.dag)]), aes(d=D, m=M)) + geom_roc(n.cuts=0)+ xlab("false positive rate") + ylab("true positive rate")+theme(axis.text=element_text(size=12),axis.title=element_text(size=18))
dev.off()

# ROC for variable selection
pdf(paste("roc.z",n,p,q,sparsity,".pdf",sep = ""))
ggplot(data=data.frame(D=1*true.z, M=est.z), aes(d=D, m=M)) + geom_roc(n.cuts=0)+ xlab("false positive rate") + ylab("true positive rate")+theme(axis.text=element_text(size=12),axis.title=element_text(size=18))
dev.off()
