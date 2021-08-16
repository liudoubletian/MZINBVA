#' An approach to generate multilevel data
#'
#' @param sub.n the number of subjects
#' @param pos.n the number of clusters
#' @param vis.n the number of visits
#' @param otu.n the number of taxa

#' @return data.mat  a list of data including covariates, count data
#' @examples
#' sub.n <- 10
#' pos.n <- 3
#' vis.n <- 3
#' otu.n <- 100
#' sim_dat <- sim.data(sub.n, pos.n, vis.n, otu.n)
#' @import Matrix
#' @import MASS
#' @export

sim.data <- function(sub.n, pos.n, vis.n, otu.n){
  set.seed(1)
  ID <-  sort(rep(seq(1:sub.n),each=pos.n*vis.n))
  cluster <- rep(rep(c(1:pos.n),each=vis.n),sub.n)
  
  x=matrix(rep(c(0,1),each=sub.n/2),ncol=1)
  x=rep(x,each=(pos.n*vis.n))
  x1=x2=as.matrix(data.frame(itc=1,x))
  
  w_l=list()
  for(i in 1:sub.n){
    w_l[[i]]=matrix(1,pos.n*vis.n,1)
  }
  w <- as.matrix(bdiag(w_l))
  
  
  v_i=matrix(0,ncol=pos.n,nrow=sub.n*pos.n*vis.n)
  for(k in 1:pos.n){
    v_i[which(cluster==k),k]=1
  }
  
  group_l=list()
  for(i in 1:sub.n){
    group_l[[i]]=v_i[(1+(i-1)*(pos.n*vis.n)):(i*(pos.n*vis.n)),]
  }
  
  
  v=as.matrix(bdiag(group_l))
  
  Y_mat <- matrix(NA,sub.n*pos.n*vis.n,otu.n)
  
  
  beta.mat <- matrix(c(0,0), 2,otu.n)
  gamma.mat <- matrix(c(0,0), 2,otu.n)
  
  beta.mat[1,] <- rnorm(otu.n,2,sd=1)
  gamma.mat[1,] <- rnorm(otu.n,0.5,sd=1)
  
  beta.mat[2,] <- runif(otu.n,-2,-0.1)
  gamma.mat[2,] <- runif(otu.n,0.1,2)
  
  
  for(j in 1:otu.n){
    set.seed(j)
    bi_sd <- runif(1,0.1,1)
    bi <- as.matrix(rnorm(sub.n,mean=0,sd=bi_sd))
    hij_sd <-  runif(1,0.1,1)
    hij <- as.matrix(rnorm(sub.n*pos.n,mean=0,sd=hij_sd))
    di_sd <- runif(1,0.1,1)
    di <- as.matrix(rnorm(sub.n,mean=0,sd=di_sd))
    fij_sd <- runif(1,0.1,1)
    fij <- as.matrix(rnorm(sub.n*pos.n,mean=0,sd=fij_sd))
    
    logit.p  <- x1 %*% gamma.mat[,j]+w%*%di+v%*%fij
    p  <- 1 / (1 + exp(-logit.p))
    ind.mix <- rbinom(sub.n*pos.n*vis.n, 1, p)
    
    lamda=exp(x2%*%matrix(beta.mat[,j])+w%*%bi+v%*%hij)
    for (i in 1:(sub.n*pos.n*vis.n)) {
      if (ind.mix[i] == 0) {
        phi <- runif(1,0.1,10)
        Y_mat[i,j] = rnbinom(1, mu=lamda[i],size=phi)
      }
      if (ind.mix[i] == 1) {
        Y_mat[i,j] = 0
      }
    }
  }
  
  data.mat <- list(w=w,v=v,x1=x,x2=x,Y_mat=Y_mat,ID=ID,cluster=cluster)
  return(data.mat)
  
}



