#' A sandwich covariance matrix for model parameters
#'
#' @param res.list the list of data including the estimation of model parameters and variational parameters
#' @param nodes the number of nodes required for parallel
#' @return p_val the p-value for model parameters
#' @examples
#' sub.n <- 10
#' pos.n <- 3
#' vis.n <- 3
#' otu.n <- 100
#' sim_dat <- sim.data(sub.n, pos.n, vis.n, otu.n)
#' Y_mat <- sim_dat$Y_mat
#' data <- data.frame(Y=Y_mat[,1],ID=sim_dat$ID,cluster=sim_dat$cluster,x1=sim_dat$x1,x2=sim_dat$x2)
#' est_m <- step_alg(data)
#' res.t <- data.pro(est_m,nodes=3)
#' @export
#' @import parallel
#' @import MASS




total_elbo_split <- function(beta,data,gamma,phi,psib,psig,lamdab,lamdah,lamdad,lamdaf,pijk,nodes,offset){
  Y <- data$y;ID <-data$ID;cluster <-data$cluster;x1 <-data$x1;x2 <-data$x2;
  p1 <- length(unique(ID));p2 <- length(unique(cluster))
  mubi=psib[p1]
  muhij=psib[(1+p1):(p1+p2)]
  sigmabi=psib[(p1+p2+1)]
  sigmahij=psib[(p1+p2+2):(p1+p2*2+1)]
  
  mudi=psig[p1]
  mufij=psig[(1+p1):(p1+p2)]
  sigmadi=psig[(p1+p2+1)]
  sigmafij=psig[(p1+p2+2):(p1+p2*2+1)]
  
  new.mubi <- rep(mubi,length(ID))
  new.muhij <- muhij[cluster]
  new.sigmabi <- rep(sigmabi,length(ID))
  new.sigmahij <- sigmahij[cluster]
  
  new.mudi <- rep(mudi,length(ID))
  new.mufij <- mufij[cluster]
  new.sigmadi <-  rep(sigmadi,length(ID))
  new.sigmafij <- sigmafij[cluster]
  
  eq2=x2%*%beta+new.mubi+new.muhij+offset
  eq1=x1%*%gamma+new.mudi+new.mufij
  B1=(exp(x2%*%beta+new.mubi+new.sigmabi^2/2+new.muhij+new.sigmahij^2/2+offset))
  C1=(exp(-x2%*%beta-new.mubi+new.sigmabi^2/2-new.muhij+new.sigmahij^2/2-offset))
  F1=(exp(x1%*%gamma+new.mudi+(new.sigmadi^2)/2+new.mufij+(new.sigmafij^2)/2))
  
  eq11=-x1%*%gamma-new.mudi-new.mufij
  D1=-phi*log(1+B1/phi)
  pijk[Y==0]=(1/(1+exp(eq11+D1)))[Y==0]
  
  ini1=sum(pijk*eq1-log(1+F1))
  ini2=(1-pijk)*(lfactorial(Y+phi-1)-lfactorial(Y)-log(gamma(phi))-(Y+phi)*log(1+phi*C1)+
                 phi*log(phi)-phi*eq2)
  ini3=(1-pijk)*((-phi)*log(1+B1/phi))
  
  
  elbo_c=log(sigmabi/lamdab)-0.5*((mubi/lamdab)^2+(sigmabi/lamdab)^2)+0.5+
    sum(log(sigmahij/lamdah))-0.5*sum((muhij/lamdah)^2+(sigmahij/lamdah)^2)+0.5*length(sigmahij)+
    log(sigmadi/lamdad)-0.5*((mudi/lamdad)^2+(sigmadi/lamdad)^2)+0.5+
    sum(log(sigmafij/lamdaf))-0.5*sum((mufij/lamdaf)^2+(sigmafij/lamdaf)^2)+0.5*length(sigmahij)-
    sum(pijk[pijk!=0]*log(pijk[pijk!=0]))-sum((1-pijk[pijk!=1])*log(1-pijk[pijk!=1]))
  
  elbo=ini1+sum(ifelse(Y>0,ini2,ini3))+elbo_c
  return(-elbo)
  
}




one_elbo_theta=function(theta,data,gamma,lamdad,lamdaf,psib,psig,pijk,nodes,offset){
  Y <- data$y;ID <-data$ID;cluster <-data$cluster;x1 <-data$x1;x2 <-data$x2
  p1 <- ncol(x1);
  beta=theta[1:p1];phi=theta[(p1+1)];
  lamdab=theta[(p1+2)];lamdah=theta[(p1+3)];
  total_elbo_split(beta=beta,phi=phi,lamdab=lamdab,lamdah=lamdah,data=data,gamma=gamma,psib=psib,psig=psig,lamdad=lamdad,lamdaf=lamdaf,pijk=pijk,nodes=nodes,offset)
}



one_elbo_thetapsi=function(thetapsi,psig,data,pijk,gamma,lamdad,lamdaf,nodes,offset){
  Y <- data$y;ID <-data$ID;cluster <-data$cluster;x1 <-data$x1;x2 <-data$x2
  p1 <- ncol(x1);
  beta=thetapsi[1:p1];phi=thetapsi[(p1+1)];
  lamdab=thetapsi[(p1+2)];lamdah=thetapsi[(p1+3)];
  psib <- thetapsi[-(1:(p1+3))]
  
  total_elbo_split(beta=beta,phi=phi,lamdab=lamdab,lamdah=lamdah,data=data,gamma=gamma,psib=psib,psig=psig,lamdad=lamdad,lamdaf=lamdaf,pijk=pijk,nodes=nodes,offset)
}


sandwich_cov <- function(beta,data,phi,gamma,lamdab,lamdah,lamdad,lamdaf,psib,psig,pijk,verbose=FALSE,nodes,offset) {
  require(numDeriv)
  require(parallel)
  N <- length(data)
  theta <- c( beta,phi,lamdab,lamdah)
  theta_ln <- length(theta)
  abs <- mclapply(1:N, function(i) {
    theta <- c( beta,phi,lamdab,lamdah)
    theta_ln <- length(theta)
    thetapsi <- c(theta,psib=psib[[i]])
    psi_ln <- length(thetapsi)-theta_ln
    
    grad <- grad(one_elbo_theta,theta,data=data[[i]],gamma=gamma,lamdad=lamdad,lamdaf=lamdaf,psib=psib[[i]],psig=psig[[i]],pijk=pijk[[i]],nodes=nodes,offset=offset[[i]])
    hess <- hessian(one_elbo_thetapsi,thetapsi,data=data[[i]],gamma=gamma,lamdad=lamdad,lamdaf=lamdaf,psig=psig[[i]],pijk=pijk[[i]],nodes=nodes,offset=offset[[i]])
    
    Ai <- hess[1:theta_ln,1:theta_ln] - hess[1:theta_ln, -(1:theta_ln)] %*% solve(hess[(theta_ln + 1):(theta_ln + psi_ln), (theta_ln + 1):(theta_ln + psi_ln)]) %*% hess[-(1:theta_ln), 1:theta_ln]
    Bi <- tcrossprod(grad) 
    c(c(Ai), c(Bi))
  }, mc.cores=nodes)
  
  abs <- matrix(unlist(abs), nrow=length(abs), byrow=TRUE)
  abhat <- colMeans(abs)
  ahat <- matrix(abhat[1:theta_ln^2], nrow=theta_ln)
  bhat <- matrix(abhat[-(1:theta_ln^2)], nrow=theta_ln)
  
  return(list(ahat=ahat, bhat=bhat, sand_cov=solve(ahat) %*% bhat %*% solve(ahat) / N))
}


data.pro <- function(kk,nodes=3){
  library(parallel)
  library(MASS)
  data <- data.frame(y=kk$Y,ID=kk$ID,cluster=kk$cluster)
  N <- length(unique(data$ID))
  offset <- kk$offset
  x1 <- as.data.frame(kk$x1)
  x2 <- as.data.frame(kk$x2)
  beta <- kk$beta; gamma <- kk$gamma; phi <- kk$phi;
  lamdab <- kk$lamdab;lamdah <- kk$lamdah;lamdad <- kk$lamdad;lamdaf <- kk$lamdaf
  mubi <- kk$mubi;muhij <- kk$muhij; sigmabi <- kk$sigmabi; sigmahij <- kk$sigmahij; 
  mudi <- kk$mudi;mufij <- kk$mufij; sigmadi <- kk$sigmadi; sigmafij <- kk$sigmafij;
  pijk <- kk$pijk;numk <- kk$num_k
  uniq_id <- unique(ID)
  
  data.split <- lapply(1:N, function(j) {
    sub_data <- subset(data, data$ID == uniq_id[j])
    y <- sub_data$y
    x1 <- as.matrix(cbind(1,x1[which(data$ID == uniq_id[j]),]))
    x2 <- as.matrix(cbind(1,x2[which(data$ID == uniq_id[j]),]))
    list(x1=x1,x2=x2, y=y,ID=sub_data$ID,cluster=sub_data$cluster)
  })
  
  pijk.split <- lapply(1:N, function(j) {
    pijk[which(data$ID==uniq_id[j])]
  })
  
   off.split <- lapply(1:N, function(i) {
    offset[which(data$ID==i)]
  })
  
  
  
  cl.cu <- unlist(lapply(1:N, function(j) {
    sub_data <- subset(data, data$ID == uniq_id[j])
    length(unique(sub_data$cluster))
  }))
  
  
  
  psib.split <- list()
  len <- 0
  for(i in 1:N) {
    psib.split[[i]] <- c(mubi[i],muhij[(sum(len)+1):(sum(len)+cl.cu[i])],sigmabi[i],sigmahij[(sum(len)+1):(sum(len)+cl.cu[i])])
    len <- append(len,cl.cu[i])
  }
  
  psig.split <- list()
  len <- 0
  for(i in 1:N) {
    psig.split[[i]] <- c(mudi[i],mufij[(sum(len)+1):(sum(len)+cl.cu[i])],sigmadi[i],sigmafij[(sum(len)+1):(sum(len)+cl.cu[i])])
    len <- append(len,cl.cu[i])
  }
  
  dd <- sandwich_cov(beta,data=data.split,phi,gamma,lamdab,lamdah,lamdad,lamdaf,psib=psib.split,psig=psig.split,pijk=pijk.split,verbose=FALSE,nodes,offset=off.split))
  
  cov.inf <- dd$sand_cov[2,2]
  tes_val <- t(matrix(beta[2]))%*%ginv(cov.inf)%*%matrix(beta[2])
  b_val_san <- 1-pchisq(as.numeric(tes_val),df=1)
  
  return(b_val_san)
}
