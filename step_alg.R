#' An algrithm for variational parameters estimation of multilevel zero inflated negative binomial model
#'
#' @param data the multilevel data
#' @return beta the model paraeters for negative binomial
#' @return gamma the model paraeters for logistic regression
#' @examples
#' sub.n <- 10
#' pos.n <- 3
#' vis.n <- 3
#' otu.n <- 100
#' sim_dat <- sim.data(sub.n, pos.n, vis.n, otu.n)
#' Y_mat <- sim_dat$Y_mat
#' sim_dat["Y_mat"] <- NULL
#' data <- c(sim_dat, list(Y=Y_mat[,1]))
#' est_m <- step_alg(data)
#' @export




step_alg=function(data,trace=FALSE){
  n.init=n.i=1; eps=1e-3; max.iter=100;
  current.loglik <- -1e6; iter <- 1; err <- 10; div=1e5
  cur.mubi=cur.muhij=cur.mudi=cur.mufij= -1e6
  cur.sigmabi=cur.sigmahij=cur.sigmadi=cur.sigmafij= -1e6
  cur.bet=cur.gam=cur.phi= -1e6
#######observe data###  
  Y <- data$Y;
  x1 <- data$x1
  x2 <- data$x2
  w <- data$w
  v <- data$v
  ID <- data$ID
  cluster  <- data$cluster
  
  num_l <- length(unique(ID))
  num_k <- sum(unlist(lapply(1:num_l,function(x){length(unique(cluster[ID==x]))})))
  num_p <- length(Y)
  l <- ncol(x1)
  k <- ncol(x2)
#####parameter initialize ####  

  beta=new.beta=rep(0, ncol(x1))
  gamma=new.gamma=rep(0, ncol(x2))
  new.phi=phi=0.01
  new.mubi=mubi= rep(0.01,num_l)
  new.muhij=muhij=matrix(rep(0.01,num_k))
  new.mudi=mudi=rep(0.01,num_l)
  new.mufij=mufij=rep(0.01,num_k)
  new.sigmabi=sigmabi=rep(0.01,num_l)
  new.sigmahij=sigmahij=rep(0.01,num_k)
  new.sigmadi=sigmadi=rep(0.01,num_l)
  new.sigmafij=sigmafij=rep(0.01,num_k)
  new.lamdab=lamdab=1
  new.lamdah=lamdah=1
  new.lamdad=lamdad=1
  new.lamdaf=lamdaf=1
  new.pijk=pijk=ifelse(Y==0,mean(Y==0),0)
  
###################
  loglik=c()
  while((div> eps*(abs(current.loglik)+eps)) && iter <= max.iter) {
    
    
    delta.sigmahij.required <- 1e-3; p.iter <- 1; p.max.iter <- 10; delta.sigmahij<-sigmahij ;
    while(!all(delta.sigmahij < delta.sigmahij.required) & (p.iter < p.max.iter)){
      
      
      
      
      
      
      
      q=try(nlminb(start =as.vector(mubi), total_mubi_neg , gradient = grad_mubi_neg,lower = rep(-Inf,num_l), upper = rep(Inf,num_l),
                   Y=Y,x1=x1,x2=x2,v=v,w=w,lamdab=new.lamdab,lamdah=new.lamdah,lamdad=new.lamdad,lamdaf=new.lamdaf,
                   phi=new.phi,pijk=new.pijk,gamma=new.gamma,beta=new.beta,muhij=new.muhij,mudi=new.mudi,
                   mufij=new.mufij,sigmabi=new.sigmabi,sigmahij=new.sigmahij,sigmadi=new.sigmadi,
                   sigmafij=new.sigmafij,
                   num_l=num_l,num_k=num_k,l=l,k=k,control = list(trace = 0, iter.max = 100)), silent = TRUE)
      
      
      if(!inherits(q, "try-error")) {
        if(q$convergence != 0) { if(trace) cat("Optimization algorithm did not converge when it tried to update variational parameters on step ", iter,"\n") }
        if(p.iter > 1 && cur.mubi > -q$objective) { if(trace) cat("Optimization did not improve estimates of variational parameters.\n"); new.mubi <- mubi }
        else {
          if(trace) cat("Variational parameters updated", "\n")
          new.mubi <- q$par
        }
      }
      else { new.mubi <- mubi}
      
      cur.mubi<- -q$objective
      
      q=try(nlminb(start =as.vector(sigmabi), total_sigmabi_neg , gradient = grad_sigmabi_neg,lower = rep(-Inf,num_l), upper = rep(Inf,num_l),
                   Y=Y,x1=x1,x2=x2,v=v,w=w,lamdab=new.lamdab,lamdah=new.lamdah,lamdad=new.lamdad,lamdaf=new.lamdaf,
                   phi=new.phi,pijk=new.pijk,gamma=new.gamma,beta=new.beta,mubi=new.mubi,muhij=new.muhij,mudi=new.mudi,
                   mufij=new.mufij,sigmahij=new.sigmahij,sigmadi=new.sigmadi,
                   sigmafij=new.sigmafij,
                   num_l=num_l,num_k=num_k,l=l,k=k,control = list(trace = 0, iter.max = 100)), silent = TRUE)
      
      
      if(!inherits(q, "try-error")) {
        if(q$convergence != 0) { if(trace) cat("Optimization algorithm did not converge when it tried to update variational parameters on step ", iter,"\n") }
        if(p.iter > 1 && cur.sigmabi > -q$objective) { if(trace) cat("Optimization did not improve estimates of variational parameters.\n"); new.sigmabi <- sigmabi }
        else {
          if(trace) cat("Variational parameters updated", "\n")
          new.sigmabi <- q$par
        }
      }
      else { new.sigmabi <- sigmabi}
      
      cur.sigmabi<- -q$objective
      
      q=try(nlminb(start =as.vector(muhij), total_muhij_neg , gradient = grad_muhij_neg,lower = rep(-Inf,num_k), upper = rep(Inf,num_k),
                   Y=Y,x1=x1,x2=x2,v=v,w=w,lamdab=new.lamdab,lamdah=new.lamdah,lamdad=new.lamdad,lamdaf=new.lamdaf,
                   phi=new.phi,pijk=new.pijk,gamma=new.gamma,beta=new.beta,mubi=new.mubi,mudi=new.mudi,
                   mufij=new.mufij,sigmabi=new.sigmabi,sigmahij=new.sigmahij,sigmadi=new.sigmadi,
                   sigmafij=new.sigmafij,
                   num_l=num_l,num_k=num_k,l=l,k=k,control = list(trace = 0, iter.max = 100)), silent = TRUE)
      
      
      if(!inherits(q, "try-error")) {
        if(q$convergence != 0) { if(trace) cat("Optimization algorithm did not converge when it tried to update variational parameters on step ", iter,"\n") }
        if(p.iter > 1 && cur.muhij > -q$objective) { if(trace) cat("Optimization did not improve estimates of variational parameters.\n"); new.muhij <- muhij }
        else {
          if(trace) cat("Variational parameters updated", "\n")
          new.muhij <- q$par
        }
      }
      else { new.muhij <- muhij}
      
      cur.muhij<- -q$objective
      
      
      
      
      
      q=try(nlminb(start =as.vector(sigmahij), total_sigmahij_neg , gradient = grad_sigmahij_neg,lower = rep(-Inf,num_k), upper = rep(Inf,num_k),
                   Y=Y,x1=x1,x2=x2,v=v,w=w,lamdab=new.lamdab,lamdah=new.lamdah,lamdad=new.lamdad,lamdaf=new.lamdaf,
                   phi=new.phi,pijk=new.pijk,gamma=new.gamma,beta=new.beta,mubi=new.mubi,muhij=new.muhij,mudi=new.mudi,
                   mufij=new.mufij,sigmabi=new.sigmabi,sigmadi=new.sigmadi,
                   sigmafij=new.sigmafij,
                   num_l=num_l,num_k=num_k,l=l,k=k,control = list(trace = 0, iter.max = 100)), silent = TRUE)
      
      
      if(!inherits(q, "try-error")) {
        if(q$convergence != 0) { if(trace) cat("Optimization algorithm did not converge when it tried to update variational parameters on step ", iter,"\n") }
        if(p.iter > 1 && cur.sigmahij > -q$objective) { if(trace) cat("Optimization did not improve estimates of variational parameters.\n"); new.sigmahij <- sigmahij }
        else {
          if(trace) cat("Variational parameters updated", "\n")
          new.sigmahij <- q$par
        }
      }
      else { new.sigmahij <- sigmahij}
      cur.sigmahij<- -q$objective
      
      new.lamdab=sqrt(sum(new.sigmabi^2+new.mubi^2)/num_l)
      
      
      new.lamdah=sqrt(sum(new.sigmahij^2+new.muhij^2)/(num_k))
      
      
      
      
      
      
      q=try(nlminb(start =as.vector(beta), elbo_beta , gradient = grad_beta,lower = rep(-Inf,l), upper = rep(Inf,l),
                   Y=Y,x1=x1,x2=x2,v=v,w=w,lamdab=new.lamdab,lamdah=new.lamdah,lamdad=new.lamdad,lamdaf=new.lamdaf,
                   phi=new.phi,pijk=new.pijk,gamma=new.gamma,mubi=new.mubi,muhij=new.muhij,mudi=new.mudi,
                   mufij=new.mufij,sigmabi=new.sigmabi,sigmahij=new.sigmahij,
                   sigmadi=new.sigmadi,sigmafij=new.sigmafij,num_l=num_l,num_k=num_k,l=l,k=k,control = list(trace = 0, iter.max = 100)), silent = TRUE)
      
      
      
      if(!inherits(q, "try-error")) {
        if(q$convergence != 0) { if(trace) cat("Optimization algorithm did not converge when it tried to update variational parameters on step ", iter,"\n") }
        if(p.iter > 1 && cur.bet > -q$objective) { if(trace) cat("Optimization did not improve estimates of variational parameters.\n"); new.beta <- beta }
        else {
          if(trace) cat("Variational parameters updated", "\n")
          new.beta <- q$par
        }
      }
      else { new.beta <- beta}
      
      cur.bet<- -q$objective
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      q=try(nlminb(start =as.vector(phi), elbo_phi , gradient = grad_phi,lower = rep(-Inf,1), upper = rep(Inf,1),
                   Y=Y,x1=x1,x2=x2,v=v,w=w,lamdab=new.lamdab,lamdah=new.lamdah,lamdad=new.lamdad,lamdaf=new.lamdaf,
                   gamma=new.gamma,pijk=new.pijk,beta=new.beta,mubi=new.mubi,muhij=new.muhij,mudi=new.mudi,
                   mufij=new.mufij,sigmabi=new.sigmabi,sigmahij=new.sigmahij,
                   sigmadi=new.sigmadi,sigmafij=new.sigmafij,num_l=num_l,num_k=num_k,l=l,k=k,control = list(trace = 0, iter.max = 100)), silent = TRUE)
      
      
      
      if(!inherits(q, "try-error")) {
        if(q$convergence != 0) { if(trace) cat("Optimization algorithm did not converge when it tried to update variational parameters on step ", iter,"\n") }
        if(p.iter > 1 && cur.phi > -q$objective) { if(trace) cat("Optimization did not improve estimates of variational parameters.\n"); new.phi <- phi }
        else {
          if(trace) cat("Variational parameters updated", "\n")
          new.phi <- q$par
        }
      }
      else { new.phi <- phi}
      #
      #
      cur.phi<- -q$objective
      
      
      
      beta=new.beta;mubi=new.mubi;muhij=new.muhij;phi=new.phi;lamdab=new.lamdab;lamdah=new.lamdah;
      sigmabi=new.sigmabi;sigmahij=new.sigmahij;
      
      
      
      delta.sigmahij <- abs(new.sigmahij-sigmahij)
      sigmahij <- new.sigmahij
      p.iter <- p.iter+1
    }
    
    
    
    delta.sigmafij.required <- 1e-3; p.iter <- 1; p.max.iter <- 10; delta.sigmafij<-sigmafij ;
    while(!all(delta.sigmafij < delta.sigmafij.required) & (p.iter < p.max.iter)){
      
      q=try(nlminb(start =as.vector(mudi), total_mudi_neg , gradient = grad_mudi_neg,lower = rep(-Inf,num_l), upper = rep(Inf,num_l),
                   Y=Y,x1=x1,x2=x2,v=v,w=w,lamdab=new.lamdab,lamdah=new.lamdah,lamdad=new.lamdad,lamdaf=new.lamdaf,
                   phi=new.phi,pijk=new.pijk,gamma=new.gamma,beta=new.beta,mubi=new.mubi,muhij=new.muhij,
                   mufij=new.mufij,sigmabi=new.sigmabi,sigmahij=new.sigmahij,sigmadi=new.sigmadi,
                   sigmafij=new.sigmafij,
                   num_l=num_l,num_k=num_k,l=l,k=k,control = list(trace = 0, iter.max = 100)), silent = TRUE)
      
      
      if(!inherits(q, "try-error")) {
        if(q$convergence != 0) { if(trace) cat("Optimization algorithm did not converge when it tried to update variational parameters on step ", iter,"\n") }
        if(p.iter > 1 && cur.mudi > -q$objective) { if(trace) cat("Optimization did not improve estimates of variational parameters.\n"); new.mudi <- mudi}
        else {
          if(trace) cat("Variational parameters updated", "\n")
          new.mudi <- q$par
        }
      }
      else { new.mudi<- mudi}
      
      cur.mudi<- -q$objective
      
      
      q=try(nlminb(start =as.vector(sigmadi), total_sigmadi_neg , gradient = grad_sigmadi_neg,lower = rep(-Inf,num_l), upper = rep(Inf,num_l),
                   Y=Y,x1=x1,x2=x2,v=v,w=w,lamdab=new.lamdab,lamdah=new.lamdah,lamdad=new.lamdad,lamdaf=new.lamdaf,
                   phi=new.phi,pijk=new.pijk,gamma=new.gamma,beta=new.beta,mubi=new.mubi,muhij=new.muhij,mudi=new.mudi,
                   mufij=new.mufij,sigmabi=new.sigmabi,sigmahij=new.sigmahij,
                   sigmafij=new.sigmafij,
                   num_l=num_l,num_k=num_k,l=l,k=k,control = list(trace = 0, iter.max = 100)), silent = TRUE)
      
      
      if(!inherits(q, "try-error")) {
        if(q$convergence != 0) { if(trace) cat("Optimization algorithm did not converge when it tried to update variational parameters on step ", iter,"\n") }
        if(p.iter > 1 && cur.sigmadi > -q$objective) { if(trace) cat("Optimization did not improve estimates of variational parameters.\n"); new.sigmadi <- sigmadi }
        else {
          if(trace) cat("Variational parameters updated", "\n")
          new.sigmadi <- q$par
        }
      }
      else { new.sigmadi <- sigmadi}
      
      cur.sigmadi<- -q$objective
      
      
      q=try(nlminb(start =as.vector(mufij), total_mufij_neg , gradient = grad_mufij_neg,lower = rep(-Inf,num_k), upper = rep(Inf,num_k),
                   Y=Y,x1=x1,x2=x2,v=v,w=w,lamdab=new.lamdab,lamdah=new.lamdah,lamdad=new.lamdad,lamdaf=new.lamdaf,
                   phi=new.phi,pijk=new.pijk,gamma=new.gamma,beta=new.beta,mubi=new.mubi,mudi=new.mudi,
                   muhij=new.muhij,sigmabi=new.sigmabi,sigmahij=new.sigmahij,sigmadi=new.sigmadi,
                   sigmafij=new.sigmafij,
                   num_l=num_l,num_k=num_k,l=l,k=k,control = list(trace = 0, iter.max = 100)), silent = TRUE)
      
      
      if(!inherits(q, "try-error")) {
        if(q$convergence != 0) { if(trace) cat("Optimization algorithm did not converge when it tried to update variational parameters on step ", iter,"\n") }
        if(p.iter > 1 && cur.mufij > -q$objective) { if(trace) cat("Optimization did not improve estimates of variational parameters.\n"); new.mufij <- mufij }
        else {
          if(trace) cat("Variational parameters updated", "\n")
          new.mufij <- q$par
        }
      }
      else { new.mufij <- mufij}
      
      cur.mufij<- -q$objective
      
      
      
      q=try(nlminb(start =as.vector(sigmafij), total_sigmafij_neg , gradient = grad_sigmafij_neg,lower = rep(-Inf,num_k), upper = rep(Inf,num_k),
                   Y=Y,x1=x1,x2=x2,v=v,w=w,lamdab=new.lamdab,lamdah=new.lamdah,lamdad=new.lamdad,lamdaf=new.lamdaf,
                   phi=new.phi,pijk=new.pijk,gamma=new.gamma,beta=new.beta,mubi=new.mubi,muhij=new.muhij,mudi=new.mudi,
                   mufij=new.mufij,sigmabi=new.sigmabi,sigmahij=new.sigmahij,
                   sigmadi=new.sigmadi,
                   num_l=num_l,num_k=num_k,l=l,k=k,control = list(trace = 0, iter.max = 100)), silent = TRUE)
      
      
      if(!inherits(q, "try-error")) {
        if(q$convergence != 0) { if(trace) cat("Optimization algorithm did not converge when it tried to update variational parameters on step ", iter,"\n") }
        if(p.iter > 1 && cur.sigmafij > -q$objective) { if(trace) cat("Optimization did not improve estimates of variational parameters.\n"); new.sigmafij <- sigmafij }
        else {
          if(trace) cat("Variational parameters updated", "\n")
          new.sigmafij <- q$par
        }
      }
      else { new.sigmafij <- sigmafij}
      
      cur.sigmafij<- -q$objective
      
      new.lamdad=sqrt(sum(new.sigmadi^2+new.mudi^2)/num_l)
      
      
      new.lamdaf=sqrt(sum(new.sigmafij^2+new.mufij^2)/(num_k))
      
      
      q=try(nlminb(start =as.vector(gamma), elbo_gamma , gradient = grad_gamma,lower = rep(-Inf,k), upper = rep(Inf,k),
                   Y=Y,x1=x1,x2=x2,v=v,w=w,lamdab=new.lamdab,lamdah=new.lamdah,lamdad=new.lamdad,lamdaf=new.lamdaf,
                   phi=new.phi,pijk=new.pijk,beta=new.beta,mubi=new.mubi,muhij=new.muhij,mudi=new.mudi,
                   mufij=new.mufij,sigmabi=new.sigmabi,sigmahij=new.sigmahij,
                   sigmadi=new.sigmadi,sigmafij=new.sigmafij,num_l=num_l,num_k=num_k,l=l,k=k,control = list(trace = 0, iter.max = 100)), silent = TRUE)
      
      
      
      if(!inherits(q, "try-error")) {
        if(q$convergence != 0) { if(trace) cat("Optimization algorithm did not converge when it tried to update variational parameters on step ", iter,"\n") }
        if(p.iter > 1 && cur.gam > -q$objective) { if(trace) cat("Optimization did not improve estimates of variational parameters.\n"); new.gamma <- gamma }
        else {
          if(trace) cat("Variational parameters updated", "\n")
          new.gamma <- q$par
        }
      }
      else { new.gamma <- gamma}
      cur.gam<- -q$objective
      
      
      
      mudi=new.mudi;mufij=new.mufij;sigmadi=new.sigmadi;sigmafij=new.sigmafij;lamdad=new.lamdad;lamdaf=new.lamdaf
      gamma=new.gamma
      
      
      delta.sigmafij <- abs(new.sigmafij-sigmafij)
      sigmafij <- new.sigmafij
      p.iter <- p.iter+1
    }
    
    B1=(exp(x2%*%new.beta+w%*%new.mubi+w%*%(new.sigmabi^2)/2+v%*%new.muhij+v%*%(new.sigmahij^2)/2))
    eq1=-x1%*%new.gamma-w%*%new.mudi-v%*%new.mufij
    D1=-new.phi*log(1+B1/new.phi)
    
    new.pijk[Y==0]=(1/(1+exp(eq1+D1)))[Y==0]
    
    pijk <- new.pijk
    
    
    objective =total_elbo(Y=Y,x1=x1,x2=x2,v=v,w=w,phi=new.phi,beta=new.beta,gamma=new.gamma,pijk=new.pijk,
                            mubi=new.mubi,muhij=new.muhij,mudi=new.mudi,mufij=new.mufij,lamdab=new.lamdab,lamdah=new.lamdah,
                            sigmabi=new.sigmabi,sigmahij=new.sigmahij,sigmadi=new.sigmadi,sigmafij=new.sigmafij,
                            lamdad=new.lamdad,lamdaf=new.lamdaf,num_l=num_l,num_k=num_k,l=l,k=k)
    new.loglik <- objective
    div=abs(new.loglik-current.loglik)
    err <- abs(new.loglik/current.loglik);
    
    current.loglik <- new.loglik
    
    iter <- iter + 1
    loglik=append(loglik,current.loglik)
  }
  
  
  
  
  out.list=list(Y=Y,x1=x1,x2=x2,w=w,v=v,num_l=num_l,num_k=num_k,num_p=num_p,l=l,ID=ID,cluster=cluster,
                k=k,beta=beta,gamma=gamma,mubi=mubi,muhij=muhij,mudi=mudi,mufij=mufij,phi=phi,pijk=pijk,
                sigmabi=sigmabi,sigmahij=sigmahij,sigmadi=sigmadi,sigmafij=sigmafij,
                lamdab=lamdab,lamdah=lamdah,lamdad=lamdad,lamdaf=lamdaf,loglik_all=loglik,loglik=current.loglik)
    return(out.list)
}
