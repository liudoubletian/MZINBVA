#' @import parallel
#' @import MASS
#' @import numDeriv
grad_beta=function(Y,w,v,x1,x2,lamdab,lamdah,lamdad,lamdaf,phi,beta,gamma,mubi,muhij,mudi,mufij,sigmabi,sigmahij,sigmadi,sigmafij,pijk,num_l,num_k,num_p,l,k,offset){
  new.beta=beta
  new.gamma=gamma
  new.mubi=mubi
  new.muhij=muhij
  new.sigmabi=sigmabi
  new.sigmahij=sigmahij
  new.mudi=mudi
  new.mufij=mufij
  new.sigmadi=sigmadi
  new.sigmafij=sigmafij
  new.lamdab=lamdab
  new.lamdah=lamdah
  new.lamdad=lamdad
  new.lamdaf=lamdaf
  new.pijk=pijk
  new.phi=phi
  
  
  B1=(exp(x2%*%new.beta+w%*%new.mubi+w%*%(new.sigmabi^2)/2+v%*%new.muhij+v%*%(new.sigmahij^2)/2+offset))
  C1=(exp(-x2%*%new.beta-w%*%new.mubi+w%*%(new.sigmabi^2)/2-v%*%new.muhij+v%*%(new.sigmahij^2)/2-offset))
  
  eq2=x2%*%new.beta+w%*%(new.mubi)+v%*%(new.muhij)+offset
  eq1=x1%*%new.gamma+w%*%(new.mudi)+v%*%(new.mufij)
  F1=(exp(x1%*%new.gamma+w%*%new.mudi+w%*%(new.sigmadi^2)/2+v%*%new.mufij+v%*%(new.sigmafij^2)/2))
  
  ini1=(1-new.pijk)*as.vector((Y+new.phi)*new.phi*C1/(1+new.phi*C1)-new.phi)*x2
  ini2=(1-new.pijk)*as.vector(-B1/(1+B1/new.phi))*x2
  if(ncol(x2)==1){
    ini3=sum(ini1[Y>0]+ini2[Y==0])
  }
  else{
    Y1 <- matrix(rep(Y,ncol(x2)),ncol=ncol(x2),byrow = FALSE)
    ini3=colSums(ifelse(Y1>0,ini1,ini2))
  }
  final_beta=ini3
  
  return(-c(final_beta))
}

grad_gamma=function(Y,w,v,x1,x2,lamdab,lamdah,lamdad,lamdaf,phi,beta,gamma,mubi,muhij,mudi,mufij,sigmabi,sigmahij,sigmadi,sigmafij,pijk,num_l,num_k,num_p,l,k,offset){
  new.beta=beta
  new.gamma=gamma
  new.mubi=mubi
  new.muhij=muhij
  new.sigmabi=sigmabi
  new.sigmahij=sigmahij
  new.mudi=mudi
  new.mufij=mufij
  new.sigmadi=sigmadi
  new.sigmafij=sigmafij
  new.lamdab=lamdab
  new.lamdah=lamdah
  new.lamdad=lamdad
  new.lamdaf=lamdaf
  new.pijk=pijk
  new.phi=phi
  
  
  B1=(exp(x2%*%new.beta+w%*%new.mubi+w%*%(new.sigmabi^2)/2+v%*%new.muhij+v%*%(new.sigmahij^2)/2+offset))
  C1=(exp(-x2%*%new.beta-w%*%new.mubi+w%*%(new.sigmabi^2)/2-v%*%new.muhij+v%*%(new.sigmahij^2)/2-offset))
  
  eq2=x2%*%new.beta+w%*%(new.mubi)+v%*%(new.muhij)+offset
  eq1=x1%*%new.gamma+w%*%(new.mudi)+v%*%(new.mufij)
  F1=(exp(x1%*%new.gamma+w%*%new.mudi+w%*%(new.sigmadi^2)/2+v%*%new.mufij+v%*%(new.sigmafij^2)/2))
  
  ini1=(new.pijk)-F1/(1+F1)
  ini1=as.vector(ini1)*x1
  if(ncol(x1)==1){
    ini3=sum(ini1)
  }
  else{
    ini3=colSums(ini1)
  }
  
  
  final_gamma=ini3
  
  return(-c(final_gamma))
}

grad_phi=function(Y,w,v,x1,x2,lamdab,lamdah,lamdad,lamdaf,phi,beta,gamma,mubi,muhij,mudi,mufij,sigmabi,sigmahij,sigmadi,sigmafij,pijk,num_l,num_k,num_p,l,k,offset){
  new.beta=beta
  new.gamma=gamma
  new.mubi=mubi
  new.muhij=muhij
  new.sigmabi=sigmabi
  new.sigmahij=sigmahij
  new.mudi=mudi
  new.mufij=mufij
  new.sigmadi=sigmadi
  new.sigmafij=sigmafij
  new.lamdab=lamdab
  new.lamdah=lamdah
  new.lamdad=lamdad
  new.lamdaf=lamdaf
  new.pijk=pijk
  new.phi=phi
  
  
  B1=(exp(x2%*%new.beta+w%*%new.mubi+w%*%(new.sigmabi^2)/2+v%*%new.muhij+v%*%(new.sigmahij^2)/2+offset))
  C1=(exp(-x2%*%new.beta-w%*%new.mubi+w%*%(new.sigmabi^2)/2-v%*%new.muhij+v%*%(new.sigmahij^2)/2-offset))
  
  eq2=x2%*%new.beta+w%*%(new.mubi)+v%*%(new.muhij)+offset
  eq1=x1%*%new.gamma+w%*%(new.mudi)+v%*%(new.mufij)
  F1=(exp(x1%*%new.gamma+w%*%new.mudi+w%*%(new.sigmadi^2)/2+v%*%new.mufij+v%*%(new.sigmafij^2)/2))
  
  ini1=(1-new.pijk)*(digamma(Y+new.phi)-digamma(new.phi)-log(1+(new.phi*C1))-(Y+new.phi)*C1/(1+new.phi*C1)+
                       log(new.phi)+1-(eq2))
  
  
  ini2=(1-new.pijk)*(-log(1+(B1/new.phi))+B1/(new.phi+B1))
  
  
  final_beta=sum(ini1[Y>0])+sum(ini2[Y==0])
  
  
  return(-c(final_beta))
}

elbo_beta=function(Y,w,v,x1,x2,lamdab,lamdah,lamdad,lamdaf,phi,beta,gamma,mubi,muhij,mudi,mufij,sigmabi,sigmahij,sigmadi,sigmafij,pijk,num_l,num_k,num_p,l,k,offset){
  new.beta=beta
  new.gamma=gamma
  new.mubi=mubi
  new.muhij=muhij
  new.sigmabi=sigmabi
  new.sigmahij=sigmahij
  new.mudi=mudi
  new.mufij=mufij
  new.sigmadi=sigmadi
  new.sigmafij=sigmafij
  new.lamdab=lamdab
  new.lamdah=lamdah
  new.lamdad=lamdad
  new.lamdaf=lamdaf
  new.pijk=pijk
  new.phi=phi
  
  
  B1=(exp(x2%*%new.beta+w%*%new.mubi+w%*%(new.sigmabi^2)/2+v%*%new.muhij+v%*%(new.sigmahij^2)/2+offset))
  C1=(exp(-x2%*%new.beta-w%*%new.mubi+w%*%(new.sigmabi^2)/2-v%*%new.muhij+v%*%(new.sigmahij^2)/2-offset))
  
  eq2=x2%*%new.beta+w%*%(new.mubi)+v%*%(new.muhij)+offset
  eq1=x1%*%new.gamma+w%*%(new.mudi)+v%*%(new.mufij)
  F1=(exp(x1%*%new.gamma+w%*%new.mudi+w%*%(new.sigmadi^2)/2+v%*%new.mufij+v%*%(new.sigmafij^2)/2))
  
  ini2=(1-new.pijk)*(-(Y+new.phi)*log(1+new.phi*C1)-new.phi*eq2)
  
  ini3=(1-new.pijk)*((-new.phi)*log(1+B1/new.phi))
  
  elbo=sum(ini2[Y>0])+sum(ini3[Y==0])
  return(-elbo)
}

elbo_gamma=function(Y,w,v,x1,x2,lamdab,lamdah,lamdad,lamdaf,phi,beta,gamma,mubi,muhij,mudi,mufij,sigmabi,sigmahij,sigmadi,sigmafij,pijk,num_l,num_k,num_p,l,k,offset){
  new.beta=beta
  new.gamma=gamma
  new.mubi=mubi
  new.muhij=muhij
  new.sigmabi=sigmabi
  new.sigmahij=sigmahij
  new.mudi=mudi
  new.mufij=mufij
  new.sigmadi=sigmadi
  new.sigmafij=sigmafij
  new.lamdab=lamdab
  new.lamdah=lamdah
  new.lamdad=lamdad
  new.lamdaf=lamdaf
  new.pijk=pijk
  new.phi=phi
  
  
  B1=(exp(x2%*%new.beta+w%*%new.mubi+w%*%(new.sigmabi^2)/2+v%*%new.muhij+v%*%(new.sigmahij^2)/2+offset))
  C1=(exp(-x2%*%new.beta-w%*%new.mubi+w%*%(new.sigmabi^2)/2-v%*%new.muhij+v%*%(new.sigmahij^2)/2-offset))
  
  eq2=x2%*%new.beta+w%*%(new.mubi)+v%*%(new.muhij)+offset
  eq1=x1%*%new.gamma
  F1=(exp(x1%*%new.gamma+w%*%new.mudi+w%*%(new.sigmadi^2)/2+v%*%new.mufij+v%*%(new.sigmafij^2)/2))
  
  ini1=sum(new.pijk*eq1-log(1+F1))
  
  elbo=ini1
  return(-elbo)
}
elbo_phi=function(Y,w,v,x1,x2,lamdab,lamdah,lamdad,lamdaf,phi,beta,gamma,mubi,muhij,mudi,mufij,sigmabi,sigmahij,sigmadi,sigmafij,pijk,num_l,num_k,num_p,l,k,offset){
  new.beta=beta
  new.gamma=gamma
  new.mubi=mubi
  new.muhij=muhij
  new.sigmabi=sigmabi
  new.sigmahij=sigmahij
  new.mudi=mudi
  new.mufij=mufij
  new.sigmadi=sigmadi
  new.sigmafij=sigmafij
  new.lamdab=lamdab
  new.lamdah=lamdah
  new.lamdad=lamdad
  new.lamdaf=lamdaf
  new.pijk=pijk
  new.phi=phi
  
  
  B1=(exp(x2%*%new.beta+w%*%new.mubi+w%*%(new.sigmabi^2)/2+v%*%new.muhij+v%*%(new.sigmahij^2)/2+offset))
  C1=(exp(-x2%*%new.beta-w%*%new.mubi+w%*%(new.sigmabi^2)/2-v%*%new.muhij+v%*%(new.sigmahij^2)/2-offset))
  
  eq2=x2%*%new.beta+w%*%(new.mubi)+v%*%(new.muhij)+offset
  eq1=x1%*%new.gamma+w%*%(new.mudi)+v%*%(new.mufij)
  F1=(exp(x1%*%new.gamma+w%*%new.mudi+w%*%(new.sigmadi^2)/2+v%*%new.mufij+v%*%(new.sigmafij^2)/2))
  
  
  ini2=(1-new.pijk)*(lfactorial(Y+new.phi-1)-log(gamma(new.phi))-(Y+new.phi)*log(1+new.phi*C1)+
                       new.phi*log(new.phi)-new.phi*eq2)
  
  ini3=(1-new.pijk)*((-new.phi)*log(1+B1/new.phi))
  
  
  elbo=sum(ini2[Y>0])+sum(ini3[Y==0])
  return(-elbo)
}



total_elbo=function(Y,w,v,x1,x2,lamdab,lamdah,lamdad,lamdaf,phi,beta,gamma,mubi,muhij,mudi,mufij,sigmabi,sigmahij,sigmadi,sigmafij,pijk,num_l,num_k,num_p,l,k,offset){
  new.beta=beta
  new.gamma=gamma
  new.mubi=mubi
  new.muhij=muhij
  new.sigmabi=sigmabi
  new.sigmahij=sigmahij
  new.mudi=mudi
  new.mufij=mufij
  new.sigmadi=sigmadi
  new.sigmafij=sigmafij
  new.lamdab=lamdab
  new.lamdah=lamdah
  new.lamdad=lamdad
  new.lamdaf=lamdaf
  new.pijk=pijk
  new.phi=phi
  
  
  B1=(exp(x2%*%new.beta+w%*%new.mubi+w%*%(new.sigmabi^2)/2+v%*%new.muhij+v%*%(new.sigmahij^2)/2+offset))
  C1=(exp(-x2%*%new.beta-w%*%new.mubi+w%*%(new.sigmabi^2)/2-v%*%new.muhij+v%*%(new.sigmahij^2)/2-offset))
  
  eq2=x2%*%new.beta+w%*%(new.mubi)+v%*%(new.muhij)+offset
  eq1=x1%*%new.gamma+w%*%(new.mudi)+v%*%(new.mufij)
  F1=(exp(x1%*%new.gamma+w%*%new.mudi+w%*%(new.sigmadi^2)/2+v%*%new.mufij+v%*%(new.sigmafij^2)/2))
  
  ini1=sum(new.pijk*eq1-log(1+F1))
  ini2=(1-new.pijk)*(lfactorial(Y+new.phi-1)-lfactorial(Y)-log(gamma(new.phi))-(Y+new.phi)*log(1+new.phi*C1)+
                       new.phi*log(new.phi)-new.phi*eq2)
  
  ini3=(1-new.pijk)*((-new.phi)*log(1+B1/new.phi))
  
  elbo_c=sum(log(new.sigmabi/new.lamdab))-0.5*sum((new.mubi/new.lamdab)^2+(new.sigmabi/new.lamdab)^2)+0.5*num_l+
    sum(log(new.sigmahij/new.lamdah))-0.5*sum((new.muhij/new.lamdah)^2+(new.sigmahij/new.lamdah)^2)+0.5*num_k+
    sum(log(new.sigmadi/new.lamdad))-0.5*sum((new.mudi/new.lamdad)^2+(new.sigmadi/new.lamdad)^2)+0.5*num_l+
    sum(log(new.sigmafij/new.lamdaf))-0.5*sum((new.mufij/new.lamdaf)^2+(new.sigmafij/new.lamdaf)^2)+0.5*num_k-
    sum(new.pijk[new.pijk!=0]*log(new.pijk[new.pijk!=0]))-sum((1-new.pijk[new.pijk!=1])*log(1-new.pijk[new.pijk!=1]))
  
  elbo=ini1+sum(ini2[Y>0])+sum(ini3[Y==0])+elbo_c
  return(elbo)
  
}
total_elbo_neg=function(Y,w,v,x1,x2,lamdab,lamdah,lamdad,lamdaf,phi,beta,gamma,mubi,muhij,mudi,mufij,sigmabi,sigmahij,sigmadi,sigmafij,pijk,num_l,num_k,num_p,l,k,offset){
  new.beta=beta
  new.gamma=gamma
  new.mubi=mubi
  new.muhij=muhij
  new.sigmabi=sigmabi
  new.sigmahij=sigmahij
  new.mudi=mudi
  new.mufij=mufij
  new.sigmadi=sigmadi
  new.sigmafij=sigmafij
  new.lamdab=lamdab
  new.lamdah=lamdah
  new.lamdad=lamdad
  new.lamdaf=lamdaf
  new.pijk=pijk
  new.phi=phi
  
  
  B1=(exp(x2%*%new.beta+w%*%new.mubi+w%*%(new.sigmabi^2)/2+v%*%new.muhij+v%*%(new.sigmahij^2)/2+offset))
  C1=(exp(-x2%*%new.beta-w%*%new.mubi+w%*%(new.sigmabi^2)/2-v%*%new.muhij+v%*%(new.sigmahij^2)/2-offset))
  
  eq2=x2%*%new.beta+w%*%(new.mubi)+v%*%(new.muhij)+offset
  eq1=x1%*%new.gamma+w%*%(new.mudi)+v%*%(new.mufij)
  F1=(exp(x1%*%new.gamma+w%*%new.mudi+w%*%(new.sigmadi^2)/2+v%*%new.mufij+v%*%(new.sigmafij^2)/2))
  
  ini1=sum(new.pijk*eq1-log(1+F1))
  ini2=(1-new.pijk)*(lfactorial(Y+new.phi-1)-lfactorial(Y)-log(gamma(new.phi))-(Y+new.phi)*log(1+new.phi*C1)+
                       new.phi*log(new.phi)-new.phi*eq2)
  
  ini3=(1-new.pijk)*((-new.phi)*log(1+B1/new.phi))
  
  elbo_c=sum(log(new.sigmabi/new.lamdab))-0.5*sum((new.mubi/new.lamdab)^2+(new.sigmabi/new.lamdab)^2)+0.5*num_l+
    sum(log(new.sigmahij/new.lamdah))-0.5*sum((new.muhij/new.lamdah)^2+(new.sigmahij/new.lamdah)^2)+0.5*num_k+
    sum(log(new.sigmadi/new.lamdad))-0.5*sum((new.mudi/new.lamdad)^2+(new.sigmadi/new.lamdad)^2)+0.5*num_l+
    sum(log(new.sigmafij/new.lamdaf))-0.5*sum((new.mufij/new.lamdaf)^2+(new.sigmafij/new.lamdaf)^2)+0.5*num_k-
    sum(new.pijk[new.pijk!=0]*log(new.pijk[new.pijk!=0]))-sum((1-new.pijk[new.pijk!=1])*log(1-new.pijk[new.pijk!=1]))
  
  elbo=ini1+sum(ini2[Y>0])+sum(ini3[Y==0])+elbo_c
  return(-elbo)
  
}





total_mubi_neg=function(Y,w,v,x1,x2,lamdab,lamdah,lamdad,lamdaf,phi,beta,gamma,mubi,muhij,mudi,mufij,sigmabi,sigmahij,sigmadi,sigmafij,pijk,num_l,num_k,num_p,l,k,offset){
  new.beta=beta
  new.gamma=gamma
  new.mubi=mubi
  new.muhij=muhij
  new.sigmabi=sigmabi
  new.sigmahij=sigmahij
  new.mudi=mudi
  new.mufij=mufij
  new.sigmadi=sigmadi
  new.sigmafij=sigmafij
  new.lamdab=lamdab
  new.lamdah=lamdah
  new.lamdad=lamdad
  new.lamdaf=lamdaf
  new.pijk=pijk
  new.phi=phi
  
  
  B1=(exp(x2%*%new.beta+w%*%new.mubi+w%*%(new.sigmabi^2)/2+v%*%new.muhij+v%*%(new.sigmahij^2)/2+offset))
  C1=(exp(-x2%*%new.beta-w%*%new.mubi+w%*%(new.sigmabi^2)/2-v%*%new.muhij+v%*%(new.sigmahij^2)/2-offset))
  
  eq2=x2%*%new.beta+w%*%(new.mubi)+v%*%(new.muhij)+offset
  eq1=x1%*%new.gamma+w%*%(new.mudi)+v%*%(new.mufij)
  F1=(exp(x1%*%new.gamma+w%*%new.mudi+w%*%(new.sigmadi^2)/2+v%*%new.mufij+v%*%(new.sigmafij^2)/2))
  
  ini2=(1-new.pijk)*(-(Y+new.phi)*log(1+new.phi*C1)-new.phi*eq2)
  
  
  elbo_c=-0.5*sum((new.mubi/new.lamdab)^2)
  elbo=sum(ini2)+elbo_c
  return(-elbo)
  
}

total_muhij_neg=function(Y,w,v,x1,x2,lamdab,lamdah,lamdad,lamdaf,phi,beta,gamma,mubi,muhij,mudi,mufij,sigmabi,sigmahij,sigmadi,sigmafij,pijk,num_l,num_k,num_p,l,k,offset){
  new.beta=beta
  new.gamma=gamma
  new.mubi=mubi
  new.muhij=muhij
  new.sigmabi=sigmabi
  new.sigmahij=sigmahij
  new.mudi=mudi
  new.mufij=mufij
  new.sigmadi=sigmadi
  new.sigmafij=sigmafij
  new.lamdab=lamdab
  new.lamdah=lamdah
  new.lamdad=lamdad
  new.lamdaf=lamdaf
  new.pijk=pijk
  new.phi=phi
  
  
  B1=(exp(x2%*%new.beta+w%*%new.mubi+w%*%(new.sigmabi^2)/2+v%*%new.muhij+v%*%(new.sigmahij^2)/2+offset))
  C1=(exp(-x2%*%new.beta-w%*%new.mubi+w%*%(new.sigmabi^2)/2-v%*%new.muhij+v%*%(new.sigmahij^2)/2-offset))
  
  eq2=x2%*%new.beta+w%*%(new.mubi)+v%*%(new.muhij)+offset
  eq1=x1%*%new.gamma+w%*%(new.mudi)+v%*%(new.mufij)
  F1=(exp(x1%*%new.gamma+w%*%new.mudi+w%*%(new.sigmadi^2)/2+v%*%new.mufij+v%*%(new.sigmafij^2)/2))
  
  ini2=(1-new.pijk)*(-(Y+new.phi)*log(1+new.phi*C1)-new.phi*eq2)
  
  
  elbo_c=-0.5*sum((new.muhij/new.lamdah)^2)
  
  elbo=sum(ini2)+elbo_c
  
  return(-elbo)
  
}

total_sigmabi_neg=function(Y,w,v,x1,x2,lamdab,lamdah,lamdad,lamdaf,phi,beta,gamma,mubi,muhij,mudi,mufij,sigmabi,sigmahij,sigmadi,sigmafij,pijk,num_l,num_k,num_p,l,k,offset){
  new.beta=beta
  new.gamma=gamma
  new.mubi=mubi
  new.muhij=muhij
  new.sigmabi=sigmabi
  new.sigmahij=sigmahij
  new.mudi=mudi
  new.mufij=mufij
  new.sigmadi=sigmadi
  new.sigmafij=sigmafij
  new.lamdab=lamdab
  new.lamdah=lamdah
  new.lamdad=lamdad
  new.lamdaf=lamdaf
  new.pijk=pijk
  new.phi=phi
  
  
  B1=(exp(x2%*%new.beta+w%*%new.mubi+w%*%(new.sigmabi^2)/2+v%*%new.muhij+v%*%(new.sigmahij^2)/2+offset))
  C1=(exp(-x2%*%new.beta-w%*%new.mubi+w%*%(new.sigmabi^2)/2-v%*%new.muhij+v%*%(new.sigmahij^2)/2-offset))
  
  eq2=x2%*%new.beta+w%*%(new.mubi)+v%*%(new.muhij)+offset
  eq1=x1%*%new.gamma+w%*%(new.mudi)+v%*%(new.mufij)
  F1=(exp(x1%*%new.gamma+w%*%new.mudi+w%*%(new.sigmadi^2)/2+v%*%new.mufij+v%*%(new.sigmafij^2)/2))
  
  ini2=(1-new.pijk)*(-(Y+new.phi)*log(1+new.phi*C1))
  
  
  elbo_c=sum(log(new.sigmabi))-0.5*sum((new.sigmabi/new.lamdab)^2)
  elbo=sum(ini2)+elbo_c
  return(-elbo)
  
}

total_sigmahij_neg=function(Y,w,v,x1,x2,lamdab,lamdah,lamdad,lamdaf,phi,beta,gamma,mubi,muhij,mudi,mufij,sigmabi,sigmahij,sigmadi,sigmafij,pijk,num_l,num_k,num_p,l,k,offset){
  new.beta=beta
  new.gamma=gamma
  new.mubi=mubi
  new.muhij=muhij
  new.sigmabi=sigmabi
  new.sigmahij=sigmahij
  new.mudi=mudi
  new.mufij=mufij
  new.sigmadi=sigmadi
  new.sigmafij=sigmafij
  new.lamdab=lamdab
  new.lamdah=lamdah
  new.lamdad=lamdad
  new.lamdaf=lamdaf
  new.pijk=pijk
  new.phi=phi
  
  
  B1=(exp(x2%*%new.beta+w%*%new.mubi+w%*%(new.sigmabi^2)/2+v%*%new.muhij+v%*%(new.sigmahij^2)/2+offset))
  C1=(exp(-x2%*%new.beta-w%*%new.mubi+w%*%(new.sigmabi^2)/2-v%*%new.muhij+v%*%(new.sigmahij^2)/2-offset))
  
  eq2=x2%*%new.beta+w%*%(new.mubi)+v%*%(new.muhij)+offset
  eq1=x1%*%new.gamma+w%*%(new.mudi)+v%*%(new.mufij)
  F1=(exp(x1%*%new.gamma+w%*%new.mudi+w%*%(new.sigmadi^2)/2+v%*%new.mufij+v%*%(new.sigmafij^2)/2))
  
  ini2=(1-new.pijk)*(-(Y+new.phi)*log(1+new.phi*C1))
  
  
  elbo_c=sum(log(new.sigmahij))-0.5*sum((new.sigmahij/new.lamdah)^2)
  
  elbo=sum(ini2)+elbo_c
  return(-elbo)
  
}

grad_mubi_neg=function(Y,w,v,x1,x2,lamdab,lamdah,lamdad,lamdaf,phi,beta,gamma,mubi,muhij,mudi,mufij,sigmabi,sigmahij,sigmadi,sigmafij,pijk,num_l,num_k,num_p,l,k,offset){
  new.beta=beta
  new.gamma=gamma
  new.mubi=mubi
  new.muhij=muhij
  new.sigmabi=sigmabi
  new.sigmahij=sigmahij
  new.mudi=mudi
  new.mufij=mufij
  new.sigmadi=sigmadi
  new.sigmafij=sigmafij
  new.lamdab=lamdab
  new.lamdah=lamdah
  new.lamdad=lamdad
  new.lamdaf=lamdaf
  new.pijk=pijk
  new.phi=phi
  
  
  B1=(exp(x2%*%new.beta+w%*%new.mubi+w%*%(new.sigmabi^2)/2+v%*%new.muhij+v%*%(new.sigmahij^2)/2+offset))
  C1=(exp(-x2%*%new.beta-w%*%new.mubi+w%*%(new.sigmabi^2)/2-v%*%new.muhij+v%*%(new.sigmahij^2)/2-offset))
  
  eq2=x2%*%new.beta+w%*%(new.mubi)+v%*%(new.muhij)+offset
  eq1=x1%*%new.gamma+w%*%(new.mudi)+v%*%(new.mufij)
  F1=(exp(x1%*%new.gamma+w%*%new.mudi+w%*%(new.sigmadi^2)/2+v%*%new.mufij+v%*%(new.sigmafij^2)/2))
  
  
  ini1=as.vector((1-new.pijk)*(((Y+new.phi)*new.phi*C1/(1+new.phi*C1))-new.phi))*w
  ini2=as.vector((1-new.pijk)*as.vector(-B1/(1+B1/new.phi)))*w
  
  Y1 <- matrix(rep(Y,ncol(ini1)),ncol=ncol(ini1),byrow = FALSE)
  ini3=colSums(ifelse(Y1>0,ini1,ini2))-new.mubi/(new.lamdab^2)
  
  
  
  
  
  return(c(-ini3))
}

grad_muhij_neg=function(Y,w,v,x1,x2,lamdab,lamdah,lamdad,lamdaf,phi,beta,gamma,mubi,muhij,mudi,mufij,sigmabi,sigmahij,sigmadi,sigmafij,pijk,num_l,num_k,num_p,l,k,offset){
  new.beta=beta
  new.gamma=gamma
  new.mubi=mubi
  new.muhij=muhij
  new.sigmabi=sigmabi
  new.sigmahij=sigmahij
  new.mudi=mudi
  new.mufij=mufij
  new.sigmadi=sigmadi
  new.sigmafij=sigmafij
  new.lamdab=lamdab
  new.lamdah=lamdah
  new.lamdad=lamdad
  new.lamdaf=lamdaf
  new.pijk=pijk
  new.phi=phi
  
  
  B1=(exp(x2%*%new.beta+w%*%new.mubi+w%*%(new.sigmabi^2)/2+v%*%new.muhij+v%*%(new.sigmahij^2)/2+offset))
  C1=(exp(-x2%*%new.beta-w%*%new.mubi+w%*%(new.sigmabi^2)/2-v%*%new.muhij+v%*%(new.sigmahij^2)/2-offset))
  
  eq2=x2%*%new.beta+w%*%(new.mubi)+v%*%(new.muhij)+offset
  eq1=x1%*%new.gamma+w%*%(new.mudi)+v%*%(new.mufij)
  F1=(exp(x1%*%new.gamma+w%*%new.mudi+w%*%(new.sigmadi^2)/2+v%*%new.mufij+v%*%(new.sigmafij^2)/2))
  
  
  ini1=as.vector((1-new.pijk)*(((Y+new.phi)*new.phi*C1/(1+new.phi*C1))-new.phi))*v
  ini2=as.vector((1-new.pijk)*as.vector(-B1/(1+B1/new.phi)))*v
  
  Y1 <- matrix(rep(Y,ncol(ini1)),ncol=ncol(ini1),byrow = FALSE)
  ini3=colSums(ifelse(Y1>0,ini1,ini2))-new.muhij/(new.lamdah^2)
  
  
  
  
  
  return(c(-ini3))
}


grad_sigmabi_neg=function(Y,w,v,x1,x2,lamdab,lamdah,lamdad,lamdaf,phi,beta,gamma,mubi,muhij,mudi,mufij,sigmabi,sigmahij,sigmadi,sigmafij,pijk,num_l,num_k,num_p,l,k,offset){
  new.beta=beta
  new.gamma=gamma
  new.mubi=mubi
  new.muhij=muhij
  new.sigmabi=sigmabi
  new.sigmahij=sigmahij
  new.mudi=mudi
  new.mufij=mufij
  new.sigmadi=sigmadi
  new.sigmafij=sigmafij
  new.lamdab=lamdab
  new.lamdah=lamdah
  new.lamdad=lamdad
  new.lamdaf=lamdaf
  new.pijk=pijk
  new.phi=phi
  
  
  B1=(exp(x2%*%new.beta+w%*%new.mubi+w%*%(new.sigmabi^2)/2+v%*%new.muhij+v%*%(new.sigmahij^2)/2+offset))
  C1=(exp(-x2%*%new.beta-w%*%new.mubi+w%*%(new.sigmabi^2)/2-v%*%new.muhij+v%*%(new.sigmahij^2)/2-offset))
  
  eq2=x2%*%new.beta+w%*%(new.mubi)+v%*%(new.muhij)+offset
  eq1=x1%*%new.gamma+w%*%(new.mudi)+v%*%(new.mufij)
  F1=(exp(x1%*%new.gamma+w%*%new.mudi+w%*%(new.sigmadi^2)/2+v%*%new.mufij+v%*%(new.sigmafij^2)/2))
  
  ini1=as.vector((1-new.pijk)*((-(Y+new.phi)*new.phi*C1*w%*%new.sigmabi/(1+new.phi*C1))))*w
  ini2=as.vector((1-new.pijk)*as.vector(-B1*w%*%new.sigmabi/(1+B1/new.phi)))*w
  
  Y1 <- matrix(rep(Y,ncol(ini1)),ncol=ncol(ini1),byrow = FALSE)
  ini3=colSums(ifelse(Y1>0,ini1,ini2))+1/new.sigmabi-new.sigmabi/(new.lamdab^2)
  
  
  
  
  return(c(-ini3))
  
}

grad_sigmahij_neg=function(Y,w,v,x1,x2,lamdab,lamdah,lamdad,lamdaf,phi,beta,gamma,mubi,muhij,mudi,mufij,sigmabi,sigmahij,sigmadi,sigmafij,pijk,num_l,num_k,num_p,l,k,offset){
  new.beta=beta
  new.gamma=gamma
  new.mubi=mubi
  new.muhij=muhij
  new.sigmabi=sigmabi
  new.sigmahij=sigmahij
  new.mudi=mudi
  new.mufij=mufij
  new.sigmadi=sigmadi
  new.sigmafij=sigmafij
  new.lamdab=lamdab
  new.lamdah=lamdah
  new.lamdad=lamdad
  new.lamdaf=lamdaf
  new.pijk=pijk
  new.phi=phi
  
  
  B1=(exp(x2%*%new.beta+w%*%new.mubi+w%*%(new.sigmabi^2)/2+v%*%new.muhij+v%*%(new.sigmahij^2)/2+offset))
  C1=(exp(-x2%*%new.beta-w%*%new.mubi+w%*%(new.sigmabi^2)/2-v%*%new.muhij+v%*%(new.sigmahij^2)/2-offset))
  
  eq2=x2%*%new.beta+w%*%(new.mubi)+v%*%(new.muhij)+offset
  eq1=x1%*%new.gamma+w%*%(new.mudi)+v%*%(new.mufij)
  F1=(exp(x1%*%new.gamma+w%*%new.mudi+w%*%(new.sigmadi^2)/2+v%*%new.mufij+v%*%(new.sigmafij^2)/2))
  
  
  ini1=as.vector((1-new.pijk)*((-(Y+new.phi)*new.phi*C1*v%*%new.sigmahij/(1+new.phi*C1))))*v
  ini2=as.vector((1-new.pijk)*as.vector(-B1*v%*%new.sigmahij/(1+B1/new.phi)))*v
  
  Y1 <- matrix(rep(Y,ncol(ini1)),ncol=ncol(ini1),byrow = FALSE)
  ini3=colSums(ifelse(Y1>0,ini1,ini2))+1/new.sigmahij-new.sigmahij/(new.lamdah^2)
  
  
  
  
  
  
  return(c(-ini3))
}

total_mudi_neg=function(Y,w,v,x1,x2,lamdab,lamdah,lamdad,lamdaf,phi,beta,gamma,mubi,muhij,mudi,mufij,sigmabi,sigmahij,sigmadi,sigmafij,pijk,num_l,num_k,num_p,l,k,offset){
  new.beta=beta
  new.gamma=gamma
  new.mubi=mubi
  new.muhij=muhij
  new.sigmabi=sigmabi
  new.sigmahij=sigmahij
  new.mudi=mudi
  new.mufij=mufij
  new.sigmadi=sigmadi
  new.sigmafij=sigmafij
  new.lamdab=lamdab
  new.lamdah=lamdah
  new.lamdad=lamdad
  new.lamdaf=lamdaf
  new.pijk=pijk
  new.phi=phi
  
  
  B1=(exp(x2%*%new.beta+w%*%new.mubi+w%*%(new.sigmabi^2)/2+v%*%new.muhij+v%*%(new.sigmahij^2)/2+offset))
  C1=(exp(-x2%*%new.beta-w%*%new.mubi+w%*%(new.sigmabi^2)/2-v%*%new.muhij+v%*%(new.sigmahij^2)/2-offset))
  
  eq2=x2%*%new.beta+w%*%(new.mubi)+v%*%(new.muhij)+offset
  eq1=w%*%(new.mudi)
  F1=(exp(x1%*%new.gamma+w%*%new.mudi+w%*%(new.sigmadi^2)/2+v%*%new.mufij+v%*%(new.sigmafij^2)/2))
  
  ini1=sum(new.pijk*eq1-log(1+F1))
  
  elbo_c=-0.5*sum((new.mudi/new.lamdad)^2)
  
  elbo=ini1+elbo_c
  return(-elbo)
  
}

total_mufij_neg=function(Y,w,v,x1,x2,lamdab,lamdah,lamdad,lamdaf,phi,beta,gamma,mubi,muhij,mudi,mufij,sigmabi,sigmahij,sigmadi,sigmafij,pijk,num_l,num_k,num_p,l,k,offset){
  new.beta=beta
  new.gamma=gamma
  new.mubi=mubi
  new.muhij=muhij
  new.sigmabi=sigmabi
  new.sigmahij=sigmahij
  new.mudi=mudi
  new.mufij=mufij
  new.sigmadi=sigmadi
  new.sigmafij=sigmafij
  new.lamdab=lamdab
  new.lamdah=lamdah
  new.lamdad=lamdad
  new.lamdaf=lamdaf
  new.pijk=pijk
  new.phi=phi
  
  
  B1=(exp(x2%*%new.beta+w%*%new.mubi+w%*%(new.sigmabi^2)/2+v%*%new.muhij+v%*%(new.sigmahij^2)/2+offset))
  C1=(exp(-x2%*%new.beta-w%*%new.mubi+w%*%(new.sigmabi^2)/2-v%*%new.muhij+v%*%(new.sigmahij^2)/2-offset))
  
  eq2=x2%*%new.beta+w%*%(new.mubi)+v%*%(new.muhij)+offset
  eq1=v%*%(new.mufij)
  F1=(exp(x1%*%new.gamma+w%*%new.mudi+w%*%(new.sigmadi^2)/2+v%*%new.mufij+v%*%(new.sigmafij^2)/2))
  
  ini1=sum(new.pijk*eq1-log(1+F1))
  
  elbo_c=-0.5*sum((new.mufij/new.lamdaf)^2)
  
  elbo=ini1+elbo_c
  return(-elbo)
  
}

total_sigmadi_neg=function(Y,w,v,x1,x2,lamdab,lamdah,lamdad,lamdaf,phi,beta,gamma,mubi,muhij,mudi,mufij,sigmabi,sigmahij,sigmadi,sigmafij,pijk,num_l,num_k,num_p,l,k,offset){
  new.beta=beta
  new.gamma=gamma
  new.mubi=mubi
  new.muhij=muhij
  new.sigmabi=sigmabi
  new.sigmahij=sigmahij
  new.mudi=mudi
  new.mufij=mufij
  new.sigmadi=sigmadi
  new.sigmafij=sigmafij
  new.lamdab=lamdab
  new.lamdah=lamdah
  new.lamdad=lamdad
  new.lamdaf=lamdaf
  new.pijk=pijk
  new.phi=phi
  
  
  B1=(exp(x2%*%new.beta+w%*%new.mubi+w%*%(new.sigmabi^2)/2+v%*%new.muhij+v%*%(new.sigmahij^2)/2+offset))
  C1=(exp(-x2%*%new.beta-w%*%new.mubi+w%*%(new.sigmabi^2)/2-v%*%new.muhij+v%*%(new.sigmahij^2)/2-offset))
  
  eq2=x2%*%new.beta+w%*%(new.mubi)+v%*%(new.muhij)+offset
  eq1=x1%*%new.gamma+w%*%(new.mudi)+v%*%(new.mufij)
  F1=(exp(x1%*%new.gamma+w%*%new.mudi+w%*%(new.sigmadi^2)/2+v%*%new.mufij+v%*%(new.sigmafij^2)/2))
  
  ini1=sum(-log(1+F1))
  
  elbo_c=sum(log(new.sigmadi))-0.5*sum((new.sigmadi/new.lamdad)^2)
  
  elbo=ini1+elbo_c
  return(-elbo)
  
}

total_sigmafij_neg=function(Y,w,v,x1,x2,lamdab,lamdah,lamdad,lamdaf,phi,beta,gamma,mubi,muhij,mudi,mufij,sigmabi,sigmahij,sigmadi,sigmafij,pijk,num_l,num_k,num_p,l,k,offset){
  new.beta=beta
  new.gamma=gamma
  new.mubi=mubi
  new.muhij=muhij
  new.sigmabi=sigmabi
  new.sigmahij=sigmahij
  new.mudi=mudi
  new.mufij=mufij
  new.sigmadi=sigmadi
  new.sigmafij=sigmafij
  new.lamdab=lamdab
  new.lamdah=lamdah
  new.lamdad=lamdad
  new.lamdaf=lamdaf
  new.pijk=pijk
  new.phi=phi
  
  
  B1=(exp(x2%*%new.beta+w%*%new.mubi+w%*%(new.sigmabi^2)/2+v%*%new.muhij+v%*%(new.sigmahij^2)/2+offset))
  C1=(exp(-x2%*%new.beta-w%*%new.mubi+w%*%(new.sigmabi^2)/2-v%*%new.muhij+v%*%(new.sigmahij^2)/2-offset))
  
  eq2=x2%*%new.beta+w%*%(new.mubi)+v%*%(new.muhij)+offset
  eq1=x1%*%new.gamma+w%*%(new.mudi)+v%*%(new.mufij)
  F1=(exp(x1%*%new.gamma+w%*%new.mudi+w%*%(new.sigmadi^2)/2+v%*%new.mufij+v%*%(new.sigmafij^2)/2))
  
  ini1=sum(-log(1+F1))
  
  elbo_c=sum(log(new.sigmafij))-0.5*sum((new.sigmafij/new.lamdaf)^2)
  
  elbo=ini1+elbo_c
  return(-elbo)
  
}

grad_mudi_neg=function(Y,w,v,x1,x2,lamdab,lamdah,lamdad,lamdaf,phi,beta,gamma,mubi,muhij,mudi,mufij,sigmabi,sigmahij,sigmadi,sigmafij,pijk,num_l,num_k,num_p,l,k,offset){
  new.beta=beta
  new.gamma=gamma
  new.mubi=mubi
  new.muhij=muhij
  new.sigmabi=sigmabi
  new.sigmahij=sigmahij
  new.mudi=mudi
  new.mufij=mufij
  new.sigmadi=sigmadi
  new.sigmafij=sigmafij
  new.lamdab=lamdab
  new.lamdah=lamdah
  new.lamdad=lamdad
  new.lamdaf=lamdaf
  new.pijk=pijk
  new.phi=phi
  
  
  B1=(exp(x2%*%new.beta+w%*%new.mubi+w%*%(new.sigmabi^2)/2+v%*%new.muhij+v%*%(new.sigmahij^2)/2+offset))
  C1=(exp(-x2%*%new.beta-w%*%new.mubi+w%*%(new.sigmabi^2)/2-v%*%new.muhij+v%*%(new.sigmahij^2)/2-offset))
  
  eq2=x2%*%new.beta+w%*%(new.mubi)+v%*%(new.muhij)+offset
  eq1=x1%*%new.gamma+w%*%(new.mudi)+v%*%(new.mufij)
  F1=(exp(x1%*%new.gamma+w%*%new.mudi+w%*%(new.sigmadi^2)/2+v%*%new.mufij+v%*%(new.sigmafij^2)/2))
  
  
  ini1=as.vector(new.pijk-F1/(1+F1))*w
  ini3=colSums(ini1)-new.mudi/(new.lamdad^2)
  
  
  
  
  return(c(-ini3))
}

grad_mufij_neg=function(Y,w,v,x1,x2,lamdab,lamdah,lamdad,lamdaf,phi,beta,gamma,mubi,muhij,mudi,mufij,sigmabi,sigmahij,sigmadi,sigmafij,pijk,num_l,num_k,num_p,l,k,offset){
  new.beta=beta
  new.gamma=gamma
  new.mubi=mubi
  new.muhij=muhij
  new.sigmabi=sigmabi
  new.sigmahij=sigmahij
  new.mudi=mudi
  new.mufij=mufij
  new.sigmadi=sigmadi
  new.sigmafij=sigmafij
  new.lamdab=lamdab
  new.lamdah=lamdah
  new.lamdad=lamdad
  new.lamdaf=lamdaf
  new.pijk=pijk
  new.phi=phi
  
  
  B1=(exp(x2%*%new.beta+w%*%new.mubi+w%*%(new.sigmabi^2)/2+v%*%new.muhij+v%*%(new.sigmahij^2)/2+offset))
  C1=(exp(-x2%*%new.beta-w%*%new.mubi+w%*%(new.sigmabi^2)/2-v%*%new.muhij+v%*%(new.sigmahij^2)/2-offset))
  
  eq2=x2%*%new.beta+w%*%(new.mubi)+v%*%(new.muhij)+offset
  eq1=x1%*%new.gamma+w%*%(new.mudi)+v%*%(new.mufij)
  F1=(exp(x1%*%new.gamma+w%*%new.mudi+w%*%(new.sigmadi^2)/2+v%*%new.mufij+v%*%(new.sigmafij^2)/2))
  
  
  ini1=as.vector(new.pijk-F1/(1+F1))*v
  ini3=colSums(ini1)-new.mufij/(new.lamdaf^2)
  
  
  
  
  
  
  return(c(-ini3))
}


grad_sigmadi_neg=function(Y,w,v,x1,x2,lamdab,lamdah,lamdad,lamdaf,phi,beta,gamma,mubi,muhij,mudi,mufij,sigmabi,sigmahij,sigmadi,sigmafij,pijk,num_l,num_k,num_p,l,k,offset){
  new.beta=beta
  new.gamma=gamma
  new.mubi=mubi
  new.muhij=muhij
  new.sigmabi=sigmabi
  new.sigmahij=sigmahij
  new.mudi=mudi
  new.mufij=mufij
  new.sigmadi=sigmadi
  new.sigmafij=sigmafij
  new.lamdab=lamdab
  new.lamdah=lamdah
  new.lamdad=lamdad
  new.lamdaf=lamdaf
  new.pijk=pijk
  new.phi=phi
  
  
  B1=(exp(x2%*%new.beta+w%*%new.mubi+w%*%(new.sigmabi^2)/2+v%*%new.muhij+v%*%(new.sigmahij^2)/2+offset))
  C1=(exp(-x2%*%new.beta-w%*%new.mubi+w%*%(new.sigmabi^2)/2-v%*%new.muhij+v%*%(new.sigmahij^2)/2-offset))
  
  eq2=x2%*%new.beta+w%*%(new.mubi)+v%*%(new.muhij)+offset
  eq1=x1%*%new.gamma+w%*%(new.mudi)+v%*%(new.mufij)
  F1=(exp(x1%*%new.gamma+w%*%new.mudi+w%*%(new.sigmadi^2)/2+v%*%new.mufij+v%*%(new.sigmafij^2)/2))
  
  ini1=as.vector(-F1*w%*%new.sigmadi/(1+F1))*w
  
  ini3=(colSums(ini1)+1/new.sigmadi-new.sigmadi/(new.lamdad^2))
  
  
  
  
  
  return(c(-ini3))
  
}

grad_sigmafij_neg=function(Y,w,v,x1,x2,lamdab,lamdah,lamdad,lamdaf,phi,beta,gamma,mubi,muhij,mudi,mufij,sigmabi,sigmahij,sigmadi,sigmafij,pijk,num_l,num_k,num_p,l,k,offset){
  new.beta=beta
  new.gamma=gamma
  new.mubi=mubi
  new.muhij=muhij
  new.sigmabi=sigmabi
  new.sigmahij=sigmahij
  new.mudi=mudi
  new.mufij=mufij
  new.sigmadi=sigmadi
  new.sigmafij=sigmafij
  new.lamdab=lamdab
  new.lamdah=lamdah
  new.lamdad=lamdad
  new.lamdaf=lamdaf
  new.pijk=pijk
  new.phi=phi
  
  
  B1=(exp(x2%*%new.beta+w%*%new.mubi+w%*%(new.sigmabi^2)/2+v%*%new.muhij+v%*%(new.sigmahij^2)/2+offset))
  C1=(exp(-x2%*%new.beta-w%*%new.mubi+w%*%(new.sigmabi^2)/2-v%*%new.muhij+v%*%(new.sigmahij^2)/2-offset))
  
  eq2=x2%*%new.beta+w%*%(new.mubi)+v%*%(new.muhij)+offset
  eq1=x1%*%new.gamma+w%*%(new.mudi)+v%*%(new.mufij)
  F1=(exp(x1%*%new.gamma+w%*%new.mudi+w%*%(new.sigmadi^2)/2+v%*%new.mufij+v%*%(new.sigmafij^2)/2))
  
  
  ini1=as.vector(-F1*v%*%new.sigmafij/(1+F1))*v
  
  ini3=(colSums(ini1)+1/new.sigmafij-new.sigmafij/(new.lamdaf^2))
  
  
  
  
  
  
  return(c(-ini3))
}

