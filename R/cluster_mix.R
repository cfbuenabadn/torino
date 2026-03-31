#'@title Xing's smoothed topic model
#'@export
#'
cluster.mix=function(y,smooth=TRUE,pi0=NULL,phi0=NULL,K,tol,maxit,nugget=FALSE,nug_control=list()){
  t0 = Sys.time()
  n=dim(y)[1]
  B=dim(y)[2]
  #if(is.null(pseudocounts)) pseudocounts=10^(round(log10(1/B/100000)))

  if(is.null(pi0)|is.null(phi0)){
    kmeans.init=kmeans(y,K,nstart=5)
  }

  if(is.null(pi0)) pi0=rep(1,n)%o%normalize(as.vector(table(kmeans.init$cluster)))

  if(is.null(phi0)){
    phi0=kmeans.init$centers
    #phi0[phi0==0]=pseudocounts
    #phi0=t(apply(phi0,1,normalize))
    row.names(phi0)=NULL
  }

  out=EMproc.mix(y,smooth,pi0,phi0,n,K,B,tol,maxit,nugget,nug_control)
  t1 = Sys.time()
  return(list(pi=out$pi,phi=out$phi,lambda=out$lambda,gamma=out$gamma,loglik=out$loglik,run_time = difftime(t1,t0)))
}

# normalize=function(x){
#   #if(sum(abs(x))!=0){
#   return(x/sum(x))
#   #}else{
#   #  return(rep(0,length(x)))
#   #}
# }

smooth.lambda = function(lambda,nugget,nug_control){
  if(nugget){
    control0 = nug_control_default()
    control1 = modifyList(control0,nug_control,keep.null = TRUE)
    return(t(apply(lambda,1,function(z){
      fit = smash.gen.poiss(z,filter.number=control1$filter.number,
      family=control1$family,nugget=control1$nugget,
      robust=control1$robust,
      robust.q = control1$robust.q,
      transformation = control1$transformation,
      method = control1$method,
      nug.init = control1$nug.init,
      ash.pm = control1$ash.pm,
      eps = control1$eps,
      maxiter = control1$maxiter,
      tol = control1$tol)
      return(fit$lambda.est)
    })))
  }else{
    return(t(apply(lambda,1,smash.poiss,cxx = FALSE)))
    #return(t(apply(lambda,1,function(z){BMSM(z)$E})))
  }

}

###mixed membership
EMupd.mix=function(y,smooth,pi,phi,n,K,B,nugget,nug_control){
  #gamma is nB*K, pi is n*K, phi is K*B, y is n*B
  gamma=pi[rep(1:n,each=B),]*t(phi)[rep(1:B,n),]
  gamma=gamma/rowSums(gamma)
  gamma[is.na(gamma)]=1/K
  gammab=(as.vector(t(y))%o%rep(1,K))*gamma
  pi.num=t(apply(array(gammab,dim=c(B,n,K)),2,colSums))
  #pi.num=(diag(1,n)[,rep(1:n,each=B)])%*%gammab
  pi=pi.num/(rowSums(y)%o%rep(1,K))
  ybt=t(apply(array(gammab,dim=c(B,n,K)),1,colSums))
  #ybt=(diag(1,B)[,rep(1:B,n)])%*%gammab
  #ybw=(diag(1,B)[,rep(1:B,n)])%*%gamma
  phi=t(ybt/(rep(1,B)%o%colSums(gammab)))
  #phi[phi==0]=pseudocounts
  #phi=t(apply(phi,1,normalize))
  #ykb=ybt/ybw
  #ykb[is.na(ykb)]=0
  #ykt=colSums(ykb)

  #lscale=((colSums(ybt)/colSums(pi))%o%rep(1,B))
  lscale=((colSums(ybt))%o%rep(1,B))
  lambda=phi*lscale
  #save(lambda,file="D:/Grad School/projects/sequence_clustering/results/analysis_k562ctcf/debug_lambda.Robj")
  phi.unsmoothed=NULL
  if(smooth==TRUE){
    phi.unsmoothed=phi
    lambda.unsmoothed=lambda
    lambda=smooth.lambda(lambda,nugget,nug_control)
    lambda[is.na(lambda)]=lambda.unsmoothed[is.na(lambda)]
    phi=lambda/lscale
  }


  return(list(pi=pi,phi=phi,phi.unsmoothed=phi.unsmoothed,lambda=lambda,gamma=gamma))
}


negloglik.mix=function(y,pi,phi,n,K,B){
  loglik.ini=log(pi%*%phi)
  yloglik=y*loglik.ini
  yloglik[is.na(yloglik)]=0
  loglik.tot=-sum(yloglik)
  return(loglik.tot)
}



# EMproc.mix=function(y,smooth,pi,phi,n,K,B,tol,maxit){
#   loglik.old=Inf
#   loglik=negloglik.mix(y,pi,phi,n,K,B)
#   cyc=0
#   while(abs(loglik-loglik.old)>tol&cyc<maxit){
#     loglik.old=loglik
#     res=EMupd.mix(y,smooth,pi,phi,n,K,B)
#     pi=res$pi
#     phi=res$phi
#     phi.unsmoothed=res$phi.unsmoothed
#     gamma=res$gamma
#     lambda=res$lambda
#     if(smooth==TRUE){
#       loglik=negloglik.mix(y,pi,phi.unsmoothed,n,K,B)
#     }else{
#       loglik=negloglik.mix(y,pi,phi,n,K,B)
#     }
#     cyc=cyc+1
# #print(cyc)
# #print(pi)
# print(loglik)
#   }
#   return(list(pi=pi,phi=phi,phi.unsmoothed,lambda=lambda,gamma=gamma,loglik=loglik))
# }


rowquantiles = function(x, q) apply(x, 1, quantile, probs = q)

normalized.norm = function(x, y){
  return(rowquantiles(abs(x - y), 0.5)*dim(y)[2])
  #  return(rowquantiles(abs(x - y), 0.95))
}

tol.criterion = function(phi.old, phi, pi.old, pi){
  return(max(max(normalized.norm(phi.old, phi)), mean(normalized.norm(pi.old, pi))))
}

EMproc.mix=function(y,smooth,pi,phi,n,K,B,tol,maxit,nugget,nug_control){
  pi.old=matrix(Inf,nrow=n,ncol=K)
  phi.old=matrix(Inf,nrow=K,ncol=B)
  cyc=0
  while(tol.criterion(phi.old,phi,pi.old,pi)>tol&cyc<maxit){
    pi.old=pi
    phi.old=phi
    res=EMupd.mix(y,smooth,pi,phi,n,K,B,nugget,nug_control)
    pi=res$pi
    phi=res$phi
    phi.unsmoothed=res$phi.unsmoothed
    gamma=res$gamma
    lambda=res$lambda
    if(smooth==TRUE){
      loglik=negloglik.mix(y,pi,phi.unsmoothed,n,K,B)
    }else{
      loglik=negloglik.mix(y,pi,phi,n,K,B)
    }
    cyc=cyc+1
    print("iteration")
    print(cyc)
    #print(pi)
    print("phi difference")
    print(max(normalized.norm(phi.old, phi)))
    print("pi difference")
    print(mean(normalized.norm(pi.old, pi)))
    print("negative loglikelihood")
    print(loglik)
  }
  return(list(pi=pi,phi=phi,phi.unsmoothed,lambda=lambda,gamma=gamma,loglik=loglik))
}


