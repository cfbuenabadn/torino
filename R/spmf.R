#'@title Poisson MF but with smoothed factors using smash.pois
#'@description
#'This is similar to Tom's sgom, but based on Poisson likelihood.
#'
#'@importFrom Rfast rowsums
#'@import Matrix
spmf = function(X,K,
                L_init = NULL,
                F_init = NULL,
                smooth_F =TRUE,
                tol=1e-3,
                maxiter=50,
                verbose=TRUE,
                printevery=10){
  # remove columns that are all 0, and are at the start or end of the matrices
  while(sum(X[,1])==0){
    cat('Removed first column that are all 0')
    cat('\n')
    X = X[,-1]
  }
  while(sum(X[,ncol(X)])==0){
    cat('Removed last column that are all 0')
    cat('\n')
    X = X[,-ncol(X)]
  }

  start_time = Sys.time()
  #browser()
  n = dim(X)[1]
  p = dim(X)[2]
  n_points = n*p

  if(is.null(L_init)|is.null(F_init)){
    kmeans.init=kmeans(X,K,nstart=5)
  }
  if(is.null(L_init)) L_init=rep(1,n)%o%normalize(as.vector(table(kmeans.init$cluster)))

  if(is.null(F_init)){
    F_init=t(kmeans.init$centers)
    colnames(F_init)=NULL
  }

  X = Matrix(X,sparse = T)
  x = summary(X)
  non0_idx = cbind(x$i,x$j)

  res = list(l=L_init,f = F_init)
  # alpha = log(res$l[x$i,] + 0.000) + log(res$f[x$j,] + 0.000)
  # exp_offset = rowMaxs(alpha)
  # alpha = alpha - outer(exp_offset,rep(1,K),FUN='*')
  # alpha = exp(alpha)
  alpha = res$l[x$i,] * res$f[x$j,]
  alpha = alpha/rowsums(alpha)
  #alpha = pmax(alpha,1e-8)
  obj = spmf_obj(x$x,tcrossprod(res$l,res$f)[non0_idx])
  for(iter in 1:maxiter){
    for(k in 1:K){
      Ez = calc_EZ(x, alpha[,k])
      res = spmf_update_rank1(Ez$rs,Ez$cs,k,res,smooth_F)
    }
    obj[iter + 1] = spmf_obj(x$x,tcrossprod(res$l,res$f)[non0_idx])
    if(abs(obj[iter+1]-obj[iter])<tol){
      break
    }
    if(verbose){
      if(iter%%printevery==0){
        cat(sprintf('At iter %d, loglik = %f',iter,obj[iter+1]))
        cat('\n')
      }
    }
    # alpha = log(res$l[x$i,] + 0.000) + log(res$f[x$j,] + 0.000)
    # exp_offset = rowMaxs(alpha)
    # alpha = alpha - outer(exp_offset,rep(1,K),FUN='*')
    # alpha = exp(alpha)
    alpha = res$l[x$i,] * res$f[x$j,]
    alpha = alpha/rowsums(alpha)
    #alpha = pmax(alpha,1e-8)
  }
  return(list(l = res$l,f=res$f,obj=obj))
}


#'@title update spmf model, rank 1
#'@importFrom smashr smash.poiss
spmf_update_rank1 = function(l_seq,f_seq,k,res,smooth_F){
  l_scale = sum(res$f[,k])
  res$l[,k] = l_seq/l_scale

  f_scale = sum(res$l[,k])
  if(smooth_F){
    res$f[,k] = smash.poiss(f_seq)/f_scale
  }else{
    res$f[,k] = f_seq/f_scale
  }
  res
}

spmf_obj = function(x,lf){
  #mean(x*log(lf+0.000) - lf)
  mean(dpois(x,lf,log = T))
}








