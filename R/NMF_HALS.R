#'@title HALS algorithm for NMF
#'@param X N by p matrix
#'@param k number of topics
#'@param smooth_method 'trendfiltering', 'wavelet' or 'runmed'
#'@param pos_u,pos_v position of observations
#'@export
#'@import wavethresh

#import genlasso. The igraph package fail to load.

NMF_HALS = function(X,K,smooth_u=F,smooth_v=T,ord=1,
                    smooth_method='trendfiltering',
                    filter.number=1,family='DaubExPhase',
                    pos_u,pos_v,maxiter=100,
                    printevery=10,tol=1e-2){


  quiet = function(x) {
    sink(tempfile())
    on.exit(sink())
    invisible(force(x))
  }

  n = dim(X)[1]
  p = dim(X)[2]
  U = matrix(runif(n*K),nrow=n,ncol=K)
  V = matrix(runif(K*p),nrow=K,ncol=p)

  loss=norm(X-U%*%V)

  for(iter in 1:maxiter){

    for(k in 1:K){
      R_k = X - U[,-k]%*%V[-k,]

      Ru = pmax(t(R_k)%*%U[,k,drop=F],0)

      if(all(Ru==0)){
        V[k,] = rep(0,p)
      }else{
        if(smooth_v){

          V[k,] = smooth_func(Ru/sum((U[,k])^2),smooth_method,pos_v,ord,filter.number,family)

        }else{
          V[k,] = Ru/sum((U[,k])^2)
        }
      }

      Rv = pmax(R_k%*%V[k,],0)

      if(all(Rv==0)){
        U[,k] = rep(0,n)
      }else{
        if(smooth_u){
          U[,k] = smooth_func(Rv/sum((V[k,])^2),smooth_method,pos_u,ord,filter.number,family)
        }else{
          U[,k] = Rv/sum((V[k,])^2)
        }
      }

    }

    loss[iter+1]=norm(X-U%*%V)

    if(iter%%printevery==0){
      print(sprintf('At iter %d, loss: %f',iter,loss[iter+1]))
    }

    if(abs(loss[iter+1]-loss[iter])<=tol){
      break
    }

  }

  return(list(U=U,V=V,loss=loss))

}


smooth_func = function(y,smooth_method,pos,ord,filter.number,family){
  n = length(y)

  if(smooth_method=='trendfiltering'){
    v_k = trendfilter(y,pos,ord=ord)
    v_k.cv = (cv.trendfilter(v_k))
    out = coef(v_k,lambda = v_k.cv$lambda.min)$beta
  }
  if(smooth_method=='wavelet'){

    if(missing(pos)){
      if(!ispowerof2(length(y))){
        reflect.x = reflect(y)
        y = reflect.x$x
        idx = reflect.x$idx
      }else{
        idx = 1:length(y)
      }

      v_k.wd = wd(y,filter.number,family)
      v_k.thresh = threshold(v_k.wd,policy = 'universal')
      out = wr(v_k.thresh)[idx]
    }else{
      grids = makegrid(pos,y)
      v_k.wd = irregwd(grids,filter.number,family)
      v_k.thresh = threshold(v_k.wd,policy = 'universal')
      out = wr(v_k.thresh)
    }
  }
  if(smooth_method=='runmed'){
    win.size = round(sqrt(n)/2)*2+1
    out = runmed(y,win.size)
  }
  out
}

