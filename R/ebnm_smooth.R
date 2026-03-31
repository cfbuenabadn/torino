#'@title Empirical Bayes wavelet smoothing DWT wrapper for flashier
#'@description Smooth homogeneous Gaussian data.
#'@param x data
#'@param s known standard error
#'@param g_init a list of priors for each scale
#'@param fix_g whether fix prior
#'@return a list of
#'  \item{posterior:}{posterior mean and 2nd moment}
#'  \item{fitted_g:}{estimated prior}
#'  \item{log_likelihood:}{log likelihood}
#'@import wavethresh
#'@importFrom ebnm ebnm
#'@importFrom smashr reflect
#'@export
ebnm_dwt = function(x, s, g_init, fix_g,
                    filter.number=1,
                    family="DaubExPhase"){


  ebnm_params=list()
  W=NULL

  n = length(x)
  J = log(n,2)
  n_orig = n
  if(ceiling(J)!=floor(J)){
    #stop('Length of x must be power of 2')
    # reflect
    x = reflect(x)
    idx = x$idx
    x = x$x
    n = length(x)
    J = log(n,2)
  }else{
    idx = 1:n
  }

  if(filter.number==1&family=='DaubExPhase'){
    wavelet_name = "haar"
  }else{
    wavelet_name = 'non-haar'
  }
  ebnm_params = modifyList(ebnm_params_default_smooth(),ebnm_params,keep.null =  TRUE)
  tsum = sum(x)/sqrt(n)
  x.w = wd(x, filter.number = filter.number,
           family = family, type = "wavelet")

  data.var = s^2
  if(length(data.var==1)){
    data.var = rep(data.var,n)
  }else if(length(unique(data.var))==1){
    data.var = rep(data.var[1],n)
  }else{
    stop('sigma must be constant for all observations.')
  }

  if(wavelet_name!='haar'){
    if(is.null(W)){
      W = (t(GenW(n,filter.number,family)))[-1,]
    }
  }



  x.w.v =  data.var
  tsum.var = x.w.v[1]
  x.w.v = x.w.v[-1]


  dKL = 0
  loglik.scale = c()
  x.w.v.s = rep(0, 2^J-1)
  fitted_g = list()
  for (j in 0:(J - 1)) {
    x.pm = rep(0, 2^j)
    #index = (((J - 1) - j) * n + 1):((J - j) * n)
    index = (n-2^(j+1)+1):(n-2^j)
    x.w.j = accessD(x.w, j)
    x.w.v.j = x.w.v[index]
    ind.nnull = (x.w.v.j != 0)

    a = ebnm(x.w.j[ind.nnull],sqrt(x.w.v.j[ind.nnull]),
             prior_family=ebnm_params$prior_family,
             g_init = g_init[[j+1]],
             fix_g = fix_g,
             optmethod = ebnm_params$optmethod,
             control = ebnm_params$control)
    fitted_g[[j+1]] = a$fitted_g
    dKL = dKL + a$log_likelihood - Eloglik(x.w.j[ind.nnull], sqrt(x.w.v.j[ind.nnull]),a$posterior$mean, a$posterior$mean^2+a$posterior$sd^2)
    x.pm[ind.nnull] = a$posterior$mean
    x.pm[!ind.nnull] = 0
    x.w = putD(x.w, j, x.pm)
    loglik.scale[j + 1] = a$log_likelihood
    x.w.v.s[index[ind.nnull]] = a$posterior$sd^2
    x.w.v.s[index[!ind.nnull]] = 0
  }
  mu.est = wr(x.w)
  loglik = sum(loglik.scale)
  #x.w.v.s = c(tsum.var,x.w.v.s)
  if(wavelet_name=='haar'){
    mu.est.var = haar_inv_var(c(x.w.v.s,0))
  }else{
    mu.est.var = colSums(W^2*x.w.v.s)
  }

  return(list(posterior=data.frame(mean = mu.est[idx],second_moment=(mu.est^2+mu.est.var)[idx]),
              fitted_g=fitted_g,
              log_likelihood = loglik/n*n_orig))
}


#'@title Empirical Bayes wavelet smoothing haar wavelet for flashier
#'@export
ebnm_dwt_haar = function(x, s, g_init=NULL, fix_g=FALSE, output){
  ebnm_dwt(x, s, g_init, fix_g,filter.number=1,family="DaubExPhase")
}

#'@title Empirical Bayes wavelet smoothing symlet wavelet for flashier
#'@export
ebnm_dwt_symlet = function(x, s, g_init=NULL, fix_g=FALSE, output){
  ebnm_dwt(x, s, g_init, fix_g,filter.number=10,family="DaubLeAsymm")
}


ebnm_params_default_smooth = function(){
  return(list(prior_family='normal_scale_mixture',
              optmethod = NULL))
}

Eloglik = function(x, s, Et, Et2) {
  # Deal with infinite SEs:
  idx = is.finite(s)
  x = x[idx]
  s = s[idx]
  Et = Et[idx]
  Et2 = Et2[idx]
  return(-0.5 * sum(log(2*pi*s^2) + (1/s^2) * (Et2 - 2*x*Et + x^2)))
}

haar = function(x,scale= sqrt(2)){
  if(length(x)==1){
    return(x)
  }
  else{
    x = matrix(x,nrow=2)
    diff = (x[1,]-x[2,])/scale
    sum = (x[1,]+x[2,])/scale
    return(c(diff, haar(sum)))
  }
}

haar_inv = function(x,scale=sqrt(2)){
  n=length(x)
  if(n==1){
    return(x)
  }
  x = matrix(scale*x,nrow=2,byrow=TRUE)
  smoothed = haar_inv(x[2,])
  return(as.vector(rbind(smoothed+x[1,], smoothed-x[1,]))/2)
}

haar_inv_var = function(v,scale=sqrt(2)){
  n=length(v)
  if(n==1){
    return(v)
  }
  v = matrix(scale^2*v,nrow=2,byrow=TRUE)
  smoothed = haar_inv_var(v[2,])
  return(as.vector(rbind(smoothed+v[1,], smoothed+v[1,]))/4)
  #return(rep((smoothed+v[1,])/4,each=2))
}



#'@title This function is an EBNM function that fits ndwt method
#'@importFrom smashr smash.gaus
#'@export
ebnm_ndwt = function(x, s, g_init, fix_g, output){


  fit = smashr::smash.gaus(x,sigma=s,post.var = T,return.loglr = T,reflect = F)


  return(list(posterior=data.frame(mean = fit$mu.est,second_moment=(fit$mu.est^2+fit$mu.est.var)),
              fitted_g=list(),
              log_likelihood = fit$loglik))
}





