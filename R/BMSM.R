#'@title Bayesian Multiscale Model for Smoothing Poisson sequence
#'@param x: a vector of observed counts
#'@param sc: a scale factor/offset, must be a scalar
#'@param bmsm_control: control parameters of BMSM, see bmsm_control_default()
#'@return estimated poisson sequence
#'@importFrom smashr reflect
#'@export


BMSM = function(x,sc=1,bmsm_control=list()){
  if(min(x) < 0){stop ("negative values in x not permitted")}

  if (!is.list(bmsm_control))
    stop("Argument \"bmsm_control\" should be a list")

  control0 = bmsm_control_default()
  if (any(!is.element(names(bmsm_control),names(control0))))
    stop("Argument \"bmsm_control\" contains unknown parameter names")

  control1 = modifyList(control0,bmsm_control,keep.null = TRUE)

  reflect = control1$reflect
  shape = control1$shape
  point_mass=control1$point_mass
  nullweight=control1$nullweight
  g_init = control1$g_init
  fix_g = control1$fix_g
  control =  control1$control
  Elogl = control1$Elogl


  # whether reflect data so that it has a length of powers of 2

  if(!smashr:::ispowerof2(length(x))){
    reflect=TRUE
  }

  if(reflect){
    reflect.x = reflect(x)
    x = reflect.x$x
    idx = reflect.x$idx
  }

  # construct Translation-Invariant table

  titable=ParentTItable(x)

  tit=titable$TItable
  ptit=titable$parent
  n=dim(tit)[2]
  J=dim(tit)[1]-1
  # ns: left child
  # nf: right child
  nt=tit[-1,]
  ns=ptit[,((1:(2*n))%%2==1)]
  nf=nt-ns
  loglik = 0

  # for each scale, perform EB shrinkage.

  post_mean = matrix(0,dim(nt)[1],dim(nt)[2])
  post_mean_log = matrix(0,dim(nt)[1],dim(nt)[2])
  post_mean_log1_p = matrix(0,dim(nt)[1],dim(nt)[2])
  pi_weights = matrix(0,dim(nt)[1],length(shape)+point_mass)

  for(s in 1:dim(ns)[1]){
    spins = 2^(s)
    fit = ebbp_beta_mixture(ns[s,],nt[s,],shape,point_mass,
                            nullweight,weight=rep(1,dim(nt)[2]),
                            g_init, fix_g, control)
    loglik = loglik + (fit$log_likelihood)/spins
    post_mean[s,] = fit$posterior$mean
    post_mean_log[s,] = fit$posterior$mean_log
    post_mean_log1_p[s,] = fit$posterior$mean_log1_p
    if(point_mass){
      pi_weights[s,] = c(fit$fitted_g$pi0,fit$fitted_g$pi)
    }else{
      pi_weights[s,] = fit$fitted_g$pi
    }
  }

  # get estimate of lambda
  est = reverse.pwave(log(tit[J+1,]/sc),log(post_mean))
  est = exp(est)
  if(Elogl){
    est_log = reverse.pwave(log(tit[J+1,]/sc),post_mean_log,post_mean_log1_p)
  }else{
    est_log=NULL
  }

  results = list("E" = est[idx],
                 'Elog' = est_log[idx],
                 "pi_weights" = pi_weights,
                 "loglik" = loglik)

  return(results)

}




#############  Shift operator functions  ##########################

rshift = function(x){L=length(x); return(c(x[L],x[-L]))}

lshift = function(x){return(c(x[-1],x[1]))}


##############  Parent TI table maker (R version)  ######################

ParentTItable=function(sig){
  n = length(sig)
  J = log2(n)

  # Create decomposition table of signal, using pairwise sums,
  # keeping just the values that are *not* redundant under the
  # shift-invariant scheme.  This is very similar to TI-tables
  # in Donoho and Coifman's TI-denoising framework.
  dmat = matrix(0, nrow=J+1, ncol=n)
  dmat[1,] = sig
  #dmat[1,] = as.matrix(sig)
  dmat2 = matrix(0, nrow=J, ncol=2*n) #the parent table

  for(D in 0:(J-1)){
    nD = 2^(J-D);
    nDo2 = nD/2;
    twonD = 2*nD;
    for(l in 0:(2^D-1)){
      ind = (l*nD+1):((l+1)*nD)
      ind2 = (l*twonD+1):((l+1)*twonD)
      x = dmat[D+1,ind]
      lsumx = x[seq(from=1,to=nD-1, by=2)] + x[seq(from=2,to=nD,by=2)]
      rx = rshift(x);
      rsumx = rx[seq(from=1,to=nD-1, by=2)] + rx[seq(from=2,to=nD,by=2)]
      dmat[D+2,ind] = c(lsumx,rsumx)
      dmat2[D+1,ind2] = c(x,rx)
    }
  }
  return(list(TItable=dmat,parent=dmat2))
}


reverse.pwave = function (est, lp, lq = NULL) {
  if (is.null(lq))
    lq = log(1 - exp(lp))
  if (length(est) == 1)
    est = rep(est, ncol(lp))

  J = nrow(lp)

  for (D in J:1) {
    nD = 2^(J - D + 1)
    nDo2 = nD/2
    for (l in 0:(2^(D - 1) - 1)) {
      ind = (l * nD + 1):((l + 1) * nD)

      estvec = est[ind]
      lpvec = lp[D, ind]
      lqvec = lq[D, ind]

      estl = estvec[1:nDo2]
      lpl = lpvec[1:nDo2]
      lql = lqvec[1:nDo2]
      nestl = interleave(estl + lpl, estl + lql)

      estr = estvec[(nDo2 + 1):nD]
      lpr = lpvec[(nDo2 + 1):nD]
      lqr = lqvec[(nDo2 + 1):nD]
      nestr = interleave(estr + lpr, estr + lqr)
      nestr = lshift(nestr)

      est[ind] = 0.5 * (nestl + nestr)
    }
  }
  return(est)
}

interleave=function(x,y){
  return(as.vector(rbind(x,y)))
}


#'@title Default control list of BMSM
#'@param reflect: whether reflect x to make it have length of power of 2.
#'@param shape: shape of Beta prior parameters, default to be 100, 50, 20, 10, 5, 2, 1.
#'@param point_mass: whether put a point mass at 1/2 when performing EB estimate of binomial probability.
#'@param nullweight: nullweight to induce more smoothness.
#'@param control: controls of mixSQP
#'@param Elogl: whether calculate and return E(log(lambda))
#'@export

bmsm_control_default = function(){
  list(reflect = T,
       shape = c(100, 50, 20, 10, 5, 1),
       point_mass=T,
       nullweight=1000,
       g_init = NULL,
       fix_g = FALSE,
       control =  NULL,
       Elogl = TRUE)
}

#
# reflect_kushal = function(x){
#   extended_len <- 2^{ceiling(log(length(x), base=2))}
#   if(extended_len > length(x)){
#     pos_to_fill <- extended_len - length(x)
#     pos_to_fill_1 <- floor(pos_to_fill/2)
#     pos_to_fill_2 <- pos_to_fill - pos_to_fill_1
#     if(pos_to_fill_1 >= 1){
#       x_ext <- c(rev(head(x, pos_to_fill_1)), x, rev(tail(x, pos_to_fill_2)))
#     }else{
#       x_ext <- c(x, rev(tail(x, pos_to_fill_2)))
#     }
#   }else if(extended_len == length(x)){
#     pos_to_fill_1 <- 0
#     x_ext <- x
#   }else{
#     stop("error in extending the vector to make its size a power of 2")
#   }
#   return(list(x_ext=x_ext,idx=(pos_to_fill_1+1):(pos_to_fill_1 + length(x))))
# }
#
#
#
#




