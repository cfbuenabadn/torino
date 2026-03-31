#'@title Empirical Bayes for estimating Binomial probability
#'@import mixsqp
#'@param x: a vector of the number of successes
#'@param nx: a vector of the number of trials.
#'@param shape: parameter of mixture Beta component, Beta(shape,shape).
#'@param point_mass: whether include delta(1/2) in the mixture component.
#'@param nullweight: extra amount of weight on the first mixing proportion.
#'@param weight: a vector of weights
#'@param control: controls of mixSQP
#'@return prior, posterior, likelihood
#'@export
#'

ebbp_beta_mixture = function(x,nx,shape=NULL,point_mass=TRUE,
                             nullweight=10, weight = rep(1,length(x)),
                             g_init = NULL, fix_g = FALSE, control =  NULL){
  n=length(x)
  if(is.null(control)){control = mixsqp_control_defaults()}
  if(is.null(shape)){
    shape = c(100, 50, 20, 10, 5, 2, 1)
  }
  if(is.null(g_init)){
    fix_g = FALSE ## then automatically unfix g if specified so
    if(point_mass){
      pi_init=rep(1/(length(shape)+1),length(shape)+1)
      g_init = list(pi=pi_init[-1],shape=shape,pi0=pi_init[1])
    }else{
      pi_init = rep(1/(length(shape)),length(shape))
      g_init = list(pi=pi_init,shape=shape)
    }
  }
  if(!fix_g){
    tmp =  compute_L_bbinom(x,nx,g_init$shape,point_mass)
    if(point_mass){x0 = c(g_init$pi0,g_init$pi)}else{x0 = g_init$pi}
    if(!is.null(nullweight)){
      Lnull = rbind(c(1,rep(0,ncol(tmp$L)-1)),tmp$L)
      weight = c(nullweight-1,weight)
      fit = mixsqp(Lnull, weight,x0 = x0, control = control)
    }else{
      fit = mixsqp(tmp$L, weight,x0 = x0, control = control)
    }
    pi = fit$x
    pi = pi/sum(pi)
  }else{
    if(point_mass){
      pi = c(g_init$pi0,g_init$pi)
    }else{
      pi = g_init$pi
    }
    tmp = compute_L_bbinom(x,nx,g_init$shape,point_mass)
  }
  fitted_g = betabinommix(pi,shape,point_mass)
  log_likelihood = sum(log(exp(tmp$l_rowmax) * tmp$L %*%  pi))
  # calculate posterior mean, posterior expectation of log
  ## posterior mean of each mixture component
  pos_alpha = outer(x,shape,'+')
  pos_beta = outer(nx-x,shape,'+')
  pm_mixture_comp = pos_alpha/(pos_alpha+pos_beta)
  if(point_mass){pm_mixture_comp = cbind(1/2,pm_mixture_comp)}
  ## posterior mixing proportion
  Pi_tilde = t(t(tmp$L) * pi)
  Pi_tilde = Pi_tilde/rowSums(Pi_tilde)

  ## posterior mean
  pm = rowSums(Pi_tilde * pm_mixture_comp)

  ## posterior expectation of log(p)

  dab = digamma(pos_alpha+pos_beta)

  pm_mixture_comp_log = digamma(pos_alpha)-dab

  ## posterior expectation of log(1-p)

  pm_mixture_comp_log2 = digamma(pos_beta)-dab

  if(point_mass){
    pm_mixture_comp_log = cbind(rep(log(1/2),n),pm_mixture_comp_log)
    pm_mixture_comp_log2 = cbind(rep(log(1/2),n),pm_mixture_comp_log2)
  }
  Elogp = rowSums(Pi_tilde * pm_mixture_comp_log)
  Elog1_p = rowSums(Pi_tilde * pm_mixture_comp_log2)
  posterior = data.frame(mean = pm,mean_log=Elogp,mean_log1_p=Elog1_p)

  return(list(fitted_g=fitted_g,posterior=posterior,
              log_likelihood=log_likelihood,Pi_tilde=Pi_tilde))

}

betabinommix = function(pi,shape,point_mass){
  if(point_mass){
    structure(list(pi = pi[-1], shape = shape, pi0 = pi[1]), class="betabinomial")
  }else{
    structure(list(pi = pi, shape = shape,pi0=NULL), class="betabinomial")
  }
}

compute_L_bbinom = function(x,nx,shape,point_mass){
  n=length(x)
  shapes = outer(rep(1,n),shape,'*')
  l = dbetabinom(outer(x,rep(1,length(shape))),outer(nx,rep(1,length(shape))),shapes,shapes)
  l_rowmax  = apply(l,1,max)
  if(point_mass){
    l0 = cbind(dbinom_log(x,nx,1/2),l)
    L = exp(l0 -  l_rowmax)
  }else{
    L = exp(l -  l_rowmax)
  }
  return(list(L = L, l_rowmax = l_rowmax))
}

# returns log density of binomial
dbinom_log = function(x,nx,p){
  lchoose(nx,round(x))+x*log(p)+(nx-x)*log(1-p)
}

# returns log density of beta-binomial distribution
dbetabinom = function(x,size,a,b){
  lchoose(size,round(x))+lbeta(x+a,size-x+b)-lbeta(a,b)
}

mixsqp_control_defaults <- function() {
  return(list(verbose = F))
}

