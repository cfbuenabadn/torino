
#'@title simulate Poisson matrix
#'@details The DGP is
#'\deqn{Y\sim Poisson(\exp(l_0f_0'+LF'+E)),}
#'\deqn{l_{0i}\sim Uniform(),f_{0j}\sim Uniform(),}
#'\deqn{l_{ki}\sim N(0,d_k),f_{kj}\sim \pi_0\delta_0 + (1-\pi_0)N(0,\sigma^2_1),}
#'\deqn{e_{ij}\sim N(0,\sigma^2).}
#'@export
sim_data_log = function(n,p,K=3,d=NULL,
                        add_background = TRUE,
                        background_unif_range = c(1,2),
                        pi0=0.8,var1 = 1,
                        var_e = 0,
                        n_simu = 10,
                        seed=12345){
  set.seed(seed)
  if(is.null(d)){
    d = rep(1,K)
  }

  Y = array(dim=c(n,p,n_simu))
  Loading = array(dim=c(n,K,n_simu))
  Factor = array(dim=c(p,K,n_simu))
  L0 = matrix(nrow=n,ncol=n_simu)
  F0 = matrix(nrow=p,ncol=n_simu)

  for(i in 1:n_simu){
    # f is sparse
    FF = matrix(0,nrow=p,ncol=K)
    for(k in 1:K){
      FF[,k] = rnorm(p,0,sd=sqrt(var1))*rbinom(p,1,1-pi0)
    }
    # L controls the norms
    L = matrix(0,nrow=n,ncol=K)
    for(k in 1:K){
      L[,k] = rnorm(n,0,sqrt(d[k]))
    }
    if(add_background){
      l0 = runif(n,min=background_unif_range[1],background_unif_range[2])
      f0 = runif(p,min=background_unif_range[1],background_unif_range[2])
      L0[,i] = l0
      F0[,i] = f0
      S0 = tcrossprod(l0,f0)
    }else{
      S0 = 1
    }
    Factor[,,i] = FF
    Loading[,,i] = L
    Lambda = S0*exp(tcrossprod(L,FF)+matrix(rnorm(n*p,0,sqrt(var_e)),nrow=n,ncol=p))
    Y[,,i] = matrix(rpois(n*p,Lambda),nrow=n,ncol=p)

  }

  return(list(Y=Y,Factor=Factor,Loading=Loading,L0=L0,F0=F0,params = list(n=n,p=p,n_simu=n_simu,
                                                              d=d,
                                                              add_background = add_background,
                                                              background_unif_range = background_unif_range,
                                                              pi0=pi0,var1 = var1,
                                                              var_e = var_e)))

}

#'@title simulate simple Poisson matrix
#'@export
sim_data_log_simple = function(N,p,K=2,l = 20,d=NULL,var_e = 0,phi=1.3){
  if(is.null(d)){
    d = rep(1,K)
  }
  p = max(p,K*l)
  Ftrue = matrix(0,nrow=p,ncol=K)
  for(k in 1:K){
    Ftrue[((k-1)*l+1):((k-1)*l+l),k] = 1
  }

  Ltrue = matrix(rnorm(N*K), ncol=K)
  Ltrue = Ltrue%*%diag(sqrt(d))

  Mu = tcrossprod(Ltrue,Ftrue)
  if(!is.null(phi)){
    var_e = calc_var_PLN_mat(Mu,phi)
  }
  Lambda = exp(Mu+matrix(rnorm(N*p,0,sqrt(var_e)),nrow=N,ncol=p))
  Y = matrix(rpois(N*p,Lambda),nrow=N,ncol=p)
  return(list(Y=Y,Ftrue = Ftrue,Ltrue=Ltrue,var_e = var_e))
}

#'@title simulate Poisson matrix based on real data output
#'@export
#'@importFrom Matrix Matrix
#'@details The DGP is
#'\deqn{Y\sim Poisson(S*\exp(LF'+E))}
sim_data_real = function(L,FF,Sigma2,n_simu = 10,seed=12345){
  set.seed(seed)
  n = nrow(L)
  p = nrow(FF)
  Y = list()
  LF = tcrossprod(L,FF)
  for(i in 1:n_simu){
    Lambda = exp(LF+matrix(rnorm(n*p,0,sqrt(Sigma2)),nrow=n,ncol=p))
    Y[[i]] = Matrix(matrix(rpois(n*p,Lambda),nrow=n,ncol=p),sparse = TRUE)
  }
  return(list(Y=Y,Factor=FF,Loading=L,n_simu=n_simu))
}

calc_var_PLN_mat = function(Mu,phi){
  n = nrow(Mu)
  p = ncol(Mu)
  v = matrix(nrow=n,ncol=p)
  for(i in 1:n){
    for(j in 1:p){
      v[i,j] = calc_var_PLN(Mu[i,j],phi)
    }
  }
  return(v)
}

calc_var_PLN = function(mu,phi){
  a = ((phi-1)*exp(-mu))^2
  return(log(1+Re(polyroot(c(-a,0,1,1))[1])))
}














