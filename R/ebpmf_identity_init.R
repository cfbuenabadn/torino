#'@title initialize the ebpmf identity model
#'@param X input data matrix
#'@param K number of topics
#'@param init init methods, or a list of init L and F
#'@importFrom fastTopics fit_poisson_nmf
#'@importFrom fastTopics init_poisson_nmf
#'@importFrom Rfast rowsums
#'@importFrom Rfast rowMaxs
#'@export
ebpmf_identity_init = function(X,
                               K,
                               init,
                               maxiter_init = 50,
                               lib_size){

  n = nrow(X)
  p = ncol(X)
  if(is.list(init)){
    L_init = init$L_init
    F_init = init$F_init
    cat('F_init should be the EF from ebpmf_identity fit...')
    cat('\n')

    if(is.null(L_init)&!is.null(F_init)){
      init_fasttopic = init_poisson_nmf(X,F=F_init,init.method='random')
      init_fit = fit_poisson_nmf(X,numiter = maxiter_init,fit0 = init_fasttopic)
      L_init = init_fit$L
    }

  }else{

    if(init == 'uniform'){
      L_init = matrix(runif(n*K),nrow=n,ncol=K)
      F_init = matrix(runif(K*p),nrow=p,ncol=K)
      ratio = median(X)/(median(L_init)*median(F_init))
      L_init = L_init*sqrt(ratio)
      F_init = F_init*sqrt(ratio)
    }
    # if(init%in%c('scd','lee')){
    #   X_init_fit = nnmf(as.matrix(X),K,method=init,loss='mkl',show.warning = F,verbose = F,max.iter = maxiter_init)
    #   L_init = X_init_fit$W
    #   F_init = t(X_init_fit$H)
    # }
    if(init == 'kmeans'){
      kmeans.init=kmeans(as.matrix(X),K,nstart=5)
      L_init = rep(1,n)%o%normalize(as.vector(table(kmeans.init$cluster)))
      F_init = t(kmeans.init$centers)
      row.names(F_init)=NULL
    }
    if(init == 'fasttopics'){
      init_fasttopic = fit_poisson_nmf(X,K,numiter = maxiter_init,init.method = 'topicscore')
      L_init = init_fasttopic$L
      F_init = init_fasttopic$F
    }
    if(init == 'fasttopics_random'){
      init_fasttopic = fit_poisson_nmf(X,K,numiter = maxiter_init,init.method = 'random')
      L_init = init_fasttopic$L
      F_init = init_fasttopic$F
    }
    if(init == 'fasttopics_em'){
      init_fasttopic = fit_poisson_nmf(X,K,numiter = maxiter_init,method = 'em')
      L_init = init_fasttopic$L
      F_init = init_fasttopic$F
    }
    if(init == 'fasttopics_ccd'){
      init_fasttopic = fit_poisson_nmf(X,K,numiter = maxiter_init,method = 'ccd')
      L_init = init_fasttopic$L
      F_init = init_fasttopic$F
    }
    if(init == 'fasttopics_mu'){
      init_fasttopic = fit_poisson_nmf(X,K,numiter = maxiter_init,method = 'mu')
      L_init = init_fasttopic$L
      F_init = init_fasttopic$F
    }
  }

  # adjust scale of L and F, mainly for stability.
  # ratio = adjLF(L_init,F_init)
  # L_init = ratio$L_init
  # F_init = ratio$F_init
  ratio = poisson_to_libsize(F_init,L_init,lib_size)
  L_init = ratio$L
  F_init = ratio$FF

  gl = list()
  gf = list()

  ql = list(El = L_init, Elogl = log(L_init+1e-8))
  qf = list(Ef = F_init, Elogf = log(F_init+1e-8))



  return(list(ql=ql,
              qf=qf,
              gl=gl,
              gf=gf,
              Hl = rep(0,K),
              Hf = rep(0,K),
              lib_size=lib_size))

}



# additional init function for smoothing F
ebpmf_identity_init_smooth = function(res,K,p,x,alpha,ebps_control,maxiter_init,ebpm.fn.f){
  res$gf$sigma2 = rep(0,K)
  res$gf$sigma2_trace = rep(0,K)
  res$qf$Ef_smooth= matrix(nrow=p,ncol=K)
  res$qf$Elogf_smooth  = matrix(nrow=p,ncol=K)

  #browser()

  for(k in 1:K){
    Ez = calc_EZ(x, alpha[,k])
    sk = ebpm.fn.f(Ez$cs,sum(res$lib_size*res$ql$El[,k]),
                   init_control=list(m_init_method = ebps_control$m_init_method_for_init),
                   general_control = list(maxiter=maxiter_init,
                                          #convergence_criteria = 'nugabs',
                                          maxiter_vga = ebps_control$maxiter_vga,
                                          make_power_of_2=ebps_control$make_power_of_2,
                                          vga_tol=ebps_control$vga_tol,
                                          tol = ebps_control$tol
                                          ),
                   smooth_control = list(wave_trans=ebps_control$wave_trans,
                                         ndwt_method = ebps_control$ndwt_method,
                                         filter.number = ebps_control$filter.number,
                                         family = ebps_control$family,
                                         ebnm_params=ebps_control$ebnm_params,
                                         warmstart=ebps_control$warmstart))
    res$gf$sigma2[k] = sk$fitted_g$sigma2
    res$gf$sigma2_trace[k] = sk$fitted_g$sigma2
    res$qf$Ef_smooth[,k] = sk$posterior$mean_smooth
    res$qf$Elogf_smooth[,k] = sk$posterior$mean_log_smooth
  }
  res
}


normalize=function(x) x/sum(x)


#'
#' #'@title initialize the ebpmf identity model
#' #'@param X input data matrix
#' #'@param K number of topics
#' #'@param init init methods, or a list of init L and F
#' #'@param init_loss mkl or mse
#' #'@importFrom NNLM nnmf
#' #'@export
#' ebpmf_identity_init = function(X,K,init,init_loss,maxiter_init = 50){
#'
#'   n = nrow(X)
#'   if(is.list(init)){
#'     L_init = init$L_init
#'     F_init = init$F_init
#'
#'     if(is.null(L_init)){
#'       X_init_fit = nnmf(as.matrix(X),K,method='lee',
#'                         loss='mse',show.warning = F,
#'                         init = list(H=t(F_init)),
#'                         verbose = F,max.iter = maxiter_init)
#'       L_init = X_init_fit$W
#'     }
#'
#'   }else{
#'
#'     if(init%in%c('scd','lee')){
#'       X_init_fit = NNLM::nnmf(as.matrix(X),K,method=init,loss=init_loss,show.warning = F,verbose = F,max.iter = maxiter_init)
#'       L_init = X_init_fit$W
#'       F_init = t(X_init_fit$H)
#'     }
#'     if(init == 'uniform'){
#'       L_init = matrix(runif(n*K),nrow=n,ncol=K)
#'       F_init = matrix(runif(K*p),nrow=p,ncol=K)
#'       ratio = median(X)/(median(L_init)*median(F_init))
#'       L_init = L_init*sqrt(ratio)
#'       F_init = F_init*sqrt(ratio)
#'     }
#'     if(init == 'kmeans'){
#'       kmeans.init=kmeans(as.matrix(X),K,nstart=5)
#'       L_init = rep(1,n)%o%normalize(as.vector(table(kmeans.init$cluster)))
#'       F_init = t(kmeans.init$centers)
#'       row.names(F_init)=NULL
#'     }
#'   }
#'
#'   # adjust scale of L and F, mainly for stability.
#'   ratio = adjLF(L_init,F_init)
#'   L_init = ratio$L_init
#'   F_init = ratio$F_init
#'
#'   gl = list(sigma2=NULL)
#'   gf = list(sigma2=NULL)
#'
#'   ql = list(El = L_init, Elogl = log(L_init+1e-10), Esmooth_l=L_init)
#'   qf = list(Ef = F_init, Elogf = log(F_init+1e-10), Ef_smooth=F_init)
#'
#'   return(list(ql=ql,
#'               qf=qf,
#'               gl=gl,
#'               gf=gf,
#'               Hl = rep(0,K),
#'               Hf = rep(0,K)))
#'
#' }
#'
