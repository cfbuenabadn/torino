#'@title this function fits flash during ebpmf log iterations
#'@param fit_flash An object from running flash()
ebpmf_log_flash_update = function(fit_flash,sigma2,ones_n,ones_p,iter,loadings_sign,factors_sign,ebnm_fn,ebnm_fn.offset,
                                  S_dim,verbose_flash,fix_l0,fix_f0,Kmax,add_greedy_extrapolate,maxiter_backfitting,
                                  add_greedy_every,add_greedy_Kmax,add_greedy_warmstart,
                                  backfit_extrapolate,backfit_warmstart,
                                  init_fn.flash,no_backfit_kset){

  ## create an init flash_fit obj for flash_init.factor
  ## use flash_fit for init can also init the posterior second moment
  ## Also init the prior g
  flash_fit_init = fit_flash$flash_fit[c('EF','EF2','g')]
  class(flash_fit_init) = c("flash_fit","list")


  l0 = flash_fit_init$EF[[1]][,1,drop=FALSE]
  f0 = flash_fit_init$EF[[2]][,2,drop=FALSE]
  flash_fit_init$EF[[1]] = flash_fit_init$EF[[1]][,-c(1,2),drop=FALSE]
  flash_fit_init$EF2[[1]] = flash_fit_init$EF2[[1]][,-c(1,2),drop=FALSE]
  flash_fit_init$EF[[2]] = flash_fit_init$EF[[2]][,-c(1,2),drop=FALSE]
  flash_fit_init$EF2[[2]] = flash_fit_init$EF2[[2]][,-c(1,2),drop=FALSE]

  ebnm.fixed.l0 = function(x,s,g_init,fix_g,output){
    return(list(posterior=list(mean=l0,second_moment = l0^2),
                fitted_g = NULL,
                log_likelihood=sum(dnorm(x,l0,s,log=T))))
  }

  ebnm.fixed.f0 = function(x,s,g_init,fix_g,output){
    return(list(posterior=list(mean=f0,second_moment = f0^2),
                fitted_g = NULL,
                log_likelihood=sum(dnorm(x,f0,s,log=T))))
  }

  fit_flash = flash_init(fit_flash$flash_fit$Y, S = sqrt(sigma2), var_type = NULL, S_dim = S_dim) %>%
              flash_set_verbose(verbose_flash)

  if(fix_l0){
    fit_flash = flash_factors_init(fit_flash,list(l0, ones_p),ebnm_fn = ebnm.fixed.l0) %>%
                flash_factors_fix(kset = 1, which_dim = "factors")
  }else{
    fit_flash = flash_factors_init(fit_flash,list(l0, ones_p),ebnm_fn = ebnm_fn.offset) %>%
                flash_factors_fix(kset = 1, which_dim = "factors")
  }

  if(fix_f0){
    fit_flash = flash_factors_init(fit_flash,list(ones_n, f0),ebnm_fn = ebnm.fixed.f0) %>%
                flash_factors_fix(kset = 2, which_dim = "loadings")
  }else{
    fit_flash = flash_factors_init(fit_flash,list(ones_n, f0),ebnm_fn = ebnm_fn.offset) %>%
                flash_factors_fix(kset = 2, which_dim = "loadings")
  }
  fit_flash = flash_factors_init(fit_flash,flash_fit_init,ebnm_fn = ebnm_fn)
  fit_flash$flash_fit$g = flash_fit_init$g
  if(iter%%add_greedy_every==0 & fit_flash$n_factors < Kmax){
    fit_flash = flash_greedy(fit_flash, Kmax = add_greedy_Kmax,
                                 ebnm_fn=ebnm_fn,init_fn = init_fn.flash,
                                 warmstart = add_greedy_warmstart,
                                 extrapolate = add_greedy_extrapolate)
  }

  kset_backfit = (1:fit_flash$n_factors)[!(1:fit_flash$n_factors)%in%no_backfit_kset]
  fit_flash = suppressWarnings(flash_backfit(fit_flash, kset = kset_backfit,maxiter = maxiter_backfitting,extrapolate = backfit_extrapolate,warmstart = backfit_warmstart)%>%
                                 flash_nullcheck(kset=kset_backfit))


  fit_flash


}
