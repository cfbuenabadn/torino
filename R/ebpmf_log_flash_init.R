#'@title the init fit of flash in ebpmf log function
ebpmf_log_flash_init = function(M,sigma2,l0,f0,ones_n,ones_p,loadings_sign,factors_sign,ebnm_fn,ebnm_fn.offset,
                                S_dim,verbose_flash,fix_l0,fix_f0,Kmax,add_greedy_extrapolate,maxiter_backfitting,
                                backfit_extrapolate,backfit_warmstart,
                                init_fn.flash,no_backfit_kset,n_refit_max,deal_with_no_init_factor,var_type,flash_est_sigma2){

  n = nrow(M)
  p = ncol(M)

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

  if(flash_est_sigma2){
    fit_flash = flash_init(M,var_type = var_type, S=NULL)%>%
      flash_set_verbose(verbose_flash)
  }else{
    fit_flash = flash_init(M,S=sqrt(sigma2),var_type = NULL, S_dim=S_dim)%>%
      flash_set_verbose(verbose_flash)
  }


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

  fit_flash = flash_greedy(fit_flash, Kmax = Kmax,ebnm_fn = ebnm_fn,init_fn=init_fn.flash,extrapolate = add_greedy_extrapolate)
  kset_backfit = (1:fit_flash$n_factors)[!(1:fit_flash$n_factors)%in%no_backfit_kset]
  fit_flash = suppressWarnings(flash_backfit(fit_flash,kset = kset_backfit,maxiter = maxiter_backfitting,extrapolate=backfit_extrapolate,warmstart = backfit_warmstart)%>%
                                 flash_nullcheck(kset=kset_backfit))

  n_refit = 0
  while(fit_flash$n_factors<=2&n_refit<=n_refit_max){
    n_refit = n_refit + 1
    cat(paste('No structure found yet. Re-trying...',n_refit))
    cat('\n')
    init_fn.flash = function(f){flash_greedy_init_default(f, sign_constraints = c(loadings_sign, factors_sign),seed = n_refit)}
    fit_flash = flash_greedy(fit_flash, Kmax = Kmax,ebnm_fn = ebnm_fn,init_fn=init_fn.flash,extrapolate = add_greedy_extrapolate)
    kset_backfit = (1:fit_flash$n_factors)[!(1:fit_flash$n_factors)%in%no_backfit_kset]
    if(n_refit==n_refit_max){
      fit_flash = suppressWarnings(flash_backfit(fit_flash,kset = kset_backfit,maxiter = maxiter_backfitting,extrapolate=backfit_extrapolate,warmstart = backfit_warmstart))
    }else{
      fit_flash = suppressWarnings(flash_backfit(fit_flash,kset = kset_backfit,maxiter = maxiter_backfitting,extrapolate=backfit_extrapolate,warmstart = backfit_warmstart)%>%
                                     flash_nullcheck(kset=kset_backfit))
    }
  }

  if(fit_flash$n_factors<=2){
    cat('No structure found in default initialization.')
    cat('\n')
    if(deal_with_no_init_factor=='reduce_var'){
      cat('Reducing initialization sigma2')
      cat('\n')
      while(fit_flash$n_factors<=2){
        # repeat
        sigma2 = sigma2/2
        fit_flash = flash_init(M,S=sqrt(sigma2),var_type = NULL, S_dim=S_dim)%>%
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
        fit_flash = flash_greedy(fit_flash, Kmax = Kmax,ebnm_fn = ebnm_fn,init_fn=init_fn.flash,extrapolate = add_greedy_extrapolate)
        kset_backfit = (1:fit_flash$n_factors)[!(1:fit_flash$n_factors)%in%no_backfit_kset]
        fit_flash = suppressWarnings(flash_backfit(fit_flash,kset = kset_backfit,maxiter = maxiter_backfitting,extrapolate=backfit_extrapolate,warmstart = backfit_warmstart)%>%
                                       flash_nullcheck(kset=kset_backfit))
      }
    }
    if(deal_with_no_init_factor=='flash_dryrun'){
      cat('Running flash with S=NULL')
      cat('\n')
      fit_flash = flash_init(M,S=NULL,var_type = var_type)%>%
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

      fit_flash = flash_greedy(fit_flash, Kmax = Kmax,ebnm_fn = ebnm_fn,init_fn=init_fn.flash,extrapolate = add_greedy_extrapolate)
      #kset_backfit = (1:fit_flash$n_factors)[!(1:fit_flash$n_factors)%in%no_backfit_kset]
      fit_flash = suppressWarnings(flash_backfit(fit_flash,kset = kset_backfit,maxiter = maxiter_backfitting,extrapolate=backfit_extrapolate,warmstart = backfit_warmstart)%>%
                                     flash_nullcheck(kset=kset_backfit))
    }
  }

  return(fit_flash)
}


