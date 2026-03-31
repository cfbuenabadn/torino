ebpmf_log_init_control_default = function(){
  return(list(n_cores=1,
              init_tol = 1e-5,
              init_maxiter = 100,
              conv_type='elbo',
              verbose = TRUE,
              printevery = 100,
              ebpm_init = FALSE,
              sigma2_init = NULL,
              M_init = NULL,
              n_refit_flash_init=5,
              flash_est_sigma2=FALSE,
              log_init_for_non0y=FALSE,
              deal_with_no_init_factor = 'flash_dryrun'
              ))
}


ebpmf_log_general_control_default = function(){
  return(list(batch_size = Inf,
              maxiter=50,
              conv_tol=1e-5,
              printevery=10,
              verbose=TRUE,
              save_init_val = FALSE,
              save_latent_M = FALSE,
              garbage_collection_every = 10,
              save_fit_every = Inf,
              save_fit_path = NULL,
              save_fit_name = NULL))
}

ebpmf_log_vga_control_default = function(){
  return(list(maxiter_vga = 10,
         vga_tol = 1e-5))
}

ebpmf_log_sigma2_control_default = function(){
  return(list(est_sigma2 = TRUE,
              a0 = 1,
              b0 = 1,
              cap_var_mean_ratio = 0,
              return_sigma2_trace = FALSE))
}

#'This function returns init.fn.flash and no_backfit_kset
flash_extra_control = function(loadings_sign,factors_sign,fix_l0,fix_f0){
  if(loadings_sign==0&factors_sign==0){
    # this is faster than the default init method in flash
    init.fn.flash = function(f){flash_greedy_init_irlba(f)}
  }else{
    init.fn.flash = function(f){flash_greedy_init_default(f, sign_constraints = c(loadings_sign, factors_sign))}
  }
  no_backfit_kset = NULL
  if(fix_l0){
    no_backfit_kset = c(no_backfit_kset,1)
  }
  if(fix_f0){
    no_backfit_kset = c(no_backfit_kset,2)
  }
  return(list(init.fn.flash=init.fn.flash,no_backfit_kset=no_backfit_kset))
}

#'@title default flash parameters
#'@importFrom ebnm ebnm_point_normal
#'@importFrom ebnm ebnm_normal
ebpmf_log_flash_control_default = function(){
  return(list(ebnm.fn = ebnm_point_normal,
              ebnm.fn.offset = ebnm_normal,
              loadings_sign = 0,
              factors_sign = 0,
              fix_l0=TRUE,
              fix_f0=FALSE,
              Kmax=30,
              add_greedy_Kmax = 1,
              add_greedy_warmstart = TRUE,
              add_greedy_extrapolate = FALSE,
              add_greedy_every = 1,
              backfit_extrapolate = TRUE,
              backfit_warmstart = TRUE,
              maxiter_backfitting = 1,
              verbose_flash=0
  ))
}

#'@title Check loadings and factors signs, to make sure they match ebnm.fn
check_flash_signs = function(flash_control){
  ebnm.fn = flash_control$ebnm.fn
  loadings_sign = flash_control$loadings_sign
  factors_sign = flash_control$factors_sign
  if(length(ebnm.fn)==1){
    ebnm.fn.l = ebnm.fn
    ebnm.fn.f = ebnm.fn
  }else{
    ebnm.fn.l = ebnm.fn[[1]]
    ebnm.fn.f = ebnm.fn[[2]]
  }
  loadings_sign_from_ebnm.fn = test_ebnm_fn_sign(ebnm.fn.l)
  if(loadings_sign_from_ebnm.fn!=loadings_sign){
    warning(paste('Detected umatched loadings_sign from ebnm.fn, setting loadings_sign to', loadings_sign_from_ebnm.fn))
    flash_control$loadings_sign = loadings_sign_from_ebnm.fn
  }
  factors_sign_from_ebnm.fn = test_ebnm_fn_sign(ebnm.fn.f)
  if(factors_sign_from_ebnm.fn!=factors_sign){
    warning(paste('Detected umatched factors_sign from ebnm.fn, setting factors_sign to', factors_sign_from_ebnm.fn))
    flash_control$factors_sign = factors_sign_from_ebnm.fn
  }
  return(flash_control)
}

#'@title test ebnm.fn 's sign
test_ebnm_fn_sign = function(ebnm.fn){
  temp = ebnm.fn(c(-10,0,10),c(1e-3,1e-3,1e-3))
  if(all(temp$posterior$mean>=0)){
    return(1)
  }else if(all(temp$posterior$mean<=0)){
    return(-1)
  }else{
    return(0)
  }
}



