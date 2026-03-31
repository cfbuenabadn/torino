#'@title Fit empirical Bayes Poisson matrix factorization with log link function
#'@param Y count data matrix, can be sparse format
#'@param l0,f0 The background loadings and factors, see the model in ‘Details’.
#'@param var_type variance type, "by_row", "by_col" or "constant", see the model in ‘Details’
#'@param general_control A list of parameters controlling the behavior of the algorithm. See ‘Details’.
#'@param vga_control A list of parameters controlling the behavior of the VGA step. See ‘Details’.
#'@param sigma2_control A list of parameters controlling the behavior of updating variance. See ‘Details’.
#'@param flash_control A list of parameters controlling the behavior of the flash step. See ‘Details’.
#'@param verbose TRUE to print the model fitting progress
#'@return A list of:
#'  \item{fit_flash:}{fitted flash object}
#'  \item{elbo:}{evidence lower bound value}
#'  \item{K_trace:}{trace of number of factors}
#'  \item{elbo_trace:}{trace of elbo}
#'  \item{sigma2:}{the variance estimates}
#'  \item{run_time:}{run time of the algorithm}
#'
#'@details The model is
#'\deqn{y_{ij}\sim \text{Poisson}(\exp(\mu_{ij})),}
#'\deqn{\mu_{ij} = l_{i0} + f_{j0} + \sum_k l_{ik}f_{jk} + \epsilon_{ij},}
#'\deqn{l_{i0}\sim g_{l_0}(\cdot), f_{j0}\sim g_{f_0}(\cdot),}
#'\deqn{l_{ik}\sim g_{l_k}(\cdot),f_{jk}\sim g_{f_k}(\cdot),}
#'\deqn{\epsilon_{ij}\sim N(0,\sigma^2_{ij}).}
#'
#'
#'The \code{init_control} argument is a list in which any of the following
#'named components will override the default algorithm settings (as
#'defined by \code{ebpmf_log_init_control_default}):
#'
#'\describe{
#'\item{\code{sigma2_init}}{The init value of sigma2}
#'\item{\code{M_init}}{the initial value for latent M}
#'\item{\code{init_tol}}{tolerance for initialization}
#'\item{\code{init_maxiter}}{max iteration for initialization}
#'\item{\code{verbose}}{TRUE to print initialization progress}
#'\item{\code{printevery}}{Set a number to determine how often to print progress}
#'\item{\code{ebpm_init}}{whether use ash_pois for single gene model, as init for vga}
#'\item{\code{conv_type}}{for init vga fit, use either 'elbo' or 'sigma2abs' for convergence criteria}
#'\item{\code{n_cores}}{Can utilize more than 1 core to perform initialization, using `mclapply` function.}
#'\item{\code{flash_est_sigma2}}{TRUE to use flash for initializing sigma2}
#'\item{\code{log_init_for_non0y}{If TRUE, then for non-0 counts, use log(Y/exp(offset)) as init values}}
#'\item{\code{n_refit_flash_init}}{The times to refit flash using another seed if no structure was found in initialization}
#'\item{\code{deal_with_no_init_factor}}{If no factor found in initialization, use 'reduce_var' to reduce init var for flash, or 'flash_dryrun' for not providing the variance}
#'}
#'
#'The \code{general_control} argument is a list in which any of the following
#'named components will override the default algorithm settings (as
#'defined by \code{ebpmf_log_general_control_default}):
#'\describe{
#'\item{\code{batch_size}}{Set this number to 1000 or 10000 or similar to reduce memory usage for vga step by looping subsets of dataset.}
#'\item{\code{maxiter}}{max iteration allowed.}
#'\item{\code{conv_tol}}{tolerance for convergence}
#'\item{\code{printevery}}{How often to print progress over iterations}
#'\item{\code{verbose}}{TRUE to print progress}
#'\item{\code{garbage_collection_every}}{How often to perform `gc()` to reduce memory usage}
#'\item{\code{save_init_val}}{TRUE to return initailization values of latent mu and sigma2}
#'\item{\code{save_latent_M}}{TRUE to return latent M, its size can be very large}
#'\item{\code{save_fit_every}}{How often to save intermediate results?}
#'\item{\code{save_fit_path}}{Where to save intermediate results path}
#'\item{\code{save_fit_name}}{What is the name to save intermediate results}
#'}
#'
#'The \code{flash_control} argument is a list in which any of the following
#'named components will override the default algorithm settings (as
#'defined by \code{ebpmf_log_flash_control_default}):
#'
#'\describe{
#'
#'\item{\code{ebnm.fn}}{see `?flash`, `ebnm_fn`.}
#'\item{\code{ebnm.fn.offset}}{The prior for \eqn{l_0}, \eqn{f_0}, if not fixing them.}
#'\item{\code{loadings_sign}}{see `?flash_greedy_init_default` sign_constraints, must match ebnm.fn}
#'\item{\code{factors_sign}}{see `?flash_greedy_init_default` sign_constraints, must match ebnm.fn}
#'\item{\code{fix_l0}}{fix  \eqn{l_0}?}
#'\item{\code{fix_f0}}{fix  \eqn{f_0}?}
#'\item{\code{Kmax}}{see `?flash`, `greedy_Kmax`.}
#'\item{\code{add_greedy_Kmax}}{The Kmax in add_greedy in iterations}
#'\item{\code{add_greedy_warmstart}}{see `?flash_greedy`}
#'\item{\code{add_greedy_extrapolate}}{see `?flash_greedy`}
#'\item{\code{add_greedy_every}}{perform flash_greedy every `add_greedy_every` iterations.}
#'\item{\code{maxiter_backfitting}}{max iterations for the flash backfitting,see `?flash_backfit`}
#'\item{\code{backfit_extrapolate}}{see `?flash_backfit`}
#'\item{\code{backfit_warmstart}}{see `?flash_backfit`}
#'\item{\code{verbose_flash}}{whether print flash updates}
#'}
#'
#'The \code{vga_control} argument is a list in which any of the following
#'named components will override the default algorithm settings (as
#'defined by \code{ebpmf_log_vga_control_default}):
#'
#'\describe{
#'\item{\code{maxiter_vga}}{max iterations for vga step Newton's method}
#'\item{\code{vga_tol}}{tolerance for stopping the optimization.}
#'}
#'
#'The \code{sigma2_control} argument is a list in which any of the following
#'named components will override the default algorithm settings (as
#'defined by \code{ebpmf_log_sigma2_control_default}):
#'
#'\describe{
#'\item{\code{est_sigma2}}{whether estimate the variance term or fix it at sigma2_init}
#'\item{\code{a0,b0}}{Inverse-Gamma(a0,b0) prior on sigma2 for regularization.}
#'\item{\code{cap_var_mean_ratio}}{only update sigma2 when if var/mean > (1+cap_var_mean_ratio). i.e. when overdispersion is low enough, stop updating sigma2 to boost convergence.}
#'\item{\code{return_sigma2_trace}}{TRUE to return the sigma2 values along the iterations. internal usage only}
#'}
#'
#'
#'
#'
#'@import flashier
#'@import magrittr
#'@importFrom parallel mclapply
#'@importFrom Matrix Diagonal
#'@importFrom matrixStats colMaxs
#'@importFrom Rfast rowMaxs
#'@export
ebpmf_log = function(Y,l0=NULL,f0=NULL,
                     var_type='by_col',
                     general_control = list(),
                     vga_control = list(),
                     flash_control = list(),
                     sigma2_control = list(),
                     init_control = list(),
                     verbose=TRUE
                     ){

  start_time = Sys.time()

  n = nrow(Y)
  p = ncol(Y)
  num_points = n*p

  if(is.null(l0)){
    l0 = log(cbind(rowMeans(Y)))
  }
  if(is.null(f0)){
    f0 = log(cbind(colSums(Y)/sum(exp(l0))))
  }
  if(length(l0)==1){
    l0 = cbind(rep(l0,n))
  }
  if(length(f0)==1){
    f0 = cbind(rep(f0,p))
  }

  general_control = modifyList(ebpmf_log_general_control_default(),general_control,keep.null = TRUE)
  vga_control = modifyList(ebpmf_log_vga_control_default(),vga_control,keep.null = TRUE)
  flash_control = modifyList(ebpmf_log_flash_control_default(),flash_control,keep.null = TRUE)
  flash_control = check_flash_signs(flash_control)
  flash_control = c(flash_control,flash_extra_control(flash_control$loadings_sign,flash_control$factors_sign,flash_control$fix_l0,flash_control$fix_f0))
  sigma2_control = modifyList(ebpmf_log_sigma2_control_default(),sigma2_control,keep.null = TRUE)
  init_control = modifyList(ebpmf_log_init_control_default(),init_control,keep.null = TRUE)
  init_val = ebpmf_log_init(Y,l0,f0,init_control$sigma2_init,
                            var_type,
                            init_control$M_init,
                            init_control$verbose,
                            init_control$n_cores,
                            init_control$init_tol,
                            init_control$printevery,
                            init_control$ebpm_init,
                            init_control$conv_type,
                            init_control$init_maxiter,
                            init_control$log_init_for_non0y)
  run_time_vga_init = difftime(Sys.time(),start_time,units = 'secs')

  sigma2 = init_val$sigma2_init
  M = init_val$M_init
  if(!general_control$save_init_val){
    init_val = NULL
  }
  gc()

  if(is(Y,'sparseMatrix')){
    const = - sum_lfactorial_sparseMat(Y)
  }else{
    const = - sum(lfactorial(Y))
  }

  if(var_type=='by_row'){
    S.dim = 1
    var_offset_for_obj = p
    var.type = 1
  }else if(var_type=='by_col'){
    S.dim = 2
    var_offset_for_obj = n
    var.type = 2
  }else if(var_type=='constant'){
    S.dim = NULL
    var_offset_for_obj = p*n
    var.type = 0
  }else{
    stop('Non-supported variance type')
  }

  ## split data
  n_batch = ceiling(n/general_control$batch_size)
  if(n_batch>1){
    batches = split(1:n, cut(seq_along(1:n),n_batch , labels = FALSE))
    # transform Y to a list of sub-Y's
    Y = lapply(batches,function(b){Y[b,]})
  }else{
    # To speed up calculation when Y is small and dense Y can be fitted into memory
    Y = as.matrix(Y)
  }

  if(verbose){
    cat('\n')
    cat('Running initial EBMF fit')
    cat('\n')
  }


  t0 = Sys.time()
  ones_n = cbind(rep(1,n))
  ones_p = cbind(rep(1,p))

  ###########################
  # need to change f0 and M if using both non-negative loadings and factors
  # otherwise it is likely no factor can be founded!
  # basically to change baseline f0
  if(flash_control$loadings_sign ==1 & flash_control$factors_sign == 1 & !flash_control$fix_f0){
    f0 = cbind(apply(M-tcrossprod(l0,ones_p),2,min))
  }
  ###########################
  fit_flash = ebpmf_log_flash_init(M,sigma2,l0,f0,ones_n,ones_p,flash_control$loadings_sign,flash_control$factors_sign,
                                   flash_control$ebnm.fn,flash_control$ebnm.fn.offset,
                                   S.dim,flash_control$verbose_flash,flash_control$fix_l0,flash_control$fix_f0,flash_control$Kmax,
                                   flash_control$add_greedy_extrapolate,flash_control$maxiter_backfitting,
                                   flash_control$backfit_extrapolate,flash_control$backfit_warmstart,
                                   flash_control$init.fn.flash,flash_control$no_backfit_kset,
                                   init_control$n_refit_flash_init,
                                   init_control$deal_with_no_init_factor,var.type,init_control$flash_est_sigma2)
  rm(M)
  gc()
  run_time_flash_init =  difftime(Sys.time(),t0,units = 'secs')
  sigma2 = fit_flash$residuals_sd^2
  if(general_control$save_init_val){
    init_val$sigma2_after_flash = sigma2
  }
  K_trace = fit_flash$n_factors
  sigma2_trace = sigma2
  obj = -Inf

  run_time_vga = c()
  run_time_flash = c()

  if(verbose){
    cat('Running iterations...')
    cat('\n')
  }

  for(iter in 1:general_control$maxiter){

    t0 = Sys.time()
    if(n_batch >1){
      v_sum=0
      sym = 0
      ssexp = 0
      slogv = 0
      for(i_b in 1:n_batch){
        b = batches[[i_b]]
        ## this is for speeding up
        y_b = as.matrix(Y[[i_b]])
        res = vga_pois_solver_mat_newton(fit_flash$flash_fit$Y[b,],
                                         y_b,
                                         1,
                                         tcrossprod(fit_flash$flash_fit$EF[[1]][b,],fit_flash$flash_fit$EF[[2]]),
                                         my_ifelse(var_type=='by_row',sigma2[b],adjust_var_shape(sigma2,var_type,length(b),p)),
                                         maxiter=vga_control$maxiter_vga,tol=vga_control$vga_tol,return_V=TRUE)

        fit_flash$flash_fit$Y[b,] = res$M

        ### These are for ELBO calculation later ###
        ############################################################
        sym = sym + sum(y_b*res$M)
        ssexp = ssexp + sum(exp(res$M+res$V/2))
        slogv = slogv + sum(log(res$V)/2+0.9189385)
        ############################################################

        ### This is for updating sigma2, and elbo calculation ###
        ############################################################
        if(var_type=='by_row'){
          v_sum = c(v_sum,rowSums(res$V))
        }else if(var_type=='by_col'){
          v_sum=v_sum+colSums(res$V)
        }else if(var_type=='constant'){
          v_sum=v_sum+sum(res$V)
        }
        ############################################################
      }
      rm(y_b)
      rm(res)
      if(var_type=='by_row'){
        v_sum = v_sum[-1]
      }
      if(sigma2_control$est_sigma2){
        sigma2 = ebpmf_log_update_sigma2(fit_flash,sigma2,v_sum,var_type,
                                         sigma2_control$cap_var_mean_ratio,sigma2_control$a0,sigma2_control$b0,n,p)
      }


    }else{
      res = vga_pois_solver_mat_newton(fit_flash$flash_fit$Y,
                                       Y,
                                       1,
                                       fitted(fit_flash),
                                       adjust_var_shape(sigma2,var_type,n,p),
                                       maxiter = vga_control$maxiter_vga,
                                       tol=vga_control$vga_tol,return_V = TRUE)
      fit_flash$flash_fit$Y = res$M
      if(general_control$save_latent_M){
        V = res$V
      }

      ### These are for ELBO calculation later ###
      ############################################################
      sym = sum(Y*res$M)
      ssexp = sum(exp(res$M+res$V/2))
      slogv = sum(log(res$V)/2+0.9189385)

      ### This is for estimating sigma2
      if(var_type=='constant'){
        v_sum =sum(res$V)
      }else if(var_type=='by_col'){
        v_sum =colSums(res$V)
      }else if(var_type=='by_row'){
        v_sum =rowSums(res$V)
      }
      ############################################################
      if(sigma2_control$est_sigma2){
        sigma2=ebpmf_log_update_sigma2(fit_flash,sigma2,v_sum,var_type,
                                       sigma2_control$cap_var_mean_ratio,sigma2_control$a0,sigma2_control$b0,n,p)
      }
      rm(res)
    }

    run_time_vga[iter] = difftime(Sys.time(),t0,units='secs')
    if(sigma2_control$return_sigma2_trace){
      sigma2_trace = rbind(sigma2_trace,sigma2)
    }

    ## fit flash
    t0 = Sys.time()
    fit_flash = ebpmf_log_flash_update(fit_flash,sigma2,ones_n,ones_p,iter,flash_control$loadings_sign,flash_control$factors_sign,
                                       flash_control$ebnm.fn,flash_control$ebnm.fn.offset,
                                       S.dim,flash_control$verbose_flash,flash_control$fix_l0,flash_control$fix_f0,flash_control$Kmax,
                                       flash_control$add_greedy_extrapolate,flash_control$maxiter_backfitting,flash_control$add_greedy_every,
                                       flash_control$add_greedy_Kmax,flash_control$add_greedy_warmstart,
                                       flash_control$backfit_extrapolate,flash_control$backfit_warmstart,
                                       flash_control$init.fn.flash,flash_control$no_backfit_kset)
    run_time_flash[iter] = difftime(Sys.time(),t0,units='secs')
    K_trace[iter+1] = fit_flash$n_factors
    KL_LF = sum(flash_fit_get_KL(fit_flash$flash_fit,1)) + sum(flash_fit_get_KL(fit_flash$flash_fit,2))

    # check convergence
    obj[iter + 1] = calc_ebpmf_log_obj(n,p,sym,ssexp,slogv,v_sum,sigma2,fit_flash$flash_fit$R2,KL_LF,const,var_offset_for_obj,sigma2_control$a0,sigma2_control$b0)
    if((obj[iter + 1]-obj[iter])/num_points< general_control$conv_tol){
      break
    }
    if(verbose){
      if(iter%%general_control$printevery==0){
        cat(paste('iter ',iter, ', avg elbo=',round(obj[iter+1]/num_points,log10(1/general_control$conv_tol)),", K=",fit_flash$n_factors,sep = ''))
        cat('\n')
      }
    }


    if(iter%%general_control$save_fit_every==0){
      saveRDS(list(fit_flash=list(L_pm=fit_flash$L_pm,F_pm = fit_flash$F_pm,pve = fit_flash$pve),
                   elbo=obj[length(obj)],
                   K_trace=K_trace,
                   elbo_trace=obj,
                   sigma2 = sigma2,
                   sigma2_trace = sigma2_trace,
                   run_time = difftime(Sys.time(),start_time,units='auto'),
                   run_time_break_down = list(run_time_vga_init = run_time_vga_init,
                                              run_time_flash_init = run_time_flash_init,
                                              run_time_vga = run_time_vga,
                                              run_time_flash = run_time_flash)),
              file=paste(general_control$save_fit_path,general_control$save_fit_name,'_iter',iter,'.rds',sep=''))
    }

    if(iter%%general_control$garbage_collection_every==0){gc()}


  }

  end_time = Sys.time()
  if(!general_control$save_latent_M){fit_flash$flash_fit$Y = NULL}
  if(general_control$save_latent_M){fit_flash$flash_fit$V = V}

  return(list(fit_flash=fit_flash,
              elbo=obj[length(obj)],
              K_trace=K_trace,
              elbo_trace=obj,
              sigma2 = sigma2,
              sigma2_trace = sigma2_trace,
              init_val=init_val,
              run_time = difftime(end_time,start_time,units='auto'),
              run_time_break_down = list(run_time_vga_init = run_time_vga_init,
                                         run_time_flash_init = run_time_flash_init,
                                         run_time_vga = run_time_vga,
                                         run_time_flash = run_time_flash)))
}

#'@title Calc elbo
calc_ebpmf_log_obj = function(n,p,sym,ssexp,slogv,sv,sigma2,R2,KL_LF,const,ss,a0,b0){
  val = sym - ssexp + slogv + 0.5*n*p - sum(ss*log(2*pi*sigma2)/2)- sum(sv/2/sigma2) - sum(R2/2/sigma2) + const+ KL_LF - sum((a0+1)*log(sigma2)) - sum(b0/sigma2)
  # val = sum(Y*M - S*exp(M+V/2)   + log(2*pi*V)/2 + 0.5 ) - sum(ss*log(2*pi*sigma2)/2)- sum(sv/2/sigma2) - sum(R2/2/sigma2) + const+ KL_LF
  return(val)
}

#'@title update sigma2
ebpmf_log_update_sigma2 = function(fit_flash,sigma2,v_sum,var_type,cap_var_mean_ratio,a0,b0,n,p){
  if(var_type=='constant'){
    if(cap_var_mean_ratio>0){
      # ((exp(sigma2)-1)*exp(max(fitted(fit_flash))))>cap_var_mean_ratio
      if(max(fitted(fit_flash)) > (log(cap_var_mean_ratio)-log((exp(sigma2)-1))-sigma2/2)){
        sigma2 = ((v_sum +sum(fit_flash$flash_fit$R2))/2+b0)/(n*p/2+a0+1)
      }
    }else{
      sigma2 = ((v_sum +sum(fit_flash$flash_fit$R2))/2+b0)/(n*p/2+a0+1)
    }
  }else if(var_type=='by_row'){
    if(cap_var_mean_ratio>0){
      update_idx = which(rowMaxs(fitted(fit_flash)) > (log(cap_var_mean_ratio)-log((exp(sigma2)-1))-sigma2/2))
      #update_idx = which(((exp(sigma2)-1)*exp(rowMaxs(fitted(fit_flash))))>cap_var_mean_ratio)
      if(!is.null(update_idx)){
        sigma2[update_idx] = (((v_sum+fit_flash$flash_fit$R2)/2+b0)/(p/2+a0+1))[update_idx]
      }
    }else{
      sigma2 = ((v_sum+fit_flash$flash_fit$R2)/2+b0)/(p/2+a0+1)
    }
  }else if(var_type=='by_col'){
    if(cap_var_mean_ratio>0){
      update_idx = which(colMaxs(fitted(fit_flash)) > (log(cap_var_mean_ratio)-log((exp(sigma2)-1))-sigma2/2))
      #update_idx = which(((exp(sigma2)-1)*exp(colMaxs(fitted(fit_flash))))>cap_var_mean_ratio)
      if(!is.null(update_idx)){
        sigma2[update_idx] = (((v_sum+fit_flash$flash_fit$R2)/2+b0)/(n/2+a0+1))[update_idx]
      }
    }else{
      sigma2 = ((v_sum+fit_flash$flash_fit$R2)/2+b0)/(n/2+a0+1)
    }
  }else{
    stop('Non-supported var type')
  }
  return(sigma2)
}

# V_M = function(M,Y,EL,EF,sigma2,var_type='by_col'){
#   if(var_type=='by_col'){
#     p = ncol(M)
#     return(colSums((1/(Y%*%Diagonal(p,sigma2)+tcrossprod(EL,EF)-M+1))%*%Diagonal(p,sigma2)))
#   }
#   if(var_type=='by_row'|var_type=='constant'){
#     return(rowSums(sigma2/(simga2*Y-M+tcrossprod(EL,EF)+1)))
#   }
# }

#'@title Calc V given M
V_M = function(M,Y,Beta,Sigma2,var_type='by_col'){
  if(var_type=='by_col'){
    return(colSums(Sigma2/(Y*Sigma2+Beta-M+1)))
  }
  if(var_type=='by_row'|var_type=='constant'){
    return(rowSums(Sigma2/(Sigma2*Y-M+Beta+1)))
  }
}

#'@title calculate sum of log factorial of all elements of a sparse Matrix
sum_lfactorial_sparseMat = function(Y){
  sum(lfactorial(Y@x))
}

#'@title adjust var shape for vga. For computation purpose.
adjust_var_shape = function(sigma2,var_type,n,p){
  if(var_type=='constant'){
    sigma2 = matrix(sigma2,nrow=n,ncol=p)
  }else if(var_type=='by_row'){
    sigma2 = matrix(sigma2,nrow=n,ncol=p,byrow = F)
  }else if(var_type=='by_col'){
    sigma2 = matrix(sigma2,nrow=n,ncol=p,byrow = T)
  }else{
    stop('Non-supported var type')
  }
  sigma2
}

my_ifelse = function(test,yes,no){
  if(test){
    return(yes)
  }else{
    return(no)
  }
}

#'@title get posterior variance of b_ij = sum_j l_ik * f_jk
get_var_b = function(fit_flash){
  n = nrow(fit_flash$flash_fit$Y)
  p = ncol(fit_flash$flash_fit$Y)
  Vb = matrix(nrow=n,ncol=p)
  for(i in 1:n){
    for(j in 1:p){
      temp = fit_flash$L_pm[i,]^2*(fit_flash$F_psd[j,]^2) + fit_flash$F_pm[j,]^2*(fit_flash$L_psd[i,]^2)
      Vb[i,j] = sum(temp)
    }
  }
  Vb
}

