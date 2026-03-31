#'@title this function initializes the ebpmf with log link
#'@description sigma2, M needs initialization.
#'@importFrom vebpm ebpm_normal
#'@importFrom ashr ash_pois

ebpmf_log_init = function(Y,l0,f0,sigma2,
                          var_type,
                          M_init,
                          verbose,
                          n_cores,
                          init_tol,
                          printevery,
                          ebpm_init=TRUE,
                          conv_type,
                          init_maxiter,
                          log_init_for_non0y,
                          L_init=NULL,
                          F_init=NULL){
  n = nrow(Y)
  p = ncol(Y)
  # we need a. M; b. sigma2
  # case 1: both are given in input.
  if(!is.null(sigma2) & !is.null(M_init)){
    return(list(sigma2_init=sigma2,M_init=M_init))
  }

  # then we fit poisGG with mean = 0, var =sigma2
  if(is.null(M_init)&is.null(L_init)&is.null(F_init)){
    if(verbose){
      cat('Initializing')
      cat('\n')
    }
    # pre-estimate sigma2, assuming LF = 0?.
    if(var_type=='constant'){
      if(verbose){
        cat('Solving VGA constant...For large matrix this may require large memory usage')
      }
      if(ebpm_init){
        #init_val = ebpm_exponential_mixture(as.vector(Y),s = exp(as.vector(outer(l0,f0,FUN='+'))))
        #init_var_vga = mean(init_val$posterior$mean_log^2)
        init_val = ash_pois(as.vector(Y),scale = exp(as.vector(outer(l0,f0,FUN='+'))),link='identity',mode=1,method='shrink')
        #init_var_vga = mean(log(init_val$result$PosteriorMean)^2)
        init_var_vga = NULL
        init_val = suppressWarnings(ebpm_normal(as.vector(Y),
                                                g_init = list(mean=as.vector(outer(l0,f0,FUN='+')),var=init_var_vga),
                                                q_init = list(m_init = log(init_val$result$PosteriorMean) + as.vector(outer(l0,f0,FUN='+')),v_init=init_var_vga),
                                                fix_g = c(TRUE,FALSE),tol=init_tol,conv_type = conv_type,
                                                maxiter = init_maxiter))
      }else{
        init_val = suppressWarnings(ebpm_normal(as.vector(Y),
                                                g_init = list(mean=as.vector(outer(l0,f0,FUN='+')),var=NULL),
                                                fix_g = c(TRUE,FALSE),tol=init_tol,conv_type = conv_type,
                                                maxiter = init_maxiter))
      }

      M = matrix(init_val$posterior$mean_log,nrow=n,ncol=p)
      sigma2_init = init_val$fitted_g$var
      rm(init_val)
    }
    if(var_type=='by_row'){
      if(verbose){
        cat('Solving VGA for row 1...')
      }
      init_val = mclapply(1:n,function(i){
        if(verbose){
          if(i%%printevery==0){
            cat(paste(i,'...'))
          }
        }
        if(ebpm_init){
          #fit = ebpm_exponential_mixture(Y[i,],s = drop(exp(l0[i]+f0)))
          #init_var_vga = mean(fit$posterior$mean_log^2)
          fit = ash_pois(Y[i,],scale = drop(exp(l0[i]+f0)),link='identity',mode=1,method='shrink')
          #init_var_vga = mean(log(fit$result$PosteriorMean)^2)
          init_var_vga = NULL
          fit = suppressWarnings(ebpm_normal(Y[i,],
                                             g_init = list(mean=drop(l0[i]+f0),var=init_var_vga),
                                             q_init = list(m_init=log(fit$result$PosteriorMean) + drop(l0[i]+f0),v_init = init_var_vga),
                                             fix_g = c(TRUE,FALSE),tol=init_tol,conv_type = conv_type,
                                             maxiter = init_maxiter))
        }else{
          fit = suppressWarnings(ebpm_normal(Y[i,],
                                             g_init = list(mean=drop(l0[i]+f0),var=NULL),
                                             fix_g = c(TRUE,FALSE),tol=init_tol,conv_type = conv_type,
                                             maxiter = init_maxiter))
        }

        return(list(sigma2 = fit$fitted_g$var,mean_log = fit$posterior$mean_log))
      },mc.cores = n_cores)
      sigma2_init = unlist(lapply(init_val,function(fit){fit$sigma2}))
      M = do.call(rbind,lapply(init_val,function(fit){fit$mean_log}))
      rm(init_val)
    }
    if(var_type=='by_col'){
      if(verbose){
        cat('Solving VGA for column 1...')
      }

      if(n_cores>1){
        init_val = mclapply(1:p,function(j){
          if(verbose){
            if(j%%printevery==0){
              cat(paste(j,'...'))
            }
          }
          if(ebpm_init){
            fit = ash_pois(Y[,j],scale = drop(exp(l0+f0[j])),link='identity',mode=1,method='shrink')
            #init_var_vga = mean(log(fit$result$PosteriorMean)^2)
            init_var_vga = NULL
            fit = suppressWarnings(ebpm_normal(Y[,j],
                                               g_init = list(mean=drop(l0+f0[j]),var=init_var_vga),
                                               q_init = list(m_init=log(fit$result$PosteriorMean) + drop(l0+f0[j]),v_init = init_var_vga),
                                               fix_g = c(TRUE,FALSE),tol=init_tol,conv_type = conv_type,
                                               maxiter = init_maxiter))
          }else{
            fit = suppressWarnings(ebpm_normal(Y[,j],
                                               g_init = list(mean=drop(l0+f0[j]),var=NULL),
                                               fix_g = c(TRUE,FALSE),tol=init_tol,conv_type = conv_type,
                                               maxiter = init_maxiter))
          }

          return(list(sigma2 = fit$fitted_g$var,mean_log = fit$posterior$mean_log))
        },mc.cores = n_cores)
        sigma2_init = unlist(lapply(init_val,function(fit){fit$sigma2}))
        M = do.call(cbind,lapply(init_val,function(fit){fit$mean_log}))
        rm(init_val)
      }else{
        M = matrix(nrow=n,ncol=p)
        sigma2_init = rep(0,p)
        for(j in 1:p){
          if(verbose){
            if(j%%printevery==0){
              cat(paste(j,'...'))
            }
          }
          if(ebpm_init){
            fit = ash_pois(Y[,j],scale = drop(exp(l0+f0[j])),link='identity',mode=1,method='shrink')
            #init_var_vga = mean(log(fit$result$PosteriorMean)^2)
            init_var_vga = NULL
            fit = suppressWarnings(ebpm_normal(Y[,j],
                                               g_init = list(mean=drop(l0+f0[j]),var=init_var_vga),
                                               q_init = list(m_init=log(fit$result$PosteriorMean),v_init = init_var_vga),
                                               fix_g = c(TRUE,FALSE),tol=init_tol,conv_type = conv_type,
                                               maxiter = init_maxiter))
          }else{
            fit = suppressWarnings(ebpm_normal(Y[,j],
                                               g_init = list(mean=drop(l0+f0[j]),var=NULL),
                                               fix_g = c(TRUE,FALSE),tol=init_tol,conv_type = conv_type,
                                               maxiter = init_maxiter))
          }
          M[,j] = fit$posterior$mean_log
          sigma2_init[j] = fit$fitted_g$var
        }
        rm(fit)
      }
    }
    if(log_init_for_non0y){
      idx0 = which(Y!=0)
      alpha0 = (tcrossprod(l0,rep(1,p))+tcrossprod(rep(1,n),f0))[idx0]
      M[idx0] = log(Y[idx0]/exp(alpha0)) + alpha0
    }
    return(list(sigma2_init=sigma2_init,M_init=M))
  }
  ## case 3: M_init is given, sigma2 is unknown, this case is useful when using glmpca for initialization.
  ## In this case we will run flash to init sigma2
  if(!is.null(M_init)&is.null(sigma2)){
    return(list(sigma2_init=NULL,M_init=M_init))
  }
  ## case 4: if L_init and F_init are given
  if(!is.null(L_init)&!is.null(F_init)){
    # do we need this case?
  }
  # if(!is.null(M_init)){
  #   if(verbose){
  #     cat('Initializing sigma2...')
  #   }
  #
  #   if(var_type=='constant'){
  #     if(verbose){
  #       cat('Solving VGA constant...For large matrix this may require large memory usage')
  #     }
  #     init_val = suppressWarnings(pois_mean_GG(as.vector(Y),as.vector(tcrossprod(l0,f0)),prior_mean = as.vector(M_init),prior_var = sigma2,tol=init_tol))
  #     sigma2_init = init_val$fitted_g$var
  #     rm(init_val)
  #   }
  #   if(var_type=='by_row'){
  #     if(verbose){
  #       cat('Solving VGA for row 1...')
  #     }
  #     init_val = mclapply(1:n,function(i){
  #       if(verbose){
  #         if(i%%printevery==0){
  #           cat(paste(i,'...'))
  #         }
  #       }
  #       fit = suppressWarnings(pois_mean_GG(Y[i,],l0[i]*f0,prior_mean = M_init[i,],prior_var = sigma2[i],tol=init_tol))
  #       return(fit$fitted_g$var)
  #     },mc.cores = n_cores)
  #     sigma2_init = unlist(init_val)
  #     rm(init_val)
  #   }
  #   if(var_type=='by_col'){
  #     if(verbose){
  #       cat('Solving VGA for column 1...')
  #     }
  #
  #     if(n_cores>1){
  #       init_val = mclapply(1:p,function(i){
  #         if(verbose){
  #           if(i%%printevery==0){
  #             cat(paste(i,'...'))
  #           }
  #         }
  #         fit = suppressWarnings(pois_mean_GG(Y[,i],l0*f0[i],prior_mean = M_init[,i],prior_var = sigma2[i],tol=init_tol))
  #         return(fit$fitted_g$var)
  #       },mc.cores = n_cores)
  #       sigma2_init = unlist(init_val)
  #       rm(init_val)
  #     }else{
  #       M = matrix(nrow=n,ncol=p)
  #       sigma2_init = rep(0,p)
  #       for(j in 1:p){
  #         if(verbose){
  #           if(j%%printevery==0){
  #             cat(paste(j,'...'))
  #           }
  #         }
  #         fit = suppressWarnings(pois_mean_GG(Y[,j],l0*f0[j],prior_mean = M_init[,j],prior_var = sigma2[j],tol=init_tol))
  #         sigma2_init[j] = fit$fitted_g$var
  #       }
  #       rm(fit)
  #     }
  #   }
  #   return(list(sigma2_init=sigma2_init,M_init=M_init))
  # }
}

