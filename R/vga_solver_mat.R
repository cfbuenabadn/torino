

#'@title a matrix version of the vga solver using Newton's method
#'@importFrom vebpm vga_pois_solver_bisection
#'@importFrom Rfast Pmin
vga_pois_solver_mat_newton = function(M,X,S,Beta,Sigma2,maxiter=1000,tol=1e-8,return_V = TRUE){

  #X = as.matrix(X)
  const0 = (Sigma2*X+Beta + 1)
  const1 = 1/Sigma2
  const2 = Sigma2/2
  const3 = Beta/Sigma2


  # make sure m < sigma2*x+beta
  M = Pmin(M,const0-1)
  # idx = (M>(const0-1))
  # if(sum(idx)>0){
  #   M[idx] =suppressWarnings(vga_pois_solver_bisection(c(X[idx]),c(S[idx]),c(Beta[idx]),c(Sigma2[idx]),maxiter = 10)$m)
  # }
  for(i in 1:maxiter){
    temp = (const0-M)
    sexp = S*exp(M+const2/temp)
    # f = X - sexp - (M-Beta)/Sigma2
    f = X - sexp - M*const1 + const3
    if(max(abs(f))<tol){
      break
    }
    # f_grad = -sexp*(1+const2/temp^2)-const1
    # direction = (X - sexp - (M-Beta)/Sigma2)/(-sexp*(1+const2/temp^2)-const1)
    M = M - f/(-sexp*(1+const2/temp^2)-const1)

  }
  #gc()
  if(return_V){
    return(list(M=M,V=Sigma2/temp))
  }else{
    return(M)
  }

}

#'@title vga solver for fixed iteration, not necessary to convergence
vga_pois_solver_mat_newton_fixed_iter = function(M,V,Y,S,Beta,Sigma2,maxiter=1000){

  #M = pmin(M,Sigma2*Y+Beta)
  if(is.null(V)){
    V = Sigma2/2
  }
  #V = pmin(V,Sigma2)

  for(i in 1:maxiter){
    # newton for M
    temp = S*exp(M+V/2)
    # M = M - (Y-temp-(M-Beta)/Sigma2)/(-temp-1/Sigma2)
    M = M - (Y/temp-1-(M-Beta)/Sigma2/temp)/(-1-1/Sigma2/temp)
    # newton for V
    temp = S*exp(M+V/2)*V/2
    #num_V = -temp-V/2/Sigma2 + 0.5
    #denom_V = -temp*(1 + V/2) - V/2/Sigma2
    #V = V/exp((-temp-V/2/Sigma2 + 0.5)/(-temp*(1 + V/2) - V/2/Sigma2))
    V = V/exp((-1-V/2/Sigma2/temp + 0.5/temp)/(-(1 + V/2) - V/2/Sigma2/temp))
  }
  return(list(M=M,V=V))
}


#'@title a matrix version of the vga solver using Newton's method, low memory mode by partitioning the rows
#'@importFrom vebpm vga_pois_solver_bisection
#'@importFrom Rfast Pmin
vga_pois_solver_mat_newton_low_memory = function(M,Y,l0,f0,EL,EF,
                                                 sigma2,var_type,
                                                 maxiter=1000,tol=1e-8){

  n = nrow(M)
  p = ncol(M)
  Beta = tcrossprod(EL,EF)
  M = Pmin(M,Beta)

  if(var_type=='by_col'){
    const0 = Y%*%Diagonal(p,sigma2)+Beta+ 1
    for(i in 1:maxiter){
      sexp = l0*exp(M+(1/(const0-M))%*%Diagonal(p,sigma2/2))%*%Diagonal(p,f0)
      f = Y-sexp + (Beta-M)%*%Diagonal(p,1/sigma2)
      if(max(abs(f))<tol){
        break
      }
      # f_grad = - sexp*(1+0.5*(const0-M)^(-2)) - matrix(1/sigma2,nrow=n,ncol=p,byrow = TRUE)
      # direction = (X - sexp - (M-Beta)/Sigma2)/(-sexp*(1+const2/temp^2)-const1)
      M = M - f/(- sexp*(1+0.5*(const0-M)^(-2)) - matrix(1/sigma2,nrow=n,ncol=p,byrow = TRUE))
    }
  }

  if(var_type=='by_row' | var_type=='constant'){
    const0 = sigma2*Y+Beta+ 1
    for(i in 1:maxiter){
      sexp = l0*exp(M+(1/(const0-M))*sigma2/2)%*%Diagonal(p,f0)
      f = Y-sexp + (Beta-M)/sigma2
      print(range(f))
      if(max(abs(f))<tol){
        break
      }
      # f_grad = - sexp*(1+0.5*(const0-M)^(-2)) - matrix(1/sigma2,nrow=n,ncol=p,byrow = TRUE)
      # direction = (X - sexp - (M-Beta)/Sigma2)/(-sexp*(1+const2/temp^2)-const1)
      M = M - f/(- sexp*(1+0.5*(const0-M)^(-2)) - 1/sigma2)
    }
  }
  return(as.matrix(M))
}




#' #'@title a matrix version of the vga solver using Newton's method. log version Not optimal because log(const0-M-1) can be NaN
#' #'@importFrom vebpm vga_pois_solver_bisection
#' vga_pois_solver_mat_newton = function(M,X,S,Beta,Sigma2,maxiter=1,tol=1e-8){
#'
#'   const0 = Sigma2*X+Beta + 1
#'   const1 = log(S*Sigma2)
#'
#'   # make sure m < sigma2*x+beta
#'   idx = (M>(const0-1))
#'   if(sum(idx)>0){
#'     M[idx] =suppressWarnings(vga_pois_solver_bisection(c(X[idx]),c(S[idx]),c(Beta[idx]),c(Sigma2[idx]),maxiter = 10)$m)
#'   }
#'   for(i in 1:maxiter){
#'     f = M + Sigma2/(const0-M)/2 - log(const0-M-1) + const1
#'     if(max(abs(f))<tol){
#'       break
#'     }
#'     M = M - f/(1 + Sigma2/(const0-M)^2/2 + 1/(const0-M-1))
#'   }
#'
#'   return(list(M=M,V=Sigma2/(const0-M)))
#' }



#' #'@title a matrix version of the solver. Transform to vector version.
#' #'@importFrom vebpm vga_pois_solver
#' vga_pois_solver_mat = function(init_Val,X,S,Beta,Sigma2,maxiter=1000,tol=1e-8){
#'
#'   n = nrow(X)
#'   p = ncol(X)
#'   # transform them to vector
#'   x = as.vector(X)
#'   s = as.vector(S)
#'   beta = as.vector(Beta)
#'   sigma2 = as.vector(Sigma2)
#'   init_val = as.vector(init_Val)
#'
#'   res = vga_pois_solver(init_val,x,s,beta,sigma2,maxiter=maxiter,tol=tol,method='newton')
#'
#'   return(list(M = matrix(res$m,nrow=n,ncol=p,byrow = F),V = matrix(res$v,nrow=n,ncol=p,byrow = F)))
#'
#' }
#'
#'

# vga_pois_solver_mat_newton_fixed_iter_debug = function(M,V,Y,S,Beta,Sigma2,KL_LF,const,sigma2,fit_flash,var_type,maxiter=1000,tol=1e-8){
#
#   f_mv = function(v,m,y,s,theta,sigma2){
#     y*m - s*exp(m + v/2) - (m^2 + v -2*m*theta)/(2*sigma2) + log(v)/2
#   }
#
#   if(!is.null(V)){
#     print(paste('init0 vga obj vga=',sum(f_mv(V,M,Y,S,Beta,Sigma2))))
#     print(paste('init0,elbo is',round(calc_split_PMF_obj_flashier(Y,S,sigma2,M,V,fit_flash,KL_LF,const,var_type),3)))
#   }
#
#   M = pmin(M,Sigma2*Y+Beta)
#   if(is.null(V)){
#     V = Sigma2/(Sigma2*Y-M+Beta+1)
#   }
#   V = pmin(V,Sigma2)
#
#
#   print(paste('init obj vga=',sum(f_mv(V,M,Y,S,Beta,Sigma2))))
#   print(paste('init,elbo is',round(calc_split_PMF_obj_flashier(Y,S,sigma2,M,V,fit_flash,KL_LF,const,var_type),3)))
#
#   for(i in 1:maxiter){
#     # newton for M
#     M = M - (Y-S*exp(M+V/2)-(M-Beta)/Sigma2)/(-S*exp(M+V/2)-1/Sigma2)
#     print(paste('M obj vga=',sum(f_mv(V,M,Y,S,Beta,Sigma2))))
#     print(paste('M,elbo is',round(calc_split_PMF_obj_flashier(Y,S,sigma2,M,V,fit_flash,KL_LF,const,var_type),3)))
#     # newton for V
#     num_V = -S*exp(M+V/2)*V/2-V/2/Sigma2 + 0.5
#     denom_V = -V/2*S*exp(M+V/2)*(1 + V/2) - V/2/Sigma2
#     #V = exp(log(V)-(num_V)/(denom_V))
#     V = V/exp(num_V/denom_V)
#     print(paste('V obj vga=',sum(f_mv(V,M,Y,S,Beta,Sigma2))))
#     print(paste('V,elbo is',round(calc_split_PMF_obj_flashier(Y,S,sigma2,M,V,fit_flash,KL_LF,const,var_type),3)))
#   }
#   print(range(V))
#
#   return(list(M=M,V=V))
# }
