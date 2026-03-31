#'@title Summary of ebmpmf output
#'@export
#'

plot_stm = function(EL,EF){

}

# summary_stm = function(fit_obj,lambda,makeplot=T,typef='l'){
#   K = fit_obj$input$K
#   X = fit_obj$input$X
#   print(paste0('Mean KL:',round(mKL(lambda,fit_obj$Lambda_hat),5)))
#   print(paste0('Log likelihood:',round(fit_obj$loglik[length(fit_obj$loglik)]),5))
#   #print(paste0("Sum of EZ:",sum(fit_obj$EZ)))
#   if(makeplot){
#
#     par(mfrow=c(1,2))
#     image(lambda,main='True lambda')
#     image(fit_obj$Lambda_hat,main='Estimated lambda')
#     nr = floor(sqrt(K))
#     nc = ceiling(K/nr)
#
#     lf = poisson2multinom(t(fit_obj$qf$Ef),fit_obj$ql$El)
#
#     par(mfrow=c(nr,nc))
#     for(i in 1:K){
#       barplot(lf$L[,i],main=paste0('ebmpmf:estimate of loading',i),ylab = '')
#     }
#     for(i in 1:K){
#       plot(lf$FF[,i],main=paste0('ebmpmf:estimate of factor',i),ylab = '',type=typef)
#     }
#   }
#   par(mfrow=c(1,1))
# }
