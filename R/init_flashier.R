#'@title init flash add greedy
#'@description this init function extracts previous init values so no need to perform rank1 svd again.
init.fn.fix = function(f){f$init_vals}

#' #'@title init flash add greedy for dense matrix
#' #'@description use irlba for dense matrix, faster.
#' init.fn.dense = function(f,dim.signs){
#'   res = init.fn.irlba(f)
#'   if(all(dim.signs==0) | is.null(dim.signs)){
#'     return(res)
#'   }else{
#'     if(dim.signs[1]==1){
#'       res[[1]] = pmax(res[[1]],0)
#'     }
#'     if(dim.signs[1]==-1){
#'       res[[1]] = pmin(res[[1]],0)
#'     }
#'     if(dim.signs[2]==1){
#'       res[[2]] = pmax(res[[2]],0)
#'     }
#'     if(dim.signs[2]==-1){
#'       res[[2]] = pmin(res[[2]],0)
#'     }
#'     return(res)
#'   }
#' }
