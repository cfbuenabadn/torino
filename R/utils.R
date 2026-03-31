#'@title Make Init L and F on similar scale
#'@param L N by K loading matrix
#'@param FF p by K factor matrix
#'@export
adjLF = function(L,FF){
  gammaL = colSums(L)
  gammaF = colSums(FF)
  adjScale = sqrt(gammaL*gammaF)
  L = t(t(L) * (adjScale/gammaL))
  FF = t(t(FF) * (adjScale/gammaF))
  return(list(L_init = L, F_init=FF))
}

#'@title Log transformation of scRNA-seq count matrix for EBMF
#'@export
log_for_ebmf = function(Y){
  log(1+median(rowSums(Y))/0.5*Y/rowSums(Y))
}


#'@title calculate mean KL divergence of 2 nonnegative matrices
#'@export
#'

mKL = function(A,B){
  D = A*log(A/B)-A+B
  mean(as.matrix(D),na.rm=T)
}


#'@title standard loadings and factors from Poisson matrix Factorization
#'@param L: n by k matrix
#'@param FF: k by p matrix
#'@export
#'

poisson_to_multinom <- function (FF, L) {
  L <- t(t(L) * colSums(FF))
  s <- rowSums(L)
  L <- L / s
  FF <- scale_cols(FF)
  return(list(FF = FF,L = L,s = s))
}

poisson_to_libsize <- function (FF, L, lib_size) {
  res = poisson_to_multinom(FF,L)
  size = res$s/lib_size
  multinom_to_poisson(res$FF,res$L,size)
}

multinom_to_poisson <- function (FF, L,size) {
  L = L * size
  res = adjLF(L,FF)
  return(list(FF = res$F_init,L = res$L_init))
}

scale_cols <- function (A)
  apply(A,2,function (x) x/sum(x))


calc_EZ = function(x, prob){
  Ez = sparseMatrix(i = x$i, j = x$j, x = x$x * prob)
  return(list(rs = Matrix::rowSums(Ez), cs = Matrix::colSums(Ez)))
}

softmax3d=function(x){
  #x = x - array(apply(x,c(1,2),max),dim=dim(x))
  x = exp(x)
  p=as.vector(x)/as.vector(rowSums(x,dims=2))
  p = pmax(p,1e-10)
  dim(p) <- dim(x)
  return(p)
}

rowMax = function(X){
  do.call(pmax.int, c(na.rm = TRUE, as.data.frame(X)))
}
