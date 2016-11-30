#' Discretize
#'
#' Converts a landscape of continuous-valued points into a discrete set of \emph{N} points
#' 
#' @param Input vector or matrix of numeric values representing some continuously-valued landscape parameter
#' @param N The number of points to be created on the landscape
#' 
#' @export
#' @return Vector or matrix of the same dimensions as \code{Input}, but consisting of \emph{N} 1s and the remainder 0s
#' @examples 
#' X <- matrix(rnorm(10),nrow=10,ncol=1)
#' discretize(Input=X,N=5)
#' 
#' @seealso \code{\link{K.obj}},\code{\link{random.landscape}}
#' 

discretize <- function(Input,N){
  D = dim(Input)
  H <- sort(Input,decreasing=TRUE,na.last=TRUE)
  Hmin <- H[N]
  Output <- Input >= Hmin
  Output = as.numeric(Output)
  Output = matrix(Output,nrow=D[1],ncol=D[2])
}