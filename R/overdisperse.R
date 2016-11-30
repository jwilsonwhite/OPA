#' Overdisperse
#' 
#' Cause a landscape of \emph{N} points to become overdispersed, with a minimum distance between them
#'
#' @param Input vector or matrix of binary values (1 and 0) representing a regularly spaced point pattern
#' @param Mindist the minimum distance allowed between points
#' 
#' @export
#' @return A vector or matrix of the same dimensions and number of points as \code{Input}, but rearranged so that no points are closer than \code{Mindist}
#' @examples 
#' X = round(runif(20))
#' overdisperse(Input=X,Mindist=2)
#' 
#' @seealso \code{\link{K.obj}},\code{\link{random.landscape}},\code{\link{discretize}}
#' 


overdisperse <- function(Input,Mindist=0){
  N <- sum(Input) # total number of points in Input landscape
  Y <- matrix(rep(1:dim(Input)[2],times=dim(Input)[1]),nrow=dim(Input)[1],ncol=dim(Input)[2])
  X <- t(Y)
  # Distance matrix:
  Q <- dist( cbind(matrix(X,nrow=prod(dim(X))), matrix(Y,nrow=prod(dim(Y)))),diag=TRUE,upper=TRUE)
  Q = as.matrix(Q)
  Q[Q<1e-10] <- NA # eliminate 'self' distances

  Index <- 1:dim(Q)[1] # Point indices
  OKrows <- Index[Input>0] # Values that are points
  Other <- Index[Input==0] # # indices that are not points
  Qtmp <- Q[OKrows,OKrows]
  OK <- min(Qtmp >= Mindist,na.rm=TRUE) # Are all of the distances larger than Mindist?

  x <- 1
  while (!OK){
    # find column of Qtmp with smallest mean distance
    Minvec <- apply(Qtmp,1,min,na.rm=TRUE) 
    Smallrows <- match(min(Minvec),Minvec) 
    OKtmp <- rep(TRUE,length(OKrows))
    OKtmp[Smallrows] = FALSE
    Remove = OKrows[!OKtmp] # which one removed?
    OKrows = OKrows[OKtmp]
    Try = sample(Other,1) # which new point to try
    Other = Other[Other!=Try] # take it out of the list
    OKrows[length(OKrows)+1] <- Try # Add additional, randomly sampled new point to replace the one removed
    Other[length(Other)+1]=Remove # add removed point to the list of unused points to choose
    
    Qtmp <- Q[OKrows,OKrows] # reset the points in the temp distance matrix
    x <- x + 1
    OK <- min(Qtmp >= Mindist,na.rm=TRUE)
    if (x > 2*dim(Q)[1]){ 
      OK <- TRUE} # Stopper if we have been at this for a while
  } # end while
  
  # format output matrix
  Out <- rep(0,dim(Q)[1])
  Out[OKrows] = 1
  Out = matrix(Out,dim(Input)[1],dim(Input)[2])
  Out
} # end function
  
