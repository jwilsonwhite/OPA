#' No overlap
#' 
#' Filters an object-pattern landscape to eliminate overlapping objects
#' 
#' @param Input vector or matrix of binary values (1 and 0) representing the centroids of regularly-spaced object pattern.
#' @param R vector or matrix of the same dimension as Input, containing the radii of each object.
#' @param Oa an optional value specifying the proportional overlap allowed. Must be between 0 (default) and 1, inclusive.
#' 
#' @export
#' @return Vector or matrix of the same dimensions and number of points as Input, but with overlapping objects removed.
#' @examples 
#' X = matrix(round(runif(10)),nrow=5,ncol=2) ## create a simple object pattern
#' R = matrix(rnorm(10),nrow=5,ncol=2) ## Associated radii
#' no.overlap(Input=X,R=R,Oa=0)
#' 
#' @seealso \code{\link{K.obj}},\code{\link{random.landscape}}
#' 

no.overlap <- function(Input,R,Oa=0){
  
  # error checking:
  try(if( Oa < 0 | Oa > 1) stop("Oa cannot be < 0 or > 1",call.=FALSE))
Oa <- -Oa # make negative for later use
  
# Distance matrix for each point:
N <- sum(Input) # total number of points in Input landscape
Y <- matrix(rep(1:dim(Input)[2],times=dim(Input)[1]),nrow=dim(Input)[1],ncol=dim(Input)[2])
X <- t(Y)
Q <- dist( cbind(matrix(X,nrow=prod(dim(X))), matrix(Y,nrow=prod(dim(Y)))),diag=TRUE,upper=TRUE)
Q <- as.matrix(Q)
Q[Q<1e-10] <- NA # eliminate 'self' distances

Index = 1:dim(Q)[1] # point indices

# Get pairwise distances including radii
R2 <- matrix(rep(R,times=dim(Q)[1]),nrow=dim(Q)[1],ncol=dim(Q)[2])
Rmat <- R2 + t(R2) # all pairwise radius sums

QR <- (Q-Rmat)/Q; # proportional difference between distance & summed radii.  If negative, points overlap

# Extract the actual points from QR
isP <- Input > 0;
OKP <- Index[isP]
QRtmp <- QR[isP,isP];
OK <- !any(QRtmp<Oa,na.rm=TRUE) # apply overlap criterion

while (!OK){

# Find the point with the most failures 
Tmp <- apply(QRtmp<Oa,1,sum,na.rm=TRUE)
BadP <- OKP[Tmp == max(Tmp)]
BadP = BadP[1] # in case there is a tie
  
isP[BadP]=FALSE # remove that point
OKP <- Index[isP] # remove from list of points

QRtmp <- QR[isP,isP]
OK <- !any(QRtmp<Oa,na.rm=TRUE) # apply overlap criterion

if (length(OKP) == 1){ # in case we are down to one point
OK <- TRUE} # end if
} # end while

Out <- Input*0
Out[OKP] <- 1
Out

} # end function