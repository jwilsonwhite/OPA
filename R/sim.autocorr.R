#' Simulate autocorrelated landscape
#' 
#' Simulate a 1- or 2-dimensional landscape of autocorrelated values.
#' Uses FFT method based on Keitt (2000).
#'
#' @param n length of spatial domain (points will be spaced equally). For 1-D this is the length; for 2-D this is the dimension of one side of the domain.
#' @param mean the spatial mean value for the landscape parameter
#' @param sd the standard deviation of landscape values
#' @param range the desired range of the landscape variogram. That is, the decorrelation scale.
#' @param dimension either 1 (for 1-D landscape) or 2 (for 2-D)
#' 
#' @export
#' @return Vector or matrix of dimension \code{n} containing the landscape values for each point.
#' @examples 
#' sim.autocorr()
#' 
#' @seealso \code{\link{random.landscape}}
#' @references \href{http://link.springer.com/article/10.1023/A:1008193015770}{Keitt (2000)}
#' 
#' 
sim.autocorr <- function(n=10,mean=0,sd=1,range=1,dimension=1){

  try(if(range<=0) stop("range parameter must be >0",call.=FALSE))
  try(if(!any(dimension==c(1,2))) stop("dimension must be 1 or 2",call.=FALSE))
  try(if(length(mean)!=1) stop("mean must be length 1",call.=FALSE))
  try(if(length(sd)!=1) stop("sd must be length 1",call.=FALSE))
  
  # Determine dimension & calculate total # points
  if (dimension == 1){nn = n}
  else{nn = n^2}
      
  # rename variables for convenience
  Beta <- mean
  Sigma <- sd

# calculate distance matrix
  if (dimension == 1){
  Dist = 0:(n-1)    
  }else{ # 2D
  dist1 <- matrix(rep(0:(n-1),n),n,n)
  dist2 <- t(dist1)
  Dist = sqrt(dist1^2 + dist2^2)}

# spherical variogram
Vgm = 1.5*Dist/range - 0.5*(Dist/range)^3
Vgm[Vgm>1]=1
Vgm[Vgm<0]=0

if (range > 0){
Vgm[Dist > range] = 1}

# scaled covariogram is 1-variogram
Covgm = 1 - Vgm;

# perform fast Fourier transform on covariogram
Z <- fft(Covgm);

# Theta will randomize phase (imaginary part) of FFT spectrum
Theta <- runif(nn)*2*pi;

# Now do inverse FFT with real part of covariogram, but with randomized phase
Z <- Re(Z)*exp(1i*Theta);
Z <- fft(Z,inverse=TRUE)

# Now take modulus of inverse FFT result (and rescale by number of points)
Z2 <- Mod(Z)/(nn);
# Bolker suggests it would be best to confirm that transform is
# Hermitian (i.e. Z2 == conj(Z2') ?
                          
                            
# Normalize autocorrelated noise
if (sd(Z2,na.rm=TRUE) > 0){
  Z2 <- (Z2 - mean(Z2))/sd(Z2,na.rm=TRUE)
}else{
  Z2 <- Z2 - mean(Z2)}

# Apply Gaussian noise
D <- Z2 + rnorm(nn,0,1)

# Rescale to have mean Beta and SD Sigma
D = D/sd(D)*Sigma + Beta
} # end function

                            