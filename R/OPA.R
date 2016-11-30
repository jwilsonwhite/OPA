#' The Object Pattern Analysis Package
#' 
#' \pkg{OPA} allows you to perform object pattern analysis. This is an extension of the Ripley's \eqn{K} summary statistic for point pattern analysis, but
#' for objects that have substantial two-dimensional area that cannot be appropriately approximated as points.
#' 
#' The function you will primarily use in \pkg{OPA} is \code{\link{K.obj}}, which performs the analysis. Refer to vignettes for examples 
#' of how to format data and interpret results. \pkg{OPA} requires \pkg{spatstat}; data objects in both packages are of the class \code{ppp} and 
#' function such as \code{\link[spatstat]{plot.ppp}}, \code{\link[spatstat]{plot.fv}}  and \code{\link[spatstat]{envelope}} can be used to process
#' data and output of \code{\link{K.obj}}.
#' 
#' The other useful function is \code{\link{random.landscape}}. This function can be used to generate randomly generated object-pattern
#' landscapes to use for deviations from randomness in an observed object-pattern landscape (see vignettes). It can also be used independently to create
#' point-pattern or object-pattern landscapes that are overdispersed, underdispersed, or spatially random.
#' 
#' @name OPA
#' @import spatstat
#' @importFrom stats dist fft quantile rnorm runif
#' 
#' @examples 
#' \donttest{
#' X = random.landscape(maxR = 100,n=50,range.pts = 10,range.radii = 1,radius.mean=1,radius.sd=1) 
#' # A 100 m x 100 m landscape of objects that are spatially clustered at a range of 10 m
#' K = K.obj(X,Type='rect')
#' # plot the K function versus distance (r)
#' plot(K) 
#' # plot the scaled L statistic instead (expectation = 1 under complete spatial randomness)
#' plot(K,L~r) 
#' 
#' # Use spatstat::envelope to calculate a significance band around the K statistic 
#' # (note that this is computationally intensive and may take a long time!):
#' KE = envelope(X,fun=K.obj,funargs=list(Type='rect'),
#' nsim=1e3, simulate=expression(random.landscape(100,50,radius.mean=1,radius.sd=1)))
#' 
#' # Use spatstat::envelope to calculate a significance band around the L statistic, 
#' # which is a scaled transformation of the L statistic. 
#' # To perform this scaling you must directly specify the \code{Eps} parameter:
#' Eps = 0.01
#' KE = envelope(X,fun=K.obj,funargs=list(Type='rect',Eps=Eps),
#' nsim=1e3, simulate=expression(random.landscape(100,50,radius.mean=1,radius.sd=1)),
#' transform = expression(./(2*pi*r*Eps)))
#' plot(KE,legendmath=FALSE)
#' 
#' # Perform the same analysis on the sponge dataset:
#' data(spongemap)
#' units = abs(diff(spongemap$window$xrange)) # diameter of the landscape
#' n = length(spongemap$marks) # number of objects
#' rm = mean(spongemap$marks) # mean object radius
#' rs = sd(spongemap$marks) # sd of object radius
#' Eps = 0.01
#' SE = envelope(spongemap,fun=K.obj,funargs=list(Type='circle',Eps = Eps),
#' nsim=1e3, simulate=expression(random.landscape(units=units,n=n,radius.mean=rm,
#' radius.sd=rs,Type='circle')),transform = expression(./(2*pi*r*Eps)))
#' plot(SE,legendmath=FALSE)
#'  
#' }
"_PACKAGE"


