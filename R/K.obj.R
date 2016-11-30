#' Object pattern analysis
#' 
#' Calculate \eqn{K_OR}: object-pattern equivalents to the point-pattern Ripley's \eqn{K} summary statistic
#' 
#' @author Will White, \email{whitejw@@uncw.edu}
#' 
#'@param X an object of class \code{"ppp"} containing \code{x} and \code{y} coordinates & object radii in the \code{marks} field.
#'@param Type character string indicating the shape of the landscape. Must be either \code{'rect'} or \code{'circle'}.
#'@param r distances at which to calculate \eqn{K}. The default will choose reasonable values.
#'@param Eps the width of the ring used in the ring pattern calculation. It works best if it is quite small relative to the radii of the data and the sample radii \code{r}. It will be automatically calculated if not provided.
#'       Specifying \code{Eps = 0} will coerce the code to perform regular point-pattern Ripley's \eqn{K} calculation, ignoring the object radii.
#'@param which.marks optional argument; if the \code{ppp} object has a \code{marks} field with multiple columns, this indicates which column contains the object radii
#'
#'@export
#'
#'@return An object of class \code{"fv"}, similar to that created by \code{\link[spatstat]{Kest}}. 
#'        It is essentially a data frame with columns
#'        \itemize{
#'        \item \code{r}. the vector of values of the argument \code{r} at which \eqn{K_OR} has been estimated
#'        \item \code{K}. the vector of \eqn{K} or \eqn{K_OR} values (depending on the options selected) at each value of \code{r}
#'        \item \code{theo} the theoretical expected value of \eqn{K} ot \eqn{K_OR} under complete spatial randomnes
#'        \item \code{L} the value of \eqn{K} ot \eqn{K_OR} rescaled by \code{theo} to have expectation of 1 under complete spatial randomness.
#'        }
#'        
#'        Values of \code{L} > 1 indicate clustering; values < 1 indicate regularity. Note that the value of \code{L} corresonding to the smallest 
#'        value of \code{r} is often >> 1 just because it lies interior to some of the data objects.

#' @details This function uses input and output objects in a way similar to package \code{spatstat}, and is specifically analogous to 
#'          \code{\link[spatstat]{Kest}}. To compute simulation envelopes for the \eqn{K} or \eqn{L} function, use \code{\link[spatstat]{envelope}} in the spatstat package.
#'          If \code{Eps = 0}, the function will compute the standard point-pattern analysis (PPA) Ripley's \eqn{K}. This is the same calculation as in 
#'          \code{\link[spatstat]{Kest}}, although that function has a wider range of edge-correction options and is much faster for PPA.
#'          If all of the radii in \code{X} are equal to zero, this function will also calculate the PPA version of \eqn{K}, but if \code{Eps > 0} or
#'          is \code{NULL}, the 'ring function' version of \eqn{K} will be calculated, as suggested by Wiegand & Molony (2004).
#'          
#' @examples 
#' X = random.landscape()
#' K.obj(X,Type='rect',Eps=0.1)
#' 
#' @seealso \code{\link{sim.autocorr}},  \code{\link{random.landscape}}, \code{\link[spatstat]{Kest}}, \code{\link[spatstat]{ppp.object}}
#' 
#' @references 
#' \href{http://onlinelibrary.wiley.com/doi/10.1111/j.0030-1299.2004.12497.x/full}{Wiegand and Moloney 2004}, Diegnan LK, Pawlik JR, Gleason ACR, White JW. "Object Pattern Analysis: an extension of Ripley's \eqn{K} to two-dimensional objects." Submitted, Methods in Ecology & Evolution 



K.obj <- function(X,Type,Eps=NULL,r=NULL,which.marks=1){
  
  try( if(!any(Type==c('circle','rect'))) stop("Type must be either 'circle' or 'rectangle'",call.=FALSE))
  
 # Error checks:
  try( if(class(X) != "ppp") stop("X must be an object of class 'ppp'",call.=FALSE))
  
  try( if(!any(class(X$marks) == c("data.frame","numeric"))) stop("marks in X must either be a numeric vector or a data frame",call.=FALSE))
  
  try( if(!any(X$window$type == c("rectangle","polygonal"))) stop("window must be rectangular or circular; polygons and masks not supported",call.=FALSE)) 
  # When owin() creates a circular window (using disc() ) it really just creates a high-order polygonal approximation. This means that the associated type
  # stored in window$type is 'polygonal'. However this code is not designed to work with polygons, so we have to assume the user has really created a circle here.
  
  # Extract the radii
  if (class(X$marks)=="numeric"){ # vector
    Rad = X$marks
  }else{ # data frame
    Rad = X$marks[,which.marks]}
  try(if(class(Rad)!="numeric") stop("radii must be numeric values",call.=FALSE))
  
  # Extract x,y
  x = X$x
  y = X$y
  
  #If the data are all points, calculate point-pattern Ripley's K
  isPPA = all(Rad==0) # if all radii are 0, do PPA
  
  # Can only have Eps = 0 if doing PPA (standard area-based PPA, not ring calculation)
 # if (!is.null(Eps)){
#  try( if(!isRad & (Eps == 0)) stop("Eps must be nonzero if calculating object pattern analysis (i.e., if one or more object radii is nonzero)",call.=FALSE))
#  }
  
  # sample size
  n = length(Rad)
  
  # Center the data
  MuX <- mean(X$window$xrange) # center x coordinate
  MuY <- mean(X$window$yrange) # center y coordinate
  x = x - MuX
  y = y - MuY
  
  # Operations on the window - determine radius, etc
  if (Type == 'circle'){
    # ppp objects have a window stored as a polygon (default size = 128 vertices)
    # Back out to obtain radius & center by taking half the diameter
    maxR <- diff( X$window$xrange)/2} # note this is a scalar value

  if (Type == 'rect'){
  maxR <- min(c( diff( X$window$xrange)/2, diff( X$window$yrange)/2 ))  # min dimension ('radius') in either direction
  
  # Vertex points of the rectangle. X$window gives the southwest corner
  Vertices <- cbind( c(X$window$xrange[1], X$window$xrange[2],X$window$xrange[2],X$window$xrange[1] ),
                     c(X$window$yrange[1], X$window$xrange[1],X$window$xrange[2],X$window$xrange[2] ))
  # Center the vertices
  Vertices[,1]=Vertices[,1]-MuX
  Vertices[,2]=Vertices[,2]-MuY
  } # end if circle or rect
  
  # Define scales of influence
  if (is.null(r)){
  r <- seq(from=quantile(Rad,0.75),to=maxR/2,length.out=512)} # Ripley's rule of thumb, as in Kest
  
  # Width of ring 
  if (is.null(Eps)){
    Eps <- 1e-3*min(r)}
  
  # Calculate mu, the proportional area occupied by the objects in X
  if (Type == 'circle'){
  # circle areas:
    Ca <- sum(Rad^2) # area of each object
    totalA <- pi*maxR^2 # total area of the study plot (used later for PPA if needed)
    Mu <- Ca/(maxR^2) # proportion occupied (the pi terms cancel out)
  }else{ # rectangle
    Ca <- pi*sum(Rad^2) # area of all objects
    totalA <- (2*maxR)^2
    Mu <- Ca/totalA}
    
  # pre-allocate some variables
  A <- array(NA,dim=c(n,n,length(r)))
  Ai <- array(NA,dim=c(n,length(r)))
  wi <- rep(NA,n)
  wo <- rep(NA,n)
  
  
# Loops to calculate K.obj
for (rr in 1:length(r)){ # loop over each scale of clustering
  for (i in 1:n){ # loop over each data point

    Disc.i <- c(x[i],y[i],Rad[i]) # define the focal circle
    
# Calculate weighting factors wo & wi
# (proportion of the study area that is within the circle of influence)
   
   if (Type == 'circle'){
   Disc.SA <- c(0,0,maxR) # disc defining the study area
   # If circle lies entirely within study area
   OL1 <- is.overlap(Disc.i,Disc.SA,r[rr]+Eps/2);
   OL2 <- is.overlap(Disc.i,Disc.SA,r[rr]-Eps/2);
   
   if (OL1){
   wo[i] = 1}
   else{
   wo[i] = Afxn(Disc.i,Disc.SA,r[rr]+Eps/2)/(pi*(r[rr]+Eps/2)^2)}
   
   if (OL2){
   wi[i] = 1}
   else{
   wi[i] = Afxn(Disc.i,Disc.SA,r[rr]-Eps/2)/(pi*(r[rr]-Eps/2)^2)}

   }else{ # rectangle:
   # If circle lies entirely within rectangular study area
   OL1 = is.overlap.rect(c(x[i],y[i],Rad[i]),Vertices,r[rr]+Eps/2)
   OL2 = is.overlap.rect(c(x[i],y[i],Rad[i]),Vertices,r[rr]-Eps/2) 
   
   if (OL1){
   wo[i] = 1;
   }else{
   wo[i] = Afxn.rect(Disc.i,Vertices,r[rr]+Eps/2)/(pi*(r[rr]+Eps/2)^2)}
   
   if (OL2){
   wi[i] = 1
   }else{
   wi[i] = Afxn.rect(Disc.i,Vertices,r[rr]-Eps/2)/(pi*(r[rr]-Eps/2)^2)}
  } # end if/else for Type
  # End calculation of wo & wi

  # Loop over each potential neighbor point:
  for (j in 1:n){
  if (j == i){
  A[i,j,rr] = 0} # if it's the same point
  else{
    d <- ( (x[i] - x[j])^2 + (y[i] - y[j])^2)^0.5 # Euclidean distance
    rj <- Rad[j]
    Disc.j <- c(x[j],y[j],rj) # disc for object j
   
   # If just doing regular Ripley's K:
   if (isPPA){
     
    # If doing PPA Ripley's K, can either use standard approach (count all points within Rad[r]) or use ring method:
    # Standard (non-ring) method:
     if (Eps == 0){

        if (r[rr] > d){
        A[i,j,rr] = 1/wo[i]} # if the two points are within the radius.  Scale by proportional overlap, wo (ignore wi in this case)
        else{
        A[i,j,rr] = 0}
     }else{
       # Ring method
        if (r[rr]+Eps/2 > d & r[rr]-Eps/2 <= d){ # if inside the ring
        A[i,j,rr] = 1/mean(c(wi[i],wo[i])) # scale by average proportional overlap of the ring 
        }else{
        A[i,j,rr]=0}
       }
     # end standard Ripley's PPA
   
   }else{ # otherwise do Ripley's Kc:
     
     # If i is totally inside j:
     if (rj >= d + r[rr] + Eps/2){ 
     A[i,j,rr] = pi*( (r[rr]+Eps/2)^2 - (r[rr]-Eps/2)^2 )} # total area of object i
     
     if ( (rj < d + r[rr] + Eps/2) & (rj > d + r[rr] - Eps/2)){ # only the inner one is totally enveloped
     Outer =  Afxn(Disc.i,Disc.j,r[rr]+Eps/2)
     Inner = pi*(r[rr]-Eps/2)^2 # total area of i
     A[i,j,rr] = Outer/wo[i] - Inner/wi[i]}
     
     if (rj < d + r[rr] - Eps/2){ # some opportunity for overlap with both rings
      if (r[rr] + Eps/2 > d - rj){ # if they actually intersect (otherwise don't calculate)
        if (r[rr] + Eps/2 > rj + d){ # complete envelopment
        Outer = pi*rj^2}
        else{
        Outer = Afxn(Disc.i,Disc.j,r[rr]+Eps/2)}
        if (r[rr] - Eps/2 > d - rj){ # if the inner ring also intersect   
        if (r[rr] - Eps/2 > rj + d){ # complete envelopment
        Inner = pi*rj^2}
        else{
        Inner = Afxn(Disc.i,Disc.j,r[rr]-Eps/2)}
        }else{ # if inner ring does not intersect
        Inner = 0}
   
      A[i,j,rr] = Outer/wo[i] - Inner/wi[i]
      }else{ # else if they don't intersect
      A[i,j,rr] = 0}
     } # end if there is some opportunity to overlap
   
   }# end if isRad
   }# end if j == i
   }# end j loop
   
   Ai[i,rr] = sum(A[i,,rr],na.rm=TRUE) # sum of all j overlaps
   
   
   } # end i loop
   } # end r loop
   
   if (!isPPA){
   # normalized ring object function, as in Deignan et al. Eq. 7
   K = colSums(Ai)/n/Mu
   Theo = 2*pi*r*Eps
   L = K/Theo
   }else{
   if (Eps == 0){   
  #   browser()
   K = colSums(Ai)/(n)/(n/totalA) # K_or(r)
   Theo = (pi*r^2) # the theoretical value expected under CSR
   L = K/Theo # K rescaled to L_or(r) so that CSR has a value of 1
   }else{
   K = colSums(Ai)/n/(n/totalA) 
   Theo = 2*pi*r*Eps
   L = K/Theo}
   } # end if isRad
   

# Align output with that of Kest in spatstat
Out <- data.frame(r=r,theo=Theo,K=K,L=L)
Out <- fv(Out,argu='r',valu='K')
Out
} # end function K.obj

# Helper functions: 

#----------------------------------------------------------------
# Afxn
# Overlap between two circles 
Afxn <- function(Di,Dj,r){

# radii
rj <- Dj[3] # radius of target circle

# distance between circle centers
d <- ( (Di[1]-Dj[1])^2 + (Di[2]-Dj[2])^2)^0.5

if (d == 0){ # if they share a center
A = pi*( min(r,rj))^2}
else{
  if (d - rj > r){
    warning(c('Di=',as.character(Di)))
    warning(c('Dj=',as.character(Dj)))
    warning(c('r=',as.character(r)))
warning(c('d=',as.character(d)))
warning("Objects do not intersect, overlap is a complex number")}

# distance to radical line
x <- abs((d^2 - rj^2 + r^2)/(2*d)) 

# if the radical line is between the two circle centers
if (x <= d){ 
# Proceed with the formula proposed by Protazio:
  A <- intersect.Aij(d,rj,r)}
else{
# If not, define a new circle, cX.
if (r > rj){  #cX has radius r and distance from c2 of 2x - d
# Area is total c2 - overlap between cX & c2 + overlap between c1 & cX
c2.area <- pi*rj^2
cXc2.overlap <- intersect.Aij(2*x-d,r,rj)
c1cX.overlap <- intersect.Aij(2*x,r,r)
A = c2.area - cXc2.overlap + c1cX.overlap}
else{  # if r < rj, then the intersection happens behind the center of r
c1.area = pi*r^2
cXc1.overlap = intersect.Aij(2*x-d,rj,r)
c2cX.overlap = intersect.Aij(2*x,rj,rj)
A = c1.area - cXc1.overlap + c2cX.overlap}
} # end if/else x <= d
} # end if/else d == 0
A
} # end Afxn
#----------------------------------------------------------------  

#----------------------------------------------------------------  
# Afxn.rect 
# overlap between circle Di of radius r and rectangle defined by Vertices (SW,SE,NE, NW)

# Based on solution at http://stackoverflow.com/questions/622287/area-of-intersection-between-circle-and-rectangle
Afxn.rect <- function(Di,Vertices,r){ 

# Start with full circle area; subtract areas in each of the four external half-planes

Atotal <- pi*r^2

# signed (negative = outside) distances to each of the rectangle edges
d <- rep(0,4)
d[1] = (Di[1] - Vertices[1,1]) # distance to W side
d[2] = (Vertices[3,2] - Di[2]) # distance to N side
d[3] = (Vertices[2,1] - Di[1]) # distance to E side
d[4] = (Di[2] - Vertices[1,2]) # distance to S side

# should have already tested for circle being interior
if (any(d < 0)){
warning("Circle vertex is outside of study area")}

# Calculate overlap with each edge of the rectangle:
A.edge <- d*0
A.edge[d <= r] <- r * (r * acos( d[d <= r]/r ) - d[d <= r]*sin(acos( d[d <= r]/r)) ) # If d > r, acos returns NaN and a warning.


# Special case: overlaps a vertex
d2 <- d 
d2[5] <- d[1]; # for use in the following loop
A.c = rep(0,4)
Corners = c('SW','SE','NE','NW')
for (i in 1:4){ # rotate through each pair of corners
if (sqrt( (Di[1]-Vertices[i,1])^2 + (Di[2]-Vertices[i,2])^2) < r & d2[i] < r & d2[i+1] < r)  # if it intersects both sides & crosses the vertex
A.c[i] = A.rect.corner(Di,Vertices[i,],r,Corners[i])
} # end for-loop

# Area is total area, minus overlaps, plus corner overlap (to avoid double-counting)
A = Atotal - sum(A.edge) + sum(A.c)
A
} # END Afxn_rect
#----------------------------------------------------------------  

#----------------------------------------------------------------  
# A.rect.corner
# Overlap between circle (center Di, radius r) and quadrant defined by rectangle vertex V
A.rect.corner <- function(Di,V,r,Corner){
  
# 1) get intersection points with each of the two lines:
  IntX <- c(0,0)
  IntX[1] = (r^2 - (V[2] - Di[2])^2)^0.5 + Di[1]
  IntX[2] = -(r^2 - (V[2] - Di[2])^2)^0.5 + Di[1]
  
  IntY <- c(0,0)
  IntY[1] = (r^2 - (V[1] - Di[1])^2)^0.5 + Di[2]
  IntY[2] = -(r^2 - (V[1] - Di[1])^2)^0.5 + Di[2]

# 2) choose intersection points that are in the correct quadrant
Int = matrix(rep(0,4),2,2)
if (Corner == 'NW'){
Int[1,] = c(V[1], max(IntY)) # northmost along Y line
Int[2,] = c(min(IntX), V[2]) # westmost along X line
}
if (Corner == 'NE'){
Int[1,] = c(V[1], max[IntY]) # northmost along Y line
Int[2,] = c(max(IntX), V[2]) # eastmost along X line
}
if (Corner == 'SE'){
Int[1,] = c(V[1], min(IntY)) # southmost along Y line
Int[2,] = c(max(IntX), V[2]) # eastmost along X line
}
if (Corner == 'SW'){
Int[1,] = c(V[1], min(IntY)) # southmost along Y line
Int[2,] = c(min(IntX), V[2]) # westmost along X line
}

# 3) Find the length of the chord connecting the two intersection points
Z = sqrt( (Int[1,1] - Int[2,1])^2 + (Int[1,2] - Int[2,2])^2)

# 4) Find the angle of the segment
Th = asin( 0.5*Z/r)

# 5) Area of the segment
A.seg = (r^2) * Th

# 6) Area of the large isoceles triangle
A.tri.1 = r^2 * cos(Th) * sin(Th)

# 7) Area of the smaller triangle
A.tri.2 = 0.5 * ( abs(diff(Int[,1])) * abs(diff(Int[,2])) );

A = A.seg - A.tri.1 + A.tri.2
A
} # end A.rect.corner
#----------------------------------------------------------------

#----------------------------------------------------------------
# intersect.Aij 
# Area of intersection between two circles with radii rj and r and distance d separating centers
intersect.Aij <- function(d,rj,r){

if (d == 0){ # if they overlap (same center)
if (r > rj){ # overlap determined by the smaller circle
Aij = pi*rj^2}
else{
Aij = pi*r^2}}
else{
# Formula proposed by Protazio (and found correctly on Wolfram Mathworld) for an intersection forming a lens
Aij = (rj^2*acos( (d^2 + rj^2 - r^2) / (2*d*rj) ) + 
       r^2*acos( (d^2 + r^2 - rj^2) / (2*d*r) ) ) -
      (0.5 * ( (-d + rj + r)*(d + rj - r)*(d - rj + r)*(d + rj + r) )^0.5)
# (Note typo in Protazio's (2007) Eq. 4.3.  This is the correct version)
} # end if/else
Aij
} # end intersect.Aij
#----------------------------------------------------------------

#----------------------------------------------------------------
# is.overlap
# Is circle at i with radius r completely inside circle j? (returns a logical)
is.overlap <-function(Di,Dj,r){
  
  d <- ( (Di[1]-Dj[1])^2 + (Di[2]-Dj[2])^2)^0.5
  rj <- Dj[3]
  
  OL <- d + r < rj
  OL
} #end is.overlap
#----------------------------------------------------------------


#----------------------------------------------------------------
# is.overlap.rect
# Is circle at i with radius r completely inside rectangle with Vertices?
# Vertices listed as NW, NE, SE, SW
is.overlap.rect <-function(Di,Vertices,r){

OL <- abs(Di[1]) + r <= abs(Vertices[1,1]) &
      abs(Di[2]) + r <= abs(Vertices[1,2])
OL
} #end is.overlap.rect
#----------------------------------------------------------------
  
