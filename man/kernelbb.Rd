\name{kernelbb}
\alias{kernelbb}
\alias{print.liker}
\alias{liker}
\title{Estimation of Kernel Brownian Bridge Home-Range}
\description{
  \code{kernelbb} is used to estimate the utilization distribution of an
  animal using the brownian bridge approach of the kernel method (for
  autocorrelated relocations; Bullard 1991, Horne et al. 2007). \cr
  \code{liker} can be used to find the maximum likelihood
  estimation of the parameter sig1, using the approach defined in Horne
  et al. 2007 (see Details).
}
\usage{
kernelbb(ltr, sig1, sig2, grid = 40, same4all = FALSE, byburst = FALSE,
         extent = 0.5, nalpha = 25)

liker(tr, rangesig1, sig2, le = 1000,
      byburst = FALSE, plotit = TRUE)

\method{print}{liker}(x, ...)
}

\arguments{
  \item{ltr,tr}{an object of class \code{ltraj} of type II (time recorded),
    regular or not (see \code{help(as.ltraj)}).}
  \item{sig1}{first smoothing parameter for the brownian bridge method
    (related to the speed of the animals; it can be estimated by the
    function \code{liker}).}
  \item{sig2}{second smoothing parameter for the brownian bridge method
    (related to the imprecision of the relocations, supposed known).}
  \item{grid}{a number giving the size of the grid on
    which the UD should be estimated.  Alternatively, this parameter may
    be an object of class \code{SpatialPixels}, or a list of objects of
  class \code{SpatialPixels}, with named elements corresponding to each
  level of the factor id}
  \item{same4all}{logical. If \code{TRUE}, the same grid
    is used for all animals.  If \code{FALSE}, one grid per
    animal is used}
  \item{byburst}{logical. Whether the brownian bridge estimation should
    be done by burst.}
  \item{extent}{a value indicating the extent of the grid used for the
    estimation (the extent of the grid on the abscissa is equal to
    \code{(min(xy[,1]) + extent * diff(range(xy[,1])))}).}
  \item{nalpha}{a parameter used internally to compute the integral
    of the Brownian bridge.  The integral is computed by cutting each
    step built by two relocations into \code{nalpha} sub-intervals. }
  \item{rangesig1}{the range of possible values of sig1 within which
    the likelihood should be maximized.}
  \item{le}{The number of values of sig1 tested within the specified
    range.}
  \item{plotit}{logical. Whether the results of the function should be
    plotted.}
  \item{x}{an object of class \code{khr} returned by \code{kernelbb}.}
  \item{\dots}{additionnal parameters to be passed to the generic
    functions \code{print}}
}
\details{
  The function \code{kernelbb} uses the brownian bridge approach to
  estimate the Utilization Distribution of an animal with serial
  autocorrelation of the relocations (Bullard 1991, Horne et al. 2007).
  Instead of simply smoothing the relocation pattern (which is the case
  for the function \code{kernelUD}), it takes into account the fact that
  between two successive relocations r1 and r2, the animal has moved
  through a continuous path, which is not necessarily linear.  A
  brownian bridge estimates the density of probability that this path
  passed through any point of the study area, given that the animal was
  located at the point r1 at time t1 and at the point r2 at time t2,
  with a certain amount of inaccuracy (controled by the parameter sig2,
  see Examples).  Brownian bridges are placed over the different
  sections of the trajectory, and these functions are then summed over the
  area.  The brownian bridge approach therefore smoothes a trajectory.\cr
  
  The brownian bridge estimation relies on two smoothing parameters,
  \code{sig1} and \code{sig2}.  The parameter \code{sig1} is related to
  the speed of the animal, and describes how far from the line joining
  two successive relocations the animal can go during one time unit
  (here the time is measured in second). The function \code{liker} can
  be used to estimate this value using the maximum likelihood approach
  described in Horne et al. (2007). The larger this parameter is,
  and the more wiggly the trajectory is likely to be.
  The parameter \code{sig2} is equivalent to the parameter 
  \code{h} of the classical kernel method: it is related to the
  inaccuracy of the relocations, and is supposed known (See examples for
  an illustration of 
  the smoothing parameters).\cr

  The functions \code{getvolumeUD} and \code{getverticeshr} can then be
  used to conpute the home ranges (see \code{kernelbb}).  More
  generally, more details on the generic parameters of \code{kernelbb}
  can be found on the help page of \code{kernelUD}.
  
}
\value{
  An object of class \code{estUDm}

  \code{liker} returns an object of class \code{liker}, with one
  component per animal (or per burst, depending on the value of
  the parameter \code{perburst}), containing the value of (i) optimized
  sig1, (ii) sig2, and (iii) a data frame named "cv" with the tested
  values of sig1 and the corresponding log-likelihood.
}
\references{

  Bullard, F. (1991) \emph{Estimating the home range of an animal: a Brownian
    bridge approach}. Master of Science, University of North Carolina,
  Chapel Hill.

  Horne, J.S., Garton, E.O., Krone, S.M. and Lewis, J.S. (2007)
  Analyzing animal movements using brownian bridge. \emph{Ecology},
  \bold{in press}.
}
\author{ Clement Calenge \email{clement.calenge@oncfs.gouv.fr} }

\seealso{\code{\link{as.ltraj}} for further information concerning
  objects of class \code{ltraj}. \code{\link{kernelUD}} for the
  classical kernel estimation. ,  \code{\link{mcp}} for
  estimation of home ranges using the minimum convex polygon, and for
  help on the function \code{plot.hrsize}. 
}
\examples{

\dontrun{

#########################################################
#########################################################
#########################################################
###
###        Example of a typical case study
###       with the brownian bridge approach
###

## Load the data
data(puechcirc)
x <- puechcirc[1]

## Field studies have shown that the mean standard deviation (relocations
## as a sample of the actual position of the animal) is equal to 58
## meters on these data (Maillard, 1996, p. 63). Therefore
sig2 <- 58

## Find the maximum likelihood estimation of the parameter sig1
## First, try to find it between 10 and 100.
liker(x, sig2 = 58, rangesig1 = c(10, 100))

## Wow! we expected a too large standard deviation! Try again between
## 1 and 10:
liker(x, sig2 = 58, rangesig1 = c(1, 10))

## So that sig1 = 6.23

## Now, estimate the brownian bridge
tata <- kernelbb(x, sig1 = 6.23, sig2 = 58, grid = 100)
image(tata)


## OK, now look at the home range
image(tata)
plot(getverticeshr(tata, 95), add=TRUE, lwd=2)





#########################################################
#########################################################
#########################################################
###
###      Comparison of the brownian bridge approach
###            with the classical approach
###


## Take an illustrative example: we simulate a trajectory
set.seed(2098)
pts1 <- data.frame(x = rnorm(25, mean = 4.5, sd = 0.05),
                   y = rnorm(25, mean = 4.5, sd = 0.05))
pts1b <- data.frame(x = rnorm(25, mean = 4.5, sd = 0.05),
                    y = rnorm(25, mean = 4.5, sd = 0.05))
pts2 <- data.frame(x = rnorm(25, mean = 4, sd = 0.05),
                   y = rnorm(25, mean = 4, sd = 0.05))
pts3 <- data.frame(x = rnorm(25, mean = 5, sd = 0.05),
                   y = rnorm(25, mean = 4, sd = 0.05))
pts3b <- data.frame(x = rnorm(25, mean = 5, sd = 0.05),
                    y = rnorm(25, mean = 4, sd = 0.05))
pts2b <- data.frame(x = rnorm(25, mean = 4, sd = 0.05),
                    y = rnorm(25, mean = 4, sd = 0.05))
pts <- do.call("rbind", lapply(1:25, function(i) {
          rbind(pts1[i,], pts1b[i,], pts2[i,], pts3[i,],
                pts3b[i,], pts2b[i,])
}))
dat <- 1:150
class(dat) <- c("POSIXct","POSIXt")
x <- as.ltraj(pts, date=dat, id = rep("A", 150))

## See the trajectory:
plot(x)


## Now, we suppose that there is a precision of 0.05
## on the relocations
sig2 <- 0.05
## and that sig1=0.1
sig1 <- 0.1

## Now fits the brownian bridge home range
(kbb <- kernelbb(x, sig1 = sig1,
                 sig2 = sig2))


## Now fits the classical kernel home range
coordinates(pts) <- c("x","y")
(kud <- kernelUD(pts))


###### The results

opar <- par(mfrow=c(2,2), mar=c(0.1,0.1,2,0.1))
plot(pts, pch=16)
title(main="The relocation pattern")
box()
plot(x, axes=FALSE, main="The trajectory")
box()
image(kud)
title(main="Classical kernel home range")
plot(getverticeshr(kud, 95), add=TRUE)
box()
image(kbb)
title(main="Brownian bridge kernel home range")
plot(getverticeshr(kbb, 95), add=TRUE)
box()
par(opar)



###############################################
###############################################
###############################################
###
###         Image of a brownian bridge.
###          Fit with two relocations
###               


xx <- c(0,1)
yy <- c(0,1)
date <- c(0,1)
class(date) <- c("POSIXt", "POSIXct")
tr <- as.ltraj(data.frame(x = xx,y = yy), date, id="a")

## Use of different smoothing parameters
sig1 <- c(0.05, 0.1, 0.2, 0.4, 0.6)
sig2 <- c(0.05, 0.1, 0.2, 0.5, 0.7)

y <- list()
for (i in 1:5) {
  for (j in 1:5) {
     k <- paste("s1=", sig1[i], ", s2=", sig2[j], sep = "")
     y[[k]]<-kernelbb(tr, sig1[i], sig2[j])
   }
 }

## Displays the results
opar <- par(mar = c(0,0,2,0), mfrow = c(5,5))
foo <- function(x)
   {
     image(y[[x]])
     title(main = names(y)[x])
     points(tr[[1]][,c("x","y")], pch = 16)
   }
lapply(1:length(y), foo)

par(opar)


}
}
\keyword{spatial}
