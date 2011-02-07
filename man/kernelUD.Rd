\name{kernelUD}
\alias{kernelUD}
\alias{print.estUDm}
\alias{image.estUD}
\alias{image.estUDm}
\alias{plotLSCV}
\alias{getvolumeUD}
\alias{kernel.area}
\alias{estUDm2spixdf}
\alias{as.data.frame.estUD}
\title{Estimation of Kernel Home-Range}
\description{
  The function \code{kernelUD} estimates the UD of one or several
  animals.\cr
  \code{plotLSCV} allows to explore the results of the least-square
  cross-validation algorithm used to find the best smoothing value.\cr
  \code{image} allows a graphical display of the estimates. \cr
  \code{getvolumeUD} and \code{kernel.area} provide utilities
  for home range and home-range size estimation.\cr
  \code{getverticeshr} stores the home range contour as an object of
  class \code{SpatialPolygonsDataFrame} (package sp), with one row per
  animal.\cr
  \code{estUDm2spixdf} can be used to convert the result into an object
  of class \code{SpatialPixelsDataFrame}\cr
  \code{as.data.frame.estUD} can be used to convert an object of class
  \code{estUD} as a data frame.
}
\usage{

kernelUD(xy, h = "href", grid = 60,
         same4all = FALSE, hlim = c(0.1, 1.5),
         kern = c("bivnorm", "epa"), extent = 0.5,
         boundary = NULL)

\method{print}{estUDm}(x, \dots)

\method{image}{estUD}(x, \dots)

\method{image}{estUDm}(x, \dots)

\method{as.data.frame}{estUD}(x, row.names, optional, \dots)

plotLSCV(x)

getvolumeUD(x, standardize = FALSE)

kernel.area(x, percent = seq(20, 95, by = 5),
            unin = c("m", "km"),
            unout = c("ha", "km2", "m2"), standardize = FALSE)

estUDm2spixdf(x)

}

\arguments{
  \item{xy}{An object inheriting the class \code{SpatialPoints}
    containing the x and y relocations of the animal.  If \code{xy}
    inherits the class \code{SpatialPointsDataFrame}, it should contain
    only one column (factor) corresponding to the identity of the
    animals for each relocation.}
  \item{h}{a character string or a number.  If \code{h} is set to
    \code{"href"}, the ad hoc method is used for the smoothing parameter
    (see details).  If \code{h} is set to 
    \code{"LSCV"}, the least-square cross validation method is
    used. Note that \code{"LSCV"} is not available if \code{kern = "epa"}. 
    Alternatively, \code{h} may be set to any given numeric value}
  \item{grid}{a number giving the size of the grid on
    which the UD should be estimated.  Alternatively, this parameter may
    be an object inheriting the class \code{SpatialPixels}.  For the
    function \code{kernelUD}, it may in addition be a list of objects of
    class \code{SpatialPixels}, with named elements corresponding to each
    level of the factor id.}
  \item{hlim}{a numeric vector of length two.  If \code{h = "LSCV"},
    the function minimizes the cross-validation criterion for values of
    h ranging from \code{hlim[1]*href} to \code{hlim[2]*href}, where
    \code{href} is the smoothing parameter computed with the ad hoc
    method (see below)}
  \item{kern}{a character string.  If \code{"bivnorm"}, a bivariate
    normal kernel is used.  If \code{"epa"}, an Epanechnikov kernel is
    used.}
  \item{extent}{a value controlling the extent of the grid used for the
    estimation (the extent of the grid on the abscissa is equal to
    \code{(min(abscissa.relocations) + extent *
      diff(range(abscissa.relocations)))}, and similarly for the ordinate).}
  \item{same4all}{logical. If \code{TRUE}, the same grid
    is used for all animals.  If \code{FALSE}, one grid per
    animal is used}
  \item{boundary}{If, not \code{NULL}, an object inheriting the class
    \code{SpatialLines} defining a barrier that cannot be crossed by the
    animals.  There are constraints on the shape of the barrier that
    depend on the smoothing parameter \code{h} (***see details***)}
  \item{x}{an object of class \code{estUD} (UD for one animal) or
    \code{estUDm} (UD for several animals).  For the 
    function \code{estUDm2spixdf}, an object of class \code{estUDm}
    only.  For the function \code{as.data.frame.estUD}, an object of
    class \code{estUD} only.}
  \item{percent}{for \code{kernel.area}, a vector of percentage levels
    for home-range size estimation.  For \code{getverticeshr}, a single
    value giving the percentage level for home-range estimation.}
  \item{standardize}{a logical value indicating whether the UD should be
    standardized over the area of interest, so that the volume under the
    UD and *over the area* is equal to 1.}
  \item{unin}{the units of the relocations coordinates.  Either \code{"m"} for
    meters (default) or \code{"km"} for kilometers}
  \item{unout}{the units of the output areas.  Either \code{"m2"} for
    square meters, \code{"km2"} for square kilometers or \code{"ha"} for
    hectares (default)}
  \item{row.names}{unused argument here}
  \item{optional}{unused argument here}
  \item{\dots}{additionnal parameters to be passed to the generic
    functions \code{print} and \code{image}}
}
\details{
  The Utilization Distribution (UD) is the bivariate function giving the
  probability density that an animal is found at a point according to
  its geographical coordinates.  Using this model, one can define 
  the home range as the minimum area in which an animal has some
  specified probability of being located. The
  functions used here correspond to the approach described in Worton
  (1995).\cr

  The kernel method has been recommended by many authors for the
  estimation of the utilization distribution (e.g. Worton, 1989, 1995).
  The default method for the estimation of the
  smoothing parameter is the \emph{ad hoc} method, i.e. for a bivariate
  normal kernel
  \deqn{h = \sigma n^{- \frac{1}{6}}}{h = Sigma*n^(-1/6)}
  where \deqn{\sigma = 0.5 (\sigma(x)+\sigma(y))}{Sigma =
    0.5*(sd(x)+sd(y))} which supposes that the UD is
  bivariate normal.  If an Epanechnikov kernel is used, this value is
  multiplied by 1.77 (Silverman, 1986, p. 86).
  Alternatively, the smoothing parameter h may be
  computed by Least Square Cross Validation (LSCV).  The estimated value
  then minimizes the Mean Integrated Square Error (MISE), i.e. the
  difference in volume between the true UD and the estimated UD.  Note
  that the cross-validation criterion cannot be minimized in some
  cases.  According to Seaman and Powell (1998) \emph{"This is a
  difficult problem that has not been worked out by statistical
  theoreticians, so no definitive response is available at this
  time"} (see Seaman and Powell, 1998 for further details and tricky
  solutions).  \code{plotLSCV} allows to have a diagnostic of the
  success of minimization of the cross validation criterion (i.e. to
  know whether the minimum of the CV criterion occurs within the scanned
  range).  Finally, the UD is then estimated over a grid.\cr
  
  The default kernel is the bivariate normal kernel, but the
  Epanechnikov kernel, which requires less computer time is also
  available for the estimation of the UD. \cr
  
  The function \code{getvolumeUD} modifies the UD component of the
  object passed as argument: that the pixel values of the resulting
  object are equal to the percentage of the smallest home range
  containing this pixel.  This function is used in the function
  \code{kernel.area}, to compute the home-range size.  Note, that the
  function \code{plot.hrsize} (see the help page of this  function) can
  be used to display the home-range size estimated at various levels.\cr

  The parameter \code{boundary} allows to define a barrier that cannot
  be crossed by the animals.  When this parameter is
  set, the method described by Benhamou and Cornelis (2010) for
  correcting boundary biases is used.  The boundary can possibly be
  defined by several nonconnected lines, each one being built by several
  connected segments. Note that there are constraints on these segments
  (not all kinds of boundary can be defined): (i) each segment length
  should at least be equal to \code{3*h} (the size of "internal lane"
  according to the terminology of Benhamou and Cornelis), (ii) the angle
  between two line segments should be greater that \code{pi/2} or lower
  that \code{-pi/2}.  The UD of all the pixels located within a
  band defined by the boundary and with a width equal to \code{6*h}
  ("external lane") is set to zero.
}
\value{
  
  The function \code{kernelUD} returns either: (i) an object belonging
  to the S4 class \code{estUD} (see \code{?estUD-class}) when the object
  \code{xy} passed as argument contains the relocations of only one
  animal (i.e., belong to the class \code{SpatialPoints}), or (ii) a
  list of elements of class \code{estUD} when the object
  \code{xy} passed as argument contains the relocations of several
  animals (i.e., belong to the class \code{SpatialPointsDataFrame}).\cr
  
  The function \code{getvolumeUD} returns an object of the same class as
  the object passed as argument (\code{estUD} or \code{estUDm}).  \cr
  
  \code{kernel.area} returns a data frame of subclass \code{hrsize},
  with one column per animal and one row per level of
  estimation of the home range.\cr

  \code{getverticeshr} returns an object of class
  \code{SpatialPolygonsDataFrame}.\cr

  \code{estUDm2spixdf} returns an object of class
  \code{SpatialPixelsDataFrame}.
  
}
\references{

  Silverman, B.W. (1986)
  \emph{Density estimation for statistics and data
    analysis}. London: Chapman \& Hall.
  
  Worton, B.J. (1989) Kernel methods for estimating the utilization
  dirstibution in home-range studies. \emph{Ecology}, \bold{70}, 164--168.

  Worton, B.J. (1995) Using Monte Carlo simulation to evaluate
  kernel-based home range estimators. \emph{Journal of Wildlife Management},
  \bold{59},794--800.

  Seaman, D.E. and Powell, R.A. (1998) \emph{Kernel home range estimation
    program (kernelhr)}. Documentation of the
  program. \url{ftp://ftp.im.nbs.gov/pub/software/CSE/wsb2695/KERNELHR.ZIP}.

  Benhamou, S. and Cornelis, D. (2010) Incorporating Movement Behavior
  and Barriers to Improve Biological Relevance of Kernel Home Range
  Space Use Estimates. \emph{Journal of Wildlife Management}, \bold{74},
  1353--1360.

}
\author{ Clement Calenge \email{clement.calenge@oncfs.gouv.fr} }

\seealso{\code{\link{mcp}} for help on the function \code{plot.hrsize}. 
}

\examples{

## Load the data
data(puechabonsp)
loc <- puechabonsp$relocs

## have a look at the data
head(as.data.frame(loc))
## the first column of this data frame is the ID


## Estimation of UD for the four animals
(ud <- kernelUD(loc[,1]))

## The UD of the four animals
image(ud)

## Calculation of the 95 percent home range
ver <- getverticeshr(ud, 95)

## and display on an elevation map:
elev <- puechabonsp$map
image(elev, 1)
plot(ver, add=TRUE, col=rainbow(4))
legend(699000, 3165000, legend = names(ud), fill = rainbow(4))

## Example of estimation using LSCV
udbis <- kernelUD(loc[,1], h = "LSCV")
image(udbis)


## Compare the estimation with ad hoc and LSCV method
## for the smoothing parameter
(cuicui1 <- kernel.area(ud)) ## ad hoc
plot(cuicui1)
(cuicui2 <- kernel.area(udbis)) ## LSCV
plot(cuicui2)

## Diagnostic of the cross-validation
plotLSCV(udbis)



## Use of the same4all argument: the same grid
## is used for all animals
## BTW, we indicate a grid with a fine resolution:
udbis <- kernelUD(loc[,1], same4all = TRUE, grid = 100)
image(udbis)


## Estimation of the UD on a map
## (e.g. for subsequent analyses on habitat selection)
## Measures the UD in each pixel of the map
udbis <- kernelUD(loc[,1], grid = elev)
image(udbis)


##########################################
##
## Estimating the UD with the presence of a barrier
## The boars are located on the plateau of Puechabon (near
## Montpellier, France), and their movements are limited by the
## Herault river.

## We first map the elevation:
image(elev)

## Then, we used the function locator() to identify the limits of the
## segments of this barrier. BEWARE! The boundary should satisfy the two
## constraints: (i) segment length > 3*h, (ii) no angle lower than pi/2
## between successive segments. We choose a smoothing parameter of 100
## m, so that no segment length should be less than 300 m.
## Because the resolution of the map is 100 m, this means that no
## segment should cover less than 3 pixels. We have used the function
## locator() to digitize this barrier and then the function dput to
## have the following limits:

bound <- structure(list(x = c(701751.385381925, 701019.24105475,
                        700739.303517889,
                        700071.760160759, 699522.651915378,
                        698887.40904327, 698510.570051342,
                        698262.932999504, 697843.026694212,
                        698058.363261028),
                        y = c(3161824.03387414,
                        3161824.03387414, 3161446.96718494,
                        3161770.16720425, 3161479.28718687,
                        3161231.50050539, 3161037.5804938,
                        3160294.22044937, 3159389.26039528,
                        3157482.3802813)), .Names = c("x", "y"))

lines(bound, lwd=3)

## We convert bound to SpatialLines:
bound <- do.call("cbind",bound)
Slo1 <- Line(bound)
Sli1 <- Lines(list(Slo1), ID="frontier1")
barrier <- SpatialLines(list(Sli1))

## estimation of the UD
kud <- kernelUD(loc[,1], h=100, grid=100, boundary=barrier)

## Result:
image(kud)

## Have a closer look to Calou:
kud2 <- kud[[2]]
image(kud2, col=grey(seq(1,0,length=15)))
title(main="Home range of Calou")
points(loc[slot(loc,"data")[,1]=="Calou",], pch=3, col="blue")
plot(getverticeshr(kud2, 95), add=TRUE, lwd=2)
lines(barrier, col="red", lwd=3)



}
\keyword{spatial}

