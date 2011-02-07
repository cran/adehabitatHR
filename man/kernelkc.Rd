\name{kernelkc}
\alias{kernelkc}
\alias{kernelkcbase}
\alias{getvolumeUDk}
\alias{getvolumeUDs}
\alias{getverticeshrk}
\alias{getverticeshrs}
\alias{exwc}
\title{Kernel Smoothing in Space and Time of the Animals' Use of Space}
\description{
  These functions estimate the utilization distribution (UD) in space
  and time of animals monitored using radio-telemetry, using the product
  kernel estimator advocated by Keating and Cherry (2009).\cr

  Note that this approach has also been useful for the analysis of
  recoveries in programs involving ringed birds (Calenge et al. in
  prep., see section examples below).\cr
  
  \code{kernelkc} estimate the UD of several animals from an object of
  class \code{ltraj}.\cr  
  \code{kernelkcbase} estimate one UD from a data frame with three
  columns indicating the spatial coordinates and associated timing.\cr  
  \code{exwc} allows to search for the best value of the
  time smoothing parameter in the case where the time is considered as a
  circular variable (see details).
}
\usage{
kernelkc(tr, h, tcalc, t0, grid = 40, circular = FALSE,
         cycle = 24 * 3600, same4all = FALSE,
         byburst = FALSE, extent = 0.5)

kernelkcbase(xyt, h, tcalc, t0, grid=40, circular=FALSE,
             cycle=24*3600, extent=0.5)

exwc(hv)

}
\arguments{
  \item{tr}{an object of class \code{ltraj}}
  \item{xyt}{a data frame with three columns indicating the x and y
    coordinates, as well as the timing of the relocations. }
  \item{h}{a numeric vector with three elements indicating the value of
    the smoothing parameters: the first and second elements are
    the smoothing parameters of the X and Y coordinates respectively,
    the third element is the smoothing parameter for the time
    dimension.  If \code{circular=TRUE} it should be a smoothing
    parameter in the interval 0-1 (see details).  If
    \code{circular=FALSE} this smoothing parameter should be given in
    seconds.}
  \item{tcalc}{the time at which the UD is to be estimated}
  \item{t0}{if \code{circular=TRUE}, this parameter indicates the time
    at which the time cycle begins (see examples).}
  \item{grid}{a number giving the size of the grid on which the UD should
    be estimated.  Alternatively, this parameter may be an object
    of class \code{asc}, or a list of objects of class \code{asc}, with
    named elements corresponding to each level of the burst/id}
  \item{circular}{logical.  Indicates whether the time should be
    considered as a circular variable (e.g., the 31th december 2007 is
    considered to be one day before the 1st january 2007) or not (e.g.,
    the 31th december 2007 is considered to be one year after the 1st
    january 2007).}
  \item{cycle}{if \code{circular=TRUE}, the duration of the time cycle.
    for \code{kernelkc}, it should be given in seconds, and for
    \code{kernelkcbase}, in the units of the data (the units of the
    third column of \code{xyt}).}
  \item{same4all}{logical.  If \code{TRUE}, the same grid is used for all
    levels of id/burst.  If \code{FALSE}, one grid per id/burst is used.}
  \item{byburst}{logical.  Indicates whether one UD should be estimated
    by burst of \code{tr}, or whether the data should be pooled across
    all bursts of each value of id in \code{tr}}
  \item{extent}{a value indicating the extent of the grid used for the
    estimation (the extent of the grid on the abscissa is equal
    to \code{(min(xy[,1]) + extent * diff(range(xy[,1])))}). }
  \item{hv}{a value of smoothing parameter for the time dimension.}
  \item{\dots}{additional arguments to be passed to the function
    \code{contour}.}
}
\details{
  Keating and Cherry (2009) advocated the estimation of the UD in
  time and space using the product kernel estimator.  These functions
  implement exactly this methodology.\cr

  For the spatial coordinates, the
  implemented kernel function is the biweight kernel.\cr

  Two possible
  approaches are possible to manage the time in the estimation process:
  (i) the time may be considered as a linear variable (e.g., the 31th
  december 2007 is considered to be one day before the 1st january
  2007), or (ii) the time may be considered as a circular variable
  (e.g., the 31th december 2007 is considered to be one year after the
  1st january 2007).\cr

  If the time is considered as a linear variable, the kernel function
  used in the estimation process is the biweight kernel.  If the time is
  considered as a circular variable, the implemented kernel is the
  wrapped Cauchy distribution (as in the article of Keating and Cherry).
  In this latter case, the smoothing parameter should be chosen in the
  interval 0-1, with a value of 1 corresponding to a stronger
  smoothing.\cr

  These functions can only be used on objects of class "ltraj", but
  the estimation of the UD in time and space is also possible with other
  types of data (see the help page of \code{kernelkcbase}).  Note that
  both \code{kernelkc} and \code{kernelkcbase} return conditional
  probability density function (pdf), i.e. the pdf to relocate an animal
  at a place, given that it has been relocated at time \code{tcalc}
  (i.e. the volume under the UD estimated at time \code{tcalc} is equal
  to 1 whatever \code{tcalc}).\cr 

  The function \code{exwc} draws a graph of the wrapped
  Cauchy distribution for the chosen \code{h} parameter (for circular
  time), so that it is possible to make one's mind concerning the weight
  that can be given to the neighbouring points of a given time point.
  Note that although Keating and Cherry (2009) advocated the use of
  an automatic algorithm to select "optimal" values for the smoothing
  parameter, it is not implemented in adehabitatHR.  Indeed, different
  smoothing parameters may allow to identify patterns at different
  scales, and we encourage the user to try several values before
  subjectively choosing the value which allows to more clearly identify
  the patterns of the UD.


}
\value{
  \code{kernelkc} returns a list of class "\code{estUDm}" containing
  objects of class \code{estUD}, mapping one estimate of the UD per burst
  or id (depending on the value of the parameter \code{byburst}).\cr

  \code{kernelkcbase} returns an object of class "\code{estUD}" mapping
  the  estimated UD.
}
\references{
  Keating, K. and Cherry, S. (2009) Modeling utilization
  distributions in space and time. \emph{Ecology}, \bold{90}:
  1971--1980.

  Calenge, C., Guillemain, M., Gauthier-Clerc, M. and Simon, G. in
  prep.  A new exploratory approach to the study of the spatio-temporal
  distribution of ring recoveries - the example of Teal (Anas crecca)
  ringed in Camargue, Southern France.

}
\author{Clement Calenge \email{clement.calenge@oncfs.gouv.fr}}
\seealso{\code{\link{as.ltraj}} for additional information on objects of
  class \code{ltraj}, \code{\link{kernelUD}} for the "classical" kernel
  home range estimates. }
\examples{

\dontrun{


################################################
##
## Illustrates the analysis of recoveries of
## ringed data

data(teal)
head(teal)

## compute the sequence of dates at which the
## probability density function (pdf) of recoveries is to be estimated

vv <- seq(min(teal$date), max(teal$date), length=50)
head(vv)

## The package "maps" should be installed for the example below
library(maps)


re <- lapply(1:length(vv), function(i) {

             ## Estimate the pdf. We choose a smoothing parameter of
             ## 2 degrees of lat-long for X and Y coordinates,
             ## and of 2 months for the time
             uu <- kernelkcbase(teal, c(2.5,2.5,2*30*24*3600), tcalc =
                                vv[i], grid=100, extent=0.1)

             ## now, we show the result
             ## potentially, we could type
             ##
             ## jpeg(paste("prdefu", i, ".jpg", sep=""))
             ##
             ## to store the figures in a file, and then to build a
             ## movie with the resulting files:
             ##

             image(uu, col=grey(seq(1,0, length=8)))
             title(main=vv[i])

             ## highlight the area on which there is a probability
             ## equal to 0.95 to recover a bird
             ## ****warning! The argument standardize=TRUE should
             ## be passed, because the UD is defined in space and
             ## time, and because we estimate the UD just in space
             plot(getverticeshr(uu, 95, standardize=TRUE), add=TRUE,
                  border="red", lwd=2)

             ## The map:
             map(xlim=c(-20,70), ylim=c(30,80), add=TRUE)

             ## and if we had typed jpeg(...) before, we have to type
             ## dev.off()
             ## to close the device. When we have finished this loop
             ## We could combine the resulting files with imagemagick
             ## (windows) or mencoder (linux)
             })




################################################
##
## Illustrates how to explore the UD in time and
## space with the bear data

data(bear)

## compute the sequence of dates at which the UD is to be
## estimated
vv <- seq(min(bear[[1]]$date), max(bear[[1]]$date), length=50)
head(vv)

## estimates the UD at each time point
re <- lapply(1:length(vv), function(i) {

             ## estimate the UD. We choose a smoothing parameter of
             ## 1000 meters for X and Y coordinates, and of 72 hours
             ## for the time (after a visual exploration)
             uu <- kernelkc(bear, h = c(1000,1000,72*3600),
                            tcalc= vv[i], grid=100)

             ## now, we show the result
             ## potentially, we could type
             ##
             ## jpeg(paste("UD", i, ".jpg", sep=""))
             ##
             ## to store the figures in a file, and then to build a
             ## movie with the resulting files:
             ##
             image(uu, col=grey(seq(1,0,length=10)))
             title(main=vv[i])

             ## highlight the 95 percent home range
             ## we set standardize = TRUE because we want to estimate
             ## the home range in space from a UD estimated in space and
             ## time
             plot(getverticeshr(uu, 95, standardize=TRUE), lwd=2,
                  border="red", add=TRUE)

             ## and if we had typed jpeg(...) before, we have to type
             ## dev.off()
             ## to close the device. When we have finished this loop
             ## We could combine the resulting files with imagemagick
             ## (windows) or mencoder (linux)
             })


## Or, just show the home range:
re <- lapply(1:length(vv), function(i) {

             uu <- kernelkc(bear, h = c(1000,1000,72*3600),
                            tcalc= vv[i])

             pc <- getverticeshr(uu, 95, standardize=TRUE)
             plot(pc, xlim=c(510000, 530000),
                  ylim=c(6810000, 6825000))
             title(main=vv[i])
             })








##################################################
##
## Example with several wild boars (linear time)

## load wild boar data
data(puechcirc)


## keep only the first two circuits:
puechc <- puechcirc[1:2]


## Now load the map of the elevation
data(puech)


## compute the time point at which the UD is to be estimated
vv <- seq(min(puechcirc[[2]]$date), max(puechcirc[[2]]$date),
          length=50)


## The estimate the UD
re <- lapply(1:length(vv),
             function(i) {

                 ## We choose a smoothing parameter of 300 meters for
                 ## the x and y coordinates and of one hour for the time
                 ## (but try to play with these smoothing parameters)

                 uu <- kernelkc(puechcirc, h=c(300,300,3600),
                                tcalc = vv[i], same4all=TRUE,
                                extent=0.1)

                 ## show the elevation
                 image(puech$maps, 
                       xlim=c(698000,704000),
                       ylim=c(3156000,3160000))
                 title(main=vv[i])

                 ## and the UD, with contour lines
                 colo <- c("green","blue")
                 lapply(1:length(uu), function(i) {
                        contour(as(uu[[i]],"SpatialPixelsDataFrame"),
                                add=TRUE, col=colo[i])
                 })

                 ## the blue contour lines show the UD of the mother and
                 ## the red ones correspond to her son. Adult wild boars
                 ## are known to be more "shy" that the youger ones.
                 ## Here, the low elevation corresponds to crop area
                 ## (vineyards). The young boar is the first and the
                 ## last in the crops
             })








##################################################
##
## Example with the bear, to illustrate (circular time)

data(bear)

## We consider a time cycle of 24 hours.
## the following vector contains the time points on the
## time circle at which the UD is to be estimated (note that
## the time is given in seconds)
vv <- seq(0, 24*3600-1, length=40)


## for each time point:
re <- lapply(1:length(vv),
             function(i) {

                 ## Estimation of the UD for the bear. We choose
                 ## a smoothing parameter of 1000 meters for the spatial
                 ## coordinates and a smoothing parameter equal to 0.2
                 ## for the time. We set the beginning of the time
                 ## cycle at midnight (no particular reason, just to
                 ## illustrate the function). So we pass, as t0, any
                 ## object of class POSIXct corresponding t a date at
                 ## this hour, for example the 12/25/2012 at 00H00
                 t0 <- as.POSIXct("2012-12-25 00:00")
                 uu <- kernelkc(bear, h=c(1000,1000,0.2), cycle=24*3600,
                                 tcalc=vv[i], t0=t0, circular=TRUE)

                 ## shows the results
                 ## first compute the hour for the title
                 hour <- paste(floor(vv[i]/3600), "hours",
                               floor((vv[i]\%\%3600)/60), "min")

                 ## compute the 95\% home range
                 pc <- getverticeshr(uu, 95, standardize=TRUE)
                 plot(pc, xlim=c(510000, 530000),
                      ylim=c(6810000, 6825000))
                 title(main=hour)

                 ## compute the 50\% home range
                 pc <- getverticeshr(uu, 50, standardize=TRUE)
                 plot(pc, add=TRUE, col="blue")
                 
             })
 ## Now, each home range computed at a given time point corresponds to
 ## the area used by the animal at this time period.  We may for example
 ## try to identify the main difference in habitat composition of the
 ## home-range between different time, to identify differences in
 ## habitat use between different time of the day. We do not do it here
 ## (lack of example data)










##################################################
##
## Example of the use of the function kernelkcbase and
## related functions


## load the data
data(puechabon)
locs <- puechabon$locs

## keeps only the wild boar Jean
locs <- locs[slot(locs, "data")[,1]=="Jean",]

## compute the number of days since the beginning
## of the monitoring
dd <- cumsum(c(0, diff(strptime(slot(locs, "data")[,4], "\%y\%m\%d"))))
dd

## compute xyt. Note that t is here the number of
## days since the beginning of the monitoring (it
## is not an object of class POSIXt, but it may be)
xyt <- data.frame(as.data.frame(coordinates(locs)), dd)

## Now compute the time points at which the UD is to be estimated:
vv <- 1:61

## and finally, show the UD changed with time:
re <- lapply(1:length(vv),
             function(i) {
                 ud <- kernelkcbase(xyt, h=c(300,300,20), 
                                    tcalc=vv[i], grid=100)
                 image(ud, main=vv[i])
                 plot(getverticeshr(ud, 95, standardize=TRUE),
                       border="red", lwd=2, add=TRUE)

                 ## Just to slow down the process
                 Sys.sleep(0.2)
                 })



}



}
\keyword{spatial}
