\name{LoCoH}
\alias{LoCoH}
\alias{LoCoH.k}
\alias{LoCoH.r}
\alias{LoCoH.a}
\alias{LoCoH.k.area}
\alias{LoCoH.r.area}
\alias{LoCoH.a.area}
\title{Estimating LoCoH home ranges}
\description{
  The functions computes the home range of one or several animals using
  the LoCoH family of methods.\cr
  The functions \code{LoCoH.k}, \code{LoCoH.r}, and \code{LoCoH.a}
  implement the k-LoCoH, r-LoCoH, and a-LoCoH respectively (Getz et
  al. 2007).\cr
  The functions \code{LoCoH.k.area}, \code{LoCoH.r.area}, and
  \code{LoCoH.a.area} compute the curve showing the relationships
  between the home-range size (computed to a specified percent) and the
  k, r or a parameters respectively.\cr
}
\usage{
LoCoH.k(xy, k=5, unin = c("m", "km"),
        unout = c("ha", "m2", "km2"),
        duplicates=c("random","remove"), amount = NULL)

LoCoH.r(xy, r, unin = c("m", "km"),
        unout = c("ha", "m2", "km2"),
        duplicates=c("random","remove"), amount = NULL)

LoCoH.a(xy, a, unin = c("m", "km"),
        unout = c("ha", "m2", "km2"),
        duplicates=c("random","remove"), amount = NULL)



LoCoH.k.area(xy, krange, percent=100, unin = c("m", "km"),
             unout = c("ha", "m2", "km2"),
             duplicates=c("random","remove"), amount = NULL)

LoCoH.r.area(xy, rrange, percent=100, unin = c("m", "km"),
             unout = c("ha", "m2", "km2"),
             duplicates=c("random","remove"), amount = NULL)

LoCoH.a.area(xy, arange, percent=100, unin = c("m", "km"),
             unout = c("ha", "m2", "km2"),
             duplicates=c("random","remove"), amount = NULL)

}
\arguments{
  \item{xy}{An object inheriting the class \code{SpatialPoints}
    containing the x and y relocations of the animal.  If \code{xy}
    inherits the class \code{SpatialPointsDataFrame}, it should contain
    only one column (a factor) corresponding to the identity of the
    animals for each relocation.}
  \item{k}{numeric.  The number of nearest neighbors minus one out of
    which to create convex hulls}
  \item{r}{numeric.  The convex hulls are created out of all points
    within r distance from the root points}
 \item{a}{numeric. Create convex hulls from the maximum number of
   nearest neighbors such that the sum of their distances is less than
   or equal to this parameter}
 \item{unin}{the units of the relocations coordinates.  Either
   \code{"m"} for meters or \code{"km"} for kilometers}
 \item{unout}{the units of the output areas.  Either \code{"m2"} for
   square meters, \code{"km2"} for square kilometers or \code{"ha"} for
   hectares}
 \item{duplicates}{a setting to determine how duplicated points are
   handled.  If "\code{random}" the duplicated points are slightly moved
   randomly. If "\code{remove}" the duplicated points are removed.}
 \item{amount}{if \code{duplicates == random}, this parameter controls the
   amount of noise added to the data (see the help page of \code{jitter}
   for additional information on this parameter).}
 \item{krange}{a vector containing the values of k for which the home
   range size is to be estimated.}
 \item{arange}{a vector containing the values of k for which the home
   range size is to be estimated.}
 \item{rrange}{a vector containing the values of k for which the home
   range size is to be estimated.}
 \item{percent}{the percentage level of the home range.  For the
   function \code{plot.LoCoH}, this value could also be the character
   string \code{"all"}, indicating that all the polygons are to be
   displayed.} 
}
\value{
  The functions \code{LoCoH.*} return either objects of class
  \code{SpatialPolygonsDataFrame} (if the relocations of only one
  animals are passed as the \code{xy} argument) or a list of
  \code{SpatialPolygonsDataFrame} (if the relocations of several
  animals are passed as the \code{xy} argument).\cr

  The functions \code{LoCoH.*.area} return invisibly either a vector (if
  the  relocations of only one animals are passed as the \code{xy}
  argument) or a data frame containing the home-range sizes for various
  values of k, r (rows) for the different animals (columns).\cr
}
\note{
 The LoCoH family of methods for locating Utilization Distributions
 consists of three algorithms: Fixed k LoCoH, Fixed r LoCoH, and
 Adaptive LoCoH.  All the algorithms work by constructing a small convex
 hull for each relocation, and then incrementally merging the hulls
 together from smallest to largest into isopleths.  The 10\% isopleth
 contains 10\% of the points and represents a higher utilization than the
 100\% isopleth that contains all the points.
 \cr
 \cr
 \bold{Fixed k LoCoH:} Also known as k-NNCH, Fixed k LoCoH is described
 in Getz and Willmers (2004).  The convex hull for each point is
 constructed from the (k-1) nearest neighbors to that point.  Hulls are
 merged together from smallest to largest based on the area of the
 hull.
 \cr
 \bold{Fixed r LoCoH:} In this case, hulls are created from all points
 within \code{r} distance of the root point.  When merging hulls, the
 hulls are primarily sorted by the value of k generated for each hull
 (the number of points contained in the hull), and secondly by the area
 of the hull.
 \cr
 \bold{Adaptive LoCoH:}  Here, hulls are created out of the maximum
 nearest neighbors such that the sum of the distances from the nearest
 neighbors is less than or equal to d.  Use the same hull sorting as
 Fixed r LoCoH.
 \cr
 Fixed r LoCoH and Adaptive LoCoH are discussed in Getz et al (2007).
 \cr
 All of these algorithms can take a significant amount of time.  Time
 taken increases exponentially with the size of the data set.
}
\author{Clement Calenge \email{clement.calenge@oncfs.gouv.fr}\cr
  with contributions from Scott Fortmann-Roe \email{scottfr@gmail.com}}
\seealso{\code{\link{MCHu}}, \code{\link{getverticeshr}}.}
\references{
Getz, W.M. & Wilmers, C.C. (2004). A local nearest-neighbor convex-hull
construction of home ranges and utilization
distributions. \emph{Ecography}, \bold{27}, 489--505.
\cr
Getz, W.M., Fortmann-Roe, S.B, Lyons, A., Ryan, S., Cross,
P. (2007). LoCoH methods for the construction of home ranges and
utilization distributions. \emph{PLoS ONE}, \bold{2}: 1--11. 
}
\examples{
\dontrun{

## Load the data
data(puechabonsp)

## The relocations:
locs <- puechabonsp$relocs
locsdf <- as.data.frame(locs)
head(locsdf)

## Shows the relocations
plot(locs, col=as.numeric(locsdf[,1]))


## Examinates the changes in home-range size for various values of k
## Be patient! the algorithm can be very long
ar <- LoCoH.k.area(locs[,1], k=c(8:13))

## 12 points seems to be a good choice (rough asymptote for all animals)
## the k-LoCoH method:
nn <- LoCoH.k(locs[,1], k=12)

## Graphical display of the results
plot(nn, border=NA)

## the object nn is a list of objects of class
## SpatialPolygonsDataFrame
length(nn)
names(nn)
class(nn[[1]])

## shows the content of the object for the first animal
as.data.frame(nn[[1]])

## The 95\% home range is the smallest area for which the
## proportion of relocations included is larger or equal
## to 95\% In this case, it is the 22th row of the
## SpatialPolygonsDataFrame.
## The area covered by the home range is for this first animal
## equal to 22.87 ha.

## shows this area:
plot(nn[[1]][11,])


## rasterization of the home ranges:
## use the map of the area:
image(puechabonsp$map)
ras <- MCHu.rast(nn, puechabonsp$map, percent=100)  
opar <- par(mfrow=c(2,2))
lapply(1:4, function(i) { image(ras,i); box()})
par(opar)


## r-LoCoH and a-LoCoH can be applied similarly

}

}
\keyword{spatial}
