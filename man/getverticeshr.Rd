\name{getverticeshr}
\alias{getverticeshr}
\alias{getverticeshr.estUD}
\alias{getverticeshr.estUDm}
\alias{getverticeshr.MCHu}
\alias{getverticeshr.default}
\title{
  Extract the home-range contour of one or several animals
}
\description{
  These functions allow the extraction of the home-range contours
  computed using various methods (kernel home range, cluster home range,
  etc.)
}
\usage{
getverticeshr(x, percent = 95, \dots)
\method{getverticeshr}{estUD}(x, percent = 95, ida = NULL, unin = c("m", "km"),
                              unout = c("ha", "km2", "m2"),
                              standardize = FALSE, \dots)
\method{getverticeshr}{estUDm}(x, percent = 95, whi = names(x),
                               unin = c("m", "km"),
                               unout = c("ha", "km2", "m2"),
                               standardize = FALSE, \dots)
\method{getverticeshr}{MCHu}(x, percent = 95, whi = names(x), \dots)
\method{getverticeshr}{default}(x, percent = 95, \dots)
}
\arguments{
  \item{x}{For \code{getverticeshr.estUD}, an object of class
    \code{estUD}.  For \code{getverticeshr.estUDm}, an object of class
    \code{estUDm}.  For \code{getverticeshr.MCHu}, an object of class
    \code{MCHu}.}
  \item{percent}{a single value giving the percentage level for
    home-range estimation}
  \item{ida}{a character string indicating the id of the polygons
    corresponding to the home range in the resulting
    \code{SpatialPolygonsDataFrame} (see the help page of
    \code{SpatialPolygonsDataFrame}).  By default it is set to
    \code{"homerange"}}
  \item{unin}{the units of the relocations coordinates.  Either "m" for
    meters (default) or "km" for kilometers}
  \item{unout}{the units of the output areas.  Either "m2" for square
    meters, "km2" for square kilometers or "ha" for hectares
    (default)}
  \item{whi}{a vector of character strings indicating which animals
    should be returned.}
  \item{standardize}{a logical value indicating whether the UD should be
    standardized over the area of interest, so that the volume
    under the UD and *over the area* is equal to 1..}
  \item{\dots}{
    Additional arguments to be passed to and from other methods
  }
}
\value{
  An object of class \code{SpatialPolygonsDataFrame} containing the
  selected home range contours of the animals.
}
\author{
Clement Calenge \email{clement.calenge@oncfs.gouv.fr}
}
\note{
  The function \code{getverticeshr.default} is present for compatibility
  purposes.  Its use generates an error.
}
\seealso{
  \code{\link{kernelUD}}, \code{\link{kernelbb}} or
  \code{\link{kernelkc}} for methods generating objects of classes
  \code{estUD} and \code{estUDm}, 
  \code{\link{clusthr}}, \code{\link{LoCoH.a}} and
  \code{\link{CharHull}} for methods generating objects of class
  \code{\link{MCHu}}. 
}
\examples{

### Example with a kernel home range
data(puechabonsp)
loc <- puechabonsp$relocs
     
## have a look at the data
head(as.data.frame(loc))
## the first column of this data frame is the ID
     
## Estimation of UD for the four animals
(ud <- kernelUD(loc[,1]))

## Calculates the home range contour
ver <- getverticeshr(ud, percent=95)
ver
plot(ver)

## Example with a cluster home range
clu <- clusthr(loc[,1])
ver2 <- getverticeshr(clu, percent=95)
ver2
plot(ver2)

}
\keyword{spatial}
