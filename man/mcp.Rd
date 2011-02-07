\name{mcp}
\alias{mcp}
\alias{hr.rast}
\alias{plot.hrsize}
\alias{mcp.area}
\title{Estimation of the Home Range Using the Minimum Convex Polygon
  Estimator}
\description{
  \code{mcp} computes the home range of several
  animals using the Minimum Convex Polygon estimator.\cr
  \code{mcp.area} is used for home-range size estimation.\cr
  \code{hr.rast} is used to rasterize a minimum convex polygon.\cr
  \code{plot.hrsize} is used to display the home-range size estimated at
  various levels.
}
\usage{
mcp(xy, percent=95, unin = c("m", "km"),
    unout = c("ha", "km2", "m2"))

mcp.area(xy, percent = seq(20,100, by = 5),
         unin = c("m", "km"),
         unout = c("ha", "km2", "m2"), plotit = TRUE)

hr.rast(mcp, w)

\method{plot}{hrsize}(x, \dots)
}
\arguments{
  \item{xy}{An object inheriting the class \code{SpatialPoints}
    containing the x and y relocations of the animal.  If \code{xy}
    inherits the class \code{SpatialPointsDataFrame}, it should contain
    only one column (a factor) corresponding to the identity of the
    animals for each relocation.}
  \item{percent}{A single number for the function \code{mcp} and a
    vector for the function \code{mcp.area}:  100 minus the proportion
    of outliers to be excluded from the computation.}
  \item{unin}{the units of the relocations coordinates.  Either
    \code{"m"} (default) for meters or \code{"km"} for kilometers}
  \item{unout}{the units of the output areas.  Either \code{"m2"} for
    square meters, \code{"km2"} for square kilometers or \code{"ha"} for
    hectares (default)}
  \item{plotit}{logical.  Whether the plot should be drawn.}
  \item{x}{an objet of class \code{hrsize} returned by the function
    \code{mcp.area}, or \code{kernel.area} (see \code{kernelUD()})}
  \item{mcp}{an objet of class \code{SpatialPolygons} returned by the
    function \code{mcp}.}
  \item{w}{an objet of class \code{SpatialPixelsDataFrame} used as a
    reference for the rasterization.}
  \item{\dots}{additional arguments to be passed to the function
    \code{plot}.} 
}
\details{
  This function computes the Minimum Convex Polygon estimation after the
  removal of \code{(100 minus percent)} percent of the relocations the
  farthest away from the centroid of the home range (computed by the
  arithmetic mean of the coordinates of the relocations for each
  animal).\cr
}
\value{
  \code{mcp} returns an object of class \code{SpatialPolygonsDataFrame},
  in which the first column contains the ID of the animals, and the
  second contains the home range size.\cr
  \code{mcp.area} returns a data frame of class \code{hrsize}, 
  with one column per animal and one row per level of
  estimation of the home range.\cr
  \code{hr.rast} returns an object of class
  \code{SpatialPixelsDataFrame}.
}
\references{
  Mohr, C.O. (1947) Table of equivalent populations of north
  american small mammals. \emph{The American Midland Naturalist},
  \bold{37}, 223-249.}
\author{ Clement Calenge \email{clement.calenge@oncfs.gouv.fr} }

\seealso{ \code{\link[grDevices]{chull}},
  \code{\link[sp]{SpatialPolygonsDataFrame-class}} for additionnal
  information on the class \code{SpatialPolygonsDataFrame}.}

\examples{

data(puechabonsp)
rel <-  puechabonsp$relocs


## estimates the MCP
cp <- mcp(rel[,1])

## The home-range size
as.data.frame(cp)

## Plot the home ranges
plot(cp)

## ... And the relocations
plot(rel, col=as.data.frame(rel)[,1], add=TRUE)


## Computation of the home-range size:
cuicui1 <- mcp.area(rel[,1])

## Rasterization
ii <- hr.rast(cp, puechabonsp$map)

opar <- par(mfrow=c(2,2))
lapply(1:4, function(i) {image(ii, i); box()})
par(opar)


}
\keyword{hplot}
\keyword{spatial}
