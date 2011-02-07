\name{MCHu}
\alias{MCHu}
\alias{print.MCHu}
\alias{plot.MCHu}
\alias{MCHu.rast}
\alias{MCHu2hrsize}
\title{The Class "MCHu": Managing Home Ranges Built by Multiple Convex
  Hulls}
\description{
  The class \code{"MCHu"} is designed to store home ranges built by
  multiple convex hulls, for example built using the single-linkage
  cluster algorithm (function \code{clusterhr}) or the LoCoH
  (e.g. function \code{LoCoH.k}).\cr
  The function \code{plot.MCHu} allows to graphically display the
  home-ranges.\cr
  \code{MCHu.rast} allows to compute a raster map of the home ranges.\cr
  \code{MCHu2hrsize} allows to compute the home range size for specified
  percentage levels for the home range (see \code{help(plot.hrsize)}).  
}
\usage{
\method{print}{MCHu}(x, \dots)

\method{plot}{MCHu}(x, percent="all", points=NULL, \dots)

MCHu.rast(x, spdf, percent=100)

MCHu2hrsize(x, percent=seq(20,100, by=10), plotit=TRUE)
}
\arguments{

  \item{x}{an object of class \code{MCHu}.}
  \item{spdf}{an object of class \code{SpatialPixelsDataFrame}.}
  \item{points}{an object of class \code{SpatialPoints} or
   \code{SpatialPointsDataFrame} with one column (a factor storing the
   identity of the animal for each relocation), containing the
   relocations of the animal(s).}
 \item{percent}{the percentage level of the home range.  For the
   function \code{plot.MCHu}, this value could also be the character
   string \code{"all"}, indicating that all the polygons are to be
   displayed.} 
 \item{plotit}{a logical value indicating whether the results should be
   plotted.}
 \item{\dots}{additional arguments to be passed to the functions
   \code{print} and \code{plot}.}

}
\details{
  The class \code{"MCHu"} is basically a list of objects of class
  \code{SpatialPolygonsDataFrame}, with one data frame per animal. 
}
\value{
  The function \code{MCHu.rast} returns an object of class
  \code{SpatialPixelsDataFrame}.\cr
  The function \code{MCHu2hrsize} returns an object of class
  \code{hrsize} (see \code{?mcp.area}).
}
\author{Clement Calenge \email{clement.calenge@oncfs.gouv.fr}}
\seealso{\code{\link{clusthr}} and \code{\link{LoCoH}} for home range
  estimation methods returning this class of objects.}
\examples{

data(puechabonsp)

## The relocations:
locs <- puechabonsp$relocs
locsdf <- as.data.frame(locs)
head(locsdf)

## Shows the relocations
plot(locs, col=as.numeric(locsdf[,1]))


## 12 points seems to be a good choice (rough asymptote for all animals)
## the k-LoCoH method:
nn <- LoCoH.k(locs[,1], k=12)

## Graphical display of the results
plot(nn, border=NA)

## Rasterize the home range on the elevation map:
image(puechabonsp$map)
(oo <- MCHu.rast(nn, puechabonsp$map))
image(oo)


}
\keyword{spatial}
\keyword{hplot}
