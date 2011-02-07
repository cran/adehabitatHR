\name{CharHull}
\alias{CharHull}
\title{
  Estimation of the Home Range by Delaunay Triangulation method
}
\description{
  The function \code{CharHull} implements the method developed by Downs
  and Horner (2009) for the home range estimation.
}
\usage{
CharHull(xy, unin = c("m", "km"),
         unout = c("ha", "m2", "km2"),
         duplicates = c("random", "remove"), amount = NULL)
}
\arguments{
  \item{xy}{an object inheriting the class \code{SpatialPoints}
    containing the x and y coordinates of the relocations of the animal.
    If \code{xy} inherits the class \code{SpatialPointsDataFrame}, it
    should contain only one column (factor) corresponding to the
    identity of the animals for each relocation.}
  \item{unin}{the units of the relocations coordinates.  Either "m"
    (default) for meters or "km" for kilometers}
  \item{unout}{the units of the output areas.  Either "m2" for square
    meters, "km2" for square kilometers or "ha" for hectares
    (default)}
  \item{duplicates}{a setting to determine how duplicated points are handled.
    If "random" the duplicated points are slightly moved
    randomly. If "remove" the duplicated points are removed.}
  \item{amount}{if \code{duplicates == random}, this parameter controls
    the amount of noise added to the data (see the help page of
    \code{jitter} for additional information on this parameter).}
}
\details{
  This method consists in the computation of the Delaunay triangulation
  of the set of relocations.  Then, the triangles are ordered from the
  smallest to the largest. It is possible to select a given percentage
  of the smallest triangles (measured by their area) as the home-range
  estimation.  The contour can be extracted with the function
  \code{getverticeshr}
}
\value{
  an object of the class \code{MCHu}
}
\references{
  Downs J.A. and Horner, M.W. (2009) A Characteristic-Hull Based Method
  for  Home Range Estimation. \emph{Transactions in GIS}, \bold{13},
  527--537. 
}
\author{
Clement Calenge \email{clement.calenge@oncfs.gouv.fr}
}
\note{
  This function relies on the function \code{tri.mesh} and
  \code{triangles} from the package \code{tripack}.
}
\seealso{
  \code{\link{MCHu}} for further information on the class
  \code{MCHu}, and \code{\link[sp]{SpatialPolygonsDataFrame-class}} for
  additional information on this class. See \code{\link{getverticeshr}}
  to extract a given home range contour.
}
\examples{
data(puechabonsp)
lo<-puechabonsp$relocs[,1]
     
## Home Range Estimation
res <- CharHull(lo)

## Displays the home range
plot(res)
 
## Computes the home range size
MCHu2hrsize(res)

## Computes the 95 percent home range
ver <- getverticeshr(res)
ver
plot(ver)

}
\keyword{spatial}
