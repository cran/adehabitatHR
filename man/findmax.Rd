\name{findmax}
\alias{findmax}
\title{Find Local Maxima on a Map of Class 'SpatialPixelsDataFrame'}
\description{
  \code{findmax} finds the local maxima on a map of class
  \code{SpatialPixelsDataFrame}.
}
\usage{
  findmax(x)
}
\arguments{
  \item{x}{a map of class \code{SpatialPixelsDataFrame} with one column}
}
\details{
  This function may be useful, among other things, to identify the local
  modes of the utilization distribution of an animal estimated
  using \code{kernelUD}. 
}
\value{
  an object of class \code{SpatialPoints} containing the x and y
  coordinates of the local maxima.
}
\author{ Clement Calenge \email{clement.calenge@oncfs.gouv.fr} }
\seealso{ \code{\link[sp]{SpatialPixelsDataFrame-class}} for additionnal
  information on objects of class \code{SpatialPixelsDataFrame}.}
\examples{

data(puechabonsp)

## estimates the UD
kud <- kernelUD(puechabonsp$relocs[,1])

## displays the maximum
image(kud[[1]])
points(findmax(kud[[1]]))

}
\keyword{spatial}
