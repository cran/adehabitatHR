\name{kver2spol}
\alias{kver2spol}
\alias{khr2estUDm}
\title{Conversion of old classes from adehabitat to classes from adehabitatHR}
\description{
  These functions convert home ranges available in adehabitat toward
  classes available in the package \code{adehabitatHR}.\cr
  \cr
  \code{kver2spol} converts an object of class \code{kver} into an
  object of class \code{SpatialPolygons}.
  \cr
  \code{khr2estUDm} converts an object of class \code{khr} (kernel UD)
  into an object of class \code{estUDm}.
  
}
\usage{
kver2spol(kv)
khr2estUDm(x)
}
\arguments{
  \item{kv}{an object of class \code{kver}.}
  \item{x}{an object of class \code{khr}.}
}
\author{ Clement Calenge \email{clement.calenge@oncfs.gouv.fr}}
\keyword{hplot}
