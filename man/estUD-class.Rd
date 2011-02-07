\name{estUD-class}
\docType{class}
\alias{estUD-class}
\alias{coerce,estUD,data.frame-method}
\alias{show,estUD-method}

\title{Class "estUD": Storing Utilization Distributions in R}
\description{This class is an extension of the class
  \code{SpatialPixelsDataFrame} of the package sp, and is designed to
  store the utilization distribution of an animal}
\section{Objects from the Class}{
  Objects of class "\code{estUD}" can be created using the functions
  \code{kernelUD} and \code{getvolumeUD}.
}
\section{Slots}{
  \describe{
    \item{\code{h}:}{Object of class \code{"list"} containing all
      information concerning the smoothing parameters used in the
      estimation process}
    \item{\code{vol}:}{Object of class \code{"logical"}
      indicating whether the mapped values coorespond to the UD
      or to the volume under the UD (see \code{?kernelUD})}
    \item{\code{data}:}{Object of class \code{"data.frame"}
      containing the values of the UD}
  }
}
\section{Extends}{
  Class \code{"\linkS4class{SpatialPixelsDataFrame}"}, directly.
}
\section{Methods}{
  \describe{
    \item{coerce}{\code{signature(from = "estUD", to = "data.frame")}:
      converts the object into a data frame}
    \item{show}{\code{signature(object = "estUD")}: printing method of
      the object }
  }
}
\author{ Clement Calenge \email{clement.calenge@oncfs.gouv.fr} }

\seealso{ \code{\linkS4class{SpatialPixelsDataFrame}} for additional
  information about this class, and \code{\link{kernelUD}} for
  additional information about the methods generating such objects.
}
\examples{

## load the data
data(puechabonsp)

## estimate one UD for each animal
jj <- kernelUD(puechabonsp$relocs[,1])
image(jj)
jj


## Consider the first animal
jj[[1]]
class(jj[[1]])
image(jj[[1]])


}
\keyword{classes}
