\name{kerneloverlap}
\alias{kerneloverlap}
\title{Spatial Interaction between Animals Monitored Using Radio-Tracking}
\description{
  This function implements all the indices of kernel home-range overlap
  reviewed by Fieberg and Kochanny (2005).
}
\usage{
kerneloverlap(xy, method = c("HR", "PHR", "VI", "BA", "UDOI",
              "HD"), percent = 95, conditional = FALSE, \dots)
}
\arguments{
  \item{xy}{an object of class \code{SpatialPointsDataFrame} containing
    only one column (which is a factor indicating the identity
    associated to the relocations))}
  \item{method}{the desired method for the estimation of overlap (see
    details)}
  \item{percent}{the percentage level of the home range estimation}
  \item{conditional}{logical.  If \code{TRUE}, the function sets to 0
    the pixels of the grid over which the UD is estimated,
    outside the home range of the animal estimated at a level of
    probability equal to \code{percent}.  Note that this argument has no
    effect when \code{meth="HR"}.}
  \item{\dots}{additional arguments to be passed to the function
    \code{kernelUD} for the kernel estimation of the utilization
    distribution.}
}
\details{
  Fieberg and Kochanny (2005) made an extensive review of the indices of
  overlap between utilization distributions (UD) of two animals.  The
  function \code{kerneloverlap} implements these indices.  The argument
  \code{method} allows to choose an index\cr

  The choice \code{method="HR"} computes the proportion of the home
  range of one animal covered by the home range of another one, i.e.:
  \deqn{HR_{i,j} = A_{i,j} / A_i}{HR[i,j] = A[i,j] / A[i]},
  where \eqn{A_{i,j}}{A[i,j]} is the area of the intersection between
  the two home ranges and \eqn{A_i}{A[i]} is the area of the home range
  of the animal i.\cr

  The choice \code{method="PHR"} computes the volume under the UD of the
  animal j that is inside the home range of the animal i (i.e., the
  probability to find the animal j in the home range of i). That is:
  \deqn{PHR_{i,j} = \int \int_{A_i} UD_j(x,y) dxdy}{PHR[i,j] =
  integral_over_x integral_over_y A[i] UD[j](x,y) dxdy} where
  \eqn{UD_j(x,y)}{UD[j](x,y)} is the value of the utilization
  distribution of the animal j at the point x,y.\cr

  The choice \code{method="VI"} computes the volume of the intersection
  between the two UD, i.e.:
  \deqn{VI = \int_x \int_y min(UD_i(x,y),UD_j(x,y)) dxdy}{VI =
  integral_over_x integral_over_y min(UD[i](x,y),UD[j](x,y)) dxdy}\cr

  Other choices rely on the computation of the joint distribution of the
  two animals under the hypothesis of independence UD[i](x,y) *
  UD[j](x,y). 
  
  The choice \code{method="BA"} computes the Bhattacharyya's affinity
  \deqn{BA = \int_x \int_y \sqrt{UD_i(x,y)} \times \sqrt{UD_j(x,y)}}{BA =
  integral_over_x integral_over_y sqrt(UD[i](x,y)) * sqrt(UD[j](x,y))}\cr

  The choice \code{method="UDOI"} computes a measure similar to the
  Hurlbert index of niche overlap:
  \deqn{UDOI = A_{i,j} \int_x \int_y UD_i(x,y) \times
  UD_j(x,y)}{UDOI = A[i,j] integral_over_x integral_over_y UD[i](x,y) *
  UD[j](x,y)}\cr
  
  The choice \code{method="HD"} computes the Hellinger's distance:
  \deqn{HD = \int_x \int_y ((\sqrt UD_i(x,y) - \sqrt UD_j(x,y))^2
  dxdy)^{1/2}}{HD = integral_over_x integral_over_y ((sqrt(UD[i](x,y) -
  sqrt(UD[j](x,y)))^2 dxdy)^(1/2))}\cr 

}
\value{
  A matrix giving the value of indices of overlap for all pairs of
  animals.
}
\references{
Fieberg, J. and Kochanny, C.O. (2005) Quantifying home-range overlap: the
importance of the utilization distribution. \emph{Journal of Wildlife
  Management}, \bold{69}, 1346--1359.
}
\author{Clement Calenge \email{clement.calenge@oncfs.gouv.fr},
  based on a work of John Fieberg}
\seealso{\code{\link{kernelUD}} for additional information on kernel
  estimation of home ranges}
\examples{

data(puechabonsp)

kerneloverlap(puechabonsp$relocs[,1],
              grid=200, meth="VI", conditional=TRUE)

kerneloverlap(puechabonsp$relocs[,1],
              grid=200, meth="HR")

kerneloverlap(puechabonsp$relocs[,1],
              grid=200, meth="PHR")

kerneloverlap(puechabonsp$relocs[,1],
              grid=200, meth="BA")

kerneloverlap(puechabonsp$relocs[,1],
              grid=200, meth="UDOI")

kerneloverlap(puechabonsp$relocs[,1],
              grid=200, meth="HD")

}
\keyword{spatial}
