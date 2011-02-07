\name{clusthr}
\alias{clusthr}
\title{Estimation of the Home Range by Single-Linkage Cluster Analysis}
\description{
  \code{clusthr} allows the estimation of the home range by
  single-linkage cluster analysis (see details).
}
\usage{
clusthr(xy, unin = c("m", "km"),
        unout = c("ha", "m2", "km2"),
        duplicates=c("random","remove"), amount = NULL)
}
\arguments{
  \item{xy}{an object inheriting the class \code{SpatialPoints}
    containing the x and y relocations of the animal.  If \code{xy}
    inherits the class \code{SpatialPointsDataFrame}, it should contain
    only one column (factor) corresponding to the identity of the
    animals for each relocation.}
  \item{unin}{the units of the relocations coordinates.  Either
    \code{"m"} (default) for
    meters or \code{"km"} for kilometers}
  \item{unout}{the units of the output areas.  Either \code{"m2"} for
    square meters, \code{"km2"} for square kilometers or \code{"ha"} for
    hectares (default)}
  \item{duplicates}{a setting to determine how duplicated points are
    handled.  If "\code{random}" the duplicated points are slightly moved
    randomly. If "\code{remove}" the duplicated points are removed.}
  \item{amount}{if \code{duplicates == random}, this parameter controls the
    amount of noise added to the data (see the help page of \code{jitter}
    for additional information on this parameter).}
}
\details{
  This method estimates home range using the single-linkage cluster
  analysis modified by Kenward et
  al. (2001). The clustering process is described hereafter: the three
  locations with the minimum mean of nearest-neighbour joining distances
  (NNJD) form the first cluster. At each step, two distances are
  computed: (i) the minimum mean NNJD between three locations (which
  corresponds to the next potential cluster) and (ii) the minimum of the
  NNJD between a cluster "c" and the closest location.  If (i) is smaller
  that (ii), another cluster is defined with these three locations.
  If (ii) is smaller than (i), the cluster "c" gains a new location.  If
  this new location belong to another cluster, the two cluster
  fuses. The process stop when all relocations are assigned to the same
  cluster.\cr

  At each step of the clustering process, the proportion of all
  relocations which are assigned to a cluster is computed (so that the
  home range can be defined to enclose a given proportion of the
  relocations at hand, i.e. to an uncomplete process).
  At a given step, the home range is defined as
  the set of minimum convex polygon enclosing the relocations in the
  clusters.\cr

  Note that a given home-range contour can be extracted using the
  function \code{getverticeshr}.
}
\value{
  The function \code{clusthr} returns either objects of class
  \code{SpatialPolygonsDataFrame} (if the relocations of only one
  animals are passed as the \code{xy} argument) or a list of
  \code{SpatialPolygonsDataFrame} of class \code{MCHu} -- Multiple
  Convex Hull (if the relocations of several animals are passed as the
  \code{xy} argument).
}
\references{
  Kenwward R.E., Clarke R.T., Hodder K.H. and Walls S.S. (2001) Density
  and linkage estimators of homre range: nearest neighbor clustering
  defines multinuclear cores. \emph{Ecology}, \bold{82}, 1905--1920.
}
\author{Clement Calenge \email{clement.calenge@oncfs.gouv.fr}}
\seealso{\code{\link{MCHu}} for further information on the class
  \code{MCHu}, and \code{\link[sp]{SpatialPolygonsDataFrame-class}} for
  additional information on this class. See \code{\link{getverticeshr}}
  to extract a given home range contour.}
\examples{

data(puechabonsp)
lo<-puechabonsp$relocs[,1]

## Home Range Estimation
res <- clusthr(lo)

## Displays the home range
plot(res)

## Computes the home range size
MCHu2hrsize(res)

## get the 95 percent home range:
plot(getverticeshr(res, percent=95))


}
\keyword{spatial}
\keyword{hplot}
