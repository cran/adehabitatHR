\name{BRB}
\alias{BRB}
\alias{BRB.D}
\alias{BRB.likD}
\title{
  Utilization Distribution of an Animal Based on Biased Random Bridges
}
\description{
  This function estimates the utilization distribution of one/several
  animals using a biased random bridge approach (Benhamou and Cornelis
  2010, Benhamou 2011).  This function also allows the decomposition of
  the utilization distribution into (i) an intensity distribution
  reflecting the average time spent by the animal in the habitat
  patches, and (ii) a recursion distribution reflecting the number of
  visits of the animal in the habitat patches (Benhamou and
  Riotte-Lambert, 2012).
}
\usage{
BRB(ltr, D, Tmax, Lmin, hmin, type=c("UD","ID", "RD"), radius = NULL,
    maxt = NULL, filtershort = TRUE, habitat = NULL, activity = NULL,
    grid = 200, b=FALSE, same4all=FALSE, extent=0.5, tau = NULL,
    boundary=NULL) 

BRB.D(ltr, Tmax = NULL, Lmin = NULL, habitat = NULL, activity = NULL)

BRB.likD(ltr, Dr=c(0.1,100),
         Tmax = NULL, Lmin = NULL,
         habitat = NULL, activity = NULL)

}
\arguments{
  \item{ltr}{
    an object of class \code{ltraj}.
  }
  \item{D}{
    a number corresponding to the diffusion parameter (in squared
    "units" per second, where "units" denote the units of the relocation
    coordinates) used for the estimation.  Alternatively this parameter
    may be an object of class \code{DBRB} returned by the function
    \code{BRB.D}.
  }
  \item{Dr}{
    a vector of length two giving the lower and upper limits of the
    diffusion coefficient, within which the maximum likelihood could be
    found.
  }
  \item{Tmax}{
    the maximum duration (in seconds) allowed for a step built by successive
    relocations.  All steps characterized by a duration dt greater than
    \code{Tmax} are not taken into account in the calculations.
  }
  \item{Lmin}{
    the minimum distance (in units of the coordinates) between successive 
    relocations, defining intensive use or resting (See details).  The
    distance should be specified in units of the relocation coordinates
    (i.e. in metres if they are specified in metres).
  }
  \item{hmin}{
    The minimum smoothing parameter (in units of the relocations
    coordinates), applied to all recorded relocations. See details for a
    description of this parameter.
  }
  \item{type}{
    The type of distribution expected by the user: \code{"UD"} returns
    the utilization distribution using the approach described by
    Benhamou and Cornelis (2010).  \code{"ID"} returns the intensity
    distribution described in Benhamou and Riotte-Lambert (2012), i.e. a
    distribution reflecting the average time spent in habitat patches.
    \code{"RD"} returns the recursion distribution described in Benhamou
    and Riotte-Lambert (2012), i.e. a distribution reflecting the number
    of visits in habitat patches.  
  }
  \item{radius}{
    If \code{type = "ID"} or \code{"RD"}, the radius of the patches (in
    units of the relocation coordinates) used in the calculation of the
    residence time or the number of visits.  If \code{NULL}, the radius
    is set to \code{3*hmin}.
  }
  \item{maxt}{
    If \code{type = "ID"} or \code{"RD"}, maximum time threshold (in
    seconds) that the animal is allowed to spend outside the patch
    before that we consider that the animal actually left the patch (see
    \code{?residenceTime}).
  }
  \item{filtershort}{
    logical indicating the behaviour of the function
    when the length of a step is lower than \code{Lmin} (see details).
    It must be set to TRUE if track segments shorter than
    Lmin are assumed to correspond to resting periods, which thereby
    will be filtered out systematically, or to FALSE to take these short
    track segments into account when not associated to resting (animals
    active without moving more than Lmin). They then will be given a
    null diffusion coefficient.
  }
  \item{habitat}{
    optionally, an object of class \code{SpatialPixelsDataFrame} with
    one column describing the habitat type on the area.
  }
  \item{activity}{
    optionally, a character indicating the name of the variable in the
    infolocs component of \code{ltr} indicating the proportion of time
    between relocation i and relocation i+1 during which the animal was
    active (Users of adehabitatLT prior to version 0.3 should read the
    section Warning below).
  }
  \item{grid}{
    a number giving the size of the grid on which the UD should
    be estimated.  Alternatively, this parameter may be an object
    of class \code{SpatialPixels}, or a list of objects of class
    \code{SpatialPixels} with as many elements as there are bursts in
    \code{ltr}.
  }
  \item{b}{
    logical specifying how relocation and movement variances are
    combined.  If \code{TRUE}, the relocation variance progressively
    merges with the movement component; if \code{FALSE}, the relocation
    variance has a constant weight (see Benhamou, 2011).
  }
  \item{same4all}{
    logical.  if \code{TRUE}, the same grid is used for the estimation
    of all bursts.  If \code{FALSE}, one grid is used per burst.
  }
  \item{extent}{
    a value indicating the extent of the grid used for the
    estimation (the extent of the grid on the abscissa is equal
    to \code{(min(xy[,1]) + extent * diff(range(xy[,1])))}).
  }
  \item{tau}{
    interpolation time (tau, in seconds).  Defaults to \code{tmin/10},
    where \code{tmin} is the minimum duration of a step in \code{ltr}. 
  }
  \item{boundary}{
    If, not \code{NULL}, an object inheriting the class
    \code{SpatialLines} defining a barrier that cannot be crossed by the
    animals. There are constraints on the shape of the barrier that
    depend on the smoothing parameter (***see details***)
  }
}
\details{

  The function \code{BRB} uses the biased random bridge approach 
  to estimate the Utilization Distribution of an animal with serial
  autocorrelation of the relocations.  This approach is similar to the
  Brownian bridge approach (see \code{?kernelbb}), with several
  noticeable improvements.  Actually, the Brownian bridge approach
  supposes that the animal movement is random and purely diffusive
  between two successive relocations: it is supposed that the animal
  moves in a purely random fashion from the starting relocation and
  reaches the next relocation randomly. The BRB approach goes further by
  adding an advection component (i.e., a "drift") to the purely
  diffusive movement: is is supposed that the animal movement is
  governed by a drift component (a general tendency to move in the
  direction of the next relocation) and a diffusion component (tendency
  to move in other directions than the direction of the drift).
  
  The BRB approach is based on the biased random walk model.  This model
  is the following: at a given time t, the speed of the animal is drawn
  from a probability density function (pdf) and the angle between the
  step and the east direction is drawn from a circular pdf with given
  mean angle and concentration parameters.  A biased random walk occurs
  when this angular distribution is not uniform (i.e. there is a
  preferred direction of movement).  Now, consider two successive
  relocations r1 = (x1, y1) and r2 = (x2, y2) collected respectively at times t1
  and t2.  The aim of the Biased Random Bridges approach is to estimate
  the pdf that the animal is located at a given place r = (x,y) at time
  ti (with t1 < ti < t2), given that it is located at r1,r2 at times
  t1,t2, and given that the animal moves according a biased random walk
  with an advection component determined by r1 and r2.
  
  Benhamou (2011) proposed an approximation for this pdf, noting
  that it can be approximated by a circular bivariate normal
  distribution with mean location corresponding to  (x1 + pi*(x2-x1), y1
  + pi*(y2-y1)), where pi = (ti-t1)/(t2-t1).  The variance-covariance
  matrix of this distribution is diagonal, with both diagonal elements
  corresponding to the diffusion coefficient D.  This coefficient D can
  be estimated using the plug-in method, using the function
  \code{BRB.D} (for details, see Benhamou, 2011).  Note that the
  diffusion parameter \code{D} can be estimated for each habitat type is
  a habitat map is available.  Note that the function \code{BRB.likD}
  can be used alternatively to estimate the diffusion coefficient using
  the maximum likelihood method.

  An important aspect of the BRB approach is that the drift component is
  allowed to change in direction and strength from one step to the
  other, but should remain constant during each of them.  For this
  reason, it is required to set an upper time threshold \code{Tmax}.
  Steps characterized by a longer duration are not taken into account
  into the estimation of the pdf. This upper threshold should be based
  on biological grounds.

  As for the Brownian bridge approach, this conditional pdf based on
  biased random walks takes an infinite value at times ti = t1 and ti =
  t2 (because, at these times, the relocation of the animal is known
  exactly).  Benhamou proposed to circumvent this drawback by
  considering that the true relocation of the animal at times t1 and t2
  is not known exactly.  He noted: "a GPS fix should be
  considered a punctual sample of the possible locations at which the
  animal may be observed at that time, given its current motivational
  state and history. Even if the recording noise is low, the relocation
  variance should therefore be large enough to encompass potential
  locations occurring in the same habitat patch as the recorded
  location". He proposed two ways to include this "relocation
  uncertainty" component in the pdf: (i) either the relocation variance
  progressively merges with the movement component, (ii) or the
  relocation variance has a constant weight. This is controlled by the
  parameter \code{b} of the function. In both cases, the minimum
  uncertainty over the relocation of an animal is observed for ti = t1
  or t2.  This minimum standard deviation corresponds to the parameter
  \code{hmin}. According to Benhamou and Cornelis, "\code{hmin} must be
  at least equal to the standard deviation of the localization errors
  and also must integrate uncertainty of the habitat map when UDs are
  computed for habitat preference analyses. Beyond these technical
  constraints, \code{hmin} also should incorporate a random component
  inherent to animal behavior because any recorded location, even if
  accurately recorded and plotted on a reliable map, is just a punctual
  sample of possible locations at which the animal may be found at that
  time, given its current motivational state and history. Consequently,
  hmin should be large enough to encompass potential locations occurring
  in the same habitat patch as the recorded location".

  Practically, the BRB approach can be carried out with the help of the
  movement-based kernel density estimation (MKDE) developed by Benhamou
  and Cornelis (2010).  This method consists in dividing each step i
  into Ti/tau intervals, where Ti is the duration of the step (in
  seconds) and \code{tau} is the interpolation time (in seconds).  A
  kernel density estimation is then used to estimate the required pdf,
  with a smoothing parameter varying with each interpolated location ri
  and corresponding to: hi^2 = hmin^2 + 4pi(1-pi)(hmax^2 -
  hmin^2)Ti/Tmax.  In this equation, \code{hmax^2} corresponds to
  \code{hmin^2+ D*Tmax/2} if \code{b} is FALSE and to \code{D*Tmax/2}
  otherwise.  Note that this smoothing parameter may be a 
  function of the habitat type where the interpolated relocation occurs
  if the diffusion parameters are available for each habitat types.

  The special case where a given step covers a distance lower than
  \code{Lmin} merits further details.  When the parameter
  \code{filtershort = TRUE}, it is always assumed that the animal was
  resting at this time, and this step is filtered out before the
  estimation.  When the parameter \code{filtershort = FALSE}, this
  assumption is not made.  In this case, the behaviour of the function
  depends on the availability of a variable measuring the activity of
  the animal (when the name of the variable containing the ai in the
  \code{infolocs} component is passed as the parameter \code{activity};
  see \code{?infolocs} for additionnal information on this component).
  If the animal was active during the step, the smoothing parameter hi
  is set to \code{hmin} for this step. This procedure allows to give
  more weight to the immediate surroundings of this relocation
  (indicating an intensive use of these immediate surroundings).  If the
  animal was inactive, then the animal was resting and the step is
  filtered out before the estimation.  Note however that activity value
  may sometimes be relatively high while the animal is resting, e.g. if
  disturbed by flies, possibly requiring manual correction of activity
  values based on the distance moved between relocations).

  If no activity variable is available and \code{filtershort = FALSE},
  it is always suppose that the animal was active between the two
  relocations, the step is not filtered out and the
  smoothing parameter hi is set to \code{hmin} for this step.
  
  The parameter \code{boundary} allows to define a barrier that cannot
  be crossed by the animals.  When this parameter is
  set, the method described by Benhamou and Cornelis (2010) for
  correcting boundary biases is used.  The boundary can possibly be
  defined by several nonconnected lines, each one being built by several
  connected segments. Note that there are constraints on these segments
  (not all kinds of boundary can be defined): (i) each segment length
  should at least be equal to \code{3*h} (the size of "internal lane"
  according to the terminology of Benhamou and Cornelis), (ii) the angle
  between two line segments should be greater that \code{pi/2} or lower
  that \code{-pi/2}.  The UD of all the pixels located within a
  band defined by the boundary and with a width equal to \code{6*h}
  ("external lane") is set to zero.

  Benhamou and Riotte-Lambert (2012) showed that the space use at any
  given location, as estimated by this approach, can be seen as the
  product between the mean residence time per visit times the number of
  visits of the location.  They proposed an approach allowing the
  decomposition of the UD into two components: (i) the intensity
  distribution reflecting this average residence time and (ii) a
  recursion distribution reflecting the number of visits.  This function
  allows to estimate these two components, by setting the argument
  \code{type} to \code{"ID"} and \code{"RD"} respectively.  
  
  Note that all the methods available to deal with objects of class
  \code{estUDm} are available to deal with the results of the function
  \code{BRB} (see \code{?kernelUD}).
  
}
\value{
  \code{BRB} returns an object of class \code{estUDm} when the UD is
  estimated for several animals, and \code{estUD} when only one animal
  is studied.
  
  \code{BRB.D} and \code{BRB.likD} returns a list of class \code{DBRB},
  with one component per burst containing a data frame with the
  diffusion parameters.
}
\references{
  Benhamou, S. (2011) Dynamic approach to space and habitat use
  based on biased random bridges \emph{PLOS One}, \bold{6}, 1--8.


  Benhamou, S. and Cornelis, D. (2010) Incorporating Movement Behavior
  and Barriers to Improve Biological Relevance of Kernel Home Range
  Space Use Estimates. \emph{Journal of Wildlife Management}, \bold{74},
  1353--1360.

  Benhamou, S. and Riotte-Lambert, L. (2012) Beyond the Utilization
  Distribution: Identifying home range areas that are intensively
  exploited or repeatedly visited. \emph{Ecological Modelling},
  \bold{227},  112--116.
}
\author{
  Clement Calenge \email{clement.calenge@ofb.gouv.fr}, based on a C
  translation of the Pascal source code of the program provided by Simon
  Benhamou.
}
\section{Warning}{
  Users of the version 0.2 of adehabitatHR should be careful that there
  was a slight inconsistency in the package design: whereas all the
  parameters characterizing the steps in an object of class
  \code{"ltraj"} (e.g. \code{dist, dx, dy}) describe the step between
  relocations i and i+1, it was expected for \code{BRB} that the
  activity described the proportion of activity time between relocation
  i-1 and i.  This inconsistency has now been corrected since version
  0.3.
}
\seealso{
  \code{\link{kernelbb}} for the Brownian bridge kernel estimation,
    \code{\link{kernelUD}} and \code{\link{estUD-class}} for additional
    information about objects of class \code{estUDm} and \code{estUD},
    \code{\link[adehabitatLT]{infolocs}} for additional information about the
    \code{infolocs} component, \code{\link[adehabitatLT]{as.ltraj}} for additional
    information about the class \code{ltraj}.
}
\examples{

## Example dataset used by Benhamou (2011)
data(buffalo)

## The trajectory:
buffalo$traj

## The habitat map:
buffalo$habitat

## Show the dataset
plot(buffalo$traj, spixdf = buffalo$habitat)

## Estimate the diffusion component for each habitat type
## Using the plug-in method
vv <- BRB.D(buffalo$traj, Tmax = 180*60, Lmin = 50,
            habitat = buffalo$habitat, activ = "act")

vv

## Note that the values are given here as m^2/s, whereas
## they are given as m^2/min in Benhamou (2011). The
## values in m^2 per min are:
vv[[1]][,2]*60

## Approximately the same values, with slight differences due to
## differences in the way the program of Benhamou (2011) and the present
## one deal with the relocations occurring on the boundary between two
## different habitat types
## Note that an alternative estimation of the Diffusion coefficient
## could be found using maximum likelihood
vv2 <- BRB.likD(buffalo$traj, Tmax = 180*60, Lmin = 50,
                habitat = buffalo$habitat, activ = "act")
vv2
vv[[1]][,2]*60


## Estimation of the UD with the same parameters as those chosen by
## Benhamou (2011)
ud <- BRB(buffalo$traj, D = vv, Tmax = 180*60, tau = 300, Lmin = 50, hmin=100,
          habitat = buffalo$habitat, activity = "act", grid = 50, b=0,
          same4all=FALSE, extent=0.5)
ud

## Show the UD.
image(ud)

\dontrun{
## Example of the decomposition of the UD into a recursion distribution
## and a intensity distribution (Benhamou and Riotte-Lambert 2012).
##
## 1. Intensity Distribution using the same parameters as Benhamou and
## Riotte-Lambert (2012)

id <- BRB(buffalo$traj, D = 440/60, Tmax = 3*3600, Lmin = 50, type = "ID",
          hmin=100, radius = 300, maxt = 2*3600, activity="act", filtershort=FALSE,
          grid = 200, extent=0.1)

rd <- BRB(buffalo$traj, D = 440/60, Tmax = 3*3600, Lmin = 50, type = "RD",
          hmin=100, radius = 300, maxt = 2*3600, activity="act", filtershort=FALSE,
          grid = 200, extent=0.1)

ud <- BRB(buffalo$traj, D = 440/60, Tmax = 3*3600, Lmin = 50, 
          hmin=100, radius = 300, maxt = 2*3600, activity="act", filtershort=FALSE,
          grid = 200, extent=0.1)

par(mfrow = c(2,2), mar=c(0,0,2,0))
image(getvolumeUD(id))
title("ID")
image(getvolumeUD(rd))
title("RD")
image(getvolumeUD(ud))
title("UD")
}




}
\keyword{spatial}
