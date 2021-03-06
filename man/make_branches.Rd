% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/make_branches.R
\name{make_branches}
\alias{make_branches}
\title{Make point clouds that represent branches}
\usage{
make_branches(
  ttops,
  crown_length_fun = function(z) -1.3034187 + z * 0.6937502,
  crown_diameter_fun = function(z) 0.8787043 + z * 0.252169,
  branchs_per_whorl = 5,
  branch_angle = 80,
  interval = 1,
  branch_diameter_fun = function(branch_length) 0.002942 + branch_length * 0.0103,
  spacing = c(0.05, 0.3),
  parallel = TRUE
)
}
\arguments{
\item{ttops}{output from \code{link[lidR]{find_trees}}.}

\item{crown_length_fun}{function. It should take
\emph{height} in m and return \emph{crown lenght} in m.}

\item{crown_diameter_fun}{function. It should take
\emph{height} in m and return \emph{crown diameter} in m.}

\item{branchs_per_whorl}{numeric vector. Single value.}

\item{branch_angle}{numeric. Single value. Angle in degrees.}

\item{interval}{numeric vector. single value. It is
the mean length between successive whorls.}

\item{branch_diameter_fun}{function. It should take
\emph{branch length} in m and return
\emph{branch diameter} in m.}

\item{spacing}{numeric vector. Range that control spacing between points.
The first number are for points near the center and the last for the
farthest points.}

\item{parallel}{logical. Allows parallel processing.}
}
\value{
\code{data.frame} with columns x, y, and z.
}
\description{
Make point clouds that represent branches
}
\details{
\code{crown_length_fun} and \code{crown_diameter_fun} were
fitted by the author of this package as part of a work in progress.
\code{branch_diameter_fun} is from Fernandez and Norero (2006).
}
\examples{
require(lidR)
path <- system.file("external/745_r75m_d1m_local.las", package = "hpsynt")
las <- readLAS(path)

chm <- grid_canopy(las, 1, p2r(0.05, knnidw()))
ttops <- find_trees(las, lmf(3.7, shape = "circular"))
ttops <- ttops[1:10,] # to make the example faster

branches <- make_branches(ttops, parallel = FALSE)

plot(LAS_from_xyz(branches))





}
\references{
\preformatted{* Paulina Fernandez M, Norero A (2006). Relation between length
and diameter of Pinus radiata branches.
Scandinavian Journal of Forest Research.
21: 124-129. - doi: 10.1080/02827580500533177
* Webster C, Mazzotti G, Essery R, Jonas T (2020). Enhancing
airborne LiDAR data for improved forest structure representation in
shortwave transmission models.
Remote Sensing of Environment.
249: 112017. - doi: 10.1016/j.rse.2020.112017
}
}
\seealso{
\code{\link{LAS_from_xyz}}, \code{\link{make_trunks}}
}
