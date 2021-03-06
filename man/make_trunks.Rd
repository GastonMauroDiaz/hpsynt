% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/make_trunks.R
\name{make_trunks}
\alias{make_trunks}
\title{Make point clouds that represent trunks}
\usage{
make_trunks(
  ttops,
  dbh_fun = function(x) 8.867643 + x * 1.349327,
  spacing = c(0.05, 0.3),
  cylinder_height = 1.5,
  parallel = TRUE
)
}
\arguments{
\item{ttops}{output from \code{link[lidR]{find_trees}}.}

\item{dbh_fun}{function. It should take \emph{height} in m and return \emph{dbh} in cm.}

\item{spacing}{numeric vector. Range that control spacing between points.
The first number are for points near the center and the last for the
farthest points.}

\item{cylinder_height}{The height of a cylinder that represents the lower
part of the trunk. See reference.}

\item{parallel}{logical. Allows parallel processing.}
}
\value{
\code{data.frame} with columns x, y, and z.
}
\description{
Make point clouds that represent trunks
}
\examples{
require(lidR)
path <- system.file("external/745_r75m_d1m_local.las", package = "hpsynt")
las <- readLAS(path)

chm <- grid_canopy(las, 1, p2r(0.05, knnidw()))
ttops <- find_trees(las, lmf(3.7, shape = "circular"))
ttops <- ttops[1:10,] # to make the example faster

trunks <- make_trunks(ttops, parallel = FALSE)

plot(LAS_from_xyz(trunks))
}
\references{
Webster C, Mazzotti G, Essery R, Jonas T (2020). Enhancing
airborne LiDAR data for improved forest structure representation in
shortwave transmission models.
Remote Sensing of Environment. 249: 112017. - doi: 10.1016/j.rse.2020.112017
}
\seealso{
\code{\link{LAS_from_xyz}}, \code{\link{make_branches}}
}
