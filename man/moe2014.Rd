% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/moe2014.R
\name{moe2014}
\alias{moe2014}
\title{Moese et al. (2014)'s method}
\usage{
moe2014(
  X,
  Y,
  Z,
  station_point,
  lens_height = 1,
  size_weight = 0.5,
  image_size = 1000,
  filename,
  col = NULL
)
}
\arguments{
\item{X}{numeric vector. X coordinates of 3D points (a point cloud).}

\item{Y}{numeric vector. Y coordinates of 3D points (a point cloud).}

\item{Z}{numeric vector. Z coordinates of 3D points (a point cloud).}

\item{lens_height}{numeric. Single value. Height above ground of the virtual lens.}

\item{size_weight}{numeric. Single value. It controls the printed sizes. This parameter
was introduced since the original method uses Matlab to print (see reference), which
seems to manage point size different than R.}

\item{image_size}{numeric. Single value. Width in pixels}

\item{filename}{numeric vector. Argument to pass to \code{\link[grDevices]{jpeg}}.
It must have the same length than X.}
}
\value{
None. It writes a JPEG on \code{filename}.
}
\description{
Moese et al. (2014)'s method
}
\examples{
require(lidR)
path <- system.file("external/745_r75m_d1m_local.las", package = "hpsynt")
las <- readLAS(path)

filename <- tempfile("moe2014", fileext = ".jpg")
moe2014(las@data$X, las@data$Y, las@data$Z,
        size_weight = 1,
        image_size = 1000,
        filename = filename)
canopy <- brick(filename)
plotRGB(canopy)

\dontrun{
        # demonstrating how to combine functions
        chm <- grid_canopy(las, 1, p2r(0.05, knnidw()))
        ttops <- find_trees(las, lmf(3.7, shape = "circular"))
        trunks <- make_trunks(ttops, spacing = c(0.05, 2))
        branches <- make_branches(ttops, spacing = c(0.05, 2))

        las2 <- LAS_from_xyz(rbind(trunks, branches))
        plot(las2)

        filename2 <- tempfile("moe2014", fileext = ".jpg")
        moe2014(las2@data$X, las2@data$Y, las2@data$Z,
                size_weight = 0.5,
                image_size = 1000,
                filename = filename2)
        wood <- brick(filename2)
        plotRGB(wood)

        # combine data
        plotRGB(wood * canopy)

        # cleaning

        unlink(filename2)

}
unlink(filename)

}
\references{
Moeser D, Roubinek J, Schleppi P, Morsdorf F, Jonas T (2014). Canopy closure,
LAI and radiation transfer from airborne LiDAR synthetic images. Agricultural and
Forest Meteorology. 197: 158-168. - doi: 10.1016/j.agrformet.2014.06.008
}
\seealso{
\code{\link{var2012}}
}
