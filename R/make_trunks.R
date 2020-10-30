
.getCircleEdgeXY <-
  function(center = c(0, 0),
           diameter = 1,
           npoints = 100) {
    radio = diameter / 2
    t <- seq(0, 2 * pi, length.out = npoints)
    x <- center[1] + radio * cos(t)
    y <- center[2] + radio * sin(t)
    return(data.frame(x, y))
  }


.make_cylinder <- function(x, y, dbh, spacing, height) { #single tree coordinates
   dbh <- dbh/100 # cm to m
   circle <- .getCircleEdgeXY(c(x, y), diameter = dbh, pi * dbh / spacing)
   fun <- function(z) cbind(circle, z)
   cylinder <- Map(fun, seq(0, height, spacing))
   x <- unlist(Map(function(i) cylinder[[i]][,1], seq_along(cylinder)))
   y <- unlist(Map(function(i) cylinder[[i]][,2], seq_along(cylinder)))
   z <- unlist(Map(function(i) cylinder[[i]][,3], seq_along(cylinder)))
   data.frame(x = x, y = y, z = z)
}

.make_cone <- function(x, y, #single tree coordinates
                       dbh,
                       spacing,
                       cylinder_height,
                       height
                       ) {
  dbh <- dbh/100 # cm to m
  steps <- seq(cylinder_height + spacing, height, spacing)
  diameters <- seq(dbh, 0.01,length.out = length(steps))

  fun1 <- function(diameter, z) {
    circle <- .getCircleEdgeXY(c(x, y),
                               diameter,
                               pi * diameter / spacing)
    cbind(circle, z)
  }

  cone <- Map(fun1, diameters, steps)
  x <- unlist(Map(function(i) cone[[i]][,1], seq_along(cone)))
  y <- unlist(Map(function(i) cone[[i]][,2], seq_along(cone)))
  z <- unlist(Map(function(i) cone[[i]][,3], seq_along(cone)))
  data.frame(x = x, y = y, z = z)
}

.calc_spacing_per_tree <- function(ttops, spacing) {
  # to spherical coordinates to calculates distance from center for
  # calculating spacing with decay
  xyz <- cbind(ttops@coords[,1] - mean(ttops@coords[,1]),
               ttops@coords[,2] - mean(ttops@coords[,2]),
               0)
  sph <- geometry::cart2sph(xyz)
  distance_from_center <- sph[,3]
  spacing <- seq(spacing[1], spacing[2],
                 length.out = length(distance_from_center))
  indices <- match(distance_from_center, sort(distance_from_center))
  spacing <- spacing[indices]
  spacing

}

#' Make point clouds that represent trunks
#'
#' @param ttops output from \code{link[lidR]{find_trees}}.
#' @param dbh_fun function. It should take *height* in m and return *dbh* in cm.
#' @param spacing numeric vector. Range that control spacing between points.
#'     The first number are for points near the center and the last for the
#'     farthest points.
#' @param cylinder_height The height of a cylinder that represents the lower
#'     part of the trunk. See reference.
#' @param parallel logical. Allows parallel processing.
#'
#'
#' @references Webster C, Mazzotti G, Essery R, Jonas T (2020). Enhancing
#'     airborne LiDAR data for improved forest structure representation in
#'     shortwave transmission models.
#'     Remote Sensing of Environment. 249: 112017. - doi: 10.1016/j.rse.2020.112017
#'
#' @return \code{data.frame} with columns x, y, and z.
#' @export make_trunks
#'
#' @seealso \code{\link{LAS_from_xyz}}, \code{\link{make_branches}}
#' @example /inst/examples/make_trunks_example.R
#' @import lidR
make_trunks <- function(ttops,
                        dbh_fun = function(x) 8.867643 + x * 1.349327,
                        spacing = c(0.05, 0.3),
                        cylinder_height = 1.5,
                        parallel = TRUE
                        ) {

  spacing <- .calc_spacing_per_tree(ttops, spacing)

  dbhs <- dbh_fun(ttops@data$Z)

  # cylinders
  fun2 <- function(x, y, dbh, spacing) {
    .make_cylinder(x, y, dbh, spacing, cylinder_height)
  }
  cylinders <- Map(fun2, ttops@coords[,1], ttops@coords[,2], dbhs, spacing)
  x <- unlist(Map(function(i) cylinders[[i]][,1], seq_along(cylinders)))
  y <- unlist(Map(function(i) cylinders[[i]][,2], seq_along(cylinders)))
  z <- unlist(Map(function(i) cylinders[[i]][,3], seq_along(cylinders)))
  cylinders <- data.frame(x = x, y = y, z = z)

  # cones
  fun3 <- function(x, y, dbh, spacing, height) {
    .make_cone(x, y, dbh, spacing, cylinder_height, height)
  }

  if(parallel) {
    # go parallel
    ## Initiate cluster
    no_threads <- parallel::detectCores() - 1
    cl <- parallel::makeCluster(no_threads)
    parallel::clusterExport(cl,
                            c(".make_cylinder",
                              ".make_cone",
                              ".getCircleEdgeXY"),
                            environment())
    cones <- parallel::clusterMap(cl, fun3,
                                     ttops@coords[,1], ttops@coords[,2],
                                     dbhs, spacing, ttops@data$Z)

    ## finish
    parallel::stopCluster(cl)
  } else {
    cones <- Map(fun3, ttops@coords[,1], ttops@coords[,2],
                 dbhs, spacing, ttops@data$Z)
  }


  x <- unlist(Map(function(i) cones[[i]][,1], seq_along(cones)))
  y <- unlist(Map(function(i) cones[[i]][,2], seq_along(cones)))
  z <- unlist(Map(function(i) cones[[i]][,3], seq_along(cones)))
  cones <- data.frame(x = x, y = y, z = z)

  # cones on top of cylinders
  rbind(cylinders, cones)

}
