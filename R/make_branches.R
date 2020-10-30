
.make_whorl <- function(branch_lenght,
                        branch_diameter,
                        branch_angle,
                        branchs_per_whorl,
                        spacing) {


  #make one branch as if it were a vertical cone (pich)
  if (branch_lenght <= spacing) branch_lenght <- spacing * 1.1
  branch <- .make_cone(0,0, branch_diameter * 100, spacing, 0, branch_lenght)
  # branch <- seq(0, branch_lenght, spacing)
  # branch <- data.frame(x = 0, y = 0, z = branch)

  #incline the branch
  rotmat <- orientlib::eulerzyx(0,0, -branch_angle * pi / 180)
  rotmat <- orientlib::rotmatrix(rotmat)

  for (i in 1:nrow(branch)) {
    branch[i,] <- rotmat@x[,,1] %*% matrix(as.numeric(branch[i,]) , ncol = 1)
  }

  # duplicate and move the inclined branch around trunk axis (yaw)
  fun <- function(branch, yaw) {
    rotmat <- orientlib::eulerzyx(yaw,0,0)
    rotmat <- orientlib::rotmatrix(rotmat)
    for (i in 1:nrow(branch)) {
      branch[i,] <- rotmat@x[,,1] %*% matrix(as.numeric(branch[i,]) , ncol = 1)
    }
    branch
  }

  yaw <- 2 * pi / branchs_per_whorl

  randomness <- 0.3
  branches <- list()
  branches[[1]] <- fun(branch, rnorm(1, yaw, yaw * randomness) - yaw)
  for (i in 2:branchs_per_whorl) {
    branches[[i]] <- fun(branches[[i-1]], rnorm(1, yaw, yaw * randomness))
  }

  x <- unlist(Map(function(i) branches[[i]][,1], seq_along(branches)))
  y <- unlist(Map(function(i) branches[[i]][,2], seq_along(branches)))
  z <- unlist(Map(function(i) branches[[i]][,3], seq_along(branches)))

  data.frame(x, y, z)
}

.make_branches_per_tree <- function(x, y, height,
                                    crown_length, crown_diameter,
                                    branchs_per_whorl,
                                    branch_angle,
                                    interval,
                                    branch_diameter_fun,
                                    spacing
                                    ) {


  whorls_heights <- seq(height - crown_length,
                        height,
                        interval)

  .calc_crown_diameter_for_any_section <- function(section_height) {
    model <- stats::line(c(height - crown_length, height), c(crown_diameter, 0))
    model$coefficients[1] + section_height * model$coefficients[2]
  }

  branch_lenghts <- .calc_crown_diameter_for_any_section(whorls_heights) / 2
  branch_diameters <- branch_diameter_fun(branch_lenghts)

  fun <- function(branch_lenght, branch_diameter) {
    .make_whorl(branch_lenght, branch_diameter,
                branch_angle = branch_angle,
                branchs_per_whorl = branchs_per_whorl,
                spacing = spacing)
  }

  indices <- branch_lenghts > branch_diameters

  #make all needed branches, but with c(0,0) origin
  branches <- Map(fun, branch_lenghts[indices], branch_diameters[indices])

  #georreference branches
  x <- unlist(Map(function(i) branches[[i]][,1] + x, seq_along(branches)))
  y <- unlist(Map(function(i) branches[[i]][,2] + y, seq_along(branches)))
  z <- unlist(Map(function(i) branches[[i]][,3] + whorls_heights[i], seq_along(branches)))

  return(data.frame(x, y, z))

}


#' Make point clouds that represent branches
#'
#' @inheritParams make_trunks
#' @param crown_length_fun function. It should take
#'     *height* in m and return *crown lenght* in m.
#' @param crown_diameter_fun function. It should take
#'     *height* in m and return *crown diameter* in m.
#' @param branchs_per_whorl numeric vector. Single value.
#' @param branch_angle numeric. Single value. Angle in degrees.
#' @param interval numeric vector. single value. It is
#'     the mean length between successive whorls.
#' @param branch_diameter_fun function. It should take
#'     *branch length* in m and return
#'     *branch diameter* in m.
#'
#'
#' @return \code{data.frame} with columns x, y, and z.
#' @export make_branches
#'
#' @details \code{crown_length_fun} and \code{crown_diameter_fun} were
#'     fitted by the author of this package as part of a work in progress.
#'     \code{branch_diameter_fun} is from Fernandez and Norero (2006).
#'
#' @references
#'     * Paulina Fernandez M, Norero A (2006). Relation between length
#'     and diameter of Pinus radiata branches.
#'     Scandinavian Journal of Forest Research.
#'     21: 124-129. - doi: 10.1080/02827580500533177
#'     * Webster C, Mazzotti G, Essery R, Jonas T (2020). Enhancing
#'     airborne LiDAR data for improved forest structure representation in
#'     shortwave transmission models.
#'     Remote Sensing of Environment.
#'     249: 112017. - doi: 10.1016/j.rse.2020.112017
#'
#' @seealso \code{\link{LAS_from_xyz}}, \code{\link{make_trunks}}
#' @example /inst/examples/make_branches_example.R
#' @importFrom stats rnorm
make_branches <- function(ttops,
                          crown_length_fun = function(z) -1.3034187 + z * 0.6937502,
                          crown_diameter_fun = function(z) 0.8787043 + z * 0.2521690,
                          branchs_per_whorl = 5,
                          branch_angle = 80,
                          interval = 1,
                          branch_diameter_fun = function(branch_length) 0.002942 + branch_length * 0.0103,
                          spacing = c(0.05, 0.3),
                          parallel = TRUE
                          ){

  spacing <- .calc_spacing_per_tree(ttops, spacing)

  crown_length <- crown_length_fun(ttops@data$Z)
  crown_diameter <- crown_diameter_fun(ttops@data$Z)

  fun <- function(x, y, height, crown_length, crown_diameter, spacing) {
     .make_branches_per_tree(x, y, height,
                             crown_length, crown_diameter,
                             branchs_per_whorl = branchs_per_whorl,
                             branch_angle = branch_angle,
                             interval = interval,
                             branch_diameter_fun = branch_diameter_fun,
                             spacing)
  }


  if(parallel) {
    # go parallel
    ## Initiate cluster
    no_threads <- parallel::detectCores() - 1
    cl <- parallel::makeCluster(no_threads)
    parallel::clusterExport(cl,
                            c(".make_branches_per_tree",
                              ".make_whorl",
                              ".make_cone",
                              ".getCircleEdgeXY"),
                            environment())
    branches <- parallel::clusterMap(cl, fun,
                                     ttops@coords[,1], ttops@coords[,2], ttops@data$Z,
                                     crown_length, crown_diameter, spacing)
    ## finish
    parallel::stopCluster(cl)
  } else {
    branches <- Map(fun,
                    ttops@coords[,1], ttops@coords[,2], ttops@data$Z,
                    crown_length, crown_diameter, spacing)
  }

  x <- unlist(Map(function(i) branches[[i]][,1], seq_along(branches)))
  y <- unlist(Map(function(i) branches[[i]][,2], seq_along(branches)))
  z <- unlist(Map(function(i) branches[[i]][,3], seq_along(branches)))

  data.frame(x = x, y = y, z = z)
}

