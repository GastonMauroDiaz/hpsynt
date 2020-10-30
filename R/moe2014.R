#' Moese et al. (2014)'s method
#'
#' @inheritParams var2012
#' @param size_weight numeric. Single value. It controls the printed sizes. This parameter
#'    was introduced since the original method uses Matlab to print (see reference), which
#'    seems to manage point size different than R.
#'
#' @return None. It writes a JPEG on \code{filename}.
#' @export moe2014
#'
#' @references Moeser D, Roubinek J, Schleppi P, Morsdorf F, Jonas T (2014). Canopy closure,
#' LAI and radiation transfer from airborne LiDAR synthetic images. Agricultural and
#' Forest Meteorology. 197: 158-168. - doi: 10.1016/j.agrformet.2014.06.008
#'
#' @seealso \code{\link{var2012}}
#' @example /inst/examples/moe2014_example.R
#' importFrom grDevice dev.off
moe2014 <- function(X, Y, Z,
                    lens_height = 1,
                    size_weight = 0.5,
                    image_size = 1000,
                    filename) {

  X <- X -mean(X)
  Y <- Y - mean(Y)
  X <- X[Z > lens_height]
  Y <- Y[Z > lens_height]
  Z <- Z[Z > lens_height]
  Z <- Z - lens_height

  # to spherical coordinates
  xyz <- cbind(X, Y, Z)
  sph <- geometry::cart2sph(xyz)

  # build hp
  # the trick is to think the spherical as polar
  # theta is r
  prz <- cbind(sph[,1], pi/2 - sph[,2], Z)
  xyz <- geometry::pol2cart(prz)

  xy <- data.frame(xyz[,1], xyz[,2])
  xy[,1] <- -xy[,1] # flip in east-west direction

  model <- stats::line(c(0,pi/2), c(7, 0.5)) #parameters from reference
  fun <- function(x) model$coefficients[1] + x * model$coefficients[2]
  cex <- fun(prz[,2]) * size_weight
  xlim = c(-pi/2, pi/2)

  grDevices::jpeg(filename, image_size, image_size, "px", quality = 100)
  graphics::par(mar = c(0,0,0,0), xaxs = "i", yaxs = "i")
  plot(xy,
       xlim = xlim, ylim = xlim,
       pch = 16, cex = cex)

  grDevices::dev.off()

}
