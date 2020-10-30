.inverse_distance_weight <- function(distance,
                                     sphere_radius = 0.15/2, #parameter from paper
                                     picture_plane_distance = 0.008) {
  # picture_plane_distance is from linear perspective theory.
  # In normal lenses, this value is the focal length if you
  # want a canvas (picture_plane) of 24*36 mm size (35 mm film).
  # In HP the canvas should be a square of 24*24 mm

  tan_alpha <- sphere_radius*0.5/distance
  radius_24mm <- tan_alpha*picture_plane_distance
  #depicted sphere in cm for 20cm canvas
  radius_20cm <- (20 * radius_24mm / 2.4) * 100
  radius_20cm * 2

}

# sphere_radius based on GSD
.inverse_distance_weight(0.75)
.inverse_distance_weight(0.75, 0.07/2)

.radial_distortion <- function(x, theta) {
  theta <- theta * 180/pi
  x + x * theta/90  * (pi/2 - 1)
}
# check. According to the authors, below must be equal to 1.57...
.radial_distortion(1, pi/2)

#' Varhola et al. (2012)'s method
#'
#' @param X numeric vector. X coordinates of 3D points (a point cloud).
#' @param Y numeric vector. Y coordinates of 3D points (a point cloud).
#' @param Z numeric vector. Z coordinates of 3D points (a point cloud).
#' @param lens_height numeric. Single value. Height above ground of the virtual lens.
#' @param sphere_diameter numeric. Single value. Diameter in meters.
#' @param minimum_base_constant numeric. Single value. It controls the minimum
#'     size to be printed. See reference for details.
#' @param distance_thr numeric. Single value. Minimum allowed distance to virtual lens.
#' @param image_size numeric. Single value. Width in pixels
#' @param filename character. Argument to pass to \code{\link[grDevices]{jpeg}}.
#'
#' @return None. It writes a JPEG on \code{filename}.
#' @export var2012
#'
#' @references Varhola A, Frazer GW, Teti P, Coops NC (2012). Estimation of forest
#' structure metrics relevant to hydrologic modeling using coordinate
#' transformation of airborne laser scanning data. Hydrology and Earth System Sciences
#' Discussions. 9: 5531-5573. - doi: 10.5194/hessd-9-5531-2012
#'
#' @seealso \code{\link{moe2014}}
#' @example /inst/examples/var2012_example.R
#' importFrom grDevice dev.off
var2012 <- function(X, Y, Z,
                    lens_height = 1,
                    sphere_diameter = 0.15,
                    minimum_base_constant = 0.215,
                    distance_thr = 0.75,
                    image_size = 750,
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
  indices <- prz[,3] > distance_thr
  prz <- prz[indices,]
  distance_from_lens <- sph[indices, 3]
  xyz <- geometry::pol2cart(prz)

  xy <- data.frame(xyz[,1], xyz[,2])
  xy[,1] <- -xy[,1] # flip in east-west direction

  cex <- .inverse_distance_weight(distance_from_lens, sphere_diameter/2)
  cex <- .radial_distortion(cex, prz[,2] )
  cex[cex < minimum_base_constant] <- minimum_base_constant

  xlim = c(-pi/2, pi/2)

  # files <- dir(output_path, pattern = "jpg", ignore.case = TRUE)
  # name <- paste0(output_path, "/", length(files) + 1, ".jpg")

  grDevices::jpeg(filename, 20, 20, "cm", res = 2.54*image_size/20,
                  quality = 100, pointsize = 64)
  # I calibrate pointsize to produce the correct
  # point size for 20 cm images so replicate the original research
    graphics::par(mar = c(0,0,0,0), xaxs = "i", yaxs = "i")
    plot(xy,
         xlim = xlim, ylim = xlim,
         pch = 16, cex = cex)

  grDevices::dev.off()

}
