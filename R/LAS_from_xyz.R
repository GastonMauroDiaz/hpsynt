
#' Transform raw coordinates into an object from the LAS class
#'
#' @param xyz data.frame or matrix with 3D coordinates (xyz).
#'
#' @return An object from the \code{\linkS4class{LAS}}.
#' @export LAS_from_xyz
#'
#' @seealso \code{\link[lidR]{writeLAS}}
#' @example /inst/examples/LAS_from_xyz_example.R
#' importFrom methods as new (madness!!)
LAS_from_xyz <- function(xyz) {
  .make_fake_las <- function(X, Y, Z){
    data_template_names <- c("X", "Y", "Z", "gpstime","Intensity", "ReturnNumber", "NumberOfReturns",
                             "ScanDirectionFlag", "EdgeOfFlightline",  "Classification",
                             "Synthetic_flag", "Keypoint_flag", "Withheld_flag", "ScanAngleRank", "UserData",
                             "PointSourceID",  "R",  "G",  "B")
    data_template <- matrix(ncol = length(data_template_names), nrow = length(X))
    data_template <- data.frame(data_template)
    names(data_template) <- data_template_names
    data_template[] <- 0

    fake_las <- new("LAS")
    fake_las@data <- as(data_template, "data.table")
    fake_las@data$X <- X
    fake_las@data$Y <- Y
    fake_las@data$Z <- Z
    fake_las
  }
  .make_fake_las(xyz[,1],xyz[,2],xyz[,3])
}
