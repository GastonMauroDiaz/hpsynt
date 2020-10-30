require(lidR)
path <- system.file("external/745_r75m_d1m_local.las", package = "hpsynt")
las <- readLAS(path)

filename <- tempfile("var2012", fileext = ".jpg")
var2012(las@data$X, las@data$Y, las@data$Z,
       lens_height = 1,
       sphere_diameter = 0.5,
       minimum_base_constant = 0.0215,
       distance_thr = 0.75,
       image_size = 2000,
       filename = filename
       )
plotRGB(brick(filename))
unlink(filename)
