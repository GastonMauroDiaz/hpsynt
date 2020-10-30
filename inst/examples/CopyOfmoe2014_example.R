require(lidR)
las <- readLAS("inst/external/745_r75m_d1m_local.las")

fun <- function(X,Y,Z) moe2014(X,Y,Z,
                               lens_height = 1,
                               size_weight = 0.5,
                               image_size = 1000,
                               filename = "temp.jpg")

cloud_metrics(my_las, fun(X,Y,Z))
