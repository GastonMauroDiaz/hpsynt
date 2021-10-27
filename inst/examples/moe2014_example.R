require(lidR)
path <- system.file("external/745_r75m_d1m_local.las", package = "hpsynt")
las <- readLAS(path)

las <- clip_circle(las, 0, 0, 10)

chm <- grid_canopy(las, 1, p2r(0.05, knnidw()))
ttops <- find_trees(las, lmf(3.7, shape = "circular"))
las <- segment_trees(las, silva2016(chm, ttops))

las <- filter_poi(las, treeID == 5 | treeID == 10)

station_point <- c(0,0,0)

filename <- tempfile("moe2014", fileext = ".jpg")
moe2014(las@data$X, las@data$Y, las@data$Z,
        station_point,
        size_weight = 1,
        image_size = 1000,
        filename = filename,
        col = las@data$treeID)
canopy <- brick(filename)
plotRGB(canopy)


# \dontrun{
#         # demonstrating how to combine functions
#         chm <- grid_canopy(las, 1, p2r(0.05, knnidw()))
#         ttops <- find_trees(las, lmf(3.7, shape = "circular"))
#         trunks <- make_trunks(ttops, spacing = c(0.05, 2))
#         branches <- make_branches(ttops, spacing = c(0.05, 2))
#
#         las2 <- LAS_from_xyz(rbind(trunks, branches))
#         plot(las2)
#
#         filename2 <- tempfile("moe2014", fileext = ".jpg")
#         moe2014(las2@data$X, las2@data$Y, las2@data$Z,
#                 station_point,
#                 size_weight = 0.5,
#                 image_size = 1000,
#                 filename = filename2)
#         wood <- brick(filename2)
#         plotRGB(wood)
#
#         # combine data
#         plotRGB(wood * canopy)
#
#         # cleaning
#
#         unlink(filename2)
#
# }
unlink(filename)

