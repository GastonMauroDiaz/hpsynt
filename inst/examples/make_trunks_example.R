require(lidR)
path <- system.file("external/745_r75m_d1m_local.las", package = "hpsynt")
las <- readLAS(path)

chm <- grid_canopy(las, 1, p2r(0.05, knnidw()))
ttops <- find_trees(las, lmf(3.7, shape = "circular"))
ttops <- ttops[1:10,] # to make the example faster

trunks <- make_trunks(ttops, parallel = FALSE)

plot(LAS_from_xyz(trunks))
