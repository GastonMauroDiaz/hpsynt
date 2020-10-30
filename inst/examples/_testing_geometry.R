require(geometry)
require(lidR)

tree <- readLAS("inst/external/single_tree.las")
tree <- filter_poi(tree, Z > 1)
plot(tree)
xyz <- as.matrix(tree@data[,1:3])
xyz[,1] <- xyz[,1] - mean(xyz[,1])
xyz[,2] <- xyz[,2] - mean(xyz[,2])
head(xyz)
foo <- geometry::convhulln(xyz, "FA")
foo$area
plot(foo)

m <- geometry::to.mesh3d(foo)

n <- 1000
rp_x <- runif(n, min(xyz[,1]), max(xyz[,1]))
rp_y <- runif(n, min(xyz[,2]), max(xyz[,2]))
rp_z <- runif(n, min(xyz[,3]), max(xyz[,3]))
rp <- cbind(rp_x, rp_y, rp_z)
head(rp)

indices <- inhulln(foo, rp)
indices

.writePC <- function(x, path) {

  if (ncol(x) == 3) cbind(x, rep(0, nrow(x)))
  write.table(x, path, quote = FALSE, row.names = FALSE, col.names = FALSE)
}

.writePC(rp[indices,], "tree_enhaced.xyz")
