require(lidR)

x <- runif(1000, 0, 50)
y <- runif(1000, 0, 50)
z <- runif(1000, 0, 20)

las <- LAS_from_xyz(data.frame(x, y, z))
las
plot(las)




