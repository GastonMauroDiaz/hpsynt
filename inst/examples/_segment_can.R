require(lidR)

my_las <- readLAS("inst/external/point_cloud.las")

SegmentCan <- function(my_las){
  chm <- lidR::grid_canopy(my_las, 0.5, p2r(0.05, knnidw()))

  loc <- lidR::find_trees(my_las, lmf(3.7, shape = "circular"))
  loc <- data.frame(coordinates(loc), loc@data$Z)
  colnames(loc) <- c("x", "y", "height")

  canopy <- rLiDAR::ForestCAS(chm, loc)

  # remove small tree crowns
  canopy[[2]] <- canopy[[2]][!(canopy[[2]]$ca < 4),]

  canopy
}

ClassTrees <- function(canopy, #esto es lo que devuelve SegmentCam
                       my_las, #parece ser la ruta a un archivo las
                       species = 1,
                       biome = 10,
                       branches = 1){

  lasdat <- my_las

  # calculate crown diameter
  canopy[[2]]$cd_m <- 2*sqrt(canopy[[2]]$ca/pi)

  treedim <- itcSegment::dbh(H = canopy[[2]]$z, CA = canopy[[2]]$cd_m, biome = biome)

  boundaryTrees <- canopy[[1]]


  # empty matrix for classifying branch points (ltc)
  lastc         <- matrix(data=-NaN,nrow=length(lasdat$X),ncol=1)
  treex         <- matrix(data=-NaN,nrow=length(lasdat$X),ncol=1)
  treey         <- matrix(data=-NaN,nrow=length(lasdat$X),ncol=1)
  horzdist      <- matrix(data=-NaN,nrow=length(lasdat$X),ncol=1)
  ba            <- matrix(data=-NaN,nrow=length(lasdat$X),ncol=1)

  ltcdf <- data.frame(lasdat$X, lasdat$Y, lasdat$Z, treex, treey, horzdist, lastc, ba)

  colnames(ltcdf)[1] <- "lasdat_X"
  colnames(ltcdf)[2] <- "lasdat_Y"
  colnames(ltcdf)[3] <- "lasdat_Z"

  ltcdf$lasdat_Z[ltcdf$lasdat_Z < 1.0] = NaN # don't bother with branches within 1 metre of the ground
  ltcdf <- ltcdf[complete.cases(ltcdf$lasdat_Z), ]


  # empty matrix for tree trunk and crown statistics (dbh)
  treeptsX      <- matrix(data=NaN,nrow=length(boundaryTrees$Trees),ncol=1)
  treeptsY      <- matrix(data=NaN,nrow=length(boundaryTrees$Trees),ncol=1)
  treeptsH_m    <- matrix(data=NaN,nrow=length(boundaryTrees$Trees),ncol=1)
  treeptsCA_m2  <- matrix(data=NaN,nrow=length(boundaryTrees$Trees),ncol=1)
  treeptsdbh_cm <- matrix(data=NaN,nrow=length(boundaryTrees$Trees),ncol=1)

  dbhdf <- data.frame(treeptsX, treeptsY, treeptsH_m, treeptsdbh_cm, treeptsCA_m2)


  # classify points in lidar cloud based on what tree boundary it is in
  for (i in seq_along(1:length(boundaryTrees$Trees))){

    pix <- boundaryTrees@polygons[[i]]@plotOrder[1]
    outline <- boundaryTrees@polygons[[i]]@Polygons[[pix]]@coords

    # find which chm maxima is in that polygon because loc and boundaryTrees don't match
    tpts_in <- point.in.polygon(canopy[[2]]$x,canopy[[2]]$y,outline[,1],outline[,2], mode.checked = FALSE)

    if (length(which(tpts_in == 1) > 0)) {
      dbhdf$treeptsX[i]      <- canopy[[2]]$x[which(tpts_in==1)]
      dbhdf$treeptsY[i]      <- canopy[[2]]$y[which(tpts_in==1)]
      dbhdf$treeptsH_m[i]    <- canopy[[2]]$z[which(tpts_in==1)]
      dbhdf$treeptsCA_m2[i]  <- canopy[[2]]$ca[which(tpts_in==1)]
      dbhdf$treeptsdbh_cm[i] <- treedim[which(tpts_in==1)]
    }

    if (branches == 1){ # only bother looping through and classifying branches if the user wants them
      # now classify the lidar by which polygon it falls in
      if (length(which(tpts_in > 0) > 0)) {
        pts_in <- point.in.polygon(ltcdf$lasdat_X,
                                   ltcdf$lasdat_Y,
                                   outline[,1],
                                   outline[,2],
                                   mode.checked = FALSE)

        ltcdf$lastc[which(pts_in>0)] <- i

        ltcdf$treex[which(pts_in>0)]    <- boundaryTrees@polygons[[i]]@Polygons[[pix]]@labpt[1]
        ltcdf$treey[which(pts_in>0)]    <- boundaryTrees@polygons[[i]]@Polygons[[pix]]@labpt[2]
        ltcdf$horzdist[which(pts_in>0)] <- sqrt(((ltcdf$lasdat_X[which(pts_in>0)]-
                                                    ltcdf$treex[which(pts_in>0)])**2)
                                                + ((ltcdf$lasdat_Y[which(pts_in>0)]-
                                                      ltcdf$treey[which(pts_in>0)])**2))

        # assign a branch angle
        if (species == 1){
          ang <- runif(1, 60, 100)
        } else if (species == 2){
          ang <- runif(1, 50, 100)
        } else if (species == 3){
          ang <- runif(1, 45, 75)
        } else if (species == 4){
          ang <- runif(1, 32, 50)
        } else {
          ang <- 90
        }

        ltcdf$ba[which(pts_in>0)] <- ang
        rm(ang)

      }

    }

  }

  ltcdf$horzdist[ltcdf$horzdist < 0.2] = NaN

  return(list(ltcdf=ltcdf,
              dbhdf=dbhdf,
              canopy=canopy,
              treedim=treedim))

}

GenGrid <- function(dbhdf, #classTrees devuelve una lista con esto
                    canopy, #esto es lo que devuelve SegmentCam. classTrees devuelve una lista con esto
                    treedim, #classTrees devuelve una lista con esto. Es esto itcSegment::dbh(H = canopy[[2]]$z, CA = canopy[[2]]$cd_m, biome = biome)
                    xlim,
                    ylim,
                    spacing,
                    epsgstr){

  Xmat <- seq(from = xlim[1], to = xlim[2], by = spacing)
  Ymat <- seq(from = ylim[1], to = ylim[2], by = spacing)
  cnew <- meshgrid(Xmat,Ymat)
  xnew <- c(cnew$X)
  ynew <- c(cnew$Y)

  rad <- ((treedim/100) / 2)

  df <- as.data.frame(canopy[[2]])
  xy <- df[,c(1,2)]

  if (epsgstr == '0000'){ # doesn't remove points in trunks

    cat("\n \n Not processing grid points for those in trunks. \n")


  } else { # remove points within trunks

    projstr = strcat('+init=EPSG:',epsgstr)

    spdf <- SpatialPointsDataFrame(coords = xy,data=df,proj4string = CRS(projstr))

    sptA <- gBuffer(spdf,byid=TRUE,id=NULL,width=rad,quadsegs = 10,
                    capStyle = "ROUND",joinStyle = "ROUND")

    lastcA <- matrix(data=-NaN,nrow=length(xnew),ncol=1)

    for (i in seq_along(1:length(sptA$X))){

      outline <- sptA@polygons[[i]]@Polygons[[1]]@coords

      pts_in <- point.in.polygon(xnew,ynew,outline[,1],outline[,2], mode.checked = FALSE)

      lastcA[which(pts_in==1)] <- 1

    }

    ptsdf = data.frame(xnew,ynew)
    ptsdf$xnew[which(lastcA>0)] <- NaN
    ptsdf$ynew[which(lastcA>0)] <- NaN

    ptsdf <- ptsdf[complete.cases(ptsdf),]

    return(ptsdf)

  }

}

canopy <- SegmentCan(my_las)
foo <- ClassTrees(canopy, my_las)


chm <- lidR::grid_canopy(my_las, 0.5, p2r(0.05, knnidw()))

require(pracma)
require(rgeos)
baz <- GenGrid(foo$dbhdf, canopy, foo$treedim,
               c(xmin(chm), xmax(chm)),
               c(ymin(chm), ymax(chm)),
               1,
               "22181")


head(foo$treedim)
head(foo$dbhdf)

writePC(foo$ltcdf, "temp.xyz")

stop()
chm <- grid_canopy(my_las, 0.5, p2r(0.05, knnidw()))

p <- find_trees(my_las, lmf(3.7, shape = "circular"))
crowns <- silva2016(chm, p)()
