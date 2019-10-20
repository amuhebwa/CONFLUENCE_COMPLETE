#.rs.restartR()
rm(list=ls())

library(ncdf4)  

library(tidyverse)
library(tidyr)
library(dplyr)
library(schoolmath)
library(hydroTSM)
library(quantmod)
library(matrixStats)

# setwd("/Users/amuhebwa/Documents/Python/BAM_PREDIAGNOSTICS")
filepath = paste(getwd(), "/original_data/", sep="")
riverreach_files <-list.files(filepath, pattern = "*.nc", full.names = FALSE)

for (files in riverreach_files){
  Reach_in=nc_open(paste(filepath,files,sep=""))
  QReach <- ncvar_get(Reach_in,"XS_Timseries/Q")
  S1k <- t(ncvar_get(Reach_in,"Reach_Timeseries/S_1km"))
  S90 <- t(ncvar_get(Reach_in,"Reach_Timeseries/S_90m"))
  H1k <- t(ncvar_get(Reach_in,"XS_Timseries/H_1km"))
  H90 <- t(ncvar_get(Reach_in,"XS_Timseries/H_90m"))
  Wxs <- t(ncvar_get(Reach_in,"XS_Timseries/W"))
  GID <- ncvar_get(Reach_in,"Reach_Timeseries/GridID")
  HID <- ncvar_get(Reach_in,"Reach_Timeseries/HydroID")
  NDID <- ncvar_get(Reach_in,"Reach_Timeseries/NextDownID") 
  
  ## Filtering the over bank flows using the 2-year flood criteria     
  
  H_Reach=rowMeans(H1k)
  H_Reach[is.negative(H_Reach)] <- 0
  
  # Create dates as a Date class object starting from 2016-01-01
  dates <- seq(as.Date("1990-01-01"), length = 9131, by = "days")
  # Use xts() to create smith
  H_Reach_mean <- xts(x = H_Reach,order.by=dates)
  H_Reach_ann <- to.yearly(H_Reach_mean)
  WL_Fdc=fdc(H_Reach_ann)
  f <- splinefun(WL_Fdc[,2],H_Reach_ann[,2])
  OB_WL=f(c(.5))
  H_OBF_s= H_Reach > OB_WL
  Timestamp_OBF= which(H_OBF_s, arr.ind = FALSE, useNames = TRUE)
  
  ## Filtering the high water surface slope  
  
  Slope_threshold_up <- 0.003
  Slope_threshold_down <- 0
  
  # Mountainous area check
  SReach_mountainous= S1k > Slope_threshold_up
  Timestamp_Slopeup= which(SReach_mountainous, arr.ind = FALSE, useNames = TRUE)
  
  # Flowing upstream check
  
  SReach_flowupstream= S1k < Slope_threshold_down
  Timestamp_flowupstream= which(SReach_flowupstream, arr.ind = FALSE, useNames = TRUE)
  
  Qmean=colMeans(QReach)
  #Qmean=colMins(QReach)
  
  # Negative Discharge check
  
  QReachmean=colMeans(QReach)
  QReach_negative= QReachmean < 0
  Timestamp_Qneg= which(QReach_negative, arr.ind = FALSE, useNames = TRUE)
  # Height profile check
  
  dH=apply(H1k,1,diff)
  dH_max=colMaxs(dH)
  
  Height_profile = dH_max > 0
  Timestamp_profile_wrong= which(Height_profile, arr.ind = FALSE, useNames = TRUE)
  
  ##--------------------Zero width flags----------------------------------------------------
  Wreach=rowMeans(Wxs)
  Wxs_zerowidths= Wreach == 0
  Timestamp_zerowidths= which(Wxs_zerowidths, arr.ind = FALSE, useNames = TRUE)
  
  
  # FLAGS LEGEND #
  # -11111 is a flag for over-bank flow
  # -22222 is a flag for high slopes that indicate mountainous areas which is not correct 
  # -33333 is a flag for wrong flow direction (flow going upstream)
  # -44444 is a flag for negative LISFLOOD model discharge
  # -55555 is a flag for a water depth profile that is messed up (depth is not decreasing downstream)
  H1k[Timestamp_OBF,]=-11111
  H90[Timestamp_OBF,]=-11111
  S1k[Timestamp_OBF]=-11111
  S90[Timestamp_OBF]=-11111
  QReach[,Timestamp_OBF]=-11111
  Wxs[Timestamp_OBF,]=-11111
  
  H1k[Timestamp_Slopeup,]=-22222
  H90[Timestamp_Slopeup,]=-22222
  S1k[Timestamp_Slopeup]=-22222
  S90[Timestamp_Slopeup]=-22222
  QReach[,Timestamp_Slopeup]=-22222
  Wxs[Timestamp_Slopeup,]=-22222
  
  H1k[Timestamp_flowupstream,]=-33333
  H90[Timestamp_flowupstream,]=-33333
  S1k[Timestamp_flowupstream]=-33333
  S90[Timestamp_flowupstream]=-33333
  QReach[,Timestamp_flowupstream]=-33333
  Wxs[Timestamp_flowupstream,]=-33333
  
  H1k[Timestamp_Qneg,]=-44444
  H90[Timestamp_Qneg,]=-44444
  S1k[Timestamp_Qneg]=-44444
  S90[Timestamp_Qneg]=-44444
  QReach[,Timestamp_Qneg]=-44444
  Wxs[Timestamp_Qneg,]=-44444
  
  H1k[Timestamp_profile_wrong,]=-55555
  H90[Timestamp_profile_wrong,]=-55555
  S1k[Timestamp_profile_wrong]=-55555
  S90[Timestamp_profile_wrong]=-55555
  QReach[,Timestamp_profile_wrong]=-55555
  Wxs[Timestamp_profile_wrong,]=-55555
  
  H1k[Timestamp_zerowidths,]=-66666
  H90[Timestamp_zerowidths,]=-66666
  S1k[Timestamp_zerowidths]=-66666
  S90[Timestamp_zerowidths]=-66666
  QReach[,Timestamp_zerowidths]=-66666
  Wxs[Timestamp_zerowidths,]=-66666
  
  H1k=t(H1k)
  H90=t(H90)
  Wxs=t(Wxs)
  
  nc_close(Reach_in)
  rm(Reach_in)
  rm(H_Reach)
  
  filename=substr(files,1,(nchar(files)-3)) 
  reachname=substr(files,19,(nchar(files)-3))
  
  #-- writing the modified netcdf files -----------------
  #---------- Creating the Netcdf file ------------------
  
  ncpath = paste(getwd(), "/data/", sep="")
  ncname <- paste(filename,sep="")  
  ncfname <- paste(ncpath, ncname, ".nc", sep="")
  
  #---------- define variables-------------------------------
  
  Xsec90m <- as.array(seq(1,dim(H90)[1],1)) #90 m
  Xsec1k <- as.array(seq(1,dim(H90)[1],1)) # 1 km
  time <- as.array(seq(1,length(S90),1))
  Reach <- as.array(seq(1,1,1))
  
  Xsec90mdim <- ncdim_def("XS_90m","orthogonals",as.integer(Xsec90m))
  Xsec1kdim <- ncdim_def("Stageloc_1km","orthogonals",as.integer(Xsec1k))
  XSTimestepsdim <- ncdim_def("Time steps","day",as.integer(time))
  Reachdim <- ncdim_def("Reach","RiverReaches",as.integer(Reach))
  dimnchar <- ncdim_def("nchar", "", 1:12, create_dimvar=FALSE )
  
  #--------- Adding the priors---------------------------------
  
  #Prior information includes prior distrbution of mean annual flow (Q),
  #roughness coefficient (n), bankfull depth (z) and width (w)
  
  fillvalue = -999999
  dlname <- "Width_LISFLOOD_derived"
  XS_W <- ncvar_def("XS_Timseries/W","meters",list(Xsec90mdim,XSTimestepsdim),1.e30,dlname,prec="double")
  dlname <- "Water SUrface Elevation_90m"# 
  XS_H_90m <- ncvar_def("XS_Timseries/H_90m","meters",list(Xsec90mdim,XSTimestepsdim),1.e30,dlname,prec="double")
  dlname <- "Water SUrface Elevation_1km"
  XS_H_1km <- ncvar_def("XS_Timseries/H_1km","meters",list(Xsec1kdim,XSTimestepsdim),1.e30,dlname,prec="double")
  dlname <- "Discharge"
  XS_Q <- ncvar_def("XS_Timseries/Q","cubic meters per second",list(Xsec1kdim,XSTimestepsdim),1.e30,dlname,prec="double")
  dlname <- "Mean Discharge"
  Reach_Q <- ncvar_def("Reach_Timseries/Q","cubic meters per second",list(Reachdim,XSTimestepsdim),1.e30,dlname,prec="double")
  dlname <- "Slope_1km"
  ReachSlope_1km <- ncvar_def("Reach_Timeseries/S_1km","m/m",list(Reachdim,XSTimestepsdim),1.e30,dlname,prec="double")
  dlname <- "Slope_90m"
  ReachSlope_90m <- ncvar_def("Reach_Timeseries/S_90m","m/m",list(Reachdim,XSTimestepsdim),1.e30,dlname,prec="double")
  dlname <- "GridID"
  GridID <- ncvar_def("Reach_Timeseries/GridID","",list(dimnchar,Reachdim),prec="char")
  dlname <- "HydroID"
  HydroID <- ncvar_def("Reach_Timeseries/HydroID","dimensionless",list(Reachdim,Reachdim),fillvalue,dlname,prec="double")
  dlname <- "NextDownID"
  NextDownID <- ncvar_def("Reach_Timeseries/NextDownID","dimensionless",list(Reachdim,Reachdim),fillvalue,dlname,prec="double")
  
  
  #--------- create netCDF file and put arrays---------------------
  ncout <- nc_create(ncfname,list(XS_W,XS_H_90m,XS_H_1km,XS_Q,Reach_Q,ReachSlope_90m,ReachSlope_1km,GridID,HydroID,NextDownID),force_v4=TRUE)
  
  #---------------- Insert the variables 
  ncvar_put(ncout,XS_W,Wxs)
  ncvar_put(ncout,XS_H_90m,H90)
  ncvar_put(ncout,XS_H_1km,H1k)
  ncvar_put(ncout,XS_Q,QReach)
  ncvar_put(ncout,Reach_Q,Qmean)
  ncvar_put(ncout,ReachSlope_90m,S90)
  ncvar_put(ncout,ReachSlope_1km,S1k)
  ncvar_put(ncout,GridID,GID)
  ncvar_put(ncout,HydroID,HID)
  ncvar_put(ncout,NextDownID,NDID)
  
  
  #---------------- add global attributes
  title=paste("Severn SWOT-like data for reach", reachname)
  ncatt_put(ncout,0,"title",title)
  
  nc_close(ncout)
  
  rm(QReach)
  rm(S1k)
  rm(S90)
  rm(H1k)
  rm(H90)
  rm(Wxs)
}
