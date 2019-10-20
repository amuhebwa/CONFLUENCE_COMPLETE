rm(list=ls())


library(ncdf4)  
library(tidyverse)
library(tidyr)
library(dplyr)
library(readxl)
library(schoolmath)

#setwd("/Users/amuhebwa/Documents/Python/CONFLUENCE_POSTDIAGNOSTICS")
filepath <- paste(getwd(), "/integrator.nc", sep="")
routing_table <- read_excel(paste(getwd(),"/Routing_table.xlsx", sep = ""))

Reach_in = nc_open(filepath)

# NOTE THAT THESE VARIABLE NAMES IN THE NETCDF FILES ARE GOING TO CHANGE
Qstar <- t(ncvar_get(Reach_in,"Qstar"))
Qpred <- t(ncvar_get(Reach_in,"Qpred"))
Q_Integrator <- t(ncvar_get(Reach_in,"Q_Integrator"))
########################################################################
S_Reach <- t(ncvar_get(Reach_in,"S"))
n_Reach <- t(ncvar_get(Reach_in,"n_Integrator"))  ## doesn't change with time
Ao_Reach <- t(ncvar_get(Reach_in,"Ab_Intgrator"))## doesn't change with time
n_BAM_Reach <- t(ncvar_get(Reach_in,"n"))  ## doesn't change with time
Ao_BAM_Reach <- t(ncvar_get(Reach_in,"Abase"))## doesn't change with time
W_Reach <- t(ncvar_get(Reach_in,"W"))
dA_Reach <- t(ncvar_get(Reach_in,"dA"))
H_Reach <- t(ncvar_get(Reach_in,"H")) 
HydroID_Reach<- ncvar_get(Reach_in,"HydroID")

nc_close(Reach_in)

nReaches=dim(W_Reach)[2]

for (i in 1:nReaches){
  RelInitialerror=abs(Q_Integrator[,i]-Qpred[,i])/Qpred[,i]
  Rel_Error_threshold <- 0.75
  
  
  # Initial error check
  QReach_Rel_errorhigh= RelInitialerror > Rel_Error_threshold
  Timestamp_REL_error= which(QReach_Rel_errorhigh, arr.ind = FALSE, useNames = TRUE)
  
  ## Comparing pre and post integrator flows #################################################################
  
  #calculating the inital error ~ DIFFERENCE BETWEEN aLGORITHM dISCHARGE AND iNTEGRATOR dISCharge
  
  Priors_error=abs(Q_Integrator[,i]-Qstar[,i])/Qstar[,i]
  Prior_Error_threshold <- 0.75
  
  # Comparing to the priors
  
  QReach_Prior_errorhigh= Priors_error > Prior_Error_threshold
  Timestamp_PRIOR_error= which(QReach_Prior_errorhigh, arr.ind = FALSE, useNames = TRUE)
  
  #Qpred[Timestamp_REL_error]= -77777
  Q_Integrator[Timestamp_REL_error,i]= -77777
  
  Q_Integrator[Timestamp_PRIOR_error,i]= -88888
  
  GridID_Reach=routing_table[which(routing_table$HydroID == HydroID_Reach[i]),c(1)]
  reachname=GridID_Reach
  filename=paste("SevernRiver_Reach_",reachname,sep="")
  
  
  #-- writing the modified netcdf files -----------------
  ## We need to follow the structure of the output file we get tomorrow
  #---------- Creating the Netcdf file ------------------
  ncpath <- paste(getwd(), '/output_files/', sep = "")
  ncname <- paste(filename,sep="")  
  ncfname <- paste(ncpath, ncname, ".nc", sep="")
  
  #---------- define variables-------------------------------
 
  time <- as.array(seq(1,dim(W_Reach)[1],1))
  Reach <- as.array(seq(1,1,1))
  
 
  Timestepsdim <- ncdim_def("Time steps","day",as.integer(time))
  Reachdim <- ncdim_def("Reach","RiverReaches",as.integer(Reach))
  
  fillvalue <- -9999
 
  dlname <- "Manning's Roughness coefficient_BAM"
  n_BAM <- ncvar_def("n_BAM","dimensionless",list(Reachdim,Reachdim),fillvalue,dlname,prec="double")
  dlname <- "Manning's Roughness coefficient_Integrator"
  n_Int <- ncvar_def("n_Int","dimensionless",list(Reachdim,Reachdim),fillvalue,dlname,prec="double")
  dlname <- "Average Water SUrface Area_BAM"
  Ao_BAM <- ncvar_def("Ao_BAM","meters squared",list(Reachdim,Reachdim),fillvalue,dlname,prec="double")
  dlname <- "Average Water SUrface Area_Integrator"
  Ao_Int <- ncvar_def("Ao_Int","meters squared",list(Reachdim,Reachdim),fillvalue,dlname,prec="double")
  dlname <- "Water Surface Elevation"# should we average per reach?
  H <- ncvar_def("H","meters",list(Reachdim,Timestepsdim),fillvalue,dlname,prec="double")
  dlname <- "River Width"
  w <- ncvar_def("w","meters",list(Reachdim,Timestepsdim),fillvalue,dlname,prec="double")
  dlname <- "Reach Slope"
  S <- ncvar_def("Reach_Timeseries/S_1km","m/m",list(Reachdim,Timestepsdim),fillvalue,dlname,prec="double")
  dlname <- "Difference from average Water surface area"
  dA <- ncvar_def("dA","meters squared",list(Reachdim,Timestepsdim),fillvalue,dlname,prec="double")
  
  dlname <- "Reach Discharge_Prior_mean"
  Q_mean <- ncvar_def("Q_mean","cubic meters per second",list(Reachdim,Timestepsdim),fillvalue,dlname,prec="double")
  dlname <- "Reach Discharge_BAM"
  Q_BAM <- ncvar_def("Q_BAM","cubic meters per second",list(Reachdim,Timestepsdim),fillvalue,dlname,prec="double")
  dlname <- "Reach Discharge_Integrator"
  Q_Int <- ncvar_def("Q_Int","cubic meters per second",list(Reachdim,Timestepsdim),fillvalue,dlname,prec="double")
  
  #--------- create netCDF file and put arrays---------------------
  ncout <- nc_create(ncfname,list(n_BAM,n_Int,Ao_BAM,Ao_Int,H,w,S,dA,Q_mean,Q_BAM,Q_Int),force_v4=TRUE)
  
  #---------------- Insert the variables 
  ncvar_put(ncout,n_BAM,n_BAM_Reach[,i])
  ncvar_put(ncout,n_Int,n_Reach[,i])
  ncvar_put(ncout,Ao_BAM,Ao_BAM_Reach[,i])
  ncvar_put(ncout,Ao_Int,Ao_Reach[,i])
  ncvar_put(ncout,H,H_Reach[,i])
  ncvar_put(ncout,w,W_Reach[,i])
  ncvar_put(ncout,S,S_Reach[,i])
  ncvar_put(ncout,dA,dA_Reach[,i])
  ncvar_put(ncout,Q_mean,Qstar[,i])
  ncvar_put(ncout,Q_BAM,Qpred[,i])
  ncvar_put(ncout,Q_Int,Q_Integrator[,i])
  
  #---------------- add global attributes
  title=paste("Severn SWOT-like data for reach", reachname)
  ncatt_put(ncout,0,"title",title)
  
  nc_close(ncout)
  
  
}