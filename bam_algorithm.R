library(doParallel)
library(ncdf4)
library(bamr)
library(dplyr)
library(doParallel)
library(foreach)
library(easyNCDF)
computation_function <- function(filename) {
  library(doParallel)
  library(ncdf4)
  library(bamr)
  library(dplyr)
  library(doParallel)
  library(foreach)
  library(easyNCDF)

  calcdA_mat <- function(w, h) {
    stopifnot(all(dim(w) == dim(h)))
    dA <- w
    for (i in 1:nrow(dA)) {
      dA[i, ] <- calcdA_vec(w[i, ], h[i, ])
    }
    dA
  }

  calcdA_vec <- function(w, h) {
    words <- order(w)
    warr <- w[words]
    harr <- h[words]
    delh <- c(0, diff(harr))
    delA <- cumsum(warr * delh)
    dA <- 1:length(w)
    dA[words] <- delA
    dA
  }

  data_in <- nc_open(filename)
  W_obs <- ncvar_get(data_in, "XS_Timseries/W")
  Q_obs <- ncvar_get(data_in, "XS_Timseries/Q")
  H_obs <- ncvar_get(data_in, "XS_Timseries/H_1km")
  S_obs <- ncvar_get(data_in, "Reach_Timeseries/S_1km")
  HydroID = ncvar_get(data_in, "Reach_Timeseries/HydroID")
  NextDownID = ncvar_get(data_in, "Reach_Timeseries/NextDownID")
  GridID = ncvar_get(data_in, "Reach_Timeseries/GridID")


  number_of_days = 690
  W_obs <- W_obs[ , 320:number_of_days, drop=FALSE]
  H_obs <- H_obs[ , 320:number_of_days , drop=FALSE]
  Q_obs <- Q_obs[ , 320:number_of_days , drop=FALSE]
  S_obs <- S_obs[320:number_of_days]

  # make a copy of W_obs, H_obs and S_obs
  S_obs[S_obs==0]=1e-6
  H_obs1 = H_obs
  W_obs1 = W_obs
  S_obs1 = S_obs
  Q_obs1 = Q_obs
  S_obs1 = matrix(unlist(S_obs1), ncol = length(S_obs1), byrow = TRUE)


  S_obs[S_obs==0]=1e-6
  S_obs[S_obs<0]=NA
  W_obs[W_obs<0]=NA
  H_obs[H_obs<0]=NA
  Q_obs[Q_obs<0]=NA



  W_obs= W_obs[, colSums(is.na(W_obs) ) < nrow(W_obs)]
  W_obs[W_obs==0]=1e-6

  # Mariam Edits : 10/03/2019
  H_obs= H_obs[, colSums(is.na(H_obs) ) < nrow(H_obs)]
  H_obs[H_obs==0]=1e-6

  Q_obs= Q_obs[, colSums(is.na(Q_obs) ) < nrow(Q_obs)]
  Q_obs[Q_obs==0]=1e-6

  # convert S_obs to matrix
  S_obs <- matrix(unlist(S_obs),ncol = length(S_obs), byrow = TRUE)
  S_obs[S_obs==0]=1e-6
  S_obs <- S_obs[rep(1:nrow(S_obs), times=nrow(W_obs)), ]
  S_obs= S_obs[, colSums(is.na(S_obs) ) < nrow(S_obs)]

  dA=calcdA_mat(W_obs,H_obs)

  Q_prior <- abs(mean(Q_obs,na.rm=TRUE))

  bamdata=bamr::bam_data(max_xs=nrow(W_obs), w=W_obs,dA=dA,s=S_obs,Qhat=as.vector(Q_prior))

  priors=bam_priors(bamdata,logQ_sd = cv2sigma(1))

  priors$lowerbound_logQc=0.01
  priors$upperbound_logn=log(0.05)

  # Mariam edits: 10/03/2019
  priors$lowerbound_logQ=min(na.rm=TRUE, log(Q_obs))
  priors$upperbound_logQ=max(na.rm=TRUE ,log(Q_obs))

  if (max(dA,na.rm=TRUE)<2000){
    priors$lowerbound_A0= 1
    priors$upperbound_A0= 1e3
  }

  #==================================================
  #Qbounds
  if (priors$logWc_hat == -Inf){
    priors$logWc_hat=mean(log(W_obs),na.rm=TRUE) #in case of error
  }

  if (priors$lowerbound_logQ == -Inf){
    priors$lowerbound_logQ=  priors$logQc_hat- (priors$upperbound_logQ-priors$logQc_hat) #in case of error
  }

  if (  priors$upperbound_logQ == -Inf){
    priors$upperbound_logQ=  priors$logQc_hat+ (priors$logQc_hat-priors$lowerbound_logQ) #in case of error
  }

  estimates <- bam_estimate(bamdata=bamdata,bampriors=priors, variant = 'manning_amhg', iter=5)
  required_data <- as.data.frame(estimates)

  A0 <-mean(rowMeans(unname(required_data[, grep('^A0', names(required_data))])))
  n <- mean(as.numeric(required_data[, grep('^logn', names(required_data))]))
  n <- exp(n)

  dA <- as.matrix(colMeans(dA))
  S1 <- t(S_obs1)
  W1 <- as.matrix(colMeans(W_obs1))
  H1 <- as.matrix(colMeans(H_obs1))

  empty_indices_of_H1 <- which(is.na(H1))
  not_non_indices_of_H1 <- which(!is.na(H1))
  # create an empty matrix withh same dimensions as H
  empty_matrix = t(matrix(0L, nrow = nrow(H1), ncol = ncol(H1)))
  partially_filled_mat <- replace(empty_matrix, empty_indices_of_H1, NA)

  predicted_q <- required_data[, grep('^logQ_amhg', names(required_data))]
  predicted_q <- colMeans(predicted_q)
  predicted_q <- matrix(data=predicted_q, ncol=ncol(Q_obs), nrow=nrow(Q_obs))
  #predicted_q <- as.matrix(colMeans(predicted_q))
  predicted_q <- mean(colMeans(predicted_q))
  Q_obs1 = colMeans(Q_obs1)
  Q_obs1 = t(as.matrix(Q_obs1))
  #Q_obs1 = mean(Q_obs1, na.rm = TRUE)
  predicted_q <- (replace(partially_filled_mat, not_non_indices_of_H1, predicted_q))
  #predicted_q <- mean(predicted_q, na.rm = TRUE)
  dA <- t(replace(partially_filled_mat, not_non_indices_of_H1, dA))
  current_result = list("A0"=A0, "n"=n, "S"=S1, "H"=H1, "dA"=dA, "W"=W1, "PredQ"=predicted_q, 'Q_obs'=Q_obs1, "HydroID"=HydroID, "NextDownID"=NextDownID, "GridID"=GridID)
  return(current_result)
}

# prepare the results for storing into a netcdf file.
format_estimated_results <- function(results, output_filename) {

  output_filename <- paste(getwd(), "/bam_output.nc", sep="")
  routing_table_name <- paste(getwd(), "/routing_table.csv", sep="")

  Abase <- as.array(as.numeric(results[1,]))
  n <- as.array(as.numeric(results[2,]))
  S <- results[3, ]
  H <- results[4, ]
  dA <- results[5, ]
  W <- results[6, ]
  pred_q <- results[7, ]
  Q_star <- results[8, ]
  HydroID <- as.array(as.numeric(results[9, ]))
  NextDownId <- as.array(as.numeric(results[10, ]))
  GridId <- as.array(results[11, ])

  dat = list(GridId, HydroID, NextDownId)
  dat = unlist(dat)
  dat = matrix(dat, ncol=3)
  colnames(dat) = c("GridID", "HydroID", "NextDownID")
  dat <- as.data.frame(dat)
  dat$HydroID <- as.numeric(as.character(dat$HydroID))
  dat$NextDownID <- as.numeric(as.character(dat$NextDownID))

  write.csv(dat, routing_table_name, row.names = FALSE)

  pred_q <- as.data.frame(pred_q)
  row.names(pred_q)<- NULL ; colnames(pred_q)<- NULL
  #pred_q = t(pred_q)
  #pred_q <- array(pred_q, dim = c(nrow(pred_q), ncol(pred_q)))
  pred_q <- array(unlist(pred_q), dim = c(length(n), dim(pred_q)[2]/length(n)))
  names(dim(pred_q)) <- c('predqRow','predqCol')

  qStar <- as.data.frame(Q_star)
  row.names(qStar)<- NULL ; colnames(qStar)<- NULL
  #qStar = t(qStar)
  #qStar <- array(qStar, dim = c(nrow(qStar), ncol(qStar)))
  qStar <- array(unlist(qStar), dim = c(length(n), dim(qStar)[2]/length(n)))
  names(dim(qStar)) <- c('qStarRow','qStarCol')

  H <- as.data.frame(H)
  row.names(H)<- NULL ; colnames(H)<- NULL
  H = t(H)
  H<- array(H, dim = c(nrow(H), ncol(H)))
  names(dim(H)) <- c('HRow','HCol')

  S <- as.data.frame(S)
  row.names(S)<- NULL ; colnames(S)<- NULL
  S <- t(S)
  S <- array(S, dim = c(nrow(S), ncol(S)))
  names(dim(S)) <- c('SRow','SCol')

  dA <- as.data.frame(dA)
  row.names(dA) <- NULL ; colnames(dA) <- NULL
  dA <- t(dA)
  dA <- array(dA, dim = c(nrow(dA), ncol(dA)))
  names(dim(dA)) <- c('dARow','dACol')


  W <- as.data.frame(W)
  row.names(W) <- NULL ; colnames(W) <- NULL
  W <- t(W)
  W <- array(W, dim=c(nrow(W), ncol(W)))
  names(dim(W)) <- c('WRow','WCol')

  names(dim(Abase)) <- c("Abase")
  names(dim(n)) <- c("n")
  names(dim(HydroID)) <- c("HydroID")
  names(dim(GridId)) <- c("GridId")

  names(dim(S)) < c(c("S"), dim(S)[2])
  names(dim(dA)) < c(c("dA"), dim(dA)[2])
  names(dim(W)) < c(c("W"), dim(W)[2])
  names(dim(H)) < c(c("H"), dim(H)[2])
  names(dim(pred_q)) < c(c("Qpred"), dim(pred_q)[2])
  names(dim(qStar)) < c(c("QStar"), dim(qStar)[2])
  ArrayToNc(list("Abase"=Abase, "n"=n, "S"=S, "dA"=dA, "W"=W, "H"=H, "Qpred"=pred_q, "Qstar"=qStar, "HydroID"=HydroID), output_filename)
}

#setwd("/Users/amuhebwa/Documents/Python/BAM_ALGORITHM_FINAL")
filepath <- paste(getwd(), "/data", sep="")
all_files <- list.files(filepath, pattern="*.nc", full.names=TRUE)
cl <- makeCluster(future::availableCores())
registerDoParallel()
doFuture::registerDoFuture()

start_time= Sys.time()
results <- foreach(filename=all_files, .combine = cbind, .errorhandling = 'remove', .packages = c("ncdf4", "doParallel","foreach", "bamr", "dplyr", "parallel", "easyNCDF")) %dopar% computation_function(filename)
stopCluster(cl)
end_time = Sys.time()
elapsed_time=end_time-start_time
print(elapsed_time)
# drop NaNs resulting from errors in some of the files!
format_estimated_results(results, output_filename)
