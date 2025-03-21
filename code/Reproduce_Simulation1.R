##################################################################
#### Code for reproduce the numerical results of Simulation 1 ####
##################################################################
#### Preliminary ####
# Load functions
sapply(c("code/EM_MLE.R","code/EM_CHMM.R","code/EM_FHMM.R","code/EM_Oracle.R"), source)

# Required library to load data 
library(httr)
library(parallel) # Use multicore.
core_num <- 25 # Set number of cores.

# Base directory of data
home_dir <- "https://raw.githubusercontent.com/JS-Bok/Clustered-hidden-Markov-models/main/data/"

#### Model 1 ####
# Raw URL of the file
file_url <- paste0(home_dir,"Sim1_Model1.Rdata")

# Load data
temp_file <- tempfile(fileext = ".Rdata")
response <- GET(file_url)
writeBin(content(response, as = "raw"), temp_file)
load(temp_file)
rm(response)
unlink(temp_file)

# Define a function generating numerical results across OE, MLE, MPLE
model_measurements <- function(data, dist_class, A_true, phi_true){
  OE_result <- EM_Oracle_Sim1_A1(data,dist_class = dist_class, m=10)
  MLE_result <- EM_MLE(data,dist_class = dist_class, m=10)
  MPLE_result <- EM_MPLE(data,dist_class = dist_class, m=10)
  OE_val <- data.frame(TE=n2_sqrt(A_true-OE_result$transition),
                       EE=n2_sqrt(unlist(phi_true)-unlist(OE_result$par)),
                       TP=TFPN(cf(A_true)$index,
                               cf(OE_result$transition)$index)$TP,
                       FP=TFPN(cf(A_true)$index,
                               cf(OE_result$transition)$index)$FP)
  MLE_val <- data.frame(TE=n2_sqrt(A_true-MLE_result$transition),
                       EE=n2_sqrt(unlist(phi_true)-unlist(MLE_result$par)),
                       TP=TFPN(cf(A_true)$index,
                               cf(MLE_result$transition)$index)$TP,
                       FP=TFPN(cf(A_true)$index,
                               cf(MLE_result$transition)$index)$FP)
  MPLE_val <- data.frame(TE=n2_sqrt(A_true-MPLE_result$transition),
                       EE=n2_sqrt(unlist(phi_true)-unlist(MPLE_result$par)),
                       TP=TFPN(cf(A_true)$index,
                               cf(MPLE_result$transition)$index)$TP,
                       FP=TFPN(cf(A_true)$index,
                               cf(MPLE_result$transition)$index)$FP)
  return(list(OE=OE_val, MLE=MLE_val, MPLE=MPLE_val))
}

# Define true parameters
A_true <- matrix(c(1/5, 1/5, 1/5, 1/5, 1/30, 1/30, 1/30, 1/30, 1/30, 1/30, 
                   1/5, 1/5, 1/5, 1/5, 1/30, 1/30, 1/30, 1/30, 1/30, 1/30, 
                   1/5, 1/5, 1/5, 1/5, 1/30, 1/30, 1/30, 1/30, 1/30, 1/30, 
                   1/30, 1/30, 1/30, 1/5, 1/5, 1/5, 1/5, 1/30, 1/30, 1/30, 
                   1/30, 1/30, 1/30, 1/5, 1/5, 1/5, 1/5, 1/30, 1/30, 1/30, 
                   1/30, 1/30, 1/30, 1/5, 1/5, 1/5, 1/5, 1/30, 1/30, 1/30, 
                   1/30, 1/30, 1/30, 1/5, 1/5, 1/5, 1/5, 1/30, 1/30, 1/30, 
                   1/30, 1/30, 1/30, 1/30, 1/30, 1/30, 1/5, 1/5, 1/5, 1/5, 
                   1/30, 1/30, 1/30, 1/30, 1/30, 1/30, 1/5, 1/5, 1/5, 1/5, 
                   1/30, 1/30, 1/30, 1/30, 1/30, 1/30, 1/5, 1/5, 1/5, 1/5 
),nrow = 10, byrow = TRUE)

dist_class <- "norm"
mean <- seq(from = -400, by = 100, length.out = 10)
sd <- rep(25, times=10)
phi_true <- mapply(function(xx,yy){list(mean=xx, sigma=matrix(yy))}, xx=mean, yy=sd, SIMPLIFY = FALSE)

# Calculate model measurements
Model1_result <- mclapply(Sim1_Model1, function(x){
  result_4000 <- model_measurements(x[1:4000], dist_class=dist_class, A_true=A_true, phi_true=phi_true)
  result_5000 <- model_measurements(x[1:5000], dist_class=dist_class, A_true=A_true, phi_true=phi_true)
  result_6000 <- model_measurements(x[1:6000], dist_class=dist_class, A_true=A_true, phi_true=phi_true)
  return(list(n_4000=result_4000, n_5000=result_5000,n_6000=result_6000))
}, mc.cores = core_num)











#### Model 2 ####
# Raw URL of the file
file_url <- paste0(home_dir,"Sim1_Model2.Rdata")

# Load data
temp_file <- tempfile(fileext = ".Rdata")
response <- GET(file_url)
writeBin(content(response, as = "raw"), temp_file)
load(temp_file)
rm(response)
unlink(temp_file)

# Define a function generating numerical results across OE, MLE, MPLE
model_measurements <- function(data, dist_class, A_true, phi_true){
  OE_result <- EM_Oracle_Sim1_A2(data,dist_class = dist_class, m=10)
  MLE_result <- EM_MLE(data,dist_class = dist_class, m=10)
  MPLE_result <- EM_MPLE(data,dist_class = dist_class, m=10)
  OE_val <- data.frame(TE=n2_sqrt(A_true-OE_result$transition),
                       EE=n2_sqrt(unlist(phi_true)-unlist(OE_result$par)),
                       TP=TFPN(cf(A_true)$index,
                               cf(OE_result$transition)$index)$TP,
                       FP=TFPN(cf(A_true)$index,
                               cf(OE_result$transition)$index)$FP)
  MLE_val <- data.frame(TE=n2_sqrt(A_true-MLE_result$transition),
                        EE=n2_sqrt(unlist(phi_true)-unlist(MLE_result$par)),
                        TP=TFPN(cf(A_true)$index,
                                cf(MLE_result$transition)$index)$TP,
                        FP=TFPN(cf(A_true)$index,
                                cf(MLE_result$transition)$index)$FP)
  MPLE_val <- data.frame(TE=n2_sqrt(A_true-MPLE_result$transition),
                         EE=n2_sqrt(unlist(phi_true)-unlist(MPLE_result$par)),
                         TP=TFPN(cf(A_true)$index,
                                 cf(MPLE_result$transition)$index)$TP,
                         FP=TFPN(cf(A_true)$index,
                                 cf(MPLE_result$transition)$index)$FP)
  return(list(OE=OE_val, MLE=MLE_val, MPLE=MPLE_val))
}

# Define true parameters
A_true <- matrix(c(1/3, 1/3, 1/24, 1/24, 1/24, 1/24, 1/24, 1/24, 1/24, 1/24,
                   1/3, 1/3, 1/24, 1/24, 1/24, 1/24, 1/24, 1/24, 1/24, 1/24,
                   1/24, 1/24, 1/3, 1/3, 1/24, 1/24, 1/24, 1/24, 1/24, 1/24,
                   1/24, 1/24, 1/3, 1/3, 1/24, 1/24, 1/24, 1/24, 1/24, 1/24,
                   1/24, 1/24, 1/24, 1/24, 1/3, 1/3, 1/24, 1/24, 1/24, 1/24,
                   1/24, 1/24, 1/24, 1/24, 1/3, 1/3, 1/24, 1/24, 1/24, 1/24,
                   1/24, 1/24, 1/24, 1/24, 1/24, 1/24, 1/3, 1/3, 1/24, 1/24,
                   1/24, 1/24, 1/24, 1/24, 1/24, 1/24, 1/3, 1/3, 1/24, 1/24,
                   1/24, 1/24, 1/24, 1/24, 1/24, 1/24, 1/24, 1/24, 1/3, 1/3,
                   1/24, 1/24, 1/24, 1/24, 1/24, 1/24, 1/24, 1/24, 1/3, 1/3
),nrow = 10, byrow = TRUE)

dist_class <- "norm"
mean <- seq(from = -400, by = 100, length.out = 10)
sd <- rep(25, times=10)
phi_true <- mapply(function(xx,yy){list(mean=xx, sigma=matrix(yy))}, xx=mean, yy=sd, SIMPLIFY = FALSE)


# Calculate model measurements
Model2_result <- mclapply(Sim1_Model2, function(x){
  result_4000 <- model_measurements(x[1:4000], dist_class=dist_class, A_true=A_true, phi_true=phi_true)
  result_5000 <- model_measurements(x[1:5000], dist_class=dist_class, A_true=A_true, phi_true=phi_true)
  result_6000 <- model_measurements(x[1:6000], dist_class=dist_class, A_true=A_true, phi_true=phi_true)
  return(list(n_4000=result_4000, n_5000=result_5000,n_6000=result_6000))
}, mc.cores = core_num)


#### Model 3 ####
# Raw URL of the file
file_url <- paste0(home_dir,"Sim1_Model3.Rdata")

# Load data
temp_file <- tempfile(fileext = ".Rdata")
response <- GET(file_url)
writeBin(content(response, as = "raw"), temp_file)
load(temp_file)
rm(response)
unlink(temp_file)

# Define a function generating numerical results across OE, MLE, MPLE
model_measurements <- function(data, dist_class, A_true, phi_true){
  OE_result <- EM_Oracle_Sim1_A1(data,dist_class = dist_class, m=10)
  MLE_result <- EM_MLE(data,dist_class = dist_class, m=10)
  MPLE_result <- EM_MPLE(data,dist_class = dist_class, m=10)
  OE_val <- data.frame(TE=n2_sqrt(A_true-OE_result$transition),
                       EE=n2_sqrt(unlist(phi_true)-unlist(OE_result$par)),
                       TP=TFPN(cf(A_true)$index,
                               cf(OE_result$transition)$index)$TP,
                       FP=TFPN(cf(A_true)$index,
                               cf(OE_result$transition)$index)$FP)
  MLE_val <- data.frame(TE=n2_sqrt(A_true-MLE_result$transition),
                        EE=n2_sqrt(unlist(phi_true)-unlist(MLE_result$par)),
                        TP=TFPN(cf(A_true)$index,
                                cf(MLE_result$transition)$index)$TP,
                        FP=TFPN(cf(A_true)$index,
                                cf(MLE_result$transition)$index)$FP)
  MPLE_val <- data.frame(TE=n2_sqrt(A_true-MPLE_result$transition),
                         EE=n2_sqrt(unlist(phi_true)-unlist(MPLE_result$par)),
                         TP=TFPN(cf(A_true)$index,
                                 cf(MPLE_result$transition)$index)$TP,
                         FP=TFPN(cf(A_true)$index,
                                 cf(MPLE_result$transition)$index)$FP)
  return(list(OE=OE_val, MLE=MLE_val, MPLE=MPLE_val))
}

# Define true parameters
A_true <- matrix(c(1/5, 1/5, 1/5, 1/5, 1/30, 1/30, 1/30, 1/30, 1/30, 1/30, 
                   1/5, 1/5, 1/5, 1/5, 1/30, 1/30, 1/30, 1/30, 1/30, 1/30, 
                   1/5, 1/5, 1/5, 1/5, 1/30, 1/30, 1/30, 1/30, 1/30, 1/30, 
                   1/30, 1/30, 1/30, 1/5, 1/5, 1/5, 1/5, 1/30, 1/30, 1/30, 
                   1/30, 1/30, 1/30, 1/5, 1/5, 1/5, 1/5, 1/30, 1/30, 1/30, 
                   1/30, 1/30, 1/30, 1/5, 1/5, 1/5, 1/5, 1/30, 1/30, 1/30, 
                   1/30, 1/30, 1/30, 1/5, 1/5, 1/5, 1/5, 1/30, 1/30, 1/30, 
                   1/30, 1/30, 1/30, 1/30, 1/30, 1/30, 1/5, 1/5, 1/5, 1/5, 
                   1/30, 1/30, 1/30, 1/30, 1/30, 1/30, 1/5, 1/5, 1/5, 1/5, 
                   1/30, 1/30, 1/30, 1/30, 1/30, 1/30, 1/5, 1/5, 1/5, 1/5 
),nrow = 10, byrow = TRUE)

dist_class <- "mvnorm"

phi_true=list()
mean_mat <- matrix(c(0, 5, 10, 15, 20, 25, 30, 35, 40, 45,
                     0, 10, 20, 30, 40, 50, 60, 70, 80, 90,
                     0, 15, 30, 45, 60, 75, 90, 105, 120, 135), ncol = 3 )
for (i in c(1:10)){
  phi_true[[i]] <- list(mean = as.vector(mean_mat[i,]), sigma=matrix(c(1, 0.1, 0.1,
                                                                       0.1, 2^2, 0.1,
                                                                       0.1, 0.1, 3^2), ncol = 3))
}
# Calculate model measurements
Model3_result <- mclapply(Sim1_Model3, function(x){
  result_4000 <- model_measurements(x[1:4000,], dist_class=dist_class, A_true=A_true, phi_true=phi_true)
  result_5000 <- model_measurements(x[1:5000,], dist_class=dist_class, A_true=A_true, phi_true=phi_true)
  result_6000 <- model_measurements(x[1:6000,], dist_class=dist_class, A_true=A_true, phi_true=phi_true)
  return(list(n_4000=result_4000, n_5000=result_5000,n_6000=result_6000))
}, mc.cores = core_num)



#### Model 4 ####
# Raw URL of the file
file_url <- paste0(home_dir,"Sim1_Model4.Rdata")

# Load data
temp_file <- tempfile(fileext = ".Rdata")
response <- GET(file_url)
writeBin(content(response, as = "raw"), temp_file)
load(temp_file)
rm(response)
unlink(temp_file)

# Define a function generating numerical results across OE, MLE, MPLE
model_measurements <- function(data, dist_class, A_true, phi_true){
  OE_result <- EM_Oracle_Sim1_A2(data,dist_class = dist_class, m=10)
  MLE_result <- EM_MLE(data,dist_class = dist_class, m=10)
  MPLE_result <- EM_MPLE(data,dist_class = dist_class, m=10)
  OE_val <- data.frame(TE=n2_sqrt(A_true-OE_result$transition),
                       EE=n2_sqrt(unlist(phi_true)-unlist(OE_result$par)),
                       TP=TFPN(cf(A_true)$index,
                               cf(OE_result$transition)$index)$TP,
                       FP=TFPN(cf(A_true)$index,
                               cf(OE_result$transition)$index)$FP)
  MLE_val <- data.frame(TE=n2_sqrt(A_true-MLE_result$transition),
                        EE=n2_sqrt(unlist(phi_true)-unlist(MLE_result$par)),
                        TP=TFPN(cf(A_true)$index,
                                cf(MLE_result$transition)$index)$TP,
                        FP=TFPN(cf(A_true)$index,
                                cf(MLE_result$transition)$index)$FP)
  MPLE_val <- data.frame(TE=n2_sqrt(A_true-MPLE_result$transition),
                         EE=n2_sqrt(unlist(phi_true)-unlist(MPLE_result$par)),
                         TP=TFPN(cf(A_true)$index,
                                 cf(MPLE_result$transition)$index)$TP,
                         FP=TFPN(cf(A_true)$index,
                                 cf(MPLE_result$transition)$index)$FP)
  return(list(OE=OE_val, MLE=MLE_val, MPLE=MPLE_val))
}

# Define true parameters
A_true <- matrix(c(1/3, 1/3, 1/24, 1/24, 1/24, 1/24, 1/24, 1/24, 1/24, 1/24,
                   1/3, 1/3, 1/24, 1/24, 1/24, 1/24, 1/24, 1/24, 1/24, 1/24,
                   1/24, 1/24, 1/3, 1/3, 1/24, 1/24, 1/24, 1/24, 1/24, 1/24,
                   1/24, 1/24, 1/3, 1/3, 1/24, 1/24, 1/24, 1/24, 1/24, 1/24,
                   1/24, 1/24, 1/24, 1/24, 1/3, 1/3, 1/24, 1/24, 1/24, 1/24,
                   1/24, 1/24, 1/24, 1/24, 1/3, 1/3, 1/24, 1/24, 1/24, 1/24,
                   1/24, 1/24, 1/24, 1/24, 1/24, 1/24, 1/3, 1/3, 1/24, 1/24,
                   1/24, 1/24, 1/24, 1/24, 1/24, 1/24, 1/3, 1/3, 1/24, 1/24,
                   1/24, 1/24, 1/24, 1/24, 1/24, 1/24, 1/24, 1/24, 1/3, 1/3,
                   1/24, 1/24, 1/24, 1/24, 1/24, 1/24, 1/24, 1/24, 1/3, 1/3
),nrow = 10, byrow = TRUE)

dist_class <- "mvnorm"

phi_true=list()
mean_mat <- matrix(c(0, 5, 10, 15, 20, 25, 30, 35, 40, 45,
                     0, 10, 20, 30, 40, 50, 60, 70, 80, 90,
                     0, 15, 30, 45, 60, 75, 90, 105, 120, 135), ncol = 3 )
for (i in c(1:10)){
  phi_true[[i]] <- list(mean = as.vector(mean_mat[i,]), sigma=matrix(c(1, 0.1, 0.1,
                                                                       0.1, 2^2, 0.1,
                                                                       0.1, 0.1, 3^2), ncol = 3))
}
# Calculate model measurements
Model4_result <- mclapply(Sim1_Model4[1:2], function(x){
  result_4000 <- model_measurements(x[1:4000,], dist_class=dist_class, A_true=A_true, phi_true=phi_true)
  result_5000 <- model_measurements(x[1:5000,], dist_class=dist_class, A_true=A_true, phi_true=phi_true)
  result_6000 <- model_measurements(x[1:6000,], dist_class=dist_class, A_true=A_true, phi_true=phi_true)
  return(list(n_4000=result_4000, n_5000=result_5000,n_6000=result_6000))
}, mc.cores = core_num)










