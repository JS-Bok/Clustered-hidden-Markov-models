##################################################################
#### Code for reproduce the numerical results of Simulation 2 ####
##################################################################
#### Preliminary ####
# Load functions
sapply(c("code/EM_MLE.R","code/EM_MPLE.R","code/EM_FHMM.R","code/EM_Oracle.R", "code/EM_MPLE-GroupLASSO.R"), source)

# Required library to load data 
library(parallel) # Use multicore.
core_num <- 100 # Set number of cores.

# Seeds for data generation
set.seed(2026, kind = "L'Ecuyer-CMRG")
Model1_seed <- sample(1:1e8,500)
Model2_seed <- sample(1:1e8,500)
Model3_seed <- sample(1:1e8,500)
Model4_seed <- sample(1:1e8,500)

#### Model 1 ####
# Data generation
m <- 9

data_dim <- 1

# Transition matrix
A1_real <- diag(3)+0.2
A1_real <- A1_real / rowSums(A1_real)
A2_real <- matrix(1, ncol=3, nrow=3)+0.2
A2_real <- A2_real / rowSums(A2_real)

A_true <- kronecker(A1_real, A2_real)

# Emission parameters
dist_class <- "norm"
mean <- seq(from = -400, by = 100, length.out = 9)
sd <- rep(25, times=9)
phi_true <- mapply(function(xx,yy){list(mean=xx, sigma=matrix(yy^2))}, xx=mean, yy=sd, SIMPLIFY = FALSE)


data <- mclapply(Model1_seed,function(x){
  set.seed(x, kind = "L'Ecuyer-CMRG")
  HMM_simulation(dimension = data_dim,
                 transition = A_true, 
                 emission = phi_true, 
                 size = 10000)
},mc.cores = core_num)

dist_class <- "norm"
mean <- seq(from = -400, by = 100, length.out = 9)
sd <- rep(25, times=9)
phi_true <- mapply(function(xx,yy){list(mean=xx, sigma=matrix(yy))}, xx=mean, yy=sd, SIMPLIFY = FALSE)


# Define a function generating numerical results across OE, MLE, MPLE
model_measurements <- function(data, dist_class, A_true, phi_true){
  OE_result <- EM_Oracle_Sim2_A1(data,dist_class = dist_class, m=9)
  MLE_result <- EM_MLE(data,dist_class = dist_class, m=9)
  MPLE_GLASSO_result <- EM_MPLE_GLASSO(data,dist_class = dist_class, m=9)
  MPLE_result <- EM_MPLE(data,dist_class = dist_class, m=9)
  FHMM_result <- fhmm(data, M=2, K=3, cyc=100)
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
  MPLE_GLASSO_val <- data.frame(TE=n2_sqrt(A_true-MPLE_GLASSO_result$transition),
                                EE=n2_sqrt(unlist(phi_true)-unlist(MPLE_GLASSO_result$par)),
                                TP=TFPN(cf(A_true)$index,
                                        cf(MPLE_GLASSO_result$transition)$index)$TP,
                                FP=TFPN(cf(A_true)$index,
                                        cf(MPLE_GLASSO_result$transition)$index)$FP)
  MPLE_val <- data.frame(TE=n2_sqrt(A_true-MPLE_result$transition),
                         EE=n2_sqrt(unlist(phi_true)-unlist(MPLE_result$par)),
                         TP=TFPN(cf(A_true)$index,
                                 cf(MPLE_result$transition)$index)$TP,
                         FP=TFPN(cf(A_true)$index,
                                 cf(MPLE_result$transition)$index)$FP)
  
  fhmm_par <- FHMM_result$par
  
  fhmm_A <- kronecker(FHMM_result$P[1:3,], FHMM_result$P[4:6,])
  
  FHMM_val <- data.frame(TE=n2_sqrt(A_true-fhmm_A),
                         EE=n2_sqrt(unlist(phi_true)-unlist(fhmm_par)),
                         TP=TFPN(cf(A_true)$index,
                                 cf(fhmm_A)$index)$TP,
                         FP=TFPN(cf(A_true)$index,
                                 cf(fhmm_A)$index)$FP)
  
  return(list(OE=OE_val, MLE=MLE_val, MPLE_GLASSO=MPLE_GLASSO_val, MPLE=MPLE_val, FHMM=FHMM_val))
}

# Calculate model measurements
Model1_result <- mclapply(data, function(x){
  result_5000 <- model_measurements(x[1:5000,], dist_class=dist_class, A_true=A_true, phi_true=phi_true)
  result_7500 <- model_measurements(x[1:7500,], dist_class=dist_class, A_true=A_true, phi_true=phi_true)
  result_10000 <- model_measurements(x[1:10000,], dist_class=dist_class, A_true=A_true, phi_true=phi_true)
  return(list(n_5000=result_5000, n_7500=result_7500,n_10000=result_10000))
}, mc.cores = core_num)





#### Model 2 ####
# Data generation
m <- 9

data_dim <- 1

# Define true parameters
A1_real <- matrix(c(1,0,0,
                    0,1,0,
                    1,0,0), ncol = 3, byrow = TRUE)+0.2
A1_real <- A1_real / rowSums(A1_real)
A2_real <- matrix(c(0,0,1,
                    0,1,0,
                    0,0,1), ncol = 3, byrow = TRUE)+0.2
A2_real <- A2_real / rowSums(A2_real)

A_true <- kronecker(A1_real, A2_real)

dist_class <- "norm"
mean <- seq(from = -400, by = 100, length.out = 9)
sd <- rep(25, times=9)
phi_true <- mapply(function(xx,yy){list(mean=xx, sigma=matrix(yy^2))}, xx=mean, yy=sd, SIMPLIFY = FALSE)

data <- mclapply(Model2_seed,function(x){
  set.seed(x, kind = "L'Ecuyer-CMRG")
  HMM_simulation(dimension = data_dim,
                 transition = A_true, 
                 emission = phi_true, 
                 size = 10000)
},mc.cores = core_num)

dist_class <- "norm"
mean <- seq(from = -400, by = 100, length.out = 9)
sd <- rep(25, times=9)
phi_true <- mapply(function(xx,yy){list(mean=xx, sigma=matrix(yy))}, xx=mean, yy=sd, SIMPLIFY = FALSE)


# Define a function generating numerical results across OE, MLE, MPLE
model_measurements <- function(data, dist_class, A_true, phi_true){
  OE_result <- EM_Oracle_Sim2_A2(data,dist_class = dist_class, m=9)
  MLE_result <- EM_MLE(data,dist_class = dist_class, m=9)
  MPLE_GLASSO_result <- EM_MPLE_GLASSO(data,dist_class = dist_class, m=9)
  MPLE_result <- EM_MPLE(data,dist_class = dist_class, m=9)
  FHMM_result <- fhmm(data, M=2, K=3, cyc=100)
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
  MPLE_GLASSO_val <- data.frame(TE=n2_sqrt(A_true-MPLE_GLASSO_result$transition),
                                EE=n2_sqrt(unlist(phi_true)-unlist(MPLE_GLASSO_result$par)),
                                TP=TFPN(cf(A_true)$index,
                                        cf(MPLE_GLASSO_result$transition)$index)$TP,
                                FP=TFPN(cf(A_true)$index,
                                        cf(MPLE_GLASSO_result$transition)$index)$FP)
  MPLE_val <- data.frame(TE=n2_sqrt(A_true-MPLE_result$transition),
                         EE=n2_sqrt(unlist(phi_true)-unlist(MPLE_result$par)),
                         TP=TFPN(cf(A_true)$index,
                                 cf(MPLE_result$transition)$index)$TP,
                         FP=TFPN(cf(A_true)$index,
                                 cf(MPLE_result$transition)$index)$FP)
  
  fhmm_par <- FHMM_result$par
  
  fhmm_A <- kronecker(FHMM_result$P[1:3,], FHMM_result$P[4:6,])
  
  FHMM_val <- data.frame(TE=n2_sqrt(A_true-fhmm_A),
                         EE=n2_sqrt(unlist(phi_true)-unlist(fhmm_par)),
                         TP=TFPN(cf(A_true)$index,
                                 cf(fhmm_A)$index)$TP,
                         FP=TFPN(cf(A_true)$index,
                                 cf(fhmm_A)$index)$FP)
  
  return(list(OE=OE_val, MLE=MLE_val, MPLE_GLASSO=MPLE_GLASSO_val, MPLE=MPLE_val, FHMM=FHMM_val))
}


# Calculate model measurements
Model2_result <- mclapply(data, function(x){
  result_5000 <- model_measurements(x[1:5000], dist_class=dist_class, A_true=A_true, phi_true=phi_true)
  result_7500 <- model_measurements(x[1:7500], dist_class=dist_class, A_true=A_true, phi_true=phi_true)
  result_10000 <- model_measurements(x[1:10000], dist_class=dist_class, A_true=A_true, phi_true=phi_true)
  return(list(n_5000=result_5000, n_7500=result_7500,n_10000=result_10000))
}, mc.cores = core_num)










#### Model 3 ####
# Data generation
m <- 9

data_dim <- 3

# Define true parameters
A1_real <- diag(3)+0.2
A1_real <- A1_real / rowSums(A1_real)
A2_real <- matrix(1, ncol=3, nrow=3)+0.2
A2_real <- A2_real / rowSums(A2_real)

A_true <- kronecker(A1_real, A2_real)

dist_class <- "mvnorm"

phi_true=list()
mean_mat <- matrix(c(0, 5, 10, 15, 20, 25, 30, 35, 40, 
                     0, 10, 20, 30, 40, 50, 60, 70, 80, 
                     0, 15, 30, 45, 60, 75, 90, 105, 120), ncol = 3 )
for (i in c(1:9)){
  phi_true[[i]] <- list(mean = as.vector(mean_mat[i,]), sigma=matrix(c(1, 0.1, 0.1,
                                                                       0.1, 2^2, 0.1,
                                                                       0.1, 0.1, 3^2), ncol = 3))
}


# Data generation
data <- mclapply(Model3_seed,function(x){
  set.seed(x, kind = "L'Ecuyer-CMRG")
  HMM_simulation(dimension = data_dim,
                 transition = A_true, 
                 emission = phi_true, 
                 size = 10000)
},mc.cores = core_num)


# Define a function generating numerical results across OE, MLE, MPLE
model_measurements <- function(data, dist_class, A_true, phi_true){
  OE_result <- EM_Oracle_Sim2_A1(data,dist_class = dist_class, m=9)
  MLE_result <- EM_MLE(data,dist_class = dist_class, m=9)
  MPLE_GLASSO_result <- EM_MPLE_GLASSO(data,dist_class = dist_class, m=9)
  MPLE_result <- EM_MPLE(data,dist_class = dist_class, m=9)
  FHMM_result <- fhmm(data, M=2, K=3, cyc=100)
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
  MPLE_GLASSO_val <- data.frame(TE=n2_sqrt(A_true-MPLE_GLASSO_result$transition),
                                EE=n2_sqrt(unlist(phi_true)-unlist(MPLE_GLASSO_result$par)),
                                TP=TFPN(cf(A_true)$index,
                                        cf(MPLE_GLASSO_result$transition)$index)$TP,
                                FP=TFPN(cf(A_true)$index,
                                        cf(MPLE_GLASSO_result$transition)$index)$FP)
  MPLE_val <- data.frame(TE=n2_sqrt(A_true-MPLE_result$transition),
                         EE=n2_sqrt(unlist(phi_true)-unlist(MPLE_result$par)),
                         TP=TFPN(cf(A_true)$index,
                                 cf(MPLE_result$transition)$index)$TP,
                         FP=TFPN(cf(A_true)$index,
                                 cf(MPLE_result$transition)$index)$FP)
  
  fhmm_par <- FHMM_result$par
  
  fhmm_A <- kronecker(FHMM_result$P[1:3,], FHMM_result$P[4:6,])
  
  FHMM_val <- data.frame(TE=n2_sqrt(A_true-fhmm_A),
                         EE=n2_sqrt(unlist(phi_true)-unlist(fhmm_par)),
                         TP=TFPN(cf(A_true)$index,
                                 cf(fhmm_A)$index)$TP,
                         FP=TFPN(cf(A_true)$index,
                                 cf(fhmm_A)$index)$FP)
  
  return(list(OE=OE_val, MLE=MLE_val, MPLE_GLASSO=MPLE_GLASSO_val, MPLE=MPLE_val, FHMM=FHMM_val))
}

# Calculate model measurements
Model3_result <- mclapply(data, function(x){
  result_5000 <- model_measurements(x[1:5000,], dist_class=dist_class, A_true=A_true, phi_true=phi_true)
  result_7500 <- model_measurements(x[1:7500,], dist_class=dist_class, A_true=A_true, phi_true=phi_true)
  result_10000 <- model_measurements(x[1:10000,], dist_class=dist_class, A_true=A_true, phi_true=phi_true)
  return(list(n_5000=result_5000, n_7500=result_7500,n_10000=result_10000))
}, mc.cores = core_num)





#### Model 4 ####
# Data generation
m <- 9

data_dim <- 3

# Define true parameters
A1_real <- matrix(c(1,0,0,
                    0,1,0,
                    1,0,0), ncol = 3, byrow = TRUE)+0.2
A1_real <- A1_real / rowSums(A1_real)
A2_real <- matrix(c(0,0,1,
                    0,1,0,
                    0,0,1), ncol = 3, byrow = TRUE)+0.2
A2_real <- A2_real / rowSums(A2_real)

A_true <- kronecker(A1_real, A2_real)

dist_class <- "mvnorm"

phi_true=list()
mean_mat <- matrix(c(0, 5, 10, 15, 20, 25, 30, 35, 40, 
                     0, 10, 20, 30, 40, 50, 60, 70, 80, 
                     0, 15, 30, 45, 60, 75, 90, 105, 120), ncol = 3 )
for (i in c(1:9)){
  phi_true[[i]] <- list(mean = as.vector(mean_mat[i,]), sigma=matrix(c(1, 0.1, 0.1,
                                                                       0.1, 2^2, 0.1,
                                                                       0.1, 0.1, 3^2), ncol = 3))
}

# Data generation
data <- mclapply(Model4_seed,function(x){
  set.seed(x, kind = "L'Ecuyer-CMRG")
  HMM_simulation(dimension = data_dim,
                 transition = A_true, 
                 emission = phi_true, 
                 size = 10000)
},mc.cores = core_num)


# Define a function generating numerical results across OE, MLE, MPLE
model_measurements <- function(data, dist_class, A_true, phi_true){
  OE_result <- EM_Oracle_Sim2_A2(data,dist_class = dist_class, m=9)
  MLE_result <- EM_MLE(data,dist_class = dist_class, m=9)
  MPLE_GLASSO_result <- EM_MPLE_GLASSO(data,dist_class = dist_class, m=9)
  MPLE_result <- EM_MPLE(data,dist_class = dist_class, m=9)
  FHMM_result <- fhmm(data, M=2, K=3, cyc=100)
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
  MPLE_GLASSO_val <- data.frame(TE=n2_sqrt(A_true-MPLE_GLASSO_result$transition),
                                EE=n2_sqrt(unlist(phi_true)-unlist(MPLE_GLASSO_result$par)),
                                TP=TFPN(cf(A_true)$index,
                                        cf(MPLE_GLASSO_result$transition)$index)$TP,
                                FP=TFPN(cf(A_true)$index,
                                        cf(MPLE_GLASSO_result$transition)$index)$FP)
  MPLE_val <- data.frame(TE=n2_sqrt(A_true-MPLE_result$transition),
                         EE=n2_sqrt(unlist(phi_true)-unlist(MPLE_result$par)),
                         TP=TFPN(cf(A_true)$index,
                                 cf(MPLE_result$transition)$index)$TP,
                         FP=TFPN(cf(A_true)$index,
                                 cf(MPLE_result$transition)$index)$FP)
  
  fhmm_par <- FHMM_result$par
  
  fhmm_A <- kronecker(FHMM_result$P[1:3,], FHMM_result$P[4:6,])
  
  FHMM_val <- data.frame(TE=n2_sqrt(A_true-fhmm_A),
                         EE=n2_sqrt(unlist(phi_true)-unlist(fhmm_par)),
                         TP=TFPN(cf(A_true)$index,
                                 cf(fhmm_A)$index)$TP,
                         FP=TFPN(cf(A_true)$index,
                                 cf(fhmm_A)$index)$FP)
  
  return(list(OE=OE_val, MLE=MLE_val, MPLE_GLASSO=MPLE_GLASSO_val, MPLE=MPLE_val, FHMM=FHMM_val))
}


# Calculate model measurements
Model4_result <- mclapply(data, function(x){
  result_5000 <- model_measurements(x[1:5000,], dist_class=dist_class, A_true=A_true, phi_true=phi_true)
  result_7500 <- model_measurements(x[1:7500,], dist_class=dist_class, A_true=A_true, phi_true=phi_true)
  result_10000 <- model_measurements(x[1:10000,], dist_class=dist_class, A_true=A_true, phi_true=phi_true)
  return(list(n_5000=result_5000, n_7500=result_7500,n_10000=result_10000))
}, mc.cores = core_num)



save(Model1_result, Model2_result, Model3_result, Model4_result, file = "results/Simulation2.Rdata")








