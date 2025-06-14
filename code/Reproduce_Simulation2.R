##################################################################
#### Code for reproduce the numerical results of Simulation 2 ####
##################################################################
#### Preliminary ####
# Load functions
sapply(c("code/EM_MLE.R","code/EM_MPLE.R","code/EM_FHMM.R","code/EM_Oracle.R"), source)

# Required library to load data 
library(parallel) # Use multicore.
core_num <- 25 # Set number of cores.

#### Model 1 ####
load("data/Sim2_Model1.Rdata")

# Define a function generating numerical results across OE, MLE, MPLE
model_measurements <- function(data, dist_class, A_true, phi_true){
  OE_result <- EM_Oracle_Sim2_A1(data,dist_class = dist_class, m=9)
  MLE_result <- EM_MLE(data,dist_class = dist_class, m=9)
  MPLE_result <- EM_MPLE(data,dist_class = dist_class, m=9)
  FHMM_result <- fhmm(data, M=2, K=3)
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
  # Transfrom FHMM to HMM 
  hmm2fhmm_par <- function(x, M, K){
    index_mat <- sapply(1:M, function(x){rep(1:K, K^(x-1), each=K^(M-x))})
    fhmm_mean <- apply(index_mat, 1, function(xx){
      x_new <- 0
      for (i in 1:M) {
        x_new <- x_new + x[K*(i-1)+xx[i]]
      }
      return(x_new)
    })
    return(fhmm_mean)
  }
  Mu_hmm <- apply(X=FHMM_result$Mu, MARGIN=2, hmm2fhmm_par, M=2, K=3)
  if(dist_class=="norm"){
    fhmm_par <- apply(Mu_hmm, 1,function(x){list(mean=x, sd = sqrt(FHMM_result$Cov))}, simplify = TRUE)
  }else{
    fhmm_par <- apply(Mu_hmm, 1,function(x){list(mean=x, sigma = FHMM_result$Cov)}, simplify = TRUE)
  }
  fhmm_A <- kronecker(FHMM_result$P[1:3,], FHMM_result$P[4:6,])
  
  FHMM_val <- data.frame(TE=n2_sqrt(A_true-fhmm_A),
                         EE=n2_sqrt(unlist(phi_true)-unlist(fhmm_par)),
                         TP=TFPN(cf(A_true)$index,
                                 cf(fhmm_A)$index)$TP,
                         FP=TFPN(cf(A_true)$index,
                                 cf(fhmm_A)$index)$FP)
  
  return(list(OE=OE_val, MLE=MLE_val, MPLE=MPLE_val, FHMM=FHMM_val))
}

# Define true parameters
# A1_real과 A2_real 정의
A1_real <- diag(3)+0.2
A1_real <- A1_real / rowSums(A1_real)
A2_real <- matrix(1, ncol=3, nrow=3)+0.2
A2_real <- A2_real / rowSums(A2_real)

A_true <- kronecker(A1_real, A2_real)

dist_class <- "norm"
mean <- seq(from = -400, by = 100, length.out = 9)
sd <- rep(25, times=9)
phi_true <- mapply(function(xx,yy){list(mean=xx, sigma=matrix(yy))}, xx=mean, yy=sd, SIMPLIFY = FALSE)

# Calculate model measurements
Model1_result <- mclapply(Sim2_Model1, function(x){
  result_4000 <- model_measurements(x[1:4000], dist_class=dist_class, A_true=A_true, phi_true=phi_true)
  result_5000 <- model_measurements(x[1:5000], dist_class=dist_class, A_true=A_true, phi_true=phi_true)
  result_6000 <- model_measurements(x[1:6000], dist_class=dist_class, A_true=A_true, phi_true=phi_true)
  return(list(n_4000=result_4000, n_5000=result_5000,n_6000=result_6000))
}, mc.cores = core_num)





#### Model 2 ####
load("data/Sim2_Model2.Rdata")

# Define a function generating numerical results across OE, MLE, MPLE
model_measurements <- function(data, dist_class, A_true, phi_true){
  OE_result <- EM_Oracle_Sim2_A2(data,dist_class = dist_class, m=9)
  MLE_result <- EM_MLE(data,dist_class = dist_class, m=9)
  MPLE_result <- EM_MPLE(data,dist_class = dist_class, m=9)
  FHMM_result <- fhmm(data, M=2, K=3)
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
  # Transfrom FHMM to HMM 
  hmm2fhmm_par <- function(x, M, K){
    index_mat <- sapply(1:M, function(x){rep(1:K, K^(x-1), each=K^(M-x))})
    fhmm_mean <- apply(index_mat, 1, function(xx){
      x_new <- 0
      for (i in 1:M) {
        x_new <- x_new + x[K*(i-1)+xx[i]]
      }
      return(x_new)
    })
    return(fhmm_mean)
  }
  Mu_hmm <- apply(X=FHMM_result$Mu, MARGIN=2, hmm2fhmm_par, M=2, K=3)
  if(dist_class=="norm"){
    fhmm_par <- apply(Mu_hmm, 1,function(x){list(mean=x, sd = sqrt(FHMM_result$Cov))}, simplify = TRUE)
  }else{
    fhmm_par <- apply(Mu_hmm, 1,function(x){list(mean=x, sigma = FHMM_result$Cov)}, simplify = TRUE)
  }
  fhmm_A <- kronecker(FHMM_result$P[1:3,], FHMM_result$P[4:6,])
  
  FHMM_val <- data.frame(TE=n2_sqrt(A_true-fhmm_A),
                         EE=n2_sqrt(unlist(phi_true)-unlist(fhmm_par)),
                         TP=TFPN(cf(A_true)$index,
                                 cf(fhmm_A)$index)$TP,
                         FP=TFPN(cf(A_true)$index,
                                 cf(fhmm_A)$index)$FP)
  
  return(list(OE=OE_val, MLE=MLE_val, MPLE=MPLE_val, FHMM=FHMM_val))
}

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
phi_true <- mapply(function(xx,yy){list(mean=xx, sigma=matrix(yy))}, xx=mean, yy=sd, SIMPLIFY = FALSE)


# Calculate model measurements
Model2_result <- mclapply(Sim2_Model2, function(x){
  result_4000 <- model_measurements(x[1:4000], dist_class=dist_class, A_true=A_true, phi_true=phi_true)
  result_5000 <- model_measurements(x[1:5000], dist_class=dist_class, A_true=A_true, phi_true=phi_true)
  result_6000 <- model_measurements(x[1:6000], dist_class=dist_class, A_true=A_true, phi_true=phi_true)
  return(list(n_4000=result_4000, n_5000=result_5000,n_6000=result_6000))
}, mc.cores = core_num)










#### Model 3 ####
load("data/Sim2_Model3.Rdata")

# Define a function generating numerical results across OE, MLE, MPLE
model_measurements <- function(data, dist_class, A_true, phi_true){
  OE_result <- EM_Oracle_Sim2_A1(data,dist_class = dist_class, m=9)
  MLE_result <- EM_MLE(data,dist_class = dist_class, m=9)
  MPLE_result <- EM_MPLE(data,dist_class = dist_class, m=9)
  FHMM_result <- fhmm(data, M=2, K=3)
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
  # Transfrom FHMM to HMM 
  hmm2fhmm_par <- function(x, M, K){
    index_mat <- sapply(1:M, function(x){rep(1:K, K^(x-1), each=K^(M-x))})
    fhmm_mean <- apply(index_mat, 1, function(xx){
      x_new <- 0
      for (i in 1:M) {
        x_new <- x_new + x[K*(i-1)+xx[i]]
      }
      return(x_new)
    })
    return(fhmm_mean)
  }
  Mu_hmm <- apply(X=FHMM_result$Mu, MARGIN=2, hmm2fhmm_par, M=2, K=3)
  if(dist_class=="norm"){
    fhmm_par <- apply(Mu_hmm, 1,function(x){list(mean=x, sd = sqrt(FHMM_result$Cov))}, simplify = TRUE)
  }else{
    fhmm_par <- apply(Mu_hmm, 1,function(x){list(mean=x, sigma = FHMM_result$Cov)}, simplify = TRUE)
  }
  fhmm_A <- kronecker(FHMM_result$P[1:3,], FHMM_result$P[4:6,])
  
  FHMM_val <- data.frame(TE=n2_sqrt(A_true-fhmm_A),
                         EE=n2_sqrt(unlist(phi_true)-unlist(fhmm_par)),
                         TP=TFPN(cf(A_true)$index,
                                 cf(fhmm_A)$index)$TP,
                         FP=TFPN(cf(A_true)$index,
                                 cf(fhmm_A)$index)$FP)
  
  return(list(OE=OE_val, MLE=MLE_val, MPLE=MPLE_val, FHMM=FHMM_val))
}

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
# Calculate model measurements
Model3_result <- mclapply(Sim2_Model3, function(x){
  result_4000 <- model_measurements(x[1:4000,], dist_class=dist_class, A_true=A_true, phi_true=phi_true)
  result_5000 <- model_measurements(x[1:5000,], dist_class=dist_class, A_true=A_true, phi_true=phi_true)
  result_6000 <- model_measurements(x[1:6000,], dist_class=dist_class, A_true=A_true, phi_true=phi_true)
  return(list(n_4000=result_4000, n_5000=result_5000,n_6000=result_6000))
}, mc.cores = core_num)





#### Model 4 ####
load("data/Sim2_Model4.Rdata")

# Define a function generating numerical results across OE, MLE, MPLE
model_measurements <- function(data, dist_class, A_true, phi_true){
  OE_result <- EM_Oracle_Sim2_A2(data,dist_class = dist_class, m=9)
  MLE_result <- EM_MLE(data,dist_class = dist_class, m=9)
  MPLE_result <- EM_MPLE(data,dist_class = dist_class, m=9)
  FHMM_result <- fhmm(data, M=2, K=3)
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
  # Transfrom FHMM to HMM 
  hmm2fhmm_par <- function(x, M, K){
    index_mat <- sapply(1:M, function(x){rep(1:K, K^(x-1), each=K^(M-x))})
    fhmm_mean <- apply(index_mat, 1, function(xx){
      x_new <- 0
      for (i in 1:M) {
        x_new <- x_new + x[K*(i-1)+xx[i]]
      }
      return(x_new)
    })
    return(fhmm_mean)
  }
  Mu_hmm <- apply(X=FHMM_result$Mu, MARGIN=2, hmm2fhmm_par, M=2, K=3)
  if(dist_class=="norm"){
    fhmm_par <- apply(Mu_hmm, 1,function(x){list(mean=x, sd = sqrt(FHMM_result$Cov))}, simplify = TRUE)
  }else{
    fhmm_par <- apply(Mu_hmm, 1,function(x){list(mean=x, sigma = FHMM_result$Cov)}, simplify = TRUE)
  }
  fhmm_A <- kronecker(FHMM_result$P[1:3,], FHMM_result$P[4:6,])
  
  FHMM_val <- data.frame(TE=n2_sqrt(A_true-fhmm_A),
                         EE=n2_sqrt(unlist(phi_true)-unlist(fhmm_par)),
                         TP=TFPN(cf(A_true)$index,
                                 cf(fhmm_A)$index)$TP,
                         FP=TFPN(cf(A_true)$index,
                                 cf(fhmm_A)$index)$FP)
  
  return(list(OE=OE_val, MLE=MLE_val, MPLE=MPLE_val, FHMM=FHMM_val))
}

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
# Calculate model measurements
Model4_result <- mclapply(Sim2_Model4[1:2], function(x){
  result_4000 <- model_measurements(x[1:4000,], dist_class=dist_class, A_true=A_true, phi_true=phi_true)
  result_5000 <- model_measurements(x[1:5000,], dist_class=dist_class, A_true=A_true, phi_true=phi_true)
  result_6000 <- model_measurements(x[1:6000,], dist_class=dist_class, A_true=A_true, phi_true=phi_true)
  return(list(n_4000=result_4000, n_5000=result_5000,n_6000=result_6000))
}, mc.cores = core_num)



save(Model1_result, Model2_result, Model3_result, Model4_result, file = "results/Simulation2.Rdata")








