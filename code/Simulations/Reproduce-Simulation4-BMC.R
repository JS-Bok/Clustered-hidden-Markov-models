##################################################################
#### Code for reproduce the numerical results of Simulation 4 ####
##################################################################
#### Preliminary ####
# Load functions
sapply(c("code/Functions/BMC.R",
         "code/Functions/CMC.R"), source)


dir.create("results/Simulation", recursive = TRUE, showWarnings = FALSE)

# Required library to load data 
library(parallel) # Use multicore.
core_num <- 1 # Package 'Parallel' only support single core for Windows.

# Seeds for data generation
set.seed(2026, kind = "L'Ecuyer-CMRG")
Model1_seed <- sample(1:1e8,500)
Model2_seed <- sample(1:1e8,500)
Model3_seed <- sample(1:1e8,500)
Model4_seed <- sample(1:1e8,500)

data_dim   <- 3
dist_class <- "mvnorm"

#### 1) Well-separated transition matrix generator (A)   ####

make_A_WS <- function(M, K, zeta){
  # (1) Partition rule for well-separated rows
  # - K=3: M=10 -> (3,4,3), M=20 -> (6,8,6)
  make_groups_WS <- function(M, K){
    if (K == 3){
      if (M == 10) sizes <- c(3, 4, 3)
      else if (M == 20) sizes <- c(6, 8, 6)
    } else if (K == 5){
      if (M == 10) sizes <- c(2, 2, 2, 2, 2)
      else if (M == 20) sizes <- c(4, 4, 4, 4, 4)
    } else if (K == M){
      sizes <- rep(1,times=M)
    }
    
    ends   <- cumsum(sizes)
    starts <- c(1, head(ends, -1) + 1)
    groups <- Map(seq, starts, ends)  # list of index vectors
    return(groups)
  }
  
  # (2) Build block-uniform B
  make_block_uniform_B <- function(M, groups){
    B <- matrix(0, M, M)
    for (g in groups){
      mk <- length(g)
      B[g, g] <- matrix(1 / mk, nrow = mk, ncol = mk)
    }
    return(B)
  }
  
  # (3) Mix: A = zeta B + (1-zeta) U
  make_A_mix <- function(B, zeta){
    M <- nrow(B)
    U <- matrix(1 / M, M, M)
    A <- zeta * B + (1 - zeta) * U
    A <- A / rowSums(A)
    return(A)
  }
  
  groups <- make_groups_WS(M, K)
  B      <- make_block_uniform_B(M, groups)
  A      <- make_A_mix(B, zeta)
  return(list(A = A, groups = groups))
}


#### 2) Emission parameter generator (phi_true)           ####
make_phi_true <- function(m, data_dim = 3){
  phi_true <- vector("list", m)
  
  # same mean pattern you used (just generalized to m)
  mean_mat <- matrix(c(seq(from = 4,  by = 4,  length.out = m),
                       seq(from = 8,  by = 8,  length.out = m),
                       seq(from = 12, by = 12, length.out = m)), ncol = data_dim)
  
  Sigma <- matrix(c(1,   0.1, 0.1,
                    0.1, 2^2, 0.1,
                    0.1, 0.1, 3^2), ncol = data_dim)
  
  for (i in 1:m){
    phi_true[[i]] <- list(mean = as.vector(mean_mat[i, ]), sigma = Sigma)
  }
  return(phi_true)
}


#### 3) Model fitting + measurement ####
model_measurements <- function(data, dist_class, A_true, phi_true, ord, prefix){
  t <- length(data)
  m <- nrow(A_true)
  
  emb <- embed(data,2)
  
  pair_mat <- table(emb[,2],emb[,1])
  
  BMC_result <- run_bmc_autoK(pair_mat)
  #saveRDS(BMC_result, file = paste0("results/Simulation/", prefix, "-", t, "-BMC-", ord, ".rds"))
  
  BMC_index <- split(seq_along(BMC_result$clustering_final), factor(BMC_result$clustering_final, levels = unique(BMC_result$clustering_final)))
  
  BMC_val        <- data.frame(
    TPC = TFPN(cf(A_true)$index, BMC_index)$TP,
    FPC = TFPN(cf(A_true)$index, BMC_index)$FP
  )
  
  return(list(
    BMC = BMC_val
  ))
}


#### 4) One WS model runner (500 reps; n=5000 & 10000)     ####
run_one_WS_model <- function(model_name, seed_vec, m, K, zeta,
                             data_dim = 3, dist_class = "mvnorm",
                             n_list = c(5000, 10000), size_full = 10000){
  
  # Build A_true
  A_obj  <- make_A_WS(M = m, K = K, zeta = zeta)
  A_true <- A_obj$A
  
  # Build emission
  phi_true <- make_phi_true(m = m, data_dim = data_dim)
  
  # Prefix for saved files
  prefix <- paste0("S1-WS-", model_name, "-M", m, "-K", K, "-zeta", zeta)
  
  # Main loop: for each replicate, simulate once then fit for n=5000/10000
  result_list <- mclapply(seq_along(seed_vec), function(ord){
    
    set.seed(seed_vec[ord], kind = "L'Ecuyer-CMRG")
    
    dat_full <- HMM_simulation(
      dimension  = data_dim,
      transition = A_true,
      emission   = phi_true,
      size       = size_full,
      return_X = T
    )
    
    # Fit for each n in n_list
    out <- list()
    for (n in n_list){
      out[[paste0("n_", n)]] <- model_measurements(
        data       = dat_full$X[1:n],
        dist_class = dist_class,
        A_true     = A_true,
        phi_true   = phi_true,
        ord        = ord,
        prefix     = prefix
      )
      
      # saveRDS(out[[paste0("n_", n)]],
      #         file = paste0("results/Simulation/", prefix, "-", n, "-", ord, ".rds"))
    }
    
    return(out)
    
  }, mc.cores = core_num, mc.preschedule = FALSE)
  
  save(result_list, file = paste0("results/Simulation/", prefix, "-BMC.Rdata"))
  
  return(list(A_true = A_true, phi_true = phi_true, results = result_list))
}



#### 5) Run simulations ####

# Well-separated setting (zeta=0.8)
Sim4_WS_M1_BMC <- run_one_WS_model(model_name = "Sim4-M1", seed_vec = Model1_seed, m = 10, K = 3, zeta = 0.8,
                               data_dim = data_dim, dist_class = dist_class)

Sim4_WS_M2_BMC <- run_one_WS_model(model_name = "Sim4-M2", seed_vec = Model2_seed, m = 10, K = 5, zeta = 0.8,
                               data_dim = data_dim, dist_class = dist_class)

Sim4_WS_M3_BMC <- run_one_WS_model(model_name = "Sim4-M3", seed_vec = Model3_seed, m = 20, K = 3, zeta = 0.8,
                               data_dim = data_dim, dist_class = dist_class)

Sim4_WS_M4_BMC <- run_one_WS_model(model_name = "Sim4-M4", seed_vec = Model4_seed, m = 20, K = 5, zeta = 0.8,
                               data_dim = data_dim, dist_class = dist_class)

# Nearly-tied setting (zeta=0.2)
Sim4_NT_M1_BMC <- run_one_WS_model(model_name = "Sim4-M1", seed_vec = Model1_seed, m = 10, K = 3, zeta = 0.2,
                               data_dim = data_dim, dist_class = dist_class)

Sim4_NT_M2_BMC <- run_one_WS_model(model_name = "Sim4-M2", seed_vec = Model2_seed, m = 10, K = 5, zeta = 0.2,
                               data_dim = data_dim, dist_class = dist_class)

Sim4_NT_M3_BMC <- run_one_WS_model(model_name = "Sim4-M3", seed_vec = Model3_seed, m = 20, K = 3, zeta = 0.2,
                               data_dim = data_dim, dist_class = dist_class)

Sim4_NT_M4_BMC <- run_one_WS_model(model_name = "Sim4-M4", seed_vec = Model4_seed, m = 20, K = 5, zeta = 0.2,
                               data_dim = data_dim, dist_class = dist_class)

save(Sim4_WS_M1_BMC, Sim4_WS_M2_BMC, Sim4_WS_M3_BMC, Sim4_WS_M4_BMC, 
     Sim4_NT_M1_BMC, Sim4_NT_M2_BMC, Sim4_NT_M3_BMC, Sim4_NT_M4_BMC, 
     file = "results/Simulation/Simulation4_results-BMC.Rdata")