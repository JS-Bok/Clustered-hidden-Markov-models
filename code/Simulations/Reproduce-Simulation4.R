##################################################################
#### Code for reproduce the numerical results of Simulation 4 ####
##################################################################
#### Preliminary ####
# Load functions
sapply(c("code/Functions/EM_MLE.R",
         "code/Functions/EM_MPLE.R",
         "code/Functions/EM_Oracle.R",
         "code/Functions/CMC.R"), source)


dir.create("results/Simulation", recursive = TRUE, showWarnings = FALSE)

# Required library to load data 
library(parallel) # Use multicore.
core_num <- 50 # Set number of cores.

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

model_measurements <- function(data, dist_class, A_true, phi_true,
                               ord, prefix){
  t <- nrow(data)
  m <- nrow(A_true)
  
  OE_result <- CMC_Oracle(data, A_true = A_true, m = m)
  #saveRDS(OE_result, file = paste0("results/Simulation/", prefix, "-", t, "-OE-", ord, ".rds"))
  
  MLE_result <- CMC_MLE(data, m = m)
  #saveRDS(MLE_result, file = paste0("results/Simulation/", prefix, "-", t, "-MLE-", ord, ".rds"))
  
  MPLE_result_all <- CMC_MPLE_SaveAllLambda(data)
  bic_vec   <- sapply(MPLE_result_all, `[[`, "BIC")
  MPLE_result <- MPLE_result_all[[ which.min(bic_vec) ]]
  #saveRDS(MPLE_result, file = paste0("results/Simulation/", prefix, "-", t, "-MPLE-", ord, ".rds"))
  
  MPLE_FSCAD_result_all <- CMC_MPLE_SaveAllLambda(data, penalty_type = "fused_scad")
  bic_vec   <- sapply(MPLE_FSCAD_result_all, `[[`, "BIC")
  MPLE_FSCAD_result <- MPLE_FSCAD_result_all[[ which.min(bic_vec) ]]
  #saveRDS(MPLE_FSCAD_result, file = paste0("results/Simulation/", prefix, "-", t, "-FSCAD-", ord, ".rds"))
  
  MPLE_GLASSO_result_all <- CMC_MPLE_SaveAllLambda(data, penalty_type = "group_fused_lasso")
  bic_vec   <- sapply(MPLE_GLASSO_result_all, `[[`, "BIC")
  MPLE_GLASSO_result <- MPLE_GLASSO_result_all[[ which.min(bic_vec) ]]
  #saveRDS(MPLE_GLASSO_result, file = paste0("results/Simulation/", prefix, "-", t, "-GLASSO-", ord, ".rds"))
  
  MPLE_GSCAD_result_all <- CMC_MPLE_SaveAllLambda(data, penalty_type = "group_scad")
  bic_vec   <- sapply(MPLE_GSCAD_result_all, `[[`, "BIC")
  MPLE_GSCAD_result <- MPLE_GSCAD_result_all[[ which.min(bic_vec) ]]
  #saveRDS(MPLE_GSCAD_result, file = paste0("results/Simulation/", prefix, "-", t, "-GSCAD-", ord, ".rds"))
  
  ## ------------------------------------------------------------
  ## measurement extractor including predictive likelihood error
  ## ------------------------------------------------------------
  extract_meas <- function(x){
    
    tfpn_x <- TFPN(cf(A_true)$index, cf(x$transition)$index)
    out <- data.frame(
      TPC = tfpn_x$TP,
      FPC = tfpn_x$FP
    )
    return(out)
  }
  
  OE_val           <- extract_meas(OE_result)
  MLE_val          <- extract_meas(MLE_result)
  MPLE_GLASSO_val  <- extract_meas(MPLE_GLASSO_result)
  MPLE_GSCAD_val   <- extract_meas(MPLE_GSCAD_result)
  MPLE_FSCAD_val   <- extract_meas(MPLE_FSCAD_result)
  MPLE_val         <- extract_meas(MPLE_result)
  
  return(list(
    OE          = OE_val,
    MLE         = MLE_val,
    MPLE_GLASSO = MPLE_GLASSO_val,
    MPLE_GSCAD  = MPLE_GSCAD_val,
    MPLE_FSCAD  = MPLE_FSCAD_val,
    MPLE        = MPLE_val
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
  prefix <- paste0(model_name, "-M", m, "-K", K, "-zeta", zeta)
  
  # Main loop (parallel): for each replicate, simulate once then fit for n=5000/10000
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
  
  save(result_list, file = paste0("results/Simulation/", prefix, ".Rdata"))
  
  return(list(A_true = A_true, phi_true = phi_true, results = result_list))
}



#### 5) Run simulations ####

# Well-separated setting (zeta=0.8)
Sim4_WS_M1 <- run_one_WS_model(model_name = "Sim4-M1", seed_vec = Model1_seed, m = 10, K = 3, zeta = 0.8,
                          data_dim = data_dim, dist_class = dist_class)

Sim4_WS_M2 <- run_one_WS_model(model_name = "Sim4-M2", seed_vec = Model2_seed, m = 10, K = 5, zeta = 0.8,
                          data_dim = data_dim, dist_class = dist_class)

Sim4_WS_M3 <- run_one_WS_model(model_name = "Sim4-M3", seed_vec = Model3_seed, m = 20, K = 3, zeta = 0.8,
                          data_dim = data_dim, dist_class = dist_class)

Sim4_WS_M4 <- run_one_WS_model(model_name = "Sim4-M4", seed_vec = Model4_seed, m = 20, K = 5, zeta = 0.8,
                          data_dim = data_dim, dist_class = dist_class)

# Nearly-tied setting (zeta=0.2)
Sim4_NT_M1 <- run_one_WS_model(model_name = "Sim4-M1", seed_vec = Model1_seed, m = 10, K = 3, zeta = 0.2,
                               data_dim = data_dim, dist_class = dist_class)

Sim4_NT_M2 <- run_one_WS_model(model_name = "Sim4-M2", seed_vec = Model2_seed, m = 10, K = 5, zeta = 0.2,
                               data_dim = data_dim, dist_class = dist_class)

Sim4_NT_M3 <- run_one_WS_model(model_name = "Sim4-M3", seed_vec = Model3_seed, m = 20, K = 3, zeta = 0.2,
                               data_dim = data_dim, dist_class = dist_class)

Sim4_NT_M4 <- run_one_WS_model(model_name = "Sim4-M4", seed_vec = Model4_seed, m = 20, K = 5, zeta = 0.2,
                               data_dim = data_dim, dist_class = dist_class)

save(Sim4_WS_M1, Sim4_WS_M2, Sim4_WS_M3, Sim4_WS_M4, 
     Sim4_NT_M1, Sim4_NT_M2, Sim4_NT_M3, Sim4_NT_M4, 
     file = "results/Simulation/Simulation4_results.Rdata")