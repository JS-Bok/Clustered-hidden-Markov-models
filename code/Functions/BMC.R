library(reticulate)
# # 1) Generate Python environment (You must generate it first!!!)
# conda_create("bmc-py", packages = "python=3.11")
# 
# use_condaenv("bmc-py", required = TRUE)
# 
# py_install(c("numpy<2", "BMCToolkit"), envname = "bmc-py", pip = TRUE)



np  <- import("numpy", convert = FALSE)
bmc <- import("BMCToolkit", convert = FALSE)

run_bmc_autoK <- function(freq_mat,
                          trim_level = 0L) {
  freq_mat <- as.matrix(freq_mat)
  
  if (nrow(freq_mat) != ncol(freq_mat)) {
    stop("freq_mat must be a square matrix.")
  }
  if (any(freq_mat < 0)) {
    stop("freq_mat must be nonnegative.")
  }
  
  n_states <- nrow(freq_mat)
  traj_len <- sum(freq_mat)
  
  if (traj_len <= n_states) {
    stop("sum(freq_mat) must be larger than nrow(freq_mat).")
  }
  
  # Python numpy array로 명시적으로 변환
  freq_py <- np$array(freq_mat, dtype = "float64", order = "C")
  
  # 1) trimming
  trimmed_py <- bmc$trim_count_matrix(freq_py, as.integer(trim_level))
  
  # PyPI example parameter choices
  improve_iter <- ceiling(log(n_states))
  ratio <- traj_len / n_states
  radius <- sqrt(ratio^(1+0.9)/n_states)#ratio / sqrt(n_states * log(ratio))
  neighborhood_size_threshold <- n_states / (ratio)^(0.9-0.1) #n_states * log(ratio) / ratio
  singular_value_threshold <- ratio^0.75 #sqrt(ratio) * log(ratio)
  num_indices <- n_states #as.integer(ceiling(log(n_states)))
  
  # 2) K 자동추정
  k_hat_py <- bmc$compute_num_clusters(
    trimmed_py,
    radius,
    neighborhood_size_threshold,
    singular_value_threshold,
    num_indices
  )
  k_hat <- py_to_r(k_hat_py)
  
  # 3) spectral clustering
  cl_init_py <- bmc$compute_spectral_clustering(trimmed_py, as.integer(k_hat))
  
  # 4) cluster improvement
  cl_final_py <- bmc$compute_cluster_improvement(
    trimmed_py,
    cl_init_py,
    as.integer(improve_iter)
  )
  
  # 5) optional: BMC parameter estimate
  pars_py <- bmc$compute_bmcs_parameters(trimmed_py, cl_final_py)
  
  list(
    k_hat = k_hat,
    ratio = ratio,
    radius = radius,
    neighborhood_size_threshold = neighborhood_size_threshold,
    singular_value_threshold = singular_value_threshold,
    num_indices = num_indices,
    trimmed_matrix = py_to_r(trimmed_py),
    clustering_init = py_to_r(cl_init_py),
    clustering_final = py_to_r(cl_final_py),
    parameters = py_to_r(pars_py)
  )
}
