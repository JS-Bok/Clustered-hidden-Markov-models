######################################################
#### EM_algorithms for Oracle Estimations of CHMM ####
######################################################

### Input elements###
# y is a list of data, the matrix with a size (length of data) * (dimension of data).
# A_true is a ture transition matrix (for identical row structure)
# A is a initial transition matrix.
# dist_class is a distribution class of emission probability.
# parameter is a parameter of the emission distribution.
# tol is a tolerance of the relative error of log-likelihood.
# m is a number of hidden state.

### Output elements ###
# observation = data.
# transition = transition matrix
# par = parameter of the emission probabilities
# LH = log-likelihood
# dist_class = distribution class of the emission distribution


source("code/functions/Auxiliary_Functions.R")


EM_Oracle <- function(
    y,
    A_true,            
    A = NULL,
    dist_class,
    parameter = NULL,
    tol = 1e-6,
    m = NULL,
    max_iter = 1000
){
  
  if (!is.list(y)) y <- list(y)
  
  ## ---- 1) Extract row-tying structure from A_true ----
  oracle_extract_groups <- function(A_true){
    m <- nrow(A_true)
    A_cmp <- A_true
    
    used <- rep(FALSE, m)
    groups <- list()
    
    for (i in seq_len(m)) {
      if (used[i]) next
      diffs <- apply(A_cmp, 1, function(r) max(abs(r - A_cmp[i, ])))
      idx <- which(diffs == 0)
      used[idx] <- TRUE
      groups[[length(groups) + 1]] <- sort(idx)
    }
    
    # group order by first appearance
    groups <- groups[order(sapply(groups, `[`, 1))]
    groups
  }
  
  groups <- oracle_extract_groups(A_true)
  G <- length(groups)
  p <- m - 1
  
  row_to_group <- integer(m)
  for (k in seq_along(groups)) row_to_group[groups[[k]]] <- k
  if (any(row_to_group == 0)) stop("Failed to build row_to_group mapping (some rows not assigned).")
  
  ## ---- 2) Build mapping for grad/hess aggregation (S) ----
  # grad_oracle = S %*% grad_full
  # Hess_oracle = S %*% Hess_full %*% t(S)
  S <- matrix(0, nrow = G * p, ncol = m * p)
  I_p <- diag(p)
  
  for (k in seq_along(groups)) {
    rr <- ((k - 1) * p + 1):(k * p)
    for (i in groups[[k]]) {
      cc <- ((i - 1) * p + 1):(i * p)
      S[rr, cc] <- S[rr, cc] + I_p
    }
  }
  
  ## ---- 3) Define automatic origin<->oracle transforms ----
  # "pick" oracle vector from full A_vec (rep row or mean row within each group)
  origin2oracle_pick <- function(vec_full){
    full_mat <- matrix(vec_full, nrow = m, ncol = p, byrow = TRUE)
  
    reps <- vapply(groups, function(g) g[1], integer(1))
    oracle_mat <- full_mat[reps, , drop = FALSE]
  
    as.vector(t(oracle_mat))  # row-wise vec
  }
  
  # expand oracle vector back to full vector with tied rows replicated
  oracle2origin <- function(vec_oracle){
    oracle_mat <- matrix(vec_oracle, nrow = G, ncol = p, byrow = TRUE)
    full_mat <- oracle_mat[row_to_group, , drop = FALSE]
    as.vector(t(full_mat))  # row-wise vec
  }
  
  origin2oracle_sum <- function(vec_full){
    as.vector(S %*% vec_full)
  }
  
  origin2oracle_hess <- function(mat_full){
    S %*% mat_full %*% t(S)
  }
  
  ## ---- 4) Initialize A / emission parameters if needed ----
  if (is.null(A)) {
    # generic safe init (positive + row-stochastic)
    A <- matrix(1/m, nrow = m, ncol = m)
  } else {
    if (!is.matrix(A) || nrow(A) != m || ncol(A) != m) stop("A must be an m x m matrix.")
  }
  
  if (is.null(parameter)) {
    if (dist_class == "mvnorm") {
      parameter <- init_par_mv(y, m, dist_class)$emission
    } else {
      parameter <- init_par(y, m, dist_class)$emission
    }
  }
  
  ## ---- 5) Oracle EM step (automatic structure) ----
  Oracle_step_auto <- function(y, A, dist_class, parameter, prob_prev = NULL){
    pb <- txtProgressBar(min = 0, max = 6, style = 3, width = 50, char = "=")
    
    ## E-step
    x <- stationary_dist(A)
    m_loc <- ncol(A)
    
    U <- list()
    V_sum <- list()
    
    tot_time <- length(unlist(y))
    ydim <- 1
    if (dist_class == "mvnorm"){
      ydim <- ncol(y[[1]])
      tot_time <- length(unlist(y)) / ydim
    }
    
    setTxtProgressBar(pb, 1)
    cat(" E_step: Calculating Expected Emission Probabilities")
    
    if (is.null(prob_prev)){
      em_mat_list <- lapply(y, em_prob_mat, x = x, dist_class = dist_class, parameter = parameter)
    } else {
      em_mat_list <- prob_prev
    }
    
    setTxtProgressBar(pb, 2)
    cat(" E-step: Forward Backward Algotithm                         ")
    
    f_mat_list <- mapply(function(yy, em) { Log_f_algorithm(x, yy, A, em) },
                         yy = y, em = em_mat_list, SIMPLIFY = FALSE)
    b_mat_list <- mapply(function(yy, em) { Log_b_algorithm(x, yy, A, em) },
                         yy = y, em = em_mat_list, SIMPLIFY = FALSE)
    
    for (obs in seq_along(y)){
      sub_time <- length(y[[obs]]) / ydim
      em_mat <- em_mat_list[[obs]]
      f_mat <- f_mat_list[[obs]]
      b_mat <- b_mat_list[[obs]]
      
      Uj <- matrix(NA, nrow = m_loc, ncol = sub_time)
      fb_mat1 <- f_mat + b_mat
      
      for (t in 1:(sub_time - 1)){
        Uj[, t] <- exp(fb_mat1[, t] - logsum(fb_mat1[, t]))
        
        junk <- exp(
          matrix(rep(f_mat[, t], m_loc), byrow = FALSE, nrow = m_loc) +
            matrix(rep(b_mat[, t + 1] + log(em_mat[, t + 1]), m_loc), byrow = TRUE, nrow = m_loc) +
            log(A) - logsum(fb_mat1[, t])
        )
        
        V_sum[[length(V_sum) + 1]] <- junk
      }
      
      Uj[, sub_time] <- exp(fb_mat1[, sub_time] - logsum(fb_mat1[, sub_time]))
      U[[length(U) + 1]] <- Uj
    }
    
    V_sum <- Reduce("+", V_sum)
    
    ## M-step: emission parameters
    setTxtProgressBar(pb, 3)
    cat(" M-step: Maximizing emission parameters                        ")
    
    if (dist_class == "pois"){
      new_par <- list()
      for (i in 1:m_loc){
        U_sum <- Reduce("+", lapply(U, function(X) sum(X[i, ])))
        UY_sum <- sum(mapply(function(X, Y) sum(X * Y[i, ]), X = y, Y = U))
        new_par[[i]] <- list(mean = as.numeric(UY_sum / U_sum))
      }
      parameter <- new_par
      
      # keep only initial gamma for pi-estimation usage downstream (same as your code)
      U_junk <- list(matrix(0, nrow = m_loc, ncol = 1))
      for (el in U) U_junk[[1]] <- U_junk[[1]] + el[, 1]
      U <- U_junk
    }
    
    if (dist_class == "norm"){
      new_par <- list()
      for (i in 1:m_loc){
        U_sum <- Reduce("+", lapply(U, function(X) sum(X[i, ])))
        UY_sum <- sum(mapply(function(X, Y) sum(X * Y[i, ]), X = y, Y = U))
        junk_mu <- UY_sum / U_sum
        
        UY2_sum <- Reduce("+",
                          mapply(function(X, Y){
                            Reduce("+", mapply("*", as.list(Y[i, ]),
                                               lapply(as.list(data.frame(t(X))), function(xx) ((junk_mu - matrix(xx))^2)),
                                               SIMPLIFY = FALSE))
                          }, X = y, Y = U, SIMPLIFY = FALSE)
        )
        
        new_par[[i]] <- list(
          mean = as.numeric(UY_sum / U_sum),
          sd   = as.numeric(sqrt(UY2_sum / U_sum))
        )
      }
      parameter <- new_par
      
      U_junk <- list(matrix(0, nrow = m_loc, ncol = 1))
      for (el in U) U_junk[[1]] <- U_junk[[1]] + el[, 1]
      U <- U_junk
    }
    
    if (dist_class == "mvnorm"){
      new_par <- list()
      for (i in 1:m_loc){
        U_sum <- Reduce("+", lapply(U, function(X) sum(X[i, ])))
        
        UY_sum <- rowSums(mapply(function(X, Y){
          rowSums(t(X) * matrix(rep(Y[i, ], time = ydim), nrow = ydim, byrow = TRUE))
        }, X = y, Y = U))
        
        junk_mu <- UY_sum / U_sum
        
        UY2_sum <- Reduce("+",
                          mapply(function(X, Y){
                            Reduce("+", mapply("*", as.list(Y[i, ]),
                                               lapply(as.list(data.frame(t(X))), function(xx){
                                                 (junk_mu - matrix(xx)) %*% t(junk_mu - matrix(xx))
                                               }),
                                               SIMPLIFY = FALSE))
                          }, X = y, Y = U, SIMPLIFY = FALSE)
        )
        
        new_par[[i]] <- list(mean = UY_sum / U_sum, sigma = UY2_sum / U_sum)
      }
      parameter <- new_par
      
      U_junk <- list(matrix(0, nrow = m_loc, ncol = 1))
      for (el in U) U_junk[[1]] <- U_junk[[1]] + el[, 1]
      U <- U_junk
    }
    
    ## M-step: transition matrix A (automatic oracle row-tying)
    setTxtProgressBar(pb, 4)
    cat(" M-step: Maximizing transition matrix               ")
    
    A_vec <- as.vector(t(A[, 1:(m_loc - 1)]))
    
    beta <- 0.5
    LB_par <- 1e6
    t_0 <- 1
    d <- G * (m_loc - 1)
    
    while (d / LB_par > 1e-12){
      PG_step <- 1
      e2 <- 1
      
      while (e2 > 1e-15){
        PG_step <- PG_step + 1
        
        if (PG_step == 2){
          grad_value <- origin2oracle_sum(grad_LB(A_vec, U, V_sum, tot_time, LB_par))
          hess <- origin2oracle_hess(hessian_nonpen(A_vec, V_sum, LB_par, m_loc, tot_time))
          hess_grad <- tryCatch({
            R <- chol(hess)
            backsolve(R, forwardsolve(t(R), grad_value))
          }, error = function(e) {
            tryCatch({
              solve(hess, grad_value)
            }, error = function(e) {
              hess_inv <- tryCatch({chol2inv(chol(hess))},
                                   error = function(e) { tryCatch({armaInv(hess)},
                                                                  error = function(e) { return(inv_svd(hess)) }
                                   ) }
              )
              hess_inv %*% grad_value
            }
            )
            
          })
        }
        
        if (sum(grad_value * hess_grad) / 2 < 1e-15) break
        
        pre_gfunc <- g_func_ori(A_vec, U, V_sum, tot_time) + LB_term(A_vec, LB_par, m_loc)
        
        BLS <- 1
        t <- t_0
        
        while (BLS == 1){
          A_oracle_now <- origin2oracle_pick(A_vec)
          A_oracle_new <- A_oracle_now - t * as.vector(hess_grad)
          A_vec_new <- oracle2origin(A_oracle_new)
          
          if (is.na(sum(log(vec2mat(A_vec_new, m_loc))))){
            t <- t * beta
            next
          }
          
          LHS <- g_func_ori(A_vec_new, U, V_sum, tot_time) + LB_term(A_vec_new, LB_par, m_loc)
          RHS <- pre_gfunc - t * sum(grad_value * hess_grad) * 0.5
          
          if (LHS <= RHS) BLS <- 0
          t <- t * beta
        }
        
        grad_new <- origin2oracle_sum(grad_LB(A_vec_new, U, V_sum, tot_time, LB_par))
        hess <- origin2oracle_hess(hessian_nonpen(A_vec_new, V_sum, LB_par, m_loc, tot_time))
        hess_grad <- tryCatch({
          R <- chol(hess)
          backsolve(R, forwardsolve(t(R), grad_new))
        }, error = function(e) {
          tryCatch({
            solve(hess, grad_new)
          }, error = function(e) {
            hess_inv <- tryCatch({chol2inv(chol(hess))},
                                 error = function(e) { tryCatch({armaInv(hess)},
                                                                error = function(e) { return(inv_svd(hess)) }
                                 ) }
            )
            hess_inv %*% grad_new
          }
          )
          
        })
        e2 <- sum(grad_new * hess_grad) / 2
        
        A_vec <- A_vec_new
        grad_value <- grad_new
        
        if (PG_step > 1e2) break
      }
      
      LB_par <- LB_par * 100
    }
    
    ## rebuild A
    A_sub <- matrix(A_vec, ncol = (m_loc - 1), byrow = TRUE)
    B <- 1 - matrix(apply(A_sub, 1, sum), nrow = m_loc)
    A <- cbind(A_sub, B)
    
    ## likelihood
    x <- stationary_dist(A)
    
    setTxtProgressBar(pb, 5)
    cat(" Calculating Likelihood                          ")
    
    prob_prev <- lapply(y, em_prob_mat, x = x, dist_class = dist_class, parameter = parameter)
    LH <- mLHf(y, A, dist_class, parameter, prob_prev) / tot_time
    
    setTxtProgressBar(pb, 6)
    cat(paste(" Scaled Likelihood: ", LH, "                        ", sep = " "))
    close(pb)
    
    list(prior = x, par = parameter, transition = A, Likelihood = LH,
         dist_c = dist_class, observation = y, prob_prev = prob_prev)
  }
  
  ## ---- 6) EM iterations ----
  Start_time <- Sys.time()
  
  for (z in 1:max_iter){
    cat(paste("EM_step:", z, "\n", sep = " "))
    
    if (z == 1){
      BWHMM_new <- Oracle_step_auto(y, A, dist_class, parameter)
      BWHMM_old <- BWHMM_new
    } else {
      BWHMM_new <- Oracle_step_auto(y, A, dist_class, parameter, BWHMM_old$prob_prev)
    }
    
    parameter <- BWHMM_new$par
    A <- BWHMM_new$transition
    
    if (z > 1){
      if (abs(BWHMM_new$Likelihood - BWHMM_old$Likelihood) / abs(BWHMM_old$Likelihood) <= tol) break
    }
    BWHMM_old <- BWHMM_new
  }
  
  End_time <- Sys.time()
  
  output <- list(
    observation = y,
    transition = cf(BWHMM_new$transition)$rc_trans,
    par = BWHMM_new$par,
    LH = mLHf(y, cf(BWHMM_new$transition)$rc_trans, BWHMM_new$dist_c, BWHMM_new$par),
    dist_class = BWHMM_new$dist_c,
    time = End_time - Start_time
  )
  
  print(End_time - Start_time)
  return(output)
}