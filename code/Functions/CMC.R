#######################################
#### EM_algorithms for MPLE of CMC ####
#######################################

source("code/functions/Auxiliary_Functions.R")
source("code/functions/EM_MLE.R")

#### MPLE of CMC ####
CMC_MPLE_singleLambda <- function(y,m = NULL, return_y=TRUE, penalty_type='group_fused_scad', Lambda){
  
  if (!is.list(y)) { 
    y <- list(y)
  }
  
  tot_time <- length(unlist(y))
  
  if (is.null(m)){
    m <- max(unlist(y))
  }
  
  
  N_trans <- lapply(y, function(yy){
    emb <- embed(yy,2)
    table(emb[,2],emb[,1])
  }) %>% Reduce(f='+')
    
  N_init <- lapply(y, table) %>% Reduce(f='+') %>% list
  
  A <- N_trans/rowSums(N_trans)
  x <- stationary_dist(A)
  
    
  ### maximizing A ###
  # A_vec is a m*(m-1) vector consist of the rows of the transition matrix except last column.
  
  A_vec <- as.vector(t(A[,c(1:(m-1))]))
  
  if (penalty_type %in% c("group_fused_scad", "fused_scad", "group_fused_lasso")) {
    
    # Pairwise-Differential matrix
    D <- diff_mat(m, (m-1))
    tDD <- t(D) %*% D
    
    # initial alpha and beta in ADMM
    a <- D %*% A_vec
    b <- rep(0, length = nrow(D))
    
  } else if (penalty_type %in% c("group_scad", "scad")) {
    
    D <- diag(m*(m-1))
    tDD <- diag(m*(m-1))
    
    # initial alpha and beta in ADMM
    a <- A_vec
    b <- rep(0, length = m*(m-1))
    
  } else {
    stop("Unknown penalty_type")
  }
  
  
  # Step for ADMM
  ADMM_step <- 0
  
  # SCAD parameter
  gamma <- 3.7
  
  # Step size for ADMM
  rho <- max(5*Lambda,1.1/(gamma-1))
  NM_step_seq <- seq()
  A_vec_bADMM <- A_vec
  
  # Perform ADMM
  for (ad in c(1:100)) {
    ADMM_step <- ADMM_step+1
    NM_step <- 1 
    e2 <- 1
    
    A_vec_previous <- A_vec
    A_vec_old <- A_vec
    
    # Learning rate of the backtracking in newton's method 
    beta <- 0.5
    # Step size of Log-Barrier method
    LB_par <- 1e10
    # Initial step size of Log-Barrier method
    t_0 <- 1
    # Pre-calculated value for efficient calculation.
    Dab <- t(D) %*% (a-b)
    
    ### Update A_vec ###
    if(ad==1){
      DA <- as.numeric(D %*% A_vec)   # cache D %*% A_vec
    }
    
    
    while (m*(m-1)/LB_par > 1e-14) {
      
      NM_step <- 1
      e2 <- 1
      
      while (e2 > 1e-15) {
        NM_step <- NM_step + 1
        Update_ind <- which(A_vec > 1e-6)
        
        if (NM_step == 2) {
          # Smoothing
          A_smooth <- matrix(A_vec,ncol=(m-1), byrow = TRUE)
          for (i in 1:m) {
            if (sum(A_smooth[i,])>=1){
              A_smooth[i,] <- pmax(A_smooth[i,]-(1e-6),0)
            }
          }
          A_vec <- as.vector(t(A_smooth))
          
          grad_value <- grad_LB(A_vec, N_init, N_trans, tot_time, LB_par) + rho*(tDD %*% A_vec - Dab)
          hess <- hessian(A_vec, N_trans, rho, LB_par, m, tot_time)  
          
          grad_value <- grad_value[Update_ind]
          hess <- hess[Update_ind, Update_ind, drop = FALSE]
          
          hess_grad <- tryCatch({
            R <- chol(hess)
            backsolve(R, forwardsolve(t(R), grad_value))
          }, error = function(e) {
            tryCatch({
              solve(hess, grad_value)
            }, error = function(e) {
              hess_inv <- tryCatch({chol2inv(chol(hess))},
                                   error = function(e) { tryCatch({inv_svd(hess)},
                                                                  error = function(e) { return(armaInv(hess)) }
                                   ) }
              )
              hess_inv %*% grad_value
            }
            )
            
          })
        }
        
        step <- as.numeric(hess_grad)               # Newton direction on active set
        gtd  <- 0.5 * sum(grad_value * step)        # RHS
        
        if (gtd < 1e-15) break
        
        # full step vector (inactive=0) : D %*% step_full one time
        step_full <- numeric(length(A_vec))
        step_full[Update_ind] <- step
        
        # residual r = D A_vec - a + b
        r <- DA - a + b
        r2 <- drop(crossprod(r))                    # ||r||^2
        
        # u = D step_full
        u <- as.numeric(D %*% step_full)
        ru <- drop(crossprod(r, u))                 # r^T u
        u2 <- drop(crossprod(u))                    # ||u||^2
        
        A_vec <- pmax(A_vec,1e-300)
        pre_gfunc <- g_func_ori(A_vec, N_init, N_trans, tot_time) +
          (rho/2) * r2 +
          LB_term(A_vec, LB_par, m)
        
        accepted <- FALSE
        t <- t_0
        beta <- 0.5
        
        repeat {
          A_vec_new <- A_vec
          A_vec_new[Update_ind] <- A_vec[Update_ind] - t * step
          
          pen2 <- r2 - 2*t*ru + (t*t)*u2
          if (pen2 < 0) pen2 <- 0
          
          LHS <- g_func_ori(A_vec_new, N_init, N_trans, tot_time) +
            (rho/2) * pen2 +
            LB_term(A_vec_new, LB_par, m)
          
          RHS <- pre_gfunc - t * gtd
          
          if (is.finite(LHS) && (LHS <= RHS)) {
            accepted <- TRUE
            break
          }
          
          t <- t * beta
          if (t < 1e-12) break
        }
        
        if (!accepted) break   
        
        # accept step 
        A_vec <- A_vec_new
        DA <- DA - t * u
        
        # Smoothing
        A_smooth <- matrix(A_vec,ncol=(m-1), byrow = TRUE)
        for (i in 1:m) {
          if (sum(A_smooth[i,])>=1){
            A_smooth[i,] <- pmax(A_smooth[i,]-(1e-6),0)
          }
        }
        A_vec <- as.vector(t(A_smooth))
        
        grad_new <- grad_LB(A_vec, N_init, N_trans, tot_time, LB_par) + rho*(tDD %*% A_vec - Dab)
        hess <- hessian(A_vec, N_trans, rho, LB_par, m, tot_time)
        
        grad_new <- grad_new[Update_ind]
        hess <- hess[Update_ind, Update_ind, drop = FALSE]
        
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
        
        e2 <- 0.5 * sum(grad_new * as.numeric(hess_grad))
        grad_value <- grad_new
        
        if (NM_step > 1e1) break
      }
      
      LB_par <- LB_par * 100
    }
    
    NM_step_seq <- append(NM_step_seq, (NM_step-1))
    
    ### Update a ###
    DA <- as.numeric(D %*% A_vec)
    a_old <- a
    GTO_variable <- (DA + b)
    
    if (penalty_type == "group_fused_scad") {
      for (i in c(1:(m*(m-1)/2))) {
        a[(1+(m-1)*(i-1)):((m-1)*i)] <-
          update_a(GTO_variable[(1+(m-1)*(i-1)):((m-1)*i)],
                   rho, gamma, Lambda)
      }
    } else if (penalty_type == "fused_scad") {
      ## fused SCAD (entrywise)
      a <- vapply(GTO_variable,
                  FUN = update_a,
                  FUN.VALUE = numeric(1),
                  rho = rho, gamma = gamma, Lambda = Lambda)
    } else if (penalty_type == "group_scad") {
      for (r in 1:m) {
        idx <- (1+(m-1)*(r-1)):((m-1)*r)
        a[idx] <- update_a(GTO_variable[idx], rho, gamma, Lambda)
      }
      
    } else if (penalty_type == "scad") {
      a <- vapply(GTO_variable, update_a, numeric(1), rho=rho, gamma=gamma, Lambda=Lambda)
      
    } else if (penalty_type == "group_fused_lasso") {
      for (i in c(1:(m*(m-1)/2))) {
        a[(1+(m-1)*(i-1)):((m-1)*i)] <- GTO(GTO_variable[(1+(m-1)*(i-1)):((m-1)*i)], Lambda/rho)
      }
    }
    
    ### Update b ###
    b <- b + DA - a
    
    # Values for Stopping criteria
    r <- DA - a
    s <- rho * t(D) %*% (a-a_old)
    
    e1_pri <- sqrt(nrow(D))*1e-6 + 1e-3 * max(n2_sqrt(DA), n2_sqrt(a))
    e1_dual <- sqrt(ncol(D))*1e-6 + 1e-3 * n2_sqrt(t(D) %*% (rho*b))
    er <- n2_sqrt(r)
    es <-n2_sqrt(s)
    
    # #Adaptive Penalty Parameter(B. Wohlberg 2017)
    if (ADMM_step < 50){
      if( er > es*10){
        rho <- rho*1.5
      }
      if (es > er*10){
        rho <- max(rho/1.5,1.1/(gamma-1))
      }
    }
    # After 50 steps, rho is fixed.
    
    #Stopping criteria
    if ((er < e1_pri && es < e1_dual)) {break}
  }
  # Obatin updated A
  A <- matrix(A_vec,ncol=(m-1), byrow = TRUE)
  B <- 1-matrix(apply(A, MARGIN = 1, sum), nrow=m)
  A <- cbind(A, B)
  
  
  ### calculate pi ###
  x <- stationary_dist(A)
  
  ### calculate Likelihood ###
  SCAD_variable <- D %*% A_vec
  if (penalty_type == "group_fused_scad") {
    Pen <- 0
    for (i in c(1:(m*(m-1)/2))) {
      Pen <- Pen + SCAD_pen(
        n2_sqrt(SCAD_variable[(1+(m-1)*(i-1)):((m-1)*i)]),
        Lambda, gamma
      )
    }
  } else if (penalty_type == "fused_scad") {
    Pen <- sum(vapply(abs(SCAD_variable),
                      FUN = SCAD_pen,
                      FUN.VALUE = numeric(1),
                      Lambda = Lambda, gamma = gamma))
  } else if (penalty_type == "group_scad") {
    
    Pen <- 0
    for (i in 1:m) {
      idx <- (1+(m-1)*(i-1)):((m-1)*i)
      Pen <- Pen + SCAD_pen(n2_sqrt(SCAD_variable[idx]), Lambda, gamma)
    }
    
  } else if (penalty_type == "scad") {
    
    Pen <- sum(vapply(abs(SCAD_variable), SCAD_pen, numeric(1), Lambda=Lambda, gamma=gamma))
    
  } else if (penalty_type == "group_fused_lasso") {
    Pen <- 0
    for (i in c(1:(m*(m-1)/2))) {
      Pen <- Pen + n2_sqrt(SCAD_variable[(1+(m-1)*(i-1)):((m-1)*i)])*Lambda
    }
    
  }
  
  A_safe <- pmax(A, 1e-300)
  
  LH <- sum(N_trans * log(A_safe))
  
  pLH <- LH - Pen
  
  Start_time <- Sys.time()
  print(Lambda)
  if (penalty_type %in% c("group_fused_scad", "fused_scad", "group_fused_lasso")) {
    groups <- cluster_from_a(a, m = nrow(A))
    A <- refit_A_by_groups_mean(A, groups)
  }
  
  
  BIC <- -2*sum(N_trans * log(A)) + (length(cf(A)$index) * (m-1)) * log(tot_time)
  
  End_time <- Sys.time()
  output <- list(
    observation = y,
    transition = A,
    LH = sum(N_trans * log(A)),
    BIC = BIC,
    Lambda = Lambda,
    time = End_time - Start_time
  )
  
  if(!(return_y)){
    output <- list(
      transition = A,
      LH = sum(N_trans * log(A)),
      BIC = BIC,
      Lambda = Lambda,
      time = End_time - Start_time
    )
  }
  
  return(output)
}

CMC_MPLE_SaveAllLambda <- function(y,
                                  Lambda_vec = 2^seq(-7.9,-2,0.1),
                                  penalty_type = "group_fused_scad"
){
  k=0
  CMC_result <- list()
  for (lam in Lambda_vec) {
    tmp <- tryCatch({
      CMC_MPLE_singleLambda(y=y, 
                           Lambda=lam,
                           return_y=FALSE,
                           penalty_type = penalty_type)
    }, error = function(e) { NULL }) 
    if (!is.null(tmp)) {
      k = k+1
      CMC_result[[k]] <- tmp
    }
  }
  if (length(CMC_result)==0) {
    print("CHMM fitting failed for all lambda values.")
  }
  return(CMC_result)
}





######################################
#### EM_algorithms for MLE of CMC ####
######################################

#### EM for MLE of CMC ####
CMC_MLE <- function(y, m=NULL){
  
  if (!is.list(y)) { 
    y <- list(y)
  }
  
  tot_time <- length(unlist(y))
  
  if (is.null(m)){
    m <- max(unlist(y))
  }
  
  
  N_trans <- lapply(y, function(yy){
    emb <- embed(yy,2)
    table(emb[,2],emb[,1])
  }) %>% Reduce(f='+')
  
  N_init <- lapply(y, table) %>% Reduce(f='+') %>% list
  
  A <- N_trans/rowSums(N_trans)
  x <- stationary_dist(A)
  
  Start_time <- Sys.time()

  ############# E-step ############
  ydim <- 1
  
  ########### M-step ###########
  
  ### maximizing parameter ###
  
  ### maximizing A ###
  # A_vec is a m*(m-1) vector consist of the rows of the transition matrix except last column.
  A_vec <- as.vector(t(A[,c(1:(m-1))]))
  A_vec_previous <- A_vec
  A_vec_old <- A_vec
  
  # Learning rate of the backtracking in newton's method 
  beta <- 0.5
  # Step size of Log-Barrier method
  LB_par <- 1e6
  # Initial step size of Log-Barrier method
  t_0 <- 1
  while (m*(m-1)/LB_par > 1e-12) {
    NM_step <- 1
    e2 <- 1
    # Newton's method
    while (e2 > 1e-15) {
      NM_step <- NM_step+1
      
      if(NM_step == 2){
        # Smoothing last column
        A_smooth <- matrix(A_vec,ncol=(m-1), byrow = TRUE)
        for (i in 1:m) {
          if (sum(A_smooth[i,])>=1){
            A_smooth[i,] <- pmax(A_smooth[i,]-(1e-6),0)
          }
        }
        A_vec <- as.vector(t(A_smooth))
        
        
        Update_ind <- which(A_vec > 1e-6)
        grad_value <- grad_LB(A_vec, N_init,N_trans, tot_time,LB_par)[Update_ind]
        hess <- hessian_nonpen(A_vec, N_trans, LB_par,m,tot_time)[Update_ind,Update_ind]
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
      if (sum(grad_value * hess_grad)/2 < 1e-15) { break }
      BLS <- 1
      t <- t_0
      BLSTEP <- 0
      A_vec <- pmax(A_vec,1e-300)
      pre_gfunc <- g_func_ori(A_vec, N_init,N_trans,tot_time)+ LB_term(A_vec,LB_par,m)
      while (BLS == 1) {
        BLSTEP <- BLSTEP+1
        A_vec_new <- A_vec
        A_vec_new[Update_ind] <- A_vec[Update_ind] - t * as.vector(hess_grad)
        if (t < 1e-12) break
        if (is.na(sum(log(vec2mat(A_vec_new,m))))) {
          t <- t*beta
          BLS <- 1
          next
        }
        # Backtracking line search
        LHS <- g_func_ori(A_vec_new, N_init,N_trans,tot_time) + LB_term(A_vec_new,LB_par,m)
        RHS <- pre_gfunc-t*sum(grad_value*hess_grad)*0.5
        
        if(LHS<=RHS){
          BLS <- 0
        }
        t <- t*beta
      }
      
      # Smoothing
      A_smooth <- matrix(A_vec_new,ncol=(m-1), byrow = TRUE)
      for (i in 1:m) {
        if (sum(A_smooth[i,])>=1){
          A_smooth[i,] <- pmax(A_smooth[i,]-(1e-6),0)
        }
      }
      A_vec_new <- as.vector(t(A_smooth))
      
      grad_new <- grad_LB(A_vec_new, N_init,N_trans, tot_time, LB_par)[Update_ind]
      hess <- hessian_nonpen(A_vec_new, N_trans, LB_par,m,tot_time)[Update_ind,Update_ind]
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
      #Stopping criteria for newton's method
      e2 <- sum(grad_new * hess_grad)/2
      
      A_vec_old <- A_vec
      A_vec <- A_vec_new
      grad_value <- grad_new
      if (NM_step > 1e3) {
        break
      }
    }
    LB_par <- LB_par * 100
  }
  
  
  
  # Obatin updated A
  A <- matrix(A_vec,ncol=(m-1), byrow = TRUE)
  B <- 1-matrix(apply(A, MARGIN = 1, sum), nrow=m)
  A <- cbind(A, B)
  
  
  ### calculate pi ###
  x <- stationary_dist(A)
  
  ### calculate Likelihood ###

  A_safe <- pmax(A, 1e-300)
  
  LH <- sum(N_trans * log(A_safe))
  
  
  End_time <- Sys.time()
  
  output <- list(
    observation = y,
    transition = A,
    LH = LH,
    time = End_time - Start_time
  )
  print(End_time - Start_time)
  return(output)
}

######################################
#### EM_algorithms for OE of CMC ####
######################################

#### EM for OE ####
CMC_Oracle <- function(y,A_true, m=NULL){
  
  ## Extract row-tying structure from A_true
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
  
  ## Build mapping for grad/hess aggregation (S)
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
  
  ## Define automatic origin<->oracle transforms
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
  
  if (!is.list(y)) { 
    y <- list(y)
  }
  
  tot_time <- length(unlist(y))
  
  if (is.null(m)){
    m <- max(unlist(y))
  }
  m_loc <- ncol(A_true)
  
  N_trans <- lapply(y, function(yy){
    emb <- embed(yy,2)
    table(emb[,2],emb[,1])
  }) %>% Reduce(f='+')
  
  N_init <- lapply(y, table) %>% Reduce(f='+') %>% list
  
  A <- N_trans/rowSums(N_trans)
  
  x <- stationary_dist(A)
  
  Start_time <- Sys.time()
  
  ############# E-step ############
  ydim <- 1
  
  ########### M-step ###########
  
  ### maximizing parameter ###
  
  ### maximizing A ###
  # A_vec is a m*(m-1) vector consist of the rows of the transition matrix except last column.
  A_raw <- N_trans / rowSums(N_trans)
  A_vec_raw <- as.vector(t(A_raw[, 1:(m_loc - 1)]))
  A_vec <- oracle2origin(origin2oracle_pick(A_vec_raw))
  
  beta <- 0.5
  LB_par <- 1e6
  t_0 <- 1
  d <- G * (m_loc - 1)
  
  while (d / LB_par > 1e-12){
    PG_step <- 1
    e2 <- 1
    
    while (e2 > 1e-15){
      PG_step <- PG_step+1
      
      if (PG_step == 2){
        # Smoothing
        A_smooth <- matrix(A_vec,ncol=(m-1), byrow = TRUE)
        for (i in 1:m) {
          A_smooth[i,] <- pmax(A_smooth[i,],1e-8)
          if (sum(A_smooth[i,])>=1){
            A_smooth[i,] <- pmax(A_smooth[i,]-(1e-6),1e-8)
          }
        }
        A_vec <- as.vector(t(A_smooth))
        
        
        grad_value <- origin2oracle_sum(grad_LB(A_vec, N_init, N_trans, tot_time, LB_par))
        hess <- origin2oracle_hess(hessian_nonpen(A_vec, N_trans, LB_par, m_loc, tot_time))
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
      
      A_vec <- pmax(A_vec,1e-300)
      pre_gfunc <- g_func_ori(A_vec, N_init, N_trans, tot_time) + LB_term(A_vec, LB_par, m_loc)
      
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
        
        LHS <- g_func_ori(A_vec_new, N_init, N_trans, tot_time) + LB_term(A_vec_new, LB_par, m_loc)
        RHS <- pre_gfunc - t * sum(grad_value * hess_grad) * 0.5
        
        if (LHS <= RHS) BLS <- 0
        if (t < 1e-12) break
        t <- t * beta
      }
      # Smoothing
      A_smooth <- matrix(A_vec_new,ncol=(m-1), byrow = TRUE)
      for (i in 1:m) {
        if (sum(A_smooth[i,])>=1){
          A_smooth[i,] <- pmax(A_smooth[i,]-(1e-6),1e-8)
        }
      }
      A_vec_new <- as.vector(t(A_smooth))
      
      grad_new <- origin2oracle_sum(grad_LB(A_vec_new, N_init, N_trans, tot_time, LB_par))
      hess <- origin2oracle_hess(hessian_nonpen(A_vec_new, N_trans, LB_par, m_loc, tot_time))
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
  
  
  ### calculate pi ###
  x <- stationary_dist(A)
  
  ### calculate Likelihood ###
  
  A_safe <- pmax(A, 1e-300)
  
  LH <- sum(N_trans * log(A_safe))
  
  
  End_time <- Sys.time()
  
  output <- list(
    observation = y,
    transition = A,
    LH = LH,
    time = End_time - Start_time
  )
  print(End_time - Start_time)
  return(output)
}