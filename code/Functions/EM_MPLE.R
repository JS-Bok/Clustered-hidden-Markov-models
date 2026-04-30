################################
#### EM_algorithms for MPLE ####
################################
### Input elements###
# y is a list of data, the matrix with a size (length of data) * (dimension of data).
# A is a initial transision matrix.
# dist_class is a distribution class of emission probability.
# parameter is a parameter of the emission distribution.
# tol is a tolerance of the relative error of log-likelihood.
# m is a number of hidden state.
# Lambda_vec is a ordered vector of penalty parameters that will be compared with BIC.

### Output elements ###
# observation = data.
# transition = transition matrix
# par = parameter of the emission probabilities
# LH = log-likelihood
# BIC = BIC score
# dist_class = distribution class of the emission distribution
# Lambda = selected penalty parameter

source("code/functions/Auxiliary_Functions.R")
source("code/functions/EM_MLE.R")

#### EM for MPLE ####
EM_MPLE_singleLambda <- function(y,A=NULL,dist_class,parameter=NULL,Lambda, tol=1e-5, max_step=30,m=NULL, return_y=TRUE, penalty_type){
  
  if (!is.list(y)) { 
    y <- list(y)
  }
  if(is.null(A) || is.null(parameter)){
    if(is.null(m)){
      print("Initial parameters or Number of states is necessary.")
      break
    }
    if(dist_class=="mvnorm"){
      A <- init_par_mv(y, m, dist_class)$transition
      parameter <- init_par_mv(y, m, dist_class)$emission
    }else{
      A <- init_par(y, m, dist_class)$transition
      parameter <- init_par(y, m, dist_class)$emission
    }
  }
  A_hat <- A
  PSHMM <- list()
  x <- stationary_dist(A)
  
  ##### EM Algorithm with penalized Stationary HMM ####
  EM_MPLE_step <- function(x, y, A, dist_class, parameter, Lambda, A_hat, prob_prev=NULL, penalty_type){
    
    pb <- txtProgressBar(min = 0,
                         max = 6,
                         style = 3,
                         width = 50,
                         char = "="
    )
    ############# E-step ############
    x <- stationary_dist(A)
    m <- ncol(A)
    U <- list()
    V_sum <- list()
    tot_time <- length(unlist(y))
    ydim <- 1
    if (dist_class == "mvnorm"){
      ydim <- ncol(y[[1]])
      tot_time <- length(unlist(y))/ydim
    }
    
    
    setTxtProgressBar(pb, 1)
    cat(" E_step: Calculating Expected Emission Probabilities")
    
    if(is.null(prob_prev)){
      em_mat_list <- lapply(y, em_prob_mat, x=x, dist_class=dist_class, parameter=parameter)
    } else {
      em_mat_list <- prob_prev
    }
    
    setTxtProgressBar(pb, 2)
    cat(" E-step: Forward Backward Algotithm                         ")
    f_mat_list <- mapply(function(yy, em) {Log_f_algorithm(x, yy, A, em)}, yy=y, em=em_mat_list,SIMPLIFY = FALSE)
    b_mat_list <- mapply(function(yy, em) {Log_b_algorithm(x, yy, A, em)}, yy=y, em=em_mat_list,SIMPLIFY = FALSE)
    for (obs in c(1:length(y))) {
      sub_time <- length(unlist(y[[obs]]))/ydim
      em_mat <- em_mat_list[[obs]]
      f_mat <- f_mat_list[[obs]]
      b_mat <- b_mat_list[[obs]]
      
      
      Uj <- matrix(NA, nrow = m, ncol = sub_time)
      fb_mat1 <- f_mat+b_mat
      for (t in c(1:{sub_time-1})) {
        Uj[,t] <- exp(fb_mat1[,t]-logsum(fb_mat1[,t]))
        
        junk <- exp(matrix(rep(f_mat[,t],m), byrow = FALSE, nrow=m) 
                    + matrix(rep(b_mat[,{t+1}]+log(em_mat[,{t+1}]),m),byrow = TRUE, nrow=m)
                    +log(A)-logsum(fb_mat1[,t]))
        V_sum[[length(V_sum)+1]] <- junk
      }
      Uj[,sub_time] <- exp(fb_mat1[,sub_time]-logsum(fb_mat1[,sub_time]))
      U[[length(U)+1]] <- Uj
    }
    V_sum <- Reduce('+',V_sum)
    
    
    
    ########### M-step ###########
    
    setTxtProgressBar(pb, 3)
    cat(" M-step: Maximizing emission parameters                        ")
    ### maximizing parameter ###
    
    if(dist_class == "pois"){
      new_par <- list()
      for (i in c(1:m)){
        UY_sum <- matrix(0,nrow = ydim, ncol = 1)
        
        U_sum <- Reduce('+', lapply(U, function(X){return(sum(X[i,]))}))
        
        UY_sum <- UY_sum+sum((mapply(function(X,Y) { return(sum(X* Y[i,])) }, X=y, Y=U)))
        
        
        new_par[[i]] <- list(mean=as.numeric(UY_sum/U_sum))
      }
      parameter <- new_par
      U_junk <- list(matrix(0,nrow=m, ncol=1))
      for (el in U) {
        U_junk[[1]] <- U_junk[[1]]+el[,1]
      }
      U <- U_junk
    }
    
    
    if(dist_class == "norm"){
      new_par <- list()
      for (i in c(1:m)){
        UY_sum <- matrix(0,nrow = ydim, ncol = 1)
        UY2_sum <- matrix(0,nrow = ydim, ncol = ydim)
        
        U_sum <- Reduce('+', lapply(U, function(X){return(sum(X[i,]))}))
        
        UY_sum <- UY_sum+sum((mapply(function(X,Y) { return(sum(X* Y[i,])) }, X=y, Y=U)))
        
        junk_mu <- UY_sum/U_sum
        
        UY2_sum <- Reduce('+', 
                          mapply(function(X,Y) {
                            return(Reduce('+',mapply('*', as.list(Y[i,]), 
                                                     lapply(as.list(data.frame(t(X))), function(xx) {return(((junk_mu-matrix(xx))^2))}), 
                                                     SIMPLIFY = FALSE)))
                          }, X=y, Y=U, SIMPLIFY = FALSE)
        )
        new_par[[i]] <- list(mean=as.numeric(UY_sum/U_sum), sd=as.numeric(sqrt(UY2_sum/U_sum)))
      }
      parameter <- new_par
      U_junk <- list(matrix(0,nrow=m, ncol=1))
      for (el in U) {
        U_junk[[1]] <- U_junk[[1]]+el[,1]
      }
      U <- U_junk
    }
    
    if(dist_class == "mvnorm"){
      new_par <- list()
      for (i in c(1:m)){
        UY_sum <- matrix(0,nrow = ydim, ncol = 1)
        UY2_sum <- matrix(0,nrow = ydim, ncol = ydim)
        
        U_sum <- Reduce('+', lapply(U, function(X){return(sum(X[i,]))}))
        
        UY_sum <- UY_sum+rowSums((mapply(function(X,Y) { return(rowSums(t(X)* matrix( rep(Y[i,], time=ydim), nrow = ydim, byrow = TRUE))) }, X=y, Y=U)))
        
        junk_mu <- UY_sum/U_sum
        
        UY2_sum <- Reduce('+', 
                          mapply(function(X,Y) {
                            return(Reduce('+',mapply('*', as.list(Y[i,]), 
                                                     lapply(as.list(data.frame(t(X))), function(xx) {return(((junk_mu-matrix(xx)) %*% t(junk_mu-matrix(xx))))}), 
                                                     SIMPLIFY = FALSE)))
                          }, X=y, Y=U, SIMPLIFY = FALSE)
        )
        new_par[[i]] <- list(mean=UY_sum/U_sum, sigma=UY2_sum/U_sum)
      }
      parameter <- new_par
      U_junk <- list(matrix(0,nrow=m, ncol=1))
      for (el in U) {
        U_junk[[1]] <- U_junk[[1]]+el[,1]
      }
      U <- U_junk
    }
    
    
    ### maximizing A ###
    setTxtProgressBar(pb, 4)
    cat(" M-step: Maximizing transition matrix               ")
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
            grad_value <- grad_LB(A_vec, U, V_sum, tot_time, LB_par) + rho*(tDD %*% A_vec - Dab)
            hess <- hessian(A_vec, V_sum, rho, LB_par, m, tot_time)  
            
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
          
          pre_gfunc <- g_func_ori(A_vec, U, V_sum, tot_time) +
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
            
            LHS <- g_func_ori(A_vec_new, U, V_sum, tot_time) +
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
          
          grad_new <- grad_LB(A_vec, U, V_sum, tot_time, LB_par) + rho*(tDD %*% A_vec - Dab)
          hess <- hessian(A_vec, V_sum, rho, LB_par, m, tot_time)
          
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
    
    setTxtProgressBar(pb, 5)
    cat(" Calculating Likelihood                          ")
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
    prob_prev <- lapply(y, em_prob_mat, x=x, dist_class=dist_class, parameter=parameter)
    LH <- mLHf(y, A, dist_class, parameter, prob_prev)/tot_time
    pLH <- LH - Pen
    setTxtProgressBar(pb, 6)
    cat(paste(" Scaled Likelihood: ", LH, ", Penalized Likelihood: ", pLH, sep = " "))
    L <- list(prior=x, par=parameter, transition=A, 
              pLH=pLH, LH=LH, ADMM_step=ADMM_step, 
              NM_step=NM_step_seq, dist_c=dist_class, 
              prob_prev = prob_prev,
              a=a)
    return(L)
  }
  
  
  par_mat <- data.frame(matrix(nrow = length(x)))
  LH_mat <- data.frame(matrix(nrow = length(x)))
  Start_time <- Sys.time()
  for (z in c(1:max_step)) {
    
    
    if (z == 1){
      PSHMM_new <- EM_MPLE_step(x,y,A,dist_class,parameter, Lambda, A_hat, penalty_type=penalty_type)
      End_time <- Sys.time()
      print(End_time - Start_time)
      PSHMM_new$RT <- (End_time - Start_time)
      PSHMM_old <- PSHMM_new
    } else {
      PSHMM_new <- EM_MPLE_step(x,y,A,dist_class,parameter, Lambda, A_hat, PSHMM_old$prob_prev, penalty_type=penalty_type)
      End_time <- Sys.time()
      print(End_time - Start_time)
      PSHMM_new$RT <- (End_time - Start_time)
    }
    PSHMM <- append(PSHMM, list(PSHMM_new))
    prob_prev <- PSHMM_new$prob_prev
    x <- PSHMM_new$prior
    parameter <- PSHMM_new$par
    A <- PSHMM_new$transition
    if (z == 1){
      next
    }
    if( abs(PSHMM_new$pLH - PSHMM_old$pLH)/abs(PSHMM_new$pLH) <= tol) {
      break
    }
    PSHMM_old <- PSHMM_new
  }
  print(Lambda)
  if (penalty_type %in% c("group_fused_scad", "fused_scad", "group_fused_lasso")) {
    groups <- cluster_from_a(PSHMM_new$a, m = nrow(PSHMM_new$transition))
    A <- refit_A_by_groups_mean(PSHMM_new$transition, groups)
  } else {
    A <- PSHMM_new$transition
  }
  PSHMM_new$BIC <- BICscore(y, A, 
                            PSHMM_new$dist_c, 
                            PSHMM_new$par)
  End_time <- Sys.time()
  output <- list(
    observation = y,
    transition = A,
    par = PSHMM_new$par,
    LH = mLHf(y, A, PSHMM_new$dist_c, PSHMM_new$par),
    BIC = PSHMM_new$BIC,
    dist_class = PSHMM_new$dist_c,
    Lambda = Lambda,
    time = End_time - Start_time
  )
  
  if(!(return_y)){
    output <- list(
      transition = A,
      par = PSHMM_new$par,
      LH = mLHf(y, A, PSHMM_new$dist_c, PSHMM_new$par),
      BIC = PSHMM_new$BIC,
      dist_class = PSHMM_new$dist_c,
      Lambda = Lambda,
      time = End_time - Start_time
    )
  }
  
  return(output)
}

EM_MPLE <- function(y,
                    A=NULL,
                    dist_class,
                    parameter=NULL, 
                    tol=1e-5,
                    m=NULL,
                    Lambda_vec = 2^seq(-7.9,-2,0.1),
                    penalty_type = "group_fused_scad",
                    initial_hmm = FALSE){
  if(!all(c("LH", "par", "transition") %in% names(initial_hmm))){
    initial_hmm <- EM_MLE(y,A,dist_class,parameter, tol=1e-6, m)
  }

  # Initialize with HMM
  HMM_LH <- initial_hmm$LH
  m <- nrow(initial_hmm$transition)
  if (!is.list(y)) { 
    y <- list(y)
  }
  tot_time <- length(unlist(y))
  if (dist_class == "mvnorm"){
    tot_time <- tot_time/ncol(y[[1]])
  }
  
  initial_hmm$BIC <- BICscore(y, initial_hmm$transition, 
                              initial_hmm$dist_class, 
                              initial_hmm$par)
  
  initial_hmm$Lambda <- 0
  CHMM_result <- NULL
  threshold <- 0.5 * log(tot_time) * (m-1)^2 
  
  Start_time <- Sys.time()
    for (lam in Lambda_vec) {
      tmp <- tryCatch({
        EM_MPLE_singleLambda(y=y, 
                             A=initial_hmm$transition, 
                             dist_class=dist_class, 
                             parameter=initial_hmm$par, 
                             Lambda = lam, tol = tol,
                             penalty_type = penalty_type)
      }, error = function(e) { NULL }) 
      
      if (!is.null(tmp)) {
        LH_diff <- abs(tmp$LH - HMM_LH)  
        if (LH_diff > threshold) {  
          break
        }
        
        if (is.null(CHMM_result) || CHMM_result$BIC > tmp$BIC) {
          CHMM_result <- tmp
        }
      }
    }
  End_time <- Sys.time()
  
  
  
  if (is.null(CHMM_result)) {
    print("CHMM fitting failed for all lambda values. return MLE.")
    CHMM_result <- initial_hmm
  }
  
  CHMM_result$time <- End_time - Start_time
  return(CHMM_result)
}


EM_MPLE_SaveAllLambda <- function(y,
                            A=NULL,
                            dist_class,
                            parameter=NULL, 
                            tol=1e-5,
                            m=NULL,
                            Lambda_vec = 2^seq(-7.9,-2,0.1),
                            penalty_type = "group_fused_scad",
                            initial_hmm = FALSE
){
  
  if(!all(c("LH", "par", "transition") %in% names(initial_hmm))){
    initial_hmm <- EM_MLE(y,A,dist_class,parameter, tol=1e-6, m)
  }
  
  # Initialize with HMM
  HMM_LH <- initial_hmm$LH
  m <- nrow(initial_hmm$transition)
  if (!is.list(y)) { 
    y <- list(y)
  }
  tot_time <- length(unlist(y))
  if (dist_class == "mvnorm"){
    tot_time <- tot_time/ncol(y[[1]])
  }
  
  CHMM_result <- list()
  k=0
  
  initial_hmm$BIC <- BICscore(y, initial_hmm$transition, 
                              initial_hmm$dist_class, 
                              initial_hmm$par)
  
  initial_hmm$Lambda <- 0
  for (lam in Lambda_vec) {
    tmp <- tryCatch({
      EM_MPLE_singleLambda(y=y, 
                           A=initial_hmm$transition, 
                           dist_class=dist_class, 
                           parameter=initial_hmm$par, 
                           lam, tol = tol,
                           return_y=FALSE,
                           penalty_type = "group_fused_scad")
    }, error = function(e) { NULL }) 
    if (!is.null(tmp)) {
      k = k+1
      CHMM_result[[k]] <- tmp
    }
  }
  if (length(CHMM_result)==0) {
    print("CHMM fitting failed for all lambda values. return MLE.")
    CHMM_result <- initial_hmm
  }
  return(CHMM_result)
}


