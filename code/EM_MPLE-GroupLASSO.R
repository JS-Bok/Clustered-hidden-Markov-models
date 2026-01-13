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

source("code/Auxiliary_Functions.R")
source("code/EM_MLE.R")

#### EM for MPLE ####
EM_MPLE_GLASSO_singleLambda <- function(y,A=NULL,dist_class,parameter=NULL,Lambda, tol=1e-5, m=NULL){
  
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
  EM_MPLE_step <- function(x, y, A, dist_class, parameter, Lambda, A_hat, prob_prev=NULL){
    
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
    
    # Pairwise-Differential matrix
    D <- diff_mat(m, (m-1))
    tDD <- t(D) %*% D
    
    # initial alpha and beta in ADMM
    a <- D %*% A_vec
    b <- rep(0, times = (m-1)*m*(m-1)/2)
    
    # Step for ADMM
    ADMM_step <- 0
    
    # SCAD parameter
    gamma <- 3.7
    
    # Step size for ADMM
    rho <- max(50*Lambda,1.1/(gamma-1))
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
      while (m*(m-1)/LB_par > 1e-14) {
        NM_step <- 1
        e2 <- 1
        
        while (e2 > 1e-15) {
          NM_step <- NM_step+1
          Update_ind <- which(A_vec > 1e-6)
          # A_vec[A_vec<1e-6] <- 1e-30
          if(NM_step == 2){
            grad_value <- grad_LB(A_vec, U,V_sum, tot_time,LB_par)+rho*(tDD%*%A_vec - Dab)
            hess <- hessian(A_vec, V_sum, rho, LB_par,m,tot_time)
            grad_value <- grad_value[Update_ind]
            hess <- hess[Update_ind, Update_ind]
            hess_inv <- tryCatch({chol2inv(chol(hess))},
                                 error = function(e) { tryCatch({armaInv(hess)},
                                                                error = function(e) { return(inv_svd(hess)) }
                                 ) }
            )
            hess_grad <- hess_inv %*% grad_value
          }
          if (sum(grad_value * hess_grad)/2 < 1e-15) { break }
          
          BLS <- 1
          t <- t_0
          BLSTEP <- 0
          pre_gfunc <- g_func_ori(A_vec, U,V_sum,tot_time) +rho/2*n2_square(D %*% A_vec-a+b) + LB_term(A_vec,LB_par,m)
          while (BLS == 1) {
            BLSTEP <- BLSTEP+1
            A_vec_new <- A_vec
            A_vec_new[Update_ind] <-  A_vec[Update_ind] - t * as.vector(hess_grad)
            if (is.na(sum(log(vec2mat(A_vec_new,m))))) {
              t <- t*beta
              BLS <- 1
              
              next
            }
            # Backtracking line search
            LHS <- g_func_ori(A_vec_new, U,V_sum,tot_time)+rho/2*n2_square(D %*% A_vec_new-a+b) + LB_term(A_vec_new,LB_par,m)
            RHS <-  pre_gfunc-t*sum(grad_value*hess_grad)*0.5
            
            if(LHS<=RHS){
              BLS <- 0
              break
            }
            
            t <- t*beta
            
          }
          
          
          grad_new <- grad_LB(A_vec_new, U,V_sum, tot_time, LB_par)+rho*(tDD%*%A_vec_new - Dab)
          hess <- hessian(A_vec_new, V_sum, rho, LB_par,m,tot_time)
          grad_new <- grad_new[Update_ind]
          hess <- hess[Update_ind, Update_ind]
          # hess_inv <- tryCatch({chol2inv(chol(hess))},
          #                      error = function(e) { return(inv_svd(hess)) }
          # )
          hess_inv <- tryCatch({chol2inv(chol(hess))},
                               error = function(e) { tryCatch({armaInv(hess)},
                                                              error = function(e) { return(inv_svd(hess)) }
                               ) }
          )
          hess_grad <- hess_inv %*% grad_new
          # Stopping criteria for Newton's method
          e2 <- sum(grad_new * hess_grad)/2
          if(t < 1e-8){
            break
          }
          
          
          e2_old <- e2
          A_vec_old <- A_vec
          A_vec <- A_vec_new
          grad_value <- grad_new
          if (NM_step > 1e1) {
            err <- 0
            # print("Larger than 1e1 step")
            break}
        }
        LB_par <- LB_par * 100
      }
      
      NM_step_seq <- append(NM_step_seq, (NM_step-1))
      
      ### Update a ###
      a_old <- a
      GTO_variable <- (D %*% A_vec + b)
      
      for (i in c(1:(m*(m-1)/2))) {
        a[(1+(m-1)*(i-1)):((m-1)*i)] <- GTO(GTO_variable[(1+(m-1)*(i-1)):((m-1)*i)], Lambda/rho)
      }
      
      ### Update b ###
      b <- b + D %*% A_vec - a
      
      # Values for Stopping criteria
      r <- D %*% A_vec - a
      s <- rho * t(D) %*% (a-a_old)
      
      e1_pri <- sqrt(nrow(D))*1e-6 + 1e-3 * max(n2_sqrt(D %*% A_vec), n2_sqrt(a))
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
    Pen <- 0
    SCAD_variable <- D %*% A_vec
    for (i in c(1:(m*(m-1)/2))) {
      Pen <- Pen + n2_sqrt(SCAD_variable[(1+(m-1)*(i-1)):((m-1)*i)])*Lambda
    }
    prob_prev <- lapply(y, em_prob_mat, x=x, dist_class=dist_class, parameter=parameter)
    LH <- mLHf(y, A, dist_class, parameter, prob_prev)/tot_time
    pLH <- LH - Pen
    setTxtProgressBar(pb, 6)
    cat(paste(" Scaled Likelihood: ", LH, ", Penalized Likelihood: ", pLH, sep = " "))
    L <- list(prior=x, par=parameter, transition=A, 
              pLH=pLH, LH=LH, ADMM_step=ADMM_step, 
              NM_step=NM_step_seq, dist_c=dist_class, 
              prob_prev = prob_prev)
    return(L)
  }
  
  
  par_mat <- data.frame(matrix(nrow = length(x)))
  LH_mat <- data.frame(matrix(nrow = length(x)))
  Start_time <- Sys.time()
  for (z in c(1:30)) {
    
    
    if (z == 1){
      PSHMM_new <- EM_MPLE_step(x,y,A,dist_class,parameter, Lambda, A_hat)
      End_time <- Sys.time()
      print(End_time - Start_time)
      PSHMM_new$RT <- (End_time - Start_time)
      PSHMM_old <- PSHMM_new
    } else {
      PSHMM_new <- EM_MPLE_step(x,y,A,dist_class,parameter, Lambda, A_hat, PSHMM_old$prob_prev)
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
  PSHMM_new$BIC <- BICscore(y, PSHMM_new$transition, 
                            PSHMM_new$dist_c, 
                            PSHMM_new$par)
  output <- list(
    observation = y,
    transition = cf(PSHMM_new$transition)$rc_trans,
    par = PSHMM_new$par,
    LH = mLHf(y, cf(PSHMM_new$transition)$rc_trans, PSHMM_new$dist_c, PSHMM_new$par),
    BIC = PSHMM_new$BIC,
    dist_class = PSHMM_new$dist_c,
    Lambda = Lambda
  )
  
  return(output)
}

EM_MPLE_GLASSO <- function(y,
                    A=NULL,
                    dist_class,
                    parameter=NULL, 
                    tol=1e-5,
                    m=NULL,
                    Lambda_vec = 2^seq(-7.8,-3,0.1) ){
  
  initial_hmm <- EM_MLE(y,A,dist_class,parameter, tol=1e-6, m)
  
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

  CHMM_result <- NULL
  threshold <- 0.5 * log(tot_time) * (m-1)^2  # 임계값 계산
  

    for (lam in Lambda_vec) {
      tmp <- tryCatch({
        EM_MPLE_GLASSO_singleLambda(y=y, A=initial_hmm$transition, dist_class=dist_class, parameter=initial_hmm$par, lam, tol = tol)
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

  
  if (is.null(CHMM_result)) {
    stop("CHMM_GLASSO fitting failed for all lambda values")
  }
  
  return(CHMM_result)
}


