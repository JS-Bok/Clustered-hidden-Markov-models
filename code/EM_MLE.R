###############################
#### EM_algorithms for MLE ####
###############################
### Input elements###
# y is a list of data, the matrix with a size (length of data) * (dimension of data).
# A is a initial transition matrix.
# dist_class is a distribution class of emission probability.
# parameter is a parameter of the emission distribution.
# tol is a tolerance of the relative error of log-likelihood.

### Output elements ###
# observation = data.
# transition = transition matrix
# par = parameter of the emission probabilities
# LH = log-likelihood
# dist_class = distribution class of the emission distribution


source("code/Auxiliary_Functions.R")

#### EM for MLE ####
EM_MLE <- function(y,A=NULL,dist_class,parameter=NULL, tol=1e-6, m=NULL){
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
  
  
  x <- stationary_dist(A)
  Start_time <- Sys.time()
  
  EM_MLE_step <- function(y, A, dist_class, parameter, prob_prev=NULL){
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
          grad_value <- grad_LB(A_vec, U,V_sum, tot_time,LB_par)
          hess <- hessian_nonpen(A_vec, V_sum, LB_par,m,tot_time)
          hess_inv <- tryCatch({chol2inv(chol(hess))},
                               error = function(e) { return(inv_svd(hess)) }
          )
          hess_grad <- hess_inv %*% grad_value
        }
        if (sum(grad_value * hess_grad)/2 < 1e-15) { break }
        BLS <- 1
        t <- t_0
        BLSTEP <- 0
        pre_gfunc <- g_func_ori(A_vec, U,V_sum,tot_time)+ LB_term(A_vec,LB_par,m)
        while (BLS == 1) {
          BLSTEP <- BLSTEP+1
          A_vec_new <- A_vec - t * as.vector(hess_grad)
          if (is.na(sum(log(vec2mat(A_vec_new,m))))) {
            t <- t*beta
            BLS <- 1
            next
          }
          # Backtracking line search
          LHS <- g_func_ori(A_vec_new, U,V_sum,tot_time) + LB_term(A_vec_new,LB_par,m)
          RHS <- pre_gfunc-t*sum(grad_value*hess_grad)*0.5
          
          if(LHS<=RHS){
            BLS <- 0
          }
          t <- t*beta
        }
        
        grad_new <- grad_LB(A_vec_new, U,V_sum, tot_time, LB_par)
        hess <- hessian_nonpen(A_vec_new, V_sum, LB_par,m,tot_time)
        hess_inv <- tryCatch({chol2inv(chol(hess))},
                             error = function(e) { return(inv_svd(hess)) }
        )
        hess_grad <- hess_inv %*% grad_new
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
    
    setTxtProgressBar(pb, 5)
    cat(" Calculating Likelihood                          ")
    ### calculate Likelihood ###
    
    
    prob_prev <- lapply(y, em_prob_mat, x=x, dist_class=dist_class, parameter=parameter)
    LH <- mLHf(y, A, dist_class, parameter, prob_prev)/tot_time
    setTxtProgressBar(pb, 6)
    cat(paste(" Scaled Likelihood: ", LH, "                        ", sep = " "))
    close(pb)
    
    L <- list(prior=x, par=parameter, transition=A, Likelihood=LH, dist_c=dist_class, observation=y, prob_prev = prob_prev)
    return(L)
    
  }
  
  for (z in c(1:1000)) {
    cat(paste('EM_step:', z,'\n', sep = " "))
    
    if (z == 1){
      BWHMM_new <- EM_MLE_step(y,A,dist_class,parameter)
      BWHMM_old <- BWHMM_new
    } else {
      BWHMM_new <- EM_MLE_step(y,A,dist_class,parameter, BWHMM_old$prob_prev)
    }
    
    
    x <- BWHMM_new$prior
    parameter <- BWHMM_new$par
    A <- BWHMM_new$transition
    
    
    if (z == 1){
      next
    }
    if(abs(BWHMM_new$Likelihood - BWHMM_old$Likelihood)/abs(BWHMM_old$Likelihood) <= tol) {
      
      break
    }
    BWHMM_old <- BWHMM_new
  }
  
  End_time <- Sys.time()
  
  output <- list(
    observation = y,
    transition = cf(BWHMM_new$transition)$rc_trans,
    par = BWHMM_new$par,
    LH = mLHf(y, cf(BWHMM_new$transition)$rc_trans, BWHMM_new$dist_c, BWHMM_new$par),
    dist_class = BWHMM_new$dist_c
  )
  print(End_time - Start_time)
  return(output)
}



