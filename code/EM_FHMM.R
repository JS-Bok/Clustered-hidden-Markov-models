#######################################
#### EM_algorithms for MLE of FHMM ####
#######################################
# This algorithm follows the EM-algorithm with exact E-step, introduced by (Ghahramani & Jordan, 1995)
# fhmm is a EM algorithm for MLE of FHMM.

### Input elements###
# X is a list of data, the matrix with a size (length of data) * (dimension of data).
# M is a number of chains.
# K is a number of states per chain
# cyc is a maximum number of cycles of EM steps
# tol is a tolerance of the error of log-likelihood.

### Output elements ###
# Mu = mean vectors
# Cov = output covariance matrix
# P = transition matrices
# LL = log-likelihood curve

# Import required libraries
library(MASS)
library(GMCM)
source("code/Auxiliary_Functions.R")

#### Main Function ####


fhmm <- function(X, M = 2, K = 3, cyc = 1000, tol = 1e-5) {
  #### Auxiliary functions for FHMM ####
  fhmm2hmm <- function(x, M, K){
    index_mat <- sapply(1:M, function(x){rep(1:K, K^(x-1), each=K^(M-x))})
    index_x <- cbind(index_mat, x)
    hmm_x <- as.vector(sapply(1:M, function(x){
      sapply(1:K, function(y){sum(index_x[index_x[,x]==y,(M+1)])})
    }))
    return(hmm_x)
  }
  
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
  
  # Generate MK*MK Pairwise expectation from K^M U-vector
  PwExpectation <- function(U_vec,M,K){
    index_mat <- sapply(1:M, function(x){rep(1:K, K^(x-1), each=K^(M-x))})
    U_index <- cbind(index_mat, U_vec)
    PE <- diag(as.vector(sapply(1:M, function(x){
      sapply(1:K, function(y){sum(U_index[U_index[,x]==y,(M+1)])})
    })))
    local_mat <- matrix(0,ncol = K, nrow = K)
    for (i in 2:M) {
      for (j in 1:(i-1)) {
        
        for ( ik in 1:K) {
          for (jk in 1:K) {
            local_mat[ik,jk] <- sum(U_index[U_index[,i]==ik & U_index[,j]==jk,(M+1)])
          }
        }
        PE[ (K*(i-1)+1):(K*i), (K*(j-1)+1):(K*j) ] <- local_mat
        PE[ (K*(j-1)+1):(K*j), (K*(i-1)+1):(K*i) ] <- t(local_mat)
      }
    }
    return(PE)
  }
  
  # EM algorithm ####
  if (!is.matrix(X)) { 
    X <- matrix(X)
  }
  if (!is.list(X)) { 
    X <- list(X)
  }
  y <- X
  p <- ncol(X[[1]])
  ydim <- p
  N <- length(unlist(X))/p
  m <- K^M
  if(p==1){
    dist_class <- "norm"
  }else{
    dist_class <- "mvnorm"
  }
  
  
  # if (is.null(T)) T <- N
  # 
  # if (N %% T != 0) {
  #   stop("Error: Data matrix length must be multiple of sequence length T")
  # }
  ### Initialize
  XX <- Reduce(rbind, X)
  if(p==1){
    dist_class <- "norm"
  }else{
    dist_class <- "mvnorm"
  }
  
  if(dist_class=="mvnorm"){
    parameter <- init_par_mv(X, m, dist_class)$emission
  }else{
    parameter <- init_par(X, m, dist_class)$emission
  }
  
  Pi <- matrix(1, nrow = M, ncol = K)
  Pi <- t(Pi / rowSums(Pi))
  
  P_list <- list()
  # P_list[[length(P_list)+1]] <- diag(K)+0.1
  # P_list[[length(P_list)+1]] <- matrix(1,ncol = K,nrow = K)
  # if(K>2){
  # for (i in 3:K) {
  #   P_list[[length(P_list)+1]] <- matrix(runif(K^2,min = 0,max = 1),nrow = K)
  # }
  # }
  
  for (i in 1:M) {
    P_list[[length(P_list)+1]] <- 0.8 * diag(K) * M + rep(0.2 / K, K)
  }
  # for (i in 1:M) {
  #   P_list[[length(P_list)+1]] <- matrix(runif(K^2,min = 0,max = 1),nrow = K)
  # }
  P <- Reduce(rbind,P_list)
  P <- P / rowSums(P)
  
  LL <- numeric()
  lik <- 0
  
  
  
  for (cycle in 1:cyc) {
    #### FORWARD-BACKWARD ### EXACT E STEP
    
    ############################
    Pi_hmm <- exp(hmm2fhmm_par(log(as.vector(Pi)), M, K))
    
    A <- P[(K*(M-1)+1):(K*M),]
    for (i in (M-1):1) {
      A <- kronecker(P[(K*(i-1)+1):(K*i),], A)
    }
    Pi_hmm[Pi_hmm==0] <- 1e-10 # Prevent Machine Zero Occurrence
    
    em_mat_list <- lapply(X, em_prob_mat, x=Pi_hmm, dist_class=dist_class, parameter=parameter)
    ### Calculate expectations ###
    
    f_mat_list <- mapply(function(yy, em) {Log_f_algorithm(Pi_hmm, yy, A, em)}, yy=X, em=em_mat_list,SIMPLIFY = FALSE)
    b_mat_list <- mapply(function(yy, em) {Log_b_algorithm(Pi_hmm, yy, A, em)}, yy=X, em=em_mat_list,SIMPLIFY = FALSE)
    
    U <- list()
    V_sum <- list()
    tot_time <- N
    for (obs in c(1:length(X))) {
      sub_time <- length(unlist(X[[obs]]))/p
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
    
    
    
    
    
    
    ############################
    
    
    #### Calculate Likelihood and determine convergence
    
    oldlik <- lik
    lik <- Reduce('+',lapply(f_mat_list, function(x){logsum(x[,ncol(x)])}))/tot_time
    LL <- c(LL, lik)
    message(sprintf("cycle %i log likelihood = %f ", cycle, lik))
    
    if (cycle <= 2) {
      likbase <- lik
    } else if (is.nan(lik)){
      break
    } else if (abs(lik-oldlik)/abs(oldlik)<tol) {
      break
    }
    
    #### M STEP 
    
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
    
    # transition matrix 
    V_sum_new <- apply(V_sum,MARGIN=2,fhmm2hmm,M=M,K=K)
    V_sum_new <- apply(V_sum_new,MARGIN=1,fhmm2hmm,M=M,K=K)
    
    
    Pi_tmp <- Reduce('+',lapply(U,function(x){x[,1]}))
    for (i in 1:M) {
      P[(K*(i-1)+1):(K*i),] <- V_sum_new[(K*(i-1)+1):(K*i),(K*(i-1)+1):(K*i)]
      Pi[,i] <- Pi_tmp[(K*(i-1)+1):(K*i)]
    }
    
    P <- P/rowSums(P)
    
    Pi <- apply(Pi, 2, function(x){exp(log(x)-logsum(log(x)))})
  }
  
  return(list(par=parameter, P = P, Pi = Pi, LL = LL))
}


