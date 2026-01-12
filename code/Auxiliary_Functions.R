##### Required packages ####
library('Rcpp')
library('RcppEigen')
library('tibble')
library('RcppHMM')
library('RcppArmadillo')
library('SimDesign')

#### Matrix inversion ####
# Matrix inversion using svd
inv_svd <- function(mat) {
  mat_svd <- svd(mat)
  return(mat_svd$v %*% diag(1/(mat_svd$d)) %*% t(mat_svd$u))
}

# Matrix inversion using Rcpp
Rcpp::cppFunction("arma::mat armaInv(arma::mat x) { return arma::inv(x); }", depends="RcppArmadillo")

# Matrix inversion using SVD
inv_svd <- function(mat) {
  mat_svd <- svd(mat)
  return(mat_svd$v %*% diag(1/(mat_svd$d)) %*% t(mat_svd$u))
}

#### Norm function ####
n2_square <- function(x) {sum(x^2)}
n2_sqrt <- function(x) {sqrt(sum(x^2))}


#### Log-Sum function ############
# logxpy(log(x),log(y))=log(x+y)
logxpy <- function(lx,ly) max(lx,ly) + log1p(exp(-abs(lx-ly)))

logsum <- function(lc) max(lc) + log(sum(exp(lc-max(lc))))


#### Pairwise difference matrix generating function ####
diff_mat <- function(m,n){
  D1 <- matrix(0, nrow = n*m*(m-1)/2, ncol = m*n)
  D2 <- matrix(0, nrow = n*m*(m-1)/2, ncol = m*n)
  for (i in c(0:{m-2})) {
    for (j in c(0:{m-i-2})) {
      D1[c({n*j+n*i*(m-(i+1)/2)+1}:{n*j+n*i*(m-(i+1)/2)+n}),c({i*n+1}:{i*n+n})] <- diag(n)
    }
  }
  for (i in c(0:{m-2})) {
    for (j in c(0:{m-i-2})) {
      D2[c({n*j+n*i*(m-(i+1)/2)+1}:{n*j+n*i*(m-(i+1)/2)+n}),c({n*(i+j+1)+1}:{n*(i+j+1)+n})] <- diag(n)
    }
  }
  return(D1-D2)
}


#### Functions for optimization ####
# A_vec is a old point,
# U and V_sum are matrix calculated in E-step,
# LB_par is a "t" in Log-barrier method,
# rho is a penalty parameter of ADMM,
# m is a number of states,
# tot_time is a total length of observation sequence.



# Gradient function for the non-penalized estimation
grad <- function(A_vec, U, V_sum, t){
  m <- nrow(V_sum)
  A1 <- vec2mat(A_vec, m)
  junk <- -V_sum/A1+V_sum[,m]/(A1[,m])
  gd2 <- as.vector(t(junk[c(1:m),c(1:(m-1))]))
  
  
  Pi <- stationary_dist(A1)
  gd1 <- rep(0, length = (m*(m-1)))
  junk1 <- Matrix::solve(diag(m)-t(A1)+1)
  
  for (i in c(1:m)) {
    for (j in c(1:(m-1))) {
      junk2 <- rep(0, length = m)
      junk2[i] <- (Pi[j]-Pi[m])
      gd1[(m-1)*(i-1)+j] <- t(U[[1]]) %*% (-( junk1 %*% junk2 )/Pi)
    }
  }
  return((gd1+gd2)/t)
}

# Gradient function for the non-penalized estimation
grad_LB<- function(A_vec, U, V_sum, t, LB_par){
  m <- nrow(V_sum)
  A1 <- vec2mat(A_vec, m)
  junk <- -V_sum/A1+V_sum[,m]/(A1[,m])
  gd2 <- as.vector(t(junk[c(1:m),c(1:(m-1))]))
  LB_junk <- -1/A1+1/(A1[,m])
  gd3 <- as.vector(t(LB_junk[c(1:m),c(1:(m-1))]))
  
  gd1 <- rep(0, length = (m*(m-1)))
  junk1 <- Matrix::solve(diag(m)-t(A1)+1)
  Pi <- junk1 %*% matrix(1, nrow=ncol(junk1), ncol=1)
  for (i in c(1:m)) {
    for (j in c(1:(m-1))) {
      junk2 <- rep(0, length = m)
      junk2[j] <- -Pi[i]
      junk2[m] <- Pi[i]
      gd1[(m-1)*(i-1)+j] <- -t(U[[1]]) %*% (-( junk1 %*% junk2 )/Pi)
    }
  }
  return((gd1+gd2)/t+gd3/LB_par)
}

# Calculating the value of objective function, used for the backtracking line search for newton's method
g_func_ori <- function(A_vec, U, V_sum, t){
  m <- nrow(V_sum)
  A1 <- vec2mat(A_vec, m)
  g2 <- -sum(V_sum * log(A1))
  
  Pi <- stationary_dist(A1)
  g1 <- -t(U[[1]]) %*% log(Pi)
  
  return((g1+g2)/t)
}

# Calculating value of Log-barrier term.
LB_term <- function(A_vec, LB_par,m){
  A1 <- vec2mat(A_vec, m)
  return(-sum(log(A1))/LB_par)
}

# Calculating the hessian matrix for the newton's method(non-penalized estimation)
hessian_nonpen <- function(A_vec, V_sum, LB_par,m,tot_time) {
  hess <- matrix(0, ncol = m*(m-1), nrow = m*(m-1))
  A1 <- vec2mat(A_vec,m)
  for (i1 in c(1:m)) {
    for (j1 in c(1:(m-1))) {
      for (i2 in c(1:m)) {
        for (j2 in c(1:(m-1))) {
          if ( i1 == i2 && j1 == j2){
            hess[((i1-1)*(m-1)+j1),((i2-1)*(m-1)+j2)] <- (V_sum[i1,j1]/tot_time+1/LB_par)/(A1[i1,j1]^2)+(V_sum[i1,m]/tot_time+1/LB_par)/(A1[i1,m]^2)
          } else if (i1 == i2){
            hess[((i1-1)*(m-1)+j1),((i2-1)*(m-1)+j2)] <- (V_sum[i1,m]/tot_time+1/LB_par)/(A1[i1,m]^2)
          }
        }
      }
    }
  }
  return(hess)
}

# Calculating the hessian matrix for the newton's method(penalized estimation)
hessian <- function(A_vec, V_sum, rho, LB_par,m,tot_time) {
  hess <- matrix(0, ncol = m*(m-1), nrow = m*(m-1))
  A1 <- vec2mat(A_vec,m)
  for (i1 in c(1:m)) {
    for (j1 in c(1:(m-1))) {
      for (i2 in c(1:m)) {
        for (j2 in c(1:(m-1))) {
          if ( i1 == i2 && j1 == j2){
            hess[((i1-1)*(m-1)+j1),((i2-1)*(m-1)+j2)] <- (V_sum[i1,j1]/tot_time+1/LB_par)/(A1[i1,j1]^2)+(V_sum[i1,m]/tot_time+1/LB_par)/(A1[i1,m]^2)+rho*(m-1)
          } else if (i1 == i2){
            hess[((i1-1)*(m-1)+j1),((i2-1)*(m-1)+j2)] <- (V_sum[i1,m]/tot_time+1/LB_par)/(A1[i1,m]^2)
          } else if (j1 == j2){
            hess[((i1-1)*(m-1)+j1),((i2-1)*(m-1)+j2)] <- -rho
          }
          
        }
      }
    }
  }
  return(hess)
}




#### Functions for ADMM ####
# Group Thresholding Operater
GTO <- function(x,t) {
  if (norm(x, type="2") > 0){
    return(x * max(0, 1-(t/norm(x, type="2"))))
  }
  else{
    return(x)
  }
}

# Updating alpha in ADMM
update_a <- function(K,rho,gamma,Lambda){
  K_norm <- n2_sqrt(K)
  if (Lambda/rho > K_norm ){
    a_new <- 0
  }else if (K_norm<= Lambda*(1+1/rho)) {
    a_new <- (1-(Lambda/rho)/K_norm) * K 
  }else if ( K_norm > gamma*Lambda ) {
    a_new <- K
  } else {
    a_new <- (rho*(gamma-1)*K_norm - gamma*Lambda)/(rho*(gamma-1)-1) * K / K_norm
  }
  return(a_new)
}

# Calculating the value of SCAD function
SCAD_pen <- function(var_norm, Lambda, gamma){
  #var_norm <- norm(D%*%A_vec[group],type="2")
  if ( var_norm <= Lambda){
    pen_value <- Lambda * var_norm
  } else if (var_norm <= gamma*Lambda){
    pen_value <- (var_norm^2-2*gamma*Lambda*var_norm+Lambda^2)/(2*(1-gamma))
  } else {
    pen_value <- Lambda^2 * (gamma+1) /2
  }
  return(pen_value)
}

#### Stationary distribution ####
# Given a transition matrix, it returns the induced stationary distribution.
stationary_dist <- function(x) Matrix::solve( -t(x)+1+diag(nrow(x)) ) %*% matrix(1, nrow=nrow(x), ncol=1)

#### Inverse of Vectorization ####
# Given a vector x, vectorization of some transition matrix with last column removed, it returns the original transition matrix.
vec2mat <- function(x,m){
  A_junk <- matrix(x,nrow=m,ncol=(m-1), byrow = TRUE)
  B_junk <- 1-matrix(apply(A_junk, MARGIN = 1, sum), nrow=m)
  return(cbind(A_junk, B_junk))
}
#### Multivariate Normal Distribution ####
# X is a input vector, mean is a mean vector, and sigma_inv is a inverse of the covariance matrix
mvnorm_inv <- function(X,mean,sigma_inv) {
  k <- length(X)
  junk <- as.numeric(X-mean)
  return(exp(-(k/2)*log(2*pi) -1/2*(t(junk) %*% sigma_inv %*% (junk)) +1/2*log(det(sigma_inv))))
}

#### Emission probability matrix ####
# x is a initial distribution, 
# y is a data, 
# dist_class is a distribution class of the emission probability, 
# and parameter is a parameter of the emission probabilities.
em_prob_mat <- function(x, y, dist_class, parameter){
  if (dist_class == "pois"){
    y <- matrix(y, ncol = 1)
    junk <- function(li) { return(apply(y, MARGIN = 1, FUN = dpois, mean=li$mean)) }
    K <- t(sapply(parameter, junk))
  }
  
  
  
  if (dist_class == "norm"){
    y <- matrix(y, ncol = 1)
    junk <- function(li) { return(apply(y, MARGIN = 1, FUN = dnorm, mean=li$mean, sd=li$sd)) }
    K <- t(sapply(parameter, junk))
  }
  
  
  if (dist_class == "mvnorm"){
    par_inv <- lapply(parameter, function(li) {
      sigma_inv <- tryCatch(solve(li$sigma), error = function(e){inv_svd(li$sigma)})
      return(list(mean=li$mean, sigma_inv=sigma_inv))
    })
    junk <- function(li) { return(apply(y, MARGIN = 1, FUN = mvnorm_inv, mean=li$mean, sigma_inv=li$sigma_inv)) }
    K <- t(sapply(par_inv, junk))
  }
  
  K[K==0] <- .Machine$double.xmin*1e10
  K[which(is.na(K))] <- .Machine$double.xmin*1e10
  K[which(is.infinite(K) & K>0)] <- .Machine$double.xmax*1e-10 
  K[which(is.infinite(K) & K<0)] <- .Machine$double.xmin*1e10
  return(K)
  
}



#### Forward probability matrix ####
# The algorithm follows the book "Hidden Markov Models for Time Series"-W. Zucchini, et al.(2016)
# x is a initial distribution, 
# y is a data, 
# A is a transition matrix
# p_mat is a emission probability matrix
Log_f_algorithm <- function(x, y, A, p_mat){
  ### x is a initial probability of hidden states, and A is transition probability A[i,j]=P(X(t+1)=j \ X(t)=i) ###
  m <- nrow(A)
  n <- ncol(p_mat)
  f_m <- matrix(NA,m,n)
  foo <- t(x)*p_mat[,1]
  sumfoo <- sum(foo)
  f_scale <- log(sumfoo)
  foo <- foo/sumfoo
  f_m[,1] <- f_scale+log(foo)
  for (i in c(2:n)) {
    foo <- foo %*% A * p_mat[,i]
    sumfoo <- sum(foo)
    f_scale <- f_scale+log(sumfoo)
    foo <- foo/sumfoo
    f_m[,i] <- log(foo)+f_scale
  }
  return(f_m)
}




#### Backward probability matrix ####
# The algorithm follows the book "Hidden Markov Models for Time Series"-W. Zucchini, et al.(2016)
# x is a initial probability, 
# y is a data, 
# A is a transition matrix
# p_mat is a emission probability matrix
Log_b_algorithm <- function(x, y, A, p_mat){
  m <- nrow(A)
  n <- ncol(p_mat)
  b_m <- matrix(NA,m,n)
  b_m[,n] <- rep(0,m)
  foo <- rep(1/m,m)
  b_scale <- log(m)
  for (i in c((n-1):1)) {
    foo <- A %*% (p_mat[,(i+1)]*foo)
    b_m[,i] <- log(foo)+b_scale
    sumfoo <- sum(foo)
    foo <- foo/sumfoo
    b_scale <- b_scale+log(sumfoo)
  }
  return(b_m)
}









#### Initialize EM for MLE  ####
# y is a data,
# m is a number of staets,
# dist_c is a cistribution class
# Available input of dist_c in init_par is "norm" and "pois", and init_par_mv is "mvnorm".
init_par <- function(y, m, dist_c){
  if(is.list(y)){
    y <- unlist(y)
  }
  trans <- 0.8 * diag(m) + rep(0.2 / m, m)
  mean <- quantile(quantile(y, c(0.0001,0.9999), na.rm = FALSE,names = F, type = 7), probs = seq(0, 1, 1/(m+1)), na.rm = FALSE,names = F, type = 7)[2:(m+1)]
  para <- list()
  if (dist_c == "pois"){
    para <- list(mean=mean)
    for (i in c(1:m)) {
      para[[i]] <- list(mean = mean[i])
    }
  }
  
  
  if (dist_c == "norm"){
    sd<- rep(sd(y) / m, times = m)
    for (i in c(1:m)) {
      para[[i]] <- list(mean = mean[i], sd = sd[i])
    }
  }
  
  
  return(list(transition=trans,emission=para))
}


init_par_mv <- function(y, m, dist_c){
  if(is.list(y)){
    y <- Reduce(rbind,y)
  }
  d <- ncol(y)
  trans <- 0.8 * diag(m) + rep(0.2 / m, m)
  
  #mean
  mean <- matrix(NA, nrow = m, ncol = d)
  for (i in c(1:d)){
    mean[,i] <- quantile(quantile(as.vector(y[,i]), c(0.0001,0.9999), na.rm = FALSE,names = F, type = 7), probs = seq(0, 1, 1/(m+1)), na.rm = FALSE,names = F, type = 7)[2:(m+1)]
  }
  para <- list()
  
  
  if (dist_c == "mvnorm"){
    #sigma
    sigma <- diag(diag(var(y)))/m^2
    
    for (i in c(1:m)) {
      para[[i]] <- list(mean = mean[i,], sigma = sigma)
    }
  }
  
  return(list(transition=trans,emission=para))
}


#### Clustering function ####
# A is a transition matrix.
# It returns the indices of states of each clusters (ind), 
# the reconstructed transition matrix after clustering step (rc_trans),
# the transition matrix of clusters (transition),
# and the matrix of clustering probabilities (cluster).
cf <- function(A, th=1e-2) {
  m <- nrow(A)
  cr <- list()
  ind <- list()
  cr[[1]] <- matrix(A[1,], nrow = 1)
  ind[[1]] <- c(1)
  for (i in c(2:m)) {
    co <- 0
    for (j in c(1:length(cr))) {
      if (all(abs(cr[[j]][1,]-A[i,]) < th )) {
        cr[[j]] <- rbind(cr[[j]], A[i,])
        ind[[j]][(length(ind[[j]])+1)] = i
        co <- 1
      }
      if (co == 1){
        break
      }
    }
    if (co == 0){
      cr[[(length(cr)+1)]] <- matrix(A[i,], nrow = 1)
      ind[[(length(ind)+1)]] <- c(i)
    }
  }
  
  cm <- matrix(0, nrow = length(cr), ncol = m)
  for (i in c(1:length(cr))) {
    cm[i,] <- apply(cr[[i]],2,mean)
  }
  
  trans <- matrix(0, nrow = length(cr), ncol = length(cr))
  for (i in (1:length(cr))) {
    for (j in ind[[i]]) {
      trans[,i] <- trans[,i]+cm[,j]
    }
  }
  
  clust <- cm
  for (i in c(1:length(cr))) {
    clust[,ind[[i]]] <- clust[,ind[[i]]]/trans[,i]
  }
  
  rc_matrix <- matrix(0, ncol=m, nrow=m)
  for (i in c(1:length(ind))) {
    rc_matrix[ind[[i]],] <- matrix(rep(cm[i,], times = length(ind[[i]])), byrow = TRUE, ncol = m)
  }
  
  L = list(index=ind, rc_trans = rc_matrix,transition=trans, cluster=clust)
  return(L)
}

#### Function for TPR, FPR, TP, FP ####
# TI = cf(True transition matrix)$index
# EI = cf(Estimated transition matrix)$index
TFPN<- function(TI,EI){
  m <- max(unlist(TI))
  pa1 <- rep(0, length = m*(m-1)/2)
  pa2 <- rep(0, length = m*(m-1)/2)
  for (k in TI) {
    for (i in c(1:(m-1))) {
      for (j in c((i+1):m)) {
        if({i %in% k}&&{j %in% k}){
          pa1[(i-1)*m-i*(i+1)/2+j] <- 1
        }
      }
      
    }
  }
  for (k in EI) {
    for (i in c(1:(m-1))) {
      for (j in c((i+1):m)) {
        if({i %in% k}&&{j %in% k}){
          pa2[(i-1)*m-i*(i+1)/2+j] <- 1
        }
      }
      
    }
  }
  TP <- 0
  FP <- 0
  FN <- 0
  TN <- 0
  
  for (l in c(1:(m*(m-1)/2))) {
    if({pa1[l] == 0} && {pa2[l] == 0}){
      TP <- TP+1
    }
    if({pa1[l] == 1} && {pa2[l] == 0}){
      FP <- FP+1
    }
    if({pa1[l] == 0} && {pa2[l] == 1}){
      FN <- FN+1
    }
    if({pa1[l] == 1} && {pa2[l] == 1}){
      TN <- TN+1
    }
  }
  L = list(FPR=FP/(FP+TN), TPR=TP/(TP+FN), TP=TP, FP=FP)
  
  return(L)
  
}



#### Log-likelihood function ####
# y is a data,
# A is a transition matrix,
# dist_class is a distribution class,
# parameter is a parameter for emission probability, 
# and if there is a prob_prev, a pre-calculated probability distribution, we use them.
LHf <- function(y, A, dist_class, parameter, prob_prev=NULL, x=NULL) {
  if (is.null(x)){
    x <- stationary_dist(A)
  }
  
  if (is.null(prob_prev)){
    f_prob <- Log_f_algorithm(x, y, A, em_prob_mat(x, y, dist_class, parameter))
  } else {
    f_prob <- Log_f_algorithm(x, y, A, prob_prev)
  }
  LH <- logsum(as.vector(f_prob[,ncol(f_prob)]))
  return(LH)
}

# mLHf is a Likelihood function for multiple data
mLHf <- function(y, A, dist_class, parameter, prob_prev=NULL, x=NULL){
  LH <- 0
  
  junk_func <- function(junk) {LHf(junk, A, dist_class, parameter,x=x)}
  if (is.null(prob_prev)){
    junk_func <- function(junk) {LHf(junk, A, dist_class, parameter,x=x)}
    mLH <- Reduce('+', lapply(y, junk_func))
  } else {
    junk_func <- function(junk1,junk2) {LHf(junk1, A, dist_class, parameter , junk2,x=x)}
    mLH <- Reduce('+', mapply(junk_func, junk1=y, junk2=prob_prev, SIMPLIFY = FALSE))
  }
  
  return( mLH )
}



#### BIC score function ####
# iLH is a pre-calculated log-likelihood
BICscore <- function(y, A, dist_class, parameter, iLH=NULL){
  if (dist_class == 'mvnorm'){
    x <- stationary_dist(A)
    n <- length(unlist(y))/ncol(y[[1]])
    cf_A <- cf(A)
    if (is.null(iLH)){  LH <- mLHf(y, cf_A$rc_trans, dist_class, parameter) }
    else{LH <- iLH}
    clust_ind <- length(cf_A$index)
    df <- length(unlist(parameter)) + clust_ind * (ncol(A) - 1)
    return(-2*LH + df*log(n))
  }
  cf_A <- cf(A)
  x <- stationary_dist(A)
  n <- length(unlist(y))
  if (is.null(iLH)){  LH <- mLHf(y, cf_A$rc_trans, dist_class, parameter) }
  else{LH <- iLH}
  clust_ind <- length(cf_A$index)
  df <-length(unlist(parameter)) + clust_ind * (ncol(A) - 1)
  return(-2*LH + df*log(n))
}
#### HMM Simulation for Gaussian Distribution ####
# dimension = a dimension of the data.
# init_dist = initial distribution of the markov chain. (As a default, stationary distribution will be calculated.)
# transition = transition matrix.
# emission = list of parameters of emission probabilities, (mean vector, covariance matrix).
# size = length of data.
# Use SimDesign package for the sampling of multivariate gaussian distribution. 
HMM_simulation <- function(dimension, init_dist=NULL, transition, emission, size){
  m <- nrow(transition)
  
  # If there is no initial distribution, set as stationary distribution.
  if(is.null(init_dist)){
    init_dist <- Matrix::solve( -t(transition)+1+diag(m) ) %*% matrix(1, nrow=m, ncol=1) 
  }
  
  # Generate underlying markov sequence
  markov_chain <- sample(x = seq(1, m, by = 1), 1, prob = init_dist)
  for (i in 2:size) 
  {
    last_state <- markov_chain[i - 1]
    markov_chain <- c(markov_chain, sample(x = seq(1, m, by=1), 1, prob = transition[last_state, ]))
  }
  observation <- matrix(NA, nrow = size, ncol = dimension)
  
  # Generate observation sequence
  for (i in 1:size) 
  { 
    observation[i,] <- rmvnorm(n = 1, mean = emission[[markov_chain[i]]]$mean, sigma = emission[[markov_chain[i]]]$sigma)
  }			
  return(observation)
}

#### Cluster ordering function ####
Clust_Ordering <- function(A) {
  m <- nrow(A)
  cr <- list()
  gsf_order <- c(1)
  while (length(gsf_order) < m) {
    tmp <- apply(A, 1, function(x) { return(norm( A[gsf_order[length(gsf_order)],] - x,  type = "2"))})
    tmp[gsf_order] <- 100
    tmp_order <- which.min(tmp)
    gsf_order <- append(gsf_order, tmp_order)
  }
  A_new <- A[gsf_order, gsf_order]
  
  L = list(order=gsf_order, transition=A_new)
  return(L)
}
