################################################################################
#### Code for reproduce the numerical results of Protein structure analysis ####
################################################################################
#### Preliminary ####
# Load functions
sapply(c("code/functions/EM_MLE.R","code/functions/EM_MPLE.R","code/functions/EM_FHMM.R"), source)


# Load data
load("data/Protein_data.Rdata")
Protein_descriptor <- lapply(Protein_data, function(x) x$descriptor_vec)

#### MLE ####
# This code might take several hours to complete.
### Initialization ###
y <- Reduce(rbind,Protein_descriptor)
d <- ncol(y)
m=27

A_initial <- 0.8 * diag(m) + rep(0.2 / m, m)

#mean
emission <- list()
mean <- matrix(NA, nrow = m, ncol = d)
for (i in c(1:d)){
  mean[,i] <- quantile(quantile(as.vector(y[,i]), c(0.0001,0.9999), na.rm = FALSE,names = F, type = 7), probs = seq(0, 1, 1/(m+1)), na.rm = FALSE,names = F, type = 7)[2:(m+1)]
}

#sigma
sigma <- var(y)/m^2

for (i in c(1:m)) {
  emission[[i]] <- list(mean = mean[i,], sigma = sigma)
}

Protein_MLE <- EM_MLE(Protein_descriptor, 
                      A = A_initial,
                      parameter = emission,
                      dist_class = "mvnorm", 
                      m=27)

#### MPLE ####
# !!! This code might take dozens of hours to complete. !!! 
Protein_MPLE <- EM_MPLE(Protein_descriptor,
                        A = Protein_MLE$transition,
                        dist_class = "mvnorm",
                        parameter = Protein_MLE$par,
                        m=27)

### If parallel computing is available, 
### uncomment the line below and run it to compute the results for each lambda in parallel.
# library(parallel) # Use multicore.
# core_num <- 60 # Set number of cores.
# 
# run_one_lambda <- function(lam) {
#     tmp <- EM_MPLE_singleLambda(
#       y = Protein_descriptor,
#       A = Protein_MLE$transition,
#       dist_class = "mvnorm",
#       parameter = Protein_MLE$par,
#       Lambda = lam,
#       m = 27,
#       penalty_type = "group_fused_scad"
#     )
#     tmp
# }
# res_list <- mclapply(2^seq(-7.8, -2, 0.1), run_one_lambda, mc.cores = core_num)
# Protein_MPLE <- res_list[[
#   which.min(
#     unlist(lapply(res_list, `[[`, 'BIC'))
#     )
#   ]]

#### FHMM ####
Protein_FHMM <- fhmm(Protein_descriptor, M=3, K=3,parameter = emission, tol=1e-5)


save(Protein_MLE, Protein_MPLE, Protein_FHMM, file = "results/Protein_structure_analysis.Rdata")


