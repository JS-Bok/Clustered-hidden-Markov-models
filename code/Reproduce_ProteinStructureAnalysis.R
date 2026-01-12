################################################################################
#### Code for reproduce the numerical results of Protein structure analysis ####
################################################################################
#### Preliminary ####
# Load functions
sapply(c("code/EM_MLE.R","code/EM_MPLE.R","code/EM_FHMM.R","code/EM_Oracle.R"), source)

# Required library to load data 
library(parallel) # Use multicore.
core_num <- 25 # Set number of cores.

# Load data
load("data/Protein_data.Rdata")
Protein_descriptor <- lapply(Protein_data, function(x) x$descriptor_vec)

#### MLE ####
# This code might take several hours to complete.
### Initialization ###
y <- Reduce(rbind,y)
d <- ncol(y)

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

phi_initial = list(transition=trans,emission=emission)

Protein_MLE <- EM_MLE(Protein_descriptor, 
                      A = A_initial,
                      parameter = phi_initial,
                      dist_class = "mvnorm", 
                      m=27)

#### MPLE ####
# This code might take several hours to complete. 
Protein_MPLE <- EM_MPLE(Protein_descriptor,
                        A = Protein_MLE$transition,
                        dist_class = "mvnorm",
                        parameter = Protein_MLE$par,
                        m=27)

#### FHMM ####
Protein_FHMM <- fhmm(Protein_descriptor, M=3, K=3, tol=1e-5)


save(Protein_MLE, Protein_MPLE, Protein_FHMM, file = "results/Protein_structure_analysis.Rdata")


