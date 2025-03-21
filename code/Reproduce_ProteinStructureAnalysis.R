################################################################################
#### Code for reproduce the numerical results of Protein structure analysis ####
################################################################################
#### Preliminary ####
# Load functions
sapply(c("code/EM_MLE.R","code/EM_CHMM.R","code/EM_FHMM.R","code/EM_Oracle.R"), source)

# Required library to load data 
library(httr)
library(parallel) # Use multicore.
core_num <- 25 # Set number of cores.

# Base directory of data
home_dir <- "https://raw.githubusercontent.com/JS-Bok/Clustered-hidden-Markov-models/main/data/"

# Raw URL of the file
file_url <- paste0(home_dir,"Protein_data.Rdata")

# Load data
temp_file <- tempfile(fileext = ".Rdata")
response <- GET(file_url)
writeBin(content(response, as = "raw"), temp_file)
load(temp_file)
rm(response)
unlink(temp_file)
Protein_descriptor <- lapply(Protein_data, function(x) x$descriptor_vec)

#### MLE ####
# This code might take several hours to complete.
Protein_MLE <- EM_MLE(Protein_descriptor, dist_class = "mvnorm", m=27)

#### MPLE ####
# This code might take several hours to complete. 
Protein_MPLE <- EM_MPLE(Protein_descriptor,
                        A = Protein_MLE$transition,
                        dist_class = "mvnorm",
                        parameter = Protein_MLE$par,
                        m=27)

#### FHMM ####
Protein_FHMM <- fhmm(Protein_descriptor, M=3, K=3, tol=1e-5)
