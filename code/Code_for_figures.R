###################################################################
#### Code for generating figures in protein structure analysis ####
###################################################################
# Required libraries
library('tidyverse')
library('gridExtra')
library('magick')
library('ggplot2')
library('grid')
library('parallel') # Use multicore.
core_num <- 25 # Set number of cores.

# Load data
source("code/Auxiliary_Functions.R")
load("data/Protein_data.Rdata")
load("results/Protein_structure_analysis.Rdata")

#### Reconstruct hidden sequences and match with secondary structure ####
# Data load
SS_seq <- lapply(Protein_data,function(x){x$secondary_structure})

y <- lapply(Protein_data,function(x){x$descriptor_vec})

A <- Clust_Ordering(Protein_MPLE$transition)$transition # This is the 27 by 27 transition matrix in Figure S1.

ind <- cf(A)$index
par <- Protein_MPLE$par[Clust_Ordering(Protein_MPLE$transition)$order]
dist_class <- "mvnorm"

### Define Viterbi Algorithm for reconstruction of hidden sequences ###
Viterbi <- function(y, A, em_mat){
  n <- nrow(y)
  m <- nrow(A)
  x <- stationary_dist(A)
  xi <- matrix(0, n, m)
  foo <- x*em_mat[,1]
  xi[1,] <- foo/sum(foo)
  for (i in 2:n) {
    foo <- apply(xi[i-1,]*A,2,max)*em_mat[,i]
    xi[i,] <- foo/sum(foo)
  }
  iv <- numeric(n)
  iv[n] <- which.max(xi[n,])
  for (i in (n-1):1) {
    iv[i] <- which.max(A[,iv[i+1]]*xi[i,])
  }
  return(iv)
}

RC_seq <- function(y, A, par, dist_class, ind){
  x <- stationary_dist(A)
  em_mat_list <- lapply(y, em_prob_mat, x=x, dist_class=dist_class, parameter=par)
  Vit_tmp <- function(y,em_mat){return(Viterbi(y,A,em_mat))}
  rc_s <- mapply(Vit_tmp, y=y, em_mat=em_mat_list)
  for (i in 1:length(ind)) {
    if(i==1){
      ind_mat<- as.data.frame(cbind(i,ind[[i]]))
      next
    }
    ind_mat<- rbind(ind_mat,cbind(i,ind[[i]]))
  }
  rc_s_clust <- rc_s
  for (i in 1:nrow(ind_mat)) {
    tmp2 <- ind_mat[i,2]
    tmp1 <- ind_mat[i,1]
    rc_s_clust <- sapply(rc_s_clust, function(x,tmp1,tmp2) {x[x==tmp2] <- tmp1; return(x)},tmp1=tmp1,tmp2=tmp2)
  }
  return(list(rc_s=rc_s, rc_s_clust=rc_s_clust))
}

## Reconstruction of hidden sequences in CHMM ##
RC_seq_tmp <- RC_seq(y, A=A, par=par, dist_class=dist_class, ind=ind)
RC_seq_tmp2 <- mapply(function(original,clustered){data.frame(original=original,clustered=clustered)}, 
                      original=RC_seq_tmp$rc_s, clustered=RC_seq_tmp$rc_s_clust, 
                      SIMPLIFY = FALSE)

## Reconstruction of hidden sequences in FHMM ##
# Mu = Protein_FHMM$Mu
# Cov = Protein_FHMM$Cov
P = Protein_FHMM$P
Pi = Protein_FHMM$Pi
M=3
K=3

A_fhmm <- P[(K*(M-1)+1):(K*M),]
for (i in (M-1):1) {
  A_fhmm <- kronecker(P[(K*(i-1)+1):(K*i),], A_fhmm)
}

parameter <- Protein_FHMM$par

em_mat_list <- lapply(y, em_prob_mat, x=Pi_hmm, dist_class=dist_class, parameter=parameter)
Vit_tmp <- function(y,em_mat){return(Viterbi(y,A_fhmm,em_mat))}
fhmm_rc <- mapply(Vit_tmp, y=y, em_mat=em_mat_list)


SS_match <- mapply(function(x,y,z){data.frame(original=x$original, clustered=x$clustered, fhmm=y,SS=z)},
                   x=RC_seq_tmp2, y=fhmm_rc, z=SS_seq, SIMPLIFY = FALSE)



#### Code generating Figure S1 ####
figS1 <- round(A,3)*100
colnames(figS1) <- 1:27
write.csv(figS1, "results/Figure_S1.csv")



#### Code generating Figure 3 ####
SS_match_combined <- Reduce(rbind, SS_match) 

SS_match_combined[SS_match_combined=="H"] <- "Alpha Helix"
SS_match_combined[SS_match_combined=="I"] <- "Pi Helix"
SS_match_combined[SS_match_combined=="G"] <- "3-10 Helix"
SS_match_combined[SS_match_combined=="B"] <- "Isolated beta-bridge"
SS_match_combined[SS_match_combined=="E"] <- "Extended strand"
SS_match_combined[SS_match_combined=="T"] <- "Hydrogen bonded turn"
SS_match_combined[SS_match_combined=="S"] <- "Bend"
SS_match_combined[SS_match_combined=="-"] <- "Other or Coil"
SS_match_combined[SS_match_combined=="n"] <- "Other or Coil"


SS_match_combined$original <- as.factor(SS_match_combined$original)



SS_match_combined$clustered <- as.factor(SS_match_combined$clustered)
SS_match_combined$clustered <- sapply(SS_match_combined$clustered, function(x){return(paste("C",x,sep =""))})
SS_match_combined$clustered <- as.factor(SS_match_combined$clustered)
SS_match_combined$clustered <- factor(SS_match_combined$clustered, 
                             levels=c("C1",  "C2", "C3", "C4", "C5", "C6", "C7", "C8", "C9", "C10", 
                                      "C11",  "C12",  "C13",  "C14",  "C15",  "C16",  "C17",  "C18"))



max_y <- 5250


stacked_plot_original <- ggplot(SS_match_combined, aes(x = original, fill = SS)) +
  geom_bar(position = "stack") +
  scale_y_continuous(limits = c(0, max_y)) +
  theme_minimal() +
  theme(legend.position = "none", 
        axis.title.x = element_blank(), 
        axis.title.y = element_blank())

ggsave("results/Figure_3-original.pdf", stacked_plot_original, width = 7, height = 6, units = "in", dpi = 600)


stacked_plot_clustered <- ggplot(SS_match_combined, aes(x = clustered, fill = SS)) +
  geom_bar(position = "stack") +
  scale_y_continuous(limits = c(0, max_y)) +
  theme_minimal() +
  theme(axis.title.x = element_blank(), 
        axis.title.y = element_blank())+
  labs(fill = "Secondary structure") 

ggsave("results/Figure_3-clustered.pdf", stacked_plot_clustered, width = 7, height = 6, units = "in", dpi = 600)


#### Code generating Figure S2 ####

K=3; M=3
index_mat <- sapply(1:M, function(x){rep(1:K, K^(x-1), each=K^(M-x))})

fhmm_x1 <- SS_match_combined$fhmm
fhmm_x1[fhmm_x1 %in% which(index_mat[,1]==1)] <- 1
fhmm_x1[fhmm_x1 %in% which(index_mat[,1]==2)] <- 2
fhmm_x1[fhmm_x1 %in% which(index_mat[,1]==3)] <- 3

SS_match_combined$fhmm_x1 <- fhmm_x1

fhmm_x2 <- SS_match_combined$fhmm
fhmm_x2[fhmm_x2 %in% which(index_mat[,2]==1)] <- 1
fhmm_x2[fhmm_x2 %in% which(index_mat[,2]==2)] <- 2
fhmm_x2[fhmm_x2 %in% which(index_mat[,2]==3)] <- 3

SS_match_combined$fhmm_x2 <- fhmm_x2


fhmm_x3 <- SS_match_combined$fhmm
fhmm_x3[fhmm_x3 %in% which(index_mat[,3]==1)] <- 1
fhmm_x3[fhmm_x3 %in% which(index_mat[,3]==2)] <- 2
fhmm_x3[fhmm_x3 %in% which(index_mat[,3]==3)] <- 3


SS_match_combined$fhmm_x3 <- fhmm_x3



max_y <- 15000


stacked_plot_fhmm_x1 <- ggplot(SS_match_combined, aes(x = fhmm_x1, fill = SS)) +
  geom_bar(position = "stack") +
  scale_y_continuous(limits = c(0, max_y)) +
  theme_minimal() +
  theme(legend.position = "none", 
        axis.title.x = element_blank(), 
        axis.title.y = element_blank())

ggsave("results/Figure_S2-x1.pdf", stacked_plot_fhmm_x1, width = 1.5, height = 6, units = "in", dpi = 600)


stacked_plot_fhmm_x2 <- ggplot(SS_match_combined, aes(x = fhmm_x2, fill = SS)) +
  geom_bar(position = "stack") +
  scale_y_continuous(limits = c(0, max_y)) +
  theme_minimal() +
  theme(legend.position = "none", 
        axis.title.x = element_blank(), 
        axis.title.y = element_blank())

ggsave("results/Figure_S2-x2.pdf", stacked_plot_fhmm_x2, width = 1.5, height = 6, units = "in", dpi = 600)

stacked_plot_fhmm_x3 <- ggplot(SS_match_combined, aes(x = fhmm_x3, fill = SS)) +
  geom_bar(position = "stack") +
  scale_y_continuous(limits = c(0, max_y)) +
  theme_minimal() +
  theme(axis.title.x = element_blank(), 
        axis.title.y = element_blank())+
  labs(fill = "Secondary structure") 

ggsave("results/Figure_S2-x3.pdf", stacked_plot_fhmm_x3, width = 3.36, height = 6, units = "in", dpi = 600)

#### Code generating Figure 4 ####
library('tibble')
library('NGLVieweR')


# Load 2PNE data
pne <- SS_match[[136]]


## Clusted

c_1 <- which(pne$clustered==1)+1
c_1_1 <- c_1[c_1<59]
c_1_2 <- c_1[c_1>=59]

c_1_1 <- sapply(c_1_1, function(x){paste(x,x+1,sep = "-")})
c_1_1 <- paste(c_1_1, collapse = " or ")
c_1_2 <- sapply(c_1_2, function(x){paste(x,x+1,sep = "-")})
c_1_2 <- paste(c_1_2, collapse = " or ")

NGLVieweR("2PNE") %>%
  stageParameters(backgroundColor = "white") %>%
  addRepresentation("backbone", param = list(
    sele = "1-81",
    colorScheme = "uniform",
    radiusScale = 0.5  
  )) %>%
  addRepresentation("backbone", param = list(
    sele = c_1_1,
    colorValue = "red",
    colorScheme = "uniform",
    radiusScale = 2.0  
  )) %>%
  addRepresentation("backbone", param = list(
    sele = c_1_2,
    colorValue = "red",
    colorScheme = "uniform",
    radiusScale = 2.0  
  )) 



## Raw

s_1 <- which(pne$original==1)+1
s_1 <- sapply(s_1, function(x){paste(x,x+1,sep = "-")})
s_3 <- which(pne$original==3)+1
s_3 <- sapply(s_3, function(x){paste(x,x+1,sep = "-")})
s_4 <- which(pne$original==4)+1
s_4 <- sapply(s_4, function(x){paste(x,x+1,sep = "-")})
s_1 <- paste(s_1, collapse = " or ")
s_3 <- paste(s_3, collapse = " or ")
s_4 <- paste(s_4, collapse = " or ")

NGLVieweR("2PNE") %>%
  stageParameters(backgroundColor = "white") %>%
  addRepresentation("backbone", param = list(
    sele = "1-81",
    colorScheme = "uniform",
    radiusScale = 0.5
  )) %>%
  addRepresentation("backbone", param = list(
    sele = s_1,
    colorValue = "purple", 
    radiusScale = 2.0
  )) %>%
  addRepresentation("backbone", param = list(
    sele = s_3,
    colorValue = "lightblue", 
    radiusScale = 2.0
  )) %>%
  addRepresentation("backbone", param = list(
    sele = s_4,
    colorValue = "#ff7f0e", # 연한 초록 (밝은 명암)
    radiusScale = 2.0
  ))

## Secondary Structure
oc <- which(pne$SS=="-")+1
oc <- sapply(oc, function(x){paste(x,x+1,sep = "-")})
oc <- paste(oc, collapse = " or ")

NGLVieweR("2PNE") %>%
  stageParameters(backgroundColor = "white") %>%
  addRepresentation("backbone",param = list(
    sele = "1-81",
    colorScheme = "uniform",
    radiusScale = 0.5
  ) )%>%
  addRepresentation("backbone",
                    param = list(
                      sele = oc,
                      colorValue = "#2ca02c",
                      colorScheme = "uniform",
                      radiusScale = 2.0
                    ) 
  )

###------------------------------------------------------------------###
###* We saved the four figures in result folder with size 350 X 550 *###
###------------------------------------------------------------------###
library(magick)  
library(ggplot2)
library(grid)
library(patchwork)  

# image path
clustered_path <- "results/Figure_4-clustered.png"
original_path <- "results/Figure_4-original.png"
ss_path <- "results/Figure_4-SS.png"

crop_image <- function(image_path, crop_width = 50, crop_height = 75) {
  img <- image_read(image_path)
  info <- image_info(img)
  new_width <- info$width - (2 * crop_width)  # 양쪽에서 crop_width 만큼씩 자름
  new_height <- info$height - crop_height
  image_crop(img, geometry = paste0(new_width, "x", new_height, "+", crop_width, "+0"))
}

ss_img <- crop_image(ss_path, crop_width = 10, crop_height = 50)
original_img <- crop_image(original_path, crop_width = 10, crop_height = 50)
clustered_img <- crop_image(clustered_path, crop_width = 10, crop_height = 50)

scale_height <- function(image, scale_factor = 1.1) {
  info <- image_info(image)
  new_height <- as.integer(info$height * scale_factor)  # 높이만 scale_factor 배로 조정
  image_resize(image, paste0(info$width, "x", new_height, "!"))  # 강제로 높이만 조정
}

# original_img <- scale_height(original_img, scale_factor = 1.05)  # 10% 증가

#--- 2. raster 이미지를 ggplot 객체로 변환 ----------------------------------------
ss_plot <- ggplot() +
  annotation_custom(rasterGrob(as.raster(ss_img), interpolate = TRUE),
                    xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf) +
  theme_void()

original_plot <- ggplot() +
  annotation_custom(rasterGrob(as.raster(original_img), interpolate = TRUE),
                    xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf) +
  theme_void()

clustered_plot <- ggplot() +
  annotation_custom(rasterGrob(as.raster(clustered_img), interpolate = TRUE),
                    xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf) +
  theme_void()

#--- 3. names -------------------------------------------------
ss_plot        <- ss_plot        + ggtitle("(a) Secondary structure")
original_plot  <- original_plot  + ggtitle("(b) Original state")
clustered_plot <- clustered_plot + ggtitle("(c) Cluster")


title_theme <- theme(
  plot.title = element_text(hjust = .1),
  plot.margin = margin(t = 5, r = 5, b = 5, l=0)   # 틈 조금
)

ss_plot        <- ss_plot        + title_theme
original_plot  <- original_plot  + title_theme
clustered_plot <- clustered_plot + title_theme

#--- 4.legend -----------------------------------------------------------
legend_plot <- ggplot() +
  geom_point(aes(x = 0,   y = 1.15), colour = "#2ca02c", size = 5) +
  annotate("text", x = 0.1, y = 1.15, label = "Other or Coil", hjust = 0, size = 5) +
  geom_point(aes(x = 1.1, y = 1.15), colour = "purple", size = 5) +
  annotate("text", x = 1.2, y = 1.15, label = "State 1", hjust = 0, size = 5) +
  geom_point(aes(x = 1.1, y = 0.95), colour = "lightblue", size = 5) +
  annotate("text", x = 1.2, y = 0.95, label = "State 3", hjust = 0, size = 5) +
  geom_point(aes(x = 1.1, y = 0.75), colour = "#ff7f0e", size = 5) +
  annotate("text", x = 1.2, y = 0.75, label = "State 4", hjust = 0, size = 5) +
  geom_point(aes(x = 2.2, y = 1.15), colour = "red", size = 5) +
  annotate("text", x = 2.3, y = 1.15,
           label = expression(Cluster~C[1]), hjust = 0, size = 5) +
  xlim(0, 3) + ylim(0.5, 1.2) +
  theme_void()

#--- 5. patchwork -------------------------------------------------


row_plots <- (ss_plot | original_plot | clustered_plot) 

final_plot <- row_plots / legend_plot +
  plot_layout(heights = c(3, 1))

final_plot


ggsave("results/Figure_4.pdf", final_plot, width = 7, height = 5.7, dpi = 600)


