
rm(list=ls())

# load the data and the functions
source("/home/lc3521/supMDS/NM_method_reg_multi.R")
load("/home/lc3521/datafile/all_seqs_17countires_common.RData")
load("/home/piaac_data/background/BQ_17countries/bg_vars_demo_occu_scores.RData")
load("/home/lc3521/datafile/multi_targets.RData") 

args <- commandArgs(trailingOnly = T)
ncores <- as.numeric(args[1])
item_code <- c("U01a", "U01b", "U02", "U03a", "U04a", "U06a", "U06b", "U07", "U11b", "U16", "U19a", "U19b", "U23")
chosen_distance <- c("f1", "g1", "LCS", "lendiff")
n <- length(idx_match)  # total sample size
subset_size <- 200
K <- 10
lambda <- 1e-1
tol <- 1e-4
lambda2 <- 1

set.seed(123)
init_obj <- sample(1:n, subset_size)

# calculate theta's
theta_SingleMatrix_list <- list()
for(i_item in 1:length(item_code)){
  theta_SingleMatrix_list[[i_item]] <- list()
  cat('i_item=', i_item, "\n")
  item <- item_code[i_item]
  
  for(i_distance in 1:4){
    D <- get(load(paste0("/home/lc3521/dist_matrices/", item, "_D_", chosen_distance[i_distance], "_17countries_common.RData")))[idx_match, idx_match]
    theta_SingleMatrix_list[[i_item]][[i_distance]] <- mds_large(D/median(D), K, init_obj, ncores)
    rm(D)
  }
  D_f <- get(load(paste0("/home/lc3521/dist_matrices/",item,"_D_f1_17countries_common.RData")))[idx_match, idx_match]
  D_g <- get(load(paste0("/home/lc3521/dist_matrices/",item,"_D_g1_17countries_common.RData")))[idx_match, idx_match]
  D_oss <- D_f+D_g
  theta_SingleMatrix_list[[i_item]][[5]] <- mds_large(D_oss/median(D_oss), K, init_obj, ncores)
  D_oss <- D_f/median(D_f) + D_g/median(D_g)
  theta_SingleMatrix_list[[i_item]][[6]] <- mds_large(D_oss/median(D_oss), K, init_obj, ncores)
  
  names(theta_SingleMatrix_list[[i_item]]) <- c(chosen_distance, "OSS", "mOSS")
  rm(D_f, D_g, D_oss)
}

save(theta_SingleMatrix_list, file=paste0("/home/lc3521/supMDS/real_data/data/Theta_SingleMatrix_K", K, ".RData"))




