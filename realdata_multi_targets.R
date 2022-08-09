
rm(list=ls())
library(mirt)

# load the data and the functions
source("/home/lc3521/supMDS/NM_fcns.R")
load("/home/lc3521/multi_targets/data/bg_vars_17countries_common.RData")  # chosen targets
load("/home/lc3521/datafile/common_IDs_noskip_17ctys.RData")  # common IDs & no skip
chosen_targets <- chosen_targets[-(1:4)] # only consider the other targets

# parameters
args <- commandArgs(trailingOnly = T)
ncores <- as.numeric(args[1])
d_num <- as.numeric(args[2])
minit <- 3
maxit <- 50 
l_beta <- 1e-2
subset_size <- 200
K <- 10

# initialization
item_code <- c("U01a", "U01b", "U02", "U03a", "U04a", "U06a", "U06b", "U07", "U11b", "U16", "U19a", "U19b", "U23")
if (d_num == 4){
  chosen_distance <- c("f1", "g1", "LCS", "lendiff")
}
if (d_num == 2){
  chosen_distance <- c("f1", "g1")
}
idx_match <- match(common_IDs_noskip, IDs_common)  # IDs_common are IDs that answered all items

y_list <- list()
for(i_target in 1:length(chosen_targets)){
  y_list[[i_target]] <- bgs_targets_common[idx_match, chosen_targets[i_target]]  # the chosen targets
}
idx_na <- lapply(y_list, function(x) which(is.na(x)))
idx_na <- Reduce(intersect, idx_na)  # NA indices intersection
idx_match <- idx_match[-idx_na]
#print(length(idx_match))
y_list <- list()
for(i_target in 1:length(chosen_targets)){
  y_list[[i_target]] <- bgs_targets_common[idx_match, chosen_targets[i_target]]  # the chosen targets
}

n <- length(idx_match)  # total sample size
set.seed(123)
data_split <- sample(1:5, n, replace = T)

set.seed(123)
res_save <- list()
for(i_cv in 1:5){
  res_save[[i_cv]] <- list()
  idx_test <- which(data_split == i_cv)
  idx_train <- setdiff(1:n, idx_test)
  y_list_train <- lapply(y_list, function(y) as.vector(scale(y[idx_train])))
  y_list_test <- lapply(y_list, function(y) as.vector(scale(y[idx_test])))
  
  for(i_item in 1:length(item_code)){
    item <- item_code[i_item]
    D_list <- list()
    for(i in 1:length(chosen_distance)) { 
      D_list[[i]] <- get(load(paste0("/home/lc3521/dist_matrices/", item, "_D_", chosen_distance[i], "_17countries_common.RData")))
      D_list[[i]] <-  D_list[[i]][idx_match, idx_match] 
    }
    D_list <- lapply(D_list, function(x) x/median(x))
    
    # optimization
    if (d_num == 4){
      beta_simplex <- matrix(abs(c(rnorm(3, 0, 0.2), rnorm(1, 10, 1), rnorm(3, 0, 0.2), 
                                   rnorm(1, 10, 1), 
                                   rnorm(3, 0, 0.2), rnorm(1, 10, 1))), nr=3)
      beta_simplex <- split(t(beta_simplex), rep(1:ncol(beta_simplex)))
      tol <- 3e-1
    }
    if (d_num == 2) {
      beta_simplex <- list(abs(rnorm(1, 0, 0.2)), abs(rnorm(1, 10, 1)))
      tol <- 1e-1
    }
    init_obj <- sample(1:length(idx_train), subset_size)
    t1 <- Sys.time()
    res <- main_large_NM_multi(lapply(D_list, function(x) x[idx_train, idx_train]), y_list_train, K, 
                               beta_simplex, init_obj, l_beta, reg=T, 
                               tol, alpha=1, gamma=2, rho=1/2, sigma=1/2, 
                               minit, maxit, ncores, output=F)
    res$theta_full <- mds_large(Reduce("+", Map("*", D_list, c(1,res$beta_est))), K, init_obj, ncores)
    t2 <- Sys.time(); print(t2-t1)
    res_save[[i_cv]][[i_item]] <- res
    
    cat("i_cv=", i_cv, "item=", item, "iter=", res$iter, "beta est=", res$beta_est, "\n")
    
  }
  save(res_save, file=paste0("/home/lc3521/supMDS/real_data/multi_", d_num, "d_cv_res_0809.RData"))
}


