library(glmnet)

source("/home/lc3521/supMDS/NM_fcns.R")
load("/home/lc3521/datafile/multi_targets.RData") 

# parameters
args <- commandArgs(trailingOnly = T)
ncores <- 4
subset_size <- 200
K <- 10

item_code <- c("U01a", "U01b", "U02", "U03a", "U04a", "U06a", "U06b", "U07", "U11b", "U16", "U19a", "U19b", "U23")
chosen_distance <- c("f1", "g1", "LCS", "lendiff")
n <- length(idx_match)  # total sample size

set.seed(123)
data_split <- sample(1:5, n, replace = T)
init_obj <- sample(1:n, subset_size)

#### Single Items ####
load(paste0("/home/lc3521/supMDS/real_data/data/Theta_SingleMatrix_K10.RData"))
load("/home/lc3521/supMDS/real_data/multi_4d_cv_res_0809.RData")

res_cv <- array(NA, dim=c(5, 13, length(y_list), 2, 2, 7), 
                dimnames=list(i_cv=1:5, n_item=1:13, target=1:length(y_list), 
                              IS_OS=c("IS", "OS"), reg_method=c("lm", "ridge"),
                              method=c("supMDS", "CAD", "UAD", "LCS", "SLD", "OSS", "mOSS")))
beta_cv <- array(NA, dim=c(5, 13, 3), 
                 dimnames=list(i_cv=1:5, n_item =1:13, dim=1:3))

for(i_cv in 1:5){
  idx_test <- which(data_split == i_cv)
  idx_train <- setdiff(1:n, idx_test)
  y_list_train <- lapply(y_list, function(y) as.vector(scale(y[idx_train])))
  y_list_test <- lapply(y_list, function(y) as.vector(scale(y[idx_test])))
  for(i_item in 1:13){
    beta_cv[i_cv, i_item, ] <- res_save[[i_cv]][[i_item]]$beta_est
    
    theta <- res_save[[i_cv]][[i_item]]$theta_est
    pred_res <- pred(theta, y_list, idx_train, idx_test)
    res_cv[i_cv, i_item, , "OS", "ridge", 1] <- pred_res$r2_OS
    res_cv[i_cv, i_item, , "OS", "lm", 1] <- pred_res$r2_lm_OS
    res_cv[i_cv, i_item, , "IS", "ridge", 1] <- pred_res$r2_IS
    res_cv[i_cv, i_item, , "IS", "lm", 1] <- pred_res$r2_lm_IS
    
    for(i_method in 1:6){
      theta <- theta_SingleMatrix_list[[i_item]][[i_method]]
      pred_res <- pred(theta, y_list, idx_train, idx_test)
      res_cv[i_cv, i_item, , "OS", "ridge", i_method+1] <- pred_res$r2_OS
      res_cv[i_cv, i_item, , "OS", "lm", i_method+1] <-  pred_res$r2_lm_OS
      res_cv[i_cv, i_item, , "IS", "ridge", i_method+1] <- pred_res$r2_IS
      res_cv[i_cv, i_item, , "IS", "lm", i_method+1] <- pred_res$r2_lm_IS
    }
    cat("i_cv=", i_cv, "i_item=", i_item, "\n")
  }
}
save(res_cv, file="/home/lc3521/supMDS/real_data/multi_4d_SingleItems_cv_res.RData")
save(beta_cv, file="/home/lc3521/supMDS/real_data/multi_4d_SingleItems_cv_beta_res.RData")

#### All Items ####
load(paste0("/home/lc3521/supMDS/real_data/data/Theta_SingleMatrix_K10.RData"))
load("/home/lc3521/supMDS/real_data/multi_4d_cv_res.RData")

res_data <- array(NA, dim=c(5, length(y_list), 2, 7), 
                  dimnames=list(i_cv=1:5, Target=1:length(y_list),
                                IS_OS=c("IS", "OS"),
                                Method=c("supMDS", "CAD", "UAD", 'LCS', 'SLD', "OSS", "mOSS")))

for(i_cv in 1:5){
  n <- length(data_split)
  idx_test <- which(data_split == i_cv)
  idx_train <- setdiff(1:n, idx_test)
  
  theta_opt = theta_f1 = theta_g1 = theta_LCS = theta_lendif = theta_oss = theta_moss <- matrix(NA, n, 0)
  for(i_item in 1:13){
    theta_opt <- cbind(theta_opt, res_save[[i_cv]][[i_item]]$theta_est)
    theta_f1 <- cbind(theta_f1, theta_SingleMatrix_list[[i_item]]$f1)
    theta_g1 <- cbind(theta_g1, theta_SingleMatrix_list[[i_item]]$g1)
    theta_LCS <- cbind(theta_LCS, theta_SingleMatrix_list[[i_item]]$LCS)
    theta_lendif <- cbind(theta_lendif, theta_SingleMatrix_list[[i_item]]$lendiff)
    theta_oss <- cbind(theta_oss, theta_SingleMatrix_list[[i_item]]$OSS)
    theta_moss <- cbind(theta_moss, theta_SingleMatrix_list[[i_item]]$mOSS)
  }
  pred1 <- pred(theta_opt, y_list, idx_train, idx_test)
  pred2 <- pred(theta_f1, y_list, idx_train, idx_test)
  pred3 <- pred(theta_g1, y_list, idx_train, idx_test)
  pred4 <- pred(theta_LCS, y_list, idx_train, idx_test)
  pred5 <- pred(theta_lendif, y_list, idx_train, idx_test)
  pred6 <- pred(theta_oss, y_list, idx_train, idx_test)
  pred7 <- pred(theta_moss, y_list, idx_train, idx_test)
  res_data[i_cv, , "IS", 1] <- pred1$r2_IS
  res_data[i_cv, , "OS", 1] <- pred1$r2_OS
  res_data[i_cv, , "IS", 2] <- pred2$r2_IS
  res_data[i_cv, , "OS", 2] <- pred2$r2_OS
  res_data[i_cv, , "IS", 3] <- pred3$r2_IS
  res_data[i_cv, , "OS", 3] <- pred3$r2_OS
  res_data[i_cv, , "IS", 4] <- pred4$r2_IS
  res_data[i_cv, , "OS", 4] <- pred4$r2_OS
  res_data[i_cv, , "IS", 5] <- pred5$r2_IS
  res_data[i_cv, , "OS", 5] <- pred5$r2_OS
  res_data[i_cv, , "IS", 6] <- pred6$r2_IS
  res_data[i_cv, , "OS", 6] <- pred6$r2_OS
  res_data[i_cv, , "IS", 7] <- pred7$r2_IS
  res_data[i_cv, , "OS", 7] <- pred7$r2_OS
}
save(res_data, file="/home/lc3521/supMDS/real_data/multi_4d_AllItems_cv_res.RData")


#### plot ####
# single items
library(ggplot2)
library(latex2exp)

res_cv <- get(load("/Users/lingchen/Desktop/Supervised MDS/SupMDS codes/data/realdata/multi_4d_SingleItems_cv_res.RData"))
res_cv <- apply(res_cv, 2:6, function(x) mean(x, na.rm=T)) # average over cv's
res_cv <- apply(res_cv, c(2:5), function(x) mean(x, na.rm=T))  # average over items
res_cv_mat <- as.data.frame.table(res_cv[, "OS", "ridge", 1:5])

ggplot(data=res_cv_mat, aes(x=method, y=Freq)) + 
  geom_boxplot(outlier.size=.5) +
  ylab(TeX("$OSR^2$")) + xlab("Method") + 
  theme(axis.text.x = element_text(angle = 45, vjust = 0.7), text = element_text(size = 9)) + 
  ggtitle("Single Items. Avg over items.")
ggsave("/Users/lingchen/Desktop/Supervised MDS/SupMDS codes/images/real data/multi_4d_SingleItems.pdf", 
       width = 4, height = 4)

beta_data <- get(load("/Users/lingchen/Desktop/Supervised MDS/SupMDS codes/data/realdata/multi_4d_SingleItems_cv_beta_res.RData"))
#beta_data <- apply(beta_data, 2:4, function(x) mean(x, na.rm=T))
dimnames(beta_data)$dim <- c("UAD", "LCS", "SLD")
beta_data_mat <- as.data.frame.table(beta_data)

ggplot(beta_data_mat, aes(x=dim, y=Freq)) + 
  geom_boxplot(outlier.size=.5) + 
  ylab(TeX("Estimated $\\omega$")) + 
  xlab('Method') + 
  theme(text = element_text(size = 8)) # + ggtitle(TeX("Estimated $\\omega$ with $D_0=d_{CAD},\ D_1=d_{UAD},\ D_2=d_{LCS},\ D_3=d_{SLD}$"))
ggsave("/Users/lingchen/Desktop/Supervised MDS/writeup/images/multi_4d_BetaEst.pdf", 
       width = 4, height = 4)

# all items
load("/Users/lingchen/Desktop/Supervised MDS/SupMDS codes/data/realdata/multi_4d_AllItems_cv_res.RData")

res_data <- apply(res_data, 2:4, function(x) mean(x, na.rm=T))
res_data_mat <- as.data.frame.table(res_data[,"OS", c(1, 6, 2:5)])
ggplot(data=res_data_mat, aes(x=method, y=Freq))  + 
  geom_boxplot(outlier.size=.5) +
  ylab(TeX("$OSR^2$")) + xlab("Method") + 
  theme(axis.text.x = element_text(vjust = 0.7), text = element_text(size = 9)) #+ 
#ggtitle("All Items")
ggsave("/Users/lingchen/Desktop/Supervised MDS/writeup/images/multi_4d_AllItems.pdf", 
       width = 4, height = 4)

res_data_mat <- as.data.frame.table(res_data[,"OS", c(1, 6)])
ggplot(data=res_data_mat, aes(x=method, y=Freq))  + 
  geom_boxplot(outlier.size=.5) +
  ylab(TeX("$OSR^2$")) + xlab("Method") + 
  theme(axis.text.x = element_text(vjust = 0.7), text = element_text(size = 9)) #+ 
#ggtitle("All Items")
ggsave("/Users/lingchen/Desktop/Supervised MDS/writeup/images/multi_4d_AllItems_OSS.pdf", 
       width = 4, height = 4)



