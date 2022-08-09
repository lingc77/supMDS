
#### process data ####
library(ProcData)

item_code <- c("U01a", "U01b", "U02", "U03a", "U04a", "U06a", "U06b", "U07", "U11b", "U16", "U19a", "U19b", "U23")
# no data for U21

#### seqs of the common IDs of the 17 countries ####
item_seqs <- list()
ID_list <- list()
for(i_item in 1:length(item_code)){
  item <- item_code[i_item]
  file_path <- paste0("/data/home/piaac_data/data_Britt_full/", item, "_data_cleaned.csv")
  item_seqs[[i_item]] <- read.seqs(file_path, style = "multiple", id_var = "ID", action_var = "Coding")
  ID_list[[i_item]] <- names(item_seqs[[i_item]]$action_seqs)
}
IDs_all <- Reduce(union, ID_list)
length(IDs_all) # 40230

IDs_common <- Reduce(intersect, ID_list) # IDs that answered all the items
length(IDs_common)

item_seqs_common <- list()
for(i_item in 1:length(item_code)){
  idx_common <- match(IDs_common, ID_list[[i_item]])
  item_seqs_common[[i_item]] <- sub_seqs(item_seqs[[i_item]], idx_common)
}
library(stringr)
IDs_common <- sapply(IDs_common, function(x) str_replace(x, pattern = '_', replacement = ''))
save(IDs_common, item_code, item_seqs_common, file="/home/lc3521d/datafile/all_seqs_17countires_common.RData")

# summary data
#summ <- list()
#for(i_item in 1:13){
#  summ[[i_item]] <- summary(item_seqs[[i_item]])
#}
item_seqs_noskip <- list()
for(i_item in 1:13){
  len <- sapply(item_seqs[[i_item]]$action_seqs, length)
  idx_skip <- which(len == 3)
  item_seqs_noskip[[i_item]] <- sub_seqs(item_seqs[[i_item]], 
                                         (1:length(item_seqs[[i_item]]$action_seqs))[-idx_skip])
  item_seqs_noskip[[i_item]] <- remove_repeat(item_seqs_noskip[[i_item]])
}

ID_list_noskip <- list()
summ_noskip <- list()
for(i_item in 1:13){
  ID_list_noskip[[i_item]] <- names(item_seqs_noskip[[i_item]]$action_seqs)
  #summ_noskip[[i_item]] <- summary(item_seqs_noskip[[i_item]])
}
IDs_noskip_common <- Reduce(intersect, ID_list_noskip) 
length(IDs_noskip_common)

#### distance matrices ####
library(parallel)
library(qualV)
load("/home/lc3521/datafile/all_seqs_17countires_common.RData")
source("/home/lc3521/supMDS/functions_dist.R")

#methods <- c("f1", "g1", "g2", "LCS", "lendiff", "edit")
methods <- c("edit")
for(method in methods){
  for(i_item in 1:13){
    item <- item_code[i_item]
    cat(method, item, "\n")
    t1 <- Sys.time()
    D_res <- calculate_dist_mat(item_seqs_common[[i_item]]$action_seqs, method=method, ncores=4)
    t2 <- Sys.time(); cat("D cal time="); print(t2-t1)
    
    save(D_res, file = paste0("/home/lc3521/dist_matrices/", item, "_D_", method, "_17_countries_common.RData"))
  }
}
# modify f1, g1 and g2.
methods <- c("f1", "g1", "g2")
for(method in methods){
  for(i_item in 1:13){
    item <- item_code[i_item]
    cat(method, item, "\n")
    
    D_res0 <- get(load(paste0("/home/lc3521/dist_matrices/", item, "_D_", method, "_17_countries_common.RData")))
    n <- nrow(D_res0)
    D_res <- matrix(0, n, n)
    t1 <- Sys.time()
    for(i in 1:(n-1)){
      for(j in (i+1):n){
        D_res[i,j] <- D_res0[i, j] / (length(item_seqs_common[[i_item]]$action_seqs[[i]]) + length(item_seqs_common[[i_item]]$action_seqs[[j]]))
      }
    }
    D_res <- D_res + t(D_res)
    t2 <- Sys.time(); cat("D cal time="); print(t2-t1)
    
    save(D_res, file = paste0("/home/lc3521/dist_matrices/", item, "_D_", method, "_17countries_common.RData"))
  }
}

#### the four targets for 17ctys ####
load("/home/piaac_data/background/BQ_17countries/bg_vars_demo_occu_scores.RData")
load("/home/lc3521/datafile/common_IDs_noskip_17ctys.RData")

IDs <- all_bgs_skills_scores[, 1]
idx_match <- match(common_IDs_noskip, IDs)
age <- all_bgs_skills_scores[idx_match, "age"]
lit <- all_bgs_skills_scores[idx_match, "Literacy"]
num <- all_bgs_skills_scores[idx_match, "Numeracy"]
PSTRE <- all_bgs_skills_scores[idx_match, "PSTRE"]

targets <- data.frame(common_IDs_noskip = common_IDs_noskip, age=age, lit=lit, num=num, PSTRE=PSTRE)

save(targets, file="/home/lc3521/datafile/targets.RData")


#### target summaries #####
load("/home/lc3521/multi_targets/data/bg_vars_17countries_common.RData") 
load("/home/lc3521/datafile/common_IDs_noskip_17ctys.RData") 

chosen_targets <- chosen_targets[-(1:4)] # only consider the other targets
idx_match <- match(common_IDs_noskip, IDs_common)  # IDs_common are IDs that answered all items

y_list <- list()
for(i_target in 1:length(chosen_targets)){
  y_list[[i_target]] <- bgs_targets_common[idx_match, chosen_targets[i_target]]  # the chosen targets
}
idx_na <- lapply(y_list, function(x) which(is.na(x)))
idx_na <- Reduce(intersect, idx_na)  # NA indices intersection
idx_match <- idx_match[-idx_na]
y_list <- list()
for(i_target in 1:length(chosen_targets)){
  y_list[[i_target]] <- bgs_targets_common[idx_match, chosen_targets[i_target]]  # the chosen targets
}
save(idx_match, y_list, idx_match, file='/home/lc3521/datafile/multi_targets.RData')



