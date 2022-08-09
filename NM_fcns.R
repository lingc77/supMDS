
library(parallel)
library(glmnet)
library(MASS)
library('Rcpp')
library('inline')

rcpp_inc <- '
using namespace Rcpp;
using namespace arma;
'
src_matmult <- '
mat m1 = as<mat>(m1in);
mat m2 = as<mat>(m2in);
mat cp = m1 * m2;
return(wrap(cp));
'

#' \code{matmult_cpp} calculates the multiplication of two matrices using cpp.
#' \code{mds_large} extracts \code{K} features from a large distance matrix by
#' multidimensional scaling. The algorithm first selects a relatively small subset of 
#' objects to perform the classical MDS. Then the coordinate of each of the 
#' other objects are obtained by minimizing the loss function related to the target'
#' distance to those in the subset through BFGS.
#' \code{mds_l2_small} extracts \code{K} features from a full distance matrix by l2-regularized
#' multidimensional scaling.
#' \code{mds_l2_large} extracts \code{K} features from a large distance matrix by l2-regularized
#' multidimensional scaling. The algorithm first selects a relatively small subset of 
#' objects to perform the l2-MDS. Then the coordinate of each of the other objects 
#' are obtained by minimizing the loss function related to the target'
#' distance to those in the subset through BFGS.
#' \code{ylist_pred} gives in-sample prediction results for a list of y.
#' \code{ylist_pred_OS} gives extracted features and  in-sample and out-of-sample prediction 
#' results for a list of y, given the  training and testing indices.

matmult_cpp <- cxxfunction(signature(m1in="numeric", m2in="numeric"), src_matmult, 
                           plugin='RcppArmadillo', rcpp_inc)

mds_large <- function(D_full, K, init_obj, ncores, pca = FALSE) {
  #' D_full: distance matrix
  #' K: feature dimension
  #' init_obj: the indices of fully-connected objects
  #' ncores: number of cores to be used
  #' pca: logical. Whether to use PCA for the final features
  #' return: extracted features
  #' required package: parallel
  n <- nrow(D_full)
  if (K > n) {
    stop("Feature dimension larger than distance matrix dimension!\n")
  }
  if (length(init_obj) >= n){
    stop("Length of initial indices larger than distance matrix dimension!\n")
  }
  
  theta <- matrix(0, n, K)
  subset_size <- length(init_obj)
  remn_obj <- setdiff(1:n, init_obj)
  
  theta[init_obj, ] <- cmdscale(D_full[init_obj, init_obj], k = K)
  
  obj_fun <- function (theta_j, theta_m_mat, d_vec) {
    theta_j_mat <- matmult_cpp(cbind(rep(1, subset_size)), t(theta_j))
    theta_diff <- theta_m_mat - theta_j_mat
    res <- sum((d_vec - sqrt(rowSums((theta_diff)^2)))^2)
    res
  }
  grad_fun <- function (theta_j, theta_m_mat, d_vec) {
    theta_j_mat <- matmult_cpp(cbind(rep(1, subset_size)), t(theta_j))
    theta_diff <- theta_m_mat - theta_j_mat
    dhat_vec <- sqrt(rowSums((theta_diff)^2))
    res <- 2 * colSums((d_vec / dhat_vec - 1) * theta_diff)
    res
  }
  opt_fcn <- function (d_vec) {
    opt_res <- optim(rnorm(K), fn = obj_fun, gr = grad_fun, method = "BFGS", 
                     theta_m_mat = theta[init_obj, ], d_vec = d_vec)
    opt_res$par
  }
  
  # split into lists of columns
  dist_list <- split(D_full[init_obj, remn_obj], rep(1:length(remn_obj), each = subset_size))
  theta[remn_obj, ] <- t(mcmapply(opt_fcn, dist_list, mc.cores = ncores))
  
  if (pca) theta <- prcomp(theta, scale.=T)$x
  
  return(theta)
}

mds_l2_small <- function(D_full, K, lambda=1e-1, tol, pca = FALSE){
  #' D_full: distance matrix
  #' K: feature dimension
  #' lambda: l2 penalization parameter for feature extraction
  #' tol: convergence tolerance for feature extraction
  #' pca: logical. Whether to use PCA for the final features
  #' return: extracted regularized features and the target function value of each iteration
  #' required package: MASS
  n <- nrow(D_full)
  if (n > 5000) warning("Using the small dataset method for a large dataset!\n")
  
  d2 <- sum(D_full^2)
  V_reg <- diag(n + lambda, n) - matmult_cpp(matrix(rep(1, n), n, 1), matrix(rep(1, n), 1, n))
  V_reg_inv <- ginv(V_reg)
  B_fcn <- function(X) { # 
    B <- - D_full / as.matrix(dist(X))
    B[is.na(B)] <- 0
    B[is.infinite(B)] <- 0
    diag(B) <- - apply(B, 1, sum)
    
    h <- sum(diag(matmult_cpp(matmult_cpp(t(X), (V_reg - 2 * B)), X))) + d2
    return(list(B = B, h = h))
  }
  
  # initialization
  X_tmp <- cmdscale(D_full, k=K)
  res <- B_fcn(X_tmp)
  B_tmp <- res$B; h_tmp <- res$h
  obj_save <- c()
  num <- 0
  
  while (TRUE) {
    num <- num + 1
    X_update <- matmult_cpp(matmult_cpp(V_reg_inv, B_tmp), X_tmp)
    res <- B_fcn(X_update)
    B_update <- res$B; h_update <- res$h
    obj_save <- c(obj_save, h_tmp)
    if (abs((h_tmp - h_update)) < tol) {
      break
    }
    X_tmp <- X_update
    B_tmp <- B_update
    h_tmp <- h_update
  }
  
  if (pca) X_update <- prcomp(X_update, scale.=T)$x
  
  return(list(theta = X_update, obj_save = obj_save))
}

mds_l2_large <- function(D_full, K, init_obj, lambda=1e-1, tol, ncores, pca = FALSE){
  #' D_full: distance matrix
  #' K: feature dimension
  #' init_obj: the indices of fully-connected objects
  #' lambda: l2 penalization parameter for feature extraction. 84 v.s. 5
  #' tol: convergence tolerance for feature extraction
  #' ncores: number of cores to be used
  #' pca: logical. Wheter to use PCA for the final features
  #' return: extracted regularized features and the target function values of each iteration
  #' required package: MASS
  n <- nrow(D_full)
  theta <- matrix(0, n, K)
  subset_size <- length(init_obj)
  remn_obj <- setdiff(1:n, init_obj)
  
  res_init <- mds_l2_small(D_full[init_obj, init_obj], K, lambda, tol, pca=pca)
  theta[init_obj, ] <- res_init$theta
  obj_save <- res_init$obj_save
  
  obj_fun <- function(theta_j, theta_m_mat, d_vec) {
    theta_j_mat <- matmult_cpp(cbind(rep(1, subset_size)), t(theta_j))
    theta_diff <- theta_m_mat - theta_j_mat
    res <- sum((d_vec - sqrt(rowSums((theta_diff)^2)))^2) + lambda * sum(theta_j^2)
    res
  }
  grad_fun <- function(theta_j, theta_m_mat, d_vec) {
    theta_j_mat <- matmult_cpp(cbind(rep(1, subset_size)), t(theta_j))
    theta_diff <- theta_m_mat - theta_j_mat
    dhat_vec <- sqrt(rowSums((theta_diff)^2))
    res <- 2 * colSums((d_vec / dhat_vec - 1) * theta_diff) + 2 * lambda * theta_j
    res
  }
  opt_fcn <- function(d_vec){
    opt_res <- optim(rnorm(K), fn = obj_fun, gr = grad_fun, method = "BFGS", 
                     theta_m_mat = theta[init_obj, ], d_vec = d_vec)
    opt_res$par
  }
  
  dist_list <- split(D_full[init_obj, remn_obj], rep(1:length(remn_obj), each = subset_size))
  theta[remn_obj,] <- t(mcmapply(opt_fcn, dist_list, mc.cores = ncores))
  
  if (pca) theta <- prcomp(theta, scale.=T)$x
  
  return(list(theta = theta, obj_save = obj_save))
}

ylist_pred <- function(X, y_list, reg=T){
  #' X: feature matrix
  #' y_list: list of targets. Can contain NA but should be the same length as nrow(X)
  #' reg: logical. Whether to use l2 regularization for prediction
  #' return: a list of in-sample R2, SSE and reg for each y
  n <- nrow(X)
  p <- ncol(X)
  M <- length(y_list)
  
  X <- as.matrix(scale(X))
  
  pred <- function(y){
    idx_na <- which(is.na(y))
    if(length(idx_na) > 1) y <- y[-idx_na]; X <- X[-idx_na, ]  # only consider non-NA
    mod <- lm(y ~ ., data = data.frame(y = y, X))
    if(reg) {
      beta_est <- summary(mod)$coefficients[-1, 1]
      sig_est <- summary(mod)$sigma
      lambda <- p * sig_est^2 / sum(beta_est^2)  # estimated lambda that minimizes the out-of-sample SSE
      y_hat <- n / (n + lambda) * matmult_cpp(X, matrix(beta_est, nc = 1))
    } else {
      y_hat <- predict(mod, newdata = data.frame(X))
    }
    return(list(r.squared=cor(y, y_hat)^2, sse=sum((y-y_hat)^2)))
  }
  lapply(y_list, pred)
}

ylist_pred_OS <- function(X, y_list, idx_train, idx_test) {
  #' D_full: distance matrix including both training and testing indices
  #' y_list: list of targets. Can contain NA but should be the same length as nrow(D_full)
  #' K: feature dimension
  #' init_obj: the indices of fully-connected objects for MDS
  #' ncores: number of cores to be used
  #' idx_train: indices for the training objects
  #' idx_test: indices for the testing objects
  #' return: a list of in-sample and out-of-sample R2, SSE and reg for each y
  n <- nrow(X)
  if (any(sapply(y_list, length) != n)) {
    stop("Unmatched dimension for distance matrix and at least one target!\n")
  }
  r2_ridge_OS <- r2_ridge_IS <- r2_OLS_OS <- r2_OLS_IS <- rep(NA, length(y_list))
  for(i in 1:length(y_list)){
    idx_train_nonna <- idx_train; idx_test_nonna <- idx_test
    y <- y_list[[i]]
    idx_na <- which(is.na(y))
    if(length(idx_na) > 1) {
      idx_train_nonna <- intersect(idx_train_nonna, (1:n)[-idx_na])
      idx_test_nonna <- intersect(idx_test_nonna, (1:n)[-idx_na])
    }
    y_train <- y[idx_train_nonna]; y_test <- y[idx_test_nonna]
    X_train <- X[idx_train_nonna, ]; X_test <- X[idx_test_nonna, ]
    
    # ridge
    cv_res <- cv.glmnet(X_train, y_train)
    mod <- glmnet(X, y, lambda=cv_res$lambda.min)
    y_ridge_IS <- predict(mod, newx = X_train)
    y_ridge_OS <- predict(mod, newx = X_test)
    # OLS
    mod <- lm(y ~ ., data.frame(y_train, X_train))
    y_OLS_IS <- predict(mod, newdata = data.frame(X_train))
    y_OLS_OS <- predict(mod, newdata = data.frame(X_test))
    
    r2_ridge_IS[i] <- cor(y_ridge_IS, y_train)^2
    r2_ridge_OS[i] <- cor(y_ridge_OS, y_test)^2
    r2_OLS_IS[i] <- cor(y_OLS_IS, y_train)^2
    r2_OLS_OS[i] <- cor(y_OLS_OS, y_test)^2
  }
  
  return(list(r2_OLS_IS=r2_OLS_IS, r2_OLS_OS=r2_OLS_OS, 
              r2_ridge_IS=r2_ridge_IS, r2_ridge_OS=r2_ridge_OS))
}
 
NM_beta_update_multi <- function(D_list, y_list, K, init_obj, X_simplex, beta_simplex, loss_simplex, 
                                 l_beta, reg=T, alpha=1, gamma=2, rho=1/2, sigma=1/2, 
                                 ncores=detectCores()-1) {
  #' D_list: list of distance candidate matrices
  #' y_list: list of y
  #' K: feature dimension
  #' init_obj: indices of fully-connected objects
  #' X_simplex: list of X simplex: ordered
  #' beta_simplex: list of beta simplex: ordered
  #' loss_simplex: vector of loss simplex: ordered
  #' l_beta: l2 regularization parameter for beta 
  #' reg: logical. Whether to use l2 regularization for prediction
  #' ncores: number of cores to use
  
  dim_beta <- length(beta_simplex[[1]])
  beta_o <- colMeans(do.call(rbind, beta_simplex[1:dim_beta])) # mean of the first dim_beta vectors
  
  # Reflection
  type <- "R"
  beta_r <- beta_o + alpha * (beta_o - beta_simplex[[dim_beta+1]])
  beta_r <- sapply(beta_r, function(x) max(x, 0))
  D_sum <- Reduce("+", Map("*", D_list, c(1, beta_r)))
  X_r <- mds_large(D_sum/median(D_sum), K, init_obj, ncores)
  loss_r <- mean(sapply(ylist_pred(X_r, y_list, reg), "[[", "sse")) + 
    l_beta*mean(as.vector(beta_r)^2)
  if (loss_r < loss_simplex[dim_beta] & loss_r >= loss_simplex[1]) {
    beta_simplex[[dim_beta+1]] <- beta_r
    X_simplex[[dim_beta+1]] <- X_r
    loss_simplex[dim_beta+1] <- loss_r
  } else if (loss_r < loss_simplex[1]) {
    # Expansion
    type <- "E"
    beta_e <- beta_o + gamma * (beta_r - beta_o)
    beta_e <- sapply(beta_e, function(x) {max(x,0)})
    D_sum <- Reduce("+", Map("*", D_list, c(1, beta_e)))
    X_e <- mds_large(D_sum/median(D_sum), K, init_obj, ncores)
    loss_e <- mean(sapply(ylist_pred(X_e, y_list, reg), "[[", "sse")) + 
      l_beta*mean(as.vector(beta_e)^2)
    if (loss_e < loss_r){
      beta_simplex[[dim_beta+1]] <- beta_e
      X_simplex[[dim_beta+1]] <- X_e
      loss_simplex[dim_beta+1] <- loss_e
    } else {
      beta_simplex[[dim_beta+1]] <- beta_r
      X_simplex[[dim_beta+1]] <- X_r
      loss_simplex[dim_beta+1] <- loss_r
    }
  }
  else { # h_r >= loss_simplex[dim_beta]
    # Contraction
    type <- "C"
    if (loss_r < loss_simplex[dim_beta+1]) {
      beta_c <- beta_o + rho * (beta_r - beta_o)
      D_sum <- Reduce("+", Map("*", D_list, c(1, beta_c)))
      X_c <- mds_large(D_sum/median(D_sum), K, init_obj, ncores)
      loss_c <- mean(sapply(ylist_pred(X_c, y_list, reg), "[[", "sse")) + 
        l_beta*mean(as.vector(beta_c)^2)
      if(loss_c < loss_r){
        beta_simplex[[dim_beta+1]] <- beta_c
        X_simplex[[dim_beta+1]] <- X_c
        loss_simplex[dim_beta+1] <- loss_c
      } else {
        # shrink
        type <- "S"
        for(i in 2:(dim_beta+1)){
          beta_simplex[[i]] <- beta_simplex[[1]] + sigma * (beta_simplex[[i]] - beta_simplex[[1]])
          D_sum <- Reduce("+", Map("*", D_list, c(1, beta_simplex[[i]])))
          X_simplex[[i]] <- mds_large(D_sum/median(D_sum), K, init_obj, ncores)
          loss_simplex[i] <- mean(sapply(ylist_pred(X_simplex[[i]], y_list, reg), "[[", "sse")) + 
            l_beta*mean(as.vector(beta_simplex[[i]])^2)
        }
      }
    } else { # loss_r >= loss_n+1
      beta_c <- beta_o + rho * (beta_simplex[[dim_beta+1]] - beta_o)
      D_sum <- Reduce("+", Map("*", D_list, c(1, beta_c)))
      X_c <- mds_large(D_sum/median(D_sum), K, init_obj, ncores)
      loss_c <- mean(sapply(ylist_pred(X_c, y_list, reg), "[[", "sse")) + 
        l_beta*mean(as.vector(beta_c)^2)
      
      if (loss_c < loss_simplex[dim_beta+1]) {
        beta_simplex[[dim_beta+1]] <- beta_c
        X_simplex[[dim_beta+1]] <- X_c
        loss_simplex[dim_beta+1] <- loss_c
      } else {
        # shrink
        type <- "S"
        for(i in 2:(dim_beta+1)){
          beta_simplex[[i]] <- beta_simplex[[1]] + sigma * (beta_simplex[[i]] - beta_simplex[[1]])
          D_sum <- Reduce("+", Map("*", D_list, c(1, beta_simplex[[i]])))
          X_simplex[[i]] <- mds_large(D_sum/median(D_sum), K, init_obj, ncores)
          loss_simplex[i] <- mean(sapply(ylist_pred(X_simplex[[i]], y_list, reg), "[[", "sse")) + 
            l_beta*mean(as.vector(beta_simplex[[i]])^2)
        }
      }
    }
  }
  
  # Reorder by the value of h(X)
  order_idx <- sort(loss_simplex, index.return=T)$ix
  beta_simplex <- beta_simplex[order_idx]
  loss_simplex <- loss_simplex[order_idx]
  X_simplex <- X_simplex[order_idx]
  
  return(list(X_simplex=X_simplex, beta_simplex=beta_simplex, loss_simplex=loss_simplex, type=type))
}

main_large_NM_multi <- function(D_list, y_list, K, beta_simplex, init_obj, l_beta, 
                                reg=T, tol=1e-1, alpha=1, gamma=2, rho=1/2, sigma=1/2, 
                                minit=3, maxit=50, ncores=detectCores()-1, output=T){
  #' D_list: list of distance candidate matrices
  #' y_list: list of y 
  #' K: feature dimension
  #' beta_simplex: list of beta simplex initialization
  #' init_obj: indices of fully-connected objects
  #' l_beta: l2 regularity parameter for beta
  #' reg: logical. Whether to use l2 regularization for prediction
  #' tol: tolerance for objective value variance to terminate
  #' minit: minimum number of iterations
  #' maxit: maximum number of iterations
  #' return: 
  #' r2_cp: matrix of which each row is the R2 of the targets at each iteration
  if (var(sapply(beta_simplex, length)) > 0) {
    stop("Length of beta's are different!\n")
  }
  dim_beta <- length(beta_simplex[[1]])
  if (length(D_list) != dim_beta+1){
    stop("Number of distance matrices does not match beta's length")
  }
  n <- nrow(D_list[[1]])
  if(length(y_list) > 1) {
    if (var(sapply(y_list, length)) > 0 ) {
      stop("Ununiform lengths of targets!\n")
    }
  }
  if(length(y_list[[1]]) != n){
    stop("Unmatched target length and distance matrix dimension!\n")
  }
  if (any(sapply(y_list, function(y) all(is.na(y))))) {
    stop("At least one target vector is all NA's!\n")
  }
  if (K > n) {
    stop("Feature dimension larger than distance matrix dimension!\n")
  }
  
  X_simplex <- list()
  loss_simplex <- rep(NA, dim_beta+1)
  sse_IS_cp <- r2_IS_cp <- matrix(NA, nr=0, nc=length(y_list))
  
  # Initialize the simplex
  for(i in 1:(1+dim_beta)){
    D_sum <- Reduce("+", Map("*", D_list, c(1, beta_simplex[[i]])))
    X_simplex[[i]] <- mds_large(D_sum/median(D_sum), K, init_obj, ncores)
    pred_res <- ylist_pred(X_simplex[[i]], y_list, reg)
    loss_simplex[i] <- mean(sapply(pred_res, "[[", "sse")) + l_beta*mean(beta_simplex[[i]]^2)
    r2_IS_cp <- rbind(r2_IS_cp, sapply(pred_res, "[[", "r.squared"))
    sse_IS_cp <- rbind(sse_IS_cp, sapply(pred_res, "[[", "sse"))
  }
  order_idx <- sort(loss_simplex, index.return=T)$ix
  loss_simplex <- loss_simplex[order_idx]
  beta_simplex <- beta_simplex[order_idx]
  X_simplex <- X_simplex[order_idx]
  sse_IS_cp <- sse_IS_cp[order_idx, ]; r2_IS_cp <- r2_IS_cp[order_idx, ]
  beta_cp <- list(beta_simplex); loss_cp <- list(loss_simplex) 
  S <- sqrt(sum((loss_simplex - mean(loss_simplex))^2) / dim_beta)
  S_cp <- S
  
  # iteration
  iter <- 0
  while( ((S > tol) & (iter < maxit)) | (iter < minit)){
    t1 <- Sys.time()
    iter <- iter + 1
    res <- NM_beta_update_multi(D_list, y_list, K, init_obj, X_simplex, beta_simplex, loss_simplex,
                                l_beta, reg, alpha, gamma, rho, sigma, ncores)
    beta_simplex <- res$beta_simplex; beta_cp <- append(beta_cp, beta_simplex)
    loss_simplex <- res$loss_simplex; loss_cp <- rbind(loss_cp, loss_simplex)
    X_simplex <- res$X_simplex
    S <- sqrt(sum((loss_simplex-mean(loss_simplex))^2) / dim_beta)
    pred_res <- ylist_pred(X_simplex[[1]], y_list, reg)
    r2_IS_cp <- rbind(r2_IS_cp, sapply(pred_res, "[[", "r.squared"))
    sse_IS_cp <- rbind(sse_IS_cp, sapply(pred_res, "[[", "sse"))
    S_cp <- c(S_cp, S)
    
    # printout
    if(output){
      if(dim_beta == 3) {
        cat("iter=", iter, "beta1=", sprintf("%.2f", beta_simplex[[1]]), 
            "loss=", sprintf("%.2f", loss_simplex[1:2]), "S=", round(S, 2), 
            "r2=", round(mean(r2_IS_cp[nrow(r2_IS_cp), ]), 3), res$type, "\n")
      }
      if(dim_beta == 2) {
        cat("iter=", iter, "beta1=", sprintf("%.2f", beta_simplex[[1]]), 
            "loss=", sprintf("%.2f", loss_simplex[1:2]), "S=", round(S, 2), 
            "r2=", round(mean(r2_IS_cp[nrow(r2_IS_cp), ]), 3), res$type, "\n")
      }
      if(dim_beta == 1) {
        cat("iter=", iter, "beta1=", sprintf("%.2f", beta_simplex[[1]]), 
            "loss=", sprintf("%.2f", loss_simplex[1:2]), "S=", round(S, 2), 
            "r2=", round(mean(r2_IS_cp[nrow(r2_IS_cp), ]), 3), res$type, "\n")
      }
      cat("t_iter="); print(Sys.time()-t1)
    }
  }
  cat("Final mean IS r2 =", round(mean(r2_IS_cp[nrow(r2_IS_cp), ]), 2), "\n")
  
  res <- list(beta_est=beta_simplex[[1]], theta=X_simplex[[1]], y_list=y_list, S_cp=S_cp,
              beta_cp=beta_cp, loss_cp=loss_cp, sse_IS_cp=sse_IS_cp, r2_IS_cp=r2_IS_cp, iter=iter)
  return(res)
}







