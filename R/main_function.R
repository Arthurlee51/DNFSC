#' Stability-Based Factor Number Estimation
#'
#' This is the core function of the package. It computes the stability-based criteria—SC1, SC2, and SC3—and provides corresponding estimates for the number of factors, as outlined in Lee & Chen (2024+).
#'
#'
#' @param X A numeric data matrix matrix of dimensions \eqn{n \times p}.
#' @param K.cand An integer vector specifying the candidate set for the number of factors to be evaluated. The default is \code{1:10}.
#' @param averaging A logical parameter indicating whether to perform multiple data splits and average the resulting instability measures. If \code{TRUE} (default), the function conducts multiple splits.
#' @return A list containing the results of the estimation, including:
#' \item{SC1}{Computed values of SC1 for each value in \code{K.cand}.}
#' \item{SC2}{Computed values of SC2 for each value in \code{K.cand}.}
#' \item{SC3}{Computed values of SC3 for each value in \code{K.cand}.}
#' \item{K_hat_SC1}{Estimated number of factors under SC1.}
#' \item{K_hat_SC2}{Estimated number of factors under SC2.}
#' \item{K_hat_SC3}{Estimated number of factors under SC3.}
#' \item{instability}{Vector of loading instabilities for each value in \code{K.cand}.}
#' @examples
#' # Generate example dataset
#' set.seed(123)
#' X <- example_data.func(scenario =1, signal = 1, n = 1500, p=500)
#' result <- dnfsc_sc(X, K.cand = 1:10,averaging = TRUE )
#' @export
dnfsc_sc <-function(X, K.cand = 1:10, averaging = TRUE) {
  #Get n and p.
  n <- nrow(X)
  p <- ncol(X)

  #n_split: number of splits.
  if(averaging==TRUE){
    n_split = 10
  }else{
    n_split = 1
  }


  #Get first terms of SC1, SC2 and SC3.
  #c_values for SC1 using deterministic c_k
  c_values <- 1 - K.cand/length(K.cand)

  #fit_SC2: First term of SC2
  fit_SC2 <- Get_fit_SC2.func(X,K.cand,p,n)

  #fit_SC3: First term of SC3
  fit_SC3 <- Get_fit_SC3.func(X,K.cand,p,n)
  #=======================================================================================
  #instability: matrix storing the instability
  instab_mat <- matrix(0, nrow = n_split, ncol = length(K.cand))

  #Splitting procedure
  for( s in 1:n_split){
    #Split the data in half to compute stability measure
    X1_ind <- sample(1:n, round(n/2),  replace = FALSE)
    X1 <- X[X1_ind,]
    X2 <- X[setdiff(1:n, X1_ind), ]


    #Perform SVD to X1 and X2
    #Svd to X1/sqrt{p*n} to get \tilde{v}_j. Scaled by sqrt(p*n) to ensure numerical stability.
    svd_1 = svd(X1/sqrt(p*n))
    svd_2 = svd(X2/sqrt(p*n))

    #Compute stability measure
    for ( i in 1:length(K.cand)){
      k = K.cand[i]
      instab_mat[s,i]<- subspace(svd_1$v[,1:k],svd_2$v[,1:k])
    }
  }
  instability <- colMeans(instab_mat)
  #Compute SC1, SC2 and SC3
  SC1<- c_values + instability
  SC2 <- fit_SC2 + instability
  SC3 <- fit_SC3 + instability
  K_hat_SC1 <- K.cand[which.min(SC1)]
  K_hat_SC2 <- K.cand[which.min(SC2)]
  K_hat_SC3 <- K.cand[which.min(SC3)]

  return(list(SC1 = SC1, SC2 = SC2, SC3 = SC3, K_hat_SC1 = K_hat_SC1, K_hat_SC2 = K_hat_SC2, K_hat_SC3 = K_hat_SC3, instability = instability ))
}


