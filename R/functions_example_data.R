#' Example Data Generator
#'
#' This function generates example data according to the settings in Section 4 of Lee & Chen (2024+).
#'
#'
#' @param scenario An integer specifying the error term distribution: 1 for homogeneous error and 2 for heterogeneous error.
#' @param signal An integer indicating the signal strength: 1 for strong signal, 2 for weak signal, and 3 for a mix of equal and mixed-order signals.
#' @param n An integer specifying the sample size.
#' @param p An integer specifying the number of features.
#' @return A numeric matrix of dimensions \eqn{n \times p} representing the generated data.
#' @examples
#' X <- example_data.func(scenario =1, signal = 1, n = 1500, p=500)
#' @export
example_data.func <- function(scenario,signal, n, p) {
  #K = 4 across all situations
  K=4
  #Signal matrix, where the rows correspond to conditions (i), (ii) and (iii)
  Signal_matrx <- rbind(sqrt(p)*c(6,5,4,3),p^{1/6}*c(6,5,3,3), c(3*p^{1/3},3*p^{1/3},3*p^{1/6}, 3*p^{1/6}) )

  if (scenario == 1 ) {
    Epsilon <- matrix(rnorm(n*p) , nrow = n, ncol = p)
    Gamma <- matrix(runif(n*K,-0.5,0.5), nrow = n, ncol = K)
  }

  if (scenario ==2 ) {
      Q_epsilon <- matrix(rnorm(n*n), nrow = n, ncol = n)
      Q_epsilon <-qr.Q(qr(Q_epsilon))
      Epsilon <-  diag(sqrt((1:n)/n ))%*%Q_epsilon%*%matrix(rt(n*p,10), nrow = n, ncol = p)
      Gamma <- matrix(runif(n*K,-0.5,0.5), nrow = n, ncol = K)
    }

  #Get D: The singular values of \Lambda
  D = diag(Signal_matrx[signal,])
  #Get Q: normalised matrix from iid normal
  Z <- matrix(rnorm(p*K), nrow = p, ncol = K)
  qr_decomp <- qr(Z)
  Q <- qr.Q(qr_decomp)
  Lambda = Q%*%D
  X <-Gamma%*%t(Lambda) + Epsilon
  return(X)
}

