#Script containing all supporting functions.
#Function to compute the sin(angle) between two vectors
sin_vec <- function(vector1, vector2) {
  # Ensure vectors are numeric and have the same length
  if(!is.numeric(vector1) | !is.numeric(vector2) | length(vector1) != length(vector2)) {
    stop("Both vectors must be numeric and of the same length")
  }
  
  # Compute the dot product of the two vectors
  dot_product <- sum(vector1 * vector2)
  
  # Compute the magnitudes of the vectors
  magnitude1 <- sqrt(sum(vector1^2))
  magnitude2 <- sqrt(sum(vector2^2))
  
  # Compute the cosine of the angle
  cos_angle <- dot_product / (magnitude1 * magnitude2)
  
  # Compute the angle in radians
  angle_radians <- acos(cos_angle)
  
  sin(angle_radians)
  
  return(sin(angle_radians))
}


#Function to compare the sin angle between two matrices.
sin_matrix <- function(A, B){
  if (nrow(A) != nrow(B)) {
    stop("Row dimensions of A and B must be the same.")
  }
  
  # Compute orthonormal bases using SVD for stability
  A <- svd(A)$u  # Use qr.Q for orthonormal basis
  B <- svd(B)$u
  
  # Define threshold for small angles
  threshold <- sqrt(2) / 2
  
  # Most accurate method (adapted from reference)
  s <- svd(t(A) %*% B)$d  # Use t() for transpose
  costheta <- min(s)  # Cosine of the angle
  
  # Check for small angles
  if (costheta < threshold) {
    theta <- acos(min(1, costheta))  # Angle from cosine
  } else {
    # Ignore angles due to different sizes (adapted from reference)
    if (ncol(A) < ncol(B)) {
      sintheta <- norm(t(A) - (t(A) %*% B) %*% t(B))  # Sine from difference due to B
    } else {
      sintheta <- norm(t(B) - (t(B) %*% A) %*% t(A))  # Sine from difference due to A
    }
    theta <- asin(min(1, sintheta))  # Angle from sine (recompute)
  }
  sin(theta)
}

subspace<- function(A,B){
  if(is.vector(A)){
    out <- sin_vec(A,B)
  }
  
  if(is.matrix(A)){
    out <- sin_matrix(A,B)
  }
  return(out)
}

#recover the matrix based on the first k eigenvectors
svd_recover.func <- function(svd_object,k){
  if(k==1){
    out <- svd_object$d[1]*svd_object$u[,1]%*% t(svd_object$v[,1])
  }else{
    out <- matrix(svd_object$u[,1:k], ncol=k) %*% matrix(diag(svd_object$d[1:k]),ncol=k) %*% t(matrix(svd_object$v[,1:k],ncol=k))
  }
  return(out)
}

#Function to get the first term of SC2
Get_fit_SC2.func = function(X,K.cand,p,n){
  log_full_eigen_values <- log(p*svd( t(X)%*%X/(p*n) )$d[1:length(K.cand)] +1  )
  #Compute l(0), l(1), ... l(K_{max}). l_values(k) is equivalent to l(k-1) in the paper.
  l_values <- rep(0, length(K.cand))
  l_values[1] <-  sum(log_full_eigen_values)
  for( k in 2:length(K.cand)){
    l_values[k] <- l_values[k-1] - log_full_eigen_values[k-1] 
  }
  #Return fit: First term of SC2
  return(c(l_values[-1]/l_values[1],0))
}

Get_fit_SC3.func = function(X,K.cand,p,n){
  svd_object = svd(X)
  denom = log(sum(X^2)/(n*p) +1)
  fit <- rep(0, length(K.cand))
  for ( i in 1: length(K.cand)){
    k = K.cand[i]
    #Get C, the common component estimated by the kth factor 
    C <- svd_recover.func(svd_object,k)
    V = sum( (X -C)^2 )/(n*p)
    fit[i] = log(V+1)/denom 
  }
  return(fit)
}



#Information criterion proposed in Bai and Ng (2002). Only relevant for reproduction of result
IC.func = function(X,K.cand){
  n <- nrow(X)
  p <- ncol(X)
  #Perform svd first
  svd_object = svd(X)
  #penalty term
  pen= ( (n+p)/(n*p) )* log( (n*p)/(n+p) )
  IC = rep(0, length(K.cand))
  #recover k=1, and compute fit measure
  denom = log(sum(X^2)/(n*p) +1)
  for ( i in 1: length(K.cand)){
    k = K.cand[i]
    #Get C, the common component estimated by the kth factor 
    C <- svd_recover.func(svd_object,k)
    V = sum( (X -C)^2 )/(n*p)
    IC[i] = log(V) + k*pen
  }
  Khat_IC <- K.cand[which.min(IC)]
  return( list(Khat_IC =Khat_IC, IC = IC  ) )
}

#Functions to output AIC and BIC criteria proposed in Bai et al (2018)
BPY.func = function(X, K.cand){
  n <- nrow(X)
  p <- ncol(X)
  #Perform svd first
  full_eigen_values <- p*svd( t(X)%*%X/(p*(n-1) )) $d
  #Initiate AIC and BIC
  AIC = rep(0, length(K.cand))
  BIC = rep(0, length(K.cand))
  
  for ( i in 1: length(K.cand)){
    k = K.cand[i]
    #Get C, the common component estimated by the kth factor 
    bar_sigma <- sum(full_eigen_values[(i+1):p ] )/(p-i)
    ICprep <- (p-i)*log(bar_sigma) - sum(log(full_eigen_values[(i+1):p]) )
    pen1= (p-k-1)*(p-k+2)/n
    pen2 = ( (p-k-1)*(p-k+2)/(2*n) )*log(n)
    AIC[i] = ICprep  - pen1
    BIC[i] = ICprep - pen2
  }
  Khat_A <- K.cand[which.min(AIC)]
  Khat_B <- K.cand[which.min(BIC)]
  return(list(Khat_A = Khat_A, Khat_B = Khat_B, AIC = AIC, BIC = BIC))
}
