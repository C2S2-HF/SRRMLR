#' Generate data for simulation
#' @param n number of observations
#' @param p dimension of predictor vector
#' @param q dimension of response vector
#' @param s sparsity level of the coefficient matrix
#' @param r the specified rank
#' @param type 'AR' for Auto Regressive, 'CS' for Compound Symmetries
#' @param C the true coefficient matrix
#' @return Y generated response matrix
#' @return X generated covariate matrix
#' @return C generated true coefficient matrix
#' @return B generated B matrix
#' @return V generated V matrix
#' @importFrom stats runif
#' @importFrom stats rmultinom
#' @importFrom MASS mvrnorm
#' @export
gendata <- function(n, p, q = 8, s = 5, r = 5, type, C = NULL){
  Sigma <- array(NA, dim = c(p, p))
  if(type == "CS"){
    Sigma <- 0.5*diag(p)+matrix(0.5,nrow=p,ncol=p)
  }
  else{
    for(i in 1:p){
      for(j in 1:p){
        Sigma[i,j] <- 0.5^(abs(i-j))
      }
    }
  }
  X <- mvrnorm(n, rep(0,p), Sigma)
  B <- NULL
  V <- NULL
  if(is.null(C)) {
    B1 <- matrix(ncol = r, nrow = s, sample(c(-1,1), r*s, replace=TRUE) * runif(r*s,1/q,2/q))
    B2 <- matrix(ncol = r, nrow = p-s, 0)
    B <- rbind(B1, B2)
    V <- matrix(ncol = r, nrow = q, runif(q*r, 0, q))
    C <- B %*% t(V)
  }
  z1 <- as.matrix(exp(X%*%C)[,-1])
  d <- apply(z1,1,sum)+1
  p1 <- as.matrix(1/d)
  p2 <- z1/d
  temp <- cbind(p1,p2)
  Y <- t(apply(temp, 1, rmultinom, n = 1, size = 1))
  list(Y = Y, X = X, C = C, B = B, V = V)
}
