require(Rcpp)
require(RcppArmadillo)
require(SRRMLR)

# sourceCpp("srrmlr.cpp")
# gendata <- function(n, p, q = 8, s = 5, r = 5, type, C = NULL) {
#   ## Generate the covariance matrix of predictors
#   Sigma <- array(NA, dim = c(p, p))
#   if(type == "CS"){
#     Sigma <- 0.5*diag(p)+matrix(0.5,nrow=p,ncol=p)
#   }
#   else{
#     for(i in 1:p){
#       for(j in 1:p){
#         Sigma[i,j] <- 0.5^(abs(i-j))
#       }
#     }
#   }
#   X <- mvrnorm(n, rep(0,p), Sigma)
#   X <- cbind(rep(1, n), X)
#   B <- NULL
#   V <- NULL
#   if(is.null(C)) {
#     B1 <- matrix(ncol = r, nrow = s+1, runif((s+1) * r, -2/q, 2/q))
#     B2 <- matrix(ncol = r, nrow = p-s, 0)
#     B <- rbind(B1, B2)
#     V <- matrix(ncol = r, nrow = q, runif(q*r, 0, q))
#     C <- B %*% t(V)
#   }
#   z1 <- as.matrix(exp(X%*%C)[,-1])
#   d <- apply(z1,1,sum)+1
#   p1 <- as.matrix(1/d)
#   p2 <- z1/d
#   temp <- cbind(p1,p2)
#   Y <- t(apply(temp, 1, rmultinom, n = 1, size = 1))
#   list(Y = Y, X = X, C = C, B = B, V = V)
# }

type <- "CS"
n <- 400
p <- 200
q <- 8
s <- r <- 5
data <- gendata(n = n, p = p, q = q, s = s, r = r, type = type)
valdata <- gendata(n = 1000, p = p, q = q, s = s, r = r, type = type, C = data$C)
testdata <- gendata(n= 1000, p = p, q = q, s = s, r = r, type = type, C = data$C)

s.max = ceiling(5*(n/(q*log(p*q)))^(1/4))
r.max = min(s.max,q-1,ceiling(5*(n/(q*log(p*q)))^(1/3)))
s_list <- 4:s.max
r_list <- 4:r.max
#srrmlr_one(valdata$X, valdata$Y, s.max, r.max)
srrmlr_simu(data$X, data$Y, valdata$X, valdata$Y, testdata$X, testdata$Y, s_list, r_list, data$C, s, r)

s_list <- 4:6
r_list <- 4:6
#srrmlr_real(data$X, data$Y, testdata$X, testdata$Y, s_list, r_list, cv=5)


