source("Function.R")
require(doParallel)
require(glmnet)
require(nnet)
require(SRRMLR)
# require(Rcpp)
# require(RcppArmadillo)
# sourceCpp("srrmlr.cpp")

data.all <- read.csv("ham_rgb.csv")

set.seed(123)

n <- dim(data.all)[1]    
p <- dim(data.all)[2] - 1   
X <- as.matrix(data.all[,1:p])
names <- colnames(X)
X <- scale(X)  
y_label = as.factor(data.all$label)
q <- length(levels(y_label))  
m <- q - 1
Y <- array(0, dim = c(n, q))   
colnames(Y) <- levels(y_label)
for(lab in levels(y_label)){
  set <- which(y_label==lab)
  idx = as.numeric(lab)+1
  Y[set,idx] <- 1
}

real <- function(X, Y){
  n <- dim(X)[1]
  p <- dim(X)[2]
  q <- dim(Y)[2]
  m <- q - 1
  n.train = floor(0.7*n)
  train <- sample(1:n, n.train)
  X.train <- X[train,]
  Y.train <- Y[train,]
  X.test <- X[-train,]
  Y.test <- Y[-train,]
  
  Pred <- rep(NA, 3)
  r.hat <- rep(NA, 3)
  s.hat <- rep(NA, 3)
  
  s.max = ceiling(10*(n/(q*log(p*q)))^(1/4))
  r.max = min(s.max,q-1,ceiling(10*(n/(q*log(p*q)))^(1/3)))
 
  print('srrmlr begins')
  r_list <- 3:r.max
  s_list <- 10:s.max
  bess <- srrmlr_real(X.train, Y.train, X.test, Y.test, s_list, r_list, cv=5)

  Pred[1] <- bess$Pred
  r.hat[1] <- bess$r
  s.hat[1] <- bess$s

  print('mlasso begins')
  try({
    unglm.fit <- cv.glmnet(X.train, Y.train, family = "multinomial", type.measure = "deviance", nfolds = 5, type.multinomial = "ungrouped")
    c1 <- coef(unglm.fit, s = "lambda.min")
    c2 <- NULL
    for(j in 2:q){
      z <- as.matrix(c1[[j]])[-1]
      colnames(z) <- rownames(z) <- NULL
      c2 <- cbind(c2, z)
    }
    unglm.C <- c2
    unglm.aset <- which(apply(unglm.C, 1, sum) != 0)
    unglm.s <- length(unglm.aset)
    unglm.Pred <- f(Y = Y.test, X = X.test, C = unglm.C)
    if(is.null(unglm.Pred)){}else{
      if(unglm.Pred == Inf){}else{
        Pred[2] <- unglm.Pred
        s.hat[2] <- unglm.s
      }
    }
   }, silent=TRUE)
  
  print('multinom begins')
  try({
    mul.fit <- multinom(Y.train ~ X.train, MaxNWts = 5000, trace = FALSE)
    mul.C <- t(coef(mul.fit))[-1,]
    mul.aset <- which(apply(mul.C, 1, sum) != 0)
    mul.s <- length(mul.aset)
    mul.Pred <- f(Y = Y.test, X = X.test, C = mul.C)
    if(is.null(mul.Pred)){}else{
      if(mul.Pred == Inf){}else{
        Pred[3] <- mul.Pred
        s.hat[3] <- mul.s
      }
    }
  },silent=TRUE)
  
  return(cbind(Pred, r.hat, s.hat))
}

set.seed(123)
n.rep <- 100

cl = makeCluster(10)
registerDoParallel(cl)
fit <- foreach(i = 1:n.rep, .packages = c("glmnet", "SRRMLR", "nnet")) %dopar% real(X=X,Y=Y)
stopImplicitCluster()
stopCluster(cl)

save(fit,file="ham.RData")


