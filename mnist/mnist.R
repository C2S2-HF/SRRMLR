require(Rcpp)
require(RcppArmadillo)
require(SRRMLR)
library(nnet)
# library(dslabs)
# library(corpcor)
# mnist <- read_mnist(destdir="./")

set.seed(1234)
# load("raw.Rdata")
# data.all <- rbind(mnist$train$images,mnist$test$images)
# labels <- c(mnist$train$labels,mnist$test$labels)
# pca <- prcomp(data.all[, 1:784], retx = TRUE)
# save(pca, file="PCA.RData")
# save(labels, file="labels.RData")


load("PCA.RData")
load("labels.RData")
pca_s <- 100
var_per <- sum((pca$sdev^2/sum(pca$sdev^2))[1:pca_s])
var_per
data.all <- pca$x[,1:pca_s]

n <- dim(data.all)[1]
p <- dim(data.all)[2]
q <- 10
X <- as.matrix(data.all)
names <- colnames(X)
X <- scale(X)
Y <- array(0, dim = c(n, q))
for(i in 1:10){
  set <- which(labels==(i-1))
  Y[set,i] <- 1
}

f <- function(Y, X, C) {
  z <- X%*%C
  z <- exp(z)
  z <- log(apply(z,1,sum)+1)
  loglike <- sum(diag(X%*%C%*%t(Y[,-1])))-sum(z)
  return(-loglike)
}


real <- function(X, Y){
  n <- dim(X)[1]
  p <- dim(X)[2]
  q <- dim(Y)[2]
  m <- q - 1
  n.train = 1:60000
  X.train <- X[n.train,]
  Y.train <- Y[n.train,]
  X.test <- X[-n.train,]
  Y.test <- Y[-n.train,]
  Pred <- r.hat <- s.hat <- time <- rep(NA, 3)
  r.max = 9
  
  indices <- sample(1:60000, 60000)
  X.train <- X.train[indices,]
  Y.train <- Y.train[indices,]
  print("srrmlr begins")
  r_list <- 9
  s_list <- 99
  bess <- srrmlr_real(X.train, Y.train, X.test, Y.test, s_list, r_list, cv=5)

  Pred[1] <- bess$Pred
  r.hat[1] <- bess$r
  s.hat[1] <- bess$s
  time[1] <- bess$time
  
  print('mlasso begins')
  try({
    t1 <- Sys.time()
    unglm.fit <- cv.glmnet(X.train, Y.train, family = "multinomial", type.measure = "deviance", nfolds = 5, type.multinomial = "ungrouped")
    t2 <- Sys.time()
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
    
    Pred[2] <- unglm.Pred
    s.hat[2] <- unglm.s
    time[2] <- t2-t1
  }, silent=TRUE)
  
  print('multinom begins')
  try({
    t1 <- Sys.time()
    mul.fit <- multinom(Y.train ~ X.train, MaxNWts = 5000, trace = FALSE)
    t2 <- Sys.time()
    mul.C <- t(coef(mul.fit))[-1,]
    mul.aset <- which(apply(mul.C, 1, sum) != 0)
    mul.s <- length(mul.aset)
    mul.Pred <- f(Y = Y.test, X = X.test, C = mul.C)
    if(is.null(mul.Pred)){}else{
      if(mul.Pred == Inf){}else{
        Pred[3] <- mul.Pred
        s.hat[3] <- mul.s
        time[3] <- t2-t1
      }
    }
  },silent=TRUE)
  
  return(cbind(Pred, r.hat, s.hat, time))
}

fit <- real(X=X, Y=Y)

save(fit,file="mnist.RData")
