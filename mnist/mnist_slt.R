require(Rcpp)
library(SRRMLR)
library(glmnet)
library(nnet)
library(VGAM)
library(dslabs)
library(corpcor)
library(doParallel)
# mnist <- read_mnist(destdir="./")
# data.all <- rbind(mnist$train$images,mnist$test$images)
# labels <- c(mnist$train$labels,mnist$test$labels)
# 
# set.seed(1234)
# pca <- prcomp(data.all[, 1:784], retx = TRUE)

load(file="PCA.RData")
load(file="labels.RData")
pca_s <- 100
var_per <- sum((pca$sdev^2/sum(pca$sdev^2))[1:pca_s])
var_per
data.all <- pca$x[,1:pca_s]

f <- function(Y, X, C) {
  z <- X%*%C
  z <- exp(z)
  z <- log(apply(z,1,sum)+1)
  loglike <- sum(diag(X%*%C%*%t(Y[,-1])))-sum(z)
  return(-loglike)
}


real <- function(data.all, labels, n.train=500){
  
  n <- 1000
  p <- 100
  q <- 10
  
  obj <-array(0,dim=c(n, p+q))
  for(i in 1:q){
    set <- sample(which(labels==(i-1)), 100)
    obj[(100*(i-1)+1):(100*i),1:100] <- data.all[set,]
    obj[(100*(i-1)+1):(100*i),(100+i)] <- 1
  }
  obj[,1:100] <-scale(obj[,1:100])
  
  train <- sample(1:n, n.train)
  X.train <- obj[train,1:100]
  Y.train <- obj[train,101:110]
  X.test <- obj[-train,1:100]
  Y.test <- obj[-train,101:110]
  
  m <- q - 1
  size <- round(n.train/5)
  
  Pred <- r.hat <- s.hat <- rep(NA, 4)
  
  s.max = ceiling(20*(n/(q*log(p*q)))^(1/4))
  r.max = min(s.max,q-1,ceiling(20*(n/(q*log(p*q)))^(1/3)))

  print("srrmlr begins")
  r_list <- 5:r.max
  s_list <- 15:s.max
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
    
    Pred[2] <- unglm.Pred
    s.hat[2] <- unglm.s
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
  
  print('rrvglm begins')
  err.vglm <- rep(NA, r.max)
  for(i in 1:r.max){
    sum_rr <- rep(NA, 5)
    for(k in 1:5){
      try({
        set <- ((k-1)*size+1):(k*size)
        fit.vglm <- rrvglm(Y.train[-set,] ~ X.train[-set,], multinomial(refLevel = 1), Rank = i)
        sum_rr[k] <- f(Y = Y.train[set,], X = X.train[set,], C = coef(fit.vglm, matrix = TRUE)[-1,])
      }, silent = TRUE)
    }
    err.vglm[i] <- mean(sum_rr,na.rm=TRUE)
  }
  if(all(is.na(err.vglm))){}else{
    vglm.r <- which.min(err.vglm)
    try({
      vglm.fit <- rrvglm(Y.train ~ X.train, multinomial(refLevel = 1), Rank = vglm.r)
      vglm.C <- coef(vglm.fit, matrix = TRUE)[-1,]
      vglm.aset <- which(apply(vglm.C, 1, sum) != 0)
      vglm.s <- length(vglm.aset)
      vglm.Pred <- f(Y = Y.test, X = X.test, C = vglm.C)
      if(vglm.Pred == Inf){} else {
        Pred[4] <- vglm.Pred
        r.hat[4] <- vglm.r
        s.hat[4] <- vglm.s
      }
    }, silent = TRUE)
  }
  return(cbind(Pred, r.hat, s.hat))
}

n.rep <- 100
cl = makeCluster(10)
registerDoParallel(cl)
fit <- foreach(i = 1:n.rep, .packages = c("glmnet", "SRRMLR", "nnet", "VGAM")) %dopar% real(data.all=data.all, labels=labels, n.train=500)
stopImplicitCluster()
stopCluster(cl)
save(fit,file="mnistslt5.RData")

