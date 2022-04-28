require(doParallel)
require(MASS)
require(nnet)
require(glmnet)
require(VGAM)
require(SRRMLR)

## Compute negative log likelihood value
f <- function(Y, X, C) {
  z <- exp(X%*%C)
  z <- log(apply(z,1,sum)+1)
  loglike <- sum(diag(X%*%C%*%t(Y[,-1])))-sum(z)
  return(-loglike)
}

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


longwait <- function(n, p, type){
  q <- 8
  s <- r <- 5
  data <- gendata(n = n, p = p, q = q, s = s, r = r, type = type)
  valdata <- gendata(n = 1000, p = p, q = q, s = s, r = r, type = type, C = data$C)
  testdata <- gendata(n= 1000, p = p, q = q, s = s, r = r, type = type, C = data$C)
  Y <- data$Y
  X <- data$X
  C <- data$C
  m <- q - 1
  P.set <- 1:s
  N.set <- setdiff(1:p, P.set)
  C.hat <- array(NA, dim = c(5, p, m))
  Pred <- Est <- r.hat <- s.hat <- Sen <- Spe <- time <- rep(NA, 5)
  s.max = ceiling(5*(n/(q*log(p*q)))^(1/4))
  r.max = min(s.max,q-1,ceiling(5*(n/(q*log(p*q)))^(1/3)))
  
  print("srrmlr begins")
  t1 <- Sys.time()
  s_list <- 1:s.max
  r_list <- 1:r.max
  bess <- srrmlr_simu(data$X, data$Y, valdata$X, valdata$Y, testdata$X, testdata$Y, s_list, r_list, data$C, s, r)
  t2 <- Sys.time()
  Pred[1] <- bess$Pred
  Est[1] <- bess$Est
  s.hat[1] <- bess$s
  r.hat[1] <- bess$r
  Sen[1] <- bess$Sen
  Spe[1] <- bess$Spe
  time[1] <- difftime(t2,t1,units="s")

  print("mlasso begins")
  t1 <- Sys.time()
  try({
    fit.unglm <- glmnet(X, Y, family = "multinomial", type.multinomial = "ungrouped")
    C.unglm <- array(0, dim = c(length(fit.unglm$lambda), p, m))
    err.unglm <- rep(Inf, length(fit.unglm$lambda))
    for(i in 1:length(fit.unglm$lambda)){
      c1 <- coef(fit.unglm, fit.unglm$lambda[i])
      c2 <- NULL
      for(j in 2:q){
        z <- as.matrix(c1[[j]])[-1]
        colnames(z) <- rownames(z) <- NULL
        c2 <- cbind(c2, z)
      }
      C.unglm[i,,] <- c2
      err.unglm[i] <- f(Y = valdata$Y, X = valdata$X, C = c2)
    }
    id.unglm <- which.min(err.unglm)
    unglm.C <- C.unglm[id.unglm,,]
    unglm.aset <- which(apply(unglm.C, 1, sum)!= 0)
    unglm.inaset <- setdiff(1:p, unglm.aset)
    unglm.s <- length(unglm.aset)
    unglm.Est <- sum((unglm.C-C[, -1])^2)/(p*m)
    unglm.Pred <- f(Y = testdata$Y, X = testdata$X, C = unglm.C)
    unglm.Sen <- length(intersect(unglm.aset, P.set))/length(P.set)
    unglm.Spe <- length(intersect(unglm.inaset, N.set))/length(N.set)
   t2 <- Sys.time()
    if(is.null(unglm.Pred)){}else{
      if(unglm.Pred == Inf){}else{
        Pred[2] <- unglm.Pred
        Est[2] <- unglm.Est
        s.hat[2] <- unglm.s
        Sen[2] <- unglm.Sen
        Spe[2] <- unglm.Spe
        time[2] <- difftime(t2,t1,units="s")
      }
    }
  },silent=TRUE)
  
  # print("glmnet begins")
  # t1 <- Sys.time()
  # try({
  #   fit.glm <- glmnet(X, Y, family = "multinomial", type.multinomial = "grouped")
  #   C.glm <- array(0, dim = c(length(fit.glm$lambda), p, m))
  #   err.glm <- rep(Inf, length(fit.glm$lambda))
  #   for(i in 1:length(fit.glm$lambda)){
  #     c1 <- coef(fit.glm, fit.glm$lambda[i])
  #     c2 <- NULL
  #     for(j in 2:q){
  #       z <- as.matrix(c1[[j]])[-1]
  #       colnames(z) <- rownames(z) <- NULL
  #       c2 <- cbind(c2, z)
  #     }
  #     C.glm[i,,] <- c2
  #     err.glm[i] <- f(Y = valdata$Y, X = valdata$X, C = c2)
  #   }
  #   id.glm <- which.min(err.glm)
  #   glm.C <- C.glm[id.glm,,]
  #   glm.aset <- which(apply(glm.C, 1, sum) != 0)
  #   glm.inaset <- setdiff(1:p, glm.aset)
  #   glm.s <- length(glm.aset)
  #   glm.Est <- sum((glm.C - C[, -1])^2)/(p*m)
  #   glm.Pred <- f(Y = testdata$Y, X = testdata$X, C = glm.C)
  #   glm.Sen <- length(intersect(glm.aset, P.set))/length(P.set)
  #   glm.Spe <- length(intersect(glm.inaset, N.set))/length(N.set)
  #  t2 <- Sys.time()
  #   if(is.null(glm.Pred)){}else{
  #     if(glm.Pred == Inf){}else{
  #       Pred[3] <- glm.Pred
  #       Est[3] <- glm.Est
  #       s.hat[3] <- glm.s
  #       Sen[3] <- glm.Sen
  #       Spe[3] <- glm.Spe
  #       time[3] <- difftime(t2,t1,units="s")
  #       }
  #     }
  #   },silent=TRUE)

  print("nnet begins")
  
  try({
    t1 <- Sys.time()
    fit.mul <- multinom(Y ~ X[, -1], MaxNWts = 6000, trace = FALSE)
    mul.C <- t(coef(fit.mul))
    mul.aset <- which(apply(mul.C, 1, sum) != 0)
    mul.inaset <- setdiff(1:p, mul.aset)
    mul.s <- length(mul.aset)
    mul.Est <- sum((mul.C - C[, -1])^2)/(p*m)
    mul.Pred <- f(Y = testdata$Y, X = testdata$X, C = mul.C)
    mul.Sen <- length(intersect(mul.aset, P.set))/length(P.set)
    mul.Spe <- length(intersect(mul.inaset, N.set))/length(N.set)
    t2 <- Sys.time()
    if(is.null(mul.Pred)){}else{
      if(mul.Pred == Inf){}else{
        Pred[4] <- mul.Pred
        Est[4] <- mul.Est
        s.hat[4] <- mul.s
        Sen[4] <- mul.Sen
        Spe[4] <- mul.Spe
        time[4] <- difftime(t2,t1,units="s")
        }
      }
  },silent=TRUE)
  
  print("rrvglm begins")
  if(p==50){
    t1 <- Sys.time()
    C.vglm <- array(NA, dim = c(r.max, p, m))
    err.vglm <- rep(Inf, r.max)
    for(i in 1:r.max){
      print(i)
      try({
        fit.vglm <- rrvglm(Y ~ X, multinomial(refLevel = 1), Rank = i)
        C.vglm[i,,] <- coef(fit.vglm, matrix = TRUE)[-1,]
        err.vglm[i] <- f(Y = valdata$Y, X = valdata$X, C = C.vglm[i,,])
      },silent=TRUE)
    }
    
    try({
      id.vglm <- which.min(err.vglm)
      vglm.C <- C.vglm[id.vglm,,]
      vglm.r <- id.vglm
      vglm.aset <- which(apply(vglm.C, 1, sum) != 0)
      vglm.inaset <- setdiff(1:p, vglm.aset)
      vglm.s <- length(vglm.aset)
      vglm.Est <- sum((vglm.C - C[, -1])^2)/(p*m)
      vglm.Pred <- f(Y = testdata$Y, X = testdata$X, C = vglm.C)
      vglm.Sen <- length(intersect(vglm.aset, P.set))/length(P.set)
      vglm.Spe <- length(intersect(vglm.inaset, N.set))/length(N.set)
      t2 <- Sys.time()
      if(is.null(vglm.Pred)){}else{
        if(vglm.Pred == Inf){}else{
          Pred[5] <- vglm.Pred
          Est[5] <- vglm.Est
          r.hat[5] <- vglm.r
          s.hat[5] <- vglm.s
          Sen[5] <- vglm.Sen
          Spe[5] <- vglm.Spe
          time[5] <- difftime(t2,t1,units="s")
          }
        }
    },silent=TRUE)
  }
  return(cbind(Pred, Est, r.hat, s.hat, Sen, Spe, time))
}

get_table <- function(obj,n.rep){
  set.bg <- set.mul <- set.vglm <- NULL
  tab.hat <- array(0, dim = c(n.rep, 5, 7))
  for(i in 1:n.rep){
    tab <- obj[[i]]
    tab[2:4,3] <- 0
    if(any(is.na(tab[1:2,])) == F){set.bg <- c(set.bg, i)}
    if(any(is.na(tab[4,])) == F){set.mul <- c(set.mul, i)}
    if(any(is.na(tab[5,])) == F){set.vglm <- c(set.vglm, i)}
    tab.hat[i,,] <- tab
  }
  ave.bg <- apply(tab.hat[set.bg,1:2,], c(2,3), mean)
  ave.mul <- apply(tab.hat[set.mul,4,], 2, mean)
  ave.vglm <- apply(tab.hat[set.vglm,5,], 2, mean)
  est.ave <- rbind(ave.bg, ave.mul, ave.vglm)
  est.ave[2:3,3] <- NA
  
  var.bg <- apply(tab.hat[set.bg,1:2,], c(2,3), var)
  var.mul <- apply(tab.hat[set.mul,4,], 2, var)
  var.vglm <- apply(tab.hat[set.vglm,5,], 2, var)
  est.var <- rbind(var.bg, var.mul, var.vglm)
  est.var[2:3,3] <- NA
  colnames(est.ave) <- colnames(est.var) <- c("Pred", "Est", "r.hat", "s.hat", "Sen", "Spe", "Time")
  rownames(est.ave) <- rownames(est.var) <- c("SRRMLR", "GLMNET(mlasso)", "NNET", "RRVGLM")
  list(est.ave=est.ave,est.var=est.var)
}

