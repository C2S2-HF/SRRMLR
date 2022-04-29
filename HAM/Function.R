## Compute negative log likelihood value
f <- function(Y, X, C) {
  z <- X%*%C
  z <- exp(z)
  z <- log(apply(z,1,sum)+1)
  loglike <- sum(diag(X%*%C%*%t(Y[,-1])))-sum(z)
  return(-loglike)
}


get_table <- function(obj, n.rep){
  set.bg <- set.mul <- set.vglm <- NULL
  tab.hat <- array(0, dim = c(n.rep, 3, 3))
  for(i in 1:n.rep){
    tab <- obj[[i]]
    tab[2:3,2] <- 0
    if(any(is.na(tab[1:2,])) == F){set.bg <- c(set.bg, i)}
    if(any(is.na(tab[3,])) == F){set.mul <- c(set.mul, i)}
    tab.hat[i,,] = tab
  }
  ave.bg <- apply(tab.hat[set.bg,1:2,], c(2,3), mean)
  var.bg <- apply(tab.hat[set.bg,1:2,], c(2,3), sd)
  if(length(set.mul)==1){
    ave.mul <- tab.hat[set.mul,3,]
    var.mul <- 0
  }else{
    ave.mul <- apply(tab.hat[set.mul,3,], 2, mean)
    var.mul <- apply(tab.hat[set.mul,3,], 2, sd)
  }
  est.ave <- rbind(ave.bg, ave.mul)
  est.ave[2:3,2] <- NA
  est.var <- rbind(var.bg, var.mul)
  est.var[2:3,2] <- NA
  colnames(est.ave) <- colnames(est.var) <- c("Pred", "r.hat", "s.hat")
  rownames(est.ave) <- rownames(est.var) <- c("SRRMLR", "mlasso", "NNET")
  list(rea=est.ave, rev=est.var)
}

