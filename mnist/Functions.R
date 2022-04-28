get_table <- function(obj,n.rep){
  set.bg <- set.mul <- set.vglm <- NULL
  tab.hat <- array(0, dim = c(n.rep, 4, 3))
  for(i in 1:n.rep){
    tab <- obj[[i]]
    tab[2:3,2] <- 0
    if(any(is.na(tab[1:2,])) == F){set.bg <- c(set.bg, i)}
    if(any(is.na(tab[3,])) == F){set.mul <- c(set.mul, i)}
    if(any(is.na(tab[4,])) == F){set.vglm <- c(set.vglm, i)}
    tab.hat[i,,] <- tab
  }
  ave.bg <- apply(tab.hat[set.bg,1:2,], c(2,3), mean)
  ave.mul <- apply(tab.hat[set.mul,3,], 2, mean)
  ave.vglm <- apply(tab.hat[set.vglm,4,], 2, mean)
  est.ave <- rbind(ave.bg, ave.mul, ave.vglm)
  est.ave[2:3,2] <- NA
  
  var.bg <- apply(tab.hat[set.bg,1:2,], c(2,3), sd)
  var.mul <- apply(tab.hat[set.mul,3,], 2, sd)
  var.vglm <- apply(tab.hat[set.vglm,4,], 2, sd)
  est.var <- rbind(var.bg, var.mul, var.vglm)
  est.var[2:3,2] <- NA
  colnames(est.ave) <- colnames(est.var) <- c("Pred", "r.hat", "s.hat")
  rownames(est.ave) <- rownames(est.var) <- c("SRRMLR", "GLMNET(mlasso)", "NNET", "RRVGLM")
  list(est.ave=est.ave,est.var=est.var)
}

