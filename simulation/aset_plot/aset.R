library(MASS)
library(nnet)
library(ggplot2)
library(ggpubr)
library(Rcpp)
library(RcppArmadillo)
setwd("~/Desktop/SRRMLR/simulation/aset_plot/")
sourceCpp("srrmlr.cpp")

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

grph <- function(n, s, type, ylabel=FALSE){
  p <- 50
  q <- 8
  r <- 5
  list_ <- matrix(0, 100, 1000)
  for(i_ in 1:100){
    print(i_)
    data <- gendata(n = n, p = p, q = q, s = s, r = r, type = type)
    best <- srrmlr_once(X = data$X, Y = data$Y, s = s, r = r)
    obj <- best$Aset[,1:best$iters]
    for(i in 1:best$iters){
      list_[i_,i] <- length(intersect(best$Aset[, i], 0:(s-1)))
    }
  }
  for(i in 1:100){
    for(j in 2:200){
      if(list_[i,j]==0){
        list_[i,j] = list_[i,j-1]
      }
    }
  }
  lis <- apply(list_[,1:200],2,mean)
  Iteration <- 1:200
  ActiveSet <- lis
  df <- data.frame(cbind(Iteration,ActiveSet))
  plt <- ggplot(data=df, aes(x=Iteration,y=ActiveSet))+geom_point()+
                geom_smooth(method="loess",color="black",fill="gray")+
                ggtitle(paste0("n:",n,"  mode:",type))+
                theme(axis.text.x = element_text(size = 16,color="black"),
                      axis.text.y = element_text(size = 16,color="black"),
                      plot.title = element_text(hjust = 0.5))+ylim(s-4,s)
  if(ylabel){
    plt <- plt + ylab(expression(paste("intersect(",hat(s),", s*)")))
  }else{
    plt <- plt + ylab(NULL)
  }
  plt
}

p1 <- grph(n=200, s=10, type="CS", ylabel=TRUE)
p2 <- grph(n=400, s=10, type="CS")
p3 <- grph(n=800, s=10, type="CS")
p4 <- grph(n=200, s=10, type="AR", ylabel=TRUE)
p5 <- grph(n=400, s=10, type="AR")
p6 <- grph(n=800, s=10, type="AR")
plot = ggarrange(p1, p2, p3, p4, p5, p6, ncol=3, nrow=2)
ggsave(plot, filename="s10.png", width=12, height=8)



