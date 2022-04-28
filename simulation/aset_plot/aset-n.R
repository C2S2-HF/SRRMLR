library(MASS)
library(ggplot2)
library(ggpubr)
library(latex2exp)
library(Rcpp)
library(RcppArmadillo)
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

grph <- function(type, p, s, ylabel=FALSE){
  q <- 8
  r <- 5
  true_s = 0:(s-1)
  n = 400:800
  ActiveSet = c()
  for(i in n){
    print(i)
    data <- gendata(n = i, p = p, q = q, s = s, r = r, type = type)
    best <- srrmlr_once(X = data$X, Y = data$Y, r = r, s = s)
    ActiveSet <- c(ActiveSet, length(intersect(true_s, best$Aset[,best$iters])))
  }
  df <- data.frame(cbind(n,ActiveSet))
  p <- ggplot(data=df, aes(x=n,y=ActiveSet))+geom_point()+
    geom_smooth(method="loess",color="black",fill="gray")+
    coord_cartesian(xlim = c(min(n), max(n)), ylim = c(s-5, s))+
    ggtitle(paste0("p:", p, "  mode:", type))+
    theme(axis.text.x = element_text(size = 16,color="black"),
          axis.text.y = element_text(size = 16,color="black"),
          plot.title = element_text(hjust = 0.5))
  if(ylabel){
    p <- p + ylab(expression(paste("intersect(",hat(s),", s*)")))
  }else{
    p <- p + ylab(NULL)
  }
  p
}

set.seed(67)
s = 5
s = 10
p1 <- grph(type="CS", p=50, s=s, ylabel=TRUE)
p2 <- grph(type="CS", p=200, s=s)
p3 <- grph(type="CS", p=500, s=s)
p4 <- grph(type="AR", p=50, s=s, ylabel=TRUE)
p5 <- grph(type="AR", p=200, s=s)
p6 <- grph(type="AR", p=500, s=s)
plot = ggarrange(p1, p2, p3, p4, p5, p6, ncol=3, nrow=2)
ggsave(plot, filename="csar5_3.png", width=12, height=8)
ggsave(plot, filename="csar10_3.png", width=12, height=8)
