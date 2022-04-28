library(MASS)
library(ggplot2)
library(Rcpp)
library(ggpubr)
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

grph <- function(n, s, type, ylabel=FALSE){
  q <- 8
  r <- 5
  p <- 50
  df = matrix(0, 100, 1000)
  for(i in 1:100){
    print(i)
    data <- gendata(n = n, p = p, q = q, s = s, r = r, type = type)
    temp <- srrmlr_one(X = data$X, Y = data$Y, r = r, s = s)
    for(j in 1:temp$iters){
      df[i, j] = temp$ilp[j]
    }
  }
  df_ = matrix(0, 20000, 2)
  for(i in 1:200){
    for(j in 1:100){
      df_[100*(i-1)+j, 1] = df[j,i]
      df_[100*(i-1)+j, 2] = i
    }
  }
  df_ = data.frame(df_)
  df_$X2 = as.factor(df_$X2)
  
  df.m <- c()
  df.s <- c()
  for(i in unique(df_$X2)){
    df.m <- c(df.m, mean(df_[df_$X2==i,]$X1))
    df.s <- c(df.s, var(df_[df_$X2==i,]$X1))
  }
  dat<-data.frame(iter=unique(df_$X2), df.m=df.m, df.s=df.s)
  
  plt <- ggplot(dat, aes(x=iter, y = df.m)) +
    geom_point(position = position_dodge(0)) +
    geom_errorbar(aes(ymin = df.m - df.s, ymax = df.m + df.s), 
                  width = 0.8, position = position_dodge(0))+
    ggtitle(paste0("n:",n,"  mode:",type))+
    labs(x='outer iteration (m)') +
    scale_x_discrete(breaks=seq(0,190,30))+
    theme(axis.text.x = element_text(size = 16,color="black"),
          axis.text.y = element_text(size = 16,color="black"),
          plot.title = element_text(hjust = 0.5))
  if(ylabel){
    plt <- plt + ylab("inner iteration (k)")
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
ggsave(plot, filename="ioloop10.png", width=12, height=8)



type = "CS"
n <- 400
p <- 50
for(i in 1:100){
  data <- gendata(n = n, p = p, q = q, s = s, r = r, type = type)
  temp <- srrmlr_one(X = data$X, Y = data$Y, r = r, s = s)
  for(j in 1:temp$iters){
    cat(temp$msel[j], " ", file='conv_cs50.txt', append=TRUE)
  }
  cat("\n", file='conv_cs50.txt', append=TRUE)
}
p <- 200
for(i in 1:100){
  data <- gendata(n = n, p = p, q = q, s = s, r = r, type = type)
  temp <- srrmlr_one(X = data$X, Y = data$Y, r = r, s = s)
  for(j in 1:temp$iters){
    cat(temp$msel[j], " ", file='conv_cs200.txt', append=TRUE)
  }
  cat("\n", file='conv_cs200.txt', append=TRUE)
}
p <- 500
for(i in 1:100){
  data <- gendata(n = n, p = p, q = q, s = s, r = r, type = type)
  temp <- srrmlr_one(X = data$X, Y = data$Y, r = r, s = s)
  for(j in 1:temp$iters){
    cat(temp$msel[j], " ", file='conv_cs500.txt', append=TRUE)
  }
  cat("\n", file='conv_cs500.txt', append=TRUE)
}


type = "AR"
p <- 50
for(i in 1:100){
  data <- gendata(n = n, p = p, q = q, s = s, r = r, type = type)
  temp <- srrmlr_one(X = data$X, Y = data$Y, r = r, s = s)
  for(j in 1:temp$iters){
    cat(temp$msel[j], " ", file='conv_ar50.txt', append=TRUE)
  }
  cat("\n", file='conv_ar50.txt', append=TRUE)
}
p <- 200
for(i in 1:100){
  data <- gendata(n = n, p = p, q = q, s = s, r = r, type = type)
  temp <- srrmlr_one(X = data$X, Y = data$Y, r = r, s = s)
  for(j in 1:temp$iters){
    cat(temp$msel[j], " ", file='conv_ar200.txt', append=TRUE)
  }
  cat("\n", file='conv_ar200.txt', append=TRUE)
}
p <- 500
for(i in 1:100){
  data <- gendata(n = n, p = p, q = q, s = s, r = r, type = type)
  temp <- srrmlr_one(X = data$X, Y = data$Y, r = r, s = s)
  for(j in 1:temp$iters){
    cat(temp$msel[j], " ", file='conv_ar500.txt', append=TRUE)
  }
  cat("\n", file='conv_ar500.txt', append=TRUE)
}
