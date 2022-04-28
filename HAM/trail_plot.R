library(ggplot2)
library(reshape2)

get_variables <- function(fit){
  pred <- s <- c()
  for(i in fit){
    pred <- rbind(pred, i[c(1,2,3),1])
    s <- rbind(s, i[c(1,2,3),3])
  }
  colnames(pred) <- colnames(s) <- c("SRRMLR","mLasso","NNET")
  return(list(pred=melt(pred,varnames=c("Trail","Method")),s=melt(s,varnames=c("Trail","Method"))))
}

getted <- get_variables(fit)

data = na.omit(rbind(cbind(getted$pred,title='Error'),cbind(getted$s,title='Sparsity')))
ggplot(data=data,aes(x=Trail,y=value,color=Method,shape=Method))+geom_point()+geom_line()+
  facet_wrap(~title,scales='free')+xlab("Trail")+ylab(NULL)

