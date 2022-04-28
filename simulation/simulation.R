source("Function.R")

n.rep <- 100
n <- 400
type = "AR"

cl_num = 20
# Simulation 1 with parallel
p <- 50
cl = makeCluster(cl_num)
registerDoParallel(cl)
fit <- foreach(i = 1:n.rep, .packages = c("MASS","glmnet","nnet","VGAM","SRRMLR")) %dopar% longwait(n = n, p = p, type = type)
stopImplicitCluster()
stopCluster(cl)

save(fit,file="fit50AR.Rdata")
p50 = get_table(obj=fit,n.rep=n.rep)
save(p50,file='p50ar.RData')


# Simulation 2 with parallel
p <- 200
cl = makeCluster(cl_num)
registerDoParallel(cl)
fit <- foreach(i = 1:n.rep, .packages = c("MASS","glmnet","nnet","VGAM","SRRMLR")) %dopar% longwait(n = n, p = p, type = type)
stopImplicitCluster()
stopCluster(cl)

save(fit,file="fit200AR.Rdata")
p200 = get_table(obj=fit,n.rep=n.rep)
save(p200,file='p200ar.RData')

# Simulation 3 with parallel
p <- 500
cl = makeCluster(cl_num)
registerDoParallel(cl)
fit <- foreach(i = 1:n.rep, .packages = c("MASS","glmnet","nnet","VGAM","SRRMLR")) %dopar% longwait(n = n, p = p, type = type)
stopImplicitCluster()
stopCluster(cl)

save(fit,file="fit500AR.Rdata")
p500 = get_table(obj=fit,n.rep=n.rep)
save(p500,file='p500ar.RData')

type = "CS"

# test
# longwait(n = n, p = 50, type = type)
# stop("Done")


# Simulation 1 with parallel
p <- 50
cl = makeCluster(cl_num)
registerDoParallel(cl)
fit <- foreach(i = 1:n.rep, .packages = c("MASS","glmnet","nnet","VGAM","SRRMLR")) %dopar% longwait(n = n, p = p, type = type)
stopImplicitCluster()
stopCluster(cl)

save(fit,file="fit50CS.Rdata")
p50 = get_table(obj=fit,n.rep=n.rep)
save(p50,file='p50cs.RData')


# Simulation 2 with parallel
p <- 200
cl = makeCluster(cl_num)
registerDoParallel(cl)
fit <- foreach(i = 1:n.rep, .packages = c("MASS","glmnet","nnet","VGAM","SRRMLR")) %dopar% longwait(n = n, p = p, type = type)
stopImplicitCluster()
stopCluster(cl)

save(fit,file="fit200CS.Rdata")
p200 = get_table(obj=fit,n.rep=n.rep)
save(p200,file='p200cs.RData')

# Simulation 3 with parallel
p <- 500
cl = makeCluster(cl_num)
registerDoParallel(cl)
fit <- foreach(i = 1:n.rep, .packages = c("MASS","glmnet","nnet","VGAM","SRRMLR")) %dopar% longwait(n = n, p = p, type = type)
stopImplicitCluster()
stopCluster(cl)

save(fit,file="fit500CS.Rdata")
p500 = get_table(obj=fit,n.rep=n.rep)
save(p500,file='p500cs.RData')

stop("Done!")
