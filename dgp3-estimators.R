source("~/these_doctorat/Simulations/fctVARXcoint.R")
setwd("~/these_doctorat/Simulations")

library(parallel)

#bivariate VARX(1,0)  Model, bivariate VAR(1) for Xt

phi1 = matrix(c(0.6,.12,1,.7),2,2)
v0 = matrix(c(2,4,1.5,1.25),2,2)
#v1 = matrix(c(-.3,.5,.17,-.44),2,2)
sigma.e = matrix(c(25,5.4,5.4,9.0),2,2)
phix.1 = matrix(c(.8,-.21,-.19,.5),2,2)
sigma.b = matrix(c(4,-.8,-.8,2),2,2)

sim3 = function(n) {
  out = genVARX(n, phi1, v0, sigma.e, p=1, s=0, phix.1, sigma.b, pX=1)
  Y = out$y
  X = out$x
  
  W = diff(Y)
  Ylag1 = Y[1:(n-1),]
  Xlag0 = X[2:n,]
  Z = cbind(Ylag1, Xlag0)
  
  #Estimateurs des moindres carrés de rang plein
  
  para.ls = t(solve(t(Z)%*%Z, t(Z)%*%W)) 
  resi = W - Z %*% t(para.ls)
  nprime = nrow(resi)
  cov.ls = (t(resi) %*% resi) / nrow(resi)
  
  #Estimation de rang réduit
  
  C.ls = para.ls[,1:2]
  out2 = reducedRank(C.ls, cov.ls)
  
  #Estimateurs de vraisemblance maximale
  
  f = function(par) {
    
    AB = cbind(par[1:2], par[1:2] * par[3])
    V0 = matrix(par[4:7],2,2)
    Epsilon = W - Ylag1 %*% t(AB) - Xlag0 %*% t(V0)
    cov = crossprod(Epsilon) / nprime
    
    log(cov[1,1] * cov[2,2] - cov[1,2] * cov[2,1])
  }
  
  #départ avec les estimés des moindres carrés
  
  start.ls = c(out2$A,out2$B0, para.ls[,-(1:2)])
  
  # estimés du maximum de vraisemblance
  
  out3 = nlminb(start = start.ls, objective = f)
  par.mle = out3$par
  
  AB_est = cbind(par.mle[1:2], par.mle[1:2] * par.mle[3])
  V0_est = matrix(par.mle[4:7],2,2)
  Epsilon = W - Ylag1 %*% t(AB_est) - Xlag0 %*% t(V0_est)
  cov.mle = crossprod(Epsilon) / nprime
 
  #gives the standardized estimates of B0
  scale_factor = sqrt(sum(Ylag1[, 2]^2))
  true_b0 = -2.5
  stat.lse = scale_factor * (out2$B0 - true_b0)
  stat.mle = scale_factor * (par.mle[3] - true_b0)
  
  return(unname(c(start.ls, cov.ls, par.mle, cov.mle, stat.lse, stat.mle)))
}

cl = makeCluster(detectCores()-1)
clusterExport(cl, varlist = c("phi1", "v0", "sigma.e", "phix.1", 
                              "sigma.b", "genVARX", "reducedRank", "sim3"))

# It computes the means and sample standard deviations of the LS and ML estimates 
# from all the 10,000 replications 

sink(file = "dgp3.txt", append = F)
bigtic = Sys.time()
for (n in c(50,100,200,400,1000)) {
  clusterExport(cl,"n")
  clusterSetRNGStream(cl, iseed = 1)
  tic = Sys.time()
  nsim = 10000
  res = t(parSapply(cl, 1:nsim, function(i) {sim3(n)}))
  cat("number of observations :", n, "\n")
  tac = Sys.time() - tic
  print(tac); cat("\n")
  
  cat("mean and standard deviation of parameters estimates of VECM\n")
  cat(round(colMeans(res[,1:11]),3),"\n")
  cat(round(colMeans(res[,12:22]),3),"\n\n")
  cat(round(colSD(res[,1:11]),3),"\n")
  cat(round(colSD(res[,12:22]),3),"\n\n")
  
  cat("mean and variance of the standardized estimates of B0\n")
  cat("LSE, MLE\n")
  cat(round(colMeans(res[,23:24]),3),"\n")
  cat(round(var(res[,23:24]),3),"\n\n")
}
bigtac = Sys.time() - bigtic
print("total time")
print(bigtac)
sink()

#close the connexions of clusters and free memory
stopCluster(cl)

#Results for standardises B0 estimators

#means
0.455; 0.284 #n=50
0.104; 0.039 #n=100
0.181; 0.153 #n=200
0.046; 0.032 #n=400
-0.067; -0.074 #n=1000
#variance
141.091; 140.933; 140.933; 140.968 #n=50
108.327; 108.338; 108.338; 108.372 #n=100
95.014; 95.009; 95.009; 95.007   #n=200
87.870; 87.866; 87.866; 87.863   #n=400
85.174; 85.175; 85.175; 85.177   #n=1000