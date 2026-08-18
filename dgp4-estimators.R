source("~/these_doctorat/Simulations/fctVARXcoint.R")
setwd("~/these_doctorat/Simulations")

library(parallel)

#bivariate VARX(1,1)  Model, bivariate VAR(1) for Xt

phi1 = matrix(c(0.6,.12,1,.7),2,2)
v0 = matrix(c(2,4,1.5,1.25),2,2)
v1 = matrix(c(-.3,.5,.17,-.44),2,2)
sigma.e = matrix(c(25,5.4,5.4,9.0),2,2)
phix.1 = matrix(c(.8,-.21,-.19,.5),2,2)
sigma.b = matrix(c(4,-.8,-.8,2),2,2)

sim4 = function(n) {
  out = genVARX(n, phi1, list(v0,v1), sigma.e, p=1, s=1, phix.1, sigma.b, pX=1)
  Y = out$y
  X = out$x
  
  W = diff(Y)
  Ylag1 = Y[1:(n-1),]
  Xlag0 = X[2:n,]
  Xlag1 = X[1:(n-1),]
  Z = cbind(Ylag1, Xlag0, Xlag1)
  
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
    V1 = matrix(par[8:11],2,2)
    Epsilon = W - Ylag1 %*% t(AB) - Xlag0 %*% t(V0) - Xlag1 %*% t(V1)
    cov = crossprod(Epsilon) / nprime
    
    log(cov[1,1] * cov[2,2] - cov[1,2] * cov[2,1])
  }
  
  #départ avec les estimés des moindres carrés
  
  start.ls = c(out2$A,out2$B0, para.ls[,-(1:2)])
  
  # estimés du maximum de vraisemblance
  
  out3 = nlminb(start = start.ls, objective = f)
  par.mle = out3$par
  
  # covariance par maximum de vraisemblance
  
  AB_est = cbind(par.mle[1:2], par.mle[1:2] * par.mle[3])
  V0.est = matrix(par.mle[4:7],2,2)
  V1.est = matrix(par.mle[8:11],2,2)
  
  #mêmes donnees, seuls les paramètres changent  
  Epsilon =  W - Ylag1 %*% t(AB_est) - Xlag0 %*% t(V0.est) - Xlag1 %*% t(V1.est)
  cov.mle = crossprod(Epsilon) / nprime
  
  #gives the standardized estimates of B0
  scale_factor = sqrt(sum(Ylag1[, 2]^2))
  true_b0 = -2.5
  stat.lse = scale_factor * (out2$B0 - true_b0)
  stat.mle = scale_factor * (par.mle[3] - true_b0)
  
  return(unname(c(start.ls, cov.ls, par.mle, cov.mle, stat.lse, stat.mle)))
}

cl = makeCluster(detectCores()-1)
clusterExport(cl, varlist = c("phi1", "v0", "v1", "sigma.e", "phix.1", 
                              "sigma.b", "genVARX", "reducedRank", "sim4"))

# It computes the means and sample standard deviations of the LS and ML estimates 
# from all the 10,000 replications 

sink(file = "dgp4.txt", append = F)
bigtic = Sys.time()
for (n in c(50,100,200,400,1000)) {
  clusterExport(cl,"n")
  clusterSetRNGStream(cl, iseed = 1)
  tic = Sys.time()
  nsim = 10000
  res = t(parSapply(cl, 1:nsim, function(i) {sim4(n)}))
  cat("number of observations :", n, "\n")
  tac = Sys.time() - tic
  print(tac); cat("\n")
  
  cat("mean and standard deviation of parameters estimates of VECM\n")
  cat(round(colMeans(res[,1:15]),3),"\n")
  cat(round(colMeans(res[,16:30]),3),"\n\n")
  cat(round(colSD(res[,1:15]),3),"\n")
  cat(round(colSD(res[,16:30]),3),"\n\n")
  
  cat("mean and variance of the standardized estimates of B0\n")
  cat("LSE, MLE\n")
  cat(round(colMeans(res[,31:32]),3),"\n")
  cat(round(cov(res[,31:32]),3),"\n\n")
}
bigtac = Sys.time() - bigtic
print("total time")
print(bigtac)
sink()

#Results of estimates of kappa:

#means
#lse      mle
0.682; 0.419 #n=50
0.128; 0.032 #n=100
0.216; 0.175 #n=200
0.066; 0.045 #n=400
-0.065; -0.075 #n=1000

#variances
152.984; 152.951; 152.951; 153.477 #n=50
111.757; 111.788; 111.788; 111.875 #n=100
96.686; 96.675; 96.675; 96.672 #n=200
88.459; 88.453; 88.453; 88.449 #n=400
85.161; 85.164; 85.164; 85.168 #n=1000


#close the connexions of clusters and free memory
stopCluster(cl)