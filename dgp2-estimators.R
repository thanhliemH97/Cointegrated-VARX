source("~/these_doctorat/Simulations/fctVARXcoint.R")
setwd("~/these_doctorat/Simulations")

library(parallel)

#bivariate VARX(1,1)  Model, univariate AR(1) for Xt

phi1 = matrix(c(0.6,.12,1,.7),2,2)
v0 = matrix(c(2,4),2,1)
v1 = matrix(c(-.3,.5),2,1)
sigma.e = matrix(c(25,5.4,5.4,9.0),2,2)
phix.1 = matrix(0.8)
sigma.b = matrix(4)

sim2 = function(n) {
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
  cov.ls = (t(resi) %*% resi) / nprime
  
  #Estimation de rang réduit
  
  C.ls = para.ls[,1:2]
  out2 = reducedRank(C.ls, cov.ls)
  
  #Estimateurs de vraisemblance maximale
  
  f = function(par) {
    
    AB = cbind(par[1:2], par[1:2] * par[3])
    Epsilon = W - Ylag1 %*% t(AB) - Xlag0 %*% t(par[4:5]) - Xlag1 %*% t(par[6:7])
    cov = crossprod(Epsilon) / nprime
    
    log(cov[1,1] * cov[2,2] - cov[1,2] * cov[2,1])
  }
  
  #départ avec les estimés des moindres carrés
  
  start.ls = c(out2$A,out2$B0, para.ls[,-(1:2)])
  
  # estimés du maximum de vraisemblance
  
  out3 = nlminb(start = start.ls, objective = f)
  par.mle = out3$par

  AB_est = cbind(par.mle[1:2], par.mle[1:2] * par.mle[3])
  Epsilon = W - Ylag1 %*% t(AB_est) - Xlag0 %*% t(par.mle[4:5]) -
    Xlag1 %*% t(par.mle[6:7])
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
                              "sigma.b", "genVARX", "reducedRank", "sim2"))

# It computes the means and sample standard deviations of the LS and ML estimates 
# from all the 10,000 replications 
#sink(file = "dgp2.txt", append = F)
bigtic = Sys.time()
for (n in c(50,100,200,400,1000)) {
  clusterExport(cl,"n")
  clusterSetRNGStream(cl, iseed = 1)
  tic = Sys.time()

  res = t(parSapply(cl, 1:10000, function(i) {sim2(n)}))
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
#sink()

#close the connexions of clusters and free memory
stopCluster(cl)

#Results of estimates of kappa:

#means
0.490; 0.271 #n=50
0.191; 0.097 #n=100
0.222; 0.178 #n=200
-0.095; -0.116 #n=400
-0.039; -0.048 #n=1000
#variance
114.427; 114.671; 114.671; 115.277  #n=50
98.859; 98.878; 98.878; 98.946   #n=100
92.116; 92.133; 92.133; 92.160   #n=200
85.447; 85.448; 85.448; 85.452    #n=400
84.373; 84.376; 84.376; 84.380  #1000