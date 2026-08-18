source("~/these_doctorat/Simulations/fctVARXcoint.R")

library(parallel)

#bivariate VARX(1,0)  Model, univarite VAR(1) for Xt

phi1 = matrix(c(0.6,.12,1,.7),2,2)
v0 = matrix(c(2,4),2,1)
sigma.e = matrix(c(25,5.4,5.4,9.0),2,2)
phix.1 = matrix(0.8)
sigma.b = matrix(4)

#verification of eigen values of companion matrix

#valeur propre 1 <=> racine 1 pour det{\Phi(z)} = 0
#valeur propre 0.3 <=> racine  10/3 pour det{\Phi(z)} = 0

eigen(phi1)$values 

# determining the rank of matrices

qr(phi1)$rank #rang 2
qr(diag(1,2)-phi1)$rank #rang 1

#Checking the convergence of the standardised form of \hat{B}_0 and \tilde{B}_0 

#Its asymptotic variance: 84.472

A = matrix(c(-.4,.12), 2,1)
sigma.e = matrix(c(25,5.4,5.4,9.0), 2,2)
solve(t(A)%*%solve(sigma.e)%*%A) [1]

################################################################################

#This function generates VARX data, and gives the LS and ML estimates

sim1 = function(n) {
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
  cov.ls = (t(resi) %*% resi) / nprime
  
  #Estimation de rang réduit
  
  C.ls = para.ls[,1:2]
  out2 = reducedRank(C.ls, cov.ls)
  
  #Estimateurs de vraisemblance maximale
  
  f = function(par) {
    
    AB = cbind(par[1:2], par[1:2] * par[3])
    Epsilon = W - Ylag1 %*% t(AB) - Xlag0 %*% t(par[4:5])
    cov = crossprod(Epsilon) / nprime
    
    log(cov[1,1] * cov[2,2] - cov[1,2] * cov[2,1])
  }
  
  #départ avec les estimés des moindres carrés
  
  start.ls = c(out2$A,out2$B0, para.ls[,3])
  
  # estimés du maximum de vraisemblance
  
  out3 = nlminb(start = start.ls, objective = f)
  par.mle = out3$par
  
  AB_est = cbind(par.mle[1:2], par.mle[1:2] * par.mle[3])
  Epsilon = W - Ylag1 %*% t(AB_est) - Xlag0 %*% t(par.mle[4:5])
  cov.mle = crossprod(Epsilon) / nprime
  
  #gives the standardized estimates of B0
  scale_factor = sqrt(sum(Ylag1[, 2]^2))
  true_b0 = -2.5
  stat.lse = scale_factor * (out2$B0 - true_b0)
  stat.mle = scale_factor * (par.mle[3] - true_b0)
    
  return(unname(c(start.ls, cov.ls, par.mle, cov.mle, stat.lse, stat.mle)))
}

cl = makeCluster(detectCores()-1)
clusterExport(cl, varlist = c("phi1", "v0", "sigma.e", "phix.1","sigma.b",
                              "genVARX","reducedRank","sim1"))

# It computes the means and sample standard deviations of the LS and ML estimates 
# from all the 10,000 replications 

#sink(file = "dgp2.txt", append = F)
bigtic = Sys.time()
for (n in c(50,100,200,400,1000)) {
  clusterExport(cl,"n")
  clusterSetRNGStream(cl, iseed = 1)
  tic = Sys.time()
  # change n for 50, 100, 200, 400
  res = t(parSapply(cl, 1:10000, function(i) {sim1(n)}))
  cat("number of observations :", n, "\n")
  tac = Sys.time() - tic
  print(tac); cat("\n")
  
  cat("mean and standard deviation of parameters estimates of VECM\n")
  cat(round(colMeans(res[,1:9]),3),"\n")
  cat(round(colMeans(res[,10:18]),3),"\n\n")
  cat(round(colSD(res[,1:9]),3),"\n")
  cat(round(colSD(res[,10:18]),3),"\n\n")
  
  cat("mean and variance of the standardized estimates of B0\n")
  cat("LSE, MLE\n")
  cat(round(colMeans(res[,19:20]),3),"\n")
  cat(round(diag(Var(res[,19:20]),3)),"\n\n")
}
bigtac = Sys.time() - bigtic
print("total time")
print(bigtac)
#sink()

#Results of estimates of kappa:

#means
#lse      mle
0.364;  0.212 #n=50
0.129;  0.064 #n=100
0.182;  0.152 #n=200
-0.046; -0.066 #n=300
-0.106; -0.121 #n=400
0.096;  0.084 #n=500
0.036;  0.026 #n=600
-0.011; -0.019 #n=700
-0.034; -0.041 #n=800
-0.040; -0.046 #n=900
-0.042; -0.049 #n=1000
0.079;  0.074 #n=1500
0.015;  0.008 #n=2000

#variances
#lse      mle
109.966;110.297 #n=50
97.383; 97.408 #n=100
91.517; 91.532 #n=200
88.866; 88.870 #n=300
84.980; 84.984 #n=400
85.729; 85.727 #n=500
85.250; 85.249 #n=600
83.941; 83.938 #n=700
86.375; 86.375 #n=800
86.418; 86.420 #n=900
84.186; 84.190 #n=1000
84.074; 84.072 #n=1500
83.660; 83.658 #n=2000


#close the connexions of clusters and free memory
stopCluster(cl)
