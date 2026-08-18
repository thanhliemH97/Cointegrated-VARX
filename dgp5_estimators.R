source("~/these_doctorat/Simulations/fctVARXcoint.R")

library(parallel)
cl = makeCluster(detectCores()-1)
clusterEvalQ(cl, library(MTS))

D = diag(c(1, 0.5, 0.4)) #1 unit eigenvalues, others are less than 1

P = matrix(c(.5, 1, 2,
              1, 1, 1,
              2, .5, 0), 3, 3, byrow=T)

(C = P %*% (D-diag(1,3)) %*% solve(P)); -0.25*C[,1]-0.5*C[,2]

det(P) # not = 0
(sumofPhi = P %*% D %*% solve(P))

sumofPhi = round(sumofPhi,3)

eigen(sumofPhi)$values

(Phi1 = matrix(c(0.2, 0.25, -0.2,
                 -0.2, 0.4, 0.1,
                 0.3, -0.4, 0.8),3,3, byrow=T))

(Phi2 = round(sumofPhi - Phi1, 3))

bigPhi = rbind(cbind(Phi1,Phi2),
               cbind(diag(1,3),0,0,0))

round(eigen(bigPhi)$values,4)
round(Mod(eigen(bigPhi)$values),4)

(C = round(sumofPhi - diag(1,3), 3))
round(eigen(C)$values,4) # une valeur propre 0
qr(C)$rank #rang 2 = 3-1

A = C[,1:2]
B = cbind(diag(1,2), matrix(c(-.25,-.5),2,1))
round(A%*%B,3) == C

# parameters in the exogenous part of VARX

V02 = matrix(c(0.4,0.6,0,0.3,0.6,-0.2,0,-0.1,0.7),3,3)
sigma.e2 = matrix(c(25,5,3,5,9,4,3,4,16),3,3)

#covariance matrix of B0
solve(t(A)%*%solve(sigma.e2)%*%A)
#     [,1]     [,2]
#[1,] 53.19967 20.49935
#[2,] 20.49935 18.79869

# VAR(1) parameters of X_t
phix.2 = matrix(c( 0.2,-.6,.1,.4,.3,.7,-.4,.3,.2),3,3) 
sigma.b2 = matrix(c(4,-.8,-.6,-.8, 2,-.1,-.6,-.1,1),3,3)

Mod(eigen(phix.2)$values) #valeurs propres inférieurs à 1

#VECM parameters
cbind(C,-Phi2,V02)

# this function generates VARX data, and gives the LS and ML estimates

sim5 = function(n) {
  out = genVARX(n, list(Phi1, Phi2), V02, sigma.e2, p=2, s=0, phix.2, sigma.b2,
                pX=1)
  Y = out$y
  X = out$x
  W = diff(Y)
  Wlag0 = W[2:(n-1),]; Wlag1 = W[1:(n-2),]
  Ylag1 = Y[2:(n-1),]
  Xlag0 = X[3:n,]
  Z = cbind(Ylag1, Wlag1, Xlag0)
  
  #Estimateurs des moindres carrés de rang plein
  
  para.ls = t(solve(t(Z)%*%Z, t(Z)%*%Wlag0)) 
  resi = Wlag0 - Z %*% t(para.ls)
  nprime = nrow(resi)
  cov.ls = (t(resi) %*% resi) / nrow(Wlag0)
  
  #Estimation de rang réduit
  C.ls = para.ls[,1:3]
  out2 = reducedRank(C.ls, cov.ls, r=2)
  
  #Estimateurs de vraisemblance maximale
  
  f = function(par) {
    
    A1 = par[1:3]
    A2 = par[4:6]
    AB = cbind(A1, A2, A1 * par[7] + A2 * par[8])
    Phi1star = matrix(par[9:17],3,3)
    V0 = matrix(par[18:26],3,3)
    
    Epsilon =  Wlag0 - Ylag1 %*% t(AB) - Wlag1 %*% t(Phi1star) - 
      Xlag0 %*% t(V0)
    
    cov = crossprod(Epsilon) / nprime
    
    loglike = log(det(cov))
    return(loglike)
  }
  
  #départ avec les estimés des moindres carrés
  start.ls = c(out2$A,out2$B0, para.ls[,-(1:3)])
  
  # estimés du maximum de vraisemblance
  out3 = nlminb(start = start.ls, objective = f)
  par.mle = out3$par
  
  # covariance par maximum de vraisemblance
  A1.est = par.mle[1:3]
  A2.est = par.mle[4:6]
  AB.est = cbind(A1.est, A2.est, A1.est * par.mle[7] + A2.est * par.mle[8])
  Phi1star.est = matrix(par.mle[9:17],3,3)
  V0.est = matrix(par.mle[18:26],3,3)
  
  #mêmes donnees, seuls les paramètres changent  
  Epsilon =  Wlag0 - Ylag1 %*% t(AB.est) - Wlag1 %*% t(Phi1star.est) - 
    Xlag0 %*% t(V0.est)
  
  cov.mle = crossprod(Epsilon) / nprime
  
  # vecteur aléatoire du théorème
  B0.lse = c(out2$B0)
  B0.mle = par.mle[7:8]
  B0 = c(-.25,-.5)
  scale_factor = sqrt(sum(Ylag1[,3]^2))
  stat.lse = scale_factor*(B0.lse - B0)
  stat.mle = scale_factor*(B0.mle - B0)
  return(unname(c(start.ls, cov.ls, par.mle, cov.mle, stat.lse, stat.mle)))
}

clusterExport(cl, 
              varlist = c("sim5", "Phi1", "Phi2", "V02", "sigma.e2", 
                          "phix.2", "sigma.b2", "genVARX", "reducedRank"))

clusterSetRNGStream(cl, iseed = 1)
tic = Sys.time()
nsim = 10000
res = t(parSapply(cl, 1:nsim, function(i) {sim5(n)}))
(tac = Sys.time() - tic)
round(colMeans(res),3)
round(colSD(res),3)

bigtic = Sys.time()
for (n in c(50,100,200,400,1000)) {
  clusterExport(cl,"n")
  clusterSetRNGStream(cl, iseed = 1)
  tic = Sys.time()
  res = t(parSapply(cl, 1:10000, function(i) {sim5(n)}))
  cat("number of observations :", n, "\n")
  tac = Sys.time() - tic
  print(tac); cat("\n")
  
  cat("mean and standard deviation of parameters estimates of VECM\n")
  cat(round(colMeans(res[,1:35]),3),"\n")
  cat(round(colMeans(res[,36:70]),3),"\n\n")
  cat(round(colSD(res[,1:35]),3),"\n")
  cat(round(colSD(res[,36:70]),3),"\n\n")
  
  cat("mean and variance of the standardized estimates of B0\n")
  cat("LSE, MLE\n")
  cat(round(colMeans(res[,71:74]),3),"\n")
  cat(round(cov(res[,71:72]),3),"\n\n")
  cat(round(cov(res[,73:74]),3),"\n\n")
}
bigtac = Sys.time() - bigtic
print("total time")
print(bigtac)
