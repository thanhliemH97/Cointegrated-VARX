source("~/these_doctorat/Simulations/fctVARXcoint.R")

library(parallel)
cl = makeCluster(detectCores()-1)
clusterEvalQ(cl, {library(MASS) 
  library(MTS)})

sim3.1 = function(n, eigenvalues, r) {
  
  ## --- Constantes fixes ---
  Phi1star = matrix(c(.15,.19,.11,
                      -.2,-.08,-.05,
                      .45,.30,.32), 3,3)
  SigmaEps = matrix(c(.47,.2,.18,
                      .2,.32,.27,
                      .18,.27,.30), 3,3)
  phix.1 = 0.8
  sigma.b = 4
  V0 = matrix(c(.08,.13,-.05),3,1)
  Q = matrix(c(-.29,-.01,-.75,
               -.47,-.85,1.39,
               -.57,1,-.55), 3,3)
  P = solve(Q)
  
  ## --- Matrices autorégressives à varier ---
  
  D = diag(eigenvalues,3)
  C = P %*% D %*% Q - diag(1,3)
  Phi2 = -Phi1star
  Phi1 = C + diag(1,3) + Phi1star
  
  ## --- Simulation VARX(2,0) ---
  out = genVARX(n, list(Phi1,Phi2), V0, SigmaEps, p=2, s=0, phix.1, sigma.b, pX=1)
  Y = out$y
  X = out$x
  
  ## --- Construction des variables retardées ---
  W = diff(Y)
  Ylag1 = Y[2:(n-1),]
  Wlag0 = W[2:(n-1),]
  Xlag0 = X[3:n,]
  nobs = nrow(Wlag0)
  Z = cbind(Ylag1, Xlag0)
  
  ## --- Estimation OLS ---
  out2 = ols(Z, Wlag0)
  para.ls = t(out2$parameters)
  cov.ls = out2$covariance
  C.ls = para.ls[,1:3]
  Theta2.ls = para.ls[,-(1:3)]
  
  ## --- Point de départ pour la MV ---
  if (r > 0) {
    out3 = reducedRank(C.ls, cov.ls, r)
    start.ls = c(out3$A, out3$B0, Theta2.ls)
  } else {
    start.ls = c(Theta2.ls)
  }
  
  ## --- Fonction de log-vraisemblance modifiée à optimiser ---
  loglikelihood = function(para) {
    if (r == 2) {
      A = matrix(para[1:6], 3, 2)
      B = matrix(cbind(diag(1,2), para[7:8]), 2, 3)
      V0 = matrix(para[9:11],3,1)
      terme.AB = Ylag1 %*% t(A %*% B)
    } else if (r == 1) {
      A = matrix(para[1:3], 3, 1)
      B = matrix(c(1, para[4:5]), 1, 3)
      V0 = matrix(para[6:8],3,1)
      terme.AB = Ylag1 %*% t(A %*% B)
    } else {
      V0 = matrix(para[1:3], 3, 1)
      terme.AB = 0
    }
    
    Eps = Wlag0 - terme.AB - Xlag0 %*% t(V0)
    cov = crossprod(Eps) / nobs
    log(det(cov))
  }
  
  # --- Minimisation ---
  opt = nlminb(start = start.ls, objective = loglikelihood)
  para.mle = opt$par
  
  # --- Retour : statistique du test ---
  nobs * (loglikelihood(para.mle) - log(det(cov.ls)))
}
sim3.2 = function(n, eigenvalues, r) {
  
  ## --- Constantes fixes ---
  Phi1star = matrix(c(.15,.19,.11,
                      -.2,-.08,-.05,
                      .45,.30,.32), 3,3)
  SigmaEps = matrix(c(.47,.2,.18,
                      .2,.32,.27,
                      .18,.27,.30), 3,3)
  phix.1 = 0.8
  sigma.b = 4
  V0 = matrix(c(.08,.13,-.05),3,1)
  Q = matrix(c(-.29,-.01,-.75,
               -.47,-.85,1.39,
               -.57,1,-.55), 3,3)
  P = solve(Q)
  
  ## --- Matrices autorégressives à varier ---
  
  D = diag(eigenvalues,3)
  C = P %*% D %*% Q - diag(1,3)
  Phi2 = -Phi1star
  Phi1 = C + diag(1,3) + Phi1star
  
  ## --- Simulation VARX(2,0) ---
  out = genVARX(n, list(Phi1,Phi2), V0, SigmaEps, p=2, s=0, phix.1, sigma.b, pX=1)
  Y = out$y
  X = out$x
  
  ## --- Construction des variables retardées ---
  W = diff(Y)
  Ylag1 = Y[2:(n-1),]
  Wlag0 = W[2:(n-1),]
  Wlag1 = W[1:(n-2),]
  Xlag0 = X[3:n,]
  nobs = nrow(Wlag0)
  
  Z = cbind(Ylag1, Wlag1, Xlag0)
  
  ## --- Estimation OLS ---
  out2 = ols(Z, Wlag0)
  para.ls = t(out2$parameters)
  cov.ls = out2$covariance
  C.ls = para.ls[,1:3]
  Theta2.ls = para.ls[,-(1:3)]
  
  ## --- Point de départ pour la MV ---
  if (r > 0) {
    out3 = reducedRank(C.ls, cov.ls, r)
    start.ls = c(out3$A, out3$B0, Theta2.ls)
  } else {
    start.ls = c(Theta2.ls)
  }
  
  ## --- Fonction de log-vraisemblance à optimiser ---
  loglikelihood =  function(para) {
    if (r == 2) {
      A = matrix(para[1:6], 3, 2)
      B = matrix(cbind(diag(1,2), para[7:8]), 2, 3)
      Phi1s = matrix(para[9:17], 3, 3)
      V0 = matrix(para[18:20], 3, 1)
      AB = A %*% B
      terme.AB = Ylag1 %*% t(AB)
    } else if (r == 1) {
      A = matrix(para[1:3], 3, 1)
      B = matrix(c(1, para[4:5]), 1, 3)
      Phi1s = matrix(para[6:14], 3, 3)
      V0 = matrix(para[15:17], 3, 1)
      terme.AB = Ylag1 %*% t(A %*% B)
    } else {
      Phi1s = matrix(para[1:9], 3, 3)
      V0 = matrix(para[10:12], 3, 1)
      terme.AB = 0
    }
    Eps = Wlag0 - terme.AB - Wlag1 %*% t(Phi1s) - Xlag0 %*% t(V0)
    cov = crossprod(Eps) / nobs
    log(det(cov))
  }
  
  ## --- Maximisation ---
  opt = nlminb(start = start.ls, objective = loglikelihood)
  para.mle = opt$par
  
  ## --- Retour : statistique du test ---
  nobs * (loglikelihood(para.mle) - log(det(cov.ls)))
}
sim3.3 = function(n, eigenvalues, r) {
  
  ## --- Constantes fixes ---
  Phi1star = matrix(c(.15,.19,.11,
                      -.2,-.08,-.05,
                      .45,.30,.32), 3,3)
  SigmaEps = matrix(c(.47,.2,.18,
                      .2,.32,.27,
                      .18,.27,.30), 3,3)
  phix.1 = 0.8
  sigma.b = 4
  V0 = matrix(c(.08,.13,-.05),3,1)
  Q = matrix(c(-.29,-.01,-.75,
               -.47,-.85,1.39,
               -.57,1,-.55), 3,3)
  P = solve(Q)
  
  ## --- Matrices autorégressives à varier ---
  
  D = diag(eigenvalues,3)
  C = P %*% D %*% Q - diag(1,3)
  Phi2 = -Phi1star
  Phi1 = C + diag(1,3) + Phi1star
  
  ## --- Simulation VARX(2,0) ---
  out = genVARX(n, list(Phi1,Phi2), V0, SigmaEps, p=2, s=0, phix.1, sigma.b, pX=1)
  Y = out$y
  X = out$x
  
  ## --- Construction des variables retardées ---
  W = diff(Y)
  Ylag1 = Y[3:(n-1),]
  Wlag0 = W[3:(n-1),] 
  Wlag1 = W[2:(n-2),]
  Wlag2 = W[1:(n-3),]
  Xlag0 = X[4:n,]
  
  nobs = nrow(Wlag0)
  Z = cbind(Ylag1, Wlag1, Wlag2, Xlag0)
  
  ## --- Estimation OLS ---
  out2 = ols(Z, Wlag0)
  para.ls = t(out2$parameters)
  cov.ls = out2$covariance
  C.ls = para.ls[,1:3]
  Theta2.ls = para.ls[,-(1:3)]
  
  ## --- Point de départ pour la MV ---
  if (r > 0) {
    out3 = reducedRank(C.ls, cov.ls, r)
    start.ls = c(out3$A, out3$B0, Theta2.ls)
  } else {
    start.ls = c(Theta2.ls)
  }
  
  ## --- Fonction de log-vraisemblance à optimiser ---
  loglikelihood = function(para) {
    if (r == 2) {
      A = matrix(para[1:6], 3, 2)
      B = matrix(cbind(diag(1,2), para[7:8]), 2, 3)
      Phi1star = matrix(para[9:17], 3, 3)
      Phi2star = matrix(para[18:26], 3, 3)
      V0 = matrix(para[27:29], 3, 1)
      AB = A %*% B
      terme.AB = Ylag1 %*% t(AB)
    } else if (r == 1) {
      A = matrix(para[1:3], 3, 1)
      B = matrix(c(1,para[4:5]), 1, 3)
      Phi1star = matrix(para[6:14], 3, 3)
      Phi2star = matrix(para[15:23], 3, 3)
      V0 = matrix(para[24:26], 3, 1)
      terme.AB = Ylag1 %*% t(A %*% B)
    } else {
      Phi1star = matrix(para[1:9], 3, 3)
      Phi2star = matrix(para[10:18], 3, 3)
      V0 = matrix(para[19:21], 3, 1)
      terme.AB = 0
    }
    
    Eps = Wlag0 - terme.AB - Wlag1 %*% t(Phi1star) - Wlag2 %*% t(Phi2star) -
      Xlag0 %*% t(V0)
    cov = crossprod(Eps) / nobs
    log(det(cov))
  }
  
  ## --- Maximisation ---
  opt = nlminb(start = start.ls, objective = loglikelihood)
  para.mle = opt$par
  
  ## --- Retour : statistique du test ---
  nobs * (loglikelihood(para.mle) - log(det(cov.ls)))
}
sim3.4 = function(n, eigenvalues, r) {
  
  ## --- Constantes fixes ---
  Phi1star = matrix(c(.15,.19,.11,
                      -.2,-.08,-.05,
                      .45,.30,.32), 3,3)
  SigmaEps = matrix(c(.47,.2,.18,
                      .2,.32,.27,
                      .18,.27,.30), 3,3)
  phix.1 = 0.8
  sigma.b = 4
  V0 = matrix(c(.08,.13,-.05),3,1)
  Q = matrix(c(-.29,-.01,-.75,
               -.47,-.85,1.39,
               -.57,1,-.55), 3,3)
  P = solve(Q)
  
  ## --- Matrices autorégressives à varier ---
  D = diag(eigenvalues,3)
  C = P %*% D %*% Q - diag(1,3)
  Phi2 = -Phi1star
  Phi1 = C + diag(1,3) + Phi1star
  
  ## --- Simulation VARX(2,0) ---
  out = genVARX(n, list(Phi1,Phi2), V0, SigmaEps, p=2, s=0, phix.1, sigma.b, pX=1)
  Y = out$y
  X = out$x
  
  ## --- Construction des variables retardées ---
  W = diff(Y)
  Ylag1 = Y[4:(n-1),]
  Wlag0 = W[4:(n-1),] 
  Wlag1 = W[3:(n-2),]
  Wlag2 = W[2:(n-3),]
  Wlag3 = W[1:(n-4),]
  Xlag0 = X[5:n,]
  
  nobs = nrow(Wlag0)
  Z = cbind(Ylag1, Wlag1, Wlag2, Wlag3, Xlag0)
  
  ## --- Estimation MCO ---
  out2 = ols(Z, Wlag0)
  para.ls = t(out2$parameters)
  cov.ls = out2$covariance
  C.ls = para.ls[,1:3]
  Theta2.ls = para.ls[,-(1:3)]
  
  ## --- Point de départ pour la MV ---
  if (r > 0) {
    out3 = reducedRank(C.ls, cov.ls, r)
    start.ls = c(out3$A, out3$B0, Theta2.ls)
  } else {
    start.ls = c(Theta2.ls)
  }
  
  ## --- Fonction de log-vraisemblance à optimiser ---
  loglikelihood = function(para) {
    if (r == 2) {
      A = matrix(para[1:6], 3, 2)
      B = matrix(cbind(diag(1,2), para[7:8]), 2, 3)
      Phi1star = matrix(para[9:17], 3, 3)
      Phi2star = matrix(para[18:26], 3, 3)
      Phi3star = matrix(para[27:35], 3, 3)
      V0 = matrix(para[36:38], 3, 1)
      AB = A %*% B
      terme.AB = Ylag1 %*% t(AB)
    } else if (r == 1) {
      A = matrix(para[1:3], 3, 1)
      B = matrix(c(1,para[4:5]), 1, 3)
      Phi1star = matrix(para[6:14], 3, 3)
      Phi2star = matrix(para[15:23], 3, 3)
      Phi3star = matrix(para[24:32], 3, 3)
      V0 = matrix(para[33:35], 3, 1)
      terme.AB = Ylag1 %*% t(A %*% B)
    } else {
      Phi1star = matrix(para[1:9], 3, 3)
      Phi2star = matrix(para[10:18], 3, 3)
      Phi3star = matrix(para[19:27], 3, 3)
      V0 = matrix(para[28:30], 3, 1)
      terme.AB = 0
    }
    
    Eps = Wlag0 - terme.AB - Wlag1 %*% t(Phi1star) - Wlag2 %*% t(Phi2star) -
      Wlag3 %*% t(Phi3star) - Xlag0 %*% t(V0)
    cov = crossprod(Eps) / nobs
    log(det(cov))
  }
  
  ## --- Maximisation ---
  opt = nlminb(start = start.ls, objective = loglikelihood)
  para.mle = opt$par
  
  ## --- Retour : statistique du test ---
  nobs * (loglikelihood(para.mle) - log(det(cov.ls)))
}
sim3.5 = function(n, eigenvalues, r) {
  
  ## --- Constantes fixes ---
  Phi1star = matrix(c(.15,.19,.11,
                      -.2,-.08,-.05,
                      .45,.30,.32), 3,3)
  SigmaEps = matrix(c(.47,.2,.18,
                      .2,.32,.27,
                      .18,.27,.30), 3,3)
  phix.1 = 0.8
  sigma.b = 4
  V0 = matrix(c(.08,.13,-.05),3,1)
  Q = matrix(c(-.29,-.01,-.75,
               -.47,-.85,1.39,
               -.57,1,-.55), 3,3)
  P = solve(Q)
  
  ## --- Matrices autorégressives à varier ---
  
  D = diag(eigenvalues,3)
  C = P %*% D %*% Q - diag(1,3)
  Phi2 = -Phi1star
  Phi1 = C + diag(1,3) + Phi1star
  
  ## --- Simulation VARX(2,0) ---
  out = genVARX(n, list(Phi1,Phi2), V0, SigmaEps, p=2, s=0, phix.1, sigma.b, pX=1)
  Y = out$y
  X = out$x
  
  ## --- Construction des variables retardées ---
  W = diff(Y)
  Ylag1 = Y[2:(n-1),]
  Wlag0 = W[2:(n-1),]
  Wlag1 = W[1:(n-2),]
  nobs = nrow(Wlag0)
  Z = cbind(Ylag1, Wlag1)
  
  ## --- Estimation OLS ---
  out2 = ols(Z, Wlag0)
  para.ls = t(out2$parameters)
  cov.ls = out2$covariance
  C.ls = para.ls[,1:3]
  Theta2.ls = para.ls[,-(1:3)]
  
  ## --- Point de départ pour la MV ---
  if (r > 0) {
    out3 = reducedRank(C.ls, cov.ls, r)
    start.ls = c(out3$A, out3$B0, Theta2.ls)
  } else {
    start.ls = c(Theta2.ls)
  }
  
  ## --- Fonction de log-vraisemblance à optimiser ---
  loglikelihood =  function(para) {
    if (r == 2) {
      A = matrix(para[1:6], 3, 2)
      B = matrix(cbind(diag(1,2), para[7:8]), 2, 3)
      Phi1s = matrix(para[9:17], 3, 3)
      terme.AB = Ylag1 %*% t(A %*% B)
    } else if (r == 1) {
      A = matrix(para[1:3], 3, 1)
      B = matrix(c(1, para[4:5]), 1, 3)
      Phi1s = matrix(para[6:14], 3, 3)
      terme.AB = Ylag1 %*% t(A %*% B)
    } else {
      Phi1s = matrix(para[1:9], 3, 3)
      terme.AB = 0
    }
    Eps = Wlag0 - terme.AB - Wlag1 %*% t(Phi1s)
    cov = crossprod(Eps) / nobs
    log(det(cov))
  }
  
  ## --- Maximisation ---
  opt = nlminb(start = start.ls, objective = loglikelihood)
  para.mle = opt$par
  
  ## --- Retour : statistique du test ---
  nobs * (loglikelihood(para.mle) - log(det(cov.ls)))
}

sim3Parallel = function(nsim, n, eigenvalues, r, p, exo=T) {
  clusterSetRNGStream(cl, iseed = 1)
  clusterExport(cl,varlist=c("n","sim3.1", "sim3.2", "sim3.3", "sim3.4", "sim3.5",
                             "genVARX","ols","reducedRank"))
  tic = Sys.time()
  if (p==1) {
    testStat = parSapply(cl,1:nsim, function(i) sim3.1(n,eigenvalues,r))
  } else if (p==2) {
    if (exo==T) {
      testStat = parSapply(cl,1:nsim, function(i) sim3.2(n,eigenvalues,r))
    } else {
      testStat = parSapply(cl,1:nsim, function(i) sim3.5(n,eigenvalues,r))
    }
    
  } else if (p==3) {
    testStat = parSapply(cl,1:nsim, function(i) sim3.3(n,eigenvalues,r))
  } else if (p==4) {
    testStat = parSapply(cl,1:nsim, function(i) sim3.4(n,eigenvalues,r))
  }
  
  tac = Sys.time() - tic
  print(tac)
  
  if (r==2) {
    crit = c(4.122, 2.974) #case 1(a)
  } else if (r==1) {
    crit = c(12.336, 10.481) #case 2(a)
  } else if (r==0) {
    crit = c(24.29, 21.787) #case 3(a)
  }
  #hist(testStat, seq(min(testStat), max(testStat) + 0.2, 0.2), probability = T,
  #     main = "", xlab = "", ylab = "")
  #print(round(quantile(testStat, c(.005,.01,.025,.05,.1,.5,.9,.95,.975,.99,.995)),3))
  return(c(mean(testStat>crit[1]),mean(testStat>crit[2])))
}

nsim=10000; n=50; p=1

sim3Parallel(nsim,n,c(1,.87,.73),r=2, p)
sim3Parallel(nsim,n,c(1,1,.73),r=1, p)
sim3Parallel(nsim,n,c(1,1,1),r=0, p)
sim3Parallel(nsim,n,c(.96,.87,.73),r=2, p)
sim3Parallel(nsim,n,c(1,.96,.73),r=1, p)
sim3Parallel(nsim,n,c(1,1,.96),r=0, p)
sim3Parallel(nsim,n,c(.92,.87,.73),r=2, p)
sim3Parallel(nsim,n,c(1,.92,.73),r=1, p)
sim3Parallel(nsim,n,c(1,1,.92),r=0, p)

#p = 2 fixed, the VAR(2) is ajusted
sim3Parallel(nsim,n,c(1,.83,.73),r=2, 2,F)
sim3Parallel(nsim,n,c(1,1,.73),r=1, 2,F)
sim3Parallel(nsim,n,c(1,1,1),r=0, 2,F)
sim3Parallel(nsim,n,c(.96,.87,.73),r=2, 2,F)
sim3Parallel(nsim,n,c(1,.96,.73),r=1, 2,F)
sim3Parallel(nsim,n,c(1,1,.96),r=0, 2,F)
sim3Parallel(nsim,n,c(.92,.87,.73),r=2, 2,F)
sim3Parallel(nsim,n,c(1,.92,.73),r=1, 2,F)
sim3Parallel(nsim,n,c(1,1,.92),r=0, 2,F)