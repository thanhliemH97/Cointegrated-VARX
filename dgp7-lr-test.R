source("~/these_doctorat/Simulations/fctVARXcoint.R")

library(parallel)
cl = makeCluster(detectCores()-1)
clusterEvalQ(cl, {library(MASS) 
  library(MTS)})

#Le processus générateur de données 4 est un VARX(2,1) avec
#3 variables endogènes et 3 variables exogènes.

sim4.1 = function(n, eigenvalues, r) {
  
  ## --- Constantes fixes ---
  Phi1star = matrix(c(.15,.19,.11,
                      -.2,-.08,-.05,
                      .45,.30,.32), 3,3)
  SigmaEps = matrix(c(.47,.2,.18,
                      .2,.32,.27,
                      .18,.27,.30), 3,3)
  phix.1 = matrix(c(.2,-.6,.1,
                    .4,.3,.7,
                    -.4,.3,.2),3,3)
  sigma.b = matrix(c(4,.8,.6,
                     .8,2,.1,
                     .6,.1,1),3,3)
  V0 = matrix(c(.08,.13,-.05,
                .19,-.26,.11,
                -.4,.32,.22),3,3)
  V1 = diag(.1,3)
  Q = matrix(c(-.29,-.01,-.75,
               -.47,-.85,1.39,
               -.57,1,-.55), 3,3)
  P = solve(Q)
  
  ## --- Matrices autorégressives à varier ---
  
  D = diag(eigenvalues,3)
  C = P %*% D %*% Q - diag(1,3)
  Phi2 = -Phi1star
  Phi1 = C + diag(1,3) + Phi1star
  
  ## --- Simulation VARX(2,1) ---
  out = genVARX(n, list(Phi1,Phi2), list(V0,V1), SigmaEps, p=2, s=1, 
                phix.1, sigma.b, pX=1)  
  Y = out$y
  X = out$x
  
  ## --- Construction des variables retardées ---
  W = diff(Y)
  Ylag1 = Y[2:(n-1),]
  Wlag0 = W[2:(n-1),]
  Xlag0 = X[3:n,]
  Xlag1 = X[2:(n-1),]
  nobs = nrow(Wlag0)
  Z = cbind(Ylag1, Xlag0, Xlag1)
  
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
      V0 = matrix(para[9:17],3,3)
      V1 = matrix(para[18:26],3,3)
      terme.AB = Ylag1 %*% t(A %*% B)
    } else if (r == 1) {
      A = matrix(para[1:3], 3, 1)
      B = matrix(c(1, para[4:5]), 1, 3)
      V0 = matrix(para[6:14], 3, 3)
      V1 = matrix(para[15:23], 3, 3)
      terme.AB = Ylag1 %*% t(A %*% B)
    } else {
      V0 = matrix(para[1:9], 3, 3)
      V1 = matrix(para[10:18], 3, 3)
      terme.AB = 0
    }
    
    Eps = Wlag0 - terme.AB - Xlag0 %*% t(V0) - Xlag1 %*% t(V1)
    cov = crossprod(Eps) / nobs
    log(det(cov))
  }
  
  # --- Minimisation de la log-vraisemblance modifiée ---
  opt = nlminb(start = start.ls, objective = loglikelihood)
  para.mle = opt$par
  
  # --- Retour : statistique du test ---
  nobs * (loglikelihood(para.mle) - log(det(cov.ls)))
}
sim4.2 = function(n, eigenvalues, r) {
  
  ## --- Constantes fixes ---
  Phi1star = matrix(c(.15,.19,.11,
                      -.2,-.08,-.05,
                      .45,.30,.32), 3,3)
  SigmaEps = matrix(c(.47,.2,.18,
                      .2,.32,.27,
                      .18,.27,.30), 3,3)
  phix.1 = matrix(c(.2,-.6,.1,
                    .4,.3,.7,
                    -.4,.3,.2),3,3)
  sigma.b = matrix(c(4,.8,.6,
                     .8,2,.1,
                     .6,.1,1),3,3)
  V0 = matrix(c(.08,.13,-.05,
                .19,-.26,.11,
                -.4,.32,.22),3,3)
  V1 = diag(.1,3)
  Q = matrix(c(-.29,-.01,-.75,
               -.47,-.85,1.39,
               -.57,1,-.55), 3,3)
  P = solve(Q)
  
  ## --- Matrices autorégressives à varier ---
  
  D = diag(eigenvalues,3)
  C = P %*% D %*% Q - diag(1,3)
  Phi2 = -Phi1star
  Phi1 = C + diag(1,3) + Phi1star
  
  ## --- Simulation VARX(2,1) ---
  out = genVARX(n, list(Phi1,Phi2), list(V0,V1), SigmaEps, p=2, s=1, 
                phix.1, sigma.b, pX=1)  
  Y = out$y
  X = out$x
  
  ## --- Construction des variables retardées ---
  W = diff(Y)
  Ylag1 = Y[2:(n-1),]
  Wlag0 = W[2:(n-1),]
  Wlag1 = W[1:(n-2),]
  Xlag0 = X[3:n,]
  Xlag1 = X[2:(n-1),]
  nobs = nrow(Wlag0)
  
  Z = cbind(Ylag1, Wlag1, Xlag0, Xlag1)
  
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
      V0 = matrix(para[18:26], 3, 3)
      V1 = matrix(para[27:35], 3, 3)
      terme.AB = Ylag1 %*% t(A %*% B)
    } else if (r == 1) {
      A = matrix(para[1:3], 3, 1)
      B = matrix(c(1, para[4:5]), 1, 3)
      Phi1s = matrix(para[6:14], 3, 3)
      V0 = matrix(para[15:23], 3, 3)
      V1 = matrix(para[24:32], 3, 3)
      terme.AB = Ylag1 %*% t(A %*% B)
    } else {
      Phi1s = matrix(para[1:9], 3, 3)
      V0 = matrix(para[10:18], 3, 3)
      V1 = matrix(para[19:27], 3, 3)
      terme.AB = 0
    }
    Eps = Wlag0 - terme.AB - Wlag1 %*% t(Phi1s) - Xlag0 %*% t(V0) - Xlag1 %*% t(V1)
    cov = crossprod(Eps) / nobs
    log(det(cov))
  }
  
  ## --- Maximisation ---
  opt = nlminb(start = start.ls, objective = loglikelihood)
  para.mle = opt$par
  
  ## --- Retour : statistique du test ---
  nobs * (loglikelihood(para.mle) - log(det(cov.ls)))
}
sim4.3 = function(n, eigenvalues, r) {
  
  ## --- Constantes fixes ---
  Phi1star = matrix(c(.15,.19,.11,
                      -.2,-.08,-.05,
                      .45,.30,.32), 3,3)
  SigmaEps = matrix(c(.47,.2,.18,
                      .2,.32,.27,
                      .18,.27,.30), 3,3)
  phix.1 = matrix(c(.2,-.6,.1,
                    .4,.3,.7,
                    -.4,.3,.2),3,3)
  sigma.b = matrix(c(4,.8,.6,
                     .8,2,.1,
                     .6,.1,1),3,3)
  V0 = matrix(c(.08,.13,-.05,
                .19,-.26,.11,
                -.4,.32,.22),3,3)
  V1 = diag(.1,3)
  Q = matrix(c(-.29,-.01,-.75,
               -.47,-.85,1.39,
               -.57,1,-.55), 3,3)
  P = solve(Q)
  
  ## --- Matrices autorégressives à varier ---
  
  D = diag(eigenvalues,3)
  C = P %*% D %*% Q - diag(1,3)
  Phi2 = -Phi1star
  Phi1 = C + diag(1,3) + Phi1star
  
  ## --- Simulation VARX(2,1) ---
  out = genVARX(n, list(Phi1,Phi2), list(V0,V1), SigmaEps, p=2, s=1, 
                phix.1, sigma.b, pX=1)  
  Y = out$y
  X = out$x
  
  ## --- Construction des variables retardées ---
  W = diff(Y)
  Ylag1 = Y[3:(n-1),]
  Wlag0 = W[3:(n-1),] 
  Wlag1 = W[2:(n-2),]
  Wlag2 = W[1:(n-3),]
  Xlag0 = X[4:n,]
  Xlag1 = X[3:(n-1),]
  
  nobs = nrow(Wlag0)
  Z = cbind(Ylag1, Wlag1, Wlag2, Xlag0, Xlag1)
  
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
      V0 = matrix(para[27:35], 3, 3)
      V1 = matrix(para[36:44],3,3)
      terme.AB = Ylag1 %*% t(A %*% B)
    } else if (r == 1) {
      A = matrix(para[1:3], 3, 1)
      B = matrix(c(1,para[4:5]), 1, 3)
      Phi1star = matrix(para[6:14], 3, 3)
      Phi2star = matrix(para[15:23], 3, 3)
      V0 = matrix(para[24:32], 3, 3)
      V1 = matrix(para[33:41], 3,3)
      terme.AB = Ylag1 %*% t(A %*% B)
    } else {
      Phi1star = matrix(para[1:9], 3, 3)
      Phi2star = matrix(para[10:18], 3, 3)
      V0 = matrix(para[19:27], 3, 3)
      V1 = matrix(para[28:36], 3, 3)
      terme.AB = 0
    }
    
    Eps = Wlag0 - terme.AB - Wlag1 %*% t(Phi1star) - Wlag2 %*% t(Phi2star) -
      Xlag0 %*% t(V0) - Xlag1 %*% t(V1)
    cov = crossprod(Eps) / nobs
    log(det(cov))
  }
  
  ## --- Maximisation ---
  opt = nlminb(start = start.ls, objective = loglikelihood)
  para.mle = opt$par
  
  ## --- Retour : statistique du test ---
  nobs * (loglikelihood(para.mle) - log(det(cov.ls)))
}
sim4.4 = function(n, eigenvalues, r) {
  
  ## --- Constantes fixes ---
  Phi1star = matrix(c(.15,.19,.11,
                      -.2,-.08,-.05,
                      .45,.30,.32), 3,3)
  SigmaEps = matrix(c(.47,.2,.18,
                      .2,.32,.27,
                      .18,.27,.30), 3,3)
  phix.1 = matrix(c(.2,-.6,.1,
                    .4,.3,.7,
                    -.4,.3,.2),3,3)
  sigma.b = matrix(c(4,.8,.6,
                     .8,2,.1,
                     .6,.1,1),3,3)
  V0 = matrix(c(.08,.13,-.05,
                .19,-.26,.11,
                -.4,.32,.22),3,3)
  V1 = diag(.1,3)
  Q = matrix(c(-.29,-.01,-.75,
               -.47,-.85,1.39,
               -.57,1,-.55), 3,3)
  P = solve(Q)
  
  ## --- Matrices autorégressives à varier ---
  
  D = diag(eigenvalues,3)
  C = P %*% D %*% Q - diag(1,3)
  Phi2 = -Phi1star
  Phi1 = C + diag(1,3) + Phi1star
  
  ## --- Simulation VARX(2,1) ---
  out = genVARX(n, list(Phi1,Phi2), list(V0,V1), SigmaEps, p=2, s=1, 
                phix.1, sigma.b, pX=1)  
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
  Xlag1 = X[4:(n-1),]
  
  nobs = nrow(Wlag0)
  Z = cbind(Ylag1, Wlag1, Wlag2, Wlag3, Xlag0, Xlag1)
  
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
      V0 = matrix(para[36:44], 3, 3)
      V1 = matrix(para[45:53], 3, 3)
      terme.AB = Ylag1 %*% t(A %*% B)
    } else if (r == 1) {
      A = matrix(para[1:3], 3, 1)
      B = matrix(c(1,para[4:5]), 1, 3)
      Phi1star = matrix(para[6:14], 3, 3)
      Phi2star = matrix(para[15:23], 3, 3)
      Phi3star = matrix(para[24:32], 3, 3)
      V0 = matrix(para[33:41], 3, 3)
      V1 = matrix(para[42:50], 3, 3)
      terme.AB = Ylag1 %*% t(A %*% B)
    } else {
      Phi1star = matrix(para[1:9], 3, 3)
      Phi2star = matrix(para[10:18], 3, 3)
      Phi3star = matrix(para[19:27], 3, 3)
      V0 = matrix(para[28:36], 3, 3)
      V1 = matrix(para[37:45], 3, 3)
      terme.AB = 0
    }
    
    Eps = Wlag0 - terme.AB - Wlag1 %*% t(Phi1star) - Wlag2 %*% t(Phi2star) -
      Wlag3 %*% t(Phi3star) - Xlag0 %*% t(V0) - Xlag1 %*% t(V1)
    cov = crossprod(Eps) / nobs
    log(det(cov))
  }
  
  ## --- Maximisation ---
  opt = nlminb(start = start.ls, objective = loglikelihood)
  para.mle = opt$par
  
  ## --- Retour : statistique du test ---
  nobs * (loglikelihood(para.mle) - log(det(cov.ls)))
}
sim4.5 = function(n, eigenvalues, r) {
  
  ## --- Constantes fixes ---
  Phi1star = matrix(c(.15,.19,.11,
                      -.2,-.08,-.05,
                      .45,.30,.32), 3,3)
  SigmaEps = matrix(c(.47,.2,.18,
                      .2,.32,.27,
                      .18,.27,.30), 3,3)
  phix.1 = matrix(c(.2,-.6,.1,
                    .4,.3,.7,
                    -.4,.3,.2),3,3)
  sigma.b = matrix(c(4,.8,.6,
                     .8,2,.1,
                     .6,.1,1),3,3)
  V0 = matrix(c(.08,.13,-.05,
                .19,-.26,.11,
                -.4,.32,.22),3,3)
  V1 = diag(.1,3)
  Q = matrix(c(-.29,-.01,-.75,
               -.47,-.85,1.39,
               -.57,1,-.55), 3,3)
  P = solve(Q)
  
  ## --- Matrices autorégressives à varier ---
  
  D = diag(eigenvalues,3)
  C = P %*% D %*% Q - diag(1,3)
  Phi2 = -Phi1star
  Phi1 = C + diag(1,3) + Phi1star
  
  ## --- Simulation VARX(2,1) ---
  out = genVARX(n, list(Phi1,Phi2), list(V0,V1), SigmaEps, p=2, s=1, 
                phix.1, sigma.b, pX=1)  
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
sim4.6 = function(n, eigenvalues, r) {
  
  ## --- Constantes fixes ---
  Phi1star = matrix(c(.15,.19,.11,
                      -.2,-.08,-.05,
                      .45,.30,.32), 3,3)
  SigmaEps = matrix(c(.47,.2,.18,
                      .2,.32,.27,
                      .18,.27,.30), 3,3)
  phix.1 = matrix(c(.2,-.6,.1,
                    .4,.3,.7,
                    -.4,.3,.2),3,3)
  sigma.b = matrix(c(4,.8,.6,
                     .8,2,.1,
                     .6,.1,1),3,3)
  V0 = matrix(c(.08,.13,-.05,
                .19,-.26,.11,
                -.4,.32,.22),3,3)
  V1 = diag(.1,3)
  Q = matrix(c(-.29,-.01,-.75,
               -.47,-.85,1.39,
               -.57,1,-.55), 3,3)
  P = solve(Q)
  
  ## --- Matrices autorégressives à varier ---
  
  D = diag(eigenvalues,3)
  C = P %*% D %*% Q - diag(1,3)
  Phi2 = -Phi1star
  Phi1 = C + diag(1,3) + Phi1star
  
  ## --- Simulation VARX(2,1) ---
  out = genVARX(n, list(Phi1,Phi2), list(V0,V1), SigmaEps, p=2, s=1, 
                phix.1, sigma.b, pX=1)  
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
      V0 = matrix(para[18:26], 3, 3)
      terme.AB = Ylag1 %*% t(A %*% B)
    } else if (r == 1) {
      A = matrix(para[1:3], 3, 1)
      B = matrix(c(1, para[4:5]), 1, 3)
      Phi1s = matrix(para[6:14], 3, 3)
      V0 = matrix(para[15:23], 3, 3)
      terme.AB = Ylag1 %*% t(A %*% B)
    } else {
      Phi1s = matrix(para[1:9], 3, 3)
      V0 = matrix(para[10:18], 3, 3)
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

sim4Parallel = function(nsim, n, eigenvalues, r, p, s=1) {
  clusterSetRNGStream(cl, iseed = 1)
  clusterExport(cl,varlist=c("n","genVARX","ols","reducedRank",
                             "sim4.1", "sim4.2", "sim4.3", "sim4.4", 
                             "sim4.5", "sim4.6"))
  tic = Sys.time()
  if (s==1) {
    if (p==1) {
      testStat = parSapply(cl,1:nsim, function(i) sim4.1(n,eigenvalues,r))
    } else if (p==2) {
      testStat = parSapply(cl,1:nsim, function(i) sim4.2(n,eigenvalues,r))
    } else if (p==3) {
      testStat = parSapply(cl,1:nsim, function(i) sim4.3(n,eigenvalues,r))
    } else if (p==4) {
      testStat = parSapply(cl,1:nsim, function(i) sim4.4(n,eigenvalues,r))
    }
  } else {
    if (s==-1) {
      testStat = parSapply(cl,1:nsim, function(i) sim4.5(n,eigenvalues,r))
    } else if (s==0){
      testStat = parSapply(cl,1:nsim, function(i) sim4.6(n,eigenvalues,r))
    }
  }
  
  tac = Sys.time() - tic
  print(tac)
  
  if (r==2) {
    crit = c(4.122, 2.974) #Case 1 (a)
  } else if (r==1) {
    crit = c(12.336, 10.481) #Case 2(a)
  } else if (r==0) {
    crit = c(24.29, 21.787) #Case 3(a)
  }
  #hist(testStat, seq(min(testStat), max(testStat) + 0.2, 0.2), probability = T,
  #     main = "", xlab = "", ylab = "")
  #print(round(quantile(testStat, c(.005,.01,.025,.05,.1,.5,.9,.95,.975,.99,.995)),3))
  return(c(mean(testStat>crit[1]),mean(testStat>crit[2])))
}

#Ajuster VARX(1,1), VARX(2,1), VARX(3,1) et VARX(4,1) aux données VARX(2,1)
nsim = 10000; n = 400; p=4 #varier p
sim4Parallel(nsim, n, c(1,.87,.73),   r=2, p)
sim4Parallel(nsim, n, c(1,1,.73),     r=1, p)
sim4Parallel(nsim, n, c(1,1,1),       r=0, p)
sim4Parallel(nsim, n, c(.96,.87,.73), r=2, p)
sim4Parallel(nsim, n, c(1,.96,.73),   r=1, p)
sim4Parallel(nsim, n, c(1,1,.96),     r=0, p)
sim4Parallel(nsim, n, c(.92,.87,.73), r=2, p)
sim4Parallel(nsim, n, c(1,.92,.73),   r=1, p)
sim4Parallel(nsim, n, c(1,1,.92),     r=0, p)

s=1 #varier s
sim4Parallel(nsim, n, c(1,.87,.73),   r=2, 2, s)
sim4Parallel(nsim, n, c(1,1,.73),     r=1, 2, s)
sim4Parallel(nsim, n, c(1,1,1),       r=0, 2, s)
sim4Parallel(nsim, n, c(.96,.87,.73), r=2, 2, s)
sim4Parallel(nsim, n, c(1,.96,.73),   r=1, 2, s)
sim4Parallel(nsim, n, c(1,1,.96),     r=0, 2, s)
sim4Parallel(nsim, n, c(.92,.87,.73), r=2, 2, s)
sim4Parallel(nsim, n, c(1,.92,.73),   r=1, 2, s)
sim4Parallel(nsim, n, c(1,1,.92),     r=0, 2, s)
