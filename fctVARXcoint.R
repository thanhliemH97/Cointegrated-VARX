library(MTS)

phi = matrix(c(0.6,.12,1,.7),2,2)
ve = matrix(c(2,4),2,1)
sigma = matrix(c(25,5.4,5.4,9.0),2,2)
phiX = matrix(0.8)
sigmaX = matrix(4)
p=1; s=0; pX=1
phi0=NULL
genVARX = function(n, phi, ve, sigma, p, s, phiX, sigmaX, pX, phi0=NULL){
  
  #Le modèle est dans le sens Y = XB + e plutôt que Y = BX + e
  
  #n : nombre d'observations
  #phi : liste de matrices de coefficients autorégressifs
  #ve : liste de matrices de coefficients régressifs exogènes
  #sigma : matrice de covariance de bruits blancs
  #phiX : array de VAR(pX) pour la variable exogène
  #p : ordre autorégressif maximal
  #s : ordre exogène maximal. s=0 signifie que X_t non décalé est présent
  #pX : ordre autorégressif maximal de la variable exogène
  
  #Si avoir une seule matrice phi, changer sa classe pour une liste.
  if (is.matrix(phi)) phi = list(phi)
  if (is.matrix(ve)) ve = list(ve)
  if (is.numeric(phiX)) phiX = as.matrix(phiX)
  
  k = ncol(phi[[1]]) # dimension de la série chronologique endogène
  l = ncol(phiX) # dimension de la série chronologique exogène
  if (is.null(phi0)) phi0 = rep(0,k)
  
  sigma  = if (missing(sigma))  diag(1, k) else as.matrix(sigma)
  sigmaX = if (missing(sigmaX)) diag(1, l) else as.matrix(sigmaX)
  
  chauffe = n+1 #temps de chauffe de "n+1" 
  nn = n + chauffe 
  m = max(p,s) #max pour éviter des erreurs de récursion dans la simulation
  w = MASS::mvrnorm(nn, mu=base::rep(0,k), Sigma=sigma) #bruit blanc vectoriel
  x = MTS::VARMAsim(nn, arlags = 1:pX, phi = phiX, sigma = sigmaX)$series
  x = as.matrix(x)
  y = matrix(0, ncol=k, nrow=nn) #variable endogène initialisée à 0
  
  PhiMat  = t(do.call(cbind, phi))
  VeMat = if (s >= 0) t(do.call(cbind, ve)) else NULL
  
  #Récursion 
  for (t in (m+1):nn) {
    Ylag = as.vector(t(y[(t-1):(t-p), , drop = FALSE]))
    Xlag = as.vector(t(x[t:(t-s), , drop = FALSE]))
    y[t,] = c(Ylag, Xlag) %*% rbind(PhiMat, VeMat) + t(phi0) + w[t,]
  }
  #enlever les "chauffe" premiers nombres aléatoires
  list(y = y[-(1:chauffe), , drop = F], x = x[-(1:chauffe), , drop = F])
}

ols = function(X,Y){ 
  # X: variables explicatives, Y: variables expliquées

  param = qr.solve(X, Y)
  resi = Y - X %*% param
  cov = (t(resi) %*% resi) / nrow(Y)
  return(list(parameters=param, covariance=cov))
}

VARXorder2 = function(y, x, max_p = 13, max_s = 3, output = T, intercept = T) {
  y = as.matrix(y)
  x = as.matrix(x)
  n = dim(y)[1]
  k = dim(y)[2]
  n_x = dim(x)[1]
  l = dim(x)[2]
  d = max(max_p,max_s)
  
  if (n_x != n) {
    cat("Adjustment made for different nobs:", c(n, n_x), "\n")
    n = min(n, n_x)
  }
  n_prime = n - d
  y_trunc = y[(d+1):n, ]
  
  aic = matrix(0, max_p + 1, max_s + 2)
  rownames(aic) = paste0("p=", 0:max_p)
  colnames(aic) = paste0("s=", -1:max_s)
  bic = aic
  hqc = aic
  
  loglike = function(y,z) {
    y = as.matrix(y)
    n = nrow(y)
    ztz = t(z) %*% z
    zty = t(z) %*% y
    beta = solve(ztz, zty)
    resi = y - z %*% beta
    sigma = t(resi) %*% resi / n
    if (class(sigma)[1] == "matrix") log(det(sigma))
    else log(sigma)
  }
  
  ln_ds = loglike(y_trunc, rep(1, n_prime))
  
  aic[1, 1] = ln_ds
  bic[1, 1] = ln_ds
  hqc[1, 1] = ln_ds
  
  for (p in 1:max_p) {
    if (intercept) {z = rep(1, n_prime)}
    else {z = NULL}
    
    for (i in 1:p) z = cbind(z, y[(d+1 - i):(n - i), ])
    
    ln_ds = loglike(y_trunc, z)
    npar = k * k * p
    
    aic[p + 1, 1] = ln_ds + 2 * npar / n_prime
    bic[p + 1, 1] = ln_ds + log(n_prime) * npar / n_prime
    hqc[p + 1, 1] = ln_ds + 2 * log(log(n_prime)) * npar / n_prime
  }
  
  for (s in 0:max_s) {
    
    if (intercept) {z = rep(1, n_prime)}
    else {z = NULL}
    
    for (j in 0:s) z = cbind(z, x[(d+1 - j):(n - j), ])
    
    ln_ds = loglike(y_trunc, z)
    npar = k * l*(s+1)
    
    aic[1, s + 2] = ln_ds + 2 * npar / n_prime
    bic[1, s + 2] = ln_ds + log(n_prime) * npar / n_prime
    hqc[1, s + 2] = ln_ds + 2 * log(log(n_prime)) * npar / n_prime
    
    for (p in 1:max_p) {
      if (intercept) {z = rep(1, n_prime)}
      else {z = NULL}
      
      for (i in 1:p) z = cbind(z, y[(d+1 - i):(n - i), ])
      
      for (j in 0:s) z = cbind(z, x[(d+1 - j):(n - j), ])
      
      ln_ds = loglike(y_trunc, z)
      npar = k * (p*k + l*(s+1) )
      
      aic[p + 1, s + 2] = ln_ds + 2 * npar / n_prime
      bic[p + 1, s + 2] = ln_ds + log(n_prime) * npar / n_prime
      hqc[p + 1, s + 2] = ln_ds + 2 * log(log(n_prime)) * npar / n_prime
    }
  }
  ind.min = function(A) {
    index = which(A == min(A), arr.ind = TRUE)
    p = index[1,1] - 1
    s = index[1,2] - 2 # because s = -1, 0,1,...
    return(c(p,s))
  }
  
  aic_order = ind.min(aic)
  bic_order = ind.min(bic)
  hqc_order = ind.min(hqc)
  
  if (output) {
    cat("selected order(p,s): aic = ", aic_order, "\n")
    cat("selected order(p,s): hqc = ", hqc_order, "\n")
    cat("selected order(p,s): bic = ", bic_order, "\n")
  }
  VARXorder2 = list(aic = aic, aicor = aic_order,  
                    hqc = hqc, hqcor = hqc_order,
                    bic = bic, bicor = bic_order)
}

VARX.fit = function(y,x,p,s,intercept=T){
  y = as.matrix(y); x = as.matrix(x)
  n = nrow(y); k = ncol(y); l = ncol(x)
  maxlag = max(p,s)
  nobs = n - maxlag  # nombre d'observations retenues
  ylag0 = y[(maxlag+1):n,]
  z = NULL
  namesofz = NULL
  if (intercept) {z=1; namesofz = "intercept"}
  for (i in 1:p) {
    z = cbind(z, y[(maxlag+1 - i):(n - i), ])
    namesofz = c(namesofz, rep(paste0("p=",i),k))
    }
  for (j in 0:s) {
    z = cbind(z, x[(maxlag+1 - j):(n - j), ])
    namesofz = c(namesofz, rep(paste0("s=",j),l))
    }
  ztz = crossprod(z)
  zty = crossprod(z, ylag0)
  para = solve(ztz, zty)
  rownames(para) = namesofz
  resi = ylag0 - z %*% para
  sigma = crossprod(resi) / nobs
  npara = nrow(zty)
  varCoef = sigma %x% solve(ztz) #car on a empilé de sorte que le modèle de
  #régression s'écrive Y = XB+E, donc hat{B} = (X'X)^{-1}X'Y
  
  se = matrix(sqrt(diag(varCoef)), byrow=F, nrow=npara, ncol=k)
  rownames(se) = namesofz; colnames(se) = colnames(y)
  return(list(parameters=para,covariance=sigma,residuals=resi,se=se))
}

ecvarx.fit = function(y,x,p,s){
  y = as.matrix(y)
  x = as.matrix(x)
  n = nrow(y)
  w = diff(y)
  k = ncol(y); l = ncol(x)
  maxlag = max(p,s)
  ylag1 = y[maxlag:(n-1),]
  wlag0 = w[maxlag:(n-1),]

  z = ylag1
  namesofz = rep("C",k)
  if (p>1) {
    for (i in 1:(p-1)) {
      z = cbind(z, w[(maxlag - i):(n - 1 - i), ])
      namesofz = c(namesofz, rep(paste0("p=",i),k))}
  }
  if (s>=0) {
    for (j in 0:s) {
      z = cbind(z, x[(maxlag + 1 - j):(n - j), ])
      namesofz = c(namesofz, rep(paste0("s=",j),l))
      }
  }
  
  para = solve(t(z) %*% z, t(z) %*% wlag0)
  rownames(para) = namesofz
  resi = wlag0 - z %*% para
  sigma = (t(resi) %*% resi) / (n-maxlag)
  varCoef = sigma %x% solve(t(z)%*%z) #car on a empilé de sorte que le modèle de
  #regression s'écrive Y = XB+E, donc hat{B} = (X'X)^{-1}X'Y
  npara = nrow(para)
  se = matrix(sqrt(diag(varCoef)), byrow=F, nrow=npara, ncol=k)
  rownames(se) = namesofz; colnames(se) = colnames(y)
  return(list(parameters=para,covariance=sigma,residuals=resi,se=se))
}

reducedRank = function(C, Sigma, r=1) {
  C = as.matrix(C)
  A = C[,1:r]
  C2 = C[,-(1:r)]
  
  part1 = t(A) %*% qr.solve(Sigma, A)  
  part2 = t(A) %*% qr.solve(Sigma, C2)
  
  B0 = solve(part1, part2)
  
  return(list(A=A, B0=B0))
}

colSD = function(matrix){ #écart-type de chaque colonne d'une matrice
  apply(matrix,2,function(y) sd(y)) #"2" indique que la fonction s'applique aux colonnes
}
colVar = function(matrix){ #variance de chaque colonne d'une matrice
  apply(matrix,2,function(y) var(y)) 
}
colDiff = function(matrix, lag=1, differences=1){
  apply(matrix,2,function(y) diff(y, lag=1, differences=1))
}

VAR.fit = function(y,p,intercept=T){
  y = as.matrix(y)
  n = nrow(y); k = ncol(y)
  ylag0 = y[(p+1):n,]
  z = NULL; namesofz = NULL
  if (intercept) {z=1; namesofz = "intercept"}
  for (i in 1:p) {
    z = cbind(z, y[(p+1 - i):(n - i), ])
    namesofz = c(namesofz, rep(paste0("p=",i),k))
  }
  ztz = t(z) %*% z
  zty = t(z) %*% ylag0
  para = solve(ztz, zty)
  rownames(para) = namesofz
  resi = ylag0 - z %*% para
  sigma = (t(resi) %*% resi) / (n-p)
  k = ncol(y) 
  npara = nrow(zty)
  varCoef = sigma %x% solve(ztz)
  se = matrix(sqrt(diag(varCoef)), byrow=F, nrow=npara, ncol=k)
  rownames(se) = namesofz; colnames(se) = colnames(y)
  return(list(parameters=para,covariance=sigma,residuals=resi,se=se))
}

ecvar.fit = function(y,p,cnst=F){
  y = as.matrix(y)
  n = nrow(y)
  w = diff(y)
  ylag1 = y[p:(n-1),]
  wlag0 = w[p:(n-1),]
  maxlag = p
  if (cnst==T) z = cbind(1,ylag1)
  else z=ylag1
  
  if (p>1) {
    for (i in 1:(p-1)) z = cbind(z, w[(p - i):(n - 1 - i), ])
  }
  para = solve(t(z) %*% z, t(z) %*% wlag0)
  resi = wlag0 - z %*% para
  sigma = (t(resi) %*% resi) / (n-p)
  
  return(list(parameters=para,covariance=sigma,residuals=resi))
}

qs <- function(x) {
  cte <- sqrt(5/3)
  ind0 <- (x == 0)
  ind1 <- (x != 0)
  out <- rep(1, length(x))
  out[ind0] <- 1
  out[ind1] <- 9/(5 * x[ind1]^2 * pi^2) * (sin(cte * pi * x[ind1])/(cte * pi * x[ind1]) - cos(cte * pi * x[ind1]))
  out
}
daniell <- function(x) {
  ind0 <- (x == 0)
  ind1 <- (x != 0)
  out <- x
  out[ind0] <- 1
  out[ind1] <- sin(pi * x[ind1])/(pi * x[ind1])
  out
}

parzen <- function(x) {
  out <- rep(0, length(x))
  A <- (abs(x) > 6/pi)
  B <- (abs(x) > 3/pi) & (abs(x) <= 6/pi)
  C <- (abs(x) <= 3/pi)
  out[A] <- 0
  out[B] <- 2 * (1 - abs((pi * x[B])/6))^3
  out[C] <- 1 - 6 * ((pi * x[C])/6)^2 + 6 * abs((pi * x[C])/6)^3
  out
}


bartlett <- function(x) {
  out <- rep(0, length(x) )
  A <- (abs(x) <= 1)
  B <- (abs(x) > 1)
  out[A] <- 1 - abs(x[A])
  out[B] <- 0
  out
}


tronque <- function(x) {
  out <- rep(0, length(x))
  A <- (abs(x) <= 1)
  B <- (abs(x) > 1)
  out[A] <- 1
  out[B] <- 0
  out
}
##fonction de calcul de l'autocovariance d'ordre k###
acf.lagk.mult <- function(k,X){
  # la matrice des donnees X est de dimension n x d
  n <- nrow(X)
  d <- ncol(X)
  meanX <- apply(X,2,mean)
  output <- matrix(0, nrow=d, ncol=d)
  for (t in (k+1): n){
    mat <- (X[t,] - meanX) %*% t(X[t-k, ] - meanX)
    output <- output+mat/n      
  }
  output
}

##########la fonction de calcul de la statistique#########
trmat <- function(A) {
  sum(diag(A))
}
statTn <- function(res, noyau, Pn){
  n <- nrow(res)
  nm1 <- n-1
  d <- ncol(res)
  Mn <- 0
  Vn <- 0
  Qn <- 0
  sigmam1 <- solve(acf.lagk.mult(0,res))
  for (j in 1:nm1){
    acf.lagj <- acf.lagk.mult(j, res)
    mat <- t( acf.lagj ) %*%  sigmam1  %*%  acf.lagj  %*%  sigmam1
    Mn <- Mn + (1-j/n)*( noyau(j/Pn) )^2
    Vn <- Vn + (1- j/n)*(1-(j+1)/n)*( noyau(j/Pn) )^4
    Qn <- Qn+trmat(mat)*( noyau(j/Pn) )^2
  }
  Tn <- ( n*Qn-(d^2)*Mn )/sqrt(2*d*d*Vn)
  return(Tn)
}

# newtonRaphson1 = function(par){
# 
#   condition = Inf; i = 1
# 
#   while(condition>1e-10){
# 
#     A = matrix(par[1:2],2,1)
#     B = matrix(c(1,par[3]),1,2)
#     V0 = matrix(par[4:5],2,1)
# 
#     Epsilon =  W - Ylag1 %*% t(A%*%B) - Xtrunc %*% t(V0)
#     cov = (t(Epsilon) %*% Epsilon) / nrow(W)
# 
#     sum1  = matrix(0, 5, 5); sum2 = matrix(0, 5, 1)
# 
#     for (t in 2:n) {
# 
#       bloc1 = Y[t-1,2]  %x% A
#       Ut = cbind(t(B%*%Y[t-1,]), t(X[t,]))
#       bloc2 = Ut %x% diag(1,2)
#       Zstar_t = cbind(bloc1,bloc2)
#       sum1 = sum1 + t(Zstar_t) %*% solve(cov, Zstar_t)
#       sum2 = sum2 + t(Zstar_t) %*% solve(cov, Epsilon[t-1,])
#     }
#     remain = solve(sum1, sum2)
#     newpar = par + remain[c(2,3,1,4,5),]
#     par = newpar
#     condition = sqrt(sum(remain^2))
#     i = i+1
#   }
#   print(i)
#   return(par)
# }