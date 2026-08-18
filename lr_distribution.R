library(parallel)
cl = makeCluster(detectCores()-2)
clusterEvalQ(cl, library(MASS))

#generates a value of the asymptotic distribution of the log of LR test

sim1 = function(n, d) {
  if (d==1) {
    a = rnorm(n)
    z = cumsum(a)
    
    m1 = sum(a[-1] * z[-n])      
    m2 = sum(z[-n] ** 2)
    
    dist = m1^2 / m2
    return(dist)
    
  } else{
    mu = rep(0, d)
    Sigma = diag(1, d)
    
    a = t(mvrnorm(n, mu, Sigma))
    z = t(apply(a, 1, cumsum))
    
    m1 = a[, -1] %*% t(z[, -n])
    m2 = z[, -n] %*% t(z[, -n])
    
    dist = m1 %*% solve(m2, t(m1))
    traceDist = sum(diag(dist))
    return(traceDist)
  }
}


#gives the histogram and the quantile of the distribution

distLogLike = function(n,nsim,d) {
  clusterSetRNGStream(cl, iseed = 1)  
  clusterExport(cl,varlist=c("sim1","n","d"))
  tic = Sys.time()
  dist = c(parSapply(cl, 1:nsim, function(i) {sim1(n,d)}))
  tac = Sys.time() - tic
  print(tac)
  hist(dist, seq(min(dist), max(dist)+.2, .2), probability = T,
       main = "", xlab = "", ylab = "")
  print(round(quantile(dist, c(0.005,.01,.025,.05,.1,.5,.9,.95,.975,.99,.995)),3))
}

png("figure10_2.png", width = 1440, height = 480)
par(mfrow = c(1, 3)) 
n=10000; d=1; nsim = 1000000
#pour démarrer parSapply, il faut initialiser n et d sinon ne marche pas
distLogLike(n, 1000000,1)
title(main = "d=1", xlab = "Values", ylab = "Density")
distLogLike(n, 1000000,2)
title(main = "d=2", xlab = "Values", ylab = "Density")
distLogLike(n, 1000000,3)
title(main = "d=3", xlab = "Values", ylab = "Density")
par(mfrow=c(1,1))
dev.off()

#close the connexions of clusters and free memory
stopCluster(cl)

###############################################################################
#different trials with different n for convergence investigation

n=100; d=1 #initialize arguments, if not, don't work
distLogLike(100,1000000,1)
distLogLike(1000,1000000,1)
distLogLike(2000,1000000,1)
distLogLike(5000,1000000,1)
distLogLike(10000,1000000,1)

distLogLike(100,1000000,2)
distLogLike(1000,1000000,2)
distLogLike(2000,1000000,2)
distLogLike(5000,1000000,2)
distLogLike(10000,1000000,2)

distLogLike(100,1000000,3)
distLogLike(1000,1000000,3)
distLogLike(2000,1000000,3)
distLogLike(5000,1000000,3)
distLogLike(10000,1000000,3)

distLogLike(100,1000000,4)
distLogLike(1000,1000000,4)
distLogLike(2000,1000000,4)
distLogLike(5000,1000000,4)
distLogLike(10000,1000000,4)

distLogLike(100,1000000,5)
distLogLike(1000,1000000,5)
distLogLike(2000,1000000,5)
distLogLike(5000,1000000,5)
distLogLike(10000,1000000,5)

######################################################################
#generates a value of the asymptotic distribution of the log of LR test when
#the model includes a constant term

sim2 = function(n, d) {
  if (d==1) {
    a = rnorm(n)
    z = cumsum(a)
    
    m1 = sum(a[-1] * (z[-n]-mean(z[-n]))  )      
    m2 = sum((z[-n]-mean(z[-n])) ** 2)
    
    dist = m1^2 / m2
    return(dist)
    
  } else{
    mu = rep(0, d)
    Sigma = diag(1, d)
    
    a = t(mvrnorm(n, mu, Sigma))
    z = t(apply(a, 1, cumsum))
    
    m1 = a[, -1] %*% t(z[, -n] - rowMeans(z[, -n]) ) 
    m2 = (z[, -n] - rowMeans(z[, -n])) %*% t(z[, -n] - rowMeans(z[, -n])) 
    
    dist = m1 %*% solve(m2, t(m1))
    traceDist = sum(diag(dist))
    return(traceDist)
  }
}

distLogLike2 = function(n,nsim,d) {
  clusterSetRNGStream(cl, iseed = 1)  
  clusterExport(cl,varlist=c("sim2","n","d"))
  tic = Sys.time()
  dist = c(parSapply(cl, 1:nsim, function(i) {sim2(n,d)}))
  tac = Sys.time() - tic
  print(tac)
  hist(dist, seq(min(dist), max(dist)+1, 1), probability = T,
       main = "", xlab = "", ylab = "")
  print(round(quantile(dist, c(0.005,.01,.025,.05,.1,.5,.9,.95,.975,.99,.995)),3))
}

n=10000; d=1 #initialize arguments, if not, don't work
distLogLike2(10000,1000000,1)
distLogLike2(10000,1000000,2)
distLogLike2(10000,1000000,3)
distLogLike2(10000,1000000,4)
distLogLike2(10000,1000000,5)
stopCluster(cl)

png("figure11_1.png", width = 1440, height = 480)
par(mfrow = c(1, 3)) 
n=10000; d=1; nsim = 1000000
#pour démarrer parSapply, il faut initialiser n et d sinon ne marche pas
distLogLike2(n, 1000000,1)
title(main = "d=1", xlab = "Values", ylab = "Density")
distLogLike2(n, 1000000,2)
title(main = "d=2", xlab = "Values", ylab = "Density")
distLogLike2(n, 1000000,3)
title(main = "d=3", xlab = "Values", ylab = "Density")
par(mfrow=c(1,1))
dev.off()

#> distLogLike2(n, 1000000,1)
#Time difference of 4.163034 mins
#0.5%     1%   2.5%     5%    10%    50%    90%    95%  97.5%    99%  99.5% 
#0.001  0.003  0.019  0.073  0.257  2.450  6.590  8.188  9.755 11.788 13.331 
#Il y a eu 28 avis (utilisez warnings() pour les visionner)
#> title(main = "d=1", xlab = "Values", ylab = "Density")
#> distLogLike2(n, 1000000,2)
#Time difference of 14.20776 mins
#0.5%     1%   2.5%     5%    10%    50%    90%    95%  97.5%    99%  99.5% 
#2.474  2.845  3.507  4.174  5.071  9.433 15.871 18.131 20.228 22.852 24.716
#> distLogLike2(n, 1000000,3)
#Time difference of 20.04201 mins
#0.5%     1%   2.5%     5%    10%    50%    90%    95%  97.5%    99%  99.5% 
#9.188  9.959 11.206 12.395 13.902 20.418 29.013 31.878 34.512 37.757 40.137 