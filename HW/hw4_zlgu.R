setwd("/Users/gzl/Dropbox/STAT547/hw4")

#########  #########  #########    ######### 
######### Q2. Bivariate smoothing  ######### 
#########  #########   #########   ######### 

library(SemiPar) 
data("scallop") 
library(locfit)
dv <- c(1,2,3)
hv <- c(1,2,3)
n <- nrow(scallop) 
m <- 100
png("q2.png",width = 1500, height = 1500, units = "px", pointsize = 30)
par (mfrow = c(3,3))
for (deg in dv) {
  for(h in hv){
    res <- locfit(tot.catch~ lp(latitude, longitude, deg=1, h=h), scallop, ev=lfgrid(mg=m))
    plot(res, type = 'persp', main=paste0('h = ',h, ' degree = ', deg), theta=180, phi=50 ) 
  }
}
dev.off()

######### ######### ######### 
#######  Q3. FPCA   ######### 
######### ######### ######### 

# 3.1 Presmoothing
data <- read.table("yeast.txt", skip  = 4)
X <- as.matrix(data)
T <- as.numeric(substr(names(data), 6, 8))

png("q3a.png",width = 1500, height = 1500, units = "px", pointsize = 30)
par (mfrow = c(3,2))
# original data & smoothing
matplot(T, t(X), type = 'l', main = "Original data")
m <- 30
smooth.X = apply(X, 1,function(y){
  res=locfit(y ~ lp(T, deg = 2, h = 10), ev=lfgrid(mg=m)) 
  predict(res, newdata=T)
})
matplot(T, (smooth.X), type = 'l', main = "Smoothed data, deg = 2, h = 10")


# fpca
fpca <- function(T, X) { 
  mu <- colMeans(X, na.rm = TRUE) # In case missing data points in original datasets
  G <- cov(X, use = 'complete.obs') # In case missing data points in original datasets
  eig <- eigen(G) 
  n <- nrow(X) 
  m <- ncol(X) 
  lam <- eig$values * diff(range(T)) / m 
  phi <- eig$vectors / sqrt(diff(range(T) / m)) 
  Xcenter <- X - matrix(mu, nrow=n, ncol=m, byrow=TRUE) 
  xi <- Xcenter %*% phi * diff(range(T)) / m 
  list(mu=mu, G=G, phi=phi, lam=lam, xi=xi) 
} 

orig.res <- fpca(T, (X))
smooth.res <- fpca(T, t(smooth.X))

orig.mu <- orig.res$mu 
smooth.mu<- smooth.res$mu

plot(T, orig.mu, type = 'l', main="Mean function of original data") 
plot(T, smooth.mu, type = 'l', main="Mean function with presmoothing")
matplot(T, orig.res$phi[, 1:3], type='l', main='The first three eigenfunctions of original data')
matplot(T, smooth.res$phi[, 1:3], type='l', main='The first three eigenfunctions with presmoothing')
dev.off()

# 3.2 first derivative
png("q3b.png",width = 1500, height = 1500, units = "px", pointsize = 30)
par(mfrow = c(3,1))
deriv.X = apply(X, 1, function(y){
  res=locfit(y ~ lp(T, deg = 2, h = 10), deriv=1, ev=lfgrid(mg=m)) 
  predict(res, newdata=T)
})
matplot(T, deriv.X, type = 'l', main= 'Estimated derivative curve, deg = 2, h = 10')
deriv.res= fpca(T, t(deriv.X))
deriv.mu=deriv.res$mu
plot(T, deriv.mu, type = 'l', main='mean function of derivatives')
matplot(T, deriv.res$phi[, 1:3], type='l', main='The first three eigenfunctions of derivative curves')
dev.off()



