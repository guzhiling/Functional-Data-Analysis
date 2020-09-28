library(fds) 
library(fda)
library(fdapace)
setwd("~/Dropbox/STAT547/hw2")

######### ######### ######### 
#########     Q1    ######### 
######### ######### ######### 
# Bivariate PCA
df <- fda ::gait
# Hip angle  39 * 20, rows are individuals
X1 <- t(df[,,1])
#  Knee angle 39 * 20, rows are individuals
X2 <- t(df[,,2])
T <- as.numeric(colnames(X1))

fpca <- function(T, X) { 
  mu <- colMeans(X) 
  G <- cov(X) 
  eig <- eigen(G) 
  n <- nrow(X) 
  m <- ncol(X) 
  lam <- eig$values * diff(range(T)) / m 
  phi <- eig$vectors / sqrt(diff(range(T) / m)) 
  Xcenter <- X - matrix(mu, nrow=n, ncol=m, byrow=TRUE) 
  xi <- Xcenter %*% phi * diff(range(T)) / m 
  list(mu=mu, G=G, phi=phi, lam=lam, xi=xi) 
} 



fpca.bi <- function(T, X1,X2) { 
  mu1 <- colMeans(X1) 
  mu2 <- colMeans(X2)

  n <- nrow(X1) 
  m <- ncol(X1) 
  
  X1center <- X1 - matrix(mu1, nrow=n, ncol=m, byrow=TRUE) 
  X2center <- X2 - matrix(mu2, nrow=n, ncol=m, byrow=TRUE) 
  G12 <- t(X1center) %*% (X2center)
  G21 <- t(X2center) %*% (X1center)
  G11 <- t(X1center) %*% (X1center)
  G22 <- t(X2center) %*% (X2center)
  G <- rbind(cbind(G11,G12),cbind(G21,G22))/(n-1)
  
  eig <- eigen(G) 
  lam <- eig$values * diff(range(T)) / m 
  phi <- eig$vectors / sqrt(diff(range(T) / m)) 
  phi1 <- phi[1:m,]
  phi2 <- phi[(m+1):(2*m),]
  xi1 <- (X1center %*% phi1 )* diff(range(T)) / m 
  xi2 <- X2center %*% phi2 * diff(range(T)) / m 
  list(mu1=mu1, mu2 = mu2, G=G, phi=phi, lam=lam, xi1=xi1,xi2=xi2, phi1 = phi1, phi2 = phi2) 
} 


res <- fpca.bi(T, X1,X2) 
par(mfrow = c(2,2))
plot(T, res$mu1) 
plot(T, res$mu2) 
c=4
png("q1a.png",width = 2000, height = 1200, units = "px", pointsize = 30)
par(mfrow = c(2,2))
matplot(T, res$mu1, type='l', xlab = "Proportion of gait cycle", ylab = "Hip angle", main = "Hip curve for PC1") 
lines(T, res$mu1+c*res$phi1[,1], type='l', col = 2, lty = 2) 
lines(T, res$mu1-c*res$phi1[,1], type='l', col = 4, lty = 2) 
matplot(T, res$mu2, type='l', xlab = "Proportion of gait cycle", ylab = "Knee angle", main = "Knee curve for PC1") 
lines(T, res$mu2+c*res$phi2[,1], type='l', col = 2, lty = 2) 
lines(T, res$mu2-c*res$phi2[,1], type='l', col = 4, lty = 2)
matplot(T, res$mu1, type='l', xlab = "Proportion of gait cycle", ylab = "Hip angle", main = "Hip curve for PC2") 
lines(T, res$mu1+c*res$phi1[,2], type='l', col = 2, lty = 2) 
lines(T, res$mu1-c*res$phi1[,2], type='l', col = 4, lty = 2)
matplot(T, res$mu2, type='l', xlab = "Proportion of gait cycle", ylab = "Knee angle", main = "Knee curve for PC2") 
lines(T, res$mu2+c*res$phi2[,2], type='l', col = 2, lty = 2) 
lines(T, res$mu2-c*res$phi2[,2], type='l', col = 4, lty = 2)
dev.off()

png("q1b.png",width = 2000, height = 1200, units = "px", pointsize = 30)
par(mfrow = c(2,2))
plot(res$mu1,res$mu2, type='l', col = "white", xlab = "hip angle", ylab = "knee angle", main = paste0("Mean Curve"))
text(res$mu1,res$mu2, label=as.character(seq_len(nrow(phi1)))) 

plot(res$mu1,res$mu2, type='l', col = "white", xlab = "hip angle", ylab = "knee angle", main = paste0("PC1 (", round(lam[1]/sum(lam)*100,2), "% of Variation)"))
text(res$mu1,res$mu2, label=as.character(seq_len(nrow(phi1)))) 
arrows(x0=res$mu1,y0=res$mu2, x1 = res$mu1+c*res$phi1[,1], y1 = res$mu2+c*res$phi2[,1], length = 0.1, angle = 30,
       code = 2, col = par("fg"), lty = par("lty"),
       lwd = par("lwd"))
plot(res$mu1,res$mu2, type='l', col = "white", xlab = "hip angle", ylab = "knee angle", main = paste0("PC2 (", round(lam[2]/sum(lam)*100,2), "% of Variation)"))
text(res$mu1,res$mu2, label=as.character(seq_len(nrow(phi1)))) 
arrows(x0=res$mu1,y0=res$mu2, x1 = res$mu1+c*res$phi1[,2], y1 = res$mu2+c*res$phi2[,2], length = 0.1, angle = 30,
       code = 2, col = par("fg"), lty = par("lty"),
       lwd = par("lwd"))

plot(res$mu1,res$mu2, type='l', col = "white", xlab = "hip angle", ylab = "knee angle", main = paste0("PC3 (", round(lam[3]/sum(lam)*100,2), "% of Variation)"))
text(res$mu1,res$mu2, label=as.character(seq_len(nrow(phi1)))) 
arrows(x0=res$mu1,y0=res$mu2, x1 = res$mu1+c*res$phi1[,3], y1 = res$mu2+c*res$phi2[,3], length = 0.1, angle = 30,
       code = 2, col = par("fg"), lty = par("lty"),
       lwd = par("lwd"))
dev.off()



######### ######### ######### 
#########     Q2    ######### 
######### ######### ######### 

## Simulation of functional data 
library(fdapace) 
K <- 20 
n <- 20
m <- 101
T <- seq(0, 1, length.out = m) 
mu <- -(T - 0.5)^2 + 1 
plot(T, mu, ylim=c(0, 1)) 

# ?CreateBasis 
phiSim <- CreateBasis(K=K, T, type='fourier')  # 20 columns of eigen vectors
matplot(T, phiSim[, 1:5], type='l') 
# lamSim <- seq_len(K)^(-2) 
lamSim <- exp(-seq_len(K)) # lam= var explained for each phi (20*1)


# Gaussian Process
# rows are for individuals, col are for different FPC 
png("gauss.png",width = 2000, height = 600, units = "px", pointsize = 30)
par(mfrow = c(1,2))
n=20
xi <- matrix(rnorm(n * K), n, K) %*% diag(sqrt(lamSim))  # variance = lamSim  200*20
X <- matrix(mu, nrow=n, ncol=m, byrow=TRUE) + xi %*% t(phiSim) 
matplot(T, t(X), type='l',ylab =  "Gaussian Processes") 
n=200
xi <- matrix(rnorm(n * K), n, K) %*% diag(sqrt(lamSim))  # variance = lamSim  200*20
X <- matrix(mu, nrow=n, ncol=m, byrow=TRUE) + xi %*% t(phiSim)
X.adj <- apply(X,1,mean)
qqnorm(X.adj);qqline(X.adj, col = 2) # non normal obviously 
legend("topleft", paste("Shapiro P-value =", round(shapiro.test(X.adj)[[2]],3) ), col=2)
dev.off()


# Non-Gaussian Process
png("nongauss.png",width = 2000, height = 600, units = "px", pointsize = 30)
par(mfrow = c(1,2))
n=20
xi <-  matrix(runif(n * K,-1,1), n, K ) %*% diag(sqrt(3*lamSim))
X <- matrix(mu, nrow=n, ncol=m, byrow=TRUE) + xi %*% t(phiSim) 
matplot(T, t(X), type='l',ylab =  "Non-Gaussian Processes") 
n=200
xi <-  matrix(runif(n * K,-1,1), n, K ) %*% diag(sqrt(3*lamSim))
X <- matrix(mu, nrow=n, ncol=m, byrow=TRUE) + xi %*% t(phiSim) 
X.adj1 <- apply(X,1,mean)
qqnorm(X.adj1);qqline(X.adj1, col = 2) 
legend("topleft", paste("Shapiro P-value =", round(shapiro.test(X.adj1)[[2]],3) ), col=2)
dev.off()


######### ######### ######### 
#########     Q3    ######### 
######### ######### ######### 
# The first column contains the name of the genes. 
# The rest of the columns corresponds to 18 equally spaced time points from 0 to 119 minutes
# The first 44 genes are related to G1 phase regulation, and the rest 45 genes are not.
data <- read.table("yeast.txt")
# note there is missing data points
X <- data[, is.na(colSums(data))==0]
T <- seq(0,119, 7)[is.na(colSums(data))==0]


## Non-equally spaced grid 
mu <- colMeans(X) 
G <- cov(X) 
plot(mu) 
plot(T, mu) 
eqGrid <- seq(min(T), max(T), length.out=18) 

# ?ConvertSupport: Perform linear interpolation 
mureg <- fdapace::ConvertSupport(T, eqGrid, mu=mu) 
points(eqGrid, mureg, pch='x') 
Greg <- fdapace::ConvertSupport(T, eqGrid, Cov=G) # Note the isCrossCov argument 
eig <- eigen(Greg) 
lam <- eig$values * diff(range(eqGrid)) / length(eqGrid) 
phi <- eig$vectors / sqrt(diff(range(eqGrid)) / length(eqGrid)) 
matplot(eqGrid, phi[, 1:2], type='l') # first: average height, second: growth spurt 

Xreg <- t(fdapace::ConvertSupport(T, eqGrid, phi=t(X))) 
xi <- (Xreg - matrix(mureg, nrow=nrow(X), ncol=ncol(Xreg), byrow=TRUE)) %*%  
  phi *  diff(range(T)) / length(eqGrid) 

FVE <- cumsum(abs(lam)/ sum(abs(lam)))

png("q3.png",width = 2000, height = 1200, units = "px", pointsize = 30) 
par(mfrow = c(2,2))
plot(xi[, 1], xi[, 2]) 
plot(FVE, xlab = "i-th eigen functions")

plot(eqGrid, mureg, type='l', lwd=2, ylim=c(-1,1)) 
lines(eqGrid, mureg + sqrt(lam[1]) * phi[, 1], type='l', lwd=2, col='red') 
lines(eqGrid, mureg - sqrt(lam[1]) * phi[, 1], type='l', lwd=2, col='blue') 
title('1st') 



plot(eqGrid, mureg, type='l', lwd=2, ylim=c(-1,1)) 
lines(eqGrid, mureg + sqrt(lam[2]) * phi[, 2], type='l', lwd=2, col='red') 
lines(eqGrid, mureg - sqrt(lam[2]) * phi[, 2], type='l', lwd=2, col='blue') 
title('2nd') 
dev.off()
