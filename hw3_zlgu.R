library(fds) 
library(fda)
library(fdapace)
library(SemiPar) 
library (gss)
library(locfit)
setwd("~/Dropbox/STAT547/hw3")

######### ######### ######### 
#########     Q2    ######### 
######### ######### ######### 


n.time <- 1000 # numer of points on each trajectory
n.real <- 50 # number of trajectories
n.base <- 10 # number of basis
T <- seq(0,1, length.out = n.time)
Z <- matrix(rnorm(n.time * n.real), ncol = n.real, nrow = n.base) 

base.X<- function(t){
    1/ (((1:n.base)-.5)*pi) * sqrt(2) * sin(pi * t * ((1:n.base)-.5))
}

base.B<- function(t){
  1/ ((1:n.base)*pi) * sqrt(2) * sin(pi * t * (1:n.base))
}



# ---------------------------
# Brownian motion simulation
# ---------------------------
##  every column is a scaled base on [0,1]
png("q2a.png",width = 2000, height = 1000, units = "px", pointsize = 30)
par(mfrow = c(1,2))
X.base <- t(sapply(T, base.X))
matplot(T, X.base, type='l', xlab = "time", ylab ="Scaled Base of Brownian Motion")
X <-  X.base %*% Z
matplot(T, X, type='l', xlab = "time", ylab ="50 Realizations Brownian Process")
dev.off()

# ---------------------------
# Brownian bridge simulation
# ---------------------------
##  every column is a scaled base on [0,1]
png("q2b.png",width = 2000, height = 1000, units = "px", pointsize = 30)
par(mfrow = c(1,2))
B.base <- t(sapply(T, base.B))
matplot(T, B.base, type='l', xlab = "time", ylab ="Scaled Base of Brownian Bridge")
B <-  B.base %*% Z
matplot(T, B, type='l', xlab = "time", ylab ="50 Realizations Brownian Bridge")
dev.off()


######### ######### ######### 
#########     Q4    ######### 
######### ######### ######### 
# install.packages("gss")
data ( LakeAcidity )
str(LakeAcidity)


######### ######### ######### 
#########     Q5   ######### 
######### ######### ######### 

# scale the original data for better approximation
cal <- LakeAcidity$cal
ph <- LakeAcidity$ph  
plot(ph,cal , xlab = "Surface PH level", ylab = "Calcium concentration level", main = "Plot of Calcium against PH level")
png("q5a.png",width = 2000, height = 1000, units = "px", pointsize = 30)
par(mfrow = c(1,2))
# 1-D design plot
stripchart(cal, pch='|',xlab = "Surface PH level",  main = "Dot plot of PH levels")
# Kernel density estimate
plot(density(cal), main = "KDE of PH levels")
dev.off()


# ------------------------------------------------------
#     Comparison between different smoothing methods
# ------------------------------------------------------

#------------------#
# Local polynomial 
#------------------#
# Loader, C. (1999) Local Regression and Likelihood

df <- data.frame(X=ph,Y = cal)
n <- length(cal)
png("q5b1.png",width = 2000, height = 800, units = "px", pointsize = 30)
par(mfrow = c(1,3))
deg <- 1; h = .05
res <- locfit(Y ~ lp(X, deg= deg, h = h), df)
plot(df,  xlab = "Surface PH level", ylab = "Calcium concentration level", main =paste0( "p:",deg, " h: ", h))
lines(res)
deg <- 3
res <- locfit(Y ~ lp(X, deg= deg), df)
plot(df,  xlab = "Surface PH level", ylab = "Calcium concentration level", main =paste0( "p:",deg, " h: ", h))
lines(res)
deg <- 5
res <- locfit(Y ~ lp(X, deg= deg), df)
plot(df,  xlab = "Surface PH level", ylab = "Calcium concentration level", main =paste0( "p:",deg, " h: ", h))
lines(res)
dev.off()


png("q5b2.png",width = 2000, height = 800, units = "px", pointsize = 30)
par(mfrow = c(1,3))
deg <- 2; h = 1
res <- locfit(Y ~ lp(X, deg= deg, h = h), df)
plot(df,  xlab = "Surface PH level", ylab = "Calcium concentration level", main =paste0( "p:",deg, " h: ", h))
lines(res)
h = 2
res <- locfit(Y ~ lp(X, deg= deg), df)
plot(df,  xlab = "Surface PH level", ylab = "Calcium concentration level", main =paste0( "p:",deg, " h: ", h))
lines(res)
h = 3
res <- locfit(Y ~ lp(X, deg= deg), df)
plot(df,  xlab = "Surface PH level", ylab = "Calcium concentration level", main =paste0( "p:",deg, " h: ", h))
lines(res)
dev.off()


#----------------------------#
# Smoothing splines
#----------------------------#
png("q5b3.png",width = 1000, height = 800, units = "px", pointsize = 30)
plot(df,  xlab = "Surface PH level", ylab = "Calcium concentration level", main =paste0( "Smoothing spline with different lambdas"), pch =16)
lines(smooth.spline( LakeAcidity$ph  , LakeAcidity$cal, lambda=1e-1), lty=2, col=1,lwd=2) 
lines(smooth.spline( LakeAcidity$ph  , LakeAcidity$cal, lambda=1e-3), lty=3, col=2,lwd=2) 
lines(smooth.spline( LakeAcidity$ph  , LakeAcidity$cal, lambda=1e-4), lty=4, col=3,lwd=2)
lines(smooth.spline( LakeAcidity$ph  , LakeAcidity$cal, lambda=1e-5), lty=5, col=4,lwd=2)
legend("topleft",legend = c(1e-1,1e-3,1e-4,1e-5), col = 1:4, lty = 2:5, lwd =2)
dev.off()

#----------------------------#
# Smoothing splines
#----------------------------#

# B-splines basis 
library(splines) 
res <- bs(x, intercept=TRUE, 
          knots=interior,
          degree = d - 1,
          Boundary.knots = boundary) 
newGrid <- seq(min(ph), max(ph), length=100) 
phNew <- (ph - min(ph)) / diff(range(ph)) 


png("q5b4.png",width = 1000, height = 800, units = "px", pointsize = 30)
plot(df,  xlab = "Surface PH level", ylab = "Calcium concentration level", main =paste0( "B-spline with different number of interior knots"), pch =16)
d <- 4 

for(nKnots in c(1,2,3,4,5)){
  # nKnots <- 5 
  knots01 <- seq(0, 1, length.out = nKnots + 2) 
  knots01 <- knots01[-c(1, length(knots01))] 
  knots <- knots01 * diff(range(ph)) + min(ph) 
  X <- bs(ph, degree=d-1, knots=knots, intercept = TRUE) 
  dat2 <- data.frame(Y = cal, X = X) 
  m2 <- lm(Y ~ . - 1, dat2) 
  xNew2 <- bs(newGrid, degree=d - 1, knots=knots, intercept=TRUE) 
  pred2 <- predict(m2, newdata=data.frame(X=xNew2)) 
  lines(newGrid, pred2, lty=2, lwd = 2,col = nKnots) 
}

legend("topleft",legend = 1:5, col = 1:5, lty = 2, lwd =2)
dev.off()



#----------------------------#
# Penalized splines 
#----------------------------#

png("q5b5.png",width = 1000, height = 800, units = "px", pointsize = 30)
plot(df,  xlab = "Surface PH level", ylab = "Calcium concentration level", main =paste0( "Penalized spline with different  J, lambda = 1e-4"), pch =16)
for(nknots in c(5,10,20,40,50)){
  lines(smooth.spline(LakeAcidity$ph  , LakeAcidity$cal, lambda=1e-4, nknots = nknots),  lty=2, lwd=2, col=nknots) 
}
legend("topleft",legend = c(5,10,20,40,50), col =c(5,10,20,40,50), lty = 2, lwd =2)
dev.off()


png("q5b6.png",width = 1000, height = 800, units = "px", pointsize = 30)
plot(df,  xlab = "Surface PH level", ylab = "Calcium concentration level", main =paste0( "Penalized spline with different lambdas, J =20"), pch =16)
for(l in c(1e-1,1e-2,1e-3,1e-4,1e-5)){
  lines(smooth.spline(LakeAcidity$ph  , LakeAcidity$cal, lambda=l, nknots = 20),  lty=2, lwd=2, col=l*1e5) 
}
legend("topleft",legend = c(1e-1,1e-2,1e-3,1e-4,1e-5), col =c(1e-1,1e-2,1e-3,1e-4,1e-5)*1e5, lty = 2, lwd =2)
dev.off()





