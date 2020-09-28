setwd("/Users/gzl/Dropbox/STAT547/hw5/latex")
library(fdapace)

###  Q1   ###
n <- 50  # number of samples
M <- 100 # number of obs per sample
c0 <- 1  # coeff for eigen values
K <- 3  # number of eigen components
##   Simulate dense functional data
df.d <- MakeGPFunctionalData(n, M =M, mu = rep(0, M), K =K,
                     lambda = sapply(rep(1, K), function(x) {c0 * x^(-2)}), sigma = 0, basisType = "cos")

##   Simulate sparse functional data
df.s <- MakeSparseGP(n, rdist = runif, sparsity = 2:9, muFun = function(x)
                    rep(0, length(x)), K = K, lambda = rep(1, K), sigma = 0,
                    basisType = "cos", CovFun = NULL)

## (a) Compare component number between different selection methods
#### dense case
df.d.Ly <- split(df.d$Y, row(df.d$Y)) # change the original matrix form into list
df.d.Lt <- rep(list(df.d$pts),n)
res.d1 <- FPCA(df.d.Ly,df.d.Lt, optns = list(methodSelectK = 'FVE'))
res.d2 <- FPCA(df.d.Ly,df.d.Lt, optns = list(methodSelectK = 'AIC'))
res.d3 <- FPCA(df.d.Ly,df.d.Lt, optns = list(methodSelectK = 'BIC'))
res.compare.d <- data.frame("TRUE" = K,"FVE" = res.d1$selectK,"AIC" = res.d2$selectK,"BIC" = res.d3$selectK )
rownames(res.compare.d) <- "Number of Component: dense"
res.compare.d

#### sparse case
res.s1 <- FPCA(df.s$Ly,df.s$Lt, optns = list(methodSelectK = 'FVE'))
res.s2 <- FPCA(df.s$Ly,df.s$Lt, optns = list(methodSelectK = 'AIC'))
res.s3 <- FPCA(df.s$Ly,df.s$Lt, optns = list(methodSelectK = 'BIC'))
res.compare.s <- data.frame("TRUE" = K,"FVE" = res.s1$selectK,"AIC" = res.s2$selectK,"BIC" = res.s3$selectK )
rownames(res.compare.s) <- "Number of Component: sparse"
res.compare.s


## (b) & (c) Compare estimated FPC scores between integration method and BLUP
n <- 50  # number of samples
M <- 20 # number of obs per sample
c0 <- 1  # coeff for eigen values
K <- 3  # number of eigen components
##   Simulate dense functional data
df.d <- MakeGPFunctionalData(n, M =M, mu = rep(0, M), K =K,
                             lambda = sapply(rep(1, K), function(x) {c0 * x^(-2)}), sigma = 0, basisType = "cos")

##   Simulate sparse functional data
df.s <- MakeSparseGP(n, rdist = runif, sparsity = 2:3, muFun = function(x)
  rep(0, length(x)), K = K, lambda = rep(1, K), sigma = 0,
  basisType = "cos", CovFun = NULL)


#### dense case
res.d1 <- FPCA(df.d.Ly,df.d.Lt, optns = list(methodSelectK = 'FVE', methodXi = 'CE'))
res.d2 <- FPCA(df.d.Ly,df.d.Lt, optns = list(methodSelectK = 'FVE', methodXi = 'IN'))
xiEst.compare.d <- abs(res.d1$xiEst-res.d2$xiEst)
#### sparse case
res.s1 <- FPCA(df.s$Ly,df.s$Lt, optns = list(methodSelectK = 'FVE',maxK = K, methodXi = 'CE'))
res.s2 <- FPCA(df.s$Ly,df.s$Lt, optns = list(methodSelectK = 'FVE',maxK = K, methodXi = 'IN'))
xiEst.compare.s <- abs(res.s1$xiEst-res.s2$xiEst)

png("q1bc.png",width = 2200, height = 800, units = "px", pointsize = 30)
par(mfrow= c(1,2))
matplot(xiEst.compare.d, type = 'l', ylim = c(0,3), xlab = "i-th sample", main = "Dense: difference in xiEst between BLUP and Integration")
matplot(xiEst.compare.s, type = 'l', ylim = c(0,3), xlab = "i-th sample", main = "Sparse: difference in xiEst between BLUP and Integration")
dev.off()

res.compare <- data.frame(rbind(res.d1$lambda,res.d2$lambda,res.s1$lambda,res.s2$lambda))
rownames(res.compare) <- c("dense, CE", "dense, IN", "sparse, CE", "sparse, IN")
res.compare


###  Q2   ###
## Extend the PACE method for bivariate sparsely observed longitudinal data
dat <- read.delim('http://www.statsci.org/data/oz/wallaby.txt') %>%
  select(Age, Weight, Leng, Anim) %>%  # select takes the columns
  mutate(Age = Age / 365.24, Weight=Weight / 1e4, Leng = Leng / 1e3) %>%
  na.omit %>%
  filter(Age >= truncAge[1] & Age <= truncAge[2]) %>%
  dplyr::rename(t = Age, y1 = Weight, y2 = Leng,   id = Anim) 
ngrid <- 50 
hmu <- 0.1 

resMu1 <- locfit(y1 ~ lp(t, deg=1, h=hmu), dat, ev=lfgrid(mg=ngrid)) 
resMu2 <- locfit(y2 ~ lp(t, deg=1, h=hmu), dat, ev=lfgrid(mg=ngrid)) 
tGrid <- seq(min(df$t), max(df$t), length.out=ngrid)  # working grid
muHat1 <- predict(resMu1) 
muHat2 <- predict(resMu2) 

par(mfrow=c(1,2))
plot(dat$t, dat$y1) 
lines(tGrid, muHat1) 
plot(dat$t, dat$y2) 
lines(tGrid, muHat2) 
persp(dat$t, dat$y1, dat$y2)

## Data Normalizing
muObs1 <- predict(resMu1, dat$t) # interpolate 
muObs2 <- predict(resMu2, dat$t) # interpolate 
dat$y1.norm <- (dat$y1 - min(dat$y1))/diff(range(dat$y1))
dat$y2.norm <- (dat$y2 - min(dat$y2))/diff(range(dat$y2))
# dat$y1.norm <- scale((dat$y1))
# dat$y2.norm <- scale(dat$y2)

## From the paper, we can stack the dataset and retrieve the principal components for each random variable
## by looking for the corresponding entries.

## Dataframe stacking
#### Note the time for y_2 is shifted, otherwise error of duplicate time t would occur
#### Also note the shift in time t would not affect the estimation of mean function, covariance function or eigen function. 
#### It would only affect the FPC scores.
df <- data.frame('t' = c(dat$t, dat$t+diff(range(dat$t))), 'y' =  c(dat$y1, dat$y2), 'id' =c(dat$id,dat$id))
samp <- MakeFPCAInputs(df$id, df$t, df$y) 
res <- FPCA(samp$Ly, samp$Lt, list(userBwMu=0.1, userBwCov=0.2, FVEthreshold=1)) 
class(res) 
plot(res) 
res$lambda
res$xiEst




## Covariance estimation
#### Get raw cov 
muObs <- predict(resMu, df$t) # interpolate 
df$yCenter <- df$y - muObs 
rcov <- plyr::ddply(df, 'id', function(d) { 
  raw <- c(tcrossprod(d$yCenter)) 
  TT <- expand.grid(t1=d$t, t2=d$t) 
  cbind(TT, raw=raw) 
}) 

#### fit raw covariance
library(rgl) 
plot3d(rcov$t1, rcov$t2, rcov$raw) 
#### bandwidth for covariance estimation
hCov <- 0.2 
resCov <- locfit(raw ~ lp(t1, t2, deg=1, h=hCov), rcov %>% filter(t1 != t2), ev=lfgrid(mg=ngrid), kern='epan') 
covMat <- matrix(predict(resCov), ngrid, ngrid) 
persp3d(tGrid, tGrid, covMat, add=FALSE, col='white') 
#### find the eigen value and functions of fitted raw covariances
eig <- eigen(covMat) 
eig$values 
#### remove the negative eigen values
rmInd <- eig$values <= 0 
eig$values <- eig$values[!rmInd] 
eig$vectors <- eig$vectors[, !rmInd, drop=FALSE] 
lam <- eig$values * diff(range(tGrid)) / ngrid 
phi <- eig$vectors / sqrt(diff(range(tGrid)) / ngrid) 
matplot(tGrid, phi[, 1:3], type='l') 
#### diagnonal elemenet 
diagRes <- locfit(raw ~ lp(t1, deg=1, h=hCov),  
                  rcov %>% filter(t1 == t2),  # %>% passes the result to the outside function, 
                  ev=lfgrid(mg=ngrid)) 
Vhat <- predict(diagRes) 
lines3d(tGrid, tGrid, Vhat) 
sig2 <- mean(Vhat - diag(covMat)) 

## BLUP step 
K <- 3 
LambdaK <- diag(lam[seq_len(K)], nrow=K) 
xi <- plyr::ddply(df, 'id', function(d) { 
  #browser() 
  str(d) 
  # Linear interpolation 
  Phii <- matrix(fdapace::ConvertSupport(tGrid, d$t, phi=phi[, seq_len(K), drop=FALSE]), ncol=K) 
  SigYi <- fdapace::ConvertSupport(tGrid, d$t, Cov=covMat) + diag(sig2, nrow=nrow(d)) 
  res <- c(LambdaK %*% t(Phii) %*% solve(SigYi, d$yCenter)) 
  names(res) <- paste0('xi', seq_len(K)) 
  res 
}) 


plot(xi$xi1, xi$xi2) 

##### X_K 
XK <- as.matrix(xi[, -1]) %*% t(phi[, seq_len(K), drop=FALSE]) + matrix(muHat, nrow=nrow(xi), ncol=ngrid, byrow=TRUE) 
plotInd <- c(1:3, 35) 
matplot(tGrid, t(XK[plotInd, ]), type='l') 
lines(tGrid, muHat, lwd=3) 
points(dat$t[dat$id == 53], dat$y[dat$id == 53], cex=2, col='green') 
points(dat$t[dat$id == 127], dat$y[dat$id == 127], cex=2, col='cyan') 

plot(xi$xi1, xi$xi2) 
points(xi$xi1[plotInd], xi$xi2[plotInd], type='p', col=plotInd, cex=4) 


## Derivative of mean function
FPC <- FPCAder(res, list(method='FPC1', bwMu=0.2, bwCov=0.2)) 
DPC <- FPCAder(res, list(method='DPC', bwMu=0.2, bwCov=0.2)) 
DPCFVE <- cumsum(DPC$lambdaDer) / sum(DPC$lambdaDer) 
xGrid <- DPC$workGrid 
plot(xGrid, DPC$mu, type='l', main='mu') 
plot(xGrid, DPC$muDer, type='l', main='mu\'') 

phiShow <- res$phi[, seq_len(K)] 
matplot(xGrid, phiShow, type='l', main='phi_k') 
phiPrime <- FPC$phiDer[, seq_len(K)] 
matplot(xGrid, phiPrime, main='phi\'', type='l') 
phi1 <- DPC$phiDer[, seq_len(K)] 
matplot(xGrid, phi1, main='phi1', type='l') 

allID <- unique(trajDat$id) 
plotInd <- c(6, 11, 13, 29) 
trajDat1 <- filter(trajDat, id %in% allID[plotInd]) 
fitFPC <- fdapace:::fitted.FPCA(FPC, K=K) 
fitFPCder <- fitted(FPC, K=K) # S3 function 
fitDPC <- fitted(DPC, K=K) 
plotDat <- rbind( 
  cbind(type='FPC', reshape2::melt(fitFPC[plotInd, ], varnames=c('i', 'tInd'))),  
  cbind(type='FPCder', reshape2::melt(fitFPCder[plotInd, ], varnames=c('i', 'tInd'))),  
  cbind(type='DPC', reshape2::melt(fitDPC[plotInd, ], varnames=c('i', 'tInd')))  
) %>% mutate(t = FPC$workGrid[tInd], id = allID[plotInd[i]]) 



p1 <- ggplot(trajDat1, aes(x=t, y=y)) + geom_point() + geom_line() + geom_line(aes(x=t, y=value), data=plotDat %>% filter(type == 'FPC')) + facet_wrap(~id, ncol=1) + ggtitle('X(t)') 
p2 <- ggplot(plotDat %>% filter(type %in% c('FPCder', 'DPC')), aes(x=t, y=value, group=type, color=type)) + geom_line() + facet_wrap(~id, ncol=1) + ggtitle('X\'(t)') 



p <- cowplot::plot_grid(p1, p2, rel_widths=c(1, 1.5), rows=1) 

print(p) # DPC removes redundant information; shrinkage may be better 
