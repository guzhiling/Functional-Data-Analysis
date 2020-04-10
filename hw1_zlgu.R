 # install.packages('tidyverse')    
# install.packages('rgl')
# Sys.setenv(TZ="America/Chicago")
# wd <- "/Users/gzl/Dropbox/STAT547/hw1"
# setwd(wd)
getwd()
library(fda)
library(tidyverse)
library(reshape2)
library(rgl)


########################
########   Q1   ########
########################


png('Q1.png',width = 2000, height = 1200, units = "px", pointsize = 40)

# (a) 
# Input a vector, return a fd object
nobs <- nrow(pinch)
x <- ((1:nobs)-1)*2
nrep <- ncol(pinch)

bspline4 <- function(y) { 
  basis <- create.bspline.basis(rangeval =c(min(x), max(x)), nbasis = 15, norder = 4)
  fdobj <- smooth.basis (x, y,  fdParobj=basis)
}


# Store the fd objects in fd_list
fd_list <- list()
# Color the curves gradiantly
col <- sapply(1:nrep, function(i){rgb (1,0.8*i/nrep, 0.8*i/nrep, alpha =0.7)})

  par(mar = c(5, 4, 4, 4) + 0.3)  # Leave space for z axis
    # PLOT
    plot(1, type="n", xlab="", ylab="", xlim=range(x), ylim=range(pinch))
    for (i in 1:nrep){
      y <- pinch[,i]
      fd_i <- bspline4(y)
      fd_list <- append(fd_list, fd_i)
      lines(fd_i,lwd=2,col = col[i])
    }


# (b) 
y.mean <- apply(pinch, 1, mean)
# y.mean <- colMeans(pinch)
y.sd <- apply(pinch , 1, sd)
mean_col <- 'black'
sd_col <- 'green'

    # PLOT
    lines(x,y.mean, bg = mean_col, col = mean_col, cex =0.6, axes = FALSE, lwd=3)
    par(new = TRUE)
    plot(x, y.sd, type = "l",  bg = sd_col, col = sd_col, cex =0.5, axes = FALSE, xlab = '', ylab = '',xlim=range(x), ylim=range(y.sd), lwd=2)
    axis(side = 4, at = pretty(range(y.sd)))
    axis(side = 2, at = pretty(range(y.mean)))  
    # axis(side = 1, at = pretty(range(x))) 
    mtext("Pinch force &",  side=2, line=3, col = 'black',   font =2 )
    mtext("Pointwise Mean",side=2,line=2,  col = mean_col, font =2)
    mtext("Pointwise SD", side=4, line=2.5, col = sd_col,   font =2)
    mtext("Milliseconds", side=1, line=2.5, col = 'black',   font =2)
    legend('topright', legend = c(sapply ((1:ncol(pinch)), function(x){ paste("",x)})), col = c(col), text.col =c(col), lwd=2,cex = 0.6)
    legend("topleft", legend =c("Mean","SD"), 
           col=c(mean_col, sd_col),text.col =c(mean_col, sd_col), lty =1, lwd=4)
    title(main = "Pinch Force through time for 20 replications")
# title(ylab= "Pinch force", col.lab=rgb(0,0,0))
# title(xlab= "Milliseconds", col.lab=rgb(0,0,0),   font =2)
# box()
dev.off()



# (c)
t.pinch <- t(pinch)  # each row is an observation
df <- data.frame(t.pinch)
colnames(df) <- x   # change the column units to milliseconds
head(df)
Ghat <- cov(df)
matplot(x = x, y = Ghat, type='l', xlab="Milliseconds", ylab= "covariance columns") 
rho <- seq(0, 1, length.out=ncol(df))
persp3d(rho, rho, Ghat, col='white', xlab='s', ylab='t' , zlab = 'Sample Cov',main= "Perspective of Sample Covariance")
contour(rho, rho, Ghat, xlab='s', ylab='t',  main= "Contour Plot of Sample Covariance")

# (d) 
Sig <- Ghat
eig <- eigen(Sig)
FVE <- cumsum(eig$values) / sum(eig$values)
plot(FVE[1:4], type='b', axes = FALSE, ann = FALSE, sex = 0.6) # scree plot
  axis(side =1 , at = 1:4)
  axis(side =2 , at = pretty(range(FVE)))
  title (xlab = 'i-th Eigen Value')
  title (ylab = 'Variance explained')
  box()
  abline(0.9,0,lwd =1, col = 'red', lty=2)
  text(c(0.91,1),labels = "0.9 Cutoff", col = 'red', pos = 4)

for (i in 1:length(FVE)){
  if(FVE[i] >= 0.9){
    flag = i
    break
  }
}
cat(paste('The', flag, "-th eigen value explains 90% of variations" ))



########################
########   Q2   ########
########################
getwd()
d <- read.csv("./DataSets/Dow_companies_data.csv")
data = d[,2:31]

# (a)
xom <- d$XOM
str(xom)
r <- sapply (xom, function(x){(x-xom[1])/xom[1]})*100
cat("Exxon–Mobil stock price increase in 2013 has increased by ", r[length(r)], " percent in return.")

# (b) 
# Create a 252 by 30 matrix cr that contains the cumulative returns on all stocks in DJIA. 
# Plot all the cumulative return functions in one plot. Add the mean and median functions in color.

cr <- apply(d[,2:ncol(d)], 2, function(x){(x - x[1])/x[1]}) *100
rownames(cr) <-  1:252  # d[1:252,1]
colnames(cr) <-  1:30   # colnames(d)[2:31]

median <- apply (cr,1, median); mean_col = 'black'
mean <- apply(cr,1, mean); median_col = rgb(0,0,1,alpha =0.5)

  png('Q2b.png',width = 2000, height = 1200, units = "px", pointsize = 40)
  plot(1, type="n", xlab="", ylab="", axes = FALSE, xlim = c(1,252), ylim =c(-20,90))
  for (i in 1:ncol(cr)){
    lines(rownames(cr), cr[,i], col = rgb(1, 1*i/ncol(cr),1*i/ncol(cr)))
  }
  lines(rownames(cr), mean, col = mean_col,lwd =3)
  lines(rownames(cr), median, col = median_col,lwd =3,lty =1)
  
  axis(side =1, at = seq(1,252,by=10))
  axis(side =2 , at = seq(-20,90,by=10))
  box()
  legend("topleft", legend =c("Mean","Median"), 
         col=c(mean_col, median_col),text.col =c(mean_col, median_col), lty =1, lwd=4)
  title(main = "Cumulative return plot for 30 stocks")
  dev.off();

# (c)
prob <- c(0.9, 0.6, 0.5, 0.3)
prob_col <- sapply(prob, function(x) {rgb(1,x,x )})

png('Q2c.png',width = 2000, height = 1200, units = "px", pointsize = 40)
  fbplot(fit = cr, ylim = c(-20,90), xlab="Trading day", ylab = "Cumulative return",prob=prob, 
                          col=prob_col, axes=FALSE)
  axis(side =1, at = seq(1,252, by =10))
  axis(side =2 , at = seq(-20,90,by=10))
  title(main = "Functional Boxplot for 30 stocks")
  box()
  legend("topleft", legend =c(0.9, 0.6, 0.5, 0.3), 
         col=prob_col,text.col = prob_col, lty =1, lwd=4)
dev.off()



########################
########   Q3   ########
########################



