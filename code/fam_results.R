#######################################
### Create plots of results for the ###
### (simulated) Foot and Mouth data ###
########################################

# install any missing packages
list.of.packages <- c("RColorBrewer", "classInt")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

library(RColorBrewer)
library(classInt)

###############################
## set up 

# set preferred seed
set.seed(12345)

# read in the data
data<-read.table("../data/FMD/data_rd.txt") # removal data
n_c<-read.table("../data/FMD/n_c.txt") # number of cows per farm
n_s<-read.table("../data/FMD/n_s.txt") # number of sheep per farm
xc<-read.table("../data/FMD/x_coord_rd.txt") # x coordinate of farm location
yc<-read.table("../data/FMD/y_coord_rd.txt") # y coordinate of farm location
N<-5378 # population size
n<-1021 # outbreak (final) size

# order the numbers of cows and sheep/ coordinates such that infected farms are first
n_c<-c(n_c[4358:5378,1], n_c[1:4357,1])
n_s<-c(n_s[4358:5378,1], n_s[1:4357,1])
xc<-c(xc[4358:5378,1], xc[1:4357,1])
yc<-c(yc[4358:5378,1], yc[1:4357,1])

total_animals<-n_c+n_s


###############################
## Plot the locations of farms

# by infected/non infected farms
plot(3:4,5:6,type="n", xlab="", ylab="", xaxt="n", yaxt="n", main="Location of farms in 2001 FMD Cumbrian outbreak dataset")
brks<<-classIntervals(c(rep(0,N-n), rep(1,n)), n=2, style="quantile")
brks<- brks$brks
cols<-colors[findInterval(c(rep(0,N-n), rep(1,n)), brks, all.inside=TRUE)]
points(xc+0.05,yc+0.05 , pch=21, bg=cols, cex=0.7, lwd=.4)


###############################
## Plot profile likelihoods

beta0<-read.table("beta0prof.txt")
gamma<-read.table("gammaprof.txt")
v<-read.table("vprof.txt")
ep<-read.table("epprof.txt")
xi<-read.table("xiprof.txt")
ze<-read.table("zeprof.txt")


par(mfrow=c(3,2), mai = c(0.6, 0.6, 0.3, 0.3))


plot(seq(0, 1000/100000000, length.out=1001),beta0[,1], main="", type="l", ylab="log likelihood", xaxt="n", xlab="")
axis(1,cex.axis=1.2)
mtext(expression(beta[0]), side=1, line=2.2, cex=1.1)
abline(v=0.000000607, col="black", lty=1)
abline(v=0.000000710, col="black", lty=2)

plot(seq(1/100, 1, length.out=100),gamma[,1], main="", type="l",  ylab="log likelihood", xaxt="n", xlab="")
axis(1,cex.axis=1.2)
mtext(expression(gamma), side=1, line=2.2, cex=1.1)
abline(v=0.52, col="black", lty=1)
abline(v=0.54, col="black", lty=2)

plot(seq(0, 100/10000, length.out=101)[2:101],v[2:101,1], main="", type="l",  ylab="log likelihood", xaxt="n", xlab="")
axis(1,cex.axis=1.2)
mtext("v", side=1, line=2.2, cex=1.1)
abline(v=0.0065, col="black", lty=1)
abline(v=0.0049, col="black", lty=2)

plot(seq(0.1, 10, length.out=100),ep[,1], main="", type="l",  ylab="log likelihood", xaxt="n", xlab="")
axis(1,cex.axis=1.2)
mtext(expression(epsilon), side=1, line=2.2, cex=1.1)
abline(v=1.45, col="black", lty=1)
abline(v=2.32, col="black", lty=2)

plot(seq(0, 10, length.out=101),xi[,1], main="", type="l", ylab="log likelihood", xaxt="n", xlab="")
axis(1,cex.axis=1.2)
mtext(expression(xi), side=1, line=2.2, cex=1.1)
abline(v=2.32, col="black", lty=1)
abline(v=2.20, col="black", lty=2)

plot(seq(0.01, 1, length.out=100),ze[,1], main="", type="l", ylab="log likelihood", xaxt="n", xlab="")
axis(1,cex.axis=1.2)
mtext(expression(zeta), side=1, line=2.2, cex=1.1)
abline(v=0.32, col="black", lty=1)
abline(v=0.31, col="black", lty=2)




