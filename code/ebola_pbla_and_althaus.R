#######################################
### Analysis of CDC Ebola data with ###
###    PBLA and Althaus methods     ###
#######################################

# Install any missing and required packages
list.of.packages <- c("chron", "deSolve", "bbmle", "RColorBrewer", "ggplot2", "reshape2")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

library(chron)
library(deSolve)
library(bbmle)
library(RColorBrewer)
library(ggplot2)
library(reshape2)

myPalette <- colorRampPalette(rev(brewer.pal(11, "Spectral")), space="Lab")


##################################################

# Read in and set up the data
fdata<- read.csv("../data/Ebola/ebola_cdc_deathsdata.csv")
fdata[,1] <- chron(as.character(fdata[,1]),format=c(dates = "d/m/y"))

##set start dates
gbegin <- chron("2 Dec 2013", format=c(dates = "day mon year")) 
gtimes <- as.numeric(fdata[,1] - gbegin)
slbegin <- chron("27 May 2014", format=c(dates = "day mon year")) - 58
sltimes <- as.numeric(fdata[,1] - slbegin)
lbegin <- chron("27 Mar 2014", format=c(dates = "day mon year")) - 0
ltimes <- as.numeric(fdata[,1] - lbegin)

## transform data into individual based form
gdeaths<-0
j<-1
for (i in 1:length(fdata[,2])){
  if (fdata[i,2]>0){
    while (j<=fdata[i,2]){
      gdeaths[j]<-gtimes[i]
      j=j+1
    }
  }
}
sldeaths<-0
j<-1
for (i in 1:length(fdata[,4])){
  if (fdata[i,4]>0){
    while (j<=fdata[i,4]){
      sldeaths[j]<-sltimes[i]
      j=j+1
    }
  }
}
ldeaths<-0
j<-1
for (i in 1:length(fdata[,3])){
  if (fdata[i,3]>0){
    while (j<=fdata[i,3]){
      ldeaths[j]<-ltimes[i]
      j=j+1
    }
  }
}

## Spread the data out - over the whole time period since the last observation
for(i in 58:1){
  gdeaths[i]<-gdeaths[i+1]-0.1
}
i<-60
while(i<length(gdeaths)){
  is<-i
  while (gdeaths[i]==gdeaths[i+1]){
    i=i+1
  }
  am<-(gdeaths[i]-gdeaths[is-1])/(i-is+1)
  for (j in (is):i){
    gdeaths[j]<-gdeaths[j-1]+am
  }
  i=i+1
}
for(i in 3:1){
  sldeaths[i]<-sldeaths[i+1]-0.1
}
i<-5
while(i<length(sldeaths)){
  is<-i
  while (sldeaths[i]==sldeaths[i+1]){
    i=i+1
  }
  am<-(sldeaths[i]-sldeaths[is-1])/(i-is+1)
  for (j in (is):i){
    sldeaths[j]<-sldeaths[j-1]+am
  }
  i=i+1
}
for(i in 5:1){
  ldeaths[i]<-ldeaths[i+1]-0.1
}
i<-7
while(i<length(ldeaths)){
  is<-i
  while (ldeaths[i]==ldeaths[i+1]){
    i=i+1
  }
  am<-(ldeaths[i]-ldeaths[is-1])/(i-is+1)
  for (j in (is):i){
    ldeaths[j]<-ldeaths[j-1]+am
  }
  i=i+1
}

# Scale data - helps numerically (will be scaled back up after)
gdeaths<-gdeaths/1000
sldeaths<-sldeaths/1000
ldeaths<-ldeaths/1000


##################################################
## Log likelihood 
PBLA_SEIR<-function(x,gt,N,r,c){
  
  # parameters
  b0 <- x[1]
  k <- x[2]
  nI <- length(r) # number infected
  g <- gt
  
  b<-matrix(NA, nI, nI)
  for (j in 1:nI){
    for (i in 1:nI){
      if (r[j]<(r[i]-1/c)){
        t<- r[j]-1/(2*g)-(3/(4*g))*exp(-g*(r[i]-r[j]-1/c))
      }
      else{t<-r[i]-1/c-1/g-(1/(4*g))*exp(-g*(r[j]-r[i]+1/c)) }
      #t<-min(r[j],r[k]-1/g-1/c)
      b[j,i] <- b0*exp(-k*(t))/N
    }
  } 
  
  
  #only need i>nI for B calculation below:
  B<-c(rep(0,nI))
  for (j in 1:nI){
    #for (i in (nI+1):N){
    B[j] <- (N-nI)*b0*exp(-k*(r[j]-1/(2*g)))/N     
    #}
  }
  d <- g + B
  
  if (b<0 || g<0){like <- -1000000}
  else{
    
    # Prior for initial infective
    # pri_initial <- rep(1/N,nI)
    pri_initial <- rep(0,nI)
    pri_initial[1] <- 1
    chi_phi <- rep(0,nI)
    psi <- rep(1,nI)
    like <- 0
    # Product term, product of chi_phi[j] psi[j]
    for (j in 1:nI){
      for  (i in 1:nI){      
        if (i!=j){        
          if ((r[j]-(1/c))<r[i]){          
            denom <- 1 - (d[j]*b[i,j]/( (b[i,j]+d[i])*(d[j]+d[i]) ))*exp(-d[i]*(r[i]-(r[j]-(1/c))))
            chi_phi[j] <- chi_phi[j] + b[i,j]*exp(-d[i]*(r[i]-(r[j]-(1/c))))*(d[j]*d[i]/( (b[i,j]+d[i])*(d[j]+d[i])*denom))
          }
          else{
            denom <- ( d[i]/((d[i])+b[i,j]) ) + (d[i]*b[i,j]/( (b[i,j]+d[i])*(d[j]+d[i]) ))*exp(-d[j]*((r[j]-(1/c))-r[i]))
            chi_phi[j] <- chi_phi[j] + b[i,j]*exp(-d[j]*((r[j]-(1/c))-r[i]))*(d[j]*d[i]/( (b[i,j]+d[i])*(d[j]+d[i])*denom))
          }
          
        }  
      }
      #chi_phi[j] <- chi_phi[j]
      
      for (L in 1:nI){
        if (L!=j){ 
          if ((r[j]-(1/c))<r[L]){
            psi[j] <- psi[j]*(1 - (d[j]*b[L,j]/( (b[L,j]+d[L])*(d[j]+d[L]) ))*exp(-d[L]*(r[L]-(r[j]-(1/c)))))
          }
          else{
            psi[j] <- psi[j]*( ( d[L]/((d[L])+b[L,j]) ) + (d[L]*b[L,j]/( (b[L,j]+d[L])*(d[j]+d[L]) ))*exp(-d[j]*((r[j]-(1/c))-r[L])) )
          }
        }
      }
      like <- like + log(chi_phi[j]) + log(psi[j]) 
    }
    
    # sum term
    # Assumes prior for initial infection time is uniform on (-infinity, r[1])
    s <- 0 
    for (L in 1:nI){
      if (pri_initial[L] > 0){    
        
        s <- s + exp( -d[L]*((r[L])-r[1]) + log(pri_initial[L]) - log(chi_phi[L]) - log(psi[L]) )
        #s <- s + exp( -d*(((r[k]-c)-r[1])*d) + log(pri_initial[k]) - log(chi_phi[k]) - log(psi[k]) )
      }
    }
    
    like <- like + log(s) + sum(log(g) - log(d))
  }
  return(like)
}



##################################################
# Analyses 


#################################################
# Althaus analysis - deaths only 
# (code adapted from https://github.com/calthaus/Ebola)

# Definition of the SEIR model
SEIR <- function(t, x, parms) {
  with(as.list(c(parms,x)),{
    if(t < tau1) beta <- beta0
    else beta <- beta0*exp(-k*(t-tau1))
    N <- S + E + I + R
    dS <- - beta*S*I/N
    dE <- beta*S*I/N - sigma*E
    dI <- sigma*E - gamma*I
    dR <- 0 #(1-f)*gamma*I
    dD <- f*gamma*I
    dC <- sigma*E
    der <- c(dS,dE,dI,dR,dD,dC)
    list(der)
  })
}

# Negative log-likelihood
nll <- function(beta0,k,f,tau0,tau1,sigma,gamma) {
  pars <- c(beta0=beta0,k=k,f=f,tau0=tau0,tau1=tau1,sigma=sigma,gamma=gamma)
  pars <- trans(pars)
  times <- c(0,data$times+pars["tau0"])
  simulation <- as.data.frame(ode(init,times,SEIR,parms=pars))
  simulation <- simulation[-1,]
  ll <- sum(dpois(data$deaths,simulation$D,log=TRUE))
  return(-ll)
}

# Parameter transformation
trans <- function(pars) {
  pars["beta0"] <- exp(pars["beta0"])
  pars["k"] <- exp(pars["k"])
  pars["f"] <- plogis(pars["f"])
  pars["tau0"] <- exp(pars["tau0"])
  pars["tau1"] <- exp(pars["tau1"])
  return(pars)
}

# GUINEA: Prepare the data, set the initial values and fit the model
data <- fdata[25:156,1:2]
names(data) <- c("times","deaths")
begin <- chron("2 Dec 2013", format=c(dates = "day mon year")) 
delay <- as.numeric(data$times[1] - begin)
data$times <- data$times - data$times[1]
N <- 1e6		
init <- c(S = N - 1, E = 0, I = 1, R = 0, D = 0, C = 1)
fixed <- c(tau0 = log(delay), tau1 = -Inf, sigma = 1/5.3, gamma = 1/5.61, f=20)
free <- c(beta0 = log(0.2), k = log(0.001))
fit <- mle2(nll,start=as.list(free),fixed=as.list(fixed),method="Nelder-Mead",control=list(maxit=1e3))
trans(coef(fit))


# SIERRA LEONE: Prepare the data, set the initial values and fit the model
data <- fdata[40:156,c(1,4)]
names(data) <- c("times","deaths")
begin <- min(data$times) 
data$times <- data$times - data$times[1]
N <- 1e6		
init <- c(S = N - 1, E = 0, I = 1, R = 0, D = 0, C = 1)
fixed <- c(tau1 = -Inf, sigma = 1/5.3, gamma = 1/5.61,f = 30)
free <- c(beta0 = log(0.2), k = log(0.001), tau0 = log(60))
fit <- mle2(nll,start=as.list(free),fixed=as.list(fixed),method="Nelder-Mead",control=list(maxit=1e3))
trans(coef(fit))


# LIBERIA: Prepare the data, set the initial values and fit the model - free k
data <- fdata[27:156,c(1,3)]
names(data) <- c("times","deaths")
begin <- min(data$times) 
data$times <- data$times - data$times[1]
N <- 1e6		
init <- c(S = N - 1, E = 0, I = 1, R = 0, D = 0, C = 1)
fixed <- c(tau1 = -Inf, sigma = 1/5.3, gamma = 1/5.61, f = 30)
free <- c(beta0 = log(0.2), k = log(0.001),  tau0 = log(60))
fit <- mle2(nll,start=as.list(free),fixed=as.list(fixed),method="Nelder-Mead",control=list(maxit=1e3))
trans(coef(fit))


#################################################
# Now perform PBLA analysis

# GUINEA

# Althaus method results, with CDC data, scaled back up
b0<-0.2306*10^(3)
k<-0.0007118*10^(3)
g<-1/(5.61*10^(-3))
c<-1/(5.3*10^(-3))
N<-10^6
x<-c(b0,k)

#  get PBLA MLEs
optg<-optim(par=c(200,1), fn=PBLA_SEIR, gt=g, N=N, r=gdeaths, c=c, control=list("fnscale"=-1, trace=3), hessian=FALSE)$par
# check the result:
optg

# and PBLA profile log-likelihoods (with other parameters fixed to MLE values)
#b0#
marggb<-0
for (j in 1:50){
  marggb[j]<-PBLA_SEIR(c(j*10,optg[2]),g,N,gdeaths,c)
}

#k#
marggk<-0
for (j in 0:50){
  marggk[j]<-PBLA_SEIR(c(optg[1],j/10),g,N,gdeaths,c)
}

## contour plots for b0, k
no<-30
matg<-matrix(NA, no,no)
for (i in 0:no){
  for (j in 0:no){
    print(j)
    matg[i,j]<-PBLA_SEIR(c(10*i+100,j/6),g,N,gdeaths,c)
  }
}
filled.contour(x = seq(100, 400, length.out = no)/1000, y = seq(0, 5, length.out = no)/1000,matg,col=myPalette(20),
               plot.axes={axis(1,seq(100,400,length.out=11)/1000)
                 axis(2, seq(0,5, length.out=11)/1000) 
              points(b0/1000,k/1000, pch=16, cex=1.5)}, main="Contour plot for Guinea deaths CDC data")
title( xlab=expression('b'[0]), ylab="k")


# SIERRA LEONE

# Althaus method results, with CDC data, scaled back up
b0<-0.2588*10^(3)
k<-0.001310*10^(3)
g<-1/(5.61*10^(-3))
c<-1/(5.3*10^(-3))
N<-10^6
x<-c(b0,k)

#  get PBLA MLEs
opts<-optim(par=c(200,10), PBLA_SEIR, gt=g, N=N, r=sldeaths, c=c, control=list("fnscale"=-1, trace=3), hessian=FALSE)$par
# check result:
opts

# and PBLA profile log-likelihoods (with other parameters fixed to MLE values)
#b0#
margsb<-0
for (j in 1:50){
  margsb[j]<-PBLA_SEIR(c(j*20,opts[2]),g,N,sldeaths,c)
}

#k#
margsk<-0
for (j in 0:50){
  margsk[j]<-PBLA_SEIR(c(opts[1],j/5),g,N,sldeaths,c)
}

# contour plot
no<-30
mats<-matrix(NA, no,no)
for (i in 0:no){
  print(i)
  for (j in 0:no){
    mats[i,j]<-PBLA_SEIR(c(50*i,j/1.5),g,N,sldeaths,c)
  }
}
filled.contour(x = seq(0, 1500, length.out = no)/1000, y = seq(0, 20, length.out = no)/1000,mats,col=myPalette(20),plot.axes={axis(1,seq(0,1500,length.out=16)/1000)
  axis(2, seq(0,20, length.out=11)/1000) 
  points(b0/1000,k/1000)}, main="Contour plot for Sierra Leone deaths CDC data")
title( xlab=expression('b'[0]), ylab="k")


# LIBERIA

# Althaus method results, with CDC data, scaled back up
b0<-0.2807*10^(3)
k<-0.001219*10^(3)
g<-1/(5.61*10^(-3))
c<-1/(5.3*10^(-3))
N<-10^6
x<-c(b0,k)

#  get PBLA MLEs
optl<-optim(par=c(100,0.0000001), PBLA_SEIR, gt=g, N=N, r=ldeaths, c=c, control=list("fnscale"=-1, trace=3), hessian=TRUE)$par
# check result:
optl

# and PBLA profile log-likelihoods (with other parameters fixed to MLE values)
#b0#
marglb<-0
for (j in 1:50){
  marglb[j]<-PBLA_SEIR(c(j*20,optl[2]),g,N,ldeaths,c)
}
#k#
marglk<-0
for (j in 0:50){
  marglk[j]<-PBLA_SEIR(c(optl[1],j/5),g,N,ldeaths,c)
}

# contour plot
no<-30
matl<-matrix(NA, no,no)
for (i in 0:no){
  print(i)
  for (j in 0:no){
    matl[i,j]<-PBLA_SEIR(c(i*100,j),g,N,ldeaths,c)
  }
}
filled.contour(x = seq(0, 3000, length.out = no)/1000, y = seq(0, 30, length.out = no)/1000,matl,col=myPalette(30),plot.axes={axis(1,seq(0,3000,length.out=7)/1000)
  axis(2, seq(0,30, length.out=16)/1000) 
  points(b0/1000,k/1000)}, main="Contour plot for Liberia deaths CDC data")
title( xlab=expression('b'[0]), ylab="k")


#################################################
# Write/read necessary files (so we don't need to re-run the optimisation)
write(file="marggb.txt", marggb, ncolumns=1)
write(file="marggk.txt", marggk, ncolumns=1)
write(file="margsb.txt", margsb, ncolumns=1)
write(file="marglb.txt", marglb, ncolumns=1)
write(file="marglk.txt", marglk, ncolumns=1)

marggb<-read.table(file="marggb.txt")[,1]
marggk<-read.table(file="marggk.txt")[,1]
margsb<-read.table(file="margsb.txt")[,1]
margsk<-read.table(file="margsk.txt")[,1]
marglb<-read.table(file="marglb.txt")[,1]
marglk<-read.table(file="marglk.txt")[,1]



#################################################
# Finally, prepare plots of results

# (Profile likelihoods, with parameters scaled back, 
# and overlaid Althaus MLEs)

#GUINEA

# Althaus method results, with CDC data, scaled back up
b0<-0.2306*10^(3)
k<-0.0007118*10^(3)
g<-1/(5.61*10^(-3))
c<-1/(5.3*10^(-3))
N<-10^6
x<-c(b0,k)

par(mfrow=c(1,2),oma = c(0, 0, 1, 0))
#b0#
x<-seq(10,500,length.out=50)
plot(x/1000,marggb, type="l", ylab="log likelihood", xaxt="n", xlab="")
axis(1)
mtext(expression('b'[0]), side=1, line=2.8, cex=1.2)
abline(v=b0/1000, col="black", lty=2)

#k#
x<-seq(0,5,length.out=50)
plot(x/1000,marggk, type="l", ylab="log likelihood", xaxt="n", xlab="")
axis(1)
mtext("k", side=1, line=2.8, cex=1.2)
abline(v=k/1000, col="black", lty=2)

mtext("Guinea", outer = TRUE, line=-2, cex = 1.5)


# SIERRA LEONE

# Althaus method results, with CDC data, scaled back up
b0<-0.2588*10^(3)
k<-0.001310*10^(3)
g<-1/(5.61*10^(-3))
c<-1/(5.3*10^(-3))
N<-10^6
x<-c(b0,k)

#b0#
x<-seq(20,1000,length.out=50)
plot(x/1000,margsb, type="l", ylab="log likelihood", xaxt="n", xlab="")
axis(1)
mtext(expression('b'[0]), side=1, line=2.8, cex=1.2)
abline(v=b0/1000, col="black", lty=2)

#k#
x<-seq(0,10,length.out=50)
plot(x/1000,margsk, type="l", ylab="log likelihood", xaxt="n", xlab="")
axis(1)
mtext("k", side=1, line=2.8, cex=1.2)
abline(v=k/1000, col="black", lty=2)

mtext("Sierra Leone", outer = TRUE, line=-2, cex = 1.5)

# LIBERIA

# Althaus method results, with CDC data, scaled back up
b0<-0.2807*10^(3)
k<-0.001219*10^(3)
g<-1/(5.61*10^(-3))
c<-1/(5.3*10^(-3))
N<-10^6
x<-c(b0,k)

#b0#
x<-seq(20,1000,length.out=50)
plot(x/1000,marglb, type="l", ylab="log likelihood", xaxt="n", xlab="")
axis(1)
mtext(expression('b'[0]), side=1, line=2.8, cex=1.2)
abline(v=b0/1000, col="black", lty=2)

#k#
x<-seq(0,10,length.out=50)
plot(x/1000,marglk, type="l", ylab="log likelihood", xaxt="n", xlab="")
axis(1)
mtext("k", side=1, line=2.8, cex=1.2)
abline(v=k/1000, col="black", lty=2)

mtext("Liberia", outer = TRUE, line=-2, cex = 1.5)






