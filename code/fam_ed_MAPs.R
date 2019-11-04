
##############################################
###   Code to find Eichner Dietz MLEs/MAPs  ###
### for the (simulated) Foot and Mouth data ###
##############################################

###############################
## set up 

# set seed if you wish
set.seed(12345)

print("Read in the data...")

# Load in the data and define necessary quantities
data<-read.table("../data/FMD/data_rd.txt")
n_c<-read.table("../data/FMD/n_c.txt") # number of cows
n_s<-read.table("../data/FMD/n_s.txt") # number of sheep
xc<-read.table("../data/FMD/x_coord_rd.txt") # x co-ordinate of farm
yc<-read.table("../data/FMD/y_coord_rd.txt") # y co-ordinate of farm
N<-5378 # population size
n<-1021 # final outbreak size

# Jitter the data, in a random order for all equal sets of removals
nd<-data[,2]
nd[2]<-nd[2]+0.01

c<-1
for (i in 3:(length(data[,2])+1)){
  if (i!=(length(data[,2])+1)){
    if (data[i,2]==data[i-1,2]){ #add one to counter
      c<-c+1
    }else{ 
      if (c>1){#Jitter previous set, and reset counter
        sam<-sample(c(1:c), size = c, replace = FALSE)
        for (j in 1:c){
          nd[i-c-1+sam[j]]<-data[i-1,2]+0.01*(j-1)
        }
        c<-1}
    }
  }else{
    if (c>1){#Jitter last ones if needed
      sam<-sample(c(1:c), size = c, replace = FALSE)
      for (j in 1:c){
        nd[i-c-1+sam[j]]<-data[i-1,2]+0.01*j
      }
      c<-1}}
}
data[,2]<-nd

# re-order numbers of cows and sheep, so infected farms are first
n_c<-c(n_c[4358:5378,1], n_c[1:4357,1])
n_s<-c(n_s[4358:5378,1], n_s[1:4357,1])
xc<-c(xc[4358:5378,1], xc[1:4357,1])
yc<-c(yc[4358:5378,1], yc[1:4357,1])
total_animals<-n_c+n_s

###############################
## Define required functions/quantities
print("Define required functions/quantities...")

# Euclidean distance matrix
pmat<-matrix(NA,N,N)
for (i in 1:N){
  for (j in 1:N){
    pmat[i,j]<-(xc[i]-xc[j])^2 + (yc[i]-yc[j])^2
  }
}

# trapezium rule
trapz <- function (x, y) {
  idx = 2:length(x)
  return(as.double((x[idx] - x[idx - 1]) %*% (y[idx] + y[idx-1]))/2)
}

# ED negative log posterior, for gamma infectious periods (since NLM minimises)
post_o_min<-function(theta,  m=4, N, data, n_c, n_s, xc, yc, pmat, K = 1000, init.value = -10){
  # k = number of trapezia, init.value = lower limit of numerical integrals, m= shape parameter of gamma distribution
  
  # removal times
  r<-data[,2]
  nI <- length(r)
  
  # parameters to be estimated
  g<-theta[2]
  v<-theta[3]
  ep<-theta[4]
  xi<-theta[5]
  ze<-theta[6]
  
  b<-matrix(NA, N, N)
  b<-theta[1]*(v/(pmat +v^2))
  b<-b*(ep*((n_c)^ze) + (n_s)^ze)
  b<-t(t(b)*(xi*((n_c)^ze) + (n_s)^ze))
  
  # begin log likelihood calculation:
  if (any(theta<0)){ # need non-negative parameter values
    return(10000000.00)
  }
  
  #contribution from infectives
  lhc <- 0
  for (k in 2:length(r)){
    
    x <- seq(init.value, max(r), length.out=K)
    sumC <- rep(0,K)   
    for (i in 1:length(x)){
      sumM<-0
      for (l in 0:(m-1)){
        sumM <- sumM + (1/factorial(l))*(m-l)*(g*(r-pmin(x[i],r)))[-k]^l
      }
      sumC[i] <- sum((b[1:length(r),k]*exp(-g*(r-pmin(x[i],r))))[-k]*(sumM)) # use pmin instead of min for element by element min  
    }
    
    y <- rep(0,K)
    for (j in 1:length(r)){
      
      
      if (j!=k){
        choose.index <- which(x < min(r[k], r[j])) # integral limit
        
        sumM2<-0
        for (l in 0:(m-1)){
          sumM2 <- sumM2 + ((r[j] - x[choose.index])^l)*(g^l)/factorial(l) 
        }
        
        y[choose.index ]<- y[choose.index] + exp(-g*(r[j]+r[k]) + 2*g*x[choose.index] - (1/g)*sumC[choose.index])*b[j,k]*
          ((r[k]-x[choose.index])^(m-1))*(sumM2)
      }
    }  
    
    lhc <- lhc + log(g^m) - log(factorial(m-1)) + log(trapz(x,y))
  }  
  
  #contribution from non-infectives
  nilhc <- -(m/g)*sum(b[1:nI,(nI+1):N])
  
  #overall (negative log) posterior 
  loglh <- lhc + nilhc
  loglh <- loglh + log(pgamma(theta[1], shape=0.001, scale = 1/0.001)) +
    log(pgamma(g, shape=0.001, scale = 1/0.001)) +
    log(pgamma(v, shape=1.0, scale = 1/0.1)) +
    log(pgamma(ep, shape=1.0, scale = 1/0.001)) +
    log(pgamma(xi, shape=1.0, scale = 1/0.001)) +
    log(pgamma(ze, shape=1.0, scale = 1/0.001))
  return(-loglh)  
  
}
  
###############################
## Now begin analysis

print("Now run NLM...")  
  
## Use non-linear minimization (found best for high-dimensional lhs) to find MLEs of 6 model parameters

init.values <- c(1.404e-06, 8.900e-01, 9.540e-03, 3.152e+00, 4.778e+00, 6.360e-01)
opt2<-nlm(post_o_min, init.values,pmat=pmat, m=4, N=N, data=data, n_c=n_c, n_s=n_s,xc=xc,yc=yc,
          K=500,init.value = -10, typsize = c(0.0000001,0.1,0.001,1.0,1.0,0.1), 
          fscale = 6000, print.level=2 )


  
