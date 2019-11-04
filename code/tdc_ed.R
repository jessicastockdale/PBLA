##################################################
### Code to generate Eichner Dietz (ED) method ###
###   results for the Tristan da Cunha data    ###
##################################################

# Install any missing and required packages
list.of.packages <- c("ggplot2")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

library("ggplot2")

##################################################
## Read in the data (files provided in data folder of github repository)
data<-read.table("../data/TristanDaCunha/tdc_jitteredtimes.txt", header=FALSE)
c<-read.table("../data/TristanDaCunha/tdc_agegroups.txt", header=FALSE)
data <- cbind(data,c)

##################################################
## Set up required variables and functions
datar<-data[1:40,1]
N<-254 # population size
nI <-40 # number infected
g<-0.371 # gamma (removal rate)
beta <- c(0.00451, 0.00181, 0.00131) # beta (infection rate)
theta <- c(0.00451, 0.00181, 0.00131, g) # set of parameters to be estimated

# Trapezium rule
trapz <- function (x, y) {
    idx = 2:length(x)
    return(as.double((x[idx] - x[idx - 1]) %*% (y[idx] + y[idx-1]))/2)
  }


# ED log likelihood (to be minimised in NLM, so made negative)
likelihood_o_min<-function(theta, N, nI, data, K = 1000, init.value = -10){
    
    # all parameters must be non-negative
    if (any(theta<0)){
      return(100000.00)
    }
    
    # set up
    r<-data[1:40,1] # removal times
    type<-data[,2] # age groups
    g<-theta[4] # gamma removal rate
    # number never infected in each group
    m1 <- 16
    m2 <- 30
    m3 <- 168
    
    #contribution from infectives
    lhc <- 0
    for (k in 2:length(r)){
      
      # which beta, from age groupings
      if (type[k]==1){b = theta[1]
      }else if (type[k]==2){b = theta[2]
      }else{b = theta[3]}
      
      # time to integrate over
      x <- seq(init.value, max(r), length.out=K)
      
      # now calculate the contribution...
      sumA <- rep(0,K)   
      for (i in 1:length(x)){
        sumA[i] <- sum(exp(-g*(r-pmin(x[i],r)))[-k]) 
        }
      
      y <- rep(0,K)
      for (j in 1:length(r)){
        if (j!=k){
          choose.index <- which(x < min(r[k], r[j])) # integral limits
          y[choose.index ]<- y[choose.index] + exp(-g*(r[j]+r[k]) + 2*g*x[choose.index] - (b/g)*sumA[choose.index])
        }
      }  
      
      lhc <- lhc + log(b) + log(g) + log(trapz(x,y))
    }  
  
  #contribution from non-infectives
  nilhc <- -(nI/g)*(theta[1]*m1 + theta[2]*m2 + theta[3]*m3)
  
  #overall negative log likelihood
  loglh <- lhc + nilhc
  return(-loglh)  
}

# Since we would like MAPs, we also need to include priors
post<-function(theta, N, nI, data, K = 1000, init.value = -10){
  likelihood_o_min(theta, N, nI, data, K = 1000, init.value = -10) + 
    log(pgamma(theta[1], shape=0.0000001, scale = 1/0.00001)) +
    log(pgamma(theta[2], shape=0.0000001, scale = 1/0.00001)) +
    log(pgamma(theta[3], shape=0.0000001, scale = 1/0.00001)) +
    log(pgamma(theta[4], shape=0.0001, scale = 1/0.001))
}



##################################################
## Now, perform the optimisation

# note: warning 'NA/Inf replaced by maximum positive value' from NLM may occur - unacceptable parameter values just being tested
opt3<-nlm(post, c(0.1,0.1,0.1,0.01), N=N, nI=nI, data=data, typsize = c(0.001, 0.001, 0.001, 0.1), 
          fscale = -225, print.level=2 )
write(opt3$estimate, file ="ed_tdc_opt.txt", sep = ",", ncolumns=1, append = TRUE)
#$estimate
# [1] 0.01040.00408 0.00289 0.879



##################################################
## Create plots comparing PBLA, DA-MCMC and ED 

# Read in the results from PBLA:
beta1<-read.table("beta1.txt")
beta2<-read.table("beta2.txt")
beta3<-read.table("beta3.txt")
gamma<-read.table("gamma.txt")

# and for ED (in case we ran the optimiser in another session):
opt3$estimate <- t(read.table("ed_tdc_opt.txt"))

# We add the Hayakawa et al DA-MCMC results in manually

# Create figures:
pl<-ggplot(data=data.frame(beta1[,1]), aes(beta1[,1])) + 
  geom_histogram(col="gray2", fill="gray2", alpha = .2) +
  ggtitle("") + theme_light() +
  theme(plot.title = element_text(hjust = 0.5, size = 30), axis.text=element_text(size=18,colour="black"),
        axis.title=element_text(size=25), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  labs(x=expression(beta[1]), y="Count") + 
  xlim(c(0,0.03))
pl + geom_vline(xintercept = 0.00451, linetype="solid", 
                color = "black", size=1) + geom_vline(xintercept = opt3$estimate[1], linetype="dashed", 
                                                      color = "black", size=1)


pl<-ggplot(data=data.frame(beta2[,1]), aes(beta2[,1])) + 
  geom_histogram(col="gray2", fill="gray2", alpha = .2) +
  ggtitle("") + theme_light() +
  theme(plot.title = element_text(hjust = 0.5, size = 30), axis.text=element_text(size=18,colour="black"),
        axis.title=element_text(size=25), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  labs(y="Count") + 
  scale_x_continuous(name=expression(beta[2]),breaks=c(0.0,0.005,0.01), labels=c(0.0,0.005,0.01), limits=c(0,0.01))
pl + geom_vline(xintercept = 0.00181, linetype="solid", 
                color = "black", size=1) + geom_vline(xintercept = opt3$estimate[2], linetype="dashed", 
                                                      color = "black", size=1)

pl<-ggplot(data=data.frame(beta3[,1]), aes(beta3[,1])) + 
  geom_histogram(col="gray2",  fill="gray2",  alpha = .2) +
  ggtitle("") + theme_light() +
  theme(plot.title = element_text(hjust = 0.5, size = 30), axis.text=element_text(size=18,colour="black"),
        axis.title=element_text(size=25), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  labs(x=expression(beta[3]), y="Count") + 
  xlim(c(0,0.006))
pl + geom_vline(xintercept = 0.00131, linetype="solid", 
                color = "black", size=1) + geom_vline(xintercept = opt3$estimate[3], linetype="dashed", 
                                                      color = "black", size=1)

pl<-ggplot(data=data.frame(gamma[,1]), aes(gamma[,1])) + 
  geom_histogram(col="gray2", fill="gray2",  alpha = .2) +
  ggtitle("") + theme_light() +
  theme(plot.title = element_text(hjust = 0.5, size = 30), axis.text=element_text(size=18,colour="black"),
        axis.title=element_text(size=25), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  labs(x=expression(gamma), y="Count") + 
  xlim(c(0,2))
pl + geom_vline(xintercept = 0.371, linetype="solid", 
                color = "black", size=1) + geom_vline(xintercept = opt3$estimate[4], linetype="dashed", 
                                                      color = "black", size=1)

## R0
R0_ed<- (opt3$estimate[1]*25 + opt3$estimate[2]*36 + opt3$estimate[3]*193)/opt3$estimate[4]
R0<- (beta1*25 + beta2*36 + beta3*193)/gamma

pl<-ggplot(data=data.frame(R0[,1]), aes(R0[,1])) + 
  geom_histogram(col="gray2", fill="gray2", alpha = .2) +
  ggtitle("") + theme_light() +
  theme(plot.title = element_text(hjust = 0.5, size = 30), axis.text=element_text(size=18,colour="black"),
        axis.title=element_text(size=25), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  labs(x=expression(R[0]), y="Count") + 
  xlim(c(0,2.5))
pl + geom_vline(xintercept = 1.2, linetype="solid", 
                color = "black", size=1) + geom_vline(xintercept = R0_ed, linetype="dashed", 
                                                      color = "black", size=1)











