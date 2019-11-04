rm(list=ls())

args=(commandArgs(TRUE))

############################################################
###   Code to generate PBLA MAPs and profile likelihoods ###
###        for the (simulated) Foot and Mouth data       ###
############################################################

# (generating both MAPs and profile likelihoods for each parameter can be quite
# time consuming, so in reality we might want to parallelize these)

# install any missing packages
list.of.packages <- c("Rcpp", "RcppGSL")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

library(Rcpp)
library(RcppGSL)

###############################
## set up 

# To improve computation time, we use Rcpp and the likelihood/posterior coded in C
sourceCpp('lh_foropt_post.cpp')
sourceCpp('lh_foropt.cpp')

print("Read in the data...")

set.seed(12345)

# read in data
data<-read.table("../data/FMD/data_rd.txt") # removal times
n_c<-read.table("../data/FMD/n_c.txt") #number of cows per farm
n_s<-read.table("../data/FMD/n_s.txt") #number of sheep per farm
xc<-read.table("../data/FMD/x_coord_rd.txt") # x coordinate of farm
yc<-read.table("../data/FMD/y_coord_rd.txt")  # y coordinate of farm
N<-5378 # population size
n<-1021 # final outbreak size
m<-4 # shape parameter of gamma distributed infectious periods

print("Prepare to run...")

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

# set up and order cow/sheep vectors correctly
n_c<-n_c[,1]
n_c<-c(n_c[4358:5378], n_c[1:4357])
n_s<-n_s[,1]
n_s<-c(n_s[4358:5378], n_s[1:4357])
xc<-xc[,1]
xc<-c(xc[4358:5378], xc[1:4357])
yc<-yc[,1]
yc<-c(yc[4358:5378], yc[1:4357])

total_animals<-n_c+n_s

print("Set up distance matrix")
# euclidean distance matrix
pmat<-matrix(NA,N,N)
for (i in 1:N){
  for (j in 1:N){
    pmat[i,j]<-(xc[i]-xc[j])^2 + (yc[i]-yc[j])^2
  }
}

r<-data[,2]


###############################
## NLM optimisation
print("Ready to optimise")
# note: warning 'NA/Inf replaced by maximum positive value' from NLM may occur - unacceptable parameter values just being tested

# first transform to negative log posterior
nlog_post<-function(x,R2=r,pmat=pmat, N=N,n=n, m=4, n_c=n_c, n_s=n_s){
  return(-log_post(x, R2=r, pmat=pmat, N=N,n=n, m=4, n_c=n_c, n_s=n_s))
}
# then minimise
nlm(f=nlog_post, p=c(0.00001, 0.8, 0.0001, 2.0, 0.5, 0.8),  print.level = 2, R2=r, pmat=pmat, N=N,n=n, m=4, n_c=n_c, n_s=n_s)
# this provides the MAPs

###############################
## Profile log likelihoods

 beta0<-0
 params=c(0.000000703,0.445,0.00477,1.576,2.389,0.318)
 for (i in 0:1000) {
    params[1]<-i/100000000
     print(i)
     beta0[i+1] = log_lh(params=params,  R2=r, pmat=pmat,  N=N,  n=n,  m=4, n_c=n_c, n_s=n_s)
 }
 print("Completed beta0")
 write.table(beta0, "beta0prof.txt")


 gamma<-0
 params=c(0.000000703,0.445,0.00477,1.576,2.389,0.318)
 for (i in 1:100) {
     params[2]<-i/100
     print(i)
     gamma[i] = log_lh(params=params,  R2=r,  pmat=pmat,  N=N,  n=n,  m=4, n_c=n_c, n_s=n_s)
 }
 print("Completed gamma")
 write.table(gamma, "gammaprof.txt")

v<-0
params=c(0.000000703,0.445,0.00477,1.576,2.389,0.318)
for (i in 1:100) {
  params[3]<-i/10000
  print(i)
  v[i+1] = log_lh(params=params,  R2=r, pmat=pmat,  N=N,  n=n,  m=4, n_c=n_c, n_s=n_s)
}
print("Completed v")
write.table(v, "vprof.txt")

ep<-0
params=c(0.000000703,0.445,0.00477,1.576,2.389,0.318)
for (i in 1:100) {
  params[4]<-i/10
  print(i)
  ep[i] = log_lh(params=params,  R2=r,  pmat=pmat,  N=N,  n=n,  m=4, n_c=n_c, n_s=n_s)
}
print("Completed ep")
write.table(ep, "epprof.txt")
gc()

xi<-0
params=c(0.000000703,0.445,0.00477,1.576,2.389,0.318)
for (i in 0:100) {
  params[5]<-i/10
  print(i)
  xi[i+1] = log_lh(params=params,  R2=r, pmat=pmat,  N=N,  n=n,  m=4, n_c=n_c, n_s=n_s)
}
print("Completed xi")
write.table(xi, "xiprof.txt")


ze<-0
params=c(0.000000703,0.445,0.00477,1.576,2.389,0.318)
for (i in 1:100) {
  params[6]<-i/100
  print(i)
  ze[i] = log_lh(params=params,  R2=r,  pmat=pmat,  N=N,  n=n,  m=4, n_c=n_c, n_s=n_s)
}
print("Completed ze")
write.table(ze, "zeprof.txt")
gc()




