library(deSolve)
## HETEROGENEOUS MODEL 
HSL1L2IT_model=function(current_timepoint, state_values, params)
{
  #creat state variables (local variables)
  S=state_values[1] # susceptible 
  L1=state_values[2] #early latency
  L2=state_values[3] #late latency
  I=state_values[4] #infectious
  T=state_values[5] #Under treatment
  
  with(
    as.list(params), 
    {
      #compute derivative
      dS=(mu*N+phi*T+mui*I+mut*T)-(k*log(1+beta*I/k/N))*S-mu*S
      dL1=(k*log(1+beta*I/k/N))*S+(r*k*log(1+beta*I/k/N)*L2)-(epsilon+mu+kappa)*L1
      dL2=(kappa*L1+gamma*I)-(r*k*log(1+beta*I/k/N))*L2-(nu+mu)*L2
      dI=(epsilon*L1+nu*L2+omega*T)-(mu+mui+gamma+delta)*I
      dT=(delta*I)-(mut+mu+omega+phi)*T
      dinc=(epsilon*L1+nu*L2)
      dnotif=delta*I
      #combine results
      results=c(dS, dL1, dL2, dI, dT,dinc,dnotif)
      list(results)
    }
  )
}

##Initial values for sub population: 
A=1 #Fully susceptible hosts
B=0        #Early latent hosts
C=0        #Late latent hosts
D=1e-3        #Infectious hosts
E=0        #Treated hosts




#=====================================================================
#SENSITIVITY ANALYSIS
#Latin hypercube sampling
#install.packages('lhs')
require(lhs) #add the lhs library
h <- 1000 #choose number of points to simulate
set.seed(6242015)#random number generator
lhs<-maximinLHS(h,13) #simulate h= number of simulations, 13=number of parameters
#To map these points in the unit cube to our parameters, we need minimum and maximum values for each.
mu.min=0.015
mu.max=0.015
mut.min=0.052
mut.max=0.052
mui.min=0.12
mui.max=0.12
r.min=0.21
r.max=0.21
k.min<-0
k.max<-1
beta.min<-4
beta.max<-24 
epsilon.min<-0.01
epsilon.max<-0.5
kappa.min<-0.1
kappa.max<-0.9
gamma.min<-0.03
gamma.max<-0.5
nu.min<-0
nu.max<-0.1
omega.min<-0
omega.max<- 1
phi.min<-0.1
phi.max<- 2
delta.min<-0.5
delta.max<-6
#Now we can generate a “parameter set” by rescaling our simulated latin hypercube sample
params.set <- cbind(
  k = lhs[,1]*(k.max-k.min)+k.min,
  beta = lhs[,2]*(beta.max-beta.min)+beta.min,
  epsilon =lhs[,3]*(epsilon.max-epsilon.min)+epsilon.min,
  kappa = lhs[,4]*(kappa.max-kappa.min)+kappa.min,
  gamma = lhs[,5]*(gamma.max-gamma.min)+gamma.min,
  nu = lhs[,6]*(nu.max-nu.min)+nu.min,
  omega = lhs[,7]*(omega.max-omega.min)+omega.min,
  phi = lhs[,8]*(phi.max-phi.min)+phi.min,
  delta = lhs[,9]*(delta.max-delta.min)+delta.min,
  mu = lhs[,10]*(mu.max-mu.min)+mu.min,
  mui = lhs[,11]*(mui.max-mui.min)+mui.min,
  mut = lhs[,12]*(mut.max-mut.min)+mut.min,
  r = lhs[,13]*(r.max-r.min)+r.min)
View(params.set)
#creat matrix to save whole info
output_matrix = data.frame(params.set)
incidence = data.frame('incidence'=rep(NA,h))
output_matrix = cbind(output_matrix, incidence) #add incidence column

#(#THese are all the simulated parameters we need to exlore how these parameters affect the incidence


#we set up a nested loop, first to cycle through different values of INCIDENCE, then to cycle through
#different simulated parameter sets. Note the pre-allocated data frame and use of the counter j, also to
#speed up evaluation



for(i in 1:h){
  
  #initial_values=c(S=0.5418456958, L1=0.0197116035, L2=0.4348963153, I=0.0026741092, T=0.0008722763,inc=0.0026741092,notif=0.0008722763)#equilibrium
  initial_values=c(S=A-D, L1=B, L2=C, I=D, T=E,inc=D,notif=E)
  
  params <- as.list(c(params.set[i,]))
  times=seq(0, 1000, by = 1)
  out <- as.data.frame(lsoda(initial_values, times, HSL1L2IT_model, params))
  
  Incidence=diff(out$inc)*100000
  output_matrix$incidence[i] = Incidence
  
}
View(output_matrix) #now we have incidence and each parameters in one matrix




#Figure 1 of Appendix 6
#We can plot each simulated value as a point
par(mfrow=c(3,3))
plot(output_matrix$epsilon, output_matrix$incidence, type = 'p' ,lwd=2,pch=19, cex=0.9, col='blue',
     xlab=expression(epsilon),
     ylab='Incidence',main = 'PRCC= 0.91')

plot(output_matrix$beta, output_matrix$incidence, type = 'p' ,lwd=2,pch=19, cex=0.9, col='blue',
     xlab=expression(beta),
     ylab='Incidence',main = 'PRCC= 0.85')
plot(output_matrix$delta, output_matrix$incidence, type = 'p' ,lwd=2,pch=19, cex=0.9, col='blue',
     xlab=expression(delta),
     ylab='Incidence',main = 'PRCC= -0.66')
plot(output_matrix$omega, output_matrix$incidence, type = 'p' ,lwd=2,pch=19, cex=0.9, col='blue',
     xlab=expression(omega),
     ylab='Incidence',main = 'PRCC= 0.14')
plot(output_matrix$kappa, output_matrix$incidence, type = 'p' ,lwd=2,pch=19, cex=0.9, col='blue',
     xlab=expression(kappa),
     ylab='Incidence',main = 'PRCC= -0.14')



plot(output_matrix$k, output_matrix$incidence, type = 'p' ,lwd=2,pch=19, cex=0.9, col='blue',
     xlab='dispersion parameter',
     ylab='Incidence',main = 'PRCC= 0.12')

plot(output_matrix$nu, output_matrix$incidence, type = 'p' ,lwd=2,pch=19, cex=0.9, col='blue',
     xlab=expression(nu),
     ylab='Incidence',main = 'PRCC= 0.10')

plot(output_matrix$gamma, output_matrix$incidence, type = 'p' ,lwd=2,pch=19, cex=0.9, col='coral4',
     xlab=expression(gamma),
     ylab='Incidence',main = 'PRCC= -0.06')

plot(output_matrix$phi, output_matrix$incidence, type = 'p' ,lwd=2,pch=19, cex=0.9, col='coral4',
     xlab=expression(phi),
     ylab='Incidence',main = 'PRCC= -0.008')




#which parameter affects more,inspecting the partial correlation
#Partial rank correlations can be computed using the pcc function in the R package sensitivity 
#install.packages('sensitivity')
library(sensitivity)
bonferroni.alpha <- 0.05/9
prcc <- pcc(output_matrix[,1:9], output_matrix[,14], nboot = 1000, rank=TRUE, conf=1-bonferroni.alpha)
save(prcc, file='prcc.Rdata')

#We can view a table of the resulting partial correlation coefficients. if none of the (penalized)
#confidence intervals contains zero, we conclude that all are significant and produce a plot showing their
#relative magnitudes.
library(sensitivity)
load('prcc.Rdata')
summary <- print(prcc)
Corl=data.frame(summary)
View(Corl)
#plote the partial corelation coeeficient

par(mar=c(9,4,4,2)+0.1)
plot(Corl$original, main='Partial rank correlation coefficients', ylim=c(-1,1),
     xlab='', ylab='Coefficient',
     axes=FALSE)
axis(2)
axis(1, at=seq(1:9), labels=row.names(Corl), las=2)
mtext(text='Parameter', side=1, line=4.5)
box()
for(i in 1:9) lines(c(i,i),c(Corl[i,4], Corl[i,5]))
abline(h=0)


#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++