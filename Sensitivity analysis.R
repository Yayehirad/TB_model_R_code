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
# Set parameters 
risk_after_reinfection=0.79
life_expectancy_years=65
proportion_of_progress_to_activeTB_from_early_latency_over_2yeras=0.075
proportion_of_progress_to_late_latency_from_early_latency_over_2yeras=0.925
proportion_of_progress_to_activeTB_from_late_latency_over_20yeras=0.07
TB_fatality_during_active_disea_over_3years=0.37
proporton_spontaneous_recovery_from_active_disease_over_3years=0.63
proportion_of_unsuccessful_treatment_over_0.5year=0.231
proportion_TBdeath_during_teatment_over_0.5year=0.026
proportion_successful_treatment_over_0.5year=0.728
case_detection_over_1year=0.95
dispersion_parameter=0.1
transmission_events_per_infectious_person_year=6
total_population=1
r=1-risk_after_reinfection 
N=total_population
mu=1/life_expectancy_years #Non-TB specific death rate
beta=transmission_events_per_infectious_person_year #Effective contact rate 

phi=proportion_successful_treatment_over_0.5year*2 
kappa=proportion_of_progress_to_late_latency_from_early_latency_over_2yeras/2 
epsilon=proportion_of_progress_to_activeTB_from_early_latency_over_2yeras/2 

gamma=proporton_spontaneous_recovery_from_active_disease_over_3years/3 
nu=proportion_of_progress_to_activeTB_from_late_latency_over_20yeras/20 



omega=proportion_of_unsuccessful_treatment_over_0.5year*2
mui=TB_fatality_during_active_disea_over_3years/3 
mut=proportion_TBdeath_during_teatment_over_0.5year*2
delta=case_detection_over_1year*(gamma+mui+mu)/(1-case_detection_over_1year) 
k=dispersion_parameter
#force of infection
##Initial values for sub population: 
A=1 #Fully susceptible hosts
B=0        #Early latent hosts
C=0        #Late latent hosts
D=1e-3        #Infectious hosts
E=0        #Treated hosts


params=c(N=N,mu=mu, k=k, beta=beta, epsilon=epsilon, kappa=kappa,gamma=gamma,nu=nu,omega=omega, 
                 mut=mut,mui=mui,phi=phi, delta=delta)


#Initial state values for differential equasions
initial_values=c(S=A-D, L1=B, L2=C, I=D, T=E,inc=D,notif=E)#full dynamics
#equilibrium
#initial_values=c(S=0.5418456958, L1=0.0197116035, L2=0.4348963153, I=0.0026741092, T=0.0008722763,inc=0.0026741092,notif=0.0008722763)#equilbrium values
## Output Time points
times=seq(0, 2000, by = 1)
#Simulate the HS1S2L1L2IT transmission 
Houtput=as.data.frame(lsoda(initial_values, times, HSL1L2IT_model, params))
#View(Houtput)
### fill vectors with results
HvS=Houtput$S
HvL1=Houtput$L1
HvL2=Houtput$L2
HvI=Houtput$I
HvT=Houtput$T
Hvinc=Houtput$inc
Hvnotif=Houtput$notif
Hvtime=Houtput$time
HvN=HvS+HvL1+HvL2+HvI+HvT
incidence.baseline=diff(Hvinc)
plot(incidence.baseline)
max(incidence.baseline)
#################################
Incidence=diff(Houtput$inc) #incidence

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
#To facilitate later plotting, we decide how many different levels of I to consider.
#levels <- 15
#Even though we have simulated 100 points, to speed up computing we will use only a fraction of these
#for this demonstration.
#h2 <-30

#we set up a nested loop, first to cycle through different values of I, then to cycle through
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
write.csv(x=output_matrix,'..//Documents/sensitivity.csv')


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

write.csv(x=Corl,file='..//Documents/corl.csv')
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
#tornado plot


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


 
 


#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++