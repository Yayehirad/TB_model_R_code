#HOMOGENEOUS ASSUMPTION (MODEL 1)
library(deSolve)

SL1L2IT_model=function(current_timepoint, state_values, parameters)
{
  #creat state variables (local variables)
  S=state_values[1] # susceptible 
  L1=state_values[2] #early latenyy
  L2=state_values[3] #late latency
  I=state_values[4] #infectious
  T=state_values[5] #Under treatment
  
  with(
    as.list(parameters), #variable names within parameters can be used
    {
      #compute derivative
      
      dS=mu*N+mui*I+mut*T+phi*T - (beta*I/N)*S-mu*S
      dL1=(beta*I/N)*S+(r*beta*I/N)*L2  - (epsilon+mu+kappa)*L1
      dL2=kappa*L1+gamma*I -  nu*L2-(r*beta*I/N)*L2-mu*L2
      dI=epsilon*L1+nu*L2+omega*T - (mu+mui+gamma+delta)*I
      dT=delta*I - (mu+mut+omega+phi)*T
      
      dinc=epsilon*L1+nu*L2
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
proportion_of_unsuccessful_treatment_over_0.5year=0.145
proportion_TBdeath_during_teatment_over_0.5year=0.025
proportion_successful_treatment_over_0.5year=0.83
case_detection_over_1year=0.65

transmission_events_per_infectious_person_year=7
total_population=1
r=1-risk_after_reinfection #risk of progression after re-infection
N=total_population
mu=1/life_expectancy_years #Non-TB specific death rate
beta=transmission_events_per_infectious_person_year #Effective contact rate 

phi=proportion_successful_treatment_over_0.5year*2 # Cure rate after treatment(0.769 over 1/2 yrs)

kappa=proportion_of_progress_to_late_latency_from_early_latency_over_2yeras/2 #Rate of progression to late latency after infection(0.925 over 2yrs)
epsilon=proportion_of_progress_to_activeTB_from_early_latency_over_2yeras/2 #Rate progression to active TB from early latency(0.075)

gamma=proporton_spontaneous_recovery_from_active_disease_over_3years/3 #Rate of spontaneous recovery to latency(0.63 over 3yrs)   
nu=proportion_of_progress_to_activeTB_from_late_latency_over_20yeras/20 #Rate of progression to active TB from late LTBI(0.1 over 20yrs)



omega=proportion_of_unsuccessful_treatment_over_0.5year*2 #Unsuccessful treatment or non-cure rate (0.231 over 1/2 yrs)
mui=TB_fatality_during_active_disea_over_3years/3 #TB specific death during Infectious period
mut=proportion_TBdeath_during_teatment_over_0.5year*2# TB related death during Treatment
delta=case_detection_over_1year*(gamma+mui+mu)/(1-case_detection_over_1year) #Case detection rate()
#Disease dynamic parameters
parameter_list=c(N=N,mu=mu,mut=mut,mui=mui, r=r,beta=beta, epsilon=epsilon, kappa=kappa,gamma=gamma,nu=nu,
                 omega=omega, phi=phi, delta=delta)
##Initial values for sub population: 

A=1 #Fully susceptible hosts
B=0        #Early latent hosts
C=0        #Late latent hosts
D=1e-6        #Infectious hosts
E=0        #Treated hosts
#Compute total population

#Initial state values for differential equasions
initial_values=c(S=A-D, L1=B, L2=C, I=D, T=E,inc=D,notif=E)

times=seq(0, 2000, by = 1)
#Simulate the SL1L2IT transmission 
output=as.data.frame(lsoda(initial_values, times, SL1L2IT_model, parameter_list))
#View(output)
### fill vectors with results
vS=output$S
vL1=output$L1
vL2=output$L2
vI=output$I
vT=output$T
vinc=output$inc
vnotif=output$notif
vtime=output$time
vN=vS+vL1+vL2+vI+vT
plot(vtime,vN,ylim=c(0,1))

#prevalence 
Prev_homo=((vI+vT)/vN)*100000
#incidence 
Incidence_homo=(diff(vinc)/vN)*100000
max(Incidence_homo)
#notification rate per 100,000
Notif_homo=(diff(vnotif)/vN)*100000

#====================================================================================
#Hetrogeneous assumption (MODEL TWO)
#
library(deSolve)
## creat SLALBIT function ## let LA=L1, LB=L2
HSL1L2IT_model=function(current_timepoint, state_values, parameters)
{
  #creat state variables (local variables)
  S=state_values[1] #fully susceptible 
  L1=state_values[2] #early latency
  L2=state_values[3] #late latency
  I=state_values[4] #infectious
  T=state_values[5] #Under treatment
  
  with(
    as.list(parameters), #variable names within parameters can be used
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
proportion_of_unsuccessful_treatment_over_0.5year=0.145
proportion_TBdeath_during_teatment_over_0.5year=0.025
proportion_successful_treatment_over_0.5year=0.83
case_detection_over_1year=0.65


total_population=1
r=1-risk_after_reinfection #risk of progression after re-infection
N=total_population
mu=1/life_expectancy_years #Non-TB specific death rate
 

phi=proportion_successful_treatment_over_0.5year*2 # Cure rate after treatment(0.769 over 1/2 yrs)
kappa=proportion_of_progress_to_late_latency_from_early_latency_over_2yeras/2 #Rate of progression to late latency after infection(0.925 over 2yrs)
epsilon=proportion_of_progress_to_activeTB_from_early_latency_over_2yeras/2 #Rate progression to active TB from early latency(0.075)

gamma=proporton_spontaneous_recovery_from_active_disease_over_3years/3 #Rate of spontaneous recovery to latency(0.63 over 3yrs)   
nu=proportion_of_progress_to_activeTB_from_late_latency_over_20yeras/20 #Rate of progression to active TB from late LTBI(0.1 over 20yrs)



omega=proportion_of_unsuccessful_treatment_over_0.5year*2 #Unsuccessful treatment or non-cure rate (0.231 over 1/2 yrs)
mui=TB_fatality_during_active_disea_over_3years/3 #TB specific death during Infectious period
mut=proportion_TBdeath_during_teatment_over_0.5year*2# TB related death during Treatment
delta=case_detection_over_1year*(gamma+mui+mu)/(1-case_detection_over_1year) #Case detection rate()

dispersion_parameter=0.05
transmission_events_per_infectious_person_year=7
beta=transmission_events_per_infectious_person_year #Effective contact rate
k=dispersion_parameter
##Initial values for sub population: 
A=1 #Fully susceptible hosts
B=0        #Early latent hosts
C=0        #Late latent hosts
D=1e-6        #Infectious hosts
E=0        #Treated hosts
parameter_list=c(N=N,mu=mu, k=k, beta=beta, epsilon=epsilon, kappa=kappa,gamma=gamma,nu=nu,omega=omega, 
                 mut=mut,mui=mui,phi=phi, delta=delta)
#Initial state values for differential equasions
initial_values=c(S=A-D, L1=B, L2=C, I=D, T=E,inc=D,notif=E)

## Output Time points
times=seq(0, 2000, by = 1)
#Simulate the S1S2L1L2IT transmission 
Houtput=as.data.frame(lsoda(initial_values, times, HSL1L2IT_model, parameter_list))
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
#incidence  
Incidence_hetero=(diff(Hvinc)/HvN)*100000
max(Incidence_hetero)

HIncidence0.2=(diff(Hvinc)/HvN)*100000

max(HIncidence0.1)


HIncidence1=(diff(Hvinc)/HvN)*100000
HIncidence0.5=(diff(Hvinc)/HvN)*100000
HIncidence0.1=(diff(Hvinc)/HvN)*100000
HIncidence0.05=(diff(Hvinc)/HvN)*100000






#prevalence 
Prev_hetero=((HvI+HvT)/HvN)*100000


#notification rate per 100,000
Notif=(diff(Hvnotif)/HvN)*100000




#=====================================================================

++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
##### Force of infection

k=0.05
I=c(0,0.001,0.002,0.003,0.004,0.005,0.006,0.007,0.008,0.009,0.01,0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.1,0.2)
b=7
N=1
lamb_0.05=k*log(1+b*I/k/N)
lambd0=b*I/N
print(lamb)
print(lambd0)
lamb_1=as.data.frame(lamb)

View(lambd)
View(lambd0)
plot(lambd)
write.csv(x=lamb_0.05,file='..//Documents/lamb0.05.csv')
write.csv(x=lambd0,file='..//Documents/lambds00.csv')
lambdas1=read.csv(file.choose())
View(lambdas1)
plot(lambdas1$I,lambdas1$lambda_o,ylim=c(0,0.05),xlim =c(0,0.006))
plot(lambdas1$I,lambdas1$lambda_o,type = 'n',ylim=c(0,0.05),xlim =c(0,0.006),ylab = "Force of infection per year",
     xlab="TB prevalence")
lines(lambdas1$I,lambdas1$lambda_o,type = 'l',col='black',lwd=3)

lines(lambdas1$I,lambdas1$lambda_k_1, type = 'l',col='blue',lwd=3)
  lines(lambdas1$I,lambdas1$lambda_k_0.5, type = 'l',col='purple',lwd=3)
lines(lambdas1$I,lambdas1$lambda_k_0.1, type = 'l',col='magenta',lwd=3)
lines(lambdas1$I,lambdas1$lambda_k_0.05, type = 'l',col='red',lwd=3)
#lines(lambdas1$I,lambdas1$lambda_k_0.01, type = 'l',col='darkorchid',lwd=3)

legend("topleft", c("homogeneous","k=1","k=0.5","k=0.1","k=0.05"),lty=c(1,1,1,1,1,1),lwd=c(3,3,3,3,3,3),
       col=c('black','blue','purple','magenta','red'))
#FIGURE 3
library(reshape2)
library(ggplot2)
View(lambdas1)
lambdam=melt(lambdas1,id='I')
View(lambdam)
ggplot()+
  geom_line(data = lambdam,aes(x=I,y=value,colour=variable),size=1.05)+
  ylab("Force of infection per year")+
  xlab("TB prevalence")+
  theme(legend.position = c(0.8,0.2))+

  scale_colour_manual(values=c("black","green","purple","magenta","red","pink"),name="dispersion",
                      labels=c("Homogeneous","k=10","k=1","k=0.5","k=0.1","k=0.05"))
  
#them(legend.text=element_text("Homogeneous","k=10","k=1","k=0.5","k=0.1","k=0.05"))

# 
  #scale_fill_discrete(
                     # breaks=c("lambda_o","lambda_k_10","lambda_k_1","lambda_k_0.5","lambda_k_0.1","lambda_k_0.05"),
                      #)
ggplot()+
  geom_line(data = lambdas1,aes(x=I,y=lambda_o),size=1.5,colour='black')+
  geom_line(data = lambdas1,aes(x=I,y=lambda_k_10),size=1.05,colour='green')+
  geom_line(data = lambdas1,aes(x=I,y=lambda_k_1),size=1.05,colour='purple')+
  geom_line(data = lambdas1,aes(x=I,y=lambda_k_0.5),size=1.05,colour='magenta')+
  geom_line(data = lambdas1,aes(x=I,y=lambda_k_0.1),size=1.05,colour='red')+
  geom_line(data = lambdas1,aes(x=I,y=lambda_k_0.05),size=1.05,colour='pink')+
  ylab("Force of infection per year")+
  xlab("TB prevalence")+
  theme(legend.position = c(0.8,0.2))+
  
  scale_colour_manual(values=c("black","green","purple","magenta","red","pink"),name="dispersion",
                      labels=c("Homogeneous","k=10","k=1","k=0.5","k=0.1","k=0.05"))
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
kbetai=read.csv(file.choose())
View(kbetai)
plot(kbetai$K,kbetai$X24,log="x",type = 'n',ylim=c(0,400),xlim=c(0.01,1),ylab = "Incidence per 100,000 popn.",
     xlab="Dispersion parameter, k")

#lines(kbetai$K,kbetai$X24,type = 'l',col='blue',lwd=3,lty=1)
#lines(kbetai$K,kbetai$X20,type = 'l',col='green',lwd=3,lty=1)
#lines(kbetai$K,kbetai$X15,type = 'l',col='goldenrod',lwd=3,lty=1)
lines(kbetai$K,kbetai$X10,type = 'l',col='black',lwd=3,lty=1)
lines(kbetai$K,kbetai$X9,type = 'l',col='blue',lwd=3,lty=1)
lines(kbetai$K,kbetai$X8,type = 'l',col='purple',lwd=3,lty=1)
lines(kbetai$K,kbetai$X7,type = 'l',col='violet',lwd=3,lty=1)

lines(kbetai$K,kbetai$X6,type = 'l',col='pink',lwd=3,lty=1)
lines(kbetai$K,kbetai$X5,type = 'l',col='orange',lwd=3,lty=1)
lines(kbetai$K,kbetai$X4,type = 'l',col='palegoldenrod',lwd=3,lty=1)
#lines(kbetap$K,kbetap$X3,type = 'l',col='gray8',lwd=3)


#legend("topleft",legend=c("24","20","15","10","5","4","3"), title = "beta",bg='lightblue',cex=0.8,
# bty = "n",lty=c(1,1,1,1,1,1,1),lwd=c(3,3,3,3,3,3,3),col=c("blue",'green','goldenrod','chocolate3','darkorchid','firebrick1','gray8'))



legend("topleft",legend=c("24","20","15","10","8","7","6","5","4"), title = "beta",bg='lightblue',cex=0.8,
       bty = "n",lty=c(1,1,1,1,1,1,1,1,1),lwd=c(3,3,3,3,3,3,3,3,3),col=c("blue",'green','goldenrod','chocolate3','chocolate4',
                                                                         'cornflowerblue','darkorchid4','darkorchid','firebrick1'))
legend("topleft",legend=c("10","9","8","7","6","5","4"), title = "ß",bg='lightblue',cex=0.8,
       bty = "n",lty=c(1,1,1,1,1,1),lwd=c(3,3,3,3,3,3),col=c('black','blue','purple',
                                                                         'violet','pink','orange','palegoldenrod'))

#legend("topright", inset=c(-0.2,0),xpd = T,legend=c("24","20","15","10","5","4","3"), title = "beta", 
#bty = "n",lty=c(1,1,1,1,1,1,1),lwd=c(3,3,3,3,3,3,3),col=c("blue",'green','goldenrod','chocolate3','darkorchid','firebrick1','gray8'))
#FIGURE 5
library(reshape2)
library(ggplot2)
kbetagg=read.csv(file.choose())
kbetaM=melt(kbetagg,id='K')
View(kbetaM)
ggplot(data = kbetaM,aes(x=K,y=value,colour=variable))+
  geom_line(size=1.05)+
  ylim(0,350)+
  xlim(0.01,1)+
  ylab("Incidence per 100,000 population")+
  xlab("Dispersion parameter, k")+
  scale_x_log10()+
  scale_colour_manual(values=c("gold1","goldenrod","tan3","chocolate3","chocolate4","coral4","darkorchid4","firebrick1"),
                      name="ß",
                    breaks=c("X10", "X9", "X8","X7","X6","X5","X4"),
                    labels=c("10", "9", "8","7","6","5","4"))+
  theme(legend.justification = "top")
  #theme(legend.title = element_text(colour="black", size=10, face="bold"))+
 # theme(legend.text = element_text(colour="black", size=10, face="bold"))+


#+
  
  
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++






#================================================================================================================
#MODEL THREE
library(deSolve)
#SUBCLASS FOR INFECTIOUS GROUP
## creat SLALBIT function ## let LA=L1, LB=L2
subclass_model=function(current_timepoint, state_values, parameters)
{
  #creat state variables (local variables)
  S=state_values[1] #fully susceptible 
  L1=state_values[2] #early latency
  L2=state_values[3] #late latency
  I0=state_values[4] #infectious non-spreaders
  I1=state_values[5] #infectious spreaders
  I2=state_values[6] #infectious Super-spreaders
  T=state_values[7] #Under treatment
  
  with(
    as.list(parameters), #variable names within parameters can be used
    {
      #compute derivative
      dS=(mu*N+phi*T+mui*I0+mui*I1+mui*I2+mut*T)-((beta0*I0+beta1*I1+beta2*I2)/N)*S-mu*S
      dL1=((beta0*I0+beta1*I1+beta2*I2)/N)*S+(r*(beta0*I0+beta1*I1+beta2*I2)/N)*L2-(epsilon0+epsilon1+epsilon2+mu+kappa)*L1
      dL2=kappa*L1+(gamma0*I0+gamma1*I1+gamma2*I2)-(r*(beta0*I0+beta1*I1+beta2*I2)/N)*L2-(nu0+nu1+nu2+mu)*L2
      dI0=(epsilon0*L1+nu0*L2+omega0*T)-(mui+mu+gamma0+delta0)*I0
      dI1=(epsilon1*L1+nu1*L2+omega1*T)-(mui+mu+gamma1+delta1)*I1
      dI2=(epsilon2*L1+nu2*L2+omega2*T)-(mui+mu+gamma2+delta2)*I2
      dT=delta0*I0+delta1*I1+delta2*I2-(mut+mu+omega0+omega1+omega2+phi)*T
      
      dinc=(epsilon0*L1+nu0*L2)+(epsilon1*L1+nu1*L2)+(epsilon2*L1+nu2*L2)
      dnotif=delta0*I0+delta1*I1+delta2*I2
      #combine results
      results=c(dS, dL1, dL2, dI0,dI1,dI2, dT,dinc,dnotif)
      list(results)
    }
  )
}
# Set parameters 
a=0.642 #proportion of Io
b=0.226#proportion of I1
c=0.132##proportion of I2
risk_after_reinfection=0.79
life_expectancy_years=65
proportion_of_progress_to_activeTB_from_early_latency_over_2yeras=0.075
proportion_of_progress_to_late_latency_from_early_latency_over_2yeras=0.925
proportion_of_progress_to_activeTB_from_late_latency_over_20yeras=0.07
TB_fatality_during_active_disea_over_3years=0.37
proporton_spontaneous_recovery_from_active_disease_over_3years=0.63
proportion_of_unsuccessful_treatment_over_0.5year=0.145
proportion_TBdeath_during_teatment_over_0.5year=0.025
proportion_successful_treatment_over_0.5year=0.83

case_detection_non_infectious_over_1year=0.65
case_detection_infectious_nonsuper_spreaders_over_1year=0.65
case_detection_super_spreaders_over_1year=0.65


transmission_events_per_non_infectious_person_year=0
transmission_events_per_infectious_non_super_spreader_person_year=4
transmission_events_per_infectious_super_spreader_person_year=46
total_population=1

r=1-risk_reduction_after_reinfection #risk of progression after re-infection
N=total_population
mu=1/life_expectancy_years #Non-TB specific death rate

beta0=transmission_events_per_non_infectious_person_year 
beta1=transmission_events_per_infectious_non_super_spreader_person_year
beta2=transmission_events_per_infectious_super_spreader_person_year

phi=proportion_successful_treatment_over_0.5year*2 # Cure rate after treatment(0.769 over 1/2 yrs)
kappa=proportion_of_progress_to_late_latency_from_early_latency_over_2yeras/2 #Rate of progression to late latency after infection(0.925 over 2yrs)
epsilon=proportion_of_progress_to_activeTB_from_early_latency_over_2yeras/2 #Rate progression to active TB from early latency(0.075)

gamma=proporton_spontaneous_recovery_from_active_disease_over_3years/3 #Rate of spontaneous recovery to latency(0.63 over 3yrs)   
nu=proportion_of_progress_to_activeTB_from_late_latency_over_20yeras/20 #Rate of progression to active TB from late LTBI(0.1 over 20yrs)



omega=proportion_of_unsuccessful_treatment_over_0.5year*2 #Unsuccessful treatment or non-cure rate (0.231 over 1/2 yrs)
mui=TB_fatality_during_active_disea_over_3years/3 #TB specific death during Infectious period
mut=proportion_TBdeath_during_teatment_over_0.5year*2# TB related death during Treatment


gamma0=gamma
gamma1=gamma
gamma2=gamma

nu0=nu*a
nu1=nu*b
nu2=nu*c

epsilon0= epsilon*a
epsilon1= epsilon*b
epsilon2= epsilon*c


omega0=omega*a
omega1=omega*b
omega2=omega*c

cdr0=case_detection_non_infectious_over_1year
cdr1=case_detection_infectious_nonsuper_spreaders_over_1year
cdr2=case_detection_super_spreaders_over_1year

delta0=cdr0*(gamma+mui+mu)/(1-cdr0)
delta1=cdr1*(gamma+mui+mu)/(1-cdr1)
delta2=cdr2*(gamma+mui+mu)/(1-cdr2)
##Initial values for sub population: 
A=1 #Fully susceptible hosts
B=0        #Early latent hosts
C=0        #Late latent hosts
D=1e-6*a     #Infectious hosts proportion of extrapulmonary 15% of all infectious 
E=1e-6*b # Normal spreaders
F=1e-6*c # Super-spreader 10% of TB patients are super-spreaders
G=0 #Treated hosts
H=0 #incidence
I=0 # diagnosed 

parameter_list=c(N=N,mu=mu,beta0=beta0, beta1=beta1, epsilon0=epsilon0,epsilon1=epsilon1,epsilon2=epsilon2, 
                 kappa=kappa,gamma0=gamma0,gamma1=gamma1,gamma2=gamma2,nu0=nu0,nu1=nu1,nu2=nu2,
                 omega0=omega0, omega1=omega1,omega2=omega2,phi=phi, 
                 delta0=delta0,delta1=delta1,delta2=delta2)

#Initial state values for differential equasions
initial_values=c(S=A-(D+E+F), L1=B, L2=C, I0=D,I1=E,I2=F,T=G,inc=D+E+F,notif=I)
#at Equilibrium
#initial_values=c(S=0.4915928963, L1=0.0227268406, L2=0.4819654218, I0= 0.0017543290,I1=0.0006859876,I2=0.0003711080,T=0.0009034167,inc= 0.0017543290 +0.0006859876+ 0.0003711080,notif=0.0009034167)
## Output Time points
times=seq(0, 2000, by = 1)
#Simulate the S1S2L1L2IT transmission 
soutput=as.data.frame(lsoda(initial_values, times, subclass_model, parameter_list))
#View(soutput)
### fill vectors with results
svS=soutput$S
svL1=soutput$L1
svL2=soutput$L2
svI0=soutput$I0
svI1=soutput$I1
svI2=soutput$I2
svT=soutput$T
svinc=soutput$inc
svnotif=soutput$notif
svtime=soutput$time
svN=svS+svL1+svL2+svI0+svI1+svI2+svT
plot(svtime,svN,ylim = c(0,1))


##########

##########
#prevalence 
svPrev=((svI0+svI1+svI2+svT)/svN)*100000
svPrev
max(svPrev)
#equi
esvPrev=((svI0+svI1+svI2+svT)/svN)*100000
esvPrev
max(esvPrev)
#plote prevalence
plot(svtime,esvPrev,type='l', col = 'blue',lwd=3,main="Prevalence",
     ylab="fraction of infectious and treated",xlab="Time in years",xlim=c(0,1000),ylim = c(0,350)) 

#incidence the fraction of new cases at each time step 
sIncidence=(diff(svinc)/svN)*100000
max(sIncidence)
#equi
esIncidence=diff(svinc)*100000
max(esIncidence)
#plot incidence
plot (svtime,sIncidence, type='l', col = 'red',lwd=3,main="Incidence",
      ylab="incidence per 100,000 population",xlab="Time in years",xlim = c(0,300),ylim = c(0,300))

#notification rate per 100,000
Notif=(diff(svnotif)/svN)*100000
max(Notif)
plot (svtime, Notif, type='l', col = 'black',lwd=3,main="Notification",
      ylab="TB notification per 100,000 population",xlab="Time in years",xlim = c(0,1600),ylim = c(0,200))

#Intervention
#targetting 99%of SSs
#prevalence 
svPrev=((svI0+svI1+svI2+svT)/svN)*100000
svPrev
max(svPrev)
isvPrev=((svI0+svI1+svI2+svT)/svN)*100000
max(isvPrev)
#plote prevalence
plot(svtime,svPrev,type='n', col = 'blue',lwd=3,main="Prevalence",
     ylab="Prevalence per 100,000 Popn",xlab="Time in years",xlim=c(0,1000),ylim = c(0,300))

lines(svtime,svPrev,type='l', col = 'blue',lwd=3)
lines(svtime,isvPrev,type='l', col = 'red',lwd=3)
legend('right',c('No intervention','Targeting Super-spreaders'),col=c('blue','red'),
       lty = c(1,1),lwd=c(3,3))
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#incidence the fraction of new cases at each time step 
BsIncidence=(diff(svinc)/svN)*100000 #baseline 65%
max(esIncidence)
M1sIncidence=(diff(svinc)/svN)*100000 # 66%
T1sIncidence=(diff(svinc)/svN)*100000 #75%

M2sIncidence=(diff(svinc)/svN)*100000 #untargeted 67%
T2sIncidence=(diff(svinc)/svN)*100000 #baseline 80%

M3sIncidence=(diff(svinc)/svN)*100000 #baseline 68%
T3sIncidence=(diff(svinc)/svN)*100000 #baseline 90%

M4sIncidence=(diff(svinc)/svN)*100000 #baseline 69%
T4sIncidence=(diff(svinc)/svN)*100000 #baseline 99%


plot (svtime, sIncidence, type='n', col = 'red',lwd=3,
      ylab="Incidence per 100,000 population",xlab="Time in years",xlim = c(0,50),ylim = c(0,300), 
      panel.first = grid(lty = 1, lwd = 1))
lines(svtime,esIncidence,type='l', col = 'red',lwd=3)
lines(svtime,isIncidence,type='l', col = 'darkviolet',lwd=3)

lines(svtime,TisIncidence,type='l', col = 'blue',lwd=3)
lines(svtime,tisIncidence,type='l', col = 'darkgreen',lwd=3)

legend(0,80,c('initial','mass 75%','90% superspreaders','99% superspreaders'),col=c('red','darkviolet','blue','darkgreen'),
       lty = c(1,1,1,1),lwd=c(3,3,3,3),cex = 0.75)


all_interventions<-cbind(svtime,BsIncidence,M1sIncidence,T1sIncidence,M2sIncidence,T2sIncidence,M3sIncidence,T3sIncidence,M4sIncidence,T4sIncidence)

View(all_interventions)
write.csv(x=all_interventions,file='..//Documents/intervention.csv')

data_intervention=read.csv(file.choose())
View(data_intervention)
par(mfrow=c(2,2))
plot (data_intervention$svtime, data_intervention$BsIncidence,type='n', col = 'red',lwd=3,
      ylab="Incidence per 100,000 population",xlab="Time in years",xlim = c(0,50),ylim = c(0,300), 
      panel.first = grid(lty = 1, lwd = 1),main="a) increase by 1%",cex = 0.1)
lines(data_intervention$svtime,data_intervention$BsIncidence,type='l', col = 'red',lwd=3)
lines(data_intervention$svtime,data_intervention$M1sIncidence,type='l', col = 'darkviolet',lwd=3)
lines(data_intervention$svtime,data_intervention$T1sIncidence,type='l', col = 'blue',lwd=3)
#legend(0,80,c('baseline','mass increase',' targeted increase '), col=c('red','darkviolet','blue'),   lty = c(1,1,1),lwd=c(3,3,3),cex = 0.75)
plot (data_intervention$svtime, data_intervention$BsIncidence, type='n', col = 'red',lwd=3,
      ylab="Incidence per 100,000 population",xlab="Time in years",xlim = c(0,50),ylim = c(0,300), 
      panel.first = grid(lty = 1, lwd = 1),main="b) increase by 2%")
lines(data_intervention$svtime,data_intervention$BsIncidence,type='l', col = 'red',lwd=3)
lines(data_intervention$svtime,data_intervention$M2sIncidence,type='l', col = 'darkviolet',lwd=3)
lines(data_intervention$svtime,data_intervention$T2sIncidence,type='l', col = 'blue',lwd=3)
#legend(0,80,c('baseline 65%','mass 67%','80% superspreaders'),   col=c('red','darkviolet','blue'),   lty = c(1,1,1),lwd=c(3,3,3),cex = 0.75)
plot (data_intervention$svtime, data_intervention$BsIncidence, type='n', col = 'red',lwd=3,
      ylab="Incidence per 100,000 population",xlab="Time in years",xlim = c(0,50),ylim = c(0,300), 
      panel.first = grid(lty = 1, lwd = 1),main="c) increase by 3%")
lines(data_intervention$svtime,data_intervention$BsIncidence,type='l', col = 'red',lwd=3)
lines(data_intervention$svtime,data_intervention$M3sIncidence,type='l', col = 'darkviolet',lwd=3)
lines(data_intervention$svtime,data_intervention$T3sIncidence,type='l', col = 'blue',lwd=3)
#legend(0,80,c('baseline 65%','mass 68%','90% superspreaders'), col=c('red','darkviolet','blue'),lty = c(1,1,1),lwd=c(3,3,3),cex = 0.75)
plot (data_intervention$svtime, data_intervention$BsIncidence, type='n', col = 'red',lwd=3,
      ylab="Incidence per 100,000 population",xlab="Time in years",xlim = c(0,50),ylim = c(0,300), 
      panel.first = grid(lty = 1, lwd = 1),main="d) increase by 4%")
lines(data_intervention$svtime,data_intervention$BsIncidence,type='l', col = 'red',lwd=3)
lines(data_intervention$svtime,data_intervention$M4sIncidence,type='l', col = 'darkviolet',lwd=3)
lines(data_intervention$svtime,data_intervention$T4sIncidence,type='l', col = 'blue',lwd=3)
#legend(0,80,c('baseline 65%','mass 69%','99% superspreaders'),col=c('red','darkviolet','blue'), lty = c(1,1,1),lwd=c(3,3,3),cex = 0.75)
#Each intervention senarios
only1=read.csv(file.choose()) #intervention_01_nt
only2=read.csv(file.choose())#intervention_02_nt
only3=read.csv(file.choose())#intervention_03_nt
only4=read.csv(file.choose())#intervention_04_nt
View(only2)
#FIGURE 6
library(reshape2)
library(ggplot2)
m_data_intervention=melt(data_intervention,id="svtime")
View(only1)
M_only1=melt(only1,id='svtime')
M_only2=melt(only2,id='svtime')
M_only3=melt(only3,id='svtime')
M_only4=melt(only4,id='svtime')
par(mfrow=c(2,2))
ggplot()+
  geom_line(data = M_only1,aes(x=svtime,y=value,colour=variable,linetype=variable),size=1.06)+
  geom_line(data = M_only2,aes(x=svtime,y=value,colour=variable,linetype=variable),size=1.06)+
  geom_line(data = M_only3,aes(x=svtime,y=value,colour=variable,linetype=variable),size=1.06)+
  geom_line(data = M_only4,aes(x=svtime,y=value,colour=variable,linetype=variable),size=1.06)+
  ylab("Incidence per 100,000 popn")+
  xlab("Time in years")+
  
  xlim(0,50)+
  ylim(0,300)+
 scale_colour_manual(values=c("black","yellow","red","blue","green",   "yellow","red","blue","green"),
                     name='CDR',labels=c("Baseline 65%","Mass 66%","Mass 67%","Mass 68%","Mass 69%","Targeted 75%","Targeted 80%","Targeted 90%","Targeted 99%"))+
  scale_linetype_manual(values = c("solid","dotted","dotted","dotted","dotted","solid","solid","solid","solid"),
                        name='CDR',labels=c("Baseline 65%","Mass 66%","Mass 67%","Mass 68%","Mass 69%","Targeted 75%","Targeted 80%","Targeted 90%","Targeted 99%"))
  
  
        

 
  
 

ggplot()+
  geom_line(data = M_only1,aes(x=svtime,y=value,linetype=variable),size=1.06,colour="black")+
  geom_line(data = M_only2,aes(x=svtime,y=value,linetype=variable),size=1.06,colour="blue")+
  geom_line(data = M_only3,aes(x=svtime,y=value,linetype=variable),size=1.06,colour="magenta")+
  geom_line(data = M_only4,aes(x=svtime,y=value,linetype=variable),size=1.06,colour="red")+
  ylab("Incidence per 100,000 popn")+
  xlab("Time in years")+
  
  xlim(0,50)+
  ylim(0,300)+
  theme(legend.justification = "top")+
  

  scale_fill_manual(values=c("red","black","blue","magenta","red",   "red","black","blue","magenta","red"),name='CDR',
                      labels=c("Baseline 65%","Mass 66%","Mass 67%","Mass 68%","Mass 69%","Targeted 75%","Targeted 80%","Targeted 90%","Targeted 99%"))


p1=ggplot(data = M_only1,aes(x=svtime,y=value,colour=variable))+
  geom_line(size=1.06)+
  ylab("Incidence per 100,000 popn")+
  xlab("Time in years")+
  ggtitle("a) 1% increase")+
  xlim(0,50)+
   ylim(0,300)+
  theme(legend.position = "bottom")+
  scale_colour_manual(values=c("red","purple","blue"),name='',
                      labels=c("Baseline","Mass","Targeted"))

  
p2= ggplot(data = M_only2,aes(x=svtime,y=value,colour=variable))+
    geom_line(size=1.06)+
    ylab("Incidence per 100,000 popn")+
    xlab("Time in years")+
    ggtitle("b) 2% increase")+
    xlim(0,50)+
    ylim(0,300)+
    theme(legend.position = "bottom")+
   scale_colour_manual(values=c("red","purple","blue"),name='',
                      labels=c("Baseline","Mass","Targeted"))
p3= ggplot(data = M_only3,aes(x=svtime,y=value,colour=variable))+
    geom_line(size=1.06)+
    ylab("Incidence per 100,000 popn")+
    xlab("Time in years")+
    ggtitle("c) 3% increase")+
    xlim(0,50)+
    ylim(0,300)+
    theme(legend.position = "bottom")+
    scale_colour_manual(values=c("red","purple","blue"),name='',
                      labels=c("Baseline","Mass","Targeted"))

p4=ggplot(data = M_only4,aes(x=svtime,y=value,colour=variable))+
    geom_line(size=1.06)+
    ylab("Incidence per 100,000 popn")+
    xlab("Time in years")+
    ggtitle("d) 4% increase")+
    xlim(0,50)+
    ylim(0,300)+
    theme(legend.position = "bottom")+
    scale_colour_manual(values=c("red","purple","blue"),name='',
                      labels=c("Baseline","Mass","Targeted"))
  
library(ggpubr) 
ggarrange(p1, p2, p3, p4, ncol=2, nrow=2, common.legend = TRUE, legend="bottom")

 scale_colour_manual(values=c("red","green","blue"),name="Intervention",
     
                                      labels=c("Baseline","Mass","Targeted"))

 #one graph
 ggplot()+
   geom_line(data = M_only4,aes(x=svtime,y=value,colour=variable))+
   geom_line(data = M_only3,aes(x=svtime,y=value,colour=variable))+
   geom_line(data = M_only2,aes(x=svtime,y=value,colour=variable))+
   geom_line(data = M_only1,aes(x=svtime,y=value,colour=variable))+
   ylab("Incidence per 100,000 popn")+
   xlab("Time in years")+

   xlim(0,50)+
   ylim(0,300)+
   theme(legend.justification =  "top")+
   scale_colour_manual(values=c("red","purple","blue"),name='',
                       labels=c("Baseline","Mass","Targeted"))
#+++++++++++++++++++++++++++++++++++++++++++
#Incidence comparison homogenous and different level of heterogeneity and Subclass model
par(bg="grey")
#FIGURE 4
plot ( svtime,sIncidence,type='n',ylab = 'Incidence per 100,000 Popn.', xlab='Time in years'
       , ylim=c(0,3e+02) ,xlim=c(0,700))

lines(svtime,sIncidence,type='l', col = 'green',lwd=3)
lines(Incidence_homo~ vtime, data = output, type='l',lwd=3,  col = 'blue')


lines(HIncidence1 ~ Hvtime,lwd=3,   col = 'magenta')
lines(HIncidence0.5 ~ Hvtime,lwd=3,   col = 'red')
lines(HIncidence0.1 ~ Hvtime,lwd=3,   col = 'brown3')
lines(HIncidence0.05 ~ Hvtime,lwd=3,   col = 'chocolate4')


legend(10,250, box.lty=0,c("Stratified model","Homogeneous model","Heterogeneous,k=1","Heterogeneous,k=0.5","Heterogeneous,k=0.1","Heterogeneous,k=0.05"),
       lty=c(1,1,1,1,1,1),lwd=c(3,3,3,3,3,3),
       col=c("green","blue","magenta",'red','brown3','chocolate4')) 

aii_inc<-cbind(svtime,sIncidence,HIncidence1,HIncidence0.5,HIncidence0.1,HIncidence0.05)
View(aii_inc)



#==============================================================================================================


#=================================================================
##### Force of infection



k=10
I=c(0,1,2,3,4,5,6,7,8,9,10,15,20,30,40,50,60,70,80,90,100)
b=7
N=1000
lamb=k*log(1+b*I/k/N)
lambd0=b*I/N
print(lamb)
lamb_10=as.data.frame(lamb)

View(lambd)
View(lambd0)
plot(lambd)
write.csv(x=lamb_10,file='..//Documents/lamb_10.csv')
write.csv(x=lambd0,file='..//Documents/lambd0.csv')
lambdas=read.csv(file.choose())
View(lambdas)
plot(lambdas$I,lambdas$lambda_o,type = 'n',ylim=c(0,0.8),xlim =c(0,100),ylab = "force of infection",
     xlab="number of infectious, I", main = "comparison between forces of infection")
lines(lambdas$I,lambdas$lambda_o,type = 'l',col='blue',lwd=3)
#lines(lambdas$I,lambdas$lambda_k_10, type = 'l',col='red',lwd=3)
#lines(lambdas$I,lambdas$lambda_k_1, type = 'l',col='darkred',lwd=3)
lines(lambdas$I,lambdas$lambda_k_0.1, type = 'l',col='chocolate3',lwd=3)
lines(lambdas$I,lambdas$lambda_k_0.01, type = 'l',col='darkorchid',lwd=3)

legend("topleft", c("??_o","??_k(k=0.1)","??_k(k=0.01)"),lty=c(1,1,1),lwd=c(3,3,3),
       col=c("blue",'chocolate3','darkorchid'))



#par(mfrow=c(2,2))

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


#plot
plot (vS ~ vtime,type='n',ylab = 'Population size', main = 'TB infectious size',xlab='Time in years'
      , ylim=c(0,3e+02),xlim=c(0,700) )
lines(vI ~ vtime, data = output, type='l',lwd=2.5,  col = 'red')
lines(HvI ~ Hvtime, data = Houtput, type='l',lwd=2.5,  col = 'Blue')
legend("bottomright", c("Homogeneous","Heterogeneous,k=0.1"),lty=c(1,1),lwd=c(2.5,2.5),
       col=c("red","blue"))

#&** K for sub-clss
sub_beta<-read.csv(file.choose())
View(sub_beta)
library(MASS)
fitdistr(sub_beta$no_inf_index,'negative binomial')->NB.fit
NB.fit$estimate #give us 0.2 which is higher
confint(NB.fit)
#let us produce a negative binomial random numbe with k=0.1 and mu = Ro of the homogeneous population (Ro=1.26)
Randum_sub<-data.frame(rnegbin (1000,mu=7,theta=0.1))
View(Randum_sub)
write.csv(x=Randum_sub,file='..//Desktop/Randum_suby.csv')
table(Randum_sub)#give us 739 (73.9%)are non infectious
#cut point for superspreader
pois=rpois(100000,lambda=7)
quantile(pois, c(.25, .5, .75, .9, .95, .99))#4 is the cut point for beta of super spreading and there are 12.1% SS, thus the mean number after 4 is
#the mean beta for super spreaders is 10.75 and the mean number of beta for enfectious but non-superspreaders is 1.55 and there are 12.4% infectious non-SS I(0.755,0.124,0.121),(0,2,11)
