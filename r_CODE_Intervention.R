library(deSolve)
#SUBCLASS FOR INFECTIOUS GROUP

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
a=0.624 #proportion of Io
b=0.244#proportion of I1
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
case_detection_super_spreaders_over_1year=0.99


transmission_events_per_non_infectious_person_year=0
transmission_events_per_infectious_non_super_spreader_person_year=5
transmission_events_per_infectious_super_spreader_person_year=44
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
D=1e-3*a     #Infectious hosts proportion of extrapulmonary 15% of all infectious 
E=1e-3*b # Normal spreaders
F=1e-3*c # Super-spreader 10% of TB patients are super-spreaders
G=0 #Treated hosts
H=0 #incidence
I=0 # diagnosed 

parameter_list=c(N=N,mu=mu,beta0=beta0, beta1=beta1, epsilon0=epsilon0,epsilon1=epsilon1,epsilon2=epsilon2, 
                 kappa=kappa,gamma0=gamma0,gamma1=gamma1,gamma2=gamma2,nu0=nu0,nu1=nu1,nu2=nu2,
                 omega0=omega0, omega1=omega1,omega2=omega2,phi=phi, 
                 delta0=delta0,delta1=delta1,delta2=delta2)

#Initial state values for differential equasions
#initial_values=c(S=A-(D+E+F), L1=B, L2=C, I0=D,I1=E,I2=F,T=G,inc=D+E+F,notif=I)
#at Equilibrium
initial_values=c(S=0.4915928963, L1=0.0227268406, L2=0.4819654218, I0= 0.0017543290,I1=0.0006859876,I2=0.0003711080,T=0.0009034167,inc= 0.0017543290 +0.0006859876+ 0.0003711080,notif=0.0009034167)
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

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#estimate the proportionof I_0, I_1,I_2
#the mean number of infections per index (m)
d=1/(delta+gamma+mui+mu)
m=beta*d
k=0.1
print(m)

#--------------

library(MASS)
fitdistr(sub_beta$no_inf_index,'negative binomial')->NB.fit
NB.fit$estimate #give us 0.2 which is higher
confint(NB.fit)
#let us produce a negative binomial random numbe with k=0.1 and mu = Ro of the homogeneous population (Ro=1.26)
library(MASS)
Randum_sub<-data.frame(rnegbin (1000,mu=m,theta=k))
View(Randum_sub)
write.csv(x=Randum_sub,file='..//Desktop/Randum_subz.csv')
table(Randum_sub)#give us 739 (73.9%)are non infectious
#cut point for superspreader
pois=rpois(100000,lambda=m)
quantile(pois, c(.25, .5, .75, .9, .95, .99))#4 is the cut point for beta of super spreading and there are 12.1% SS, thus the mean number after 4 is
#the mean beta for super spreaders is 10.75 and the mean number of beta for enfectious but non-superspreaders is 1.55 and there are 12.4% infectious non-SS I(0.755,0.124,0.121),(0,2,11)

a=0.642 #proportion of Io
b=0.226#proportion of I1
c=0.132##proportion of I2
# the average number of infections per index(m) for each subcompartment
m_0=0
m_1=4.34
m_2=51.045
beta0=m_0/d#0
beta1=m_1/d#4
beta2=m_2/d#51


