
#Numerical solution for endemic equilibrium size of each compartments for The three Models

#x=[1,2,3,4,5]i.e x[S,La,Lb,I,T]
library(nleqslv)
#Heterogeneous Model
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

dispersion_parameter=0.1
transmission_events_per_infectious_person_year=7
beta=transmission_events_per_infectious_person_year #Effective contact rate 
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
k=dispersion_parameter
#zero finder function
f <- function(x){
  
  lambda=k*log(1+beta*x[4]/k/N)
  lambdad=0.21*lambda
  s=lambda+mu
  a=epsilon+mu+kappa
  b=nu+lambdad+mu
  i=delta+mui+mu+gamma
  t=mut+omega+phi+mu
  
  y<-numeric(5) 
  y[1]<- x[1]-((mu*N+mui*x[4]+mut*x[5])+phi*x[5])/s
  y[2]<-x[2]-(lambda*x[1]+lambdad*x[3])/a
  y[3]<-x[3]-(kappa*x[2]+gamma*x[4])/b
  y[4]<- x[4]-(epsilon*x[2]+omega*x[5]+nu*x[3])/i
  y[5]<-x[5]-delta*x[4]/t
  return(y)
}
x_start=c(0.3,0.2962,0.0018,0.002,0.4)
hsolution = nleqslv(x_start, f,control=list(ftol=1e-20))
#View(hsolution)
print(hsolution$x)*100000

sum(hsolution$x)
#+++++++++++++++++++
#Homogenoues Model
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
r=1-risk_reduction_after_reinfection #risk of progression after re-infection
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
f <- function(x){
  
  lambda=beta*x[4]/N
  lambdad=0.21*lambda
  s=lambda+mu
  a=epsilon+mu+kappa
  b=nu+lambdad+mu
  i=delta+mui+mu+gamma
  t=mut+omega+phi+mu
 # Pie=mu*N+mui*x[4]+mut*x[3]
  y<-numeric(5) 
  y[1]<- x[1]-((mu*N+mui*x[4]+mut*x[5])+phi*x[5])/s
   y[2]<-x[2]-(lambda*x[1]+lambdad*x[3])/a
   y[3]<-x[3]-(kappa*x[2]+gamma*x[4])/b
   y[4]<- x[4]-(epsilon*x[2]+omega*x[5]+nu*x[3])/i
  y[5]<-x[5]-delta*x[4]/t
 
  
  return(y)
}
x_start=c(0.3,0.2962,0.0018,0.002,0.4)

solution = nleqslv(x_start, f, control=list(ftol=1e-20))
#View(solution)

print(solution$x)
sum(solution$x)
print(solution$fvec)
sum=sum(solution$x)
print(sum)
############################################ Stratified Infectious compartments
library(nleqslv)
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
case_detection_super_spreaders_over_1year=0.65


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
delta=case_detection_over_1year*(gamma+mui+mu)/(1-case_detection_over_1year) #Case detection rate()

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

f <- function(x){
  
  lambda=(beta0*x[4]+beta1*x[5]+beta2*x[6])/N
  lambdad=0.21*lambda
  s=lambda+mu
  a=epsilon0+epsilon1+epsilon2+mu+kappa
  b=nu0+nu1+nu2+lambdad+mu
  i=delta0+mui+mu+gamma0
  j=delta1+mui+mu+gamma1
  h=delta2+mui+mu+gamma2
  t=mut+omega0+omega1+omega2+phi+mu
  
  y<-numeric(7) 
  y[1]<- x[1]-((mu*N+mui*(x[4]+x[5]+x[6])+mut*x[7])+phi*x[7])/s
  
  y[2]<-x[2]-(lambda*x[1]+lambdad*x[3])/a
  
  y[3]<-x[3]-(kappa*x[2]+gamma0*x[4]+gamma1*x[5]+gamma2*x[6])/b
  
  y[4]<- x[4]-(epsilon0*x[2]+omega0*x[7]+nu0*x[3])/i
  
  y[5]<- x[5]-(epsilon1*x[2]+omega1*x[7]+nu1*x[3])/j
  
  y[6]<- x[6]-(epsilon2*x[2]+omega2*x[7]+nu2*x[3])/h
    
  y[7]<-x[7]-(delta0*x[4]+delta1*x[5]+delta2*x[6])/t
  return(y)
}
x_start=c(0.537947,0.019486,0.4397222,0.001076628,0.0006459771,0.0004306514,0.00069151)
Ssolution = nleqslv(x_start, f,control=list(ftol=1e-20))
#View(Ssolution)
print(Ssolution$x)
sum(Ssolution$x)
