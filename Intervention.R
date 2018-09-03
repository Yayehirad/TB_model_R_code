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

case_detection_non_infectious_over_1year=0.65
case_detection_infectious_nonsuper_spreaders_over_1year=0.65
case_detection_super_spreaders_over_1year=0.99


transmission_events_per_non_infectious_person_year=0
transmission_events_per_infectious_non_super_spreader_person_year=6
transmission_events_per_infectious_super_spreader_person_year=20
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

nu0=nu*0.5
nu1=nu*0.3
nu2=nu*0.2

epsilon0= epsilon*0.5
epsilon1= epsilon*0.3
epsilon2= epsilon*0.2


omega0=omega*0.5
omega1=omega*0.3
omega2=omega*0.2

cdr0=case_detection_non_infectious_over_1year
cdr1=case_detection_infectious_nonsuper_spreaders_over_1year
cdr2=case_detection_super_spreaders_over_1year

delta0=cdr0*(gamma+mui+mu)/(1-cdr0)
delta1=cdr1*(gamma+mui+mu)/(1-cdr1)
delta2=cdr2*(gamma+mui+mu)/(1-cdr2)


parameter_list=c(N=N,mu=mu,beta0=beta0, beta1=beta1, epsilon0=epsilon0,epsilon1=epsilon1,epsilon2=epsilon2, 
                 kappa=kappa,gamma0=gamma0,gamma1=gamma1,gamma2=gamma2,nu0=nu0,nu1=nu1,nu2=nu2,
                 omega0=omega0, omega1=omega1,omega2=omega2,phi=phi, 
                 delta0=delta0,delta1=delta1,delta2=delta2)

#Initial state values for differential equasions

#at Equilibrium
initial_values=c(S=0.5621359845, L1=0.0186135753, L2=0.4158730485, I0= 0.0012733407,I1=0.0007640044,I2=0.0005093363,T=0.0008307103,inc= 0.0012733407+ 0.0007640044 +0.0005093363,notif=0.0008307103)
## Output Time points
times=seq(0, 200, by = 1)
#Simulate the S1S2L1L2IT transmission 
soutput=as.data.frame(lsoda(initial_values, times, subclass_model, parameter_list))
View(soutput)
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
#incidence with different CDR 
esIncidence=(diff(svinc)/svN)*100000#baseline 65%
View(esIncidence)
isIncidence=(diff(svinc)/svN)*100000 #untargeted 75%
View(isIncidence)
TisIncidence=(diff(svinc)/svN)*100000 #90% SS
View(TisIncidence)
tisIncidence=(diff(svinc)/svN)*100000#target 99% SS
View(tisIncidence)
#plote change in incidence through time
plot (svtime, esIncidence, type='n', col = 'red',lwd=3,main="Intervention comparison",
      ylab="Incidence per 100,000 population",xlab="Time in years",xlim = c(0,50),ylim = c(0,250),panel.first = grid(lty = 3, lwd = 0.5))
lines(svtime,esIncidence,type='l', col = 'red',lwd=3)
lines(svtime,isIncidence,type='l', col = 'darkviolet',lwd=3)
lines(svtime,TisIncidence,type='l', col = 'blue',lwd=3)
lines(svtime,tisIncidence,type='l', col = 'darkgreen',lwd=3)

legend(0,80,c('initial CDR of 65%','mass CDR increase to 75%','90% CDR towards superspreaders','99% CDR towards superspreaders'),col=c('red','darkviolet','blue','darkgreen'),
       lty = c(1,1,1,1),lwd=c(3,3,3,3),cex = 0.75)

