#basic reproductive number\R_o

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


#R_o=beta*nu(epsilon+kappa)/(delta+mui+mu+gamma)(epsilon+kappa+mu)(nu+mu)
p=(epsilon/(epsilon+kappa+mu)) +((kappa/(epsilon+kappa+mu))*(nu/(nu+mu))) 
D=1/(delta+mui+mu+gamma)


R_oh=(beta*D*p)
print(R_oh)

#for startified model
#proportions 
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
cdr=case_detection_non_infectious_over_1year
cdr0=case_detection_non_infectious_over_1year
cdr1=case_detection_infectious_nonsuper_spreaders_over_1year
cdr2=case_detection_super_spreaders_over_1year
delta=cdr*(gamma+mui+mu)/(1-cdr)
delta0=cdr0*(gamma+mui+mu)/(1-cdr0)
delta1=cdr1*(gamma+mui+mu)/(1-cdr1)
delta2=cdr2*(gamma+mui+mu)/(1-cdr2)
#Ro =ß_0 xP_(0 ) xD+ß_1 xP_(1 ) xD+ß_2 xP_(2 ) xD
D=1/(delta+mui+mu+gamma)
p_0=(epsilon0/(epsilon0+epsilon1+epsilon2+kappa+mu))+(kappa/(epsilon0+epsilon1+epsilon2+kappa+mu))*(nu0/(nu0+nu1+nu2+mu))
p_1=(epsilon1/(epsilon0+epsilon1+epsilon2+kappa+mu))+(kappa/(epsilon0+epsilon1+epsilon2+kappa+mu))*(nu1/(nu0+nu1+nu2+mu))
p_2=(epsilon2/(epsilon0+epsilon1+epsilon2+kappa+mu))+(kappa/(epsilon0+epsilon1+epsilon2+kappa+mu))*(nu2/(nu0+nu1+nu2+mu))


R_os=(beta0*p_0*D)+(beta1*p_1*D)+(beta2*p_2*D)
print(R_os)
#the mean number of infections per index (m)
d=1/(delta+gamma+mui+mu)
m=beta*d
print(m)
