#####################################################
# DETERMINISTIC MODEL OF SMA IN GENERAL POPULATION AND POST-DISCHARGE COHORT
#####################################################

run_repeat_mod<-function(inc_sma_genpop,eir,pmc,
                         output="summary",ft_um_pd0=0.5,ft_sev_pd0=0.5,
                         prop_hosp=0.5, prop_hosp_paton=0.5,
                         pmc_covs=c(0.765, 0.879,0.894), r_eir,
                         prob_um_plac, prob_um_pmc, max_beta_phi, 
                         w_scale_risk, w_shape_risk,
                         adm) {
  
  al_eff<-0.98
  ft_um_pd<- ft_um_pd0*al_eff
  ft_sev_pd<-ft_sev_pd0*al_eff
    
  ##### DHA PIP DRUG 2
  w_scale2<-28.1
  w_slope2<-4.4
  dp_eff<-0.99
  max_proph_time=50
  time<-seq(from=0,to=max_proph_time,by=1)
  p_protect2<- dp_eff*exp(-((time/w_scale2)^w_slope2))

  pmc_times<-c(2,6,10)*7-13  ## time from first dose
  
  n_pd<-26*7-14-7  ## number of pd compartments allowing for 2 weeks of AL, and cutting off last week.
  
  pp1<-pp2<-pp3<-rep(0,n_pd)
  pp1[pmc_times[1]:(pmc_times[1]+max_proph_time)]<-p_protect2
  pp2[pmc_times[2]:(pmc_times[2]+max_proph_time)]<-p_protect2
  pp3[pmc_times[3]:(pmc_times[3]+max_proph_time)]<-p_protect2

  cfr_sma_h<-0.074
  cfr_sma_c<-2*cfr_sma_h
  cfr_sev_other_h<-0.01  ## (5921+472)/(20988+621373) from WMR 2015 Kenya & Uganda
  cfr_sev_c<-2*cfr_sev_other_h
  
  
  
  if(!pmc)    pmc_covs<-rep(0,3)
  dt<-1
  ## daily probability of SMA and other severe disease in the community, non post discharge
  # total: inc sma genpop is the total incidence from the G compartment, incl unhospitalised.
  # calculate hospitalised incidence at the end to get estimates to match to hospitalised data.
  prob_sma_gen_d<-1-exp((-inc_sma_genpop/365)*dt)

  prob_sm_plac<-1-prob_um_plac
  prob_sm_pmc<-1-prob_um_pmc
  prob_sma_of_sm<-43.9/101.9  ## SMA out of SM readmissions in Kwambai et al
  dur_p<-15 # 2 days treatment seeking + 13 days of prophylaxis provided by uncomplicated malaria treatment.
  hosp_stay<-3 # days in hospital per admission
  dur_pd_al<-13 # for simplicity, and PMC starts on day 14
  
  
  # cumulative risk of readmission (severe) within first 2 weeks is 3% (Kwambai)
  # rearrange 1-e(-r*14)=0.03
  # assume prob SMA is same as in trial prob_sma_of_sm
  inc_sev_2week<- -log(1-0.03)/14
  prob_sev_2week_d<- 1-exp(-inc_sev_2week*dt)
  
  # risk over time, scaled by EIR and post discharge time. Risk over time params fitted from discharge.
  # Here PD tracking starts after 14 days.
  #beta_phi<-max_beta_phi*exp(-eir*r_eir)
  ### start at day 15 - was fitted for days 1-14 but here we start PD after 14 days.
  
  if(eir<=30) {  # no post discharge data>30 eir
    beta_phi<-max_beta_phi*exp(-eir*r_eir)*exp(-((15:(n_pd+14)/w_scale_risk)^w_shape_risk))
  } else {
    beta_phi<-max_beta_phi*exp(-30*r_eir)*exp(-((15:(n_pd+14)/w_scale_risk)^w_shape_risk))
  }
  
  ########### Used Kwambai prop SMA out of SM. See opoka also.
  ## prob of sma over time since discharge from last episode
  if(eir<=30) {
    inc_sma_plac<-beta_phi*eir*prob_sm_plac*prob_sma_of_sm  ## pppy
    inc_sma_pmc<-beta_phi*eir*prob_sm_pmc*prob_sma_of_sm  ## pppy. WIll only use when on PMC.
  } else {
    inc_sma_plac<-beta_phi*30*prob_sm_plac*prob_sma_of_sm  ## pppy
    inc_sma_pmc<-beta_phi*30*prob_sm_pmc*prob_sma_of_sm  ## pppy. WIll only use when on PMC.
  }
  ## prob of sma over time since discharge from last episode
  # placebo
  prob_sma_pd_time_plac<-1-exp((-inc_sma_plac/365)*dt)
  ## PMC arms. Here prob sm is lower when drug protection>1%
  prob_sma_pd_time_pmc_temp<-1-exp((-inc_sma_pmc/365)*dt)
  for(i in 1:length(pmc_times)) {
    curr_pp<-eval(parse(text=paste0("pp",i)))
    assign(paste0("prob_sma_pd_time_pmc",i),c(prob_sma_pd_time_plac[which(curr_pp<=0.01 & 1:length(curr_pp)<pmc_times[i])],
                             prob_sma_pd_time_pmc_temp[which(curr_pp>0.01 & 1:length(curr_pp)>=pmc_times[i])],
                             prob_sma_pd_time_plac[which(curr_pp<=0.01 & 1:length(curr_pp)>pmc_times[i])]))
  }
  
  ## prob of other type of sev malaria over time since discharge from last episode
  if(eir<=30) {
    inc_sev_other_plac<-beta_phi*eir*prob_sm_plac*(1-prob_sma_of_sm)  ## pppy
    inc_sev_other_pmc<-beta_phi*eir*prob_sm_pmc*(1-prob_sma_of_sm)  ## pppy
  } else {
    inc_sev_other_plac<-beta_phi*30*prob_sm_plac*(1-prob_sma_of_sm)  ## pppy
    inc_sev_other_pmc<-beta_phi*30*prob_sm_pmc*(1-prob_sma_of_sm)  ## pppy
  }
  ## prob of other type of sev malaria over time since discharge from last episode
  prob_sev_other_pd_time_plac<-1-exp((-inc_sev_other_plac/365)*dt)
  prob_sev_other_pd_time_pmc_temp<-1-exp((-inc_sev_other_pmc/365)*dt)
  
  for(i in 1:length(pmc_times)) {
    curr_pp<-eval(parse(text=paste0("pp",i)))
    assign(paste0("prob_sev_other_pd_time_pmc",i),c(prob_sev_other_pd_time_plac[which(curr_pp<=0.01 & 1:length(curr_pp)<pmc_times[i])],
                                                  prob_sev_other_pd_time_pmc_temp[which(curr_pp>0.01 & 1:length(curr_pp)>=pmc_times[i])],
                                                  prob_sev_other_pd_time_plac[which(curr_pp<=0.01 & 1:length(curr_pp)>pmc_times[i])]))
  }
  
  #prob_st_pd_time<-prob_sev_other_pd_time*ft_sev_pd 
  
  ## prob of uncomplicated malaria over time since discharge from last episode
  if(eir<=30) {
    inc_um_high_plac<-beta_phi*eir*prob_um_plac  ## pppy
    inc_um_high_pmc<-beta_phi*eir*prob_um_pmc  ## pppy
  } else {
    inc_um_high_plac<-beta_phi*30*prob_um_plac  ## pppy
    inc_um_high_pmc<-beta_phi*30*prob_um_pmc  ## pppy
  }
  ## prob of uncomplicated malaria over time since discharge from last episode
  prob_um_pd_time_plac<-1-exp((-inc_um_high_plac/365)*dt)  ## daily prob
  prob_um_pd_time_pmc_temp<-1-exp((-inc_um_high_pmc/365)*dt)
  for(i in 1:length(pmc_times)) {
    curr_pp<-eval(parse(text=paste0("pp",i)))
    assign(paste0("prob_um_pd_time_pmc",i),c(prob_um_pd_time_plac[which(curr_pp<=0.01 & 1:length(curr_pp)<pmc_times[i])],
                                             prob_um_pd_time_pmc_temp[which(curr_pp>0.01 & 1:length(curr_pp)>=pmc_times[i])],
                                             prob_um_pd_time_plac[which(curr_pp<=0.01 & 1:length(curr_pp)>pmc_times[i])]))
  }
  
  
  
  runTime<-1000
  N0<-1
  G<-rep(0,runTime)  ## gen pop
  G[1]<-N0
  H<-rep(0,runTime)  ## in hospital
  #D<-rep(0,runTime)  ## severe disease state, not hospitalised. Don't need to explicitly model as the process is identical to H and can follow these people through given the fixed time period.
  PDal<-rep(0,runTime)  ## on AL after hospital.
  inc_sma_all<-rep(0,runTime)  ## track inc sma all
  ### NB code only works if duration in each pd compartment same as dt
  PD2_c<-matrix(nrow=runTime,ncol=n_pd) # post SMA, no hospital (not eligible for PMC).
  PD<-matrix(nrow=runTime,ncol=n_pd) # post SMA, post hospital (eligible for PMC)
  Ut<-matrix(nrow=runTime,ncol=n_pd) # uncomplicated treated
  Svt<-matrix(nrow=runTime,ncol=n_pd) # other severe, treated
  ## track those in Pmc doses 1,2,3.
  ## If don't get dose 2, stay in 1, etc.
  ## For simplicity track all of these for the same amount of time up to max effect (114 days) even though people cannot enter e.g. Pmc2 until later on.
  Pmc1<-Pmc2<-Pmc3<-inc_sev_other_all<-matrix(nrow=runTime,ncol=n_pd)
  PD[1,]<-0
  PD2_c[1,]<-0
  Ut[1,]<-0
  Svt[1,]<-0
  Pmc1[1,]<-0
  Pmc2[1,]<-0
  Pmc3[1,]<-0
  inc_sev_other_all[1,]<-0
  inc_sev_other_all[,1]<-0
  N_check<-rep(NA,runTime)

  deaths<-0
  
  ##########################################
  # RUN MODEL OVER TIME
  ##########################################
  for(i in 2:runTime) {

    # move everyone in the postdischarge state along a day in time (i) and in postdischarge time (j). 
    for(j in 2:n_pd) {
      ## other types of severe malaria / malaria requiring hospitalisation - track in postdischarge groups. No need to track in the community because it will not change any states.
      inc_sev_other_all[i,j]<-PD[i-1,j-1]*prob_sev_other_pd_time_plac[j-1] +
        PD2_c[i-1,j-1]*prob_sev_other_pd_time_plac[j-1] +
        Pmc1[i-1,j-1]*prob_sev_other_pd_time_pmc1[j-1]*(1-pp1[j-1]) +
        Pmc2[i-1,j-1]*prob_sev_other_pd_time_pmc2[j-1]*(1-pp2[j-1]) +
        Pmc3[i-1,j-1]*prob_sev_other_pd_time_pmc3[j-1]*(1-pp3[j-1])
      
      # post SMA state, was not hospitalised.
      PD2_c[i,j]<- PD2_c[i-1,j-1]  - PD2_c[i-1,j-1]*prob_sma_pd_time_plac[j-1] - 
        PD2_c[i-1,j-1]*prob_um_pd_time_plac[j-1]*ft_um_pd -
        ((1-prop_hosp)*(1-cfr_sev_c)*ft_sev_pd + prop_hosp*(1-cfr_sev_other_h))*PD2_c[i-1,j-1]*prob_sev_other_pd_time_plac[j-1] -  # non fatal treated sev.
        ((1-prop_hosp)*cfr_sev_c +prop_hosp*cfr_sev_other_h)*PD2_c[i-1,j-1]*prob_sev_other_pd_time_plac[j-1] + # fatal hosp sev.
        ifelse((i>dur_p & j>dur_p),PD2_c[i-dur_p,j-dur_p]*prob_um_pd_time_plac[j-dur_p]*ft_um_pd,0) +
        ifelse((i>(dur_p+hosp_stay) & j>(dur_p+hosp_stay)),((1-prop_hosp)*(1-cfr_sev_c)*ft_sev_pd + prop_hosp*(1-cfr_sev_other_h))*PD2_c[i-(dur_p+hosp_stay),j-(dur_p+hosp_stay)]*prob_sev_other_pd_time_plac[j-(dur_p+hosp_stay)],0) 

      
      # post discharge. Have received PMC with probability acc to coverage.
      PD[i,j]<- PD[i-1,j-1] - PD[i-1,j-1]*prob_sma_pd_time_plac[j-1] - 
        PD[i-1,j-1]*prob_um_pd_time_plac[j-1]*ft_um_pd -
        ((1-prop_hosp)*(1-cfr_sev_c)*ft_sev_pd + prop_hosp*(1-cfr_sev_other_h))*PD[i-1,j-1]*prob_sev_other_pd_time_plac[j-1] -  # non fatal treated sev.
        ((1-prop_hosp)*cfr_sev_c +prop_hosp*cfr_sev_other_h)*PD[i-1,j-1]*prob_sev_other_pd_time_plac[j-1] + # fatal sev.
        ifelse((i>dur_p & j>dur_p),PD[i-dur_p,j-dur_p]*prob_um_pd_time_plac[j-dur_p]*ft_um_pd,0) +
        ifelse((i>(dur_p+hosp_stay) & j>(dur_p+hosp_stay)),((1-prop_hosp)*(1-cfr_sev_c)*ft_sev_pd + prop_hosp*(1-cfr_sev_other_h))*PD[i-(dur_p+hosp_stay),j-(dur_p+hosp_stay)]*prob_sev_other_pd_time_plac[j-(dur_p+hosp_stay)],0) -
        ifelse(j==pmc_times[2], pmc_covs[2]*PD[i-1,j-1],0) -
        ifelse(j==pmc_times[3], pmc_covs[3]*PD[i-1,j-1],0)
      
      # assume complete protection for a fixed 14 days in  treated cases.
      #  track treated uncomplicated + sev malaria as they will be protected.
      Pmc1[i,j]<- Pmc1[i-1,j-1]  - 
        Pmc1[i-1,j-1]*prob_sma_pd_time_pmc1[j-1]*(1-pp1[j-1]) -
        Pmc1[i-1,j-1]*prob_um_pd_time_pmc1[j-1]*ft_um_pd*(1-pp1[j-1]) -
        ((1-prop_hosp)*(1-cfr_sev_c)*ft_sev_pd + prop_hosp*(1-cfr_sev_other_h))*Pmc1[i-1,j-1]*prob_sev_other_pd_time_pmc1[j-1]*(1-pp1[j-1]) -  # non fatal treated sev.
        ((1-prop_hosp)*cfr_sev_c +prop_hosp*cfr_sev_other_h)*Pmc1[i-1,j-1]*prob_sev_other_pd_time_pmc1[j-1]*(1-pp1[j-1]) -
        ifelse(j==pmc_times[2], pmc_covs[2]*Pmc1[i-1,j-1],0) -
        ifelse(j==pmc_times[3], pmc_covs[3]*Pmc1[i-1,j-1],0) +
        ifelse((i>dur_p & j>dur_p),Pmc1[i-dur_p,j-dur_p]*prob_um_pd_time_pmc1[j-dur_p]*ft_um_pd*(1-pp1[j-dur_p]),0) +
        ifelse((i>(dur_p+hosp_stay) & j>(dur_p+hosp_stay)),((1-prop_hosp)*(1-cfr_sev_c)*ft_sev_pd + prop_hosp*(1-cfr_sev_other_h))*Pmc1[i-(dur_p+hosp_stay),j-(dur_p+hosp_stay)]*prob_sev_other_pd_time_pmc1[j-(dur_p+hosp_stay)]*(1-pp1[j-(dur_p+hosp_stay)]),0)
      
      Pmc2[i,j]<- Pmc2[i-1,j-1] + ifelse(j==pmc_times[2], pmc_covs[2]*PD[i-1,j-1],0) -
        Pmc2[i-1,j-1]*prob_sma_pd_time_pmc2[j-1]*(1-pp2[j-1]) -
        Pmc2[i-1,j-1]*prob_um_pd_time_pmc2[j-1]*ft_um_pd*(1-pp2[j-1]) -
        ifelse(j==pmc_times[3], pmc_covs[3]*Pmc2[i-1,j-1],0) -
        ((1-prop_hosp)*(1-cfr_sev_c)*ft_sev_pd + prop_hosp*(1-cfr_sev_other_h))*Pmc2[i-1,j-1]*prob_sev_other_pd_time_pmc2[j-1]*(1-pp2[j-1]) -  # non fatal treated sev.
        ((1-prop_hosp)*cfr_sev_c +prop_hosp*cfr_sev_other_h)*Pmc2[i-1,j-1]*prob_sev_other_pd_time_pmc2[j-1]*(1-pp2[j-1]) +
        ifelse(j==pmc_times[2], pmc_covs[2]*Pmc1[i-1,j-1],0) +
        ifelse((i>dur_p & j>dur_p),Pmc2[i-dur_p,j-dur_p]*prob_um_pd_time_pmc2[j-dur_p]*ft_um_pd*(1-pp2[j-dur_p]),0) +
        ifelse((i>(dur_p+hosp_stay) & j>(dur_p+hosp_stay)),((1-prop_hosp)*(1-cfr_sev_c)*ft_sev_pd + prop_hosp*(1-cfr_sev_other_h))*Pmc2[i-(dur_p+hosp_stay),j-(dur_p+hosp_stay)]*prob_sev_other_pd_time_pmc2[j-(dur_p+hosp_stay)]*(1-pp2[j-(dur_p+hosp_stay)]),0)
      
      Pmc3[i,j]<- Pmc3[i-1,j-1] + ifelse(j==pmc_times[3], pmc_covs[3]*PD[i-1,j-1],0) -
        Pmc3[i-1,j-1]*prob_sma_pd_time_pmc3[j-1]*(1-pp3[j-1]) -
        Pmc3[i-1,j-1]*prob_um_pd_time_pmc3[j-1]*ft_um_pd*(1-pp3[j-1]) -
        ((1-prop_hosp)*(1-cfr_sev_c)*ft_sev_pd + prop_hosp*(1-cfr_sev_other_h))*Pmc3[i-1,j-1]*prob_sev_other_pd_time_pmc3[j-1]*(1-pp3[j-1])  -  # non fatal treated sev.
        ((1-prop_hosp)*cfr_sev_c +prop_hosp*cfr_sev_other_h)*Pmc3[i-1,j-1]*prob_sev_other_pd_time_pmc3[j-1]*(1-pp3[j-1]) +
        ifelse(j==pmc_times[3], pmc_covs[3]*Pmc1[i-1,j-1],0) +
        ifelse(j==pmc_times[3], pmc_covs[3]*Pmc2[i-1,j-1],0) +
        ifelse((i>dur_p & j>dur_p),Pmc3[i-dur_p,j-dur_p]*prob_um_pd_time_pmc3[j-dur_p]*ft_um_pd*(1-pp3[j-dur_p]),0) +
        ifelse((i>(dur_p+hosp_stay) & j>(dur_p+hosp_stay)),((1-prop_hosp)*(1-cfr_sev_c)*ft_sev_pd + prop_hosp*(1-cfr_sev_other_h))*Pmc3[i-(dur_p+hosp_stay),j-(dur_p+hosp_stay)]*prob_sev_other_pd_time_pmc3[j-(dur_p+hosp_stay)]*(1-pp3[j-(dur_p+hosp_stay)]),0)
      
      
      Ut[i,j]<- Ut[i-1,j-1] + PD[i-1,j-1]*prob_um_pd_time_plac[j-1]*ft_um_pd +
        PD2_c[i-1,j-1]*prob_um_pd_time_plac[j-1]*ft_um_pd +
        Pmc1[i-1,j-1]*prob_um_pd_time_pmc1[j-1]*ft_um_pd*(1-pp1[j-1]) +
        Pmc2[i-1,j-1]*prob_um_pd_time_pmc2[j-1]*ft_um_pd*(1-pp2[j-1]) +
        Pmc3[i-1,j-1]*prob_um_pd_time_pmc3[j-1]*ft_um_pd*(1-pp3[j-1]) - 
        ifelse((i>dur_p & j>dur_p),
               PD[i-dur_p,j-dur_p]*prob_um_pd_time_plac[j-dur_p]*ft_um_pd +
                 PD2_c[i-dur_p,j-dur_p]*prob_um_pd_time_plac[j-dur_p]*ft_um_pd +
                 Pmc1[i-dur_p,j-dur_p]*prob_um_pd_time_pmc1[j-dur_p]*ft_um_pd*(1-pp1[j-dur_p]) +
                 Pmc2[i-dur_p,j-dur_p]*prob_um_pd_time_pmc2[j-dur_p]*ft_um_pd*(1-pp2[j-dur_p]) +
                 Pmc3[i-dur_p,j-dur_p]*prob_um_pd_time_pmc3[j-dur_p]*ft_um_pd*(1-pp3[j-dur_p]) ,0)
      
      
      # Also need to track other sev malaria amongst post-SMA group, as treatment will provide protection.
      Svt[i,j]<- Svt[i-1,j-1] + ((1-prop_hosp)*(1-cfr_sev_c)*ft_sev_pd + prop_hosp*(1-cfr_sev_other_h))*inc_sev_other_all[i,j]  - 
        ifelse((i>(dur_p+hosp_stay) & j>(dur_p+hosp_stay)),
               ((1-prop_hosp)*(1-cfr_sev_c)*ft_sev_pd + prop_hosp*(1-cfr_sev_other_h))*
                 inc_sev_other_all[i+1-(dur_p+hosp_stay),j+1-(dur_p+hosp_stay)],0)

    }  # end of j loop over post discharge time.
    
    ## all incident SMA cases
    inc_sma_all[i]<-G[i-1]*prob_sma_gen_d +
      sum(PD[i-1,]*prob_sma_pd_time_plac) +
      sum(PD2_c[i-1,]*prob_sma_pd_time_plac) +
      sum(Pmc1[i-1,]*prob_sma_pd_time_pmc1*(1-pp1)) +
      sum(Pmc2[i-1,]*prob_sma_pd_time_pmc2*(1-pp2)) +
      sum(Pmc3[i-1,]*prob_sma_pd_time_pmc3*(1-pp3))
    
    G[i]<-G[i-1] - G[i-1]*prob_sma_gen_d + PD[i-1,n_pd]*(1-prob_sma_pd_time_plac[n_pd]) +
      PD2_c[i-1,n_pd]*(1-prob_sma_pd_time_plac[n_pd]) +
      Ut[i-1,n_pd] + Svt[i-1,n_pd] + Pmc1[i-1,n_pd]*(1-prob_sma_pd_time_pmc1[n_pd]*(1-pp1[n_pd])) + 
      Pmc2[i-1,n_pd]*(1-prob_sma_pd_time_pmc2[n_pd]*(1-pp2[n_pd])) + 
      Pmc3[i-1,n_pd]*(1-prob_sma_pd_time_pmc3[n_pd]*(1-pp3[n_pd])) +
      #new births to replace fatal cases and maintain constant population
      ((1-prop_hosp)*cfr_sma_c +prop_hosp*cfr_sma_h)*inc_sma_all[i] +
      ((1-prop_hosp)*cfr_sev_c +prop_hosp*cfr_sev_other_h)*sum(inc_sev_other_all[i,2:n_pd]) 
    
    
    # Fixed days in hosp or of illness.
    H[i]<- H[i-1] + 
      ## SMA, hospitalised and survived OR not hospitalised and survived in the communuty
      ((1-prop_hosp)*(1-cfr_sma_c)+prop_hosp*(1-cfr_sma_h))*inc_sma_all[i] -
      ifelse(i<=hosp_stay,0,
             ((1-prop_hosp)*(1-cfr_sma_c)+prop_hosp*(1-cfr_sma_h))*inc_sma_all[i-hosp_stay])
    
    # those leaving hospital and surviving enter PDal. 
    # + severe cases in the community treated with AL.
    PDal[i]<-PDal[i-1] +  ifelse(i<=hosp_stay,0,
                                 ((1-prop_hosp)*(1-cfr_sma_c)+prop_hosp*(1-cfr_sma_h))*inc_sma_all[i-hosp_stay] #+
    ) -
      ifelse(i<=(dur_pd_al+hosp_stay),0,
             ((1-prop_hosp)*(1-cfr_sma_c)+prop_hosp*(1-cfr_sma_h))*inc_sma_all[i-hosp_stay-dur_pd_al] #+
      ) 
      
    
    
    #############
    # PD starts 14 days after discharge.
    # hospitalised for SMA, got PMC but did not take the first dose
    PD[i,1] <- ifelse(i<=(dur_pd_al+hosp_stay),0,
                      prop_hosp*(1-cfr_sma_h)*(1-pmc_covs[1])*inc_sma_all[i-hosp_stay-dur_pd_al])
    
    # not hospitalised for SMA, will not get any PMC doses
    PD2_c[i,1] <-ifelse(i<=(dur_pd_al+hosp_stay),0,
                          (1-prop_hosp)*(1-cfr_sma_c)*inc_sma_all[i-hosp_stay-dur_pd_al])  # never hospitalised. Postdischarge, no PMC.
    
    # PMC1 starts 14 days after discharge. Go directly from PDal (G with delay)
    Pmc1[i,1] <-ifelse(i<=(dur_pd_al+hosp_stay),0,
                       prop_hosp*(1-cfr_sma_h)*pmc_covs[1]*inc_sma_all[i-hosp_stay-dur_pd_al])
    # dose 2, 3 of PMC. No one enters these states immediately after discharge.
    Pmc2[i,1]<-0
    Pmc3[i,1]<-0
    
    # ppl cannot enter PD and get symptomatic malaria on same day, so is always zero.
    Ut[i,1]<-0
    Svt[i,1]<-0
    N_check[i]<-G[i]+ H[i] + PDal[i] + sum(PD[i,]) + sum(PD2_c[i,]) + 
      sum(Ut[i,]) + sum(Svt[i,]) +
      sum(Pmc1[i,]) +sum(Pmc2[i,]) + sum(Pmc3[i,]) 
    
  }

    deaths<- -log(1-(((1-prop_hosp)*cfr_sma_c +prop_hosp*cfr_sma_h)*inc_sma_all[i] +
    ((1-prop_hosp)*cfr_sev_c +prop_hosp*cfr_sev_other_h)*sum(inc_sev_other_all[i,2:n_pd])))
    
    
  ## inc sma at equilib. Convert from probability back to inc
  inc_sma_all_eq<- -log(1-inc_sma_all[runTime-1])
  
  ## inc sma hosp at equilib.
  inc_sma_all_hosp_eq<- -log(1-(prop_hosp*(G[runTime-1]*prob_sma_gen_d +
    sum(PD[runTime-1,]*prob_sma_pd_time_plac) +
    sum(PD2_c[runTime-1,]*prob_sma_pd_time_plac) +
    sum(Pmc1[runTime-1,]*prob_sma_pd_time_pmc1*(1-pp1)) +
    sum(Pmc2[runTime-1,]*prob_sma_pd_time_pmc2*(1-pp2)) +
    sum(Pmc3[runTime-1,]*prob_sma_pd_time_pmc3*(1-pp3)))))
  
  ## inc sma hosp in Paton at equilib. If prop hosp is different from Paton.
  # What would hospitalised incidence have been if using Paton hospitalisation probability
  inc_sma_all_paton_eq<- -log(1-(prop_hosp_paton*(G[runTime-1]*prob_sma_gen_d +
                                                sum(PD[runTime-1,]*prob_sma_pd_time_plac) +
                                                sum(PD2_c[runTime-1,]*prob_sma_pd_time_plac) +
                                                sum(Pmc1[runTime-1,]*prob_sma_pd_time_pmc1*(1-pp1)) +
                                                sum(Pmc2[runTime-1,]*prob_sma_pd_time_pmc2*(1-pp2)) +
                                                sum(Pmc3[runTime-1,]*prob_sma_pd_time_pmc3*(1-pp3)))))
  
  ## Separate how many episodes are occurring post-SMA versus in general popn
  num_sma_gen_pppd <- -log(1- (G[runTime-1]*prob_sma_gen_d))
  num_sma_pd_pppd <- -log(1-(sum(PD[runTime-1,]*prob_sma_pd_time_plac) +
                               sum(PD2_c[runTime-1,]*prob_sma_pd_time_plac) +
                               sum(Pmc1[runTime-1,]*prob_sma_pd_time_pmc1*(1-pp1)) +
                               sum(Pmc2[runTime-1,]*prob_sma_pd_time_pmc2*(1-pp2)) +
                               sum(Pmc3[runTime-1,]*prob_sma_pd_time_pmc3*(1-pp3))))
  
  num_sev_oth_pd_pppd <- -log(1-(sum(PD[runTime-1,]*prob_sev_other_pd_time_plac) +
                                   sum(PD2_c[runTime-1,]*prob_sev_other_pd_time_plac) +
                                   sum(Pmc1[runTime-1,]*prob_sev_other_pd_time_pmc1*(1-pp1)) +
                               sum(Pmc2[runTime-1,]*prob_sev_other_pd_time_pmc2*(1-pp2)) +
                               sum(Pmc3[runTime-1,]*prob_sev_other_pd_time_pmc3*(1-pp3))))

    ## Incidence post SMA per person (among post SMA group)
  inc_sma_pd_eq<- -log(1-((sum(PD[runTime-1,]*prob_sma_pd_time_plac) +
                             sum(PD2_c[runTime-1,]*prob_sma_pd_time_plac) +
                             sum(Pmc1[runTime-1,]*prob_sma_pd_time_pmc1*(1-pp1)) +
                             sum(Pmc2[runTime-1,]*prob_sma_pd_time_pmc2*(1-pp2)) +
                             sum(Pmc3[runTime-1,]*prob_sma_pd_time_pmc3*(1-pp3))))/
                         (sum(PD[runTime-1,]) + sum(PD2_c[runTime-1,]) + sum(Pmc1[runTime-1,]) +
                            sum(Pmc2[runTime-1,]) + sum(Pmc3[runTime-1,]) +sum(Ut[runTime-1]) +
                                                          sum(Svt[runTime-1])))

  inc_sev_oth_pd_eq<- -log(1-((sum(PD[runTime-1,]*prob_sev_other_pd_time_plac) +
                                 sum(PD2_c[runTime-1,]*prob_sev_other_pd_time_plac) +
                                 sum(Pmc1[runTime-1,]*prob_sev_other_pd_time_pmc1*(1-pp1)) +
                             sum(Pmc2[runTime-1,]*prob_sev_other_pd_time_pmc2*(1-pp2)) +
                             sum(Pmc3[runTime-1,]*prob_sev_other_pd_time_pmc3*(1-pp3))))/
                         (sum(PD[runTime-1,]) + sum(PD2_c[runTime-1,]) + sum(Pmc1[runTime-1,]) +
                            sum(Pmc2[runTime-1,]) + sum(Pmc3[runTime-1,])+sum(Ut[runTime-1]) +
                                                          sum(Svt[runTime-1])))
  
  inc_sma_1_3_pd_eq<- -log(1-((sum(PD[runTime-1,1:84]*prob_sma_pd_time_plac[1:84]) +
                                 sum(PD2_c[runTime-1,1:84]*prob_sma_pd_time_plac[1:84]) +
                                 sum(Pmc1[runTime-1,1:84]*prob_sma_pd_time_pmc1[1:84]*(1-pp1[1:84])) +
                                 sum(Pmc2[runTime-1,1:84]*prob_sma_pd_time_pmc2[1:84]*(1-pp2[1:84])) +
                                 sum(Pmc3[runTime-1,1:84]*prob_sma_pd_time_pmc3[1:84]*(1-pp3[1:84]))))/
                             (sum(PD[runTime-1,1:84])+ sum(PD2_c[runTime-1,1:84])+sum(Pmc1[runTime-1,1:84])+
                                sum(Pmc2[runTime-1,1:84])+sum(Pmc3[runTime-1,1:84])+
                                sum(Ut[runTime-1,1:84]) + sum(Svt[runTime-1,1:84])))
  
  inc_sev_oth_1_3_pd_eq<- -log(1-((sum(PD[runTime-1,1:84]*prob_sev_other_pd_time_plac[1:84]) +
                                     sum(PD2_c[runTime-1,1:84]*prob_sev_other_pd_time_plac[1:84]) +
                                     sum(Pmc1[runTime-1,1:84]*prob_sev_other_pd_time_pmc1[1:84]*(1-pp1[1:84])) +
                                 sum(Pmc2[runTime-1,1:84]*prob_sev_other_pd_time_pmc2[1:84]*(1-pp2[1:84])) +
                                 sum(Pmc3[runTime-1,1:84]*prob_sev_other_pd_time_pmc3[1:84]*(1-pp3[1:84]))))/
                             (sum(PD[runTime-1,1:84])+ sum(PD2_c[runTime-1,1:84])+
                                sum(Pmc1[runTime-1,1:84])+
                                sum(Pmc2[runTime-1,1:84])+sum(Pmc3[runTime-1,1:84])+
                                sum(Ut[runTime-1,1:84]) + sum(Svt[runTime-1,1:84])))
  ### number of PMC given out per day if perfect coverage.
  ### Include those who don't adhere to PMC (as they will still be given it)
  num_pmc_eq<-PD[runTime,1] + Pmc1[runTime,1] + Pmc2[runTime,1] + Pmc3[runTime,1]
  
  
  ### Total annual incidence sev
  if(output=="rds" | output=="summary2R") {
    summ<-c(inc_sma_all_eq=365*inc_sma_all_eq, 
            inc_sma_all_hosp_eq=365*inc_sma_all_hosp_eq,
            inc_sma_all_paton_eq=365*inc_sma_all_paton_eq,
            inc_sma_pd_eq=365*inc_sma_pd_eq,
            inc_sev_oth_pd_eq=365*inc_sev_oth_pd_eq,
            inc_sma_1_3_pd_eq=365*inc_sma_1_3_pd_eq,
            inc_sev_oth_1_3_pd_eq=365*inc_sev_oth_1_3_pd_eq,
            num_pmc_eq=365*num_pmc_eq, 
            num_sma_gen_pppy=365*num_sma_gen_pppd,
            num_sma_pd_pppy=365*num_sma_pd_pppd,
            num_sev_oth_pd_pppy=365*num_sev_oth_pd_pppd,
            deaths_pppy=365*deaths)
    if(output=="rds")
      saveRDS(summ,
            file=paste0("out_summ_adm",adm,".rds"))
    if(output=="summary2R")
      return(summ)
    } else {  
    return(data.frame(time=1:runTime,N=N_check,G=G, H=H, PDal=PDal, PD=rowSums(PD),
                      PD2_c=rowSums(PD2_c), Ut=rowSums(Ut),   
                      Svt=rowSums(Svt),Pmc1=rowSums(Pmc1),Pmc2=rowSums(Pmc2),
                      Pmc3=rowSums(Pmc3)))
  }
}




