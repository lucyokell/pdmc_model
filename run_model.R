
run_pmc_model<-function(n_cohort=50000,N=1000,n_years=4,adm=999, EIR_d_annual=NULL,pmc_cov=0, pmc_drug="dp", 
                        pmc_weeks=c(2,6,10), 
                        start_measure= 3*7, # time start measuring outcome, since discharge.
                        end_measure=14*7,EIR_fully_specified=NULL,save_sev_track=FALSE,save_age_inc=FALSE,
                        save_summary_over_time=FALSE, save_rel_foi=FALSE,
                        burnin_d=50,
                        output_root="output/",min_age_cohort=4*30.5,max_age_cohort=5*365)
  {
  ######## SIMULATE PMC
  ## to do update demog.
  ## PMC coverage is random per dose.
  
  ##############################################################
  ## set parameters
  ##############################################################
  
  # infection blocking immunity parameters
  bh	=	0.590076		
  bmin	=	.5
  db	=	3650
  kb	=	2.15506
  ub	=	7.19919
  IB0	<- 43.8787
  
  ## exposure by age
  rho	<-	0.85
  a0	<-	2920
  eta <- 1/(21*365) # death rate in humans
  max_age<-60  # set max age to avoid unrealistically old people.
  
  # clinical immunity parameters
  phi0	<-	0.791666		
  phi1	<-	0.000737
  dc	<-	10950
  kc	<-	2.36949
  uc	<-	6.06349
  P_IC_M	<-	0.774368
  dm	<-	67.6952
  
  IC0	<- 18.02366
  
  # severe immunity
  IV0	=	1.054525
  theta0	=	0.0749265
  theta1	=	0.000119
  dv	=	10950
  kv	=	1.9944
  uv	=	11.3096
  fv0	=	0.140045
  av00	=	2484.49
  gammav	=	2.91793
  P_IV_M	=	0.196823
  dvm	=	76.5871
  
  prob_sev_fatal<-0.06
  
  q0_sma <-0.548   ##q0 determines starting prevalence for age=0 (theoretically)
  q1_sma <-0.178   ##q1 determines the ending prevalence at oldest age
  r_sma <-0.252    ##r determines the exponential growth- how fast does it reach the ending prevalence
  
  
  
  # incubation period
  dur_E	=	12			#{the duration that you stay in the Exposed compartment}
  
  # daily non malaria death rate
  eta<-1/(21*365)
  
  # PMC parameters
  hospital_stay<-4 ## days. https://www.thelancet.com/journals/langlo/article/PIIS2214-109X(19)30345-6/fulltext Take quinine or artesunate during this time.
  pmc_times<-pmc_weeks*7 + hospital_stay  # time since admission.
  start_pmc<-burnin_d

  delay_sev_boost<-FALSE
  delay_boost_time<-6*30.5
  max_fup<-6*30.5  # max follow up time after PMC
  # Ignore delay to hospital at the moment.
  prob_hosp<-1  ## probability get to hospital if have severe malaria
  # Take AL at end of hospital stay.
  # Take DP week 2,6,10 after discharge. 100% coverage at the moment.
  # check if beginning or end of week?
  
  # Mechanisms of recurrence.
  fixed_theta<-FALSE
  
  
  
  #######################################################################
  # IMMUNITY BURN-IN
  # FIRST RUN A COHORT FROM AGE 0-20, TRACK THEIR IMMUNITY THROUGHOUT
  #######################################################################
  #IB_test5d. Add severe immunity 
  # allow multiple bites per person. 
  days<-20*365  ## max age
  dt_5d<-0.5   ## dt in age, days.
  
  #### EIR per time step for the duration of the run:
  if(!is.null(EIR_d_annual)) {
    EIR_dt_annual<-rep(EIR_d_annual,each=1/dt_5d)
  } 
  if(!is.null(EIR_fully_specified)) {  ### take the first year of EIR for immunity burnin_d
    EIR_dt_annual<-rep(EIR_fully_specified[1:365],each=1/dt_5d)
  }
  EIR_test<-rep(EIR_dt_annual,20)
  
  ### Set up matrices to store immunity over time per person
  IB_test5d_i<-ICA_test5d_i<-IVA_test5d_i<-IVM_test5d_i<-matrix(nrow=N,ncol=days/dt_5d)
  ### Initialise immunity at 0 at age 0
  IB_test5d_i[,1]<-ICA_test5d_i[,1]<-IVA_test5d_i[,1]<-rep(0,N)
  ### Set up matrices to record time each type of immunity boosted
  time_ib_boost<-time_ica_boost<-time_sev_boost<-matrix(nrow=N,ncol=days/dt_5d)
  ### Initialise immune boost time at large negative number at age=0
  time_ib_boost[,1]<-time_ica_boost[,1]<-time_sev_boost[,1]<-rep(-1000,N)
  ### Precompute severe immunity variable which depends only on age.
  fva<- 1- (1-fv0)/(1+(dt_5d*(1:(days/dt_5d))/av00)^gammav) 
  #plot(1:(days/dt_5d),fva,type="l") # check sensible shape
  
  ##### Heterogeneity
  # rlnorm params are the mean and sd of the log distribution 
  sigma2	<-	1.67 
  temp_rel_foi<- rlnorm(N, meanlog = -sigma2/2, sdlog = sqrt(sigma2))   # mean of log normal is exp(mu+sigma^2/2), so rearrange to get mean=1.
  #mean(temp_rel_foi)
  rel_foi<-rep(1,N)
  het<-T
  if(het) rel_foi<-temp_rel_foi/mean(temp_rel_foi)  # normalise.
  #mean(rel_foi)   # check mean should be 1
  
  
  #### Run the cohort age 0-20
  for(j in 2:(days/dt_5d)) {  # days are age. 
    
    ### Update variables
    IB_test5d_i[,j]<-IB_test5d_i[,j-1]*exp(-1/db*dt_5d)  #update daily decline
    ICA_test5d_i[,j]<-ICA_test5d_i[,j-1]*exp(-1/dc*dt_5d)  #update daily decline
    IVA_test5d_i[,j]<-IVA_test5d_i[,j-1]*exp(-1/dv*dt_5d)  #update daily decline
    time_ib_boost[,j]<-time_ib_boost[,j-1]  ## carry over boost times
    time_ica_boost[,j]<-time_ica_boost[,j-1]  ## carry over boost times
    time_sev_boost[,j]<-time_sev_boost[,j-1]  ## carry over boost times
    ## just track b_i as a vector per person. Not a matrix over time right now.
    b_i<- bh*((1-bmin)/(1+(IB_test5d_i[,j]/IB0)^kb) + bmin) # update b. Minimum b is bmin*bh (not bmin alone.)
    foi_age<-1-rho*exp(-(j*dt_5d)/a0) # update foi_age with current age
    
    ### Allocate bites, boost infection blocking immunity
    nbites<-rpois(N,EIR_test[j]*foi_age*rel_foi*dt_5d) # draw daily infectious bites per person
    if(anyNA(nbites)) {
      print(j)
      print(EIR_test[j])
    }
    is.bite<-which(nbites>0)  ## index of who gets infectious bite
    is.bite.boost<-is.bite[which((j*dt_5d-time_ib_boost[is.bite,j])>ub)] # which current time since last boost is more than ub
    time_ib_boost[is.bite.boost,j] <-j*dt_5d  ### record time of boosting infection blocking immunity
    IB_test5d_i[is.bite.boost,j]<-IB_test5d_i[is.bite.boost,j] + 1  ## increase infection blocking immunity by one for those who are bitten and boosted
    
    ### Draw infections, boost clinical immunity
    is.inf<-is.bite[which(runif(length(is.bite))<b_i[is.bite])]
    is.inf.boost<-is.inf[which((j*dt_5d-time_ica_boost[is.inf,j])>uc)]
    ICA_test5d_i[is.inf.boost,j]<-ICA_test5d_i[is.inf.boost,j] + 1
    time_ica_boost[is.inf.boost,j] <-j*dt_5d
    
    ### Boost severe malaria immunity
    is.sev.boost<-is.inf[which((j*dt_5d-time_sev_boost[is.inf,j])>uv)]
    IVA_test5d_i[is.sev.boost,j]<-IVA_test5d_i[is.sev.boost,j] + 1
    time_sev_boost[is.sev.boost,j] <-j*dt_5d
  }
  ## add maternal severe immunity + decay.
  IVM_test5d_i[,1]<-IVA_test5d_i[,20*365/dt_5d] * P_IV_M  ## each person aged 20.
  ivm_decay<-function(x) {x*exp(-1/dvm*2:ncol(IVM_test5d_i)*dt_5d)}
  for(i in 1:nrow(IVM_test5d_i)) {
    IVM_test5d_i[i,2:ncol(IVM_test5d_i)] <- ivm_decay(IVM_test5d_i[i,1])
  }
  ## Compute infection blocking and severe immunity over time.
  b_i<-theta_i<-matrix(nrow=N,ncol=days/dt_5d)
  b_i<- bh*((1-bmin)/(1+(IB_test5d_i/IB0)^kb) + bmin) # update b. Minimum b is bmin*bh (not bmin alone.)
  computeTheta<-function(iva,ivm) {
    theta0*(theta1 + (1-theta1)/(1+fva*((iva+ivm)/IV0)^gammav ))  ## Check brackets
  }
  for(i in 1:nrow(theta_i)) theta_i[i,]<-computeTheta(IVA_test5d_i[i,],IVM_test5d_i[i,])

  ######################################################################################
  # RUN PMC SIMULATION OVER TIME
  ######################################################################################
  
  ## sample some individuals at random from the above simulation - randomly by individual and weighted by exponential age structure:
  inds_indiv<-round(runif(n_cohort,min=1,max=N))
  possible_inds_age<-((min_age_cohort+1)/dt_5d):(max_age_cohort/dt_5d)
  inds_age<-sample(possible_inds_age,size=n_cohort,replace=T,prob=exp(-eta*possible_inds_age*dt_5d))
  inds<-cbind(inds_indiv,inds_age)
  cohort<-data.frame(age_d=inds_age*dt_5d,
                     IB=IB_test5d_i[inds], 
                     IVA=IVA_test5d_i[inds],
                     IVM=IVM_test5d_i[inds],
                     b=b_i[inds], 
                     theta=theta_i[inds],
                     age_ib_boost=time_ib_boost[inds],
                     age_sev_boost=time_sev_boost[inds],
                     rel_foi=rel_foi[inds_indiv]/mean(rel_foi[inds_indiv]),  ##re-normalise for this population
                     last_sev_time=-1000,  #time last had a severe episode
                     last_sma_time=-1000,  #time last had a SMA episode
                     hosp_time=-1000, # time of first going to hospital.
                     last_drug_time=-1000, #timing of AL or DP at home or at discharge
                     last_drug_type=NA)   # type of PMC
  cohort<-cohort %>%
    dplyr::mutate(fva= 1- (1-fv0)/(1+(cohort$age_d/av00)^gammav),
                  foi_age=1-rho*exp(-(age_d)/a0),
                  prob_sma=q0_sma+ (q1_sma-q0_sma) * (1-exp(-r_sma*(age_d/365))))
  
  dt<-1
  if(!is.null(EIR_d_annual)) {
    EIR_dt_annual<-rep(EIR_d_annual,each=1/dt)
    eir_d<-rep(EIR_dt_annual,n_years)
  }
  if(!is.null(EIR_fully_specified)) {
    eir_d<-rep(EIR_fully_specified,each=1/dt)
  }
  
  ### Timings. Infection becomes severe a fixed 12 days after bite (if it's going to).
  
  ##### AL DRUG 1. 
  # max proph time
  max_proph_time<-50
  w_scale1<-93.5*0.139
  w_slope1<-93.5
  time<-seq(from=0,to=max_proph_time,by=dt)
  p_protect1<- exp(-((time/w_scale1)^w_slope1))
  #plot(time,p_protect1)
  
  ##### DHA PIP DRUG 2
  w_scale2<-28.1
  w_slope2<-4.4
  time<-seq(from=0,to=max_proph_time,by=dt)
  p_protect2<- exp(-((time/w_scale2)^w_slope2))
  #plot(time,p_protect2)
  
  
  ## create empty matrix of people x days to track severe episodes.
  # Record hospital stays and PMC in here too for now.
  sev_track<-matrix(rep(-1,nrow(cohort)*length(eir_d)),nrow=nrow(cohort),ncol=length(eir_d))
  # track totals for different episodes
  tot_bites<-0
  tot_sev<-0 # incidence of severe
  tot_sma<-0 # incidence of SMA
  tot_inf_no_pmc<-0 # incidence of infection that would have happened without PMC
  tot_bites_postd<-0  # total infections during postdischarge period between start_measure and end measure
  tot_inf_postd<-0  # total infections during postdischarge period between start_measure and end measure
  tot_sev_postd<-0  # total severe disease during postdischarge period between start_measure and end measure
  tot_protect_start_end<-0 # number of infections prevented by PMC.
  person_time_postd<-0
  sev_inc_age<-NULL
  rel_foi_sev<-NULL
  rel_foi_sma<-NULL
  
  res_over_time<-NULL
  if(save_summary_over_time) {
    tot_weeks=floor(length(eir_d)*dt/7)
    res_over_time<-data.frame(weeks=1:tot_weeks,n_sma=0,persondays_all=rep(n_cohort*7),n_reinf_pd=0,persondays_pd=0)
  }
  
  for(j in 1:length(eir_d)) {  # days are TIME 
  #for(j in 101:200) {  # days are TIME 
    
    ### Update age, immunity
    cohort<-cohort %>%  
      dplyr::mutate(age_d=age_d+dt,
                    IB=IB*exp(-1/db*dt),
                    IVA=IVA*exp(-1/dv*dt),
                    IVM=IVM*exp(-1/dvm*dt),
                    b=bh*((1-bmin)/(1+(IB/IB0)^kb) + bmin),
                    fva= 1- (1-fv0)/(1+(age_d/av00)^gammav),
                    theta= theta0*(theta1 + (1-theta1)/(1+fva*((IVA+IVM)/IV0)^gammav)),  ## add IVM. Check brackets
                    foi_age=1-rho*exp(-(age_d)/a0))
    
    # who is currently in postdischarge state
    is.postd<-NULL
    is.postd<-which(j*dt - (cohort$last_sma_time +hospital_stay) < end_measure &
                                     j*dt - (cohort$last_sma_time +hospital_stay) > start_measure)
    
    ### Modify severe immunity acc to assumption
    if(fixed_theta) {
      inds<-which(j*dt - cohort$last_sev_time <6*30.5)
      cohort$theta[inds]<-theta0
    }
    
    ### New infectious bites, infection blocking immunity
    nbites<-rpois(n_cohort,eir_d[j]*cohort$foi_age*cohort$rel_foi*dt) # infectious bites per person this time step
    is.bite<-which(nbites>0)
    if(j/dt>burnin_d) tot_bites<-tot_bites + length(is.bite)
    is.bite.boost<-is.bite[which((cohort$age_d[is.bite]-cohort$age_ib_boost[is.bite])>ub)] # which current time since last boost is more than ub
    cohort$age_ib_boost[is.bite.boost] <-cohort$age_d[is.bite.boost]
    cohort$IB[is.bite.boost]<-cohort$IB[is.bite.boost] + 1
    
    is.latent<-is.bite[which(runif(length(is.bite))<cohort$b[is.bite])]
    sev_track[is.latent,j]<-0  ## bite that will lead to infection,without intervention.
    
    ### New infections emerging from the liver
    is.inf<-NULL
    if(j>dur_E/dt) {
      is.inf<-which(sev_track[,(j-dur_E/dt)]==0 & sev_track[,j]!=2)
    }
    if(j/dt>burnin_d) tot_inf_no_pmc<-tot_inf_no_pmc+length(is.inf)
    
      ### PMC effect in reducing infections. Who is still potentially protected - all.
    is.inf.on.pmc<-is.inf[which((j*dt - cohort$last_drug_time[is.inf]) < max_proph_time)]
    
    # track number of infectious bites among postdischarge group (not necessarily same as population)
    is.bite.postd.sma<-is.bite[which(is.bite %in% is.postd)]
    # among infected, who was discharged recently 
    is.inf.postd.sma<-is.inf[which(is.inf %in% is.postd)]
    
    is.protected<-NULL
    if(length(is.inf.on.pmc)>0) {
      # extract probabilities of protection for each person.
      curr_p_protect<-vector(length=length(is.inf.on.pmc))
      time_since_pmc<-j*dt - cohort$last_drug_time[is.inf.on.pmc]
      for(i in 1:length(is.inf.on.pmc)) {
        if(cohort$last_drug_type[is.inf.on.pmc[i]]=="al")
          curr_p_protect[i]<-p_protect1[time_since_pmc[i]/dt]
        if(cohort$last_drug_type[is.inf.on.pmc[i]]=="dp")
          curr_p_protect[i]<-p_protect2[time_since_pmc[i]/dt]
      }
      is.protected<-is.inf.on.pmc[which(runif(length(is.inf.on.pmc))<curr_p_protect)]
      # remove those who are protected
      is.inf<-is.inf[!is.inf %in% is.protected]
    }
    
    #### Count up bites, infections and protections during the post discharge follow up period only (start_measure to end_measure)
    if(j/dt>burnin_d) tot_bites_postd<-tot_bites_postd+length(is.bite.postd.sma)
    is.protected.postd<-NULL
    if(length(is.inf.postd.sma)>0) {
      if(j/dt>burnin_d) tot_inf_postd<-tot_inf_postd+length(is.inf.postd.sma)
      if(length(is.protected)>0) is.protected.postd<-is.protected[is.protected %in% is.inf.postd.sma]
      if(j/dt>burnin_d) tot_protect_start_end<-tot_protect_start_end+length(is.protected.postd)
    }
    ### Severe episodes
    # Now at the start of the severe episode, record hospital stay and PMC.
    # code 0=bite which will lead to sev if not stopped.
    # -1=nothing, 1=start of episode, 2=in hospital, 3=al at discharge, 4=AL at discharge. 5=home pmc course started at home.
    # Can start a new severe episode if severe bite durE days ago and not in hospital now.
    # TODO add drug effects.
    is.start.sev<-is.inf[which(runif(length(is.inf))<cohort$theta[is.inf])]
    if(j/dt>burnin_d) {
      tot_sev<-tot_sev + length(is.start.sev)
      rel_foi_sev<-c(rel_foi_sev,cohort$rel_foi[is.start.sev])
    }
    is.start.sma<-is.start.sev[which(runif(length(is.start.sev))<cohort$prob_sma[is.start.sev])]
    cohort$last_sma_time[is.start.sma]<-j*dt
    if(j/dt>burnin_d) {
      tot_sma<-tot_sma + length(is.start.sma)
      rel_foi_sma<-c(rel_foi_sma,cohort$rel_foi[is.start.sma])
    }
    if(save_age_inc & j/dt>burnin_d) {
      sev_inc_age<-c(sev_inc_age,cohort$age_d[is.start.sev])
    }
    
    is.sev.postd.sma<-is.start.sev[which(is.start.sev %in% is.postd)]
    
    if(length(is.sev.postd.sma)>0 & j/dt>burnin_d) tot_sev_postd<-tot_sev_postd+length(is.sev.postd.sma)
      
    curr_person_time_postd<-dt*length(is.postd)
    
    if(j/dt>burnin_d) {
      person_time_postd<-person_time_postd+  curr_person_time_postd
    }
    
    sev_track[is.start.sev,j]<-1
    cohort$last_sev_time[is.start.sev]<-j*dt
    
    ## severe immunity boosting.
    # normal
    if(!delay_sev_boost) {
      is.sev.boost<-is.inf[which((cohort$age_d[is.inf] -cohort$age_sev_boost[is.inf])>uv)]
    }
    # delayed (no boosting allowed within x months of having a severe episode)
    if(delay_sev_boost) {
      is.sev.boost<-is.inf[which((cohort$age_d[is.inf] -cohort$age_sev_boost[is.inf])>uv
                                 & (j*dt - cohort$last_sev_time[is.inf])>delay_boost_time)]
    }
    cohort$IVA[is.sev.boost]<-cohort$IVA[is.sev.boost] + 1
    cohort$age_sev_boost[is.sev.boost] <-cohort$age_d[is.sev.boost]
    
    ### Hospitalization
    # eligible to start hospital (no delay at the moment)
    is.start.hosp<-which(sev_track[,j-1]==1)
    ## gets to go to hospital
    is.start.hosp<-is.start.hosp[which(runif(length(is.start.hosp))<prob_hosp)]
    cohort$hosp_time[is.start.hosp]<-j*dt
    ## still in hospital?
    is.hosp<-which((j*dt - cohort$hosp_time)<hospital_stay)
    sev_track[is.hosp,j]<-2
    
    ## who was hospitalised 1 hospital stay ago.
    is.discharged<-which((j*dt-cohort$hosp_time)==hospital_stay)
    
    ## case fatality
    is.died.mal<-is.discharged.alive<-is.al.discharge<-NULL
    if(length(is.discharged)>0) {
      survive<-runif(length(is.discharged))>prob_sev_fatal
      is.discharged.alive<-is.discharged[which(survive)]
      sev_track[is.discharged.alive,j]<-3
      is.died.mal<-is.discharged[which(!survive)]
      ## reset characteristics for those who died:
      if(j*dt>start_pmc) is.al.discharge<-is.discharged.alive[which(runif(length(is.discharged.alive))<pmc_cov)]
      sev_track[is.al.discharge,j]<-4
      cohort$last_drug_time[is.al.discharge]<-j*dt
      cohort$last_drug_type[is.al.discharge]<-"al"
    } 
    
    ### died or aged out of age band (but keep in those on PMC to the end of their follow up).
    is.aged<- which(cohort$age_d>max_age_cohort & 
                      (j*dt - cohort$last_sev_time)>max_fup)
    is.died<-which(runif(n_cohort)<eta)
    is.removed<-c(is.died.mal,is.aged,is.died)
    is.removed<-is.removed[which(!duplicated(is.removed))]
    if(length(is.removed)>0) {
      sev_track[is.removed,j]<-9  # record removal.
      # reset characteristics for those who are removed.
      cohort$last_drug_time[is.removed]<- -1000
      cohort$last_sev_time[is.removed]<- -1000
      cohort$last_sma_time[is.removed]<- -1000
      cohort$hosp_time[is.removed]<- -1000
      cohort$last_drug_type[is.removed]<- NA
      cohort$age_d[is.removed]<-min_age_cohort  
      # reset immunity to what the person would have had at the age of inclusion into study
      # rel foi stays the same.
      cohort$IB[is.removed]<-IB_test5d_i[inds_indiv[is.removed],min_age_cohort/dt_5d]
      cohort$IVA[is.removed]<-IVA_test5d_i[inds_indiv[is.removed],min_age_cohort/dt_5d]
      cohort$IVM[is.removed]<-IVM_test5d_i[inds_indiv[is.removed],min_age_cohort/dt_5d]
      cohort$b[is.removed]<-b_i[inds_indiv[is.removed],min_age_cohort/dt_5d]
      cohort$theta[is.removed]<-theta_i[inds_indiv[is.removed],min_age_cohort/dt_5d]
      cohort$age_ib_boost[is.removed]<-time_ib_boost[inds_indiv[is.removed],min_age_cohort/dt_5d]
      cohort$age_sev_boost[is.removed]<-time_sev_boost[inds_indiv[is.removed],min_age_cohort/dt_5d]
      
    }
    
        ### Start PMC
    # eligible for PMC
    is.pmc<-NULL
    if(j*dt>start_pmc) {
      is.pmc<-which((j*dt - cohort$hosp_time) %in% pmc_times)
      # PMC coverage - actually get PMC
      is.pmc<-is.pmc[which(runif(length(is.pmc))<pmc_cov)]
      sev_track[is.pmc,j]<-5
      cohort$last_drug_time[is.pmc]<-j*dt
      cohort$last_drug_type[is.pmc]<-pmc_drug
    }
    
    # summ up per week
    if(save_summary_over_time) {
      week_ind<-which(res_over_time$weeks==ceiling((j*dt-1)/7))  ## minus 1 so that 1-7 are week 1, 8-14 are week 2 etc.
      res_over_time$n_sma[week_ind]<-res_over_time$n_sma[week_ind]+length(is.start.sma)
      res_over_time$n_reinf_pd[week_ind]<-res_over_time$n_reinf_pd[week_ind]+length(is.inf.postd.sma)
      res_over_time$persondays_pd[week_ind]<-res_over_time$persondays_pd[week_ind]+curr_person_time_postd
    }
  } ### end of pmc simulation
  
  if(!is.null(EIR_d_annual)) {
    mean_EIR0<-mean(EIR_d_annual)*365
  } else mean_EIR0<-NA
  mean_EIR_post_burnin_d<-365*mean(eir_d[(burnin_d/dt):length(eir_d)])
  res<-data.frame(adm=adm, tot_bites=tot_bites, tot_inf=tot_inf_no_pmc,tot_sev=tot_sev,tot_sma=tot_sma, 
                  inc_inf_ppy=tot_inf_no_pmc/(n_cohort*((length(eir_d)-burnin_d)/365)), 
                  inc_sev_ppy=tot_sev/(n_cohort*((length(eir_d)-burnin_d)/365)), 
                  inc_sma_ppy=tot_sma/(n_cohort*((length(eir_d)-burnin_d)/365)), 
                  mean_EIR=mean_EIR0,
                  mean_EIR_post_burnin_d=mean_EIR_post_burnin_d,
                  tot_inf_postd=tot_inf_postd,
                  tot_bites_postd,
                  tot_sev_postd=tot_sev_postd,
                  tot_protect_pmc=tot_protect_start_end,
                  person_time_postd=person_time_postd)
  
  saveRDS(res,paste0(output_root,"adm",adm))
  if(save_sev_track) {
    # get rid of the burnin_d period
    sev_track<-sev_track[(burnin_d/dt+1):ncol(sev_track)]
    saveRDS(sev_track,paste0(output_root,"sev_track",adm))
  }
  if(save_age_inc) {
    saveRDS(sev_inc_age,paste0(output_root,"sev_inc_age",adm))
    saveRDS(cohort$age_d,paste0(output_root,"ages_all",adm))
  }
  if(save_rel_foi) {
    saveRDS(cohort$rel_foi,paste0(output_root,"rel_foi_all",adm))
    saveRDS(rel_foi_sev,paste0(output_root,"rel_foi_sev",adm))
    saveRDS(rel_foi_sma,paste0(output_root,"rel_foi_sma",adm))
  }
  if(save_summary_over_time) {
    saveRDS(res_over_time,paste0(output_root,"res_over_time",adm))
  }
}
        
