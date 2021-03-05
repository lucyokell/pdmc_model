####################### RUN PMC ON THE CLUSTER ##################################################
### THIS IS A PROJECT FILE. NB if running other R code at same time, and have setwd to something else, need to change back to
# this project drive (Q, Lucy, PMC_post_discharge etc) for it to work.
setwd("Q:/Lucy/PMC_post_discharge")

#################################################################################################
# ### TODO LIST PMC MODEL CODE
#################################################################################################
# check each step.
# Faster burnin function for immunity from drug resistance code.
# accurate age structure of under fives.
# track last severe episode in the burn in part of the code, since kids join studies at age>0
# Add Amani's mechanism.


##################################################################################################
##  SET UP  
##################################################################################################

#install.packages("drat")
drat:::add("mrc-ide")
#install.packages("didehpc")
#install.packages(c("storr", "provisionr", "context", "queuer", "didehpc"))
#options(didehpc.username="lokell")

library(didehpc)
library(storr)
library(provisionr)
library(context)
library(queuer)

#options(didehpc.cluster = "fi--didemrchnb")
options(didehpc.cluster = "fi--dideclusthn")
didehpc::web_login()

packages<-c("dplyr")
## SOURCES MUST BE FILES WITH FUNCTIONS IN
sources<-c("run_model.R")

#config <- didehpc::didehpc_config(rtools = TRUE, home = "Q:/", cluster = "fi--didemrchnb") 
config <- didehpc::didehpc_config(rtools = TRUE, home = "Q:/", cluster = "fi--dideclusthn") 
#config <- didehpc::didehpc_config(rtools = TRUE, home = "Q:/", cluster = "fi--didemrchnb", template = "24Core", cores = 24) 
#config <- didehpc::didehpc_config(home = "Q:/", cluster = "fi--didemrchnb", template = "16Core", cores = 16) 
#config <- didehpc::didehpc_config(home = "Q:/", cluster = "fi--didemrchnb", template = "20Core", cores = 20) 
#config <- didehpc::didehpc_config(home = "Q:/", cluster = "fi--didemrchnb", template = "24Core", cores = 8) 
#config <- didehpc::didehpc_config(home = "Q:/", cluster = "fi--didemrchnb", template = "32Core", cores = 32) 

ctx <- context::context_save("contexts", packages = packages, sources = sources)
obj <- didehpc::queue_didehpc(ctx,config=config)


t <- obj$enqueue(sessionInfo())

###########################################################################################
### READ IN DATA
###########################################################################################

# EIR - either set constant or read in seasonal patterns from admin units
EIR_test<-100/365 # constant

eir_adms<-read.csv("adm1_eir_annual.csv")
eir_var_inds<-grep("X",names(eir_adms))
# plot(1:length(EIR_test),EIR_test)


#####################################################################################
# Check the model is working locally and on cluster with test runs.
#####################################################################################

########### Run model once locally for a test
adm0<-2
source("run_model.R")
EIR_d_annual0<-as.numeric(eir_adms[adm0,eir_var_inds])/365
run_pmc_model(adm=adm0,N=1000,n_cohort=1000,n_years=1,burnin=0,EIR_d_annual=EIR_d_annual0,output_root = "trial_output/",save_sev_track=T)

readRDS(paste0("trial_output/adm",adm0))
x<-readRDS(paste0("trial_output/sev_track",adm0))
x
########### send a test job to the cluster
set.seed(1)

adm0<-3
EIR_d_annual0<-as.numeric(eir_adms[adm0,eir_var_inds])/365
t1<-obj$enqueue(run_pmc_model(adm=adm0,N=1000,n_cohort=1000,n_years=1,EIR_d_annual = EIR_d_annual0)) 
t1$status()
t1$log()
readRDS("output/adm3")


#####################################################################################
# Set up multiple cluster runs
#####################################################################################

##### RUN PMC MODEL FOR ALL ADMINS
for(i in 1:nrow(eir_adms)) {
#for(i in 401:600) {
  adm0<-eir_adms$DIDE_CODE[i]
  print(adm0)
  EIR_d_annual0<-as.numeric(eir_adms[i,eir_var_inds])/365  # run using the most recent EIR.
  # model is run for n_years, but outputs only recorded after burnin
  t1<-obj$enqueue(run_pmc_model(adm=adm0,N=1000,n_cohort=50000,EIR_d_annual = EIR_d_annual0,output_root="output/no_pmc_",
                                n_years=10,burnin=720,pmc_cov = 0, pmc_drug="dp", pmc_weeks=c(2,6,10),start_measure= 3*7, # time start measuring outcome, since discharge.
                                end_measure=14*7,min_age_cohort = 4*30.5,max_age_cohort = 5*365))  ## pd outcomes from start measure to end measure. 
}
t1$status()
t1$log()

readRDS("output/no_pmc_adm100041")

#####################################################################################
# Combine results
#####################################################################################
##### Somewhere to output results.
reinfection_risk<-dplyr::select(eir_adms,DIDE_CODE:NAME_1)
reinfection_risk<-reinfection_risk %>%
  dplyr::mutate(tot_sev=NA,
                tot_sma=NA,
                inc_inf_ppy=NA,
                inc_sev_ppy=NA,
                inc_sma_ppy=NA,
                EIR=NA,
                tot_inf_postd=NA,
                tot_bites_postd=NA,
                tot_sev_postd=NA,
                person_time_postd=NA  )

for(i in 1:nrow(eir_adms)) {
  #x<-readRDS(paste0("output/adm",i))
  adm<-eir_adms$DIDE_CODE[i]
  x<-readRDS(paste0("output/no_pmc_adm",adm))
  reinfection_risk$tot_sev[i]<-x$tot_sev
  reinfection_risk$tot_sma[i]<-x$tot_sma
  if(length(x$inc_inf_ppy)>0) reinfection_risk$inc_inf_ppy[i]<-x$inc_inf_ppy
  reinfection_risk$inc_sev_ppy[i]<-x$inc_sev_ppy
  reinfection_risk$inc_sma_ppy[i]<-x$inc_sma_ppy
  reinfection_risk$EIR[i]<-x$mean_EIR
  reinfection_risk$tot_inf_postd[i]<-x$tot_inf_postd
  reinfection_risk$tot_bites_postd[i]<-x$tot_bites_postd
  reinfection_risk$person_time_postd[i]<-x$person_time_postd
}

reinfection_risk<-reinfection_risk %>%
  dplyr::mutate(inc_bites_postd_ppy=tot_bites_postd/(person_time_postd/365))

# reinfection_risk<-reinfection_risk %>%
#   mutate(risk_reinf=any_inf_recent_sev3/tot_recent_sev3,
#          inc_sev_recent_sev3=tot_sev_recent_sev3/(person_time_sev3/(30.5*3)),  ## incidence during follow up.
#          inc_inf_recent_sev3=tot_inf_recent_sev3/(person_time_sev3/(30.5*3)),
#          inc_reinf_pd_100py=100*tot_inf_postd_start_end/(person_time_postd_start_end/365),
#          inf_avert_pd_100py=100*tot_protect_pmc/(person_time_postd_start_end/365),
#          pmc_eff=tot_protect_pmc/tot_inf_postd_start_end
#   )  

#write.csv(reinfection_risk,"reinfection_risk.csv")
#write.csv(reinfection_risk,"reinfection_risk_10y_50000.csv")
write.csv(reinfection_risk,"reinfection_risk_10y_50000_2.csv")




##### RUN PMC MODEL FOR ALL ADMINS. WITH PMC TURNED ON.
#### ASSUME CONSTANT RISK OF SEVERE MALARIA FOR FIRST 3 MONTHS.
for(i in 1:nrow(eir_adms)) {
  #for(i in 562:nrow(eir_adms)) {
  adm0<-eir_adms$DIDE_CODE[i]
  print(adm0)
  EIR_d_annual0<-as.numeric(eir_adms[i,eir_var_inds])/365  # run using the most recent EIR.
  # model is run for n_years, but outputs only recorded after burnin
  t1<-obj$enqueue(run_pmc_model(adm=adm0,N=1000,n_cohort=50000,EIR_d_annual = EIR_d_annual0,output_root="output/pmc_",
                                n_years=10,burnin=720,pmc_cov = 1, pmc_drug="dp", pmc_weeks=c(2,6,10),start_measure= 3*7, # time start measuring outcome, since discharge.
                                end_measure=14*7,min_age_cohort = 4*30.5,max_age_cohort = 5*365))  ## pd outcomes from start measure to end measure. 
}
t1$status()
t1$log()

readRDS("output/no_pmc_adm100035")


