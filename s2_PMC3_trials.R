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
library(dplyr)

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


#####################################################################################
# Trial validation params
#####################################################################################

# Trial most recent
hospital_stay<-4
pmc_drug<-"dp"
pmc_times<-c(2,6,10)*7 + hospital_stay  # time since admission.
start_measure<- 3*7 # time start measuring outcome, since discharge.
end_measure<-14*7

# Trial Phiri et al 2012
pmc_drug<-"al"
pmc_times<-c(4,8)*7 + hospital_stay  # time since admission.
start_measure<- 4*7 # time start measuring outcome, since discharge.
end_measure<-12*7


#####################################################################################
# Check the model is working locally and on cluster with test runs.
#####################################################################################

########### Run model once locally for a test
# pmc_cov<-0
# fixed_theta<-TRUE
# EIR_test<-5/365
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


###################################################
# RUN TRIAL SITES
###################################################

####################################################
# USING ADMIN LEVEL EIR
########### Malawi Phiri 2012 ######################
# EIRS MADE IN s1_output_EIR_epi_trials.R
eir_phiri<-read.csv("C:/Users/lokell/Dropbox (SPH Imperial College)/PMC post discharge/res_phiri2012_eir.csv",header=T)
eir_var_inds_phiri<-grep("X",names(eir_phiri))
# dates = June 2006 - August 2009 for RECRUITMENT
# dates in model launcher years:
start_trial<-2006+5/12 - 2017
end_trial<-2009+8/12 +3/12 - 2017

## mean EIR (not population weighted):
trial_eir<-eir_phiri %>% 
  dplyr::select(X1:X3620)
year<- ((-17*365):round(end_trial*365))/365
inds<-which(year>start_trial)
inds<-inds[inds<=ncol(trial_eir)]  ## slight fudge - should work, check this. only a few days out (5)
mean(unlist(trial_eir[,inds]))
trial_eir<-as.matrix(trial_eir[,inds])
for(i in 1:nrow(eir_phiri)) print(mean(trial_eir[i,]),na.rm=T)

#######
# PLOT EIR
par(mar=c(4,4,5,2))
plot((1:length(eir_var_inds_phiri))/365+2000,eir_phiri[1,eir_var_inds_phiri],type="l",ylab="EIR annual",xlab="year",
     ylim=c(0,max(eir_phiri[,eir_var_inds_phiri])))
for(i in 2:nrow(eir_phiri)) lines((1:length(eir_var_inds_phiri))/365+2000,eir_phiri[i,eir_var_inds_phiri],col=i)
lines(rep(start_trial+2017,2),c(0,1000),lty=2)
lines(rep(end_trial+2017,2),c(0,1000),lty=2)
legend(2002,205,c(eir_phiri$NAME_1),col=c("black",2:nrow(eir_phiri)),lty=1,xpd=T,bty='n')

malawi_admins<-c("Blantyre","Chikwawa","Thyolo","Zomba")
malawi_inds<-which(eir_adms$NAME_0=="Malawi" & eir_adms$NAME_1 %in% malawi_admins)
malawi_dide_codes<-eir_adms$DIDE_CODE[malawi_inds]

for(i in 1:length(malawi_inds)) {
  adm0<-malawi_inds[i]
  EIR_d_annual0<-as.numeric(eir_adms[adm0,eir_var_inds])/365
  t1<-obj$enqueue(run_pmc_model(adm=adm0,N=1000,n_cohort=50000,EIR_d_annual = EIR_d_annual0,output_root = "trial_output/",
                                save_sev_track = T)) 
}
t1$status()
t1$log()


########### Kwambai 2020 ######################
eir_kwambai<-read.csv("C:/Users/lokell/Dropbox (SPH Imperial College)/PMC post discharge/res_kwambai_eir.csv",header=T)
eir_var_inds_kwambai<-grep("X",names(eir_kwambai))
# dates = May 2016 - May 2018 for RECRUITMENT
# dates in model launcher years:
start_trial<-2016+4/12 - 2017
end_trial<-2018+5/12 +3/12 - 2017

## mean EIR (not population weighted):
trial_eir<-eir_kwambai %>% 
  dplyr::select(X1:X6566) %>%
  dplyr::filter(!is.na(X1))
year<- ((-17*365):round(end_trial*365))/365
inds<-which(year>start_trial & year<end_trial)
inds<-inds[which(inds<=ncol(trial_eir))]  ## slight fudge - should work, check this. only a few days out (5)
trial_eir<-as.matrix(trial_eir[,inds])
mean(trial_eir,na.rm=T)
eir_kwambai<-dplyr::filter(eir_kwambai,!is.na(X1))
for(i in 1:nrow(eir_kwambai)) print(mean(trial_eir[i,]),na.rm=T)
trial_eir_kwambai<-dplyr::select(eir_kwambai,DIDE_CODE:NAME_1)
trial_eir_kwambai$trial_eir<-NA
for(i in 1:nrow(trial_eir_kwambai)) trial_eir_kwambai$trial_eir[i]<-mean(trial_eir[i,],na.rm=T)
write.csv(trial_eir_kwambai, file="C:/Users/lokell/Dropbox (SPH Imperial College)/PMC post discharge/trial_eir_kwambai.csv",row.names = F)


#######
# PLOT EIR
par(mar=c(4,4,5,2))
plot((1:length(eir_var_inds_kwambai))/365+2000,eir_kwambai[1,eir_var_inds_kwambai],type="l",ylab="EIR annual",xlab="year",
     ylim=c(0,max(eir_kwambai[,eir_var_inds_kwambai],na.rm=T)))
for(i in 2:nrow(eir_kwambai)) lines((1:length(eir_var_inds_kwambai))/365+2000,eir_kwambai[i,eir_var_inds_kwambai],col=i)
lines(rep(start_trial+2017,2),c(0,1000),lty=2)
lines(rep(end_trial+2017,2),c(0,1000),lty=2)
legend(2002,205,c(eir_kwambai$NAME_1),col=c("black",2:nrow(eir_kwambai)),lty=1,xpd=T,bty='n')


for(i in 1:nrow(eir_adms)) {
  #for(i in 1:50) {
  adm0<-i
  EIR_d_annual0<-as.numeric(eir_adms[adm0,eir_var_inds])/365
  t1<-obj$enqueue(run_pmc_model(adm=adm0,N=1000,n_cohort=50000,EIR_d_annual = EIR_d_annual0)) 
}
t1$status()
t1$log()



######################################################
# SITE SPECIFIC TRIAL PREDICTIONS USING MAP PREV 20KM AROUND THE HOSPITALS
######################################################

########### All hospitals ######################
# EIRS MADE IN output_EIR.R
eir_trial<-read.csv("C:/Users/lokell/Dropbox (SPH Imperial College)/PMC post discharge/res_trial_eir2.csv",header=T)
eir_var_inds_trial<-grep("X",names(eir_trial))
### trial characteristics.
hosp<-read.csv("C:/Users/lokell/Dropbox (SPH Imperial College)/PMC post discharge/hospital_locations.csv")
hosp<-filter(hosp,!(fup_end_month>4 & study=="Kwambai"))

#######
# PLOT EIR
par(mar=c(4,4,5,2))
plot((1:length(eir_var_inds_trial))/365+2000,eir_trial[1,eir_var_inds_trial],type="l",ylab="EIR annual",xlab="year",
     ylim=c(0,max(eir_trial[,eir_var_inds_trial],na.rm=T)))
for(i in 2:nrow(eir_trial)) lines((1:length(eir_var_inds_trial))/365+2000,eir_trial[i,eir_var_inds_trial],col=i)

#####################
# RUN MODEL FOR EACH TRIAL SITE. 
# Run from 2000 to allow burnin for maternal immunity.
# NB initial immunity set up in the model automatically uses the first year of the EIR and runs a cohort up to age 20.
# age specific and follow up specific.
for(i in 1:nrow(eir_trial)) {
  #for(i in 1:14) {
  adm0<-paste0(eir_trial$study[i],"_",eir_trial$DIDE_CODE[i])
  print(adm0)
  EIR_d_2000_2018<-as.numeric(eir_trial[i,eir_var_inds_trial])/365  ## but never use 2018 as do not have intervention coverage yet.
  days_start_eir<-max(hosp$start_year_trial[i]-5-2000, 0)*365  ## start max of 5 years before trial to run in immunity for 5 yr olds.
  days_end_eir<-(min(hosp$end_year_trial[i],2017)-2000)*365  
  EIR_d_trial<-EIR_d_2000_2018[days_start_eir:days_end_eir]   ### only run to the end of the trial period.
  start_fup<-hosp$fup_start_month[i]*30.5
  end_fup<-hosp$fup_end_month[i]*30.5
  days_start_trial<-max(hosp$start_year_trial[i]-2000, 0)*365 - days_start_eir
  #test<-EIR_d_trial[days_start_trial:length(EIR_d_trial)]
  
  age_min<-hosp$age_low[i]*365 + 1  ## add one day so as to deal with studies which have age=0 (cannot have index=0)
  age_max<-hosp$age_high[i]*365
  t2<-obj$enqueue(run_pmc_model(adm=adm0,N=1000,
                                n_cohort=100000,
                                EIR_fully_specified = EIR_d_trial,
                                pmc_cov = 0, 
                                start_measure= start_fup, # time start measuring outcome, since discharge.
                                end_measure=end_fup,
                                burnin_d=days_start_trial,  ### RELATIVE TO BEGINNING OF EIR_D_TRIAL
                                min_age_cohort=age_min,
                                max_age_cohort=age_max,
                                save_summary_over_time = T,
                                output_root="trial_output/"))
  
}
t2$status()
t2$log()

y<-readRDS("trial_output/admOpoka_2016_103262")
x<-readRDS("trial_output/res_over_timeOpoka_2016_103262")
plot(x$weeks,x$n_sma/(x$persondays_all/365),type="l")
plot(x$weeks,x$n_reinf_pd/(x$persondays_pd/365),type="l")
## check output is correct.
temp<-x %>% dplyr::filter(weeks>=days_start_trial/7)
sum(temp$n_reinf_pd)
sum(temp$persondays_pd)
sum(temp$n_reinf_pd)/ (sum(temp$persondays_pd)/365)
y$tot_inf_postd/(y$person_time_postd/365)
y$tot_bites_postd/(y$person_time_postd/365)


#### add modelled reinfections to trial data
trial_hospital_epi<-read.csv("C:/Users/lokell/Dropbox (SPH Imperial College)/PMC post discharge/trial_hospital_epi.csv")
trial_hospital_epi<-trial_hospital_epi %>%
  dplyr::mutate(eir_pmc_model=NA,
                n_bites_pd=NA,
                n_reinf_pd=NA,
                person_days_postd=NA,
                inc_reinf_pd_ppy=NA,
                inc_bites_pd_ppy=NA)
for(i in 1:nrow(trial_hospital_epi)) {
  adm0<-paste0(trial_hospital_epi$study[i],"_",trial_hospital_epi$DIDE_CODE[i])
  y<-readRDS(paste0("trial_output/adm",adm0))
  trial_hospital_epi$eir_pmc_model[i]<-y$mean_EIR_post_burnin_d
  trial_hospital_epi$n_bites_pd[i]<-y$tot_bites_postd
  trial_hospital_epi$n_reinf_pd[i]<-y$tot_inf_postd
  trial_hospital_epi$person_days_postd[i]<-y$person_time_postd
  trial_hospital_epi$inc_bites_pd_ppy[i]<-y$tot_bites_postd/(y$person_time_postd/365)
  trial_hospital_epi$inc_reinf_pd_ppy[i]<-y$tot_inf_postd/(y$person_time_postd/365)
}

write.csv(trial_hospital_epi,"C:/Users/lokell/Dropbox (SPH Imperial College)/PMC post discharge/trial_hospital_epi.csv")


############################
# get age specific sev inc for certain fixed EIRs
####### fixed eirs.

fixed_eirs<-c(100,10,1)
for(i in 1:length(fixed_eirs)) {
  EIR_d_annual0<-rep(fixed_eirs[i]/365,365)
  n_years0<-50
  burnin0<-365
  t1<-obj$enqueue(run_pmc_model(adm=paste0("eir",fixed_eirs[i]),N=1000,n_cohort=50000,n_years=n_years0,burnin=burnin0,
                                EIR_d_annual=EIR_d_annual0,output_root = "output/",save_age_inc=T,save_rel_foi=T)) 
}
t1$status()
t1$log()
readRDS(paste0("output/admeir",fixed_eirs[i]))

for(i in 1:length(fixed_eirs)) {
  sev_age<-readRDS(paste0("output/sev_inc_ageeir",fixed_eirs[i])) / 365
  all_ages<-readRDS(paste0("output/ages_alleir",fixed_eirs[i]))/ 365
  brks<-seq(0,5,by=1/12)
  
  sev_age_summ<-cut(sev_age,breaks = brks,include.lowest = T)
  sev_age_summ<-as.numeric(table(sev_age_summ))
  py_pop<-cut(all_ages,breaks = brks,include.lowest = T)
  py_pop<-as.numeric(table(py_pop))
  py_pop<-py_pop*(n_years0 - burnin0/365)
  sev_inc<-sev_age_summ/py_pop
  sev_inc_age<-data.frame(min_age=brks[1:(length(brks)-1)],max_age=brks[2:length(brks)],sev_inc=sev_inc)
  saveRDS(sev_inc_age,paste0("output/sev_inc_age_eir",fixed_eirs[i]))
}


