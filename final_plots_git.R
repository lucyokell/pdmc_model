#### PLOTS

library(dplyr)
library(gplots)
library(vioplot)
library(rstan)  ## read in model fits
library(epitools)
library(rgdal)
library(readxl)
source("color.bar.R")
source("reduce_shape_file.R")

########
# DATA / MODEL OUTPUTS
#########

dat3d<-readRDS("daily_kwambai_3D.RDS") # not publicly available
# read in model fit
fit<-readRDS("pmc14g_cut8_fit.rds") ## too large to upload to github (run smaller posterior sample?)
params<-rstan::extract(fit,pars=c("EIR_y","prob_um_plac","prob_um_pmc","max_beta_phi","r",
                                  "w_scale_risk","w_shape_risk",
                                  "res"))
params_for_mod<-rstan::extract(fit,pars=c("EIR_y","prob_um_plac","prob_um_pmc","max_beta_phi","r",
                                          "w_scale_risk","w_shape_risk"))
saveRDS(params_for_mod,file="model_fit_stan_output_paramsonly_14g.rds")
rm(fit)


f2<-readRDS("pmc14g_cut8_fit_gqs.rds") ## generated quantities from stan. 
params_gqs<-rstan::extract(f2,pars=c("res2","events_tot"))
saveRDS(params_gqs,file="model_fit_stan_gqs.rds") # too large for github


rm(f2) # to save on memory
imperial_mod<-read.csv("adm1_epi_2019_sma2_paton_cal_repmod.csv")
shp1 <- readOGR(dsn="shape_files", layer="SSA_shapefile_adm1")
shp0 <- readOGR(dsn="shape_files", layer="SSA_shapefile_adm0")
shp2 <- readOGR(dsn="shape_files", layer="Africa_Boundaries")
pop<-read.csv("population_admin.csv")
imperial_mod<-left_join(imperial_mod,pop,by=c("ISO","NAME_1"))


#############################
# Parameter summary - post discharge stan model
#############################
runTime<-182-7  ## 6 months minus 7 days to cut off last week.
params_for_mod<-readRDS("model_fit_stan_output_paramsonly_14g.rds")
round(quantile(params_for_mod$max_beta_phi,probs = c(0.5,0.025,0.975)),digits=2)
round(quantile(params_for_mod$r,probs = c(0.5,0.025,0.975)),digits=3)
round(quantile(params_for_mod$w_shape_risk,probs = c(0.5,0.025,0.975)),digits=2)
round(quantile(params_for_mod$w_scale_risk,probs = c(0.5,0.025,0.975)),digits=1)
round(quantile(1-params_for_mod$prob_um_plac,probs = c(0.5,0.025,0.975)),digits=2)
round(quantile(1-params_for_mod$prob_um_pmc,probs = c(0.5,0.025,0.975)),digits=2)

## output as EIR, LCI, UCI, rounded
for(i in 1:9) print(
round(quantile(params_for_mod$EIR_y[,i],probs = c(0.5,0.025,0.975)),digits=1))


#################
# Process data
#################
## trial data:
trial_epi<-read.csv("trial_hospital_epi.csv")
trial_epi<-dplyr::filter(trial_epi,study %in% c("Kwambai", "Opoka_2016","Opoka_2019"))
trial_epi<-trial_epi %>%
  dplyr::mutate(plac_sm_lci3=100*(pois.exact(placebo_sm_events3,pt=placebo_sm_py3)$lower),
                plac_sm_uci3=100*(pois.exact(placebo_sm_events3,pt=placebo_sm_py3)$upper),
                plac_um_lci3=ifelse(!is.na(placebo_um_events3), 100*(pois.exact(placebo_um_events3,pt=placebo_um_py3)$lower),NA),
                plac_um_uci3=100*(pois.exact(placebo_um_events3,pt=placebo_um_py3)$upper),
                plac_sm_inc_p100py16=100*(placebo_sm_events3+placebo_sm_events6)/(placebo_sm_py3+placebo_sm_py6),
                plac_um_inc_p100py16=100*(placebo_um_events3+placebo_um_events6)/(placebo_um_py3+placebo_um_py6),
                pmc_sm_inc_p100py16=100*(pmc_sm_events3+pmc_sm_events6)/(pmc_sm_py3+pmc_sm_py6),
                pmc_um_inc_p100py16=100*(pmc_um_events3+pmc_um_events6)/(pmc_um_py3+pmc_um_py6),
                plac_um_inc_p100py3=100*placebo_um_events3/placebo_um_py3,
                plac_sm_inc_p100py3=100*placebo_sm_events3/placebo_sm_py3,
                plac_sm_inc_p100py6=100*placebo_sm_events6/placebo_sm_py6,
                pmc_um_inc_p100py3=100*pmc_um_events3/pmc_um_py3,
                sm_lci16=100*(pois.exact((placebo_sm_events3+placebo_sm_events6),pt=(placebo_sm_py3+placebo_sm_py6))$lower),
                sm_uci16=100*(pois.exact((placebo_sm_events3+placebo_sm_events6),pt=(placebo_sm_py3+placebo_sm_py6))$upper),
                tot_episodes_plac_um_sm=placebo_um_events3+placebo_um_events6+placebo_sm_events3+placebo_sm_events6,
                tot_episodes_pmc_um_sm=pmc_um_events3+pmc_um_events6+pmc_sm_events3+pmc_sm_events6,
                tot_inc_plac_um_sm_p100py=100*tot_episodes_plac_um_sm/(placebo_um_py3+placebo_um_py6),
                tot_inc_pmc_um_sm_p100py=100*tot_episodes_pmc_um_sm/(pmc_um_py3+pmc_um_py6),
                tot_inc_plac_um_sm_lci=100*(pois.exact((tot_episodes_plac_um_sm),pt=(placebo_sm_py3+placebo_sm_py6))$lower),
                tot_inc_plac_um_sm_uci=100*(pois.exact((tot_episodes_plac_um_sm),pt=(placebo_sm_py3+placebo_sm_py6))$upper),
                
  )  

### CI for pmc um events.
inds<-which(!is.na(trial_epi$pmc_um_events3))
trial_epi$pmc_um_lci3<-trial_epi$pmc_um_uci3<-NA
trial_epi$pmc_um_lci3[inds]<-100*(pois.exact(trial_epi$pmc_um_events3[inds],pt=trial_epi$pmc_um_py3[inds])$lower)
trial_epi$pmc_um_uci3[inds]<-100*(pois.exact(trial_epi$pmc_um_events3[inds],pt=trial_epi$pmc_um_py3[inds])$upper)


##########################
## Generate drug protection curves for plotting
##### AL DRUG 1. 
max_proph_time<-50
w_scale1<-93.5*0.139
w_slope1<-93.5
dt<-0.01
time<-seq(from=0,to=max_proph_time,by=dt)
p_protect1<- exp(-((time/w_scale1)^w_slope1))

##### DHA PIP DRUG 2
w_scale2<-28.1
w_slope2<-4.4
time<-seq(from=0,to=max_proph_time,by=dt)
p_protect2<- exp(-((time/w_scale2)^w_slope2))

#treat weeks 2,6,10 after discharge, as well as AL at time of 1st episode
pmc_times<-c(2,6,10)*7
PE<-data.frame(time=seq(from=0,to=26*7,by=dt),pp_pmc=0,pp_plac=0,pp_pmc_ad=0)
## allow for adherence.
PE$pp_plac[1:length(p_protect1)]<-100*p_protect1
PE$pp_pmc[1:length(p_protect1)]<-100*p_protect1
PE$pp_pmc_ad[1:length(p_protect1)]<-100*p_protect1
pmc_ad<-c(0.765,0.879,0.894)  #adherence

for(i in 1:length(pmc_times)) {
  start<-which(PE$time==pmc_times[i])
  inds<-start:(start+length(p_protect2)-1)
  inds<-inds[which(inds<=nrow(PE))]
  PE$pp_pmc[inds]<-100*p_protect2[1:length(inds)]
  PE$pp_pmc_ad[inds]<-100*p_protect2[1:length(inds)]*pmc_ad[i]
}

##########################
## Process shape files
shp2<-getSmallPolys(shp2,minarea = 10)  ## remove small islands
shp2<-subset(shp2, NAME_0!="Mauritius")
shp2<-subset(shp2, NAME_0!="Seychelles")
shp2<-subset(shp2, NAME_0!="Sao Tome and Principe")
shp2<-subset(shp2, NAME_0!="Reunion")
shp2<-subset(shp2, NAME_0!="Mayotte")
shp2<-subset(shp2, NAME_0!="French Southern Territories")
shp2<-subset(shp2, NAME_0!="Cape Verde")


africa<-subset(shp1,CONTINENT=="Africa")
africa@data<-africa@data %>%
  mutate(NAME_0=ifelse(NAME_0=="Côte d'Ivoire","Cote d'Ivoire",NAME_0))
africa<-subset(africa,NAME_0!="Mauritius")
africa<-subset(africa,NAME_0!="Seychelles")


## Process model fit to get model predictions and 95% CI
### process stan output
# res2 is 4D number of cases over the time period, dimensions:
# samples, days, var, site(EIR)
# was generated per 50 children.

EIRpred<-exp(log(seq(0.1,30,length.out = 100))) # the EIRs used in generated quantities in stan

fit_line<-data.frame(eir=EIRpred,
                     pred_pmc_um_inc3=NA,
                     pred_pmc_um_inc6=NA,
                     pred_pmc_sm_inc3=NA,
                     pred_pmc_sm_inc6=NA,
                     pred_plac_um_inc3=NA,
                     pred_plac_um_inc6=NA,
                     pred_plac_sm_inc3=NA,
                     pred_plac_sm_inc6=NA,
                     pred_plac_um_inc3_lci=NA,
                     pred_plac_um_inc3_uci=NA,
                     pred_plac_sm_inc3_lci=NA,
                     pred_plac_sm_inc3_uci=NA
)
                     
# sum the incidence over each day(rowSums)
# then take the mean over the samples.
# *100 -> per 100 per day
# person days per person = length(15:104) or 105:
# person years for each person = length(15:runTime)/365
# for each EIR level (100 EIR values used for predictions)
for(i in 1:100) {
  
  fit_line$pred_plac_sm_inc3[i]<-100*mean(rowSums(params_gqs$res2[,15:104,1,i]))/(100*length(15:104)/365)  
  fit_line$pred_plac_sm_inc6[i]<-100*mean(rowSums(params_gqs$res2[,105:runTime,1,i]))/(100*length(105:runTime)/365)  
  fit_line$pred_plac_um_inc3[i]<-100*mean(rowSums(params_gqs$res2[,15:104,3,i]))/(100*length(15:104)/365)  
  fit_line$pred_plac_um_inc6[i]<-100*mean(rowSums(params_gqs$res2[,105:runTime,3,i]))/(100*length(105:runTime)/365)  

  fit_line$pred_pmc_sm_inc3[i]<-100*mean(rowSums(params_gqs$res2[,15:104,2,i]))/(100*length(15:104)/365)  
  fit_line$pred_pmc_sm_inc6[i]<-100*mean(rowSums(params_gqs$res2[,105:runTime,2,i]))/(100*length(105:runTime)/365)  
  fit_line$pred_pmc_um_inc3[i]<-100*mean(rowSums(params_gqs$res2[,15:104,4,i]))/(100*length(15:104)/365)  
  fit_line$pred_pmc_um_inc6[i]<-100*mean(rowSums(params_gqs$res2[,105:runTime,4,i]))/(100*length(105:runTime)/365)  
  fit_line$pred_plac_um_inc3_lci[i]=100*quantile(rowSums(params_gqs$res2[,15:104,3,i]),probs=0.025)/(100*length(15:104)/365)
  fit_line$pred_plac_um_inc3_uci[i]=100*quantile(rowSums(params_gqs$res2[,15:104,3,i]),probs=0.975)/(100*length(15:104)/365)
  fit_line$pred_plac_sm_inc3_lci[i]=100*quantile(rowSums(params_gqs$res2[,15:104,1,i]),probs=0.025)/(100*length(15:104)/365)
  fit_line$pred_plac_sm_inc3_uci[i]=100*quantile(rowSums(params_gqs$res2[,15:104,1,i]),probs=0.975)/(100*length(15:104)/365)
  
}


# add a row of zeros to cover zero EIR settings
fit_line<-rbind(rep(0,ncol(fit_line)),fit_line)

if(write2file) saveRDS(fit_line,"fit_line_stan.RDS")

######## Daily fit
###################
#### Get data back to 2D
dat2d<-dat3d[,,1]
for(i in 2:9) dat2d<-dat2d+dat3d[,,i]
### remove last 7 days of data:
dat2d<-dat2d[1:(168-7),]
# order: SM plac, SM PMC, UM plac, UM PMC
fit_daily<-data.frame(days=1:nrow(dat2d)+14, 
                      dat_plac_sm=cumsum(dat2d[,1]),
                      dat_pmc_sm=cumsum(dat2d[,2]),
                      dat_plac_um=cumsum(dat2d[,3]),
                      dat_pmc_um=cumsum(dat2d[,4]),
                      # sum over sites, average over samples.
                      # samples, days, outcome,site
                      model_plac_sm=cumsum(rowSums(apply(params$res[,15:runTime,1,],c(2,3),mean))),
                      model_pmc_sm=cumsum(rowSums(apply(params$res[,15:runTime,2,],c(2,3),mean))),
                      model_plac_um=cumsum(rowSums(apply(params$res[,15:runTime,3,],c(2,3),mean))),
                      model_pmc_um=cumsum(rowSums(apply(params$res[,15:runTime,4,],c(2,3),mean))),
                      model_plac_sm_lci=apply(apply(apply(params$res[,15:runTime,1,],c(1,2),sum),1,cumsum), 1,quantile,probs=0.025),
                      model_plac_sm_uci=apply(apply(apply(params$res[,15:runTime,1,],c(1,2),sum),1,cumsum), 1,quantile,probs=0.975),
                      model_pmc_sm_lci=apply(apply(apply(params$res[,15:runTime,2,],c(1,2),sum),1,cumsum), 1,quantile,probs=0.025),
                      model_pmc_sm_uci=apply(apply(apply(params$res[,15:runTime,2,],c(1,2),sum),1,cumsum), 1,quantile,probs=0.975),
                      model_plac_um_lci=apply(apply(apply(params$res[,15:runTime,3,],c(1,2),sum),1,cumsum), 1,quantile,probs=0.025),
                      model_plac_um_uci=apply(apply(apply(params$res[,15:runTime,3,],c(1,2),sum),1,cumsum), 1,quantile,probs=0.975),
                      model_pmc_um_lci=apply(apply(apply(params$res[,15:runTime,4,],c(1,2),sum),1,cumsum), 1,quantile,probs=0.025),
                      model_pmc_um_uci=apply(apply(apply(params$res[,15:runTime,4,],c(1,2),sum),1,cumsum), 1,quantile,probs=0.975)
)



# Total fit numbers
tot_fits_daily<-data.frame(trial_arm=rep(c("pmc","plac"),each=4),
                           var=rep(c("tot_um3","tot_um6","tot_sm3","tot_sm6"),2),
                           dat=NA,
                           model=0,
                           model_lci=NA,
                           model_uci=NA)
### add up total events of each type across sites and days
# order in dat3d: days x SM plac, SM PMC, UM plac, UM PMC x site
tot_fits_daily$dat[1]<-sum(dat3d[1:90,4,])
tot_fits_daily$dat[2]<-sum(dat3d[91:(runTime-14),4,])
tot_fits_daily$dat[3]<-sum(dat3d[1:90,2,])
tot_fits_daily$dat[4]<-sum(dat3d[91:(runTime-14),2,])
tot_fits_daily$dat[5]<-sum(dat3d[1:90,3,])
tot_fits_daily$dat[6]<-sum(dat3d[91:(runTime-14),3,])
tot_fits_daily$dat[7]<-sum(dat3d[1:90,1,])
tot_fits_daily$dat[8]<-sum(dat3d[91:(runTime-14),1,])



### add up total events of each type across sites and days
# events_tot order of results: SM plac3, SM PMC3, UM plac3, UM PMC3, SM plac6, SM PMC6, UM plac6, UM PMC6.

tot_fits_daily$model[1]<-mean(params_gqs$events_tot[,4])
tot_fits_daily$model[2]<-mean(params_gqs$events_tot[,8])
tot_fits_daily$model[3]<-mean(params_gqs$events_tot[,2])
tot_fits_daily$model[4]<-mean(params_gqs$events_tot[,6])
tot_fits_daily$model[5]<-mean(params_gqs$events_tot[,3])
tot_fits_daily$model[6]<-mean(params_gqs$events_tot[,7])
tot_fits_daily$model[7]<-mean(params_gqs$events_tot[,1])
tot_fits_daily$model[8]<-mean(params_gqs$events_tot[,5])

inds<-grep("ci",names(tot_fits_daily))
tot_fits_daily[1,inds]<-quantile(params_gqs$events_tot[,4],probs=c(0.025,0.975))
tot_fits_daily[2,inds]<-quantile(params_gqs$events_tot[,8],probs=c(0.025,0.975))
tot_fits_daily[3,inds]<-quantile(params_gqs$events_tot[,2],probs=c(0.025,0.975))
tot_fits_daily[4,inds]<-quantile(params_gqs$events_tot[,6],probs=c(0.025,0.975))
tot_fits_daily[5,inds]<-quantile(params_gqs$events_tot[,3],probs=c(0.025,0.975))
tot_fits_daily[6,inds]<-quantile(params_gqs$events_tot[,7],probs=c(0.025,0.975))
tot_fits_daily[7,inds]<-quantile(params_gqs$events_tot[,1],probs=c(0.025,0.975))
tot_fits_daily[8,inds]<-quantile(params_gqs$events_tot[,5],probs=c(0.025,0.975))


##########################################
# Relative incidence over time since discharge
##########################################
## Calculate credible intervals
pd_risk<-matrix(nrow=length(params_for_mod$w_scale_risk),ncol=runTime+1)
for(i in 1:nrow(pd_risk)) pd_risk[i,]<-
  params_for_mod$max_beta_phi[i]*
  exp(-((0:runTime/params_for_mod$w_scale_risk[i])^params_for_mod$w_shape_risk[i]))
pd_risk_med<-apply(pd_risk,2,quantile,probs=0.5)
pd_risk_lci<-apply(pd_risk,2,quantile,probs=0.025)
pd_risk_uci<-apply(pd_risk,2,quantile,probs=0.975)

###############################
# EXPECTED EIR IN 0-5 VERSUS SYMPTOMATIC INCIDENCE
###############################
## exposure by age
rho	<-	0.85
a0	<-	2920
eta <- 1/(21*365) # death rate in humans to give approx age structure in under 5's.
foi_age<-1-rho*exp(-(1:(5*365))/a0) # update foi_age with current age
den<-exp(-eta*(1:(5*365)))
# weighted average (weighted by population)
rel_eir_0_5<-sum(den*foi_age)/(sum(den))
# estimate EIR in 0-5 year olds using EIR in adults.
trial_epi$eir_0_5<-trial_epi$eir*rel_eir_0_5


###############################
## ADD MODEL RESULTS TO SHAPE FILE DATA
###############################
## categorise EIR
## put all the highest EIR into one category - beyond the max data EIR, we assume the same as EIR=30
fit_line<-readRDS("fit_line_stan.RDS")
imperial_mod<-imperial_mod %>%
  dplyr::mutate(eir_cat=fit_line$eir[findInterval(eir,fit_line$eir,all.inside = T)])
imperial_mod<-left_join(imperial_mod,fit_line,by=c("eir_cat"="eir"))
## calculate outputs
hosp_name<-c("03","05","07")
hosp_params<-c(0.3,0.5,0.7)
paton_hosp_name<-hosp_name
paton_hosp_params<-hosp_params
pmc_name<-c("no_pmc","pmc")
pmc_params<-c(FALSE,TRUE)

for(h in 1:length(hosp_name)) {
  if(h==1) { 
    ph<-3
  } else if(h==2) {
    ph<-2
  } else {
    ph<-1
  } 
  imperial_mod<- imperial_mod %>% 
      dplyr::mutate(## number of severe episodes averted by PMC in the whole population per year
                !!paste0("repmod_num_sev_avert_pppy_all_h",hosp_name[h],"_ph",paton_hosp_name[ph]):= 
                  !!as.name(paste0("repmod_no_pmc_h",hosp_name[h],"__ph",paton_hosp_name[ph],"_num_sma_pd_pppy")) +
                  !!as.name(paste0("repmod_no_pmc_h",hosp_name[h],"__ph",paton_hosp_name[ph],"_num_sev_oth_pd_pppy")) -
                  !!as.name(paste0("repmod_pmc_h",hosp_name[h],"__ph",paton_hosp_name[ph],"_num_sma_pd_pppy")) -
                  !!as.name(paste0("repmod_pmc_h",hosp_name[h],"__ph",paton_hosp_name[ph],"_num_sev_oth_pd_pppy")),
                  ## number of deaths averted by PMC in the whole population per year
                !!paste0("repmod_deaths_avert_pppy_h",hosp_name[h],"_ph",paton_hosp_name[ph]):= 
                  !!as.name(paste0("repmod_no_pmc_h",hosp_name[h],"__ph",paton_hosp_name[ph],"_deaths_pppy")) -
                  !!as.name(paste0("repmod_pmc_h",hosp_name[h],"__ph",paton_hosp_name[ph],"_deaths_pppy")),
                  ## number of severe episodes averted by 1 PMC
                !!paste0("repmod_num_sev_avert_pp_perPMC_h",hosp_name[h],"_ph",paton_hosp_name[ph]):= 
                  !!as.name(paste0("repmod_num_sev_avert_pppy_all_h",hosp_name[h],"_ph",paton_hosp_name[ph])) /
                  !!as.name(paste0("repmod_pmc_h",hosp_name[h],"__ph",paton_hosp_name[ph],"_num_pmc_eq")),
                ## number needed to treat to avert one severe episode.
      !!paste0("trts_per_sev_avert_h",hosp_name[h],"_ph",paton_hosp_name[ph]):= 
        1 / !!as.name(paste0("repmod_num_sev_avert_pp_perPMC_h",hosp_name[h],"_ph",paton_hosp_name[ph])),
      ## number of deaths averted by 1 PMC
      !!paste0("repmod_num_deaths_avert_pp_perPMC_h",hosp_name[h],"_ph",paton_hosp_name[ph]):= 
        !!as.name(paste0("repmod_deaths_avert_pppy_h",hosp_name[h],"_ph",paton_hosp_name[ph])) /
        !!as.name(paste0("repmod_pmc_h",hosp_name[h],"__ph",paton_hosp_name[ph],"_num_pmc_eq")),
      ## number needed to treat to avert one death
      !!paste0("trts_per_death_avert_h",hosp_name[h],"_ph",paton_hosp_name[ph]):= 
        1 / !!as.name(paste0("repmod_num_deaths_avert_pp_perPMC_h",hosp_name[h],"_ph",paton_hosp_name[ph])))
}

africa@data <- left_join(africa@data, imperial_mod, by="DIDE_CODE")


##### Generate some categories for mapping
breaks_pmc_demand_eq<- seq(0,0.6,by=0.02)
col_palette_sev<-colorRampPalette(rev(c("blue1","cyan")))(length(breaks_pmc_demand_eq)-1)
africa_pmc_demand_eq_cat<-cut(100*africa@data$repmod_pmc_h05__ph05_num_pmc_eq,breaks_pmc_demand_eq,include.lowest=T)
levels(africa_pmc_demand_eq_cat)<-col_palette_sev
africa_pmc_demand_eq_cat<-as.character(africa_pmc_demand_eq_cat)

africa_pmc_demand_eq_cat_low<-cut(100*africa@data$repmod_pmc_h03__ph07_num_pmc_eq,breaks_pmc_demand_eq,include.lowest=T)
levels(africa_pmc_demand_eq_cat_low)<-col_palette_sev
africa_pmc_demand_eq_cat_low<-as.character(africa_pmc_demand_eq_cat_low)

africa_pmc_demand_eq_cat_high<-cut(100*africa@data$repmod_pmc_h07__ph03_num_pmc_eq,breaks_pmc_demand_eq,include.lowest=T)
levels(africa_pmc_demand_eq_cat_high)<-col_palette_sev
africa_pmc_demand_eq_cat_high<-as.character(africa_pmc_demand_eq_cat_high)


### severe prevented per person per PMC.
breaks_sev_episodes_prevent_pp_pd<-seq(0,0.6,0.05)
col_palette_sev_prevent_pd<-colorRampPalette(c("yellow", "orange", "red3", "brown4","maroon4"))(length(breaks_sev_episodes_prevent_pp_pd)-1)
africa_sev_avert_cat<-cut(africa@data$repmod_num_sev_avert_pp_perPMC_h05_ph05,breaks_sev_episodes_prevent_pp_pd,include.lowest=T)
levels(africa_sev_avert_cat)<-col_palette_sev_prevent_pd
africa_sev_avert_cat<-as.character(africa_sev_avert_cat)

### treatments preventing 1 sev episode
breaks_trts_per_sev_avert<-c(seq(0,100,10),seq(150,250,50),350)
col_palette_trts_per_sev<-colorRampPalette(c("yellow", "chartreuse3", "blue","purple"))(length(breaks_trts_per_sev_avert)-1)
africa_trts_per_sev_cat<-cut(africa@data$trts_per_sev_avert_h05_ph05,breaks_trts_per_sev_avert,include.lowest=T)
inds<-which(africa@data$trts_per_sev_avert_h05_ph05>350) ## set above this to max category (will later label '350+')
africa_trts_per_sev_cat[inds]<-levels(africa_trts_per_sev_cat)[length(levels(africa_trts_per_sev_cat))]
levels(africa_trts_per_sev_cat)<-col_palette_trts_per_sev
africa_trts_per_sev_cat<-as.character(africa_trts_per_sev_cat)



#############################
# FIGURES 
#############################


###########################################
# MODEL FIT TO TRIAL DATA, TOTAL COUNTS & BY DAY
## save underlying model fits for upload to github
saveRDS(tot_fits_daily %>% dplyr::select(-dat), "tot_fits_mod_git.rds")
tot_fits_daily<-readRDS("tot_fits_mod_git.rds")
saveRDS(fit_daily %>% dplyr::select(days,model_plac_sm:model_pmc_um_uci), "daily_fits_mod_git.rds")
fit_daily<-readRDS("daily_fits_mod_git.rds")

if(plot2file) {
  ####################
  # Fit to total counts
  ####################
  tiff(file="trial_fits_weib.tiff",width = 1200,height = 2700,res=300,
       compression = "lzw")
  par(mar=c(4.5,4,4,2),mfrow=c(3,1))
  plotCI((1:8)+0.2,tot_fits_daily$model,ui=tot_fits_daily$model_uci,li=tot_fits_daily$model_lci,
         col="black",pch=1,gap=0.23,ylab="Number of events",xlab="",xlim=c(0,9),
         sfrac=0,xaxt="n",ylim=c(0,300),las=1)
  plotCI((1:4)+0.2,tot_fits_daily$model[1:4],ui=tot_fits_daily$model_uci[1:4],li=tot_fits_daily$model_lci[1:4],
         col="blue",pch=1,gap=0.23,add=T,sfrac=0)
  points(1:8,tot_fits_daily$dat,col="black",pch=19)
  points(1:4,tot_fits_daily$dat[1:4],col="blue",pch=19)
  lines(c(4.5,4.5),c(-5,350),lty=2)
  mtext("PMC",at=3.5,cex=0.7)
  mtext("placebo",at=5.5,cex=0.7)
  mtext("A",at=-0.9,line=0.5,cex=1)
  legend(7.5,370,c("data","model"),pch=c(19,1),col="black",lty=c(NA,1),bty='n',xpd=T)
  axis(1,at=1:8,labels=c("UM 3-14","UM 15-25","SM 3-14","SM 15-25","UM 3-14","UM 15-25","SM 3-14","SM 15-25"),
       tot_fits_daily$var,las=2)
  
  ###################################
  # Daily fit to events.
  ###################################
  par(mar=c(4,4,2,2))
  plot(fit_daily$days,fit_daily$dat_plac_sm,type="l",xlab="days since discharge",
       ylab="cumulative number of hospialized cases")
  polygon(c(fit_daily$days,rev(fit_daily$days)),
          c(fit_daily$model_plac_sm_lci,rev(fit_daily$model_plac_sm_uci)), 
          col="gray90",border=NA)
  lines(fit_daily$days,fit_daily$model_plac_sm,col="black",lty=2)
  lines(fit_daily$days,fit_daily$dat_plac_sm,col="black",lwd=2)
  polygon(c(fit_daily$days,rev(fit_daily$days)),
          c(fit_daily$model_pmc_sm_lci,rev(fit_daily$model_pmc_sm_uci)), 
          col="lightblue",border=NA)
  lines(fit_daily$days,fit_daily$model_pmc_sm,col="dodgerblue",lty=2)
  lines(fit_daily$days,fit_daily$dat_pmc_sm,col="blue",lwd=2)
  legend("topleft",c("data placebo","data PMC","model placebo","model PMC"),
         lty=c(1,1,2,2),lwd=c(2,2,1,1), bty='n',
         col=c("black","blue","black","dodgerblue"))
  mtext("B",at=-1.7,line=0.5,cex=1)
  
  
  plot(fit_daily$days,fit_daily$dat_plac_um,type="l",xlab="days since discharge",
       ylab="cumulative number of uncomplicated cases")
  polygon(c(fit_daily$days,rev(fit_daily$days)),
          c(fit_daily$model_plac_um_lci,rev(fit_daily$model_plac_um_uci)), 
          col="gray90",border=NA)
  lines(fit_daily$days,fit_daily$model_plac_um,col="black",lty=2)
  lines(fit_daily$days,fit_daily$dat_plac_um,col="black",lwd=2)
  polygon(c(fit_daily$days,rev(fit_daily$days)),
          c(fit_daily$model_pmc_um_lci,rev(fit_daily$model_pmc_um_uci)), 
          col="lightblue",border=NA)
  lines(fit_daily$days,fit_daily$model_pmc_um,col="dodgerblue",lty=2)
  lines(fit_daily$days,fit_daily$dat_pmc_um,col="blue",lwd=2)
  mtext("C",at=-1.7,line=0.5,cex=1)
  
  dev.off()
}


############### 
# MODEL FIT vs CASES by EIR
###############
fit_line<-readRDS("fit_line_stan.RDS")
if(plot2file) {
  tiff(file="trial_fit_by_eir.tiff",width = 2200,height = 2500,res=300,
       compression = "lzw")
  par(mfrow=c(2,1))
  par(mar=c(4,4,2,13))
  ## Uncomplicated malaria
  inds<-which(trial_epi$study=="Kwambai")
  plotCI(trial_epi$eir[inds],trial_epi$plac_um_inc_p100py3[inds], ui=trial_epi$plac_um_uci3[inds], 
         li=trial_epi$plac_um_lci3[inds], xlab="Prior annual EIR", xlim=c(0,30),
         ylab="Uncomplicated incidence (per 100 p-y)",pch=19,col="firebrick",gap=0,
         sfrac=0,ylim=c(0,max(trial_epi$plac_um_uci3[inds])))
  polygon(c(fit_line$eir,rev(fit_line$eir)),
          c(fit_line$pred_plac_um_inc3_lci,rev(fit_line$pred_plac_um_inc3_uci)),
          col="lightblue",border="lightblue")
  points(imperial_mod$eir,imperial_mod$clin_inc_0_5*100,col="black",pch=19)
  plotCI(trial_epi$eir[inds],trial_epi$plac_um_inc_p100py3[inds], ui=trial_epi$plac_um_uci3[inds], 
         li=trial_epi$plac_um_lci3[inds],pch=19,col="firebrick",gap=0,
         sfrac=0,add=T)
  lines(fit_line$eir,fit_line$pred_plac_um_inc3,col="blue")
  mtext("A",at=-5.5,line=0.5,cex=1)
  legend(32,150,c("data post-discharge (trial)","data post-discharge\n (comparison)",
                  "model post-discharge","model general popn 0-5yr"),
         xpd=T,pch=c(19,19,NA,19),lty=c(1,1,1,NA),col=c("firebrick","purple","blue","black"),
         bty='n',cex=1)
  
  # hospitalised malaria
  par(mar=c(4,4,1,13))
  plotCI(trial_epi$eir,pmax(trial_epi$plac_sm_inc_p100py3,0.01), 
         ui=trial_epi$plac_sm_uci3, 
         li=pmax(trial_epi$plac_sm_lci3,0.01), xlab="Prior annual EIR", xlim=c(0,30),
         ylab="Hospitalized incidence (per 100 p-y)",pch=19,col="firebrick",gap=0,
         sfrac=0,log='y',ylim=c(0.01,500),yaxt='n')
  ats<-c(0.01,0.1,1,10,100,500)
  axis(2,at=ats,labels=ats,las=1)
  polygon(c(fit_line$eir,rev(fit_line$eir)),
          c(fit_line$pred_plac_sm_inc3_lci,rev(fit_line$pred_plac_sm_inc3_uci)),
          col="lightblue",border="lightblue")
  points(imperial_mod$eir,imperial_mod$sev_inc_0_5*100,col="black",pch=19)
  plotCI(trial_epi$eir,pmax(trial_epi$plac_sm_inc_p100py3,0.01), 
         ui=trial_epi$plac_sm_uci3, 
         li=pmax(trial_epi$plac_sm_lci3,0.01),pch=19,col="firebrick",gap=0,
         sfrac=0,add=T)
  inds<-which(trial_epi$study!="Kwambai")
  plotCI(trial_epi$eir[inds],
         trial_epi$plac_sm_inc_p100py3[inds], ui=trial_epi$plac_sm_uci3[inds], 
         li=trial_epi$plac_sm_lci3[inds], pch=19,col="purple",gap=0,
         sfrac=0,add=T)
  plotCI(11.7, 47.8, ui=62.5, 
         li=35.8, pch=19,col="purple",gap=0,
         sfrac=0,add=T)
  lines(fit_line$eir,fit_line$pred_plac_sm_inc3,col="blue")
  mtext("B",at=-5.5,line=0.5,cex=1)
  dev.off()
  
}

#########################################
# Combined map and contribution of post discharge to total morbidity
#########################################
# Model predictions are in imperial_mod (processed prior to plotting)
if(plot2file) {
  tiff(file="pmc_impact_maps2.tiff", width=5000,height=1300,res=300,compression="lzw")
  layout(matrix(c(1,2,3,4,5),1,5), widths=c(7,1,7,1,5.5), heights=c(1,1,1,1,1))
  par(mar=c(0,0,0,0))
  ## Map of hospitalised episodes averted per PMC
  plot(shp2, lwd=0.02,xlim=c(-25,70),ylim=c(-35,35))
  plot(africa, add=T, col=africa_sev_avert_cat, border=africa_sev_avert_cat)
  plot(shp2, lwd=0.02,xlim=c(-25,70),ylim=c(-35,35),add=T)
  text(-24,34,"A",cex=3)
  # colour scale
  labels<-as.character(round(breaks_sev_episodes_prevent_pp_pd,digits=2))
  labels[1:(floor(length(labels)/2))*2]<-""
  incr<- 400 / (length(breaks_sev_episodes_prevent_pp_pd)-1)
  ticks=seq(from= -400, to=0, by=incr)
  par(mar=c(3,0,3,4))
  color.bar(col_palette_sev_prevent_pd, min(ticks),max=max(max(ticks)), 
            ticks=ticks, labels=labels,cex.axis0 = 1.5)
  
  # Number needed to treat map
  par(mar=c(0,0,0,0))
  plot(shp2, lwd=0.02,xlim=c(-25,70),ylim=c(-35,35))
  plot(africa, add=T, col=africa_trts_per_sev_cat, border=africa_trts_per_sev_cat)
  plot(shp2, lwd=0.02,xlim=c(-25,70),ylim=c(-35,35),add=T)
  text(-24,34,"B",cex=3)
  # colour scale
  labels<-as.character(round(breaks_trts_per_sev_avert,digits=2))
  labels[length(labels)]<-paste0(labels[length(labels)],"+")
  incr<- 400 / (length(breaks_trts_per_sev_avert)-1)
  ticks=seq(from= -400, to=0, by=incr)
  par(mar=c(3,0,3,4))
  color.bar(col_palette_trts_per_sev, min(ticks),max=max(max(ticks)), 
            ticks=ticks, labels=labels,cex.axis0 = 1.5)

  # Contribution of recurrent SMA episodes to total SMA burden
  par(mar=c(5,4.5,7,2))
  plot(imperial_mod$eir,100*imperial_mod$repmod_no_pmc_h05__ph05_inc_sma_ppy,col="black",pch=19,
       xlab="annual EIR",
       ylab="SMA incidence per 100 p-y",las=1,cex=0.9, cex.axis=1.4,cex.lab=1.5)
  points(imperial_mod$eir,100*(imperial_mod$repmod_no_pmc_h05__ph05_num_sma_pd_pppy),col="black",
         pch=21,bg="purple",cex=0.9)
  points(imperial_mod$eir,100*(imperial_mod$repmod_pmc_h05__ph05_inc_sma_ppy),col="black",pch=21,
         bg="white",cex=0.9)
  points(imperial_mod$eir,100*(imperial_mod$repmod_pmc_h05__ph05_num_sma_pd_pppy),col="purple",pch=21,
         bg="white",cex=0.9)
  legend(10,0.85,c("Total episodes, no PMC","Recurrent episodes<25 wks, no PMC",
                   "Total episodes, 100% PMC coverage","Recurrent episodes<25 wks, 100% PMC coverage"),
         col=c("black","purple","black","purple"),
         pch=c(19,19,1,1),xpd=T,bty='n',cex=1.2)
  mtext("C",at=-25, line=3.3,cex=2)
  dev.off()
  
}


########################################
### Protection over time in PMC vs placebo
########################################
tiff(file="p_protect_by_arm.tiff", width=3500,height=1300,res=300,compression="lzw")
par(mfrow=c(1,2))
plot(PE$time,PE$pp_plac,xlab="days since discharge",ylab="% protected",
     col="blue",type="l")
mtext("AL",at=0)
plot(PE$time,PE$pp_pmc,xlab="days since discharge",ylab="% protected",type="l",
     col="blue")
mtext("AL",at=0)
mtext("DP",at=pmc_times)
dev.off()


###########################################
# Decline in risk over time since PMC
###########################################
tiff(file="risk_over_time_pd.tiff", width=2000,height=1300,res=300,compression="lzw")
plot(15:runTime,pd_risk_med[15:runTime], type="l",
     ylab="relative incidence, placebo",
     xlab="days since hospital discharge", ylim=c(0,1.1),las=1)
polygon(c(15:runTime,runTime:15),c(pd_risk_lci[15:runTime],rev(pd_risk_uci[15:runTime])), col="lightblue", border = NA )
lines(15:runTime, pd_risk_med[15:runTime])
dev.off()


######################
# RELATIVE INCIDENCE PER BITE AS EIR CHANGES
######################
if(plot2file) {
  tiff(file="inc_symptoms_eir_supp.tiff", width=1500,height=1300,res=300,compression="lzw")
  plot(0.1:30,median(params_for_mod$max_beta_phi)*exp(-median(params_for_mod$r)*0.1:30),
       type="l",ylab="relative incidence, placebo",xlab="EIR",ylim=c(0,1.2))
  dev.off()
}


#########################################
# EIR IN 0-5 YR OLDS VS TOTAL INCIDENCE
#########################################
tiff(file="supp_eir_0_5_vs_inc.tiff", width=1500,height=1400,res=300,compression="lzw")
plotCI(trial_epi$eir_0_5, trial_epi$tot_inc_plac_um_sm_p100py/100,
       ui=trial_epi$tot_inc_plac_um_sm_uci/100,
       li=trial_epi$tot_inc_plac_um_sm_lci/100,
       sfrac=0,gap=0,
       ylim=c(0,7),xlab="Prior annual EIR, 0-5 yr olds", ylab="Incidence of symptomatic malaria (UM + SM) pppy",
       xlim=c(0,11),pch=19,col="firebrick",las=1)
inds<-grep("Opoka", trial_epi$study)
plotCI(trial_epi$eir_0_5[inds], trial_epi$tot_inc_plac_um_sm_p100py[inds]/100,
       ui=trial_epi$tot_inc_plac_um_sm_uci[inds]/100,
       li=trial_epi$tot_inc_plac_um_sm_lci[inds]/100,
       sfrac=0,gap=0,
       pch=19,col="purple",add=T)
lines(c(0,11),c(0,11),lty=2)
legend(0,9,c("data post-discharge (K)","data post-discharge (O)"),
       xpd=T,pch=c(19,19),lty=c(1,1),col=c("firebrick","purple"),
       bty='n',cex=1)
dev.off()

#####################################
# WHY EXCLUDED LAST WEEK OF FOLLOW UP
#####################################
#### Plot full data
dat2d<-dat3d[,,1]
for(i in 2:9) dat2d<-dat2d+dat3d[,,i]
dat2d<-dat2d[1:168,]
fit_daily<-data.frame(days=1:nrow(dat2d)+14, 
                      dat_plac_sm=cumsum(dat2d[,1]),
                      dat_pmc_sm=cumsum(dat2d[,2]),
                      dat_plac_um=cumsum(dat2d[,3]),
                      dat_pmc_um=cumsum(dat2d[,4]))
if(plot2file) {
  tiff(file="supp_cut_last_week.tiff", width=2300,height=1400,res=300,compression="lzw")
  par(mfrow=c(1,2))
  plot(fit_daily$days,fit_daily$dat_plac_um,type="l",xlab="days since discharge",
       ylab="cumulative number of uncomplicated cases")
  polygon(c(175,182,182,175),c(0,0,500,500),col="orange1",border=NA)
  lines(fit_daily$days,fit_daily$dat_pmc_um)
  lines(fit_daily$days,fit_daily$dat_plac_um)
  plot(fit_daily$days,fit_daily$dat_plac_sm,type="l",xlab="days since discharge",
       ylab="cumulative number of hospitalised cases")
  polygon(c(175,182,182,175),c(0,0,500,500),col="orange1",border=NA)
  lines(fit_daily$days,fit_daily$dat_pmc_sm)
  lines(fit_daily$days,fit_daily$dat_plac_sm)
  dev.off()
  
}


### TOTAL VERSUS RECURRENT SMA - sensitivity analysis: hospitalisation rate

par(mfrow=c(1,3),mar=c(5,4,7,2))
plot(imperial_mod$eir,100*imperial_mod$repmod_no_pmc_h07__ph03_inc_sma_ppy,col="black",pch=19,
     xlab="annual EIR",
     ylab="SMA incidence per 100 p-y",las=1,cex=0.9)
points(imperial_mod$eir,100*(imperial_mod$repmod_no_pmc_h07__ph03_num_sma_pd_pppy),col="black",
       pch=21,bg="purple",cex=0.9)
points(imperial_mod$eir,100*(imperial_mod$repmod_pmc_h07__ph03_inc_sma_ppy),col="black",pch=21,
       bg="white",cex=0.9)
points(imperial_mod$eir,100*(imperial_mod$repmod_pmc_h07__ph03_num_sma_pd_pppy),col="purple",pch=21,
       bg="white",cex=0.9)
mtext("A",at = 0)
plot(imperial_mod$eir,100*imperial_mod$repmod_no_pmc_h05__ph05_inc_sma_ppy,col="black",pch=19,
     xlab="annual EIR",
     ylab="SMA incidence per 100 p-y",las=1,cex=0.9)
points(imperial_mod$eir,100*(imperial_mod$repmod_no_pmc_h05__ph05_num_sma_pd_pppy),col="black",
       pch=21,bg="purple",cex=0.9)
points(imperial_mod$eir,100*(imperial_mod$repmod_pmc_h05__ph05_inc_sma_ppy),col="black",pch=21,
       bg="white",cex=0.9)
points(imperial_mod$eir,100*(imperial_mod$repmod_pmc_h05__ph05_num_sma_pd_pppy),col="purple",pch=21,
       bg="white",cex=0.9)
mtext("B",at = 0)
legend(0,0.85,c("Total episodes, no PMC","Recurrent episodes<25 wks, no PMC",
                "Total episodes, 100% PMC coverage","Recurrent episodes<25 wks,\n 100% PMC coverage"),
       col=c("black","purple","black","purple"),
       pch=c(19,19,1,1),xpd=T,bty='n')
plot(imperial_mod$eir,100*imperial_mod$repmod_no_pmc_h03__ph07_inc_sma_ppy,col="black",pch=19,
     xlab="annual EIR",
     ylab="SMA incidence per 100 p-y",las=1,cex=0.9)
points(imperial_mod$eir,100*(imperial_mod$repmod_no_pmc_h03__ph07_num_sma_pd_pppy),col="black",
       pch=21,bg="purple",cex=0.9)
points(imperial_mod$eir,100*(imperial_mod$repmod_pmc_h03__ph07_inc_sma_ppy),col="black",pch=21,
       bg="white",cex=0.9)
points(imperial_mod$eir,100*(imperial_mod$repmod_pmc_h03__ph07_num_sma_pd_pppy),col="purple",pch=21,
       bg="white",cex=0.9)
mtext("C",at = 0)


###########################################
# PMC DEMAND - SMA INCIDENCE IN HOSPITALS (minus mortality, with and without PMC)
###########################################
if(plot2file) {
  tiff(file="pmc_demand_maps.tiff", width=3800,height=1300,res=300,compression="lzw")
  layout(matrix(c(1,2,3,4),1,4), widths=c(2,5,5,5), heights=c(1,1,1,1))
  # colour scale
  labels<-round(breaks_pmc_demand_eq,digits=2)
  incr<- 400 / (length(breaks_pmc_demand_eq)-1)
  ticks=seq(from= -400, to=0, by=incr)
  par(mar=c(3,5.5,3,4))
  color.bar(col_palette_sev, min(ticks),max=max(max(ticks)), 
            ticks=ticks, labels=labels,cex.axis0 = 1.5)
  print(2)
  par(mar=c(0,0,0,0))
  plot(shp2, lwd=0.02,xlim=c(-20,55),ylim=c(-35,35))
  plot(africa, add=T, col=africa_pmc_demand_eq_cat_low, border=africa_pmc_demand_eq_cat_low)
  plot(shp2, lwd=0.02,add=T)
  text(-20,35,"A",cex=2)
  print(1)
  
  par(mar=c(0,0,0,0))
  plot(shp2, lwd=0.02,xlim=c(-20,55),ylim=c(-35,35))
  plot(africa, add=T, col=africa_pmc_demand_eq_cat, border=africa_pmc_demand_eq_cat)
  plot(shp2, lwd=0.02,add=T)
  text(-20,35,"B",cex=2)
  print(3)
  
  par(mar=c(0,0,0,0))
  plot(shp2, lwd=0.02,xlim=c(-20,55),ylim=c(-35,35))
  plot(africa, add=T, col=africa_pmc_demand_eq_cat_high, border=africa_pmc_demand_eq_cat_high)
  plot(shp2, lwd=0.02,add=T)
  text(-20,35,"C",cex=2)
  print(4)
  dev.off()

}


################################################
# NUMERIC RESULTS BY COUNTRY FOR TABLE
################################################
country_summ<-imperial_mod %>%  ## take weighted averages of the adm units (wgt by population)
  dplyr::group_by(NAME_0) %>%
  dplyr::summarise(pop_0_5=round(sum(population_0_5)),
                   inc_sma_paton_p100py=round(100*sum(inc_sma_ppy_0_5_paton*population_0_5)/sum(population_0_5),2),
                   inc_sma_p100py=round(100*sum(repmod_pmc_h05__ph05_inc_sma_ppy*population_0_5)/sum(population_0_5),2),
                   inc_sma_p100py_low=round(100*sum(repmod_pmc_h03__ph07_inc_sma_ppy*population_0_5)/sum(population_0_5),2),
                   inc_sma_p100py_high=round(100*sum(repmod_pmc_h07__ph03_inc_sma_ppy*population_0_5)/sum(population_0_5),2),
                   num_pmc_y_adm_init=(1-0.074)*round(sum(repmod_no_pmc_h05__ph05_inc_sma_all_hosp_eq*population_0_5)),
                   num_pmc_y_adm_init_low=(1-0.074)*round(sum(repmod_no_pmc_h03__ph07_inc_sma_all_hosp_eq*population_0_5)),
                   num_pmc_y_adm_init_high=(1-0.074)*round(sum(repmod_no_pmc_h07__ph03_inc_sma_all_hosp_eq*population_0_5)),
                   num_pmc_y_adm=round(sum(repmod_pmc_h05__ph05_num_pmc_eq*population_0_5)),
                   num_pmc_y_adm_low=round(sum(repmod_pmc_h03__ph07_num_pmc_eq*population_0_5)),
                   num_pmc_y_adm_high=round(sum(repmod_pmc_h07__ph03_num_pmc_eq*population_0_5)),
                   sev_avert_adm_py_all=sum(repmod_num_sev_avert_pppy_all_h05_ph05*population_0_5),
                   sev_avert_adm_py_all_low=sum(repmod_num_sev_avert_pppy_all_h03_ph07*population_0_5),
                   sev_avert_adm_py_all_high=sum(repmod_num_sev_avert_pppy_all_h07_ph03*population_0_5),
                   sev_avert_per_100_pmc=round(100*sev_avert_adm_py_all/num_pmc_y_adm,1),
                   pmc_per_sev_avert=round(num_pmc_y_adm/sev_avert_adm_py_all,1),
                   deaths_avert_adm_py_h05_ph05 = sum(repmod_deaths_avert_pppy_h05_ph05*population_0_5),
                   pmc_per_deaths_avert_h05_ph05=round(num_pmc_y_adm/deaths_avert_adm_py_h05_ph05,0),
                   deaths_avert_adm_py_h07_ph03 = sum(repmod_deaths_avert_pppy_h07_ph03*population_0_5),
                   pmc_per_deaths_avert_h07_ph03=round(num_pmc_y_adm_high/deaths_avert_adm_py_h07_ph03,0),
                   deaths_avert_adm_py_h03_ph07 = sum(repmod_deaths_avert_pppy_h03_ph07*population_0_5),
                   pmc_per_deaths_avert_h03_ph07=round(num_pmc_y_adm_low/deaths_avert_adm_py_h03_ph07,0)
                   
  )

## reorder for table
country_summ_table<-select(country_summ,NAME_0,pop_0_5,inc_sma_p100py,
                           inc_sma_p100py_low,inc_sma_p100py_high,
                           sev_avert_per_100_pmc,
                     pmc_per_sev_avert, pmc_per_deaths_avert_h05_ph05,
                     pmc_per_deaths_avert_h03_ph07, pmc_per_deaths_avert_h07_ph03,
                     num_pmc_y_adm,num_pmc_y_adm_low,num_pmc_y_adm_high
                     )
hbhi_inds<-which(country_summ_table$NAME_0 %in% c("Burkina Faso","Cameroon","Nigeria",
                                            "Democratic Republic of the Congo",
                                            "Mozambique","Uganda","Ghana",
                                            "Niger","Mali","Tanzania"))
country_summ_table<- rbind(country_summ_table[hbhi_inds,],country_summ_table[!(1:nrow(country_summ_table) %in% hbhi_inds),])

write.csv(country_summ_table,"table_country_summ_2019_full.csv")

### Format for table in the manuscript
country_summ_table<-country_summ_table %>%
  mutate(pop_0_5_text=format(pop_0_5, nsmall=0, big.mark=","),
    inc_sma_p100py_text=paste0(inc_sma_p100py," (",inc_sma_p100py_low,"-",inc_sma_p100py_high,")"),
         pmc_per_deaths_avert_text=paste0(pmc_per_deaths_avert_h05_ph05," (",pmc_per_deaths_avert_h03_ph07,"-",pmc_per_deaths_avert_h07_ph03,")"),
    num_pmc_y_adm_text=format(num_pmc_y_adm, nsmall=0, big.mark=","),
    num_pmc_y_adm_low_text=format(num_pmc_y_adm_low, nsmall=0, big.mark=","),
    num_pmc_y_adm_high_text=format(num_pmc_y_adm_high, nsmall=0, big.mark=","),
    num_pmc_y_text=paste0(num_pmc_y_adm_text,
                               " (",num_pmc_y_adm_low_text, "-",
                               num_pmc_y_adm_high_text,
                               ")")
  )
country_summ_table_text<-select(country_summ_table,NAME_0,pop_0_5_text,inc_sma_p100py_text,
                           sev_avert_per_100_pmc,
                           pmc_per_sev_avert, pmc_per_deaths_avert_text,
                           num_pmc_y_text
)

write.csv(country_summ_table_text,"table_country_summ_2019.csv")


