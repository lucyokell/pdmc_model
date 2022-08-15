
##############################################
# RUN PMC REPEAT MODEL ON THE CLUSTER
##############################################

###########################################################################################
### SOURCE CODE
###########################################################################################
library(dplyr)
library(reshape2)
library(ggplot2)
source("repeat_model_pmc6c_function.R")
source("fit_inc_sma_function.R")

##################################
# READ IN IMPERIAL MODEL OUTPUTS
##################################
## to do update input after running Paton.R
imperial_mod<-read.csv("adm1_epi_2019_sma3_paton.csv")

##################################
# READ IN FITTED PARAMS
##################################
fit<-readRDS("pmc14g_cut8_fit_paramsonly_correct_negbin.rds")
params<-fit
r_eir<-median(params$r)
prob_um_plac<-median(params$prob_um_plac)
prob_um_pmc<-median(params$prob_um_pmc)
max_beta_phi<-median(params$max_beta_phi)
w_scale_risk<-median(params$w_scale_risk)
w_shape_risk<-median(params$w_shape_risk)

#########################################
# CHECK THE MODEL RUNS LOCALLY WITHOUT ERRORS AND GIVES SENSIBLE OUTPUT
#########################################
source("repeat_model_pmc6c_function.R")
## summary output of incidence
i<-52  # pick an admin-1 area
run_repeat_mod(eir=imperial_mod$eir[i],inc_sma_genpop = imperial_mod$inc_sma_ppy_0_5_paton[i],
               pmc=F,ft_um_pd0=1,ft_sev_pd0=1,prop_hosp=0.5,output="summary2R",
               r_eir=r_eir,
               prob_um_plac=prob_um_plac,
               prob_um_pmc=prob_um_pmc,
               max_beta_phi=max_beta_phi,
               w_scale_risk=w_scale_risk,
               w_shape_risk=w_shape_risk,
               adm="test")

# now with output of state variables over time
i<-52
x<-run_repeat_mod(eir=imperial_mod$eir[i],inc_sma_genpop = imperial_mod$inc_sma_ppy[i],
               pmc=TRUE,ft_sev_pd0=0.5,prop_hosp=0.5,prop_hosp_paton=0.5,output="x",
               r_eir=r_eir,
               prob_um_plac=prob_um_plac,
               prob_um_pmc=prob_um_pmc,
               max_beta_phi=max_beta_phi,
               w_scale_risk=w_scale_risk,
               w_shape_risk=w_shape_risk,
               adm="test")

# Plot example output 
x<-melt(x,id.vars = "time")
ggplot(data=x,aes(time,value)) +
  geom_line(color="blue",size=1) +
  facet_wrap(~variable, scales= "free_y") +
  theme_bw()


#########################################################
# RECALIBRATE THE MODEL to include post-SMA incidence from PMC study placebo group
### match the total SMA incidence in the population when this post-SMA incidence is included
### do this by having lower risks of SMA in those not experiencing SMA in the past 6 months
#########################################################
i<-52  # choose an example admin-1 subnational area.
inc_sma_genpop0<-imperial_mod$inc_sma_ppy_0_5_paton[i]
eir0<-imperial_mod$eir[i]

## check the fitting function works.
fit_inc_sma(inc_sma_genpop=0.002154645,
            target_tot_inc_sma=inc_sma_genpop0,
            eir=eir0,
            prop_hosp0=0.7,
            prop_hosp_paton0=0.7,
            r_eir=r_eir,
            prob_um_plac=prob_um_plac,
            prob_um_pmc=prob_um_pmc,
            max_beta_phi=max_beta_phi,
            w_scale_risk=w_scale_risk,
            w_shape_risk=w_shape_risk
)

## now find the best fit
prop_hosp_paton0<-0.5
t2<-optimise(f=fit_inc_sma, interval=c(0,inc_sma_genpop0/prop_hosp_paton0),
                         target_tot_inc_sma = inc_sma_genpop0,eir=eir0,
                         prop_hosp0=0.5,
                         prop_hosp_paton0=prop_hosp_paton0,
                         r_eir=r_eir,
                         prob_um_plac=prob_um_plac,
                         prob_um_pmc=prob_um_pmc,
                         max_beta_phi=max_beta_phi,
                         w_scale_risk=w_scale_risk,
                         w_shape_risk=w_shape_risk, 
                         tol=1e-4)

t2


##### Check model calibration worked
## Input the fitted SMA incidence in the general population to check the total outputted
    # matches the total in the IC model
run_repeat_mod(inc_sma_genpop=t2$minimum,
               eir=eir0,pmc=F,ft_um_pd0=0.5,ft_sev_pd0 = 0.5,prop_hosp=0.5,output="summary2R",
               prop_hosp_paton=0.5,
               r_eir=r_eir,
               prob_um_plac=prob_um_plac,
               prob_um_pmc=prob_um_pmc,
               max_beta_phi=max_beta_phi,
               w_scale_risk=w_scale_risk,
               w_shape_risk=w_shape_risk, 
               adm="test")
#target total inc - compare with 'inc_sma_all_hosp_eq'
imperial_mod$inc_sma_ppy_0_5_paton[i]

###################################
# Recalibrate SMA all admin areas
###################################
##### 70% hospitalization
prop_hosp_paton0<-0.7
### this loop will take some hours to run - consider using a cluster or read in the results (below)
for(i in 1:nrow(imperial_mod)) {
    print(i)
    inc_sma_genpop0<-imperial_mod$inc_sma_ppy_0_5_paton[i]
    eir0<-imperial_mod$eir[i]
    assign(paste0("w",i),optimise(f=fit_inc_sma, interval=c(0,inc_sma_genpop0/prop_hosp_paton0),
                                              target_tot_inc_sma = inc_sma_genpop0,eir=eir0,
                                              prop_hosp0=0.7,prop_hosp_paton=prop_hosp_paton0,
                                              r_eir=r_eir,
                                              prob_um_plac=prob_um_plac,
                                              prob_um_pmc=prob_um_pmc,
                                              max_beta_phi=max_beta_phi,
                                              w_scale_risk=w_scale_risk,
                                              w_shape_risk=w_shape_risk, 
                                              tol=0.00001))
}

for(j in 1:nrow(imperial_mod)) saveRDS(eval(parse(text=paste0("w",j,"$minimum"))), file=paste0("w",j,".rds"))

# check all have run
try_read<-function(i) { try(readRDS(paste0("w",i,".rds"))) }
i<-1:nrow(imperial_mod)
check<-unlist(purrr::map(i,try_read))
check

imperial_mod<-read.csv("adm1_epi_2019_sma3_paton.csv")
imperial_mod$inc_sma_ppy_genpop_prophosp_paton07<-NA
for(i in 1:nrow(imperial_mod)) {
  if(imperial_mod$inc_sma_ppy[i]>0) {
    temp<-readRDS(paste0("output/w",i,".rds"))
    imperial_mod$inc_sma_ppy_genpop_prophosp_paton07[i]<-as.numeric(temp)
  } else {
    imperial_mod$inc_sma_ppy_genpop_prophosp_paton07[i]<-0
  }
}


##### 50% hospitalization
prop_hosp_paton0<-0.5
### this loop will take some hours to run - consider using a cluster or read in the results (below)
for(i in 1:nrow(imperial_mod)) {
  print(i)
  inc_sma_genpop0<-imperial_mod$inc_sma_ppy_0_5_paton[i]
  eir0<-imperial_mod$eir[i]
  assign(paste0("z",i),optimise(f=fit_inc_sma, interval=c(0,inc_sma_genpop0/prop_hosp_paton0),
                                target_tot_inc_sma = inc_sma_genpop0,eir=eir0,
                                prop_hosp0=0.5,prop_hosp_paton=prop_hosp_paton0,
                                r_eir=r_eir,
                                prob_um_plac=prob_um_plac,
                                prob_um_pmc=prob_um_pmc,
                                max_beta_phi=max_beta_phi,
                                w_scale_risk=w_scale_risk,
                                w_shape_risk=w_shape_risk, 
                                tol=0.00001))
}

for(j in 1:nrow(imperial_mod)) saveRDS(eval(parse(text=paste0("z",j,"$minimum"))), file=paste0("z",j,".rds"))

imperial_mod<-read.csv("adm1_epi_2019_sma3_paton.csv")
imperial_mod$inc_sma_ppy_genpop_prophosp_paton05<-NA
for(i in 1:nrow(imperial_mod)) {
  if(imperial_mod$inc_sma_ppy[i]>0) {
    temp<-readRDS(paste0("output/w",i,".rds"))
    imperial_mod$inc_sma_ppy_genpop_prophosp_paton05[i]<-as.numeric(temp)
  } else {
    imperial_mod$inc_sma_ppy_genpop_prophosp_paton05[i]<-0
  }
}


##### 30% hospitalization
prop_hosp_paton0<-0.3
### this loop will take some hours to run - consider using a cluster or read in the results (below)
for(i in 1:nrow(imperial_mod)) {
  print(i)
  inc_sma_genpop0<-imperial_mod$inc_sma_ppy_0_5_paton[i]
  eir0<-imperial_mod$eir[i]
  assign(paste0("y",i),optimise(f=fit_inc_sma, interval=c(0,inc_sma_genpop0/prop_hosp_paton0),
                                target_tot_inc_sma = inc_sma_genpop0,eir=eir0,
                                prop_hosp0=0.3,prop_hosp_paton=prop_hosp_paton0,
                                r_eir=r_eir,
                                prob_um_plac=prob_um_plac,
                                prob_um_pmc=prob_um_pmc,
                                max_beta_phi=max_beta_phi,
                                w_scale_risk=w_scale_risk,
                                w_shape_risk=w_shape_risk, 
                                tol=0.00001))
}

for(j in 1:nrow(imperial_mod)) saveRDS(eval(parse(text=paste0("y",j,"$minimum"))), file=paste0("y",j,".rds"))

imperial_mod<-read.csv("adm1_epi_2019_sma3_paton.csv")
imperial_mod$inc_sma_ppy_genpop_prophosp_paton03<-NA
for(i in 1:nrow(imperial_mod)) {
  if(imperial_mod$inc_sma_ppy[i]>0) {
    temp<-readRDS(paste0("output/y",i,".rds"))
    imperial_mod$inc_sma_ppy_genpop_prophosp_paton03[i]<-as.numeric(temp)
  } else {
    imperial_mod$inc_sma_ppy_genpop_prophosp_paton03[i]<-0
  }
}

write.csv(imperial_mod, "adm1_epi_2019_sma3_paton_cal.csv",row.names = F)


####################
## Now run calibrated model with and without PMC, with different prop hosp (x4 sets of runs)
####################
imperial_mod<-read.csv("adm1_epi_2019_sma3_paton_cal.csv")

## scenarios = 30% hospitalisation in Paton et al, 70% in other settings
## scenarios = 50% hospitalisation in Paton et al, 50% in other settings
## scenarios = 70% hospitalisation in Paton et al, 30% in other settings
hosp_name<-c("03","05","07")
hosp_params<-c(0.3,0.5,0.7)
paton_hosp_name<-hosp_name
paton_hosp_params<-hosp_params
pmc_name<-c("no_pmc","pmc")
pmc_params<-c(FALSE,TRUE)
h<-3
ph<-1
p<-1
for(h in 1:length(hosp_name)) {
   for(p in 1:length(pmc_name)) {
     if(h==1) { 
       ph<-3
     } else if(h==2) {
       ph<-2
     } else {
       ph<-1
     }
     for(i in 1:nrow(imperial_mod)) {
       run_repeat_mod(inc_sma_genpop = eval(parse(text=paste0("imperial_mod$inc_sma_ppy_genpop_prophosp_paton",paton_hosp_name[ph],"[i]"))),
                      eir=imperial_mod$eir[i],
                      pmc=pmc_params[p],
                      output="rds",
                      ft_um_pd0 = 0.5,
                      ft_sev_pd0 = 0.5,
                      prop_hosp=hosp_params[h] ,
                      prop_hosp_paton=paton_hosp_params[ph] ,
                      r_eir=r_eir,
                      prob_um_plac=prob_um_plac,
                      prob_um_pmc=prob_um_pmc,
                      max_beta_phi=max_beta_phi,
                      w_scale_risk=w_scale_risk,
                      w_shape_risk=w_shape_risk,
                      adm=paste0(pmc_name[p],"_prophosp",hosp_name[h],"_", "_prophosppaton",paton_hosp_name[ph],"_",i))
     }
  }
}

grp$status()
grp$wait(timeout=1200)

## save results
for(h in 1:length(hosp_name)) {
  for(p in 1:length(pmc_name)) {
    if(h==1) { 
      ph<-3
    } else if(h==2) {
      ph<-2
    } else {
      ph<-1
    }
    imperial_mod[[paste0("repmod_",pmc_name[p],"_h",hosp_name[h],"_", "_ph",paton_hosp_name[ph], "_inc_sma_ppy")]]<-NA
      imperial_mod[[paste0("repmod_",pmc_name[p],"_h",hosp_name[h],"_", "_ph",paton_hosp_name[ph], "_inc_sma_all_paton_eq_check")]]<-NA
      imperial_mod[[paste0("repmod_",pmc_name[p],"_h",hosp_name[h],"_", "_ph",paton_hosp_name[ph], "_inc_sma_all_hosp_eq")]]<-NA
      imperial_mod[[paste0("repmod_",pmc_name[p],"_h",hosp_name[h],"_", "_ph",paton_hosp_name[ph], "_inc_sma_pd_eq")]]<-NA
      imperial_mod[[paste0("repmod_",pmc_name[p],"_h",hosp_name[h],"_", "_ph",paton_hosp_name[ph], "_inc_sev_oth_pd_eq")]]<-NA
      imperial_mod[[paste0("repmod_",pmc_name[p],"_h",hosp_name[h],"_", "_ph",paton_hosp_name[ph], "_inc_sma_1_3_pd_eq")]]<-NA
      imperial_mod[[paste0("repmod_",pmc_name[p],"_h",hosp_name[h],"_", "_ph",paton_hosp_name[ph], "_inc_sev_oth_1_3_pd_eq")]]<-NA
      imperial_mod[[paste0("repmod_",pmc_name[p],"_h",hosp_name[h],"_", "_ph",paton_hosp_name[ph], "_num_pmc_eq")]]<-NA
      imperial_mod[[paste0("repmod_",pmc_name[p],"_h",hosp_name[h],"_", "_ph",paton_hosp_name[ph], "_num_sma_gen_pppy")]]<-NA
      imperial_mod[[paste0("repmod_",pmc_name[p],"_h",hosp_name[h],"_", "_ph",paton_hosp_name[ph], "_num_sma_pd_pppy")]]<-NA
      imperial_mod[[paste0("repmod_",pmc_name[p],"_h",hosp_name[h],"_", "_ph",paton_hosp_name[ph], "_num_sev_oth_pd_pppy")]]<-NA
      imperial_mod[[paste0("repmod_",pmc_name[p],"_h",hosp_name[h],"_", "_ph",paton_hosp_name[ph], "_deaths_pppy")]]<-NA
  }
}


                
for(h in 1:length(hosp_name)) {
  for(p in 1:length(pmc_name)) {
    if(h==1) { 
      ph<-3
    } else if(h==2) {
      ph<-2
    } else {
      ph<-1
    }
    for(i in 1:nrow(imperial_mod)) {
        curr_file<-readRDS(paste0("out_summ_adm",pmc_name[p],"_prophosp",hosp_name[h],"_", "_prophosppaton",paton_hosp_name[ph],"_",i,".rds"))
        vars<-names(curr_file)
        imperial_mod[[paste0("repmod_",pmc_name[p],"_h",hosp_name[h],"_", "_ph",paton_hosp_name[ph], "_inc_sma_ppy")]][i]<-as.numeric(curr_file[which(vars=="inc_sma_all_eq")])
        imperial_mod[[paste0("repmod_",pmc_name[p],"_h",hosp_name[h],"_", "_ph",paton_hosp_name[ph], "_inc_sma_all_paton_eq_check")]][i]<-as.numeric(curr_file[which(vars=="inc_sma_all_paton_eq")])
        imperial_mod[[paste0("repmod_",pmc_name[p],"_h",hosp_name[h],"_", "_ph",paton_hosp_name[ph], "_inc_sma_all_hosp_eq")]][i]<-as.numeric(curr_file[which(vars=="inc_sma_all_hosp_eq")])
        imperial_mod[[paste0("repmod_",pmc_name[p],"_h",hosp_name[h],"_", "_ph",paton_hosp_name[ph], "_inc_sma_pd_eq")]][i]<-as.numeric(curr_file[which(vars=="inc_sma_pd_eq")])
        imperial_mod[[paste0("repmod_",pmc_name[p],"_h",hosp_name[h],"_", "_ph",paton_hosp_name[ph], "_inc_sev_oth_pd_eq")]][i]<-as.numeric(curr_file[which(vars=="inc_sev_oth_pd_eq")])
        imperial_mod[[paste0("repmod_",pmc_name[p],"_h",hosp_name[h],"_", "_ph",paton_hosp_name[ph], "_inc_sma_1_3_pd_eq")]][i]<-as.numeric(curr_file[which(vars=="inc_sma_1_3_pd_eq")])
        imperial_mod[[paste0("repmod_",pmc_name[p],"_h",hosp_name[h],"_", "_ph",paton_hosp_name[ph], "_inc_sev_oth_1_3_pd_eq")]][i]<-as.numeric(curr_file[which(vars=="inc_sev_oth_1_3_pd_eq")])
        imperial_mod[[paste0("repmod_",pmc_name[p],"_h",hosp_name[h],"_", "_ph",paton_hosp_name[ph], "_num_pmc_eq")]][i]<-as.numeric(curr_file[which(vars=="num_pmc_eq")])
        imperial_mod[[paste0("repmod_",pmc_name[p],"_h",hosp_name[h],"_", "_ph",paton_hosp_name[ph], "_num_sma_gen_pppy")]][i]<-as.numeric(curr_file[which(vars=="num_sma_gen_pppy")])
        imperial_mod[[paste0("repmod_",pmc_name[p],"_h",hosp_name[h],"_", "_ph",paton_hosp_name[ph], "_num_sma_pd_pppy")]][i]<-as.numeric(curr_file[which(vars=="num_sma_pd_pppy")])
        imperial_mod[[paste0("repmod_",pmc_name[p],"_h",hosp_name[h],"_", "_ph",paton_hosp_name[ph], "_num_sev_oth_pd_pppy")]][i]<-as.numeric(curr_file[which(vars=="num_sev_oth_pd_pppy")])
        imperial_mod[[paste0("repmod_",pmc_name[p],"_h",hosp_name[h],"_", "_ph",paton_hosp_name[ph], "_deaths_pppy")]][i]<-as.numeric(curr_file[which(vars=="deaths_pppy")])
      }
    #}
  }
}

write.csv(imperial_mod,"adm1_epi_2019_sma2_paton_cal_repmod3.csv",row.names = F)

