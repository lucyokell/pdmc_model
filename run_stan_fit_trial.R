################## 
# RUN THE POST DISCHARGE COHORT MODEL AND FIT TO PMC TRIAL DATA
##################

## Packages
library(rstan)
library(dplyr)

## Set to working directory / project.
setwd("C:/Users/lokell/OneDrive - Imperial College London/Paper drafts/PMC/submission/code")

################ READ IN DATA
trial_epi<-read.csv("trial_hospital_epi_git.csv")
dat3d<-readRDS("daily_kwambai_3D.RDS")
n_fup<-read.csv("n_fup_time2.csv")

## cut the last 7 days
dat3d<-dat3d[1:(dim(dat3d)[1]-7),,]

### MODEL PARAMETERS
runTime<-26*7 - 7   # in days
time<-0:runTime

# Timing of PMC post discharge, days
pmc_times<-c(2*7,6*7,10*7)

# prob of leaving AL protection by day since treatment
shape_al<-93.5
scale_al<-0.139
prob_rp<-pgamma(1:runTime, shape = shape_al, scale =scale_al,lower.tail = T) - 
  pgamma(0:(runTime-1), shape = shape_al, scale =scale_al,lower.tail = T)

# prob of DP protection by day since treatment
# Weibull distribution
w_scale2<-28.1
w_slope2<-4.4
dp_eff<-0.99
p_protect_dp<- dp_eff*exp(-((time/w_scale2)^w_slope2))


################ PROCESS TRIAL DATA
site<-1:9   ## sites fitted to.

## N at each fup time
## re-order to be the same as the trial epi order
site_names<-tolower(gsub(" ","",trial_epi$label[site]))
n_fup_plac<-dplyr::select(n_fup, paste0(site_names,"_Placebo"))
n_fup_pmc<-dplyr::select(n_fup, paste0(site_names,"_DP"))

## loss to follow up
prop_fup_plac<-as.matrix(n_fup_plac[2:162,]) ## day 2 onwards.
## calculate relative proportion remaining each day, as proportion of those there the day before.
for(i in 1:nrow(prop_fup_plac)) prop_fup_plac[i,]<-
  as.numeric(n_fup_plac[i+2,]/n_fup_plac[i+1,])

# do not model loss to follow up prior to day 14
filler<-matrix(rep(1,length(site)*14),nrow=14,ncol=length(site))
prop_fup_plac<-rbind(filler,prop_fup_plac)

prop_fup_pmc<-as.matrix(n_fup_pmc[2:162,]) ## day 2 onwards.
for(i in 1:nrow(prop_fup_pmc)) prop_fup_pmc[i,]<-
  as.numeric(n_fup_pmc[i+2,]/n_fup_pmc[i+1,])
prop_fup_pmc<-rbind(filler,prop_fup_pmc)


#### prior EIR
priorEIR_y<-trial_epi$eir[site]
# mean of prior is priorEIR_y (from MAP)
CV<-0.1  ## coefficient of variation
shape_EIR_y=1/(CV^2)
rate_EIR_y<-shape_EIR_y / priorEIR_y

## Prior CI's for table:
for(i in 1:length(rate_EIR_y)) print(round(qgamma(p=c(0.5,0.025,0.975),shape=shape_EIR_y,rate=rate_EIR_y[i]),digits=1))

## Number of EIRs for predictions in posterior plots
nPred<-100


####### set-up the stan model######
model14g<-stan_model("stan_code_fit_trial.stan")

## set up data & parameter list to pass to stan.
data_stan<-list(nSites=length(site),
          runTime=runTime,
          N_plac_init=as.numeric(n_fup_plac[1,]),
          N_pmc_init=as.numeric(n_fup_pmc[1,]),
          prop_fup_plac=prop_fup_plac,
          prop_fup_pmc=prop_fup_pmc,
          daily_events=dat3d,
          shape_EIR_y=shape_EIR_y,
          rate_EIR_y=rate_EIR_y,
          prob_rp=prob_rp,
          last_pmc_time= -170,
          pmc_times=pmc_times,
          length_pp_dp=length(p_protect_dp),
          p_protect_dp=p_protect_dp,
          nPred=nPred,
          EIRpred=exp(log(seq(0.1,30,length.out = nPred))))


#### Optional: ask stan to run on multiple cores
options(mc.cores = 4)

#### Run stan model. This will take a number of hours for 17500 iterations.
fit<-sampling(model14g, 
              data=data_stan,
              iter=200, #17500,  ## 200 iter takes 1 minute.
              warmup= 100, #5000,
              chains=1, #4,
              control = list(adapt_delta=0.95))
params<-rstan::extract(fit)
saveRDS(fit,file="model_fit_stan_output_14g.rds")


##################
## generated quantities
##################
# Run separately to avoid memory issues

fit<-readRDS("model_fit_stan_output_14g.rds")
model14g_gqs<-stan_model("stan_code_fit_trial_gqs.stan")
test<-as.matrix(fit)
rm(fit)  # to avoid memory overload

## set no loss to follow up to make calculations of py easier.
prop_fup_pred<-rep(1,nrow(prop_fup_plac))

data_stan_pred<-list(nSites=length(site),
                runTime=runTime,
                N_plac_init=as.numeric(n_fup_plac[1,]),
                N_pmc_init=as.numeric(n_fup_pmc[1,]),
                prop_fup_plac=prop_fup_plac,
                prop_fup_pmc=prop_fup_pmc,
                daily_events=dat3d,
                shape_EIR_y=shape_EIR_y,
                rate_EIR_y=rate_EIR_y,
                prob_rp=prob_rp,
                last_pmc_time= -170,
                pmc_times=pmc_times,
                length_pp_dp=length(p_protect_dp),
                p_protect_dp=p_protect_dp,
                nPred=100,
                EIRpred=exp(log(seq(0.1,30,length.out = 100))),
                prop_fup_pred=prop_fup_pred
)

## Sample posterior fits
f2 <- gqs(model14g_gqs, data=data_stan_pred, draws = test[sample(1:nrow(test),
                                                                size=100, replace = F),])
saveRDS(f2,file="posterior_fit_gqs.rds")


