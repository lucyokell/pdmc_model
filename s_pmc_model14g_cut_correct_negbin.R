#####################################
# Run stan model to fit trial data
#####################################

##Packages
library(rstan)
library(dplyr)

################ TRIAL DATA
trial_epi<-read.csv("trial_hospital_epi.csv")
dat3d<-readRDS("daily_kwambai_3D.RDS")
## cut the last 7 days
dim(dat3d)
dat3d<-dat3d[1:(dim(dat3d)[1]-7),,]

###Model params
dt<-1
runTime<-26*7 - 7
time<-0:runTime
dur_inc<-12  ##incubation period in days
dur_rec<-10  ##recovery period in days
dur_time_to_treat<-2 ##time from symptoms until patient takes first dose
shape_al<-93.5
scale_al<-0.139
shape_al*scale_al


pmc_times<-c(14,6*7,10*7)

## EIR in PMC age children, adjusting for body surface area.
rho	<-	0.85
a0	<-	2920

prob_treat_today <-1-exp(-(1/dur_time_to_treat)*dt)
# prob of leaving AL protection by day since treatment
prob_rp<-pgamma(1:runTime, shape = shape_al, scale =scale_al,lower.tail = T) - 
  pgamma(0:(runTime-1), shape = shape_al, scale =scale_al,lower.tail = T)
plot(prob_rp)
sum(prob_rp) ## check sums to 1
shape_al*scale_al


# prob of DP protection by day since treatment
# Weibull distribution
w_scale2<-28.1
w_slope2<-4.4
dp_eff<-0.99
p_protect_dp<- dp_eff*exp(-((time/w_scale2)^w_slope2))
plot(time,p_protect_dp)


################ TRIAL DATA
#View(trial_epi)
site<-1:9 ## 

## N at each fup time
n_fup<-read.csv("n_fup_time2.csv")
## re-order to be the same as the trial epi order
site_names<-tolower(gsub(" ","",trial_epi$label[site]))
## Lira = Kamuli
names(n_fup)<-gsub("lira","kamuli", names(n_fup))
n_fup_plac<-dplyr::select(n_fup, paste0(site_names,"_Placebo"))
n_fup_pmc<-dplyr::select(n_fup, paste0(site_names,"_DP"))

prop_fup_plac<-as.matrix(n_fup_plac[2:162,]) ## day 2 onwards.
for(i in 1:nrow(prop_fup_plac)) prop_fup_plac[i,]<-
  as.numeric(n_fup_plac[i+2,]/n_fup_plac[i+1,])

# do not model loss prior to day 14
filler<-matrix(rep(1,length(site)*14),nrow=14,ncol=length(site))
prop_fup_plac<-rbind(filler,prop_fup_plac)

prop_fup_pmc<-as.matrix(n_fup_pmc[2:162,]) ## day 2 onwards.
for(i in 1:nrow(prop_fup_pmc)) prop_fup_pmc[i,]<-
  as.numeric(n_fup_pmc[i+2,]/n_fup_pmc[i+1,])
prop_fup_pmc<-rbind(filler,prop_fup_pmc)

priorEIR_y<-trial_epi$eir[site]

#### prior EIR
# mean of prior is priorEIR_y (from MAP)
CV<-0.1  ## coefficient of variation
shape_EIR_y=1/(CV^2)
rate_EIR_y<-shape_EIR_y / priorEIR_y
## check prior looks sensible
#hist(rgamma(10000,shape=shape_EIR_y,rate=rate_EIR_y[4]))

## Prior CI's for table:
for(i in 1:length(rate_EIR_y)) print(round(qgamma(p=c(0.5,0.025,0.975),shape=shape_EIR_y,rate=rate_EIR_y[i]),digits=1))




####### set-up the stan model######
model14g<-stan_model("pmc_model14g_cut8_correct_negbin.stan")

nPred<-100
data_stan<-list(nSites=length(as.numeric(n_fup_plac[1,])),
          runTime=runTime,
          N_plac_init=as.numeric(n_fup_plac[1,]),
          N_pmc_init=as.numeric(n_fup_pmc[1,]),
          prop_fup_plac=prop_fup_plac,
          prop_fup_pmc=prop_fup_pmc,
          daily_events=dat3d,
          shape_EIR_y=shape_EIR_y,
          rate_EIR_y=rate_EIR_y,
          prob_rp=prob_rp,
          prob_treat_today=prob_treat_today,
          last_pmc_time= -170,
          pmc_times=pmc_times,
          length_pp_dp=length(p_protect_dp),
          p_protect_dp=p_protect_dp,
          nPred=nPred,
          EIRpred=exp(log(seq(0.1,30,length.out = nPred))))
options(mc.cores = 2)


fit<-sampling(model14g, 
              data=data_stan,
              iter=17500,  ## 200 iter takes 1 minute.
              warmup=5000,
              chains=4,
              control = list(adapt_delta=0.95))
params<-rstan::extract(fit)
saveRDS(fit,file="pmc14g_cut8_fit_correct_negbin3.rds")


##################
## generated quantities
##################

fit<-readRDS("pmc14g_cut8_fit_correct_negbin3.rds")
model14g_gqs<-stan_model("pmc_model14g_cut8_gqs_correct_negbin.stan")
test<-as.matrix(fit)
rm(fit)

## set no loss to follow up to make calculations of py easier.
prop_fup_pred<-rep(1,nrow(prop_fup_plac))

data_stan_pred<-list(nSites=length(as.numeric(n_fup_plac[1,])),
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
                prop_fup_pred=prop_fup_pred  ,
                py3_plac= colSums(n_fup_plac[1:(15*7-1-14),])/365
)

f2 <- gqs(model14g_gqs, data=data_stan_pred, draws = test[sample(1:nrow(test),
                                                                size=1000, replace = F),])

saveRDS(f2,file="pmc14g_cut8_fit_gqs_correct_negbin_loglik.rds")


plot(params$EIR_y[,2],type="l")
plot(params$max_beta_phi,type="l")
plot(params$r,type="l")
plot(params$prob_um_plac,type="l")
plot(params$prob_um_pmc,type="l")
plot(params$w_scale_risk ,type="l")
plot(params$w_shape_risk,type="l")
plot(params$k,type="l")
plot(params$lp__,type="l")
summary(params$lp__)

plot(density(params$EIR_y[,2]))
plot(density(params$max_beta_phi))
plot(density(params$k))


colMeans(params$EIR_y)

eir<-1:100
max_beta_phi<-mean(params$max_beta_phi)
r<-mean(params$r)
plot(eir,max_beta_phi*exp(-r*eir))
rw_scale_risk<-mean(params$w_scale_risk)
w_shape_risk<-mean(params$w_shape_risk)
plot(eir,eir*max_beta_phi*exp(-r*eir)*exp(-(((1)/w_scale_risk)^w_shape_risk)))
eir<-30.9
prob_inf_d= 1 - exp(-eir/365);
plot(exp(-((1:182/w_scale_risk)^w_shape_risk)))

test<-vector()
params_gqs<-rstan::extract(f2)
for(i in 1:100) test[i]<-mean(rowSums(params$res2[,,1,i]))
plot(test)


