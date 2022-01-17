
library(dplyr)

setwd("C:/Users/lokell/Dropbox (SPH Imperial College)/PMC post discharge")

## read in Imperial College model predictions.
# NB in this file, 'sev_inc' variables = hospitalised malaria
imperial_mod<-read.csv("adm1_epi_2019_sma2.csv")

#SMA - prevalence fit from Paton et al. Science 2021.
alpha<- -13.51
beta<- 7.2
gamma<- 7.94
r<-0.95

sma<-data.frame(prev_2_10=seq(0.0001,0.67,0.001)) #min and max of current MAP admin-1 prev2-10
sma$prev_2_10<-round(sma$prev_2_10,digits = 3) # round exactly for matching later

sma$log_lambda<-alpha + beta/(1+exp(-gamma*sma$prev_2_10))
sma$paton_sma_0_10<-exp(sma$log_lambda)

#########################
# read in SMA inc from Imperial model to get ratio of 0-10 inc to 0-5 inc
#########################
### First smooth over IC model predictions.
# 0-5
imod<-imperial_mod[order(imperial_mod$prev_2_10),]
lmod1<-loess(inc_sma_ppy ~ prev_2_10, data=imod, span=0.80) # 10% smoothing span
imod$sma_smooth_0_5<-predict(lmod1)
sma$smooth_0_5 <- predict(lmod1, newdata = sma)
# set smoothed line>0
sma$smooth_0_5[which(sma$smooth_0_5<0)]<-min(sma$smooth_0_5[which(sma$smooth_0_5>0)])/2
ind<-which(imod$prev_2_10>0.66)
sma$smooth_0_5[which(sma$prev_2_10>0.66)]<-imod$sma_smooth_0_5[ind]
# 0-10
lmod1<-loess(inc_sma_ppy_0_10 ~ prev_2_10, data=imod, span=0.80) # 10% smoothing span
imod$sma_smooth_0_10<-predict(lmod1)
sma$smooth_0_10 <- predict(lmod1, newdata = sma) 
ind<-which(imod$prev_2_10>0.66)
sma$smooth_0_10[which(sma$prev_2_10>0.66)]<-imod$sma_smooth_0_10[ind]

plot(imperial_mod$prev_2_10,imperial_mod$inc_sma_ppy)
points(imperial_mod$prev_2_10,imperial_mod$inc_sma_ppy_0_10,col="blue")
lines(sma$prev_2_10, sma$smooth_0_5,col="red")
lines(sma$prev_2_10, sma$smooth_0_10,col="green")

plot(sma$prev_2_10, sma$smooth_0_5, type="l")
lines(sma$prev_2_10, sma$smooth_0_10)

####### calculate ratio of 0_5 to 0_10
sma$ratio_5_10<-sma$smooth_0_5/sma$smooth_0_10
#assume 0_5 always > 0_10
sma$ratio_5_10[which(sma$ratio_5_10<1)]<-1

#############
# recalculate Paton et al in other groups
#############
sma$paton_sma_0_5<-sma$paton_sma_0_10*sma$ratio_5_10

###############
# Match Imperial mod settings by prevalence and generate Paton estimates.
###############
imperial_mod<-imperial_mod %>%
        dplyr::mutate(prev_2_10_round=round(map_prev_210_2019,digits=3),
                      ind_sma=match(prev_2_10_round,sma$prev_2_10),
                      prev_2_10_sma=sma$prev_2_10[ind_sma],
                      inc_sma_ppy_0_5_paton=sma$paton_sma_0_5[ind_sma] #,
                      #sev_inc_0_5_paton=sma$paton_sev_inc_0_5[ind_sma]
        )

write.csv(imperial_mod,"adm1_epi_2019_sma3_paton.csv")
write.csv(sma,"sma_paton.csv")

###############################
# plot comparison
###############################

tiff(file="paton_sma_age.tiff", width=1700,height=1500,res=300,compression="lzw")
with(sma,
     {
       plot(100*prev_2_10,1000*paton_sma_0_5,type="l",col="blue",  #smooth_0_10,type="l",
            xlab="PfPR2-10 (%)",ylab="SMA incidence per 1000 person years")
       lines(100*prev_2_10, 1000*paton_sma_0_10,col="blue",lty=2)
     }
)
legend("topleft", c("Paton 0-5yr (estimated)","Paton 3mo-9y"),lty=c(1,2),col=c("blue","blue"),bty='n')
dev.off()