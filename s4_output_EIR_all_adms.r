########### TO GET SEVERE MALARIA INCIDENCE, EIR, PREV for admin units (and trial areas)
#### NOW CAN USE CLUSTER FOR THIS
write2file<-F

#### To check current malariaLaunchR files, parameters etc go to: C:/Users/lokell/Documents/R/win-library/4.0/MalariaLaunchR/model_files

# Install the DRAT package------------------------------------------------------------------
#install.packages("drat")

#First you must add a new repository, here called malaria. 
# Set the W: drive to your mapping of the malaria drive (projects.dide.ic.ac.uk/malaria)
#drat::addRepo("malaria", "file:///w:/drat")

# Now you should be able to install the package directly from the new repo:
#install.packages("MalariaLaunchR")

# Have a look at the vignettes --------------------------------------------------------------

# vignette("Introduction", package = "MalariaLaunchR")
# vignette("Flexible", package = "MalariaLaunchR")
#vignette("Fitting", package = "MalariaLaunchR")
# ?Model_launcher  ### For information on arguments 

# Load the package---------------------------------------------------------------------------
library(MalariaLaunchR)
# load other packages needed:
library(magrittr)
library(dplyr)
library(raster)
library(rgdal)

##########################################################################
# READ IN DATA AND INITIAL PROCESSING
##########################################################################
#Load the input site-specific parameters -----------------------------------------------------
#source: malaria drive
#W:/Global_site_file_update/PROJECT_20180511_GADM/final_fitted_and_key_inputs/
## copied to local folder:

fulldata<-read.csv("C:/Users/lokell/Dropbox (SPH Imperial College)/Programming/malariaLaunchR examples/GADM_site_file_fitted.csv")  
source("reformat_site_adm_data.R")


################################################################
# OBTAIN EIR FOR ALL ADMINS 2000-2017 daily.
################################################################
# Run locally as running on cluster involves saving all output files - too large.
##### Make results output
res_eir_2017<-res_epi_2017<-dplyr::select(data_draw0, DIDE_CODE,CONTINENT, ISO, NAME_0, NAME_1)
### set up year timings using a test run
adm<-1
subset<- data_test[c(1,1+adm)]     ##subset = each time you keep 1 column as the "value" column
names(subset)[2] <- "value"
subset$value[which(subset$par=="num_people")]<-1000
Run1 <- Model_launcher(OutputName = "Run1",
                       Options = subset,
                       Parms="Default_parms2.txt")
##for now year is set to 0 (2017)
#plot(Run1$Y$year,Run1$Y$EIR,type="l")
Run1<-Run1$Y %>%
  dplyr::filter(year>= -1 & year<0) %>%
  dplyr::select(year, EIR)
### format is: admin1 details, horizontal EIR format:
# make extra cols to store EIR
res_eir_2017<-cbind(res_eir_2017,data.frame(matrix(ncol = nrow(Run1), nrow = nrow(data_draw0))))
## setup to record sev_inc, clin_inc etc.
res_epi_2017<-res_epi_2017 %>%
  mutate(prev_2_10=NA,
         prev_0_5=NA,
         inf_inc_0_5=NA,
         clin_inc_0_5=NA,
         sev_inc_0_5=NA,
         sev_inc_0_10=NA,
         sev_inc_0_1=NA,
         sev_inc_1_2=NA,
         sev_inc_2_3=NA,
         sev_inc_3_4=NA,
         sev_inc_4_5=NA)

eir_var_inds<-grep("X",names(res_eir_2017))
for(i in 1:nrow(data_draw0)) {  ##loop through 605 adm1 units
  #for(i in 1:2) {  ##loop through 605 adm1 units
  adm<-data_draw0$DIDE_CODE[i]
  print(adm)
  subset<- data_test[c(1,1+i)]     ##subset = each time you keep 1 column as the "value" column
  names(subset)[2] <- "value"
  
  #I could specify year0 and add that to the list
  ###the default is year0=2017 (Pete)-- the site says 2014 but is outdated
  Run1 <- Model_launcher(OutputName = "Run1",
                         OutputRoot = "//fi--san03/homes/lokell/Lucy/PMC_post_discharge/output",
                         Options = subset,
                         Parms = "Default_parms2.txt")
  
  ## on cluster
  # t1<-obj$enqueue(Model_launcher(OutputName = paste0("Run",adm),
  #                                OutputRoot = "//fi--san03/homes/lokell/Lucy/PMC_post_discharge/output",
  #                                Options = subset, 
  #                                Parms = "Default_parms2.txt"))
  # 
  ##for now year is set to 0 (2017)
  #plot(Run1$Y$year,Run1$Y$EIR,type="l")
  Run1_2017<-Run1$Y %>%
    dplyr::filter(year>= -1 & year<0)
  res_eir_2017[i,eir_var_inds] <-Run1_2017$EIR
  res_epi_2017[i,"prev_0_5"]<-mean(Run1_2017$prev_0_5)
  res_epi_2017[i,"prev_2_10"]<-mean(Run1_2017$prev_2_10)
  res_epi_2017[i,"inf_inc_0_5"]<-mean(Run1_2017$inf_inc_0_5)
  res_epi_2017[i,"clin_inc_0_5"]<-mean(Run1_2017$clin_inc_0_5)
  res_epi_2017[i,"sev_inc_0_5"]<-mean(Run1_2017$sev_inc_0_5)
  res_epi_2017[i,"sev_inc_0_10"]<-mean(Run1_2017$sev_inc_0_10)
  res_epi_2017[i,"sev_inc_0_1"]<-mean(Run1_2017$sev_inc_0_1)
  res_epi_2017[i,"sev_inc_1_2"]<-mean(Run1_2017$sev_inc_1_2)
  res_epi_2017[i,"sev_inc_2_3"]<-mean(Run1_2017$sev_inc_2_3)
  res_epi_2017[i,"sev_inc_3_4"]<-mean(Run1_2017$sev_inc_3_4)
  res_epi_2017[i,"sev_inc_4_5"]<-mean(Run1_2017$sev_inc_4_5)
}


if(write2file) write.csv(res_eir_2017, file="C:/Users/lokell/Dropbox (SPH Imperial College)/PMC post discharge/adm1_eir_annual.csv",row.names = F)
if(write2file) write.csv(res_epi_2017, file="C:/Users/lokell/Dropbox (SPH Imperial College)/PMC post discharge/adm1_epi_2017.csv",row.names = F)

