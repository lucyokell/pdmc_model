########### TO GET SEVERE MALARIA INCIDENCE, EIR, PREV for and trial areas
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

hosp<-read.csv("C:/Users/lokell/Dropbox (SPH Imperial College)/PMC post discharge/hospital_locations.csv")
shp0 <- readOGR(dsn="C:/Users/lokell/Dropbox (SPH Imperial College)/PMC post discharge/shape_files", layer="SSA_shapefile_adm0")

#Load the input site-specific parameters -----------------------------------------------------
#source: malaria drive
#W:/Global_site_file_update/PROJECT_20180511_GADM/final_fitted_and_key_inputs/
## copied to local folder:

fulldata<-read.csv("GADM_site_file_fitted.csv")  
source("reformat_site_adm_data.R")  ## put model params into easier format to use.
##### Makes data_test (fully subsetted to africa and single draw, reformatted) and data_draw0 (subsetted to Africa and single draw)



##################################################################################################
##  SET UP CLUSTER
##################################################################################################

#setwd("Q:/Lucy/Programming/MalariaLaunchR")
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

context::context_log_start()
root <- "context" #"context_workers"


#options(didehpc.cluster = "fi--didemrchnb")
options(didehpc.cluster = "fi--dideclusthn")
didehpc::web_login()

# shares 
#share <- didehpc::path_mapping("malaria", "M:", "//fi--didef3/malaria", "M:")
share <- didehpc::path_mapping("malaria", "Z:", "//fi--didef3/malaria", "Z:")

# load malria launchR function. Hayley specific function so should not need it?
#sources <- c("s3_malaria_function.R")

# create a context  
ctx <- context::context_save(root,
                             #sources=sources,
                             package_sources=provisionr::package_sources(repos="file:///Z:/drat"),
                             packages=c("MalariaLaunchR"))

#config the cluster 
config <- didehpc::didehpc_config(shares=share,
                                  #cluster="fi--dideclusthn", 
                                  cluster="fi--didemrchnb", 
                                  use_workers = FALSE)

# create a queue - interface to the cluster queue 
obj <- didehpc::queue_didehpc(ctx, config=config)


obj$enqueue(sessionInfo())

#################################
# TEST RUNS
#################################
# local test
Model_launcher(OutputName = "Run1", OutputRoot = "//fi--san03/homes/lokell/Lucy/PMC_post_discharge/output",
               Parms = "Default_parms2.txt")
# cluster test
### NB to change output vars, default parms files etc for use on the cluster:
# need to add them to 'working directory'\context\lib\windows\4.0\MalariaLaunchR\model_files
t1<-obj$enqueue(Model_launcher(OutputName = "Run2",OutputRoot = "//fi--san03/homes/lokell/Lucy/PMC_post_discharge/output",
                               Parms = "Default_parms2.txt", Options = "num_people 1000")) 
t1$status()
t1$log()

names(read.table("//fi--san03/homes/lokell/Lucy/PMC_post_discharge/output/Run2.txt",header=T))




##############################################
# RUN MALARIA LAUNCHR IN TRIAL SITES, EXTRACT LOCAL EIR, PREV ETC AT ADMIN LEVEL
##############################################

##### Phiri 2012 ###########################
malawi_admins<-c("Blantyre","Chikwawa","Thyolo","Zomba")
malawi_inds<-which(data_draw0$NAME_0=="Malawi" & data_draw0$NAME_1 %in% malawi_admins)
malawi_dide_codes<-data_draw0$DIDE_CODE[malawi_inds]
# dates = June 2006 - August 2009 for RECRUITMENT
# dates in model launcher years:
# start_trial<-2006+5/12 - 2017
# end_trial<-2009+8/12 +3/12 - 2017
res_phiri2012_eir<-res_phiri2012_epi<-data_draw0 %>%
  dplyr::select(DIDE_CODE,CONTINENT, ISO, NAME_0, NAME_1) %>%
  dplyr::filter(DIDE_CODE %in% malawi_dide_codes)

### do a test run to check size of output and set up results file
adm<-malawi_inds[1]
subset<- data_test[c(1,1+adm)]     
names(subset)[2] <- "value"
subset$value[which(subset$par=="num_people")]<-1000 # do a shorter run
Run1 <- Model_launcher(OutputName = "Run1",
                       Options = subset)
Run1<-Run1$Y %>%
  dplyr::filter(year<0) %>%
  dplyr::select(year, EIR)

# make extra cols to store EIR
res_phiri2012_eir<-cbind(res_phiri2012_eir,data.frame(matrix(ncol = nrow(Run1), nrow = nrow(res_phiri2012_eir))))
eir_var_inds_malawi<-grep("X",names(res_phiri2012_eir))

start_time = Sys.time()
for(i in 1:length(malawi_inds)) {  ##loop through malawi adm1 units
  adm<-malawi_inds[i]
  print(adm)
  subset<- data_test[c(1,1+adm)]     ##subset = each time you keep 1 column as the "value" column
  names(subset)[2] <- "value"
  
  Run1 <- Model_launcher(OutputName = "Run1",
                         Options = subset)
  
  ##for now year is set to 0 (2017)
  #plot(Run1$Y$year,Run1$Y$EIR,type="l")
  res_phiri2012_eir[i,eir_var_inds_malawi] <-Run1$Y %>%
    dplyr::filter(year<0) %>%
    dplyr::pull(EIR)
  
}

write.csv(res_phiri2012_eir, file="C:/Users/lokell/Dropbox (SPH Imperial College)/PMC post discharge/res_phiri2012_eir.csv",row.names = F)


############# Kwambai 2020 ################
# In western Kenya, we will recruit participants from
# hospitals located in areas around Lake Victoria with
# well-documented malaria transmission intensity, including
# the Jaramogi Oginga Odinga Teaching & Referral -> in kisumu city
# Hospital (JOOTRH) and Siaya, Kisumu, Homa Bay
# and Migori County referral hospitals. In Uganda, we
# will recruit from Jinja, Hoima, Masaka and Mubende
# regional referral hospitals as well as Kamuli Mission
# Hospital (Fig. 1).
kwambai_admins<-c("Siaya","Kisumu","Homa Bay","Migori")  # in Kenya
kwambai_inds1<-which(data_draw0$NAME_0=="Kenya" & data_draw0$NAME_1 %in% kwambai_admins)
kwambai_admins_ug<-c("Jinja","Hoima","Masaka","Mubende","Kamuli")  # in Uganda
kwambai_inds2<-which(data_draw0$NAME_0=="Uganda" & data_draw0$NAME_1 %in% kwambai_admins_ug)
kwambai_inds<-c(kwambai_inds1,kwambai_inds2)
kwambai_dide_codes<-data_draw0$DIDE_CODE[kwambai_inds]


# dates = May 2016 - May 2018 for RECRUITMENT
# dates in model launcher years:
start_trial<-2016+4/12 - 2017
end_trial<-2018+5/12 +3/12 - 2017
end_trial<-0  ## end trial at time 0 (2017) as do not have intervetion coverage afterwards at present.
res_kwambai_eir<-res_kwambai_epi<-data_draw0 %>%
  dplyr::select(DIDE_CODE,CONTINENT, ISO, NAME_0, NAME_1) %>%
  dplyr::filter(DIDE_CODE %in% kwambai_dide_codes)
res_kwambai_epi<-res_kwambai_epi %>%
  mutate(sev_inc_0_5=NA,
         sev_inc_0_10=NA,
         clin_inc_0_5=NA,
         prev_2_10=NA,
         prev_0_5=NA)

### do a test run to check size of output and set up results file
adm<-kwambai_inds[1]
subset<- data_test[c(1,1+adm)]     
names(subset)[2] <- "value"
subset$value[which(subset$par=="num_people")]<-10000 # do a shorter run
Run1 <- Model_launcher(OutputName = "Run1",
                       Options = subset, 
                       Parms = "Default_parms2.txt")
test<-Run1$Y %>%
  dplyr::filter(year<0) %>%
  dplyr::select(year, EIR)

###### get length
length<-Run1$Y %>%
  dplyr::filter(year>start_trial & year<0) 

# make extra cols to store EIR
res_kwambai_eir<-cbind(res_kwambai_eir,data.frame(matrix(ncol = nrow(length), nrow = length(kwambai_inds))))
eir_var_inds_kwambai<-grep("X",names(res_kwambai_eir))

start_time = Sys.time()
for(i in 1:length(kwambai_inds)) {  ##loop through adm1 units
  adm<-kwambai_inds[i]
  print(adm)
  subset<- data_test[c(1,1+adm)]     ##subset = each time you keep 1 column as the "value" column
  names(subset)[2] <- "value"
  
  # Run1 <- Model_launcher(OutputName = "Run1",
  #                        Options = subset,
  #                         Parms = "Default_parms2.txt")
  # 
  # write.csv(Run1$Y, file=paste0("C:/Users/lokell/Dropbox (SPH Imperial College)/PMC post discharge/trial_output/full_output_adm",
  #     adm,".csv"),row.names = F)
  Run1$Y<-read.csv(paste0("C:/Users/lokell/Dropbox (SPH Imperial College)/PMC post discharge/trial_output/full_output_adm",
                          adm,".csv"))
  
  ##for now year is set to 0 (2017)
  plot(Run1$Y$year,Run1$Y$EIR,type="l")
  res_kwambai_eir[i,eir_var_inds_kwambai] <-Run1$Y %>%
    dplyr::filter(year<0) %>%
    dplyr::pull(EIR)
  
  trial<-Run1$Y %>%
    dplyr::filter(year>start_trial & year<0) 
  res_kwambai_epi$sev_inc_0_5[i]<-mean(trial$sev_inc_0_5)
  res_kwambai_epi$sev_inc_0_10[i]<-mean(trial$sev_inc_0_10)
  res_kwambai_epi$clin_inc_0_5[i]<-mean(trial$clin_inc_0_5)
  res_kwambai_epi$prev_2_10[i]<-mean(trial$prev_2_10)
  res_kwambai_epi$prev_0_5[i]<-mean(trial$prev_0_5)
  
}
end_time<-Sys.time()
end_time-start_time

write.csv(res_kwambai_eir, file="C:/Users/lokell/Dropbox (SPH Imperial College)/PMC post discharge/res_kwambai_eir.csv",row.names = F)
write.csv(res_kwambai_epi, file="C:/Users/lokell/Dropbox (SPH Imperial College)/PMC post discharge/res_kwambai_epi.csv",row.names = F)


##############################
# RECALIBRATE TO HOSPITAL SPECIFIC LOCATION
# NOW EXTRACT TRANSMISSION FROM MAP 20 KM AROUND HOSPITAL, RECALIBRATE MODEL AND EXTRACT EIR
##############################
LL<-data.frame(x=hosp$lon,y=hosp$lat)
hosp$map_prev_trial<-NA
for(i in 1:nrow(LL)) {
  years<-hosp$start_year_trial[i]:hosp$end_year_trial[i]
  prevs<-rep(NA,length(years))
  for(j in 1:length(years)) {
    raster_filename<-
      paste0("C:/Users/lokell/Dropbox (SPH Imperial College)/MAP/2020_Global_PfPR_2000/2020_GBD2019_Global_PfPR_",years[j],".tif") 
    temp<-raster(raster_filename)   
    prevs[j]<-raster::extract(temp,LL[i,],buffer=20000,fun=mean)
  }
  hosp$map_prev_trial[i]<-mean(prevs)
  print(i)
}

## visualise hospital locations on prevalence map
# plot(temp,xlim=c(-20,60),ylim=c(-40,40))
plot(shp0,xlim=c(-20,60),ylim=c(-40,40),axes=T)
points(hosp$lon,hosp$lat,col="red",pch=19)
plot(temp,xlim=c(30,40),ylim=c(-5,5),axes=F)
points(hosp$lon,hosp$lat,col="red")
text(hosp$lon,hosp$lat,hosp$short_name,cex=0.5,pos=2)



##################################################
# Refit mosquito density as per Malaria Launch R vignette "Fitting"
##################################################
# first make local site files. Number is now DIDE CODE to avoid mistakes.
for(i in 1:nrow(hosp)) {  ##loop through adm1 units
  adm<-hosp$DIDE_CODE[i]
  print(adm)
  subset<- data_draw0 %>% 
    dplyr::filter(DIDE_CODE==adm) %>%    ##subset = each time you keep 1 column as the "value" column
    dplyr::select(-DIDE_CODE:-NAME_1)
  subset<-data.frame(par=names(subset),value=as.numeric(subset))
  write.table(subset,file=paste0("C:/Users/lokell/Documents/R/win-library/4.0/MalariaLaunchR/model_files/sites/site_",adm,".txt"),
              sep="\t",col.names = F,row.names = F,quote=F)
}

#### fit to the mean prevalence over the trial years and the mean year.
# Fit_M not happy when fitting to years before year 0 for some reason. 
# for(i in 1:nrow(hosp)) {  ##loop through adm1 units
#   adm<-hosp$DIDE_CODE[i]
#   print(adm)
#   year_match<-round(mean(hosp$start_year_trial[i]:hosp$end_year_trial[i]))-2017
#   Fit_M(variable = 'prev_2_10', target=hosp$map_prev_trial[i], 
#                   year=year_match,
#                   Options="init_output 17",
#                   Root="C:/Users/lokell/Documents/R/win-library/4.0/MalariaLaunchR/model_files", 
#                   Site=paste0("site_",adm,".txt"), Options=NULL,
#                   overwrite = TRUE,
#                   tolerance = 0.01, interval = c(0.1, 1000), maxiter = 50)
# }

# Fit with optimise function instead.
fit_m_optim<-function(m,target) {
  Run1 <- Model_launcher(OutputName = "Run1",
                         Options = paste0("output_type 0 final_run 1 init_output 17 num_people 10000 total_M ",m),
                         Site=paste0("site_",adm,".txt"))
  model_prev<-mean(dplyr::filter(Run1$Y, year>=hosp$start_year_trial[i]-2017 & year <=min(0,hosp$end_year_trial[i]-2017)) %>%
                     dplyr::pull(prev_2_10))
  return((model_prev-target)^2)
}

for(i in 1:nrow(hosp)) {
  adm<-hosp$DIDE_CODE[i]
  print(adm)
  x<-optimise(fit_m_optim,interval=c(0.1,100),target=hosp$map_prev_trial[i],tol=0.005)
  m_fitted<-x$minimum
  ### rewrite site file
  site_file<-read.delim(paste0("C:/Users/lokell/Documents/R/win-library/4.0/MalariaLaunchR/model_files/sites/site_",adm,".txt"),
                        header=F)
  print("original total M")
  print(site_file$V2[which(site_file$V1=="total_M")])
  site_file<-site_file %>% dplyr::mutate(V2=ifelse(par=="total_M",m_fitted,V2))
  print("fitted total M")
  print(site_file$V2[which(site_file$V1=="total_M")])
  write.table(site_file,file=paste0("C:/Users/lokell/Documents/R/win-library/4.0/MalariaLaunchR/model_files/sites/site_refitted_",hosp$study[i],"_",adm,".txt"),
              sep="\t",col.names = F,row.names = F,quote=F)
}

####################################################
##### Rerun and check prevalence is right. Save EIR.
####################################################
res_trial_eir<-res_trial_epi<-data.frame(study=hosp$study, DIDE_CODE=hosp$DIDE_CODE)
res_trial_epi<-res_trial_epi %>%
  mutate(eir=NA,
         sev_inc_0_5=NA,
         sev_inc_0_10=NA,
         clin_inc_0_5=NA,
         prev_2_10=NA,
         prev_0_5=NA)

### do a test run to check size of output and set up results file
subset<- data_test[c(1,1+1)]     
names(subset)[2] <- "value"
subset$value[which(subset$par=="num_people")]<-1000 # do a shorter run
Run1 <- Model_launcher(OutputName = "Run1",
                       Options = subset, 
                       Parms = "Default_parms2.txt")
###### get length of output
length<-Run1$Y %>%
  dplyr::filter(year<0)

# make extra cols to store EIR. From 2000-2017
res_trial_eir<-cbind(res_trial_eir,data.frame(matrix(ncol = nrow(length), nrow = nrow(res_trial_eir))))
eir_var_inds_trial<-grep("X",names(res_trial_eir))

########## RUN FOR EACH TRIAL SITE using recalibrated site files
#### remember to copy refitted site files to local 'context' folder if running on the cluster.
for(i in 1:nrow(hosp)) {  ##loop through adm1 units
  adm<-hosp$DIDE_CODE[i]
  study<-hosp$study[i]
  print(adm)
  t1<-obj$enqueue(Model_launcher(OutputName = paste0("full_output_recalibrateMAP_", study,"_",  adm),
                                 OutputRoot = "//fi--san03/homes/lokell/Lucy/PMC_post_discharge/trial_output",
                                 Parms = "Default_parms2.txt",
                                 Options = "output_type 0 final_run 1 output_per_yr 365 init_output 17 num_people 100000",
                                 Site=paste0("site_refitted_",study,"_",adm,".txt")))

  #write.csv(Run1$Y, file=paste0("C:/Users/lokell/Dropbox (SPH Imperial College)/PMC post discharge/trial_output/,".csv"),row.names = F)
}
t1$status()

# Loop through results and summarise.
for(i in 1:nrow(hosp)) {  ##loop through adm1 units
  adm<-hosp$DIDE_CODE[i]
  study<-hosp$study[i]
  print(adm)
  Run1<-read.table(paste0("trial_output/full_output_recalibrateMAP_",
                           hosp$study[i],"_",  adm,".txt"),header=T)
  # # 
  ##for now year is set to 0 (2017)
  plot(Run1$year,Run1$EIR,type="l")
  res_trial_eir[i,eir_var_inds_trial] <-Run1 %>%
    dplyr::filter(year<0) %>%
    dplyr::pull(EIR)
  
  trial<-Run1 %>%
    dplyr::filter(year>hosp$start_year_trial[i]-2017 & year<min(0,hosp$end_year_trial[i]-2017)) 
  res_trial_epi$eir[i]<-mean(trial$EIR)
  res_trial_epi$sev_inc_0_5[i]<-mean(trial$sev_inc_0_5)
  res_trial_epi$sev_inc_0_10[i]<-mean(trial$sev_inc_0_10)
  res_trial_epi$clin_inc_0_5[i]<-mean(trial$clin_inc_0_5)
  res_trial_epi$prev_2_10[i]<-mean(trial$prev_2_10)
  res_trial_epi$prev_0_5[i]<-mean(trial$prev_0_5)
  
}

write.csv(res_trial_eir, file="C:/Users/lokell/Dropbox (SPH Imperial College)/PMC post discharge/res_trial_eir2.csv",row.names = F)
write.csv(res_trial_epi, file="C:/Users/lokell/Dropbox (SPH Imperial College)/PMC post discharge/res_trial_epi2.csv",row.names = F)

hosp<-read.csv("C:/Users/lokell/Dropbox (SPH Imperial College)/PMC post discharge/hospital_locations.csv")
hosp<-full_join(hosp,res_trial_epi,by=c("study","DIDE_CODE"))
LL<-data.frame(x=hosp$lon,y=hosp$lat)
hosp$map_prev_trial<-NA
for(i in 1:nrow(LL)) {
  years<-hosp$start_year_trial[i]:hosp$end_year_trial[i]
  prevs<-rep(NA,length(years))
  for(j in 1:length(years)) {
    raster_filename<-
      paste0("C:/Users/lokell/Dropbox (SPH Imperial College)/MAP/2020_Global_PfPR_2000/2020_GBD2019_Global_PfPR_",years[j],".tif") 
    temp<-raster(raster_filename)   
    prevs[j]<-raster::extract(temp,LL[i,],buffer=20000,fun=mean)
  }
  hosp$map_prev_trial[i]<-mean(prevs)
  print(i)
}
write.csv(hosp,file="C:/Users/lokell/Dropbox (SPH Imperial College)/PMC post discharge/trial_hospital_epi.csv",row.names = F)
