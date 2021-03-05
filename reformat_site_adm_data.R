###################################
# site file data preparation and reformatting
###################################


#### subset to draw==0 to select point estimate only (median). 
### Also subsetting to Africa only.

data_draw0 <- fulldata %>% 
  filter(draw==0) %>%
  filter(CONTINENT=="Africa")


#drop columns that aren't parameter inputs (5 first columns)
# remove a few variables that aren't being recognised as input in the model (in object "drop")
drop <- c("draw","ppr", "ppu", "prop_gamb_ss")  ## ppr and ppu are proportion of treatment in the public(ppu)/private(ppr) sector
data_draw0 <- data_draw0[,!(names(data_draw0) %in% drop)]
data <- dplyr::select(data_draw0, -DIDE_CODE:-NAME_1)

########## REFORMAT DATA TO BE MODEL INPUT
####another way : trying with a dataframe of two columns (par and value)
### need to transpose columns and rows and reformat as required
matrix<-data.matrix(data)
t_matrix<-t(matrix)
data_test<-as.data.frame(t_matrix)
par <- rownames(data_test)
rownames(data_test) <- NULL
data_test <- cbind(par,data_test)
ncol(data_test) ##605 admin1 units in total
data_test<-rbind(data_test, c("output_type",rep(0,ncol(data_test)-1)))
data_test<-rbind(data_test, c("final_run",rep(1,ncol(data_test)-1)))
data_test<-rbind(data_test, c("output_per_yr",rep(365,ncol(data_test)-1)))
data_test<-rbind(data_test, c("init_output",rep(17,ncol(data_test)-1)))
data_test<-rbind(data_test, c("num_people",rep(100000,ncol(data_test)-1)))
