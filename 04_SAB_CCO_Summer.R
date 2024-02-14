
####################################################
### Conditional Logistic Regression ~ SAB Analysis
### Restricted to summer events

#Written by: Emma Gause
#Date: 06/20/23
#updated: 08/28/23

#load in libraries
library(dplyr)
library(tidyverse)
library(tidyr)
library(lubridate)
library(survival)

#create path to directory
datadir <- "[/path to your data directory/]"
dataexp <- "[/path to your export directory/]"

#read in SAB PRESTO data
sab <- readRDS(paste0(dataexp, "CONUS_cohort_082423.rds"))

#read in temp data 
temp <- readRDS(paste0(dataexp, "AllTemps_2013_2022_7dayPctls.rds"))

#### NOTE: this analysis is repeated using absolute temperatures and using a 3-day exposure window. Code not shown. 

##--------------------------------------------------------------##

#remove spatial info from data
sab <- st_drop_geometry(sab)
str(sab)

#divide into summer/extended summer cases
sab <- sab %>% mutate(month = month(SAB_DATE))
sabsum <- sab %>% filter(month>5&month<9) #just June/July/August
sabext <- sab %>% filter(month>4&month<10) #includes May and September (West coast particularly, summer months extend into September)


##--------------------------------------------------------------##

#create control dates!!

# Create data structure for cco
ds <- data.frame(SAB_DATE=as.Date(character()),
                 control_date1=as.Date(character()), control_date2=as.Date(character()),
                 control_date3=as.Date(character()), control_date4=as.Date(character()),
                 control_date5=as.Date(character()), control_date6=as.Date(character()),
                 control_date7=as.Date(character()), control_date8=as.Date(character()),
                 control_date9=as.Date(character()), stringsAsFactors = FALSE)
series <- data.frame(SAB_DATE=seq(as.Date('2013-01-01'), as.Date('2022-12-31'),'day')) #sequence everyday in study years

#this gets the date for every possible same DOW in the month for every date in our sequence
#i.e. every possible combination of control dates for every date
for(i in seq(nrow(series))){
  cat(i, '')
  ds[i, 1] <- series$SAB_DATE[i]
  ds[i, 2] <- series$SAB_DATE[i]-28
  ds[i, 3] <- series$SAB_DATE[i]-21
  ds[i, 4] <- series$SAB_DATE[i]-14
  ds[i, 5] <- series$SAB_DATE[i]-7
  ds[i, 6] <- series$SAB_DATE[i]
  ds[i, 7] <- series$SAB_DATE[i]+7
  ds[i, 8] <- series$SAB_DATE[i]+14
  ds[i, 9] <- series$SAB_DATE[i]+21
  ds[i, 10] <- series$SAB_DATE[i]+28
}

#keep only the control dates that are in the same month, and flag the case data ==1
ds <- ds %>% gather('var','control_date', control_date1:control_date9) %>% 
  filter(month(SAB_DATE)==month(control_date)) %>% select(-var) %>% 
  mutate(case=case_when(SAB_DATE==control_date~1, TRUE~0))
#this leaves us with a df with every date from 2013-2021 and every single control date for each


#create an indicator for extreme heat
str(temp)
temp <- temp %>% 
  group_by(StCoFIPS) %>% mutate(Text7day = if_else(Tmax_7day>95, 1, 0),
                                Wext7day = if_else(Wmax_7day>95, 1, 0)) %>% ungroup()

table(temp$Text7day)
table(temp$Wext7day)

#change the temp units so 1 = a 10 percentile difference in temp -- better interpretation
  #In absolute temperature snesitivity analysis we used a 5 degree C difference 
  #both of these are just linear transformations to aid in the interpretation 
temp <- temp %>% mutate(Tmean_7day = Tmean_7day/10,
                        Tmax_7day = Tmax_7day/10,
                        Wmean_7day = Wmean_7day/10,
                        Wmax_7day = Wmax_7day/10)


#merge cases to control dates and to temp data
str(ds)
str(sabext) 
str(temp)
dat <- left_join(sabext, ds, by="SAB_DATE")

#data too large for join -- remove unnecessary columns
rm(ds, sab, series)
str(dat)
dat <- dat %>% select(STUDYID, SAB_DATE, SAB_WKS, LMP_ALL, conceive_date,
                      sab_year, fips, control_date, case)
datt <- left_join(dat, temp, by= c("control_date"="date", "fips"="fips"))


#get stratification of early and late SAB
table(datt$SAB_WKS, useNA = "ifany") #lots of missing!

#calculate estimated SAB weeks
datt <- datt %>% mutate(est_SABwks = as.numeric(round(difftime(SAB_DATE, LMP_ALL, units = "weeks"))))
table(datt$est_SABwks, useNA = "ifany")

#supplement missingness
datt <- datt %>% mutate(SAB_WKS_i = if_else(is.na(SAB_WKS), est_SABwks, SAB_WKS))
table(datt$SAB_WKS_i, useNA = "ifany")

#split into early and late SAB events
datt <- datt %>% mutate(early = if_else(SAB_WKS_i<8, 1, 0))
table(datt$early, datt$case, useNA = "ifany", deparse.level = 2)
early <- datt %>% filter(early==1)
late <- datt %>% filter(early==0)

table(early$case)
table(late$case)
table(datt$case)
##--------------------------------------------------------------##

table(datt$case)

## TEMPERATURE ##

# Temperature: 7-day Average
t7avg <- clogit(case ~ Tmean_7day + strata(STUDYID), data=datt, method = 'breslow') 

# Temperature: 7-day Maximum
t7max <- clogit(case ~ Tmax_7day + strata(STUDYID), data=datt, method = 'breslow') 

# Temperature: 7-day Any Extreme Heat
t7ext <- clogit(case ~ Text7day + strata(STUDYID), data=datt, method = 'breslow') 

#stratified by early and late now
#T 7-day Average
t7avg_e <- clogit(case ~ Tmean_7day + strata(STUDYID), data=early, method = 'breslow') 
t7avg_l <- clogit(case ~ Tmean_7day + strata(STUDYID), data=late, method = 'breslow') 

#T 7-day Maximum
t7max_e <- clogit(case ~ Tmax_7day + strata(STUDYID), data=early, method = 'breslow') 
t7max_l <- clogit(case ~ Tmax_7day + strata(STUDYID), data=late, method = 'breslow') 

#T 7-day Any Extreme Heat
t7ext_e <- clogit(case ~ Text7day + strata(STUDYID), data=early, method = 'breslow') 
t7ext_l <- clogit(case ~ Text7day + strata(STUDYID), data=late, method = 'breslow') 



## WBGT ##

# WBGT: 7-day Average
w7avg <- clogit(case ~ Wmean_7day + strata(STUDYID), data=datt, method = 'breslow') 

# WBGT: 7-day Maximum
w7max <- clogit(case ~ Wmax_7day + strata(STUDYID), data=datt, method = 'breslow') 

# WBGT: 7-day Any Extreme Heat
w7ext <- clogit(case ~ Wext7day + strata(STUDYID), data=datt, method = 'breslow') 

#stratified by early and late now
#W 7-day Average
w7avg_e <- clogit(case ~ Wmean_7day + strata(STUDYID), data=early, method = 'breslow') 
w7avg_l <- clogit(case ~ Wmean_7day + strata(STUDYID), data=late, method = 'breslow') 

#W 7-day Maximum
w7max_e <- clogit(case ~ Wmax_7day + strata(STUDYID), data=early, method = 'breslow') 
w7max_l <- clogit(case ~ Wmax_7day + strata(STUDYID), data=late, method = 'breslow') 

#W 7-day Any Extreme Heat
w7ext_e <- clogit(case ~ Wext7day + strata(STUDYID), data=early, method = 'breslow') 
w7ext_l <- clogit(case ~ Wext7day + strata(STUDYID), data=late, method = 'breslow') 


## GATHER COEFFICIENTS
ls()

#there's probably a recursive way to do this quicker...
t7avg_or <- cbind(exp(coef(t7avg)), exp(confint(t7avg)))
t7max_or <- cbind(exp(coef(t7max)), exp(confint(t7max)))
t7ext_or <- cbind(exp(coef(t7ext)), exp(confint(t7ext)))
w7avg_or <- cbind(exp(coef(w7avg)), exp(confint(w7avg)))
w7max_or <- cbind(exp(coef(w7max)), exp(confint(w7max)))
w7ext_or <- cbind(exp(coef(w7ext)), exp(confint(w7ext)))

t7avg_or_e <- cbind(exp(coef(t7avg_e)), exp(confint(t7avg_e)))
t7avg_or_l <- cbind(exp(coef(t7avg_l)), exp(confint(t7avg_l)))

t7max_or_e <- cbind(exp(coef(t7max_e)), exp(confint(t7max_e)))
t7max_or_l <- cbind(exp(coef(t7max_l)), exp(confint(t7max_l)))

t7ext_or_e <- cbind(exp(coef(t7ext_e)), exp(confint(t7ext_e)))
t7ext_or_l <- cbind(exp(coef(t7ext_l)), exp(confint(t7ext_l)))

w7avg_or_e <- cbind(exp(coef(w7avg_e)), exp(confint(w7avg_e)))
w7avg_or_l <- cbind(exp(coef(w7avg_l)), exp(confint(w7avg_l)))

w7max_or_e <- cbind(exp(coef(w7max_e)), exp(confint(w7max_e)))
w7max_or_l <- cbind(exp(coef(w7max_l)), exp(confint(w7max_l)))

w7ext_or_e <- cbind(exp(coef(w7ext_e)), exp(confint(w7ext_e)))
w7ext_or_l <- cbind(exp(coef(w7ext_l)), exp(confint(w7ext_l)))


extsum_res <- rbind(t7avg_or, t7max_or, t7ext_or, 
                    w7avg_or, w7max_or, w7ext_or,
                    t7avg_or_e, t7avg_or_l, 
                    t7max_or_e, t7max_or_l,
                    t7ext_or_e, t7ext_or_l,
                    w7avg_or_e, w7avg_or_l,
                    w7max_or_e, w7max_or_l,
                    w7ext_or_e, w7ext_or_l)


write.csv(extsum_res, 
          "/[directory path]/SAB_ExtSummer_CCO_Coefs_082823.csv",
          row.names = TRUE)
