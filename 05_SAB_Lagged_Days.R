
########################################################
### Analysis of temp and wbgt on SAB date vs. controls
### DLNM analysis using temperature percentiles
### Individual lagged days sensitivity analysis

#Written by: Emma Gause
#Date: 02/07/23
#updated: 01/08/24

#load in libraries
library(dplyr)
library(tidyverse)
library(tidyr)
library(lubridate)
library(sf)
library(data.table)
library(dlnm)
library(splines)
library(survival)
library(zoo)
library(ggthemes)

datadir <- "[/path to your data directory/]"
dataexp <- "[/path to your export directory/]"

#read in SAB PRESTO data
sab <- readRDS(paste0(dataexp, "CONUS_cohort_052323.rds"))

#remove spatial info from data
sab <- st_drop_geometry(sab)
str(sab)

#read in temp data 
temp <- readRDS(paste0(dataexp, "AllTemps_2013_2022_7dayAvg.rds"))

##NOTE: this analysis is repeated using absolute temperatures. Code not shown. 

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
#this leaves us with a df with every date from 2013-2022 and every single control date for each


##--------------------------------------------------------------##

#merge cases to control dates and to temp data
str(ds)
str(sab) 
str(temp)
dat <- left_join(sab, ds, by="SAB_DATE")

#data too large for join -- remove unnecessary columns
rm(ds, sab, series)
str(dat)
dat <- dat %>% select(STUDYID, SAB_DATE, SAB_WKS, LMP_ALL, conceive_date,
                      sab_year, fips, control_date, case)
datt <- left_join(dat, temp, by= c("control_date"="date", "fips"="fips"))

##--------------------------------------------------------------##

summary(datt) 
test <- datt %>% filter(is.na(Tmax_C))

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

#TEMP MODEL - 7 day temp percentiles, daily lags
str(datt)

cbtemp7 <- crossbasis(datt[,c(13,15:20)], 
                      lag=6, 
                      argvar=list(knots=c(10,50, 90)),
                      arglag=list(knots=1))

modt7 <- clogit(case ~ cbtemp7 + strata(STUDYID), data=datt)
#summary(modt7)

cbtemppredt7 <- crosspred(cbtemp7, modt7, cen=50, from = 0, to = 100, by = 5) 
#summary(cbtemppredt7)

tiff(filename = "/[directory path]/SAB_7lag_Prctls.tiff",
     width = 833, height = 633, units = "px", res = 120)
plot(cbtemppredt7,
     ptype = "slices", log = "y",
     var=95, 
     ci="b", type="p", 
     col=4, pch=19, cex=1.7,
     ylab="OR", xlab="Day Lag", main=NA,
     ylim = c(0.7, 1.6))

box(col = "black", which = "figure")

dev.off()







