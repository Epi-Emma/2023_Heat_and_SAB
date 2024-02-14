
################################################################
### Analysis of temp and wbgt on SAB date vs. controls
### DLNM Analysis using all-year events and full temp spectrum

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

#read in temperature data 
temp <- readRDS(paste0(dataexp, "temperature/Heatmetrics_CONUS_County_2012_2022.rds"))

##--------------------------------------------------------------##

#prepare temperature data 

#remove expected missingness for coastal counties: 25019 and 12087
str(temp)
temp <- temp %>% filter(StCoFIPS!="25019" & StCoFIPS!="12087")
summary(temp$Tmax_C)

#convert absolute temperatures into temperature percentiles, unique to county
#and create lagged days up to 6

# get percentiles both for absolute max temp AND max wbgt (*100 to convert to decimal to 0-100 format)
temp <- temp %>% select(StCoFIPS, date=Date, Tmax_C, WBGTmax_C) %>%
  group_by(StCoFIPS) %>%
  mutate(fips = as.character(StCoFIPS),
         tempmax_pctl=percent_rank(Tmax_C)*100,
         wbgtmax_pctl=percent_rank(WBGTmax_C)*100) %>%
  arrange(StCoFIPS, date) %>% 
  mutate(lag1_temp=lag(tempmax_pctl, 1), lag2_temp=lag(tempmax_pctl, 2),
         lag3_temp=lag(tempmax_pctl, 3), lag4_temp=lag(tempmax_pctl, 4), 
         lag5_temp=lag(tempmax_pctl, 5), lag6_temp=lag(tempmax_pctl, 6),
         lag1_wbgt=lag(wbgtmax_pctl, 1), lag2_wbgt=lag(wbgtmax_pctl, 2),
         lag3_wbgt=lag(wbgtmax_pctl, 3), lag4_wbgt=lag(wbgtmax_pctl, 4), 
         lag5_wbgt=lag(wbgtmax_pctl, 5), lag6_wbgt=lag(wbgtmax_pctl, 6)) %>% ungroup()

#now calculate moving Means/maximums across 3-days and 7-days 
temp <- temp %>% group_by(StCoFIPS) %>%
  arrange(StCoFIPS, date) %>% 
  mutate(Tmean_3day = rollapply(tempmax_pctl,3,mean,align='right',fill=NA),
         Tmax_3day = rollapply(tempmax_pctl,3,max,align='right',fill=NA),
         Tmean_7day = rollapply(tempmax_pctl,7,mean,align='right',fill=NA),
         Tmax_7day = rollapply(tempmax_pctl,7,max,align='right',fill=NA),
         Wmean_3day = rollapply(wbgtmax_pctl,3,mean,align='right',fill=NA),
         Wmax_3day = rollapply(wbgtmax_pctl,3,max,align='right',fill=NA),
         Wmean_7day = rollapply(wbgtmax_pctl,7,mean,align='right',fill=NA),
         Wmax_7day = rollapply(wbgtmax_pctl,7,max,align='right',fill=NA)) %>% ungroup()

saveRDS(temp, paste0(dataexp, "AllTemps_2013_2022_7dayPctls.rds"))

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
test <- datt %>% filter(is.na(Tmax_C)) # no missingness

#get stratification of early and late SAB
table(datt$SAB_WKS, useNA = "ifany") #lots of missing!

#calculate estimated SAB weeks
datt <- datt %>% mutate(est_SABwks = as.numeric(round(difftime(SAB_DATE, LMP_ALL, units = "weeks"))))
table(datt$est_SABwks, useNA = "ifany")

#supplement missingness in actual weeks with estimated weeks
datt <- datt %>% mutate(SAB_WKS_i = if_else(is.na(SAB_WKS), est_SABwks, SAB_WKS))
table(datt$SAB_WKS_i, useNA = "ifany")

#create strata variable for early and late SAB (less then 8 weeks vs. 8+ weeks)
datt <- datt %>% mutate(early = if_else(SAB_WKS_i<8, 1, 0))
table(datt$early, datt$case, useNA = "ifany", deparse.level = 2)
early <- datt %>% filter(early==1)
late <- datt %>% filter(early==0)

table(early$case)
table(late$case)
table(datt$case)
##--------------------------------------------------------------##

#TEMP MODEL - 7 day moving mean
cbtemp7 <- crossbasis(datt[,29], 
                      lag=0, 
                      argvar=list(knots=c(10,50,90)))

modt7 <- clogit(case ~ cbtemp7 + strata(STUDYID), data=datt)
summary(modt7)

cbtemppredt7 <- crosspred(cbtemp7, modt7, cen=50) 

#jpeg(filename = "/[directory path]/SAB_7day_avg_temp.jpg",
#     width = 480, height = 380, units = "px")
plot(cbtemppredt7, "overall", xlab="7-Day Mean Temperature (%ile)", ylab="OR", col=3,
     main=NA, ylim=c(0.5,2), log="y")
#dev.off()

#repeat in ggplot - looks prettier!
perc <- as.data.frame(as.numeric(row.names(cbtemppredt7$matRRfit)))
est <- cbtemppredt7$matRRfit
lb <- cbtemppredt7$matRRlow
ub <- cbtemppredt7$matRRhigh
df <- cbind(perc, est, lb, ub)
str(df)
colnames(df) <- c("perc", "est", "lb", "ub")


tiff(filename = "/[directory path]/SAB_7day_avg_temp.tiff",
     width = 1400, height = 1100, units = "px", res = 300)
     
tempp7a <- ggplot(df) +
  geom_hline(yintercept= 1, linetype = 'dashed', lwd = 0.9) +
  geom_line(aes(x = perc, y = est), color = 'red') +
  geom_ribbon(aes(x = perc, y = est, ymin = lb, ymax = ub), fill = "grey75", alpha=0.4) + 
  scale_y_continuous(trans = 'log10',
                     breaks = scales::log_breaks(n = 6),
                     limits = c(0.68, 3.5),
                     name = "OR (95% CI)") +
  scale_x_continuous(name = "7-day mean daily maximum temperature (percentiles)") + 
  theme_gdocs()
tempp7a
dev.off()


#TEMP MODEL - 7 day moving MAX
cbtemp7x <- crossbasis(datt[,30], 
                       lag=0, 
                       argvar=list(knots=c(10,50,90)))

modt7x <- clogit(case ~ cbtemp7x + strata(STUDYID), data=datt)
summary(modt7x)

cbtemppredt7x <- crosspred(cbtemp7x, modt7x, cen=50) 

#jpeg(filename = "/[directory path]/SAB_7day_max_temp.jpg",
#     width = 480, height = 380, units = "px")
plot(cbtemppredt7x, "overall", xlab="7-Day Max Temperature (%ile)", ylab="OR", col=3,
     main=NA, ylim=c(0.5,2.0), log="y")
#dev.off()


#in ggplot
perc <- as.data.frame(as.numeric(row.names(cbtemppredt7x$matRRfit)))
est <- cbtemppredt7x$matRRfit
lb <- cbtemppredt7x$matRRlow
ub <- cbtemppredt7x$matRRhigh
df <- cbind(perc, est, lb, ub)
str(df)
colnames(df) <- c("perc", "est", "lb", "ub")


tiff(filename = "/[directory path]/SAB_7day_max_temp.tiff",
     width = 1400, height = 1100, units = "px", res = 300)

tempp7m <- ggplot(df) +
  geom_hline(yintercept= 1, linetype = 'dashed', lwd = 0.9) +
  geom_line(aes(x = perc, y = est), color = 'red') +
  geom_ribbon(aes(x = perc, y = est, ymin = lb, ymax = ub), fill = "grey75", alpha=0.4) + 
  scale_y_continuous(trans = 'log10',
                     breaks = scales::log_breaks(n = 6),
                     limits = c(0.68, 3.5),
                     name = "OR (95% CI)") +
  scale_x_continuous(name = "7-day maximum temperature (percentiles)") + 
  theme_gdocs()
tempp7m
dev.off()


#### SENSITIVITY ANALYSIS USING DIFFERENT KNOTS
#TEMP MODEL - 7 day moving mean
cbtemp7_2k <- crossbasis(datt[,29], 
                      lag=0, 
                      argvar=list(knots=c(33,66)))

modt7_2k <- clogit(case ~ cbtemp7_2k + strata(STUDYID), data=datt)
summary(modt7_2k)

cbtemppredt7_2k <- crosspred(cbtemp7_2k, modt7_2k, cen=50) 

#jpeg(filename = "/[directory path]/SAB_7day_avg_temp.jpg",
#     width = 480, height = 380, units = "px")
plot(cbtemppredt7_2k, "overall", xlab="7-Day Mean Temperature (%ile)", ylab="OR", col=3,
     main=NA, ylim=c(0.5,2), log="y")
#dev.off()

#in ggplot
perc <- as.data.frame(as.numeric(row.names(cbtemppredt7_2k$matRRfit)))
est <- cbtemppredt7_2k$matRRfit
lb <- cbtemppredt7_2k$matRRlow
ub <- cbtemppredt7_2k$matRRhigh
df <- cbind(perc, est, lb, ub)
str(df)
colnames(df) <- c("perc", "est", "lb", "ub")


tiff(filename = "/[directory path]/SAB_7day_perc_2knot_sens.tiff",
     width = 933, height = 733, units = "px", res = 200)

tempp7a_2k <- ggplot(df) +
  geom_hline(yintercept= 1, linetype = 'dashed', lwd = 0.9) +
  geom_line(aes(x = perc, y = est), color = 'red') +
  geom_ribbon(aes(x = perc, y = est, ymin = lb, ymax = ub), fill = "grey75", alpha=0.4) + 
  scale_y_continuous(trans = 'log10',
                     breaks = scales::log_breaks(n = 6),
                     limits = c(0.68, 3.5),
                     name = "OR (95% CI)") +
  scale_x_continuous(name = "7-day mean daily maximum temperature (percentiles)") + 
  theme_gdocs()
tempp7a_2k
dev.off()


##--------------------------------------------------------------##
##--------------------------------------------------------------##

#WBGT MODEL - 7 day moving mean
cbwbgt7 <- crossbasis(datt[,33], 
                      lag=0, 
                      argvar=list(knots=c(10,50,90)))

modt7 <- clogit(case ~ cbwbgt7 + strata(STUDYID), data=datt)
summary(modt7)

cbwbgtpredt7 <- crosspred(cbwbgt7, modt7, cen=50) 

#jpeg(filename = "/[directory path]/SAB_7day_avg_wbgt.jpg",
#     width = 480, height = 380, units = "px")
plot(cbwbgtpredt7, "overall", xlab="7-Day Mean WBGT (%ile)", ylab="OR", col=3,
     main=NA, ylim=c(0.5,2))
#dev.off()

#in ggplot
perc <- as.data.frame(as.numeric(row.names(cbwbgtpredt7$matRRfit)))
est <- cbwbgtpredt7$matRRfit
lb <- cbwbgtpredt7$matRRlow
ub <- cbwbgtpredt7$matRRhigh
df <- cbind(perc, est, lb, ub)
str(df)
colnames(df) <- c("perc", "est", "lb", "ub")


tiff(filename = "/[directory path]/SAB_7day_avg_wbgt.tiff",
     width = 933, height = 733, units = "px", res = 200)

wbgtpm <- ggplot(df) +
  geom_hline(yintercept= 1, linetype = 'dashed', lwd = 0.9) +
  geom_line(aes(x = perc, y = est), color = 'red') +
  geom_ribbon(aes(x = perc, y = est, ymin = lb, ymax = ub), fill = "grey75", alpha=0.4) + 
  scale_y_continuous(trans = 'log10',
                     breaks = scales::log_breaks(n = 6),
                     limits = c(0.65, 4.0),
                     name = "OR (95% CI)") +
  scale_x_continuous(name = "7-day mean daily maximum WBGT (percentiles)") + 
  theme_gdocs()
wbgtpm
dev.off()


############
#WBGT MODEL - 7 day moving MAX
cbwbgt7x <- crossbasis(datt[,34], 
                       lag=0, 
                       argvar=list(knots=c(10,50,90)))

modt7x <- clogit(case ~ cbwbgt7x + strata(STUDYID), data=datt)
summary(modt7x)

cbwbgtpredt7x <- crosspred(cbwbgt7x, modt7x, cen=50) 

#jpeg(filename = "/[directory path]/SAB_7day_max_wbgt.jpg",
#     width = 480, height = 380, units = "px")
plot(cbwbgtpredt7x, "overall", xlab="7-Day Max WBGT (%ile)", ylab="OR", col=3,
     main=NA, ylim=c(0.5,2.0))
#dev.off()


#in ggplot
perc <- as.data.frame(as.numeric(row.names(cbwbgtpredt7x$matRRfit)))
est <- cbwbgtpredt7x$matRRfit
lb <- cbwbgtpredt7x$matRRlow
ub <- cbwbgtpredt7x$matRRhigh
df <- cbind(perc, est, lb, ub)
str(df)
colnames(df) <- c("perc", "est", "lb", "ub")


tiff(filename = "/[directory path]/SAB_7day_max_wbgt.tiff",
     width = 933, height = 733, units = "px", res = 200)

wbgtp7m <- ggplot(df) +
  geom_hline(yintercept= 1, linetype = 'dashed', lwd = 0.9) +
  geom_line(aes(x = perc, y = est), color = 'red') +
  geom_ribbon(aes(x = perc, y = est, ymin = lb, ymax = ub), fill = "grey75", alpha=0.4) + 
  scale_y_continuous(trans = 'log10',
                     breaks = scales::log_breaks(n = 6),
                     limits = c(0.65, 4.0),
                     name = "OR (95% CI)") +
  scale_x_continuous(name = "7-day maximum WBGT (percentiles)") + 
  theme_gdocs()
wbgtp7m
dev.off()



##--------------------------------------------------------------##
##--------------------------------------------------------------##

##--------------------------------------------------------------##
##--------------------------------------------------------------##

#### DO THE SAME FOR 3-DAY EXPOSURE WINDOW ####

##--------------------------------------------------------------##
##--------------------------------------------------------------##

##--------------------------------------------------------------##
##--------------------------------------------------------------##

#TEMP MODEL - 3 day moving mean
cbtemp3 <- crossbasis(datt[,27], 
                      lag=0, 
                      argvar=list(knots=c(10,50,90)))

modt3 <- clogit(case ~ cbtemp3 + strata(STUDYID), data=datt)
summary(modt3)

cbtemppredt3 <- crosspred(cbtemp3, modt3, cen=50) 

#in ggplot
perc <- as.data.frame(as.numeric(row.names(cbtemppredt3$matRRfit)))
est <- cbtemppredt3$matRRfit
lb <- cbtemppredt3$matRRlow
ub <- cbtemppredt3$matRRhigh
df <- cbind(perc, est, lb, ub)
str(df)
colnames(df) <- c("perc", "est", "lb", "ub")


tiff(filename = "/[directory path]/SAB_3day_avg_temp.tiff",
     width = 933, height = 733, units = "px", res = 200)

tempp3a <- ggplot(df) +
  geom_hline(yintercept= 1, linetype = 'dashed', lwd = 0.9) +
  geom_line(aes(x = perc, y = est), color = 'red') +
  geom_ribbon(aes(x = perc, y = est, ymin = lb, ymax = ub), fill = "grey75", alpha=0.4) + 
  scale_y_continuous(trans = 'log10',
                     breaks = scales::log_breaks(n = 6),
                     limits = c(0.7, 2.0),
                     name = "OR (95% CI)") +
  scale_x_continuous(name = "3-day mean daily maximum temperature (percentiles)") + 
  theme_gdocs()
tempp3a
dev.off()


#TEMP MODEL - 3 day moving MAX
cbtemp3x <- crossbasis(datt[,28], 
                       lag=0, 
                       argvar=list(knots=c(10,50,90)))

modt3x <- clogit(case ~ cbtemp3x + strata(STUDYID), data=datt)
summary(modt3x)

cbtemppredt3x <- crosspred(cbtemp3x, modt3x, cen=50) 

#in ggplot
perc <- as.data.frame(as.numeric(row.names(cbtemppredt3x$matRRfit)))
est <- cbtemppredt3x$matRRfit
lb <- cbtemppredt3x$matRRlow
ub <- cbtemppredt3x$matRRhigh
df <- cbind(perc, est, lb, ub)
str(df)
colnames(df) <- c("perc", "est", "lb", "ub")


tiff(filename = "/[directory path]/SAB_3day_max_temp.tiff",
     width = 933, height = 733, units = "px", res = 200)

tempp3m <- ggplot(df) +
  geom_hline(yintercept= 1, linetype = 'dashed', lwd = 0.9) +
  geom_line(aes(x = perc, y = est), color = 'red') +
  geom_ribbon(aes(x = perc, y = est, ymin = lb, ymax = ub), fill = "grey75", alpha=0.4) + 
  scale_y_continuous(trans = 'log10',
                     breaks = scales::log_breaks(n = 6),
                     limits = c(0.70, 2.0),
                     name = "OR (95% CI)") +
  scale_x_continuous(name = "3-day maximum temperature (percentiles)") + 
  theme_gdocs()
tempp3m
dev.off()


##--------------------------------------------------------------##
##--------------------------------------------------------------##

#WBGT MODEL - 3 day moving mean
cbwbgt3 <- crossbasis(datt[,31], 
                      lag=0, 
                      argvar=list(knots=c(10,50,90)))

modt3 <- clogit(case ~ cbwbgt3 + strata(STUDYID), data=datt)
summary(modt3)

cbwbgtpredt3 <- crosspred(cbwbgt3, modt3, cen=50) 

#in ggplot
perc <- as.data.frame(as.numeric(row.names(cbwbgtpredt3$matRRfit)))
est <- cbwbgtpredt3$matRRfit
lb <- cbwbgtpredt3$matRRlow
ub <- cbwbgtpredt3$matRRhigh
df <- cbind(perc, est, lb, ub)
str(df)
colnames(df) <- c("perc", "est", "lb", "ub")


tiff(filename = "/[directory path]/SAB_3day_avg_wbgt.tiff",
     width = 933, height = 733, units = "px", res = 200)

wbgtp3a <- ggplot(df) +
  geom_hline(yintercept= 1, linetype = 'dashed', lwd = 0.9) +
  geom_line(aes(x = perc, y = est), color = 'red') +
  geom_ribbon(aes(x = perc, y = est, ymin = lb, ymax = ub), fill = "grey75", alpha=0.4) + 
  scale_y_continuous(trans = 'log10',
                     breaks = scales::log_breaks(n = 6),
                     limits = c(0.65, 2.0),
                     name = "OR (95% CI)") +
  scale_x_continuous(name = "3-day mean daily maximum temperature (percentiles)") + 
  theme_gdocs()
wbgtp3a
dev.off()





#WBGT MODEL - 3 day moving MAX
cbwbgt3x <- crossbasis(datt[,32], 
                       lag=0, 
                       argvar=list(knots=c(10,50,90)))

modt3x <- clogit(case ~ cbwbgt3x + strata(STUDYID), data=datt)
summary(modt3x)

cbwbgtpredt3x <- crosspred(cbwbgt3x, modt3x, cen=50) 


#in ggplot
perc <- as.data.frame(as.numeric(row.names(cbwbgtpredt3x$matRRfit)))
est <- cbwbgtpredt3x$matRRfit
lb <- cbwbgtpredt3x$matRRlow
ub <- cbwbgtpredt3x$matRRhigh
df <- cbind(perc, est, lb, ub)
str(df)
colnames(df) <- c("perc", "est", "lb", "ub")


tiff(filename = "/[directory path]/SAB_3day_max_wbgt.tiff",
     width = 933, height = 733, units = "px", res = 200)

wbgtp3m <- ggplot(df) +
  geom_hline(yintercept= 1, linetype = 'dashed', lwd = 0.9) +
  geom_line(aes(x = perc, y = est), color = 'red') +
  geom_ribbon(aes(x = perc, y = est, ymin = lb, ymax = ub), fill = "grey75", alpha=0.4) + 
  scale_y_continuous(trans = 'log10',
                     breaks = scales::log_breaks(n = 6),
                     limits = c(0.65, 2.0),
                     name = "OR (95% CI)") +
  scale_x_continuous(name = "3-day maximum WBGT (percentiles)") + 
  theme_gdocs()
wbgtp3m
dev.off()

