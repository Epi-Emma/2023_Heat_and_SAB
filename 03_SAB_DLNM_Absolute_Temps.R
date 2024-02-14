
########################################################
### Analysis of temp and wbgt on SAB date vs. controls
### DLNM Analysis using all-year events and full temp spectrum
### ABSOLUTE TEMPERATURE SENSITIVITY ANALYSIS 

#Written by: Emma Gause
#Date: 02/07/23
#updated: 11/08/23

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

#remove expected missingness for coastal counties: 25019 and 12087
str(temp)
temp <- temp %>% filter(StCoFIPS!="25019" & StCoFIPS!="12087")
summary(temp$Tmax_C)

# get percentiles both for absolute max temp AND max wbgt (*100 to convert to decimal to 0-100 format)
temp <- temp %>% select(StCoFIPS, date=Date, Tmax_C, WBGTmax_C) %>%
  group_by(StCoFIPS) %>%
  mutate(fips = as.character(StCoFIPS)) %>%
  arrange(StCoFIPS, date) %>% 
  mutate(lag1_temp=lag(Tmax_C, 1), lag2_temp=lag(Tmax_C, 2),
         lag3_temp=lag(Tmax_C, 3), lag4_temp=lag(Tmax_C, 4), 
         lag5_temp=lag(Tmax_C, 5), lag6_temp=lag(Tmax_C, 6),
         lag1_wbgt=lag(WBGTmax_C, 1), lag2_wbgt=lag(WBGTmax_C, 2),
         lag3_wbgt=lag(WBGTmax_C, 3), lag4_wbgt=lag(WBGTmax_C, 4), 
         lag5_wbgt=lag(WBGTmax_C, 5), lag6_wbgt=lag(WBGTmax_C, 6)) %>% ungroup()


#now calculate moving Means/maximums across 3-days and 7-days 
temp <- temp %>% group_by(StCoFIPS) %>%
  arrange(StCoFIPS, date) %>% 
  mutate(Tmean_3day = rollapply(Tmax_C,3,mean,align='right',fill=NA),
         Tmax_3day = rollapply(Tmax_C,3,max,align='right',fill=NA),
         Tmean_7day = rollapply(Tmax_C,7,mean,align='right',fill=NA),
         Tmax_7day = rollapply(Tmax_C,7,max,align='right',fill=NA),
         Wmean_3day = rollapply(WBGTmax_C,3,mean,align='right',fill=NA),
         Wmax_3day = rollapply(WBGTmax_C,3,max,align='right',fill=NA),
         Wmean_7day = rollapply(WBGTmax_C,7,mean,align='right',fill=NA),
         Wmax_7day = rollapply(WBGTmax_C,7,max,align='right',fill=NA)) %>% ungroup()


#saveRDS(temp, paste0(dataexp, "AllTemps_2013_2022_7dayAbs.rds"))


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


summary(datt)
test <- datt %>% filter(is.na(Tmax_C))

##--------------------------------------------------------------##

#TEMP MODEL - 7 day moving average
cen7t <- median(datt$Tmean_7day)
p107t <- as.numeric(quantile(datt$Tmean_7day, probs = 0.10))
p907t <- as.numeric(quantile(datt$Tmean_7day, probs = 0.90))
p17t <- as.numeric(quantile(datt$Tmean_7day, probs = 0.01))
p997t <- as.numeric(quantile(datt$Tmean_7day, probs = 0.99))

cbtemp7 <- crossbasis(datt[,27], 
                      lag=0, 
                      argvar=list(knots=c(p107t, cen7t, p907t)))

modt7 <- clogit(case ~ cbtemp7 + strata(STUDYID), data=datt)
summary(modt7)

cbtemppredt7 <- crosspred(cbtemp7, modt7, cen=cen7t) 

#in ggplot
perc <- as.data.frame(as.numeric(row.names(cbtemppredt7$matRRfit)))
est <- cbtemppredt7$matRRfit
lb <- cbtemppredt7$matRRlow
ub <- cbtemppredt7$matRRhigh
df <- cbind(perc, est, lb, ub)
str(df)
colnames(df) <- c("perc", "est", "lb", "ub")


tiff(filename = "/[directory path]/SAB_ABS_7day_avg_temp.tiff",
     width = 1400, height = 1100, units = "px", res = 300)

tempa7a <- ggplot(df) +
  geom_hline(yintercept= 1, linetype = 'dashed', lwd = 0.9) +
  geom_line(aes(x = perc, y = est), color = 'red') +
  geom_ribbon(aes(x = perc, y = est, ymin = lb, ymax = ub), fill = "grey75", alpha=0.4) + 
  scale_y_continuous(trans = 'log10',
                     breaks = scales::log_breaks(n = 6),
                     limits = c(0.7, 4),
                     name = "OR (95% CI)") +
  scale_x_continuous(name = "7-day mean daily maximum temperature (Celsius)",
                     limits = c(p17t, p997t)) + 
  theme_gdocs()
tempa7a
dev.off()



#TEMP MODEL - 7 day moving MAX
cen7tx <- median(datt$Tmax_7day)
p107tx <- as.numeric(quantile(datt$Tmax_7day, probs = 0.10))
p907tx <- as.numeric(quantile(datt$Tmax_7day, probs = 0.90))
p17tx <- as.numeric(quantile(datt$Tmax_7day, probs = 0.01))
p997tx <- as.numeric(quantile(datt$Tmax_7day, probs = 0.99))

cbtemp7x <- crossbasis(datt[,28], 
                       lag=0, 
                       argvar=list(knots=c(p107tx, cen7tx, p907tx)))

modt7x <- clogit(case ~ cbtemp7x + strata(STUDYID), data=datt)
summary(modt7x)

cbtemppredt7x <- crosspred(cbtemp7x, modt7x, cen=cen7tx) 


#in ggplot
perc <- as.data.frame(as.numeric(row.names(cbtemppredt7x$matRRfit)))
est <- cbtemppredt7x$matRRfit
lb <- cbtemppredt7x$matRRlow
ub <- cbtemppredt7x$matRRhigh
df <- cbind(perc, est, lb, ub)
str(df)
colnames(df) <- c("perc", "est", "lb", "ub")


tiff(filename = "/[directory path]/SAB_ABS_7day_max_temp.tiff",
     width = 1400, height = 1100, units = "px", res = 300)

tempa7m <- ggplot(df) +
  geom_hline(yintercept= 1, linetype = 'dashed', lwd = 0.9) +
  geom_line(aes(x = perc, y = est), color = 'red') +
  geom_ribbon(aes(x = perc, y = est, ymin = lb, ymax = ub), fill = "grey75", alpha=0.4) + 
  scale_y_continuous(trans = 'log10',
                     breaks = scales::log_breaks(n = 6),
                     limits = c(0.7, 4),
                     name = "OR (95% CI)") +
  scale_x_continuous(name = "7-day maximum temperature (Celsius)",
                     limits = c(p17tx, p997tx)) + 
  theme_gdocs()
tempa7m
dev.off()



#### SENSITIVITY ANALYSIS USING 2 KNOTS
#TEMP MODEL - 7 day moving average
cen7t_2k <- median(datt$Tmean_7day)
p337t_2k <- as.numeric(quantile(datt$Tmean_7day, probs = 0.33))
p667t_2k <- as.numeric(quantile(datt$Tmean_7day, probs = 0.66))
p17t_2k <- as.numeric(quantile(datt$Tmean_7day, probs = 0.01))
p997t_2k <- as.numeric(quantile(datt$Tmean_7day, probs = 0.99))

cbtemp7_2k <- crossbasis(datt[,27], 
                      lag=0, 
                      argvar=list(knots=c(p337t_2k, p667t_2k)))

modt7_2k <- clogit(case ~ cbtemp7_2k + strata(STUDYID), data=datt)
summary(modt7_2k)

cbtemppredt7_2k <- crosspred(cbtemp7_2k, modt7_2k, cen=cen7t_2k) 

#in ggplot
perc <- as.data.frame(as.numeric(row.names(cbtemppredt7_2k$matRRfit)))
est <- cbtemppredt7_2k$matRRfit
lb <- cbtemppredt7_2k$matRRlow
ub <- cbtemppredt7_2k$matRRhigh
df <- cbind(perc, est, lb, ub)
str(df)
colnames(df) <- c("perc", "est", "lb", "ub")


tiff(filename = "/[directory path]/SAB_ABS_7day_2knot_sens.tiff",
     width = 933, height = 733, units = "px", res = 200)

tempa7a_2k <- ggplot(df) +
  geom_hline(yintercept= 1, linetype = 'dashed', lwd = 0.9) +
  geom_line(aes(x = perc, y = est), color = 'red') +
  geom_ribbon(aes(x = perc, y = est, ymin = lb, ymax = ub), fill = "grey75", alpha=0.4) + 
  scale_y_continuous(trans = 'log10',
                     breaks = scales::log_breaks(n = 6),
                     limits = c(0.7, 4),
                     name = "OR (95% CI)") +
  scale_x_continuous(name = "7-day mean daily maximum temperature (Celsius)",
                     limits = c(p17t_2k, p997t_2k)) + 
  theme_gdocs()
tempa7a_2k
dev.off()

##--------------------------------------------------------------##
##--------------------------------------------------------------##

#WBGT MODEL - 7 day moving average
cen7w <- median(datt$Wmean_7day)
p107w <- as.numeric(quantile(datt$Wmean_7day, probs = 0.10))
p907w <- as.numeric(quantile(datt$Wmean_7day, probs = 0.90))
p17w <- as.numeric(quantile(datt$Wmean_7day, probs = 0.01))
p997w <- as.numeric(quantile(datt$Wmean_7day, probs = 0.99))

cbwbgt7 <- crossbasis(datt[,31], 
                      lag=0, 
                      argvar=list(knots=c(p107w, cen7w, p907w)))

modt7 <- clogit(case ~ cbwbgt7 + strata(STUDYID), data=datt)
summary(modt7)

cbwbgtpredt7 <- crosspred(cbwbgt7, modt7, cen=cen7w) 

#in ggplot
perc <- as.data.frame(as.numeric(row.names(cbwbgtpredt7$matRRfit)))
est <- cbwbgtpredt7$matRRfit
lb <- cbwbgtpredt7$matRRlow
ub <- cbwbgtpredt7$matRRhigh
df <- cbind(perc, est, lb, ub)
str(df)
colnames(df) <- c("perc", "est", "lb", "ub")


tiff(filename = "/[directory path]/AB_ABS_7day_avg_wbgt.tiff",
     width = 933, height = 733, units = "px", res = 200)

wbgta7a <- ggplot(df) +
  geom_hline(yintercept= 1, linetype = 'dashed', lwd = 0.9) +
  geom_line(aes(x = perc, y = est), color = 'red') +
  geom_ribbon(aes(x = perc, y = est, ymin = lb, ymax = ub), fill = "grey75", alpha=0.4) + 
  scale_y_continuous(trans = 'log10',
                     breaks = scales::log_breaks(n = 6),
                     limits = c(0.62, 3),
                     name = "OR (95% CI)") +
  scale_x_continuous(name = "7-day mean daily maximum WBGT (absolute index)",
                     limits = c(p17w, p997w)) + 
  theme_gdocs()
wbgta7a
dev.off()


#WBGT MODEL - 7 day moving MAX
cen7wx <- median(datt$Wmax_7day)
p107wx <- as.numeric(quantile(datt$Wmax_7day, probs = 0.10))
p907wx <- as.numeric(quantile(datt$Wmax_7day, probs = 0.90))
p17wx <- as.numeric(quantile(datt$Wmax_7day, probs = 0.01))
p997wx <- as.numeric(quantile(datt$Wmax_7day, probs = 0.99))

cbwbgt7x <- crossbasis(datt[,32], 
                       lag=0, 
                       argvar=list(knots=c(p107wx, cen7wx, p907wx)))

modt7x <- clogit(case ~ cbwbgt7x + strata(STUDYID), data=datt)
summary(modt7x)

cbwbgtpredt7x <- crosspred(cbwbgt7x, modt7x, cen=cen7wx) 

#in ggplot
perc <- as.data.frame(as.numeric(row.names(cbwbgtpredt7x$matRRfit)))
est <- cbwbgtpredt7x$matRRfit
lb <- cbwbgtpredt7x$matRRlow
ub <- cbwbgtpredt7x$matRRhigh
df <- cbind(perc, est, lb, ub)
str(df)
colnames(df) <- c("perc", "est", "lb", "ub")


tiff(filename = "/[directory path]/SAB_ABS_7day_max_wbgt.tiff",
     width = 933, height = 733, units = "px", res = 200)

wbgta7m <- ggplot(df) +
  geom_hline(yintercept= 1, linetype = 'dashed', lwd = 0.9) +
  geom_line(aes(x = perc, y = est), color = 'red') +
  geom_ribbon(aes(x = perc, y = est, ymin = lb, ymax = ub), fill = "grey75", alpha=0.4) + 
  scale_y_continuous(trans = 'log10',
                     breaks = scales::log_breaks(n = 6),
                     limits = c(0.62, 3),
                     name = "OR (95% CI)") +
  scale_x_continuous(name = "7-day maximum WBGT (absolute index)",
                     limits = c(p17wx, p997wx)) + 
  theme_gdocs()
wbgta7m
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

#TEMP MODEL - 3 day moving average
colnames(datt)

cen3t <- median(datt$Tmean_3day)
p103t <- as.numeric(quantile(datt$Tmean_3day, probs = 0.10))
p903t <- as.numeric(quantile(datt$Tmean_3day, probs = 0.90))
p13t <- as.numeric(quantile(datt$Tmean_3day, probs = 0.01))
p993t <- as.numeric(quantile(datt$Tmean_3day, probs = 0.99))

cbtemp3 <- crossbasis(datt[,25], 
                      lag=0, 
                      argvar=list(knots=c(p103t, cen3t, p903t)))

modt3 <- clogit(case ~ cbtemp3 + strata(STUDYID), data=datt)
summary(modt3)

cbtemppredt3 <- crosspred(cbtemp3, modt3, cen=cen3t) 

#in ggplot
perc <- as.data.frame(as.numeric(row.names(cbtemppredt3$matRRfit)))
est <- cbtemppredt3$matRRfit
lb <- cbtemppredt3$matRRlow
ub <- cbtemppredt3$matRRhigh
df <- cbind(perc, est, lb, ub)
str(df)
colnames(df) <- c("perc", "est", "lb", "ub")


tiff(filename = "/[directory path]/SAB_ABS_3day_avg_temp.tiff",
     width = 933, height = 733, units = "px", res = 200)

tempa3a <- ggplot(df) +
  geom_hline(yintercept= 1, linetype = 'dashed', lwd = 0.9) +
  geom_line(aes(x = perc, y = est), color = 'red') +
  geom_ribbon(aes(x = perc, y = est, ymin = lb, ymax = ub), fill = "grey75", alpha=0.4) + 
  scale_y_continuous(trans = 'log10',
                     breaks = scales::log_breaks(n = 6),
                     limits = c(0.65, 3),
                     name = "OR (95% CI)") +
  scale_x_continuous(name = "3-day mean daily maximum temperature (Celsius)",
                     limits = c(p13t, p993t)) + 
  theme_gdocs()
tempa3a
dev.off()



#TEMP MODEL - 3 day moving MAX
cen3tx <- median(datt$Tmax_3day)
p103tx <- as.numeric(quantile(datt$Tmax_3day, probs = 0.10))
p903tx <- as.numeric(quantile(datt$Tmax_3day, probs = 0.90))
p13tx <- as.numeric(quantile(datt$Tmax_3day, probs = 0.01))
p993tx <- as.numeric(quantile(datt$Tmax_3day, probs = 0.99))

cbtemp3x <- crossbasis(datt[,26], 
                       lag=0, 
                       argvar=list(knots=c(p103tx, cen3tx, p903tx)))

modt3x <- clogit(case ~ cbtemp3x + strata(STUDYID), data=datt)
summary(modt3x)

cbtemppredt3x <- crosspred(cbtemp3x, modt3x, cen=cen3tx) 

#in ggplot
perc <- as.data.frame(as.numeric(row.names(cbtemppredt3x$matRRfit)))
est <- cbtemppredt3x$matRRfit
lb <- cbtemppredt3x$matRRlow
ub <- cbtemppredt3x$matRRhigh
df <- cbind(perc, est, lb, ub)
str(df)
colnames(df) <- c("perc", "est", "lb", "ub")


tiff(filename = "/[directory path]/SAB_ABS_3day_max_temp.tiff",
     width = 933, height = 733, units = "px", res = 200)

tempa3m <- ggplot(df) +
  geom_hline(yintercept= 1, linetype = 'dashed', lwd = 0.9) +
  geom_line(aes(x = perc, y = est), color = 'red') +
  geom_ribbon(aes(x = perc, y = est, ymin = lb, ymax = ub), fill = "grey75", alpha=0.4) + 
  scale_y_continuous(trans = 'log10',
                     breaks = scales::log_breaks(n = 6),
                     limits = c(0.65, 3),
                     name = "OR (95% CI)") +
  scale_x_continuous(name = "3-day maximum temperature (Celsius)",
                     limits = c(p13tx, p993tx)) + 
  theme_gdocs()
tempa3m
dev.off()


##--------------------------------------------------------------##
##--------------------------------------------------------------##

#WBGT MODEL - 3 day moving average
colnames(datt)
cen3w <- median(datt$Wmean_3day)
p103w <- as.numeric(quantile(datt$Wmean_3day, probs = 0.10))
p903w <- as.numeric(quantile(datt$Wmean_3day, probs = 0.90))
p13w <- as.numeric(quantile(datt$Wmean_3day, probs = 0.01))
p993w <- as.numeric(quantile(datt$Wmean_3day, probs = 0.99))

cbwbgt3 <- crossbasis(datt[,29], 
                      lag=0, 
                      argvar=list(knots=c(p103w, cen3w, p903w)))

modt3 <- clogit(case ~ cbwbgt3 + strata(STUDYID), data=datt)
summary(modt3)

cbwbgtpredt3 <- crosspred(cbwbgt3, modt3, cen=cen3w) 

#in ggplot
perc <- as.data.frame(as.numeric(row.names(cbwbgtpredt3$matRRfit)))
est <- cbwbgtpredt3$matRRfit
lb <- cbwbgtpredt3$matRRlow
ub <- cbwbgtpredt3$matRRhigh
df <- cbind(perc, est, lb, ub)
str(df)
colnames(df) <- c("perc", "est", "lb", "ub")


tiff(filename = "/[directory path]/SAB_ABS_3day_avg_wbgt.tiff",
     width = 933, height = 733, units = "px", res = 200)

wbgta3a <- ggplot(df) +
  geom_hline(yintercept= 1, linetype = 'dashed', lwd = 0.9) +
  geom_line(aes(x = perc, y = est), color = 'red') +
  geom_ribbon(aes(x = perc, y = est, ymin = lb, ymax = ub), fill = "grey75", alpha=0.4) + 
  scale_y_continuous(trans = 'log10',
                     breaks = scales::log_breaks(n = 6),
                     limits = c(0.6, 3),
                     name = "OR (95% CI)") +
  scale_x_continuous(name = "3-day mean daily maximum WBGT (absolute index)",
                     limits = c(p13w, p993w)) + 
  theme_gdocs()
wbgta3a
dev.off()


#WBGT MODEL - 3 day moving MAX
cen3wx <- median(datt$Wmax_3day)
p103wx <- as.numeric(quantile(datt$Wmax_3day, probs = 0.10))
p903wx <- as.numeric(quantile(datt$Wmax_3day, probs = 0.90))
p13wx <- as.numeric(quantile(datt$Wmax_3day, probs = 0.01))
p993wx <- as.numeric(quantile(datt$Wmax_3day, probs = 0.99))

cbwbgt3x <- crossbasis(datt[,30], 
                       lag=0, 
                       argvar=list(knots=c(p103wx, cen3wx, p903wx)))

modt3x <- clogit(case ~ cbwbgt3x + strata(STUDYID), data=datt)
summary(modt3x)

cbwbgtpredt3x <- crosspred(cbwbgt3x, modt3x, cen=cen3wx) 

#in ggplot
perc <- as.data.frame(as.numeric(row.names(cbwbgtpredt3x$matRRfit)))
est <- cbwbgtpredt3x$matRRfit
lb <- cbwbgtpredt3x$matRRlow
ub <- cbwbgtpredt3x$matRRhigh
df <- cbind(perc, est, lb, ub)
str(df)
colnames(df) <- c("perc", "est", "lb", "ub")


tiff(filename = "/[directory path]/SAB_ABS_3day_max_wbgt.tiff",
     width = 933, height = 733, units = "px", res = 200)

wbgta3m <- ggplot(df) +
  geom_hline(yintercept= 1, linetype = 'dashed', lwd = 0.9) +
  geom_line(aes(x = perc, y = est), color = 'red') +
  geom_ribbon(aes(x = perc, y = est, ymin = lb, ymax = ub), fill = "grey75", alpha=0.4) + 
  scale_y_continuous(trans = 'log10',
                     breaks = scales::log_breaks(n = 6),
                     limits = c(0.6, 3),
                     name = "OR (95% CI)") +
  scale_x_continuous(name = "3-day maximum WBGT (absolute index)",
                     limits = c(p13wx, p993wx)) + 
  theme_gdocs()
wbgta3m
dev.off()

