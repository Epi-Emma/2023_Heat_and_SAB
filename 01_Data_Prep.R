
################################################
### Heat and SAB Analysis - Data Cleaning/Prep

#Written by: Emma Gause
#Date: 01/25/23
#updated: 08/24/23

#load in libraries
library(dplyr)
library(tidyverse)
library(tidyr)
library(haven)
library(lubridate)
library(sf)
library(mapview)

datadir <- "[/path to your data directory/]"
dataexp <- "[/path to your export directory/]"

#read in SAB PRESTO data
sab <- read_sas(paste0(datadir, "heat_sab_001.sas7bdat"))

#read in 2010 county shapefile for merging location information [from US Census TIGER/Line shapefiles]
cnty <- read_sf(paste0(datadir, "counties/tl_2010_us_county10.shp"))


##--------------------------------------------------------------##

str(sab)
summary(sab)

#add year as a variable and get count of events by year
sab$sab_year <- year(sab$SAB_DATE)
table(sab$sab_year)

#we have lat/long points for participants, but we need county FIPS
#make an sf object
geosab = st_as_sf(sab, coords = c("X", "Y"), remove = FALSE)

#mapview(geosab) # some are in AK or HI -- we don't have temp data for these.
#mapview(cnty)

#keep only the CONUS counties (this is where we have temp data) - also remove unnecessary columns
str(cnty)
table(cnty$STATEFP10)
conus <- cnty %>% filter(as.numeric(STATEFP10)<57 & STATEFP10!="02" & STATEFP10!="15") %>%
              select(STATEFP10, COUNTYFP10, "fips" = GEOID10, NAME10, NAMELSAD10)
#mapview(conus) #GOOD!

# set CRS for the points to be the same as county shp
st_crs(conus)
st_crs(geosab)
st_crs(geosab) <- st_crs(conus)
#mapview(conus) + mapview(geosab)

#merge SAB lat/longs to county shape to get FIPS (for temp)
fips <- st_join(geosab, conus, join = st_intersects)
str(fips)
summary(fips)

#remove the participants from HI and AK (no county match)
fips <- fips %>% filter(!is.na(fips)) #18 removed 

##--------------------------------------------------------------##

#there is one ID that occurs twice
sabx <- fips %>% group_by(STUDYID) %>% mutate(dupes = n(),
                                            row_num = row_number())
table(sabx$dupes)
test <- sabx %>% filter(dupes==2)

#duplicate has lat/long directly on county line and was therefore merged to both counties in a 1:m... 
#just set to the first match
sabx <- sabx %>% filter(row_num==1) %>% select(-dupes, -row_num) %>% ungroup()

##--------------------------------------------------------------##

#create data for the implantation and conception dates!
sabx <- sabx %>% mutate(implant_date = ymd(LMP_ALL)+days(21), #probably won't use implant date as very imprecise
                        conceive_date = ymd(LMP_ALL)+days(14))
summary(sabx)

##--------------------------------------------------------------##
#save CONUS cohort!
saveRDS(sabx, paste0(dataexp, "CONUS_cohort_082423.rds"))

