
#########################################
#
# Calculate mortality rates for DHS data
#
#########################################

# Clear environment
rm(list = ls())

# Load packages
library(demogsurv)
library(tidyverse)
library(plyr)

# Load spreadsheet with MICS/DHS survey years and region names that have been hand-matched to HDSS sites
id.regions <- read.csv("./data/keys/id-regions.csv")

#######################################
# Download DHS data using rdhs package
#######################################

# Return api requests as data tables
Sys.setenv(rdhs_DATA_TABLE = "TRUE")

# Sub-Saharan African countries
countries <- dhs_countries()
cc <- subset(countries, RegionName == "Sub-Saharan Africa")$DHS_CountryCode

# Surveys from sub-Saharan African surveys
surveys <- dhs_surveys(countryIds = cc, surveyYearStart=1995, surveyYearEnd=2018, surveyType = "DHS")
# Note: Data used in the analysis was downloaded in November 2020. 
# Later downloads may include more surveys that have been subsequently added to the DHS data repository.

# All available birth recode files
surveysbr <- dhs_datasets(fileType = "BR", fileFormat = "flat")

# Surveys from sub-Saharan Africa that have birth recode files
brd <- surveysbr[which(surveysbr$SurveyId %in% surveys$SurveyId),]

# Use rdhs to retreive datasets, downloading them from DHS website if not already in the rdhs cache.
Sys.setenv("rdhs_LOUD_DOWNLOAD" = TRUE)
brd$path <- unlist(get_datasets(brd$FileName, download_option = "both")) 

# Load all of the datasets into R as a list.
br <- list()
for(survid in brd$SurveyId){
  print(survid)
  dat <- readRDS(brd[which(brd$SurveyId == survid),]$path)
  dat <- dat[grep("caseid|^v0|^v1|^b", names(dat))]
  br[[survid]] <- dat
}

# Convert to factors
br <- lapply(br, haven::as_factor)

# Add survey-level variables
br <- Map(data.frame,
          SurveyId = brd$SurveyId,
          CountryName = brd$CountryName,
          SurveyYear = brd$SurveyYear,
          br)

# Adjust ethiopian calendar for the CMC columns
br$ET2000DHS$b3 <- br$ET2000DHS$b3 + (92)
br$ET2005DHS$b3 <- br$ET2005DHS$b3 + (92)
br$ET2011DHS$b3 <- br$ET2011DHS$b3 + (92)
br$ET2016DHS$b3 <- br$ET2016DHS$b3 + (92)
br$ET2000DHS$v008 <- br$ET2000DHS$v008 + (92)
br$ET2005DHS$v008 <- br$ET2005DHS$v008 + (92)
br$ET2011DHS$v008 <- br$ET2011DHS$v008 + (92)
br$ET2016DHS$v008 <- br$ET2016DHS$v008 + (92)
br$ET2000DHS$v011 <- br$ET2000DHS$v011 + (92)
br$ET2005DHS$v011 <- br$ET2005DHS$v011 + (92)
br$ET2011DHS$v011 <- br$ET2011DHS$v011 + (92)
br$ET2016DHS$v011 <- br$ET2016DHS$v011 + (92)

# Create indicator for deaths, and calculate date of death
br <- lapply(br, function(x){x$death <- x$b5 == "no"; x})
br <- lapply(br, function(x){x$dod <- x$b3 + x$b7 + 0.5; x})

#######################################
# National level mortality calculations
########################################

# Calculate 5-year period nqx for each indicator of interest
nat1 <- ldply(br, calc_nqx, by=~CountryName+SurveyYear, tips = c(0, 5, 10), agegr=c(0, 1)/12, .id = "SurveyId")
nat2 <- ldply(br, calc_nqx, by=~CountryName+SurveyYear, tips = c(0, 5, 10), agegr=c(1, 3, 5, 12)/12, .id = "SurveyId")
nat3 <- ldply(br, calc_nqx, by=~CountryName+SurveyYear, tips = c(0, 5, 10), agegr=c(0, 1, 3, 5, 12)/12, .id = "SurveyId")
nat4 <- ldply(br, calc_nqx, by=~CountryName+SurveyYear, tips = c(0, 5, 10), agegr=c(12, 24, 36, 48, 60)/12, .id = "SurveyId")
nat5 <- ldply(br, calc_nqx, by=~CountryName+SurveyYear, tips = c(0, 5, 10), agegr=c(0, 1, 3, 5, 12, 24, 36, 48, 60)/12, .id = "SurveyId")
nat1$nqx <- "NMR"
nat2$nqx <- "PNMR"
nat3$nqx <- "IMR"
nat4$nqx <- "CMR"
nat5$nqx <- "U5MR"
df.nat <- rbind(nat1, nat2, nat3, nat4, nat5)

# Save mortality rates as csv
write.csv(df.nat, "./outputs/dhs-nat-5yr.csv", row.names = FALSE)

####################################################
# Subnational regional level mortality calculations
####################################################

# Keep DHS regions of interest
br.regions <- br[names(br) %in% id.regions$SurveyId]
br.regions <- lapply(br.regions, function(x){x <- merge(x, id.regions[,c("SurveyId","SurveyRegion","SurveyFinalRegion")], 
                                                        by.x = c("SurveyId", "v024"), by.y = c("SurveyId","SurveyRegion")) ; x })
br.regions <- br.regions[sapply(br.regions, nrow)>0]

# Calculate 5-year period nqx for each indicator of interest
reg1 <- ldply(br.regions, calc_nqx, by=~CountryName+SurveyYear+SurveyFinalRegion, tips = c(0, 5, 10), agegr=c(0, 1)/12, .id="SurveyId")
reg2 <- ldply(br.regions, calc_nqx, by=~CountryName+SurveyYear+SurveyFinalRegion, tips = c(0, 5, 10), agegr=c(1, 3, 5, 12)/12, .id="SurveyId")
reg3 <- ldply(br.regions, calc_nqx, by=~CountryName+SurveyYear+SurveyFinalRegion, tips = c(0, 5, 10), agegr=c(0, 1, 3, 5, 12)/12, .id="SurveyId")
reg4 <- ldply(br.regions, calc_nqx, by=~CountryName+SurveyYear+SurveyFinalRegion, tips = c(0, 5, 10), agegr=c(12, 24, 36, 48, 60)/12, .id="SurveyId")
reg5 <- ldply(br.regions, calc_nqx, by=~CountryName+SurveyYear+SurveyFinalRegion, tips = c(0, 5, 10), agegr=c(0, 1, 3, 5, 12, 24, 36, 48, 60)/12, .id="SurveyId")
reg1$nqx <- "NMR"
reg2$nqx <- "PNMR"
reg3$nqx <- "IMR"
reg4$nqx <- "CMR"
reg5$nqx <- "U5MR"
df.reg <- rbind(reg1, reg2, reg3, reg4, reg5)

# Save mortality rates as csv
write.csv(df.reg, "dhs-reg-5yr.csv", row.names = FALSE)

# Calculate 10-year period nqx for each indicator of interest
reg6 <- ldply(br.regions, calc_nqx, by=~CountryName+SurveyYear+SurveyFinalRegion, tips = c(0, 10), agegr=c(0, 1)/12, .id="SurveyId")
reg7 <- ldply(br.regions, calc_nqx, by=~CountryName+SurveyYear+SurveyFinalRegion, tips = c(0, 10), agegr=c(1, 3, 5, 12)/12, .id="SurveyId")
reg8 <- ldply(br.regions, calc_nqx, by=~CountryName+SurveyYear+SurveyFinalRegion, tips = c(0, 10), agegr=c(0, 1, 3, 5, 12)/12, .id="SurveyId")
reg9 <- ldply(br.regions, calc_nqx, by=~CountryName+SurveyYear+SurveyFinalRegion, tips = c(0, 10), agegr=c(12, 24, 36, 48, 60)/12, .id="SurveyId")
reg10 <- ldply(br.regions, calc_nqx, by=~CountryName+SurveyYear+SurveyFinalRegion, tips = c(0, 10), agegr=c(0, 1, 3, 5, 12, 24, 36, 48, 60)/12, .id="SurveyId")
reg6$nqx <- "NMR"
reg7$nqx <- "PNMR"
reg8$nqx <- "IMR"
reg9$nqx <- "CMR"
reg10$nqx <- "U5MR"
df.reg10yr <- rbind(reg6, reg7, reg8, reg9, reg10)

# Save mortality rates as csv
write.csv(df.reg10yr, "./outputs/dhs-reg-10yr.csv", row.names = FALSE)

#########################################
# Residence level mortality calculations
#########################################

# Keep DHS residence of interest
br.residence <- br[names(br) %in% id.regions$SurveyId]
id.residence <- id.regions[,c("SurveyId","SurveyResidence")]
id.residence <- id.residence[!duplicated(id.residence),]
br.residence <- lapply(br.residence, function(x){x <- merge(x, id.residence, 
                                                        by.x = c("SurveyId", "v025"), by.y = c("SurveyId","SurveyResidence")) ; x })
br.residence <- br.residence[sapply(br.residence, nrow)>0]

res1 <- ldply(br.residence, calc_nqx, by=~CountryName+SurveyYear+v025, tips = c(0, 5, 10), agegr=c(0, 1)/12, .id = "SurveyId")
res2 <- ldply(br.residence, calc_nqx, by=~CountryName+SurveyYear+v025, tips = c(0, 5, 10), agegr=c(1, 3, 5, 12)/12, .id = "SurveyId")
res3 <- ldply(br.residence, calc_nqx, by=~CountryName+SurveyYear+v025, tips = c(0, 5, 10), agegr=c(0, 1, 3, 5, 12)/12, .id = "SurveyId")
res4 <- ldply(br.residence, calc_nqx, by=~CountryName+SurveyYear+v025, tips = c(0, 5, 10), agegr=c(12, 24, 36, 48, 60)/12, .id = "SurveyId")
res5 <- ldply(br.residence, calc_nqx, by=~CountryName+SurveyYear+v025, tips = c(0, 5, 10), agegr=c(0, 1, 3, 5, 12, 24, 36, 48, 60)/12, .id = "SurveyId")
res1$nqx <- "NMR"
res2$nqx <- "PNMR"
res3$nqx <- "IMR"
res4$nqx <- "CMR"
res5$nqx <- "U5MR"
df.res <- rbind(res1, res2, res3, res4, res5)
names(df.res)[which(names(df.res) == "v025")] <- "SurveyResidence"

write.csv(df.res, "./outputs/dhs-res-5yr.csv", row.names = FALSE)


