
#########################################
#
# Calculate mortality rates for MICS data
#
#########################################

# Clear environment
rm(list = ls())

# Load packages
library(foreign)
library(demogsurv)
library(survival)
library(tidyverse)
library(plyr)

# Load spreadsheet with MICS survey names, years, and countries
id.mics <- read.csv("./data/keys/id-mics.csv")
# Load spreadsheet with MICS/DHS survey years and region names that have been hand-matched to HDSS sites
id.regions <- read.csv("./data/keys/id-regions.csv")

# Set working directory to location of MICS data
setwd("./data/mics")
# Load in MICS birth history files as list
filenames <- paste(list.files(full.names=TRUE), "/bh.sav", sep="")
filenames <- filenames[file.exists(filenames)]
l.mics <- filenames %>% map(~read.spss(., to.data.frame=TRUE))
names(l.mics) <- gsub(".*[/]([^.]+)[/].*", "\\1", filenames)
# Reset working directory to original location
setwd(sub('\\/data/mics.*', '', getwd()))

#####################
# Prepare MICS data
#####################

# Remove punctuation and decapitalize variable names
l.mics <- lapply(l.mics, function(x){names(x) <- str_replace_all(names(x), "[[:punct:]]", " ") ; return(x)})
l.mics <- lapply(l.mics, function(x){names(x) <- tolower(names(x)) ; return(x)})

# Single survey adjustments

# Load in weights from Nigeria woman's questionnaire
nigeria.wm <- read.spss("./data/mics/Nigeria_MICS5_Datasets/wm.sav", to.data.frame=TRUE)
nigeria.wm <- nigeria.wm[,c("HH1","HH2","LN","wmweight")]
l.mics$Nigeria_MICS5_Datasets<- left_join(l.mics$Nigeria_MICS5_Datasets, nigeria.wm, by = c("hh1" = "HH1", "hh2" = "HH2", "ln" = "LN"))
# Remove observations missing weight in Sao Tome and Principe survey
l.mics$SaoTomeandPrincipe_2014_MICS_Datasets <- subset(l.mics$SaoTomeandPrincipe_2014_MICS_Datasets, !is.na(wmweight))
# Add residence (hh6) for Senegal survey hh6
l.mics$SenegalDakar_MICS5_Datasets$hh6 <- "urban"

# Standardize date of birth and age at death variables across surveys

# Convert factors to characters
l.mics <- lapply(l.mics, function(x){if("bh4c" %in% names(x)){x$bh4c <- as.character(x$bh4c)}; x})
l.mics <- lapply(l.mics, function(x){if("bh4m" %in% names(x)){x$bh4m <- as.character(x$bh4m)}; x})
l.mics <- lapply(l.mics, function(x){if("bh4y" %in% names(x)){x$bh4y <- as.character(x$bh4y)}; x})
l.mics <- lapply(l.mics, function(x){if("bh4yc" %in% names(x)){x$bh4yc <- as.character(x$bh4yc)}; x})
l.mics <- lapply(l.mics, function(x){if("bh4mc" %in% names(x)){x$bh4mc <- as.character(x$bh4mc)}; x})
l.mics <- lapply(l.mics, function(x){if("bh5" %in% names(x)){x$bh5 <- as.character(x$bh5)}; x})
l.mics <- lapply(l.mics, function(x){if("bh9u" %in% names(x)){x$bh9u <- as.character(x$bh9u)}; x})
l.mics <- lapply(l.mics, function(x){if("bh9n" %in% names(x)){x$bh9n <- as.character(x$bh9n)}; x})
l.mics <- lapply(l.mics, function(x){if("bh9c" %in% names(x)){x$bh9c <- as.character(x$bh9c)}; x})
l.mics <- lapply(l.mics, function(x){if("bh9uc" %in% names(x)){x$bh9uc <- as.character(x$bh9uc)}; x})
l.mics <- lapply(l.mics, function(x){if("bh9nc" %in% names(x)){x$bh9nc <- as.character(x$bh9nc)}; x})

# Replace all missing values with NA 
missingcodes <- c("missing|not|inconsistent|dk|no response|manquant|nsp|inco|non|pas|em falta")

l.mics <- lapply(l.mics, function(x){if("bh4c" %in% names(x)){x$bh4c[grep(missingcodes, x$bh4c, ignore.case=TRUE)] <- NA }; x })
l.mics <- lapply(l.mics, function(x){if("bh4m" %in% names(x)){x$bh4m[grep(missingcodes, x$bh4m, ignore.case=TRUE)] <- NA }; x })
l.mics <- lapply(l.mics, function(x){if("bh4y" %in% names(x)){x$bh4y[grep(missingcodes, x$bh4y, ignore.case=TRUE)] <- NA }; x })
l.mics <- lapply(l.mics, function(x){if("bh4yc" %in% names(x)){x$bh4yc[grep(missingcodes, x$bh4yc, ignore.case=TRUE)] <- NA }; x })
l.mics <- lapply(l.mics, function(x){if("bh4mc" %in% names(x)){x$bh4mc[grep(missingcodes, x$bh4mc, ignore.case=TRUE)] <- NA }; x })
l.mics <- lapply(l.mics, function(x){if("bh9u" %in% names(x)){x$bh9u[grep(missingcodes, x$bh9u, ignore.case=TRUE)] <- NA }; x })
l.mics <- lapply(l.mics, function(x){if("bh9n" %in% names(x)){x$bh9n[grep(missingcodes, x$bh9n, ignore.case=TRUE)] <- NA }; x })
l.mics <- lapply(l.mics, function(x){if("bh9c" %in% names(x)){x$bh9c[grep(missingcodes, x$bh9c, ignore.case=TRUE)] <- NA }; x })
l.mics <- lapply(l.mics, function(x){if("bh9uc" %in% names(x)){x$bh9uc[grep(missingcodes, x$bh9uc, ignore.case=TRUE)] <- NA }; x })
l.mics <- lapply(l.mics, function(x){if("bh9nc" %in% names(x)){x$bh9nc[grep(missingcodes, x$bh9nc, ignore.case=TRUE)] <- NA }; x })

# Convert dob and age at death back into numeric
l.mics <- lapply(l.mics, function(x){if("bh4c" %in% names(x)){x$bh4c <- as.numeric(x$bh4c) };x})
l.mics <- lapply(l.mics, function(x){if("bh9n" %in% names(x)){x$bh9n <- as.numeric(x$bh9n) };x})
l.mics <- lapply(l.mics, function(x){if("bh9nc" %in% names(x)){x$bh9nc <- as.numeric(x$bh9nc) };x})
l.mics <- lapply(l.mics, function(x){if("bh9c" %in% names(x)){x$bh9c <- as.numeric(x$bh9c) };x})

# Standardize d/m/y language
days <- c("jour|dia|days")
months <- c("m")
years <- c("an|years")
l.mics <- lapply(l.mics, function(x){if("bh9u" %in% names(x)){x$bh9u[grep(days, x$bh9u, ignore.case=TRUE)] <- "Days" }; x })
l.mics <- lapply(l.mics, function(x){if("bh9u" %in% names(x)){x$bh9u[grep(months, x$bh9u, ignore.case=TRUE)] <- "Months" }; x })
l.mics <- lapply(l.mics, function(x){if("bh9u" %in% names(x)){x$bh9u[grep(years, x$bh9u, ignore.case=TRUE)] <- "Years" }; x })
l.mics <- lapply(l.mics, function(x){if("bh9uc" %in% names(x)){x$bh9uc[grep(days, x$bh9uc, ignore.case=TRUE)] <- "Days" }; x })
l.mics <- lapply(l.mics, function(x){if("bh9uc" %in% names(x)){x$bh9uc[grep(months, x$bh9uc, ignore.case=TRUE)] <- "Months" }; x })
l.mics <- lapply(l.mics, function(x){if("bh9uc" %in% names(x)){x$bh9uc[grep(years, x$bh9uc, ignore.case=TRUE)] <- "Years" }; x })

# bh4c (dob, CMC)
# Check if any surveys don't have bh4c
no.bh4c <- lapply(l.mics, function(x) !("bh4c" %in% colnames(x)))
names(l.mics)[unlist(no.bh4c)]
# Use ccdob as bh4c for these surveys
l.mics <- lapply(l.mics, function(x){if(!("bh4c" %in% names(x))){x$bh4c <- x$ccdob}; x})
# Check if there are missing values of 9999
has.9999 <- lapply(l.mics, function(x) "bh4c" %in% names(x) & nrow(subset(x, bh4c == 9999)) > 0)
names(l.mics)[unlist(has.9999)]
# Replace 9999 values with a CMC calculation using bh4yc bh4mc
for(i in 1:length(l.mics)){
  if(names(l.mics)[i] %in% names(l.mics)[unlist(has.9999)]){
    if("bh4yc" %in% names(l.mics[[i]])){
      l.mics[[i]]$bh4mc <- as.numeric(levels(l.mics[[i]]$bh4mc))[l.mics[[i]]$bh4mc] 
      l.mics[[i]]$bh4yc <- as.numeric(levels(l.mics[[i]]$bh4yc))[l.mics[[i]]$bh4yc] 
      l.mics[[i]][l.mics[[i]]$bh4c == 9999, "bh4c"] <-  
        (l.mics[[i]][l.mics[[i]]$bh4c == 9999, "bh4yc"]-1900)*12 + l.mics[[i]][l.mics[[i]]$bh4c == 9999, "bh4mc"]+0.5
    }
  }
}

# bh5 (still alive)
# Standardize these codes in a new variable called "death"
l.mics <- lapply(l.mics, function(x){x$death <- ifelse(x$bh5 %in% c("YES","Yes", 1, "Oui", "Sim"), 0, 1); x})

# Age at death
# bh9u (age at death unit)
# bh9n (age at death number)
# bh9uc (imputed age at death unit)
# bh9nc (imputed age at death number)
# bh9c (imputed age at death in months) 

# Create standardized age at death in new variable called b7
l.mics <- lapply(l.mics, function(x){ x$b7 <- NA ; x})

# If available, use bh9u and bh9n
l.mics <- lapply(l.mics, function(x){ if("bh9u" %in% names(x) & "bh9n" %in% names(x)){ 
                                      x$b7 <- ifelse(!is.na(x$bh9u) & !is.na(x$bh9n) & x$bh9u == "Days", 0.5, NA)
                                      x$b7 <- ifelse(!is.na(x$bh9u) & !is.na(x$bh9n) & x$bh9u == "Months", x$bh9n+0.5, x$b7)
                                      x$b7 <- ifelse(!is.na(x$bh9u) & !is.na(x$bh9n) & x$bh9u == "Years", x$bh9n*12+0.5, x$b7) }; x})

# If bh9u and bh9n are missing, use bh9uc and bh9nc
l.mics <- lapply(l.mics, function(x){if("bh9uc" %in% names(x) & "bh9nc" %in% names(x) ){ 
                                      x$b7 <- ifelse(is.na(x$b7) & !is.na(x$bh9uc) & x$bh9uc == "Days" , 0.5, x$b7)
                                      x$b7 <- ifelse(is.na(x$b7) & !is.na(x$bh9uc) & x$bh9uc == "Months" & !is.na(x$bh9nc), x$bh9nc+0.5, x$b7)
                                      x$b7 <- ifelse(is.na(x$b7) & !is.na(x$bh9uc) & x$bh9uc == "Years" & !is.na(x$bh9nc), x$bh9nc*12+0.5, x$b7) }; x})
# If still missing, use bh9c
l.mics <- lapply(l.mics, function(x){if("bh9c" %in% names(x)){ 
                                      x$b7 <- ifelse(is.na(x$b7) & !is.na(x$bh9c), x$bh9c+0.5, x$b7) }; x})

# If bh9u is missing and bh9n is not, use bh9n as age at death in months
l.mics <- lapply(l.mics, function(x){ if("bh9u" %in% names(x) & "bh9n" %in% names(x)){ 
                                      x$b7 <- ifelse(is.na(x$b7) & is.na(x$bh9u) & !is.na(x$bh9n), x$bh9n+0.5, x$b7) }; x})

# Standardize hh6
l.mics <- lapply(l.mics, function(x){x$hh6 <- as.character(x$hh6) ; x})
l.mics <- lapply(l.mics, function(x){if("hh6" %in% names(x)){x$hh6[x$hh6 %in% c("Urban","Urbain","Urbano")] <- "urban" }; x })
l.mics <- lapply(l.mics, function(x){if("hh6" %in% names(x)){x$hh6[x$hh6 %in% c("Rural")] <- "rural" }; x })
l.mics <- lapply(l.mics, function(x){ x$hh6 <- factor(x$hh6) ; return(x)})

# Recode remaining variables to DHS codes
l.mics <- lapply(l.mics, function(x){x$b3 <- x$bh4c; x})
l.mics <- lapply(l.mics, function(x){x$dod <- x$b3 + x$b7; x})
l.mics <- lapply(l.mics, function(x){x$v005 <- x$wmweight ; x})
l.mics <- lapply(l.mics, function(x){x$v021 <- x$hh1 ; x})
l.mics <- lapply(l.mics, function(x){x$v008 <- x$wdoi ; x})
l.mics <- lapply(l.mics, function(x){x$v025 <- x$hh6 ; x})
l.mics <- lapply(l.mics, function(x){x$v024 <- tolower(x$hh7) ; x})

# Create columns for SurveyId and CountryName
df.mics <- plyr::ldply(l.mics, .id = "SurveyId")
l.mics <- split(df.mics, df.mics$SurveyId)

#######################################
# National level mortality calculations
########################################


# Calculate 5-year period nqx for each indicator of interest
nat1 <- ldply(l.mics, calc_nqx, by=~SurveyId, tips = c(0, 5, 10), agegr=c(0, 1)/12, .id="SurveyId")
nat2 <- ldply(l.mics, calc_nqx, by=~SurveyId, tips = c(0, 5, 10), agegr=c(1, 3, 5, 12)/12, .id="SurveyId")
nat3 <- ldply(l.mics, calc_nqx, by=~SurveyId, tips = c(0, 5, 10), agegr=c(0, 1, 3, 5, 12)/12, .id="SurveyId")
nat4 <- ldply(l.mics, calc_nqx, by=~SurveyId, tips = c(0, 5, 10), agegr=c(12, 24, 36, 48, 60)/12, .id="SurveyId")
nat5 <- ldply(l.mics, calc_nqx, by=~SurveyId, tips = c(0, 5, 10), agegr=c(0, 1, 3, 5, 12, 24, 36, 48, 60)/12, .id="SurveyId")
nat1$nqx <- "NMR"
nat2$nqx <- "PNMR"
nat3$nqx <- "IMR"
nat4$nqx <- "CMR"
nat5$nqx <- "U5MR"
df.nat <- rbind(nat1, nat2, nat3, nat4, nat5)

# Add information for MICS country names and survey years
df.nat <- merge(df.nat, id.mics, by = "SurveyId")

# Save mortality rates as csv
write.csv(df.nat, "./outputs/mics-nat-5yr.csv", row.names = FALSE)

####################################################
# Subnational regional level mortality calculations
####################################################

# Keep MICS regions of interest
mics.regions <- l.mics[names(l.mics) %in% id.regions$SurveyId]
# Add region to Nyanza dataset
if("Kenya-NyanzaProvince_2011_MICS_Datasets" %in% names(mics.regions)){
  mics.regions$`Kenya-NyanzaProvince_2011_MICS_Datasets`$v024 <- "nyanza"
}
mics.regions <- lapply(mics.regions, function(x){x <- merge(x, id.regions[,c("SurveyId","SurveyRegion","SurveyFinalRegion")], 
                                                        by.x = c("SurveyId", "v024"), by.y = c("SurveyId","SurveyRegion")) ; x })
mics.regions <-  mics.regions[sapply(mics.regions, nrow)>0]
# Remove Nigeria. No overlap with HDSS.
mics.regions$Nigeria_MICS5_Datasets <- NULL

# Calculate 5-year period nqx for each indicator of interest
reg1 <- ldply(mics.regions, calc_nqx, by=~SurveyId+SurveyFinalRegion, tips = c(0, 5, 10), agegr=c(0, 1)/12, .id="SurveyId")
reg2 <- ldply(mics.regions, calc_nqx, by=~SurveyId+SurveyFinalRegion, tips = c(0, 5, 10), agegr=c(1, 3, 5, 12)/12, .id="SurveyId")
reg3 <- ldply(mics.regions, calc_nqx, by=~SurveyId+SurveyFinalRegion, tips = c(0, 5, 10), agegr=c(0, 1, 3, 5, 12)/12, .id="SurveyId")
reg4 <- ldply(mics.regions, calc_nqx, by=~SurveyId+SurveyFinalRegion, tips = c(0, 5, 10), agegr=c(12, 24, 36, 48, 60)/12, .id="SurveyId")
reg5 <- ldply(mics.regions, calc_nqx, by=~SurveyId+SurveyFinalRegion, tips = c(0, 5, 10), agegr=c(0, 1, 3, 5, 12, 24, 36, 48, 60)/12, .id="SurveyId")
reg1$nqx <- "NMR"
reg2$nqx <- "PNMR"
reg3$nqx <- "IMR"
reg4$nqx <- "CMR"
reg5$nqx <- "U5MR"
df.reg5yr <- rbind(reg1, reg2, reg3, reg4, reg5)
# Add information for MICS country names and survey years
df.reg5yr <- merge(df.reg5yr, id.mics, by="SurveyId")

# Save mortality rates as csv
write.csv(df.reg5yr, "mics-reg-5yr.csv", row.names = FALSE)

# Calculate 10-year period nqx for each indicator of interest
reg6 <- ldply(mics.regions, calc_nqx, by=~SurveyId+SurveyFinalRegion, tips = c(0, 10), agegr=c(0, 1)/12, .id="SurveyId")
reg7 <- ldply(mics.regions, calc_nqx, by=~SurveyId+SurveyFinalRegion, tips = c(0, 10), agegr=c(1, 3, 5, 12)/12, .id="SurveyId")
reg8 <- ldply(mics.regions, calc_nqx, by=~SurveyId+SurveyFinalRegion, tips = c(0, 10), agegr=c(0, 1, 3, 5, 12)/12, .id="SurveyId")
reg9 <- ldply(mics.regions, calc_nqx, by=~SurveyId+SurveyFinalRegion, tips = c(0, 10), agegr=c(12, 24, 36, 48, 60)/12, .id="SurveyId")
reg10 <- ldply(mics.regions, calc_nqx, by=~SurveyId+SurveyFinalRegion, tips = c(0, 10), agegr=c(0, 1, 3, 5, 12, 24, 36, 48, 60)/12, .id="SurveyId")
reg6$nqx <- "NMR"
reg7$nqx <- "PNMR"
reg8$nqx <- "IMR"
reg9$nqx <- "CMR"
reg10$nqx <- "U5MR"
df.reg10yr <- rbind(reg6, reg7, reg8, reg9, reg10)
# Add information for MICS country names and survey years
df.reg10yr <- merge(df.reg10yr, id.mics, by="SurveyId")

# Save mortality rates as csv
write.csv(df.reg10yr, "./outputs/mics-reg-10yr.csv", row.names = FALSE)

#########################################
# Residence level mortality calculations
#########################################

# Keep MICS residence of interest
l.residence <- l.mics[names(l.mics) %in% id.regions$SurveyId]
id.residence <- id.regions[,c("SurveyId","SurveyResidence")]
id.residence <- id.residence[!duplicated(id.residence),]
l.residence <- lapply(l.residence, function(x){x <- merge(x, id.residence,
                                                           by.x = c("SurveyId", "v025"),  by.y = c("SurveyId","SurveyResidence")) ; x })
l.residence <- l.residence[sapply(l.residence, nrow)>0]

res1 <- ldply(l.residence, calc_nqx, by=~SurveyId+v025, tips = c(0, 5, 10), agegr=c(0, 1)/12, .id = "SurveyId")
res2 <- ldply(l.residence, calc_nqx, by=~SurveyId+v025, tips = c(0, 5, 10), agegr=c(1, 3, 5, 12)/12, .id = "SurveyId")
res3 <- ldply(l.residence, calc_nqx, by=~SurveyId+v025, tips = c(0, 5, 10), agegr=c(0, 1, 3, 5, 12)/12, .id = "SurveyId")
res4 <- ldply(l.residence, calc_nqx, by=~SurveyId+v025, tips = c(0, 5, 10), agegr=c(12, 24, 36, 48, 60)/12, .id = "SurveyId")
res5 <- ldply(l.residence, calc_nqx, by=~SurveyId+v025, tips = c(0, 5, 10), agegr=c(0, 1, 3, 5, 12, 24, 36, 48, 60)/12, .id = "SurveyId")
res1$nqx <- "NMR"
res2$nqx <- "PNMR"
res3$nqx <- "IMR"
res4$nqx <- "CMR"
res5$nqx <- "U5MR"
df.res <- rbind(res1, res2, res3, res4, res5)
names(df.res)[which(names(df.res) == "v025")] <- "SurveyResidence"
# Add information for MICS country names and survey years
df.res <- merge(df.res, id.mics, by="SurveyId")

write.csv(df.res, "./outputs/mics-res-5yr.csv", row.names = FALSE)


