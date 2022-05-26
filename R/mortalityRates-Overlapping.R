

##############################################################################
#
# Calculate mortality rates for overlapping periods of HDSS/DHS and HDSS/MICS
#
##############################################################################

# Clear environment
rm(list = ls())

# Load packages
library(tidyverse)
library(data.table)
library(lubridate)
library(plyr)
library(demogsurv)
library(survival)

# Load spreadsheet with HDSS site names, centre ids, and countries
id.hdss <- read.csv("./data/keys/id-hdss.csv")
# Load spreadsheet with MICS/DHS survey years and region names that have been hand-matched to HDSS sites
id.regions <- read.csv("./data/keys/id-regions.csv")

# Load DHS and MICS mortality rates
dhs5yr <- read.csv("./outputs/dhs-reg-5yr.csv")
mics5yr <- read.csv("./outputs/mics-reg-5yr.csv")
dhs10yr <- read.csv("./outputs/dhs-reg-10yr.csv")
mics10yr <- read.csv("./outputs/mics-reg-10yr.csv")

# Set working directory to location of iShare HDSS data
setwd("./data/hdss")
# Load in iShare csv files as list
tbl_fread <- list.files(pattern = "*.csv") %>% map_df(~fread(.))
l.hdss <- split( tbl_fread , tbl_fread$CentreId )
# Reset working directory to original location
setwd(sub('\\/data/hdss.*', '', getwd()))

# Define nested lapply function to perform operations on list of lists.
nested_lapply <- function(data, fun) {
  lapply(data, function(sublist) { lapply(sublist, fun) })  }

# Define function to construct a model matrix for aggregating mortality rates over age groups.
.mm_aggr <- function(mf, agegr){
  
  if(any(!vapply(mf, is.factor, logical(1)))){
    v <- !vapply(mf, is.factor, logical(1))
    stop(paste("Not all 'by' variables are factors:",
               paste(names(v)[v], collapse = ", ")))
  }
  
  if(!"agegr" %in% names(mf))
    stop("'agegr' variable not in mf")
  
  ## Calculate duration of each age group
  dur <- setNames(diff(agegr), levels(mf$agegr))
  dur <- dur[mf$agegr]
  
  mf$agegr <- NULL
  
  if(length(mf))
    mf$byf <- do.call(interaction, c(mf, drop=TRUE))
  
  if(!length(mf) || length(levels(mf$byf)) == 1)
    mf$byf <- 1
  
  mm <- model.matrix(~-1+byf, mf) * dur
  
  df <- mf[!duplicated(mf$byf), , drop=FALSE]
  df <- df[order(df$byf), , drop=FALSE]
  df["byf"] <- NULL
  rownames(df) <- NULL
  
  list(df = df, mm = mm)
}

#########################################################
# Identify overlapping periods between HDSS and DHS/MICS
#########################################################

# Identify duration of data collection in each HDSS site
l.hdss <- lapply(l.hdss, function(x){ x$event.date <- decimal_date(as.Date(x$EventDate, tryFormats = "%Y/%m/%d")) ; return(x)})
# Identify first year of HDSS data collection, save in data frame
HDSSstart <- lapply(l.hdss, function(x) floor(min(subset(x, EventCode == "ENU")$event.date)))
HDSSstart <- plyr::ldply(HDSSstart)
names(HDSSstart) <- c("CentreId","HDSS.start")
# Identify last year of HDSS data collection, save in data frame
HDSSend <- lapply(l.hdss, function(x) floor(max(x$event.date, na.rm = TRUE)))
HDSSend <- plyr::ldply(HDSSend)
names(HDSSend) <- c("CentreId","HDSS.end")
HDSSend$HDSS.end <- ifelse(HDSSend$HDSS.end > 2019, 2019, HDSSend$HDSS.end)
HDSSyears <- merge(HDSSstart, HDSSend)
# Add on HDSS 'CentreId' and DHS 'SurveyFinalRegion' 
# DHS region names have chagned over time. 'SurveyFinalRegion' is the name of the region in which each HDSS is located as of the most recent DHS survey.
id.hdss <- id.regions[!duplicated(id.regions[,c("CentreId","HDSS")]),][,c("CentreId","HDSS","SurveyFinalRegion")]
# HDSSyears has the duration of data collection in each HDSS
HDSSyears <- merge(HDSSyears, id.hdss)

# Combine DHS and MICS 5yr mortality rates
fbh.5yr <- rbind(dhs5yr, mics5yr)

# Join DHS/MICS rates with HDSS that are located in the survey region
overlaps5yr <- plyr::join(fbh.5yr[,c("SurveyId","SurveyYear","SurveyFinalRegion","tips")], HDSSyears, by = "SurveyFinalRegion")
# Create columns for lower and upper years of DHS/MICS TIPS rates
overlaps5yr$year_l <- ifelse(overlaps5yr$tips == "0-4", overlaps5yr$SurveyYear - 5, overlaps5yr$SurveyYear - 10)
overlaps5yr$year_u <- overlaps5yr$year_l + 4
overlaps5yr <- overlaps5yr[!duplicated(overlaps5yr),]
# Only keep rows where the HDSS has available data for the same period as the DHS/MICS rate
overlaps5yr <- subset(overlaps5yr, HDSS.start <= year_l & year_u <= HDSS.end)
# Remove any period mortality rate that begins before 1990
overlaps5yr <- subset(overlaps5yr, year_l >= 1990)

# Combine DHS and MICS 10yr mortality rates
fbh.10yr <- rbind(dhs10yr, mics10yr)

# Join DHS/MICS rates with HDSS that are located in the survey region
overlaps10yr <- plyr::join(fbh.10yr[,c("SurveyId","SurveyYear","SurveyFinalRegion","tips")], HDSSyears, by = "SurveyFinalRegion")
overlaps10yr$year_l <- overlaps10yr$SurveyYear - 10
overlaps10yr$year_u <- overlaps10yr$SurveyYear - 1
overlaps10yr <- overlaps10yr[!duplicated(overlaps10yr),]
# Ensure that each HDSS has available data for the same period as the FBH subnational region rate
overlaps10yr <- subset(overlaps10yr, HDSS.start <= year_l & year_u <= HDSS.end)
# Remove any period mortality rate that begins before 1990
overlaps10yr <- subset(overlaps10yr, year_l >= 1990)


##############################################################################
# Calculate HDSS mortality rates for overlapping 5 year periods with DHS/MICS
##############################################################################

# Only keep HDSS sites that have overlapping 5yr periods with DHS/MICS
l.hdssOverlaps <- l.hdss[names(l.hdss) %in% overlaps5yr$CentreId]
l.hdssOverlaps <-  l.hdssOverlaps[order(names(l.hdssOverlaps))]

# Create new DOB column without time information
l.hdssOverlaps <- lapply(l.hdssOverlaps, function(x){ x$dob <- decimal_date(as.Date(x$DoB)) ; return(x)})
# Remove records of deliveries
l.hdssOverlaps <- lapply(l.hdssOverlaps, function(x){ x[x$EventCode != "DLV",]})
# Create start event column
l.hdssOverlaps <- lapply(l.hdssOverlaps, function(x){ x[x$EventCode %in% c('BTH','ENU','IMG'), 'start'] <- x[x$EventCode %in% c('BTH','ENU','IMG'), 'event.date'] ; return(x)})
# Create stop event column
l.hdssOverlaps <- lapply(l.hdssOverlaps, function(x){ x[x$EventCode %in% c('DTH','OMG','OBE'), 'stop'] <- x[x$EventCode %in% c('DTH','OMG','OBE'), 'event.date'] ; return(x)})
# Create event marker
l.hdssOverlaps <- lapply(l.hdssOverlaps, function(x){ x[x$EventCode == "DTH", "event"] <- 1 ; return(x)})
l.hdssOverlaps <- lapply(l.hdssOverlaps, function(x){ x[x$EventCode != "DTH", "event"] <- 0 ; return(x)})
# Take previous rows nonmissing 'start' value and fill down
l.hdssOverlaps <- l.hdssOverlaps %>% map(~.x %>% fill(start, .direction = "down"))

# Keep nonmissing minimum stop for each individual/start event
l.hdssOverlaps <- lapply(l.hdssOverlaps, function(x){ddply(x, ~IndividualId+dob+start, plyr::summarize, 
                                           stop = min(stop, na.rm = T),
                                           event = max(event, na.rm = T) )})

# If someone has a NA for a stop event, exclude the episode. 
# This means the individual had two start events in a row, and we'll only keep the second one (i.e. ENU, IMG, OBE)
l.hdssOverlaps <- lapply(l.hdssOverlaps, function(x){x[!is.na(x$stop),]})
l.hdssOverlaps <- lapply(l.hdssOverlaps, function(x){x[!(x$stop=="Inf"),]})

# Drop cases with dob > stop
l.hdssOverlaps <- lapply(l.hdssOverlaps, function(x){x[x$dob <= x$stop,]})
# Drop cases with start > stop
l.hdssOverlaps <- lapply(l.hdssOverlaps, function(x){x[x$start <= x$stop,]})
# If stop == Inf, set to 2019.0
l.hdssOverlaps <- lapply(l.hdssOverlaps, function(x){ x$stop[is.infinite(x$stop)] <- 2019 ; return(x)})

# Count the number of times an HDSS appears in the 5yr overlaps list
# Repeat the data frame that many times
# We will conduct a separate period mortality calculation on each list element
expanded <- list()
for(i in 1:length(l.hdssOverlaps)){
  site <- names(l.hdssOverlaps[i]) 
  reps <- nrow(subset(overlaps5yr, CentreId == site))
  expanded[[i]] <- rep(list(l.hdssOverlaps[i]), reps)
}

# Coerce into data frames
expanded <- nested_lapply(expanded, function(x){ x <- ldply(x) ; return(x) })
# Unlist to create one data frame per list element
l.hdssExp5yr <- unlist(expanded, recursive = FALSE)

# Create another list that will contain the periods for each calculation
overlaps5yr$ID <- paste(overlaps5yr$CentreId, overlaps5yr$SurveyId, overlaps5yr$year_l, overlaps5yr$year_u, sep = ".")
l.periods <- split(overlaps5yr[,c("year_l","year_u")], overlaps5yr$ID)
l.periods <- lapply(l.periods, function(x){x$year_u <- x$year_u + 1 ; return(x)})

# Pair HDSS data with periods
l.calc5yr <- Map(list, l.hdssExp5yr, l.periods) 
names(l.calc5yr) <- names(l.periods)

# Create list of mortality indicators which contain traditional age breakdowns 
agegr <- list(c(0,1)/12, c(1, 3, 5, 12)/12, c(0, 1, 3, 5, 12)/12, c(12, 24, 36, 48, 60)/12, c(0, 1, 3, 5, 12, 24, 36, 48, 60)/12)
names(agegr) <- c("NMR","PNMR","IMR","CMR","U5MR")

# For all age groups, HDSS, and periods, calculate person-years
# This will create a list of lists, where the first level is the indicator (NMR, PNMR, etc.),
# and the second level is the HDSS.
l.py <- lapply(agegr, function(d){
  lapply(l.calc5yr, function(a,b){py <- demog_pyears(~1, a[[1]], agegr=b, origin = 0, scale = 1, period = a[[2]],
                                                   dob = "dob", tstart = "start", tstop = "stop", event = "event")$data 
  return(py)}, b=d)})

# Set factor combinations for aggregating person-years.
byvar <- c("agegr", "period") 

# Create column for factor combinations (age group and period)
l.aggr <- nested_lapply(l.py, function(x){ x$byf <- interaction(x[byvar], drop = TRUE) ; return(x) })

# Divide events by person-years to calculate mortality rates for each age group and period
l.aggr <- nested_lapply(l.aggr, function(x){ x <- x[order(x$byf),] ; return(x) })
l.mx <- nested_lapply(l.aggr, function(x){ x$event/x$pyears  })

# Create model matrices to aggregate piecewise-constant rates to cumulative hazards.
# Drop any duplicate factor combinations.
l.pred <- nested_lapply(l.aggr, function(x){ data.frame(x[c(byvar,"byf")])[!duplicated(x$byf),] ; return(x) })
l.pred <- nested_lapply(l.pred, function(x){ x <- x[order(x$byf),] ; return(x) })

# For each mortality indicator (NMR, PNMR, IMR, CMR, U5MR) and each HDSS,
# create model matrix with time periods in columns and age group length in rows.
l.mm <- list()
for(i in 1:length(agegr)){
  indicator <- l.pred[[i]]
  age <- agegr[[i]]
  l.mm[[i]] <- lapply(indicator, function(x) .mm_aggr(x[c("agegr", "period")], agegr = age)$mm )
}
names(l.mm) <- names(agegr)

# Flatten list of lists
ul.mx <- unlist(l.mx, recursive = FALSE)
ul.mm <- unlist(l.mm, recursive = FALSE)

# Label columns of model matrix with period names
periodnames <- unlist(l.pred, recursive = FALSE)
periodnames <- lapply(periodnames, function(x) unique(x[,"period"]))
for(i in 1:length(periodnames)){
  colnames(ul.mm[[i]]) <- periodnames[[i]]
}

# Aggregate piecewise-constant rates to cumulative hazards
# (multiply of mortality rates for each age segment by length of segment)
l.cummx <- Map('%*%', ul.mx, ul.mm)

# Calculate nqx from mx
l.nqx <- lapply(l.cummx, function(x) 1-exp(-x))

# Format nqx into data frame
l.nqx <- lapply(l.nqx, as.data.frame(t))
l.nqx <- lapply(l.nqx, function(x){ x$year <- rownames(x) ; return(x) })
df.nqx <- plyr::ldply(l.nqx, .id = "ID")

# Format identifying columns
df.nqx$ID <- as.character(df.nqx$ID)
df.nqx$nqx <- sub('\\..*', '', df.nqx$ID)
df.nqx$ID <- sub("^[^\\.]*\\.", "", df.nqx$ID)
names(df.nqx)[which(names(df.nqx) == "value")] <- "est"
df.nqx$CentreId <- substr(df.nqx$ID, 1, 5)

# Merge on additional identifying information
df.nqx <- merge(df.nqx, id.hdss, by = "CentreId")

# Save mortality rates as csv
write.csv(df.nqx, "./outputs/hdss-overlaps-5yr.csv", row.names = FALSE)

##############################################################################
# Calculate HDSS mortality rates for overlapping 10 year periods with DHS/MICS
##############################################################################

# Count the number of times an HDSS appears in the 10 yr overlaps list
# Repeat the data frame that many times
# we will conduct a separate period mortality calculation on each
expanded <- list()
for(i in 1:length(l.hdssOverlaps)){
  site <- names(l.hdssOverlaps[i]) 
  reps <- nrow(subset(overlaps10yr, CentreId == site))
  expanded[[i]] <- rep(list(l.hdssOverlaps[i]), reps)
}

# Coerce into data frames
expanded <- nested_lapply(expanded, function(x){ x <- ldply(x) ; return(x) })
# Unlist to create one data frame per list element
l.hdssExp10yr <- unlist(expanded, recursive = FALSE)

# Create another list that will contain the periods for each calculation
overlaps10yr$ID <- paste(overlaps10yr$CentreId, overlaps10yr$SurveyId, overlaps10yr$year_l, overlaps10yr$year_u, sep = ".")
l.periods <- split(overlaps10yr[,c("year_l","year_u")], overlaps10yr$ID)
l.periods <- lapply(l.periods, function(x){x$year_u <- x$year_u + 1 ; return(x)})

# Pair HDSS data with periods
l.calc10yr <- Map(list, l.hdssExp10yr, l.periods) 
names(l.calc10yr) <- names(l.periods)

# For all age groups, HDSS, and periods, calculate person-years
# This will create a list of lists, where the first level is the indicator (NMR, PNMR, etc.),
# and the second level is the HDSS.
l.py <- lapply(agegr, function(d){
  lapply(l.calc10yr, function(a,b){py <- demog_pyears(~1, a[[1]], agegr=b, origin = 0, scale = 1, period = a[[2]],
                                                     dob = "dob", tstart = "start", tstop = "stop", event = "event")$data 
  return(py)}, b=d)})

# Set factor combinations for aggregating person-years.
byvar <- c("agegr", "period") 

# Create column for factor combinations (age group and period)
l.aggr <- nested_lapply(l.py, function(x){ x$byf <- interaction(x[byvar], drop = TRUE) ; return(x) })

# Divide events by person-years to calculate mortality rates for each age group and period
l.aggr <- nested_lapply(l.aggr, function(x){ x <- x[order(x$byf),] ; return(x) })
l.mx <- nested_lapply(l.aggr, function(x){ x$event/x$pyears  })

# Create model matrices to aggregate piecewise-constant rates to cumulative hazards.
# Drop any duplicate factor combinations.
l.pred <- nested_lapply(l.aggr, function(x){ data.frame(x[c(byvar,"byf")])[!duplicated(x$byf),] ; return(x) })
l.pred <- nested_lapply(l.pred, function(x){ x <- x[order(x$byf),] ; return(x) })


# For each mortality indicator (NMR, PNMR, IMR, CMR, U5MR) and each HDSS,
# create model matrix with time periods in columns and age group length in rows.
l.mm <- list()
for(i in 1:length(agegr)){
  indicator <- l.pred[[i]]
  age <- agegr[[i]]
  l.mm[[i]] <- lapply(indicator, function(x) .mm_aggr(x[c("agegr", "period")], agegr = age)$mm )
}
names(l.mm) <- names(agegr)

# Flatten list of lists
ul.mx <- unlist(l.mx, recursive = FALSE)
ul.mm <- unlist(l.mm, recursive = FALSE)

# Label columns of model matrix with period names
periodnames <- unlist(l.pred, recursive = FALSE)
periodnames <- lapply(periodnames, function(x) unique(x[,"period"]))
for(i in 1:length(periodnames)){
  colnames(ul.mm[[i]]) <- periodnames[[i]]
}

# Aggregate piecewise-constant rates to cumulative hazards
# (multiply of mortality rates for each age segment by length of segment)
l.cummx <- Map('%*%', ul.mx, ul.mm)

# Calculate nqx from mx
l.nqx <- lapply(l.cummx, function(x) 1-exp(-x))

# Format nqx into data frame
l.nqx <- lapply(l.nqx, as.data.frame(t))
l.nqx <- lapply(l.nqx, function(x){ x$year <- rownames(x) ; return(x) })
df.nqx <- plyr::ldply(l.nqx, .id = "ID")

# Format identifying columns
df.nqx$ID <- as.character(df.nqx$ID)
df.nqx$nqx <- sub('\\..*', '', df.nqx$ID)
df.nqx$ID <- sub("^[^\\.]*\\.", "", df.nqx$ID)
names(df.nqx)[which(names(df.nqx) == "value")] <- "est"
df.nqx$CentreId <- substr(df.nqx$ID, 1, 5)

# Merge on additional identifying information
df.nqx <- merge(df.nqx, id.hdss, by = "CentreId")

# Save mortality rates as csv
write.csv(df.nqx, "./outputs/hdss-overlaps-10yr.csv", row.names = FALSE)



