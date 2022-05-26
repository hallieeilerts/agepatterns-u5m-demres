
#################################################
#
# Calculate mortality rates for iShare HDSS data
#
#################################################

# Clear environment
rm(list = ls())

# Load packages
library(demogsurv)
library(survival)
library(tidyverse)
library(data.table)
library(lubridate)
library(plyr)

# Load spreadsheet with HDSS site names, centre ids, and countries
id.hdss <- read.csv("./data/keys/id-hdss.csv", sep=",", header=TRUE)

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

#####################
# Prepare HDSS data
#####################

# Create new event date column without time information
l.hdss <- lapply(l.hdss, function(x){ x$event.date <- decimal_date(as.Date(x$EventDate, tryFormats = "%Y/%m/%d")) ; return(x)})
# Identify first year of HDSS data collection, save in data frame
HDSSstart <- lapply(l.hdss, function(x) floor(min(subset(x, EventCode == "ENU")$event.date)))
HDSSstart <- plyr::ldply(HDSSstart)
names(HDSSstart) <- c("CentreId","HDSS.start")

# Create new DOB column without time information
l.hdss <- lapply(l.hdss, function(x){ x$dob <- decimal_date(as.Date(x$DoB, tryFormats = "%Y/%m/%d")) ; return(x)})
# Remove records of deliveries
l.hdss <- lapply(l.hdss, function(x){ x[x$EventCode != "DLV",]})
# Create start event column
l.hdss <- lapply(l.hdss, function(x){ x[x$EventCode %in% c('BTH','ENU','IMG'), 'start'] <- x[x$EventCode %in% c('BTH','ENU','IMG'), 'event.date'] ; return(x)})
# Create stop event column
l.hdss <- lapply(l.hdss, function(x){ x[x$EventCode %in% c('DTH','OMG','OBE'), 'stop'] <- x[x$EventCode %in% c('DTH','OMG','OBE'), 'event.date'] ; return(x)})
# Create event marker
l.hdss <- lapply(l.hdss, function(x){ x[x$EventCode == "DTH", "event"] <- 1 ; return(x)})
l.hdss <- lapply(l.hdss, function(x){ x[x$EventCode != "DTH", "event"] <- 0 ; return(x)})
# Take previous rows nonmissing 'start' value and fill down
l.hdss <- l.hdss %>% map(~.x %>% fill(start, .direction = "down"))

# Keep nonmissing minimum stop for each individual/start event
l.hdss <- lapply(l.hdss, function(x){ddply(x, ~IndividualId+dob+start, plyr::summarize,
                                           stop = min(stop, na.rm = T),
                                           event = max(event, na.rm = T) )})

# If someone has a NA for a stop event, exclude the episode.
# This means the individual had two start events in a row, and we'll only keep the second one (i.e. ENU, IMG, OBE)
l.hdss <- lapply(l.hdss, function(x){x[!is.na(x$stop),]})
l.hdss <- lapply(l.hdss, function(x){x[!(x$stop=="Inf"),]})

# Drop cases with dob > stop
l.hdss <- lapply(l.hdss, function(x){x[x$dob <= x$stop,]})
# Drop cases with start > stop
l.hdss <- lapply(l.hdss, function(x){x[x$start <= x$stop,]})
# If stop == Inf, set to 2019.0
l.hdss <- lapply(l.hdss, function(x){ x$stop[is.infinite(x$stop)] <- 2019 ; return(x)})

##########################
# Mortality calculations
##########################

# Create list of mortality indicators which contain traditional age breakdowns
agegr <- list(c(0,1)/12, c(1, 3, 5, 12)/12, c(0, 1, 3, 5, 12)/12, c(12, 24, 36, 48, 60)/12, c(0, 1, 3, 5, 12, 24, 36, 48, 60)/12)
names(agegr) <- c("NMR","PNMR","IMR","CMR","U5MR")

# For all age groups and HDSS, calculate person-years for five-year periods, counting backwards from the most recent year of data collection.
# This will create a list of lists, where the first level is the indicator (NMR, PNMR, etc.),
# and the second level is the HDSS.
l.py <- lapply(agegr, function(d){
  lapply(l.hdss, function(a,b){intv.max <- ifelse(max(a$stop) > 2019, 2019, max(a$stop))
  intv.int <- c(intv.max - 25,intv.max - 20,intv.max - 15,intv.max - 10,intv.max - 5, intv.max)
  py <- demog_pyears(~1, a, agegr=b, origin = 0, scale = 1, period = intv.int,
                     dob = "dob", tstart = "start", tstop = "stop",event = "event")$data
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
# Order rows and add person-years column.
l.pred <- nested_lapply(l.pred, function(x){ x <- x[order(x$byf),] ; return(x) })

# For each mortality indicator (NMR, PNMR, IMR, CMR, U5MR) and each HDSS,
# create model matrix with time periods in columns and age group length in rows.
l.mm <- list()
for(i in 1:length(agegr)){
  indicator <- l.pred[[i]]
  age <- agegr[[i]]
  l.mm[[i]] <- lapply(indicator, function(x) .mm_aggr(x[c("agegr", "period")], agegr = age))
}
names(l.mm) <- names(agegr)

# Label columns of model matrix with period names
l.mm <- nested_lapply(l.mm, function(x){ colnames(x$mm) <- x$df$period ; return(x) })
l.mm <- nested_lapply(l.mm, function(x){ x$mm })

# Flatten list of lists
ul.mx <- unlist(l.mx, recursive = FALSE)
ul.mm <- unlist(l.mm, recursive = FALSE)

# Aggregate piecewise-constant rates to cumulative hazards
# (multiply of mortality rates for each age segment by length of segment)
l.cummx <- Map('%*%', ul.mx, ul.mm)

# Calculate nqx from mx
l.nqx <- lapply(l.cummx, function(x) 1-exp(-x))

# Format nqx into data frame
l.nqx <- lapply(l.nqx, as.data.frame(t))
l.nqx <- lapply(l.nqx, function(x){ x$year <- rownames(x) ; return(x) })
df.nqx <- plyr::ldply(l.nqx)

# Format identifying columns
names(df.nqx) <- c("ID","est","year")
df.nqx$nqx <- sub("\\..*", "", df.nqx$ID)
df.nqx$CentreId <- sub('.*\\.', "", df.nqx$ID)

# Merge on additional identifying information
df.nqx <- merge(df.nqx, id.hdss, by = "CentreId")

# Remove 5 year rates that start before the beginning of HDSS data collection
df.nqx <- merge(df.nqx, HDSSstart)
df.nqx$periodstart <-  as.numeric(sub("\\-.*", "", df.nqx$year))
df.nqx <- subset(df.nqx, periodstart >= HDSS.start)
df.nqx <- df.nqx[order(df.nqx$CentreId, df.nqx$nqx, df.nqx$periodstart),]
df.nqx$periodstart <- df.nqx$HDSS.start <- df.nqx$ID <- NULL

# Save mortality rates as csv
write.csv(df.nqx, "./outputs/hdss-5yr.csv", row.names = FALSE)


