
################################################################################
#
# Master script to replicate results
#
# Hallie Eilerts, hallieeilerts@gmail.com
#
################################################################################


# This is a master script that calls other scripts one-by-one to replicate the main results of the manuscript,
# "Age patterns of under-five mortality in sub-Saharan Africa during 1990-2018: A comparison of estimates from demographic surveillance with full birth histories and the historic record"

# Step 1. Set location of downloaded manuscript files
setwd("C:/yourfolder")

# Step 2. Download HDSS and MICS data and save in "data" folder ("C:/yourfolder/data")
# http://www.indepth-ishare.org/index.php/catalog/central
# http://mics.unicef.org/surveys

# Step 3. Set up account at DHS website
# Instructions are available here: https://dhsprogram.com/data/Access-Instructions.cfm
# Once you have created an account, request access to datasets for the sub-Saharan African region.
# Once you have received access to the datasets, set your rdhs credentials using the function set_rdhs_config(). 
# This will require your DHS account email and project as arguments.
install.packages("rdhs")
library(rdhs)
set_rdhs_config(email = "youremail@email.com",
                project = "Your DHS project",
                config_path = "rdhs.json",
                global = FALSE)

# Step 4. Calculate mortality rates for DHS data
source("R/mortalityRates-DHS.R")

# Step 5. Calculate mortality rates for MICS data
source("R/mortalityRates-MICS.R")

# Step 6. Calculate mortality rates for HDSS data
source("R/mortalityRates-HDSS.R")

# Step 7. Calculate HDSS mortality rates for overlapping 5-yr periods with DHS and MICS
source("R/mortalityRates-Overlapping.R")

# Step 8. Set up account at HMD website
# https://www.mortality.org/
# Once you have created an account, set your hmd credentials
HMDusername <- "youremail@email.com"
HMDpassword <- "yourpassword"

# Step 9. Analyze age pattern of mortality under-5
source("R/analysis-AgePattern.R")

# Step 10. Analyze differences in estimates from HDSS and DHS/MICS
source("R/analysis-Differences.R")


