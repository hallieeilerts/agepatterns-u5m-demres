

##################################################
#
# Analyze differences in HDSS and FBH estimates
#
#################################################

# Clear environment
rm(list = ls())

# Load packages
library(ggplot2)
library(extrafont)
library(stargazer)

# Load id keys
# Spreadsheet with MICS/DHS regions that have been hand-matched to HDSS sites
id.regions <- read.csv("./data/keys/id-regions.csv")

# Load mortality rates that were calculated for overlapping HDSS and DHS/MICS subnational regions
hdss5yr <- read.csv("./outputs/hdss-overlaps-5yr.csv")
dhs5yr <- read.csv("./outputs/dhs-reg-5yr.csv")
mics5yr <- read.csv("./outputs/mics-reg-5yr.csv")
hdss10yr <- read.csv("./outputs/hdss-overlaps-10yr.csv")
dhs10yr <- read.csv("./outputs/dhs-reg-10yr.csv")
mics10yr <- read.csv("./outputs/mics-reg-10yr.csv")

# Load mortality rates for DHS/MICS national and residence level
dhsres <- read.csv("./outputs/dhs-res-5yr.csv")
dhsnat <- read.csv("./outputs/dhs-nat-5yr.csv")
micsres <- read.csv("./outputs/mics-res-5yr.csv")
micsnat <- read.csv("./outputs/mics-nat-5yr.csv")


####################################################################
# Subnational region level differences between 5yr period estimates
###################################################################

# Combine DHS and MICS data
dhs5yr$Source <- "DHS"
mics5yr$Source <- "MICS"
fbh5yr <- rbind(dhs5yr, mics5yr)
fbh5yr$year[fbh5yr$tips == "0-4"] <- paste(fbh5yr$SurveyYear[fbh5yr$tips == "0-4"] - 5, fbh5yr$SurveyYear[fbh5yr$tips == "0-4"] - 1, sep = ".")
fbh5yr$year[fbh5yr$tips == "5-9"] <- paste(fbh5yr$SurveyYear[fbh5yr$tips == "5-9"] - 10, fbh5yr$SurveyYear[fbh5yr$tips == "5-9"] - 6, sep = ".")
# Merge on corresponding HDSS information
fbh5yr <- merge(fbh5yr, id.regions[,c("SurveyId","SurveyFinalRegion","CentreId")], by = c("SurveyId","SurveyFinalRegion"))
fbh5yr$ID <- paste(fbh5yr$CentreId, fbh5yr$SurveyId, fbh5yr$year, sep = ".")

# Merge HDSS estimates with DHS/MICS estimates, and calculate differences
dif5yr <- merge(hdss5yr, fbh5yr, by =c("ID","nqx"),
                     suffixes = c(".hdss", ".fbh"))
dif5yr$dif <- (dif5yr$est.hdss - dif5yr$est.fbh)/(dif5yr$est.hdss + dif5yr$est.fbh)/2 * 100
dif5yr$nqx <- factor(dif5yr$nqx, levels=c("NMR","PNMR","IMR","CMR","U5MR"),
                          labels=c("Neonatal","Postneonatal","Infant","Child","Under-five"), ordered=TRUE)

############
# Table A-3
############

dif5yr %>%
  group_by(nqx, Source) %>%
  dplyr::summarise(
    Median = sprintf("%.2f",median(dif, na.rm = TRUE)),
    IQR = paste("(",sprintf("%.2f",round(quantile(dif, na.rm=TRUE)[2],2)),", ", 
                sprintf("%.2f",round(quantile(dif, na.rm=TRUE)[4],2)),")", sep=""))


###########
# Figure 3
###########

ggplot(dif5yr) +
  geom_boxplot(aes(x=Source, y=dif), outlier.shape = NA) +
  geom_jitter(aes(x=Source, y=dif, col =Source), alpha=0.5) + 
  geom_hline(yintercept = 0, col = "red3") + theme_bw() +
  scale_color_manual(values = c("#440154FF","#35B779FF")) +
  scale_fill_manual(values = c("#440154FF","#35B779FF")) +
  theme(legend.position = "none", text=element_text(size=10, family="Arial")) + 
  xlab("") + ylab("Relative percentage difference") +
  facet_wrap(~nqx, nrow=1) +
  ggsave("./outputs/fig3.jpeg",width=6,height=2.5,dpi=400) 

#####################################################################
# Subnational region level differences between 10yr period estimates
####################################################################

# Combine DHS and MICS data
dhs10yr$Source <- "DHS"
mics10yr$Source <- "MICS"
fbh10yr <- rbind(dhs10yr, mics10yr)
fbh10yr$year <- paste(fbh10yr$SurveyYear - 10, fbh10yr$SurveyYear - 1, sep = ".")
# Merge on corresponding HDSS information
fbh10yr <- merge(fbh10yr, id.regions[,c("SurveyId","SurveyFinalRegion","CentreId")], by = c("SurveyId","SurveyFinalRegion"))
fbh10yr$ID <- paste(fbh10yr$CentreId, fbh10yr$SurveyId, fbh10yr$year, sep = ".")

# Merge HDSS estimates with DHS/MICS estimates, and calculate differences
dif10yr <- merge(hdss10yr, fbh10yr, by =c("ID","nqx"),
                     suffixes = c(".hdss", ".fbh"))
dif10yr$dif <- (dif10yr$est.hdss - dif10yr$est.fbh)/(dif10yr$est.hdss + dif10yr$est.fbh)/2 * 100
dif10yr$nqx <- factor(dif10yr$nqx, levels=c("NMR","PNMR","IMR","CMR","U5MR"),
                          labels=c("Neonatal","Postneonatal","Infant","Child","Under-five"), ordered=TRUE)

##############
# Figure A-4
##############

ggplot(dif10yr) +
  geom_boxplot(aes(x=Source, y=dif), outlier.shape = NA) +
  geom_jitter(aes(x=Source, y=dif, col =Source), alpha=0.5) + 
  geom_hline(yintercept = 0, col = "red3") + theme_bw() +
  scale_color_manual(values = c("#440154FF","#35B779FF")) +
  scale_fill_manual(values = c("#440154FF","#35B779FF")) +
  theme(legend.position = "none", text=element_text(size=10, family="Arial")) + 
  xlab("") + ylab("Relative percentage difference") +
  facet_wrap(~nqx, nrow=1) +
  ggsave("./outputs/figa4.jpeg",width=6,height=2.5,dpi=400) 


###########################################################
# Residence level differences between 5yr period estimates
###########################################################

# Combine DHS and MICS residence level data
dhsres$Source <- "DHS"
micsres$Source <- "MICS"
fbhres <- rbind(dhsres, micsres)
fbhres$year[fbhres$tips == "0-4"] <- paste(fbhres$SurveyYear[fbhres$tips == "0-4"] - 5, fbhres$SurveyYear[fbhres$tips == "0-4"] - 1, sep = ".")
fbhres$year[fbhres$tips == "5-9"] <- paste(fbhres$SurveyYear[fbhres$tips == "5-9"] - 10, fbhres$SurveyYear[fbhres$tips == "5-9"] - 6, sep = ".")
# Merge on corresonding HDSS information to FBH
fbhres <- merge(fbhres, id.regions[,c("SurveyId","SurveyResidence","CentreId")], by = c("SurveyId","SurveyResidence"))
fbhres$ID <- paste(fbhres$CentreId, fbhres$SurveyId, fbhres$year, sep = ".")

# Merge HDSS estimates with residence level DHS/MICS estimates, and calculate differences
difres <- merge(hdss5yr, fbhres, by =c("ID","nqx"),suffixes = c(".hdss", ".fbh"))
difres$dif <- (difres$est.hdss - difres$est.fbh)/(difres$est.hdss + difres$est.fbh)/2 * 100
difres$nqx <- factor(difres$nqx, levels=c("NMR","PNMR","IMR","CMR","U5MR"),
                     labels=c("Neonatal","Postneonatal","Infant","Child","Under-five"), ordered=TRUE)
difres$level <- "Residence"

##########################################################
# National level differences between 5yr period estimates
##########################################################

# Combine DHS and MICS national level data
dhsnat$Source <- "DHS"
micsnat$Source <- "MICS"
fbhnat <- rbind(dhsnat, micsnat)
fbhnat$year[fbhnat$tips == "0-4"] <- paste(fbhnat$SurveyYear[fbhnat$tips == "0-4"] - 5, fbhnat$SurveyYear[fbhnat$tips == "0-4"] - 1, sep = ".")
fbhnat$year[fbhnat$tips == "5-9"] <- paste(fbhnat$SurveyYear[fbhnat$tips == "5-9"] - 10, fbhnat$SurveyYear[fbhnat$tips == "5-9"] - 6, sep = ".")
# Merge on corresponding HDSS information to FBH
fbhnat <- merge(fbhnat, id.regions[,c("SurveyId","CountryName","CentreId")], by = c("SurveyId","CountryName"))
fbhnat$ID <- paste(fbhnat$CentreId, fbhnat$SurveyId, fbhnat$year, sep = ".")

# Merge HDSS estimates with national level DHS/MICS estimates, and calculate differences
difnat <- merge(hdss5yr, fbhnat, by =c("ID","nqx"),suffixes = c(".hdss", ".fbh"))
difnat$dif <- (difnat$est.hdss - difnat$est.fbh)/(difnat$est.hdss + difnat$est.fbh)/2 * 100
difnat$nqx <- factor(difnat$nqx, levels=c("NMR","PNMR","IMR","CMR","U5MR"),
                     labels=c("Neonatal","Postneonatal","Infant","Child","Under-five"), ordered=TRUE)
difnat$level <- "Nation"

# Combine national and residence level differences with subnational region differences
difreg <- dif5yr
difreg$level <- "Region"
alldif <- dplyr::bind_rows(difreg, difres, difnat)
alldif$level <- factor(alldif$level, levels = c("Region","Residence","Nation"), ordered = TRUE)
alldif$nqx <- factor(alldif$nqx, levels = c("Neonatal", "Postneonatal", "Infant", "Child", "Under-five"), labels=c("NMR","PNMR","IMR","CMR","U5MR"), ordered=TRUE)

##############
# Figure A-5
##############

ggplot(alldif, aes(x=nqx, y=dif, col = Source, group = interaction(Source, nqx))) +
  geom_point(position=position_jitterdodge(), alpha = .3) +
  geom_boxplot(fill = NA, outlier.shape = NA, show.legend=FALSE ) +
  geom_boxplot(fill = NA, outlier.shape = NA, col = "black") +
  geom_hline(yintercept = 0, col = "red3") + theme_bw() +
  scale_color_manual(values = c("#440154FF","#35B779FF")) +
  scale_fill_manual(values = c("#440154FF","#35B779FF")) +
  theme(legend.position = "bottom", legend.title = element_blank(),text=element_text(size=10, family="Arial")) + 
  xlab("") +ylab("Relative percentage difference") +
  guides(colour = guide_legend(override.aes = list(alpha = 1))) +
  facet_wrap(~level) + 
  ggsave("./outputs/figa5.jpeg", dpi = 500, width = 6, height = 4)

