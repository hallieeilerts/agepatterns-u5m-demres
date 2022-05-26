

################################################################
#
# Analyze age pattern of mortality under-5 in different sources
#
################################################################

# Load packages
library(MortalityLaws)
library(data.table)
library(ggplot2)
library(gridExtra)
library(extrafont)
library(tidyverse)
library(plyr)

# Load spreadsheet with African countries matched to UN geoscheme regions
id.africanRegions <- read.csv("./data/keys/id-africanRegions.csv")

# Load mortality rates that were calculated for HDSS, DHS, and MICS
hdss <- read.csv("./outputs/hdss-5yr.csv")
dhs <- read.csv("./outputs/dhs-nat-5yr.csv")
mics <- read.csv("./outputs/mics-nat-5yr.csv")

# Load data from HMD
age_int  <- 5
year_int <- 5 
interval <- paste0(age_int, "x", year_int)  
cntr <- c("AUS", "AUT", "BLR", "BEL", "BGR", "CAN", "CHL", "HRV", "CZE", "DNK", "EST", "FIN", "FRATNP", "DEUTNP", "DEUTE",  
          "DEUTW", "GRC", "HKG", "HUN", "ISL", "IRL", "ISR", "ITA", "JPN", "LVA", "LTU", "LUX", "NLD", "NZL_NP", "NOR",    
          "POL", "PRT", "KOR", "RUS", "SVK", "SVN", "ESP", "SWE", "CHE", "GBR_NP",  "GBRTENW", "GBR_SCO", "GBR_NIR", "USA", "UKR")
hmdLT <- ReadHMD(what="LT_t",countries=cntr,interval=interval,
                username= HMDusername,
                password= HMDpassword, save = TRUE)

###############
# Prepare data
###############

# Restrict DHS and HDSS to same time period as MICS estimates
dhs <- subset(dhs, SurveyYear >= range(mics$SurveyYear)[1])
hdss$year_l <- as.numeric(as.character(sub('\\-.*', '', hdss$year)))
hdss <- subset(hdss, year_l >= (range(mics$SurveyYear)[1]-10))
hdss$year_l <- NULL

# Calculate DHS range of 5q0
dhsu5mrange <- quantile(subset(dhs, nqx == "U5MR")$est, probs = c(.1,.9)) 

# Format HMD data
hmd <- as.data.frame(hmdLT$data)
# Only keep ages intervals of interest
hmd <- subset(hmd, Age == "0" | Age == "1-4" | Age == "5-9")
hmd <- reshape(hmd, timevar = "Age", idvar = c("country", "Year"),direction = "wide")
# Calculate 5q0
hmd$U5MR <- (hmd$lx.0 - hmd[,"lx.5-9"])/hmd$lx.0 
hmd <- hmd[,c("country", "Year", "qx.0", "qx.1-4","U5MR")]
colnames(hmd) <- c("CountryName", "year", "IMR", "CMR", "U5MR")
hmd <- subset(hmd, !is.na(IMR))

# Remove HMD estimates that have 5q0 that is outside the range (10th and 90th decile) of DHS 5q0 
hmd <- subset(hmd, U5MR >= dhsu5mrange[1] & U5MR <= dhsu5mrange[2])
# Remove GBR_NP before 1980.
hmd$year_l <- as.numeric(sub('\\-.*', '', hmd$year))
hmd <- subset(hmd, !(CountryName == "GBR_NP" & year_l < 1980))
# Delete unnecessary columns
hmd$year_l <- NULL
# Add identifying columns
hmd$Source <- "HMD"
hmd$AfricanRegion <- "All"

# Combine data from all African sources into one wide dataframe

# Format HDSS data
hdss$Source <- "HDSS"
hdss <- dcast(setDT(hdss), ...~ nqx , value.var = "est")
# Format DHS and MICS data
dhs$Source <- "DHS"
mics$Source <- "MICS"
fbh <- rbind(dhs, mics)
fbh <- dcast(setDT(fbh), CountryName+SurveyId+Source+tips ~ nqx , value.var = "est")
africaData <- dplyr::bind_rows(hdss,fbh)
# Add information on African region
africaData <- merge(africaData, id.africanRegions, by = "CountryName")

# Combine with HMD data
allData <- dplyr::bind_rows(africaData, hmd)
allData$Source <- factor(allData$Source, levels=c("HMD","DHS","MICS","HDSS"), ordered=TRUE)

# Prepare CMR:IMR ratio data for plotting
allregions <- africaData
allregions$AfricanRegion <- "All"
u5ratio <- rbind(allData, allregions) 
u5ratio$AfricanRegion <- factor(u5ratio$AfricanRegion, levels=c("All","Southern","Eastern","Central","Western"), ordered=TRUE)
u5ratio$Source <- factor(u5ratio$Source, levels=c("HMD","DHS","MICS","HDSS"), ordered=TRUE)

# Prepare PNMR:NMR ratio data for plotting
u1ratio <- rbind(subset(allData, Source != "HMD"), allregions) 
u1ratio$AfricanRegion <- factor(u1ratio$AfricanRegion, levels=c("All","Southern","Eastern","Central","Western"), ordered=TRUE)
u1ratio$Source <- factor(u1ratio$Source, levels=c("DHS","MICS","HDSS"), ordered=TRUE)


##########
# Figure 1
##########

scaleFUN <- function(x) sprintf("%.3f", x)
panelA <- ggplot() +
  geom_point(data = subset(allData, Source == 'HMD'), aes(x=IMR, y= CMR, color = Source ), alpha=0.8, size =1) +
  geom_point(data = subset(allData, Source == 'DHS'), aes(x=IMR, y= CMR, color = Source ), alpha=0.8, size =1) +
  geom_point(data = subset(allData, Source == 'MICS'), aes(x=IMR, y= CMR, color = Source ), alpha=0.8, size =1) +
  geom_point(data = subset(allData, Source == 'HDSS'), aes(x=IMR, y= CMR, color = Source ), alpha=0.8, size =1) +
  labs(title="a") +
  scale_y_log10(name = "q(12m,5y)", breaks = c( 0.0063, 0.0125, 0.0250, 0.0500, 0.1), labels=scaleFUN) + 
  scale_x_log10(name = "q(12m)" , breaks = c( 0.0063, 0.0125, 0.0250, 0.0500, 0.1), labels=scaleFUN) + 
  coord_cartesian(ylim = c(0.005,.12), xlim=c(0.005,0.12)) +
  scale_color_manual(values = c("#440154FF", "#E1AD01","gray45", "#35B779FF")) + 
  guides(color = FALSE) +   theme_bw() + 
  theme(text=element_text(size=8, family="Arial"), legend.position = "none")
p2 <- ggplot(subset(u5ratio, AfricanRegion == "All")) +
  geom_boxplot(aes(x=Source, y=CMR/IMR), alpha=0.7, outlier.shape = NA) +
  geom_jitter(aes(x=Source, y=CMR/IMR, col = Source), alpha=0.3) + 
  labs(x="",y="", title = "b", subtitle="q(12m,5y)/q(12m)") +
  scale_y_log10() + scale_color_manual(values = c("gray45","#440154FF", "#35B779FF","#E1AD01")) + guides(color=FALSE) + 
  theme_bw() + theme(plot.subtitle = element_text(hjust = 0.5), text=element_text(size=8, family="Arial"), plot.margin=unit(c(.2,.1,0,-0.5), "cm")) 
p3 <- ggplot(subset(u5ratio, AfricanRegion != "All" & Source != "HMD")) +
  geom_boxplot(aes(x=Source, y=CMR/IMR), outlier.shape = NA) +
  geom_jitter(aes(x=Source, y=CMR/IMR, col = Source), alpha=0.3) + 
  labs(x="",y="",caption="b") + scale_y_log10() +
  scale_color_manual(values = c("#440154FF", "#35B779FF", "#E1AD01")) + guides(color=FALSE) + 
  theme_bw() + theme(axis.text.x = element_text(angle = 45), text=element_text(size=8, family="Arial"), plot.margin=unit(c(0,.1,-.5,-0.5), "cm")) +
  facet_wrap(~AfricanRegion, nrow= 1)
panelB <- arrangeGrob(p2, p3, nrow=2)

fig1 <- arrangeGrob(panelA, panelB, widths=c(4,3), padding =0)
ggsave("./outputs/fig1.jpeg", fig1, width=7,height=4,dpi=400)

###########
# Figure 2
###########

panelA <- ggplot() +
  geom_point(data = subset(allData, Source == 'DHS'), aes(x=NMR, y= PNMR, color = Source ), alpha=0.8, size =1) +
  geom_point(data = subset(allData, Source == 'MICS'), aes(x=NMR, y= PNMR, color = Source ), alpha=0.8, size =1) +
  geom_point(data = subset(allData, Source == 'HDSS'), aes(x=NMR, y= PNMR, color = Source ), alpha=0.8, size =1) +
  labs(title="a") +
  scale_y_log10(name = "q(28d,12m)", breaks = c(0.0031, 0.0063, 0.0125, 0.0250, 0.0500, 0.1), labels=scaleFUN) + 
  scale_x_log10(name = "q(28d)", breaks = c(0.0015,0.0031, 0.0063, 0.0125, 0.0250, 0.0500, 0.1), labels=scaleFUN) + 
  coord_cartesian(ylim = c(0.003,.1), xlim=c(0.00125,0.055)) +
  scale_color_manual(values = c("#440154FF", "#E1AD01", "#35B779FF")) + 
  guides(color = FALSE) +   theme_bw() + 
  theme(text=element_text(size=8, family="Arial"), legend.position = "none")
p2 <- ggplot(subset(u1ratio, AfricanRegion == "All")) +
  geom_boxplot(aes(x=Source, y=PNMR/NMR), alpha=0.7, outlier.shape = NA) +
  geom_jitter(aes(x=Source, y=PNMR/NMR, col = Source), alpha=0.3) + 
  labs(x="",y="", title="b",subtitle="q(28d,12m)/q(28d)") +
  scale_y_log10() + scale_color_manual(values = c("#440154FF", "#35B779FF","#E1AD01")) + guides(color=FALSE) + 
  theme_bw() + theme(text=element_text(size=8, family="Arial"), plot.subtitle = element_text(hjust = 0.5),plot.margin=unit(c(.2,.1,.1,-0.5), "cm")) 
p3 <- ggplot(subset(u1ratio, AfricanRegion != "All" )) +
  geom_boxplot(aes(x=Source, y=PNMR/NMR), outlier.shape = NA) +
  geom_jitter(aes(x=Source, y=PNMR/NMR, col = Source), alpha=0.3) + 
  labs(x="",y="") +  scale_y_log10() + 
  scale_color_manual(values = c("#440154FF", "#35B779FF", "#E1AD01")) + guides(color=FALSE) + 
  theme_bw() + theme(axis.text.x = element_text(angle = 45), text=element_text(size=8, family="Arial"),plot.margin=unit(c(0,.1,-.2,-0.5), "cm")) +
  facet_wrap(~AfricanRegion, nrow= 1)
panelB <- arrangeGrob(p2, p3, nrow=2)

fig2 <- arrangeGrob(panelA, panelB, widths=c(4,3), padding =0) #generates g
ggsave("./outputs/fig2.jpeg", fig2, width=7,height=4,dpi=400) 

#############
# Figure A-1
#############

p1 <- ggplot(allData) +
  geom_point(aes(x = U5MR, y=CMR/IMR, col = Source), alpha=0.8) +
  labs(y="q(12m,60m)/q(12m)",x="q(60m)", title = "a") +
  scale_color_manual(values = c("gray45","#440154FF", "#35B779FF","#E1AD01")) +
  scale_x_log10() + scale_y_log10() +
  theme_bw() + theme(text=element_text(size=8, family="Arial"), legend.title = element_blank()) + guides(color=guide_legend(nrow=1,byrow=TRUE)) 
p2 <- ggplot(allData) +
  geom_point(aes(x = IMR, y=PNMR/NMR, col = Source), alpha=0.8) +
  labs(y="q(28d,12m)/q(28d)",x="q(12m)", title = "b") +
  scale_color_manual(values = c("gray45","#440154FF", "#35B779FF","#E1AD01")) +
  scale_x_log10() + scale_y_log10() +
  theme_bw() + theme(text=element_text(size=8, family="Arial"), legend.position = "none")
g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}
mylegend <- g_legend(p1)
p3 <- grid.arrange(arrangeGrob(p1 + theme(legend.position="none"),
                               p2 + theme(legend.position="none"),
                               nrow=1),
                   mylegend, nrow=2,heights=c(10, 1))
ggsave("./outputs/figa1.jpeg", p3, width=7,height=4,dpi=500) 

#############
# Figure A-2
#############

western <- subset(africaData, AfricanRegion == "Western")
western$AfricanRegion <- "Western-All"
westsubcat <- subset(africaData, AfricanRegion2 %in% c("Western-Coastal",  "Western-Sahelian"))
westsubcat$AfricanRegion <- westsubcat$AfricanRegion2
allwest <- rbind(western, westsubcat)
allwest$Region <- factor(allwest$Region, levels=c("Western-Coastal",  "Western-Sahelian","Western-All"), ordered=TRUE)
allwest$Source <- factor(allwest$Source, levels=c("DHS","MICS","HDSS"), ordered=TRUE)

p1 <- ggplot(allwest) +
  geom_boxplot(aes(x=Source, y=PNMR/NMR), outlier.shape = NA) +
  geom_jitter(aes(x=Source, y=PNMR/NMR, col = Source), alpha=0.5) + labs(title="q(28d,12m)/q(28d)",y="",x="") + #title ="Postneonatal/Neonatal mortality") +
  scale_color_manual(values = c("#440154FF", "#35B779FF", "#E1AD01")) + guides(color=FALSE) + 
  scale_y_log10() + theme_bw() + theme(text=element_text(size=10, family="Arial"),plot.title = element_text(hjust = 0.5)) +
  facet_wrap(~AfricanRegion, nrow= 1) 
p2 <- ggplot(allwest) +
  geom_boxplot(aes(x=Source, y=CMR/IMR), outlier.shape = NA) +
  geom_jitter(aes(x=Source, y=CMR/IMR, col = Source), alpha=0.5) + labs(title="q(12m,5y)/q(12m)",y="",x="") + #title ="Postneonatal/Neonatal mortality") +
  scale_color_manual(values = c("#440154FF", "#35B779FF", "#E1AD01")) + guides(color=FALSE) + 
  scale_y_log10() + theme_bw() + theme(text=element_text(size=10, family="Arial"),plot.title = element_text(hjust = 0.5)) + # axis.text.x = element_text(angle = 45)
  facet_wrap(~AfricanRegion, nrow= 1) 

g <- arrangeGrob(p2,p1, widths=c(6), nrow=2, ncol=1) #generates g
ggsave("./outputs/figa2.jpeg", g, width=6,height=5,dpi=500) 

############
# Table A-1 
############

u5ratiotab <- u5ratio
u5ratiotab$AfricanRegion <- as.character(u5ratiotab$AfricanRegion)
u5ratiotab$AfricanRegion[u5ratiotab$AfricanRegion == "Western"] <- "Western-All"
u5ratiotab <- rbind(u5ratiotab, westsubcat)
u5ratiotab$RegionSource <- paste(u5ratiotab$AfricanRegion,"-",u5ratiotab$Source, sep="")
u5ratiotab$RegionSource <- factor(u5ratiotab$RegionSource, levels = c("All-HMD","All-DHS","All-MICS","All-HDSS",
                                                                "Southern-DHS","Southern-MICS","Southern-HDSS", 
                                                                "Eastern-DHS","Eastern-MICS","Eastern-HDSS",
                                                                "Central-DHS","Central-MICS",
                                                                "Western-All-DHS", "Western-All-MICS","Western-All-HDSS",
                                                                "Western-Coastal-DHS", "Western-Coastal-MICS","Western-Coastal-HDSS",
                                                                "Western-Sahelian-DHS", "Western-Sahelian-MICS","Western-Sahelian-HDSS"), ordered = TRUE)

medIQR.u5 <- u5ratiotab %>% 
  mutate(CMR.IMR = CMR/IMR) %>%
  select(RegionSource, CMR.IMR) %>%
  melt(id.vars = c("RegionSource"), variable.name = "Ratio") %>%
  group_by(RegionSource, Ratio) %>%
  dplyr::summarise(
    Median = sprintf("%.2f",median(value, na.rm = TRUE)),
    IQR = paste("(",sprintf("%.2f",round(quantile(value, na.rm=TRUE)[2],2)),", ", sprintf("%.2f",round(quantile(value, na.rm=TRUE)[4],2)),")", sep="")) 

u5medianTest <- u5ratiotab
u5medianTest$CMR.IMR <- u5medianTest$CMR/u5medianTest$IMR
u5dhsRatios <- subset(u5medianTest, Source == "DHS" )
u5otherRatios <- subset(u5medianTest, Source != "DHS")
ids.cmrimr <- unique(u5otherRatios$RegionSource)
ids.cmrimr.dhs <- stringr::str_replace(ids.cmrimr, "HDSS", "DHS")
ids.cmrimr.dhs <- stringr::str_replace(ids.cmrimr.dhs, "MICS", "DHS")
ids.cmrimr.dhs <- stringr::str_replace(ids.cmrimr.dhs, "HMD", "DHS")
l.cmrimr <- list()
for(i in 1:length(ids.cmrimr)){
  l.cmrimr[[i]] <- wilcox.test(subset(u5dhsRatios, RegionSource == ids.cmrimr.dhs[i])$CMR.IMR, subset(u5otherRatios , RegionSource == ids.cmrimr[i])$CMR.IMR)$p.value
}
names(l.cmrimr) <- ids.cmrimr
medTest.u5 <- ldply(l.cmrimr, .id = "RegionSource")
medTest.u5$pvalue[medTest.u5$V1 <0.01] <- "<0.01"
medTest.u5$pvalue[medTest.u5$V1 >= 0.01] <- sprintf("%.2f",round(medTest.u5$V1[medTest.u5$V1 >= 0.01], 2))
medTest.u5 <- medTest.u5[,c("RegionSource","pvalue")]

u5tab <- merge(medIQR.u5, medTest.u5, by = c("RegionSource"), all = TRUE)
u5tab$pvalue[is.na(u5tab$pvalue)] <- "-"
u5tab

############
# Table A-2 
############

u1ratiotab <- u1ratio
u1ratiotab$AfricanRegion <- as.character(u1ratiotab$AfricanRegion)
u1ratiotab$AfricanRegion[u1ratiotab$AfricanRegion == "Western"] <- "Western-All"
u1ratiotab <- rbind(u1ratiotab, westsubcat)
u1ratiotab$RegionSource <- paste(u1ratiotab$AfricanRegion,"-",u1ratiotab$Source, sep="")
u1ratiotab$RegionSource <-  factor(u1ratiotab$RegionSource, levels = c("All-DHS","All-MICS","All-HDSS",
                                                            "Southern-DHS","Southern-MICS","Southern-HDSS", 
                                                            "Eastern-DHS","Eastern-MICS","Eastern-HDSS",
                                                            "Central-DHS","Central-MICS",
                                                            "Western-All-DHS", "Western-All-MICS","Western-All-HDSS",
                                                            "Western-Coastal-DHS", "Western-Coastal-MICS","Western-Coastal-HDSS",
                                                            "Western-Sahelian-DHS", "Western-Sahelian-MICS","Western-Sahelian-HDSS" ), ordered = TRUE)

medIQR.u1<- u1ratiotab %>% 
  mutate(PNMR.NMR = PNMR/NMR) %>%
  select(RegionSource, PNMR.NMR) %>%
  melt(id.vars = c("RegionSource"), variable.name = "Ratio") %>%
  group_by(RegionSource, Ratio) %>%
  dplyr::summarise(
    Median = sprintf("%.2f",median(value, na.rm = TRUE)),
    IQR = paste("(",sprintf("%.2f",round(quantile(value, na.rm=TRUE)[2],2)),", ", sprintf("%.2f",round(quantile(value, na.rm=TRUE)[4],2)),")", sep=""))


u1medianTest <- u1ratiotab
u1medianTest$PNMR.NMR <- u1medianTest$PNMR/u1medianTest$NMR
u1dhsRatios <- subset(u1medianTest, Source == "DHS" )
u1otherRatios <- subset(u1medianTest, Source != "DHS")
ids.pnmrnmr <- unique(u1otherRatios$RegionSource)
ids.pnmrnmr.dhs <- stringr::str_replace(ids.pnmrnmr, "HDSS", "DHS")
ids.pnmrnmr.dhs <- stringr::str_replace(ids.pnmrnmr.dhs, "MICS", "DHS")
l.pnmrnmr <- list()
for(i in 1:length(ids.pnmrnmr)){
  l.pnmrnmr[[i]] <- wilcox.test(subset(u1dhsRatios, RegionSource == ids.pnmrnmr.dhs[i])$PNMR.NMR, subset(u1otherRatios , RegionSource == ids.pnmrnmr[i])$PNMR.NMR)$p.value
}
names(l.pnmrnmr) <- ids.pnmrnmr
medTest.u1 <- ldply(l.pnmrnmr, .id = "RegionSource")
medTest.u1$pvalue[medTest.u1$V1 <0.01] <- "<0.01"
medTest.u1$pvalue[medTest.u1$V1 >= 0.01] <- sprintf("%.2f",round(medTest.u1$V1[medTest.u1$V1 >= 0.01], 2))
medTest.u1 <- medTest.u1[,c("RegionSource","pvalue")]

u1tab <- merge(medIQR.u1, medTest.u1, by = c("RegionSource"), all = TRUE)
u1tab$pvalue[is.na(u1tab$pvalue)] <- "-"
u1tab

