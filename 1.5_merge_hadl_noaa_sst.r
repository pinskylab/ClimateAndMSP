
## Set working directories
if(Sys.info()["nodename"] == "pinsky-macbookair"){
	setwd('~/Documents/Rutgers/Range projections/proj_ranges/')
	}
if(Sys.info()["nodename"] == "amphiprion.deenr.rutgers.edu"){
	setwd('~/Documents/range_projections/')
	}
if(Sys.info()["user"] == "lauren"){
	setwd('~/backup/NatCap/proj_ranges/')
}

#######################################################
#  Surftemp observations are missing almost entirely from three surveys (Newfoundland Fall and Spring, WC_Ann), as well as missing from individuals hauls in the other surveys. We used the global Merged Hadley-NOAA/OI Sea Surface Temperature and Sea-Ice Concentration dataset to fill in these missing observations when possible. We validated this method by comparing the merged Hadley-NOAA values to observed survey surface temperatures when possible, and used regressions to correct for any observed bias.
# Merged Hadley-NOAA/OI Sea Surface Temperature & Sea-Ice Concentration (Hurrell et al, 2008) - See more at: https://climatedataguide.ucar.edu/climate-data/merged-hadley-noaaoi-sea-surface-temperature-sea-ice-concentration-hurrell-et-al-2008#sthash.T8eqg21N.dpuf
# Citation: James W. Hurrell, James J. Hack, Dennis Shea, Julie M. Caron, and James Rosinski, 2008: A New Sea Surface Temperature and Sea Ice Boundary Dataset for the Community Atmosphere Model. J. Climate, 21, 5145–5153. doi: http://dx.doi.org/10.1175/2008JCLI2292.1
# Data obtained from here:

####################################################### 
# Validation/Calibration
# Compare existing surftemp records to Hadley-NOAA SST to determine if there is bias. Correct for this bias using an appropriate regression:
# Use WCAnn to validate WCTri survey (same region). 
# Since only 9 surftemp measurements exist from Newfoundland (fall only), this is not enough to calibrate a relationship. Thus, use Hadley-NOAA SST for all measurements in this region.
# For all other regions, use regression with sst*region to translate Hadley-NOAA SST into "surftemp".
########################################################


library(lattice)

# Read in file with HadlSST values for all hauls.
allSST<-read.csv("data/allhauls_withHadlSST_150527.csv",row.names=1)

# Load dat to later merge with new SST.
load('data/trawl_allregionsforprojections_2015-02-02.RData') # load dat data.frame. 

# Correct some surftemp measurements before regression 
# plot(allSST$surftemp,allSST$sst) #obvious issues with surftemp = 0. Many of these also miss bottomtemp, but still an issue. Maybe some of them are true zeros, but not many.
# Since these are now coded as NA, new values will be filled-in from HADL-NOAA dataset.
allSST$surftemp[allSST$surftemp==0]<-NA
#allSST$bottemp[allSST$bottemp>40 &is.na(allSST$bottemp)==F]<-NA
allSST$surftemp[allSST$surftemp>38 &is.na(allSST$surftemp)==F]<-NA

# Use complete cases for regression
allSSTc<-allSST[complete.cases(allSST[,c("sst","surftemp")]),]

########################################################
# Regress trawl surftemp on Hadl-NOAA SST
########################################################

#All regions
lmAll<-lm(surftemp~sst*region,data=allSSTc)
# For WestCoast only:
lmWC<-lm(surftemp~sst,data=allSSTc[allSSTc$region=="AFSC_WCTri",]) 

########################################################
# Use regressions to bias-correct Hadl-Noaa SST
########################################################

# Use regional/survey regressions to "correct" regional data, except for Newfoundland and WC. 
allSST$adjSST[!allSST$region %in% c("AFSC_WCTri", "NWFSC_WCAnn", "DFO_NewfoundlandFall", "DFO_NewfoundlandSpring")] <- predict(lmAll, newdata=allSST[!allSST$region %in% c("AFSC_WCTri", "NWFSC_WCAnn", "DFO_NewfoundlandFall", "DFO_NewfoundlandSpring"),])

# Use lmWC to bias-correct NWFSC_WCAnn (and WCTri) SST values
allSST$adjSST[allSST$region %in% c("AFSC_WCTri", "NWFSC_WCAnn")] <- predict(lmWC,newdata=allSST[allSST$region %in% c("AFSC_WCTri", "NWFSC_WCAnn"),])

# Use Hadley-NOAA SST directly for all Newfoundland hauls, and don't use any measured surftemp (n=7)
allSST$adjSST[allSST$region %in% c("DFO_NewfoundlandFall", "DFO_NewfoundlandSpring")] <- allSST$sst[allSST$region %in% c("DFO_NewfoundlandFall", "DFO_NewfoundlandSpring")]

########################################################
# Some simple diagnostics
########################################################

# cor(allSST$surftemp,allSST$adjSST,use="c")
# [1] 0.9759759

plot(allSST$sst,allSST$adjSST,col=rainbow(length(unique(allSST$region)))[allSST$region],xlab="Hadl-NOAA SST",ylab="Bias-corrected SST",pch=16)
plot(allSST$surftemp,allSST$adjSST,col=rainbow(length(unique(allSST$region)))[allSST$region],xlab="Survey surftemp",ylab="Bias-corrected SST",pch=1,cex=0.5)

# A few data points from SoGulf and ScotianShelf have  high SST and very low surftemp. 
anoms<-allSST[allSST$surftemp < 5 & allSST$adjSST > 11 & is.na(allSST$surftemp)==F & is.na(allSST$adjSST)==F,]
# Lots of subsetting and plotting (not shown) suggests that these are highly unlikely to be true observations (more likely transcription errors)

########################################################
# Create a new surftemp2 column to merge with dat file:
########################################################

# Default is:
allSST$surftemp2 <- allSST$surftemp
# If Newfoundland, use adjSST (HADL-NOAA)
allSST$surftemp2[allSST$region %in% c("DFO_NewfoundlandFall", "DFO_NewfoundlandSpring")] <- allSST$adjSST[allSST$region %in% c("DFO_NewfoundlandFall", "DFO_NewfoundlandSpring")]
# If WCAnn, use adjSST (HADL-NOAA)
allSST$surftemp2[allSST$region == "NWFSC_WCAnn"] <- allSST$adjSST[allSST$region == "NWFSC_WCAnn"]
# If otherwise missing, using adjSST (HADL-NOAA)
allSST$surftemp2[is.na(allSST$surftemp)==T] <- allSST$adjSST[is.na(allSST$surftemp)==T]
# If an unreasonable value (esp SoGulf and ScotianShelf), use adjSST (HADL-NOAA)
allSST$surftemp2[allSST$surftemp < 5 & allSST$adjSST > 11 & is.na(allSST$surftemp)==F & is.na(allSST$adjSST)==F] <- allSST$adjSST[allSST$surftemp < 5 & allSST$adjSST > 11 & is.na(allSST$surftemp)==F & is.na(allSST$adjSST)==F]

########################################################
# Merge allSST with dat file to create a new column
########################################################

#Haul IDs got tweaked in allSST - WCAnn were rounded. 

allhauls<-unique(dat$haulid) #114879
# Since order is preserved,
# replace haulid with allhauls - test below with full dataframes to be sure data are correctly matched
allSST$haulid<-allhauls

#change colnames in newSST (lat/lon were altered to merge with Hadl-NOAA SST)
colnames(allSST)[3:4]<-c("lon_new","lat_new")
colnames(allSST)[10:11]<-c("lat","lon")

#test<-merge(dat,allSST,by="haulid",sort=F,all=T)
#plot(test$surftemp.x,test$surftemp.y) # perfect. order was preserved.
#rm(test)

allSST2<-allSST[,c("haulid","surftemp2")]
dat2<-merge(dat,allSST2,by="haulid",sort=F,all=T)

########################################################
# Rename columns so that we don't have to re-write code.
########################################################

dat2$surftemp_orig<-dat2$surftemp
dat2$surftemp<-dat2$surftemp2
dat2<-dat2[,colnames(dat2)!="surftemp2"]
dat<-dat2
########################################################
# Save new dat file with surftemp
########################################################

save(dat,file=paste("data/trawl_allregionsforprojections_wSST_",Sys.Date(),".RData",sep=""))
