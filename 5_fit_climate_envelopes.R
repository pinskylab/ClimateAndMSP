## Set working directories
if(Sys.info()["nodename"] == "pinsky-macbookair"){
	setwd('~/Documents/Rutgers/Range projections/proj_ranges/')
	}
if(Sys.info()["nodename"] == "amphiprion.deenr.rutgers.edu"){
	setwd('~/Documents/range_projections/')
	}
# could add code for Lauren's working directory here


# Fit climate-envelope models
## Still very much in progress

#################
### Load data ###	
#################
load('data/trawl_allregionsforprojections_2015-02-02.RData') # load dat data.frame. Has all trawl observations from all regions. wtcpue has the standardized biomass estimates. They are standardized within regions, but not across regions.

source("CSquareCode.r") #taken from VMStools in googlecode

############################################
# Standardize species names across regions #
############################################
# see spptaxonomy_2015-02-09_plusManual.csv for a useful conversion table
Spptax<-read.csv("data/spptaxonomy_plusManual.csv") #note: new column in CSV file 
spptax<-apply(Spptax,2,tolower)
dat$sppl<-tolower(dat$spp)
datspp<-unique(dat$sppl)

spptax<-as.data.frame(spptax)
sum(datspp %in% spptax$taxon)# 798 of 4937 spp matched

#Match when possible and assign new genus species, otherwise keep old name
#Find taxa with matches in the taxonomy table
#In some cases, the "name" field of taxonomy has spelling errors (missing "n"s). In this case, genus+species should be used. But in some cases, "genus species" are not given if full taxonomy is missing. The spptaxonomy.csv file now has a "newname" column which merges these optimally.
matches<-datspp %in% spptax$taxon + datspp %in% spptax$name + datspp %in% spptax$common + datspp %in% spptax$genusspecies  #826 matches/4937

#For those species datspp[matches>0], replace current name with spptax$newname
newnames=character(length(datspp))
for (i in 1:length(datspp)){
	if(matches[i]==0){
		newnames[i]<-datspp[i]
	} else if(datspp[i] %in% spptax$taxon){
		ind<-match(datspp[i],spptax$taxon)
		newnames[i]<-as.character(spptax$newname[ind])
	} else if(datspp[i] %in% spptax$name){
		ind<-match(datspp[i],spptax$name)
		newnames[i]<-as.character(spptax$newname[ind])
	} else if(datspp[i] %in% spptax$genusspecies){
		ind<-match(datspp[i],spptax$genusspecies)
		newnames[i]<-as.character(spptax$newname[ind])
	} else if(datspp[i] %in% spptax$common){ # is this needed? (Malin 5/4/2015)
		ind<-match(datspp[i],spptax$common)
		newnames[i]<-as.character(spptax$newname[ind])
	}
}

length(unique(datspp))
length(unique(newnames)) # down to 4831 unique spp, from 4937

oldnewnames<-cbind(sppl=datspp,sppnew=newnames)
#Now merge new names with old names in dat file
dat<-merge(dat,oldnewnames) #dat now contains "sppnew"

######################
# Add useful columns #
######################

# Have a biomass that never goes to zero (useful for fitting log-links with a stratum effect) 
#dat$wtcpuena = dat$wtcpue
#dat$wtcpuena[dat$wtcpuena == 0] = 1e-4
#dat$wtcpuenal = log(dat$wtcpuena)

# other useful columns
dat$presfit = dat$wtcpue > 0 # indicator for where present
#dat$stratumfact = as.factor(dat$stratum)
#dat$yrfact = as.factor(dat$year)
#dat$regionfact<-as.factor(dat$region)
#dat$sppregion = paste(dat$sppnew, dat$region, sep='_')
dat$ocean[dat$region %in% c("AFSC_Aleutians", "AFSC_EBS", "AFSC_GOA", "AFSC_WCTri", "NWFSC_WCAnn")] <- "Pac"
dat$ocean[dat$region %in% c("DFO_NewfoundlandFall", "DFO_NewfoundlandSpring", "DFO_ScotianShelf","DFO_SoGulf","NEFSC_NEUSFall", "NEFSC_NEUSSpring","SEFSC_GOMex")] <- "Atl"
#dat$ocean[dat$region == "SEFSC_GOMex"] <- "Gulf" #Or should Gulf of Mex group with Altantic?
dat$sppocean = paste(dat$sppnew, dat$ocean, sep='_') 

# add rugosity
rugfile<-read.csv("data/trawl_latlons_rugosity_forMalin_2015_02_10.csv")
dat<-merge(dat,rugfile) #lose 69 instances of lumpenus lampretaeformis b/c missing lat/lon
rm(rugfile)
dat$logrugosity<-log(dat$rugosity+0.01) #log-transformed rugosity gave better model fits in initial tests of rugosity covariate

dat$cs1<-as.factor(CSquare(dat$lon,dat$lat,1)) #classify into 1 degree squares
#dat$cs6m<-as.factor(CSquare(dat$lon,dat$lat,0.1)) #classify into 6 arcminute squares

##########################
## Pick the spp to model #
##########################

# Perhaps also a minimum #years, or obs/year cutoff?

# Identify taxa caught in at least 10 survey years and at least 300 total observations > 0 in a given survey (hopefully avoiding changes in species classification through time).
# Species are identified by species+ocean for later trimming of dat.
# Some surveys have zero hauls, but not so systematically: these are very small catches (<1kg)

myregions<-unique(dat$region)
myspp<-NULL
for(r in myregions){
#	surveyyrs<-names(table(dat$year[dat$region==r]))
	regdat<-dat[dat$region==r & dat$wtcpue > 0 & !is.na(dat$wtcpue),] # trim to focal region and rows with catch > 0
	yrocc<-table(regdat$year,regdat$sppocean) #table of number of occurrances each year
	sumyrs<-apply(yrocc,2,function(x) sum(x>0)) #identify years with more than zero catch of each taxon
	sumobs<-colSums(yrocc)
	min1<-colnames(yrocc)[sumyrs>=10 & sumobs >= 300]  #identify taxa with at least 8 years of catch and at least 40 total catch records
	myspp<-c(myspp,min1) #now with myspp plus ocean
}
#table(myspp) #shows which species are selected in multiple surveys
myspp<-unique(myspp) 
length(myspp) # 663 unique taxa_ocean (up from 621 if we require presence in all years of a survey)

#remove  any taxa missing species name, egg cases, anemones, families
drop<-myspp[grep("spp",myspp)] #158 with "spp." (missing species name)
#drop<-c(drop,myspp[grep("egg",myspp)]) #7 with "egg"
#drop<-c(drop,myspp[grep("anemone",myspp)]) #2 anemones
# drop<-c(drop,myspp["teuthida","liparidinae","bathylagus sp.","lampanyctus sp.","caridea", "carinariidae","antipatharia","annelida" ,"crustacea shrimp") # also missing species name
droplist<-c("sp.", "egg","anemone","unident", "unknown", "teuthida","liparidinae","bathylagus sp.","lampanyctus sp.","caridea", "carinariidae","antipatharia","annelida" ,"crustacea shrimp", "artediellus _Atl", "asteroidea s.c._Atl", "bivalvia c._Atl", "cephalopoda c._Atl", "decapoda o._Atl", "gastropoda o._Atl", "mytilidae f._Atl", "ophiuroidea s.c._Atl", "paguridae f._Atl", "paguroidea s.f._Atl", "penaeus _Atl", "polychaeta c._Atl", "porifera p._Atl", "pycnogonida s.p._Atl", "scyphozoa c._Atl", "sepiolodae f._Atl", "anthozoa c._Atl", "beryciformes (order)_Atl", "cirripedia s.c._Atl", "clypeasteroida o._Atl", "etropus _Atl", "euphausiacea o._Atl", "fistularia _Atl", "halichondria cf. sitiens_Pac", "holothuroidea c._Atl", "macrouriformes (order)_Atl", "melanostomiidae (stomiatidae)_Atl", "octopoda o._Atl", "pandalidae f._Atl", "paralepididae _Atl", "seapen (order)_Atl", "symphurus _Atl", "tunicata s.p._Atl", "urophycis _Atl", "melanostomiidae (stomiatidae)_Atl", "gorgonocephalidae,asteronychidae f._Atl", "trash species in catch_Atl")
for (i in 1:length(droplist)) {drop<-c(drop,myspp[grep(droplist[i],myspp, fixed=TRUE)])}


# remove any taxa that lack a space (e.g., only one word, like a Family or Order name)
drop2 <- myspp[!grepl(' ', myspp)] # 92 more taxa to remove

myspp<-myspp[!(myspp %in% drop | myspp %in% drop2)]
myspp<-sort(myspp) 
length(myspp) #down to 551 spp. to model (but up from 512 if we require presence in all years). A few species are present in multiple ocean basins and thus will have multiple models.

# However, we'll need to go back to the full dat file to determine where zero hauls occur, unless we assume each haul had at least one of these species. Check this here:
#dat2<-dat[dat$sppocean %in% myspp,]
#length(unique(dat$haulid))  #114766
#length(unique(dat2$haulid)) #114766

dat<-dat[dat$sppocean %in% myspp,] 
#remove unnecessary columns:
dat<-dat[,!colnames(dat) %in% c("spp","sppl")] #could probably remove others

save(dat,file="data/dat_selectedspp.Rdata")




