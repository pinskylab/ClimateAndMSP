## Script to compare species distribution projections for the present against trawl data in MPAs

## Set working directories depending on computer
if(Sys.info()["nodename"] == "pinsky-macbookair"){
	setwd('~/Documents/Rutgers/Range projections/proj_ranges/')
	modfolder <- '../CEmodels_nonstationaritytest/'
	nthreads=2
	}
if(Sys.info()["nodename"] == "amphiprion.deenr.rutgers.edu"){
	setwd('~/Documents/range_projections/')
	.libPaths(new='~/R/library/') # so that it can find my old packages
	modfolder <- 'CEmodels_nonstationaritytest/'
	nthreads=10
	}
if(Sys.info()["user"] == "lauren"){
	setwd('~/backup/NatCap/proj_ranges/')
	modfolder <- 'output/CEmodels_nonstationaritytest/'
	nthreads=1
}

## Load packages
require(maptools)
require(rgdal)
require(rgeos)
require(data.table)
source('standardizeSppNames.r') # to fix species names across regions

## Read in data
load('data/trawl_allregionsforprojections_wSST_2015-06-02.RData') # trawl data. in 'dat'
wdpa = readShapePoly('cmsp_data/WDPA/NA_marine_MPA/mpinsky-search-1382225374362.shp', proj4string = CRS('+proj=longlat +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +no_defs')) # load the MPA data
	wdpa = wdpa[wdpa$marine == 1,] # trim to only marine (remove 20 of 3023)
	wdpa@data$wdpapolyID <- as.numeric(sapply(wdpa@polygons, FUN=slot, name='ID')) # extract polygon ID
wdpacov <- read.csv('data/wdpa_cov_by_grid0.25.csv', row.names=1) # shows wich MPAs are in which grid cells. each line is a unique grid cell-MPA combination.
	wdpacov <- as.data.table(wdpacov)
	wdpacov[wdpacov$lon<0, lon:=lon+360] # convert lon to positive

load('data/presmap_fitallreg_xreg_rcp45.RData') # loads 'presmap' with pres/abs from the ensemble projections
	presmap <- presmap[presmap$period=='2006-2020',] # trim to the 'present'
	

# Standardize species names
oldnewnames <- standardizeSppNames(dat$spp)
dat <- merge(dat, oldnewnames) # takes a couple minutes
	head(dat)
dat$ocean[dat$region %in% c("AFSC_Aleutians", "AFSC_EBS", "AFSC_GOA", "AFSC_WCTri", "NWFSC_WCAnn")] <- "Pac"
dat$ocean[dat$region %in% c("DFO_NewfoundlandFall", "DFO_NewfoundlandSpring", "DFO_ScotianShelf","DFO_SoGulf","NEFSC_NEUSFall", "NEFSC_NEUSSpring","SEFSC_GOMex")] <- "Atl"
#dat$ocean[dat$region == "SEFSC_GOMex"] <- "Gulf" #Or should Gulf of Mex group with Altantic?
dat$sppocean = paste(dat$sppnew, dat$ocean, sep='_') 


## Make points from trawl hauls
datu <- dat[!duplicated(dat$haulid) & !is.na(dat$lat) & !is.na(dat$lon),]
	dim(datu)
SP <- SpatialPointsDataFrame(datu[,c('lon', 'lat')], data=datu[,c('haulid', 'year', 'lat', 'lon')], proj4string=CRS(proj4string(wdpa)))
	# plot(SP, pch=16, cex=0.02, axes=TRUE); abline(h=c(25,65)) # all the points


## Assemble species list for each MPA from the trawl data
ov <- over(wdpa, SP, returnList=TRUE) # extract which haulids fall in which MPAs. takes a couple minutes
numhauls <- sapply(ov, FUN=nrow) # calculate number of hauls in each MPA.
inds <- which(numhauls > 0) # MPAs with >0 hauls
sppinMPAs <- data.frame(wdpapolyID=numeric(0), sppocean=character(0), presobs=logical(0)) # to hold the lists of species
for(i in 1:length(inds)){
	thesespp <- sort(unique(dat$sppocean[dat$haulid %in% ov[[inds[i]]]$haulid]))
	temp <- data.frame(wdpapolyID=wdpa@data$wdpapolyID[inds[i]], sppocean=thesespp, presobs=TRUE)
	sppinMPAs <- rbind(sppinMPAs, temp)
}
	length(unique(sppinMPAs$wdpapolyID)) # 127 MPAs with data
	sppinMPAs <- as.data.table(sppinMPAs)
	

## Assemble species list in each MPA from the projections data
	# merge MPA ID into species data (using data.tables)
setkey(wdpacov, lat, lon)
presmap <- as.data.table(presmap)
setkey(presmap, lat, lon)
presmap2 <- wdpacov[presmap,.(sppocean, pres, wdpapolyID, lat, lon), allow.cartesian=TRUE] # datatable join, using the keys set in each dt. the order specifies to keep all rows of thesesppbymod. allow.cartesian is needed because the result is larger than either parent
	dim(presmap2)
	
	# summarize by unique species in each MPA 
pressppbyMPA <- presmap2[!is.na(wdpapolyID),.(prespred = as.logical(max(pres))),by="sppocean,wdpapolyID"] # DT aggregate function
	dim(pressppbyMPA) # 137796

## Compare predictions and observations
m1 <- pressppbyMPA[wdpapolyID %in% sppinMPAs$wdpapolyID,] # trim predictions to MPAs with observations
m2 <- sppinMPAs[wdpapolyID %in% pressppbyMPA$wdpapolyID & sppocean %in% pressppbyMPA$sppocean,] # trim observations to MPAs with predictions AND to species with predictions
comp <- merge(m1, m2, by=c('wdpapolyID', 'sppocean'), all.x=TRUE, all.y=TRUE)

	# check that predictions and observations exist for each MPA to be compared
	wdpas <- sort(unique(comp$wdpapolyID))
	for(i in 1:length(wdpas)){
		inds <- comp$wdpapolyID==wdpas[i]
		if(sum(!is.na(comp$prespred[inds]))==0) print(paste('Missing predictions for wdpa', wdpas[i])) 
		if(sum(!is.na(comp$presobs[inds]))==0) print(paste('Missing observations for wdpa', wdpas[i])) 
	}
	
comp[is.na(presobs), presobs:=FALSE] # assume spp is missing if not observed
length(unique(comp$wdpapolyID)) # 121 MPAs with data

	# check that both Pac and Atl species don't appear in the same MPA
	wdpas <- sort(unique(comp$wdpapolyID))
	for(i in 1:length(wdpas)){
		inds <- comp$wdpapolyID==wdpas[i]
		if(any(grepl('Pac', comp$sppocean[inds])) & any(grepl('Atl', comp$sppocean[inds]))) print(paste('Both Atl and Pac in wdpa', wdpas[i]))
	}

	# write out
write.csv(comp, file='output/MPA_vs_trawl_comparison.csv')

## Examine comparisons
comp <- read.csv('output/MPA_vs_trawl_comparison.csv')

i=2
table(list(pred=comp[wdpapolyID==wdpas[i],prespred], obs=comp[wdpapolyID==wdpas[i],presobs])) # examine one MPA

	# commission and ommission error rates
wdpas <- sort(unique(comp$wdpapolyID))
rates <- data.frame(wdpapolyID=wdpas, truepospred=NA, falsepospred=NA, truenegpred=NA, falsenegpred=NA, comm=NA, omm=NA)
for(i in 1:length(wdpas)){
	inds <- comp$wdpapolyID==wdpas[i]
	rates$truepospred[i] <- with(comp[inds,], sum(prespred & presobs)) # misclassication of absences
	rates$falsepospred[i] <- with(comp[inds,], sum(prespred & !presobs)) # misclassication of absences
	rates$truenegpred[i] <- with(comp[inds,], sum(!prespred & !presobs)) # misclassication of absences
	rates$falsenegpred[i] <- with(comp[inds,], sum(!prespred & presobs)) # misclassication of absences
	rates$comm[i] <- with(comp[inds,], sum(prespred & !presobs)/sum(!presobs)) # misclassication of absences
	rates$omm[i] <- with(comp[inds,], sum(!prespred & presobs)/sum(presobs)) # misclassification of presences
}

	# examine error rates vs. size
rates <- merge(rates, wdpa@data[,c('wdpapolyID', 'rep_m_area', 'rep_area')], all.x=TRUE)

	with(rates, plot(comm ~ I(rep_m_area+0.1), log='x'))
	with(rates, plot(comm ~ I(rep_area+0.1), log='x'))
	with(rates, plot(omm ~ I(rep_m_area+0.1), log='x'))
	with(rates, plot(omm ~ I(rep_area+0.1), log='x'))
	
	# examine error rates vs. # trawl hauls
nhauls <- data.frame(wdpapolyID=names(numhauls), nhauls=numhauls)
rates <- merge(rates, nhauls, all.x=TRUE)
	
	with(rates, plot(comm ~ nhauls, log='x'))
	with(rates, plot(omm ~ nhauls, log='x'))
	
	# models
mod1 <- glm(cbind(falsepospred, falsepospred+truenegpred) ~ rep_m_area*nhauls, family=binomial, data=rates)
	summary(mod1)
mod2 <- glm(cbind(falsenegpred, falsenegpred+truepospred) ~ rep_m_area*nhauls, family=binomial, data=rates)
	summary(mod2)
