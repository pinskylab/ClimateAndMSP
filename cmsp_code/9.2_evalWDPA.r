# evaluate protected areas

## Set working directories
if(Sys.info()["nodename"] == "pinsky-macbookair"){
	setwd('~/Documents/Rutgers/Range projections/proj_ranges/')
	presmapbymodfolder <- '../data/'
	}
if(Sys.info()["nodename"] == "amphiprion.deenr.rutgers.edu"){
	setwd('~/Documents/range_projections/')
	.libPaths(new='~/R/x86_64-redhat-linux-gnu-library/3.1/') # so that it can find my old packages
	presmapbymodfolder <- 'data/'
	}
# could add code for Lauren's working directory here

############
## Flags
############

## choose which runs to use
## runtype refers to how the Species Distribution Models (SDMs) were fit
## projtype refers to how the SDM projections were done
#runtype <- 'test'; projtype=''
#runtype <- ''; projtype=''`
#runtype <- 'testK6noSeas'; projtype='_xreg'
runtype <- 'fitallreg'; projtype='_xreg'

# choose the rcp
rcp <- 85
otherrcp <- 45

# select initial and final timeperiod for these grids
periods <- c('2006-2020', '2081-2100')

####################
## helper functions
####################
require(Hmisc)
require(data.table)
require(lme4) # for mixed-effects models
require(car) # for testing ME models


lu <- function(x) return(length(unique(x)))

se <- function(x,na.rm=FALSE){ # standard error
	if(!na.rm){
		return(sd(x, na.rm=FALSE)/sqrt(length(x)))
	}
	if(na.rm){
		return(sd(x, na.rm=TRUE)/sqrt(sum(!is.na(x))))
	}
}

# weighted mean function to use with summarize()
wmean <- function(x){ # values in col 1, weights in col 2
	inds <- !is.na(x[,1]) & !is.na(x[,2])
	return(weighted.mean(x[inds,1], x[inds,2])) 
}

# expects x to have columns 'periodmids' and 'rich' (in that order)
calcrichtrendmids <- function(x){
	mod <- lm(x[,2] ~ x[,1])
	return(mod$coefficients[2])
}

# turnover metrics
# expects region in column 1, lat (2), lon (3), time period in column 4, spp in col 5, pres in col 6
# periods lists the periods that can appear in col4 (in order from earliest to latest)
# only present species included
turnover <- function(x, periods){
	pers <- sort(unique(x[,4]))
	if(length(periods) != 2) stop(paste('expecting 2 time periods:', paste(unique(periods, collapse=','))))
	initcomm <- x[x[,4] == periods[1] & x[,6]==TRUE,5]
	finalcomm <- x[x[,4] == periods[2] & x[,6]==TRUE,5]
	
	a <- length(intersect(initcomm, finalcomm)) # number of shared species
	lost <- length(setdiff(initcomm, finalcomm)) # number of species lost
	gained <- length(setdiff(finalcomm, initcomm))
	return(c(nstart=length(initcomm), nend=length(finalcomm), nlost=lost, ngained=gained, flost=lost/length(initcomm), fgained=gained/length(initcomm), beta_sor=2*a/(2*a+gained+lost)))
}

# same functions, but for use with data.table
nstart <- function(period, spp, pres, periods){
	pers <- sort(unique(period))
	if(length(periods) != 2) stop(paste('expecting 2 time periods:', paste(unique(periods, collapse=','))))
	initcomm <- spp[period == periods[1] & pres==TRUE]
	return(length(initcomm))
}

nend <- function(period, spp, pres, periods){
	pers <- sort(unique(period))
	if(length(periods) != 2) stop(paste('expecting 2 time periods:', paste(unique(periods, collapse=','))))
	finalcomm <- spp[period == periods[2] & pres==TRUE]
	return(length(finalcomm))
}

nlost <- function(period, spp, pres, periods){
	pers <- sort(unique(period))
	if(length(periods) != 2) stop(paste('expecting 2 time periods:', paste(unique(periods, collapse=','))))
	initcomm <- spp[period == periods[1] & pres==TRUE]
	finalcomm <- spp[period == periods[2] & pres==TRUE]
	
	lost <- length(setdiff(initcomm, finalcomm)) # number of species lost
	return(lost)
}

ngained <- function(period, spp, pres, periods){
	pers <- sort(unique(period))
	if(length(periods) != 2) stop(paste('expecting 2 time periods:', paste(unique(periods, collapse=','))))
	initcomm <- spp[period == periods[1] & pres==TRUE]
	finalcomm <- spp[period == periods[2] & pres==TRUE]
	
	gained <- length(setdiff(finalcomm, initcomm))
	return(gained)
}

flost <- function(period, spp, pres, periods){
	pers <- sort(unique(period))
	if(length(periods) != 2) stop(paste('expecting 2 time periods:', paste(unique(periods, collapse=','))))
	initcomm <- spp[period == periods[1] & pres==TRUE]
	finalcomm <- spp[period == periods[2] & pres==TRUE]
	
	lost <- length(setdiff(initcomm, finalcomm)) # number of species lost
	return(lost/length(initcomm))
}

fgained <- function(period, spp, pres, periods){ # divide by number of spp in initial community
	pers <- sort(unique(period))
	if(length(periods) != 2) stop(paste('expecting 2 time periods:', paste(unique(periods, collapse=','))))
	initcomm <- spp[period == periods[1] & pres==TRUE]
	finalcomm <- spp[period == periods[2] & pres==TRUE]
	
	gained <- length(setdiff(finalcomm, initcomm))
	return(gained/length(initcomm))
}

fgainedalt <- function(period, spp, pres, periods){ # divide by number of spp in final community
	pers <- sort(unique(period))
	if(length(periods) != 2) stop(paste('expecting 2 time periods:', paste(unique(periods, collapse=','))))
	initcomm <- spp[period == periods[1] & pres==TRUE]
	finalcomm <- spp[period == periods[2] & pres==TRUE]
	
	gained <- length(setdiff(finalcomm, initcomm))
	return(gained/length(finalcomm))
}


beta_sor <- function(period, spp, pres, periods){
	pers <- sort(unique(period))
	if(length(periods) != 2) stop(paste('expecting 2 time periods:', paste(unique(periods, collapse=','))))
	initcomm <- spp[period == periods[1] & pres==TRUE]
	finalcomm <- spp[period == periods[2] & pres==TRUE]
	
	a <- length(intersect(initcomm, finalcomm)) # number of shared species
	lost <- length(setdiff(initcomm, finalcomm)) # number of species lost
	gained <- length(setdiff(finalcomm, initcomm))
	return(2*a/(2*a+gained+lost))
}

# pick the most common item in a list, or the first if a tie
pickone <- function(x){
	tab <- table(x)
	pick <- which.max(tab)
	return(names(tab)[pick])
}

###########################
## Load and set up data
###########################

wdpa <- read.csv('data/wdpa_cov_by_grid0.25.csv', row.names=1) # shows wich MPAs are in which grid cells. each line is a unique grid cell-MPA combination.
	wdpa$lon[wdpa$lon<0] <- wdpa$lon[wdpa$lon<0] + 360 # convert lon to positive

	wdpa <- as.data.table(wdpa)

	# set up MPA networks
	wdpa[sub_loc == 'US-CA' & mang_auth == 'California Department of Fish and Game',network:='mlpa']
	wdpa[(sub_loc %in% c('US-DE', 'US-FL', 'US-GA', 'US-MA', 'US-MD', 'US-ME', 'US-NC', 'US-NJ', 'US-NY', 'US-RI', 'US-SC', 'US-VA')) & (lon > 280), network:='eastcoast'] # by choosing by state, this won't include federal water closures
	wdpa[(sub_loc %in% c('US-AK')), network:='ak'] # by choosing by state, this won't include federal water closures
	
	# choose one and only one region per MPA, assing to a new column
	wdpa[,region.one:=pickone(region),by=wdpapolyID]

	# calculate mpa extent
	wdpa[, ':=' (lat_min=min(lat), lat_max=max(lat)), by=wdpapolyID]

#	wdpanets[network=='mlpa',plot(lon,lat)]
#	wdpanets[network=='eastcoast',plot(lon,lat)]
#	wdpanets[network=='ak',plot(lon,lat)]
#	require(maps); map(database='world2',add=TRUE)


# load ensemble mean distributions (both rcps)
#load(paste('data/presmap_', runtype, projtype, '_rcp', rcp, '.RData', sep=''))
#	presmap.1 <- presmap
#	presmap.1$rcp <- rcp
#	dim(presmap.1)
#load(paste('data/presmap_', runtype, projtype, '_rcp', otherrcp, '.RData', sep=''))
#	presmap$rcp <- otherrcp
#	dim(presmap)
#	presmap <- rbind(presmap.1, presmap)
#	dim(presmap)
#	rm(presmap.1)
#
# load presmap for all model runs, a bit slow (large files)
load(paste(presmapbymodfolder, 'presmapbymod_', runtype, projtype, '_rcp', rcp, '.RData', sep='')) # loads presmap data.frame with presence/absence information from each model (slow to load).
	presmapbymod.1 <- as.data.table(presmapbymod) # seems to remove presmapbymod as well?

	setkey(presmapbymod.1,lat, lon)
		dim(presmapbymod.1)
	presmapbymod.2 <- presmapbymod.1[.(wdpa$lat, wdpa$lon),nomatch=0] # trim to lat lon rows relevant to WDPA analysis. need nomatch=0 so that latlon values not in the data.table return 0 rows instead of a row with NA
	setkey(presmapbymod.2,period,pres)
	presmapbymod.2 <- presmapbymod.2[.(periods,TRUE),nomatch=0] # trim to periods and pres=TRUE rows
		presmapbymod.2[,sum(is.na(pres))] # should be zero
		dim(presmapbymod.2)
	presmapbymod.2[,rcp:=rcp] # add label for rcp

	rm(presmapbymod.1)

load(paste(presmapbymodfolder, 'presmapbymod_', runtype, projtype, '_rcp', otherrcp, '.RData', sep='')) # loads presmap data.frame with presence/absence information from each model (slow to load) (for the other rcp, the one not used for planning)
	presmapbymod <- as.data.table(presmapbymod)

	setkey(presmapbymod,lat,lon)
		dim(presmapbymod)
	presmapbymod.3 <- presmapbymod[.(wdpa$lat, wdpa$lon),nomatch=0] # trim to rows relevant to WDPA analysis
	setkey(presmapbymod.3,period,pres)
	presmapbymod.3 <- presmapbymod.3[.(periods,TRUE),nomatch=0] # trim to rows relevant to WDPA analysis
		presmapbymod.3[,sum(is.na(pres))] # should be zero
		dim(presmapbymod.3)
	presmapbymod.3[,rcp:=otherrcp]

# bind the presmaps together
	presmapbymod <- rbind(presmapbymod.2, presmapbymod.3)
		nrow(presmapbymod.2)
		nrow(presmapbymod.3)
		dim(presmapbymod)

	rm(presmapbymod.2)
	rm(presmapbymod.3)

# trim presmap to species present, locations and time periods relevant to wdpa analysis
#	dim(presmapbymod) # 189,871,500
#	setkey(presmapbymod,period,pres) # data.tables!
#thesesppbymod <- presmapbymod[list(periods,TRUE)] # quick, using binary search
#	rm(presmapbymod) # clean up memory

#thesesppbymod <- thesesppbymod[paste(thesesppbymod$lat, thesesppbymod$lon) %in% paste(wdpa$lat, wdpa$lon),] # "slow" vector search
#	dim(thesesppbymod) # 6,010,462

	regs <- sort(unique(presmapbymod$region))
	regstokeep <- setdiff(regs, c('NEFSC_NEUSFall', 'AFSC_WCTri', 'DFO_NewfoundlandSpring')) # remove some surveys to avoid duplicates in a region
	setkey(presmapbymod, region)
	thesesppbymod <- presmapbymod[.(regstokeep)] # trim out duplicate regions
	dim(thesesppbymod) # 5,183,120? 7,601,343

	rm(presmapbymod)

##############################################
# Calculate turnover within whole MPAs
# Use ensemble mean distributions
##############################################
#
## trim presmap to species present, locations and time periods relevant to wdpa analysis
#thesespp <- presmap[paste(presmap$lat, presmap$lon) %in% paste(wdpa$lat, wdpa$lon) & presmap$period %in% periods & presmap$pres,]
#	dim(thesespp) # 
#	sort(unique(thesespp$region))
#	thesespp <- thesespp[!(thesespp$region %in% c('NEFSC_NEUSFall', 'AFSC_WCTri', 'DFO_NewfoundlandSpring')),] # trim out duplicate regions
#	dim(thesespp) # 
#	
## merge in MPA information
#thesespp <- merge(thesespp, wdpa[,c('wdpapolyID', 'lat', 'lon')], by=c('lat', 'lon'))
#	dim(thesespp) # 519,584
##	sort(unique(wdpa$wdpapolyID))
##	sort(unique(thesespp$wdpapolyID)) # lose some MPAs
#	
## summarize by unique species in each period in each MPA in each RCP (pres or not)
#thesesppbyMPA <- aggregate(list(pres=thesespp$pres), by=list(sppocean=thesespp$sppocean, period=thesespp$period, wdpapolyID=thesespp$wdpapolyID, rcp=thesespp$rcp), FUN=max)
#	thesesppbyMPA$pres <- as.logical(thesesppbyMPA$pres)
#	dim(thesesppbyMPA) # 110,588
#
## set up data.frame to hold turnover results
#wdpaturnbyMPA <- wdpa[!duplicated(wdpa$wdpapolyID) & wdpa$wdpapolyID %in% thesespp$wdpapolyID, c('wdpapolyID', 'area_wdpa', 'country', 'sub_loc', 'name', 'orig_name', 'desig', 'desig_eng', 'desig_type', 'iucn_cat', 'marine', 'rep_m_area', 'rep_area', 'status', 'status_yr', 'gov_type', 'mang_auth', 'int_crit', 'mang_plan', 'official', 'is_point', 'no_take', 'no_tk_area', 'metadata_i', 'action')] # list of MPAs to analyze
#	nrow(wdpaturnbyMPA) # 562
#
#ngrids <- aggregate(list(ngrids=wdpa$wdpapolyID), by=list(wdpapolyID=wdpa$wdpapolyID), FUN=length) # calculate number of grids per MPA
#wdpaturnbyMPA <- merge(wdpaturnbyMPA, ngrids)
#	nrow(wdpaturnbyMPA) # 562
#
#wdpaturnbyMPA$rcp <- rcp # mark the first rcp
#temp <- wdpaturnbyMPA # create a duplicate dataframe to hold the other rcp
#temp$rcp <- otherrcp #
#wdpaturnbyMPA <- rbind(wdpaturnbyMPA, temp) # merge the two together
#	nrow(wdpaturnbyMPA) # 1124
#
#	# loop through MPAs and rcps and analyze turnover (quick)
#wdpaturnbyMPA$nstart <- wdpaturnbyMPA$nend <- wdpaturnbyMPA$nlost <- wdpaturnbyMPA$ngained <- wdpaturnbyMPA$flost <- wdpaturnbyMPA$fgained <- wdpaturnbyMPA$beta_sor <- NA
#for(i in 1:nrow(wdpaturnbyMPA)){
#	if(i %% 10 == 0) print(paste(i, 'of', nrow(wdpaturnbyMPA), Sys.time()))
#	x<-thesesppbyMPA[thesesppbyMPA$wdpapolyID == wdpaturnbyMPA$wdpapolyID[i] & thesesppbyMPA$rcp == wdpaturnbyMPA$rcp[i],c('period', 'sppocean', 'pres')] # set up data.frame to hand to turnover()
#		x$region <- x$lat <- x$lon <- NA
#		x <- x[,c('region', 'lat', 'lon', 'period', 'sppocean', 'pres')] # re-order to what turnover() expects
#	wdpaturnbyMPA[i,c('nstart', 'nend', 'nlost', 'ngained', 'flost', 'fgained', 'beta_sor')] <- turnover(x, periods)
#}
#
## write out wdpaturnbyMPA
#write.csv(wdpaturnbyMPA, paste('data/wdpaturnbyMPA_', runtype, projtype, '_rcp', rcp, '&', otherrcp, '.csv', sep=''))
#

###############################
# Examine change within MPAs
# Ensemble mean results
##############################
#wdpaturnbyMPA <- read.csv(paste('data/wdpaturnbyMPA_', runtype, projtype, '_rcp', rcp, '&', otherrcp, '.csv', sep=''))
#
#ntk <- wdpaturnbyMPA$no_take %in% c('All', 'Part'); sum(ntk[wdpaturnbyMPA$rcp==rcp]) # no take reserves (index into wdpaturnbyMPA): 50 unique (43 in CA)
#
#
## Fraction of species lost
#	# all MPAs
#	mean(wdpaturnbyMPA$flost[wdpaturnbyMPA$rcp==rcp], na.rm=TRUE) # RCP 
#	sd(wdpaturnbyMPA$flost[wdpaturnbyMPA$rcp==rcp], na.rm=TRUE) # SD across MPAs
#	range(wdpaturnbyMPA$flost[wdpaturnbyMPA$rcp==rcp], na.rm=TRUE)
#
#	mean(wdpaturnbyMPA$flost[wdpaturnbyMPA$rcp==otherrcp], na.rm=TRUE) # other rcp
#	sd(wdpaturnbyMPA$flost[wdpaturnbyMPA$rcp==otherrcp], na.rm=TRUE)
#	range(wdpaturnbyMPA$flost[wdpaturnbyMPA$rcp==otherrcp], na.rm=TRUE)
#
#	# no-take
#	mean(wdpaturnbyMPA$flost[wdpaturnbyMPA$rcp==rcp & ntk], na.rm=TRUE) # RCP 
#	sd(wdpaturnbyMPA$flost[wdpaturnbyMPA$rcp==rcp & ntk], na.rm=TRUE) # SD across MPAs
#	range(wdpaturnbyMPA$flost[wdpaturnbyMPA$rcp==rcp & ntk], na.rm=TRUE)
#
#		t.test(asin(sqrt(wdpaturnbyMPA$flost[wdpaturnbyMPA$rcp==rcp & ntk])), asin(sqrt(wdpaturnbyMPA$flost[wdpaturnbyMPA$rcp==rcp & !ntk]))) # no take reserves projected to lose more species, p = 3e-7
#
#	mean(wdpaturnbyMPA$flost[wdpaturnbyMPA$rcp==otherrcp & ntk]) # other rcp
#	sd(wdpaturnbyMPA$flost[wdpaturnbyMPA$rcp==otherrcp & ntk])
#	range(wdpaturnbyMPA$flost[wdpaturnbyMPA$rcp==otherrcp & ntk])
#
#		t.test(asin(sqrt(wdpaturnbyMPA$flost[wdpaturnbyMPA$rcp==otherrcp & ntk])), asin(sqrt(wdpaturnbyMPA$flost[wdpaturnbyMPA$rcp==otherrcp & !ntk]))) # no take reserves projected to lose more species, p = 8e-11
#
#
## Fraction of species gained
#mean(wdpaturnbyMPA$fgained[is.finite(wdpaturnbyMPA$flost) & wdpaturnbyMPA$rcp==rcp])
#sd(wdpaturnbyMPA$fgained[is.finite(wdpaturnbyMPA$flost) & wdpaturnbyMPA$rcp==rcp])
#range(wdpaturnbyMPA$fgained[is.finite(wdpaturnbyMPA$fgained) & wdpaturnbyMPA$rcp==rcp])
#
#mean(wdpaturnbyMPA$fgained[is.finite(wdpaturnbyMPA$flost) & wdpaturnbyMPA$rcp==otherrcp])
#sd(wdpaturnbyMPA$fgained[is.finite(wdpaturnbyMPA$flost) & wdpaturnbyMPA$rcp==otherrcp])
#range(wdpaturnbyMPA$fgained[is.finite(wdpaturnbyMPA$fgained) & wdpaturnbyMPA$rcp==otherrcp])
#
## Sorenson similarity
#mean(wdpaturnbyMPA$beta_sor[is.finite(wdpaturnbyMPA$beta_sor) & wdpaturnbyMPA$rcp==rcp]) # RCP 
#sd(wdpaturnbyMPA$beta_sor[is.finite(wdpaturnbyMPA$beta_sor) & wdpaturnbyMPA$rcp==rcp]) # SD across MPAs
#
#mean(wdpaturnbyMPA$beta_sor[is.finite(wdpaturnbyMPA$beta_sor) & wdpaturnbyMPA$rcp==otherrcp]) # other rcp
#sd(wdpaturnbyMPA$beta_sor[is.finite(wdpaturnbyMPA$beta_sor) & wdpaturnbyMPA$rcp==otherrcp])
#
#	# no-take
#	mean(wdpaturnbyMPA$beta_sor[wdpaturnbyMPA$rcp==rcp & ntk], na.rm=TRUE) # RCP 
#	sd(wdpaturnbyMPA$beta_sor[wdpaturnbyMPA$rcp==rcp & ntk], na.rm=TRUE) # SD across MPAs
#	range(wdpaturnbyMPA$beta_sor[wdpaturnbyMPA$rcp==rcp & ntk], na.rm=TRUE)
#
#		t.test(asin(sqrt(wdpaturnbyMPA$beta_sor[wdpaturnbyMPA$rcp==rcp & ntk])), asin(sqrt(wdpaturnbyMPA$beta_sor[wdpaturnbyMPA$rcp==rcp & !ntk])))
#
#	mean(wdpaturnbyMPA$beta_sor[wdpaturnbyMPA$rcp==otherrcp & ntk]) # other rcp
#	sd(wdpaturnbyMPA$beta_sor[wdpaturnbyMPA$rcp==otherrcp & ntk])
#	range(wdpaturnbyMPA$beta_sor[wdpaturnbyMPA$rcp==otherrcp & ntk])
#
#		t.test(asin(sqrt(wdpaturnbyMPA$beta_sor[wdpaturnbyMPA$rcp==otherrcp & ntk])), asin(sqrt(wdpaturnbyMPA$beta_sor[wdpaturnbyMPA$rcp==otherrcp & !ntk])))
#
#
## Fraction lost by area
#with(wdpaturnbyMPA[wdpaturnbyMPA$rcp==rcp,], plot(ngrids, flost)) # by # grids
#	mod <- glm(cbind(nlost, nstart-nlost) ~ ngrids, family=binomial, data=wdpaturnbyMPA[wdpaturnbyMPA$rcp==rcp,])
#	summary(mod)
#with(wdpaturnbyMPA[wdpaturnbyMPA$rcp==rcp,], plot(area_wdpa, flost)) # by calculated area
#with(wdpaturnbyMPA[wdpaturnbyMPA$rcp==rcp,], plot(log(area_wdpa), flost, cex=nstart/100)) # by log(calculated area), scaled by # species
#	mod <- glm(cbind(nlost, nstart-nlost) ~ I(log(area_wdpa)), family=binomial, data=wdpaturnbyMPA[wdpaturnbyMPA$rcp==rcp,])
#	summary(mod)
#with(wdpaturnbyMPA[wdpaturnbyMPA$rcp==rcp,], plot(rep_m_area, flost)) # by marine area
#with(wdpaturnbyMPA[wdpaturnbyMPA$rcp==rcp,], plot(rep_area, flost)) # by full area
#with(wdpaturnbyMPA[wdpaturnbyMPA$rcp==rcp,], plot(pmin(rep_area, rep_m_area, na.rm=TRUE), flost)) # minimum of full or marine area
#	# no relationship
#
## Sorenson by area
#with(wdpaturnbyMPA[wdpaturnbyMPA$rcp==rcp,], plot(ngrids, beta_sor)) # by # grids
#	mod <- lm(I(asin(sqrt(beta_sor))) ~ ngrids, data=wdpaturnbyMPA[wdpaturnbyMPA$rcp==rcp,])
#	summary(mod)
#with(wdpaturnbyMPA[wdpaturnbyMPA$rcp==rcp,], plot(area_wdpa, beta_sor, cex=nstart/100)) # by calculated area
#	mod <- lm(I(asin(sqrt(beta_sor))) ~ area_wdpa, data=wdpaturnbyMPA[wdpaturnbyMPA$rcp==rcp,])
#	summary(mod)
#with(wdpaturnbyMPA[wdpaturnbyMPA$rcp==rcp,], plot(log(area_wdpa), beta_sor, cex=nstart/100)) # by log(calculated area), scaled by # species
#	mod <- lm(I(asin(sqrt(beta_sor))) ~ I(log(area_wdpa)), data=wdpaturnbyMPA[wdpaturnbyMPA$rcp==rcp,])
#	summary(mod)
#with(wdpaturnbyMPA[wdpaturnbyMPA$rcp==rcp,], plot(pmin(rep_area, rep_m_area, na.rm=TRUE), beta_sor)) # minimum of full or marine area
#	mod <- lm(I(asin(sqrt(beta_sor))) ~ I(pmin(rep_area, rep_m_area, na.rm=TRUE)), data=wdpaturnbyMPA[wdpaturnbyMPA$rcp==rcp,])
#	summary(mod)
#	# no relationship
#
## Fraction lost by # species
#with(wdpaturnbyMPA[wdpaturnbyMPA$rcp==rcp,], plot(nstart, flost, cex=area_wdpa/5)) # by # species, scaled by area
#
#
#	
## examine MLPA MPAs
#	wdpaturnave[wdpaturnave$sub_loc=='US-CA' & wdpaturnave$mang_auth=="California Department of Fish and Game",] # all rows
#	summary(wdpaturnave$beta_sor[wdpaturnave$sub_loc=='US-CA' & wdpaturnave$mang_auth=="California Department of Fish and Game"]) # Sorenson similarity
#	summary(wdpaturnave$flost[wdpaturnave$sub_loc=='US-CA' & wdpaturnave$mang_auth=="California Department of Fish and Game"]) # Fraction species lost
#	summary(wdpaturnave$fgained[wdpaturnave$sub_loc=='US-CA' & wdpaturnave$mang_auth=="California Department of Fish and Game"]) # Fraction species gained



##############################################
# Calculate turnover within whole MPAs
# Use individual climate models
##############################################
	
# merge MPA ID into species data (using data.tables)
setkey(wdpa, lat, lon)
setkey(thesesppbymod, lat, lon)
thesesppbymod2 <- wdpa[thesesppbymod,.(sppocean, period, rcp, model, pres, wdpapolyID, lat, lon, region.one), allow.cartesian=TRUE] # datatable join, using the keys set in each dt. the order specifies to keep all rows of thesesppbymod. allow.cartesian is needed because the result is larger than either parent
	dim(thesesppbymod2) # 10,015,534? 20,905,480, now 31740504 (4/23/2016)
	
# summarize by unique species in each period in each MPA in each RCP in each model (pres or not)
thesesppbyMPAbymod <- thesesppbymod2[,max(pres),by="sppocean,period,wdpapolyID,region.one,rcp,model"] # DT aggregate function
	thesesppbyMPAbymod[,pres:=as.logical(V1)]
	thesesppbyMPAbymod[,V1:=NULL] # drop V1
	dim(thesesppbyMPAbymod) # 2,214,439? now 2,326,903, then 2,214,439 (4/23/2015)

# calculate turnover by MPA, rcp, and model. each calc is pretty quick (a few seconds each)
thesesppbyMPAbymod[,nstart:=nstart(period, sppocean, pres, periods), by="wdpapolyID,rcp,model"]
thesesppbyMPAbymod[,nend:=nend(period, sppocean, pres, periods), by="wdpapolyID,rcp,model"]
thesesppbyMPAbymod[,nlost:=nlost(period, sppocean, pres, periods), by="wdpapolyID,rcp,model"]
thesesppbyMPAbymod[,ngained:=ngained(period, sppocean, pres, periods), by="wdpapolyID,rcp,model"]
thesesppbyMPAbymod[,flost:=flost(period, sppocean, pres, periods), by="wdpapolyID,rcp,model"]
thesesppbyMPAbymod[,fgained:=fgained(period, sppocean, pres, periods), by="wdpapolyID,rcp,model"]
thesesppbyMPAbymod[,fgainedalt:=fgainedalt(period, sppocean, pres, periods), by="wdpapolyID,rcp,model"]
thesesppbyMPAbymod[,beta_sor:=beta_sor(period, sppocean, pres, periods), by="wdpapolyID,rcp,model"]

# trim out duplicate rows within an mpa/rcp/model
setkey(thesesppbyMPAbymod, wdpapolyID,rcp,model)
wdpaturnbyMPAbymod <- unique(thesesppbyMPAbymod)
	dim(wdpaturnbyMPAbymod) # 14,612
	wdpaturnbyMPAbymod[,c('sppocean','period','pres'):=NULL]
	
# merge in mpa metadata
setkey(wdpa, wdpapolyID)
uwdpa <- unique(wdpa) # only the lines with unique mpas
	dim(uwdpa) # 634, 571 on 4/24/16
setkey(uwdpa, wdpapolyID)
setkey(wdpaturnbyMPAbymod, wdpapolyID,rcp,model)
wdpaturnbyMPAbymod <- uwdpa[wdpaturnbyMPAbymod, .(wdpapolyID, region.one, network, area_wdpa, lat_min, lat_max, country, sub_loc, name, orig_name, desig, desig_eng, desig_type, iucn_cat, marine, rep_m_area, rep_area, status, status_yr, gov_type, mang_auth, int_crit, mang_plan, official, is_point, no_take, no_tk_area, metadata_i, action,rcp,model,nstart,nend,nlost,ngained,flost,fgained,fgainedalt,beta_sor)]
	dim(wdpaturnbyMPAbymod) # 14,612


# write out wdpaturnbyMPA
write.csv(wdpaturnbyMPAbymod, paste('data/wdpaturnbyMPAbymod_', runtype, projtype, '_rcp', rcp, '&', otherrcp, '.csv', sep=''))



#############################################
# Examine turnover within MPAs
# Examine all climate models individually
#############################################
wdpaturnbyMPAbymod <- as.data.table(read.csv(paste('data/wdpaturnbyMPAbymod_', runtype, projtype, '_rcp', rcp, '&', otherrcp, '.csv', sep=''), row.names=1))

ntk <- wdpaturnbyMPAbymod$no_take %in% c('All', 'Part') # no take reserves (index into wdpaturnbyMPAbymod)
sum(ntk[wdpaturnbyMPAbymod$rcp==rcp & wdpaturnbyMPAbymod$model==1]) # 50 (43 in CA)


# Fraction of species lost
	# all MPAs
	wdpaturnbyMPAbymod[,.(mean=mean(flost,na.rm=TRUE)),by=c('rcp','model')][,.(mean=mean(mean), se=se(mean), min=min(mean), max=max(mean)),by='rcp'] # first mean is across MPAs within climate models, then average across climate models

	wdpaturnbyMPAbymod[,.(mean=mean(flost,na.rm=TRUE)),by=c('rcp','model','region.one')][,.(mean=mean(mean), se=se(mean), min=min(mean), max=max(mean)),by=c('rcp','region.one')] # first mean is across MPAs within climate models, then average within regions

	mod<- wdpaturnbyMPAbymod[,.(mean=mean(flost,na.rm=TRUE)),by=c('model','region.one')][,lm(mean~region.one-1)] # first mean is across MPAs within climate models within regions, then do an anova
		summary(mod)
		length(mod$fitted.values) # number of data points

	mod<- wdpaturnbyMPAbymod[,.(mean=mean(flost,na.rm=TRUE)),by=c('model','region.one')][,kruskal.test(x=mean, g=region.one)] # first mean is across MPAs within climate models within regions, then do a non-parametric Kruskal-Wallis Rank Sum Test
		mod

#	wdpaturnbyMPAbymod[,sd(flost,na.rm=TRUE),by=c('rcp','model')]
#	wdpaturnbyMPAbymod[,min(flost,na.rm=TRUE),by=c('rcp','model')]
#	wdpaturnbyMPAbymod[,max(flost,na.rm=TRUE),by=c('rcp','model')]

	# no-take
#	wdpaturnbyMPAbymod[ntk,mean(flost,na.rm=TRUE),by=c('rcp','model')]
#	wdpaturnbyMPAbymod[ntk,sd(flost,na.rm=TRUE),by=c('rcp','model')]
#	wdpaturnbyMPAbymod[ntk,min(flost,na.rm=TRUE),by=c('rcp','model')]
#	wdpaturnbyMPAbymod[ntk,max(flost,na.rm=TRUE),by=c('rcp','model')]

# Fraction of species gained (fraction of final community)
	# all MPAs
	wdpaturnbyMPAbymod[,.(mean=mean(fgainedalt,na.rm=TRUE)),by=c('rcp','model')][,.(mean=mean(mean), se=se(mean), min=min(mean), max=max(mean)),by='rcp'] # first mean is across MPAs within climate models

	wdpaturnbyMPAbymod[,.(mean=mean(fgainedalt,na.rm=TRUE)),by=c('rcp','model','region.one')][,.(mean=mean(mean), se=se(mean), min=min(mean), max=max(mean)),by=c('rcp','region.one')] # first mean is across MPAs within climate models, then average within regions

	mod<- wdpaturnbyMPAbymod[,.(mean=mean(fgainedalt,na.rm=TRUE)),by=c('model','region.one')][,lm(mean~region.one-1)] # first mean is across MPAs within climate models within regions, then do an anova
		summary(mod)
		length(mod$fitted.values) # number of data points

	mod<- wdpaturnbyMPAbymod[,.(mean=mean(fgainedalt,na.rm=TRUE)),by=c('model','region.one')][,kruskal.test(x=mean, g=region.one)] # first mean is across MPAs within climate models within regions, then do a non-parametric Kruskal-Wallis Rank Sum Test
		mod

#	b <- wdpaturnbyMPAbymod[,sd(fgainedalt,na.rm=TRUE),by=c('rcp','model')]; b # sd across MPAs within climate models (rcp45 or rcp85)
#		b[,mean(V1),by='rcp'] 


	# no-take
#	wdpaturnbyMPAbymod[ntk,mean(fgained,na.rm=TRUE),by=c('rcp','model')]
#	wdpaturnbyMPAbymod[ntk,sd(fgained,na.rm=TRUE),by=c('rcp','model')]

# Similarity
	# all MPAs
	wdpaturnbyMPAbymod[,.(mean=mean(beta_sor,na.rm=TRUE)),by=c('rcp','model')][,.(mean=mean(mean), se=se(mean), min=min(mean), max=max(mean)),by='rcp'] # first mean is across MPAs within climate models

	# examine similarity by MPA size
		wdpaturnbyMPAbymod[,cor.test(rep_area, beta_sor), by=c('rcp', 'model')]
		
			# average across rcps and climate models
		a <- wdpaturnbyMPAbymod[rep_area != 0,.(beta_sor=mean(beta_sor), abeta_sor=asin(sqrt(mean(beta_sor))), larea=log(mean(rep_area))),by=wdpapolyID]
		mod <- a[,lm(abeta_sor ~ larea)]; summary(mod)
		a[,plot(larea, beta_sor)] # plot all points
		i <- order(a$larea); a[, lines(larea[i], sin(fitted(mod)[i])^2)] # add fit line
		a[,range(exp(larea))]
		range(sin(fitted(mod))^2)
		
			# mixed-effects model
#		datsc <- as.data.frame(wdpaturnbyMPAbymod[rep_area != 0,.(rep_area, beta_sor, rcp, model)]) # remove zero area
#		datsc$lrep_area <- log(datsc$rep_area)
#		datsc$rcp <- as.factor(datsc$rcp)
#		datsc$model <- as.factor(datsc$model)
##		datsc$lrep_areasc <- scale(datsc$lrep_area)
#		datsc$y <- asin(sqrt(datsc$beta_sor))
#		mod <- lmer(y ~ lrep_area + (1|rcp/model), data=datsc, REML=FALSE)
#		mod2 <- lmer(y ~ 1 + (1|rcp/model), data=datsc, REML=FALSE)
#		summary(mod)
#		Anova(mod)
#		anova(mod, mod2)

#		plot(datsc$rep_area, datsc$beta_sor, log='x')
#		fit <- data.frame(rep_area=datsc$rep_area, lrep_area=datsc$lrep_area)
#		fit$y <- fixef(mod)[1] + fit$lrep_area * fixef(mod)[2]
#		fit$beta_sor <- sin(fit$y)^2
#		fit <- fit[order(fit$lrep_area),]
#		lines(fit$rep_area, fit$beta_sor, col='red')
#
#		head(fit)
#		tail(fit)

	# examine similarity by MPA latitudinal range
		datsc <- as.data.frame(wdpaturnbyMPAbymod[,.(lat_min, lat_max, beta_sor, rcp, model)])
		datsc$latrng <- datsc$lat_max - datsc$lat_min
		datsc$rcp <- as.factor(datsc$rcp)
		datsc$model <- as.factor(datsc$model)
#		datsc$latrngsc <- scale(datsc$latrng)
		datsc$y <- asin(sqrt(datsc$beta_sor))
		mod <- lmer(y ~ latrng + (1|rcp/model), data=datsc, REML=FALSE)
		mod2 <- lmer(y ~ 1 + (1|rcp/model), data=datsc, REML=FALSE)
		summary(mod)
		Anova(mod)
		anova(mod, mod2)

		plot(datsc$latrng, datsc$beta_sor)
		fit <- data.frame(latrng=datsc$latrng)
		fit$y <- fixef(mod)[1] + fit$latrng * fixef(mod)[2]
		fit$beta_sor <- sin(fit$y)^2
		fit <- fit[order(fit$latrng),]
		lines(fit$latrng, fit$beta_sor, col='red')

		head(fit)
		tail(fit)

	# examine similarity by MPA size and lat rng
		datsc <- as.data.frame(wdpaturnbyMPAbymod[rep_area != 0,.(rep_area, lat_min, lat_max, beta_sor, rcp, model)]) # remove zero area
		datsc$lrep_area <- log(datsc$rep_area)
		datsc$latrng <- datsc$lat_max - datsc$lat_min
		datsc$rcp <- as.factor(datsc$rcp)
		datsc$model <- as.factor(datsc$model)
		pvars <- c('lrep_area', 'latrng')
		scvars <- c('lrep_areasc', 'latrngsc')
		datsc[scvars] <- lapply(datsc[pvars], scale)
		datsc$y <- asin(sqrt(datsc$beta_sor))
		mod <- lmer(y ~ lrep_areasc*latrngsc + (1|rcp/model), data=datsc, REML=FALSE)
		mod2 <- lmer(y ~ lrep_areasc + latrngsc + (1|rcp/model), data=datsc, REML=FALSE)
		mod3.1 <- lmer(y ~ latrngsc + (1|rcp/model), data=datsc, REML=FALSE)
		mod3.2 <- lmer(y ~ lrep_areasc + (1|rcp/model), data=datsc, REML=FALSE)
		mod4 <- lmer(y ~ 1 + (1|rcp/model), data=datsc, REML=FALSE)
		summary(mod)
		Anova(mod)
		anova(mod, mod2)
		anova(mod3.1, mod3.2) 

		plot(mod, sin(fitted(.))^2 ~ lrep_areasc)
		plot(mod, sin(fitted(.))^2 ~ latrngsc)
	
		require(interplot) # to plot the interaction
		interplot(m=mod, var1='latrngsc', var2='lrep_areasc') +
			xlab('log(area) scaled') + 
			ylab('Coefficient for scaled latitudinal range')

		interplot(m=mod, var2='latrngsc', var1='lrep_areasc') +
			ylab('Coefficient for log(area) scaled') + 
			xlab('Latitudinal range scaled')

	# examine similarity by region change in temperature
		# read in change in temp by region
		trend1 <- read.csv('data/climTrendbyreg_rcp45.csv', row.names=1)
			trend1$rcp <- 45
		trend2 <- read.csv('data/climTrendbyreg_rcp85.csv', row.names=1)
			trend2$rcp <- 85
		nms <- c('region', 'rcp', 'delta_surf', 'delta_bott')
		trend <- as.data.table(rbind(trend1[,nms], trend2[,nms]))

		# calculate mean turnover by region
		turnbyreg <- wdpaturnbyMPAbymod[,.(beta_sor=mean(beta_sor)),by=.(rcp,region.one)]
		
		# merge together
		setkey(trend, region, rcp)
		setkey(turnbyreg, region.one, rcp)
		turnbyreg <- trend[turnbyreg] # merge
		
		# examine
		turnbyreg[rcp==85,]
		turnbyreg[rcp==45,]
		
		# plot
		turnbyreg[rcp==85,plot(delta_bott, beta_sor)]
		turnbyreg[rcp==45,plot(delta_bott, beta_sor)]

		# linear model
		turnbyreg[, cor.test(beta_sor, delta_bott), by=rcp]		
		turnbyreg[, cor.test(beta_sor, delta_surf), by=rcp]		

	# examine similarity by region change in temperature AND by lat rng
		# average across climate models and rcps
		turnbyMPA <- wdpaturnbyMPAbymod[,.(beta_sor=mean(beta_sor), latrng=mean(lat_max) - mean(lat_min), region=unique(region.one)),by=.(wdpapolyID)]
		
		# average temp trends across rcps
		trendave <- trend[,.(delta_surf=mean(delta_surf), delta_bott=mean(delta_bott)), by=region]
		
		# merge in regional temperature information
		setkey(trendave, region)
		setkey(turnbyMPA, region)
		turnbyMPA <- trendave[turnbyMPA] # merge
		
		# linear model
		mod <- turnbyMPA[, lm(asin(sqrt(beta_sor)) ~ latrng*delta_bott)]
		summary(mod)

		mod <- turnbyMPA[, lm(asin(sqrt(beta_sor)) ~ latrng*delta_surf)]
		mod2 <- turnbyMPA[, lm(asin(sqrt(beta_sor)) ~ latrng+delta_surf)] # used this in the paper
		summary(mod)
		summary(mod2)
		anova(mod, mod2)
		AIC(mod, mod2)

# Plot of fraction species lost as a density across MPAs within each climate model
	inds <- expand.grid(rcp=unique(wdpaturnbyMPAbymod$rcp), model=unique(wdpaturnbyMPAbymod$model))
	ds <- vector('list', nrow(inds))
	for(i in 1:nrow(inds)){
		ds[[i]] <- density(wdpaturnbyMPAbymod$flost[wdpaturnbyMPAbymod$rcp==inds$rcp[i] & wdpaturnbyMPAbymod$model==inds$model[i]], na.rm=TRUE, cut=1, kernel='gaus', adjust=2)
	}
	
	cols <- c('blue', 'green') # blue is rcp, green is otherrcp
	plot(ds[[1]], type='n', ylim=c(0,7), xlab='Fraction of species lost', main='')
	for(i in 1:nrow(inds)){	
		if(inds$rcp[i]==rcp) lines(ds[[i]], col=cols[1])
		if(inds$rcp[i]==otherrcp) lines(ds[[i]], col=cols[2])
	}


##############################################
# Calculate turnover within MPA networks
# By ensemble mean distributions
##############################################
# pick a MPA network
#length(unique(wdpa$wdpaid[wdpa$sub_loc=='US-CA' & wdpa$mang_auth=="California Department of Fish and Game"])) # 45 MPAs
#grids <- wdpa[wdpa$sub_loc=='US-CA' & wdpa$mang_auth=="California Department of Fish and Game", c('lat', 'lon')] #MLPA
#	dim(grids)
#
## select initial and final timeperiod for these grids
#thesespp <- presmap[paste(presmap$lat, presmap$lon) %in% paste(grids$lat, grids$lon) & presmap$period %in% c('2006-2020', '2081-2100'),]
#thesespp <- thesespp[thesespp$pres,] # trim to present spp
#	dim(thesespp) # 5340
#	sort(unique(thesespp$region))
##	thesespp <- thesespp[thesespp$region == 'AFSC_WCTri',] # trim to one region
#	thesespp <- thesespp[thesespp$region == 'NWFSC_WCAnn',] # trim to one region
#	dim(thesespp) # 3009
#	
## trim to unique species in each period
#thesespp <- thesespp[!duplicated(thesespp[,c('sppocean', 'period')]),]
#	dim(thesespp) # 203
#
## evaluate turnover
#turnover(thesespp[,c('region', 'lat', 'lon', 'period', 'sppocean', 'pres')])
#

##############################################
# Calculate turnover within select MPA networks
# Across all models in the ensemble
##############################################

# merge in MPAnetwork ID (using data.tables)
wdpanets <- wdpa[!is.na(network),] # trim out non-network rows
setkey(wdpanets, lat, lon)
setkey(thesesppbymod, lat, lon)
thesesppbymodnet <- wdpanets[thesesppbymod,.(sppocean, period, rcp, model, pres, network, lat, lon), allow.cartesian=TRUE, nomatch=0] # datatable join, using the keys set in each dt. the order specifies to keep all rows of thesesppbymod. allow.cartesian is needed because the result is larger than either parent. drop locations not in a network.
	dim(thesesppbymodnet) # 2,637,693? 5,816,336, now 8995646 (4/24/2016)
	
# summarize by unique species in each period in each MPA in each RCP in each model (pres or not)
thesesppbynetbymod <- thesesppbymodnet[,max(pres),by="sppocean,period,network,rcp,model"] # DT aggregate function
	thesesppbynetbymod[,pres:=as.logical(V1)]
	thesesppbynetbymod[,V1:=NULL] # drop V1
	dim(thesesppbynetbymod) # 30,077

# calculate turnover by MPA, rcp, and model. very quick.
thesesppbynetbymod[,nstart:=nstart(period, sppocean, pres, periods), by="network,rcp,model"]
thesesppbynetbymod[,nend:=nend(period, sppocean, pres, periods), by="network,rcp,model"]
thesesppbynetbymod[,nlost:=nlost(period, sppocean, pres, periods), by="network,rcp,model"]
thesesppbynetbymod[,ngained:=ngained(period, sppocean, pres, periods), by="network,rcp,model"]
thesesppbynetbymod[,flost:=flost(period, sppocean, pres, periods), by="network,rcp,model"]
thesesppbynetbymod[,fgained:=fgained(period, sppocean, pres, periods), by="network,rcp,model"]
thesesppbynetbymod[,fgainedalt:=fgainedalt(period, sppocean, pres, periods), by="network,rcp,model"]
thesesppbynetbymod[,beta_sor:=beta_sor(period, sppocean, pres, periods), by="network,rcp,model"]

# trim out duplicate rows within an mpa/rcp/model
setkey(thesesppbynetbymod, network,rcp,model)
wdpaturnbynetbymod <- unique(thesesppbynetbymod)
	dim(wdpaturnbynetbymod) # 78 (3 networks x 13 models x 2 scenarios)
	wdpaturnbynetbymod[,c('sppocean','period','pres'):=NULL]
	
# merge in mpa metadata
setkey(wdpanets, network)
uwdpanets <- unique(wdpanets) # only the lines with unique mpas
	dim(uwdpanets) # 3
setkey(uwdpanets, network)
setkey(wdpaturnbynetbymod, network,rcp,model)
wdpaturnbynetbymod <- uwdpanets[wdpaturnbynetbymod, .(network,rcp,model,nstart,nend,nlost,ngained,flost,fgained,fgainedalt,beta_sor)]
	dim(wdpaturnbynetbymod) # 78


# write out wdpaturnbyMPA
write.csv(wdpaturnbynetbymod, paste('data/wdpaturnbynetbymod_', runtype, projtype, '_rcp', rcp, '&', otherrcp, '.csv', sep=''))


###############################
# Examine change within MPA networks
# Across each model in the ensemble
##############################
wdpaturnbynetbymod <- as.data.table(read.csv(paste('data/wdpaturnbynetbymod_', runtype, projtype, '_rcp', rcp, '&', otherrcp, '.csv', sep=''), row.names=1)) # network results
wdpaturnbyMPAbymod <- as.data.table(read.csv(paste('data/wdpaturnbyMPAbymod_', runtype, projtype, '_rcp', rcp, '&', otherrcp, '.csv', sep=''), row.names=1)) # individual MPA results



# Fraction of species lost
	# all MPA networks
	a <- wdpaturnbynetbymod[,mean(flost,na.rm=TRUE),by=c('rcp','model')]; a
		a[,range(V1),by='rcp'] # 4%-19% or 11%-27% min
		a[,mean(V1),by='rcp'] # 10% or 21% mean flost (rcp45 or rcp85)
		a[,median(V1),by='rcp'] # 
		a[,se(V1),by='rcp'] # 1.3% or 1.2% se across models (rcp45 or rcp85)
	wdpaturnbyMPAbymod[,sd(flost,na.rm=TRUE),by=c('rcp','model')]
	wdpaturnbyMPAbymod[,min(flost,na.rm=TRUE),by=c('rcp','model')]
	wdpaturnbyMPAbymod[,max(flost,na.rm=TRUE),by=c('rcp','model')]

	# within the individual MPAs of these networks
	a <- wdpaturnbyMPAbymod[!is.na(network),mean(flost,na.rm=TRUE),by=c('rcp','model')]; a
		a[,range(V1),by='rcp'] # 
		a[,mean(V1),by='rcp'] # 
		a[,median(V1),by='rcp'] # 
		a[,se(V1),by='rcp'] # 
	

# Fraction of species gained (fraction of final community)
	# all MPA networks
	a <- wdpaturnbynetbymod[,mean(fgainedalt,na.rm=TRUE),by=c('rcp','model')]; a
		a[,mean(V1),by='rcp'] # 
		a[,se(V1),by='rcp'] # 
		a[,range(V1),by='rcp'] # 
		a[,median(V1),by='rcp'] # 
		a[,sd(V1),by='rcp'] # 
	b <- wdpaturnbyMPAbymod[,sd(fgainedalt,na.rm=TRUE),by=c('rcp','model')]; b
		b[,mean(V1),by='rcp'] # 

	# within the individual MPAs of these networks
	a <- wdpaturnbyMPAbymod[!is.na(network),mean(fgainedalt,na.rm=TRUE),by=c('rcp','model')]; a
		a[,mean(V1),by='rcp'] # 
		a[,se(V1),by='rcp'] # 
		a[,range(V1),by='rcp'] # 
		a[,median(V1),by='rcp'] # 

# Similarity (Sorenson)
	# all MPA networks
	a <- wdpaturnbynetbymod[,mean(beta_sor,na.rm=TRUE),by=c('rcp','model')]; a
		a[,mean(V1),by='rcp'] # 
		a[,se(V1),by='rcp'] # 
		a[,range(V1),by='rcp'] # 
		a[,median(V1),by='rcp'] # 
		a[,sd(V1),by='rcp'] # 
	b <- wdpaturnbyMPAbymod[,sd(fgainedalt,na.rm=TRUE),by=c('rcp','model')]; b
		b[,mean(V1),by='rcp'] # 

	# within the individual MPAs of these networks
	a <- wdpaturnbyMPAbymod[!is.na(network),mean(beta_sor,na.rm=TRUE),by=c('rcp','model')]; a
		a[,mean(V1),by='rcp'] # 
		a[,se(V1),by='rcp'] # 
		a[,range(V1),by='rcp'] # 
		a[,median(V1),by='rcp'] # 
