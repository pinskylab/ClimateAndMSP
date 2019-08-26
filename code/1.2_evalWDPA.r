# evaluate protected areas against shifts in species distribution
# calculate species gains, losses, turnover, etc.


############
## Flags
############

## choose which runs to use

# choose the rcp
rcp <- 85
otherrcp <- 45

# select initial and final timeperiod for these grids
periods <- c('2006-2020', '2081-2100')

####################
## helper functions
####################
#require(Hmisc)
require(data.table)
#require(lme4) # for mixed-effects models
#require(car) # for testing ME models


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

wdpa <- read.csv(gzfile('output/wdpa_cov_by_grid0.05.csv.gz')) # shows which MPAs are in which grid cells. each line is a unique grid cell-MPA combination.
	wdpa$lon[wdpa$lon<0] <- wdpa$lon[wdpa$lon<0] + 360 # convert lon to positive

	wdpa <- as.data.table(wdpa)

	# set up MPA networks
	wdpa[SUB_LOC == 'US-CA' & MANG_AUTH == 'State Fish and Wildlife',network:='mlpa']
	wdpa[(SUB_LOC %in% c('US-DE', 'US-FL', 'US-GA', 'US-MA', 'US-MD', 'US-ME', 'US-NC', 'US-NJ', 'US-NY', 'US-RI', 'US-SC', 'US-VA')) & (lon > 280), network:='eastcoast'] # by choosing by state, this won't include federal water closures
	wdpa[(SUB_LOC %in% c('US-AK')), network:='ak'] # by choosing by state, this won't include federal water closures
	
	# calculate mpa extent
	wdpa[, ':=' (lat_min=min(lat), lat_max=max(lat)), by=WDPA_PID]

#	wdpa[network=='mlpa',plot(lon,lat)]
#	wdpa[network=='eastcoast',plot(lon,lat)]
#	wdpa[network=='ak',plot(lon,lat)]
#	require(maps); map(database='world2',add=TRUE)


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

