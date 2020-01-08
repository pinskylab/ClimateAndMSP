# evaluate protected areas against shifts in species distribution
# calculate species gains, losses, turnover, etc.


############
## Flags
############

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
require(dplyr)
require(lme4) # for mixed-effects models
require(car) # for testing ME models
require(lattice)
require(RColorBrewer)


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

# pick the first non-NA item in a list, or else return NA
notNA <- function(x){
	x <- x[!is.na(x)]
	if(length(x)>0){ return(as.character(x[1]))}
	if(length(x)==0){ return('NA')}
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


# load presmap for all model runs, a bit slow (large files)
load(paste(presmapbymodfolder, 'presmapbymod_', runtype, projtype, '_rcp', rcp, '.RData', sep='')) # loads presmap data.frame with presence/absence information from each model (slow to load).
	presmapbymod.1 <- as.data.table(presmapbymod) # seems to remove presmapbymod as well?

	setkey(presmapbymod.1,period,pres)
	presmapbymod.1 <- presmapbymod.1[.(periods,TRUE),nomatch=0] # trim to periods and pres=TRUE rows
		presmapbymod.1[,sum(is.na(pres))] # should be zero
		dim(presmapbymod.1)
	presmapbymod.1[,rcp:=rcp] # add label for rcp

load(paste(presmapbymodfolder, 'presmapbymod_', runtype, projtype, '_rcp', otherrcp, '.RData', sep='')) # loads presmap data.frame with presence/absence information from each model (slow to load) (for the other rcp, the one not used for planning)
	presmapbymod <- as.data.table(presmapbymod)

	setkey(presmapbymod,period,pres)
	presmapbymod.3 <- presmapbymod[.(periods,TRUE),nomatch=0] # trim to rows (time periods) relevant to analysi
		presmapbymod.3[,sum(is.na(pres))] # should be zero
		dim(presmapbymod.3)
	presmapbymod.3[,rcp:=otherrcp]

# bind the presmaps together
	presmapbymod <- rbind(presmapbymod.1, presmapbymod.3)
		nrow(presmapbymod.1)
		nrow(presmapbymod.3)
		dim(presmapbymod)

	rm(presmapbymod.1)
	rm(presmapbymod.3)

# trim presmap to species present, locations and time periods relevant to wdpa analysis
	regs <- sort(unique(presmapbymod$region))
	regstokeep <- setdiff(regs, c('NEFSC_NEUSFall', 'AFSC_WCTri', 'DFO_NewfoundlandSpring')) # remove some surveys to avoid duplicates in a region
	setkey(presmapbymod, region)
	thesesppbymod <- presmapbymod[.(regstokeep)] # trim out duplicate regions
	dim(thesesppbymod) # 19,859,747

	rm(presmapbymod)
	
# Add a grid ID
thesesppbymod[,gridID:=paste(lat,lon)]




##################################################
# Generate random MPA networks
# Across gradients of area and temperature extent
##################################################
climGrid <- read.csv('data/climGrid.csv', row.names=1) # for temperature climatology and depth

temps <- seq(0.01,1,length.out=10) # how many steps of temperature extent (0 to 100%)
sizes <- seq(0.01,0.5,length.out=10) # which steps of area (could be 0 to 100%). Note that a network constrained to <100% of lat extent can't include all planning sites.
nreps <- 10 # how many repeats at each combination of lat and size

regs <- as.character(thesesppbymod[,sort(unique(region))])

# to store results
randMPAs <- expand.grid(region=regs, tempstep=temps, sizestep=sizes, repnum=1:nreps, temprng=NA, size=NA, meanturnIndiv=NA, netwrkturn=NA)

# a slow loop, especially on higher values of j and k
for(i in 1:length(regs)){
	print(regs[i])

	# set up the set of planning units
	setkey(thesesppbymod, region)
	punits <- thesesppbymod[regs[i],.(lat=lat,lon=lon,gridID=gridID)]
	punits <- as.data.frame(punits[!duplicated(punits),])
	oldrow <- nrow(punits)
	punits <- merge(punits, climGrid[climGrid$region==regs[i], c('lat', 'lon', 'bottemp.clim.int')], all.x=TRUE)
	if(nrow(punits) != oldrow) print(paste('i=', i, ': rows changed during merge with clim'))
	names(punits)[names(punits)=='bottemp.clim.int'] <- 'bt'

	tx <- diff(range(punits$bt, na.rm=TRUE))
	nu <- length(unique(punits$gridID))

	setkey(thesesppbymod,gridID,region)

	for(j in 1:length(temps)){
		cat(paste('\nj=',j, sep=''))
		for(k in 1:length(sizes)){
			cat(paste('\tk=',k,': ', sep=''))
			for(r in 1:nreps){
				cat(paste('r=',r,' '))
	
				# select a first MPA
				punits$sel <- FALSE
				firstind <- sample(1:nrow(punits),1)
				punits$sel[firstind] <- TRUE
				ft <- punits$bt[firstind] # pick a first MPA

				# calculate nearby planning units from first MPA (based on temp extent)
				punits$nearby <- FALSE # initialize "nearby" flag (within a certain fraction of lat extent)
				punits$nearby[abs(punits$bt-ft)/tx <= temps[j]] <- TRUE # mark some as nearby
				
				fracwholearea <- 1/nu
				fracnearbyarea <- 1/sum(punits$nearby)

				# add more MPAs until size criterion met or all nearby units selected
				while(fracwholearea <= sizes[k] & fracnearbyarea <= 1 & sum(punits$nearby & !punits$sel)>0 ) {
					newunit <- sample(which(punits$nearby & !punits$sel),1) # select a new unit
					punits$sel[newunit] <- TRUE # assign as selected
					fracwholearea <- sum(punits$sel)/nu
					fracnearbyarea <- sum(punits$sel)/sum(punits$nearby)
					punits$nearby <- ((punits$bt - min(punits$bt[punits$sel]))/tx <= temps[j] & punits$bt >= min(punits$bt[punits$sel])) | ((max(punits$bt[punits$sel]) - punits$bt)/tx < temps[j] & punits$bt <= max(punits$bt[punits$sel])) # measure from the min temp if new potential temps are > min, and measure from max temp if new potential temps are < max
				}
				
				# basic calcs for this random network
				ind <- with(randMPAs, which(region==regs[i] & tempstep==temps[j] & sizestep==sizes[k] & repnum==r))
				randMPAs$temprng[ind] <- with(punits[punits$sel,], max(bt)-min(bt))/tx
				randMPAs$size[ind] <- sum(punits$sel)/nu
				
				# evaluate ecological turnover for each MPA
				randMPAs$meanturnIndiv[ind] <- thesesppbymod[.(punits$gridID[punits$sel], regs[i]),.(beta_sor=beta_sor(period, sppocean, pres, periods)), by="gridID,rcp,model"][,mean(beta_sor)]

				# collapse to a network
				thesesppbymodnet <- thesesppbymod[.(punits$gridID[punits$sel], regs[i]),.(pres=max(pres)), by="sppocean,period,rcp,model"] #  aggregate function within the selected grids
				randMPAs$netwrkturn[ind] <- thesesppbymodnet[,.(beta_sor=beta_sor(period, sppocean, pres, periods)), by="rcp,model"][,mean(beta_sor)]

			}			
		}
	}
	cat('\n')
}

# write out
write.csv(randMPAs, file=paste('cmsp_data/randMPAs_byBT_', Sys.Date(), '.csv', sep=''))


