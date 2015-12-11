# evaluate protected areas

require(Hmisc)

## Set working directories
if(Sys.info()["nodename"] == "pinsky-macbookair"){
	setwd('~/Documents/Rutgers/Range projections/proj_ranges/')
	projfolder <- '../CEmodels_proj' # holds model projections (outside Git)
	modfolder <- '../CEModels' # holds the models (outside Git)
	}
if(Sys.info()["nodename"] == "amphiprion.deenr.rutgers.edu"){
	setwd('~/Documents/range_projections/')
	projfolder <- 'CEmodels_proj'
	modfolder <- 'CEmodels'
	}
# could add code for Lauren's working directory here

####################
## helper functions
####################
lu <- function(x) return(length(unique(x)))

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
# for use with Hmisc::summarize() but doesn't work for some reason
# expects region in column 1, lat (2), lon (3), time period in column 4, spp in col 5, pres in col 6
# only present species included
turnover <- function(x){
	pers <- sort(unique(x[,4]))
	if(length(pers) != 2) stop(paste('expecting 2 time periods:', unique(x[,1]), unique(x[,2]), unique(x[,3]), paste(unique(x[,4], collapse=','))))
	initcomm <- x[x[,4] == pers[1] & x[,6]==TRUE,5]
	finalcomm <- x[x[,4] == pers[2] & x[,6]==TRUE,5]
	
	a <- length(intersect(initcomm, finalcomm)) # number of shared species
	lost <- length(setdiff(initcomm, finalcomm)) # number of species lost
	gained <- length(setdiff(finalcomm, initcomm))
	return(c(nstart=length(initcomm), nend=length(finalcomm), nlost=lost, ngained=gained, flost=lost/length(initcomm), fgained=gained/length(initcomm), beta_sor=2*a/(2*a+gained+lost)))
}

###########################################################
## Calculate species richness change within PAs
###########################################################

load('data/rich.RData') # loads rich data.frame: richness by grid cell by time period (ensemble mean)
wdpa <- read.csv('data/wdpa_cov_by_grid0.25.csv', row.names=1)

# any overlapping richness projections? different regions, same grid cell
# average across regions
sum(duplicated(rich[,c('period', 'lat', 'lon')])) # yes: 4110
numdup <- aggregate(list(n = rich$region), by=list(period=rich$period, lat=rich$lat, lon=rich$lon), FUN=length)
	sum(numdup$n>1) # 4110, which means the aggregate worked
	sum(numdup$n>2) # only doubles
numdupreg <- aggregate(list(nreg = rich$region), by=list(period=rich$period, lat=rich$lat, lon=rich$lon), FUN=lu)
	numdupreg <- merge(numdupreg, numdup)
	sum(numdupreg$n>1 & numdupreg$nreg==1) # 0: no duplications within a region
	sum(numdupreg$n>1 & numdupreg$nreg==2) # 4110: all duplications involve 2 regions
	
richave <- aggregate(list(rich = rich$rich), by=list(period=rich$period, lat=rich$lat, lon=rich$lon), FUN=mean) # take the average across regions
	dim(richave)

# convert wdpa lon to + to match grid
wdpa$lon[wdpa$lon < 0] = wdpa$lon[wdpa$lon < 0] + 360

# merge richness with WDPA
	dim(wdpa)
	length(unique(wdpa$wdpapolyID)) # 625 PAs

wdparich <- merge(wdpa, richave, all.x=TRUE)
	dim(wdparich)
	length(unique(wdparich$wdpapolyID)) # 625 PAs

	# examine the merge
	sum(is.na(wdparich$rich)) # 169 values missing
	wdparich[is.na(wdparich$rich), c('lat', 'lon', 'country', 'sub_loc', 'period', 'rich')] # seem to be outside the richness projection. not sure why they were in the intersection with the climate grid
	
	# remove missing values
	wdparich <- wdparich[!is.na(wdparich$rich),]
	
# average richness within each PA and period
# important when PAs cross multiple grid cells
wdparichave <- Hmisc::summarize(wdparich[,c('rich', 'prop_grid')], by=llist(period=wdparich$period, wdpapolyID=wdparich$wdpapolyID), FUN=wmean, stat.name='rich')
	dim(wdparichave)
	length(unique(wdparichave$wdpapolyID)) # 611 PAs (lost a few)
wdparichave <- merge(wdparichave, wdpa[!duplicated(wdpa$wdpapolyID),]) # merge in other data
	dim(wdparichave) # 3045
	
# calculate trend in richness by PA
pds <- as.numeric(unlist(strsplit(as.character(wdparichave$period), split='-')))
dim(pds) <- c(2,nrow(wdparichave))
wdparichave$periodmids <- colMeans(pds) # midpoints of each time period

wdparichtrend <- Hmisc::summarize(X=wdparichave[,c('periodmids', 'rich')], by=llist(wdpapolyID=wdparichave$wdpapolyID, dummy=wdparichave$wdpapolyID), FUN=calcrichtrendmids, stat.name='richtrend') # calculate trend in richness (# taxa) by PA. last column is filled with NAs, so add a dummy column
	wdparichtrend <- wdparichtrend[,-grep('dummy', names(wdparichtrend))]
	dim(wdparichtrend)
wdparichtrend <- merge(wdparichtrend, wdpa[!duplicated(wdpa$wdpapolyID),]) # merge in other data
	dim(wdparichtrend)
	head(wdparichtrend)
wdparichtrend <- merge(wdparichtrend, wdparichave[wdparichave$period=='2006-2020',c('wdpapolyID', 'rich')]) # merge in original richness
	dim(wdparichtrend)
	names(wdparichtrend)[names(wdparichtrend)=='rich'] = 'rich2006_2020'
wdparichtrend <- merge(wdparichtrend, wdparichave[wdparichave$period=='2081-2100',c('wdpapolyID', 'rich')], all.x=TRUE) # merge in final richness
	dim(wdparichtrend)
	names(wdparichtrend)[names(wdparichtrend)=='rich'] = 'rich2081_2100'
	head(wdparichtrend)
	which(is.na(wdparichtrend$rich2081-2100)) # six missing values

# calculate fractional gain/loss of species in PAs
wdparichtrend$richchange = (wdparichtrend$rich2081_2100 - wdparichtrend$rich2006_2020)/wdparichtrend$rich2006_2020
wdparichtrend$richlogRR = log10(wdparichtrend$rich2081_2100/wdparichtrend$rich2006_2020)
	
# calculate trend in richness for all grid cells
pds <- as.numeric(unlist(strsplit(as.character(rich$period), split='-')))
dim(pds) <- c(2,nrow(rich))
rich$periodmids <- colMeans(pds) # midpoints of each time period
richtrend <- Hmisc::summarize(X=rich[,c('periodmids', 'rich')], by=llist(region=rich$region, lat=rich$lat, lon=rich$lon, dummy=rich$lon), FUN=calcrichtrendmids, stat.name='trend') # calculate richness (# taxa) by grid cell in each time period. for an unknown resion, the last column in the by list gets NAs, so padded with a dummy column here
	richtrend <- richtrend[,-grep('dummy', names(richtrend))]
	head(richtrend)

richtrend <- merge(richtrend, rich[rich$period=='2006-2020',c('region', 'lat', 'lon', 'rich')]) # merge in original richness
	dim(richtrend)
	names(richtrend)[names(richtrend)=='rich'] = 'rich2006_2020'
richtrend <- merge(richtrend, rich[rich$period=='2081-2100',c('region', 'lat', 'lon', 'rich')], all.x=TRUE) # merge in final richness
	dim(richtrend)
	names(richtrend)[names(richtrend)=='rich'] = 'rich2081_2100'
	head(richtrend)
	which(is.na(richtrend$rich2081-2100)) # three missing values

# calculate fractional gain/loss of species in grid cells
richtrend$richchange = (richtrend$rich2081_2100 - richtrend$rich2006_2020)/richtrend$rich2006_2020
richtrend$richlogRR = log10(richtrend$rich2081_2100/richtrend$rich2006_2020)
	
	
###########################################################
## Calculate community turnover within PAs
###########################################################

load('data/turn.RData') # loads turnover data.frame: sorenson similiarty by grid cell by time period (ensemble mean)
wdpa <- read.csv('data/wdpa_cov_by_grid0.25.csv', row.names=1)

# any overlapping richness projections? different regions, same grid cell
# average across regions
sum(duplicated(turn[,c('lat', 'lon')])) # yes: 1374
numdup <- aggregate(list(n = turn$region), by=list(lat=turn$lat, lon=turn$lon), FUN=length)
	sum(numdup$n>1) # 1265, which means the aggregate worked
	sum(numdup$n>2) # 109 triples or more... why?
numdupreg <- aggregate(list(nreg = turn$region), by=list(lat=turn$lat, lon=turn$lon), FUN=lu)
	numdupreg <- merge(numdupreg, numdup)
	sum(numdupreg$n>1 & numdupreg$nreg==1) # 0: no duplications within a region
	sum(numdupreg$n>1 & numdupreg$nreg>=2) # 1265: all duplications involve 2+ regions
	sum(numdupreg$n>1 & numdupreg$nreg>=3) # 109: some duplications involve 3 regions... why?
	
turnave <- aggregate(list(beta_sor = turn$beta_sor, fgained=turn$fgained, flost=turn$flost), by=list(lat=turn$lat, lon=turn$lon), FUN=mean) # take the average across regions
	dim(turnave)

# convert wdpa lon to + to match grid
wdpa$lon[wdpa$lon < 0] = wdpa$lon[wdpa$lon < 0] + 360

# merge turnover with WDPA
	dim(wdpa)
	length(unique(wdpa$wdpapolyID)) # 625 PAs

wdpaturn <- merge(wdpa, turnave, all.x=TRUE)
	dim(wdpaturn)
	length(unique(wdpaturn$wdpapolyID)) # 625 PAs

	# examine the merge
	sum(is.na(wdpaturn$beta_sor)) # 228 values missing
	sum(is.na(wdpaturn$flost)) # 228 values missing
	sum(is.na(wdpaturn$fgained)) # 228 values missing
	wdpaturn[is.na(wdpaturn$beta_sor), c('lat', 'lon', 'country', 'sub_loc', 'beta_sor')] # seem to be outside the richness projection. not sure why they were in the intersection with the climate grid
	
	# remove missing values
	wdpaturn <- wdpaturn[!is.na(wdpaturn$beta_sor),]
	
# average turnover within each PA and period
# important when PAs cross multiple grid cells
wdpaturnave1 <- Hmisc::summarize(wdpaturn[,c('beta_sor', 'prop_grid')], by=llist(wdpapolyID=wdpaturn$wdpapolyID), FUN=wmean, stat.name='beta_sor')
	dim(wdpaturnave1)
	length(unique(wdpaturnave1$wdpapolyID)) # 586 PAs (lost a few)
wdpaturnave2 <- Hmisc::summarize(wdpaturn[,c('fgained', 'prop_grid')], by=llist(wdpapolyID=wdpaturn$wdpapolyID), FUN=wmean, stat.name='fgained')
	dim(wdpaturnave2)
	length(unique(wdpaturnave2$wdpapolyID)) # 586 PAs (lost a few)
wdpaturnave3 <- Hmisc::summarize(wdpaturn[,c('flost', 'prop_grid')], by=llist(wdpapolyID=wdpaturn$wdpapolyID), FUN=wmean, stat.name='flost')
	dim(wdpaturnave3)
	length(unique(wdpaturnave3$wdpapolyID)) # 586 PAs (lost a few)
wdpaturnave <- merge(wdpaturnave1, wdpaturnave2)
wdpaturnave <- merge(wdpaturnave, wdpaturnave3)
wdpaturnave <- merge(wdpaturnave, wdpa[!duplicated(wdpa$wdpapolyID),]) # merge in other data
	dim(wdpaturnave) # 586
	
	# examine MLPA MPAs
	wdpaturnave[wdpaturnave$sub_loc=='US-CA' & wdpaturnave$mang_auth=="California Department of Fish and Game",] # all rows
	summary(wdpaturnave$beta_sor[wdpaturnave$sub_loc=='US-CA' & wdpaturnave$mang_auth=="California Department of Fish and Game"]) # Sorenson similarity
	summary(wdpaturnave$flost[wdpaturnave$sub_loc=='US-CA' & wdpaturnave$mang_auth=="California Department of Fish and Game"]) # Fraction species lost
	summary(wdpaturnave$fgained[wdpaturnave$sub_loc=='US-CA' & wdpaturnave$mang_auth=="California Department of Fish and Game"]) # Fraction species gained

##############################################
# Calculate turnover within MPA networks
##############################################
load('data/presmap.RData')
wdpa <- read.csv('data/wdpa_cov_by_grid0.25.csv', row.names=1)

# pick a MPA network
grids <- wdpa[wdpa$sub_loc=='US-CA' & wdpa$mang_auth=="California Department of Fish and Game", c('lat', 'lon')] #MLPA
	dim(grids)

# select initial and final timeperiod for these grids
thesespp <- presmap[paste(presmap$lat, presmap$lon) %in% paste(grids$lat, grids$lon) & presmap$period %in% c('2006-2020', '2081-2100'),]
thesespp <- thesespp[thesespp$pres,] # trim to present spp
	dim(thesespp) # 4858
	sort(unique(thesespp$region))
#	thesespp <- thesespp[thesespp$region == 'AFSC_WCTri',] # trim to one region
	thesespp <- thesespp[thesespp$region == 'NWFSC_WCAnn',] # trim to one region
	dim(thesespp) # 4858
	
# trim to unique species in each period
thesespp <- thesespp[!duplicated(thesespp[,c('sppocean', 'period')]),]
	dim(thesespp) # 180

# evaluate turnover
turnover(thesespp[,c('region', 'lat', 'lon', 'period', 'sppocean', 'pres')])

################
# Evaluate change
################

# Richness change in PAs
hist(wdparichtrend$richtrend) # nicely centered around 0. perhaps a slightly longer tail to the right (positive)
summary(wdparichtrend$richtrend) # mean gain/loss of species (spp/year): -0.04
summary(wdparichtrend$richtrend)*(2100-2006) # mean gain/loss of species (spp): -3 species, but -70 to +79
summary(wdparichtrend$richchange) # mean gain/loss of species (fraction): -3.5% (-85% to +400%)

# Richness change in grid cells
summary(richtrend$richchange) # mean gain/loss of species (fraction): -0.5% (-85% to +300%) (shouldn't upper bound match PAs?)


# Are richness trends within PAs more negative or more positive than the entire continental shelf?
t.test(richtrend$trend, wdparichtrend$richtrend) # trend in PAs is less negative than on average

# Is richness within PAs more likely to decline or increase than the entire shelf?
wdpainc = table(wdparichtrend$richtrend>0)
gridinc = table(richtrend$trend>0)
prop.test(x=rbind(wdpainc, gridinc)) # PAs are more likely to increase than on average across the shelf (46% vs. 36%)