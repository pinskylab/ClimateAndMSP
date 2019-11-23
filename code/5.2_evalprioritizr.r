# Evaluate the prioritzr solutions



############
## Flags
############

# choose the rcps (get to choose two)
rcps <- c(26, 85)

# choose the climate models to use for evaluating future planning (others for planning)
#bcc-csm1-1-m, bcc-csm1-1, CanESM2, CCSM4, CESM1-CAM5, CNRM-CM5, GFDL-CM3, GFDL-ESM2M, GFDL-ESM2G, GISS-E2-R, GISS-E2-H, IPSL-CM5A-LR, IPSL-CM5A-MR, MIROC-ESM, MPI-ESM-LR, NorESM1-ME
gcminds <- c(5, 6, 7, 11, 12, 13, 15, 16) # the complement of those in 5.1_prioritizr.r


# CMSP goals
consgoal <- 0.1 # proportion of presences to capture in conservation
energygoal <- 0.2 # proportion of NPV
fishgoal <- 0.5 # proportion of biomass

# choose regions for these runs (for reading in)
myregs <- c('ebs', 'goa', 'bc', 'wc', 'gmex', 'seus', 'neus', 'maritime', 'newf')

# periods of time
periods <- c('2007-2020', '2021-2040', '2041-2060', '2061-2080', '2081-2100')

# oceans
oceans <- c('Atl', 'Pac')

# choose names for writing out
runname1out <- 'hist_all'
runname2out <- '2per_all'


######################
## Helper functions
######################
require(RColorBrewer)
require(data.table)
require(maps)
require(ggplot2)

lu <- function(x) length(unique(x))



#######################################################
# Read in results and simple prep
#######################################################
# read in the prioritzr solutions
consplan1s <- vector('list', length(myregs))
for(i in 1:length(consplan1s)) consplan1s[[i]] <- fread(paste0('output/prioritizr_runs/solution_hist_', myregs[i], '.csv'), drop=1)
consplan2s <- vector('list', length(myregs))
for(i in 1:length(consplan2s)) consplan2s[[i]] <- fread(paste0('output/prioritizr_runs/solution_2per_', myregs[i], '.csv'), drop=1)

names(consplan1s) <- myregs
names(consplan2s) <- myregs


# loads presence/absence and biomass data for each species/model/rcp. very slow.
# not needed unless comparing plans against species distributions
if(!(length(rcps) %in% c(1,2))){
    stop('rcp must be length 1 or 2')
}
for (i in 1:length(rcps)){
    print(paste0('rcp', rcps[i]))
    
    # code if files are saved w/out timeperiod in the name
    # for(j in 1:length(oceans)){
    #     print(paste0('ocean', oceans[j]))
    #     prestemp <- fread(cmd = paste0('gunzip -c temp/presmap_', oceans[j], '_rcp', rcps[i], '.csv.gz'), drop = 1)
    #     biotemp <- fread(cmd = paste0('gunzip -c temp/biomassmap_', oceans[j], '_rcp', rcps[i], '.csv.gz'), drop = 1)
    # 
    #     # fix model #s for pres/abs projections (extras are in position 12 and 16)
    #     # DELETE THIS AFTER UPDATING 5.0.1_process_species_proj.r to fix this problem in that script
    #     prestemp <- prestemp[!(model %in% c(12, 16)), ]
    #     prestemp[model > 12 & model < 16, model := model - 1L]
    #     prestemp[model > 16, model := model - 2L]
    #     
    #     if(i == 1 & j == 1){
    #         presmap <- prestemp
    #         biomassmap <- biotemp
    #     } else {
    #         presmap <- rbind(presmap, prestemp)
    #         biomassmap <- rbind(biomassmap, biotemp)
    #     }
    # }
        
    # code if files are saved by time period
    for(j in 1:length(oceans)){
        for(k in 1:length(periods)){
            print(paste0('ocean', oceans[j], ' period', periods[k]))
            prestemp <- fread(cmd = paste0('gunzip -c temp/presmap_', oceans[j], '_rcp', rcps[i], '_', periods[k], '.csv.gz'), drop = 1)
            biotemp <- fread(cmd = paste0('gunzip -c temp/biomassmap_', oceans[j], '_rcp', rcps[i], '_', periods[k], '.csv.gz'), drop = 1)
    
            if(i == 1 & j == 1 & k == 1){
                presmap <- prestemp
                biomassmap <- biotemp
            } else {
                presmap <- rbind(presmap, prestemp)
                biomassmap <- rbind(biomassmap, biotemp)
            }
        }
    }
}
rm(prestemp, biotemp)
    
    # optional plots to examine presmap
    # thisspp <- 'gadus morhua'
    # thisspp <- 'urophycis chuss'
    # ggplot(presmap[spp == thisspp & model == 5 & rcp == 85], 
    #        aes(x = longrid, y = latgrid, color = poccur > 0.5)) +
    #     geom_point(size = 0.2) +
    #     facet_wrap(~ year_range, ncol = 2)
    # 
    # ggplot(biomassmap[spp == thisspp & model == 5 & rcp == 85], 
    #        aes(x = longrid, y = latgrid, color = biomass)) +
    #     geom_point(size = 0.2) +
    #     facet_wrap(~ year_range, ncol = 2)
    
# poccur threshold: how high does the probability of occurrence in the projections need to be to consider the species "present"?
# use the thresholds calculated during model fitting from Morley et al. 2018 PLOS ONE
poccurthresh <- fread('https://raw.githubusercontent.com/pinskylab/project_velocity/master/output/modeldiag_Nov2017_fitallreg_2017.csv', drop = 1)[, .(sppocean, thresh.kappa)]

# definition of planning goals by region
sppfiles <- list.files(path = 'output/prioritizr_runs/', pattern = 'spp_*', full.names = TRUE)
spps <- fread(sppfiles[1], drop = 1)
spps[, region := gsub('/|output|prioritizr_runs|spp_|\\.csv', '', sppfiles[1])]
for(i in 2:length(sppfiles)){
    temp <- fread(sppfiles[i], drop = 1)
    temp[, region := gsub('/|output|prioritizr_runs|spp_|\\.csv', '', sppfiles[i])]
    spps <- rbind(spps, temp)
}
rm(temp)

# region definitions
regiongrid <- fread(cmd = 'gunzip -c output/region_grid.csv.gz', drop = 1)


# add a zone indicator to consplans
for(i in 1:length(consplan1s)){
    consplan1s[[i]][ , zone := as.numeric(NA)]
    consplan1s[[i]][solution_1_conservation == 1 , zone := 1]
    consplan1s[[i]][solution_1_fishery == 1 , zone := 2]
    consplan1s[[i]][solution_1_energy == 1 , zone := 3]
    consplan1s[[i]][solution_1_conservation ==  0 & solution_1_fishery == 0 & solution_1_energy == 0 , zone := 4]
    consplan1s[[i]][, reg := myregs[i]]
}
for(i in 1:length(consplan2s)){
    consplan2s[[i]][ , zone <- as.numeric(NA)]
    consplan2s[[i]][solution_1_conservation == 1 , zone := 1]
    consplan2s[[i]][solution_1_fishery == 1 , zone := 2]
    consplan2s[[i]][solution_1_energy == 1 , zone := 3]
    consplan2s[[i]][solution_1_conservation ==  0 & solution_1_fishery == 0 & solution_1_energy == 0 , zone := 4]
    consplan2s[[i]][, reg := myregs[i]]
}


# Fix lon in regiongrid to match presmap (-360 to 0)
regiongrid[longrid > 0, longrid := longrid - 360] 

# Add region information to presmap
setkey(presmap, latgrid, longrid)
setkey(regiongrid, latgrid, longrid)
presmap <- merge(presmap, regiongrid[, .(latgrid, longrid, region)], all.x = TRUE) # add region information
if(presmap[is.na(region) & !duplicated(presmap[,.(latgrid, longrid)]), .N] != 0){ # 0 missing region: good!
    stop('presmap is missing >0 regions')
}
# presmap[is.na(region) & !duplicated(presmap[,.(latgrid, longrid)]), ]
# presmap[is.na(region) & !duplicated(presmap[,.(latgrid, longrid)]), plot(longrid, latgrid)]

# Add region information to biomassmap
setkey(biomassmap, latgrid, longrid)
setkey(regiongrid, latgrid, longrid)
biomassmap <- merge(biomassmap, regiongrid[, .(latgrid, longrid, region)], all.x = TRUE) # add region information
if(biomassmap[is.na(region) & !duplicated(biomassmap[,.(latgrid, longrid)]), .N] != 0){ # 0 missing region: good!
    stop('biomassmap is missing >0 regions')
}


# Add poccur threshold to presmap
poccurthresh[, ocean := gsub('.*_', '', sppocean)]
poccurthresh[, spp := gsub('_Atl|_Pac', '', sppocean)]

presmapPac <- merge(presmap[region %in% c('ebs', 'goa', 'bc', 'wc'), ], poccurthresh[ocean == 'Pac', .(spp, thresh.kappa)], by = 'spp') # have to do Atl and Pac separately since some species are in both regions but use different models
presmapAtl <- merge(presmap[region %in% c('gmex', 'seus', 'neus', 'maritime', 'newf'), ], poccurthresh[ocean == 'Atl', .(spp, thresh.kappa)], by = 'spp')
if(nrow(presmap) == nrow(presmapPac) + nrow(presmapAtl)){
    presmap <- rbind(presmapPac, presmapAtl)
    rm(presmapPac, presmapAtl)
} else {
    stop('merge of poccurthesh and presmap did not work')
}

# Fix a species name
presmap[spp == 'theragra chalcogramma', spp := 'gadus chalcogrammus']
biomassmap[spp == 'theragra chalcogramma', spp := 'gadus chalcogrammus']


####################################################################################
## Basic maps of solutions
####################################################################################

# plot map of selected grids (consplan #1 and #2)
cexs = 0.5 # to adjust
# quartz(width=5, height=18)
pdf(width = 5, height = 18, file = 'figures/planmaps.pdf')
par(mfrow = c(length(myregs),2), mai = c(0.5,0.5,0.3, 0.1), las = 1, mgp = c(2,1,0))
for(i in 1:length(consplan1s)){
    xlims <- range(consplan1s[[i]]$longrid)
    ylims <- range(consplan1s[[i]]$latgrid)
	j <- consplan1s[[i]]$solution_1_conservation == 1
	plot(consplan1s[[i]]$longrid[j], consplan1s[[i]]$latgrid[j], pch=16, 
	     col='red', cex=cexs, main='Historical only', xlab='', ylab='', xlim=xlims, ylim=ylims) # conservation
	j <- consplan1s[[i]]$solution_1_fishery == 1
	points(consplan1s[[i]]$longrid[j], consplan1s[[i]]$latgrid[j], pch=16, col='blue', cex=cexs) # fishery
	j <- consplan1s[[i]]$solution_1_energy == 1
	points(consplan1s[[i]]$longrid[j], consplan1s[[i]]$latgrid[j], pch=16, col='green', cex=cexs) # energy
	j <- consplan1s[[i]]$solution_1_conservation == 0 & 
	    consplan1s[[i]]$solution_1_fishery == 0 & consplan1s[[i]]$solution_1_energy == 0
	points(consplan1s[[i]]$longrid[j], consplan1s[[i]]$latgrid[j], pch=16, col='grey', cex=cexs) # all others
	
	j <- consplan2s[[i]]$solution_1_conservation == 1
	plot(consplan2s[[i]]$longrid[j], consplan2s[[i]]$latgrid[j], pch=16, col='red', cex=cexs, 
	     main='Two periods', xlab='', ylab='', xlim=xlims, ylim=ylims) # conservation
	j <- consplan2s[[i]]$solution_1_fishery == 1
	points(consplan2s[[i]]$longrid[j], consplan2s[[i]]$latgrid[j], pch=16, col='blue', cex=cexs) # fishery
	j <- consplan2s[[i]]$solution_1_energy == 1
	points(consplan2s[[i]]$longrid[j], consplan2s[[i]]$latgrid[j], pch=16, col='green', cex=cexs) # energy
	j <- consplan2s[[i]]$solution_1_conservation == 0 & consplan2s[[i]]$solution_1_fishery == 0 & 
	    consplan2s[[i]]$solution_1_energy == 0
	points(consplan2s[[i]]$longrid[j], consplan2s[[i]]$latgrid[j], pch=16, col='grey', cex=cexs) # all others
}

legend('topright', legend=c('Conservation', 'Fishery', 'Energy'), col=c('red', 'blue', 'green'), pch=rep(16,3), cex=0.7, title='Zones', bty='n')


dev.off()


############################
## Basic stats on the plans
############################
# have to read in the plans above (but don't need to read in the species distributions)
# climGrid <- read.csv('data/climGrid.csv', row.names=1) # for temperature climatology and depth

# proportion of each region in each zone
grids <- rbindlist(consplan1s)[, .(ngrid = .N), by = c('reg', 'zone')]
grids[, total := sum(ngrid), by  = 'reg'][, prop := round(ngrid/total, 2)]
setkey(grids, reg, zone)
grids
grids[ , .(meanprop = mean(prop), sdprop = sd(prop), seprop = sd(prop)/sqrt(.N)), by = 'zone']

# latitudinal range by zone in each plan
# do 2per plans span a wider latitudinal range?
#	a <- lapply(pusplan1s, FUN=function(x) aggregate(list(latrng_histonly=x$lat), by=list(zone=x$zone), FUN=function(x) max(x)-min(x)))
#	b <- lapply(pusplan2s, FUN=function(x) aggregate(list(latrng_2per=x$lat), by=list(zone=x$zone), FUN=function(x) max(x)-min(x)))
	a <- lapply(consplan1s, FUN=function(x) x[ , .(latsd_hist = sd(latgrid), reg = unique(reg)), by = zone])
	b <- lapply(consplan2s, FUN=function(x) x[ , .(latsd_2per = sd(latgrid), reg = unique(reg)), by = zone])
	a <- do.call('rbind',a)
	b <- do.call('rbind',b)
	d <- merge(a,b)
	d[d$zone==1,] # conservation
	d[d$zone==2,] # fishing
	d$diff <- d$latsd_2per - d$latsd_hist # positive if 2per spans a wider range
	summary(d$diff)
	wilcox.test(d$diff[d$zone==1])
	wilcox.test(d$diff[d$zone==2])
	wilcox.test(d$diff[d$zone==3])
		# note: zone 1 conservation, 2 fishery, 3 energy
	
# temperature and depth range by zone in each plan
	# add in climatology temperature
# pusplan1swclim <- pusplan1s
# pusplan2swclim <- pusplan2s
# for(i in 1:length(pusplan1swclim)){
# 	pusplan1swclim[[i]] <- merge(pusplan1s[[i]], climGrid[,c('lat', 'lon', 'bottemp.clim.int', 'surftemp.clim.int', 'depth')], all.x=TRUE)
# 	pusplan2swclim[[i]] <- merge(pusplan2s[[i]], climGrid[,c('lat', 'lon', 'bottemp.clim.int', 'surftemp.clim.int', 'depth')], all.x=TRUE)
# }
# 
# 	# do 2per plans span a wider surface temperature range?
# #	a <- lapply(pusplan1swclim, FUN=function(x) aggregate(list(temprng_histonly=x$surftemp.clim.int), by=list(zone=x$zone), FUN=function(x) max(x, na.rm=TRUE)-min(x, na.rm=TRUE)))
# #	b <- lapply(pusplan2swclim, FUN=function(x) aggregate(list(temprng_2per=x$surftemp.clim.int), by=list(zone=x$zone), FUN=function(x) max(x, na.rm=TRUE)-min(x, na.rm=TRUE)))
# 	a <- lapply(pusplan1swclim, FUN=function(x) aggregate(list(tempsd_histonly=x$surftemp.clim.int), by=list(zone=x$zone), FUN=function(x) sd(x, na.rm=TRUE)))
# 	b <- lapply(pusplan2swclim, FUN=function(x) aggregate(list(tempsd_2per=x$surftemp.clim.int), by=list(zone=x$zone), FUN=function(x) sd(x, na.rm=TRUE)))
# 	a <- do.call('rbind',a)
# 	b <- do.call('rbind',b)
# 	a$reg <- gsub('.[[:digit:]]', '', row.names(a))
# 	b$reg <- gsub('.[[:digit:]]', '', row.names(b))
# 	d <- merge(a,b)
# 	d[d$zone==2,] # conservation
# 	d[d$zone==3,] # fishing
# 	d$diff <- d$temprng_2per - d$temprng_histonly # positive if 2per spans a wider range
# 	summary(d$diff)
# 	wilcox.test(d$diff[d$zone==2])
# 	wilcox.test(d$diff[d$zone==3])
# 
# 	# do 2per plans span a wider bottom temperature range?
# #	a <- lapply(pusplan1swclim, FUN=function(x) aggregate(list(temprng_histonly=x$bottemp.clim.int), by=list(zone=x$zone), FUN=function(x) max(x, na.rm=TRUE)-min(x, na.rm=TRUE)))
# #	b <- lapply(pusplan2swclim, FUN=function(x) aggregate(list(temprng_2per=x$bottemp.clim.int), by=list(zone=x$zone), FUN=function(x) max(x, na.rm=TRUE)-min(x, na.rm=TRUE)))
# 	a <- lapply(pusplan1swclim, FUN=function(x) aggregate(list(tempsd_histonly=x$bottemp.clim.int), by=list(zone=x$zone), FUN=function(x) sd(x, na.rm=TRUE)))
# 	b <- lapply(pusplan2swclim, FUN=function(x) aggregate(list(tempsd_2per=x$bottemp.clim.int), by=list(zone=x$zone), FUN=function(x) sd(x, na.rm=TRUE)))
# 	a <- do.call('rbind',a)
# 	b <- do.call('rbind',b)
# 	a$reg <- gsub('.[[:digit:]]', '', row.names(a))
# 	b$reg <- gsub('.[[:digit:]]', '', row.names(b))
# 	d <- merge(a,b)
# 	d[d$zone==2,] # conservation
# 	d[d$zone==3,] # fishing
# 	d$diff <- d$temprng_2per - d$temprng_histonly # positive if 2per spans a wider range
# 	summary(d$diff)
# 	wilcox.test(d$diff[d$zone==2])
# 	wilcox.test(d$diff[d$zone==3])
# 
# 
# 	# do 2per plans span a wider depth range?
# 	a <- lapply(pusplan1swclim, FUN=function(x) aggregate(list(depthrng_histonly=x$depth), by=list(zone=x$zone), FUN=function(x) max(x, na.rm=TRUE)-min(x, na.rm=TRUE)))
# 	b <- lapply(pusplan2swclim, FUN=function(x) aggregate(list(depthrng_2per=x$depth), by=list(zone=x$zone), FUN=function(x) max(x, na.rm=TRUE)-min(x, na.rm=TRUE)))
# 	a <- do.call('rbind',a)
# 	b <- do.call('rbind',b)
# 	a$reg <- gsub('.[[:digit:]]', '', row.names(a))
# 	b$reg <- gsub('.[[:digit:]]', '', row.names(b))
# 	d <- merge(a,b)
# 	d[d$zone==2,] # conservation
# 	d[d$zone==3,] # fishing
# 	d$diff <- d$temprng_2per - d$temprng_histonly # positive if 2per spans a wider range
# 	summary(d$diff)
# 	wilcox.test(d$diff[d$zone==2])
# 	wilcox.test(d$diff[d$zone==3])
# 
	
######################################################################
# Analyze plans against each climate model projection                #
# across each rcp                                                    #
######################################################################

goalsbyzonebymod1 <- vector('list', length(myregs)) # to hold # presences and total biomass in conservation/fishing zones for each species and each projection, hist-only cmsp
goalsbyzonebymod2 <- vector('list', length(myregs)) # for 2per cmsp
goalsmetbymod1 <- vector('list', length(myregs)) # summary of # goals met by year_range/model/rcp, for hist cmsp
goalsmetbymod2 <- vector('list', length(myregs))

# set up models for evaluation
modelsp <- sort(unique(presmap$model)) # species in the pres/abs models
modelsb <- sort(unique(biomassmap$model)) # species in the biomass models
if(any(modelsp != modelsb)) stop('biomass and presmap models do not match')

# set up time periods for evaluation
periods <- sort(unique(presmap$year_range))

# evaluate # targets met in each time period in each model
for(i in 1:length(myregs)){
	print(i)
	
	# calculate total presences by species and time-period for conservation
	totalpres <- as.data.table(expand.grid(spp=spps[!grepl('fishery|energy', name) & region == myregs[i], spp], year_range = periods, model = modelsp, rcp=rcps)) # grid of all spp and time periods and models and rcps
	if(totalpres[, .(lu(spp)*lu(year_range)*lu(model)*lu(rcp))] != nrow(totalpres)) stop(paste('i=', i, 'and rows in totalpres are not right'))
	
	totprescalcs <- presmap[region == myregs[i], .(total = sum(poccur > poccurthresh)), by = c("spp", "year_range", "model", "rcp")]
        dim(totprescalcs)
    
	totalpres <- merge(totalpres, totprescalcs, all.x=TRUE) # make sure all spp and time periods and models and rcps represented
	totalpres[is.na(total), total := 0]
	totalpres[ , zone := 1] # for conservation

	# merge plan zones into presmap and calculate # presences
	temp1 <- merge(presmap[region == myregs[i], ], consplan1s[[i]][zone == 1, .(latgrid, longrid, zone)], by=c('latgrid', 'longrid')) # only keep the conserved zones for goals related to climate
		dim(temp1) # a data.table
	temp2 <- merge(presmap[region == myregs[i], ], consplan2s[[i]][zone == 1, .(latgrid, longrid, zone)], by=c('latgrid', 'longrid'))
	presbyzone1 <- temp1[, .(zonetotal = sum(poccur > poccurthresh)), by = c('spp', 'year_range', 'zone', 'model', 'rcp')]
		dim(presbyzone1) # 19392 (ebs)
	presbyzone2 <- temp2[, .(zonetotal = sum(poccur > poccurthresh)), by = c('spp', 'year_range', 'zone', 'model', 'rcp')]
		dim(presbyzone2)

	# add totals for pres across all zones
	setkey(totalpres, spp, year_range, model, rcp, zone)
	setkey(presbyzone1, spp, year_range, model, rcp, zone)
	setkey(presbyzone2, spp, year_range, model, rcp, zone)
	presbyzonebymod1 <- merge(presbyzone1, totalpres, all.y=TRUE)
	if(presbyzonebymod1[, sum(is.na(zonetotal))] > 0) presbyzonebymod1[is.na(zonetotal), zonetotal:= 0]
	    dim(presbyzonebymod1) # 5728 (ebs), 9024 (bc)
	    presbyzonebymod1[, length(unique(spp))] # 179 species (ebs), 141 (bc)

	presbyzonebymod2 <- merge(presbyzone2, totalpres, all.y=TRUE)
	if(presbyzonebymod2[, sum(is.na(zonetotal))] > 0) presbyzonebymod2[is.na(zonetotal), zonetotal:=0]
	    dim(presbyzonebymod2) # 
    	presbyzonebymod2[, length(unique(spp))] # 

	# calculate proportion of presences
	presbyzonebymod1[, prop := zonetotal/total]
	presbyzonebymod2[, prop := zonetotal/total]

	# calculate total biomass by species and time-period for fishery species
	totalbio <- as.data.table(expand.grid(spp = spps[grepl('fishery', name) & region == myregs[i], spp], year_range = periods, model = modelsp, rcp = rcps)) # grid of all spp and time periods and models and rcps
	if(totalbio[, .(lu(spp)*lu(year_range)*lu(model)*lu(rcp))] != nrow(totalbio)) stop(paste('i=', i, 'and rows in totalbio are not right'))
	
	totbiocalcs <- biomassmap[region == myregs[i], .(total = sum(biomass)), by = c("spp", "year_range", "model", "rcp")]
	    dim(totbiocalcs)
	
	totalbio <- merge(totalbio, totbiocalcs, all.x=TRUE) # make sure all spp and time periods and models and rcps represented
	totalbio[is.na(total), total := 0]
	totalbio[ , zone := 2] # for fisheries
	
	# merge plan zones into biomassmap and calculate total biomass
	# note: biobyzone* includes species that aren't fishery targets in this region
	temp1 <- merge(biomassmap[region == myregs[i], ], consplan1s[[i]][zone == 2, .(latgrid, longrid, zone)], by=c('latgrid', 'longrid')) # only keep the fishery zones for goals related to climate
    	dim(temp1) # a data.table
	temp2 <- merge(biomassmap[region == myregs[i], ], consplan2s[[i]][zone == 2, .(latgrid, longrid, zone)], by=c('latgrid', 'longrid'))
	biobyzone1 <- temp1[, .(zonetotal = sum(biomass)), by = c('spp', 'year_range', 'zone', 'model', 'rcp')]
	    dim(biobyzone1) # 19392 (ebs), 640 (bc)
	biobyzone2 <- temp2[, .(zonetotal = sum(biomass)), by = c('spp', 'year_range', 'zone', 'model', 'rcp')]
	    dim(biobyzone2)
	
	# add totals for biomass
	setkey(totalbio, spp, year_range, model, rcp, zone)
	setkey(biobyzone1, spp, year_range, model, rcp, zone)
	setkey(biobyzone2, spp, year_range, model, rcp, zone)
	biobyzonebymod1 <- merge(biobyzone1, totalbio, all.y=TRUE, by = c('spp', 'year_range', 'model', 'rcp', 'zone'))
	if(biobyzonebymod1[, sum(is.na(zonetotal))] > 0) biobyzonebymod1[is.na(zonetotal), zonetotal := 0]
	    nrow(biobyzonebymod1) # 320 (ebs w/ 2 year_ranges), 640 (bc w/ 2 year_ranges)
	    biobyzonebymod1[, length(unique(spp))] # 10 species
	
	biobyzonebymod2 <- merge(biobyzone2, totalbio, all.y=TRUE)
	if(biobyzonebymod2[, sum(is.na(zonetotal))] > 0) biobyzonebymod2[is.na(zonetotal), zonetotal := 0]
	    dim(biobyzonebymod2) # 
	biobyzonebymod2[, length(unique(spp))] # 
	
	# calculate proportion of biomass
	biobyzonebymod1[, prop := zonetotal/total]
	biobyzonebymod2[, prop := zonetotal/total]
	
	

	# count a goal as met if the species is not present
	presbyzonebymod1[total == 0, prop := 1]
	presbyzonebymod2[total == 0, prop := 1]
	biobyzonebymod1[total == 0, prop := 1]
	biobyzonebymod2[total == 0, prop := 1]

	# mark where goals met for conservation or fishery
	presbyzonebymod1[, metgoal := FALSE]
	presbyzonebymod2[, metgoal := FALSE]
	expr <- paste('zone == 1 & prop >=', consgoal) # the only way to dynamically construct a vector with which to choose rows (http://r.789695.n4.nabble.com/Idiomatic-way-of-using-expression-in-i-td4711229.html)
	presbyzonebymod1[eval(parse(text=expr)), metgoal := TRUE]
	presbyzonebymod2[eval(parse(text=expr)), metgoal := TRUE]
	
	biobyzonebymod1[, metgoal := FALSE]
	biobyzonebymod2[, metgoal := FALSE]
	expr2 <- paste('zone == 2 & prop >=', fishgoal)
	biobyzonebymod1[eval(parse(text=expr2)), metgoal := TRUE]
	biobyzonebymod2[eval(parse(text=expr2)), metgoal := TRUE]

    # combine conservation and fishery goals
	goalsbyzonebymod1[[i]] <- rbind(presbyzonebymod1, biobyzonebymod1)
	goalsbyzonebymod2[[i]] <- rbind(presbyzonebymod2, biobyzonebymod2)
	
	# calculate number goals met in each timeperiod
	numgoals <- totalpres[, length(unique(spp))] + totalbio[, length(unique(spp))] # total number of conservation and fishery goals
	goalsmetbymod1[[i]] <- goalsbyzonebymod1[[i]][, .(nmet = sum(metgoal)), by=c('year_range','model','rcp')]
	goalsmetbymod1[[i]]$mid <- sapply(strsplit(as.character(goalsmetbymod1[[i]]$year_range), split='-'), FUN=function(x) mean(as.numeric(x))) # find mid-year
	goalsmetbymod1[[i]]$pmet <- goalsmetbymod1[[i]]$nmet/numgoals

	goalsmetbymod2[[i]] <- goalsbyzonebymod2[[i]][, .(nmet=sum(metgoal)), by = c('year_range','model','rcp')]
	goalsmetbymod2[[i]]$mid <- sapply(strsplit(as.character(goalsmetbymod2[[i]]$year_range), split='-'), FUN=function(x) mean(as.numeric(x)))
	goalsmetbymod2[[i]]$pmet <- goalsmetbymod2[[i]]$nmet/numgoals

	# set region
	goalsbyzonebymod1[[i]][, region := myregs[i]]
	goalsbyzonebymod2[[i]][, region := myregs[i]]
	goalsmetbymod1[[i]][, region := myregs[i]]
	goalsmetbymod2[[i]][, region := myregs[i]]
}

# clean up
rm(temp1, temp2, presbyzone1, presbyzone2, presbyzonebymod1, presbyzonebymod2, biobyzone1, biobyzone2, biobyzonebymod1, biobyzonebymod2, totalbio, totalpres, totbiocalcs, totprescalcs)

# flatten lists into data.frames
goalsbyzonebymod1out <- rbindlist(goalsbyzonebymod1)
goalsbyzonebymod2out <- rbindlist(goalsbyzonebymod2)
goalsmetbymod1out <- rbindlist(goalsmetbymod1)
goalsmetbymod2out <- rbindlist(goalsmetbymod2)

dim(goalsbyzonebymod1out)
dim(goalsbyzonebymod2out)
dim(goalsmetbymod1out)
dim(goalsmetbymod2out)

# simple check
test <- goalsmetbymod1out[, .(min = min(pmet), max = max(pmet)), by = region]
if(test[, any(min < 0)] | test[, any(max > 1)]) stop('some regions have pmet < 0 or > 1')

# write out
write.csv(goalsbyzonebymod1out, file = gzfile(paste0('output/goalsbyzonebymod_', runname1out, '.csv.gz')))
write.csv(goalsbyzonebymod2out, file = gzfile(paste0('output/goalsbyzonebymod_', runname2out, '.csv.gz')))
write.csv(goalsmetbymod1out, file = paste0('output/goalsmetbymod_', runname1out, '.csv'))
write.csv(goalsmetbymod2out, file = paste0('output/goalsmetbymod_', runname2out, '.csv'))


######################################################################
# Compare and plot the planning approaches against all climate models
######################################################################
goalsmetbymod1out <- fread(paste0('output/goalsmetbymod_', runname1out, '.csv'), drop = 1)
goalsmetbymod2out <- fread(paste0('output/goalsmetbymod_', runname2out, '.csv'), drop = 1)


mods <- sort(unique(goalsmetbymod1out$model))
rcps <- sort(unique(goalsmetbymod1out$rcp))
regnames<-c('E. Bering Sea', 'Gulf of AK', 'British Columbia', 'West Coast US', 'Gulf of Mexico', 'Southeast US', 'Northeast US', 'Maritimes', 'Newfoundland')

# add plan identifier
goalsmetbymod1out$plan <- 'histonly'
goalsmetbymod2out$plan <- '2per'

# combine
goalsmetbymod <- rbind(goalsmetbymod1out, goalsmetbymod2out)


# compare goals met across plans and time periods
goalmetbymodag <- with(rbind(goalsmetbymod1out, goalsmetbymod2out), 
                       aggregate(list(pmet = pmet), by = list(plan = plan, model = model, rcp = rcp, period = year_range), FUN = mean))
with(goalmetbymodag, aggregate(list(ave_goals = pmet), by = list(plan = plan, period = period), FUN = mean))
with(goalmetbymodag, aggregate(list(ave_goals = pmet), by = list(plan = plan, period = period), FUN = function(x) sd(x)/length(x)))
with(goalmetbymodag, aggregate(list(ave_goals = pmet), by = list(plan = plan, period = period), FUN = length))


# test whether histonly meets significantly fewer goals than 2per
	all(goalsmetbymod1out$model == goalsmetbymod2out$model & goalsmetbymod1out$period == goalsmetbymod2out$period & goalsmetbymod1out$rcp == goalsmetbymod2out$rcp) # same order?

	# simple t-test on number met (pseudoreplicates within regions)
	t.test(goalsmetbymod1out$nmet[goalsmetbymod1out$period=='2081-2100'], goalsmetbymod2out$nmet[goalsmetbymod2out$period=='2081-2100'])

	# simple t-test on proportions (should use logistic) (also pseudoreplicates within regions)
	t.test(goalsmetbymod1out$pmet[goalsmetbymod1out$period=='2081-2100'], goalsmetbymod2out$pmet[goalsmetbymod2out$period=='2081-2100'], paired=FALSE)

	# GLMM model (since proportion data) on 2006-2020
	require(lme4)
	require(car)
	mod <- glmer(cbind(nmet, nmet/pmet - nmet) ~ plan + (1|region/rcp/model), data=goalsmetbymod, family='binomial', subset=goalsmetbymod$period=='2006-2020')
	summary(mod)
	Anova(mod)

	# GLMM model (since proportion data) on 2041-2060
	require(lme4)
	require(car)
	goalsmetbymod <- rbind(goalsmetbymod1out, goalsmetbymod2out)
	mod <- glmer(cbind(nmet, nmet/pmet - nmet) ~ plan + (1|region/rcp/model), data=goalsmetbymod, family='binomial', subset=goalsmetbymod$period=='2041-2060')
	summary(mod)
	Anova(mod)

	# GLMM model (since proportion data) on 2081-2100
	goalsmetbymod <- rbind(goalsmetbymod1out, goalsmetbymod2out)
	mod <- glmer(cbind(nmet, nmet/pmet - nmet) ~ plan + (1|region/rcp/model), data=goalsmetbymod, family='binomial', subset=goalsmetbymod$period=='2081-2100')
	summary(mod)
	Anova(mod)

	# GLMM model (since proportion data) across all time periods
	goalsmetbymod <- rbind(goalsmetbymod1out, goalsmetbymod2out)
	mod <- glmer(cbind(nmet, nmet/pmet - nmet) ~ plan + (1|region/rcp/model/period), data=goalsmetbymod, family='binomial')
	summary(mod)
	Anova(mod)


# plot %goals met (solution #1 and #2) for non-planning models
	# thesemods <- mods[gcminds]; type <- 'testing' # for models not used in planning
	thesemods <- mods[-gcminds]; type <- 'training' # for models used in planning
	
	colmat <- t(col2rgb(brewer.pal(4, 'Paired')))
#	colmat <- t(col2rgb(brewer.pal(4, 'Dark2')))
	cols <- rgb(red=colmat[,1], green=colmat[,2], blue=colmat[,3], alpha=c(255,255,255,255), maxColorValue=255)
	if(type == 'testing') outfile <- paste('figures/prioritizr_allregs_goalsmetbymod_', runname1out, '&', runname2out, '.pdf', sep='')
	if(type == 'training') outfile <- paste('temp_figures/prioritizr_allregs_goalsmetbymod_', runname1out, '&', runname2out, '_training.pdf', sep='')
	    outfile
	ylims <- c(0,1)

	# quartz(width=6, height=6)
	pdf(width=6, height=6, file=outfile)
	par(mfrow=c(3,3), mai=c(0.3, 0.35, 0.2, 0.05), omi=c(0.2,0.2,0,0), cex.main=0.8)

	for(i in 1:length(myregs)){
		# each model/rcp for hist only plans
	    inds <- goalsmetbymod1out$model == thesemods[1] & goalsmetbymod1out$rcp == rcps[1] & goalsmetbymod1out$region==myregs[i]
		plot(goalsmetbymod1out$mid[inds], goalsmetbymod1out$pmet[inds], xlab='', ylab='', ylim=ylims, type='l', pch=16, las=1, col=cols[1], main=regnames[i])
		for(k in 1:length(thesemods)){
			for(j in 1:length(rcps)){
				if(!(k==1 & j==1)){
					inds <- goalsmetbymod1out$model == thesemods[k] & goalsmetbymod1out$rcp == rcps[j] & goalsmetbymod1out$region==myregs[i]
					points(goalsmetbymod1out$mid[inds], goalsmetbymod1out$pmet[inds], type='l', pch=16, col=cols[1])	
				}
			}
		}
		
        # each model/rcp for 2per plans
		for(k in 1:length(thesemods)){
			for(j in 1:length(rcps)){
				inds <- goalsmetbymod2out$model == thesemods[k] & goalsmetbymod2out$rcp == rcps[j] & goalsmetbymod2out$region==myregs[i]
				points(goalsmetbymod2out$mid[inds], goalsmetbymod2out$pmet[inds], type='l', pch=16, col=cols[3])
			}	
		}

		# average across models/rcps for hist and 2per plans
		inds <- goalsmetbymod1out$region==myregs[i] & goalsmetbymod1out$model %in% thesemods
		ensmean <- aggregate(list(nmet=goalsmetbymod1out$nmet[inds], pmet=goalsmetbymod1out$pmet[inds]), by=list(mid=goalsmetbymod1out$mid[inds]), FUN=mean)
		lines(ensmean$mid, ensmean$pmet, col=cols[2], lwd=3)
		
		inds <- goalsmetbymod2out$region==myregs[i] & goalsmetbymod2out$model %in% thesemods
		ensmean2 <- aggregate(list(nmet=goalsmetbymod2out$nmet[inds], pmet=goalsmetbymod2out$pmet[inds]), by=list(mid=goalsmetbymod2out$mid[inds]), FUN=mean)
		lines(ensmean2$mid, ensmean2$pmet, col=cols[4], lwd=3)
		
		# legend in plot 1
		if(i == 1) legend('bottomright', legend = c('hist', '2per'), col = cols[c(2,4)], lwd = 1)
	}
	mtext(side=1,text='Year',line=0.5, outer=TRUE)
	mtext(side=2,text='Fraction goals met',line=0.3, outer=TRUE)
	
	dev.off()
	
	
## probability of < 80% goals met by region, plan, and time period
	# calcs
	u80func <- function(x) return(sum(x<0.8)/length(x))	
	u80 <- with(goalsmetbymod, aggregate(list(u80=pmet), by=list(region=region, period=year_range, plan=plan), FUN=u80func))


	# examine
	u80[u80$period=='2041-2060',] # middle of century, by region
	u80[u80$period=='2081-2100',] # end of century, by region
	
	with(u80[u80$period=='2007-2020',], aggregate(list(meanu80=u80), by=list(plan=plan), FUN=mean)) # means by plan (across regions), start
	with(u80[u80$period=='2041-2060',], aggregate(list(meanu80=u80), by=list(plan=plan), FUN=mean)) # means by plan (across regions), middle of century
	with(u80[u80$period=='2081-2100',], aggregate(list(meanu80=u80), by=list(plan=plan), FUN=mean)) # means by plan (across regions), end of century

	a<-aggregate(list(max_u80=u80$u80), by=list(region=u80$region, plan=u80$plan), FUN=max) # max probability of being under 80% of goals met (look across time periods)
		mean(a$max_u80[a$plan=='2per'])
		mean(a$max_u80[a$plan=='histonly'])
	
	# plot
	ggplot(u80, aes(x = period, y = u80, group = plan, color = plan)) +
	    geom_line() + 
	    geom_point() +
	    facet_wrap(~ region, ncol = 3)
	    
	