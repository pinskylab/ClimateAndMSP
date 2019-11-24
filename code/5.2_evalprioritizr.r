# Evaluate the prioritzr solutions against the range of species projections
# Save the results: whether or not a plan meets a goal in a given model/rcp/time



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
#require(RColorBrewer)
require(data.table)
#require(maps)
#require(ggplot2)

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

# definition of planning features by region
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
	
	totprescalcs <- presmap[region == myregs[i], .(total = sum(poccur > thresh.kappa)), by = c("spp", "year_range", "model", "rcp")]
        dim(totprescalcs)
    
	totalpres <- merge(totalpres, totprescalcs, all.x=TRUE) # make sure all spp and time periods and models and rcps represented
	totalpres[is.na(total), total := 0]
	totalpres[ , zone := 1] # for conservation

	# merge plan zones into presmap and calculate # presences
	temp1 <- merge(presmap[region == myregs[i], ], consplan1s[[i]][zone == 1, .(latgrid, longrid, zone)], by=c('latgrid', 'longrid')) # only keep the conserved zones for goals related to climate
		dim(temp1) # a data.table
	temp2 <- merge(presmap[region == myregs[i], ], consplan2s[[i]][zone == 1, .(latgrid, longrid, zone)], by=c('latgrid', 'longrid'))
	presbyzone1 <- temp1[, .(zonetotal = sum(poccur > thresh.kappa)), by = c('spp', 'year_range', 'zone', 'model', 'rcp')]
		dim(presbyzone1) # 19392 (ebs)
	presbyzone2 <- temp2[, .(zonetotal = sum(poccur > thresh.kappa)), by = c('spp', 'year_range', 'zone', 'model', 'rcp')]
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


	    
	