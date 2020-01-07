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
consplan1s <- fread(paste0('output/prioritizr_runs/solution_hist_', myregs[1], '.csv'), drop=1)[, region := myregs[1]]
for(i in 2:length(myregs)) consplan1s <- rbind(consplan1s, fread(paste0('output/prioritizr_runs/solution_hist_', myregs[i], '.csv'), drop=1)[, region := myregs[i]])
consplan2s <- fread(paste0('output/prioritizr_runs/solution_2per_', myregs[1], '.csv'), drop=1)[, region := myregs[1]]
for(i in 2:length(myregs)) consplan2s <- rbind(consplan2s, fread(paste0('output/prioritizr_runs/solution_2per_', myregs[i], '.csv'), drop=1)[, region := myregs[i]])

# add a zone indicator to consplans
consplan1s[ , zone := as.numeric(NA)]
consplan1s[solution_1_conservation == 1 , zone := 1]
consplan1s[solution_1_fishery == 1 , zone := 2]
consplan1s[solution_1_energy == 1 , zone := 3]
consplan1s[solution_1_conservation ==  0 & solution_1_fishery == 0 & solution_1_energy == 0 , zone := 4]
consplan2s[ , zone <- as.numeric(NA)]
consplan2s[solution_1_conservation == 1 , zone := 1]
consplan2s[solution_1_fishery == 1 , zone := 2]
consplan2s[solution_1_energy == 1 , zone := 3]
consplan2s[solution_1_conservation ==  0 & solution_1_fishery == 0 & solution_1_energy == 0 , zone := 4]

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

# Fix lon in regiongrid to match presmap (-360 to 0)
regiongrid[longrid > 0, longrid := longrid - 360] 


	
######################################################################
# Analyze plans against each climate model projection                #
# across each rcp                                                    #
######################################################################

# evaluate # targets met in each time period/rcp, by species/model/region
for(k in 1:length(periods)){
    print(paste('Starting', periods[k]))
    # loads presence/absence and biomass data for each species/model/rcp for this timeperiod
    for (i in 1:length(rcps)){
        for(j in 1:length(oceans)){
            cat(paste0('\tLoading rcp', rcps[i], ' ocean', oceans[j], ' period', periods[k], '\n'))
            prestemp <- fread(cmd = paste0('gunzip -c temp/presmap_', oceans[j], '_rcp', rcps[i], '_', periods[k], '.csv.gz'), drop = 1)
            biotemp <- fread(cmd = paste0('gunzip -c temp/biomassmap_', oceans[j], '_rcp', rcps[i], '_', periods[k], '.csv.gz'), drop = 1)
            
            if(i == 1 & j == 1){
                presmap <- prestemp
                biomassmap <- biotemp
            } else {
                presmap <- rbind(presmap, prestemp)
                biomassmap <- rbind(biomassmap, biotemp)
            }
        }
    }
    rm(prestemp, biotemp)
    
    cat('\tProcessing data\n')

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
    
    # set up models for evaluation
    modelsp <- sort(unique(presmap$model)) # species in the pres/abs models
    modelsb <- sort(unique(biomassmap$model)) # species in the biomass models
    if(any(modelsp != modelsb)) stop('biomass and presmap models do not match')
    
    cat('\tRunning calculations\n')
    
    # calculate total presences by species/region for conservation
	totalpres <- merge(spps[!grepl('fishery|energy', name), .(spp, region, merge = 1)], as.data.table(expand.grid(model = modelsp, rcp=rcps, merge = 1)), allow.cartesian=TRUE) # grid of all spp-regions by models by rcps. fool data.table into doing a cross join with a dummy merge column
	totalpres[, year_range := periods[k]]
	if(totalpres[, .(lu(paste(spp, region))*lu(year_range)*lu(model)*lu(rcp))] != nrow(totalpres)) stop(paste('k=', k, 'and rows in totalpres are not right'))
	totalpres[, merge := NULL] # remove the dummy merge column
	
	totprescalcs <- presmap[, .(total = sum(poccur > thresh.kappa)), by = c("spp", "year_range", "model", "rcp", "region")]
        dim(totprescalcs)
    
	totalpres <- merge(totalpres, totprescalcs, all.x=TRUE, by = c('spp', 'year_range', 'rcp', 'model', 'region')) # make sure all spp and time periods and models and rcps represented
	totalpres[is.na(total), total := 0]
	totalpres[ , zone := 1] # for conservation

	
	# merge plan zones into presmap and calculate # presences
	temp1 <- merge(presmap, consplan1s[zone == 1, .(latgrid, longrid, zone)], by=c('latgrid', 'longrid')) # only keep the conserved zones for goals related to climate
		dim(temp1) # a data.table
	temp2 <- merge(presmap, consplan2s[zone == 1, .(latgrid, longrid, zone)], by=c('latgrid', 'longrid'))
	presbyzone1 <- temp1[, .(zonetotal = sum(poccur > thresh.kappa)), by = c('spp', 'year_range', 'zone', 'model', 'rcp', 'region')]
		dim(presbyzone1)
	presbyzone2 <- temp2[, .(zonetotal = sum(poccur > thresh.kappa)), by = c('spp', 'year_range', 'zone', 'model', 'rcp', 'region')]
		dim(presbyzone2)

	# add totals for pres
	setkey(totalpres, spp, year_range, model, rcp, region, zone)
	setkey(presbyzone1, spp, year_range, model, rcp, region, zone)
	setkey(presbyzone2, spp, year_range, model, rcp, region, zone)
	presbyzonebymod1 <- merge(presbyzone1, totalpres, all.y=TRUE)
	if(presbyzonebymod1[, sum(is.na(zonetotal))] > 0) presbyzonebymod1[is.na(zonetotal), zonetotal:= 0]
	    dim(presbyzonebymod1)
	    presbyzonebymod1[, .(nspp = length(unique(spp))), by = 'region']

	presbyzonebymod2 <- merge(presbyzone2, totalpres, all.y=TRUE)
	if(presbyzonebymod2[, sum(is.na(zonetotal))] > 0) presbyzonebymod2[is.na(zonetotal), zonetotal:=0]
	    dim(presbyzonebymod2) # 
	    presbyzonebymod2[, .(nspp = length(unique(spp))), by = 'region']
	    
	# calculate proportion of presences
	presbyzonebymod1[, prop := zonetotal/total]
	presbyzonebymod2[, prop := zonetotal/total]

	# calculate total biomass by species and time-period for fishery species
	totalbio <- merge(spps[grepl('fishery', name), .(spp, region, merge = 1)], as.data.table(expand.grid(model = modelsb, rcp=rcps, merge = 1)), allow.cartesian=TRUE) # grid of all fishery spp-regions by models by rcps. fool data.table into doing a cross join with a dummy merge column
	totalbio[, year_range := periods[k]]
	if(totalbio[, .(lu(paste(spp, region))*lu(year_range)*lu(model)*lu(rcp))] != nrow(totalbio)) stop(paste('k=', k, 'and rows in totalpres are not right'))
	totalbio[, merge := NULL] # remove the dummy merge column
	
	totbiocalcs <- biomassmap[, .(total = sum(biomass)), by = c("spp", "year_range", "model", "rcp", "region")]
	    dim(totbiocalcs)
	
	totalbio <- merge(totalbio, totbiocalcs, all.x=TRUE) # make sure all spp and time periods and models and rcps represented
	totalbio[is.na(total), total := 0]
	totalbio[ , zone := 2] # for fisheries
	
	# merge plan zones into biomassmap and calculate total biomass
	# note: biobyzone* includes species that aren't fishery targets in this region
	temp1 <- merge(biomassmap, consplan1s[zone == 2, .(latgrid, longrid, zone)], by=c('latgrid', 'longrid')) # only keep the fishery zones for goals related to climate
    	dim(temp1) # a data.table
	temp2 <- merge(biomassmap, consplan2s[zone == 2, .(latgrid, longrid, zone)], by=c('latgrid', 'longrid'))
	biobyzone1 <- temp1[, .(zonetotal = sum(biomass)), by = c('spp', 'year_range', 'zone', 'model', 'rcp', 'region')]
	    dim(biobyzone1)
	biobyzone2 <- temp2[, .(zonetotal = sum(biomass)), by = c('spp', 'year_range', 'zone', 'model', 'rcp', 'region')]
	    dim(biobyzone2)
	
	# add totals for biomass
	setkey(totalbio, spp, year_range, model, rcp, region, zone)
	setkey(biobyzone1, spp, year_range, model, rcp, region, zone)
	setkey(biobyzone2, spp, year_range, model, rcp, region, zone)
	biobyzonebymod1 <- merge(biobyzone1, totalbio, all.y=TRUE, by = c('spp', 'year_range', 'model', 'rcp', 'zone', 'region'))
	if(biobyzonebymod1[, sum(is.na(zonetotal))] > 0) biobyzonebymod1[is.na(zonetotal), zonetotal := 0]
	    nrow(biobyzonebymod1)
	    biobyzonebymod1[, .(nspp = length(unique(spp))), by = 'region'] # 10 species
	
	biobyzonebymod2 <- merge(biobyzone2, totalbio, all.y=TRUE)
	if(biobyzonebymod2[, sum(is.na(zonetotal))] > 0) biobyzonebymod2[is.na(zonetotal), zonetotal := 0]
	    dim(biobyzonebymod2) # 
	    biobyzonebymod2[, .(nspp = length(unique(spp))), by = 'region'] # 10 species
	
	# calculate proportion of biomass
	biobyzonebymod1[, prop := zonetotal/total]
	biobyzonebymod2[, prop := zonetotal/total]
	
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

	# mark a goal as not counting if the species is not present
	presbyzonebymod1[total == 0, metgoal := NA]
	presbyzonebymod2[total == 0, metgoal := NA]
	biobyzonebymod1[total == 0, metgoal := NA]
	biobyzonebymod2[total == 0, metgoal := NA]
	
	# combine conservation and fishery goals
	# add on to previous time-period's results if already present
	if(k == 1){
	    goalsbyzonebymod1 <- rbind(presbyzonebymod1, biobyzonebymod1)
	    goalsbyzonebymod2 <- rbind(presbyzonebymod2, biobyzonebymod2)
	} else {
	    goalsbyzonebymod1 <- rbind(goalsbyzonebymod1, presbyzonebymod1, biobyzonebymod1)
	    goalsbyzonebymod2 <- rbind(goalsbyzonebymod2, presbyzonebymod2, biobyzonebymod2)
	}
	

}

# clean up
rm(temp1, temp2, presbyzone1, presbyzone2, presbyzonebymod1, presbyzonebymod2, biobyzone1, biobyzone2, biobyzonebymod1, biobyzonebymod2, totalbio, totalpres, totbiocalcs, totprescalcs)

# set up strict goal accounting: species not present counts against meeting the goal
goalsbyzonebymod1[, metgoalstrict := metgoal]
goalsbyzonebymod2[, metgoalstrict := metgoal]
goalsbyzonebymod1[total == 0, metgoalstrict := FALSE]
goalsbyzonebymod2[total == 0, metgoalstrict := FALSE]

# summarize by region/timeperiod/model/rcp
# no need for a separate nmetstrict column because == nmet in all cases
goalsmetbymod1 <- goalsbyzonebymod1[, .(mid = mean(as.numeric(unlist(strsplit(as.character(year_range), split='-')))), 
                                        ngoals = sum(!is.na(metgoal)), nmet = sum(metgoal, na.rm = TRUE), 
                                        pmet = sum(metgoal, na.rm = TRUE)/sum(!is.na(metgoal)),
                                        ngoalsstrict = .N, pmetstrict = sum(metgoal, na.rm = TRUE)/.N),
                                    by = c('year_range', 'model', 'rcp', 'region')]
goalsmetbymod2 <- goalsbyzonebymod2[, .(mid = mean(as.numeric(unlist(strsplit(as.character(year_range), split='-')))), 
                                        ngoals = sum(!is.na(metgoal)), nmet = sum(metgoal, na.rm = TRUE), 
                                        pmet = sum(metgoal, na.rm = TRUE)/sum(!is.na(metgoal)),
                                        ngoalsstrict = .N, pmetstrict = sum(metgoal, na.rm = TRUE)/.N),
                                    by = c('year_range', 'model', 'rcp', 'region')]

# check
dim(goalsbyzonebymod1)
dim(goalsbyzonebymod2)
dim(goalsmetbymod1)
dim(goalsmetbymod2)

# quick summary of the number of regular vs. strict goals: how many scenarios have fewer regular goals?
goalsmetbymod1[, .(ndiff = sum(ngoals < ngoalsstrict), n = .N), by = 'region']
goalsmetbymod2[, .(ndiff = sum(ngoals < ngoalsstrict), n = .N), by = 'region'] # should be the same

# simple check
test <- goalsmetbymod1[, .(min = min(pmet), max = max(pmet)), by = region]
if(test[, any(min < 0)] | test[, any(max > 1)]) stop('some regions have pmet < 0 or > 1')

# label which models are for planning and which for testing
goalsbyzonebymod1[, modeltype := 'planning']
goalsbyzonebymod1[model %in% gcminds, modeltype := 'testing']
goalsbyzonebymod2[, modeltype := 'planning']
goalsbyzonebymod2[model %in% gcminds, modeltype := 'testing']
goalsmetbymod1[, modeltype := 'planning']
goalsmetbymod1[model %in% gcminds, modeltype := 'testing']
goalsmetbymod2[, modeltype := 'planning']
goalsmetbymod2[model %in% gcminds, modeltype := 'testing']


# write out
write.csv(goalsbyzonebymod1, file = gzfile(paste0('output/goalsbyzonebymod_', runname1out, '.csv.gz')))
write.csv(goalsbyzonebymod2, file = gzfile(paste0('output/goalsbyzonebymod_', runname2out, '.csv.gz')))
write.csv(goalsmetbymod1, file = paste0('output/goalsmetbymod_', runname1out, '.csv'))
write.csv(goalsmetbymod2, file = paste0('output/goalsmetbymod_', runname2out, '.csv'))


# read back in if needed
goalsbyzonebymod1 <- fread(paste0('gunzip -c output/goalsbyzonebymod_', runname1out, '.csv.gz'), drop = 1)
goalsbyzonebymod2 <- fread(paste0('gunzip -c output/goalsbyzonebymod_', runname2out, '.csv.gz'), drop = 1)
goalsmetbymod1 <- fread(paste0('output/goalsmetbymod_', runname1out, '.csv'), drop = 1)
goalsmetbymod2 <- fread(paste0('output/goalsmetbymod_', runname2out, '.csv'), drop = 1)
