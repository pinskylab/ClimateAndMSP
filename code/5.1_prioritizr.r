# Set up and run Prioritizr with zones to simulate CMSP in one region
# set up to source from within R 3.5.3: source('code/5.1_prioritizr.r')
# May need to set R_MAX_VSIZE=60000000 or larger in .Renviron to avoid hitting memory limits (Sys.getenv('R_MAX_VSIZE') to query)


#############
## Parameters
#############

# choose the rcp (for runs using just one)
rcps <- c(26, 85)

# choose the climate models to use for future planning (save others for testing)
#bcc-csm1-1-m, bcc-csm1-1, CanESM2, CCSM4, CESM1-CAM5, CNRM-CM5, GFDL-CM3, GFDL-ESM2M, GFDL-ESM2G, GISS-E2-R, GISS-E2-H, IPSL-CM5A-LR, IPSL-CM5A-MR, MIROC-ESM, MPI-ESM-LR, NorESM1-ME
gcminds <- c(1, 2, 3, 4, 8, 9, 10, 14) # from running sample(1:16, 8)

# CMSP goals
consgoal <- 0.1 # proportion of presences to capture in conservation
energygoal <- 0.2 # proportion of NPV
fishgoal <- 0.5 # proportion of biomass
cost <- 0.01 # basic cost of including each planning unit in a given zone

# poccur threshold: how high does the probability of occurrence in the projections need to be to consider the species "present"?
poccurthresh <- 0.5

# oceans to read in
oceans <- c('Atl', 'Pac')

# choose region and name these runs
myregs <- c('ebs', 'goa', 'bc', 'wc', 'gmex', 'seus', 'neus', 'maritime', 'newf')
runname1s <- paste0('hist_', myregs)
runname2s <- paste0('2per_', myregs)

# which time periods to use in the multi-period planning
planningperiods <- c('2007-2020', '2081-2100')

# folders
prioritizrfolder <- 'output/prioritizr_runs/'

# optimality gap for gurobi solver
gap <- 0.01

######################
# Functions
######################
require(data.table)
library(prioritizr) # only runs in R 3.5.3 for now (Gurobi 8.1.1)


#####################
## Load data
#####################

# loads presence/absence and biomass data
if(!(length(rcps) %in% c(1,2))){
	stop('rcp must be length 1 or 2')
}
for (i in 1:length(rcps)){
    print(paste0('rcp', rcps[i]))
    
    # code if files are saved w/out timeperiod in the name
    for(j in 1:length(oceans)){
        print(paste0('ocean', oceans[j]))
        prestemp <- fread(cmd = paste0('gunzip -c temp/presmap_', oceans[j], '_rcp', rcps[i], '.csv.gz'), drop = 1)
        prestemp <- prestemp[model %in% c(1:11, 13:15, 17:18)[gcminds], ] # trim to focal climate models. have to account for extra models in pres/abs projections (18 instead of 16 for biomass)
        
        biotemp <- fread(cmd = paste0('gunzip -c temp/biomassmap_', oceans[j], '_rcp', rcps[i], '.csv.gz'), drop = 1)
        biotemp <- biotemp[model %in% c(1:16)[gcminds], ]
        
        if(i == 1 & j == 1){
            presmap <- prestemp
            biomassmap <- biotemp
        } else {
            presmap <- rbind(presmap, prestemp)
            biomassmap <- rbind(biomassmap, biotemp)
        }
    }
    
    # code if files are saved by time period
    # for(k in 1:length(periods)){
    #     prestemp <- fread(cmd = paste0('gunzip -c temp/presmap_Atl_rcp', rcps[i], '_', periods[k], '.csv.gz'), drop = 1)
    #     prestemp <- prestemp[model %in% c(1:11, 13:15, 17:18)[gcminds], ] # trim to focal climate models. have to account for extra models in pres/abs projections (18 instead of 16 for biomass)
    #     
    #     biotemp <- fread(cmd = paste0('gunzip -c temp/biomassmap_Atl_rcp', rcps[i], '_', periods[k], '.csv.gz'), drop = 1)
    #     biotemp <- biotemp[model %in% c(1:16)[gcminds], ]
    # 
    #     if(i == 1 & j == 1 & k == 1){
    #         presmap <- prestemp
    #         biomassmap <- biotemp
    #     } else {
    #         presmap <- rbind(presmap, prestemp)
    #         biomassmap <- rbind(biomassmap, biotemp)
    #     }
    # 
    #     rm(prestemp, biotemp)
    # }
}
rm(prestemp, biotemp)

# load NatCap calculations
windnpv <- fread(cmd = 'gunzip -c output/wind_npv.csv.gz', drop = 1)
wavenpv <- fread(cmd = 'gunzip -c output/wave_npv.csv.gz', drop = 1)

setnames(windnpv, c('lat', 'lon', 'npv'), c('latgrid', 'longrid', 'wind_npv'))
setnames(wavenpv, c('lat', 'lon', 'npv'), c('latgrid', 'longrid', 'wave_npv'))

# definition of fishery species by region
fisheryspps <- fread('output/fishery_spps.csv', drop = 1) # which spp to include in fishery goal in each region

# region definitions
regiongrid <- fread(cmd = 'gunzip -c output/region_grid.csv.gz', drop = 1)



################################
## Set up data for any region
################################

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

# Fix a species name
presmap[spp == 'theragra chalcogramma', spp := 'gadus chalcogrammus']

# zones
# id and names for each zone
zones <- data.frame(id = 1:3, name = c('conservation', 'fishery', 'energy'))

############################
# Run prioritizr
# Do hist-only and 2-period
#############################

for (i in 1:length(myregs)) {
    print(paste0('Starting region ', myregs[i]))
        
    ###############################
    # Set up data for this region
    ###############################
    
    # pus
    # planning features are each 1/4 deg square
    pus <- presmap[region == myregs[i], c('latgrid', 'longrid')]
    pus <- pus[!duplicated(pus),]
    	dim(pus) # 2195 (ebs), 795 (goa), (bc), (wc), 651 (gomex), (seus), (neus), (maritime), (newf)
    if(nrow(pus) == 0) stop('pus has length zero')
    pus <- pus[order(pus$latgrid, pus$longrid),]
    pus$id <- 1:nrow(pus)
    pus$dummycost <- rep(cost, nrow(pus)) # set the same cost in each planning unit. can add separate costs for each zone.
    
    
    ############################################
    ## Run prioritizr just on 2007-2020
    ############################################
    
    # spps
    # id and name for each species
    # fishery features entered separately from conservation features, even if same species
    	sppstokeep <- presmap[region == myregs[i] & year_range == planningperiods[1] & rcp %in% rcps & model %in% c(1:11, 13:15, 17:18), .(poccur = mean(poccur)), by = c('latgrid', 'longrid', 'spp')] # average across models
    		dim(sppstokeep)
    	sppstokeep <- sppstokeep[poccur >= poccurthresh, ]
    	sppstokeep <- merge(sppstokeep, pus[, .(latgrid, longrid, id)], by = c('latgrid', 'longrid')) # add pu id (and trim to focal pus)
    		setnames(sppstokeep, 'id', 'pu')
    		dim(sppstokeep)
    
    		ngrid <- sppstokeep[ , .(ngrid = length(unique(pu))), by = 'spp']
    		sppstokeep <- merge(sppstokeep, ngrid, by = 'spp')
    		sppstokeep[ , summary(ngrid)] #
    
    		nspps <- sppstokeep[ , .(nspp = length(unique(spp))), by = 'pu']
    		sppstokeep <- merge(sppstokeep, nspps, by = 'pu')
    		sppstokeep[, summary(nspp)] #
    
    		# sppstokeep <- sppstokeep[ngrid > (nrow(pus)*0.05),] # trim to species found at poccur > poccurthresh in at least 5% of grids
    
    		sppstokeep[ , length(unique(spp))] # 179 (ebs w/ 0.1 poccurthresh, 48 w/ 0.5), 161 (goa), (bc),  (wc), 195 (gmex), (seus), (neus), (maritime), (newf)
    
    	spps <- data.table(id = 1:length(unique(sppstokeep$spp)), name = gsub(' |_', '', sort(unique(sppstokeep$spp))), spp = sort(unique(sppstokeep$spp))) #  fill spaces in species names.
    
    	# add fishery features
    	spps <- rbind(spps, data.table(id = max(spps$id) + 1:fisheryspps[region == myregs[i], length(projname)], 
    	                               name = paste0(gsub(' |_', '', fisheryspps[region == myregs[i], projname]), '_fishery'),
    	                               spp = fisheryspps[region == myregs[i], projname]))
    	
    	# add wind and wave energy feature
    	spps <- rbind(spps, data.table(id = max(spps$id) + 1, name = c('energy'), spp = c(NA)))
    
    
    # puvsp
    # which features are in each planning unit
    	# Format conservation data
    	puvsppa <- presmap[region == myregs[i] & year_range == planningperiods[1] & rcp %in% rcps & model %in% c(1:11, 13:15, 17:18), .(poccur = mean(poccur)), by = c('latgrid', 'longrid', 'spp')] # pres/abs data. climate models are bcc-csm1-1-m, bcc-csm1-1, CanESM2, CCSM4, CESM1-CAM5, CNRM-CM5, GFDL-CM3, GFDL-ESM2M, GFDL-ESM2G, GISS-E2-R, GISS-E2-H, HadGEM2-ES, IPSL-CM5A-LR, IPSL-CM5A-MR, MIROC-ESM, MIROC5, MPI-ESM-LR, NorESM1-ME. only use those that match biomass projections
    		dim(puvsppa)
    	puvsppa[, amount := as.numeric(poccur >= poccurthresh)] # use pres/abs as conservation amount. should this instead be left as poccur?
    	    puvsppa[, summary(amount)]
    	    puvsppa[, sort(unique(amount))]
    	puvsppa[, poccur := NULL]
    	puvsppa[ , name := gsub(' |_', '', spp)] # trim out spaces from species names
    	
        # Format fishery data
        puvspbio <- biomassmap[region == myregs[i] & year_range == planningperiods[1] & rcp %in% rcps & model %in% 1:16 & spp %in% fisheryspps[region == myregs[i], projname], .(biomass = mean(biomass)), by = c('latgrid', 'longrid', 'spp')] # biomass data. models are bcc-csm1-1-m, bcc-csm1-1, CanESM2, CCSM4, CESM1-CAM5, CNRM-CM5, GFDL-CM3, GFDL-ESM2M, GFDL-ESM2G, GISS-E2-R, GISS-E2-H, IPSL-CM5A-LR, IPSL-CM5A-MR, MIROC-ESM, MPI-ESM-LR, NorESM1-ME
            dim(puvspbio)
            puvspbio[, length(unique(spp))] # should be 10
        puvspbio[, amount := biomass] # use biomass as amount for fishery targets
        puvspbio[ , name := paste0(gsub(' |_', '', spp), '_fishery')] # trim out spaces from species names
     
    	# Format wind and wave data
    	puvenergy <- merge(windnpv, wavenpv, by = c('latgrid', 'longrid'), all = TRUE)
    		head(puvenergy)
    		dim(windnpv)
    		dim(wavenpv)
    		dim(puvenergy)
    	puvenergy[wind_npv < 0 | is.na(wind_npv), wind_npv := 0] # set negative or NA NPV to 0
    	puvenergy[wave_npv < 0 | is.na(wave_npv), wave_npv := 0]
    	puvenergy[, amount := (wind_npv + wave_npv)/10000] # scale down to pass presolve checks
    	puvenergy[, name := 'energy']
    	
    	# combine
    	puvsp <- rbind(puvsppa[, .(name, latgrid, longrid, amount, zone = 1)], 
    	               puvspbio[, .(name, latgrid, longrid, amount, zone = 2)], 
    	               puvenergy[, .(name, latgrid, longrid, amount, zone = 3)])
    	
    	# Add species ids
    	nrow(puvsp)
    	puvsp <- merge(puvsp, spps[, .(id, name)], by = 'name') # merge in species IDs and trim to focal species
    	nrow(puvsp)
    	setnames(puvsp, 'id', 'species')
    		
        # Add planning units
        puvsp <- merge(puvsp, pus[, .(latgrid, longrid, id)], by = c('latgrid', 'longrid')) # add pu id (and trim to focal pus)
    	nrow(puvsp)
    	setnames(puvsp, 'id', 'pu')
    	
    	# Check fishery species for adequate biomass and scale up if needed
    	# Makes sure that no fishery species are eliminated by the next section checking for amount < 1e6
        fishtotals <- puvsp[grepl('fishery', name), .(total = sum(amount), name = unique(name)), by = 'species']
    	for(j in which(fishtotals[, total < 1])){
    	    scalar <- 1/fishtotals[j, total] # scale up so sum would be 1
    	    puvsp[species == fishtotals[j, species], amount := amount * scalar]
    	}
        
    	# Trim out values < 1e-6 (will throw error in prioritizr)
    	puvsp[amount < 1e-6, amount := 0]

    	# Sort and trim columns and rows
    	setkey(puvsp, pu, species) # order by pu then species
    	puvsp <- puvsp[amount > 0, ] # trim only to presences
    	
    	
        # checks
    	if(length(unique(puvsp$pu)) != nrow(pus)) stop(paste0('region: ', myregs[i], '. puvsp planning units do not match pus.')) # planning units for species + NatCap: 2195 (ebs), 661 (goa), 549 (neus), 1342 (newf)
    	if(!all(unique(puvsp$species) %in% spps$id)) stop(paste0('region: ', myregs[i], '. Some puvsp features are not in spps.')) # features that are species + fishery + NatCap
    	if(min(sort(unique(table(puvsp$species)))) < 1) stop(paste0('region: ', myregs[i], '. Some species are not in a planning unit (hist).')) # make sure all species show up in some planning units (shouldn't see any 0s)
    	if(min(sort(unique(table(puvsp$pu))) < 1)) stop(paste0('region: ', myregs[i], '. Some planning units do not have a species (hist).')) # make sure all planning units have some species (shouldn't see any 0s)
    	if(!all(sort(unique(table(puvsp$pu, puvsp$species))) %in% c(0,1))) stop(paste0('region: ', myregs[i], '. Some planning unit-species combinations appear more than once (hist).')) # should be all 0s and 1s
        if(puvsp[, max(amount) > 1e6]) stop(paste0('region:', myregs[i], '. Amount > 1e6 (hist).'))
    
    #zone target
    # set zone-specific targets: rows are features, columns are zones
	zonetarget <- matrix(0, nrow = nrow(spps), ncol = nrow(zones), dimnames = list(spps$name, zones$name))
	zonetarget[!grepl('energy|fishery', rownames(zonetarget)), 'conservation'] <- consgoal # set conservation zone target
    zonetarget[grepl('fishery', rownames(zonetarget)), 'fishery'] <- fishgoal # set fishing zone target
	zonetarget[grepl('energy', rownames(zonetarget)), 'energy'] <- energygoal # set energy goal target
    	
    # basic checks (automated)
    if(!all(colSums(zonetarget) > 0)) stop(paste0('region:', myregs[i], '. Some zone targets are 0 (hist).')) # reasonable targets?
    if(nrow(zonetarget) != nrow(spps)) stop(paste0('region: ', myregs[i], '. Zonetargets do not match spps (hist).'))
    if(!all(rownames(zonetarget) == spps$name))  stop(paste0('region: ', myregs[i], '. Zonetargets order does not match spps order (hist).'))
    if(sum(!(puvsp$pu %in% pus$id)) > 0) stop(paste0('region: ', myregs[i], '. Some planning units not in pus (hist).'))
    if(sum(!(puvsp$species %in% spps$id)) > 0) stop(paste0('region: ', myregs[i], '. Some species units not in spps (hist).'))
    if(sum(!(pus$id %in% puvsp$pu)) > 0) stop(paste0('region: ', myregs[i], '. Some pus units not in puvsp (hist).'))
    if(sum(!(spps$id %in% puvsp$species)) > 0) stop(paste0('region: ', myregs[i], '. Some species units not in puvsp (hist).'))
    	
    
    # Define the problem in prioritzr format
    # rij may need a zone column
    p1 <- problem(pus, spps, cost_column = c('dummycost', 'dummycost', 'dummycost'), rij = puvsp, zones = zones) %>%
        add_min_set_objective() %>%
        add_relative_targets(zonetarget) %>%
        add_binary_decisions() %>%
        add_gurobi_solver(gap = gap)
    
    # solve it
    if(presolve_check(p1)){
        s1 <- solve(p1)
    } else {
        stop(paste0('region:', myregs[i], '. Failed presolve check (hist).'))
    }
    
    # force solution (e.g., if fails presolve checks)
    # but beware: solution may be meaningless
    # s1 <- p1 %>%
    #     add_gurobi_solver(numeric_focus = TRUE) %>%
    #     solve(force = TRUE)
    
    # examine solution in various ways
    s1[,.(navail = sum(solution_1_conservation == 0 & solution_1_fishery == 0 & solution_1_energy == 0), 
          ncons = sum(solution_1_conservation), 
          nfish = sum(solution_1_fishery), 
          nenergy = sum(solution_1_energy))]  # number of cells in each zone
    r1 <- feature_representation(p1, s1[, .(solution_1_conservation, solution_1_fishery, solution_1_energy)]) # representation of each feature in each zone
    r1dt <- as.data.table(r1)
    r1dt[zone=='conservation',] # meeting conservation targets?
    r1dt[zone=='conservation' & !grepl('fishery|energy', feature), .(summary(absolute_held), summary(relative_held))] # meeting conservation targets?
    r1dt[zone=='fishery' & grepl('fishery', feature),] # fishery targets
    r1dt[zone=='energy',]
    
    
    # write out
    write.csv(s1, file = paste0(prioritizrfolder, 'solution_', runname1s[i], '.csv'))
    
    
    ##################################################################
    ## Set up a prioritizr run on 2006-2020 and half an ensemble mean 2081-2100
    ## This assumes that the historical-only code (previous section) has been run and is loaded in memory
    ##################################################################
    
    # spps2
    # document every species present
    # add future species to include
    	spps2 <- spps #  use the same species as in the historical-only run
    
    	sppinds <- !grepl('energy', spps2$name) # don't include energy in each time period
    	temp1 <- spps2[sppinds,]
    	spps2$name[sppinds] <- paste0(spps2$name[sppinds], gsub('-', '', planningperiods[1]))
    	temp1$name <- paste0(temp1$name, gsub('-', '', planningperiods[2]))
    	temp1$id = temp1$id + max(spps2$id) # make sure the ids don't overlap
    	spps2 <- rbind(spps2, temp1)
    
    # puvsp2
    # table of species by planning units
    # add future species to include
    # use only first 8 models, leaving later 8 for testing
    	# Format future conservation data
    	puvsppa2 <- presmap[region == myregs[i] & year_range == planningperiods[2] & rcp %in% rcps & model %in% c(1:11, 13:15, 17:18)[gcminds], .(poccur = mean(poccur)), by = c('latgrid', 'longrid', 'spp')] # pres/abs data. climate models are bcc-csm1-1-m, bcc-csm1-1, CanESM2, CCSM4, CESM1-CAM5, CNRM-CM5, GFDL-CM3, GFDL-ESM2M, GFDL-ESM2G, GISS-E2-R, GISS-E2-H, HadGEM2-ES, IPSL-CM5A-LR, IPSL-CM5A-MR, MIROC-ESM, MIROC5, MPI-ESM-LR, NorESM1-ME. only use those that match biomass projections
    		dim(puvsppa2)
    	puvsppa2[, amount := as.numeric(poccur >= poccurthresh)] # use pres/abs as conservation amount. should this instead be left as poccur?
    	    puvsppa2[, summary(amount)]
    	    puvsppa2[, sort(unique(amount))]
    	puvsppa2[ , name := gsub(' |_', '', spp)] # trim out spaces from species names and add future
    	puvsppa2[!grepl('energy', name), name := paste0(name, gsub('-', '', planningperiods[2]))] # append time period
    	
        # Format fishery data
        puvspbio2 <- biomassmap[region == myregs[i] & year_range == planningperiods[2] & rcp %in% rcps & model %in% (1:16)[gcminds] & spp %in% fisheryspps[region == myregs[i], projname], .(biomass = mean(biomass)), by = c('latgrid', 'longrid', 'spp')] # biomass data. models are bcc-csm1-1-m, bcc-csm1-1, CanESM2, CCSM4, CESM1-CAM5, CNRM-CM5, GFDL-CM3, GFDL-ESM2M, GFDL-ESM2G, GISS-E2-R, GISS-E2-H, IPSL-CM5A-LR, IPSL-CM5A-MR, MIROC-ESM, MPI-ESM-LR, NorESM1-ME
            dim(puvspbio2)
            puvspbio2[, length(unique(spp))] # should be 10
        puvspbio2[, amount := biomass/100] # use biomass as amount for fishery targets. Scale down to pass presolve checks.
        puvspbio2[ , name := paste0(gsub(' |_', '', spp), '_fishery')] # trim out spaces from species names
    	puvspbio2[!grepl('energy', name), name := paste0(name, gsub('-', '', planningperiods[2]))] # append time period
     	
    	# combine future data
    	puvsp2 <- rbind(puvsppa2[, .(name, latgrid, longrid, amount, zone = 1)], 
    	               puvspbio2[, .(name, latgrid, longrid, amount, zone = 2)])
    	
    	# Add species ids
    	nrow(puvsp2)
    	puvsp2 <- merge(puvsp2, spps2[, .(id, name)], by = 'name') # merge in species IDs and trim to focal species
    	nrow(puvsp2)
    	setnames(puvsp2, 'id', 'species')
    		
        # Add planning units
        puvsp2 <- merge(puvsp2, pus[, .(latgrid, longrid, id)], by = c('latgrid', 'longrid')) # add pu id (and trim to focal pus)
    	nrow(puvsp2)
    	setnames(puvsp2, 'id', 'pu')
    	
    	# Add historical data
    	temp1 <- puvsp # start from the same data as in the historical run
    	temp1[!grepl('energy', name), name := paste0(name, gsub('-', '', planningperiods[1]))] # append time period (except energy)
    	puvsp2 <- rbind(temp1, puvsp2)

    	# Trim out values < 1e-6 (will throw error in prioritizr)
    	puvsp2[amount < 1e-6, amount := 0]
    	
    	# Sort and trim columns and rows
    	setkey(puvsp2, pu, species) # order by pu then species
    	puvsp2 <- puvsp2[amount > 0, ] # trim only to presences
    	
    	# checks
    	if(length(unique(puvsp2$pu)) != nrow(pus)) stop(paste0('region: ', myregs[i], '. puvsp2 planning units do not match pus.')) # planning units for species + NatCap
    	if(!all(unique(puvsp2$species) %in% spps2$id)) stop(paste0('region: ', myregs[i], '. Some puvsp2 features are not in spps.')) # features that are species + fishery + NatCap
    	if(min(sort(unique(table(puvsp2$species)))) < 1) stop(paste0('region: ', myregs[i], '. Some species are not in a planning unit (2per).')) # make sure all species show up in some planning units (shouldn't see any 0s)
    	if(min(sort(unique(table(puvsp2$pu))) < 1)) stop(paste0('region: ', myregs[i], '. Some planning units do not have a species (2per.')) # make sure all planning units have some species (shouldn't see any 0s)
    	if(!all(sort(unique(table(puvsp2$pu, puvsp2$species))) %in% c(0,1))) stop(paste0('region: ', myregs[i], '. Some planning unit-species combinations appear more than once (2per).')) # should be all 0s and 1s
    	if(puvsp2[, max(amount) > 1e6]) stop(paste0('region:', myregs[i], '. Amount > 1e6 (2per).'))
    
    #zone target
    # set zone-specific targets: rows are features, columns are zones
	zonetarget2 <- matrix(0, nrow = nrow(spps2), ncol = nrow(zones), dimnames = list(spps2$name, zones$name))
	zonetarget2[!grepl('energy|fishery', rownames(zonetarget)), 'conservation'] <- consgoal # set conservation zone target
	zonetarget2[grepl('fishery', rownames(zonetarget)), 'fishery'] <- fishgoal # set fishing zone target
	zonetarget2[grepl('energy', rownames(zonetarget)), 'energy'] <- energygoal # set energy goal target
    
    # trim out species that aren't present
	nrow(spps2)
	spps2 <- spps2[name %in% puvsp2$name,]
	nrow(spps2)

	nrow(zonetarget2)
	zonetarget2 <- zonetarget2[rownames(zonetarget2) %in% puvsp2$name,]
	nrow(zonetarget2)
    
	# basic checks (automated)
	if(!all(colSums(zonetarget) > 0)) stop(paste0('region:', myregs[i], '. Some zone targets are 0 (2per).')) # reasonable targets?
	if(nrow(zonetarget) != nrow(spps)) stop(paste0('region: ', myregs[i], '. Zonetargets do not match spps (2per).'))
	if(!all(rownames(zonetarget) == spps$name))  stop(paste0('region: ', myregs[i], '. Zonetargets order does not match spps order (2per).'))
	if(sum(!(puvsp$pu %in% pus$id)) > 0) stop(paste0('region: ', myregs[i], '. Some planning units not in pus (2per).'))
	if(sum(!(puvsp$species %in% spps$id)) > 0) stop(paste0('region: ', myregs[i], '. Some species units not in spps (2per).'))
	if(sum(!(pus$id %in% puvsp$pu)) > 0) stop(paste0('region: ', myregs[i], '. Some pus units not in puvsp (2per).'))
	if(sum(!(spps$id %in% puvsp$species)) > 0) stop(paste0('region: ', myregs[i], '. Some species units not in puvsp (2per).'))
    	
    
    # Define the problem in prioritzr format
    # rij may need a zone column
    p2 <- problem(pus, spps2, cost_column = c('dummycost', 'dummycost', 'dummycost'), 
                  rij = puvsp2, zones = zones) %>%
        add_min_set_objective() %>%
        add_relative_targets(zonetarget2) %>%
        add_binary_decisions() %>%
        add_gurobi_solver(gap = gap)
    
    # solve it
    if(presolve_check(p1)){
        s2 <- solve(p2)
    } else {
        stop(paste0('region:', myregs[i], '. Failed presolve check (2per).'))
    }
    
    # force solution (e.g., if fails presolve checks)
    # but beware: solution may be meaningless
    # s2 <- p2 %>%
    #     add_gurobi_solver(numeric_focus = TRUE) %>%
    #     solve(force = TRUE)
    
    # examine solution in various ways
    s2[,.(navail = sum(solution_1_conservation == 0 & solution_1_fishery == 0 & solution_1_energy == 0), 
          ncons = sum(solution_1_conservation), 
          nfish = sum(solution_1_fishery), 
          nenergy = sum(solution_1_energy))]  # number of cells in each zone
    r2 <- feature_representation(p2, s2[, .(solution_1_conservation, solution_1_fishery, solution_1_energy)]) # representation of each feature in each zone
    r1dt <- as.data.table(r2)
    r1dt[zone=='conservation',] # meeting conservation targets?
    r1dt[zone=='conservation' & !grepl('fishery|energy', feature), .(summary(absolute_held), summary(relative_held))] # meeting conservation targets?
    r1dt[zone=='fishery' & grepl('fishery', feature),] # fishery targets
    r1dt[zone=='energy' & grepl('energy', feature),]
    
    
    # write out
    write.csv(s2, file = paste0(prioritizrfolder, 'solution_', runname2s[i], '.csv'))

}