# Summarize species pres/abs and biomass on the CMSP analysis grip

############
## Flags
############

# choose the analysis grid in degrees lat/lon
gridsize <- 0.25

# choose the rcp(s)
RCPS <- c(26, 85)


# number of climate models in the projections
NMODSpres <- 18 # for the pres/abs models. Models 12 and 16 are extraneous.
NMODSbio <- 16 # for the biomass models
MODSPRESTODROP <- c(12,16)

# path to the pres-abs and biomass projection results from Morley et al. 2018. At 0.05 grid size.
PRESABSPATH <- '/local/shared/pinsky_lab/projections_PlosOne2018/CEmodels_proj_PresAbs_May2018' 
BIOMASSPATH <- '/local/shared/pinsky_lab/projections_PlosOne2018/CEmodels_proj_Biomass_BCODMO' 


####################
## helper functions
####################
require(data.table)


###########################
## set up some parameters
###########################

# define the time-periods from the first file
files <- list.files(path = PRESABSPATH, pattern = paste0('*_Atl_rcp26*'), full.names = TRUE) # find relevant projection files
load(files[1]) # loads pred.agg
periods <- as.character(sort(unique(pred.agg$year_range)))
            

############################################
## Summarize species occurrence
## by grid cell for each model separately
## Write files for each rcp/ocean/time period
############################################

# step through each RCP, Ocean, and time period separately
for (i in 1:length(RCPS)) {
    for (ocean in c('Atl', 'Pac')){
        for (k in 1:length(periods)){
            if (exists('presmap')) rm(presmap) # will store results in presmap
    
            files <- list.files(path = PRESABSPATH, pattern = paste0('*_', ocean, '_rcp', RCPS[i], '*'), full.names = TRUE) # find relevant projection files
            print(paste0(Sys.time(), ' On RCP ', RCPS[i], ', Ocean ', ocean, ', Period ', periods[k], '. Species (', length(files), '): ')) # print status
    
            # load presmap for each model run and process the results
            for (j in 1:length(files)) {
                cat(paste0(' ', j))
                load(files[j]) # loads pred.agg data.frame
                thisspp <- gsub(paste0(PRESABSPATH, '/|_Atl|_Pac|_rcp26|_rcp85|_jas_PresAbs_projection_AGG.RData'), '', files[j])
                pred.agg <- as.data.table(pred.agg)
                
                # round lat lon to nearest CMSP grid cell
                pred.agg[, latgrid := floor(latitude/gridsize)*gridsize + gridsize/2]
                pred.agg[, longrid := floor(longitude/gridsize)*gridsize + gridsize/2]
                setkey(pred.agg, latgrid, longrid)
                
                # trim to focal year ranges
                pred.agg <- pred.agg[year_range %in% periods[k],]
                
                # calculate summary stats about species p(occur) by grid cell for each climate model
                for (m in 1:NMODSpres) { # for each climate model
                    if(!(m %in% MODSPRESTODROP)){ # drop two of the pres-abs models
                        if (m < 12) thismodnum <- m
                        if (m > 12 & m < 16) thismodnum <- m - 1 # adjust model numbers to match biomass models
                        if (m > 16) thismodnum <- m - 2
                        
                        # average p(occur) across projection grid cells into the larger CMSP analysis grid cells
                        # use sum/25 to implicitly consider 0s for cells on land (there are 25 climate projection grid cells per CMSP analysis grid)
                        sppbygrid <- pred.agg[, .(spp = thisspp, rcp = RCPS[i], model = thismodnum, poccur = sum(get(paste0('mean', m)), na.rm = TRUE)/(gridsize/0.05)^2), by = c('latgrid', 'longrid', 'year_range')]
                        #print(dim(sppbygrid))
                        
                        # add this species/model/rcp/time onto the output data.table
                        if (!exists('presmap')) { # create the output table if it doesn't exist
                            presmap <- sppbygrid
                        } else {
                            presmap <- rbind(presmap, sppbygrid)
                        }
                    }
                }
            }
            
            print(dim(presmap))
    
            # write out species presence/absence data 
            outfile <- paste0('temp/presmap_', ocean, '_rcp', RCPS[i], '_', periods[k], '.csv.gz')
            write.csv(presmap, file = gzfile(outfile))
            print(paste('wrote', outfile))
        }
    }
}

print(Sys.time())




############################################
## Summarize species biomass
## by grid cell for each model separately
## Write files for each rcp/ocean/time period
############################################
fisheryspp <- fread('output/fishery_spps.csv', drop = 1)

# step through each RCP, Ocean, and time period separately
for (i in 1:length(RCPS)) {
    for (ocean in c('Atl', 'Pac')) {
        for (k in 1:length(periods)){
            if (exists('biomassmap')) rm(biomassmap) # the output table. start fresh for each rcp/ocean
    
            files <- list.files(path = BIOMASSPATH, pattern = paste0('*_', ocean, '_rcp', RCPS[i], '*'), full.names = TRUE) # find relevant projection files
    
            ## trim to fishery species
            spps <- gsub(paste0(BIOMASSPATH, '/|_Atl|_Pac|_rcp26|_rcp85|_jas_prediction_AGG.RData'), '', files)
            if (ocean == 'Atl') oceanspp <- fisheryspp[region %in% c('gmex', 'seus', 'neus', 'maritime', 'newf'), projname] # projected species in the current ocean
            if (ocean == 'Pac') oceanspp <- fisheryspp[region %in% c('wc', 'goa', 'bc', 'ebs'), projname]
            torun <- spps %in% oceanspp
            if (length(setdiff(oceanspp, spps[torun])) > 0) stop(paste('missing species for i =', i, 'ocean =', ocean)) # make sure all the fishery species are in the available projections
            files <- files[torun]
        
            print(paste0(Sys.time(), ' On RCP ', RCPS[i], ', Ocean ', ocean, ', Period ', periods[k], '. Species (', length(files), '): ')) # print status

            # load biomass for each model run and process the results
            for (j in 1:length(files)) {
                cat(paste0(' ', j))
                load(files[j]) # loads pred.agg data.frame
                thisspp <- gsub(paste0(BIOMASSPATH, '/|_Atl|_Pac|_rcp26|_rcp85|_jas_prediction_AGG.RData'), '', files[j])
                pred.agg <- as.data.table(pred.agg)
                
                # round lat lon to nearest CMSP grid cell
                pred.agg[, latgrid := floor(latitude/gridsize)*gridsize + gridsize/2]
                pred.agg[, longrid := floor(longitude/gridsize)*gridsize + gridsize/2]
                setkey(pred.agg, latgrid, longrid)
                
                # trim to focal year ranges
                pred.agg <- pred.agg[year_range %in% periods[k],]
                
                # calculate summary stats about species p(occur) by grid cell for each climate model
                for (m in 1:NMODSbio) { # for each climate model
                    
                    # aggregate biomass across projection grid cells into the larger CMSP analysis grid cells
                    # use sum/25 to implicitly consider 0s for cells on land (there are 25 climate projection grid cells per CMSP analysis grid)
                    sppbygrid <- pred.agg[, .(spp = thisspp, rcp = RCPS[i], model = m, biomass = sum(get(paste0('mean', m)), na.rm = TRUE)/(gridsize/0.05)^2), by = c('latgrid', 'longrid', 'year_range')]
                    #print(dim(sppbygrid))
                    
                    # add this species/model/rcp/time onto the output data.table
                    if (!exists('biomassmap')) { # create the output table if it doesn't exist
                        biomassmap <- sppbygrid
                    } else {
                        biomassmap <- rbind(biomassmap, sppbygrid)
                    }
                    
                }
            }
    
            print(dim(biomassmap))
    
            # write out species presence/absence data 
            outfile <- paste0('temp/biomassmap_', ocean, '_rcp', RCPS[i], '_', periods[k], '.csv.gz')
            write.csv(biomassmap, file = gzfile(outfile))
            print(paste('wrote', outfile))
        }
    }
}

print(Sys.time())

