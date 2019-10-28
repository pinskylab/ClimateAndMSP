# Summarize species pres/abs and biomass on the CMSP analysis grip

############
## Flags
############

# choose the analysis grid in degrees lat/lon
gridsize <- 0.25

# choose the rcp(s)
RCPS <- c(26, 85)

# select initial and final timeperiod for these grids
PERIODS <- c('2007-2020', '2081-2100')

# number of climate models in the projections
NMODS <- 18
NMODSbio <- 16 # for the biomass models

# path to the pres-abs and biomass projection results from Morley et al. 2018. At 0.05 grid size.
PRESABSPATH <- '/local/shared/pinsky_lab/projections_PlosOne2018/CEmodels_proj_PresAbs_May2018' 
BIOMASSPATH <- '/local/shared/pinsky_lab/projections_PlosOne2018/CEmodels_proj_Biomass_BCODMO' 


####################
## helper functions
####################
require(data.table)



############################################
## Summarize species occurrence
## by grid cell for each model separately
############################################

# step through each RCP and Ocean separately
for (i in 1:length(RCPS)) {
    for (ocean in c('Atl', 'Pac')){
        if (exists('presmap')) rm(presmap) # will store results in presmap

        files <- list.files(path = PRESABSPATH, pattern = paste0('*_', ocean, '_rcp', RCPS[i], '*'), full.names = TRUE) # find relevant projection files
        print(paste0(Sys.time(), ' On RCP ', RCPS[i], ', Ocean ', ocean, '. Species (of ', length(files), '): ')) # print status

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
            pred.agg <- pred.agg[year_range %in% PERIODS,]
            
            # calculate summary stats about species p(occur) by grid cell for each climate model
            for (m in 1:NMODS) { # for each climate model
                
                # average p(occur) across projection grid cells into the larger CMSP analysis grid cells
                # use sum/25 to implicitly consider 0s for cells on land (there are 25 climate projection grid cells per CMSP analysis grid)
                sppbygrid <- pred.agg[, .(spp = thisspp, rcp = RCPS[i], model = m, poccur = sum(get(paste0('mean', m)), na.rm = TRUE)/(gridsize/0.05)^2), by = c('latgrid', 'longrid', 'year_range')]
                #print(dim(sppbygrid))
                
                # add this species/model/rcp/time onto the output data.table
                if (!exists('presmap')) { # create the output table if it doesn't exist
                    presmap <- sppbygrid
                } else {
                    presmap <- rbind(presmap, sppbygrid)
                }
                
            }
        }
        
        print(dim(presmap))

        # write out species presence/absence data 
        outfile <- paste0('temp/presmap_', ocean, '_rcp', RCPS[i], '.csv.gz')
        write.csv(presmap, file = gzfile(outfile))
        print(paste('wrote', outfile))
    }
    
}

print(Sys.time())




############################################
## Summarize species biomass
## by grid cell for each model separately
############################################
fisheryspp <- fread('output/fishery_spps.csv', drop = 1)


for (i in 1:length(RCPS)) {
    for (ocean in c('Atl', 'Pac')) {
        if (exists('biomassmap')) rm(biomassmap) # the output table. start fresh for each rcp/ocean

        files <- list.files(path = BIOMASSPATH, pattern = paste0('*_', ocean, '_rcp', RCPS[i], '*'), full.names = TRUE) # find relevant projection files

        ## trim to fishery species
        spps <- gsub(paste0(BIOMASSPATH, '/|_Atl|_Pac|_rcp26|_rcp85|_jas_prediction_AGG.RData'), '', files)
        if (ocean == 'Atl') oceanspp <- fisheryspp[region %in% c('gmex', 'seus', 'neus', 'maritime', 'newf'), projname] # projected species in the current ocean
        if (ocean == 'Pac') oceanspp <- fisheryspp[region %in% c('wc', 'goa_bc', 'ebs'), projname]
        torun <- spps %in% oceanspp
        if (length(setdiff(oceanspp, spps[torun])) > 0) stop(paste('missing species for i =', i, 'ocean =', ocean)) # make sure all the fishery species are in the available projections
        files <- files[torun]
    
        print(paste0(Sys.time(), ' On RCP ', RCPS[i], ', Ocean ', ocean, '. Species (of ', length(files), '): ')) # print status

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
            pred.agg <- pred.agg[year_range %in% PERIODS,]
            
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
        outfile <- paste0('temp/biomassmap_', ocean, '_rcp', RCPS[i], '.csv.gz')
        write.csv(biomassmap, file = gzfile(outfile))
        print(paste('wrote', outfile))

    }
}

print(Sys.time())

