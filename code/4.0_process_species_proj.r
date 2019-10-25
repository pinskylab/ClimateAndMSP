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

# path to the pres-abs projection results from Morley et al. 2018. At 0.05 grid size.
PRESABSPATH <- '/local/shared/pinsky_lab/projections_PlosOne2018/CEmodels_proj_PresAbs_May2018' 


####################
## helper functions
####################
require(data.table)



############################################
## Summarize species occurrence
## by grid cell for each model separately
############################################


projcols <- c('lat', 'lon', paste0('mean', 1:NMODS)) # names of columns in pred.agg to use (the projections)
for (i in 1:length(RCPS)) {
    print(paste0(Sys.time(), ' On RCP ', RCPS[i], '. Species: '))
    files <- list.files(path = PRESABSPATH, pattern = paste('*rcp', RCPS[i], '*', sep = ''), full.names = TRUE)
    
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
            
            # aggregate p(occur) across projection grid cells into the larger CMSP analysis grid cells
            sppbygrid <- pred.agg[, .(spp = thisspp, rcp = RCPS[i], model = m, poccur = mean(get(paste0('mean', m)), na.rm = TRUE)), by = c('latgrid', 'longrid')]
            #print(dim(sppbyMPA))
            
            # add this species/model/rcp onto the output data.table
            if(!exists('wdpa_by_spp_projnow')) { # create the output table if it doesn't exist
                presmap <- sppbygrid
            } else {
                presmap <- rbind(presmap, sppbygrid)
            }
            
        }
    }
}
print(Sys.time())

dim(presmap) #


# write out species presence/absence data 
write.csv(presmap, file = gzfile('temp/wdpa_by_spp_projnow.csv.gz'))


