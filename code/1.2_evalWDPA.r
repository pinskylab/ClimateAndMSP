# evaluate protected areas and networks against shifts in species distribution
# calculate species gains, losses, turnover, etc.
# also output table of species p(occur) by MPA in the current time-period (2007-2020)


############
## Flags
############

# choose the rcp(s)
RCPS <- c(26, 85)

# select initial and final timeperiod for these grids
PERIODS <- c('2007-2020', '2081-2100')

# climate models in the projections
#bcc-csm1-1-m, bcc-csm1-1, CanESM2, CCSM4, CESM1-CAM5, CNRM-CM5, GFDL-CM3, GFDL-ESM2M, GFDL-ESM2G, GISS-E2-R, GISS-E2-H, HadGEM2-ES, IPSL-CM5A-LR, IPSL-CM5A-MR, MIROC-ESM, MIROC5, MPI-ESM-LR, NorESM1-ME
NMODS <- 18 # for the pres/abs models. Models 12 (HadGEM2-ES) and 16 (MIROC5) are extraneous since not in the biomass projections used in rest of project
MODSPRESTODROP <- c(12,16)

# path to the pres-abs projection results from Morley et al. 2018. At 0.05 grid size.
PRESABSPATH <- '/local/shared/pinsky_lab/projections_PlosOne2018/CEmodels_proj_PresAbs_May2018' 


####################
## helper functions
####################
require(data.table)
# set default rounding (2 bytes) so that data.table will merge numeric values appropriately
# see https://rdrr.io/rforge/data.table/man/setNumericRounding.html
setNumericRounding(2)

#require(lme4) # for mixed-effects models
#require(car) # for testing ME models
#require(Hmisc)



###########################
## Load and set up WDPA data
###########################

wdpagrid <- fread('gunzip -c output/wdpa_cov_by_grid0.05.csv.gz', drop = 1) # shows which MPAs are in which grid cells. each line is a unique grid cell-MPA combination.
	
	# convert lon to -360 to 0 (instead of -180 to 180) to match the pres/abs projections and to plot more nicely
	wdpagrid[lon>0, lon := lon-360]

	# set up MPA networks
	wdpagrid[SUB_LOC == 'US-CA' & MANG_AUTH == 'State Fish and Wildlife',network := 'mlpa']
	wdpagrid[grepl('US-DE|US-FL|US-GA|US-MA|US-MD|US-ME|US-NC|US-NJ|US-NY|US-RI|US-SC|US-VA', SUB_LOC) & !grepl('US-N/A', SUB_LOC) & (lon > -100), network := 'eastcoast'] # by exluding US-N/A, this won't include federal water closures
	wdpagrid[(SUB_LOC %in% c('US-AK')), network := 'ak'] # by choosing by state, this won't include federal water closures
	
	# calculate mpa extent
	wdpagrid[, ':='(lat_min = min(lat), lat_max = max(lat), lon_min = min(lon), lon_max = max(lon)), by = WDPA_PID]

#	wdpagrid[network=='mlpa',plot(lon,lat)]
#	wdpagrid[network=='eastcoast',plot(lon,lat)]
#	wdpagrid[network=='ak',plot(lon,lat)]
#	require(maps); map(database='world2',add=TRUE)

	
##############################################################
# Load in pres/abs results for each RCP and each species
# and calculate summary statistics on species turnover by MPA
##############################################################
# the list of MPAs. Turnover data will be added as columns to this
wdpa <- wdpagrid[!duplicated(WDPA_PID), ] 
wdpa[,c('gridpolyID', 'lat', 'lon', 'area_wdpa', 'area_grid', 'prop_grid', 'prop_wdpa') := NULL] # remove gridcell-specific columns

# step through each species' projections
projcols <- c('lat', 'lon', paste0('mean', 1:NMODS)) # names of columns in pred.agg to use (the projections)
for (i in 1:length(RCPS)) {
    print(paste0(Sys.time(), ' On RCP ', RCPS[i], '. Species: '))
    files <- list.files(path = PRESABSPATH, pattern = paste('*rcp', RCPS[i], '*', sep = ''), full.names = TRUE)

    # load presmap for each model run and process the results
    # then calculate change in p(occur)
    for (j in 1:length(files)) {
        cat(paste0(' ', j))
        load(files[j]) # loads pred.agg data.frame
        pred.agg <- as.data.table(pred.agg)
        
        # round lat lon to nearest 0.05 so that merging is successful with WDPA
        pred.agg[, lat := floor(latitude*20)/20 + 0.025] # assumes 0.05 grid size
        pred.agg[, lon := floor(longitude*20)/20 + 0.025] # assumes 0.05 grid size
        setkey(pred.agg, lat, lon)
        
        # reorganize so that time 2 and time 1 are separate columns, rows are locations
        presbyloc <- merge(pred.agg[year_range == PERIODS[1], ..projcols], pred.agg[year_range == PERIODS[2], ..projcols], by = c('lat', 'lon'), suffixes = c('.t1', '.t2'))

        # merge pres/abs data with WDPA grid data and calculate summary stats about species change through time
        # summary stat column names are of form [sumstatname].[RCP].[modnumber]
        wdpagridspp <- merge(wdpagrid, presbyloc, all.x = TRUE, by = c('lat', 'lon'))
        for (m in 1:NMODS) { # for each climate model
            if(!(m %in% MODSPRESTODROP)){ # drop two of the pres-abs models
                if (m < 12) thismodnum <- m
                if (m > 12 & m < 16) thismodnum <- m - 1 # adjust model numbers to match biomass models
                if (m > 16) thismodnum <- m - 2

                wdpagridspp[is.na(wdpagridspp[[paste0('mean', m, '.t1')]]), (paste0('mean', m, '.t1')) := 0] # set NAs to 0 (species not present)
                wdpagridspp[is.na(wdpagridspp[[paste0('mean', m, '.t2')]]), (paste0('mean', m, '.t2')) := 0]

                # aggregate p(occur) across grid cells into full MPAs
                sppbyMPA <- wdpagridspp[, .(pinit = 1 - prod(1 - get(paste0('mean', m, '.t1'))*prop_wdpa), # p(occur whole MPA) at initial timepoint as 1 - prod(1-p(occur per grid))
                                     pfinal = 1 - prod(1 - get(paste0('mean', m, '.t2'))*prop_wdpa)
                                     ), by = WDPA_PID]
                #print(dim(sppbyMPA))
                
                # merge wpda with sppbyMPA (latter has results from the current species)
                wdpa <- merge(wdpa, sppbyMPA, by = 'WDPA_PID')
                
                # calculate summary stats by MPA and add this species onto results from previous species
                if (!(paste0('ninit.', RCPS[i], '.', thismodnum) %in% colnames(wdpa))) { # if output column in wdpa doesn't yet exist, create it
                    wdpa[, (paste0('ninit.', RCPS[i], '.', thismodnum)) := pinit]
                    wdpa[, (paste0('nfinal.', RCPS[i], '.', thismodnum)) := pfinal]
                    wdpa[, (paste0('nshared.', RCPS[i], '.', thismodnum)) := pinit*pfinal]
                    wdpa[pinit > pfinal, (paste0('nlost.', RCPS[i], '.', thismodnum)) := pinit - pfinal]
                    wdpa[pinit <= pfinal, (paste0('nlost.', RCPS[i], '.', thismodnum)) := 0]
                    wdpa[pinit < pfinal, (paste0('ngained.', RCPS[i], '.', thismodnum)) := pfinal - pinit]
                    wdpa[pinit >= pfinal, (paste0('ngained.', RCPS[i], '.', thismodnum)) := 0]
                } else {# if columns already exist, add the new species' values onto the existing values
                    wdpa[, (paste0('ninit.', RCPS[i], '.', thismodnum)) := get(paste0('ninit.', RCPS[i], '.', thismodnum)) + pinit]
                    wdpa[, (paste0('nfinal.', RCPS[i], '.', thismodnum)) := get(paste0('nfinal.', RCPS[i], '.', thismodnum)) + pfinal]
                    wdpa[, (paste0('nshared.', RCPS[i], '.', thismodnum)) := get(paste0('nshared.', RCPS[i], '.', thismodnum)) + pinit*pfinal]
                    wdpa[pinit > pfinal, (paste0('nlost.', RCPS[i], '.', thismodnum)) := get(paste0('nlost.', RCPS[i], '.', thismodnum)) + pinit - pfinal]
                    wdpa[pinit < pfinal, (paste0('ngained.', RCPS[i], '.', thismodnum)) := get(paste0('ngained.', RCPS[i], '.', thismodnum)) + pfinal - pinit]
                }
                wdpa[, c('pinit', 'pfinal') := NULL] # delete the species-specific columns now that we have the summaries we need
            }
        }
        # wdpa[,print(max(get(paste0('ninit.', RCPS[i], '.1'))))] # a test that cumulative sums are being calculated. Should increase through the loop
    }
}
print(Sys.time())

dim(wdpa) # 923 x 127

# how many rows are zero?
wdpa[, .(sum(ninit.26.1 == 0), sum(ninit.85.1 == 0))]
wdpa[, .(sum(ninit.26.1 > 0), sum(ninit.85.1 > 0))]
    # wdpa[ninit.26.1 > 0, plot(lon_max, lat_max, cex = 0.1, xlim=c(-190, 0))]
    # wdpa[ninit.26.1 == 0, points(lon_max, lat_max, cex = 0.5, col = 'red')]

# write out species data 
write.csv(wdpa, file = gzfile('temp/wdpaturnbyMPAbymod.csv.gz'))

# wdpa <- fread('gunzip -c temp/wdpaturnbyMPAbymod.csv.gz', drop = 1) # don't read the row numbers


################################################################################
# Load in pres/abs results for each RCP and each species
# and aggregate p(occur) for each species in each MPA for current time-period
# for comparison against observations
################################################################################

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
        
        # round lat lon to nearest 0.05 so that merging is successful with WDPA
        pred.agg[, lat := floor(latitude*20)/20 + 0.025] # assumes 0.05 grid size
        pred.agg[, lon := floor(longitude*20)/20 + 0.025] # assumes 0.05 grid size
        setkey(pred.agg, lat, lon)
        
        # trim to 2007-2020
        pred.agg <- pred.agg[year_range == '2007-2020',]

        # merge pres/abs data with WDPA grid data
        wdpagridspp <- merge(wdpagrid, pred.agg, by = c('lat', 'lon')) # don't include all MPAs, to drop MPAs for which this species wasn't projected (e.g., Atl vs Pac)
        
        # calculate summary stats about species p(occur) by MPA for each climate model
        for (m in 1:NMODS) { # for each climate model
            if(!(m %in% MODSPRESTODROP)){ # drop two of the pres-abs models
                if (m < 12) thismodnum <- m
                if (m > 12 & m < 16) thismodnum <- m - 1 # adjust model numbers to match biomass models
                if (m > 16) thismodnum <- m - 2

                # aggregate p(occur) across grid cells into full MPAs
                sppbyMPA <- wdpagridspp[, .(spp = thisspp, rcp = RCPS[i], model = thismodnum, poccur = 1 - prod(1 - get(paste0('mean', m))*prop_wdpa)), by = WDPA_PID] # p(occur whole MPA) at initial timepoint as 1 - prod(1-p(occur per grid))
                    #print(dim(sppbyMPA))
                
                # add this species/model/rcp onto the output data.table
                if(!exists('wdpa_by_spp_projnow')) { # create the output table if it doesn't exist
                    wdpa_by_spp_projnow <- sppbyMPA
                } else {
                    wdpa_by_spp_projnow <- rbind(wdpa_by_spp_projnow, sppbyMPA)
                }
            }            
        }
    }
}
print(Sys.time())

dim(wdpa_by_spp_projnow) #


# write out species data 
write.csv(wdpa_by_spp_projnow, file = gzfile('temp/wdpa_by_spp_projnow.csv.gz'))




##############################################
# Calculate turnover within select MPA networks
# For each model in the ensemble
##############################################
wdpanet <- wdpagrid[!is.na(network), .(lat = mean(lat), lon = mean(lon), area = sum(area_fullwdpa)), by = network ] # the list of MPA networks. Turnover data will be added as columns to this


projcols <- c('lat', 'lon', paste0('mean', 1:NMODS)) # names of columns in pred.agg to use (the projections)
for (i in 1:length(RCPS)) {
    print(paste0(Sys.time(), ' On RCP ', RCPS[i], '. Species: '))
    files <- list.files(path = PRESABSPATH, pattern = paste('*rcp', RCPS[i], '*', sep = ''), full.names = TRUE)

    # load presmap for each model run and process the results
    # then calculate change in p(occur)
    for (j in 1:length(files)) {
        cat(paste0(' ', j))
        load(files[j]) # loads pred.agg data.frame
        pred.agg <- as.data.table(pred.agg)
        
        # round lat lon to nearest 0.05 so that merging is successful with WDPA
        pred.agg[, lat := floor(latitude*20)/20 + 0.025] # assumes 0.05 grid size
        pred.agg[, lon := floor(longitude*20)/20 + 0.025] # assumes 0.05 grid size
        setkey(pred.agg, lat, lon)
 
        # reorganize so that time 2 and time 1 are separate columns, rows are locations
        presbyloc <- merge(pred.agg[year_range == PERIODS[1], ..projcols], pred.agg[year_range == PERIODS[2], ..projcols], by = c('lat', 'lon'), suffixes = c('.t1', '.t2'))

        # merge pres/abs data with WDPA grid data and calculate summary stats about species change through time
        # summary stat column names are of form [sumstatname].[RCP].[modnumber]
        wdpagridspp <- merge(wdpagrid, presbyloc, all.x = TRUE, by = c('lat', 'lon'))
        for (m in 1:NMODS) { # for each climate model
             if(!(m %in% MODSPRESTODROP)){ # drop two of the pres-abs models
                if (m < 12) thismodnum <- m
                if (m > 12 & m < 16) thismodnum <- m - 1 # adjust model numbers to match biomass models
                if (m > 16) thismodnum <- m - 2
                
                wdpagridspp[is.na(wdpagridspp[[paste0('mean', m, '.t1')]]), (paste0('mean', m, '.t1')) := 0] # set NAs to 0 (species not present)
                wdpagridspp[is.na(wdpagridspp[[paste0('mean', m, '.t2')]]), (paste0('mean', m, '.t2')) := 0]
    
                # aggregate p(occur) across grid cells into full MPA networks
                sppbynet <- wdpagridspp[!is.na(network), .(pinit = 1 - prod(1 - get(paste0('mean', m, '.t1'))*prop_wdpa), # p(occur whole MPA) at initial timepoint as 1 - prod(1-p(occur per grid))
                                     pfinal = 1 - prod(1 - get(paste0('mean', m, '.t2'))*prop_wdpa)
                                     ), by = network]
                #print(dim(sppbynet))
                
                # merge wpda with sppbynet (latter has results from the current species)
                wdpanet <- merge(wdpanet, sppbynet, by = 'network')
                
                # calculate summary stats and add this species onto results from previous species
                if (!(paste0('ninit.', RCPS[i], '.', thismodnum) %in% colnames(wdpanet))) { # if output column in wdpanet doesn't yet exist, create it
                    wdpanet[, (paste0('ninit.', RCPS[i], '.', thismodnum)) := pinit]
                    wdpanet[, (paste0('nfinal.', RCPS[i], '.', thismodnum)) := pfinal]
                    wdpanet[, (paste0('nshared.', RCPS[i], '.', thismodnum)) := pinit*pfinal]
                    wdpanet[pinit > pfinal, (paste0('nlost.', RCPS[i], '.', thismodnum)) := pinit - pfinal]
                    wdpanet[pinit <= pfinal, (paste0('nlost.', RCPS[i], '.', thismodnum)) := 0]
                    wdpanet[pinit < pfinal, (paste0('ngained.', RCPS[i], '.', thismodnum)) := pfinal - pinit]
                    wdpanet[pinit >= pfinal, (paste0('ngained.', RCPS[i], '.', thismodnum)) := 0]
                } else {# if columns already exist, add the new species' values onto the existing values
                    wdpanet[, (paste0('ninit.', RCPS[i], '.', thismodnum)) := get(paste0('ninit.', RCPS[i], '.', thismodnum)) + pinit]
                    wdpanet[, (paste0('nfinal.', RCPS[i], '.', thismodnum)) := get(paste0('nfinal.', RCPS[i], '.', thismodnum)) + pfinal]
                    wdpanet[, (paste0('nshared.', RCPS[i], '.', thismodnum)) := get(paste0('nshared.', RCPS[i], '.', thismodnum)) + pinit*pfinal]
                    wdpanet[pinit > pfinal, (paste0('nlost.', RCPS[i], '.', thismodnum)) := get(paste0('nlost.', RCPS[i], '.', thismodnum)) + pinit - pfinal]
                    wdpanet[pinit < pfinal, (paste0('ngained.', RCPS[i], '.', thismodnum)) := get(paste0('ngained.', RCPS[i], '.', thismodnum)) + pfinal - pinit]
                }
                wdpanet[, c('pinit', 'pfinal') := NULL] # delete the species-specific columns now that we have the summaries we need
             }
        }
        # wdpanet[,print(max(get(paste0('ninit.', RCPS[i], '.1'))))] # a test that cumulative sums are being calculated. Should increase through the loop
    }
}
print(Sys.time())

dim(wdpanet) # 923 x 127

# how many rows are zero?
wdpanet[, .(sum(ninit.26.1 == 0), sum(ninit.85.1 == 0))]
wdpanet[, .(sum(ninit.26.1 > 0), sum(ninit.85.1 > 0))]
    # wdpanet[ninit.26.1 > 0, plot(lon_max, lat_max, cex = 0.1, xlim=c(-190, 0))]
    # wdpanet[ninit.26.1 == 0, points(lon_max, lat_max, cex = 0.5, col = 'red')]

# write out species data 
write.csv(wdpanet, file = gzfile('temp/wdpaturnbynetbymod.csv.gz'), row.names = FALSE)

# wdpanet <- fread('gunzip -c temp/wdpaturnbynetbymod.csv.gz') # don't read the row numbers

