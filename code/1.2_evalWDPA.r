# evaluate protected areas against shifts in species distribution
# calculate species gains, losses, turnover, etc.


############
## Flags
############

# choose the rcp(s)
RCPS <- c(26, 85)

# select initial and final timeperiod for these grids
PERIODS <- c('2007-2020', '2081-2100')

# number of climate models in the projections
NMODS <- 18

# path to the pres-abs projection results from Morley et al. 2018. At 0.05 grid size.
PRESABSPATH <- '/local/shared/pinsky_lab/projections_PlosOne2018/CEmodels_proj_PresAbs_May2018' 

# set default rounding (2 bytes) so that data.table will merge numeric values appropriately
# see https://rdrr.io/rforge/data.table/man/setNumericRounding.html
setNumericRounding(2)

####################
## helper functions
####################
#require(Hmisc)
require(data.table)
#require(lme4) # for mixed-effects models
#require(car) # for testing ME models



###########################
## Load and set up WDPA data
###########################

wdpa <- fread('gunzip -c output/wdpa_cov_by_grid0.05.csv.gz', drop = 1) # shows which MPAs are in which grid cells. each line is a unique grid cell-MPA combination.
	
	# convert lon to -360 to 0 (instead of -180 to 180) to match the pres/abs projections and to plot more nicely
	wdpa[lon>0, lon := lon-360]

	# set up MPA networks
	wdpa[SUB_LOC == 'US-CA' & MANG_AUTH == 'State Fish and Wildlife',network := 'mlpa']
	wdpa[(SUB_LOC %in% c('US-DE', 'US-FL', 'US-GA', 'US-MA', 'US-MD', 'US-ME', 'US-NC', 'US-NJ', 'US-NY', 'US-RI', 'US-SC', 'US-VA')) & (lon > 280), network := 'eastcoast'] # by choosing by state, this won't include federal water closures
	wdpa[(SUB_LOC %in% c('US-AK')), network := 'ak'] # by choosing by state, this won't include federal water closures
	
	# calculate mpa extent
	wdpa[, ':='(lat_min = min(lat), lat_max = max(lat)), by = WDPA_PID]

#	wdpa[network=='mlpa',plot(lon,lat)]
#	wdpa[network=='eastcoast',plot(lon,lat)]
#	wdpa[network=='ak',plot(lon,lat)]
#	require(maps); map(database='world2',add=TRUE)

##############################################################
# Load in pres/abs results for each RCP and each species
# and calculate summary statistics by grid cell
##############################################################

projcols <- c('lat', 'lon', paste0('mean', 1:NMODS)) # names of columns in pred.agg to use (the projections)
for (i in 1:length(RCPS)) {
    print(paste0(Sys.time(), ' On RCP ', RCPS[i], '. Species: '))
    files <- list.files(path = PRESABSPATH, pattern = paste('*rcp', RCPS[i], '*', sep = ''), full.names = TRUE)

    # load presmap for each model run and process the results
    # modify this to read in files from /local/home/jamesm/Documents/range_projections/CEmodels_proj_PresAbs_May2018/ instead
    # then calculate change in p(occur)
    for (j in 1:length(files)) {
        cat(paste0(' ', j))
        load(files[j]) # loads pred.agg data.frame
        pred.agg <- as.data.table(pred.agg)
        
        # round lat lon to nearest 0.05 so that merging is successful with WDPA
        pred.agg[, lat := floor(latitude*20)/20 + 0.025] # assumes 0.05 grid size
            pred.agg[,max(abs(latitude - lat))] # check: should be very small (1e-14 or so)
        pred.agg[, lon := floor(longitude*20)/20 + 0.025] # assumes 0.05 grid size
            pred.agg[,max(abs(longitude - lon))] # check: should be very small
        setkey(pred.agg, lat, lon)
        
        # reorganize so that time 2 and time 1 are separate columns, rows are locations
        presbyloc <- merge(pred.agg[year_range == PERIODS[1], ..projcols], pred.agg[year_range == PERIODS[2], ..projcols], by = c('lat', 'lon'), suffixes = c('.t1', '.t2'))

        # merge pres/abs data with WDPA and calculate summary stats about species change through time
        # summary stat col names are of form [sumstatname].[RCP].[modnumber]
        wdpa <- merge(wdpa, presbyloc, all.x = TRUE, by = c('lat', 'lon'))
        for (m in 1:NMODS) { # for each climate model
            wdpa[is.na(wdpa[[paste0('mean', m, '.t1')]]), (paste0('mean', m, '.t1')) := 0] # set NAs to 0 (species not present)
            wdpa[is.na(wdpa[[paste0('mean', m, '.t2')]]), (paste0('mean', m, '.t2')) := 0]
            if (!(paste0('ninit.', RCPS[i], '.', m) %in% colnames(wdpa))) { # of summary columns don't yet exist, create them
                wdpa[, (paste0('ninit.', RCPS[i], '.', m)) := get(paste0('mean', m, '.t1'))] # initial number of species: sum of p(occur) values
                wdpa[, (paste0('nfinal.', RCPS[i], '.', m)) := get(paste0('mean', m, '.t2'))] # final number
                wdpa[, (paste0('nshared.', RCPS[i], '.', m)) := get(paste0('mean', m, '.t1'))*get(paste0('mean', m, '.t2'))] # number of shared species: sum of p(shared species)
                wdpa[wdpa[[paste0('mean', m, '.t1')]] > wdpa[[paste0('mean', m, '.t2')]], (paste0('nlost.', RCPS[i], '.', m)) := get(paste0('mean', m, '.t1')) - get(paste0('mean', m, '.t2'))] # number of species lost: sum of p(occur at time 1) - p(occur at time 2) for species whose p(occur) declines through time
                wdpa[wdpa[[paste0('mean', m, '.t1')]] <= wdpa[[paste0('mean', m, '.t2')]], (paste0('nlost.', RCPS[i], '.', m)) := 0] # 0 otherwise
                wdpa[wdpa[[paste0('mean', m, '.t1')]] < wdpa[[paste0('mean', m, '.t2')]], (paste0('ngained.', RCPS[i], '.', m)) := get(paste0('mean', m, '.t2')) - get(paste0('mean', m, '.t1'))] # number of species gained
                wdpa[wdpa[[paste0('mean', m, '.t1')]] >= wdpa[[paste0('mean', m, '.t2')]], (paste0('ngained.', RCPS[i], '.', m)) := 0]
            } else {# if columns already list, add the new species' values onto the existing values
                wdpa[, (paste0('ninit.', RCPS[i], '.', m)) := get(paste0('ninit.', RCPS[i], '.', m)) + get(paste0('mean', m, '.t1'))]
                wdpa[, (paste0('nfinal.', RCPS[i], '.', m)) := get(paste0('nfinal.', RCPS[i], '.', m)) + get(paste0('mean', m, '.t2'))]
                wdpa[, (paste0('nshared.', RCPS[i], '.', m)) := get(paste0('nshared.', RCPS[i], '.', m)) + get(paste0('mean', m, '.t1')) * get(paste0('mean', m, '.t2'))]
                wdpa[wdpa[[paste0('mean', m, '.t1')]] > wdpa[[paste0('mean', m, '.t2')]], (paste0('nlost.', RCPS[i], '.', m)) := get(paste0('nlost.', RCPS[i], '.', m)) + get(paste0('mean', m, '.t1')) - get(paste0('mean', m, '.t2'))]
                wdpa[wdpa[[paste0('mean', m, '.t1')]] < wdpa[[paste0('mean', m, '.t2')]], (paste0('ngained.', RCPS[i], '.', m)) := get(paste0('ngained.', RCPS[i], '.', m)) + get(paste0('mean', m, '.t2')) - get(paste0('mean', m, '.t1'))]
            }
            wdpa[, (paste0('mean', m, '.t1')) := NULL] # delete the new columns now that we have the summaries we need
            wdpa[, (paste0('mean', m, '.t2')) := NULL]
        }
        # wdpa[,print(max(get(paste0('ninit.', RCPS[i], '.1'))))] # a test that cumulative sums are being calculated. Should increase through the loop
    }
}
print(Sys.time())

dim(wdpa) # 73659 x 307

# how many rows are zero?
wdpa[, .(sum(ninit.26.1 == 0), sum(ninit.85.1 == 0))]
wdpa[, .(sum(ninit.26.1 > 0), sum(ninit.85.1 > 0))]
    # wdpa[ninit.26.1 > 0, plot(lon, lat, cex = 0.1, xlim=c(-180, 180))]
    # wdpa[ninit.26.1 == 0, points(lon, lat, cex = 0.5, col = 'red')]

# write out species data 
write.csv(wdpa, file = gzfile('temp/wdpaturnbyMPAbymod.csv.gz'))

# wdpa <- fread('gunzip -c temp/wdpaturnbyMPAbymod.csv.gz', drop = 1) # don't read the row numbers



##############################################
# Calculate turnover within select MPA networks
# Across all models in the ensemble
##############################################

# merge in MPAnetwork ID (using data.tables)
wdpanets <- wdpa[!is.na(network),] # trim out non-network rows
setkey(wdpanets, lat, lon)
setkey(thesesppbymod, lat, lon)
thesesppbymodnet <- wdpanets[thesesppbymod,.(sppocean, period, rcp, model, pres, network, lat, lon), allow.cartesian=TRUE, nomatch=0] # datatable join, using the keys set in each dt. the order specifies to keep all rows of thesesppbymod. allow.cartesian is needed because the result is larger than either parent. drop locations not in a network.
	dim(thesesppbymodnet) # 2,637,693? 5,816,336, now 8995646 (4/24/2016)
	
# summarize by unique species in each period in each MPA in each RCP in each model (pres or not)
thesesppbynetbymod <- thesesppbymodnet[,max(pres),by="sppocean,period,network,rcp,model"] # DT aggregate function
	thesesppbynetbymod[,pres:=as.logical(V1)]
	thesesppbynetbymod[,V1:=NULL] # drop V1
	dim(thesesppbynetbymod) # 30,077

# calculate turnover by MPA, rcp, and model. very quick.
thesesppbynetbymod[,nstart:=nstart(period, sppocean, pres, periods), by="network,rcp,model"]
thesesppbynetbymod[,nend:=nend(period, sppocean, pres, periods), by="network,rcp,model"]
thesesppbynetbymod[,nlost:=nlost(period, sppocean, pres, periods), by="network,rcp,model"]
thesesppbynetbymod[,ngained:=ngained(period, sppocean, pres, periods), by="network,rcp,model"]
thesesppbynetbymod[,flost:=flost(period, sppocean, pres, periods), by="network,rcp,model"]
thesesppbynetbymod[,fgained:=fgained(period, sppocean, pres, periods), by="network,rcp,model"]
thesesppbynetbymod[,fgainedalt:=fgainedalt(period, sppocean, pres, periods), by="network,rcp,model"]
thesesppbynetbymod[,beta_sor:=beta_sor(period, sppocean, pres, periods), by="network,rcp,model"]

# trim out duplicate rows within an mpa/rcp/model
setkey(thesesppbynetbymod, network,rcp,model)
wdpaturnbynetbymod <- unique(thesesppbynetbymod)
	dim(wdpaturnbynetbymod) # 78 (3 networks x 13 models x 2 scenarios)
	wdpaturnbynetbymod[,c('sppocean','period','pres'):=NULL]
	
# merge in mpa metadata
setkey(wdpanets, network)
uwdpanets <- unique(wdpanets) # only the lines with unique mpas
	dim(uwdpanets) # 3
setkey(uwdpanets, network)
setkey(wdpaturnbynetbymod, network,rcp,model)
wdpaturnbynetbymod <- uwdpanets[wdpaturnbynetbymod, .(network,rcp,model,nstart,nend,nlost,ngained,flost,fgained,fgainedalt,beta_sor)]
	dim(wdpaturnbynetbymod) # 78


# write out wdpaturnbyMPA
write.csv(wdpaturnbynetbymod, paste('data/wdpaturnbynetbymod_', runtype, projtype, '_rcp', rcp, '&', otherrcp, '.csv', sep=''))

