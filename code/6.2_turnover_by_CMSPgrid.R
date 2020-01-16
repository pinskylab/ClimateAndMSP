# calculate species turnover for each grid cell in the CMSP analysis
# use the ensemble mean


############
## Flags
############

# choose the rcps (get to choose two)
rcps <- c(26, 85)

# select initial and final timeperiod for these grids
periods <- c('2007-2020', '2081-2100')


####################
## helper functions
####################
require(data.table)




##################################################
# Calculate turnover
##################################################
# loads ensemble mean presence/absence data for each species/model/rcp
for(oce in c('Atl', 'Pac')){
	for (j in 1:length(rcps)){
	    for(k in 1:length(periods)){
	        cat(paste0('\tLoading rcp', rcps[j], ' ocean', oce, ' period', periods[k], '\n'))
	        prestemp <- fread(cmd = paste0('gunzip -c temp/presmap_', oce, '_rcp', rcps[j], '_', periods[k], '.csv.gz'), drop = 1)
            prestemp <- prestemp[, .(poccur = mean(poccur)), by = .(latgrid, longrid, year_range, spp)] # average across gcms
	        
	        if(oce == 'Atl' & j == 1 & k == 1){
	            presmap <- prestemp
	        } else {
	            presmap <- rbind(presmap, prestemp)
	        }
	    }
	}
}
rm(prestemp)

# average across rcps
presmap <- presmap[, .(poccur = mean(poccur)), by = .(latgrid, longrid, year_range, spp)]

# calculate change in poccur for each species/grid
presmapdelta <- dcast(presmap, spp + latgrid + longrid ~ year_range, value.var = 'poccur') # reshape to wide format
presmapdelta[, dpoccur := get(periods[2]) - get(periods[1])] # calculate the change
presmapdelta[, pshared := get(periods[2]) * get(periods[1])] # calculate the probability of being shared

# evaluate ecological turnover for each grid cell
sppdelta <- presmapdelta[, .(nshared = sum(pshared)), by = .(latgrid, longrid)] # number of spp shared
sppdelta <- merge(sppdelta, presmapdelta[dpoccur > 0, .(ngained = sum(dpoccur)), by = .(latgrid, longrid)]) # spp gained
sppdelta <- merge(sppdelta, presmapdelta[dpoccur < 0, .(nlost = - sum(dpoccur)), by = .(latgrid, longrid)]) # spp lost
sppdelta[, beta_sor := 1 - 2*nshared / (2*nshared + nlost + ngained)] # sorenson dissimilarity


# write out
write.csv(sppdelta, file='output/turnover_by_CMSPgrid.csv')


