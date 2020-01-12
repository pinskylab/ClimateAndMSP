# evaluate protected areas against shifts in species distribution
# calculate species gains, losses, turnover, etc.


############
## Flags
############

# choose the rcps (get to choose two)
rcps <- c(26, 85)

# select initial and final timeperiod for these grids
periods <- c('2007-2020', '2081-2100')

# network parameters
temps <- c(0.01, seq(0.1, 1,length.out = 10)) # how many steps of temperature extent (0 to 100%)
sizes <- c(0.01, seq(0.05, 0.5, length.out = 10)) # which steps of area (could be 0 to 100%).
nreps <- 3 # how many repeats at each combination of lat and size


####################
## helper functions
####################
require(data.table)




##################################################
# Generate random MPA networks
# Across gradients of area and temperature extent
##################################################
clim <- fread('gunzip -c output/climatology.csv.gz', drop = 1) # climatology from Morley et al. 2018
regions <- fread('gunzip -c output/region_grid.csv.gz', drop = 1) # CMSP region definitions from 5.0_define_CMSP.r

# round to analysis grids (0.25) and aggregate
clim[, lat := floor(latClimgrid/0.25)*0.25 + 0.25/2]
clim[, lon := floor(lonClimgrid/0.25)*0.25 + 0.25/2]
clim <- clim[, .(sbt = mean(sbt)), by = .(lat, lon)]

# merge climatology and region: the planning units
clim <- merge(clim, regions[, .(latgrid, longrid, region)], by.x = c('lat', 'lon'), by.y = c('latgrid', 'longrid'))

# the regions to analyze
regs <- as.character(regions[,sort(unique(region))])

# data frame to store results
randMPAs <- as.data.table(expand.grid(region = regs, tempstep = temps, sizestep = sizes, repnum = 1:nreps, temprng = NA_real_, 
                                      size = NA_real_, meanturnIndiv = NA_real_, sdturnIndiv = NA_real_, netwrkturn = NA_real_, 
                                      sdnetwrkturn = NA_real_))

# a slow loop, especially on higher values of j and k
for(i in 1:length(regs)){
	print(regs[i])

	# set up the set of planning units (lat, lon, sbt)
	punits <- clim[region == regs[i],]

	tx <- diff(range(punits$sbt, na.rm = TRUE)) # the range of temperatures in the region
	nu <- sum(!duplicated(punits[, .(lat, lon)])) # the number of planning units
	if(nu != nrow(punits)) stop(paste0('i=', i, '. Duplicated entries in punits'))

	# loads presence/absence data for each species/model/rcp
	if(regs[i] %in% c('ebs', 'goa', 'bc', 'wc')) ocean <- 'Pac'
	if(regs[i] %in% c('gmex', 'seus', 'neus', 'maritime', 'newf')) ocean <- 'Atl'
	for (j in 1:length(rcps)){
	    for(k in 1:length(periods)){
	        cat(paste0('\tLoading rcp', rcps[j], ' ocean', ocean, ' period', periods[k], '\n'))
	        prestemp <- fread(cmd = paste0('gunzip -c temp/presmap_', ocean, '_rcp', rcps[j], '_', periods[k], '.csv.gz'), drop = 1)

	        if(j == 1 & k == 1){
	            presmap <- prestemp
	        } else {
	            presmap <- rbind(presmap, prestemp)
	        }
	    }
	}
    rm(prestemp)

	# create and evaluate the random networks
	for(j in 1:length(temps)){ # for each temperature step
		cat(paste('\nj=',j, sep=''))

	    for(k in 1:length(sizes)){ # for each size step
			cat(paste('\tk=',k,': ', sep=''))
	        
			for(r in 1:nreps){ # for each repeat
				cat(paste('r=',r,' '))
	
				# randomly designate a first planning unit as MPA
				punits[, sel := FALSE]
				punits[sample(1:.N,1) , sel := TRUE] # pick a first MPA
				ft <- punits[sel == TRUE, sbt] # the temperature

				# calculate nearby planning units from first MPA (based on temp extent)
				punits[, nearby := FALSE] # initialize "nearby" flag (within a certain fraction of lat extent)
				punits[abs(sbt-ft)/tx <= temps[j], nearby := TRUE] # mark some as nearby in temperature space
				

				# add more MPAs until size criterion met or all nearby units selected
				fracwholearea <- 1/nu
				fracnearbyarea <- 1/sum(punits$nearby)
				while(fracwholearea <= sizes[k] & fracnearbyarea <= 1 & sum(punits$nearby & !punits$sel)>0 ) {
					newunit <- punits[, sample(which(nearby & !sel), 1)] # select a new unit
					punits[newunit, sel:= TRUE] # assign as selected
					fracwholearea <- sum(punits$sel)/nu
					fracnearbyarea <- sum(punits$sel)/sum(punits$nearby)
					minsbtsel <- punits[sel == TRUE, min(sbt)]
					maxsbtsel <- punits[sel == TRUE, max(sbt)]
					punits[, nearby := ((sbt - minsbtsel)/tx <= temps[j] & sbt >= minsbtsel) | ((maxsbtsel - sbt)/tx < temps[j] & sbt <= maxsbtsel)] # measure nearby from the min temp if new potential temps are > min, and measure from max temp if new potential temps are < max
				}
				
				
				# merge selected MPA list with species data
				thesesppbymod <- merge(punits[sel == TRUE, .(lat, lon)], presmap, by.x = c('lat', 'lon'), by.y = c('latgrid', 'longrid'), all.x = TRUE)

				# calculate change in poccur for each species/model/rcp/grid
				thesesppbymod <- dcast(thesesppbymod, spp + rcp + model + lat + lon ~ year_range, value.var = 'poccur') # reshape to wide format
                thesesppbymod[, dpoccur := get(periods[2]) - get(periods[1])] # calculate the change
                thesesppbymod[, pshared := get(periods[2]) * get(periods[1])] # calculate the probability of being shared

				# basic calcs for this random network
				ind <- randMPAs[, which(region == regs[i] & tempstep == temps[j] & sizestep == sizes[k] & repnum == r)]
				randMPAs[ind, temprng := punits[sel == TRUE, (max(sbt)-min(sbt))/tx]]
				randMPAs[ind, size := punits[, sum(sel)/nu]]
				
				# evaluate ecological turnover for each individual MPA
				mpatempstats <- thesesppbymod[, .(nshared = sum(pshared)), by = .(lat, lon, rcp, model)]
				mpatempstats <- merge(mpatempstats, thesesppbymod[dpoccur > 0, .(ngained = sum(dpoccur)), by = .(lat, lon, rcp, model)])
				mpatempstats <- merge(mpatempstats, thesesppbymod[dpoccur < 0, .(nlost = - sum(dpoccur)), by = .(lat, lon, rcp, model)])
				mpatempstats[, beta_sor := 1 - 2*nshared / (2*nshared + nlost + ngained)] # sorenson dissimilarity
				randMPAs$meanturnIndiv[ind] <- mpatempstats[, mean(beta_sor)]
				randMPAs$sdturnIndiv[ind] <- mpatempstats[, .(beta_sor = mean(beta_sor)), by = .(lat, lon)][, sd(beta_sor)] # sd across models and rcps
				
				# ecological turnover for the entire network
                thesesppbymodnet <- thesesppbymod[, .(pinit = 1 - prod(1 - get(periods[1])), pfinal = 1 - prod(1 - get(periods[2]))), by = .(spp, rcp, model)] # aggregate probability for each species across the full network at each time period
                thesesppbymodnet[, dpoccur := pfinal - pinit] # calculate the change
                thesesppbymodnet[, pshared := pfinal * pinit] # calculate the probability of being shared
                mpatempstatsnet <- thesesppbymodnet[, .(nshared = sum(pshared)), by = .(rcp, model)]
                mpatempstatsnet <- merge(mpatempstatsnet, thesesppbymodnet[dpoccur > 0, .(ngained = sum(dpoccur)), by = .(rcp, model)])
                mpatempstatsnet <- merge(mpatempstatsnet, thesesppbymodnet[dpoccur < 0, .(nlost = - sum(dpoccur)), by = .(rcp, model)])
                mpatempstatsnet[, beta_sor := 1 - 2*nshared / (2*nshared + nlost + ngained)] # sorenson dissimilarity
                randMPAs$netwrkturn[ind] <- mpatempstatsnet[, mean(beta_sor)]
                randMPAs$sdnetwrkturn[ind] <- mpatempstatsnet[, sd(beta_sor)] # sd across models and rcps
                
			}			
		}
	}
	cat('\n')
}

# write out
write.csv(randMPAs, file='output/randMPAs_byBT.csv')


