# stats about impacts on protected areas from shifting species


############
## Flags
############



####################
## helper functions
####################
# require(Hmisc)
require(data.table)
# require(lme4) # for mixed-effects models
# require(car) # for testing ME models


lu <- function(x) return(length(unique(x)))

se <- function(x,na.rm=FALSE){ # standard error
	if(!na.rm){
		return(sd(x, na.rm=FALSE)/sqrt(length(x)))
	}
	if(na.rm){
		return(sd(x, na.rm=TRUE)/sqrt(sum(!is.na(x))))
	}
}

# weighted mean function to use with summarize()
wmean <- function(x){ # values in col 1, weights in col 2
	inds <- !is.na(x[,1]) & !is.na(x[,2])
	return(weighted.mean(x[inds,1], x[inds,2])) 
}


###########################################
# Examine MPA size relative to grid size
###########################################
wdpagrid <- fread('gunzip -c output/wdpa_cov_by_grid0.05.csv.gz', drop = 1) # shows which MPAs are in which grid cells. each line is a unique grid cell-MPA combination.


# plot vs. grid size
wdpagrid[!duplicated(WDPA_PID), hist(log10(area_fullwdpa), xlab = 'log10(area in m2)', main = 'MPA vs. grid cell sizes')]
abline(v = log10(unique(wdpagrid$area_grid)), col = '#FF000001')

# number > grid size
wdpagrid[!duplicated(WDPA_PID), .N] # 925
wdpagrid[!duplicated(WDPA_PID), sum(area_fullwdpa >= area_grid)] # 332
wdpagrid[!duplicated(WDPA_PID), sum(area_fullwdpa < area_grid)] # 593
wdpagrid[!duplicated(WDPA_PID), sum(area_fullwdpa < area_grid)/.N] # 64%

#############################################
# Examine turnover within MPAs
# Examine all climate models individually
#############################################
# read in turnover data
wdpaturnbyMPAbymod <- fread('gunzip -c temp/wdpaturnbyMPAbymod.csv.gz', drop=1) # Turnover in each MPA. Don't read the row numbers

ntk <- wdpaturnbyMPAbymod[,NO_TAKE %in% c('All', 'Part')] # no take reserves (index into wdpaturnbyMPAbymod)
sum(ntk) # 74 no take reserves
wdpaturnbyMPAbymod[ntk, sort(SUB_LOC)]
wdpaturnbyMPAbymod[ntk, sum(grepl('US-CA', SUB_LOC))] # 49 no-take in California

# calculate Sorenson turnover and fractional change
for (r in c(26, 85)) {
    for (m in 1:18) {
        wdpaturnbyMPAbymod[, (paste0('beta_sor.', r, '.', m)) := 
                               2*get(paste0('nshared.', r, '.', m)) /
                               (2*get(paste0('nshared.', r, '.', m)) + get(paste0('ngained.', r, '.', m)) + get(paste0('nlost.', r, '.', m)))]
        wdpaturnbyMPAbymod[, (paste0('flost.', r, '.', m)) := 
                               get(paste0('nlost.', r, '.', m)) /
                               get(paste0('ninit.', r, '.', m))]
        wdpaturnbyMPAbymod[, (paste0('fgained.', r, '.', m)) := 
                               get(paste0('ngained.', r, '.', m)) /
                               get(paste0('nfinal.', r, '.', m))]
    }
}

# convert to long format
wdpalong <- melt(wdpaturnbyMPAbymod, id.vars = c('WDPA_PID', 'NAME', 'MANG_AUTH', 'SUB_LOC', 'lat_min', 'lat_max', 'lon_min', 'lon_max', 'area_fullwdpa', 'network'), 
                 measure.vars = patterns('26|85'), variable.name = 'measure', value.name = 'val')
wdpalong[, c('measure', 'rcp', 'model') := tstrsplit(measure, '.', fixed = TRUE)] # split name of turnover measure apart from RCP and model #
    dim(wdpalong) # 266400 x 14

# read in change in temp by region
# trend1 <- read.csv('data/climTrendbyreg_rcp45.csv', row.names=1)
# 	trend1$rcp <- 45
# trend2 <- read.csv('data/climTrendbyreg_rcp85.csv', row.names=1)
# 	trend2$rcp <- 85
# nms <- c('region', 'rcp', 'delta_surf', 'delta_bott')
# trend <- as.data.table(rbind(trend1[,nms], trend2[,nms]))


# Fraction of species lost (fraction of original community)
    # all MPAs
wdpalong[measure == 'flost', .(mean = mean(val,na.rm = TRUE)), 
         by = c('rcp','WDPA_PID')][, .(mean = mean(mean), se = se(mean), min = min(mean), max = max(mean)), 
                                by = 'rcp'] # first mean is across climate models with MPAs/RCPs, then average across MPAs

wdpalong[measure == 'flost', .(mean = mean(val,na.rm = TRUE), lat = mean(c(lat_min, lat_max)), pac = (lon_min< -100)[1]), 
         by = c('rcp','WDPA_PID')][, plot(lat, mean, col = ifelse(rcp == 26, ifelse(pac, 'light blue', 'pink'),
                                                       ifelse(pac, 'blue', 'red')))] # plot MPA ensemble mean vs. lat. Ocean by color (Pac is blue, Atl is red), RCPs by shade (26 is light, 85 is dark.


	# no-take
wdpalong[ntk & measure == 'flost', .(mean = mean(val,na.rm = TRUE)), 
         by = c('rcp','WDPA_PID')][, .(mean = mean(mean), se = se(mean), min = min(mean), max = max(mean)), 
                                by = 'rcp'] # first mean is across climate models with MPAs/RCPs, then average across MPAs


# Fraction of species gained (fraction of final community)
	# all MPAs
wdpalong[measure == 'fgained', .(mean = mean(val,na.rm = TRUE)), 
         by = c('rcp','WDPA_PID')][, .(mean = mean(mean), se = se(mean), min = min(mean), max = max(mean)), 
                                by = 'rcp'] # first mean is across climate models with MPAs/RCPs, then average across MPAs

wdpalong[measure == 'fgained', .(mean = mean(val,na.rm = TRUE), lat = mean(c(lat_min, lat_max)), pac = (lon_min< -100)[1]), 
         by = c('rcp','WDPA_PID')][, plot(lat, mean, 
                                          col = ifelse(rcp == 26, ifelse(pac, 'light blue', 'pink'),
                                                       ifelse(pac, 'blue', 'red')))] # plot MPA ensemble mean vs. lat. Ocean by color (Pac is blue, Atl is red), RCPs by shade (26 is light, 85 is dark.


	# no-take
wdpalong[ntk & measure == 'fgained', .(mean = mean(val,na.rm = TRUE)), 
         by = c('rcp','WDPA_PID')][, .(mean = mean(mean), se = se(mean), min = min(mean), max = max(mean)), 
                                by = 'rcp'] # first mean is across climate models with MPAs/RCPs, then average across MPAs

                                   
# Similarity
	# all MPAs
wdpalong[measure == 'beta_sor', .(mean = mean(val,na.rm = TRUE)), 
         by = c('rcp','WDPA_PID')][, .(mean = mean(mean), se = se(mean), min = min(mean), max = max(mean)), 
                                by = 'rcp'] # first mean is across climate models with MPAs/RCPs, then average across MPAs

wdpalong[measure == 'beta_sor', .(mean = mean(val,na.rm = TRUE), lat = mean(c(lat_min, lat_max)), pac = (lon_min< -100)[1]), 
         by = c('rcp','WDPA_PID')][, plot(lat, mean, 
                                          col = ifelse(rcp == 26, ifelse(pac, 'light blue', 'pink'),
                                                       ifelse(pac, 'blue', 'red')))] # plot MPA ensemble mean vs. lat. Ocean by color (Pac is blue, Atl is red), RCPs by shade (26 is light, 85 is dark.


	# no-take
wdpalong[ntk & measure == 'beta_sor', .(mean = mean(val,na.rm = TRUE)), 
         by = c('rcp','WDPA_PID')][, .(mean = mean(mean), se = se(mean), min = min(mean), max = max(mean)), 
                                by = 'rcp'] # first mean is across climate models with MPAs/RCPs, then average across MPAs

	# examine similarity vs. MPA size
wdpalong[measure == 'beta_sor', cor.test(area_fullwdpa, val), by = c('rcp', 'model')]

wdpalong[measure == 'beta_sor', .(beta_sor = mean(val,na.rm = TRUE), log10size = mean(log10(area_fullwdpa))), 
         by = c('rcp','WDPA_PID')][, plot(log10size, beta_sor, 
                                          col = ifelse(rcp == 26, 'blue', 'red'))] # plot MPA ensemble mean vs. size. RCP by color.

				
	# examine similarity by MPA latitudinal range
wdpalong[measure == 'beta_sor', cor.test(lat_max - lat_min, val), by = c('rcp', 'model')]

wdpalong[measure == 'beta_sor', .(beta_sor = mean(val,na.rm = TRUE), latsize = mean(log10(lat_max - lat_min + 0.1))), 
         by = c('rcp','WDPA_PID')][, plot(latsize, beta_sor, 
                                          col = ifelse(rcp == 26, 'blue', 'red'))] # plot MPA ensemble mean vs. size. RCP by color.

		
# examine similarity by MPA size and lat rng
wdpalong[measure == 'beta_sor', .(beta_sor = mean(val, na.rm = TRUE), latsize = mean(log10(lat_max - lat_min + 0.1)), log10size = mean(log10(area_fullwdpa))), 
         by = c('rcp','WDPA_PID')][, summary(mod <- lm(beta_sor~rcp*latsize*log10size))] # linear model with 3-way interaction


		require(interplot) # to plot the interaction
		interplot(m=mod, var1='latsize', var2='log10size') +
			xlab('log(area) scaled') + 
			ylab('Coefficient for scaled latitudinal range')

		interplot(m=mod, var2='latrngsc', var1='lrep_areasc') +
			ylab('Coefficient for log(area) scaled') + 
			xlab('Latitudinal range scaled')



################################################################
# Examine change within MPA networks and the component MPAs
# Across each model in the ensemble
################################################################
wdpaturnbynetbymod <- fread('gunzip -c temp/wdpaturnbynetbymod.csv.gz') # network results
# also need wdpalong from above section
		
# Size of networks
	wdpaturnbyMPAbymod[!is.na(network),.(num_mpas=length(unique(WDPA_PID))), by=network]

# calculate Sorenson turnover and fractional change
for (r in c(26, 85)) {
    for (m in 1:18) {
        wdpaturnbynetbymod[, (paste0('beta_sor.', r, '.', m)) := 
                               2*get(paste0('nshared.', r, '.', m)) /
                               (2*get(paste0('nshared.', r, '.', m)) + get(paste0('ngained.', r, '.', m)) + get(paste0('ngained.', r, '.', m)))]
        wdpaturnbynetbymod[, (paste0('flost.', r, '.', m)) := 
                               get(paste0('nlost.', r, '.', m)) /
                               get(paste0('ninit.', r, '.', m))]
        wdpaturnbynetbymod[, (paste0('fgained.', r, '.', m)) := 
                               get(paste0('ngained.', r, '.', m)) /
                               get(paste0('nfinal.', r, '.', m))]
    }
}

# convert to long format
wdpanetlong <- melt(wdpaturnbynetbymod, id.vars = c('network', 'lat', 'lon', 'area'), 
                 measure.vars = patterns('26|85'), variable.name = 'measure', value.name = 'val')
wdpanetlong[, c('measure', 'rcp', 'model') := tstrsplit(measure, '.', fixed = TRUE)] # split name of turnover measure apart from RCP and model #
    dim(wdpanetlong) # 864 x 8

# Fraction of species lost
	# all MPA networks
wdpanetlong[measure == 'flost', .(mean = mean(val,na.rm = TRUE)), 
         by = c('rcp','network')][, .(mean = mean(mean), se = se(mean), min = min(mean), max = max(mean)), 
                                by = 'rcp'] # first mean is across climate models within networks x RCPs, then average across networks

	# within the individual MPAs of these networks
wdpalong[measure == 'flost' & !is.na(network), .(mean = mean(val,na.rm = TRUE)), 
         by = c('rcp','WDPA_PID')][, .(mean = mean(mean), se = se(mean), min = min(mean), max = max(mean)), 
                                by = 'rcp'] # first mean is across climate models within MPAs x RCPs, then average across MPAs

# Fraction of species gained (fraction of final community)
	# all MPA networks
wdpanetlong[measure == 'fgained', .(mean = mean(val,na.rm = TRUE)), 
         by = c('rcp','network')][, .(mean = mean(mean), se = se(mean), min = min(mean), max = max(mean)), 
                                by = 'rcp'] # first mean is across climate models within networks x RCPs, then average across networks

	# within the individual MPAs of these networks
wdpalong[measure == 'fgained' & !is.na(network), .(mean = mean(val,na.rm = TRUE)), 
         by = c('rcp','WDPA_PID')][, .(mean = mean(mean), se = se(mean), min = min(mean), max = max(mean)), 
                                by = 'rcp'] # first mean is across climate models within MPAs x RCPs, then average across MPAs

# Similarity (Sorenson)
	# all MPA networks
wdpanetlong[measure == 'beta_sor', .(mean = mean(val,na.rm = TRUE)), 
         by = c('rcp','network')][, .(mean = mean(mean), se = se(mean), min = min(mean), max = max(mean)), 
                                by = 'rcp'] # first mean is across climate models within networks x RCPs, then average across networks

	# within the individual MPAs of these networks
wdpalong[measure == 'beta_sor' & !is.na(network), .(mean = mean(val,na.rm = TRUE)), 
         by = c('rcp','WDPA_PID')][, .(mean = mean(mean), se = se(mean), min = min(mean), max = max(mean)), 
                                by = 'rcp'] # first mean is across climate models within MPAs x RCPs, then average across MPAs

