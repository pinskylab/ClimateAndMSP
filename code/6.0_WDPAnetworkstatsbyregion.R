# basic stats on the area of MPAs and the area and range of temperatures spanned by MPA networks

####################
## helper functions
####################
require(data.table)




#########################################################
# Calc stats for pre-defined MPA networks
# MOVED FROM 2.0_evalRandomNetworks.r
#########################################################
require(maptools)
require(rgeos)
require(rgdal)

# read in data
clim <- fread('gunzip -c output/climatology.csv.gz', drop = 1) # climatology from Morley et al. 2018
wdpagrid <- fread('gunzip -c output/wdpa_cov_by_grid0.05.csv.gz', drop = 1) # gridded MPA data. shows which MPAs are in which grid cells. each line is a unique grid cell-MPA combination.
networks <- fread('gunzip -c temp/wdpaturnbyMPAbymod.csv.gz', drop = 1) # MPA network descriptions
regions <- fread('gunzip -c output/region_grid.csv.gz', drop = 1) # CMSP region definitions from 5.0_define_CMSP.r

# round to analysis grids (0.25)
clim[, lat := floor(latClimgrid/0.25)*0.25 + 0.25/2]
clim[, lon := floor(lonClimgrid/0.25)*0.25 + 0.25/2]
wdpagrid[, lat := floor(lat/0.25)*0.25 + 0.25/2]
wdpagrid[, lon := floor(lon/0.25)*0.25 + 0.25/2]

# average climate by analysis grid
clim2 <- clim[, .(sbt = mean(sbt, na.rm = TRUE)), by = c('lat', 'lon')]

# convert all lon to -360 to 0 so that Aleutians stay together
wdpagrid[lon > 0, lon := lon - 360]
regions[longrid > 0, longrid := longrid - 360]

# merge climate grids with wdpa names, network definitions, and region definitions. keep all grid points, whether MPA or not
thesegrids <- merge(wdpagrid[, .(lat, lon, WDPA_PID, prop_wdpa)], clim2, by = c('lat', 'lon'), all = TRUE)
thesegrids <- merge(thesegrids, networks[, .(WDPA_PID, network)], by = 'WDPA_PID', all.x = TRUE)
thesegrids <- merge(thesegrids, regions[, .(latgrid, longrid, region)], by.x = c('lat', 'lon'), by.y = c('latgrid', 'longrid'), all = TRUE)

    # check that it worked
    thesegrids[, summary(sbt)]
    thesegrids[, sort(unique(region))]
    thesegrids[, sum(is.na(region))]
    thesegrids[,plot(lon, lat, col = c('black', 'red')[1 + is.na(sbt)], cex = 0.1)] # red for grid cells without SBT
    thesegrids[,plot(lon, lat, col = c('black', 'red')[1 + !is.na(network)], cex = 0.1)] # red for grid cells in a network
    thesegrids[!is.na(network),plot(lon, lat, cex = 0.1)] # only grid cells in a network
    
    
# calc size (# grid cells) and temperature range by region
regstats <- thesegrids[,.(temprng=diff(range(sbt, na.rm = TRUE)), size = length(unique(cbind(lat,lon)))), by = region]

# calc size and temperature range of each network by region
# size calculated as fractions of a grid cell
netstats <- thesegrids[!is.na(network),.(temprngnet = diff(range(sbt, na.rm = TRUE)), sizenet = sum(prop_wdpa)), by = .(region, network)]

# merge regional and network stats
setkey(regstats, region)
setkey(netstats, region)
stats <- netstats[regstats, .(region, network, temprng, size, temprngnet, sizenet), nomatch = 0][network != 'NA',]

# calc fraction of regional ranges covered by networks
stats[, fractemp := temprngnet/temprng]
stats[, fracsize := sizenet/size]

stats[,.(region, network, fracsize, fractemp)]


# save the area calculations
write.csv(stats, file='output/MPA_network_stats.csv')


