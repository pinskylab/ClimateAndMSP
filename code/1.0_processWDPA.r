# Intersect the WDPA data with our analysis grid

require(maptools)
require(RColorBrewer)
require(rgdal)
require(rgeos)
library(sf) # for st_intersection
require(raster) # for crop

#######################
## Make gridded WDPA data ##
#######################

# read in the analysis grid from a species projection datafile
load('dataDL/morley/acanthephyra pelagica_Atl_rcp26_jas_prediction_AGG.RData') # an Atlantic projection grid. loads pred.agg
Atl<- pred.agg[!duplicated(pred.agg[, c('latitude', 'longitude')]), c('latitude', 'longitude')]
    #plot(Atl$longitude, Atl$latitude)
load('dataDL/morley/actinauge verrilli_Pac_rcp26_jas_prediction_AGG.RData') # a Pacific projection grid. loads pred.agg
Pac<- pred.agg[!duplicated(pred.agg[, c('latitude', 'longitude')]), c('latitude', 'longitude')]
    #plot(Pac$longitude, Pac$latitude, cex=0.1)

clim <- rbind(Atl, Pac)
rm(pred.agg, Atl, Pac)
clim$lon <- clim$longitude
clim$lat <- clim$latitude
clim <- clim[,c('lat', 'lon')]

clim$lon[clim$lon>180] = clim$lon[clim$lon>180] - 360 # convert lon to shp format (-180 to 180)
clim$lon[clim$lon< -180] = clim$lon[clim$lon< -180] + 360 # convert lon to shp format (-180 to 180)
gridsz = unique(round(diff(sort(unique(clim$lat))),3)) # grid size 0.05
gridsz <- gridsz[gridsz>0]


# read in the MPAs
wdpa = readOGR(dsn='dataDL/WDPA/WDPA_Aug2019_marine-shapefile', layer='WDPA_Aug2019_marine-shapefile-polygons')
# plot(wdpa) # very slow
nrow(wdpa)
wdpa = wdpa[wdpa$MARINE != 0,] # trim out 100% terrestrial (remove 185)
nrow(wdpa)

# some data exploration
#table(wdpa$MARINE)
#sort(unique(wdpa$DESIG))
#sort(unique(wdpa$DESIG_ENG))
#as.matrix(table(wdpa$IUCN_CAT))
#as.matrix(table(wdpa$NO_TAKE))


# Make the grid based on corner points
y = numeric(5*length(clim$lat)) # all the Y coords, in order
# lower right, lower left, upper left, upper right, lower right
for(i in 1:length(clim$lat)){ y[(5*(i-1)+1):(5*(i-1)+5)] <- c(clim$lat[i]-gridsz/2, clim$lat[i]-gridsz/2, clim$lat[i]+gridsz/2, clim$lat[i]+gridsz/2, clim$lat[i]-gridsz/2)}

x = numeric(5*length(clim$lon))
# lower right, lower left, upper left, upper right, lower right: clockwise so that sp sees it as an island, not a hole
for(i in 1:length(clim$lon)){ x[(5*(i-1)+1):(5*(i-1)+5)] = c(clim$lon[i]+gridsz/2, clim$lon[i]-gridsz/2, clim$lon[i]-gridsz/2, clim$lon[i]+gridsz/2, clim$lon[i]+gridsz/2) }

# Create a SpatialPolygonsDataFrame (a list of "Polygons", each of which is a list of "Polygon")
# Would be faster to create sf object directly
pgns = vector('list', length(clim$lat))
for(i in 1:length(pgns)){
    inds2 = (5*(i-1)+1):(5*(i-1)+5)
    pgns[[i]] = Polygons(list(Polygon(cbind(x[inds2],y[inds2]))), i)
}
SP <- SpatialPolygons(pgns, proj4string=CRS(proj4string(wdpa)))
SPdata <- data.frame(gridpolyID = as.numeric(sapply(slot(SP, 'polygons'), slot, 'ID')), lat = clim$lat, lon = clim$lon) # would be better to put this in SP as a SpatialPolygonsDataFrame
    length(SP)
SPsf <- st_as_sf(SP)
SPsf2 <- dplyr::bind_cols(SPsf, SPdata)
#plot(SP[1:10])
#plot(SP[1:1000]) # slow for so many
#plot(SP) # very slow

# Write out SP for use later
saveRDS(SPsf2, file='temp/SPsf2.rds')

# Trim wdpa to our analysis area
wdpacrop <- crop(wdpa, SP) # drop from 14561 to 925 elements. Slow but worth it (10 min?)
rm(wdpa)

# Find which grids intersect which PAs
# note, it would probably be faster to skip this just do the full intersection directly: out <- st_intersection(SPsf2, wdpasf). Would take a few hours, but not overnight
wdpasf <- st_as_sf(wdpacrop)        
gI <- st_intersects(SPsf2, wdpasf) # results stored as sparse binary matrix. list elements show TRUE columns for a given row: gI[[1]] is row one
    dim(gI) # matches dim SP x wdpa
    sum(sapply(gI, length)) # number of pairwise intersections 73659
ng = sum(sapply(gI, length)>0); ng # number of SP polygons that intersect at least one wdpa element
cols = which(sapply(gI, length)>0) # ids of SP polygons that intersect wdpa

# do the intersection
# takes ovenight
warnings <- character(0) # to hold warnings from the loop
messages <- character(0) # to hold messages from the loop
withCallingHandlers({ # store messages and warnings for viewing later, rather than scrolling to screen
    for(i in 1:ng){ # Steps through each grid cell
        if(i %% 1000 == 0) print(paste(i, 'of', ng))
        if(i==1){ # initialize the out object
            out <- st_intersection(SPsf2[cols[i],], wdpasf[gI[[cols[i]]],]) # only calc intersections for those polygons that actually intersect (according to st_intersects)
        } else { # add to the out object
            temp <- st_intersection(SPsf2[cols[i],], wdpasf[gI[[cols[i]]],]) # only calc intersections for those polygons that actually intersect (according to st_intersects)
            out <- rbind(out, temp)
        }
    }
}, warning=function(w){
    warnings <<- c(warnings, w$message) # catch warnings and save them
    invokeRestart("muffleWarning") # don't report the warning now
}, error=function(e){
    stop(e) # stop if there's an error
}, message=function(m){
    messages<<-c(messages,m$message)
    invokeRestart("muffleMessage")
})

dim(out)

# save the result
saveRDS(out, file=paste('temp/wdpa_by_grid', gridsz, '_intersect.rds', sep='')) # save
