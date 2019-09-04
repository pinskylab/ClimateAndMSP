# Compare species p(occur) projections to species observations within MPAs

##################
# Load functions
##################
library(sf) # for GIS functions

##################
# Read in data
##################

# read in list of focal MPAs
wdpa_by_grid <- readRDS('temp/wdpa_by_grid0.05_intersect.rds') # from 1.0_processWDPA.r
wdpalist <- sort(unique(wdpa_by_grid$WDPA_PID))
rm(wdpa_by_grid)

# read in MPA shapefile
wdpa = st_read('dataDL/WDPA/WDPA_Aug2019_marine-shapefile/WDPA_Aug2019_marine-shapefile-polygons.shp')

# trawl data from >=2016, as curated by OceanAdapt
trawl <- data.table(readRDS('dataDL/oceanadapt/all-regions-full.rds'))
trawl <- trawl[year>=2016,]
trawl <- st_sf(trawl, st_multipoint(as.matrix(cbind(trawl$lon, trawl$lat))), crs=st_crs(wdpa))

#################################
# Find species in each focal MPA
#################################

# Trim MPAs to focal MPAs
wdpa <- wdpa[wdpa$WDPA_PID %in% wdpaturnbyMPAbymod$WDPA_PID,]

# Intersect trawl data with focal MPAs
trawl_by_wdpa <- st_intersection(trawl, wdpa)