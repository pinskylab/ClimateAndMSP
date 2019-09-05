# Compare species p(occur) projections to species observations within MPAs

##################
# Load functions
##################
library(sf) # for GIS functions
library(data.table)

##################
# Read in data
##################

# read in list of focal MPAs
wdpa_by_grid <- readRDS('temp/wdpa_by_grid0.05_intersect.rds') # from 1.0_processWDPA.r
wdpalist <- sort(unique(as.character(wdpa_by_grid$WDPA_PID))) # pull out the WPDA ids we want
rm(wdpa_by_grid) # clean up

# read in MPA shapefile
wdpa <- st_read('dataDL/WDPA/WDPA_Aug2019_marine-shapefile/WDPA_Aug2019_marine-shapefile-polygons.shp')
wdpa <- wdpa[wdpa$WDPA_PID %in% wdpalist,] # Trim MPAs to focal MPAs

# trawl data from >=2016, as curated by OceanAdapt
trawl <- readRDS('dataDL/oceanadapt/all-regions-full.rds')
trawl <- trawl[trawl$year >= 2016,]
trawl <- st_as_sf(trawl, crs = st_crs(wdpa), coords = c('lon', 'lat'))
    # plot(trawl['region'], key.pos = 4, axes = TRUE, cex=0.02, key.width=lcm(5), asp=4)

#################################
# Find species in each focal MPA
# ONLY DO THIS ONCE
# THEN READ IN SAVED FILE INSTEAD SINCE TIME-CONSUMING
#################################

# Intersect trawl data with focal MPAs
trawl_by_wdpa <- st_intersection(trawl, wdpa) # 8 hours

# Save out
saveRDS(trawl_by_wdpa, file='temp/trawl_by_wdpa.rds')

#############################
# Make species list by MPA
#############################
# read in intersection of MPAs with trawl tow data
trawl_by_wdpa <- readRDS('temp/trawl_by_wdpa.rds')

# summarize by MPA
wdpaspp <- data.table(st_drop_geometry(trawl_by_wdpa))[, .(noccur=.N), by=.(WDPA_PID, NAME, SUB_LOC, spp)]

    wdpaspp[,hist(noccur, breaks = c(seq(-0.1, 10, by=1), 1000), xlim=c(0,11))] # most observations are of a species 1x in an MPA
    wdpaspp[, .(nspp=.N), by=WDPA_PID][,hist(nspp, breaks = 20)] # most MPAs have 50-100 spp

# write out
write.csv(wdpaspp, file = gzfile('output/wdpa_trawlsppobs.csv.gz'), row.names = FALSE)

