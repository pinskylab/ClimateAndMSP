# Define planning regions and choose which species count as important to fisheries


###########################
## Define planning regions
###########################
require(data.table)
require(sf)

# read in data
grid <- readRDS('temp/SPsf2.rds') # analysis grid
eez <- st_read('dataDL/marineregions/World_EEZ_v10_20180221/eez_v10.shp')
coast <- st_read('dataDL/natcap/Marine/Land/global_polyline.shp') # useful for checking success by making maps

# get ready to define regions
indsCAN <- st_intersects(eez[eez$Pol_type == '200NM' & eez$Territory1 == 'Canada',], grid, sparse = FALSE) # slow...
grid$region <- NA

# define each region (order is important, as later ones build off earlier ones)
# NEUS
grid$region[grid$lat > 35.6 & grid$lat < 45 & grid$lon > - 100 & grid$lon < -60 & !indsCAN] <- 'neus'

# Newfoundland
grid$region[(indsCAN & grid$lon > -100 & (grid$lon > -53 | (grid$lat > 49 & grid$lon > -57) | (grid$lat > 52) | (grid$lat > 47.8 & grid$lon > -54) | (grid$lat > 47 & grid$lon > -53.5))) | (grid$lon > -53 & grid$lon < -40)] <- 'newf'

# maritime Canada (including French waters around St. Pierre & Miquelon)
grid$region[(indsCAN & is.na(grid$region) & grid$lon > -100) | (grid$lat > 44 & grid$lat < 48 & grid$lon > -57.2 & grid$lon < -55)] <- 'maritime'

# SEUS
grid$region[(grid$lon > -80.6 | (grid$lon > -82.3 & grid$lat > 28)) & grid$lat <= 35.6] <- 'seus'

# GoMex
grid$region[grid$lon > -100 & grid$lat < 31 & is.na(grid$region)] <- 'gomex'

# west coast US
grid$region[grid$lon < -100 & grid$lat < 50 & !indsCAN] <- 'wc'

# British Columbia
grid$region[grid$lon < -100 & indsCAN] <- 'bc'

# GoAK
grid$region[grid$lon < -100 & grid$lon > -156 & is.na(grid$region)] <- 'goa'

# EBS and Aleutians
grid$region[grid$lon <= -156 | grid$lon > 0] <- 'ebs'

# any missing?
sum(is.na(grid$region))
    with(grid[is.na(grid$region), ], plot(lon, lat, cex=0.1))


# plot
with(grid[grid$region == 'newf', ], plot(lon, lat, cex=0.1))
with(grid[grid$region == 'maritime', ], plot(lon, lat, cex=0.1))
with(grid[grid$region == 'neus', ], plot(lon, lat, cex=0.1))
with(grid[grid$region == 'seus', ], plot(lon, lat, cex=0.1))
with(grid[grid$region == 'gomex', ], plot(lon, lat, cex=0.1))
with(grid[grid$region == 'wc', ], plot(lon, lat, cex=0.1))
with(grid[grid$region == 'bc', ], plot(lon, lat, cex=0.1))
with(grid[grid$region == 'goa', ], plot(lon, lat, cex=0.1)) # odd chunk out of SW corner? but nearly all is land
    plot(st_geometry(coast), add=TRUE)
with(grid[grid$region == 'ebs', ], plot(lon, lat, cex=0.1))

# write out
saveRDS(grid, 'output/region_grid.rds')


###########################################################
# Determine which species are of fishery importance
# from Sea Around Us fishery catch data
###########################################################
require(data.table)

## get list of projected species
wdpa_by_spp_projnow <- readRDS('temp/wdpa_by_spp_projnow.rds') # from 1.2_evalWDPA.r. A convenient summary of species.
projspps <- wdpa_by_spp_projnow[, sort(unique(as.character(spp)))]

# get tables of landings by LME from Sea Around Us
ebs <- fread('dataDL/sau/SAU LME 1 v47-1.csv')
ai <- fread('dataDL/sau/SAU LME 65 v47-1.csv') # will combine with ebs
goabc <- fread('dataDL/sau/SAU LME 2 v47-1.csv') # shared with goa and bc
wc <- fread('dataDL/sau/SAU LME 3 v47-1.csv')
gmex <- fread('dataDL/sau/SAU LME 5 v47-1.csv')
seus <- fread('dataDL/sau/SAU LME 6 v47-1.csv')
neus <- fread('dataDL/sau/SAU LME 7 v47-1.csv')
maritime <- fread('dataDL/sau/SAU LME 8 v47-1.csv')
newf <- fread('dataDL/sau/SAU LME 9 v47-1.csv')

# add region
ebs$region <- 'ebs'
ai$region <- 'ebs'
goabc$region <- 'goa_bc'
wc$region <- 'wc'
gmex$region <- 'gmex'
seus$region <- 'seus'
neus$region <- 'neus'
maritime$region <- 'maritime'
newf$region <- 'newf'

# combine and aggregate recent reported commercial landings by species
sau <- rbind(ebs, ai, goabc, wc, gmex, seus, neus, maritime, newf)
sppston <- sau[year >= 1995 & fishing_sector == 'Industrial' & catch_type == 'Landings' & reporting_status == 'Reported', .(tonnes = sum(tonnes)), by = c('region', 'scientific_name')]

# turn SAU names to lower case to match projections
sppston[, scientific_name := tolower(scientific_name)]

# set up vector hold names that match projections
sppston[, projname := as.character(NA)]

# match names from SAU to names from projections
options(warn=1)
for(j in 1:nrow(sppston)){
	ind <- agrep(sppston[j, scientific_name], projspps) # find index into spps
	if(length(ind)==1) sppston[j, projname := projspps[ind]] # only enter if exactly one match
}

# order by decreasing landings
setkey(sppston, region, tonnes)

# manually fix any that agrep missed or filled in mistakenly (checked by eye)
sppston[scientific_name == 'clupea pallasii pallasii', projname := 'clupea pallasii']


# rank the remaining species and remove those that aren't in the projections
sppstontrim <- sppston[!is.na(projname), .(sau_name = scientific_name, projname = projname, tonnes = tonnes, rank = rank(tonnes)), by = 'region']


# trim to top 10 that are in projections for each region
sppstontrim[, maxrank := max(rank), by='region']
sppstontrim <- sppstontrim[(maxrank - rank) < 10, ]
sppstontrim[ , c('rank', 'maxrank') := NULL] # remove unneeded columns

# check that it worked
sppstontrim[, .N, by='region']

# write out
write.csv(sppstontrim, file='output/fishery_spps.csv')
