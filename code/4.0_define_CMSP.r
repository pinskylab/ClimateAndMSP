# Define planning regions and choose which species count as important to fisheries


###########################
## Define planning regions
###########################
require(data.table)
require(sf)

# read in data
grid <- readRDS('temp/SPsf2.rds') # analysis grid
eez <- st_read('dataDL/marineregions/World_EEZ_v10_20180221/eez_v10.shp')
# coast <- st_read('dataDL/natcap/Marine/Land/global_polyline.shp') # useful for checking success by making maps

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
with(grid[grid$region == 'goa', ], plot(lon, lat, cex=0.1)) # odd chunk out of SW corner?
with(grid[grid$region == 'ebs', ], plot(lon, lat, cex=0.1))

# write out
saveRDS(grid, 'output/region_grid.rds')

###########################################################
# Determine which species are of fishery importance
# from Sea Around Us
###########################################################
## get list of projected species
load(paste('data/presmap_', runtype, projtype, '_rcp', rcp, '.RData', sep='')) # loads biomassavemap data.frame
head(presmap)
spps <- sort(unique(as.character(presmap$sppocean)))
sppstrim <- gsub('_Atl|_Pac', '', spps)

# get tables of landings by LME from Sea Around Us
ebs <- read.csv('cmsp_data/SAU/SAU LME 1 v1-40.csv')
al <- read.csv('cmsp_data/SAU/SAU LME 65 v1-40.csv')
goa <- read.csv('cmsp_data/SAU/SAU LME 2 v1-40.csv')
wcann <- read.csv('cmsp_data/SAU/SAU EEZ 848 v1-40.csv')
wctri <- wcann # west coast annual and trienniel: same LME
gmex <- read.csv('cmsp_data/SAU/SAU LME 5 v1-40.csv')
neuss <- read.csv('cmsp_data/SAU/SAU LME 7 v1-40.csv')
neusf <- neuss # NEUS spring and fall: same LME
scot <- read.csv('cmsp_data/SAU/SAU LME 8 v1-40.csv')
sgulf <- scot # southern Gulf of St. Lawrence, same LME as Scotian Shelf 
newff <- read.csv('cmsp_data/SAU/SAU LME 9 v1-40.csv')
newfs <- newff # newfoundland spring and fall: same LME

# add region
ebs$region <- 'AFSC_EBS'
al$region <- 'AFSC_Aleutians'
goa$region <- 'AFSC_GOA'
wcann$region <- 'NWFSC_WCAnn'
wctri$region <- 'AFSC_WCTri'
gmex$region <- 'SEFSC_GOMex'
neuss$region <- 'NEFSC_NEUSSpring'
neusf$region <- 'NEFSC_NEUSFall'
scot$region <- 'DFO_ScotianShelf'
sgulf$region <- 'DFO_SoGulf'
newff$region <- 'DFO_NewfoundlandFall'
newfs$region <- 'DFO_NewfoundlandSpring'

# combine into a list
sau <- list(ebs=ebs, al=al, goa=goa, wcann=wcann, wctri=wctri, gmex=gmex, neuss=neuss, neusf=neusf, scot=scot, sgulf=sgulf, newff=newff, newfs=newfs)

# aggregate recent reported landings by species
sppston <- sau
for(i in 1:length(sau)){
	sppston[[i]] <- with(sau[[i]][sau[[i]]$reporting_status=='Reported' & sau[[i]]$catch_type=='Landings' & sau[[i]]$year>=1990,], aggregate(list(tonnes = tonnes), by=list(region=region, scientific_name = scientific_name, common_name=common_name), FUN=sum))
}

# turn SAU names to lower case to match projections
for(i in 1:length(sppston)){
	sppston[[i]]$scientific_name <- tolower(sppston[[i]]$scientific_name)
}

# set up vector hold names that match projections
for(i in 1:length(sppston)){
	sppston[[i]]$projname <- NA
}


# match names from SAU to names from projections
options(warn=1)
for(i in 1:length(sppston)){
	for(j in 1:length(sppston[[i]]$scientific_name)){
		ind <- agrep(sppston[[i]]$scientific_name[j], sppstrim) # find index into spps and sppstrim
		if(length(ind)==1) sppston[[i]]$projname[j] <- spps[ind] # only enter if exactly one match
	}
}


# order by decreasing landings
for(i in 1:length(sppston)){
	sppston[[i]] <- sppston[[i]][order(sppston[[i]]$tonnes, decreasing=TRUE),]
}

# manually fix a few that agrep missed or filled in mistakenly
for(i in 1:length(sppston)){
	sppston[[i]]$projname[sppston[[i]]$scientific_name == 'clupea pallasii pallasii' & sppston[[i]]$region %in% atlregs] <- 'clupea pallasii_Atl'
	sppston[[i]]$projname[sppston[[i]]$scientific_name == 'clupea pallasii pallasii' & sppston[[i]]$region %in% pacregs] <- 'clupea pallasii_Pac'
	sppston[[i]]$projname[sppston[[i]]$scientific_name == 'sebastes'] <- NA
	sppston[[i]]$projname[sppston[[i]]$scientific_name == 'doryteuthis pealeii' & sppston[[i]]$region %in% atlregs] <- 'loligo pealeii_Atl'
	sppston[[i]]$projname[sppston[[i]]$scientific_name == 'doryteuthis pealeii' & sppston[[i]]$region %in% pacregs] <- 'loligo pealeii_Pac'
}

# top 10 that are in projections for each region
sppston[[1]][1:(which(cumsum(!is.na(sppston[[1]]$projname))==10)[1]),]
sppston[[2]][1:(which(cumsum(!is.na(sppston[[2]]$projname))==10)[1]),]
sppston[[3]][1:(which(cumsum(!is.na(sppston[[3]]$projname))==10)[1]),]
sppston[[4]][1:(which(cumsum(!is.na(sppston[[4]]$projname))==10)[1]),]
sppston[[5]][1:(which(cumsum(!is.na(sppston[[5]]$projname))==10)[1]),]
sppston[[6]][1:(which(cumsum(!is.na(sppston[[6]]$projname))==10)[1]),]
sppston[[7]][1:(which(cumsum(!is.na(sppston[[7]]$projname))==10)[1]),]
sppston[[8]][1:(which(cumsum(!is.na(sppston[[8]]$projname))==10)[1]),]
sppston[[9]][1:(which(cumsum(!is.na(sppston[[9]]$projname))==10)[1]),]
sppston[[10]][1:(which(cumsum(!is.na(sppston[[10]]$projname))==10)[1]),]
sppston[[11]][1:(which(cumsum(!is.na(sppston[[11]]$projname))==10)[1]),]
sppston[[12]][1:(which(cumsum(!is.na(sppston[[12]]$projname))==10)[1]),]

# combine into table to write out
temp <- sppston[[1]][1:(which(cumsum(!is.na(sppston[[1]]$projname))==10)[1]),]
temp <- temp[!is.na(temp$projname),]
out <- temp
for(i in 2:length(sppston)){
	temp <- sppston[[i]][1:(which(cumsum(!is.na(sppston[[i]]$projname))==10)[1]),]
	temp <- temp[!is.na(temp$projname),]
	out <- rbind(out,temp)
}

# write out
write.csv(out, file='cmsp_data/fishery_spps.csv')