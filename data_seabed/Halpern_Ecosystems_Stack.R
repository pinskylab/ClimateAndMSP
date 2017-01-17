setwd("~/Documents/Collaborations/Rutgers/Sediment Grain/Halpern2008_Data/")

library(data.table)
library(raster)
library(sp)
library(rgdal)
library(rasterVis)

### The original projection of the grid (discovered by opening prj.adf as a text file)
### World Mollweide projection (ESRI:54009)
proj_moll <- "+proj=moll +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"
proj_latlon <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"


### Deep Soft Bottom Presence
deep_soft <- raster("grid/d_s_bottom/w001001.adf")
### Deep Hard Bottom
deep_hard <- raster("grid/d_h_bottom/w001001.adf")

### Shelf Soft and Hard Bottom
shelf_soft <-  raster("grid/soft_shelf/w001001.adf")
shelf_hard <-  raster("grid/hard_shelf/w001001.adf")

### Slope Soft and Hard Bottom
slope_soft <-  raster("grid/soft_slope/w001001.adf")
slope_hard <-  raster("grid/hard_slope/w001001.adf")

### Rocky Reef
rocky_reef <- raster("grid/rky_reef/w001001.adf")
### Subtidal soft bottom
subt_soft <- raster("grid/s_t_s_bottom/w001001.adf")


crs(deep_soft) <- proj_moll
crs(deep_hard) <- proj_moll

crs(shelf_soft) <- proj_moll
crs(shelf_hard) <- proj_moll

crs(slope_soft) <- proj_moll
crs(slope_hard) <- proj_moll

crs(rocky_reef) <- proj_moll
crs(subt_soft) <- proj_moll

### Raster stack
hab <- stack(list(deep_soft, slope_soft,  shelf_soft, subt_soft, deep_hard, shelf_hard,  slope_hard, rocky_reef))
names(hab) <- c("deep_soft", "slope_soft",  "shelf_soft", "subt_soft", "deep_hard", "shelf_hard",  "slope_hard", "rocky_reef")

rm(deep_soft, deep_hard, shelf_soft, shelf_hard, slope_soft, slope_hard, rocky_reef, subt_soft)

# ======================
# = Get Haul Locations =
# ======================
### Will substitute with the RData file Jim sends me
# library(trawlData)
# haul.neus.raw <- clean.neus[,list(lat=mean(lat), lon=mean(lon)), by=list(haulid)]
#
# haul.neus <- clean.neus[,list(lat=mean(lat), lon=mean(lon)), by=list(haulid)]

load("haulInfo_Dec26_2016.RData")

haul.dt <- as.data.table(subset(hauls, !(is.na(lat))& !(is.na(lon)))) # 69 hauls without lat lon
### make longitudes < -180 equal to 180 - (180+lon)
haul.dt[,"lon.adj":=ifelse(lon< -180, 180+(180+lon), lon)]

coordinates(haul.dt) <- cbind(x=haul.dt$lon.adj, y=haul.dt$lat)
proj4string(haul.dt) <- proj_latlon


### Convert to World Mollweide projection (ESRI:54009)
haul_moll <- spTransform(haul.dt, CRS(proj_moll))


# =====================================
# = Classify Habitat at each location =
# =====================================

class_hab <- function(hab_stack, df_pts){
	hab_extract <- vector("list", dim(hab_stack)[3])
	dt <- data.table(as.data.frame(df_pts))
	hab_names <- names(hab_stack)
	
	for(i in 1:dim(hab)[3]){
		hab_extract[[i]] <- extract(hab_stack[[i]], df_pts, method="simple")
		name <- hab_names[i]
		dt[,(eval(name)):=hab_extract[[i]]]
	}
	dt_names <- names(dt)
	dt_melt <- melt(dt, id.vars=dt_names[!(dt_names%in%hab_names)], variable.name="habitat")
	dt_melt[,"hab_factor":=as.numeric(habitat)]

	### Subset to only values ==1 to get habitat assignment
	dt_hab_assign <- dt_melt[value==1]
	return(dt_hab_assign)
}

haul_hab_out <- class_hab(hab, haul_moll)

# neus.grid.out <- class_hab(hab, grid_moll)
#
# rug_hab_out <- class_hab(hab, rug_moll)

saveRDS(haul_hab_out, "habitat_at_haul_latlon.rds")

# ================
# = Soft vs Hard =
# ================
hab.dist <- haul_hab_out[,list(soft=sum(hab_factor %in% c(1:4)), hard=sum(hab_factor %in% c(5:8)))]
#136902 soft, 11028 hard

