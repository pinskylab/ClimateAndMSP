### Alaska Sediments
setwd("~/Documents/Collaborations/Rutgers/Sediment Grain")

library(data.table)
library(rgdal)
library(sp)
library(raster)
library(gstat)

bathy <- raster("ETOPO2v2c_f4.nc")
proj_latlon <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"



goa <- as.data.table(read.csv("GulfofAlaskaDig/GulfofAlaskaDigitizationProject_NOSSeafloorCharacter.csv"))
ebs <- as.data.table(read.csv("ebssed/EBSSED.csv"))
ebs[,"Longitude":=LONGITUDE.N.20.5]
ebs[,"Latitude":=LATITUDE.N.20.5]

ai <- readOGR("AI_sediments", "Aleutian_sediments")
ai.data <- as.data.table(ai@data)

ai.data[Longitude >0, min(Longitude)]

proj_NAD83 <- proj4string(ai)



# =====================================
# = Folk Classification to Grain Size =
# =====================================
folk.phi <- readRDS("folk_usSEABED_perc_phi.rds")

# ====================================================
# = Estimate phi for sediment classes for AI and GOA =
# ====================================================
ai.sed.codes <- as.data.table(read.csv("AI_sediments/ai.sed.codes.complete.csv"))
goa.sed.codes <- as.data.table(read.csv("GulfofAlaskaDig/goa.sed.codes.complete.csv"))

### Merge with phi estimates for each sediment code
ai.sed.codes2 <- merge(ai.sed.codes, folk.phi, by="FOLKCODE")
goa.sed.codes2 <- merge(goa.sed.codes, folk.phi, by="FOLKCODE")



### Merge with data
ai.data2 <- merge(ai.data, ai.sed.codes2, by="Sediment")
goa2 <- merge(goa, goa.sed.codes2, by=c("NOSSedCode", "Definition"))


# ============================================
# = Create master set for all of Alaska data =
# ============================================
goa.master <- goa2[,list(Longitude=Longitude, Latitude=Latitude, FOLKCODE=FOLKCODE, GRAVEL=GRAVEL, SAND=SAND, MUD=MUD, GRAINSIZE=GRAINSIZE)]
ai.master <- ai.data2[,list(Longitude=Longitude, Latitude=Latitude, FOLKCODE=FOLKCODE, GRAVEL=GRAVEL, SAND=SAND, MUD=MUD, GRAINSIZE=GRAINSIZE)]
ebs.master <- ebs[,list(Longitude=Longitude, Latitude=Latitude, FOLKCODE=HIGH_RES__.C.5, 
	GRAVEL=X__GRAVEL.N.20.5, SAND=X__SAND.N.20.5, MUD=X__MUD.N.20.5, GRAINSIZE=MEAN_PHI.N.20.5)]

ak.master <- rbind(goa.master, ai.master, ebs.master)
ak.E <- ak.master[Longitude<0, list(min.Lon=-180, max.Lon=max(Longitude), min.Lat=min(Latitude), max.Lat=max(Latitude))]
ak.W <- ak.master[Longitude>0, list(min.Lon=170, max.Lon=180, min.Lat=min(Latitude), max.Lat=max(Latitude))]

ak.E.e <- extent(c(ak.E$min.Lon, ak.E$max.Lon, ak.E$min.Lat, ak.E$max.Lat))
ak.W.e <- extent(c(ak.W$min.Lon, ak.W$max.Lon, ak.W$min.Lat, ak.W$max.Lat))

bathy.E <- crop(bathy, ak.E.e)
bathy.W <- crop(bathy, ak.W.e)

bathy.mosaic <- mosaic(x=bathy.E, y=bathy.W, fun=mean) #won't be any overlap, but have to provide a function in case there is

### Make sediment data into spatial pixels data frame
coordinates(ak.master) <- cbind(ak.master$Longitude, ak.master$Latitude)
proj4string(ak.master) <- proj_NAD83

# ==================================
# = Make bathy into spatial pixels =
# ==================================
### Convert bathy into spatial pixels data frame
raster2sed_mosaic <- function(r, sed){
	bathy.atl.df <- as.data.frame(r, xy=T)

	### Get in same format as sediment data
	bathy.atl.df$WATERDEPTH <- ifelse(bathy.atl.df$layer<0 & bathy.atl.df$layer>=-3000, -1*bathy.atl.df$layer, NA)
	bathy.atl.df$Longitude <- bathy.atl.df$x
	bathy.atl.df$Latitude <- bathy.atl.df$y
	bathy.atl.df$z <- NULL
	bathy.atl.df$x <- NULL
	bathy.atl.df$y <- NULL

	coordinates(bathy.atl.df) <- ~ Longitude + Latitude
	bathy.atl.df@proj4string <- sed@proj4string
	return(bathy.atl.df)
}

bathy.sed <- raster2sed_mosaic(bathy.mosaic, ak.master)

# =============================================
# = Inverse Distance Weighting for Grain Size =
# =============================================
ak.grain.idw <- idw(GRAINSIZE~1, subset(ak.master, !(is.na(GRAINSIZE))), subset(bathy.sed, WATERDEPTH>0 & WATERDEPTH<750))
# saveRDS(ak.grain.idw, "ak.combined.grain.idw.rds")

spplot(subset(ak.grain.idw, Longitude<0), "var1.pred")

ak.grain.df <- as.data.frame(ak.grain.idw)
ak.grain.raster <- rasterFromXYZ(ak.grain.df[,c("Longitude", "Latitude", "var1.pred")])
plot(ak.grain.raster, xlim=c(-180, -130), main="Grain Size")
# ==============================================
# = Inverse Distance Weighting for Percentages =
# ==============================================
ak.gravel.idw <- idw(GRAVEL~1, subset(ak.master, !(is.na(GRAVEL))), subset(bathy.sed, WATERDEPTH>0 & WATERDEPTH<750))
# saveRDS(ak.gravel.idw, "ak.combined.gravel.idw.rds")

ak.gravel.df <- as.data.frame(ak.gravel.idw)
ak.gravel.raster <- rasterFromXYZ(ak.gravel.df[,c("Longitude", "Latitude", "var1.pred")])
plot(ak.gravel.raster, xlim=c(-180, -130), main="Gravel")

ak.mud.idw <- idw(MUD~1, subset(ak.master, !(is.na(MUD))), subset(bathy.sed, WATERDEPTH>0 & WATERDEPTH<750))
# saveRDS(ak.mud.idw, "ak.combined.mud.idw.rds")

ak.mud.df <- as.data.frame(ak.mud.idw)
ak.mud.raster <- rasterFromXYZ(ak.mud.df[,c("Longitude", "Latitude", "var1.pred")])
plot(ak.mud.raster, xlim=c(-180, -130), main="Mud")

ak.sand.idw <- idw(SAND~1, subset(ak.master, !(is.na(SAND))), subset(bathy.sed, WATERDEPTH>0 & WATERDEPTH<750))
# saveRDS(ak.mud.idw, "ak.combined.sand.idw.rds")

ak.sand.df <- as.data.frame(ak.sand.idw)
ak.sand.raster <- rasterFromXYZ(ak.sand.df[,c("Longitude", "Latitude", "var1.pred")])
plot(ak.sand.raster, xlim=c(-180, -130), main="Sand")

par(mfrow=c(2,2))
plot(ak.grain.raster, xlim=c(-180, -130), main="Grain Size")
plot(ak.gravel.raster, xlim=c(-180, -130), main="Gravel")
plot(ak.mud.raster, xlim=c(-180, -130), main="Mud")
plot(ak.sand.raster, xlim=c(-180, -130), main="Sand")



ak.sed.stack <- stack(ak.grain.raster, ak.gravel.raster, ak.sand.raster, ak.mud.raster)
names(ak.sed.stack) <- c("GRAINSIZE", "GRAVEL","SAND",  "MUD")
plot(ak.sed.stack, xlim=c(-180, -130), ylim=c(50,70))
# saveRDS(ak.sed.stack, "sed_grid_out/ak.sed.stack.out.rds")



# ==========================
# = Extract data for hauls =
# ==========================
load("haulInfo_Dec26_2016.RData")
haul.dt <- as.data.table(subset(hauls, !(is.na(lat))& !(is.na(lon)))) # 69 hauls without lat lon
### make longitudes < -180 equal to 180 - (180+lon)
haul.dt[,"lon.adj":=ifelse(lon< -180, 180+(180+lon), lon)]

ak.hauls <- haul.dt[region %in% c("AFSC_EBS", "AFSC_Aleutians", "AFSC_GOA")]

coordinates(ak.hauls) <- cbind(ak.hauls$lon.adj, ak.hauls$lat)
proj4string(ak.hauls) <- proj_NAD83


ak.hauls.grain.idw <- idw(GRAINSIZE~1, subset(ak.master, !(is.na(GRAINSIZE))), ak.hauls)
ak.hauls.gravel.idw <- idw(GRAVEL~1, subset(ak.master, !(is.na(GRAVEL))), ak.hauls)
ak.hauls.sand.idw <- idw(SAND~1, subset(ak.master, !(is.na(SAND))), ak.hauls)
ak.hauls.mud.idw <- idw(MUD~1, subset(ak.master, !(is.na(MUD))), ak.hauls)


spplot(ak.hauls.grain.idw, "var1.pred")

ak.hauls.dt <- as.data.table(as.data.frame(ak.hauls))
ak.hauls.dt[,"GRAINSIZE":=ak.hauls.grain.idw$var1.pred]
ak.hauls.dt[,"GRAVEL":=ak.hauls.gravel.idw$var1.pred]
ak.hauls.dt[,"SAND":=ak.hauls.sand.idw$var1.pred]
ak.hauls.dt[,"MUD":=ak.hauls.mud.idw$var1.pred]

saveRDS(ak.hauls.dt, "hauls.ak.idw.rds")

ak.hauls.dt <- readRDS("haul_sp/ak.hauls.idw.rds")
ak.hauls.sp <- copy(ak.hauls.dt)
coordinates(ak.hauls.sp) <- cbind(ak.hauls.sp$coords.x1, ak.hauls.sp$coords.x2)
proj4string(ak.hauls.sp) <- proj_NAD83
writeOGR(ak.hauls.sp, "haul_sp", "haul_ak_out_sp", driver="ESRI Shapefile")



