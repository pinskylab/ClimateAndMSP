### Canadian Grain Size
setwd("~/Documents/Collaborations/Rutgers/Sediment Grain")

library(data.table)
library(gstat)
library(sp)
library(rgdal)
library(raster)

proj_latlon <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"


# =====================================
# = Folk Classification to Grain Size =
# =====================================
folk.phi <- readRDS("folk_usSEABED_perc_phi.rds")


# ====================
# = Shore Shapefiles =
# ====================
shore <- readOGR("canvec_1M_CA_Land_shp", "shoreline_1")
global_shore <- readOGR("gshhg-shp-2.3.5-1/GSHHS_shp/l", "GSHHS_l_L1") #low resolution of continental land masses from GSHHS



# ==========================
# = Canadian Sediment Data =
# ==========================
### Source: Expedition Database http://ed.gdr.nrcan.gc.ca/grainsize_detail_e.php
### Downloaded between -80 and 0 longitude and 40 and 90 latitude
### Mean= Grain Size
### X.Gravel = % Gravel
exped <- as.data.table(read.csv("CanadaExpedition/CanadaExpeditionDatabase.csv"))

### Grand Banks Sediment Thickness Database
### Source: http://geogratis.gc.ca/api/en/nrcan-rncan/ess-sst/97fc16ab-aadc-52f5-b33a-9145a78dd21c.html
gbanks<- readOGR("GrandBanks_SedThick_of_7513/ArcMapProject/Shapefiles", "GrandBanks_IsopachClasses")
gbanks.df <- as.data.frame(gbanks)
gbanks.data <- as.data.table(gbanks@data)
gbanks.data[,"Poly":=as.numeric(rownames(gbanks.df))]

gbanks.sed.type <- gbanks.data[,list(num.shape=length(unique(Shape_Area))), by=list(DominantSe)]
# write.csv(gbanks.sed.type, "GrandBanks_SedThick_of_7513/gbanks.sed.type.csv")

gbanks.sed.complete <- as.data.table(read.csv("GrandBanks_SedThick_of_7513/gbanks.sed.type.complete.csv"))

gbanks.sed.phi <- merge(gbanks.sed.complete, folk.phi, by="FOLKCODE")
gbanks.sed.phi.poly <- merge(gbanks.data, gbanks.sed.phi, by="DominantSe", all.x=T)
setorder(gbanks.sed.phi.poly, Poly)

gbanks@data <- gbanks.sed.phi.poly


dev.new()
spplot(gbanks, "GRAINSIZE")
spplot(gbanks, "MUD")
spplot(gbanks, "SAND")
spplot(gbanks, "GRAVEL")
spplot(gbanks, "FOLKCODE")

gbanks_latlon <- spTransform(gbanks, proj_latlon)
gbanks.coords <- coordinates(gbanks_latlon)

# writeOGR(obj=gbanks_latlon, dsn="GrandBanks_SedThick_of_7513", layer="gbanks_sed_latlon", driver="ESRI Shapefile")




# ### Laurentian Channel and SW Grand Banks sediments
# laur <- readOGR("Laurentian_Sed_of_6451", "LaurentianGeology")
# laur.data <- as.data.table(laur@data)
# laur.sed.types <- laur.data[,list(num.poly=length(unique(OBJECTID))), by=list(CODE)]
# setorder(laur.sed.types, num.poly)
# write.csv(laur.sed.types, "Laurentian_Sed_of_6451/laur.sed.types.csv")

# ==========================
# = Prep for Interpolation =
# ==========================

bathy <- raster("ETOPO2v2c_f4.nc")
canada.extent <- extent(-68,-40, 40, 65)
bathy.canada <- crop(bathy, canada.extent)

plot(bathy.canada, zlim=c(-1500, 0))
points(Latitude ~ Longitude, exped, col="gray", pch=20)


load("haulInfo_Dec26_2016.RData")
hauls.dt <- as.data.table(hauls)
hauls.newf <- hauls.dt[region %in% c("DFO_NewfoundlandFall","DFO_NewfoundlandSpring") & lat <=52]
hauls.lab <- hauls.dt[region %in% c("DFO_NewfoundlandFall","DFO_NewfoundlandSpring") & lat >52]
hauls.sogulf <- hauls.dt[region %in% c("DFO_SoGulf")]

hauls.scot <- hauls.dt[region %in% c("DFO_ScotianShelfSummer","DFO_ScotianShelfFall", "DFO_ScotianShelfSpring")]



plot(Latitude ~ Longitude, exped, col="gray", pch=20)
plot(global_shore, add=T)
points(lat ~ lon, hauls.newf, col="red", pch=1, cex=0.2) #pierre bank not surveyed
points(lat ~ lon, hauls.lab, col="darkgreen", pch=1, cex=0.2)
points(lat ~ lon, hauls.scot, col="blue", pch=1, cex=0.2)
points(lat ~ lon, hauls.sogulf, col="orange", pch=1, cex=0.2) # Southern Gulf of St. Lawrence not well covered
points(Latitude ~ Longitude, exped, col="gray", pch=1, cex=0.2)
points(gbanks.coords[,1], gbanks.coords[,2], col="yellow")

# ==========================
# = Interpolate grain size =
# ==========================
coordinates(exped) <- cbind(exped$Longitude, exped$Latitude)
proj4string(exped) <- proj_latlon
# writeOGR(obj=exped, dsn="CanadaExpedition", layer="exped_sp_pts", driver="ESRI Shapefile")


bathy.dt <- as.data.table(as.data.frame(bathy.canada, xy=T))
coordinates(bathy.dt) <- ~ x + y
proj4string(bathy.dt) <- proj_latlon

bathy.grain <- idw(Mean~1, exped, subset(bathy.dt, z <0 & z > - 1500))
bathy.sand <- idw(X.Sand~1, exped, subset(bathy.dt, z <0 & z > - 1500))
bathy.gravel <- idw(X.Gravel~1, exped, subset(bathy.dt, z <0 & z > - 1500))
bathy.mud <- idw(X.Mud~1, exped, subset(bathy.dt, z <0 & z > - 1500))



bathy.canada.df <- as.data.frame(bathy.grain)
bathy.canada.df$GRAINSIZE <- bathy.grain$var1.pred
bathy.canada.df$GRAVEL <- bathy.gravel$var1.pred
bathy.canada.df$SAND <- bathy.sand$var1.pred
bathy.canada.df$MUD <- bathy.mud$var1.pred


bathy.grain.raster <- rasterFromXYZ(bathy.canada.df[,c("x", "y", "GRAINSIZE")])
# plot(bathy.grain.raster, main="GRAINSIZE")
bathy.gravel.raster <- rasterFromXYZ(bathy.canada.df[,c("x", "y", "GRAVEL")])
# plot(bathy.gravel.raster, main="GRAVEL")
bathy.sand.raster <- rasterFromXYZ(bathy.canada.df[,c("x", "y", "SAND")])
# plot(bathy.sand.raster, main="SAND")
bathy.mud.raster <- rasterFromXYZ(bathy.canada.df[,c("x", "y", "MUD")])
# plot(bathy.mud.raster, main="MUD")

canada.raster.stack <- stack(bathy.grain.raster, bathy.gravel.raster, bathy.sand.raster, bathy.mud.raster)
# names(canada.raster.stack) <- c("GRAINSIZE", "GRAVEL","SAND",  "MUD")
plot(canada.raster.stack)
# saveRDS(canada.raster.stack, "sed_grid_out/canada.sed.stack.out.rds")

# ==================================
# = Interpolate for haul locations =
# ==================================
### Scotian
coordinates(hauls.scot) <- cbind(hauls.scot$lon, hauls.scot$lat)
proj4string(hauls.scot) <- proj_latlon
hauls.scot.crop <- crop(hauls.scot, extent(-67.7,-57, 41,47))
# plot(lat ~ lon, hauls.scot.crop)
# points(Latitude ~ Longitude, exped, col="red")
# abline(v=-67.7)

scotian.grain <- idw(Mean~1, exped, hauls.scot.crop)
scotian.sand <- idw(X.Sand~1, exped, hauls.scot.crop)
scotian.gravel <- idw(X.Gravel~1, exped, hauls.scot.crop)
scotian.mud <- idw(X.Mud~1, exped, hauls.scot.crop)

hauls.scot.out <- as.data.table(as.data.frame(hauls.scot.crop))
hauls.scot.out[,"GRAINSIZE":=scotian.grain$var1.pred]
hauls.scot.out[,"GRAVEL":=scotian.gravel$var1.pred]
hauls.scot.out[,"SAND":=scotian.sand$var1.pred]
hauls.scot.out[,"MUD":=scotian.mud$var1.pred]

saveRDS(hauls.scot.out, "haul_sp/hauls.scot.idw.rds")


hauls.scot.out.sp <- copy(hauls.scot.out)
coordinates(hauls.scot.out.sp) <- cbind(hauls.scot.out.sp$coords.x1, hauls.scot.out.sp$coords.x2)
proj4string(hauls.scot.out.sp) <- proj_latlon
writeOGR(hauls.scot.out.sp, "haul_sp", "haul_scot_out_sp", driver="ESRI Shapefile")

### Labrador shelf is reasonably well sampled (north of 52 degrees latitude within newfoundland survey)
coordinates(hauls.lab) <- cbind(hauls.lab$lon, hauls.lab$lat)
proj4string(hauls.lab) <- proj_latlon

lab.grain <- idw(Mean~1, exped, hauls.lab)
lab.sand <- idw(X.Sand~1, exped, hauls.lab)
lab.gravel <- idw(X.Gravel~1, exped, hauls.lab)
lab.mud <- idw(X.Mud~1, exped, hauls.lab)

hauls.lab.out <- as.data.table(as.data.frame(hauls.lab))
hauls.lab.out[,"GRAINSIZE":=lab.grain$var1.pred]
hauls.lab.out[,"GRAVEL":=lab.gravel$var1.pred]
hauls.lab.out[,"SAND":=lab.sand$var1.pred]
hauls.lab.out[,"MUD":=lab.mud$var1.pred]

saveRDS(hauls.lab.out, "haul_sp/hauls.lab.idw.rds")


hauls.lab.out.sp <- copy(hauls.lab.out)
coordinates(hauls.lab.out.sp) <- cbind(hauls.lab.out.sp$coords.x1, hauls.lab.out.sp$coords.x2)
proj4string(hauls.lab.out.sp) <- proj_latlon
writeOGR(hauls.lab.out.sp, "haul_sp", "haul_lab_out_sp", driver="ESRI Shapefile")

# =======================================================
# = Extract info for hauls for Grand Banks from polygon =
# =======================================================
coordinates(hauls.newf) <- cbind(hauls.newf$lon, hauls.newf$lat)
proj4string(hauls.newf) <- proj_latlon
hauls.newf.poly <- over(hauls.newf, gbanks_latlon)

hauls.newf.out <- as.data.table(as.data.frame(hauls.newf))
hauls.newf.out[,"GRAINSIZE":=hauls.newf.poly$GRAINSIZE]
hauls.newf.out[,"GRAVEL":=hauls.newf.poly$GRAVEL]
hauls.newf.out[,"SAND":=hauls.newf.poly$SAND]
hauls.newf.out[,"MUD":=hauls.newf.poly$MUD]

hauls.newf.GB.out <- hauls.newf.out[!(is.na(GRAINSIZE))]
saveRDS(hauls.newf.GB.out, "haul_sp/hauls.newf.poly.rds")


hauls.newf.GB.out.sp <- copy(hauls.newf.GB.out)
coordinates(hauls.newf.GB.out.sp) <- cbind(hauls.newf.GB.out.sp$coords.x1, hauls.newf.GB.out.sp$coords.x2)
proj4string(hauls.newf.GB.out.sp) <- proj_latlon
writeOGR(hauls.newf.GB.out.sp, "haul_sp", "haul_newf_GB_out_sp", driver="ESRI Shapefile")

### For hauls in newfoundland not in the GrandBanks use idw
hauls.newf.nonGB <- hauls.newf.out[is.na(GRAINSIZE)]
coordinates(hauls.newf.nonGB) <- cbind(hauls.newf.nonGB$coords.x1, hauls.newf.nonGB$coords.x2)
proj4string(hauls.newf.nonGB) <- proj_latlon

newf.nonGB.grain <- idw(Mean~1, exped, hauls.newf.nonGB)
newf.nonGB.gravel <- idw(X.Sand~1, exped, hauls.newf.nonGB)
newf.nonGB.sand <- idw(X.Gravel~1, exped, hauls.newf.nonGB)
newf.nonGB.mud <- idw(X.Mud~1, exped, hauls.newf.nonGB)

hauls.newf.nonGB.out <- as.data.table(as.data.frame(hauls.newf.nonGB))
hauls.newf.nonGB.out[,"GRAINSIZE":=newf.nonGB.grain$var1.pred]
hauls.newf.nonGB.out[,"GRAVEL":=newf.nonGB.gravel$var1.pred]
hauls.newf.nonGB.out[,"SAND":=newf.nonGB.sand$var1.pred]
hauls.newf.nonGB.out[,"MUD":=newf.nonGB.mud$var1.pred]
hauls.newf.nonGB.out[,"coords.x1.1":=NULL]
hauls.newf.nonGB.out[,"coords.x2.1":=NULL]

saveRDS(hauls.newf.nonGB.out, "haul_sp/hauls.newf.idw.rds")


hauls.newf.nonGB.out.sp <- copy(hauls.newf.nonGB.out)
coordinates(hauls.newf.nonGB.out.sp) <- cbind(hauls.newf.nonGB.out.sp$coords.x1, hauls.newf.nonGB.out.sp$coords.x2)
proj4string(hauls.newf.nonGB.out.sp) <- proj_latlon
writeOGR(hauls.newf.nonGB.out.sp, "haul_sp", "haul_newf_nonGB_out_sp", driver="ESRI Shapefile")

