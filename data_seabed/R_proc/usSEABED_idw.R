### usSEABED
setwd("~/Documents/Collaborations/Rutgers/Sediment Grain")

library(data.table)
library(rgdal)
library(raster)
library(sp)
library(gstat)

states <- readOGR(".", "states_latlon")

bathy <- raster("ETOPO2v2c_f4.nc")

# =================
# = Sediment Point Data =
# =================
### Spatial Points Dataframe from usSEABED
pac <- readOGR("pac_prs", "pac_prs")
atl <- readOGR("atl_prs", "atl_prs")
gmx <- readOGR("gmx_prs", "gmx_prs")





pac.data <- as.data.table(pac@data)
atl.data <- as.data.table(atl@data)
gmx.data <- as.data.table(gmx@data)
# atl_usgs.data <- as.data.table(atl_usgs@data)



### Extent of seabed data
pac.e <- extent(min(pac.data$LONGITUDE), max(pac.data$LONGITUDE), min(pac.data$LATITUDE), max(pac.data$LATITUDE))
atl.e <- extent(min(atl.data$LONGITUDE), max(atl.data$LONGITUDE), min(atl.data$LATITUDE), max(atl.data$LATITUDE))
gmx.e <- extent(min(gmx.data$LONGITUDE), max(gmx.data$LONGITUDE), min(gmx.data$LATITUDE), max(gmx.data$LATITUDE))

### Crop bathy to these extents
pac.bathy <- crop(bathy, pac.e)
atl.bathy <- crop(bathy, atl.e)
gmx.bathy <- crop(bathy, gmx.e)

plot(atl.bathy, zlim=c(-6000,0))

# ================================================
# = Percent Gravel, Mud, Sand for each Folk Code =
# ================================================
atl_ext <- readOGR("atl_ext", "atl_ext")
pac_ext <- readOGR("pac_ext", "pac_ext")
gmx_ext <- readOGR("gmx_ext", "gmx_ext")

atl_ext.data <- as.data.table(atl_ext@data)
pac_ext.data <- as.data.table(pac_ext@data)
gmx_ext.data <- as.data.table(gmx_ext@data)

usseabed_all <- rbind(atl_ext.data, pac_ext.data, gmx_ext.data)


folk.perc.atl <- atl_ext.data[GRAVEL != -99 & SAND != -99 & MUD !=-99 & GRAINSIZE !=-99,list(GRAVEL=mean(GRAVEL), SAND=mean(SAND), MUD=mean(MUD), GRAINSIZE=mean(GRAINSIZE)), by=list(FOLKCODE)]
folk.perc.pac <- pac_ext.data[GRAVEL != -99 & SAND != -99 & MUD !=-99 & GRAINSIZE !=-99,list(GRAVEL=mean(GRAVEL), SAND=mean(SAND), MUD=mean(MUD), GRAINSIZE=mean(GRAINSIZE)), by=list(FOLKCODE)]
folk.perc.gmx <- gmx_ext.data[GRAVEL != -99 & SAND != -99 & MUD !=-99 & GRAINSIZE !=-99,list(GRAVEL=mean(GRAVEL), SAND=mean(SAND), MUD=mean(MUD), GRAINSIZE=mean(GRAINSIZE)), by=list(FOLKCODE)]

folk.perc.us <- usseabed_all[GRAVEL != -99 & SAND != -99 & MUD !=-99 & GRAINSIZE !=-99,list(GRAVEL=mean(GRAVEL), SAND=mean(SAND), MUD=mean(MUD), GRAINSIZE=mean(GRAINSIZE)), by=list(FOLKCODE)]

# saveRDS(folk.perc.us, "folk_usSEABED_perc_phi.rds")

### phi estimates for Gravel, sand mud in order to match mean GRAINSIZE 
G.phi <- -2.5
S.phi <- 1.7
M.phi <- 7.5

folk.perc.us[,"Grain.Estimated":=0.01*GRAVEL*G.phi + 0.01*SAND*S.phi + 0.01*MUD*M.phi]

# =======================
# = Haul locations =
# =======================
load("haulInfo_Dec26_2016.RData")
haul.dt <- as.data.table(subset(hauls, !(is.na(lat))& !(is.na(lon)))) # 69 hauls without lat lon
### make longitudes < -180 equal to 180 - (180+lon)
haul.dt[,"lon.adj":=ifelse(lon< -180, 180+(180+lon), lon)]
coordinates(haul.dt) <- cbind(haul.dt$lon.adj, haul.dt$lat)

# writeOGR(haul.dt, "haul_sp", "haul_sp", driver="ESRI Shapefile")




# ============
# = Sediment Polygons =
# ============
### California Seafloor Mapping
### http://seafloor.otterlabs.org/SFMLwebDATA_SURVEYMAP.htm
# area1.hab <- readOGR("CAGeologySeries_Project_mxd/CAGeologySeries_Project_mxd/Area1", "gma1hab")
# area2.hab <- readOGR("CAGeologySeries_Project_mxd/CAGeologySeries_Project_mxd/Area2", "gma2hab")
# area3.hab <- readOGR("CAGeologySeries_Project_mxd/CAGeologySeries_Project_mxd/Area3", "GMA3HAB")
# area4.hab <- readOGR("CAGeologySeries_Project_mxd/CAGeologySeries_Project_mxd/Area4", "GMA4HAB")
# area5.hab <- readOGR("CAGeologySeries_Project_mxd/CAGeologySeries_Project_mxd/Area5", "GMA5HAB")
# area7.hab <- readOGR("CAGeologySeries_Project_mxd/CAGeologySeries_Project_mxd/Area7", "GMA7HAB")




### From dbSEABED interpolation from Gulf of Mexico Data Atlas
### https://www.ncddc.noaa.gov/website/DataAtlas/atlas.htm
### Spatial Polygons Data Frame
## v gridded values, computed with CS Interpolator with 5km search radius for rock, 20km for other sediment categories
## u gridded uncertainties
### Folk Code: position 1,2,3=mud,sand,gravel; value 0,1,2,3 indicate negligible-, slightly, x-ly, or major
### domnc: position 1,2,3,4=rock,gravel,sand,mud; value 0,2,3 indicate negligible-, subdominant-, dominant-
gmx_db <- readOGR("usSEABED_GOM_Sediments", "usSEABED_GOM_Sediments")
gmx_db$GRAINSIZE=0.01*G.phi*gmx_db$gom_gvlv+0.01*S.phi*gmx_db$gom_sndv + 0.01*M.phi*gmx_db$gom_mudv
# spplot(gmx_db)

### Atlantic Coast sediment types from CONMAP
## https://catalog.data.gov/dataset/conmapsg-continental-margin-mapping-program-conmap-sediments-grainsize-distribution-for-the-uni
### Spatial Polygons DataFrame
### Sediments
# 1=bedrock, 2=gravel, 3=gravel-sand, 4=sand, 5=clay-silt/sand, 6=sand-clay/silt, 7=clay, 8=sand-silt/clay, 9=sand/silt/clay
conmap <- readOGR("conmapsg", "conmapsg")
# spplot(conmap, "SEDIMENT")

### East Coast Sediment Texture Database
ecst <- readOGR("ecstdb2005", "ecstdb2005")

# ### Scotian Shelf Sediments
# plot(subset(ecst, LONGITUDE > -66))

# =============================
# = Get phi for CONMAP SEDNUM =
# =============================
ecst_sp <- spTransform(ecst, CRS(proj4string(conmap)))

ecst_poly <- over(ecst_sp, conmap)

ecst_data <- as.data.table(ecst@data)
ecst_data[,"SEDIMENT":=ecst_poly$SEDIMENT]
ecst_data[,"SEDNUM":=ecst_poly$SEDNUM]
ecst_data[,"MUD_PCT":=SILT_PCT + CLAY_PCT]

ecst_phi <- ecst_data[,list(GRAINSIZE=mean(MEAN, na.rm=T), GRAVEL=mean(GRAVEL_PCT, na.rm=T), SAND=mean(SAND_PCT, na.rm=T), MUD=mean(MUD_PCT, na.rm=T)), by=list(SEDIMENT, SEDNUM)]
setorder(ecst_phi, SEDNUM)
# saveRDS("ecst_phi", "ecst_conmap_phi.rds")

# =============================================
# = Get sediment data from polygons for hauls =
# =============================================

### Atlantic
proj4string(haul.dt) <- proj4string(conmap)
e.conmap <- extent(conmap)
haul.dt.atl <- crop(haul.dt, e.conmap)
haul.dt.atl2 <- subset(haul.dt.atl, !(region %in% c("SEFSC_GOMexSummer" ,"SEFSC_GOMexFall", "DFO_ScotianShelfSummer", "DFO_ScotianShelfFall"))) #DFO Scotian Shelf Spring does have hauls in conmap area
haul.dt.atl.poly <-over(haul.dt.atl2, conmap)

haul.dt.atl.out <- as.data.table(as.data.frame(haul.dt.atl2))
haul.dt.atl.out[,"SEDIMENT":=haul.dt.atl.poly$SEDIMENT]
haul.dt.atl.out[,"SEDNUM":=haul.dt.atl.poly$SEDNUM]

### Merge with ecst_phi to get GRAINSIZE and percentages
haul.dt.atl.out2 <- merge(haul.dt.atl.out, ecst_phi, by=c("SEDIMENT", "SEDNUM"))

# saveRDS(haul.dt.atl.out2, "haul_sp/hauls.atl.sed.conmap.rds")

haul.atl.sp <- as.data.frame(copy(haul.dt.atl.out2))
coordinates(haul.atl.sp) <- cbind(haul.atl.sp$coords.x1, haul.atl.sp$coords.x2)
proj4string(haul.atl.sp) <- proj4string(conmap)
writeOGR(haul.atl.sp, "haul_sp", "haul_atl_out_sp", driver="ESRI Shapefile")




### Gulf of Mexico
haul.dt.gmx <- spTransform(haul.dt, CRS(proj4string(gmx_db)))
e.gmx <- extent(gmx_db)
haul.dt.gmx2 <- crop(haul.dt.gmx, e.gmx)
haul.dt.gmx3 <- subset(haul.dt.gmx2, region%in% c("SEFSC_GOMexSummer" ,"SEFSC_GOMexFall"))

haul.dt.gmx.poly <- over(haul.dt.gmx3, gmx_db)
haul.dt.gmx.out <- as.data.table(as.data.frame(haul.dt.gmx3))
haul.dt.gmx.out[,"GRAINSIZE":=haul.dt.gmx.poly$GRAINSIZE]
haul.dt.gmx.out[,"GRAVEL":=haul.dt.gmx.poly$gom_gvlv]
haul.dt.gmx.out[,"SAND":=haul.dt.gmx.poly$gom_sndv]
haul.dt.gmx.out[,"MUD":=haul.dt.gmx.poly$gom_mudv]

# saveRDS(haul.dt.gmx.out, "haul_sp/hauls_sed_gmx_db.rds")

haul.gmx.sp <- as.data.frame(copy(haul.dt.gmx.out))
coordinates(haul.gmx.sp) <- cbind(haul.gmx.sp$coords.x1, haul.gmx.sp$coords.x2)
proj4string(haul.gmx.sp) <- proj4string(conmap)
writeOGR(haul.gmx.sp, "haul_sp", "haul_gmx_out_sp", driver="ESRI Shapefile")


# ### According to metadata for gmx_db
# # gma_folk - gridded codes for sediment Folk Codes; the 0 (or " "),1,2,3 indicate
# ### negligible-, slightly-, x-ly-, or major components of mud, sand or gravel in code positions 1,2,3 respectively;
# # they can be converted to "(x)xX" types codes;
# gmx.folk.class <- data.frame(FOLK_num=c(103,123, 130,132,203,230,232,300,302,320,322),
# 							FOLKCODE=c("G", "sG", "S", "gS", "mG", "mS", "gmS", "M", "gM", "sM", "gsM"))
#
# ### Use percent gravel, sand and mud to assign phi later
# #merge haul.dt with folk class
# haul.dt.gmx2.out2 <- merge(haul.dt.gmx2.out, gmx.folk.class, by=c("FOLK_num"), all.x=T)







# ==============================================================
# = Use Bathy Grid to interpolate usSEABED sediment point data =
# ==============================================================

# ==============================================
# = Get bathy in same formats as sediment data =
# ==============================================
### Convert raster to spatial points dataframe

raster2sed <- function(r, sed){
	bathy.atl.df <- as.data.frame(r, xy=T)

	### Get in same format as sediment data
	bathy.atl.df$WATERDEPTH <- ifelse(bathy.atl.df$z<0 & bathy.atl.df$z>=-3000, -1*bathy.atl.df$z, NA)
	bathy.atl.df$LONGITUDE <- bathy.atl.df$x
	bathy.atl.df$LATITUDE <- bathy.atl.df$y
	bathy.atl.df$z <- NULL
	bathy.atl.df$x <- NULL
	bathy.atl.df$y <- NULL

	coordinates(bathy.atl.df) <- ~ LONGITUDE + LATITUDE
	bathy.atl.df@proj4string <- sed@proj4string
	return(bathy.atl.df)
}

# bathy.atl.df <- raster2sed(atl.bathy, atl)
bathy.pac.df <- raster2sed(pac.bathy, pac)
# bathy.gmx.df <- raster2sed(gmx.bathy, gmx)



# =================
# = Interpolation =
# =================
### By inverse distance weighting
### Grain Size
pac.grain.idw <- idw(GRAINSIZE~1, subset(pac, GRAINSIZE !=-99), subset(bathy.pac.df, !(is.na(WATERDEPTH))))
# gmx.grain.idw <- idw(GRAINSIZE~1, subset(gmx, GRAINSIZE !=-99), subset(bathy.gmx.df, !(is.na(WATERDEPTH))), idp=3)
# atl.grain.idw <- idw(GRAINSIZE~1, subset(atl, GRAINSIZE !=-99), subset(bathy.atl.df, !(is.na(WATERDEPTH))), idp=3)


### Gravel percentage
pac.gravel.idw <- idw(GRAVEL ~1, subset(pac, GRAVEL !=-99), subset(bathy.pac.df, !(is.na(WATERDEPTH))))
# atl.gravel.idw <- idw(GRAVEL ~1, subset(atl, GRAVEL !=-99), subset(bathy.atl.df, !(is.na(WATERDEPTH))), idp=3)
# gmx.gravel.idw <- idw(GRAVEL ~1, subset(gmx, GRAVEL !=-99), subset(bathy.gmx.df, !(is.na(WATERDEPTH))), idp=3)

### Sand percentage
pac.sand.idw <- idw(SAND ~1, subset(pac, SAND !=-99), subset(bathy.pac.df, !(is.na(WATERDEPTH))))

### Gravel percentage
pac.mud.idw <- idw(MUD ~1, subset(pac, MUD !=-99), subset(bathy.pac.df, !(is.na(WATERDEPTH))))


# spplot(atl.grain.idw, "var1.pred")
spplot(pac.grain.idw, "var1.pred")
# spplot(gmx.grain.idw, "var1.pred")

# spplot(atl.gravel.idw, "var1.pred")
spplot(pac.gravel.idw, "var1.pred")
# spplot(gmx.gravel.idw, "var1.pred")




# atl.grain.idw.r <- rasterFromXYZ(as.data.frame(atl.grain.idw)[,c("LONGITUDE", "LATITUDE", "var1.pred")])
pac.grain.idw.r <- rasterFromXYZ(as.data.frame(pac.grain.idw)[,c("LONGITUDE", "LATITUDE", "var1.pred")])
# gmx.grain.idw.r <- rasterFromXYZ(as.data.frame(gmx.grain.idw)[,c("LONGITUDE", "LATITUDE", "var1.pred")])

# atl.gravel.idw.r <- rasterFromXYZ(as.data.frame(atl.gravel.idw)[,c("LONGITUDE", "LATITUDE", "var1.pred")])
pac.gravel.idw.r <- rasterFromXYZ(as.data.frame(pac.gravel.idw)[,c("LONGITUDE", "LATITUDE", "var1.pred")])
# gmx.gravel.idw.r <- rasterFromXYZ(as.data.frame(gmx.gravel.idw)[,c("LONGITUDE", "LATITUDE", "var1.pred")])

pac.mud.idw.r <- rasterFromXYZ(as.data.frame(pac.mud.idw)[,c("LONGITUDE", "LATITUDE", "var1.pred")])
pac.sand.idw.r <- rasterFromXYZ(as.data.frame(pac.sand.idw)[,c("LONGITUDE", "LATITUDE", "var1.pred")])

pac.sed.stack <- stack(pac.grain.idw.r, pac.gravel.idw.r,  pac.sand.idw.r, pac.mud.idw.r)
names(pac.sed.stack) <- c("GRAINSIZE", "GRAVEL","SAND",  "MUD")
plot(pac.sed.stack)
saveRDS(pac.sed.stack, "sed_grid_out/pac.sed.stack.out.rds")




# =================================
# = Interpolate at haul locations =
# =================================

haul.pac <- subset(haul.dt, region %in% c("NWFSC_WCAnn", "AFSC_WCTri") )
proj4string(haul.pac) <- proj4string(pac)

haul.pac.grain.idw <- idw(GRAINSIZE~1, subset(pac, GRAINSIZE !=-99), haul.pac)
haul.pac.gravel.idw <- idw(GRAVEL ~1, subset(pac, GRAVEL !=-99), haul.pac)
haul.pac.sand.idw <- idw(SAND ~1, subset(pac, SAND !=-99), haul.pac)
haul.pac.mud.idw <- idw(MUD ~1, subset(pac, MUD !=-99), haul.pac)

haul.pac.dt <- as.data.table(as.data.frame(haul.pac))
haul.pac.dt[,"GRAINSIZE":=haul.pac.grain.idw$var1.pred]
haul.pac.dt[,"GRAVEL":=haul.pac.gravel.idw$var1.pred]
haul.pac.dt[,"SAND":=haul.pac.sand.idw$var1.pred]
haul.pac.dt[,"MUD":=haul.pac.mud.idw$var1.pred]


# saveRDS(haul.pac.dt, "haul_sp/hauls_sed_pac_idw.rds")

haul.pac.sp <- as.data.frame(copy(haul.pac.dt))
coordinates(haul.pac.sp) <- cbind(haul.pac.sp$coords.x1, haul.pac.sp$coords.x2)
proj4string(haul.pac.sp) <- proj4string(conmap)
# writeOGR(haul.pac.sp, "haul_sp", "haul_pac_out_sp", driver="ESRI Shapefile")





# wentworth <- data.frame(mm=c(4096,256,64,4,2,1,0.5,0.25,0.125,0.0625,0.031,0.0156,0.0078, 0.0039, 0.00006),
# 	phi=c(-12,-8, -6, -2, -1, 0, 1, 2, 3, 4, 5, 6, 7, 8, 14))
# wentworth.labs <- data.frame(size.class=c("Boulder", "Cobble", "Pebble", "Granule",
# 	"Very coarse sand", "Coarse sand", "Medium sand", "Fine sand", "Very fine sand",
# 	"Coarse silt", "Medium silt", "Fine silt", "Very fine silt", "Clay"))
# library(RColorBrewer)
# wentworth.cols <-  colorRampPalette(brewer.pal(n=9, "Spectral"))(15)
#
# plot(atl.grain.idw.r, breaks=wentworth$phi, col=wentworth.cols)
# plot(atl.gravel.idw.r)

