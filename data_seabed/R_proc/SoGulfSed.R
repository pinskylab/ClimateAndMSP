setwd("~/Documents/Collaborations/Rutgers/Sediment Grain/")

library(data.table)
library(rgdal)

### Load digitized sediment data
sogulf_sed <- readOGR("Gulf of St Lawrence", "GulfofStLawrence_Sed_Digitized")
spplot(sogulf_sed, "Sed_Code")


sogulf_data <- as.data.table(sogulf_sed@data)

sogulf_sedtypes <- sogulf_data[,list(num.rec=length(unique(pkuid))), by=list(Sed_Code)]
# write.csv(sogulf_sedtypes, "Gulf of St Lawrence/sogulf_sedtypes.csv")

sogulf_sedtypes2 <- read.csv("Gulf of St Lawrence/sogulf_sedtypes_complete.csv")

# =====================================
# = Folk Classification to Grain Size from usSEABED=
# =====================================
folk.phi <- readRDS("folk_usSEABED_perc_phi.rds")

sogulf_phi <- merge(sogulf_sedtypes2, folk.phi, by="FOLKCODE", all.x=T)

# ==========================
# = Extract data for hauls =
# ==========================
load("haulInfo_Dec26_2016.RData")
haul.dt <- as.data.table(subset(hauls, !(is.na(lat))& !(is.na(lon)))) # 69 hauls without lat lon
### make longitudes < -180 equal to 180 - (180+lon)
haul.dt[,"lon.adj":=ifelse(lon< -180, 180+(180+lon), lon)]

sogulf.hauls <- haul.dt[region == "DFO_SoGulf"]
coordinates(sogulf.hauls) <- cbind(sogulf.hauls$lon, sogulf.hauls$lat)
proj4string(sogulf.hauls) <- proj4string(sogulf_sed)

sogulf.hauls$Sed_Code <- over(sogulf.hauls, sogulf_sed)$Sed_Code # points on polygon boundary and points corresponding to polygon vertex are considered inside polygon
sogulf.hauls2 <- merge(sogulf.hauls, sogulf_phi, by="Sed_Code", all.x=T)

sogulf.hauls.out <- as.data.table(as.data.frame(sogulf.hauls2)) 

sogulf.hauls.sp <- copy(sogulf.hauls.out)

coordinates(sogulf.hauls.sp) <- cbind(sogulf.hauls.sp$coords.x1, sogulf.hauls.sp$coords.x2)
proj4string(sogulf.hauls.sp) <- proj4string(sogulf_sed)
writeOGR(sogulf.hauls.sp, "haul_sp", "haul_sogulf_out_sp", driver="ESRI Shapefile")
writeOGR(subset(sogulf.hauls.sp, is.na(GRAINSIZE)), "haul_sp", "haul_sogulf_GRAIN_NA", driver="ESRI Shapefile")


sogulf.hauls.out[,"Sed_Code":=NULL]
sogulf.hauls.out[,"FOLKCODE":=NULL]
sogulf.hauls.out[,"num.rec":=NULL]
sogulf.hauls.out[,"Legend.Code":=NULL]


saveRDS(sogulf.hauls.out, "haul_sp/hauls.sogulf.poly.rds")