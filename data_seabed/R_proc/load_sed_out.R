library(raster)
library(data.table)
setwd("~/Documents/Collaborations/Rutgers/Sediment Grain/")
ak <- readRDS("sed_grid_out/ak.sed.stack.out.rds")
can <- readRDS("sed_grid_out/canada.sed.stack.out.rds")
wcan <- readRDS("sed_grid_out/canadaW.sed.stack.out.rds")
atl <- readRDS("sed_grid_out/conmap.r.stack.rds")
gmx <- readRDS("sed_grid_out/gmx_db.r.rds")
pac <- readRDS("sed_grid_out/pac.sed.stack.out.rds")

plot(ak)
plot(can)
plot(wcan)
plot(atl)
plot(gmx)
plot(pac)

### Create mosaic, with mean where overlap
rast.list <- list(ak, wcan, pac, gmx, atl, can)
rast.list$fun <- mean
sed.mosaic <- do.call(mosaic, rast.list)
names(sed.mosaic) <- names(ak)
plot(sed.mosaic, xlim=c(-180, -120), ylim=c(30,60))
plot(sed.mosaic, xlim=c(-100, -40))

saveRDS(sed.mosaic, "sed_grid_out/sed.mosaic.rds")

load("sed_grid_out/bathy_grid_course_Feb7_2017.RData")
proj.grid.bath[,"lonBathgrid.adj":=ifelse(lonBathgrid< -180, 180+(lonBathgrid+180), lonBathgrid)]
coordinates(proj.grid.bath) <- cbind(proj.grid.bath$lonBathgrid.adj, proj.grid.bath$latBathgrid)

sed.out <- extract(sed.mosaic, proj.grid.bath)
sed.out.df <- as.data.frame(sed.out)

sed.out.dt <- as.data.table(as.data.frame(proj.grid.bath))
sed.out.dt[,"GRAINSIZE":=sed.out.df$GRAINSIZE]
sed.out.dt[,"GRAVEL":=sed.out.df$GRAVEL]
sed.out.dt[,"SAND":=sed.out.df$SAND]
sed.out.dt[,"MUD":=sed.out.df$MUD]

saveRDS(sed.out.dt, "sed_grid_out/sed.projgrid.rds")

png("sed_grid_out/OverlaySedNew.png", height=5, width=5, res=300, units="in")
plot(sed.mosaic[[1]], xlim=c(-180, -40))
points(sed.out.dt$lonBathgrid.adj, sed.out.dt$latBathgrid, cex=0.25)
plot(sed.mosaic[[1]], xlim=c(-180, -40), add=T)
dev.off()

sed.out.dt[is.na(GRAINSIZE)]



# ===================
# = Halpern Habitat =
# ===================
sed.out.dt <- readRDS("sed_grid_out/sed.projgrid.rds")


halp <- readRDS("Halpern2008_Data/halpern_hab_stack.rds")
plot(halp[[2]])

proj4string(proj.grid.bath) <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"
proj.grid.moll <- spTransform(proj.grid.bath, crs(halp))



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



halp.out <- class_hab(halp, proj.grid.moll)
halp.out[,"coords.x1":=NULL]
halp.out[,"coords.x2":=NULL]
halp.out[,"value":=NULL]

benthic_hab_proj <- merge(sed.out.dt, halp.out, 
	by=c("lonClimgrid", "latClimgrid", "latBathgrid", "lonBathgrid", "depth", "rugosity", "lonBathgrid.adj"))
saveRDS(benthic_hab_proj, "sed_grid_out/benthic_hab_proj.rds")