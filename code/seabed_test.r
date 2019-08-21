setwd('/Users/mpinsky/Documents/Rutgers/Range projections/proj_ranges')
require(data.table)
require(RColorBrewer)
require(maps)

# Read in usSEABED data
# east coast
eclc <- fread('unzip -c data_seabed/usSEABED/atl_clc.zip ATL_CLC.TXT | tail -n +3') # unzip -c to read out a specific file, tail to skip first two header lines from unzip
eprs <- fread('unzip -c data_seabed/usSEABED/atl_prs.zip ATL_PRS.TXT | tail -n +3')
eext <- fread('unzip -c data_seabed/usSEABED/atl_ext.zip ATL_EXT.TXT | tail -n +3')

# gulf coast
# each file spits out some warnings, but for a column I don't need
gclc <- fread('unzip -c data_seabed/usSEABED/gmx_clc.zip GMX_CLC.TXT | tail -n +3') # unzip -c to read out a specific file, tail to skip first two header lines from unzip
gprs <- fread('unzip -c data_seabed/usSEABED/gmx_prs.zip GMX_PRS.TXT | tail -n +3')
gext <- fread('unzip -c data_seabed/usSEABED/gmx_ext.zip GMX_EXT.TXT | tail -n +3')

# pacific coast, lower 48
pclc <- fread('unzip -c data_seabed/usSEABED/pac_clc.zip PAC_CLC.txt | tail -n +3') # unzip -c to read out a specific file, tail to skip first two header lines from unzip
	setnames(pclc, c(1,2,16), c('Latitude', 'Longitude', 'Grainsize'))
pprs <- fread('unzip -c data_seabed/usSEABED/pac_prs.zip PAC_PRS.txt | tail -n +3', sep=',')
	setnames(pprs, c(1,2,16), c('Latitude', 'Longitude', 'Grainsize'))
pext <- fread('unzip -c data_seabed/usSEABED/pac_ext.zip PAC_EXT.txt | tail -n +3', sep=',')
	setnames(pext, c(1,2,16), c('Latitude', 'Longitude', 'Grainsize'))

# Add source
eclc[,src:='clc']
eprs[,src:='prs']
eext[,src:='ext']

gclc[,src:='clc']
gprs[,src:='prs']
gext[,src:='ext']

pclc[,src:='clc']
pprs[,src:='prs']
pext[,src:='ext']

# Combine datasets
cs <- c('Latitude', 'Longitude', 'Grainsize', 'src')
dat <- rbind(eclc[,cs,with=FALSE], eprs[,cs,with=FALSE], eext[,cs,with=FALSE], gclc[,cs,with=FALSE], gprs[,cs,with=FALSE], gext[,cs,with=FALSE], pclc[,cs,with=FALSE], pprs[,cs,with=FALSE], pext[,cs,with=FALSE])
	rm(eclc, eprs, eext, gclc, gprs, gext, pclc, pprs, pext)
	
# Find NAs
dat[Grainsize==-99, Grainsize:=NA]
summary(dat$Grainsize)
sum(is.na(dat$Grainsize)) # NAs
sum(!is.na(dat$Grainsize)) # useful values

# initial plot of coverage
dat[!is.na(Grainsize),plot(Longitude, Latitude, pch=16, cex=0.5, col=rgb(0.5, 0.5, 0.5, 0.1))] # good coverage

dat[,hist(Grainsize)] # -7 to 10

# rescale
dat[,gs:=(Grainsize-min(Grainsize, na.rm=TRUE))/(max(Grainsize, na.rm=TRUE)-min(Grainsize, na.rm=TRUE))]
	dat[,hist(gs)] # 0 to 1

# plot grainsize
colrmp <- colorRamp(colors=brewer.pal(9, 'OrRd'))
cols <- function(x, alpha){
	rgbs <- colrmp(x)
	return(rgb(rgbs[,1], rgbs[,2], rgbs[,3], alpha=alpha, maxColorValue=255))
}
#dat[!is.na(gs),plot(Longitude, Latitude, pch=16, cex=0.5, col=rgb(0.1, 0.1, 0.1, 0.1))] # good coverage
dat[!is.na(gs),plot(Longitude, Latitude, pch=16, cex=0.2, col=cols(gs, 40))] # good coverage



###################
## NGDC data

# Read in NGDC data
# Alaska
angdcsamp <- fread('data_seabed/ngdc_seafloor_sediment/alaska/sizesample1479922967033.txt') # reads in 26 cols: last one is blank
	angdcsamp[,V26:=NULL]
	nms <- unlist(strsplit(readLines('data_seabed/ngdc_seafloor_sediment/alaska/sizesample1479922967033.txt', n=1), split='\t')) # read col names
	setnames(angdcsamp, nms)

# Canada eastern
cngdcsamp <- fread('data_seabed/ngdc_seafloor_sediment/eastern_canada/sizesample1479924088781.txt') # reads in 26 cols: last one is blank
	cngdcsamp[,V26:=NULL]
	nms <- unlist(strsplit(readLines('data_seabed/ngdc_seafloor_sediment/eastern_canada/sizesample1479924088781.txt', n=1), split='\t')) # read col names
	setnames(cngdcsamp, nms)


# plot
angdcsamp[,plot(lon, lat, pch=16, col=rgb(0.2, 0.2, 0.2, 0.2), cex=0.5)] # alaska
cngdcsamp[,plot(lon, lat, pch=16, col=rgb(0.2, 0.2, 0.2, 0.2), cex=0.5)] # eastern canada
	map(add=TRUE)