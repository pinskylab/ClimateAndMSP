## Another script to test whether thermal envelopes change through time (non-stationarity)
## This examines summary stats on models to first decade of dataset and projected later

## Set working directories depending on computer
if(Sys.info()["nodename"] == "pinsky-macbookair"){
	setwd('~/Documents/Rutgers/Range projections/proj_ranges/')
	modfolder <- '../CEmodels_nonstationaritytest/'
	nthreads=2
	}
if(Sys.info()["nodename"] == "amphiprion.deenr.rutgers.edu"){
	setwd('~/Documents/range_projections/')
	.libPaths(new='~/R/library/') # so that it can find my old packages
	modfolder <- 'CEmodels_nonstationaritytest/'
	nthreads=10
	}
if(Sys.info()["user"] == "lauren"){
	setwd('~/backup/NatCap/proj_ranges/')
	modfolder <- 'output/CEmodels_nonstationaritytest/'
	nthreads=1
}

# Load model evaluations
modeldiag <- read.csv('output/modeldiag_nonstationarity_projection.csv')

##########################################
# Plot stats over decades for each spp
##########################################
spp <- sort(unique(modeldiag$sppocean[!is.na(modeldiag$auc)]))

# On individual points
col=rgb(0, 0, 0, 0.2)
quartz(width=12, height=9)
par(mfrow=c(3,4))
	# auc
with(modeldiag[modeldiag$sppocean==spp[1],], plot(decade, auc, type='l', col=col, ylim=c(0,1)))
for(i in 2:length(spp)) with(modeldiag[modeldiag$sppocean==spp[i],], lines(decade, auc, col=col))
	# tss
with(modeldiag[modeldiag$sppocean==spp[1],], plot(decade, tss, type='l', col=col, ylim=c(0,1)))
for(i in 2:length(spp)) with(modeldiag[modeldiag$sppocean==spp[i],], lines(decade, tss, col=col))
	# acc
with(modeldiag[modeldiag$sppocean==spp[1],], plot(decade, acc, type='l', col=col, ylim=c(0,1)))
for(i in 2:length(spp)) with(modeldiag[modeldiag$sppocean==spp[i],], lines(decade, acc, col=col))
	# sens
with(modeldiag[modeldiag$sppocean==spp[1],], plot(decade, sens, type='l', col=col, ylim=c(0,1)))
for(i in 2:length(spp)) with(modeldiag[modeldiag$sppocean==spp[i],], lines(decade, sens, col=col))
	# spec
with(modeldiag[modeldiag$sppocean==spp[1],], plot(decade, spec, type='l', col=col, ylim=c(0,1)))
for(i in 2:length(spp)) with(modeldiag[modeldiag$sppocean==spp[i],], lines(decade, spec, col=col))
	# kappa
with(modeldiag[modeldiag$sppocean==spp[1],], plot(decade, kappa, type='l', col=col, ylim=c(0,1)))
for(i in 2:length(spp)) with(modeldiag[modeldiag$sppocean==spp[i],], lines(decade, kappa, col=col))
	# tssmax
with(modeldiag[modeldiag$sppocean==spp[1],], plot(decade, tssmax, type='l', col=col, ylim=c(0,1)))
for(i in 2:length(spp)) with(modeldiag[modeldiag$sppocean==spp[i],], lines(decade, tssmax, col=col))
	# accmax
with(modeldiag[modeldiag$sppocean==spp[1],], plot(decade, accmax, type='l', col=col, ylim=c(0,1)))
for(i in 2:length(spp)) with(modeldiag[modeldiag$sppocean==spp[i],], lines(decade, accmax, col=col))
	# kappamax
with(modeldiag[modeldiag$sppocean==spp[1],], plot(decade, kappamax, type='l', col=col, ylim=c(0,1)))
for(i in 2:length(spp)) with(modeldiag[modeldiag$sppocean==spp[i],], lines(decade, kappamax, col=col))
	# rpb
with(modeldiag[modeldiag$sppocean==spp[1],], plot(decade, rpb, type='l', col=col, ylim=c(0,1)))
for(i in 2:length(spp)) with(modeldiag[modeldiag$sppocean==spp[i],], lines(decade, rpb, col=col))

# On gridded data
col=rgb(0, 0, 0, 0.2)
quartz(width=12, height=9)
par(mfrow=c(3,4))
	# auc
with(modeldiag[modeldiag$sppocean==spp[1],], plot(decade, auc.gr, type='l', col=col, ylim=c(0,1)))
for(i in 2:length(spp)) with(modeldiag[modeldiag$sppocean==spp[i],], lines(decade, auc.gr, col=col))
	# tss
with(modeldiag[modeldiag$sppocean==spp[1],], plot(decade, tss.gr, type='l', col=col, ylim=c(0,1)))
for(i in 2:length(spp)) with(modeldiag[modeldiag$sppocean==spp[i],], lines(decade, tss.gr, col=col))
	# acc
with(modeldiag[modeldiag$sppocean==spp[1],], plot(decade, acc.gr, type='l', col=col, ylim=c(0,1)))
for(i in 2:length(spp)) with(modeldiag[modeldiag$sppocean==spp[i],], lines(decade, acc.gr, col=col))
	# sens
with(modeldiag[modeldiag$sppocean==spp[1],], plot(decade, sens.gr, type='l', col=col, ylim=c(0,1)))
for(i in 2:length(spp)) with(modeldiag[modeldiag$sppocean==spp[i],], lines(decade, sens.gr, col=col))
	# spec
with(modeldiag[modeldiag$sppocean==spp[1],], plot(decade, spec.gr, type='l', col=col, ylim=c(0,1)))
for(i in 2:length(spp)) with(modeldiag[modeldiag$sppocean==spp[i],], lines(decade, spec.gr, col=col))
	# kappa
with(modeldiag[modeldiag$sppocean==spp[1],], plot(decade, kappa.gr, type='l', col=col, ylim=c(0,1)))
for(i in 2:length(spp)) with(modeldiag[modeldiag$sppocean==spp[i],], lines(decade, kappa.gr, col=col))
	# tssmax
with(modeldiag[modeldiag$sppocean==spp[1],], plot(decade, tssmax.gr, type='l', col=col, ylim=c(0,1)))
for(i in 2:length(spp)) with(modeldiag[modeldiag$sppocean==spp[i],], lines(decade, tssmax.gr, col=col))
	# accmax
with(modeldiag[modeldiag$sppocean==spp[1],], plot(decade, accmax.gr, type='l', col=col, ylim=c(0,1)))
for(i in 2:length(spp)) with(modeldiag[modeldiag$sppocean==spp[i],], lines(decade, accmax.gr, col=col))
	# kappamax
with(modeldiag[modeldiag$sppocean==spp[1],], plot(decade, kappamax.gr, type='l', col=col, ylim=c(0,1)))
for(i in 2:length(spp)) with(modeldiag[modeldiag$sppocean==spp[i],], lines(decade, kappamax.gr, col=col))
	# rpb
with(modeldiag[modeldiag$sppocean==spp[1],], plot(decade, rpb.gr, type='l', col=col, ylim=c(0,1)))
for(i in 2:length(spp)) with(modeldiag[modeldiag$sppocean==spp[i],], lines(decade, rpb.gr, col=col))