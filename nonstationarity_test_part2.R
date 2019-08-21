## Script to test whether thermal envelopes change through time (non-stationarity)
## This is part 2: analyze the model fits

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


# Load packages
library(mgcv);library(ROCR)

# Functions
# moving window mean with a given step size. window width = step size
# uses first column for determining steps, but averages all columns individually
windowmean <- function(dat, step){ 
	stind <- round(dat[,1]/step)
	out <- apply(dat, MARGIN=2, function(x) aggregate(x, by=list(stind), FUN=mean, na.rm=TRUE)[,2])
	return(as.data.frame(out))
}

############################
## Load and analyze models
############################
modfiles <- list.files(modfolder, full.names=TRUE)

sppocean <- gsub(paste(modfolder, '/CEmods_nonstationaritytest_|.RData', sep=''), '', modfiles)
nonstation <- data.frame(sppocean=sppocean, pT1=NA, pT2=NA) # hold the preferred temps from first and second half

nrow(nonstation)
for(i in 1:length(modfiles)){
	print(i)
	load(modfiles[i]) # loads mods, myregions

#	plot(mods$mygam1, pages=1)

	# predictions across temperature gradient for first and second half
	predregion <- names(which.max(table(mods$mygam1$model$regionfact))) # get the most common region in the fitting data
	if(is.null(predregion)) predregion=NA # if model only has one region
	predbio <- mean(mods$mygam1$model$biomassmean, na.rm=TRUE)
	predrug <- mean(mods$mygam1$model$logrugosity, na.rm=TRUE)
	preds <- mods$mygam1$model[,c('bottemp', 'surftemp')] # use the observed bottemp and surftemps
	preds <- windowmean(preds, step=1) # average in blocks of 1deg
	preds$half <- 1
	preds2 <- preds; preds2$half <- 2
	preds <- rbind(preds, preds2)
	preds <- merge(preds, data.frame(regionfact=predregion, biomassmean=predbio, logrugosity=predrug))
	preds$se <- preds$probpres <- NA

	preds <- preds[order(preds$bottemp),]

#	preds[,c('probpres', 'se')] <- predict(mods$mygam1, newdata=preds, se.fit=TRUE, type='response')
	preds[,c('probpres', 'se')] <- predict(mods$mygam1, newdata=preds, se.fit=TRUE, type='link')

	dim(preds)
	
	# find top quartile
	thresh1 <- quantile(preds$probpres[preds$half==1], 0.5)
	thresh2 <- quantile(preds$probpres[preds$half==2], 0.5)


	# plot the curves (set up for type='link')
	plot(NA, NA, xlim=c(-2,35), ylim=range(preds$probpres), xlab='bottemp', ylab='logit(pres)')
	polygon(x=c(preds$bottemp[preds$half==1], rev(preds$bottemp[preds$half==1])), y=c(preds$probpres[preds$half==1] + preds$se[preds$half==1], rev(preds$probpres[preds$half==1] - preds$se[preds$half==1])), border=FALSE, col='light blue')
	lines(probpres ~ bottemp, data=preds[preds$half==1,], col='blue') # for type=link
	polygon(x=c(preds$bottemp[preds$half==2], rev(preds$bottemp[preds$half==2])), y=c(preds$probpres[preds$half==2] + preds$se[preds$half==2], rev(preds$probpres[preds$half==2] - preds$se[preds$half==2])), border=FALSE, col='light grey')
	lines(probpres ~ bottemp, data=preds[preds$half==2,])
	abline(h=thresh1, lty=2, col='blue')
	abline(h=thresh2, lty=2)
	
	# preferred temperature as a weighted mean of top quartile
	inds1 <- preds$half==1 & preds$probpres > thresh1 & preds$se < 1
	inds2 <- preds$half==2 & preds$probpres > thresh2 & preds$se < 1
		sum(inds1)
		sum(inds2)
	nonstation$pT1[i] <- with(preds[inds1,], weighted.mean(bottemp, w=probpres-thresh1, na.rm=TRUE))
	nonstation$pT2[i] <- with(preds[inds2,], weighted.mean(bottemp, w=probpres-thresh2, na.rm=TRUE))
		abline(v=nonstation$pT1[i], col='blue')
		abline(v=nonstation$pT2[i], col='black')

	
}

nonstation$diff <- nonstation$pT2 - nonstation$pT1
nonstation[order(nonstation$diff),] # sorted by diff
sum(abs(nonstation$diff)<1, na.rm=TRUE)/sum(!is.na(nonstation$diff)) # fraction < 1deg
sum(is.na(nonstation$diff)) # ones that didn't work

hist(abs(nonstation$diff), breaks=seq(0, max(abs(nonstation$diff), na.rm=TRUE), length.out=40), col='grey') # histogram of abs(diff)
summary(abs(nonstation$diff)) # summary of abs(diff)

# write out
write.csv(nonstation, 'output/nonstationarity_test1.csv')