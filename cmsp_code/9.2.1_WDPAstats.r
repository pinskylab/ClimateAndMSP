# stats about impacts on protected areas from shifting species

## Set working directories
if(Sys.info()["nodename"] == "pinsky-macbookair"){
	setwd('~/Documents/Rutgers/Range projections/proj_ranges/')
	presmapbymodfolder <- '../data/'
	}
if(Sys.info()["nodename"] == "amphiprion.deenr.rutgers.edu"){
	setwd('~/Documents/range_projections/')
	.libPaths(new='~/R/x86_64-redhat-linux-gnu-library/3.1/') # so that it can find my old packages
	presmapbymodfolder <- 'data/'
	}
# could add code for Lauren's working directory here

############
## Flags
############

## choose which runs to use
## runtype refers to how the Species Distribution Models (SDMs) were fit
## projtype refers to how the SDM projections were done
#runtype <- 'test'; projtype=''
#runtype <- ''; projtype=''`
#runtype <- 'testK6noSeas'; projtype='_xreg'
runtype <- 'fitallreg'; projtype='_xreg'

# choose the rcp
rcp <- 85
otherrcp <- 45

# select initial and final timeperiod for these grids
periods <- c('2006-2020', '2081-2100')

####################
## helper functions
####################
require(Hmisc)
require(data.table)
require(lme4) # for mixed-effects models
require(car) # for testing ME models


lu <- function(x) return(length(unique(x)))

se <- function(x,na.rm=FALSE){ # standard error
	if(!na.rm){
		return(sd(x, na.rm=FALSE)/sqrt(length(x)))
	}
	if(na.rm){
		return(sd(x, na.rm=TRUE)/sqrt(sum(!is.na(x))))
	}
}

# weighted mean function to use with summarize()
wmean <- function(x){ # values in col 1, weights in col 2
	inds <- !is.na(x[,1]) & !is.na(x[,2])
	return(weighted.mean(x[inds,1], x[inds,2])) 
}


###########################################
# Examine MPA size relative to grid size
###########################################
wdpaturnbyMPAbymod <- as.data.table(read.csv(paste('data/wdpaturnbyMPAbymod_', runtype, projtype, '_rcp', rcp, '&', otherrcp, '.csv', sep=''), row.names=1)) # to get a list of MPAs to examine

wdpashp <- readShapePoly('cmsp_data/WDPA/NA_marine_MPA/mpinsky-search-1382225374362.shp', proj4string = CRS('+proj=longlat +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +no_defs')) # load the MPA data

wdpashp@data$area_deg <- sapply(slot(wdpashp, 'polygons'), slot, 'area') # extract area in degrees (since lat/lon projection)
wdpashp@data$wdpapolyID <- as.numeric(sapply(wdpashp@polygons, FUN=slot, name='ID')) # extract polygon ID

	dim(wdpashp)
wdpashp <- wdpashp[wdpashp$wdpapolyID %in% wdpaturnbyMPAbymod$wdpapolyID,] # trim shapefile to MPAs the we analyzed
	dim(wdpashp) # 562
	
summary(wdpashp@data$area_deg)

# plot vs. grid size
hist(log10(wdpashp@data$area_deg))
abline(v=log10(1/16), col='red')

# number > grid size
nrow(wdpashp)
sum(wdpashp@data$area_deg >= 1/16)
sum(wdpashp@data$area_deg >= 1/16 * 1/16)


#############################################
# Examine turnover within MPAs
# Examine all climate models individually
#############################################
# read in turnover data
wdpaturnbyMPAbymod <- as.data.table(read.csv(paste('data/wdpaturnbyMPAbymod_', runtype, projtype, '_rcp', rcp, '&', otherrcp, '.csv', sep=''), row.names=1))

ntk <- wdpaturnbyMPAbymod$no_take %in% c('All', 'Part') # no take reserves (index into wdpaturnbyMPAbymod)
sum(ntk[wdpaturnbyMPAbymod$rcp==rcp & wdpaturnbyMPAbymod$model==1]) # 50 (43 in CA)

# read in change in temp by region
trend1 <- read.csv('data/climTrendbyreg_rcp45.csv', row.names=1)
	trend1$rcp <- 45
trend2 <- read.csv('data/climTrendbyreg_rcp85.csv', row.names=1)
	trend2$rcp <- 85
nms <- c('region', 'rcp', 'delta_surf', 'delta_bott')
trend <- as.data.table(rbind(trend1[,nms], trend2[,nms]))


# Fraction of species lost
	# all MPAs
	wdpaturnbyMPAbymod[,.(mean=mean(flost,na.rm=TRUE)),by=c('rcp','model')][,.(mean=mean(mean), se=se(mean), min=min(mean), max=max(mean)),by='rcp'] # first mean is across MPAs within climate models, then average across climate models

	wdpaturnbyMPAbymod[,.(mean=mean(flost,na.rm=TRUE)),by=c('rcp','model','region.one')][,.(mean=mean(mean), se=se(mean), min=min(mean), max=max(mean)),by=c('rcp','region.one')] # first mean is across MPAs within climate models, then average within regions

	mod<- wdpaturnbyMPAbymod[,.(mean=mean(flost,na.rm=TRUE)),by=c('model','region.one')][,lm(mean~region.one-1)] # first mean is across MPAs within climate models within regions, then do an anova
		summary(mod)
		length(mod$fitted.values) # number of data points

	mod<- wdpaturnbyMPAbymod[,.(mean=mean(flost,na.rm=TRUE)),by=c('model','region.one')][,kruskal.test(x=mean, g=region.one)] # first mean is across MPAs within climate models within regions, then do a non-parametric Kruskal-Wallis Rank Sum Test
		mod

#	wdpaturnbyMPAbymod[,sd(flost,na.rm=TRUE),by=c('rcp','model')]
#	wdpaturnbyMPAbymod[,min(flost,na.rm=TRUE),by=c('rcp','model')]
#	wdpaturnbyMPAbymod[,max(flost,na.rm=TRUE),by=c('rcp','model')]

	# no-take
#	wdpaturnbyMPAbymod[ntk,mean(flost,na.rm=TRUE),by=c('rcp','model')]
#	wdpaturnbyMPAbymod[ntk,sd(flost,na.rm=TRUE),by=c('rcp','model')]
#	wdpaturnbyMPAbymod[ntk,min(flost,na.rm=TRUE),by=c('rcp','model')]
#	wdpaturnbyMPAbymod[ntk,max(flost,na.rm=TRUE),by=c('rcp','model')]

# Fraction of species gained (fraction of final community)
	# all MPAs
	wdpaturnbyMPAbymod[,.(mean=mean(fgainedalt,na.rm=TRUE)),by=c('rcp','model')][,.(mean=mean(mean), se=se(mean), min=min(mean), max=max(mean)),by='rcp'] # first mean is across MPAs within climate models

	wdpaturnbyMPAbymod[,.(mean=mean(fgainedalt,na.rm=TRUE)),by=c('rcp','model','region.one')][,.(mean=mean(mean), se=se(mean), min=min(mean), max=max(mean)),by=c('rcp','region.one')] # first mean is across MPAs within climate models, then average within regions

	mod<- wdpaturnbyMPAbymod[,.(mean=mean(fgainedalt,na.rm=TRUE)),by=c('model','region.one')][,lm(mean~region.one-1)] # first mean is across MPAs within climate models within regions, then do an anova
		summary(mod)
		length(mod$fitted.values) # number of data points

	mod<- wdpaturnbyMPAbymod[,.(mean=mean(fgainedalt,na.rm=TRUE)),by=c('model','region.one')][,kruskal.test(x=mean, g=region.one)] # first mean is across MPAs within climate models within regions, then do a non-parametric Kruskal-Wallis Rank Sum Test
		mod

#	b <- wdpaturnbyMPAbymod[,sd(fgainedalt,na.rm=TRUE),by=c('rcp','model')]; b # sd across MPAs within climate models (rcp45 or rcp85)
#		b[,mean(V1),by='rcp'] 


	# no-take
#	wdpaturnbyMPAbymod[ntk,mean(fgained,na.rm=TRUE),by=c('rcp','model')]
#	wdpaturnbyMPAbymod[ntk,sd(fgained,na.rm=TRUE),by=c('rcp','model')]

# Similarity
	# all MPAs
	wdpaturnbyMPAbymod[,.(mean=mean(beta_sor,na.rm=TRUE)),by=c('rcp','model')][,.(mean=mean(mean), se=se(mean), min=min(mean), max=max(mean)),by='rcp'] # first mean is across MPAs within climate models

	# examine similarity by MPA size
		wdpaturnbyMPAbymod[,cor.test(rep_area, beta_sor), by=c('rcp', 'model')]
		
			# average across rcps and climate models
		a <- wdpaturnbyMPAbymod[rep_area != 0,.(beta_sor=mean(beta_sor), abeta_sor=asin(sqrt(mean(beta_sor))), larea=log(mean(rep_area))),by=wdpapolyID]
		mod <- a[,lm(abeta_sor ~ larea)]; summary(mod)
		a[,plot(larea, beta_sor)] # plot all points
		i <- order(a$larea); a[, lines(larea[i], sin(fitted(mod)[i])^2)] # add fit line
		a[,range(exp(larea))]
		range(sin(fitted(mod))^2)


	# examine similarity by MPA latitudinal range
		datsc <- as.data.frame(wdpaturnbyMPAbymod[,.(lat_min, lat_max, beta_sor, rcp, model)])
		datsc$latrng <- datsc$lat_max - datsc$lat_min
		datsc$rcp <- as.factor(datsc$rcp)
		datsc$model <- as.factor(datsc$model)
#		datsc$latrngsc <- scale(datsc$latrng)
		datsc$y <- asin(sqrt(datsc$beta_sor))
		mod <- lmer(y ~ latrng + (1|rcp/model), data=datsc, REML=FALSE)
		mod2 <- lmer(y ~ 1 + (1|rcp/model), data=datsc, REML=FALSE)
		summary(mod)
		Anova(mod)
		anova(mod, mod2)

		plot(datsc$latrng, datsc$beta_sor)
		fit <- data.frame(latrng=datsc$latrng)
		fit$y <- fixef(mod)[1] + fit$latrng * fixef(mod)[2]
		fit$beta_sor <- sin(fit$y)^2
		fit <- fit[order(fit$latrng),]
		lines(fit$latrng, fit$beta_sor, col='red')

		head(fit)
		tail(fit)

	# examine similarity by MPA size and lat rng
		datsc <- as.data.frame(wdpaturnbyMPAbymod[rep_area != 0,.(rep_area, lat_min, lat_max, beta_sor, rcp, model)]) # remove zero area
		datsc$lrep_area <- log(datsc$rep_area)
		datsc$latrng <- datsc$lat_max - datsc$lat_min
		datsc$rcp <- as.factor(datsc$rcp)
		datsc$model <- as.factor(datsc$model)
		pvars <- c('lrep_area', 'latrng')
		scvars <- c('lrep_areasc', 'latrngsc')
		datsc[scvars] <- lapply(datsc[pvars], scale)
		datsc$y <- asin(sqrt(datsc$beta_sor))
		mod <- lmer(y ~ lrep_areasc*latrngsc + (1|rcp/model), data=datsc, REML=FALSE)
		mod2 <- lmer(y ~ lrep_areasc + latrngsc + (1|rcp/model), data=datsc, REML=FALSE)
		mod3.1 <- lmer(y ~ latrngsc + (1|rcp/model), data=datsc, REML=FALSE)
		mod3.2 <- lmer(y ~ lrep_areasc + (1|rcp/model), data=datsc, REML=FALSE)
		mod4 <- lmer(y ~ 1 + (1|rcp/model), data=datsc, REML=FALSE)
		summary(mod)
		Anova(mod)
		anova(mod, mod2)
		anova(mod3.1, mod3.2) 

		plot(mod, sin(fitted(.))^2 ~ lrep_areasc)
		plot(mod, sin(fitted(.))^2 ~ latrngsc)
	
		require(interplot) # to plot the interaction
		interplot(m=mod, var1='latrngsc', var2='lrep_areasc') +
			xlab('log(area) scaled') + 
			ylab('Coefficient for scaled latitudinal range')

		interplot(m=mod, var2='latrngsc', var1='lrep_areasc') +
			ylab('Coefficient for log(area) scaled') + 
			xlab('Latitudinal range scaled')

	# examine similarity by region change in temperature
		# calculate mean turnover by region
		turnbyreg <- wdpaturnbyMPAbymod[,.(beta_sor=mean(beta_sor)),by=.(rcp,region.one)]
		
		# merge together
		setkey(trend, region, rcp)
		setkey(turnbyreg, region.one, rcp)
		turnbyreg <- trend[turnbyreg] # merge
		
		# examine
		turnbyreg[rcp==85,]
		turnbyreg[rcp==45,]
		
		# plot
		turnbyreg[rcp==85,plot(delta_bott, beta_sor)]
		turnbyreg[rcp==45,plot(delta_bott, beta_sor)]

		# linear model
		turnbyreg[, cor.test(beta_sor, delta_bott), by=rcp]		
		turnbyreg[, cor.test(beta_sor, delta_surf), by=rcp]		

	# examine similarity by region change in temperature AND by lat rng
		# average across climate models and rcps
		turnbyMPA <- wdpaturnbyMPAbymod[,.(beta_sor=mean(beta_sor), latrng=mean(lat_max) - mean(lat_min), region=unique(region.one)),by=.(wdpapolyID)]
		turnbyMPA85 <- wdpaturnbyMPAbymod[rcp==85,.(beta_sor=mean(beta_sor), latrng=mean(lat_max) - mean(lat_min), region=unique(region.one)),by=.(wdpapolyID)]
		
		# average temp trends across rcps
		trendave <- trend[,.(delta_surf=mean(delta_surf), delta_bott=mean(delta_bott)), by=region]
		
		# merge in regional temperature information
		setkey(trendave, region)
		setkey(turnbyMPA, region)
		setkey(turnbyMPA85, region)
		turnbyMPA <- trendave[turnbyMPA] # merge
		turnbyMPA85 <- trendave[turnbyMPA85] # merge
		
		# linear model
		mod <- turnbyMPA[, lm(asin(sqrt(beta_sor)) ~ latrng*delta_surf)]
		mod2 <- turnbyMPA[, lm(asin(sqrt(beta_sor)) ~ latrng+delta_surf)] # used this in the paper, since not significantly worse than full model
		summary(mod)
		summary(mod2)
		anova(mod, mod2)
		AIC(mod, mod2)

			# examine data
				# rcp85 and 45 averaged
			turnbyMPA[latrng<=0.25,.(mean=mean(beta_sor), se=se(beta_sor))]
			turnbyMPA[latrng>=5,.(mean=mean(beta_sor), se=se(beta_sor))]

			turnbyMPA[delta_surf<2.5,.(mean=mean(beta_sor), se=se(beta_sor))]
			turnbyMPA[delta_surf>3.5,.(mean=mean(beta_sor), se=se(beta_sor))]

				# rcp85 only
			turnbyMPA85[latrng<=0.25,.(mean=mean(beta_sor), se=se(beta_sor))]
			turnbyMPA85[latrng>=5,.(mean=mean(beta_sor), se=se(beta_sor))]

			turnbyMPA85[delta_surf<2.5,.(mean=mean(beta_sor), se=se(beta_sor))]
			turnbyMPA85[delta_surf>3.5,.(mean=mean(beta_sor), se=se(beta_sor))]

			# effect size calculations
			preds <- expand.grid(latrng=seq(0,10,length.out=41), delta_surf=seq(2,4,length.out=41))
			preds$beta_sor_trans <- predict(mod2, newdata=preds)
			preds$beta_sor_pred <- sin(preds$beta_sor_trans)^2

				# plot
				par(mfrow=c(1,2))
				plot(preds$latrng, preds$beta_sor_pred, ylim=c(0,1))
					i <- preds$delta_surf==3
					lines(preds$latrng[i], preds$beta_sor_pred[i], col='red')
				plot(preds$delta_surf, preds$beta_sor_pred, ylim=c(0,1))
					j <- preds$latrng==0.25 # (28 km MPA)
					lines(preds$delta_surf[j], preds$beta_sor_pred[j], col='red')

				# look at predictions
				preds[i,]
				preds[j,]

# Plot of fraction species lost as a density across MPAs within each climate model
	inds <- expand.grid(rcp=unique(wdpaturnbyMPAbymod$rcp), model=unique(wdpaturnbyMPAbymod$model))
	ds <- vector('list', nrow(inds))
	for(i in 1:nrow(inds)){
		ds[[i]] <- density(wdpaturnbyMPAbymod$flost[wdpaturnbyMPAbymod$rcp==inds$rcp[i] & wdpaturnbyMPAbymod$model==inds$model[i]], na.rm=TRUE, cut=1, kernel='gaus', adjust=2)
	}
	
	cols <- c('blue', 'green') # blue is rcp, green is otherrcp
	plot(ds[[1]], type='n', ylim=c(0,7), xlab='Fraction of species lost', main='')
	for(i in 1:nrow(inds)){	
		if(inds$rcp[i]==rcp) lines(ds[[i]], col=cols[1])
		if(inds$rcp[i]==otherrcp) lines(ds[[i]], col=cols[2])
	}


################################################################
# Examine change within MPA networks and the component MPAs
# Across each model in the ensemble
################################################################
wdpaturnbynetbymod <- as.data.table(read.csv(paste('data/wdpaturnbynetbymod_', runtype, projtype, '_rcp', rcp, '&', otherrcp, '.csv', sep=''), row.names=1)) # network results
wdpaturnbyMPAbymod <- as.data.table(read.csv(paste('data/wdpaturnbyMPAbymod_', runtype, projtype, '_rcp', rcp, '&', otherrcp, '.csv', sep=''), row.names=1)) # individual MPA results

# Size of networks
	wdpaturnbyMPAbymod[!is.na(network),.(num_mpas=length(unique(wdpapolyID))), by=network]
	wdpaturnbyMPAbymod[!is.na(network) & lat_min>=35.25 & network=='eastcoast',.(num_mpas=length(unique(wdpapolyID))), by=network] # correct eastcoast to only north of Cape Hatteras (where we have data)

# Fraction of species lost
	# all MPA networks
	a <- wdpaturnbynetbymod[,mean(flost,na.rm=TRUE),by=c('rcp','model')]; a
		a[,range(V1),by='rcp'] # 4%-19% or 11%-27% min
		a[,mean(V1),by='rcp'] # 10% or 21% mean flost (rcp45 or rcp85)
		a[,median(V1),by='rcp'] # 
		a[,se(V1),by='rcp'] # 1.3% or 1.2% se across models (rcp45 or rcp85)
	wdpaturnbyMPAbymod[,sd(flost,na.rm=TRUE),by=c('rcp','model')]
	wdpaturnbyMPAbymod[,min(flost,na.rm=TRUE),by=c('rcp','model')]
	wdpaturnbyMPAbymod[,max(flost,na.rm=TRUE),by=c('rcp','model')]

	# within the individual MPAs of these networks
	a <- wdpaturnbyMPAbymod[!is.na(network),mean(flost,na.rm=TRUE),by=c('rcp','model')]; a
		a[,range(V1),by='rcp'] # 
		a[,mean(V1),by='rcp'] # 
		a[,median(V1),by='rcp'] # 
		a[,se(V1),by='rcp'] # 
	
	# average across GCMs and RCPs, keep networks separate
	a <- wdpaturnbynetbymod[,mean(flost,na.rm=TRUE),by=network]; a

# Fraction of species gained (fraction of final community)
	# average across all MPA networks, keep rcps and gcms separate
	a <- wdpaturnbynetbymod[,mean(fgainedalt,na.rm=TRUE),by=c('rcp','model')]; a
		a[,mean(V1),by='rcp'] # 
		a[,se(V1),by='rcp'] # 
		a[,range(V1),by='rcp'] # 
		a[,median(V1),by='rcp'] # 
		a[,sd(V1),by='rcp'] # 
	b <- wdpaturnbyMPAbymod[,sd(fgainedalt,na.rm=TRUE),by=c('rcp','model')]; b
		b[,mean(V1),by='rcp'] # 

	# within the individual MPAs of these networks
	a <- wdpaturnbyMPAbymod[!is.na(network),mean(fgainedalt,na.rm=TRUE),by=c('rcp','model')]; a
		a[,mean(V1),by='rcp'] # 
		a[,se(V1),by='rcp'] # 
		a[,range(V1),by='rcp'] # 
		a[,median(V1),by='rcp'] # 

	# average across GCMs and RCPs, keep networks separate
	a <- wdpaturnbynetbymod[,mean(fgainedalt,na.rm=TRUE),by=network]; a

	# average across GCMs, keep networks and RCPs separate
	a <- wdpaturnbynetbymod[,mean(fgainedalt,na.rm=TRUE),by=c('network','rcp')]; a

# Similarity (Sorenson)
	# all MPA networks
	a <- wdpaturnbynetbymod[,mean(beta_sor,na.rm=TRUE),by=c('rcp','model')]; a
		a[,mean(V1),by='rcp'] # 
		a[,se(V1),by='rcp'] # 
		a[,range(V1),by='rcp'] # 
		a[,median(V1),by='rcp'] # 
		a[,sd(V1),by='rcp'] # 
	b <- wdpaturnbyMPAbymod[,sd(fgainedalt,na.rm=TRUE),by=c('rcp','model')]; b
		b[,mean(V1),by='rcp'] # 

	# within the individual MPAs of these networks
	a <- wdpaturnbyMPAbymod[!is.na(network),mean(beta_sor,na.rm=TRUE),by=c('rcp','model')]; a
		a[,mean(V1),by='rcp'] # 
		a[,se(V1),by='rcp'] # 
		a[,range(V1),by='rcp'] # 
		a[,median(V1),by='rcp'] # 