# Climate-envelope models, using data from Trawl Data/

#############
### MODEL ###	
#############
setwd('/Users/mpinsky/Documents/Rutgers/Range projections/')
source('allsubset.glm 2012-01-30.R')
source('revgrep 2012-02-07.R')

load('Output/trawl_allregionsforprojections_2015-02-02.RData') # load dat data.frame
    dat$stratumfact = as.factor(dat$stratum)
    dat$yrfact = as.factor(dat$year)
    dat$sppregion = paste(dat$spp, dat$region, sep='_')

## Models
	require(mgcv)


## Pick the spp
spp=sort(unique(dat$spp));spp # need to adjust this to pick the species we want
sppreg = sort(unique(dat$sppregion))


## Pick the training years for each region
regs = sort(unique(dat$region))
yrs = vector('list', length(regs)); names(yrs) = regs
for(i in 1:length(yrs)){
	yrs[[i]] = sort(unique(dat$year[dat$region == regs[i]]))
}

trainingyrs = yrs; traintype = 'allyrs' # all years

# Indicator for where present, but also has at least one TRUE value per spp per stratum within the training years, for fitting stratum effects
## OR: Skip this and read it in below
	data$presfit = data$wtcpue > 0
	presfitspp = character(0)
	for(i in 1:length(regs)){
		strat = sort(unique(data$stratum[data$region == regs[i]]))
		thesespp = sort(unique(data$spp[data$region == regs[i]]))
		for(j in 1:length(thesespp)){
			for(k in 1:length(strat)){
				inds = data$region == regs[i] & data$spp == thesespp[j] & data$stratum==strat[k] & data$year %in% trainingyrs[[i]] & complete.cases(data[,c('surftemp', 'bottemp')]) # also has to have complete data, otherwise useless for model fitting
				if(sum(data$presfit[inds]) == 0){
					print(paste(regs[i], thesespp[j],strat[k]))
					m = sample(which(inds),1) # pick a random record in this stratum
					data$presfit[m] = TRUE # set it to TRUE (even though wtcpue==0)
					print(paste(data$presfit[m], data$wtcpue[m], data$wtcpuena[m], data$year[m]))
					presfitspp = c(presfitspp, as.character(data$sppregion[m]))
				}
			}
		}
	}
	sum(data$presfit & !data$pres)
	presfitspp = unique(presfitspp); presfitspp

	write.csv(data[,c('region', 'haulid', 'spp', 'presfit')], file=paste('Output/CEModels/datapresfit_', traintype, '_', Sys.Date(), '.csv', sep=''))

	# Read in presfit from before:
	datapresfitname = '../Trawl Data/Output/CEModels 2012-08-05 allyrs/datapresfit_allyrs_2012-08-05.csv'
	presfit = read.csv(datapresfitname, row.names=1, stringsAsFactors=FALSE)
	dim(data)
	dim(presfit)
	all(presfit$haulid == data$haulid & presfit$spp == data$spp & presfit$region == data$region) # make sure ALL are in the same order
 data = cbind(data, presfit[, 'presfit'])
		dim(data) # should add one column
	

options(warn=0)
for(j in 1:length(regs)){
	print(j)
	spp=sort(unique(data$spp[data$region==regs[j]]))
	sppinds = 1:length(spp)
	
	for(i in sppinds){
		# Subset to training years
		if(!(regs[j] %in% c('DFO_Newfoundland_Spring', 'DFO_Newfoundland_Fall'))){
			#inds = data$yearsurv %in% trainingyrs & complete.cases(data[,c('surftemp', 'bottemp', 'mintemp', 'maxtemp')]) & data$spp == spp[i]
#			inds = data$yearsurv %in% trainingyrs[[j]] & complete.cases(data[,c('surftemp', 'bottemp', 'depth')]) & data$spp == spp[i] & data$region == regs[j]
			inds = data$yearsurv %in% trainingyrs[[j]] & complete.cases(data[,c('surftemp', 'bottemp')]) & data$spp == spp[i] & data$region == regs[j]
		} else {
#			inds = data$yearsurv %in% trainingyrs[[j]] & complete.cases(data[,c('bottemp', 'depth')]) & data$spp == spp[i] & data$region == regs[j] # no surftemp in Newfoundland
			inds = data$yearsurv %in% trainingyrs[[j]] & complete.cases(data[,c('bottemp')]) & data$spp == spp[i] & data$region == regs[j] # no surftemp in Newfoundland
		}
		
		print(paste(regs[j], ' ', spp[i], ': ', sum(data$pres[inds]), ' presences, ', sum(inds), ' indices in training set', sep=''))
		
		if(sum(inds)>300 & sum(data$pres[inds])>=40){
				
			# GAM gaussian
#			print('GAM gaussian')		
#			if(!(regs[j] %in% c('DFO_Newfoundland_Spring', 'DFO_Newfoundland_Fall'))){
##				gam2 = gam(wtcpue ~ s(surftemp) + s(bottemp) + s(depth) + stratumfact + biomassmean, family=gaussian, data=data[inds,], select=TRUE) # mgcv has a built-in model selection function
#				gam2 = gam(wtcpue ~ s(surftemp) + s(bottemp) + stratumfact + biomassmean, family=gaussian, data=data[inds,], select=TRUE) # mgcv has a built-in model selection function
#			} else {
##				gam2 = gam(wtcpue ~ s(bottemp) + s(depth) + stratumfact + biomassmean, family=gaussian, data=data[inds,], select=TRUE) # w/out surftemp for Newfoundland
#				gam2 = gam(wtcpue ~ s(bottemp) + stratumfact + biomassmean, family=gaussian, data=data[inds,], select=TRUE) # w/out surftemp for Newfoundland
#			
#			}
			#	summary(gam2)
			#	plot(gam2, se=TRUE, shade=TRUE, pages=1, all.terms=TRUE)
		
			# GAM Hurdle gamma (pres/abs and abundance)
			print('GAM presence binomial')	
			if(!(regs[j] %in% c('DFO_Newfoundland_Spring', 'DFO_Newfoundland_Fall'))){
#				gamhp2 = gam(pres ~ s(surftemp) + s(bottemp) + s(depth) + stratumfact + biomassmean, family=binomial, data=data[inds,], select=TRUE) # mgcv has a built-in model selection function
				gamhp2 = gam(pres ~ s(surftemp) + s(bottemp) + stratumfact + biomassmean, family=binomial, data=data[inds,], select=TRUE) # mgcv has a built-in model selection function
			} else {
#				gamhp2 = gam(pres ~ s(bottemp) + s(depth) + stratumfact + biomassmean, family=binomial, data=data[inds,], select=TRUE)
				gamhp2 = gam(pres ~ s(bottemp) + stratumfact + biomassmean, family=binomial, data=data[inds,], select=TRUE)
			}

			print('GAM abundance gaussian(log)')
			try3 = try(if(!(regs[j] %in% c('DFO_Newfoundland_Spring', 'DFO_Newfoundland_Fall'))){
#				gamha2 = gam(wtcpuenal ~ s(surftemp) + s(bottemp) + s(depth) + stratumfact + biomassmean, family=gaussian, data=data[inds & data$presfit,], select=TRUE, control=list(mgcv.half=30)) # mgcv has a built-in model selection function
				gamha2 = gam(wtcpuenal ~ s(surftemp) + s(bottemp) + stratumfact + biomassmean, family=gaussian, data=data[inds & data$presfit,], select=TRUE, control=list(mgcv.half=30)) # mgcv has a built-in model selection function
			} else {
#				gamha2 = gam(wtcpuenal ~ s(bottemp) + s(depth) + stratumfact + biomassmean, family=gaussian, data=data[inds & data$presfit,], select=TRUE) # mgcv has a built-in model selection function
				gamha2 = gam(wtcpuenal ~ s(bottemp) + stratumfact + biomassmean, family=gaussian, data=data[inds & data$presfit,], select=TRUE) # mgcv has a built-in model selection function
			}
			)
			
			# Save models if model fitting worked
			if(class(try3)[1] != "try-error"){
				mods = list(gam2=gam2, gamhp2=gamhp2, gamha2 = gamha2) # without full models or glms
				save(mods, inds, file=paste('Output/CEModels/CEmods_', regs[j], '_', spp[i], '_', traintype, '_', Sys.Date(), '.RData', sep=''))
			}
		}
	}
}

##########################################
######## Plots and basic analysis ########
##########################################
setwd('/Users/mpinsky/Documents/Princeton/Trawl Data/')
require(mgcv)
data = read.csv('Output/dataCEM_all_2012-11-16.csv', row.names=1)
    data$stratumfact = as.factor(data$stratum)
    data$yrfact = as.factor(data$year)
    data$sppregion = paste(data$spp, data$region, sep='_')
	spp=sort(unique(data$spp))
	regs = sort(unique(data$region))

	# Read in presfit from before:
	datapresfitname = 'Output/CEModels 2012-08-05 allyrs/datapresfit_allyrs_2012-08-05.csv'; traintype = 'allyrs'
	presfit = read.csv(datapresfitname, row.names=1)
		dim(data)
		dim(presfit)
		all(presfit$haulid == data$haulid & presfit$spp == data$spp & presfit$region == data$region) # make sure ALL are in the same order
	data = cbind(data, presfit[, 'presfit'])
		dim(data) # should add one column
		names(data)[grep('presfit', names(data))] = 'presfit' # fix column name if concatenated on strangely

	yrs = vector('list', length(regs)); names(yrs) = regs
	for(i in 1:length(yrs)){
		j = sort(unique(data$year[data$region == regs[i]]))
		if(regs[i] == 'DFO_SoGulf') j = j[j != 2008] # because no SST in 2008
		yrs[[i]] = j
	}

	# Write out for later use
	save(yrs, regs, spp, file = paste('Output/CEM_yrs_regs_spp', Sys.Date(), '.RData', sep=''))


## Plot GAM fits for all the abundance models (multipage pdf)
	folder = 'CEModels 2012-07-26 allyr'; trainingyrs = yrs; traintype = 'allyrs'; date='2012-07-25'; vars=c('surftemp', 'bottemp') # subset of 12 years
	# folder = 'CEModels 2012-02-17 halfinter biomassmean'; trainingyrs = yrs[seq(1,length(yrs), by=2)] # subset of 12 years
	#folder = 'CEModels 2012-02-18 half biomassmean'; trainingyrs = yrs[1:12]; vars=c('bottemp', 'surftemp') # subset of 12 years
	# folder = 'CEModels 2012-02-23 2ndhalf'; trainingyrs = yrs[13:length(yrs)]; date='2012-02-23'; vars=c('surftemp', 'bottemp') # subset of 12 years

	files = list.files(path=paste('Output/', folder, sep=''), pattern='.RData')
	pdf(width=6, height=6, file=paste('Figures/GAMabund_fits_', traintype, '_', Sys.Date(), '.pdf', sep='')) # only works if only one round of CEM models are in the directory
	for(i in 1:length(files)){
		print(files[i])
		load(file=paste('Output/',folder, '/', files[i], sep=''))
		thisspp = strsplit(files[i], split=paste('AFSC_Aleutians_|AFSC_EBS_|AFSC_GOA_|DFO_ScotianShelf_|DFO_SoGulf_|NEFSCSpring_|WestCoast_Tri_|DFO_Newfoundland_Fall_|SEFSC_GOMex_|_',traintype, sep=''))[[1]][2] # extract spp name from file name
		thisreg = strsplit(files[i], split=paste('CEmods_|_', thisspp, sep=''))[[1]][2]
		theseyrs = yrs[[thisreg]]
		inds = data$year %in% trainingyrs[[thisreg]] & complete.cases(data[,vars]) & data$spp == spp[i]
	
		plot(mods$gam2, shade=TRUE, se=TRUE, pages=1, main=paste(thisspp, '\n', thisreg), all.terms=FALSE)
		sum = summary(mods$gam2)#; print(sum)
		mtext(side=3, paste(signif(sum$dev.expl,2)*100, '% dev', sep=''))	
	}
		
	dev.off()

## Plot GAM Hurdle fits for all the abundance models (multipage pdf)
	folder = 'CEModels 2012-07-26 allyr'; trainingyrs = yrs; date='2012-07-25'; vars=c('surftemp', 'bottemp') # subset of 12 years
	#folder = 'CEModels 2012-03-14 firsthalf'; trainingyrs = yrs; date='2012-03-14'; vars = c('bottemp', 'surftemp');	for(i in 1:length(yrs)) trainingyrs[[i]] = yrs[[i]][1:round(length(yrs[[i]])/2)]; traintype = 'firsthalf' # first half of years

	files=list.files(paste('Output/', folder, sep=''), pattern='.RData') # only works if only one round of CEM models are in the directory
	regs = sort(unique(data$region))

	pdf(width=6, height=6, file=paste('Figures/GAMhurdle_fits_', traintype, '_', Sys.Date(), '.pdf', sep=''))
	for(i in 1:length(files)){
		print(files[i])
		load(file=paste('Output/',folder, '/', files[i], sep=''))
		thisspp = strsplit(files[i], split=paste('AFSC_Aleutians_|AFSC_EBS_|AFSC_GOA_|DFO_ScotianShelf_|DFO_SoGulf_|NEFSCSpring_|WestCoast_Tri_|DFO_Newfoundland_Fall_|SEFSC_GOMex_|_',traintype, sep=''))[[1]][2] # extract spp name from file name
		thisreg = strsplit(files[i], split=paste('CEmods_|_', thisspp, sep=''))[[1]][2]
		theseyrs = yrs[[thisreg]]

		inds = data$year %in% trainingyrs[[thisreg]] & complete.cases(data[,vars]) & data$spp == spp[i]
		
		plot(mods$gamhp2, shade=TRUE, se=TRUE, pages=1, main=paste(thisspp, '\n', thisreg, '\npresence'), all.terms=FALSE)
		sum = summary(mods$gamhp2)#; print(sum)
		mtext(side=3, paste(signif(sum$dev.expl,2)*100, '% dev', sep=''))

		# quartz()
		plot(mods$gamha2, shade=TRUE, se=TRUE, pages=1, main=paste(thisspp, '\n', thisreg, '\nabundance'), all.terms=FALSE)
		sum = summary(mods$gamha2)#; print(sum)	
		mtext(side=3, paste(signif(sum$dev.expl,2)*100, '% dev', sep=''))
		
		rm(mods, inds)
	
	}
	
	dev.off()


# Plot predictions and reality as maps
	require(maps)
	#yrtoplot = list(one=1968:1989, two= 1990:2011)
	spp = 'Hippoglossusstenolepis'; region = 'AFSC_EBS'; yrtoplot = list(one=c(1982, 1983), two= 2005); load('Output/CEModels 2012-03-14 firsthalf/CEmods_AFSC_EBS_Hippoglossusstenolepis_firsthalf_2012-03-14.RData'); xlims = c(-180,-155); ylims = c(52,63); db='world'
	
	spp = 'Homarus americanus'; region = 'NEFSCSpring'; yrtoplot = list(one=c(1968:1979), two=2001:2009); load('Output/CEModels 2012-03-14 firsthalf/CEmods_NEFSCSpring_Homarus americanus_firsthalf_2012-03-14.RData'); xlims = c(-77,-65); ylims = c(35,46); db='usa'; wid=3.5

	spp = 'Homarus americanus'; region = 'NEFSC_Spring'; yrtoplot = list(one=c(1970), two=2005); load('Output/CEModels 2012-09-05 allyrs onlytemp/CEmods_NEFSC_Spring_Homarus americanus_allyrs_2012-09-06.RData'); xlims = c(-77,-65); ylims = c(35,46); db='usa'; wid=3.5

	spp = 'Doryteuthis pealeii'; region = 'NEFSCSpring'; yrtoplot = list(one=c(1968:1979), two=2001:2009); load('Output/CEModels 2012-03-14 firsthalf/CEmods_NEFSCSpring_Doryteuthis pealeii_firsthalf_2012-03-14.RData'); xlims = c(-77,-65); ylims = c(35,46); db='usa'; wid=3.5

	spp = 'Merluccius bilinearis'; region = 'NEFSCSpring'; yrtoplot = list(one=c(1968:1979), two=2001:2009); load('Output/CEModels 2012-03-12 allyr/CEmods_NEFSCSpring_Merluccius bilinearis_allyrs_2012-03-12.RData'); xlims = c(-77,-65); ylims = c(35,46); db='usa'; wid=3.5

	spp = 'Merluccius bilinearis'; region = 'NEFSC_Spring'; yrtoplot = list(one=c(1968:2011), two=2001:2009); load('Output/CEModels 2012-09-05 allyrs onlytemp/CEmods_NEFSC_Spring_Merluccius bilinearis_allyrs_2012-09-06.RData'); xlims = c(-77,-65); ylims = c(35,46); db='usa'; wid=3.5

	spp = 'Clupeapallasi'; region = 'AFSC_GOA'; yrtoplot = list(one=c(1984, 1987, 1990, 1993,1996,1999), two= c(2003,2005, 2007, 2009, 2011)); load('Output/CEModels 2012-03-14 firsthalf/CEmods_AFSC_GOA_Clupeapallasi_firsthalf_2012-03-14.RData'); xlims = c(-170,-130); ylims = c(53,62); db='world'; wid=6
	

	inds = data$spp == spp & data$region == region & data$year %in% as.numeric(unlist(yrtoplot))
		sum(inds)
	c = data$wtcpue[inds] # to scale pt size by biomass
	a = max(log(c), na.rm=T)/2
	cols = rev(rainbow(n=ceiling(a)+1, start=0, end=4/6))
	col = rgb(10,255,0,150,maxColorValue=255)
	prescol='black'; abscol='grey'; prespch=16; abspch=4; abscex=0.4 # colors and pchs for presences and absences
	#prescol='black'; abscol='black'; prespch=16; abspch=16; abscex=0.1 # for all points the same
	#mapcol = 'light grey'
	mapcol = 'tan'

	mods2 = list(gamhurdle = list(mods$gamhp2, mods$gamha2))
	# plot(mods$gamhp2, shade=TRUE, residuals=TRUE, pages=1)
	# plot(mods$gamha2, shade=TRUE, residuals=TRUE, pages=1)


	quartz(width=min(12,length(yrtoplot)*wid), height=4)
	# pdf(paste('Figures/Maps_CEM_', spp, '_', region, '_', Sys.Date(), '.pdf', sep=''), width=min(12,length(yrtoplot)*wid), height=4)
	par(mfcol=c(length(mods2),length(yrtoplot)), mai=c(0.3, 0.3, 0.2, 0.1)) # use mfcol to fill by columns
	if(length(yrtoplot)==2) par(mai=c(0.5, 0.5, 0.3, 0.1)) # use mfcol to fill by columns
	for(i in 1:length(yrtoplot)){
		for(j in 1:length(mods2)){
			if(length(yrtoplot[[i]])==1) title = yrtoplot[[i]]
			if(length(yrtoplot[[i]])>1) title = paste(yrtoplot[[i]][1], 'to', max(yrtoplot[[i]]))			
			#map(database='state', regions=c('massachusetts', 'rhode island', 'maine', 'connecticut', 'new york', 'new jersey', 'delaware', 'new hampshire', 'vermont', 'pennsylvania', 'maryland'), fill=TRUE, col='grey', xlim=xlims, ylim=ylims, main=title)
			map(database=db, xlim=xlims, ylim=ylims, fill=TRUE, col=mapcol)
			map.axes()
			mtext(title, side=3)
			inds = data$year %in% yrtoplot[[i]] & data$spp == spp & complete.cases(data[,c('surftemp', 'bottemp')]) & data$region == region
			inds2 = inds & data$wtcpue > 0
			inds3 = inds & data$wtcpue == 0
			c = data$wtcpue[inds2]
			#c = rep(1.1, sum(inds)); a = 10 # for all pts small
			points(data$lon[inds3], data$lat[inds3], pch=abspch, cex=abscex, col = abscol) # plot absences
			points(data$lon[inds2], data$lat[inds2], pch=prespch, cex=ceiling(log(c))/a, col = prescol) # plot observed points
			meanlat = weighted.mean(data$lat[inds2],w= data$wtcpue[inds2])
			#points(min(xlims)+0.5, meanlat, cex=2, pch=1) # plot mean observed lat
			if(i==1) baselat = meanlat
			#abline(h=meanlat, lty=2, col='black')
			#abline(h=min(data$lat[inds2]), lty=1, col='black') # min and max lat
			#abline(h=max(data$lat[inds2]), lty=1, col='black')

			if(length(mods2[[j]]) == 1){
				preds = predict(mods2[[j]][[1]], newdata = data[inds,], type='response')
			}
			if(length(mods2[[j]]) == 2){
				preds1 = predict(mods2[[j]][[1]], newdata = data[inds,], type='response')
				preds2 = exp(predict(mods2[[j]][[2]], newdata = data[inds,], type='response'))
				preds = preds1*preds2
			}
			preds[preds<0] = 0
			#points(data$lon[inds][preds>0], data$lat[inds][preds>0], pch=16, xlab='', cex=ceiling(log(preds[preds>0]))/a, ylab='', col = col) # plot predictions
			meanpredlat = weighted.mean(data$lat[inds][!is.na(preds)], w = preds[!is.na(preds)])
			#points(min(xlims)+0.5, meanpredlat, cex=2, col=col, pch=1) # plot mean predicted lat
			#if(i==1) basepredlat = meanpredlat
			#abline(h=basepredlat, lty=2, col=col)
			#abline(h=min(data$lat[inds][preds>7e-1], na.rm=TRUE), lty=1, col=col) # min and max predicted lat
			#abline(h=max(data$lat[inds][preds>7e-1], na.rm=TRUE), lty=1, col=col)
		}
	}
	
	dev.off()
	
## Model diagnostics: calculate ROC and r2 for each model
	require(ROCR)
#	folder = 'CEModels 2012-07-26 allyr'; trainingyrs = yrs; date='2012-07-25'; vars=c('surftemp', 'bottemp'); traintype = 'allyr' # all years
	#folder = 'CEModels 2012-03-14 firsthalf'; trainingyrs = yrs; date='2012-03-14'; vars = c('bottemp', 'surftemp'); traintype = 'firsthalf' # first half of years
	#for(i in 1:length(yrs)) trainingyrs[[i]] = yrs[[i]][1:round(length(yrs[[i]])/2)]; traintype = 'firsthalf'
	folder = 'CEModels 2012-09-05 allyrs onlytemp'; trainingyrs = yrs; vars=c('surftemp', 'bottemp'); traintype='allyrs'; suff='onlytemp' # all years w/ only temp (surf, bot, strat, biomass)

	files=list.files(paste('Output/', folder, sep=''), pattern='.RData') # only works if only one round of CEM models are in the directory
	regs = sort(unique(data$region))
	sppregs = sort(unique(data$sppregion))
	n = rep(NA, length(sppregs))
	
	# save ROC and r2 for each taxon in each region
	modeldiaggamh = data.frame(sppregion = sppregs, region = n, spp = n, auc = n, r2.biomass = n, r2.all = n, dev.pres = n, dev.biomass = n, stringsAsFactors=FALSE)
		modeldiaggamh = modeldiaggamh[!grepl('PORIFERA', modeldiaggamh$sppregion),] # for some reason PORIFERA gets added, but I don't have a CEM
		
	options(warn=1) # print warnings as they occur
	fileinds = 1:length(files)
	#fileinds = 221:length(files)
	#fileinds = grep('Newfoundland', files)
	for(s in fileinds){
		print(s)
		load(paste('Output/',folder, '/', files[s], sep=''))
		thisspp = strsplit(files[s], split=paste('AFSC_Aleutians_|AFSC_EBS_|AFSC_GOA_|DFO_ScotianShelf_|DFO_SoGulf_|NEFSC_Spring_|WestCoast_Tri_|DFO_Newfoundland_Fall_|SEFSC_GOMex_|_',traintype, sep=''))[[1]][2] # extract spp name from file name
		#thisreg = strsplit(files[s], split=paste('CEmods_|_', thisspp, sep=''))[[1]][2]
		thisreg = regs[sapply(regs, FUN=grepl, x=files[s])]
		theseyrs = yrs[[thisreg]]
		j = intersect(grep(thisspp, modeldiaggamh$sppregion, fixed=TRUE), grep(thisreg, modeldiaggamh$sppregion))
		if(length(j)>1 | length(j)<1) stop('Found too few or too many indices into modeldiag dataframe')
		modeldiaggamh$spp[j] = thisspp
		modeldiaggamh$region[j] = as.character(thisreg)
	
		# Find observations and predictions
		inds = complete.cases(data[,vars]) & data$spp == thisspp & data$region == thisreg
		if(thisreg == 'DFO_Newfoundland_Fall') inds = complete.cases(data[,vars[vars!='surftemp']]) & data$spp == thisspp & data$region == thisreg

		if(sum(inds)>0){
			preds1 = predict(mods$gamhp2, newdata = data[inds,], type='response')
			preds2 = exp(predict(mods$gamha2, newdata = data[inds,], type='response')) # make sure to use exp if wtcpue was log-transformed!
			#plot(preds1, preds2, main=paste(spp[s], yrs[i]), xlab='Predicted presence', ylab='Predicted wtcpue')
			preds = preds1*preds2
			preds[preds<0] = 0

			# calculate performance (in part using ROCR)
			preds1.rocr = prediction(predictions=as.numeric(preds1), labels=data$pres[inds])
				# plot the ROC curve
#				perf = performance(preds1.rocr, 'tpr', 'fpr')
#				plot(perf)
			modeldiaggamh$auc[j] = performance(preds1.rocr, 'auc')@y.values[[1]]
			modeldiaggamh$r2.biomass[j] = cor(preds2[data$presfit[inds]], data$wtcpue[inds & data$presfit])^2 # correlation of biomass where present
			modeldiaggamh$r2.all[j] = cor(preds, data$wtcpue[inds])^2 # overall biomass correlation
			modeldiaggamh$dev.pres[j] = summary(mods$gamhp)$dev.expl
			modeldiaggamh$dev.biomass[j] = summary(mods$gamha)$dev.expl
		}
	}
	
	### Save model diagnostics
	write.csv(modeldiaggamh, file = paste('Output/modeldiaggamh_', traintype, suff, '_', Sys.Date(), '.csv', sep=''))

	## Examine
	modeldiaggamh = read.csv('Output/modeldiaggamh_allyrsonlytemp_2012-09-26.csv', row.names=1)
	dropspp = c('PARALEPIDIDAE', 'MYCTOPHIDAE', 'NOTACANTHIDAE', 'LIPARIDAE', 'Actiniaria', 'Nudibranchia', 'CEPHALOPODA', 'Porifera', 'PORIFERA') # spp that aren't even at the genus level
	modeldiaggamh2 = modeldiaggamh[!(modeldiaggamh$spp %in% dropspp),]
		dim(modeldiaggamh2) # 326


	summary(modeldiaggamh2$auc)
		hist(modeldiaggamh2$auc)
		modeldiaggamh2[modeldiaggamh2$auc<0.7, c('region', 'spp', 'auc')]
	summary(modeldiaggamh2$r2.biomass)
		hist(modeldiaggamh2$r2.biomass)
	summary(modeldiaggamh2$r2.all)
		hist(modeldiaggamh2$r2.all)
	summary(modeldiaggamh2$dev.pres)
		hist(modeldiaggamh2$dev.pres)
		modeldiaggamh2[modeldiaggamh2$dev.pres<0.1, c('region', 'spp', 'dev.pres')]
	summary(modeldiaggamh2$dev.biomass)
		hist(modeldiaggamh2$dev.biomass)
		modeldiaggamh2[modeldiaggamh2$dev.biomass<0.1, c('region', 'spp', 'dev.biomass')]




