# Fit climate-envelope models
## Still very much in progress

#################
### Prep data ###	
#################
load('Output/trawl_allregionsforprojections_2015-02-02.RData') # load dat data.frame. Has all trawl observations from all regions. wtcpue has the standardized biomass estimates. They are standardized within regions, but not across regions.

# need to standardize species names across regions
# see spptaxonomy_2015-02-09_plusManual.csv for a useful conversion table
# species names aren't standardized across regions


## Pick the spp to model (WE NEED TO FIGURE OUT HOW WE WANT TO DO THIS)
## A good place to start would be to drop all taxa not identified to species? Perhaps also a minimum #years, or obs/year cutoff?
# At one point, I only used spp with at least 300 valid rows of data (present or not) and at least 40 rows of data where present
spp=sort(unique(dat$spp));spp 
sppreg = sort(unique(dat$sppregion))
# note: some species are found in both the Atlantic and the Pacific. We probably want to treat these separately? (as separate "species")
# then need to trim dat to just these species


# Add useful columns
	# Have a biomass that never goes to zero (useful for fitting log-links) 
dat$wtcpuena = dat$cpue
dat$wtcpuena[dat$wtcpuena == 0] = 1e-4
dat$wtcpuenal = log(dat$wtcpuena)

	# other useful columns
dat$stratumfact = as.factor(dat$stratum)
dat$yrfact = as.factor(dat$year)
dat$sppregion = paste(dat$spp, dat$region, sep='_')



# also need to add biomassmean column: mean biomass for a species in a region in a year




## NEED TO EXPAND DAT TO INCLUDE OBSERVATIONS OF ZERO for every species at every tow location in regions where they are found



# **IF** we want to fit stratum effects, we need this:
# Create a column that is an indicator for where a species is present, but also has at least one TRUE value per spp per stratum in each region, for fitting stratum effects. Randomly pick a row to turn into a "presence" even though not actually present.
	dat$presfit = dat$wtcpue > 0 # indicator for where present, but we'll tweak this in the code below
	presfitspp = character(0) # keeps a record of which species in which regions we had to tweak with false presences
	for(i in 1:length(regs)){ # for each region
		strat = sort(unique(dat$stratum[dat$region == regs[i]])) # get list of all strata in this region
		thesespp = sort(unique(dat$spp[dat$region == regs[i]])) # get list of all species in this region
		for(j in 1:length(thesespp)){ # for each species
			for(k in 1:length(strat)){ # for each stratum
				inds = dat$region == regs[i] & dat$spp == thesespp[j] & dat$stratum==strat[k] & complete.cases(dat[,c('surftemp', 'bottemp')]) # index into rows that we could use for fitting. it has to have complete explanatory variables, otherwise useless for model fitting. EXPAND THIS TO INCLUDE ALL EXPLANATORY VARIABLES
				if(sum(dat$presfit[inds]) == 0){ # if the rows for this species in this stratum have no presence observations...
					print(paste(regs[i], thesespp[j],strat[k]))
					m = sample(which(inds),1) # pick a random record in this stratum
					dat$presfit[m] = TRUE # set presfit to TRUE (even though wtcpue==0)
					print(paste(dat$presfit[m], dat$wtcpue[m], dat$wtcpuena[m], dat$year[m]))
					presfitspp = c(presfitspp, as.character(dat$sppregion[m]))
				}
			}
		}
	}
	sum(dat$presfit & !dat$pres) # how many rows were tweaked?
	presfitspp = unique(presfitspp); presfitspp # which species in which regions?


# write out the new dat
save()

##################
### Fit models ###	
##################

require(mgcv)

# Read in the trimmed, zero-padded dat from before
load()

spp = sort(unique(dat$spp))
for(i in 1:length(spp)){

	# Subset to rows in dat useful for fitting this model. DOESN'T YET INCLUDE ALL EXPLANATORY VARIABLES
	# also: Newfoundland and West Coast Annual don't have surftemp. How to handle this?
	inds = complete.cases(dat[,c('surftemp', 'bottemp')]) & dat$spp == spp[i]
	
	print(paste(spp[i], ': ', sum(dat$pres[inds]), ' presences, ', sum(inds), ' indices in training set', sep=''))
			
	print('fitting GAM presence binomial')	
	# note: this will fail on Newfoundland and West Coast Annual surveys, where surftemp isn't available
	# note: doesn't yet include all explanatory variables (like rugosity)
	# uses mgcv's built-in model selection function.
	gamhp2 = gam(pres ~ s(surftemp) + s(bottemp) + stratumfact + biomassmean, family=binomial, data=dat[inds,], select=TRUE)  

	print('fitting GAM abundance gaussian(log)')
	try3 = try(
		# note: this will fail on Newfoundland and West Coast Annual surveys, where surftemp isn't available
		# note: doesn't yet include all explanatory variables (like rugosity)
		# mgcv has a built-in model selection function
		gamha2 = gam(wtcpuenal ~ s(surftemp) + s(bottemp) + stratumfact + biomassmean, family=gaussian, data=dat[inds & dat$presfit,], select=TRUE, control=list(mgcv.half=30)) # includes a control parameter that seemed to help with fitting
	}
	)
	
	# Save models if model fitting worked. Currently set up to save a model for each taxon in each region (needs to be adjusted)
	if(class(try3)[1] != "try-error"){
		mods = list(gamhp2=gamhp2, gamha2 = gamha2)
		save(mods, inds, file=paste('Output/CEModels/CEmods_', regs[j], '_', spp[i], '_', traintype, '_', Sys.Date(), '.RData', sep=''))
	}
}


######################################
######## Basic model analysis ########
######################################

### CODE HERE IS FRAGMENTARY, BUT MAY BE HELPFUL?
	
## Model diagnostics: calculate ROC and r2 for each model
	require(ROCR)
	folder = 'CEModels 2012-09-05 allyrs onlytemp'; trainingyrs = yrs; vars=c('surftemp', 'bottemp'); traintype='allyrs'; suff='onlytemp' # all years w/ only temp (surf, bot, strat, biomass)

	files=list.files(paste('Output/', folder, sep=''), pattern='.RData') # only works if only one round of CEM models are in the directory
	regs = sort(unique(dat$region))
	sppregs = sort(unique(dat$sppregion))
	n = rep(NA, length(sppregs))
	
	# save ROC and r2 for each taxon in each region in a data.frame
	modeldiaggamh = data.frame(sppregion = sppregs, region = n, spp = n, auc = n, r2.biomass = n, r2.all = n, dev.pres = n, dev.biomass = n, stringsAsFactors=FALSE)
		
	options(warn=1) # print warnings as they occur
	fileinds = 1:length(files)
	for(s in fileinds){
		print(s)
		load(paste('Output/',folder, '/', files[s], sep=''))
		thisspp = strsplit(files[s], split=paste('AFSC_Aleutians_|AFSC_EBS_|AFSC_GOA_|DFO_ScotianShelf_|DFO_SoGulf_|NEFSC_Spring_|WestCoast_Tri_|DFO_Newfoundland_Fall_|SEFSC_GOMex_|_',traintype, sep=''))[[1]][2] # extract spp name from file name
		thisreg = regs[sapply(regs, FUN=grepl, x=files[s])]
		theseyrs = yrs[[thisreg]]
		j = intersect(grep(thisspp, modeldiaggamh$sppregion, fixed=TRUE), grep(thisreg, modeldiaggamh$sppregion))
		if(length(j)>1 | length(j)<1) stop('Found too few or too many indices into modeldiag dataframe')
		modeldiaggamh$spp[j] = thisspp
		modeldiaggamh$region[j] = as.character(thisreg)
	
		# Find observations and predictions
		inds = complete.cases(dat[,vars]) & dat$spp == thisspp & dat$region == thisreg
		if(thisreg == 'DFO_Newfoundland_Fall') inds = complete.cases(dat[,vars[vars!='surftemp']]) & dat$spp == thisspp & dat$region == thisreg

		if(sum(inds)>0){
			preds1 = predict(mods$gamhp2, newdata = dat[inds,], type='response')
			preds2 = exp(predict(mods$gamha2, newdata = dat[inds,], type='response')) # make sure to use exp if wtcpue was log-transformed!
			#plot(preds1, preds2, main=paste(spp[s], yrs[i]), xlab='Predicted presence', ylab='Predicted wtcpue')
			preds = preds1*preds2
			preds[preds<0] = 0

			# calculate performance (in part using ROCR)
			preds1.rocr = prediction(predictions=as.numeric(preds1), labels=dat$pres[inds])
				# plot the ROC curve
#				perf = performance(preds1.rocr, 'tpr', 'fpr')
#				plot(perf)
			modeldiaggamh$auc[j] = performance(preds1.rocr, 'auc')@y.values[[1]]
			modeldiaggamh$r2.biomass[j] = cor(preds2[dat$presfit[inds]], dat$wtcpue[inds & dat$presfit])^2 # correlation of biomass where present
			modeldiaggamh$r2.all[j] = cor(preds, dat$wtcpue[inds])^2 # overall biomass correlation
			modeldiaggamh$dev.pres[j] = summary(mods$gamhp)$dev.expl
			modeldiaggamh$dev.biomass[j] = summary(mods$gamha)$dev.expl
		}
	}
	
	### Save model diagnostics
	write.csv(modeldiaggamh, file = paste('Output/modeldiaggamh_', traintype, suff, '_', Sys.Date(), '.csv', sep=''))

