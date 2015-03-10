# Fit climate-envelope models
## Still very much in progress

#################
### Prep data ###	
#################
load('data/trawl_allregionsforprojections_2015-02-02.RData') # load dat data.frame. Has all trawl observations from all regions. wtcpue has the standardized biomass estimates. They are standardized within regions, but not across regions.


# need to standardize species names across regions
# see spptaxonomy_2015-02-09_plusManual.csv for a useful conversion table
Spptax<-read.csv("data/spptaxonomy_2015-02-09_plusManual.csv") #note: new column in CSV file 
spptax<-apply(Spptax,2,tolower)
dat$sppl<-tolower(dat$spp)
datspp<-unique(dat$sppl)

spptax<-as.data.frame(spptax)
sum(datspp %in% spptax[,1])# 792 of 4937 spp matched

#Match when possible and assign new genus species, otherwise keep old name
#Find taxa with matches in the taxonomy table
#In some cases, the "name" field of taxonomy has spelling errors (missing "n"s). In this case, genus+species should be used. But in some cases, "genus species" are not given if full taxonomy is missing. The spptaxonomy.csv file now has a "newname" column which merges these optimally.
matches<-datspp %in% spptax[,1] + datspp %in% spptax[,11] + datspp %in% spptax[,12] + datspp %in% spptax[,13]  #826 matches/4937
#For those species datspp[matches>0], replace current name with spptax$newname

newnames=NULL
for (i in 1:length(datspp)){
	if (matches[i]==0) newnames[i]<-datspp[i]
	else if(datspp[i] %in% spptax[,1]) 
	{
	ind<-match(datspp[i],spptax[,1])
	newnames[i]<-as.character(spptax[ind,14])
	} 
	else if(datspp[i] %in% spptax[,11]) 
	{
	ind<-match(datspp[i],spptax[,11])
	newnames[i]<-as.character(spptax[ind,14])
	}
	else if(datspp[i] %in% spptax[,13]) 
	{
	ind<-match(datspp[i],spptax[,13])
	newnames[i]<-as.character(spptax[ind,14])
	}
	else if(datspp[i] %in% spptax[,12]) 
	{
	ind<-match(datspp[i],spptax[,12])
	newnames[i]<-as.character(spptax[ind,14])
	}
}

oldnewnames<-cbind(sppl=datspp,sppnew=newnames)
#Now merge new names with old names in dat file
dat<-merge(dat,oldnewnames) #dat now contains "sppnew"


# Add useful columns
	# Have a biomass that never goes to zero (useful for fitting log-links) 
dat$wtcpuena = dat$wtcpue
dat$wtcpuena[dat$wtcpuena == 0] = 1e-4
dat$wtcpuenal = log(dat$wtcpuena)

	# other useful columns
dat$presfit = dat$wtcpue > 0 # indicator for where present
dat$stratumfact = as.factor(dat$stratum)
dat$yrfact = as.factor(dat$year)
dat$regionfact<-as.factor(dat$region)
dat$sppregion = paste(dat$sppnew, dat$region, sep='_')
dat$ocean[dat$region %in% c("AFSC_Aleutians", "AFSC_EBS", "AFSC_GOA", "AFSC_WCTri", "NWFSC_WCAnn")] <- "Pac"
dat$ocean[dat$region %in% c("DFO_NewfoundlandFall", "DFO_NewfoundlandSpring", "DFO_ScotianShelf","DFO_SoGulf","NEFSC_NEUSFall", "NEFSC_NEUSSpring")] <- "Atl"
dat$ocean[dat$region == "SEFSC_GOMex"] <- "Gulf" #Or should Gulf of Mex group with Altantic?
dat$sppocean = paste(dat$sppnew, dat$ocean, sep='_') 

	# add rugosity
rugfile<-read.csv("data/trawl_latlons_rugosity_forMalin_2015_02_10.csv")
dat<-merge(dat,rugfile) #lose 69 instances of lumpenus lampretaeformis b/c missing lat/lon
rm(rugfile)
dat$logrugosity<-log(dat$rugosity+0.01) #log-transformed rugosity gave better model fits in initial tests of rugosity covariate



## Pick the spp to model (WE NEED TO FIGURE OUT HOW WE WANT TO DO THIS)
## A good place to start would be to drop all taxa not identified to species? Perhaps also a minimum #years, or obs/year cutoff?

# Identify taxa caught at least once every survey year in a given survey (hopefully avoiding changes in species classification through time).
# Species are identified by species+ocean for later trimming of dat.
#Special code for DFO_SOGulf species - these include zero hauls, and some species were not caught every year. Others (chionoecetes opilio and homarus americanus) have NAs from 1971-1979.
#The other surveys have zero hauls, but not so systematically, may just be very small catches.

myregions<-unique(dat$region)
myspp<-NULL
for(r in myregions){
#	surveyyrs<-names(table(dat$year[dat$region==r]))
	regdat<-dat[dat$region==r,]
	if(r=="DFO_SoGulf") {
	 	maxcatch<-tapply(regdat$wtcpue,list(regdat$sppocean,regdat$year),max) 
	 	# check if maxcatch is always > 0	
		minofmax<-apply(maxcatch,1,min,na.rm=T) #includes spps with NA's early on. Check this.	
		min1<-names(minofmax[minofmax>0 & minofmax<Inf]) #taxa with >0 catch every year 
	}else{
	yrocc<-table(regdat$year,regdat$sppocean) #table of number of occurrances each year
	sumyrs<-apply(yrocc,2,function(x) sum(x==0)) #identify years with zero catch of taxa
	min1<-colnames(yrocc)[sumyrs==0]  #identify taxa with no years of zero catch (could relax this criterion to [sumyrs <= max0yrs], where max0yrs is maximum allowable number of years with zero catch)
	}
	myspp<-c(myspp,min1) #now with myspp plus ocean
}
#table(myspp) #shows which species are selected in multiple surveys
myspp<-unique(myspp) #length(myspp) 622 unique taxa_ocean

#remove:  any taxa missing species, egg cases, anemones missing species name

drop<-myspp[grep("spp",myspp)] #85 with "spp."
drop<-c(drop,myspp[grep("egg",myspp)]) #3 with "egg"
drop<-c(drop,myspp[grep("anemone",myspp)]) #2
drop<-c(drop,"teuthida","liparidinae","bathylagus sp.","lampanyctus sp.","caridea", "carinariidae","antipatharia","annelida" ,"crustacea shrimp")
myspp<-myspp[! myspp %in% drop]
myspp<-sort(myspp) #down to 522 spp. to model. A few species are present in multiple ocean basins and thus will have multiple models.


# At one point, I only used spp with at least 300 valid rows of data (present or not) and at least 40 rows of data where present
# spp=sort(unique(dat$spp));spp 
# sppreg = sort(unique(dat$sppregion))
# note: some species are found in both the Atlantic and the Pacific. We probably want to treat these separately? (as separate "species")
# then need to trim dat to just these species
# However, we'll need to go back to the full dat file to determine where zero hauls occur, unless we assume each haul had at least one of these species. Check this here:
#dat2<-dat[dat$sppocean %in% myspp,]
#length(unique(dat$haulid))  #114810
#length(unique(dat2$haulid)) #114640 (only 170 "zero" hauls not accounted for in dat2. meh.)

dat<-dat[dat$sppocean %in% myspp,] 


## NEED TO EXPAND DAT TO INCLUDE OBSERVATIONS OF ZERO for every species at every tow location in regions where they are found
## also need to add biomassmean column: mean biomass for a species in a region in a year (should be calculated after including zero hauls)



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

