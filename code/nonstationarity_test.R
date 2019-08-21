## Script to test whether thermal envelopes change through time (non-stationarity)
## This is part 1: fits models to beginning and end of dataset

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




# Loop through species and test models.
library(mgcv);library(ROCR)


load("data/dat_selectedspp.Rdata")
	dat$logwtcpue <- log(dat$wtcpue)
	dat$survey <- dat$region #keeps NEFSC as two separate surveys in "survey"
	dat$surveyfact <-as.factor(dat$survey)
	dat$region[dat$region %in% c("NEFSC_NEUSFall","NEFSC_NEUSSpring")] <- "NEFSC_NEUS" #NEFSC as one region in "region"
	dat$regionfact <- as.factor(dat$region) # the version ot use for model fitting

allspp = sort(unique(dat$sppocean))
n = rep(NA, length(allspp))
modeldiag = data.frame(sppocean=n, ntot=n, ntot.1=n, ntot.2=n, npres=n, npres.1=n, npres.2=n, aic1=n, aic1null=n, aic2=n, aic2null=n, anova.resid.df.1=n, anova.resid.df.1null=n, anova.resid.df.2=n, anova.resid.df.2null=n, anova.resid.dev.1=n, anova.resid.dev.1null=n, anova.resid.dev.2=n, anova.resid.dev.2null=n, anova.p.1=n, anova.F.2=n, anova.p.2=n, stringsAsFactors=FALSE) # tt is for training/testing model



######################
# Start the big loop #
######################

#Open pdf to print figures 
pdf(file=paste("figures/GAMs_nonstationarity_test.pdf",sep=""),width=8,height=6)

options(warn=1) # print warnings as they occur
allwarnings = NULL
print(paste(length(allspp), 'models to fit'))

for(i in 1:length(allspp)){ 
	fittrain = TRUE
	mygam1null <- mygam2null <- mygam1 <- mygam2 <- NULL

	sp<-allspp[i]
	print(paste(i,sp, Sys.time()))

	mydat<-dat[dat$sppocean==sp,] 
	mydat$logwtcpue[is.infinite(mydat$logwtcpue)] <- NA
	myregions<-unique(mydat$region)
	myhauls<-unique(mydat$haulid)
	myocean <- unique(mydat$ocean)




	####################################################
	# Add records for when there was a haul but no fish
	####################################################
	# without season. expand to all regions in the same ocean.
	zs<-dat[!(dat$haulid %in% myhauls) & dat$ocean %in% myocean,] #extract haulids where this species is missing

	matchpos<-match(unique(zs$haulid),zs$haulid) # Extract one record of each haulid
	zeros<-zs[matchpos,]
	#Add/change relevant columns --> zero catch for target spp.
	zeros$spp<-mydat$spp[1]
	zeros$sppl<-mydat$sppl[1]
	zeros$sppnew<-mydat$sppnew[1]
	zeros$sppocean<-mydat$sppocean[1] #may need to add "ocean"
	zeros$wtcpue<-0
	zeros$logwtcpue <- NA
	zeros$presfit<-FALSE

	mydatw0<-rbind(mydat,zeros) #combine positive hauls with zero hauls

	##########################################################
	# Calculate mean catch per haul by region for this taxon #
	##########################################################

	# (For true biomass estimate, would need to stratify the mean)
	ave.catch.wt<-tapply(mydatw0$wtcpue,list(mydatw0$year,mydatw0$region),mean,na.rm=T)
	# transform into a data.frame
	if(dim(ave.catch.wt)[2]<2) { # only one region
		avecatchyrreg<-cbind(as.data.frame(ave.catch.wt), rep(colnames(ave.catch.wt), dim(ave.catch.wt)[2]), rownames(ave.catch.wt))
		colnames(avecatchyrreg)<-c("biomassmean","region","year")
	} else { # multiple regions
		avecatchyrreg<-cbind(stack(as.data.frame(ave.catch.wt)), rep(rownames(ave.catch.wt), dim(ave.catch.wt)[2]))
		colnames(avecatchyrreg)<-c("biomassmean","region","year")
	}
	spdata<-merge(mydatw0,avecatchyrreg)
	

	####################################################
	# Trim data to complete cases
	####################################################

	spdata<-spdata[complete.cases(spdata[,c("surftemp","bottemp","rugosity","presfit")]),]
	spdata <- droplevels(spdata)


	##############################################
	# Set up data in regions with no presences
	##############################################
	spdata$logwtcpue.pad <- spdata$logwtcpue # has some zeros changed to -18 to allow abundance model fitting
	spdata$presfit.pad <- spdata$presfit # has some FALSE changed to TRUE to allow abundance model fitting


	npres <- table(spdata$regionfact[spdata$presfit])
	if(any(npres < 1)){
		mywarn <- paste('Zero presences for', i, sp, 'in', paste(names(npres)[npres<1], collapse=', '), 'so adding some as -23 (1e-10) to allow abundance model fitting')
		allwarnings <- c(allwarnings, mywarn)
		warning(mywarn)
		regstofill <- names(npres)[as.numeric(npres) == 0]

		spdata$regionfact[spdata$region %in% regstofill] <- names(npres)[which.max(npres)] # in regions with no observations, replace the region ID with that from a region with observations. this prevents a low region coefficient from explaining the zero observations.

		# if a region has no presences
		# pick some absences and force them to low abundance presences for abundance model fitting
		for(j in 1:length(regstofill)){
			theseinds <- spdata$region == regstofill[j]
			fake0s <- sample(which(theseinds), size = round(0.1 * sum(theseinds)))
			spdata$logwtcpue.pad[fake0s] <- -23
			spdata$presfit.pad[fake0s] <- TRUE
			print(paste(regstofill[j], ': Added', length(fake0s), 'fake zeros'))
		}
	}
	
	spdata <- droplevels(spdata)


	###############################################################
	#Set up data for early and late to test for non-stationarity
	###############################################################

	#Subset training and testing data by year (use first 50% and last 50%)
	spdata<-spdata[order(spdata$year,spdata$month),]

	# indices for both pres and abs
	ninds<-table(spdata$regionfact) # number of entries per region (regions as set up for fitting)
	spdata$half <- NA
	for(j in 1:length(ninds)){ # loop through each region to get first 50% and last 50%
		spdata$half[which(as.character(spdata$regionfact) == names(ninds)[j])[1:round(ninds[j]*0.5)]] <- 1 # first half
		spdata$half[which(as.character(spdata$regionfact) == names(ninds)[j])[(round(ninds[j]*0.5)+1):ninds[j]]] <- 2 # second half
	}
	spdata$half <- as.factor(spdata$half)
	if(sum(is.na(spdata$half))>0){ # should assign every row to first or second half
		mywarn <- paste('Some first/second half values undefined for', i, sp)
		allwarnings <- c(allwarnings, mywarn)
		warning(mywarn)
	}
	
	# indices for only where present (for the abundance model), including fake zeros
	firstindsp <- which(spdata$presfit.pad & spdata$half==1)
	lastindsp <- which(spdata$presfit.pad & spdata$half==2)

	# warn if too few presences overall
	if(length(firstindsp)<2){
		mywarn <- paste('Only', length(firstindsp), 'presence values in first half of dataset for', i, sp)
		allwarnings <- c(allwarnings, mywarn)
		warning(mywarn)
	}
	if(length(lastindsp)<2){
		mywarn <- paste('Only', length(lastindsp), 'presence values in second half of dataset for', i, sp)
		allwarnings <- c(allwarnings, mywarn)
		warning(mywarn)
	}

	# test if we have enough presences in first and last half sets (at least one per region as set up for model fitting)
	npresfirst <- table(spdata$regionfact[firstindsp])
	npreslast <- table(spdata$regionfact[lastindsp])
	if(any(npresfirst < 1)){
		mywarn <- paste('Zero training presences for', i, sp, 'in', paste(names(npresfirst)[npresfirst<1], collapse=', '))
		allwarnings <- c(allwarnings, mywarn)
		warning(mywarn)
		regstofill <- names(npresfirst)[as.numeric(npresfirst) == 0]
	}
	if(any(npreslast < 1)){
		mywarn <- paste('Zero testing presences for', i, sp, 'in', paste(names(npreslast)[npreslast<1], collapse=', '))
		allwarnings <- c(allwarnings, mywarn)
		warning(mywarn)
	}
	
	# make sure we have at least 6 unique levels for each variable (necessary to fit gam with 4 knots)
	levs1 <- apply(spdata[firstindsp,c('bottemp', 'surftemp', 'logrugosity')], 2, FUN=function(x) length(unique(x)))
	levs2 <- apply(spdata[lastindsp,c('bottemp', 'surftemp', 'logrugosity')], 2, FUN=function(x) length(unique(x)))
	if(any(levs1 < 6)){
		mywarn <- paste("Not enough (>=6) unique levels in first half presence set for", i, sp, ". Won't fit models")
		allwarnings <- c(allwarnings, mywarn)
		warning(mywarn)
		fittrain = FALSE
	}	
	if(any(levs2 < 6)){
		mywarn <- paste("Not enough (>=6) unique levels in second half presence set for", i, sp, ". Won't fit models")
		allwarnings <- c(allwarnings, mywarn)
		warning(mywarn)
		fittrain = FALSE
	}	
	# table(spdata$year[spdata$presfit]) # output number of presences by year
	
	####################################################
	# Figure out which model formula given data
	####################################################

	#Default models. Leave out region factor if necessary
	# null model doesn't include first/half distinction: NO NEED TO FIT THIS (also don't need to fit abunmod)
	# could simplify presmod so it doesn't have an interaction on logrugosity or biomassmean
	if(length(levels(spdata$regionfact))==1){
			mypresmod<-formula(presfit ~ s(bottemp,k=6,by=half)+s(surftemp,k=6,by=half)+s(logrugosity,k=4,by=half)+biomassmean*half)
			myabunmod<-formula(logwtcpue.pad ~ s(bottemp,k=6,by=half)+s(surftemp,k=6,by=half)+s(logrugosity,k=4,by=half)+biomassmean*half)
			mynullpresmod<-formula(presfit ~ s(bottemp,k=6)+s(surftemp,k=6)+s(logrugosity,k=4)+biomassmean) # null model not splitting by first/second half
			mynullabunmod<-formula(logwtcpue.pad ~ s(bottemp,k=6)+s(surftemp,k=6)+s(logrugosity,k=4)+biomassmean) # null model not splitting by first/second half
	} else {
			mypresmod<-formula(presfit ~ s(bottemp,k=6,by=half)+s(surftemp,k=6,by=half)+s(logrugosity,k=4,by=half)+regionfact+biomassmean*half-1)
			myabunmod<-formula(logwtcpue.pad ~ s(bottemp,k=6,by=half)+s(surftemp,k=6,by=half)+s(logrugosity,k=4,by=half)+regionfact+biomassmean*half-1)
			mynullpresmod<-formula(presfit ~ s(bottemp,k=6)+s(surftemp,k=6)+s(logrugosity,k=4)+regionfact+biomassmean-1) # null model not splitting by first/second half
			mynullabunmod<-formula(logwtcpue.pad ~ s(bottemp,k=6)+s(surftemp,k=6)+s(logrugosity,k=4)+regionfact+biomassmean-1) # null model not splitting by first/second half
	}



	####################################################
	#Fit models to All data
	####################################################

	if(fittrain){
		try2 <- tryCatch({
			mygam1<-gam(mypresmod,family="binomial",data=spdata, control=list(maxit=500, nthreads=nthreads)) # add more iterations to help convergence, plus use more threads
		}, error = function(e) {
			mywarn <- paste('Error in gam fitting for', i, sp, ':', e)
			assign('allwarnings', c(get('allwarnings', envir=.GlobalEnv), mywarn), envir=.GlobalEnv) # these assigns outside the local scope are poor form in R. But not sure how else to do it here...
			warning(mywarn)
		})
		try3 <- tryCatch({
			mygam2<-gam(myabunmod,data=spdata[spdata$presfit.pad,], control=list(maxit=500, nthreads=nthreads)) # only fit where spp is present
		}, error = function(e) {
			mywarn <- paste('Error in gam fitting for', i, sp, ':', e)
			assign('allwarnings', c(get('allwarnings', envir=.GlobalEnv), mywarn), envir=.GlobalEnv) # these assigns outside the local scope are poor form in R. But not sure how else to do it here...
			warning(mywarn)
		})
		try4 <- tryCatch({
			mygam1null<-gam(mynullpresmod,family="binomial",data=spdata, control=list(maxit=500, nthreads=nthreads))
		}, error = function(e) {
			mywarn <- paste('Error in gam fitting for', i, sp, ':', e)
			assign('allwarnings', c(get('allwarnings', envir=.GlobalEnv), mywarn), envir=.GlobalEnv) # these assigns outside the local scope are poor form in R. But not sure how else to do it here...
			warning(mywarn)
		})
		try5 <- tryCatch({
			mygam2null<-gam(mynullabunmod,data=spdata[spdata$presfit.pad,], control=list(maxit=500, nthreads=nthreads)) # only fit where spp is present
		}, error = function(e) {
			mywarn <- paste('Error in gam fitting for', i, sp, ':', e)
			assign('allwarnings', c(get('allwarnings', envir=.GlobalEnv), mywarn), envir=.GlobalEnv) # these assigns outside the local scope are poor form in R. But not sure how else to do it here...
			warning(mywarn)
		})


		####################################################
		# Plot gam smooths to check visually
		####################################################
	
		#Should write out to PDF opened before the loop
		plot(mygam1,pages=1,scale=0,all.terms=TRUE);mtext(paste(sp,"presence"),outer=T,line=-2)
		plot(mygam2,pages=1,scale=0,all.terms=TRUE);mtext(paste(sp,"abundance"),outer=T,line=-2)
		plot(mygam1null,pages=1,scale=0,all.terms=TRUE);mtext(paste(sp,"presence Null"),outer=T,line=-2)
		plot(mygam2null,pages=1,scale=0,all.terms=TRUE);mtext(paste(sp,"abundance Null"),outer=T,line=-2)


		#################
		# Compare models
		#################

		# fill in basic model data
		modeldiag$sppocean[i] = sp
		modeldiag$npres[i] = sum(spdata$presfit)
		modeldiag$npres.1[i] <- sum(spdata$presfit[spdata$half=='1'])
		modeldiag$npres.2[i] <- sum(spdata$presfit[spdata$half=='2'])

		modeldiag$ntot[i] = dim(spdata)[1]
		modeldiag$ntot.1[i] = sum(spdata$half=='1')
		modeldiag$ntot.2[i] = sum(spdata$half=='2')

		# compare AIC
		modeldiag$aic1[i] <- AIC(mygam1)
		modeldiag$aic1null[i] <- AIC(mygam1null)
		modeldiag$aic2[i] <- AIC(mygam2)
		modeldiag$aic2null[i] <- AIC(mygam2null)
	
		# anova comparison
		av1 <- anova(mygam1, mygam1null, test='Chisq') # chisq since binomial errors (see ?anova.glm)
		av2 <- anova(mygam2, mygam2null, test='F') # F since gaussian (see ?anova.glm)
	
		modeldiag$anova.resid.df.1[i] <- av1$'Resid. Df'[1]
		modeldiag$anova.resid.df.1null[i] <- av1$'Resid. Df'[2]
		modeldiag$anova.resid.df.2[i] <- av2$'Resid. Df'[1]
		modeldiag$anova.resid.df.2null[i] <- av2$'Resid. Df'[2]
	
		modeldiag$anova.resid.dev.1[i] <- av1$'Resid. Dev'[1]
		modeldiag$anova.resid.dev.1null[i] <- av1$'Resid. Dev'[2]
		modeldiag$anova.resid.dev.2[i] <- av2$'Resid. Dev'[1]
		modeldiag$anova.resid.dev.2null[i] <- av2$'Resid. Dev'[2] 
	
		modeldiag$anova.p.1[i] <- av1$'Pr(>Chi)'[2]
		modeldiag$anova.F.2[i] <- av2$F[2]
		modeldiag$anova.p.2[i] <- av2$'Pr(>F)'[2]


		####################################################
		#### Save models for later projections
		####################################################

		mods = list(mygam1=mygam1, mygam2 = mygam2, mygam1null=mygam1null, mygam2null=mygam2null)
	
		sp <- gsub('/', '', sp) # would mess up saving the file
	
		save(mods, myregions, file=paste(modfolder, 'CEmods_nonstationaritytest_', sp, '.RData', sep='')) # ~4mb file

		# write these files each time through the loop so that we can watch progress
		save(modeldiag,file="output/modeldiag_nonstationaritytest.Rdata")
		write.csv(modeldiag, file="output/modeldiag_nonstationaritytest.csv")

		write.csv(allwarnings, file='output/warnings_nonstationaritytest.csv')
	}
}


dev.off() # close the figure

