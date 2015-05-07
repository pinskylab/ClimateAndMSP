# Loop through species and fit models.
library(mgcv);library(ROCR)

load("data/dat_selectedspp.Rdata")
	dat$logwtcpue <- log(dat$wtcpue)
	dat$regionfact <- as.factor(dat$region)

	# if using season
	dat$season <- c(rep('wi', 3), rep('sp', 3), rep('su', 3), rep('fa', 3))[dat$month]
	dat$regseas <- paste(dat$region, dat$season, sep='_')

allspp = sort(unique(dat$sppocean))
n = rep(NA, length(allspp))
modeldiag = data.frame(sppocean=n, npres=n, npres.tr=n, npres.te=n, ntot=n, auc=n, auc.tt=n, r2.biomass=n, r2.biomass.tt=n, r2.all=n, r2.all.tt=n, r2.pres.1deg=n, r2.abun.1deg=n, dev.pres=n, dev.biomass=n, stringsAsFactors=FALSE) # tt is for training/testing model


## small test to see which species are only marginally present in Newfoundland or WCAnn surveys
#nosurfregions <- c("DFO_NewfoundlandSpring","DFO_NewfoundlandFall","NWFSC_WCAnn")
#nsrspp <- sort(unique(dat$sppocean[dat$region %in% nosurfregions])) # spp present in at least one of these regions
#length(nsrspp)
#length(unique(dat$sppocean))
#tab <- table(dat$sppocean[dat$wtcpue>0 & dat$sppocean %in% nsrspp], dat$region[dat$wtcpue>0 & dat$sppocean %in% nsrspp]) # counts presences in all regions
#dim(tab)
#colnames(tab) <- c('AI', 'EBS', 'GOA', 'WCTri', 'Newf_F', 'Newf_S', 'Scot', 'SoGulf', 'NEUS_F', 'NEUS_S', 'WCAnn', 'GoMex') # shorter column names for easier printing to screen
#cmax <- apply(tab, 1, which.max) # region with the most catches
#nsrspp2 <- which(!(cmax %in% c(5,6,11))) # index for species that didn't have the highest presence count in Newf or WCAnn
#tab[nsrspp2,] # spp like Clupea harengus are common in many places


######################
# Start the big loop #
######################

options(warn=1) # print warnings as they occur
allwarnings = NULL
for(i in 347:length(allspp)){ 
	fittrain = TRUE
	mygam1tt <- mygam2tt <- mygam1 <- mygam2 <- preds <- preds1 <- preds2 <- predstt <- preds1tt <- preds2tt <- NULL

	sp<-allspp[i]
	print(paste(i,sp, Sys.time()))

	mydat<-dat[dat$sppocean==sp,] 
	myregions<-unique(mydat$region)
	myhauls<-unique(mydat$haulid)
	#myocean<-mydat[1,"ocean"]
	#myregseas<-unique(mydat$regseas)

	###################################################
	# Test whether we have enough data in each region #
	###################################################

	# if ignoring season in the models
	if(any(c("DFO_NewfoundlandSpring","DFO_NewfoundlandFall","NWFSC_WCAnn") %in% myregions)){
		ninds<-table(mydat$region[mydat$presfit & complete.cases(mydat[,c("bottemp","rugosity","presfit")])]) # number of presences per region with complete data (not inc. surftemp)
	} else {
		ninds<-table(mydat$region[mydat$presfit & complete.cases(mydat[,c("bottemp","surftemp", "rugosity","presfit")])]) # number of presences per region with complete data (inc. surftemp)
	}
	myregions <- names(ninds)[ninds >= 10] # require at least 10 presences with complete data to keep a region
	mydat <- mydat[mydat$region %in% myregions,]

	# if including season in the models
#	if(any(c("DFO_NewfoundlandSpring","DFO_NewfoundlandFall","NWFSC_WCAnn") %in% myregions)){
#		ninds<-table(mydat$regseas[mydat$presfit & complete.cases(mydat[,c("bottemp","rugosity","presfit")])]) # number of presences per region with complete data (not inc. surftemp)
#	} else {
#		ninds<-table(mydat$regseas[mydat$presfit & complete.cases(mydat[,c("bottemp","surftemp", "rugosity","presfit")])]) # number of presences per region with complete data (inc. surftemp)
#	}
#	myregseas <- names(ninds)[ninds >= 50] # require at least this many presences with complete data to keep a region/season combination (not sure if this is a reasonable level)
#	mydat <- mydat[mydat$regseas %in% myregseas,]


	####################################################
	# Add records for when there was a haul but no fish
	####################################################
	# Only expand catches (pad with zeros) in regions where the species has previously been caught.

	# without season
	zs<-dat[!(dat$haulid %in% myhauls) & dat$region %in% myregions,] #extract haulids where this species is missing, but only from regions where this species has at least one record.

	# with season
#	zs<-dat[!(dat$haulid %in% myhauls) & dat$regseas %in% myregseas,] #extract haulids where this species is missing, but only from regions where this species has at least one record.

	matchpos<-match(unique(zs$haulid),zs$haulid) # Extract one record of each haulid
	zeros<-zs[matchpos,]
	#Add/change relevant columns --> zero catch for target spp.
	zeros$spp<-mydat$spp[1]
	zeros$sppl<-mydat$sppl[1]
	zeros$sppnew<-mydat$sppnew[1]
#	zeros$sppregion<-paste(zeros$sppnew,zeros$region,sep="_")
	zeros$sppocean<-mydat$sppocean[1] #may need to add "ocean"
	zeros$wtcpue<-0
	zeros$presfit<-F
#	zeros$wtcpuena<-1e-4 #for consistency, but really zero
#	zeros$wtcpuenal<-log(zeros$wtcpuena)

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

	#What if only a few records exist for catch in one of these regions - worth excuding surftemp from model? Or rather exclude region?
	#Another option would be to fill in surftemp from HadISST or another SST product
	if(any(c("DFO_NewfoundlandSpring","DFO_NewfoundlandFall","NWFSC_WCAnn") %in% myregions)){
		spdata<-spdata[complete.cases(spdata[,c("bottemp","rugosity","presfit")]),] 
	} else {
		spdata<-spdata[complete.cases(spdata[,c("surftemp","bottemp","rugosity","presfit")]),] 
	}

	####################################################
	#Set up data for training and testing to evaluate performance
	####################################################

	#Subset training and testing data by year (use first 80% from each region to predict last 20% in each region)
	spdata<-spdata[order(spdata$year,spdata$month),]

	# indices for both pres and abs
	ninds<-table(spdata$region) # number of entries per region
	traininds <- NULL; testinds <- NULL
	for(j in 1:length(ninds)){ # loop through each region to get first 80% and last 20%
		traininds <- c(traininds, which(spdata$region == names(ninds)[j])[1:round(ninds[j]*0.8)])
		testinds <- c(testinds, which(spdata$region == names(ninds)[j])[(round(ninds[j]*0.8)+1):ninds[j]])
	}
	
	# indices for only where present (for the abundance model)
	trainindsp <- intersect(traininds, which(spdata$presfit))
	testindsp <- intersect(testinds, which(spdata$presfit))

	if(length(trainindsp)<2){
		mywarn <- paste('Only', length(trainindsp), 'presence values in testing dataset for', i, sp)
		allwarnings <- c(allwarnings, mywarn)
		warning(mywarn)
	}
	if(length(testindsp)<2){
		mywarn <- paste('Only', length(testindsp), 'presence values in testing dataset for', i, sp)
		allwarnings <- c(allwarnings, mywarn)
		warning(mywarn)
	}

	# test if we have enough presences in testing and training sets (at least one per region)
	nprestrain <- table(spdata$region[trainindsp])
	nprestest <- table(spdata$region[testindsp])
	if(any(nprestrain < 1)){
		mywarn <- paste('Zero training presences for', i, sp, 'in', names(nprestrain)[nprestrain<1])
		allwarnings <- c(allwarnings, mywarn)
		warning(mywarn)
	}
	if(any(nprestest < 1)){
		mywarn <- paste('Zero testing presences for', i, sp, 'in', names(nprestest)[nprestest<1])
		allwarnings <- c(allwarnings, mywarn)
		warning(mywarn)
	}
	
	# make sure we have at least 4 unique levels for each variable (necessary to fit gam with 4 knots)
	# look at training presence indices, since the most constraining (for mygam2tt)
	# this doesn't test by season... and so wont' catch a season with very few datapoints
	levs <- apply(spdata[trainindsp,c('bottemp', 'surftemp', 'logrugosity', 'biomassmean')], 2, FUN=function(x) length(unique(x)))
	if(any(c("DFO_NewfoundlandSpring","DFO_NewfoundlandFall","NWFSC_WCAnn") %in% myregions)){
		levs <- apply(spdata[trainindsp,c('bottemp', 'logrugosity', 'biomassmean')], 2, FUN=function(x) length(unique(x)))
	}
	if(any(levs < 4)){
		mywarn <- paste("Not enough (>=4) unique levels in training presence set for", i, sp, ". Won't fit training models")
		allwarnings <- c(allwarnings, mywarn)
		warning(mywarn)
		fittrain = FALSE
	}	
	# table(spdata$year[spdata$presfit]) # output number of presences by year
	
	####################################################
	# Figure out which model formula given data
	####################################################

	# without season
	#Default models. Leave out region factor if necessary
	if(length(myregions)==1){
			mypresmod<-formula(presfit ~ s(bottemp,k=4)+s(surftemp,k=4)+s(logrugosity,k=4)+s(biomassmean,k=4))
			myabunmod<-formula(logwtcpue ~ s(bottemp,k=4)+s(surftemp,k=4)+s(logrugosity,k=4)+s(biomassmean,k=4))
	} else {
			mypresmod<-formula(presfit ~ s(bottemp,k=4)+s(surftemp,k=4)+s(logrugosity,k=4)+regionfact+s(biomassmean,k=4)-1)
			myabunmod<-formula(logwtcpue ~ s(bottemp,k=4)+s(surftemp,k=4)+s(logrugosity,k=4)+regionfact+s(biomassmean,k=4)-1)
	}
	#Replace with formula missing surftemp if necessary
	if(any(c("DFO_NewfoundlandSpring","DFO_NewfoundlandFall","NWFSC_WCAnn") %in% myregions)){
		if(length(myregions)==1){
			mypresmod<-formula(presfit~s(bottemp,k=4)+s(logrugosity,k=4)+s(biomassmean,k=4))
			myabunmod<-formula(logwtcpue~s(bottemp,k=4)+s(logrugosity,k=4)+s(biomassmean,k=4))
		} else {
			mypresmod<-formula(presfit~s(bottemp,k=4)+s(logrugosity,k=4)+regionfact+s(biomassmean,k=4)-1)
			myabunmod<-formula(logwtcpue~s(bottemp,k=4)+s(logrugosity,k=4)+regionfact+s(biomassmean,k=4)-1)
		}
	}

	# with season
	#Default models. Leave out region factor if necessary
#	if(length(myregions)==1){
#			mypresmod<-formula(presfit ~ s(bottemp,k=4,by=season) +s(surftemp,k=4,by=season) +s(logrugosity,k=4) +s(biomassmean,k=4))
#			myabunmod<-formula(logwtcpue ~ s(bottemp,k=4,by=season) +s(surftemp,k=4,by=season) +s(logrugosity,k=4) +s(biomassmean,k=4))
#	} else {
#			mypresmod<-formula(presfit ~ s(bottemp,k=4,by=season) +s(surftemp,k=4,by=season) +s(logrugosity,k=4) +regionfact +s(biomassmean,k=4)-1)
#			myabunmod<-formula(logwtcpue ~ s(bottemp,k=4,by=season) +s(surftemp,k=4,by=season) +s(logrugosity,k=4) +regionfact +s(biomassmean,k=4)-1)
#	}
#	#Replace with formula missing surftemp if necessary
#	if(any(c("DFO_NewfoundlandSpring","DFO_NewfoundlandFall","NWFSC_WCAnn") %in% myregions)){
#		if(length(myregions)==1){
#			mypresmod<-formula(presfit~s(bottemp,k=4,by=season) +s(logrugosity,k=4) +s(biomassmean,k=4))
#			myabunmod<-formula(logwtcpue~s(bottemp,k=4,by=season) +s(logrugosity,k=4)+s(biomassmean,k=4))
#		} else {
#			mypresmod<-formula(presfit~s(bottemp,k=4,by=season) +s(logrugosity,k=4)+regionfact+s(biomassmean,k=4)-1)
#			myabunmod<-formula(logwtcpue~s(bottemp,k=4,by=season) +s(logrugosity,k=4)+regionfact+s(biomassmean,k=4)-1)
#		}
#	}


	####################################	
	# Fit the training/testing models
	####################################	

	# We could use select=TRUE so that terms can be smoothed out of the model (a model selection algorithm)
	if(fittrain){
		try1 <- tryCatch({
			mygam1tt<-gam(mypresmod, family="binomial",data=spdata[traininds,]) 
			mygam2tt<-gam(myabunmod, data=spdata[trainindsp,]) # only fit where species is present
		}, error = function(e) { # ignore warnings, since no function to catch them
			mywarn <- paste('Error in training gam fitting for', i, sp, ':', e)
			allwarnings <- c(allwarnings, mywarn)
			warning(mywarn)
			assign('fittrain', FALSE, envir=.GlobalEnv) # if we hit an error in predictions, we can't calculate performance stats
		})
	}

			

	####################################################
	#Fit models to All data (no test/training split)
	####################################################

	# We could use select=TRUE so that terms can be smoothed out of the model (a model selection algorithm)
	try2 <- tryCatch({
		mygam1<-gam(mypresmod,family="binomial",data=spdata)
		mygam2<-gam(myabunmod,data=spdata[spdata$presfit,]) # only fit where spp is present
	}, error = function(e) {
		mywarn <- paste('Error in gam fitting for', i, sp, ':', e)
		allwarnings <- c(allwarnings, mywarn)
		warning(mywarn)
	})

	####################################################
	# Compare predictions to observations to assess model performance
	####################################################

	# For FULL model
	preds1 <- predict(mygam1,spdata,type="response") #can also use mygam1$fitted.values
	preds2 <- exp(predict(mygam2, newdata = spdata, type='response')) # abundance predictions
	smear = mean(exp(mygam2$residuals)) # smearing estimator for re-transformation bias (see Duan 1983, http://www.herc.research.va.gov/include/page.asp?ID=cost-regression)
	preds <- preds1*preds2*smear # adds the bias correction as well
	preds[preds<0] = 0

	# And for training/testing data set
	if(fittrain){
		try3 <- tryCatch({
			preds1tt <- predict(mygam1tt,spdata[testinds,],type="response") 
			preds2tt <- exp(predict(mygam2tt, newdata = spdata[testinds,], type='response'))
			smear = mean(exp(mygam2tt$residuals)) # smearing estimator for re-transformation bias (see Duan 1983, http://www.herc.research.va.gov/include/page.asp?ID=cost-regression)
			predstt <- preds1tt*preds2tt*smear
			predstt[predstt<0] = 0
		}, error = function(e) {
			assign('fittrain', FALSE, envir=.GlobalEnv) # if we hit an error in predictions, we can't calculate performance stats
			mywarn <- paste('Error in predicting to test data for', i, sp, ':', e)
			allwarnings <- c(allwarnings, mywarn)
			warning(mywarn)
		})
	}

	# fill in diagnostics
	modeldiag$sppocean[i] = sp
	modeldiag$npres[i] = sum(spdata$presfit)
	if(fittrain){
		modeldiag$npres.tr[i] = sum(spdata$presfit[traininds])
		modeldiag$npres.te[i] = sum(spdata$presfit[testinds])
	}
	modeldiag$ntot[i] = dim(spdata)[1]

	# calculate performance (in part using ROCR)
	preds1.rocr = prediction(predictions=as.numeric(preds1), labels=spdata$presfit)
	modeldiag$auc[i] = performance(preds1.rocr, 'auc')@y.values[[1]] # area under the ROC curve
	if(length(testindsp)>0 & fittrain){ # need presences in the test dataset
		preds1tt.rocr = prediction(predictions=as.numeric(preds1tt), labels=spdata$presfit[testinds])
		modeldiag$auc.tt[i] = performance(preds1tt.rocr, 'auc')@y.values[[1]] #
	}
	# could add true skill statistic

	modeldiag$r2.biomass[i] = cor(log(preds2[spdata$presfit]), spdata$logwtcpue[spdata$presfit])^2 # correlation of log(biomass) where present
	if(length(testindsp)>0 & fittrain) modeldiag$r2.biomass.tt[i] = cor(preds2tt[which(testinds %in% testindsp)], spdata$logwtcpue[testindsp])^2 # only if presences exist in the test dataset
	modeldiag$r2.all[i] = cor(preds, spdata$wtcpue)^2 # overall biomass correlation
	if(length(testindsp)>0 & fittrain) modeldiag$r2.all.tt[i] = cor(predstt, spdata$wtcpue[testinds])^2 # overall biomass correlation. only makes sense to do this if the species is present at least once in the testing dataset
	modeldiag$dev.pres[i] = summary(mygam1)$dev.expl
	modeldiag$dev.biomass[i] = summary(mygam2)$dev.expl

	# Some metrics at a spatially aggregated level (1x1deg square) (by year) may be more informative:
	test<-cbind(spdata,preds1,preds)
	t1<-tapply(test$preds1,list(test$year,test$cs1),mean) #average predicted p(occur)
	t2<-tapply(test$presfit,list(test$year,test$cs1),mean) #proportion of hauls with presence
	t3<-tapply(test$preds,list(test$year,test$cs1),mean) #average predicted abundance
	t4<-tapply(test$wtcpue,list(test$year,test$cs1),mean) #average observed abundance

	presr2<-round(cor(stack(as.data.frame(t2))[,1],stack(as.data.frame(t1))[,1],use="p")^2,2)
	abunr2<-round(cor(stack(as.data.frame(t4))[,1],(stack(as.data.frame(t3))[,1]),use="p")^2,2)
	#par(mfrow=c(1,2))
	#plot(stack(as.data.frame(t2))[,1],stack(as.data.frame(t1))[,1],xlab="Proportion of hauls with species present (by 1 deg square)",ylab="Mean predicted probability of occurrence", cex=0.5,main=sp)
	#mtext(paste("r^2 =",presr2))
	#plot(stack(as.data.frame(t4))[,1],(stack(as.data.frame(t3))[,1]),xlab="Average log(WTCPUE) (by 1 deg square)",ylab="Average predicted log(WTCPUE)", cex=0.5,main=sp)
	#mtext(paste("r^2 =",abunr2))

	modeldiag$r2.pres.1deg[i]<-presr2
	modeldiag$r2.abun.1deg[i]<-abunr2

	#Also save average biomassmean for each region (across all years) to use in later predictions.
	avemeanbiomass<-apply(ave.catch.wt,2,mean,na.rm=T) #ave.catch.wt from far above

	####################################################
	#### Save models for later projections
	####################################################

	#Reduce gam object size by removing unnecessary parts (leaving what's needed to predict)
	#mygam1<-stripGAM(mygam1) #only reduces by a few percent - prob not useful
	#mygam2<-stripGAM(mygam2)
	runname <- "test"
	mods = list(mygam1=mygam1, mygam2 = mygam2)
	
	sp <- gsub('/', '', sp) # would mess up saving the file
	save(mods, avemeanbiomass, myregions, file=paste('../CEmodels/CEmods_',runname, '_', sp, '_', Sys.Date(), '.RData', sep='')) # ~2mb file

	#think about figures to output - thermal response curves? spatial prediction by 1 deg square?
	#think about other data to save - number of pres/abs by region (?) 
}

save(modeldiag,file=paste("output/modeldiag_",runname,"_",Sys.Date(),".Rdata",sep=""))
write.csv(modeldiag, file=paste("output/modeldiag_",runname,"_",Sys.Date(),".csv",sep=""))

write.csv(allwarnings, file=paste('output/warnings_', runname, '_', Sys.Date(), '.csv', sep=''))