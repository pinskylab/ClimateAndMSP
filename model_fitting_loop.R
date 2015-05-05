# Loop through species and fit models.
library(mgcv);library(ROCR)

load("data/dat_selectedspp.Rdata")

allspp = sort(unique(dat$sppocean))
n = rep(NA, length(allspp))
modeldiag = data.frame(sppocean = n, npres = n, ntot = n, auc = n, auc.tt=n, r2.biomass = n, r2.biomass.tt = n, r2.all = n, r2.all.tt = n, r2.pres.1deg = n, r2.abun.1deg = n, dev.pres = n, dev.biomass = n, stringsAsFactors=FALSE)


######################
# Start the big loop #
######################

for(i in 1:length(allspp)){ 

	sp<-allspp[i]
	print(paste(i,sp, Sys.time()))

	mydat<-dat[dat$sppocean==sp,] 
	myregions<-unique(mydat$region)
	myhauls<-unique(mydat$haulid)
	#myocean<-mydat[1,"ocean"]

	####################################################
	# Add records for when there was a haul but no fish
	####################################################

	# Only expand catches (pad with zeros) in regions where the species has previously been caught. 

	zs<-dat[!dat$haulid %in% myhauls & dat$region %in% myregions,] #extract haulids where missing for regions where species has at least one record.
	matchpos<-match(unique(zs$haulid),zs$haulid) # Extract one record of each haulid
	zeros<-zs[matchpos,]
	#Add/change relevant columns --> zero catch for target spp.
	zeros$spp<-mydat$spp[1]
	zeros$sppl<-mydat$sppl[1]
	zeros$sppnew<-mydat$sppnew[1]
	zeros$sppregion<-paste(zeros$sppnew,zeros$region,sep="_")
	zeros$sppocean<-mydat$sppocean[1] #may need to add "ocean"
	zeros$wtcpue<-0
	zeros$presfit<-F
	zeros$wtcpuena<-1e-4 #for consistency, but really zero
	zeros$wtcpuenal<-log(zeros$wtcpuena)

	mydatw0<-rbind(mydat,zeros) #combine positive hauls with zero hauls

	####################################################
	# Calculate mean catch per haul by region. 
	####################################################

	# (For true biomass estimate, would need to stratify)
	ave.catch.wt<-tapply(mydatw0$wtcpue,list(mydatw0$year,mydatw0$region),mean,na.rm=T)
	if(dim(ave.catch.wt)[2]<2) {
		avecatchyrreg<-cbind(as.data.frame(ave.catch.wt), rep(colnames(ave.catch.wt), dim(ave.catch.wt)[2]), rownames(ave.catch.wt))
		colnames(avecatchyrreg)<-c("biomassmean","region","year")
		}else{
		avecatchyrreg<-cbind(stack(as.data.frame(ave.catch.wt)), rep(rownames(ave.catch.wt), dim(ave.catch.wt)[2]))
		colnames(avecatchyrreg)<-c("biomassmean","region","year")
	}
	spdata<-merge(mydatw0,avecatchyrreg)

	####################################################
	# Trim data to complete cases
	####################################################

	#What if only a few records exist for catch in one of these regions - worth excuding surftemp frmo model? Or rather exclude region?
	if(any(c("DFO_NewfoundlandSpring","DFO_NewfoundlandFall","NWFSC_WCAnn") %in% myregions)){
		spdata<-spdata[complete.cases(spdata[,c("bottemp","rugosity","presfit")]),] 
		}else{
		spdata<-spdata[complete.cases(spdata[,c("surftemp","bottemp","rugosity","presfit")]),] 
	}

	# Subset by month or season?? (should test whether it matters)
	#mymonths<-c(10,11,12)
	#spdata<-spdata[spdata$month %in% mymonths,]

	####################################################
	# Figure out which model formula given data
	####################################################

	#Default models. Leave out region factor if necessary
	if(length(myregions)==1){
			mypresmod<-formula(presfit~s(bottemp,k=4)+s(surftemp,k=4)+s(logrugosity,k=4)+s(biomassmean,k=4))
			myabunmod<-formula(wtcpuenal~s(bottemp,k=4)+s(surftemp,k=4)+s(logrugosity,k=4)+s(biomassmean,k=4))
		} else {
			mypresmod<-formula(presfit~s(bottemp,k=4)+s(surftemp,k=4)+s(logrugosity,k=4)+regionfact+s(biomassmean,k=4)-1)
			myabunmod<-formula(wtcpuenal~s(bottemp,k=4)+s(surftemp,k=4)+s(logrugosity,k=4)+regionfact+s(biomassmean,k=4)-1)
	}
	#Replace with formula missing surftemp if necessary
	if(any(c("DFO_NewfoundlandSpring","DFO_NewfoundlandFall","NWFSC_WCAnn") %in% myregions)){
		if(length(myregions)==1){
			mypresmod<-formula(presfit~s(bottemp,k=4)+s(logrugosity,k=4)+s(biomassmean,k=4))
			myabunmod<-formula(wtcpuenal~s(bottemp,k=4)+s(logrugosity,k=4)+s(biomassmean,k=4))
		} else {
			mypresmod<-formula(presfit~s(bottemp,k=4)+s(logrugosity,k=4)+regionfact+s(biomassmean,k=4)-1)
			myabunmod<-formula(wtcpuenal~s(bottemp,k=4)+s(logrugosity,k=4)+regionfact+s(biomassmean,k=4)-1)
		}
	}


	####################################################
	#Fit models for training and testing to (optionally) evaluate performance
	####################################################

	#Subset training and testing data by year (use first 80% to predict last 20%)
	# ***Should change this to be first 80% from each region
	spdata<-spdata[order(spdata$year,spdata$month),]
	ninds<-dim(spdata)[1]
	traininds<-1:round(ninds*0.8)
	testinds<- (round(ninds*0.8)+1):ninds

	mygam1tt<-gam(mypresmod, family="binomial",data=spdata[traininds,]) 
	mygam2tt<-gam(myabunmod, data=spdata[traininds,]) 


	####################################################
	#Fit models to All data (no test/training split)
	####################################################

	mygam1<-gam(mypresmod,family="binomial",data=spdata) #no surftemp
	mygam2<-gam(myabunmod,data=spdata) #no surftemp

	####################################################
	# Compare predictions to observations to assess model performance
	####################################################

	# For FULL model
	preds1 <- predict(mygam1,spdata,type="response") #can also use mygam1$fitted.values
	preds2 <- exp(predict(mygam2, newdata = spdata, type='response')) 
	preds <- preds1*preds2 #***Need to add re-transformation bias correction***
	preds[preds<0] = 0

	# And for training/testing data set
	preds1tt <- predict(mygam1tt,spdata[testinds,],type="response") 
	preds2tt <- exp(predict(mygam2tt, newdata = spdata[testinds,], type='response')) 
	predstt <- preds1tt*preds2tt #***Need to add re-transformation bias correction***
	predstt[predstt<0] = 0

	# calculate performance (in part using ROCR)
	preds1.rocr = prediction(predictions=as.numeric(preds1), labels=spdata$presfit)
	preds1tt.rocr = prediction(predictions=as.numeric(preds1tt), labels=spdata$presfit[testinds])

	modeldiag$sppocean[i] = sp
	modeldiag$npres[i] = sum(spdata$presfit)
	modeldiag$ntot[i] = dim(spdata)[1]
	modeldiag$auc[i] = performance(preds1.rocr, 'auc')@y.values[[1]] #
	modeldiag$auc.tt[i] = performance(preds1tt.rocr, 'auc')@y.values[[1]] #

	# ***Double check code below! Not sure I'm think about re-transformed log predictions correctly
	modeldiag$r2.biomass[i] = cor(log(preds2[spdata$presfit]), spdata$wtcpuenal[spdata$presfit])^2 # correlation of log(biomass) where present
	#modeldiag$r2.biomass.tt[i] = cor(preds2tt[spdata$presfit[testinds]], spdata[spdata$presfit,"wtcpuenal"][testinds])^2 # ** code not working here - how to subset?
	modeldiag$r2.all[i] = cor(preds, spdata$wtcpue)^2 # overall biomass correlation
	modeldiag$r2.all.tt[i] = cor(predstt, spdata$wtcpue[testinds])^2 # overall biomass correlation
	modeldiag$dev.pres[i] = summary(mygam1)$dev.expl
	modeldiag$dev.biomass[i] = summary(mygam2)$dev.expl


	# Some metrics at a spatially aggregated level (by year) may be more informative:
	test<-cbind(spdata,preds1,preds)
	t1<-tapply(test$preds1,list(test$year,test$cs1),mean) #average predicted p(occur)
	t2<-tapply(test$presfit,list(test$year,test$cs1),mean) #proportion of hauls with presence
	t3<-tapply(log(test$preds),list(test$year,test$cs1),mean) #average predicted abundance ***check if log-transform make sense
	t4<-tapply(test$wtcpuenal,list(test$year,test$cs1),mean) #average observed abundance *** log abundance?

	presr2<-round(cor(stack(as.data.frame(t2))[,1],stack(as.data.frame(t1))[,1],use="p")^2,2)
	abunr2<-round(cor(stack(as.data.frame(t4))[,1],(stack(as.data.frame(t3))[,1]),use="p")^2,2)
	#par(mfrow=c(1,2))
	#plot(stack(as.data.frame(t2))[,1],stack(as.data.frame(t1))[,1],xlab="Proportion of hauls with species present (by 1 deg square)",ylab="Mean predicted probability of occurrence", cex=0.5,main=sp)
	#mtext(paste("r^2 =",presr2))
	#plot(stack(as.data.frame(t4))[,1],(stack(as.data.frame(t3))[,1]),xlab="Average log(WTCPUE) (by 1 deg square)",ylab="Average predicted log(WTCPUE)", cex=0.5,main=sp)
	#mtext(paste("r^2 =",abunr2))

	modeldiag$r2.pres.1deg[i]<-presr2
	modeldiag$r2.abun.1deg[i]<-abunr2

	#Also save average biomassmean for each region to use in later predictions.
	avemeanbiomass<-apply(ave.catch.wt,2,mean,na.rm=T) #ave.catch.wt from far above

	####################################################
	#### Save models for later projections
	####################################################

	#Reduce gam object size by removing unnecessary parts (leaving what's needed to predict)
	#mygam1<-stripGAM(mygam1) #only reduces by a few percent - prob not useful
	#mygam2<-stripGAM(mygam2)
	runname <- "test"
	mods = list(mygam1=mygam1, mygam2 = mygam2)
	save(mods, avemeanbiomass, myregions, file=paste('output/CEmods_',runname, '_', sp, '_', Sys.Date(), '.RData', sep='')) #This saves a lot of info about the model - large file 22 mb. Make smaller?

	#think about figures to output - thermal response curves? spatial prediction by 1 deg square?
	#think about other data to save - number of pres/abs by region (?) 
}

save(modeldiag,file=paste("output/modeldiag_",runname,"_",Sys.Date(),".Rdata",sep=""))
