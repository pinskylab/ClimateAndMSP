# Read in temperature fields and models, then make range projections
# This could probably be sped up by switching from data.frames to data.tables

require(mgcv)
require(Hmisc)
require(parallel) # for multi-core calculations


## Set working directory
if(Sys.info()["nodename"] == "pinsky-macbookair"){
	setwd('~/Documents/Rutgers/Range projections/proj_ranges/')
	projfolder = '../CEmodels_proj'
	modfolder = '../CEModels'
	numcorestouse <- 2
	}
if(Sys.info()["nodename"] == "amphiprion.deenr.rutgers.edu"){
	setwd('~/Documents/range_projections/')
	projfolder = 'CEmodels_proj'
	modfolder = 'CEmodels'
	numcorestouse <- 20
	}
# could add code for Lauren's working directory here

##############################
# Choose the model fit to use
##############################
#runtype <- 'test'
runtype <- 'testseason'

#############################
# Choose species to project #
#############################
load(paste('output/modeldiag_', runtype, '.Rdata', sep='')) # model diagnostics

projspp <- modeldiag$sppocean[modeldiag$auc.tt >= 0.75 & !is.na(modeldiag$auc.tt)] # from Elith et al., but there may be better criteria to use
length(projspp) # number of species to project to

# find the files with these species for our chosen model fit
files <- list.files(modfolder)
files <- files[grepl(paste('_', runtype, '_', sep=''), files) & grepl(paste(gsub('/|\\(|\\)', '', projspp), collapse='|'), gsub('/|\\(|\\)', '', files))] # have to strip out parentheses and slashes from file and taxon names so that grep doesn't interpret them
length(files) # should match length of projspp

#################################
# Prep environmental data
#################################
load('data/climGrid.proj2_2015-02-10.RData') # projected temperature for each year ("clim")

# drop unneeded columns
clim <- clim[,!grepl('depthgrid', names(clim))] #  refer to GCM depth grids
clim <- clim[,!grepl('bottemp.clim|surftemp.clim|delta|latgrid|longrid', names(clim))] #  the temp climatologies, deltas, and GCM lat/lon grids (1 degree)

# fix season colum
clim$season <- as.factor(c('wi', 'sp', 'su', 'fa')[clim$season]) # convert to same format we used for model fitting

# add regionfact
clim$region<- as.factor(clim$region)
names(clim)[names(clim)=='region'] <- 'regionfact'

# add logrugosity
rugos <- read.csv('data/projectiongrid_latlons.1.16th_withRugosity_2015-05-06.csv')
	names(rugos)[names(rugos) == 'lon'] <- 'lon16th'
	names(rugos)[names(rugos) == 'lat'] <- 'lat16th'
	names(rugos)[names(rugos) == 'depth'] <- 'depth16th'

	gridsize=0.25 # size of grid of the climate data, in degrees
	rugos$lat <- floor(rugos$lat16th/gridsize)*gridsize + gridsize/2 # round to nearest grid center
	rugos$lon <- floor(rugos$lon16th/gridsize)*gridsize + gridsize/2

clim <- merge(clim, rugos) # slow
dim(clim) # 11,050,240 rows


############################################
## Project GAMS onto annual climate data  ##
############################################
	
options(warn=1) # print warnings as they occur

# VERY slow
doprojection <- function(projspp, files, clim, projfolder, modfolder, runtype){ 
	print(paste(projspp, Sys.time()))

	mod <- avemeanbiomass <- NULL

	infile <- grep(projspp, files, value=TRUE)
	load(paste(modfolder, '/', infile, sep='')) # loads mod and avemeanbiomass

	# modify the GAMs to remove the regionfact term?
	# or set all regions within the ocean of this taxon to a region name in the model? (to 'trick' the model)
	# this would allow us to project outside the regions for which we had data
	# for now, we don't do this

	# add mean biomass by region
	clim$biomassmean <- 0
	clim$biomassmean[clim$regionfact %in% names(avemeanbiomass)] <- avemeanbiomass[as.character(clim$regionfact[clim$regionfact %in% names(avemeanbiomass)])] # use region to pull the correct mean biomass values

	# smearing estimator for re-transformation bias (see Duan 1983, http://www.herc.research.va.gov/resources/faq_e02.asp)
	smear <- mean(exp(mods[['mygam2']]$residuals))

	# get seasons from this model
	myseasons <- gsub('s(bottemp):season', '', grep('bottemp', names(mods$mygam1$sp), value=TRUE), fixed=TRUE)
	if(myseasons == 's(bottemp)') myseasons <- c('wi', 'sp', 'su', 'fa') # catch if this is a non-seasonal model

	# rows of clim to project to
#	inds <- clim$regionfact %in% names(avemeanbiomass) # for non-seasonal models
	inds <- clim$regionfact %in% names(avemeanbiomass) & clim$season %in% myseasons # for models with season

	# Dataframe for this species' projections
	thisproj <- clim[inds,c('regionfact', 'lat', 'lon', 'depth16th', 'year', 'season')] # dataframe to hold projections for this taxon
	for(i in 1:13) thisproj[[paste('wtcpue.proj_', i, sep='')]] <- NA 	# Add projected biomass density columns for each climate model


	# Calculate predictions for 2020-2100 for each model
	for(i in 1:13){ # this alone takes a long time
		print(paste('model', i))
		nd <- data.frame(regionfact = clim$regionfact[inds], surftemp = clim[[paste('surftemp.proj_', i, sep='')]][inds], bottemp = clim[[paste('bottemp.proj_', i, sep='')]][inds], logrugosity = log(clim$rugosity[inds]+0.01), biomassmean = clim$biomassmean[inds], season = clim$season[inds], row.names=1:sum(inds))
		preds1 <- predict.gam(mods$mygam1, newdata = nd, type='response')
		preds2 <- exp(predict(mods$mygam2, newdata = nd, type='response'))
		preds <- preds1*preds2*smear
		preds[preds<0] <- 0
		
		pnm = paste('wtcpue.proj_', i, sep='')
		thisproj[[pnm]] = preds

		# set biomass projections on land to zero
		thisproj[[pnm]][thisproj$depth16th > 0] <- 0
	}
	
	# summarize by 1/4 degree grid
	summproj <- aggregate(thisproj[,grepl('wtcpue.proj', names(thisproj))], by=list(region=thisproj$region, lat=thisproj$lat, lon=thisproj$lon, year=thisproj$year, season=thisproj$season), FUN=mean)
	
#	print(summary(summproj))
#	print(dim(summproj))	
	
	save(summproj, file=paste(projfolder, '/summproj_', runtype, '_', projspp, '.Rdata', sep='')) # write out the projections (15MB file)
}

result <- mclapply(X= projspp, FUN=doprojection, files=files, clim=clim, projfolder=projfolder, modfolder=modfolder, runtype=runtype, mc.cores=numcorestouse) # spawn out to multiple cores. errors will be stored in results