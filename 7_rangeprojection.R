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
	climgridfolder <- '../data/'
	numcorestouse <- 2
	}
if(Sys.info()["nodename"] == "amphiprion.deenr.rutgers.edu"){
	setwd('~/Documents/range_projections/')
	projfolder = 'CEmodels_proj'
	modfolder = 'CEmodels'
	climgridfolder <- 'data/'
	numcorestouse <- 2
	}
# could add code for Lauren's working directory here

##############################
# Choose the model fit to use
##############################
#runtype <- 'test'
#runtype <- 'testseason'
runtype <- 'testK6noSeas'

#############################
# Choose species to project #
#############################
load(paste('output/modeldiag_', runtype, '.Rdata', sep='')) # model diagnostics

#projspp <- modeldiag$sppocean[modeldiag$auc.tt >= 0.75 & !is.na(modeldiag$auc.tt)] # from Elith et al., but there may be better criteria to use

#With additional criteria: ***
projspp <- modeldiag$sppocean[modeldiag$auc.tt >= 0.75 & !is.na(modeldiag$auc.tt) & ((modeldiag$dev.pres - modeldiag$dev.pres.null > 0.05) | (modeldiag$dev.biomass - modeldiag$dev.biomass.null > 0.05))] # from Elith et al., and deviance explained by temp must be > 5% for at least one model
length(projspp) # number of species to project to

# find the files with these species for our chosen model fit
files <- list.files(modfolder)
files <- files[grepl(paste('_', runtype, '_', sep=''), files) & grepl(paste(gsub('/|\\(|\\)', '', projspp), collapse='|'), gsub('/|\\(|\\)', '', files))] # have to strip out parentheses and slashes from file and taxon names so that grep doesn't interpret them
length(files) # should match length of projspp

## Remove spp from planned projections IF the projection file already exists (OPTIONAL). 
# If this step is skipped, the existing files will be overwritten.
donefiles <- list.files(projfolder, pattern=runtype) # models I made earlier
donespp <- gsub(paste('summproj_', runtype, '_', sep=''), '', gsub('.Rdata', '', donefiles))
if(length(donespp)>0){
	files <- files[!grepl(paste(gsub('/|\\(|\\)', '', donespp), collapse='|'), gsub('/|\\(|\\)', '', files))] # remove any models that we made earlier
	projspp <- projspp[!grepl(paste(gsub('/|\\(|\\)', '', donespp), collapse='|'), gsub('/|\\(|\\)', '', projspp))]
}
length(files)
length(projspp)


#################################
# Prep environmental data
#################################
load(paste(climgridfolder, 'climGrid.proj2.RData', sep='')) # projected temperature for each year ("clim")

# drop unneeded columns
clim <- clim[,!grepl('depthgrid', names(clim))] #  refer to GCM depth grids
clim <- clim[,!grepl('bottemp.clim|surftemp.clim|delta|latgrid|longrid', names(clim))] #  the temp climatologies, deltas, and GCM lat/lon grids (1 degree)

# fix season colum
clim$season <- as.factor(c('wi', 'sp', 'su', 'fa')[clim$season]) # convert to same format we used for model fitting

# change region column to combine NEUSSpring and NEUSFall
clim$region[clim$region == 'NEFSC_NEUSSpring'] = 'NEFSC_NEUS'
clim$region[clim$region == 'NEFSC_NEUSFall'] = 'NEFSC_NEUS'

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
doprojection <- function(thisprojspp, files, clim, projfolder, modfolder, runtype){ 
	mod <- avemeanbiomass <- NULL
	fileindex <- which(grepl(gsub('/|\\(|\\)', '', thisprojspp), gsub('/|\\(|\\)', '', files)))
	print(paste(fileindex, thisprojspp, Sys.time()))

	load(paste(modfolder, '/', files[fileindex], sep='')) # loads mod and avemeanbiomass

	# modify the GAMs to remove the regionfact term?
	# or set all regions within the ocean of this taxon to a region name in the model? (to 'trick' the model)
	# this would allow us to project outside the regions for which we had data
	# for now, we don't do this

	# add mean biomass by region
	clim$biomassmean <- 0
	clim$biomassmean[clim$regionfact %in% names(avemeanbiomass)] <- avemeanbiomass[as.character(clim$regionfact[clim$regionfact %in% names(avemeanbiomass)])] # use region to pull the correct mean biomass values

	# smearing estimator for re-transformation bias (see Duan 1983, http://www.herc.research.va.gov/resources/faq_e02.asp)
	smear <- mean(exp(mods[['mygam2']]$residuals))

	# get seasons from this model (all seasons since non-seasonal model)
	myseasons <- c('wi', 'sp', 'su', 'fa')

	# rows of clim to project to
	inds <- clim$regionfact %in% names(avemeanbiomass) & clim$season %in% myseasons # for models with season or without

	# Dataframe for this species' projections
	thisproj <- clim[inds,c('regionfact', 'lat', 'lon', 'depth16th', 'year', 'season')] # dataframe to hold projections for this taxon
	for(i in 1:13) thisproj[[paste('wtcpue.proj_', i, sep='')]] <- NA 	# Add projected biomass density columns for each climate model


	# Calculate predictions for 2020-2100 for each model
	for(i in 1:13){ # this alone takes a long time
		print(paste(fileindex, 'model', i))
		nd <- data.frame(regionfact = clim$regionfact[inds], surftemp = clim[[paste('surftemp.proj_', i, sep='')]][inds], bottemp = clim[[paste('bottemp.proj_', i, sep='')]][inds], logrugosity = log(clim$rugosity[inds]+0.01), biomassmean = clim$biomassmean[inds], season = clim$season[inds], row.names=1:sum(inds))
		preds1 <- predict.gam(mods$mygam1, newdata = nd, type='response')
		preds2 <- exp(predict(mods$mygam2, newdata = nd, type='response'))
		preds <- preds1*preds2*smear
		preds[preds<0] <- 0
		
		pnm = paste('wtcpue.proj_', i, sep='')
		thisproj[[pnm]] = preds # can do this because clim[inds,] and thisproj are in the same order

		# set biomass projections on land to zero
		thisproj[[pnm]][thisproj$depth16th > 0] <- 0
	}
	
	# summarize by 1/4 degree grid
	summproj <- aggregate(thisproj[,grepl('wtcpue.proj', names(thisproj))], by=list(region=thisproj$region, lat=thisproj$lat, lon=thisproj$lon, year=thisproj$year, season=thisproj$season), FUN=mean)
	
#	print(summary(summproj))
#	print(dim(summproj))	

	thisprojspp <- gsub('/', '', thisprojspp) # would mess up saving the file
	save(summproj, file=paste(projfolder, '/summproj_', runtype, '_', thisprojspp, '.Rdata', sep='')) # write out the projections (15MB file)
}

result <- mclapply(X= projspp, FUN=doprojection, files=files, clim=clim, projfolder=projfolder, modfolder=modfolder, runtype=runtype, mc.cores=numcorestouse) # spawn out to multiple cores. errors will be stored in results

print(result)