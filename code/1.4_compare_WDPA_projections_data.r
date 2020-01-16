# Compare species p(occur) projections to species observations within MPAs


########################################################
# 1. Harmonize observation species names with predictions
# ONLY DO THIS ONCE SINCE TIME-CONSUMING
# THEN SKIP to 2 AND READ IN SAVED FILE INSTEAD
########################################################
library(data.table)

# trawl data from >=2016, as curated by OceanAdapt
trawl <- data.table(readRDS('dataDL/oceanadapt/all-regions-full.rds'))
trawl <- trawl[year >= 2016,] # trim to years we want
trawl <- trawl[!duplicated(spp), .(spp)] # trim to only species names

# Read in data
wdpa_by_spp_projnow <- readRDS('temp/wdpa_by_spp_projnow.rds') # from 1.2_evalWDPA.r

# Harmonize observation species names
trawl[, spp.proj := as.character(rep(NA, .N))] # create a column to hold the new names that match the projection names
obsspp <- trawl[, tolower(sort(unique(spp)))] # the list of observed species names, in lower case
projspp <- wdpa_by_spp_projnow[, sort(unique(spp))] # the list of species names used for the projections
for (i in 1:length(obsspp)) { # step through each observed name
    if(obsspp[i] %in% projspp){
        trawl[tolower(spp) == obsspp[i], spp.proj := obsspp[i]] # if user chose one of the species names
    } else {
        inds <- agrep(obsspp[i], projspp, max.distance = 0.1) # approximate matching
        if(length(inds) > 0){
            cat(paste0('\nMatch ', obsspp[i], ' to:\n1 none of the above\n', paste(paste(2:(length(inds)+1), projspp[inds]), collapse='\n')))
            continue <- FALSE
            while(continue == FALSE) {
                n <- suppressWarnings(as.numeric(readline(prompt = 'Enter an integer: ')))
                if(!is.na(n)){
                    if(n <= (length(inds)+1)){
                        continue <- TRUE # acceptable answer. keep going
                    }
                }
            }
            n <- as.numeric(n)
            if (n == 1) {
                # if user chose 'none of the above', leave as NA
            }
            if (n > 1 & n <= (length(inds)+1)) {
                trawl[tolower(spp) == obsspp[i], spp.proj := projspp[inds][n-1]] # if user chose one of the species names, copy it in
            }
        } else {
            print(paste("No match for", obsspp[i]))
            # if no match, leave as NA
        }
        
    }
    
}

nameconv <- trawl[!is.na(spp.proj), ][!duplicated(spp), .(spp, spp.proj)] # all of the name conversions

print(nameconv, n = 600) # all of the name conversions, including tolower
nameconv[tolower(spp) != spp.proj, ] # where different (more than just tolower)

# save the species name conversions
write.csv(nameconv, file='output/name_conversions_obs_to_proj.csv', row.names = FALSE)


###################################################
# 2. Find species in each focal MPA
# ONLY DO THIS ONCE SINCE TIME-CONSUMING
# THEN SKIP to 3 AND READ IN SAVED FILE INSTEAD
###################################################
library(sf) # for GIS functions
library(data.table)

# read in list of focal MPAs (those that overlap the trawl data)
wdpa_by_grid <- readRDS('temp/wdpa_by_grid0.05_intersect.rds') # from 1.0_processWDPA.r
wdpalist <- sort(unique(as.character(wdpa_by_grid$WDPA_PID))) # pull out the WPDA ids we want
rm(wdpa_by_grid) # clean up

# read in MPA shapefile
wdpa <- st_read('dataDL/WDPA/WDPA_Aug2019_marine-shapefile/WDPA_Aug2019_marine-shapefile-polygons.shp')
wdpa <- wdpa[wdpa$WDPA_PID %in% wdpalist,] # Trim MPAs to focal MPAs

# read in name conversions from observed to projected
nameconv <- fread('output/name_conversions_obs_to_proj.csv')

# read in trawl data from >=2016, as curated by OceanAdapt
trawl <- readRDS('dataDL/oceanadapt/all-regions-full.rds')
trawl <- trawl[trawl$year >= 2016,] # trim to years we want
trawl <- merge(trawl, nameconv, by = 'spp') # merge in name conversions. also trims to just those species in the projections.
trawl <- st_as_sf(trawl, crs = st_crs(wdpa), coords = c('lon', 'lat')) # create an sf GIS object
range(trawl$year)
# plot(trawl['region'], key.pos = 4, axes = TRUE, cex=0.02, key.width=lcm(5), asp=4)

# Prep for intersections
regs <- sort(unique(trawl$region))
if(exists('wdpa_trawlsppobs_byreg')) rm(wdpa_trawlsppobs_byreg) # remove output table if it exists

# Intersect trawl data with focal MPAs by region and only add zeros within regions
for (i in 1:length(regs)){
    print(paste(regs[i], Sys.time()))
    # Intersect trawl with WDPA and find species in each MPA
    thisint <- st_intersection(trawl[trawl$region == regs[i], ], wdpa) # returns points in trawl, since trawl is the first argument. only this region. the slow step.
    thiswdpaspp <- data.table(st_drop_geometry(thisint), stringsAsFactors = FALSE)[, .(noccur=length(unique(haulid))), by=.(WDPA_PID, spp = spp.proj)] # summarize spp in each MPA
    
    # Add 0s for species not observed
    thesespp <- sort(unique(trawl$spp.proj[trawl$region == regs[i]])) # all observed species in the region
    thesewdpas <- sort(unique(thisint$WDPA_PID)) # all MPAs that overlap trawl
    thisfullwdpaspp <- CJ(spp = thesespp, WDPA_PID = thesewdpas) # make a table of all spp x WDPA combinations
    thisfullwdpaspp <- merge(thisfullwdpaspp, thiswdpaspp, all.x = TRUE, by = c('spp', 'WDPA_PID')) # merge in spp observations
    thisfullwdpaspp[is.na(noccur), noccur := 0] # set lack of obs to noccur=0

    # add number of hauls per MPA
    thiswdpahauls <- data.table(st_drop_geometry(thisint), stringsAsFactors = FALSE)[, .(nhaul=length(unique(haulid))), by=.(WDPA_PID)] # summarize number of hauls per MPA (effort)
    nrow1 <- nrow(thisfullwdpaspp)
    thisfullwdpaspp <- merge(thisfullwdpaspp, thiswdpahauls, by='WDPA_PID')
    if(nrow1 != nrow(thisfullwdpaspp)) stop('Lost rows when merging thiswdpahauls') # a simple error check
    
    # add to output table across all regions
    thisfullwdpaspp[, region := regs[i]]
    if (!exists('wdpa_trawlsppobs_byreg')) {
        wdpa_trawlsppobs_byreg <- thisfullwdpaspp # create table if it doesn't exist
    } else {
        wdpa_trawlsppobs_byreg <- rbind(wdpa_trawlsppobs_byreg, thisfullwdpaspp) # append to end if table already exists
    }
}
Sys.time()

# Summarize across surveys
wdpa_trawlsppobs <- wdpa_trawlsppobs_byreg[, .(noccur = sum(noccur), nhaul = sum(nhaul)), by=.(WDPA_PID, spp)]

# Save out
write.csv(wdpa_trawlsppobs_byreg, file=gzfile('temp/wdpa_trawlsppobs_byreg.csv.gz'))
write.csv(wdpa_trawlsppobs, file=gzfile('output/wdpa_trawlsppobs.csv.gz'))


# wdpa_trawlsppobs_byreg <- fread('gunzip -c temp/wdpa_trawlsppobs_byreg.csv.gz', drop = 1)

########################################################
# 3. Compare observed and projected species in MPAs
########################################################
library(data.table)
library(ggplot2)
library(beanplot)

# Read in data
wdpa_by_spp_obs_reg <- fread('gunzip -c temp/wdpa_trawlsppobs_byreg.csv.gz', drop = 1) # trawl survey observations of species in MPAs
wdpa_by_spp_projnow <- readRDS('temp/wdpa_by_spp_projnow.rds') # projections of current p(occur) from SDMs. From 1.2_evalWDPA.r
poccurthresh <- fread('https://raw.githubusercontent.com/pinskylab/project_velocity/master/output/modeldiag_Nov2017_fitallreg_2017.csv', drop = 1)[, .(sppocean, thresh.kappa)] # poccur threshold: how high does the probability of occurrence in the projections need to be to consider the species "present"? use the thresholds calculated during model fitting from Morley et al. 2018 PLOS ONE

# Add poccur threshold to wdpa_by_spp_obs_reg (the trawl observations)
# rather than to projections, because the trawl observations still have region data (from survey name)
# trawl spp names have already been converted to the projection names
# except one issue: aspidophoroides bartoni is a synonym in the Pacific for Aspidophoroides monopterygius. Latter is accepted, but trawl uses monopterygius, projections use bartoni.
# note that we lack thresh.kappa for some spp in trawl: those that appear in two oceans and for which we have projections in only one ocean (lebbeus groenlandicus, crangon septemspinosa, etc.)
poccurthresh[, ocean := gsub('.*_', '', sppocean)]
poccurthresh[, spp := gsub('_Atl|_Pac', '', sppocean)]
poccurthresh[spp == 'aspidophoroides bartoni', spp := 'aspidophoroides monopterygius'] # fix proj name so that it matches trawl in the right ocean
wdpa_by_spp_obs_regPac <- merge(wdpa_by_spp_obs_reg[region %in% c("Aleutian Islands", "Eastern Bering Sea", "Gulf of Alaska", "West Coast Annual"), ], 
                                poccurthresh[ocean == 'Pac', .(spp, thresh.kappa)], by = 'spp', all.x = TRUE) # have to do Atl and Pac separately since some species are in both regions but use different models
wdpa_by_spp_obs_regAtl <- merge(wdpa_by_spp_obs_reg[region %in% c("Gulf of Mexico", "Northeast US Fall", "Northeast US Spring", "Scotian Shelf", 
                                                                  "Southeast US Fall", "Southeast US Spring", "Southeast US Summer"), ], 
                                poccurthresh[ocean == 'Atl', .(spp, thresh.kappa)], by = 'spp', all.x = TRUE)
if(nrow(wdpa_by_spp_obs_reg) == nrow(wdpa_by_spp_obs_regAtl) + nrow(wdpa_by_spp_obs_regPac)){
    wdpa_by_spp_obs_reg <- rbind(wdpa_by_spp_obs_regPac, wdpa_by_spp_obs_regAtl)
    rm(wdpa_by_spp_obs_regPac, wdpa_by_spp_obs_regAtl)
} else {
    stop('merge of poccurthesh and presmap did not work')
}
wdpa_by_spp_obs_reg[is.na(thresh.kappa), table(spp, region)] # spp lacking kappa thresholds, by region

# Sumarrize trawl data by MPA (across surveys)
wdpa_by_spp_obs <- wdpa_by_spp_obs_reg[, .(noccur = sum(noccur), nhaul = sum(nhaul), nkappa = length(unique(thresh.kappa)), 
                                           thresh.kappa = unique(thresh.kappa)), by=.(WDPA_PID, spp)]
if(any(wdpa_by_spp_obs$nkappa > 1)) stop('too many kappa values for a species')

# Summarize projections by species (ensemble mean)
wdpa_by_spp_ensemble <- wdpa_by_spp_projnow[, .(poccur = mean(poccur)), by=.(WDPA_PID, spp)]

# Merge obs and proj
wdpa_by_spp_obs[, WDPA_PID := as.character(WDPA_PID)] # make sure wdpa id is a character
wdpasppobsproj <- merge(wdpa_by_spp_ensemble, wdpa_by_spp_obs[, .(spp, WDPA_PID, noccur, nhaul, thresh.kappa)], by = c('spp', 'WDPA_PID'))


# Basic stats
wdpasppobsproj[,hist(nhaul, breaks = seq(0,1000,by=10))] # most MPAs have < 10 hauls
wdpasppobsproj[,hist(nhaul, breaks = seq(0,1000,by=1), xlim = c(0, 10))] # many MPAs have < 5 hauls
wdpasppobsproj[,hist(noccur, breaks = c(seq(-0.1, 10, by=1), 1000), xlim=c(0,11))] # most observations are of a species 1x in an MPA
wdpasppobsproj[, .(nspp=sum(noccur>0)), by=WDPA_PID][,hist(nspp, breaks = 20)] # most MPAs have 30-80 spp
wdpasppobsproj[, .(nhaul=max(nhaul)), by=WDPA_PID][,hist(nhaul, breaks = 200)] # most MPAs have <10 hauls

# Compare
# note that observed absences can be because
#   1. species isn't present
#   2. species is present but wasn't observed in the few years of sampling that we've analyzed
wdpasppobsproj[,plot(noccur ~ poccur)] # biplot of # occurrences vs. projected occurrence probability
wdpasppobsproj[,plot(I(noccur/nhaul) ~ poccur)] # biplot of proportion occurrences vs. projected occurrence probability
wdpasppobsproj[,.(occur = factor(noccur > 0), poccur=poccur)][, boxplot(poccur ~ occur)] # boxplot

    # confusion matrices
conf <- wdpasppobsproj[,table(obs = noccur>0, proj = poccur>thresh.kappa)] # confusion matrix with cutoff from kappa
conf[2,2]/sum(conf[,2]) # positive predictive value
conf[1,1]/sum(conf[,1]) # negative predictive value
ao = (conf[2,2] + conf[1,1])/sum(conf); ao # overall accuracy
ae = (sum(conf[2,])*sum(conf[,2]) + sum(conf[1,])*sum(conf[,1]))/sum(conf)^2; (ao-ae)/(1-ae) # kappa
Se = conf[2,2]/sum(conf[2,]); Se # sensitivity
Sp = conf[1,1]/sum(conf[1,]); Sp # specificity
Se + Sp - 1 # true skill statistic


    # PPV and NPV vs. nhauls in an MPA
ppv <- function(occur, predoccur){
    conf <- table(occur, predoccur) # confusion matrix
    if(nrow(conf) == 2 & ncol(conf) == 2) return(conf[2,2]/sum(conf[,2])) # positive predictive value
    else return(NA_real_)
}
npv <- function(occur, predoccur){
    conf <- table(occur, predoccur) # confusion matrix
    if(nrow(conf) == 2 & ncol(conf) == 2) return(conf[1,1]/sum(conf[,1])) # positive predictive value
    else return(NA_real_)
}
tss <- function(occur, predoccur){
    conf <- table(occur, predoccur) # confusion matrix
    if(nrow(conf) == 2 & ncol(conf) == 2){
        Se = conf[2,2]/sum(conf[2,]); Se # sensitivity
        Sp = conf[1,1]/sum(conf[1,]); Sp # specificity
        return(Se + Sp - 1) # true skill statistic
    }
    else{
        return(NA_real_)
    } 
}

thresh <- data.table(min = c(1, 2, 3, 5, 10, 100), npv = NA_real_, ppv = NA_real_, npvse = NA_real_, ppvse = NA_real_)
for(i in 1:nrow(thresh)){
    thresh$npv[i] <- wdpasppobsproj[nhaul >= thresh$min[i], npv(noccur>0, poccur>thresh.kappa)]
    thresh$ppv[i] <- wdpasppobsproj[nhaul >= thresh$min[i], ppv(noccur>0, poccur>thresh.kappa)]

    # bootstrap over the values
    npvs <- rep(NA, 100)
    ppvs <- rep(NA, 100)
    tsss <- rep(NA, 100)
    for(j in 1:100){
        npvs[j] <- wdpasppobsproj[nhaul >= thresh$min[i], ][sample(.N, replace = TRUE), ][, npv(noccur>0, poccur>thresh.kappa)]
        ppvs[j] <- wdpasppobsproj[nhaul >= thresh$min[i], ][sample(.N, replace = TRUE), ][, ppv(noccur>0, poccur>thresh.kappa)]
    }
    thresh$npvse[i] <- sqrt(1/99*sum((npvs - thresh$npv[i])^2)) # bootstrap SE
    thresh$ppvse[i] <- sqrt(1/99*sum((ppvs - thresh$ppv[i])^2))
}
plot(1000, 1000, xlim = c(1, 100), ylim = c(0, 1), log = 'x', xlab = 'Minimum # hauls', ylab = 'Predictive value')
thresh[, polygon(c(min, rev(min)), c(npv+npvse, rev(npv-npvse)), col = 'orange', border = NA)]
thresh[, polygon(c(min, rev(min)), c(ppv+ppvse, rev(ppv-ppvse)), col = 'grey', border = NA)]
thresh[, lines(min, npv, type = 'o', col = 'red')]
thresh[, lines(min, ppv, type ='o')]

    # pearson and spearman correlation across all species at once
wdpasppobsproj[, cor.test(noccur, poccur)]
wdpasppobsproj[, cor.test(noccur, poccur, method = 'spearman')]

# write out NPV and PPV
write.csv(thresh, 'output/MPAvstrawl_NPV_PPV.csv')

