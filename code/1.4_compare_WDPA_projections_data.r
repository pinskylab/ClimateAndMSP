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
trawl <- st_as_sf(trawl, crs = st_crs(wdpa), coords = c('lon', 'lat'))
# plot(trawl['region'], key.pos = 4, axes = TRUE, cex=0.02, key.width=lcm(5), asp=4)

# Prep for intersections
regs <- sort(unique(trawl$region))
if(exists('wdpa_trawlsppobs')) rm(wdpa_trawlsppobs) # remove output table if it exists

# Intersect trawl data with focal MPAs by region and only add zeros within regions
for (i in 1:length(regs)){
    print(paste(regs[i], Sys.time()))
    # Intersect trawl with WDPA and find species in each MPA
    thisint <- st_intersection(trawl[trawl$region == regs[i], ], wdpa) # returns points in trawl, since trawl is the first argument. only this region
    thiswdpaspp <- data.table(st_drop_geometry(thisint), stringsAsFactors = FALSE)[, .(noccur=.N), by=.(WDPA_PID, spp = spp.proj)] # summarize spp in each MPA
    
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

# Summarize across regions
wdpa_trawlsppobs <- wdpa_trawlsppobs_byreg[, .(noccur = sum(noccur), nhaul = sum(nhaul)), by=.(WDPA_PID, spp)]

# Save out
write.csv(wdpa_trawlsppobs, file=gzfile('output/wdpa_trawlsppobs.csv.gz'))
write.csv(wdpa_trawlsppobs_byreg, file=gzfile('temp/wdpa_trawlsppobs_byreg.csv.gz'))




########################################################
# 3. Compare observed and projected species in MPAs
########################################################
library(data.table)

# Read in data
wdpa_by_spp_obs <- fread('gunzip -c output/wdpa_trawlsppobs.csv.gz') # trawl survey observations of species in MPAs
wdpa_by_spp_projnow <- readRDS('temp/wdpa_by_spp_projnow.rds') # projections of current p(occur) from SDMs. From 1.2_evalWDPA.r

# Summarize projections by species (ensemble mean)
wdpa_by_spp_ensemble <- wdpa_by_spp_projnow[, .(poccur = mean(poccur)), by=.(WDPA_PID, spp)]

# Merge obs and proj
wdpa_by_spp_obs[, WDPA_PID := as.character(WDPA_PID)] # make sure wdpa id is a character
wdpasppobsproj <- merge(wdpa_by_spp_ensemble, wdpa_by_spp_obs[, .(spp, WDPA_PID, noccur, nhaul)], by = c('spp', 'WDPA_PID'))

# Basic stats
wdpasppobsproj[,hist(noccur, breaks = c(seq(-0.1, 10, by=1), 1000), xlim=c(0,11))] # most observations are of a species 1x in an MPA
wdpasppobsproj[, .(nspp=sum(noccur>0)), by=WDPA_PID][,hist(nspp, breaks = 20)] # most MPAs have 30-80 spp
wdpasppobsproj[, .(nhaul=max(nhaul)), by=WDPA_PID][,hist(nhaul, breaks = 200)] # most MPAs have <10 hauls


# Compare
# note that observed absences can be because
#   1. species isn't present
#   2. species is present but wasn't observed in the few years of sampling that we've analyzed
wdpasppobsproj[,plot(noccur ~ poccur)] # biplot
wdpasppobsproj[,.(occur = factor(noccur > 0), poccur=poccur)][,boxplot(poccur ~ occur)] # boxplot

    # densities
cols <- c('red', 'black')
par(mfrow=c(2,1), mai = c(0.8, 0.8, 0.2, 0.1), mgp=c(2.2, 0.5, 0), las=1, tcl= -0.2)
wdpasppobsproj[noccur == 0, plot(density(poccur, from = 0, to = 1), main='Not observed', xlab='Projected p(occur)', col = cols[2])]
wdpasppobsproj[noccur > 0, plot(density(poccur, from = 0, to = 1), main='Observed present', xlab='Projected p(occur)', col = cols[2])]

    # densities by nhauls
cols <- 'black'
par(mfrow=c(3,1), mai = c(0.8, 0.8, 0.2, 0.1), mgp=c(2.2, 0.5, 0), las=1, tcl= -0.2)
wdpasppobsproj[noccur == 0 & nhaul<5, plot(density(poccur, from = 0, to = 1), main='Not observed, nhauls<10', xlab='Projected p(occur)', col = cols)]
wdpasppobsproj[noccur == 0 & nhaul>=20, plot(density(poccur, from = 0, to = 1), main='Not observed, nhauls>=20', xlab='Projected p(occur)', col = cols)]
wdpasppobsproj[noccur > 0, plot(density(poccur, from = 0, to = 1), main='Observed present', xlab='Projected p(occur)', col = cols)]

    # confusion matrices
wdpasppobsproj[,table(obs = noccur>0, proj = poccur>0.5)] # confusion matrix with cutoff of p(occur)=0.5
wdpasppobsproj[,table(obs = noccur>0, proj = poccur>0.9)] # confusion matrix with cutoff of p(occur)=0.9

    # summaries
wdpasppobsproj[,print(summary(poccur)), by = (noccur > 0)] # summary of projections with observed pres vs. abs
wdpasppobsproj[,print(summary(poccur)), by = (noccur > 0)] # summary of projections with observed pres vs. abs

    # TSS and AUC metrics?
