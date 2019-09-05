# Compare species p(occur) projections to species observations within MPAs


###################################################
# 1. Find species in each focal MPA
# ONLY DO THIS ONCE SINCE TIME-CONSUMING
# THEN SKIP to 2. AND READ IN SAVED FILE INSTEAD
###################################################
library(sf) # for GIS functions
library(data.table)

# read in list of focal MPAs
wdpa_by_grid <- readRDS('temp/wdpa_by_grid0.05_intersect.rds') # from 1.0_processWDPA.r
wdpalist <- sort(unique(as.character(wdpa_by_grid$WDPA_PID))) # pull out the WPDA ids we want
rm(wdpa_by_grid) # clean up

# read in MPA shapefile
wdpa <- st_read('dataDL/WDPA/WDPA_Aug2019_marine-shapefile/WDPA_Aug2019_marine-shapefile-polygons.shp')
wdpa <- wdpa[wdpa$WDPA_PID %in% wdpalist,] # Trim MPAs to focal MPAs

# trawl data from >=2016, as curated by OceanAdapt
trawl <- readRDS('dataDL/oceanadapt/all-regions-full.rds')
trawl <- trawl[trawl$year >= 2016,]
trawl <- st_as_sf(trawl, crs = st_crs(wdpa), coords = c('lon', 'lat'))
# plot(trawl['region'], key.pos = 4, axes = TRUE, cex=0.02, key.width=lcm(5), asp=4)

# Intersect trawl data with focal MPAs
trawl_by_wdpa <- st_intersection(trawl, wdpa) # 8 hours. returns points in trawl

# Save out
saveRDS(trawl_by_wdpa, file='temp/trawl_by_wdpa.rds')


#############################
# 2. Make observed species list by MPA
# ONLY DO THIS ONCE
# THEN SKIP to 3. AND READ IN SAVED FILE INSTEAD
#############################
library(sf) # for GIS functions
library(data.table)

# read in intersection of MPAs with trawl tow data
trawl_by_wdpa <- readRDS('temp/trawl_by_wdpa.rds')

# summarize by MPA
wdpaspp <- data.table(st_drop_geometry(trawl_by_wdpa), stringsAsFactors = FALSE)[, .(noccur=.N), by=.(WDPA_PID, NAME, SUB_LOC, spp)]

    wdpaspp[,hist(noccur, breaks = c(seq(-0.1, 10, by=1), 1000), xlim=c(0,11))] # most observations are of a species 1x in an MPA
    wdpaspp[, .(nspp=.N), by=WDPA_PID][,hist(nspp, breaks = 20)] # most MPAs have 50-100 spp

# write out
write.csv(wdpaspp, file = gzfile('output/wdpa_trawlsppobs.csv.gz'), row.names = FALSE)


########################################################
# 3. Harmonize observation species names with predictions
# ONLY DO THIS ONCE SINCE TIME-CONSUMING
# THEN SKIP to 4. AND READ IN SAVED FILE INSTEAD
########################################################
library(data.table)

# Read in data
wdpa_by_spp_obs <- fread('gunzip -c output/wdpa_trawlsppobs.csv.gz')
wdpa_by_spp_projnow <- readRDS('temp/wdpa_by_spp_projnow.rds')

# Harmonize observation species names
wdpa_by_spp_obs[, spp.proj := as.character(rep(NA, .N))] # create a column to hold the new names
obsspp <- wdpa_by_spp_obs[, tolower(sort(unique(spp)))] # the list of observed species names, in lower case
projspp <- wdpa_by_spp_projnow[, sort(unique(spp))] # the list of species names used for the projections
for (i in 1:length(obsspp)) { # step through each observed name
    if(obsspp[i] %in% projspp){
        wdpa_by_spp_obs[tolower(spp) == obsspp[i], spp.proj := obsspp[i]] # if user chose one of the species names
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
                wdpa_by_spp_obs[tolower(spp) == obsspp[i], spp.proj := projspp[inds][n-1]] # if user chose one of the species names, copy it in
            }
        } else {
            print(paste("No match for", obsspp[i]))
            # if no match, leave as NA
        }
        
    }
    
}

nameconv <- wdpa_by_spp_obs[!is.na(spp.proj), ][!duplicated(spp), .(spp, spp.proj)] # all of the name conversions

print(nameconv, n = 600) # all of the name conversions, including tolower
nameconv[tolower(spp) != spp.proj, ] # where different (more than just tolower)

# save the species name conversions
write.csv(nameconv, file='output/name_conversions_obs_to_proj.csv', row.names = FALSE)


########################################################
# 4. Compare observed and projected species in MPAs
########################################################
library(data.table)

# Read in data
wdpa_by_spp_obs <- fread('gunzip -c output/wdpa_trawlsppobs.csv.gz')
wdpa_by_spp_projnow <- readRDS('temp/wdpa_by_spp_projnow.rds')
nameconv <- fread('output/name_conversions_obs_to_proj.csv')

# Add name conversions to the observations
# Also trims to just those species that are in the projections
wdpa_by_spp_obs <- merge(wdpa_by_spp_obs, nameconv)

# Summarize projections by species (ensemble mean)
wdpa_by_spp_ensemble <- wdpa_by_spp_projnow[, .(poccur = mean(poccur)), by=.(WDPA_PID, spp)]

# Merge obs and proj
wdpa_by_spp_obs[, WDPA_PID := as.character(WDPA_PID)] # make sure wdpa id is a character
wdpasppobsproj <- merge(wdpa_by_spp_ensemble, wdpa_by_spp_obs[, .(spp.proj, WDPA_PID, noccur)], by.x = c('spp', 'WDPA_PID'), by.y = c('spp.proj', 'WDPA_PID'), all.x = TRUE) # keep all the projections

# Set lack of observations to 0
wdpasppobsproj[is.na(noccur), noccur := 0]

# Trim out the MPAs with no observations
wdpasppobsproj <- wdpasppobsproj[WDPA_PID %in% wdpa_by_spp_obs$WDPA_PID, ]

# Compare
# note that observed absences can be because
#   1. species isn't present
#   2. species is present but wasn't observed in the few years of sampling that we've analyzed
#   3. species isn't recorded by the survey (e.g., inverts in some surveys)
wdpasppobsproj[,plot(noccur ~ poccur)] # biplot
wdpasppobsproj[,boxplot(poccur ~ I(noccur>0))] # boxplot

    # densities
cols <- c('red', 'black')
par(mfrow=c(2,1), mai = c(0.8, 0.8, 0.2, 0.1), mgp=c(2.2, 0.5, 0), las=1, tcl= -0.2)
wdpasppobsproj[noccur == 0, plot(density(poccur), main='Not observed', xlab='Projected p(occur)', col = cols[2])]
wdpasppobsproj[noccur > 0, plot(density(poccur), main='Observed present', xlab='Projected p(occur)', col = cols[2])]

wdpasppobsproj[,table(obs = noccur>0, proj = poccur>0.5)] # confusion matrix with cutoff of p(occur)=0.5
wdpasppobsproj[,table(obs = noccur>0, proj = poccur>0.9)] # confusion matrix with cutoff of p(occur)=0.9
wdpasppobsproj[,print(summary(poccur)), by = (noccur > 0)] # summary of projections with observed pres vs. abs
