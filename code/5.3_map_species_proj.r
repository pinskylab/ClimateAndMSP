# Plot the summarized species distributions


############
## Flags
############

# choose the rcps (get to choose one)
rcps <- c(85)

# periods of time
periods <- c('2007-2020', '2021-2040', '2041-2060', '2061-2080', '2081-2100')

# oceans
oceans <- c('Atl', 'Pac')


######################
## Helper functions
######################
require(RColorBrewer)
require(data.table)
require(ggplot2)



#######################################################
# Read in results and simple prep
#######################################################


# loads presence/absence data for each species/model/rcp. slow (5 min?).
if(length(rcps) != 1){
    stop('rcp must be length 1')
} else {
    print(paste0('rcp', rcps))
    
    for(j in 1:length(oceans)){
        for(k in 1:length(periods)){
            print(paste0('ocean', oceans[j], ' period', periods[k]))
            prestemp <- fread(cmd = paste0('gunzip -c temp/presmap_', oceans[j], '_rcp', rcps, '_', periods[k], '.csv.gz'), drop = 1)
            #biotemp <- fread(cmd = paste0('gunzip -c temp/biomassmap_', oceans[j], '_rcp', rcps, '_', periods[k], '.csv.gz'), drop = 1)
            
            if(j == 1 & k == 1){
                presmap <- prestemp
             #   biomassmap <- biotemp
            } else {
                presmap <- rbind(presmap, prestemp)
              #  biomassmap <- rbind(biomassmap, biotemp)
            }
        }
    }
}
rm(prestemp) #, biotemp)

# poccur threshold: how high does the probability of occurrence in the projections need to be to consider the species "present"?
# use the thresholds calculated during model fitting from Morley et al. 2018 PLOS ONE
poccurthresh <- fread('https://raw.githubusercontent.com/pinskylab/project_velocity/master/output/modeldiag_Nov2017_fitallreg_2017.csv', drop = 1)[, .(sppocean, thresh.kappa)]

# region definitions
regiongrid <- fread(cmd = 'gunzip -c output/region_grid.csv.gz', drop = 1)

# Fix lon in regiongrid to match presmap (-360 to 0)
regiongrid[longrid > 0, longrid := longrid - 360] 

# definition of species planning goals by region
sppfiles <- list.files(path = 'output/prioritizr_runs/', pattern = 'spp_*', full.names = TRUE)
spps <- fread(sppfiles[1], drop = 1)
spps[, region := gsub('/|output|prioritizr_runs|spp_|\\.csv', '', sppfiles[1])]
for(i in 2:length(sppfiles)){
    temp <- fread(sppfiles[i], drop = 1)
    temp[, region := gsub('/|output|prioritizr_runs|spp_|\\.csv', '', sppfiles[i])]
    spps <- rbind(spps, temp)
}
rm(temp)

# Add region information to presmap
setkey(presmap, latgrid, longrid)
setkey(regiongrid, latgrid, longrid)
presmap <- merge(presmap, regiongrid[, .(latgrid, longrid, region)], all.x = TRUE) # add region information
if(presmap[is.na(region), .N] != 0){ # 0 missing region: good!
    stop('presmap is missing >0 regions')
}
# presmap[is.na(region) & !duplicated(presmap[,.(latgrid, longrid)]), ]
# presmap[is.na(region) & !duplicated(presmap[,.(latgrid, longrid)]), plot(longrid, latgrid)]


# Add poccur threshold to presmap
poccurthresh[, ocean := gsub('.*_', '', sppocean)]
poccurthresh[, spp := gsub('_Atl|_Pac', '', sppocean)]

presmapPac <- merge(presmap[region %in% c('ebs', 'goa', 'bc', 'wc'), ], poccurthresh[ocean == 'Pac', .(spp, thresh.kappa)], by = 'spp') # have to do Atl and Pac separately since some species are in both regions but use different models
presmapAtl <- merge(presmap[region %in% c('gmex', 'seus', 'neus', 'maritime', 'newf'), ], poccurthresh[ocean == 'Atl', .(spp, thresh.kappa)], by = 'spp')
if(nrow(presmap) == nrow(presmapPac) + nrow(presmapAtl)){
    presmap <- rbind(presmapPac, presmapAtl)
    rm(presmapPac, presmapAtl)
} else {
    stop('merge of poccurthesh and presmap did not work')
}

# Fix a species name
presmap[spp == 'theragra chalcogramma', spp := 'gadus chalcogrammus']


####################################################################################
## Basic maps of species distributions
## poccur
####################################################################################
myregs <- presmap[, sort(unique(region))]

for(j in 1:length(myregs)){
    thisreg <- myregs[j] # trim to a region first!
    thesespps <- spps[region == thisreg & !grepl('fishery|energy', name), unique(spp)]
    numspps <- length(thesespps)
    
    pdf(paste0('temp_figures/species_maps_', thisreg, '_', rcps, '.pdf'), width = 8, height = 20)
    #png(paste0('temp_figures/species_maps_', thisreg, '.png'), width = 8, height = 16, units = 'in', res = 300)
    
    print(numspps)
    for(i in 1:numspps){
        #for(i in 1:2){
        cat(i)
        thissp <- thesespps[i]
        
        p <- ggplot(presmap[spp == thissp & rcp == rcps & region == thisreg, ],
                    aes(x = longrid, y = latgrid, color = poccur > thresh.kappa)) +
            geom_point(size = 0.2) +
            facet_grid(model ~ year_range) +
            theme(strip.text.x = element_text(size = 6), 
                  strip.text.y = element_text(size = 6)) +
            labs(color = thissp)
        print(p)
    }
    
    dev.off()
    cat('\n')
    
}



####################################################################################
## Map of average change in species occurrence
## poccur
####################################################################################

# calculate change in occurrance 2013-2080
preschange2013 <- presmap[rcp == rcps & year_range == '2007-2020', .(poccur2013 = (mean(poccur) > unique(thresh.kappa))), by = c('latgrid', 'longrid', 'region', 'spp')]
preschange2090 <- presmap[rcp == rcps & year_range == '2081-2100', .(poccur2090 = (mean(poccur) > unique(thresh.kappa))), by = c('latgrid', 'longrid', 'region', 'spp')]
preschange <- merge(preschange2013, preschange2090, by = c('latgrid', 'longrid', 'region', 'spp'))
preschange[, diff := poccur2090 - poccur2013]

# plot each region
myregs <- preschange[, sort(unique(region))]

for(j in 1:length(myregs)){
    thisreg <- myregs[j] # trim to a region
    thesespps <- spps[region == thisreg & !grepl('fishery|energy', name), unique(spp)]
    numspps <- length(thesespps)

    pdf(paste0('temp_figures/species_avechange_maps_', thisreg, '_', rcps, '.pdf'), width = 8, height = 6)

    print(numspps)
    for(i in 1:numspps){
        cat(i)
        thissp <- thesespps[i]
        
        p <- ggplot(preschange[region == thisreg & spp == thissp,],
                    aes(x = longrid, y = latgrid, color = diff)) +
            geom_point(size = 1.5) +
            labs(color = thissp) +
            scale_color_gradient2(midpoint = 0, low="red", mid="white",
                                  high="blue", space ="Lab" ) +
            geom_point(data = preschange[region == thisreg & spp == thissp & poccur2013 == TRUE & poccur2090 == TRUE,], 
                   size = 1.5, color = 'black')
        print(p)
    }
    
    dev.off()
    cat('\n')
}


