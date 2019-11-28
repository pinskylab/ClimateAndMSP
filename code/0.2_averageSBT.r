# Average the bottom temperatures from Morley et al. 2018
# Useful for plotting
# To run on Amphiprion

##############
# Parameters
##############
CLIMPATH <- '/local/shared/pinsky_lab/projections_PlosOne2018/Climate_projection_PlosOne2018'

################
# Functions
################
require(data.table)
require(ggplot2)

#########################
# Run the calculations
########################

# find rcp26 and rcp85 files
files <- list.files(path = CLIMPATH, pattern = '^pred.*rcp[28]',
                    full.names = TRUE)
filesE <- files[grepl('EAST', files)]
filesW <- files[grepl('WEST', files)]

# read in east coast data
# trim to 2007-2019 and average by climate grid
print(length(filesE))
for (i in 1:length(filesE)) {
    cat(i)
    load(filesE[i])
    pred.bathE <- as.data.table(pred.bathE)
    pred <- pred.bathE[year > 2006 & year < 2020, .(SBT = mean(SBT.seasonal)), by = c('latClimgrid', 'lonClimgrid')]
    pred[, mod := gsub(paste0(CLIMPATH, '\\/|predictionEASTnoBias_rcp|85|26|_jas_|\\.RData'), '', filesE[i])]
    pred[, rcp := gsub(paste0(CLIMPATH, '\\/|predictionEASTnoBias_rcp|_jas_|\\.RData|', unique(mod)), '', filesE[i])]
    if (i == 1) aggE <- pred
    if (i != 1) aggE <- rbind(aggE, pred)
    rm(pred.bathE, pred)
}
nrow(aggE) # 196740

# read in west coast data
# trim to 2007-2019 and average by climate grid
print(length(filesW))
for (i in 1:length(filesW)) {
    cat(i)
    load(filesW[i])
    pred.bathW <- as.data.table(pred.bathW)
    pred <- pred.bathW[year > 2006 & year < 2020, .(SBT = mean(SBT.seasonal)), by = c('latClimgrid', 'lonClimgrid')]
    pred[, mod := gsub(paste0(CLIMPATH, '\\/|predictionWESTnoBias_rcp|85|26|_jas_|\\.RData'), '', filesW[i])]
    pred[, rcp := gsub(paste0(CLIMPATH, '\\/|predictionWESTnoBias_rcp|_jas_|\\.RData|', unique(mod)), '', filesW[i])]
    if (i == 1) aggW <- pred
    if (i != 1) aggW <- rbind(aggW, pred)
    rm(pred.bathW, pred)
}

# average across GCMs and rcps and combine
aveE <- aggE[, .(sbt = mean(SBT)), by = c('latClimgrid', 'lonClimgrid')]
aveW <- aggW[, .(sbt = mean(SBT)), by = c('latClimgrid', 'lonClimgrid')]
clim <- rbind(aveE, aveW)

# simple plot
ggplot(clim, aes(x = lonClimgrid, y = latClimgrid, color = sbt)) +
    geom_point(size = 0.05)

# write out
write.csv(clim, file = gzfile('output/climatology.csv.gz'))
