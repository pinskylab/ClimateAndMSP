### Stats and plots on the plans and their performance


############
## Flags
############

# choose the climate models to use for evaluating future planning (others for planning)
#bcc-csm1-1-m, bcc-csm1-1, CanESM2, CCSM4, CESM1-CAM5, CNRM-CM5, GFDL-CM3, GFDL-ESM2M, GFDL-ESM2G, GISS-E2-R, GISS-E2-H, IPSL-CM5A-LR, IPSL-CM5A-MR, MIROC-ESM, MPI-ESM-LR, NorESM1-ME
gcminds <- c(5, 6, 7, 11, 12, 13, 15, 16) # the complement of those in 5.1_prioritizr.r

# choose regions for these runs (for reading in)
myregs <- c('ebs', 'goa', 'bc', 'wc', 'gmex', 'seus', 'neus', 'maritime', 'newf')

# choose file name suffixes
runname1out <- 'hist_all'
runname2out <- '2per_all'


######################
## Helper functions
######################
require(RColorBrewer)
require(data.table)
#require(maps)
require(ggplot2)

lu <- function(x) length(unique(x))


#######################################################
# Read in results and simple prep
#######################################################
# read in the prioritzr solutions
consplan1s <- vector('list', length(myregs))
for(i in 1:length(consplan1s)) consplan1s[[i]] <- fread(paste0('output/prioritizr_runs/solution_hist_', myregs[i], '.csv'), drop=1)
consplan2s <- vector('list', length(myregs))
for(i in 1:length(consplan2s)) consplan2s[[i]] <- fread(paste0('output/prioritizr_runs/solution_2per_', myregs[i], '.csv'), drop=1)

names(consplan1s) <- myregs
names(consplan2s) <- myregs




# add a zone indicator to consplans
for(i in 1:length(consplan1s)){
    consplan1s[[i]][ , zone := as.numeric(NA)]
    consplan1s[[i]][solution_1_conservation == 1 , zone := 1]
    consplan1s[[i]][solution_1_fishery == 1 , zone := 2]
    consplan1s[[i]][solution_1_energy == 1 , zone := 3]
    consplan1s[[i]][solution_1_conservation ==  0 & solution_1_fishery == 0 & solution_1_energy == 0 , zone := 4]
    consplan1s[[i]][, reg := myregs[i]]
}
for(i in 1:length(consplan2s)){
    consplan2s[[i]][ , zone <- as.numeric(NA)]
    consplan2s[[i]][solution_1_conservation == 1 , zone := 1]
    consplan2s[[i]][solution_1_fishery == 1 , zone := 2]
    consplan2s[[i]][solution_1_energy == 1 , zone := 3]
    consplan2s[[i]][solution_1_conservation ==  0 & solution_1_fishery == 0 & solution_1_energy == 0 , zone := 4]
    consplan2s[[i]][, reg := myregs[i]]
}


####################################################################################
## Basic maps of solutions
####################################################################################

# plot map of selected grids (consplan #1 and #2)
cexs = 0.5 # to adjust
# quartz(width=5, height=18)
pdf(width = 5, height = 18, file = 'figures/planmaps.pdf')
par(mfrow = c(length(myregs),2), mai = c(0.5,0.5,0.3, 0.1), las = 1, mgp = c(2,1,0))
for(i in 1:length(consplan1s)){
    xlims <- range(consplan1s[[i]]$longrid)
    ylims <- range(consplan1s[[i]]$latgrid)
    j <- consplan1s[[i]]$solution_1_conservation == 1
    plot(consplan1s[[i]]$longrid[j], consplan1s[[i]]$latgrid[j], pch=16, 
         col='red', cex=cexs, main='Historical only', xlab='', ylab='', xlim=xlims, ylim=ylims) # conservation
    j <- consplan1s[[i]]$solution_1_fishery == 1
    points(consplan1s[[i]]$longrid[j], consplan1s[[i]]$latgrid[j], pch=16, col='blue', cex=cexs) # fishery
    j <- consplan1s[[i]]$solution_1_energy == 1
    points(consplan1s[[i]]$longrid[j], consplan1s[[i]]$latgrid[j], pch=16, col='green', cex=cexs) # energy
    j <- consplan1s[[i]]$solution_1_conservation == 0 & 
        consplan1s[[i]]$solution_1_fishery == 0 & consplan1s[[i]]$solution_1_energy == 0
    points(consplan1s[[i]]$longrid[j], consplan1s[[i]]$latgrid[j], pch=16, col='grey', cex=cexs) # all others
    
    j <- consplan2s[[i]]$solution_1_conservation == 1
    plot(consplan2s[[i]]$longrid[j], consplan2s[[i]]$latgrid[j], pch=16, col='red', cex=cexs, 
         main='Two periods', xlab='', ylab='', xlim=xlims, ylim=ylims) # conservation
    j <- consplan2s[[i]]$solution_1_fishery == 1
    points(consplan2s[[i]]$longrid[j], consplan2s[[i]]$latgrid[j], pch=16, col='blue', cex=cexs) # fishery
    j <- consplan2s[[i]]$solution_1_energy == 1
    points(consplan2s[[i]]$longrid[j], consplan2s[[i]]$latgrid[j], pch=16, col='green', cex=cexs) # energy
    j <- consplan2s[[i]]$solution_1_conservation == 0 & consplan2s[[i]]$solution_1_fishery == 0 & 
        consplan2s[[i]]$solution_1_energy == 0
    points(consplan2s[[i]]$longrid[j], consplan2s[[i]]$latgrid[j], pch=16, col='grey', cex=cexs) # all others
}

legend('topright', legend=c('Conservation', 'Fishery', 'Energy'), col=c('red', 'blue', 'green'), pch=rep(16,3), cex=0.7, title='Zones', bty='n')


dev.off()


############################
## Basic stats on the plans
############################
# have to read in the plans above (but don't need to read in the species distributions)
# climGrid <- read.csv('data/climGrid.csv', row.names=1) # for temperature climatology and depth

# proportion of each region in each zone
grids <- rbindlist(consplan1s)[, .(ngrid = .N), by = c('reg', 'zone')]
grids[, total := sum(ngrid), by  = 'reg'][, prop := round(ngrid/total, 2)]
setkey(grids, reg, zone)
grids
grids[ , .(meanprop = mean(prop), sdprop = sd(prop), seprop = sd(prop)/sqrt(.N)), by = 'zone']

# latitudinal range by zone in each plan
# do 2per plans span a wider latitudinal range?
#	a <- lapply(pusplan1s, FUN=function(x) aggregate(list(latrng_histonly=x$lat), by=list(zone=x$zone), FUN=function(x) max(x)-min(x)))
#	b <- lapply(pusplan2s, FUN=function(x) aggregate(list(latrng_2per=x$lat), by=list(zone=x$zone), FUN=function(x) max(x)-min(x)))
a <- lapply(consplan1s, FUN=function(x) x[ , .(latsd_hist = sd(latgrid), reg = unique(reg)), by = zone])
b <- lapply(consplan2s, FUN=function(x) x[ , .(latsd_2per = sd(latgrid), reg = unique(reg)), by = zone])
a <- do.call('rbind',a)
b <- do.call('rbind',b)
d <- merge(a,b)
d[d$zone==1,] # conservation
d[d$zone==2,] # fishing
d$diff <- d$latsd_2per - d$latsd_hist # positive if 2per spans a wider range
summary(d$diff)
wilcox.test(d$diff[d$zone==1])
wilcox.test(d$diff[d$zone==2])
wilcox.test(d$diff[d$zone==3])
# note: zone 1 conservation, 2 fishery, 3 energy

# temperature and depth range by zone in each plan
# add in climatology temperature
# pusplan1swclim <- pusplan1s
# pusplan2swclim <- pusplan2s
# for(i in 1:length(pusplan1swclim)){
# 	pusplan1swclim[[i]] <- merge(pusplan1s[[i]], climGrid[,c('lat', 'lon', 'bottemp.clim.int', 'surftemp.clim.int', 'depth')], all.x=TRUE)
# 	pusplan2swclim[[i]] <- merge(pusplan2s[[i]], climGrid[,c('lat', 'lon', 'bottemp.clim.int', 'surftemp.clim.int', 'depth')], all.x=TRUE)
# }
# 
# 	# do 2per plans span a wider surface temperature range?
# #	a <- lapply(pusplan1swclim, FUN=function(x) aggregate(list(temprng_histonly=x$surftemp.clim.int), by=list(zone=x$zone), FUN=function(x) max(x, na.rm=TRUE)-min(x, na.rm=TRUE)))
# #	b <- lapply(pusplan2swclim, FUN=function(x) aggregate(list(temprng_2per=x$surftemp.clim.int), by=list(zone=x$zone), FUN=function(x) max(x, na.rm=TRUE)-min(x, na.rm=TRUE)))
# 	a <- lapply(pusplan1swclim, FUN=function(x) aggregate(list(tempsd_histonly=x$surftemp.clim.int), by=list(zone=x$zone), FUN=function(x) sd(x, na.rm=TRUE)))
# 	b <- lapply(pusplan2swclim, FUN=function(x) aggregate(list(tempsd_2per=x$surftemp.clim.int), by=list(zone=x$zone), FUN=function(x) sd(x, na.rm=TRUE)))
# 	a <- do.call('rbind',a)
# 	b <- do.call('rbind',b)
# 	a$reg <- gsub('.[[:digit:]]', '', row.names(a))
# 	b$reg <- gsub('.[[:digit:]]', '', row.names(b))
# 	d <- merge(a,b)
# 	d[d$zone==2,] # conservation
# 	d[d$zone==3,] # fishing
# 	d$diff <- d$temprng_2per - d$temprng_histonly # positive if 2per spans a wider range
# 	summary(d$diff)
# 	wilcox.test(d$diff[d$zone==2])
# 	wilcox.test(d$diff[d$zone==3])
# 
# 	# do 2per plans span a wider bottom temperature range?
# #	a <- lapply(pusplan1swclim, FUN=function(x) aggregate(list(temprng_histonly=x$bottemp.clim.int), by=list(zone=x$zone), FUN=function(x) max(x, na.rm=TRUE)-min(x, na.rm=TRUE)))
# #	b <- lapply(pusplan2swclim, FUN=function(x) aggregate(list(temprng_2per=x$bottemp.clim.int), by=list(zone=x$zone), FUN=function(x) max(x, na.rm=TRUE)-min(x, na.rm=TRUE)))
# 	a <- lapply(pusplan1swclim, FUN=function(x) aggregate(list(tempsd_histonly=x$bottemp.clim.int), by=list(zone=x$zone), FUN=function(x) sd(x, na.rm=TRUE)))
# 	b <- lapply(pusplan2swclim, FUN=function(x) aggregate(list(tempsd_2per=x$bottemp.clim.int), by=list(zone=x$zone), FUN=function(x) sd(x, na.rm=TRUE)))
# 	a <- do.call('rbind',a)
# 	b <- do.call('rbind',b)
# 	a$reg <- gsub('.[[:digit:]]', '', row.names(a))
# 	b$reg <- gsub('.[[:digit:]]', '', row.names(b))
# 	d <- merge(a,b)
# 	d[d$zone==2,] # conservation
# 	d[d$zone==3,] # fishing
# 	d$diff <- d$temprng_2per - d$temprng_histonly # positive if 2per spans a wider range
# 	summary(d$diff)
# 	wilcox.test(d$diff[d$zone==2])
# 	wilcox.test(d$diff[d$zone==3])
# 
# 
# 	# do 2per plans span a wider depth range?
# 	a <- lapply(pusplan1swclim, FUN=function(x) aggregate(list(depthrng_histonly=x$depth), by=list(zone=x$zone), FUN=function(x) max(x, na.rm=TRUE)-min(x, na.rm=TRUE)))
# 	b <- lapply(pusplan2swclim, FUN=function(x) aggregate(list(depthrng_2per=x$depth), by=list(zone=x$zone), FUN=function(x) max(x, na.rm=TRUE)-min(x, na.rm=TRUE)))
# 	a <- do.call('rbind',a)
# 	b <- do.call('rbind',b)
# 	a$reg <- gsub('.[[:digit:]]', '', row.names(a))
# 	b$reg <- gsub('.[[:digit:]]', '', row.names(b))
# 	d <- merge(a,b)
# 	d[d$zone==2,] # conservation
# 	d[d$zone==3,] # fishing
# 	d$diff <- d$temprng_2per - d$temprng_histonly # positive if 2per spans a wider range
# 	summary(d$diff)
# 	wilcox.test(d$diff[d$zone==2])
# 	wilcox.test(d$diff[d$zone==3])
# 



######################################################################
# Stats comparing the planning approaches against all climate models
######################################################################
goalsmetbymod1out <- fread(paste0('output/goalsmetbymod_', runname1out, '.csv'), drop = 1)
goalsmetbymod2out <- fread(paste0('output/goalsmetbymod_', runname2out, '.csv'), drop = 1)


# add plan identifier and combine
goalsmetbymod1out$plan <- 'histonly'
goalsmetbymod2out$plan <- '2per'
goalsmetbymod <- rbind(goalsmetbymod1out, goalsmetbymod2out)


# compare goals met across plans and time periods
goalmetbymodag <- with(rbind(goalsmetbymod1out, goalsmetbymod2out), 
                       aggregate(list(pmet = pmet), by = list(plan = plan, model = model, rcp = rcp, period = year_range), FUN = mean))
with(goalmetbymodag, aggregate(list(ave_goals = pmet), by = list(plan = plan, period = period), FUN = mean))
with(goalmetbymodag, aggregate(list(ave_goals = pmet), by = list(plan = plan, period = period), FUN = function(x) sd(x)/length(x)))
with(goalmetbymodag, aggregate(list(ave_goals = pmet), by = list(plan = plan, period = period), FUN = length))


# test whether histonly meets significantly fewer goals than 2per
all(goalsmetbymod1out$model == goalsmetbymod2out$model & goalsmetbymod1out$period == goalsmetbymod2out$period & goalsmetbymod1out$rcp == goalsmetbymod2out$rcp) # same order?

# simple t-test on number met (pseudoreplicates within regions)
t.test(goalsmetbymod1out$nmet[goalsmetbymod1out$period=='2081-2100'], goalsmetbymod2out$nmet[goalsmetbymod2out$period=='2081-2100'])

# simple t-test on proportions (should use logistic) (also pseudoreplicates within regions)
t.test(goalsmetbymod1out$pmet[goalsmetbymod1out$period=='2081-2100'], goalsmetbymod2out$pmet[goalsmetbymod2out$period=='2081-2100'], paired=FALSE)

# GLMM model (since proportion data) on 2006-2020
require(lme4)
require(car)
mod <- glmer(cbind(nmet, nmet/pmet - nmet) ~ plan + (1|region/rcp/model), data=goalsmetbymod, family='binomial', subset=goalsmetbymod$period=='2006-2020')
summary(mod)
Anova(mod)

# GLMM model (since proportion data) on 2041-2060
require(lme4)
require(car)
goalsmetbymod <- rbind(goalsmetbymod1out, goalsmetbymod2out)
mod <- glmer(cbind(nmet, nmet/pmet - nmet) ~ plan + (1|region/rcp/model), data=goalsmetbymod, family='binomial', subset=goalsmetbymod$period=='2041-2060')
summary(mod)
Anova(mod)

# GLMM model (since proportion data) on 2081-2100
goalsmetbymod <- rbind(goalsmetbymod1out, goalsmetbymod2out)
mod <- glmer(cbind(nmet, nmet/pmet - nmet) ~ plan + (1|region/rcp/model), data=goalsmetbymod, family='binomial', subset=goalsmetbymod$period=='2081-2100')
summary(mod)
Anova(mod)

# GLMM model (since proportion data) across all time periods
goalsmetbymod <- rbind(goalsmetbymod1out, goalsmetbymod2out)
mod <- glmer(cbind(nmet, nmet/pmet - nmet) ~ plan + (1|region/rcp/model/period), data=goalsmetbymod, family='binomial')
summary(mod)
Anova(mod)


######################################################################
# Plots comparing the planning approaches against all climate models
######################################################################
goalsmetbymod1out <- fread(paste0('output/goalsmetbymod_', runname1out, '.csv'), drop = 1)
goalsmetbymod2out <- fread(paste0('output/goalsmetbymod_', runname2out, '.csv'), drop = 1)


mods <- sort(unique(goalsmetbymod1out$model))
rcps <- sort(unique(goalsmetbymod1out$rcp))
regnames<-c('E. Bering Sea', 'Gulf of AK', 'British Columbia', 'West Coast US', 'Gulf of Mexico', 'Southeast US', 'Northeast US', 'Maritimes', 'Newfoundland')

# plot %goals met (solution #1 and #2) for non-planning models (the test set)
    thesemods <- mods[gcminds]; type <- 'testing' # for models not used in planning
    #thesemods <- mods[-gcminds]; type <- 'training' # for models used in planning
    
    colmat <- t(col2rgb(brewer.pal(4, 'Paired')))
    #	colmat <- t(col2rgb(brewer.pal(4, 'Dark2')))
    cols <- rgb(red=colmat[,1], green=colmat[,2], blue=colmat[,3], alpha=c(255,255,255,255), maxColorValue=255)
    if(type == 'testing') outfile <- paste('figures/prioritizr_allregs_goalsmetbymod_', runname1out, '&', runname2out, '.pdf', sep='')
    if(type == 'training') outfile <- paste('temp_figures/prioritizr_allregs_goalsmetbymod_', runname1out, '&', runname2out, '_training.pdf', sep='')
    outfile
    ylims <- c(0,1)
    
    # quartz(width=6, height=6)
    pdf(width=6, height=6, file=outfile)
    par(mfrow=c(3,3), mai=c(0.3, 0.35, 0.2, 0.05), omi=c(0.2,0.2,0,0), cex.main=0.8)
    
    for(i in 1:length(myregs)){
        # each model/rcp for hist only plans
        inds <- goalsmetbymod1out$model == thesemods[1] & goalsmetbymod1out$rcp == rcps[1] & goalsmetbymod1out$region==myregs[i]
        plot(goalsmetbymod1out$mid[inds], goalsmetbymod1out$pmet[inds], xlab='', ylab='', ylim=ylims, type='l', pch=16, las=1, col=cols[1], main=regnames[i])
        for(k in 1:length(thesemods)){
            for(j in 1:length(rcps)){
                if(!(k==1 & j==1)){
                    inds <- goalsmetbymod1out$model == thesemods[k] & goalsmetbymod1out$rcp == rcps[j] & goalsmetbymod1out$region==myregs[i]
                    points(goalsmetbymod1out$mid[inds], goalsmetbymod1out$pmet[inds], type='l', pch=16, col=cols[1])	
                }
            }
        }
        
        # each model/rcp for 2per plans
        for(k in 1:length(thesemods)){
            for(j in 1:length(rcps)){
                inds <- goalsmetbymod2out$model == thesemods[k] & goalsmetbymod2out$rcp == rcps[j] & goalsmetbymod2out$region==myregs[i]
                points(goalsmetbymod2out$mid[inds], goalsmetbymod2out$pmet[inds], type='l', pch=16, col=cols[3])
            }	
        }
        
        # average across models/rcps for hist and 2per plans
        inds <- goalsmetbymod1out$region==myregs[i] & goalsmetbymod1out$model %in% thesemods
        ensmean <- aggregate(list(nmet=goalsmetbymod1out$nmet[inds], pmet=goalsmetbymod1out$pmet[inds]), by=list(mid=goalsmetbymod1out$mid[inds]), FUN=mean)
        lines(ensmean$mid, ensmean$pmet, col=cols[2], lwd=3)
        
        inds <- goalsmetbymod2out$region==myregs[i] & goalsmetbymod2out$model %in% thesemods
        ensmean2 <- aggregate(list(nmet=goalsmetbymod2out$nmet[inds], pmet=goalsmetbymod2out$pmet[inds]), by=list(mid=goalsmetbymod2out$mid[inds]), FUN=mean)
        lines(ensmean2$mid, ensmean2$pmet, col=cols[4], lwd=3)
        
        # legend in plot 1
        if(i == 1) legend('bottomright', legend = c('hist', '2per'), col = cols[c(2,4)], lwd = 1)
    }
    mtext(side=1,text='Year',line=0.5, outer=TRUE)
    mtext(side=2,text='Fraction goals met',line=0.3, outer=TRUE)
    
    dev.off()


## probability of < 80% goals met by region, plan, and time period
# calcs
u80func <- function(x) return(sum(x<0.8)/length(x))	
u80 <- with(goalsmetbymod, aggregate(list(u80=pmet), by=list(region=region, period=year_range, plan=plan), FUN=u80func))


# examine
u80[u80$period=='2041-2060',] # middle of century, by region
u80[u80$period=='2081-2100',] # end of century, by region

with(u80[u80$period=='2007-2020',], aggregate(list(meanu80=u80), by=list(plan=plan), FUN=mean)) # means by plan (across regions), start
with(u80[u80$period=='2041-2060',], aggregate(list(meanu80=u80), by=list(plan=plan), FUN=mean)) # means by plan (across regions), middle of century
with(u80[u80$period=='2081-2100',], aggregate(list(meanu80=u80), by=list(plan=plan), FUN=mean)) # means by plan (across regions), end of century

a<-aggregate(list(max_u80=u80$u80), by=list(region=u80$region, plan=u80$plan), FUN=max) # max probability of being under 80% of goals met (look across time periods)
mean(a$max_u80[a$plan=='2per'])
mean(a$max_u80[a$plan=='histonly'])

# plot
ggplot(u80, aes(x = period, y = u80, group = plan, color = plan)) +
    geom_line() + 
    geom_point() +
    facet_wrap(~ region, ncol = 3)




######################################################################
# Plots comparing the planning approaches against climate model ensembles
######################################################################
goalsmetbyens1 <- fread(paste0('output/goalsmetbyensemble_', runname1out, '.csv'), drop = 1)
goalsmetbyens2 <- fread(paste0('output/goalsmetbyensemble_', runname2out, '.csv'), drop = 1)


mods <- sort(unique(goalsmetbyens1$model))
rcps <- sort(unique(goalsmetbyens1$rcp))
regnames<-c('E. Bering Sea', 'Gulf of AK', 'British Columbia', 'West Coast US', 'Gulf of Mexico', 'Southeast US', 'Northeast US', 'Maritimes', 'Newfoundland')

# plot %goals met (solution #1 and #2) for non-planning models (the test set)
    colmat <- t(col2rgb(brewer.pal(8, 'Paired')))
    #	colmat <- t(col2rgb(brewer.pal(4, 'Dark2')))
    cols <- rgb(red=colmat[,1], green=colmat[,2], blue=colmat[,3], alpha=c(255,255,255,255), maxColorValue=255)
    outfile <- paste('temp_figures/prioritizr_allregs_goalsmetbyensemble_', runname1out, '&', runname2out, '.pdf', sep='')
    outfile
    ylims <- c(0,1)
    
    # quartz(width=6, height=6)
    pdf(width=6, height=6, file=outfile)
    par(mfrow=c(3,3), mai=c(0.3, 0.35, 0.2, 0.05), omi=c(0.2,0.2,0,0), cex.main=0.8)
    
    for(i in 1:length(myregs)){
        # each model/rcp for hist only plans
        colind <- 1
        goalsmetbyens1[model == mods[1] & rcp == rcps[1] & region==myregs[i], 
                       plot(mid, pmet, xlab='', ylab='', ylim=ylims, type='l', pch=16, las=1, 
                            col=cols[colind], main=regnames[i])]
        colnms <- paste0(mods[1], rcps[1], 'hist')
        for(k in 1:length(mods)){
            for(j in 1:length(rcps)){
                if(!(k==1 & j==1)){
                    colind <- colind + 1
                    goalsmetbyens1[model == mods[k] & rcp == rcps[j] & region==myregs[i], 
                                   points(mid, pmet, xlab='', ylab='', type='l', pch=16, col=cols[colind])]
                    colnms <- c(colnms, paste0(mods[k], rcps[j], 'hist'))
                }
            }
        }

        # each model/rcp for 2per plans
        colind <- colind + 1
        goalsmetbyens2[model == mods[1] & rcp == rcps[1] & region==myregs[i], 
                       points(mid, pmet, type='l', pch=16, col=cols[colind])]
        colnms <- c(colnms, paste0(mods[1], rcps[1], '2per'))
        for(k in 1:length(mods)){
            for(j in 1:length(rcps)){
                if(!(k==1 & j==1)){
                    colind <- colind + 1
                    goalsmetbyens2[model == mods[k] & rcp == rcps[j] & region==myregs[i], 
                                   points(mid, pmet, xlab='', ylab='', type='l', pch=16, col=cols[colind])]
                    colnms <- c(colnms, paste0(mods[k], rcps[j], '2per'))
                }
            }
        }
        
        # legend in plot 1
        if(i == 1) legend('bottomleft', legend = colnms, col = cols, lwd = 1, cex=0.5)
    }
    mtext(side=1,text='Year',line=0.5, outer=TRUE)
    mtext(side=2,text='Fraction goals met',line=0.3, outer=TRUE)
    
    dev.off()
