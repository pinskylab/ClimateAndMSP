# Plot figures for paper


#####################
## Useful functions
#####################
require(maps)
require(mapdata) # for hires world map
require(maptools)
require(rgdal)
require(rgeos)
require(raster)
require(data.table)
require(beanplot)
require(RColorBrewer)
require(lme4)
require(ggsci) # colors for Fig. 4

se <- function(x,na.rm=FALSE){ # standard error
	if(!na.rm){
		return(sd(x, na.rm=FALSE)/sqrt(length(x)))
	}
	if(na.rm){
		return(sd(x, na.rm=TRUE)/sqrt(sum(!is.na(x))))
	}
}

# Function to plot color bar
# from https://stackoverflow.com/questions/9314658/colorbar-from-custom-colorramppalette
color.bar <- function(lut, min, max=-min, nticks=11, ticks=seq(min, max, len=nticks), title='') {
    scale = (length(lut)-1)/(max-min)
    
    plot(c(0,10), c(min,max), type='n', bty='n', xaxt='n', xlab='', yaxt='n', ylab='', main=title)
    axis(4, ticks, las=1)
    for (i in 1:(length(lut)-1)) {
        y = (i-1)/scale + min
        rect(0,y,10,y+1/scale, col=lut[i], border=NA)
    }
}

##############################
# Useful quantities (for paper text)
##############################

# number of species and projections 
presmap1 <- fread(cmd = 'gunzip -c temp/presmap_Atl_rcp26_2007-2020.csv.gz', drop = 1)
presmap2 <- fread(cmd = 'gunzip -c temp/presmap_Pac_rcp26_2007-2020.csv.gz', drop = 1)
biomap1 <- fread(cmd = 'gunzip -c temp/biomassmap_Atl_rcp26_2007-2020.csv.gz', drop = 1)
biomap2 <- fread(cmd = 'gunzip -c temp/biomassmap_Pac_rcp26_2007-2020.csv.gz', drop = 1)
nspp <- presmap1[, length(unique(spp))] + presmap2[, length(unique(spp))] + biomap1[, length(unique(spp))] + biomap2[, length(unique(spp))]
nspp # number of species
8*2*4*nspp # number of projections (8 GCMs, 2 RCPs, 4 time periods for each species)

rm(presmap1, presmap2, biomap1, biomap2)

# number of mpas in each network
networks <- fread('gunzip -c temp/wdpaturnbyMPAbymod.csv.gz', drop = 1) # MPA network descriptions
networks[, length(unique(WDPA_PID)), by = network]

#########################
## Fig 1 Study regions maps
#########################
# turnover statistics
turn <- fread('output/turnover_by_CMSPgrid.csv', drop = 1)

# region definitions
regiongrid <- fread(cmd = 'gunzip -c output/region_grid.csv.gz', drop = 1)
regiongrid[longrid > 0, longrid := longrid - 360] # convert to -360 to 0 for merging with turn
regiongrid[longrid < 0, longrid2 := longrid + 360] # convert to 0 to 360 for plotting a continent-scale map

# merge regions and turnover
turn <- merge(turn, regiongrid, by = c('latgrid', 'longrid'), all.x = TRUE)


# set regions and params for each
regs <- c('ebs', 'goa', 'bc', 'wc', 'gmex', 'seus', 'neus', 'maritime', 'newf')
cexs <- c(ebs = 0.3, goa = 0.3, bc = 0.5, wc = 0.35, gmex = 0.33, seus = 0.4, neus = 0.3, maritime = 0.25, newf = 0.13)
regsnice = c('Eastern Bering Sea', 'Gulf of Alaska', 'British Columbia', 'West Coast U.S.', 'Gulf of Mexico', 'Southeast U.S.', 'Northeast U.S.', 'Maritimes', 'Newfoundland')
regsniceabbrev = c('(EBS)', '(GoA)', '(BC)', '(WC)', '(GMex)', '(SEUS)', '(NEUS)', '(Mar)', '(Newf)')
ylabs = c('Latitude (°N)', '', '', '', 'Latitude (°N)', '', '', 'Latitude (°N)', '')
xlabs = c('', '', '', 'Longitude (°E)', '', '', 'Longitude (°E)', 'Longitude (°E)', 'Longitude (°E)')
ylims = list(ebs = c(51,62.5), goa = c(54, 61), bc = c(48, 54), wc = c(32.2, 48.5), gmex = c(24,30.5), seus = c(25, 35.5), neus = c(35, 45), maritime = c(41, 52), newf = c(42, 62))
xlims = list(ebs = c(-179.5,-155), goa = c(-156, -133), bc = c(-136, -122), wc = c(-126.5, -117), gmex = c(-97.5,-81), seus = c(-82, -74), neus = c(-76.5, -66), maritime = c(-69, -53), newf = c(-68, -43))
pos <- c(ebs = 'right', goa = 'left', bc = 'left', wc = 'left', gmex = 'left', seus = 'right', neus = 'right', maritime = 'right', newf = 'left')
bcol <- 'dark grey' # background color
yfrac <- c(0.1, 0.1, 0.1, 0.05, 0.1, 0.1, 0.1, 0.1, 0.1) # fraction of plot from the bottom that the temperature ranges are written
mfgs <- list(ebs = c(1,2), goa = c(1,3), bc = c(2,1), wc = c(1,4), gmex = c(2,2), seus = c(2,3), neus = c(3,1), maritime = c(3,2), newf = c(3,3)) # plot coordinates for each region
cols <- brewer.pal(9, 'GnBu')
colfun <- colorRamp(colors = cols)
colpal <- colorRampPalette(colors = cols)

#### plot map
# quartz(width=8.7/2.54,height=6/2.54)
#png(width=8.7/2.54, height=6/2.54, units = 'in', res = 300, file='figures/Fig1_study_regions.png')
jpeg(width = 8.7/2.54, height = 6/2.54, file = 'figures/Fig1_study_regions.jpg', units = 'in', res = 1200, quality = 100)
par(mai=c(0.15, 0.08, 0.15, 0.1), omi=c(0.15, 0.2, 0, 0), tck=-0.06, mgp=c(1.2,0.4,0), las=1, cex.main=0.5, cex.axis=0.5)
layout(mat=matrix(c(1,2,3,5,4,6,7,5,8,9,10,11), byrow=TRUE, nrow=3))

# Add continent-scale map
turn[!is.na(beta_sor), plot(1, 1, xlab='', ylab='', main='', xaxt='n', xlim = c(178,310), ylim = c(24, 62))]
axis(1, mgp=c(1.2, 0.02, 0), at=c(200,250,300), labels=c(-160, -110, -60))
map('world2Hires', add=TRUE, xlim=c(170,320), col=bcol, lwd=0.2, resolution=0, fill=FALSE, wrap=TRUE) # annoying that turning fill=TRUE also draws big horizontal lines across the map
turn[!is.na(beta_sor), points(longrid2, latgrid, col=rgb(colfun(beta_sor), maxColorValue = 255), pch=15, cex=0.1)]

# Add each region
for(i in 1:length(regs)){
    cat(paste0(regs[i], ' '))
	inds <- turn[, region==regs[i] & !is.na(beta_sor)]
	turn[inds, plot(longrid, latgrid, col = rgb(colfun(beta_sor), maxColorValue = 255), pch = 15, cex = cexs[i], xlab = '', ylab = '', xlim = xlims[[i]], ylim = ylims[[i]], 
	                main = paste(regsnice[i], '\n', regsniceabbrev[i], sep=''), xaxt='n')]
	axis(1, mgp=c(1.2, 0.02, 0))
	map('worldHires',add=TRUE, col=bcol, fill=TRUE, border=FALSE, resolution=0)
}

mtext('Longitude (°E)', side = 1, outer = TRUE, cex = 0.5)
mtext('Latitude (°N)', side = 2, outer = TRUE, cex = 0.5, las = 0, line = 0.7)


# Add a color bar
par(mgp=c(2,0.5,0), mai=c(0.1, 0.3, 0.1, 0.4), tcl=-0.1)
color.bar(colpal(100), min = 0, max = 1, nticks = 5, title = 'Turnover')


dev.off()


####################################################################
# Fig 2 Compare the planning approaches against each climate projection
####################################################################
goalsmetbymod1 <- fread('output/goalsmetbymod_hist_all.csv', drop = 1)
goalsmetbymod2 <- fread('output/goalsmetbymod_2per_all.csv', drop = 1)
goalsmetbymod1 <- goalsmetbymod1[modeltype == 'testing', ] # remove the planning projections
goalsmetbymod2 <- goalsmetbymod2[modeltype == 'testing', ] # remove the planning projections
setkey(goalsmetbymod1, 'region', 'rcp', 'model', 'year_range')
setkey(goalsmetbymod2, 'region', 'rcp', 'model', 'year_range')
goalsmetbymod1[, type := 'hist']
goalsmetbymod2[, type := '2per']
goalsmetbymod <- rbind(goalsmetbymod1, goalsmetbymod2)
goalsmetbymod[, type := factor(type, levels = c('hist', '2per'))] # set order

myregs <- c('ebs', 'goa', 'bc', 'wc', 'gmex', 'seus', 'neus', 'maritime', 'newf')
regnames = c('Eastern Bering Sea', 'Gulf of Alaska', 'British Columbia', 'West Coast U.S.', 'Gulf of Mexico', 'Southeast U.S.', 'Northeast U.S.', 'Maritimes', 'Newfoundland')
mods <- sort(unique(goalsmetbymod1$model))
rcps <- sort(unique(goalsmetbymod1$rcp))

# Statistics
# goals met, mid-century or end-of-century (for text)
goalsmetbymod[year_range == '2041-2060', .(mean = mean(pmetstrict), sd = sd(pmetstrict)), by = c('type')] # mean and se across models and goals

goalsmetbymod[year_range == '2081-2100' & rcp == '26', .(mean = mean(pmetstrict), sd = sd(pmetstrict)), by = c('type')] # mean and se across models and goals
goalsmetbymod[year_range == '2081-2100' & rcp == '85', .(mean = mean(pmetstrict), sd = sd(pmetstrict)), by = c('type')] # mean and se across models and goals


# p(80% goals met), end-of-century (for text)
goalsmetbymod[year_range == '2081-2100', .(prop = 1 - sum(pmetstrict > 0.7)/.N), by = c('type')]

# goals met by region
goalsmetbymod[year_range == '2081-2100', .(mean = mean(pmet), meanstrict = mean(pmetstrict)), by= c('region', 'type')] # mean (across rcps and GCMs)
goalsmetbymod[year_range == '2081-2100', .(mean=mean(pmet), meanstrict = mean(pmetstrict)), 
              by= c('region', 'model')][, .(mean = mean(mean), se = se(mean), meanstrict = mean(meanstrict), 
                                            sestrict = se(meanstrict)), by = 'region'] # mean and se across models within regions

# statistical test (for text)
mod <- glmer(cbind(nmet, nmet/pmet - nmet) ~ type + (1|region/rcp/model/year_range), data=goalsmetbymod, family='binomial')
summary(mod)
nrow(goalsmetbymod)
cc <- confint(mod, parm="beta_")  # CIs. slow (30 sec)
ctab <- cbind(est = fixef(mod), cc) # add estimates
rtab <- exp(ctab) # to get odds ratios
print(rtab, digits = 3) # print the odds ratios and CIs

# Plot %goals met (hist and 2per solution)
summ <- goalsmetbymod[, .(mid = mean(mid), mean = mean(pmetstrict), lb = mean(pmetstrict) - sd(pmetstrict), ub = mean(pmetstrict) + sd(pmetstrict)), by = c('year_range', 'region', 'type')] # mean and 95%ci across models and goals, for plotting
colmat <- t(col2rgb(brewer.pal(6, 'PuOr'))) # dark red for histonly average, middle for lines, light red for CI. purple for 2per
cols <- rgb(red=colmat[,1], green=colmat[,2], blue=colmat[,3], alpha=c(255, 90, 200, 200, 90, 255), maxColorValue=255)
yaxts <- c('s', 'n', 'n', 's', 'n', 'n', 's', 'n', 'n')
xaxts <- c('n', 'n', 'n', 'n', 'n', 'n', 's', 's', 's')
outfile <- 'figures/Fig2_prioritizr_goalsmetbymod.pdf'; presentonly <- FALSE
# outfile <- 'figures/prioritizr_goalsmetbymod_presentonly.png'; presentonly <- TRUE
outfile
ylims <- c(0, 1)

# quartz(width=8.7/2.54, height=8.7/2.54)
#png(width=8.7/2.54, height=8.7/2.54, units = 'in', res = 300, file=outfile)
pdf(width=8.7/2.54, height=8.7/2.54, file=outfile)
par(mfrow=c(3,3), mai=c(0.05, 0.05, 0.2, 0.05), omi=c(0.4,0.4,0,0), cex.main=0.8, cex.axis=0.6, tcl=-0.15, mgp=c(1.6,0.4,0), las = 1)

for (i in 1:length(myregs)) { # for each region
    plot(0, 0, xlab='', ylab='', ylim=ylims, xlim=c(2030,2090), main=regnames[i], yaxt='n', xaxt='n')
    
    if(yaxts[i]=='s'){
        axis(2, mgp=c(2,0.4,0))
    } else {
        axis(2, labels=FALSE)
    }
    
    if(xaxts[i]=='s'){
        axis(1, mgp=c(2,0.1,0))
    } else {
        axis(1, labels=FALSE)
    }

    # plot histonly lines
    for(k in 1:length(mods)) {
        for(j in 1:length(rcps)) {
            goalsmetbymod1[model == mods[k] & rcp == rcps[j] & region==myregs[i] & year_range != '2007-2020', 
                           points(mid, pmetstrict, type='l', pch=16, lwd = 0.5, col=cols[2])]
        }
    }
    # plot 2per lines
    if(!presentonly){
        for(k in 1:length(mods)){
            for(j in 1:length(rcps)){
                goalsmetbymod2[model == mods[k] & rcp == rcps[j] & region==myregs[i]  & year_range != '2007-2020', 
                               points(mid, pmetstrict, type='l', pch=16, lwd = 0.5, col=cols[5])]
            }
        }
    }

    # plot polygons
    summ[region==myregs[i] & mid > 2014 & type == 'hist', polygon(c(mid, rev(mid)), c(lb, rev(ub)), col = cols[3], border = NA)]
    if(!presentonly) summ[region==myregs[i] & mid > 2014 & type == '2per', polygon(c(mid, rev(mid)), c(lb, rev(ub)), col = cols[4], border = NA)]

    # plot means
    summ[region==myregs[i] & mid > 2014 & type == 'hist', lines(mid, mean, col = cols[1], lwd = 2)]
    if(!presentonly) summ[region==myregs[i] & mid > 2014 & type == '2per', lines(mid, mean, col = cols[6], lwd = 2)]
    
}

mtext(side=1, text='Year', line=1.6, outer=TRUE)
mtext(side=2, text='Fraction goals met', line=1.8, outer=TRUE, las = 0)

dev.off()


	

###################################
## Fig 3 Compare area needed for each plan
###################################

# Read in plans
folder <- 'output/prioritizr_runs'
runnames <- list.files(path = folder, pattern = 'solution')

consplans <- vector('list', length(runnames))
for(i in 1:length(consplans)){
	consplans[[i]] <- fread(paste0(folder, '/', runnames[i]), drop = 1)
	consplans[[i]]$region <- gsub('solution_|2per_|hist_|.csv', '', runnames[i])
	consplans[[i]]$type <- gsub('solution_|_ebs|_goa|_bc|_wc|_gmex|_seus|_neus|_maritime|_newf|.csv', '', runnames[i])
}
consplans <- rbindlist(consplans)
regs <- c('ebs', 'goa', 'bc', 'wc', 'gmex', 'seus', 'neus', 'maritime', 'newf') # set the order

# How many planning units?
nrow(consplans)
consplans[, .N, by = region] # per region

# add a zone indicator to consplans
consplans[ , zone := as.numeric(NA)]
consplans[solution_1_conservation == 1 , zone := 1]
consplans[solution_1_fishery == 1 , zone := 2]
consplans[solution_1_energy == 1 , zone := 3]
consplans[solution_1_conservation ==  0 & solution_1_fishery == 0 & solution_1_energy == 0 , zone := 4]

# Calculate fraction grid cells that change zone
consplansw <- dcast(consplans, latgrid + longrid + region ~ type, value.var = 'zone')
consplansw[, change := hist != `2per`]

consplansw[, .(nchange = sum(change), ntot = .N, pchange = sum(change)/.N), by = region] # average in each region
consplansw[, .(nchange = sum(change), ntot = .N, pchange = sum(change)/.N), 
           by = region][, .(avepchange = mean(nchange/ntot), sd = sd(nchange/ntot))] # average proportion across regions


# Make matrix of proportion in each zone in each region, for barplot
mathist <- matrix(NA, nrow = 4, ncol=length(regs))
colnames(mathist) <- regs
rownames(mathist) <- c('conservation', 'fishery', 'energy', 'free')
mat2per <- mathist
mathistraw <- mathist
mat2perraw <- mathist
for(i in 1:length(regs)){
	for(j in 1:4){
		mathist[j,i] <- consplansw[region == regs[i], sum(hist == j)/.N]
		mat2per[j,i] <- consplansw[region == regs[i], sum(`2per` == j)/.N]
		mathistraw[j,i] <- consplansw[region == regs[i], sum(hist == j)]
		mat2perraw[j,i] <- consplansw[region == regs[i], sum(`2per` == j)]
	}
}
	colSums(mathist)
	colSums(mat2per)
	colSums(mathistraw)
	colSums(mat2perraw)
	
l <- list(mathist=mathist, mat2per=mat2per)
mat <- do.call(cbind, l)[,order(sequence(sapply(l, ncol)))]

# change in zones: 2per - hist
mat2perraw - mathistraw # number
round(mat2per - mathist,3) # fraction
range(mat2per[4,] - mathist[4,]) # fractional decrease in free space (e.g., increase in plan area)
mean(mat2per[4,] - mathist[4,]) # mean fractional decrease in free space (e.g., increase in plan area)
sd(mat2per[4,] - mathist[4,])/sqrt(9) # SE fractional decrease in free space (e.g., increase in plan area)


# modify to trick barplot into letting me use 8 colors instead of 4
# see http://r.789695.n4.nabble.com/barplot-colors-td4662538.html
matmod <- cbind(c(mat[,1], rep(0,nrow(mat))), c(rep(0,nrow(mat)), mat[,2])) 
for(i in seq(3,ncol(mat),by=2)){
	matmod <- cbind(matmod, c(mat[,i], rep(0,nrow(mat))), c(rep(0,nrow(mat)), mat[,i+1])) 
}

# plot
regnames = c('EBS', 'GoA', 'BC', 'WC', 'GMex', 'SEUS', 'NEUS', 'Mar', 'Newf')
cols <- brewer.pal(8, 'RdYlBu')
cols <- cols[c(1:4, 8:5)]

# quartz(width=8.7/2.54, height=5/2.54)
pdf(width=8.7/2.54, height=10/2.54, file='figures/Fig3_planareas.pdf')
#png(width=8.7/2.54, height=10/2.54, units = 'in', res = 300, file='figures/Fig3_planareas.png')
par(mfrow = c(2, 1), mai=c(0.2, 0.75, 0.1, 0.1), mgp=c(1.8, 0.4, 0), tcl=-0.2, las=1, cex.axis=0.8)
barplot(height=matmod, space=rep(c(1,0), ncol(mathist)), xaxt='n', col=cols, ylab='Proportion', xlim=c(1.5,40))
axis(1, at=seq(2,26,by=3), labels=NA, cex.axis=0.8, las=2, mgp=c(1.8, 0.5, 0))
mtext(side = 3, text = 'a)', adj = -0.3, line = -0.5, las = 1, cex = 1)
legend(x=28.5, y=1, fill=cols[8:5], legend=c('Free', 'Energy', 'Fishing', 'Conservation'), cex=0.5, bty='n')
legend(x=27, y=1, fill=cols[4:1], legend=rep("",4), cex=0.5, bty='n')

par(mai=c(0.8, 0.75, 0.05, 0.1), mgp = c(1.8, 0.4, 0))
consplansw[, .(pchange = sum(change)/.N), 
           by = region][, barplot(height = pchange, space = 2, xaxt = 'n', yaxt = 'n', col = 'black', 
                                  ylab = 'Proportion\nchanged', xlim = c(2,40), ylim = c(0, 0.35))]
axis(1, at = seq(2.5,26.5,by=3), labels=regnames, cex.axis=0.8, las=2, mgp=c(1.8, 0.5, 0))
axis(2, at = seq(0, 0.3, by = 0.1))
mtext(side = 3, text = 'b)', adj = -0.3, line = -0.5, las = 1, cex = 1)

dev.off()

# make a table instead
round(mathist[c('conservation', 'fishery', 'energy', 'free'), c('ebs', 'ai', 'goa', 'wc', 'gmex', 'neus', 'scot', 'sgulf', 'newf')],2)

round(mat2per[c('conservation', 'fishery', 'energy', 'free'), c('ebs', 'ai', 'goa', 'wc', 'gmex', 'neus', 'scot', 'sgulf', 'newf')],2)

round(rowMeans(mathist[c('conservation', 'fishery', 'energy', 'free'), c('ebs', 'ai', 'goa', 'wc', 'gmex', 'neus', 'scot', 'sgulf', 'newf')]),2)
round(rowMeans(mat2per[c('conservation', 'fishery', 'energy', 'free'), c('ebs', 'ai', 'goa', 'wc', 'gmex', 'neus', 'scot', 'sgulf', 'newf')]),2)



################################################
# Fig. 4 Plot the prioritizr efficiency frontier
################################################
# read in prioritizr solutions
frontierall <- fread('temp/frontierall_2019-12-31_075440.csv', drop = 1)
frontierall <- frontierall[budget == 0.75, ]
setkey(frontierall, region, budget, presweight)

# definition of planning features by region
sppfiles <- list.files(path = 'output/prioritizr_runs/', pattern = 'spp_*', full.names = TRUE)
spps <- fread(sppfiles[1], drop = 1)
spps[, region := gsub('/|output|prioritizr_runs|spp_|\\.csv', '', sppfiles[1])]
for(i in 2:length(sppfiles)){
    temp <- fread(sppfiles[i], drop = 1)
    temp[, region := gsub('/|output|prioritizr_runs|spp_|\\.csv', '', sppfiles[i])]
    spps <- rbind(spps, temp)
}
rm(temp)
ngoals <- spps[name != 'energy', .(ngoals = .N), by = region]
setkey(ngoals, region)

# set up region order and add ngoals
frontierall[, region := factor(region, levels = c('ebs', 'goa', 'bc', 'wc', 'gmex', 'seus', 'neus', 'maritime', 'newf'))] # set order
frontierall <- merge(frontierall, ngoals, by = 'region')

# how many not optimal?
frontierall[, .(notopt = sum(status != 'OPTIMAL'), total = .N), by = region]


# set up goals as a % of total
frontierall[, ':='(presperc = presgoals/ngoals, futperc = futgoals/ngoals)]

# quick plot as a check
#require(ggplot2)
#ggplot(frontierall, aes(x = presperc, y = futperc, group = budget, color = budget)) +
#    geom_path(size = 0.4) +
#    geom_point(size = 0.3) +
#    facet_wrap(~ region, nrow = 3, scales = 'free')

# Drop non-frontier points (automated)
# If two points share the same futperc, it chooses the one with the higher presperc (or vice versa)
frontierall[, todrop := 0]
for(i in 1:nrow(frontierall)){
    thisreg <- frontierall[i, region]
    thispresperc <- frontierall[i, presperc]
    k1 <- frontierall[, presperc == thispresperc & region == thisreg]
    if(length(unique(frontierall[k1, futperc])) > 1){
        mx <- frontierall[k1, max(futperc)]
        frontierall[presperc == thispresperc & futperc != mx & region == thisreg & presweight != 0 & presweight != 100, todrop := 1]
    }

    thisfutperc <- frontierall[i, futperc]
    k2 <- frontierall[, futperc == thisfutperc & region == thisreg]
    if(length(unique(frontierall[k2, presperc])) > 1){
        mx <- frontierall[k2, max(presperc)]
        frontierall[futperc == thisfutperc & presperc != mx & region == thisreg & presweight != 0 & presweight != 100, todrop := 1]
    }
}

ggplot(frontierall[todrop == 0,], aes(x = presperc, y = futperc, group = budget, color = budget)) +
    geom_path(size = 0.4) +
    geom_point(size = 0.3) +
    facet_wrap(~ region, nrow = 3, scales = 'free')

# Plot %goals met for each weighting and region
#colmat <- t(col2rgb(brewer.pal(9, 'Set1')))
#cols <- rgb(red=colmat[,1], green=colmat[,2], blue=colmat[,3], alpha=c(255,255,255,255), maxColorValue=255)
cols <- pal_lancet(alpha = 1)(9)

yaxts <- c('s', 'n', 'n', 's', 'n', 'n', 's', 'n', 'n')
xaxts <- c('n', 'n', 'n', 'n', 'n', 'n', 's', 's', 's')
myregs <- c('ebs', 'goa', 'bc', 'wc', 'gmex', 'seus', 'neus', 'maritime', 'newf')
regnames = c('Eastern Bering Sea', 'Gulf of Alaska', 'British Columbia', 'West Coast U.S.', 'Gulf of Mexico', 'Southeast U.S.', 'Northeast U.S.', 'Maritimes', 'Newfoundland')
buds <- 0.75 # the budgets to plot

#png('figures/Fig4_prioritizr_frontiers.png', height = 4, width = 6, units = 'in', res = 300)
pdf('figures/Fig4_prioritizr_frontiers.pdf', height = 4, width = 6)
layout(matrix(c(1,4,4,4,2,4,4,4,3,4,4,4), byrow = TRUE, ncol = 4))

par(mai=c(0.1, 0.1, 0.1, 0.05), cex.main = 1, cex.axis = 0.8, tcl = -0.3, mgp=c(2,0.5,0))
plot(x = c(0, 1, 1), y = c(1, 1, 0), type = 'l', bty = 'l', lwd = 2, xaxt = 'n', yaxt = 'n', xlab = '', ylab = '', xlim = c(-0.25, 1.25), ylim = c(-0.25, 1.25))
plot(x = c(0, 1), y = c(1, 0), type = 'l', bty = 'l', lwd = 2, xaxt = 'n', yaxt = 'n', xlab = '', ylab = '', xlim = c(-0.25, 1.25), ylim = c(-0.25, 1.25))

x1 <- seq(0, 1, length = 100)
plot(x = x1, y = sqrt(1 - x1^2), type = 'l', bty = 'l', lwd = 2, xaxt = 'n', yaxt = 'n', xlab = '', ylab = '', xlim = c(-0.25, 1.25), ylim = c(-0.25, 1.25))

par(mai=c(0.7, 0.7, 0.2, 0.05))
plot(0, 0, type = 'o', bty = 'l', pch = 16, col = 'white', xlab='Present goals met (proportion)', ylab='Future goals met (proportion)', ylim = c(0.2, 1), xlim = c(0.3, 1), main='')
for (i in 1:length(myregs)) { # for each region
    frontierall[region == myregs[i] & budget == buds & todrop == 0, lines(presperc, futperc, type = 'l', pch = 16, col = cols[i])]
}

legend('bottomleft', legend = regnames, col = cols, lty = 1, cex = 0.7, title = 'Regions')

dev.off()


###############################################
## Fig. 5 Management area species gains and losses plot and stats
###############################################
wdpaturnbyMPAbymod <- fread('gunzip -c temp/wdpaturnbyMPAbymod.csv.gz', drop = 1) # individual MPA results
wdpaturnbynetbymod <- fread('gunzip -c temp/wdpaturnbynetbymod.csv.gz') # network results

# Prep
    # reshape to long format
    wdpaturnbyMPAbymodl <- melt(wdpaturnbyMPAbymod, id.vars = c('WDPA_PID', 'network'), 
                                measure.vars = patterns('ninit|nfinal|nshared|nlost|ngained')) # convert to long
    wdpaturnbyMPAbymodl[, c('var', 'rcp', 'model') := tstrsplit(variable, "\\.", fixed = FALSE)][, variable := NULL] # extract rcp and GMC number from the col name
    wdpaturnbyMPAbymodl <- dcast(wdpaturnbyMPAbymodl, WDPA_PID + network + rcp + model ~ var) # group ninit -> nshared as separate columns
    
    wdpaturnbynetbymodl <- melt(wdpaturnbynetbymod, id.vars = c('network'),
                                measure.vars = patterns('ninit|nfinal|nshared|nlost|ngained')) # convert to long
    wdpaturnbynetbymodl[, c('var', 'rcp', 'model') := tstrsplit(variable, "\\.", fixed = FALSE)][, variable := NULL] # extract rcp and GMC number from the col name
    wdpaturnbynetbymodl <- dcast(wdpaturnbynetbymodl, network + rcp + model ~ var) # group ninit -> nshared as separate columns
    
    # calculate means within models/rcps (across regions)
    means <- wdpaturnbyMPAbymodl[, .(flost = mean(nlost/ninit), fgained = mean(ngained/nfinal), beta_sor=mean(2*nshared/(2*nshared + ngained + nlost))), 
                                 by=c('rcp', 'model')] 
    
    # calculate network means within climate models/rcps
    means.net <- wdpaturnbynetbymodl[, .(flost = mean(nlost/ninit), fgained = mean(ngained/nfinal), beta_sor=mean(2*nshared/(2*nshared + ngained + nlost))), 
                                     by=c('rcp', 'model')] 
    
    # calculate means within models/rcps for individual MPAs in the networks
    means.net.indiv <- wdpaturnbyMPAbymodl[!is.na(network), .(flost = mean(nlost/ninit), fgained = mean(ngained/nfinal), beta_sor=mean(2*nshared/(2*nshared + ngained + nlost))), 
                                           by=c('rcp', 'model')] 
    
    # combine network and individual MPA results
    means.net$type <- 'net'
    means.net.indiv$type <- 'ind'
    means.net <- rbind(means.net, means.net.indiv)

# Statistics for individual and networks of management zones
    # means and SE
    means.net[, .(flost = mean(flost), flost.sd = sd(flost), fgained = mean(fgained), fgained.sd = sd(fgained), 
                  beta_sor_diss = mean(1 - beta_sor), beta_sor_diss.sd = sd(1 - beta_sor)), 
              by = c('type', 'rcp')]
    
    sort(means.net[type == 'ind', beta_sor] - means.net[type == 'net', beta_sor]) # sorted differences: all negative
    
    # test
    t.test(means.net[type == 'ind', beta_sor], means.net[type == 'net', beta_sor]) # parametric
    wilcox.test(means.net[type == 'ind', beta_sor], means.net[type == 'net', beta_sor]) # non-parametric
    
    # MLPA example
    wdpaturnbyMPAbymodl[network == 'mlpa', .(flost = mean(nlost/ninit), fgained = mean(ngained/nfinal), beta_sor = 1 - mean(2*nshared/(2*nshared + ngained + nlost))), 
                        by = .(rcp, model)][, .(flost = mean(flost), fgained = mean(fgained), beta_sor = mean(beta_sor))] # individual MPAs
    wdpaturnbynetbymodl[network == 'mlpa', .(flost = mean(nlost/ninit), fgained = mean(ngained/nfinal), beta_sor = 1 - mean(2*nshared/(2*nshared + ngained + nlost)))] # individual MPAs
    
    
    
# Plot of mean MPA change and network change
cols <- list(c('#67a9cf','white','black','#2166ac'), c('#ef8a62', 'white','black','#b2182b')) # colors from Colorbrewer2 7-class RdBu
# quartz(width=8.7/2.54,height=8.7/2.54)
pdf(width=8.7/2.54, height=8.7/2.54, file='figures/Fig5_MPA_turnover.pdf')
#png(width=8.7/2.54, height=8.7/2.54, units = 'in', res = 300, file='figures/Fig5_MPA_turnover.png')

par(mfrow=c(2,2), las=2, mai=c(0.5,0.45,0.1, 0.05), omi=c(0,0,0,0), mgp=c(1.6,0.4,0), tcl=-0.2, cex.axis=0.8)

beanplot(flost ~  type + rcp, data = means.net, what = c(0,1,1,0), side = 'both', col = cols, border = NA, wd = 0.18, handlelog = FALSE, 
         names = c('RCP 2.6', 'RCP 8.5'), las = 1,
         at = c(1, 3), log = "", ylim = c(0, 0.5), xlim = c(0, 4), cut = 0.01, ylab = 'Fraction lost')
mtext(side = 3, text = 'a)', adj = -0.4, line = -0.4, las = 1, cex = 0.9)

beanplot(fgained ~ type + rcp, data = means.net, what = c(0,1,1,0), side = 'both', col = cols, border = NA, wd = 0.18, handlelog = FALSE, 
         names = c('RCP 2.6', 'RCP 8.5'), las = 1,
         at = c(1, 3.2), log = "", ylim = c(0, 0.6), xlim = c(0, 4.5), cut = 0.01, ylab = 'Fraction gained')
mtext(side = 3, text = 'b)', adj = -0.4, line = -0.4, las = 1, cex = 0.9)

beanplot(I(1-beta_sor) ~ type + rcp, data = means.net, what = c(0,1,1,0), side = 'both', col = cols, border = NA, wd = 0.18, handlelog = FALSE, 
         names = c('RCP 2.6', 'RCP 8.5'), las = 1,
         at = c(1, 3), log = "", ylim = c(0, 1), xlim = c(0, 4), cut = 0.01, ylab = 'Dissimilarity') # flost
mtext(side = 3, text = 'c)', adj = -0.4, line = -0.4, las = 1, cex = 0.9)

plot(1, 1, bty = 'n', xaxt = 'n', yaxt = 'n', col = 'white', xlab = '', ylab = '')
legend('center', legend = c('Individual', 'Network'), title = 'Type', col = c(cols[[1]][1], cols[[2]][1]), lty = 1, lwd = 10, bty = 'n')

dev.off()



###################################################
## Fig. S1 Efficiency frontiers for other budgets
###################################################
    
# read in prioritizr solutions
frontier1 <- fread('temp/frontierall_2019-12-22_071607.csv', drop = 1) # 50% budget
frontier2 <- fread('temp/frontierall_2019-12-31_075440.csv', drop = 1) # 75% and 90% budgets
frontier2 <- frontier2[budget == 0.9, ]
frontierall <- rbind(frontier1, frontier2)
setkey(frontierall, region, budget, presweight)

# definition of planning features by region
sppfiles <- list.files(path = 'output/prioritizr_runs/', pattern = 'spp_*', full.names = TRUE)
spps <- fread(sppfiles[1], drop = 1)
spps[, region := gsub('/|output|prioritizr_runs|spp_|\\.csv', '', sppfiles[1])]
for(i in 2:length(sppfiles)){
    temp <- fread(sppfiles[i], drop = 1)
    temp[, region := gsub('/|output|prioritizr_runs|spp_|\\.csv', '', sppfiles[i])]
    spps <- rbind(spps, temp)
}
rm(temp)
ngoals <- spps[name != 'energy', .(ngoals = .N), by = region]
setkey(ngoals, region)

# set up region order and add ngoals
frontierall[, region := factor(region, levels = c('ebs', 'goa', 'bc', 'wc', 'gmex', 'seus', 'neus', 'maritime', 'newf'))] # set order
frontierall <- merge(frontierall, ngoals, by = 'region')

# how many not optimal?
frontierall[, .(notopt = sum(status != 'OPTIMAL'), total = .N), by = region]


# set up goals as a % of total
frontierall[, ':='(presperc = presgoals/ngoals, futperc = futgoals/ngoals)]

# quick plot as a check
require(ggplot2)
ggplot(frontierall, aes(x = presperc, y = futperc, group = budget, color = region)) +
    geom_path(size = 0.4) +
    geom_point(size = 0.3) +
    facet_wrap(~ budget, nrow = 3, scales = 'free')

# Drop non-frontier points (automated)
# If two points share the same futperc, it chooses the one with the higher presperc (or vice versa)
frontierall[, todrop := 0]
for(i in 1:nrow(frontierall)){
    thisreg <- frontierall[i, region]
    thisbud <- frontierall[i, budget]
    thispresperc <- frontierall[i, presperc]
    k1 <- frontierall[, presperc == thispresperc & region == thisreg & budget == thisbud]
    if(length(unique(frontierall[k1, futperc])) > 1){
        mx <- frontierall[k1, max(futperc)]
        frontierall[presperc == thispresperc & futperc < mx & region == thisreg & presweight != 0 & presweight != 100, todrop := 1]
    }
    
    thisfutperc <- frontierall[i, futperc]
    k2 <- frontierall[, futperc == thisfutperc & region == thisreg & budget == thisbud]
    if(length(unique(frontierall[k2, presperc])) > 1){
        mx <- frontierall[k2, max(presperc)]
        frontierall[futperc == thisfutperc & presperc < mx & region == thisreg & presweight != 0 & presweight != 100, todrop := 1]
    }
}

ggplot(frontierall[todrop == 0, ], aes(x = presperc, y = futperc, group = budget, color = region)) +
    geom_path(size = 0.4) +
    geom_point(size = 0.3) +
    facet_wrap(~ budget, nrow = 3, scales = 'free')


# Plot %goals met for each weighting and region
#colmat <- t(col2rgb(brewer.pal(9, 'Set1')))
#cols <- rgb(red=colmat[,1], green=colmat[,2], blue=colmat[,3], alpha=c(255,255,255,255), maxColorValue=255)
cols <- pal_lancet(alpha = 1)(9)

yaxts <- c('s', 'n', 'n', 's', 'n', 'n', 's', 'n', 'n')
xaxts <- c('n', 'n', 'n', 'n', 'n', 'n', 's', 's', 's')
myregs <- c('ebs', 'goa', 'bc', 'wc', 'gmex', 'seus', 'neus', 'maritime', 'newf')
regnames = c('Eastern Bering Sea', 'Gulf of Alaska', 'British Columbia', 'West Coast U.S.', 'Gulf of Mexico', 'Southeast U.S.', 'Northeast U.S.', 'Maritimes', 'Newfoundland')
buds <- 0.5 # the budgets to plot

png('figures/FigS1_prioritizr_frontiers.png', height = 8, width = 6, units = 'in', res = 300)
par(mfrow = c(2,1), mai=c(0.7, 0.7, 0.2, 0.05), cex.main = 1, cex.axis = 0.8, tcl = -0.3, mgp=c(2,0.5,0), las = 1)

plot(0, 0, type = 'o', pch = 16, col = 'white', bty = 'l', xlab='Present goals met (proportion)', ylab='Future goals met (proportion)', ylim = c(0.2, 1), xlim = c(0.3, 1), main='')
for (i in 1:length(myregs)) { # for each region
    frontierall[region == myregs[i] & budget == 0.5 & todrop == 0, lines(presperc, futperc, type = 'l', pch = 16, col = cols[i])]
}

legend('bottomleft', legend = regnames, col = cols, lty = 1, cex = 0.7, title = 'Regions')

plot(0, 0, type = 'o', pch = 16, col = 'white', bty = 'l', xlab='Present goals met (proportion)', ylab='Future goals met (proportion)', ylim = c(0.2, 1), xlim = c(0.3, 1), main='')
for (i in 1:length(myregs)) { # for each region
    frontierall[region == myregs[i] & budget == 0.9 & todrop == 0, lines(presperc, futperc, type = 'l', pch = 16, col = cols[i])]
}

dev.off()



#################################################
## Fig. S2: Simulated management area networks
#################################################
randMPAs <- fread('output/randMPAs_byBT.csv', drop = 1) # read in the simulations
stats <- fread('output/MPA_network_stats.csv', drop = 1)

regs <- sort(unique(randMPAs$region))


# plot netwrkturn as color dots (initial plot)
colrmp <- colorRamp(brewer.pal(11, name='Spectral'))

par(mfrow=c(3,3))
for(i in 1:length(regs)){
    inds <- randMPAs$region==regs[i]
    plot(randMPAs$size[inds], randMPAs$temprng[inds], pch=16, col=rgb(colrmp(randMPAs$netwrkturn[inds]), maxColorValue=256), main=regs[i])
}


# plot netwrkturn as averages within grid squares
colrmp <- colorRamp(rev(brewer.pal(11, name='Spectral')))
colpal <- colorRampPalette(rev(brewer.pal(11, name='Spectral')))
xlabs <- c('', '', '', '', '', '', 'Proportion of area in network', 'Proportion of area in network', 'Proportion of area in network')
ylabs <- c('Proportion of thermal\nrange in network', '', '', 'Proportion of thermal\nrange in network', '', '', 'Proportion of thermal\nrange in network', '', '') 
regs <- c('ebs', 'goa', 'bc', 'wc', 'gmex', 'seus', 'neus', 'maritime', 'newf') # set plot order
regsnice = c('Eastern Bering Sea', 'Gulf of Alaska', 'British Columbia', 'West Coast U.S.', 'Gulf of Mexico', 'Southeast U.S.', 'Northeast U.S.', 'Maritimes', 'Newfoundland')

szs <- seq(min(randMPAs$size, na.rm=TRUE), max(randMPAs$size, na.rm=TRUE), length.out=10)
rngs <- seq(min(randMPAs$temprng, na.rm=TRUE), max(randMPAs$temprng, na.rm=TRUE), length.out=10)
szstep <- diff(szs)[1]
tempstep <- diff(rngs)[1]

gridave <- expand.grid(region = regs, size = szs, temprng = rngs, ave = NA)
for(i in 1:nrow(gridave)){
    inds <- randMPAs[, region == gridave$region[i] & abs(temprng - gridave$temprng[i]) < tempstep/2 & abs(size - gridave$size[i]) < szstep]
    gridave$ave[i] <- mean(randMPAs$netwrkturn[inds])		
}


png(units='in', res=300, width=7, height=6, file='figures/FigS2_randMPAs.png')
par(mgp=c(2,0.5,0), mai=c(0.2, 0.3, 0.3, 0.1), omi=c(0.25, 0.3, 0, 0), xpd=NA, tcl=-0.3, las=1)
layout(matrix(1:18, nrow = 3, byrow = TRUE), widths = c(6, 3, 6, 3, 6, 3), heights = c(1, 1, 1))
for(i in 1:length(regs)){
    # construct matrix for plotting dissimilarity and scale to 0-1
    inds <- gridave$region == regs[i] & !is.na(gridave$ave)
    thisdat <- gridave[inds,]
    minz <- min(floor(thisdat$ave*10)/10)
    maxz <- max(ceiling(thisdat$ave*10)/10)
    thisdat$newz <- thisdat$ave - minz
    thisdat$newz <- thisdat$newz/(maxz - minz)
    mat <- as.data.frame(dcast(as.data.table(thisdat), temprng ~ size, value.var='newz'))
    row.names(mat) <- mat$temprng
    mat <- t(as.matrix(mat[,2:ncol(mat)])) # transpose because of the way image handles matrices
    
    # plot
    par(mgp=c(2,0.5,0), mai=c(0.2, 0.3, 0.3, 0.05), tcl=-0.3)
    image(z=mat, x=sort(unique(thisdat$size)), y=sort(unique(thisdat$temprng)), col=colpal(100), 
          main=regsnice[i], xlab=xlabs[i], ylab=ylabs[i])
    
    # plot NA as grey
    matna <- mat
    matna[is.na(mat)] <- 1
    matna[!is.na(mat)] <- NA
    image(z=matna, x=sort(unique(thisdat$size)), y=sort(unique(thisdat$temprng)), col='grey', add=TRUE)
    
    # add dot for empirical network
    inds2 <- as.character(stats$region) == regs[i]
    if(sum(inds2) > 0){
        points(stats$fracsize[inds2], stats$fractemp[inds2], pch=10, cex=2, col='white')
    }
    
    # add color bar in next plot space
    par(mgp=c(2,0.5,0), mai=c(0.1, 0.05, 0.3, 0.5), tcl=-0.1)
    color.bar(colpal(100), min = minz, max = maxz, nticks = 5)
#    color.bar(colpal(100), min = 0, max = 0.4, nticks = 5)
    
    
}

dev.off()



#################################################################################
## Fig. S3: Compare trawl observations vs. SDM predictions in management areas
#################################################################################
wdpa_by_spp_obs_reg <- fread('gunzip -c temp/wdpa_trawlsppobs_byreg.csv.gz', drop = 1) # trawl survey observations of species in MPAs
mpatrawl <- fread('output/MPAvstrawl_NPV_PPV.csv', drop = 1) # NPV and PPV calculations

cols = brewer.pal(4, 'Paired')

png(units='in', res=300, width=8, height=4, file='figures/FigS3_MPA_trawl.png')

par(mfrow = c(1,2), mai = c(1, 1, 0.3, 0.3))
wdpa_by_spp_obs_reg[, hist(nhaul, breaks = seq(0, 450, by = 4), col = 'grey', xlab = 'Number of hauls per management area', main = '', ylab = '')]
mtext(side = 2, 'Frequency', line = 2.8, las = 0)
mtext(side = 3, 'a)', line = -0.5, at = -120, cex = 1.5)

plot(1000, 1000, xlim = c(1, 100), ylim = c(0, 1), log = 'x', xlab = 'Minimum hauls per management area', ylab = 'Predictive value', las = 1, bty = 'n')
thresh[, polygon(c(min, rev(min)), c(npv+npvse, rev(npv-npvse)), col = cols[1], border = NA)]
thresh[, polygon(c(min, rev(min)), c(ppv+ppvse, rev(ppv-ppvse)), col = cols[3], border = NA)]
thresh[, lines(min, npv, type = 'l', col = cols[2])]
thresh[, lines(min, ppv, type ='l', col = cols[4])]
mtext(side = 3, 'b)', line = -0.5, at = 0.25, cex = 1.5)

dev.off()

# overall values (for text)
thresh[min == 1, ]






##############################################################
## Table S1 problem definitions for planning in each region
## Also stats for text (# planning units total)
##############################################################
# Read in plans and species
folder <- 'output/prioritizr_runs'

runnames1 <- list.files(path = folder, pattern = 'solution')
consplans <- vector('list', length(runnames1))
for(i in 1:length(consplans)){
    consplans[[i]] <- fread(paste0(folder, '/', runnames1[i]), drop = 1)
    consplans[[i]]$region <- gsub('solution_|2per_|hist_|.csv', '', runnames[i])
    consplans[[i]]$type <- gsub('solution_|_ebs|_goa|_bc|_wc|_gmex|_seus|_neus|_maritime|_newf|.csv', '', runnames[i])
}
consplans <- rbindlist(consplans)

runnames2 <- list.files(path = folder, pattern = 'spp')
spps <- vector('list', length(runnames2))
for(i in 1:length(spps)){
    spps[[i]] <- fread(paste0(folder, '/', runnames2[i]), drop = 1)
    spps[[i]]$region <- gsub('spp_|.csv', '', runnames2[i])
}
spps <- rbindlist(spps)
spps[, type := 'conservation']
spps[grepl('fishery', name), type := 'fishery']
spps[grepl('energy', name), type := 'energy']

# How many planning units and in each region?
nrow(consplans)
area <- consplans[type == 'hist', .(area = .N), by = .(region)] # per region
area

# How many goals of each type in each region?
spps[, .N, by = .(region, type)]
cons <- spps[type == 'conservation', .(conservation = .N), by = .(region)]
fish <- spps[type == 'fishery', .(fishery = .N), by = .(region)]
spps[type == 'energy', .N, by = .(region)]

# combine into table
tables1 <- merge(area, cons, by = 'region')
tables1 <- merge(tables1, fish, by = 'region')

# order W to E
ord <- data.table(region = c('ebs', 'goa', 'bc', 'wc', 'gmex', 'seus', 'neus', 'maritime', 'newf'), order = 1:9)
tables1 <- merge(tables1, ord)
setorder(tables1, order)

# write out
write.csv(tables1[, .(region, area, conservation, fishery)], 'tables/tableS1.csv', row.names = FALSE)
