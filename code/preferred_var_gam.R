### NEUS Climate Envelope FINAL
setwd("~/Documents/Collaborations/Rutgers/RAM_Shift/")

library(data.table)
library(oddsratio) #use this to get output from plotting gam model without the plot using no.plot
library(mgcv)

### Load prepped data from neus_data_prep.R
neus.z <- readRDS("neus.z.indic.rds")

### Get rid of 2015 because no btemp
neus.z <- neus.z[year!=2015]



### Get preferred ranges of smoothed explanatory variable from binomial gam

smooth.fit <- function(dt, vars=c("btemp", "stemp")){
	num.var <- length(vars)
	rhs <- paste(paste0("s(", vars, ")"), collapse="+")
	form <- paste0("pres2 ~ ",rhs)
	out1 <- gam(as.formula(form), family=binomial, data=dt)
	mod.fit1 <- no.plot(out1)
	return(mod.fit1)
}

dt <- neus.z[spp=="Gadus morhua"]

out <- smooth.fit(dt, vars=c("btemp", "stemp", "depth"))


out.fit <- lapply(out, function(y) data.frame(x=y$x, fit=y$fit, se=y$se, var=y$xlab))
out.fit.df <- do.call("rbind", out.fit)
out.fit.dt <- as.data.table(out.fit.df)
out.fit.dt[,"fit2":=ifelse(fit<0, 0, fit)]
out.fit.dt[,"upper":=fit + 1.96*se]
out.fit.dt[,"lower":= fit - 1.96*se]

preferred.range <- out.fit.dt[,list(mean.var=weighted.mean(x, w=fit2), 
	min.var=min(x[which(fit2>0 & se <1)]), 
	max.var=max(x[which(fit2>0 & se <1)]), 
	pref.var=x[which(fit2==max(fit2) & se < 1)]), by=list(var)] 

