# I've come up with a very simple metric to measure success in management
# when looking at the time series of ewm within a lake; if we take the coeffecients from a linear regression on ewm by time and multiply them together we can obtain a measure of success whereby the more negative the metric the more successful the management on a particular lake
# success = intercept*slope (a more negative product means greater success)

library(vegan)
library(gclus)
library(cluster)
library(dtw)
library(proxy)
library(reshape2)
library(ggplot2)
library(maps)
library(lubridate)
library(AICcmodavg)
library(grDevices)
library(maptools)
library(rgdal)

# read in data and reorganize to calculate Euclidean or DTW distance on time-series trend
setwd('/Users/paulfrater/Desktop/DNR/Long_term_EWM_Study')
sf <- read.csv('ewm_master_spp_95_pct_lit_zone_frequency.csv', header=T)
sf$lake <- as.character(sf$lake)
ewm <- subset(sf, select=c(lake, year, ewm))
ewm[is.na(ewm)] <- 0
vars <- read.csv('exploratory analysis/enviroLimnoManagementData.csv', header=T)
dat <- merge(ewm, vars, by.x=c('lake'))
man <- subset(dat, managed==1)


# create year of study variable
ewm$time[ewm$year==2005] <- 1
ewm$time[ewm$year==2006] <- 2
ewm$time[ewm$year==2007] <- 3
ewm$time[ewm$year==2008] <- 4
ewm$time[ewm$year==2009] <- 5
ewm$time[ewm$year==2010] <- 6
ewm$time[ewm$year==2011] <- 7
ewm$time[ewm$year==2012] <- 8
ewm$time[ewm$year==2013] <- 9
ewm$time[ewm$year==2014] <- 10

coefs <- data.frame(NULL);
for (i in unique(ewm$lake)) {
	linmod <- lm(ewm$ewm[ewm$lake==i]~ewm$time[ewm$lake==i]);
	lake <- i;
	slope <- coef(linmod)[2];
	int <- coef(linmod)[1];
	r <- cor(ewm$ewm[ewm$lake==i], ewm$time[ewm$lake==i]);
	l <- data.frame(lake, slope, int, r)
	coefs <- rbind(coefs, l)
}
row.names(coefs) <- c(1:28)

man <- merge(man, coefs, by=c('lake'))
man$time <- ewm$time[ewm$lake %in% man$lake]

# add in treatment data in order to calculate max decrease in ewm associated with a treatment
trt <- read.csv('TreatmentData/ewm.ltt.lakes.treatment.csv', header=T, na.strings=c('n/a', 'NA'))
names(trt) <- tolower(names(trt))
trt$date.surveyed <- as.Date(as.character(trt$date.surveyed), format='%m/%d/%Y')
trt$year <- year(trt$date.surveyed)
trt <- subset(trt, select=c(lake, year, date.ewm.discovered, treated.in.year))
trt$lake <- tolower(as.character(trt$lake)); trt$lake <- gsub(' ', '.', trt$lake)
trt$trt.in.year <- 0; trt$trt.in.year[trt$treated.in.year=='Y'] <- 1
man <- merge(man, trt, by=c('lake', 'year'), all.x=T)
man <- man[!duplicated(man[c('lake', 'year')]),]
man$delta.ewm <- 0
new <- NULL;
for (i in unique(man$lake)) {
	e <- subset(man, lake==i);
	for (j in 2:nrow(e)) {
		e$delta.ewm[j] <- e$ewm[j-1] - e$ewm[j]	
	}
	new <- rbind(new, e);
}
man$delta.ewm <- new$delta.ewm


new <- NULL;
for (i in unique(man$lake)){
	e <- subset(man, lake==i);
	e$yr.last.trtd <- NULL;
	yr.last.trtd <- NULL;
	for (j in 1:nrow(e)) {
		if (e$trt.in.year[j]==1) {
			yr.last.trtd <- e$year[j];
			e$yr.last.trtd[j] <- yr.last.trtd;
		} else if (e$trt.in.year[j]==0 & is.null(yr.last.trtd)) {
			e$yr.last.trtd[j] <- NA
		} else {e$yr.last.trtd[j] <- yr.last.trtd}
	}
	new <- rbind(new, e)
}
man$yr.last.trtd <- new$yr.last.trtd
man$yrs.post.trt <- man$year - man$yr.last.trtd

# calculate delta T
trt.sub <- subset(man, trt.in.year==1)
max.d <- aggregate(trt.sub$delta.ewm, by=list(lake=trt.sub$lake), FUN=function(x) max(na.omit(x)))
max.d$x[max.d$x==-Inf] <- 1
names(max.d) <- c('lake', 'max.delta.t')
man <- merge(man, max.d, by=c('lake')) 

# calculate consistency of treatment to maintain low levels of EWM post trt
man$c <- 0;
for (i in 1:nrow(man)) {
	if (man$delta.ewm[i] >= -5) {
		man$c[i] <- man$yrs.post.trt[i]
	}
}
c <- aggregate(man$c, by=list(lake=man$lake), FUN=function(x) max(na.omit(x)))
names(c) <- c('lake', 'max.c')
man <- merge(man, c, by=c('lake'))

# calculate treatment success using the formula:
# trt.success = |intercept|*slope*max.delta.t
man$trt.success <- man$slope * man$max.delta.t - 10*man$max.c
mano <- man[!duplicated(man$lake),]

# first validate the metric by visually checking to see if it looks about right
ts.labs <- data.frame(lake=mano$lake, x=rep(2011, nrow(mano)), y=rep(50, nrow(mano)), labs=round(mano$trt.success, 2))
t <- ggplot(data=man, aes(x=year, y=ewm))
t + geom_point() + geom_smooth(method='lm') + geom_text(aes(x,y,label=labs), data=ts.labs) + facet_wrap(~lake, nrow=4) + scale_y_continuous(limits=c(-10, 100))


# plotting treatment success against lake variables
plot(trt.success~res.dens, data=mano, main='Residential Density') 
plot(trt.success~ph.11, data=mano, main='pH') 
plot(trt.success~cond.25.11, data=mano, main='Conductivity')
plot(trt.success~tds.11, data=mano, main='Total Dissolved Solids')
plot(trt.success~num.x.managed, data=mano, main='# of Times Managed')
plot(trt.success~orp.11, data=mano, main='ORP') # highly significant
plot(trt.success~hwm, data=mano, main='HWM') # the few treatments on hwm seem to be inneffective

# not significant
plot(trt.success~surf.area.acres, data=mano, main='Area') # not significant
plot(trt.success~max.lake.depth, data=mano, main='Depth') # not significant
plot(trt.success~mn.secchi.ft, data=mano, main='Clarity') # not significant
plot(trt.success~mn.totalP, data=mano, main='P') # not significant
plot(trt.success~turbidity.11, data=mano, main='Turbidity') # not significant

########################
# regressions

## nonlinear regressions
t <- mano$trt.success - min(mano$trt.success) + 1
mano$t <- t
f <- function(x, a, k) {a*(1-exp(-k*x))}
mm <- function(x, a, b) {(a*x) / (b + x)}

# res dens
rd.nl.fit <- nls(t~f(x=res.dens, a, k), data=mano, start=c(a=400, k=0.1))
rd.co <- coef(rd.nl.fit)
rd.mm.fit <- nls(t~mm(res.dens, a, b), data=mano, start=c(a=400, b=0.1))
# AIC difference for nl and mm very minimal

# conductivity
cond.nl.fit <- nls(t ~ f(x=cond.25.11, a, k), data=mano, start=c(a=400, k=1))
cond.co <- coef(cond.nl.fit)

# tds
tds.nl.fit <- nls(t~f(x=tds.11, a, k), data=mano, start=c(a=400, k=1))
tds.co <- coef(tds.nl.fit)

# num.x.managed
num.man.fit <- nls(t~f(x=num.x.managed, a, k), data=mano, start=c(a=400, k=0.1))
nm.co <- coef(num.man.fit)


## regression results for significant variables
summary(lm(trt.success~ph.11, data=mano))
summary(rd.nl.fit) # residential density nls regression
summary(cond.nl.fit) # conductivity nls regression
summary(tds.nl.fit) # total dissolved salts nls regression
summary(num.man.fit) # number of times managed nls regression
summary(lm(trt.success~orp.11, data=mano))

# perform AIC calculations for nonlinear models against linear models
AICc(rd.nl.fit) - AICc(lm(t~res.dens, data=mano))
AIC(cond.nl.fit) - AICc(lm(t~cond.25.11, data=mano))
AIC(tds.nl.fit) - AICc(lm(t~tds.11, data=mano))
AIC(num.man.fit) - AICc(lm(t~num.x.managed, data=mano))

## plotting the outputs from nonlinear regression models
# residential density
plot(mano$res.dens, t, xlab=expression(paste('Residential Density (homes mi'[shoreline]^-1,')')), ylab=expression(paste('S'[T], ' - ', 'S'[Tmin], sep='')))
curve(expr=rd.co[1]*(1-exp(-rd.co[2]*x)), add=T, col='red')
text(35, 200, expression(italic(Y==beta[1]*(1-e^paste(-beta[2]*x)))), pos=4)
text(35, 180, 'S = 70.9, d.f.=13', pos=4)

# number of times managed throughout study
plot(mano$num.x.managed, t, xlab='Number of Times Managed \n (throughout study)', ylab=expression(paste('S'[T], ' - ', 'S'[Tmin], sep='')))
curve(expr=nm.co[1]*(1-exp(-nm.co[2]*x)), add=T, col='red')
text(8, 100, expression(italic(Y==beta[1]*(1-e^paste(-beta[2]*x)))), pos=4)
text(8, 80, 'S = 85.36, d.f.=13', pos=4)



######################################
# making a four-paneled figure for the paper
######################################

pdf(file='/Users/paulfrater/Documents/Writing/Papers/EWMTreatmentSuccessPaper/figures/trtSuccessFourPaneledFig.pdf', width=10, height=10)
par(mfrow=c(2,2))
# plot S[t] against pH
plot(trt.success~ph.11, data=mano, xlab='A. pH', ylab=expression(paste('S'[T])))
abline(lm(trt.success~ph.11, data=mano), col='red')
text(8.5, -300, expression(italic(paste(R^2==0.37))), pos=4)
text(8.5, -320, expression(italic(p < 0.01)), pos=4)

#conductivity
plot(mano$cond.25.11, t, xlab=expression(paste("B. Conductivity (",mu,'g l'^'-1',')')), ylab=expression(paste('S'[T], ' - ', 'S'[Tmin], sep='')))
curve(expr=cond.co[1]*(1-exp(-cond.co[2]*x)), add=T, col='red')
text(0.35, 200, expression(italic(Y==beta[1]*(1-e^paste(-beta[2]*x)))), pos=4)
text(0.35, 180, 'S = 80.37, d.f.=13', pos=4)

# ORP
plot(trt.success~orp.11, data=mano, xlab='C. ORP (mV)', ylab=expression('S'[T]))
abline(lm(trt.success~orp.11, data=mano), col='red')
text(225, -20, expression(italic(R^2==0.62)), pos=4)
text(225, -38, expression(italic(p < 0.001)), pos=4)

# total dissolved solids
plot(mano$tds.11, t, xlab=expression(paste("D. Total Dissolved Solids (g l"^'-l',')')), ylab=expression(paste('S'[T], ' - ', 'S'[Tmin], sep='')))
curve(expr=tds.co[1]*(1-exp(-tds.co[2]*x)), add=T, col='red')
text(0.23, 100, expression(italic(Y==beta[1]*(1-e^paste(-beta[2]*x)))), pos=4)
text(0.23, 80, 'S = 78.21, d.f.=13', pos=4)
dev.off()


######################################
# create map of lakes with ecoregions
######################################
# template for creating map
# get lake and map data
wi <- map('state', region=c('wisconsin'))
vars <- read.csv('exploratory analysis/enviroLimnoManagementData.csv', header=T)
spatial <- merge(ewm, vars, by.x=c('lake'))
lakes <- subset(spatial, select=c(lake, lat, long, managed))
lakes$lake <- as.character(lakes$lake)
lakes <- lakes[!duplicated(lakes),]
man.lakes <- subset(lakes, managed==1)
# get ecoregion data
ecoreg <- readOGR('/Users/paulfrater/Documents/Writing/Papers/ewm.native.managed.unmanaged/figures', 'ecoregionsMap')
ecoreg <- spTransform(ecoreg, CRS("+proj=longlat +datum=WGS84"))
wis.eco <- subset(ecoreg, ecoreg@data$STATE_NAME=='Wisconsin')

plot(wis.eco); points(man.lakes$long, man.lakes$lat, pch=21, cex=2, bg=rgb(0,0,0,0.4))











