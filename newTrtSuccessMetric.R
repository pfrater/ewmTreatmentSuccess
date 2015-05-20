# this is code for a new metric of treatment success which simply looks at the 
# pre and post ewm treatments divided by the mean ewm abundance post-treatment

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
setwd('C:/Users/fraterp/Documents/DNR/Long_term_EWM_Study')
sf <- read.csv('ewm_master_spp_95_pct_lit_zone_frequency.csv', header=T)
sf$lake <- as.character(sf$lake)
ewm <- subset(sf, select=c(lake, year, ewm))
ewm[is.na(ewm)] <- 0
vars <- read.csv('exploratoryAnalysis/enviroLimnoManagementData.csv', header=T)
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

# calculate post treatment mean for each treatment
post.trt.mean <- aggregate(x=list(post.trt.mean=man$ewm), by=list(lake=man$lake, yr.last.trtd=man$yr.last.trtd), FUN=function(x) mean(x, na.rm=T))

man <- merge(man, post.trt.mean, by=c('lake', 'yr.last.trtd'), all=T)
man <- man[order(man$lake, man$year),]


# here is the algorithm for the new treatment success metric
# ((ewm-pre / ewm-post)*(10*c)) / mean(ewm)
man$trt.success <- NA
new.man <- data.frame(NULL)
for (i in unique(man$lake)) {
	lake <- subset(man, lake==i)
	for (j in 2:nrow(lake)) {
		if (lake$trt.in.year[j]==1) {
			if (lake$ewm[j] != 0) {
				lake$trt.success[j] <- (lake$ewm[j - 1] / lake$ewm[j]) / lake$post.trt.mean[j];
			} else if (lake$post.trt.mean[j] != 0) {
				lake$trt.success[j] <- (lake$ewm[j-1] / 0.1) / lake$post.trt.mean[j];
			} else {lake$trt.success[j] <- (lake$ewm[j-1] / 0.1) / 0.1};
		} else {lake$trt.success[j] <- NA}
	}
	new.man <- rbind(new.man, lake)
}

boxplot(trt.success ~ lake, data=new.man)






























