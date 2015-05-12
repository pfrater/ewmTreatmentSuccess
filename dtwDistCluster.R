library(vegan)
library(dtw)
library(cluster)
library(gclus)
library(maps)
library(impute)


# this is additional analysis for the long-term ewm treatment success paper
# you are going to cluster managed lakes based on dynamic time warping distance, then you will place these in an ordination (pca? pcoa? nmds maybe?) and color sites to cluster (maybe outline them), then overlay environmental variables using envfit() to find correlations with the most successful treatments

###########################################################################
# exploring dynamic time warping to create assocation matrix and cluster lakes
###########################################################################

setwd('/Users/paulfrater/Desktop/DNR/Long_term_EWM_Study')

# create the data set here and subset to just look at ewm
# the next 15 lines or so are just to assemble data in a format for clustering
pi.data <- read.csv('ewm_master_spp_lit_zone_frequency.csv', header=T)
pi.data$lake <- as.character(pi.data$lake)
lakes <- unique(pi.data$lake)
years <- unique(pi.data$year)
sample.time.grid <- expand.grid(lake=lakes, year=years)
sample.time.grid <- sample.time.grid[order(sample.time.grid$lake, sample.time.grid$year),]
lake.info <- read.csv('lake.info.csv', header=T)
managed <- read.csv('Management Data/years managed.csv', header=T)
managed <- subset(managed, select=c(lake, managed))
managed <- subset(managed, managed==1)
envLimMang <- read.csv('exploratory analysis/enviroLimnoManagementData.csv', header=T)
envLimMang$surf.area.ha <- envLimMang$surf.area.acres*0.40469
envLimMang$secchi.m <- envLimMang$mn.secchi.ft*0.3076923
envLimMang$max.depth.m <- envLimMang$max.lake.depth*0.3076923
envLimMang$shoreline.km <- envLimMang$shoreline*1.609344
water.quality <- subset(envLimMang, select=c(lake, surf.area.ha, shoreline.km, max.depth.m, bldgs, secchi.m, mn.chlorophyll, mn.totalP, mn.lake.level, mn.temp, ph.11, orp.11, turbidity.11, do.11, tds.11, cond.25.11, sand, rock, muck, mn.max.depth.plants))
site.spp <- pi.data; rownames(site.spp) <- paste(pi.data$lake, pi.data$year, sep='.'); site.spp <- subset(site.spp, select=c(ewm:length(site.spp)))
ewm <- subset(pi.data, select=c(lake, year, ewm))
ewm[is.na(ewm)] <- 0
ewm <- merge(sample.time.grid, ewm, by=c('lake', 'year'), all.x=T)
ewm <- ewm[ewm$lake %in% managed$lake,]
ewm.wide <- reshape(ewm, timevar='year', idvar='lake', direction='wide')
row.names(ewm.wide) <- ewm.wide$lake; ewm.wide <- subset(ewm.wide, select= -lake)


# there are missing years for some lakes, so first we must impute some of the data
# for lakes that did not have ewm in first sampling year, we will assume it was not present in years prior to this
ewm[11:12,3] <- 0 # berry 2005 & 2006


# now impute data to fill in missing gaps
ewm.impute <- impute.knn(as.matrix(ewm.wide), k=5)
ewm.impute <- as.data.frame(ewm.impute$data)

# plot imputed data to see if there are oddballs
imputed.matrix <- as.matrix(ewm.impute)
for (i in 1:nrow(ewm.impute)) {
	plot(imputed.matrix[i,], main=rownames(ewm.impute)[i], ylim=c(0,60));
	lines(imputed.matrix[i,]);
}

# try running dynamic time warp distance on imputed data
ewm.dist <- dist(ewm.impute, method='DTW') # it worked!

# cluster using various methods
ewm.cl.single <- hclust(ewm.dist, method='single')
ewm.cl.complete <- hclust(ewm.dist, method='complete')
ewm.cl.avg <- hclust(ewm.dist, method='average')
ewm.cl.cent <- hclust(ewm.dist, method='centroid')
ewm.cl.ward <- hclust(ewm.dist, method='ward')


# checking the cophenetic correlation for the various clustering methods
ewm.cl.single.coph <- cophenetic(ewm.cl.single)
single.coph.cor <- round(cor(ewm.dist, ewm.cl.single.coph), 3)
ewm.cl.complete.coph <- cophenetic(ewm.cl.complete)
complete.coph.cor <- round(cor(ewm.dist, ewm.cl.complete.coph), 3)
ewm.cl.avg.coph <- cophenetic(ewm.cl.avg)
avg.coph.cor <- round(cor(ewm.dist, ewm.cl.avg.coph), 3)
ewm.cl.cent.coph <- cophenetic(ewm.cl.cent)
cent.coph.cor <- round(cor(ewm.dist, ewm.cl.cent.coph), 3)
ewm.cl.ward.coph <- cophenetic(ewm.cl.ward)
ward.coph.cor <- round(cor(ewm.dist, ewm.cl.ward.coph), 3)


# based on these Shepard-like plots the average agglomerative clustering looks really solid
par(mfrow=c(2,2))
plot(ewm.dist, ewm.cl.single.coph, xlab='DTW Distance', ylab='Cophenetic Distance', asp=1, main=c('Single Linkage', paste('Cophenetic correlation', single.coph.cor)))
abline(0,1)
lines(lowess(ewm.dist, ewm.cl.single.coph), col='red')
plot(ewm.dist, ewm.cl.complete.coph, xlab='DTW Distance', ylab='Cophenetic Distance', asp=1, main=c('Complete Linkage', paste('Cophenetic correlation', complete.coph.cor)))
abline(0,1)
lines(lowess(ewm.dist, ewm.cl.complete.coph), col='red')
plot(ewm.dist, ewm.cl.avg.coph, xlab='DTW Distance', ylab='Cophenetic Distance', asp=1, main=c('UPGMA', paste('Cophenetic correlation', avg.coph.cor)))
abline(0,1)
lines(lowess(ewm.dist, ewm.cl.avg.coph), col='red')
plot(ewm.dist, ewm.cl.cent.coph, xlab='DTW Distance', ylab='Cophenetic Distance', asp=1, main=c('Centroid Average Clustering', paste('Cophenetic correlation', cent.coph.cor)))
abline(0,1)
lines(lowess(ewm.dist, ewm.cl.cent.coph), col='red')
plot(ewm.dist, ewm.cl.ward.coph, xlab='DTW Distance', ylab='Cophenetic Distance', asp=1, main=c('Ward Clustering', paste('Cophenetic correlation', ward.coph.cor)))
abline(0,1)
lines(lowess(ewm.dist, ewm.cl.ward.coph), col='red')


# trying Gower distance to see what clustering method that says is the best
# average agglomerative is best here as well
gower.single <- sum((ewm.dist - ewm.cl.single.coph)^2)
gower.complete <- sum((ewm.dist - ewm.cl.complete.coph)^2)
gower.avg <- sum((ewm.dist - ewm.cl.avg.coph)^2)
gower.cent <- sum((ewm.dist - ewm.cl.cent.coph)^2)
gower.ward <- sum((ewm.dist - ewm.cl.ward.coph)^2)
gower <- data.frame(gower.single, gower.complete, gower.avg, gower.cent, gower.ward)
which.min(gower)

# plotting the fusion level values for various groups
plot(ewm.cl.single$height, nrow(ewm.wide):2, type='S', main='Fusion levels - DTW - Single', ylab='k (number of clusters)', xlab='h (node height)', col='grey')
text(ewm.cl.single$height, nrow(ewm.wide):2, col='red', cex=0.8)
plot(ewm.cl.complete$height, nrow(ewm.wide):2, type='S', main='Fusion levels - DTW - Complete', ylab='k (number of clusters)', xlab='h (node height)', col='grey')
text(ewm.cl.complete$height, nrow(ewm.wide):2, col='red', cex=0.8)
plot(ewm.cl.avg$height, nrow(ewm.wide):2, type='S', main='Fusion levels - DTW - UPGMA', ylab='k (number of clusters)', xlab='h (node height)', col='grey')
text(ewm.cl.avg$height, nrow(ewm.wide):2, col='red', cex=0.8)
plot(ewm.cl.cent$height, nrow(ewm.wide):2, type='S', main='Fusion levels - DTW - Centroid', ylab='k (number of clusters)', xlab='h (node height)', col='grey')
text(ewm.cl.cent$height, nrow(ewm.wide):2, col='red', cex=0.8)
plot(ewm.cl.ward$height, nrow(ewm.wide):2, type='S', main='Fusion levels - DTW - Ward', ylab='k (number of clusters)', xlab='h (node height)', col='grey')
text(ewm.cl.ward$height, nrow(ewm.wide):2, col='red', cex=0.8)

# cut dendrograms to somewhat arbitrary numbers of groups
k <- 5
ewm.single.cut <- cutree(ewm.cl.single, k)
ewm.complete.cut <- cutree(ewm.cl.complete, k)
ewm.avg.cut <- cutree(ewm.cl.avg, k)
ewm.cent.cut <- cutree(ewm.cl.cent, k)
ewm.ward.cut <- cutree(ewm.cl.ward, k)

# create contingency tables to compare classifications
table(ewm.single.cut, ewm.complete.cut)
table(ewm.single.cut, ewm.avg.cut)
table(ewm.single.cut, ewm.cent.cut)
table(ewm.single.cut, ewm.ward.cut)
table(ewm.complete.cut, ewm.avg.cut)
table(ewm.complete.cut, ewm.ward.cut)
table(ewm.complete.cut, ewm.cent.cut)
table(ewm.avg.cut, ewm.cent.cut)
table(ewm.avg.cut, ewm.ward.cut)
table(ewm.cent.cut, ewm.ward.cut)


# plotting silhouette widths to look for optimal number of groups
# single linkage
asw <- numeric(nrow(ewm.wide))
for (i in 2:nrow(ewm.wide)) {
	sil <- silhouette(cutree(ewm.cl.single, k=i), ewm.dist)
	asw[i] <- summary(sil)$avg.width
}
k.best.single <- which.max(asw)
plot(1:nrow(ewm.wide), asw, type='h', main='Silhouette-optimal number of clusters-Single Linkage', xlab='k (number of groups)', ylab='Average silhouette widths')
axis(1, k.best, paste('optimum', k.best, sep='\n'), col='red', font=2, col.axis='red')
points(k.best.single, max(asw), pch=16, col='red', cex=1.5)
# complete linkage
asw <- numeric(nrow(ewm.wide))
for (i in 2:nrow(ewm.wide)) {
	sil <- silhouette(cutree(ewm.cl.complete, k=i), ewm.dist)
	asw[i] <- summary(sil)$avg.width
}
k.best.complete <- which.max(asw)
plot(1:nrow(ewm.wide), asw, type='h', main='Silhouette-optimal number of clusters-Complete Linkage', xlab='k (number of groups)', ylab='Average silhouette widths')
axis(1, k.best.complete, paste('optimum', k.best.complete, sep='\n'), col='red', font=2, col.axis='red')
points(k.best.complete, max(asw), pch=16, col='red', cex=1.5)
# UPGMA clustering
asw <- numeric(nrow(ewm.wide))
for (i in 2:nrow(ewm.wide)) {
	sil <- silhouette(cutree(ewm.cl.avg, k=i), ewm.dist)
	asw[i] <- summary(sil)$avg.width
}
k.best.avg <- which.max(asw)
plot(1:nrow(ewm.wide), asw, type='h', main='Silhouette-optimal number of clusters-UPGMA', xlab='k (number of groups)', ylab='Average silhouette widths')
axis(1, k.best.avg, paste('optimum', k.best.avg, sep='\n'), col='red', font=2, col.axis='red')
points(k.best.avg, max(asw), pch=16, col='red', cex=1.5)
# centroid clustering
asw <- numeric(nrow(ewm.wide))
for (i in 2:nrow(ewm.wide)) {
	sil <- silhouette(cutree(ewm.cl.cent, k=i), ewm.dist)
	asw[i] <- summary(sil)$avg.width
}
k.best.cent <- which.max(asw)
plot(1:nrow(ewm.wide), asw, type='h', main='Silhouette-optimal number of clusters-Centroid', xlab='k (number of groups)', ylab='Average silhouette widths')
axis(1, k.best.cent, paste('optimum', k.best.cent, sep='\n'), col='red', font=2, col.axis='red')
points(k.best.cent, max(asw), pch=16, col='red', cex=1.5)
# ward clustering
asw <- numeric(nrow(ewm.wide))
for (i in 2:nrow(ewm.wide)) {
	sil <- silhouette(cutree(ewm.cl.ward, k=i), ewm.dist)
	asw[i] <- summary(sil)$avg.width
}
k.best.ward <- which.max(asw)
plot(1:nrow(ewm.wide), asw, type='h', main='Silhouette-optimal number of clusters-Ward', xlab='k (number of groups)', ylab='Average silhouette widths')
axis(1, k.best.ward, paste('optimum', k.best.cent, sep='\n'), col='red', font=2, col.axis='red')
points(k.best.ward, max(asw), pch=16, col='red', cex=1.5)
##### all of the previous plot show that 2 is the best number of clusters; however, this is not very ecologically interesting

par(mfrow=c(2,2))
plot(ewm.cl.single)
plot(ewm.cl.complete)
plot(ewm.cl.avg)
plot(ewm.cl.cent)


# after reviewing some of the various clustering methods, I feel like the complete linkage tells a really compelling story with what is going on in these lakes. More so than any of the other methods; UPGMA is also decent

k <- 4
cutg <- cutree(ewm.cl.complete, k=k)
sil <- silhouette(cutg, ewm.dist)
silo <- sortSilhouette(sil)
rownames(silo) <- row.names(ewm.wide)[attr(silo, 'iOrd')]
plot(silo, main='Silhouette plot - DTW Distance - Complete', cex.names=0.8, col=cutg+1, nmax.lab=100)

# plot a final dendrogram of the cluster analysis
ewm.cl.reorder <- reorder.hclust(ewm.cl.complete, ewm.dist)
plot(ewm.cl.reorder, hang=-1, xlab='5 groups', sub='', ylab='Height', main='DTW Distance - Complete (reordered)')
rect.hclust(ewm.cl.reorder, k=k, border=c('blue', 'green', 'red', 'brown', 'purple'))

# map of five complete linkage clusters
map('state', 'Wisconsin')
points(latitude~longitude, data=lake.info[lake.info$lake %in% managed$lake,], xlab='longitude', ylab='latitude', pch=cutg+15, cex=cutg, col=cutg)
title(main='Four Complete Linkage Groups')

# plot the time series ewm data within each groups
for (i in 1:nrow(ewm)) {
	ewm$comp.cluster[i] <- cutg[names(cutg)==ewm$lake[i]]
}
ewm$comp.cluster[ewm$comp.cluster==1] <- 'Low Level'
ewm$comp.cluster[ewm$comp.cluster==2] <- 'Moderate Decrease'
ewm$comp.cluster[ewm$comp.cluster==3] <- 'Highly Successful'
ewm$comp.cluster[ewm$comp.cluster==4] <- 'Unsuccessful'

ewm.gg <- ggplot(data=ewm, aes(x=year, y=ewm))
ewm.gg + geom_line(aes(color=lake)) + facet_wrap(~comp.cluster) + theme_bw() + theme(legend.position='none') + labs(title='Complete Linkage Clustering on DTW Distance \n', y='% FOO EWM \n', x='\n Year') + scale_x_continuous(breaks=c(2006, 2008, 2010, 2012, 2014)) + theme(strip.background=element_rect(fill=colors))

myStripStyle <- function(which.panel, factor.levels, ...) {
    panel.rect(0, 0, 1, 1, col = colores[which.panel],border = 1);
    panel.text(x=0.5, y=0.5, lab=factor.levels[which.panel], col='black')
}
panelColors <- c('grey10', 'grey20', 'grey30', 'grey40', 'grey50', 'grey60', 'grey70', 'grey80', 'grey90')
colores <- colors()[c(124, 12, 82, 35)]
trellis.par.set(superpose.line=list(col=panelColors))
xyplot(ewm ~ year | comp.cluster, data=ewm, type='l', as.table=T, strip=myStripStyle, group=lake, scales=list(alternating=1, tck=c(1,0)), ylab='% FOO EWM', xlab='Year', main='Complete Linkage Clustering on DTW Distance')





##########################################################################
# trying some k-means partitioning
##########################################################################
ewm.kmeans.three <- kmeans(ewm.impute, centers=3, nstart=100)
km.clust.three <- ewm.kmeans.three$cluster
for (i in 1:nrow(ewm)) {
	ewm$kmeans.cluster[i] <- ewm.kmeans$cluster[names(ewm.kmeans$cluster)==ewm$lake[i]]
}
ewm.gg.kmeans <- ggplot(data=ewm, aes(x=year, y=ewm))
ewm.gg.kmeans + geom_line(aes(color=lake)) + facet_wrap(~kmeans.cluster) + theme_bw() + theme(legend.position='none') + labs(title='k-means clusters - k=4')


ewm.km.cascade <- cascadeKM(ewm.impute, inf.gr=2, sup.gr=10, iter=100, criterion='ssi')
plot(ewm.km.cascade, sortg=T)

# cascadeKM gave a k of 6 as the highest ssi; plotting output to assess visually
ewm.km.six <- kmeans(ewm.impute, centers=6, nstart=100)
km.clust.six <- ewm.km.six$cluster
for (i in 1:nrow(ewm)) {
	ewm$six.kmeans.cluster[i] <- ewm.new.km$cluster[names(ewm.new.km$cluster)==ewm$lake[i]]
}
ewm.gg.kmeans <- ggplot(data=ewm, aes(x=year, y=ewm))
ewm.gg.kmeans + geom_line(aes(color=lake)) + facet_wrap(~six.kmeans.cluster) + theme_bw() + theme(legend.position='none') + labs(title='k-means clusters - k=6')


############################################################################# here we're going to put these lakes to an ordination and see how environmental variables line up with clusters of ewm across time
###########################################################################
ewm.pcoa <- cmdscale(ewm.dist, k=2, eig=T)

# update water quality variables data.frame
wq.vars <- water.quality[water.quality$lake %in% managed$lake,]
rownames(wq.vars) <- wq.vars$lake; wq.vars <- subset(wq.vars, select= -c(lake, mn.chlorophyll, mn.totalP, mn.temp, mn.lake.level))
wq.vars$resDens <- wq.vars$bldgs / wq.vars$shoreline.km
wq.vars <- subset(wq.vars, select= -bldgs)
names(wq.vars) <- c('sa', 'shore', 'md', 'secchi', 'ph', 'orp', 'turb', 'do', 'tds', 'cond', 'sand', 'rock', 'muck', 'mdc', 'resDens')

## overlay environmental variables on top to look for correlations with our treatment groups
ewm.envfit <- envfit(ord=ewm.pcoa$points, env=wq.vars, na.rm=T)
colores <- colors()[c(12, 82, 124, 35)]
plot(ewm.pcoa$points, col=colores[cutg], pch=cutg+16, xlab='Dim1', ylab='Dim2', cex=2)
abline(h=0, v=0, lty=2, col='grey')
for (i in 1:4) {
	ordihull(ewm.pcoa$points[cutg==i,], groups=cutg[cutg==i], display='sites', draw=c('lines'), col=colores[i], lwd=1.5)
}
plot(ewm.envfit, p=0.054, cex=0.001, lwd=2)
text(locator(5), c('ResDens', 'TDS', 'Cond', 'SA', 'ORP'), cex=1.2, col='blue')



