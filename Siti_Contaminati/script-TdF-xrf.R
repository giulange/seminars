#!/usr/bin/Rscript

#_____________
#
# I N P U T S
#_____________
#
# inputs are given in this order:
#		(1) chemical element string (e.g. Chrome)
#		(2) number of geostatistical conditional simulations
# 			see this link as explanatory for Geostatistical Simulation:
#  	 			http://petrowiki.org/Geostatistical_conditional_simulation

#___________
#
# N O T E S
#___________
#
#		(1) combination of ggplot & gstat
#					https://rstudio-pubs-static.s3.amazonaws.com/46259_d328295794034414944deea60552a942.html
#		(2)	ggplot with anova by means of ad-hoc on-the-fly functions
#					http://www.cookbook-r.com/Graphs/Plotting_means_and_error_bars_(ggplot2)/

#___________________
#
# A R G U M E N T S
#___________________

# to enable arguments
args <- commandArgs(trailingOnly=TRUE)

# test passed arguments:
if(length(args)==0) {
	stop("At least one argument must be supplied (chemical element.n",call.=FALSE)
}else if (length(args)==1) {
	ChElstr = args[1]
	gcs_nsim = 2
}else {
	ChElstr = args[1]
	gcs_nsim = as.numeric( args[2] )
}

#ChElstr
print( paste("Species  :  ", ChElstr),  quote=FALSE )
print( paste("GCS nsim :  ", gcs_nsim), quote=FALSE )

# WDIR
setwd("/home/giuliano/work/Projects/terrafuochi/")

#_________________
#
# P A C K A G E S
#_________________

library("sp")
library("gstat")
library("geoR")
library("raster")
require("ggplot2")
require("gridExtra")

#___________________
#
# F U N C T I O N S
#___________________

## Gives count, mean, standard deviation, standard error of the mean, and confidence interval (default 95%).
##   data: a data frame.
##   measurevar: the name of a column that contains the variable to be summariezed
##   groupvars: a vector containing names of columns that contain grouping variables
##   na.rm: a boolean that indicates whether to ignore NA's
##   conf.interval: the percent range of the confidence interval (default is 95%)
summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
    library(plyr)

    # New version of length which can handle NA's: if na.rm==T, don't count them
    length2 <- function (x, na.rm=FALSE) {
        if (na.rm) sum(!is.na(x))
        else       length(x)
    }

    # This does the summary. For each group's data frame, return a vector with
    # N, mean, and sd
    datac <- ddply(data, groupvars, .drop=.drop,
      .fun = function(xx, col) {
        c(N    = length2(xx[[col]], na.rm=na.rm),
          mean = mean   (xx[[col]], na.rm=na.rm),
          sd   = sd     (xx[[col]], na.rm=na.rm)
        )
      },
      measurevar
    )

    # Rename the "mean" column    
    datac <- rename(datac, c("mean" = measurevar))

    datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean

    # Confidence interval multiplier for standard error
    # Calculate t-statistic for confidence interval: 
    # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
    ciMult <- qt(conf.interval/2 + .5, datac$N-1)
    datac$ci <- datac$se * ciMult

    return(datac)
}


#_________
#
# D A T A
#_________

# R E A D 
# xrf measurements
xrf 		<- read.table("/home/giuliano/work/Projects/terrafuochi/xrf_1d_ChEl.txt",header = TRUE, sep = ",",row.names=1,na.strings = "NaN")
# grid for interpolations
grd 		<- read.table("/home/giuliano/work/Projects/terrafuochi/xrf_grid.txt",header = TRUE, sep = ",")

# S Y N T H E S I S
tgc <- summarySE( xrf,measurevar="ChEl",groupvars=c("Layer","ChElthre") )

# rearrangement:
#xrf_2 	<- xrf
#idx.max <- which.max(xrf_2[,3])
#xrf_2 	<- xrf_2[-idx.max,]# delete row from data.frame
#summary(xrf_2)

# exploratory
#summary(xrf)
#length(xrf)

#_________
#
# P L O T
#_________

# RColorBrewer::display.brewer.all()

# ___basic-style___
#par(mfrow=c(1,2))
#plot( xrf, main = paste(ChElstr,"points location"))
#plot( xrf[which(xrf$ChElthre==1),], main = paste(ChElstr,":: locations exceeding Law-Threshold") )
#plot( xrf[which(xrf$ChElthre==0),], main = "below Law-Threshold" )
#plot( xrf[which(xrf$ChElthre==1),], col="red", lwd=1, main = "exceeding Law-Threshold" )

# ___advanced-style___
#ggplot(xrf[which(xrf$ChElthre==1),], aes(Easting,Northing)) + geom_point(aes(color=ChEl))
#ggplot(xrf[which(xrf$ChElthre==0),], aes(Easting,Northing)) + geom_point(aes(color=ChEl))
#qplot(Easting, Northing, data=xrf, colour=factor(ChElthre), size=ChEl)
plt1 <- qplot(Easting, Northing, data=xrf[which(xrf$Layer==1),], colour=ChEl, size=ChEl) + coord_fixed(ratio = 1) + theme_bw() + ggtitle("L1 [magnitude]") + theme(plot.title = element_text(lineheight=.8, face="bold"))
plt2 <- qplot(Easting, Northing, data=xrf[which(xrf$Layer==2),], colour=ChEl, size=ChEl) + coord_fixed(ratio = 1) + theme_bw() + ggtitle("L2 [magnitude]") + theme(plot.title = element_text(lineheight=.8, face="bold"))
plt3 <- qplot(Easting, Northing, data=xrf[which(xrf$Layer==3),], colour=ChEl, size=ChEl) + coord_fixed(ratio = 1) + theme_bw() + ggtitle("L3 [magnitude]") + theme(plot.title = element_text(lineheight=.8, face="bold"))
#plt1 <- qplot(Easting, Northing, data=xrf[which(xrf$Layer==1),], colour=ChEl) + coord_fixed(ratio = 1) + theme_bw() + ggtitle("L1 [magnitude]") + theme(plot.title = element_text(lineheight=.8, face="bold")) + scale_fill_brewer(palette="RdYlGn")
#plt2 <- qplot(Easting, Northing, data=xrf[which(xrf$Layer==2),], colour=ChEl) + coord_fixed(ratio = 1) + theme_bw() + ggtitle("L2 [magnitude]") + theme(plot.title = element_text(lineheight=.8, face="bold"))
#plt3 <- qplot(Easting, Northing, data=xrf[which(xrf$Layer==3),], colour=ChEl) + coord_fixed(ratio = 1) + theme_bw() + ggtitle("L3 [magnitude]") + theme(plot.title = element_text(lineheight=.8, face="bold"))
plt4 <- qplot(ChEl, data=xrf, geom="histogram",binwidth=25,main="All Soil Layers") + theme_bw() + xlab( paste(ChElstr,"Concentration [ppm]") )
grid.arrange(plt1,plt2,plt3,plt4,ncol=2,nrow=2)
#scale_fill_manual(values=c("#dd4050","#1aa130"))

plt1 <- qplot(Easting, Northing, data=xrf[which(xrf$Layer==1),], colour=factor(ChElthre)) + coord_fixed(ratio = 1) + theme_bw() + ggtitle("L1 [law-exceeding]") + theme(plot.title = element_text(lineheight=.8, face="bold"))
plt2 <- qplot(Easting, Northing, data=xrf[which(xrf$Layer==2),], colour=factor(ChElthre)) + coord_fixed(ratio = 1) + theme_bw() + ggtitle("L2 [law-exceeding]") + theme(plot.title = element_text(lineheight=.8, face="bold"))
plt3 <- qplot(Easting, Northing, data=xrf[which(xrf$Layer==3),], colour=factor(ChElthre)) + coord_fixed(ratio = 1) + theme_bw() + ggtitle("L3 [law-exceeding]") + theme(plot.title = element_text(lineheight=.8, face="bold"))
#plt4 <- qplot(aggregate(xrf$ChElthre, by=list(Category=xrf$Layer), FUN=sum),,geom="bar")
plt4 <- qplot(factor(Layer),data=xrf[which(xrf$ChElthre==1),],geom="bar",fill=factor(Layer)) + xlab("Soil Layers")
grid.arrange(plt1,plt2,plt3,plt4,ncol=2,nrow=2)
#grid.arrange(plt1,plt2,ncol=2,nrow=1)

# effect :: Magnitude ~ Uncertainty
ggplot(xrf, aes(x=1:length(xrf$ChEl), y=ChEl, colour=factor(ChElthre))) + 
    geom_errorbar(aes(ymin=ChEl-std, ymax=ChEl+std), width=.01) +
    xlab("Samples") +
    ylab( paste(ChElstr,"Concentration [ppm]") ) +
    scale_colour_hue(name="Law Threshold",    	# Legend label, use darker colors
                     breaks=c(0, 1),
                     labels=c("Below", "Above"),
                     l=40) +                    # Use darker colors, lightness=40
    ggtitle("The Effect of Magnitude \n on Chemical Concentration Uncertainty") +
    expand_limits(y=0) +                        # Expand y range
    #scale_y_continuous(breaks=0:20*4) +        # Set tick every 4
    theme_bw() +
    theme(legend.justification=c(1,0)
          #,legend.position=c(.5,.8)             # Position legend in bottom right
          )

# bar :: aggregating previous response
ggplot(tgc, aes(x=Layer, y=ChEl, fill=factor(ChElthre))) + 
    geom_bar(position=position_dodge(), stat="identity",
             colour="black", # Use black outlines,
             size=.3) +      # Thinner lines
    geom_errorbar(aes(ymin=ChEl-se, ymax=ChEl+se),
                  size=.3,    # Thinner lines
                  width=.2,
                  position=position_dodge(.9)) +
    xlab("Soil Layers") +
    ylab( paste(ChElstr,"Concentration [ppm]")) +
    scale_fill_hue(name="Law Threshold", # Legend label, use darker colors
                   breaks=c(0, 1),
                   labels=c("Below", "Above")) +
    ggtitle("The Effect of Magnitude and Soil Layer\non Uncertainty of Chemical Concentration") +
    #scale_y_continuous(breaks=0:20*4) +
    theme_bw()


#ggplot(tgc, aes(x=Layer, y=ChEl, colour=ChElthre)) + 
#    geom_errorbar(aes(ymin=ChEl-se, ymax=ChEl+se), width=.1) +
#    geom_point()

#ggplot(tgc, aes(x=Layer, y=ChEl, colour=ChElthre)) + 
#    geom_errorbar(aes(ymin=ChEl-ci, ymax=ChEl+ci), width=.1) +
#    geom_point()


#___________
#
# G R I D S
#___________

# data
coordinates(xrf) = ~Easting+Northing
#coordinates(xrf_2) = ~Easting+Northing
# grid
coordinates(grd) = ~Easting+Northing
gridded(grd) = TRUE
summary(xrf)

#_______________________
#
# V A R I O G R A P H Y
#_______________________

#   -omnidirectional
#			>norm
vgm.exp.n.L1 	= variogram( ChEl~1, xrf[which(xrf$Layer==1),])
vgm.exp.n.L2 	= variogram( ChEl~1, xrf[which(xrf$Layer==2),])
vgm.exp.n.L3 	= variogram( ChEl~1, xrf[which(xrf$Layer==3),])
plt1 <- ggplot(data=vgm.exp.n.L1,aes(x=dist,y=gamma)) + geom_point(color='red',size=0.5) + ggtitle("L1 [exp-vario-raw]") + theme(plot.title = element_text(lineheight=.8, face="bold"))
plt2 <- ggplot(data=vgm.exp.n.L2,aes(x=dist,y=gamma)) + geom_point(color='red',size=0.5) + ggtitle("L2 [exp-vario-raw]") + theme(plot.title = element_text(lineheight=.8, face="bold"))
plt3 <- ggplot(data=vgm.exp.n.L3,aes(x=dist,y=gamma)) + geom_point(color='red',size=0.5) + ggtitle("L3 [exp-vario-raw]") + theme(plot.title = element_text(lineheight=.8, face="bold"))
grid.arrange(plt1,plt2,plt3,ncol=2,nrow=2)
#			>log
vgm.exp.log.L1 = variogram( log(ChEl)~1, xrf[which(xrf$Layer==1),])
vgm.exp.log.L2 = variogram( log(ChEl)~1, xrf[which(xrf$Layer==2),])
vgm.exp.log.L3 = variogram( log(ChEl)~1, xrf[which(xrf$Layer==3),])
plt1 <- ggplot(data=vgm.exp.log.L1,aes(x=dist,y=gamma)) + geom_point(color='red',size=0.5) + ggtitle("L1 [exp-vario-log]") + theme(plot.title = element_text(lineheight=.8, face="bold"))
plt2 <- ggplot(data=vgm.exp.log.L2,aes(x=dist,y=gamma)) + geom_point(color='red',size=0.5) + ggtitle("L2 [exp-vario-log]") + theme(plot.title = element_text(lineheight=.8, face="bold"))
plt3 <- ggplot(data=vgm.exp.log.L3,aes(x=dist,y=gamma)) + geom_point(color='red',size=0.5) + ggtitle("L3 [exp-vario-log]") + theme(plot.title = element_text(lineheight=.8, face="bold"))
grid.arrange(plt1,plt2,plt3,ncol=2,nrow=2)
#			>log & cloud
vgm.exp.log.L1c 	= variogram( ChEl~1, xrf[which(xrf$Layer==1),], cloud=TRUE)
vgm.exp.log.L2c 	= variogram( ChEl~1, xrf[which(xrf$Layer==2),], cloud=TRUE)
vgm.exp.log.L3c 	= variogram( ChEl~1, xrf[which(xrf$Layer==3),], cloud=TRUE)
plt1 <- ggplot(data=vgm.exp.log.L1c,aes(x=dist,y=gamma)) + geom_point(color='red',size=0.2) + ggtitle("L1 [exp-cloud-log]") + theme(plot.title = element_text(lineheight=.8, face="bold"))
plt2 <- ggplot(data=vgm.exp.log.L2c,aes(x=dist,y=gamma)) + geom_point(color='red',size=0.2) + ggtitle("L2 [exp-cloud-log]") + theme(plot.title = element_text(lineheight=.8, face="bold"))
plt3 <- ggplot(data=vgm.exp.log.L3c,aes(x=dist,y=gamma)) + geom_point(color='red',size=0.2) + ggtitle("L3 [exp-cloud-log]") + theme(plot.title = element_text(lineheight=.8, face="bold"))
grid.arrange(plt1,plt2,plt3,ncol=2,nrow=2)

#par(mfrow=c(2,1))
#plot(vgm.exp_n, 	main = paste("[OMNI-DIR] Experimental Variogram",ChElstr,"[raw data]") )
#plot(vgm.exp_log, main = paste("[OMNI-DIR] Experimental Variogram",ChElstr,"[log data]") )

#   -directional
#vgm.exp_d4cl 		= variogram(log(ChEl)~1, xrf,alpha=c(0,45,90,135),cloud=TRUE)
#plot(vgm.exp_d4cl, main = paste("[4-DIR] Experimental Cloud Variogram",ChElstr,"[log data]"))
#vgm.exp_d4 			= variogram(log(ChEl)~1, xrf,alpha=c(0,45,90,135))
#plot(vgm.exp_d4, main = paste("[4-DIR] Experimental Variogram",ChElstr,"[log data]"))
#vgm.exp_d4_ang	= variogram(log(ChEl)~1, xrf,alpha=c(0,45,90,135)+6)
#plot(vgm.exp_d4_ang, main = paste("[4-DIR,6^] Experimental Variogram",ChElstr,"[log data]"))
#vgm.exp_d2 = variogram(log(ChEl)~1, xrf,alpha=c(0,90))
#plot(vgm.exp_d2, main = paste("[2-DIR] Experimental Variogram",ChElstr,"[log data]"))

#vgm.fit_d4 = fit.variogram( vgm.exp_d4, model=vgm(1.0,'Exp',100,0.0) )
#plot(vgm.exp_d4,vgm.fit_d4, main = paste("[4-DIR] Fitted Variogram",ChElstr,"[log data]"))

#   -fitting																					vgm(psill, model, range, nugget, ...)
vgm.fit.log.L1 = fit.variogram( vgm.exp.log.L1, model=vgm(0.30,'Exp',100,0.1), fit.sills = TRUE, fit.ranges = TRUE, fit.method = 7, debug.level = 1, warn.if.neg = TRUE )
vgm.fit.log.L2 = fit.variogram( vgm.exp.log.L2, model=vgm(1.00,'Exp',100,0.1), fit.sills = TRUE, fit.ranges = TRUE, fit.method = 7, debug.level = 1, warn.if.neg = TRUE )
vgm.fit.log.L3 = fit.variogram( vgm.exp.log.L3, model=vgm(1.00,'Exp',100,0.1), fit.sills = TRUE, fit.ranges = TRUE, fit.method = 7, debug.level = 1, warn.if.neg = TRUE )

#plot(vgm.exp_log,vgm.fit_log, main = paste("[OMNI-DIR] Fitted Variogram",ChElstr,"[log data]"))
# L1
maxDist <- round(max(vgm.exp.log.L1$dist))
vgm.fit.plot <- data.frame(1:maxDist,g_exp <- vgm.fit.log.L1$psill[1] + vgm.fit.log.L1$psill[2]*( 1 - exp(-(1:maxDist/vgm.fit.log.L1$range[2])) ) )
names(vgm.fit.plot) <- c("dist","gamma")
plt1 <- ggplot(data=vgm.exp.log.L1,aes(x=dist,y=gamma)) + geom_point(color='blue',size=2) + geom_line(data=vgm.fit.plot, aes(x=dist,y=gamma),color='black') + theme_bw() + 
				ggtitle("L1 [fit-vario-log]") + theme(plot.title = element_text(lineheight=.8, face="bold"))
# L2
maxDist <- round(max(vgm.exp.log.L2$dist))
vgm.fit.plot <- data.frame(1:maxDist,g_exp <- vgm.fit.log.L2$psill[1] + vgm.fit.log.L2$psill[2]*( 1 - exp(-(1:maxDist/vgm.fit.log.L2$range[2])) ) )
names(vgm.fit.plot) <- c("dist","gamma")
plt2 <- ggplot(data=vgm.exp.log.L2,aes(x=dist,y=gamma)) + geom_point(color='blue',size=2) + geom_line(data=vgm.fit.plot, aes(x=dist,y=gamma)) + theme_bw() + 
				ggtitle("L2 [fit-vario-log]") + theme(plot.title = element_text(lineheight=.8, face="bold"))
# L3
maxDist <- round(max(vgm.exp.log.L3$dist))
vgm.fit.plot <- data.frame(1:maxDist,g_exp <- vgm.fit.log.L3$psill[1] + vgm.fit.log.L3$psill[2]*( 1 - exp(-(1:maxDist/vgm.fit.log.L3$range[2])) ) )
names(vgm.fit.plot) <- c("dist","gamma")
plt3 <- ggplot(data=vgm.exp.log.L3,aes(x=dist,y=gamma)) + geom_point(color='blue',size=2) + geom_line(data=vgm.fit.plot, aes(x=dist,y=gamma)) + theme_bw() + 
				ggtitle("L3 [fit-vario-log]") + theme(plot.title = element_text(lineheight=.8, face="bold"))
grid.arrange(plt1,plt2,plt3,ncol=2,nrow=2)

#___________________________________
#
# remove low(Dist) high(Var) points
#___________________________________

# remove couples with high covariance but short distance apart:
#   <L1>
xrf_L1 			<- data.frame(	xrf[which(xrf$Layer==1),] )
coordinates(xrf_L1) = ~Easting+Northing
vgm.exp.log.L1c_rd = variogram( ChEl~1, xrf_L1, cloud=TRUE)
short_dist 	<- vgm.exp.log.L1c[which(vgm.exp.log.L1c_rd$dist<5.0),]
high_gamma 	<- short_dist[which(short_dist$gamma>100000),]
high_gamma 	<- data.frame(high_gamma)
list_remove <- sort(unique( c(high_gamma$left,high_gamma$right) ))
rem.L1 <- length(list_remove)
if(length(list_remove)==0){
	xrf_L1_rd   <- xrf_L1
}else{
	xrf_L1_rd   <- xrf_L1[-list_remove,]
}

#   <L2>
xrf_L2 			<- data.frame(	xrf[which(xrf$Layer==2),] )
coordinates(xrf_L2) = ~Easting+Northing
vgm.exp.log.L2c_rd = variogram( ChEl~1, xrf_L2, cloud=TRUE)
short_dist 	<- vgm.exp.log.L2c[which(vgm.exp.log.L2c_rd$dist<5.0),]
high_gamma 	<- short_dist[which(short_dist$gamma>100000),]
high_gamma 	<- data.frame(high_gamma)
list_remove <- sort(unique( c(high_gamma$left,high_gamma$right) ))
rem.L2 <- length(list_remove)
if(length(list_remove)==0){
	xrf_L2_rd   <- xrf_L2
}else{
	xrf_L2_rd   <- xrf_L2[-list_remove,]
}

#   <L3>
xrf_L3 			<- data.frame(	xrf[which(xrf$Layer==3),] )
coordinates(xrf_L3) = ~Easting+Northing
vgm.exp.log.L3c_rd = variogram( ChEl~1, xrf_L3, cloud=TRUE)
short_dist 	<- vgm.exp.log.L3c[which(vgm.exp.log.L3c_rd$dist<5.0),]
high_gamma 	<- short_dist[which(short_dist$gamma>100000),]
high_gamma 	<- data.frame(high_gamma)
list_remove <- sort(unique( c(high_gamma$left,high_gamma$right) ))
xrf_L3_rd   <- xrf_L3[-list_remove,]
rem.L3 <- length(list_remove)
if(length(list_remove)==0){
	xrf_L3_rd   <- xrf_L3
}else{
	xrf_L3_rd   <- xrf_L3[-list_remove,]
}

# REMOVED [fit-vario-raw]:
#coordinates(xrf_L1_rd) = ~Easting+Northing
vgm.exp.raw.L1_rd = variogram( ChEl~1, xrf_L1_rd)# vgm( psill, model, range, nugget )
psill 	= max(vgm.exp.raw.L1_rd$gamma) - min(vgm.exp.raw.L1_rd$gamma)
range 	= max(vgm.exp.raw.L1_rd$dist)/6
nugget 	= vgm.exp.raw.L1_rd$gamma[1]/2
vgm.fit.raw.L1_rd = fit.variogram( vgm.exp.raw.L1_rd, model=vgm(psill,'Exp',range,nugget), fit.sills = TRUE, fit.ranges = TRUE, fit.method = 7, debug.level = 1, warn.if.neg = TRUE )
#coordinates(xrf_L2_rd) = ~Easting+Northing
vgm.exp.raw.L2_rd = variogram( ChEl~1, xrf_L2_rd)
psill 	= max(vgm.exp.raw.L2_rd$gamma) - min(vgm.exp.raw.L2_rd$gamma)
range 	= max(vgm.exp.raw.L2_rd$dist)/6
nugget 	= vgm.exp.raw.L2_rd$gamma[1]/2
vgm.fit.raw.L2_rd = fit.variogram( vgm.exp.raw.L2_rd, model=vgm(psill,'Exp',range,nugget), fit.sills = TRUE, fit.ranges = TRUE, fit.method = 7, debug.level = 1, warn.if.neg = TRUE )
#coordinates(xrf_L3_rd) = ~Easting+Northing
vgm.exp.raw.L3_rd = variogram( ChEl~1, xrf_L3_rd)
psill 	= max(vgm.exp.raw.L3_rd$gamma) - min(vgm.exp.raw.L3_rd$gamma)
range 	= max(vgm.exp.raw.L3_rd$dist)/6
nugget 	= vgm.exp.raw.L3_rd$gamma[1]
vgm.fit.raw.L3_rd = fit.variogram( vgm.exp.raw.L3_rd, model=vgm(psill,'Exp',range,nugget), fit.sills = TRUE, fit.ranges = TRUE, fit.method = 7, debug.level = 1, warn.if.neg = TRUE )
#plot(vgm.exp_log,vgm.fit_log, main = paste("[OMNI-DIR] Fitted Variogram",ChElstr,"[raw data]"))
# L1
maxDist <- round(max(vgm.exp.raw.L1_rd$dist))
vgm.fit.plot <- data.frame(1:maxDist,g_exp <- vgm.fit.raw.L1_rd$psill[1] + vgm.fit.raw.L1_rd$psill[2]*( 1 - exp(-(1:maxDist/vgm.fit.raw.L1_rd$range[2])) ) )
names(vgm.fit.plot) <- c("dist","gamma")
plt1 <- ggplot(data=vgm.exp.raw.L1_rd,aes(x=dist,y=gamma)) + geom_point(color='blue',size=2) + geom_line(data=vgm.fit.plot, aes(x=dist,y=gamma),color='black') + theme_bw() + 
				ggtitle(paste("L1 [fit-vario-raw] rem:",rem.L1)) + theme(plot.title = element_text(lineheight=.8, face="bold"))
# L2
maxDist <- round(max(vgm.exp.raw.L2_rd$dist))
vgm.fit.plot <- data.frame(1:maxDist,g_exp <- vgm.fit.raw.L2_rd$psill[1] + vgm.fit.raw.L2_rd$psill[2]*( 1 - exp(-(1:maxDist/vgm.fit.raw.L2_rd$range[2])) ) )
names(vgm.fit.plot) <- c("dist","gamma")
plt2 <- ggplot(data=vgm.exp.raw.L2_rd,aes(x=dist,y=gamma)) + geom_point(color='blue',size=2) + geom_line(data=vgm.fit.plot, aes(x=dist,y=gamma)) + theme_bw() + 
				ggtitle(paste("L2 [fit-vario-raw] rem:",rem.L2)) + theme(plot.title = element_text(lineheight=.8, face="bold"))
# L3
maxDist <- round(max(vgm.exp.raw.L3_rd$dist))
vgm.fit.plot <- data.frame(1:maxDist,g_exp <- vgm.fit.raw.L3_rd$psill[1] + vgm.fit.raw.L3_rd$psill[2]*( 1 - exp(-(1:maxDist/vgm.fit.raw.L3_rd$range[2])) ) )
names(vgm.fit.plot) <- c("dist","gamma")
plt3 <- ggplot(data=vgm.exp.raw.L3_rd,aes(x=dist,y=gamma)) + geom_point(color='blue',size=2) + geom_line(data=vgm.fit.plot, aes(x=dist,y=gamma)) + theme_bw() + 
				ggtitle(paste("L3 [fit-vario-raw] rem:",rem.L3)) + theme(plot.title = element_text(lineheight=.8, face="bold"))
grid.arrange(plt1,plt2,plt3,ncol=2,nrow=2)

# REMOVED [fit-vario-log]:
#coordinates(xrf_L1_rd) = ~Easting+Northing
vgm.exp.log.L1_rd = variogram( log(ChEl)~1, xrf_L1_rd)
vgm.fit.log.L1_rd = fit.variogram( vgm.exp.log.L1_rd, model=vgm(0.90,'Exp',100,0.1), fit.sills = TRUE, fit.ranges = TRUE, fit.method = 7, debug.level = 1, warn.if.neg = TRUE )
#coordinates(xrf_L2_rd) = ~Easting+Northing
vgm.exp.log.L2_rd = variogram( log(ChEl)~1, xrf_L2_rd)
vgm.fit.log.L2_rd = fit.variogram( vgm.exp.log.L2_rd, model=vgm(0.90,'Exp',100,0.1), fit.sills = TRUE, fit.ranges = TRUE, fit.method = 7, debug.level = 1, warn.if.neg = TRUE )
#coordinates(xrf_L3_rd) = ~Easting+Northing
vgm.exp.log.L3_rd = variogram( log(ChEl)~1, xrf_L3_rd)
vgm.fit.log.L3_rd = fit.variogram( vgm.exp.log.L3_rd, model=vgm(0.90,'Exp',100,0.1), fit.sills = TRUE, fit.ranges = TRUE, fit.method = 7, debug.level = 1, warn.if.neg = TRUE )
#plot(vgm.exp_log,vgm.fit_log, main = paste("[OMNI-DIR] Fitted Variogram",ChElstr,"[log data]"))
# L1
maxDist <- round(max(vgm.exp.log.L1_rd$dist))
vgm.fit.plot <- data.frame(1:maxDist,g_exp <- vgm.fit.log.L1_rd$psill[1] + vgm.fit.log.L1_rd$psill[2]*( 1 - exp(-(1:maxDist/vgm.fit.log.L1_rd$range[2])) ) )
names(vgm.fit.plot) <- c("dist","gamma")
plt1 <- ggplot(data=vgm.exp.log.L1_rd,aes(x=dist,y=gamma)) + geom_point(color='blue',size=2) + geom_line(data=vgm.fit.plot, aes(x=dist,y=gamma),color='black') + theme_bw() + 
				ggtitle(paste("L1 [fit-vario-log] rem:",rem.L1)) + theme(plot.title = element_text(lineheight=.8, face="bold"))
# L2
maxDist <- round(max(vgm.exp.log.L2_rd$dist))
vgm.fit.plot <- data.frame(1:maxDist,g_exp <- vgm.fit.log.L2_rd$psill[1] + vgm.fit.log.L2_rd$psill[2]*( 1 - exp(-(1:maxDist/vgm.fit.log.L2_rd$range[2])) ) )
names(vgm.fit.plot) <- c("dist","gamma")
plt2 <- ggplot(data=vgm.exp.log.L2_rd,aes(x=dist,y=gamma)) + geom_point(color='blue',size=2) + geom_line(data=vgm.fit.plot, aes(x=dist,y=gamma)) + theme_bw() + 
				ggtitle(paste("L2 [fit-vario-log] rem:",rem.L2)) + theme(plot.title = element_text(lineheight=.8, face="bold"))
# L3
maxDist <- round(max(vgm.exp.log.L3_rd$dist))
vgm.fit.plot <- data.frame(1:maxDist,g_exp <- vgm.fit.log.L3_rd$psill[1] + vgm.fit.log.L3_rd$psill[2]*( 1 - exp(-(1:maxDist/vgm.fit.log.L3_rd$range[2])) ) )
names(vgm.fit.plot) <- c("dist","gamma")
plt3 <- ggplot(data=vgm.exp.log.L3_rd,aes(x=dist,y=gamma)) + geom_point(color='blue',size=2) + geom_line(data=vgm.fit.plot, aes(x=dist,y=gamma)) + theme_bw() + 
				ggtitle(paste("L3 [fit-vario-log] rem:",rem.L3)) + theme(plot.title = element_text(lineheight=.8, face="bold"))
grid.arrange(plt1,plt2,plt3,ncol=2,nrow=2)

# L1
#   -directional cloud
vgm.exp.log.L1_d4cl 		= variogram(log(ChEl)~1, xrf_L1_rd,alpha=c(0,45,90,135),cloud=TRUE)
plot(vgm.exp.log.L1_d4cl, main = paste("L1 [fit-vario-log-anis:4D]",ChElstr))
#plt1 <- ggplot(data=vgm.exp.log.L1_d4cl,aes(x=dist,y=gamma)) + geom_point(color='red',size=0.2) + ggtitle("L1 [exp-cloud-log]") + theme(plot.title = element_text(lineheight=.8, face="bold"))
#plt2 <- ggplot(data=vgm.exp.log.L2c,aes(x=dist,y=gamma)) + geom_point(color='red',size=0.2) + ggtitle("L2 [exp-cloud-log]") + theme(plot.title = element_text(lineheight=.8, face="bold"))
#plt3 <- ggplot(data=vgm.exp.log.L3c,aes(x=dist,y=gamma)) + geom_point(color='red',size=0.2) + ggtitle("L3 [exp-cloud-log]") + theme(plot.title = element_text(lineheight=.8, face="bold"))
#plt4 <- ...
#grid.arrange(plt1,plt2,plt3,plt4,ncol=2,nrow=2)
#   -directional
vgm.exp.log.L1_d4 			= variogram(log(ChEl)~1, xrf_L1_rd,alpha=c(0,45,90,135))
plot(vgm.exp.log.L1_d4, main = paste("L1 [fit-vario-log-anis:4D",ChElstr))

# L2
#   -directional cloud
vgm.exp.log.L2_d4cl 		= variogram(log(ChEl)~1, xrf_L2_rd,alpha=c(0,45,90,135),cloud=TRUE)
plot(vgm.exp.log.L2_d4cl, main = paste("L2 [fit-vario-log-anis:4D]",ChElstr))
#   -directional
vgm.exp.log.L2_d4 			= variogram(log(ChEl)~1, xrf_L2_rd,alpha=c(0,45,90,135))
plot(vgm.exp.log.L2_d4, main = paste("L2 [fit-vario-log-anis:4D",ChElstr))

# L3
#   -directional cloud
vgm.exp.log.L3_d4cl 		= variogram(log(ChEl)~1, xrf_L3_rd,alpha=c(0,45,90,135),cloud=TRUE)
plot(vgm.exp.log.L3_d4cl, main = paste("L3 [fit-vario-log-anis:4D]",ChElstr))
#   -directional
vgm.exp.log.L3_d4 			= variogram(log(ChEl)~1, xrf_L3_rd,alpha=c(0,45,90,135))
plot(vgm.exp.log.L3_d4, main = paste("L3 [fit-vario-log-anis:4D",ChElstr))

# L1
#   -directional
vgm.exp.log.L1_d2 = variogram(log(ChEl)~1, xrf_L1_rd,alpha=c(0,90)+6)
plot(vgm.exp.log.L1_d2, main = paste("L1 [fit-vario-log-anis:4D",ChElstr))


#___________________________
#
# I N T E R P O L A T I O N
#___________________________

#   -ok,  global
#krg.ok = krige(log(ChEl)~1, xrf, grd, model = vgm.fit_log)
#spplot(krg.ok["var1.pred"],  main = paste(names(xrf),"ok global"))

#   -ok,  local
krg.ok.log.L1 = krige(log(ChEl)~1, xrf[which(xrf$Layer==1),], grd, model = vgm.fit.log.L1, nmax=20, nmin=5, maxdist=60)
krg.ok.log.L2 = krige(log(ChEl)~1, xrf[which(xrf$Layer==2),], grd, model = vgm.fit.log.L2, nmax=20, nmin=5, maxdist=60)
krg.ok.log.L3 = krige(log(ChEl)~1, xrf[which(xrf$Layer==3),], grd, model = vgm.fit.log.L3, nmax=20, nmin=5, maxdist=60)
plt1 = spplot(krg.ok.log.L1, "var1.pred", main="L1 [pred-ok-log]")
plt2 = spplot(krg.ok.log.L2, "var1.pred", main="L2 [pred-ok-log]")
plt3 = spplot(krg.ok.log.L3, "var1.pred", main="L3 [pred-ok-log]")
print(plt1, position = c(0,.5,.5,1),more=T)		# c(xmin, ymin, xmax, ymax)
print(plt2, position = c(.5,.5,1,1),more = T) # c(xmin, ymin, xmax, ymax)
print(plt3, position = c(0,0,.5,.5))						# c(xmin, ymin, xmax, ymax)

krg.ok.log.L1_rd = krige(log(ChEl)~1, xrf_L1_rd, grd, model = vgm.fit.log.L1_rd, nmax=20, nmin=5, maxdist=60)
krg.ok.log.L2_rd = krige(log(ChEl)~1, xrf_L2_rd, grd, model = vgm.fit.log.L2_rd, nmax=20, nmin=5, maxdist=60)
krg.ok.log.L3_rd = krige(log(ChEl)~1, xrf_L3_rd, grd, model = vgm.fit.log.L3_rd, nmax=20, nmin=5, maxdist=60)
#spplot(krg.okp["var1.pred"], main = paste(ChElstr,"ok local kriging" ))
#krg_plt <- as.data.frame(krg.ok.log.L1)
#plt1<-ggplot(data=krg_plt,aes(x=Easting,y=Northing)) + c(geom_tile(data=krg_plt,aes(fill=var1.pred))) + scale_fill_gradient(low="#FEEBE2", high="#7A0177")+coord_equal()
plt1 = spplot(krg.ok.log.L1_rd, "var1.pred", main="L1 [pred-ok-log] adj")
plt2 = spplot(krg.ok.log.L2_rd, "var1.pred", main="L2 [pred-ok-log] adj")
plt3 = spplot(krg.ok.log.L3_rd, "var1.pred", main="L3 [pred-ok-log] adj")
print(plt1, position = c(0,.5,.5,1),more=T)		# c(xmin, ymin, xmax, ymax)
print(plt2, position = c(.5,.5,1,1),more = T) # c(xmin, ymin, xmax, ymax)
print(plt3, position = c(0,0,.5,.5))						# c(xmin, ymin, xmax, ymax)

plt1 = spplot(krg.ok.log.L1_rd, "var1.var", main="L1 [var-ok-log] adj")
plt2 = spplot(krg.ok.log.L2_rd, "var1.var", main="L2 [var-ok-log] adj")
plt3 = spplot(krg.ok.log.L3_rd, "var1.var", main="L3 [var-ok-log] adj")
print(plt1, position = c(0,.5,.5,1),more = T)		# c(xmin, ymin, xmax, ymax)
print(plt2, position = c(.5,.5,1,1),more = T) # c(xmin, ymin, xmax, ymax)
print(plt3, position = c(0,0,.5,.5))						# c(xmin, ymin, xmax, ymax)


#plt 		<- ggplot(xrf, aes(Easting,Northing)) + geom_point(aes(color=ChEl))#start with the base-plot 
#layer1 	<- c(geom_tile(data=idw.output,aes(fill=var1.pred)))#then create a tile layer and fill with predicted values
#layer2	<- c(geom_path(data=boroughoutline,aes(long, lat, group=group),colour = "grey40", size=1))#then create an outline layer
## now add all of the data together
#plot+layer1+layer2+scale_fill_gradient(low="#FEEBE2", high="#7A0177")+coord_equal()

#   -gls, local
#     WHY?
#      (1) to capture heterogeneity
#      (2) to quantify and assess uncertainty
#     HOW?
#      (a) variogram
#      (b) cumulative distribution function, cdf
if(gcs_nsim==0) {
    print("Conditional Gaussian simulation is skipped!")
}else{
    #krg.gcs.L1 = krige(log(ChEl)~1, xrf[which(xrf$Layer==1),], grd, model = vgm.fit.log.L1, nmax=20, nmin=5, maxdist=60, nsim=gcs_nsim)
    #krg.gcs.L2 = krige(log(ChEl)~1, xrf[which(xrf$Layer==2),], grd, model = vgm.fit.log.L2, nmax=20, nmin=5, maxdist=60, nsim=gcs_nsim)
    #krg.gcs.L3 = krige(log(ChEl)~1, xrf[which(xrf$Layer==3),], grd, model = vgm.fit.log.L3, nmax=20, nmin=5, maxdist=60, nsim=gcs_nsim)

    # rd - raw
    krg.gcs.raw.L1_rd = krige(ChEl~1, xrf_L1_rd, grd, model = vgm.fit.raw.L1_rd, nmax=20, nmin=5, maxdist=60, nsim=gcs_nsim)
    krg.gcs.raw.L2_rd = krige(ChEl~1, xrf_L2_rd, grd, model = vgm.fit.raw.L2_rd, nmax=20, nmin=5, maxdist=60, nsim=gcs_nsim)
    krg.gcs.raw.L3_rd = krige(ChEl~1, xrf_L3_rd, grd, model = vgm.fit.raw.L3_rd, nmax=20, nmin=5, maxdist=60, nsim=gcs_nsim)
    spplot(krg.gcs.raw.L1_rd, main = "L1 [pred-gcs-raw] adj" )
    spplot(krg.gcs.raw.L2_rd, main = "L2 [pred-gcs-raw] adj" )
    spplot(krg.gcs.raw.L3_rd, main = "L3 [pred-gcs-raw] adj" )

    # rd - log
    krg.gcs.log.L1_rd = krige(log(ChEl)~1, xrf_L1_rd, grd, model = vgm.fit.log.L1_rd, nmax=20, nmin=5, maxdist=60, nsim=gcs_nsim)
    krg.gcs.log.L2_rd = krige(log(ChEl)~1, xrf_L2_rd, grd, model = vgm.fit.log.L2_rd, nmax=20, nmin=5, maxdist=60, nsim=gcs_nsim)
    krg.gcs.log.L3_rd = krige(log(ChEl)~1, xrf_L3_rd, grd, model = vgm.fit.log.L3_rd, nmax=20, nmin=5, maxdist=60, nsim=gcs_nsim)
    spplot(krg.gcs.log.L1_rd, main = "L1 [pred-gcs-log] adj" )
    spplot(krg.gcs.log.L2_rd, main = "L2 [pred-gcs-log] adj" )
    spplot(krg.gcs.log.L3_rd, main = "L3 [pred-gcs-log] adj" )
}


#__________
#
# S A V E 
#__________

# experimental variogram:   paste("vgm_exp_log_L1_rd__",ChElstr,".txt",sep="")
write.table( vgm.exp.log.L1_rd, file=paste("vgm_exp_log_L1_rd__",ChElstr,".txt",sep=""),sep=",",na="NaN",col.names=TRUE,row.names=FALSE)
write.table( vgm.exp.log.L2_rd, file=paste("vgm_exp_log_L2_rd__",ChElstr,".txt",sep=""),sep=",",na="NaN",col.names=TRUE,row.names=FALSE)
write.table( vgm.exp.log.L3_rd, file=paste("vgm_exp_log_L3_rd__",ChElstr,".txt",sep=""),sep=",",na="NaN",col.names=TRUE,row.names=FALSE)
# synthetic variogram:
write.table( vgm.fit.log.L1_rd, file=paste("vgm_fit_log_L1_rd__",ChElstr,".txt",sep=""),sep=",",na="NaN",col.names=TRUE,row.names=FALSE)
write.table( vgm.fit.log.L2_rd, file=paste("vgm_fit_log_L2_rd__",ChElstr,".txt",sep=""),sep=",",na="NaN",col.names=TRUE,row.names=FALSE)
write.table( vgm.fit.log.L3_rd, file=paste("vgm_fit_log_L3_rd__",ChElstr,".txt",sep=""),sep=",",na="NaN",col.names=TRUE,row.names=FALSE)
# ok
write.table( krg.ok.log.L1_rd, file="L1_ok.txt",sep=",",na="NaN",col.names=TRUE,row.names=FALSE)
write.table( krg.ok.log.L2_rd, file="L2_ok.txt",sep=",",na="NaN",col.names=TRUE,row.names=FALSE)
write.table( krg.ok.log.L3_rd, file="L3_ok.txt",sep=",",na="NaN",col.names=TRUE,row.names=FALSE)
# gcs
#write.table( krg.gcs.raw.L1_rd, file="L1_gcs.txt",sep=",",na="NaN",col.names=TRUE,row.names=FALSE)
#write.table( krg.gcs.log.L1_rd, file="L1_gcs.txt",sep=",",na="NaN",col.names=TRUE,row.names=FALSE)