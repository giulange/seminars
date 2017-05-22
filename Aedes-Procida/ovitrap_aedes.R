#_________
#
# D A T A
#_________

# R E A D 
D1 		<- read.table("~/work/People/G.Saccone/OVI_w1.txt",header = TRUE, sep = ",",na.strings = "NaN")
D2 		<- read.table("~/work/People/G.Saccone/OVI_w2.txt",header = TRUE, sep = ",",na.strings = "NaN")
D3 		<- read.table("~/work/People/G.Saccone/OVI_w3.txt",header = TRUE, sep = ",",na.strings = "NaN")
D4 		<- read.table("~/work/People/G.Saccone/OVI_w4.txt",header = TRUE, sep = ",",na.strings = "NaN")
summary(D1)
summary(D2)
summary(D3)
summary(D4)

#_________________
#
# P A C K A G E S
#_________________

require("methods")
library("sp")
library("gstat")
library("geoR")
library("raster")
require("ggplot2")
require("gridExtra")

#_________
#
# P L O T
#_________

# ___advanced-style___
plt1 <- qplot(Easting, Northing, data=D1, colour=week, size=week) + coord_fixed(ratio = 1) + theme_bw() + ggtitle("No. Eggs [week1]") + theme(plot.title = element_text(lineheight=.8, face="bold"))
plt2 <- qplot(Easting, Northing, data=D2, colour=week, size=week) + coord_fixed(ratio = 1) + theme_bw() + ggtitle("No. Eggs [week2]") + theme(plot.title = element_text(lineheight=.8, face="bold"))
plt3 <- qplot(Easting, Northing, data=D3, colour=week, size=week) + coord_fixed(ratio = 1) + theme_bw() + ggtitle("No. Eggs [week3]") + theme(plot.title = element_text(lineheight=.8, face="bold"))
plt4 <- qplot(Easting, Northing, data=D4, colour=week, size=week) + coord_fixed(ratio = 1) + theme_bw() + ggtitle("No. Eggs [week4]") + theme(plot.title = element_text(lineheight=.8, face="bold"))
grid.arrange(plt1,plt2,plt3,plt4,ncol=2,nrow=2)

#___________
#
# G R I D S
#___________

# data
coordinates(D1) = ~Easting+Northing
coordinates(D2) = ~Easting+Northing
coordinates(D3) = ~Easting+Northing
coordinates(D4) = ~Easting+Northing

# grid
#coordinates(grd) = ~Easting+Northing
#gridded(grd) = TRUE

#_______________________
#
# V A R I O G R A P H Y
#_______________________

#   -omnidirectional
#			>norm
vgm.exp.raw.D1 	= variogram( week~1, D1 )
vgm.exp.raw.D2 	= variogram( week~1, D2 )
vgm.exp.raw.D3 	= variogram( week~1, D3 )
vgm.exp.raw.D4 	= variogram( week~1, D4 )
plt1 <- ggplot(data=vgm.exp.raw.D1,aes(x=dist,y=gamma)) + geom_point(color='red',size=0.5) + ggtitle("D1 - Experimental Variogram [raw]") + theme(plot.title = element_text(lineheight=.8, face="bold"))
plt2 <- ggplot(data=vgm.exp.raw.D2,aes(x=dist,y=gamma)) + geom_point(color='red',size=0.5) + ggtitle("D2 - Experimental Variogram [raw]") + theme(plot.title = element_text(lineheight=.8, face="bold"))
plt3 <- ggplot(data=vgm.exp.raw.D3,aes(x=dist,y=gamma)) + geom_point(color='red',size=0.5) + ggtitle("D3 - Experimental Variogram [raw]") + theme(plot.title = element_text(lineheight=.8, face="bold"))
plt4 <- ggplot(data=vgm.exp.raw.D4,aes(x=dist,y=gamma)) + geom_point(color='red',size=0.5) + ggtitle("D4 - Experimental Variogram [raw]") + theme(plot.title = element_text(lineheight=.8, face="bold"))
grid.arrange(plt1,plt2,plt3,plt4,ncol=2,nrow=2)

#___________________________________
#
# remove low(Dist) high(Var) points
#___________________________________

# remove couples with high covariance but short distance apart:
#   <D1>
vgm.exp.raw.D1c_rd = variogram( week~1, D1, cloud=TRUE)
plot(vgm.exp.raw.D1c_rd)
short_dist 	<- vgm.exp.raw.D1c_rd[which(vgm.exp.raw.D1c_rd$dist<150.0),]
high_gamma 	<- short_dist[which(short_dist$gamma>20000),]
high_gamma 	<- data.frame(high_gamma)
list_remove <- sort(unique( c(high_gamma$left,high_gamma$right) ))
rem.D1 <- length(list_remove)
if(length(list_remove)==0){
	D1_rd <- D1
}else{
	D1_rd <- D1[-list_remove,]
}
#   <D2>
vgm.exp.raw.D2c_rd = variogram( week~1, D2, cloud=TRUE)
plot(vgm.exp.raw.D2c_rd)
short_dist 	<- vgm.exp.raw.D2c_rd[which(vgm.exp.raw.D2c_rd$dist<150.0),]
high_gamma 	<- short_dist[which(short_dist$gamma>20000),]
high_gamma 	<- data.frame(high_gamma)
list_remove <- sort(unique( c(high_gamma$left,high_gamma$right) ))
rem.D2 <- length(list_remove)
if(length(list_remove)==0){
	D2_rd <- D2
}else{
	D2_rd <- D2[-list_remove,]
}
#   <D3>
vgm.exp.raw.D3c_rd = variogram( week~1, D3, cloud=TRUE)
plot(vgm.exp.raw.D3c_rd)
short_dist 	<- vgm.exp.raw.D3c_rd[which(vgm.exp.raw.D3c_rd$dist<200.0),]
high_gamma 	<- short_dist[which(short_dist$gamma>35000),]
high_gamma 	<- data.frame(high_gamma)
list_remove <- sort(unique( c(high_gamma$left,high_gamma$right) ))
rem.D3 <- length(list_remove)
if(length(list_remove)==0){
	D3_rd <- D3
}else{
	D3_rd <- D3[-list_remove,]
}




# REMOVED [fit-vario-raw]:
# <D1>
vgm.exp.raw.D1_rd = variogram( week~1, D1_rd)# vgm( psill, model, range, nugget )
plot(vgm.exp.raw.D1_rd)
psill 	= 1600#max(vgm.exp.raw.D1_rd$gamma) - min(vgm.exp.raw.D1_rd$gamma)
range 	= 500#max(vgm.exp.raw.D1_rd$dist)/6
nugget 	= 700#vgm.exp.raw.D1_rd$gamma[1]/2
vgm.fit.raw.D1_rd = fit.variogram( vgm.exp.raw.D1_rd, model=vgm(psill,'Exp',range,nugget), fit.sills = TRUE, fit.ranges = TRUE, fit.method = 7, debug.level = 1, warn.if.neg = TRUE )
plot(vgm.exp.raw.D1_rd,vgm.fit.raw.D1_rd)
# <D2>
vgm.exp.raw.D2_rd = variogram( week~1, D2_rd)# vgm( psill, model, range, nugget )
plot(vgm.exp.raw.D2_rd)
psill 	= 2000#max(vgm.exp.raw.D2_rd$gamma) - min(vgm.exp.raw.L1_rd$gamma)
range 	= 150#max(vgm.exp.raw.D2_rd$dist)/6
nugget 	= 900#vgm.exp.raw.D2_rd$gamma[1]/2
vgm.fit.raw.D2_rd = fit.variogram( vgm.exp.raw.D2_rd, model=vgm(psill,'Exp',range,nugget), fit.sills = FALSE, fit.ranges = FALSE, fit.method = 7, debug.level = 1, warn.if.neg = TRUE )
#vgm.fit.raw.D2_rd = fit.variogram( vgm.exp.raw.D2_rd, model=vgm(psill,'Gau',range,nugget), fit.sills = FALSE, fit.ranges = FALSE, fit.method = 7, debug.level = 1, warn.if.neg = TRUE )
#vgm.fit.raw.D2_rd = fit.variogram( vgm.exp.raw.D2_rd, model=vgm(psill,'Sph',range,nugget), fit.sills = FALSE, fit.ranges = FALSE, fit.method = 7, debug.level = 1, warn.if.neg = TRUE )
plot(vgm.exp.raw.D2_rd,vgm.fit.raw.D2_rd)
# <D3>
vgm.exp.raw.D3_rd = variogram( week~1, D3_rd)# vgm( psill, model, range, nugget )
plot(vgm.exp.raw.D3_rd)
psill 	= 6300#max(vgm.exp.raw.D2_rd$gamma) - min(vgm.exp.raw.L1_rd$gamma)
range 	= 150#max(vgm.exp.raw.D2_rd$dist)/6
nugget 	= 200#vgm.exp.raw.D2_rd$gamma[1]/2
vgm.fit.raw.D3_rd = fit.variogram( vgm.exp.raw.D3_rd, model=vgm(psill,'Exp',range,nugget), fit.sills = FALSE, fit.ranges = FALSE, fit.method = 7, debug.level = 1, warn.if.neg = TRUE )
#vgm.fit.raw.D2_rd = fit.variogram( vgm.exp.raw.D2_rd, model=vgm(psill,'Gau',range,nugget), fit.sills = FALSE, fit.ranges = FALSE, fit.method = 7, debug.level = 1, warn.if.neg = TRUE )
#vgm.fit.raw.D2_rd = fit.variogram( vgm.exp.raw.D2_rd, model=vgm(psill,'Sph',range,nugget), fit.sills = FALSE, fit.ranges = FALSE, fit.method = 7, debug.level = 1, warn.if.neg = TRUE )
plot(vgm.exp.raw.D3_rd,vgm.fit.raw.D3_rd)

#...