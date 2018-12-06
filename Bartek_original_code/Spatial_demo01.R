#### keywords: spatial statistics, SPARC, K-function, data courtesy of Nikki Pelot
#### Hackathon 2018

########################################################################################
#### Some useful libraries that are used later
########################################################################################

#require(caret) for machine learning 
#require(tibble) for interacting with complex tables and data frames

require(magrittr)
require(dplyr)

#require(multcomp) for ANOVA
#require(emmeans)  for ANOVA

require(spatstat) ### for spatial statistics
require(png)  ### for loading the boundary maps in PNG format
library(httr) #### for data download

library(lctools) ### for spatial autocorrelation 
library(dixon)


########################################################################################
#### Load the data describing the ROIs, and the mask defining the contours
########################################################################################

### spatial.data <- read.csv("/Users/rajwa/Dropbox/Projects/Pelot/Results-R3.csv", header=T)
### nerve.win <- readPNG("/Users/rajwa/Dropbox/Projects/Pelot/Mask 1.png")

url1 <- "https://www.dropbox.com/s/xlm5xss0cienzer/Results-R3.csv?dl=1"
GET(url1, write_disk(tf <- tempfile(fileext = ".csv")))
spatial.data <- read.csv(tf, header=T) 


url2 <- "https://www.dropbox.com/s/qgiukoom1g65yn2/Mask_1.png?dl=1"
GET(url2, write_disk(tf <- tempfile(fileext = ".png")))
nerve.win <- readPNG(tf)

nerve.win.mask <- matrix(ncol=ncol(nerve.win), nrow=nrow(nerve.win), data=!as.logical(nerve.win), byrow=F)
nerve.win.mask <- apply(t(nerve.win.mask),1,rev) 
nerve.owin <- owin(mask=nerve.win.mask)
nerve.poly <- as.polygonal(nerve.owin)

#### preview the boundary
plot(nerve.poly, box=T)

spatial.data.xy <- data.frame(X=round(spatial.data[,"XM"]), Y=round(432-spatial.data[,"YM"]))
spat.ppp <- ppp(x=spatial.data.xy[,1], y=spatial.data.xy[,2], marks=(spatial.data$GrMedian)/1000, window=nerve.poly, checkdup=F)

#### preview the points (accepted and rejected)
plot(spat.ppp)


#### preview the CIELAB b-component of the ROIs
par(mar=c(5,5,0,5))
hist(spatial.data$GrMedian, breaks=100, cex.lab=2, xlab="Median CIELAB b* component (16 bit)", lwd=2, col="grey", main="")
abline(v=33000, col="red", lwd=3)


#### define the cut-off to separate two staining hues (set on 33000)
spatial.data$B <- cut(spatial.data$GrMedian, breaks=c(12890, 33000, 49660 ))
levels(spatial.data$B) <- c("Low", "High")

#### separate the point pattern using pre-defined discretization cut-off 
spat.ppp.cut <- cut.ppp(spat.ppp, breaks=c(12890, 33000, 49660)/1000, include.lowest=T)

par(mar=c(0,0,0,0))
plot(spat.ppp.cut, cols=c("black", "black", "black","blue", "orange"), border="blue", bg=c("darkgreen", "red", "red"), 
     how="perspective", pch=21, cex=1.5, main="", legend=T, lwd=3)


########################################################################################
#### make density maps for the high- and low-valued groups
########################################################################################

levels(spat.ppp.cut$marks) <- c("Low", "High")

#quartz(width = 6, height = 6)
par(mar=c(0,0,0,0))
plot(density(split(spat.ppp.cut)[[1]], sigma=20, fractional=T, diggle=T),ribbon=F, useRaster=T, main="",legend=F, col=grey.colors(100, start = 0, end = 1, gamma=1))
plot(density(split(spat.ppp.cut)[[2]], sigma=20, fractional=T, diggle=T),ribbon=F, useRaster=T, main="",legend=F, col=grey.colors(100, start = 0, end = 1, gamma=1))

#### Split the pattern, and analyze the resulting patterns separately
spat.ppp.1 <- split(spat.ppp.cut)[[1]]
spat.ppp.2 <- split(spat.ppp.cut)[[2]]

#### Simple estimate of the Ripley's K-function

K.1 <- Kest(spat.ppp.1, correction="translate")
K.2 <- Kest(spat.ppp.2, correction="translate")
par(mar=c(5,5,5,5))
plot(K.1, . ~ r,main="Ripley's K")
plot(K.1, cbind(sqrt(trans/pi), sqrt(theo/pi)) ~ r,main="Ripley's K")
plot(K.2, . ~ r,main="Ripley's K")
plot(K.2, cbind(sqrt(trans/pi), sqrt(theo/pi)) ~ r,main="Ripley's K")


K.1.env <- envelope(spat.ppp.1, fun=Kest, funargs=list(correction="best"))
K.2.env <- envelope(spat.ppp.2, fun=Kest, funargs=list(correction="best"))
plot(K.1.env, sqrt(./pi)-r ~ r, lwd=2, cex.lab=1.5, cex.axis=1.2, main="")
plot(K.2.env, sqrt(./pi)-r ~ r, lwd=2, cex.lab=1.5, cex.axis=1.2, main="")


#### Simple estimate of the J-function

J.1 <- Jest(spat.ppp.1, correction="best")
J.2 <- Jest(spat.ppp.2, correction="best")
par(mar=c(5,5,5,5))
plot(J.1, . - 1 ~ r,main="J function")
plot(J.2, . - 1 ~ r,main="J function")

J.1.env <- envelope(spat.ppp.1, fun=Jest, funargs=list(correction="best"))
J.2.env <- envelope(spat.ppp.2, fun=Jest, funargs=list(correction="best"))
plot(J.1.env, . -1 ~ r, lwd=2, cex.lab=1.5, cex.axis=1.2, main="")
plot(J.2.env, . -1 ~ r, lwd=2, cex.lab=1.5, cex.axis=1.2, main="")

#### Compute Hopkins-Skellam statistics

hopskel.test(split(spat.ppp.cut)[[2]], nsim=9990, method="MonteCarlo")


#### Compute Compute J-Foxall function

Jfox12.env <- envelope(Y=spat.ppp.1, fun=Jfox, funargs=list(spat.ppp.2), nsim=500 )
Jfox21.env <- envelope(Y=spat.ppp.2, fun=Jfox, funargs=list(spat.ppp.1), nsim=500 )
plot(Jfox12.env, . - 1 ~ r)
plot(Jfox21.env, . - 1 ~ r)


########################################################################################
#### Compute spatial autocorrelation
########################################################################################


bw <- 6
mI <- moransI(spatial.data.xy, bw, spatial.data$GrMode)#, WType="Bi-square")
moran.table <- matrix(data=NA,nrow=1,ncol=6)
col.names <- c("Moran's I", "Expected I", "Z resampling", "P-value resampling", "Z randomization", "P-value randomization")

colnames(moran.table) <- col.names

moran.table[1,1] <- mI$Morans.I
moran.table[1,2] <- mI$Expected.I
moran.table[1,3] <- mI$z.resampling
moran.table[1,4] <- mI$p.value.resampling
moran.table[1,5] <- mI$z.randomization
moran.table[1,6] <- mI$p.value.randomization

print(moran.table)

#### local Moran
l.moran <- l.moransI(Coords,6,spatial.data$GrMedian)


########################################################################################
#### Mark correlation and Mark variogram functions
########################################################################################


mark.spat <- markcorr(spat.ppp)
mark.env <- envelope(spat.ppp, fun=markcorr)
plot(mark.spat, cbind(theo, trans) ~r, main="")
plot(mark.env, (.) ~r, main="Mark correlation function")


marvar.spat <- markvario(spat.ppp)
markvar.env <- envelope(spat.ppp, fun=markvario)
plot(markvar.env)


markcon.spat <- markconnect(spat.ppp.cut, "Low", "High")
plot(markcon.spat)
markcon.env <- envelope(spat.ppp.cut, fun=markconnect)
plot(markcon.env)

########################################################################################
#### Pair Correlation Functions
########################################################################################

pcf.spat <- pcf(spat.ppp.1)
plot(pcf.spat)
pcf.env <- envelope(spat.ppp.1, fun= pcf)
plot(pcf.env, asinh(.) ~ r)
pcf.spat <- pcf(spat.ppp.2)
plot(pcf.spat)
pcf.env <- envelope(spat.ppp.2, fun= pcf)
plot(pcf.env, asinh(.) ~ r)



pcf.spat <- pcfcross(spat.ppp.cut, "Low", "High")
plot(pcf.spat)
pcf.env <- envelope(spat.ppp.cut, fun= pcfcross)
plot(pcf.env, sqrt(.) ~ r)

########################################################################################
#### Diggle-Cressie-Loosmore-Ford and Maximum Absolute Deviation Tests
########################################################################################


mad.test(spat.ppp.1, Fest, nsim=499, rmax=2, use.theo=T)
mad.test(spat.ppp.2, Fest, nsim=499, rmax=2, use.theo=T)
dclf.test(spat.ppp.1, Kest, nsim=499, use.theo=T)
dclf.test(spat.ppp.2, Kest, nsim=499, use.theo=T)


########################################################################################
#### Fit a statistical model
########################################################################################

spat.fit.1 <- ppm(spat.ppp.cut ~ 1,  Strauss(r=10), rbord=10)
spat.fit.2 <- ppm(spat.ppp.cut ~ polynom(x,y,2),  Strauss(r=10), rbord=10)
spat.fit.3 <- ppm(spat.ppp.cut ~ marks*polynom(x,y,2),  Strauss(r=10), rbord=10)

plot(spat.fit.2)

anova(spat.fit.0, spat.fit.1, spat.fit.2, spat.fit.3, test="LR")
sapply(list(spat.fit.0,spat.fit.1,spat.fit.2, spat.fit.3), AIC)


#spat.fit.1 <- kppm(spat.ppp.1 ~ polynom(x,y,2), "Thomas", method="palm")

sim.ppp <- simulate(spat.fit.3, nsim=4)

plot(sim.ppp, cols=c("black", "black", "black","blue", "orange"), border="blue", bg=c("darkgreen", "red", "red"), 
     how="perspective", pch=21, cex=0.5, main="", legend=F, lwd=1)



