#Solitary Sandpiper z589

#Archival light-level geolocator workshop
# Michael T. Hallworth 
# Smithsonian Migratory Bird Center
# 
#Geolocation analysis with Open Source Tools 2016 North American Ornithological Congress, Washington D.C.
#Eli Bridge, Simeon Lisovski, Eldar Rakhimberdiev, and Michael Hallworth
#
# Set-up--------------------------------------------------------------------------------------------------------------
# Check to make sure the required packages are installed on your machine
# If not, they will be installed


reqPackages <- c("devtools","digest","GeoLight","geosphere","raster","fields","forecast",
                 "circular","truncnorm","parallel","bit","rgdal","CircStats","Rcpp",
                 "RcppArmadillo","ggmap","ggsn","sp","maptools","rgeos","MASS")


get.packages <- reqPackages[!(reqPackages %in% installed.packages()[,"Package"])]

if(length(get.packages)>0){
  install.packages(get.packages,repos = "https://cloud.r-project.org/", dependencies = TRUE)
}

# Load devtools library before install_github so it will work. 
library(devtools)

# Install necessary packages from Github using the devtools library #
install_github("SWotherspoon/SGAT")
install_github("SLisovski/TwGeos") 
install_github("SLisovski/GeoLight", ref = "Update_2.01", force = T)
install_github("eldarrak/FLightR")
install_github("SWotherspoon/BAStag")

# In order to use different tol values for spring and fall equinox periods you need to 
# get edits that MTHallworth made to the SGAT package - 
install_github("MTHallworth/SGAT")
install_github("MTHallworth/LLmig") # you may want to try the MigSchedule function to
                                    # determine the migration schedule. I'm also really
                                    # interested to see if the function will work for a
                                    # wide variety of species. So far so good! 


# Load required packages - note I changed where SGAT gets loaded #

library(TwGeos)
library(GeoLight)
library(TwGeos)
library(FLightR)
library(MASS)  #needed for fitting distributions
library(maptools)
library(MASS)
library(raster)
library(sp)
library(SGAT)
library(LLmig)

# read in a simple world map from the maptools package #
data(wrld_simpl)                             

setwd = ("T:/McDuffie/SOSA_Geo/SOSA_R/SOSA_589")

list.files(path = "T:/McDuffie/SOSA_GEO/SOSA_R/SOSA_589", pattern = NULL, all.files = FALSE,
           full.names = FALSE, recursive = FALSE,
           ignore.case = FALSE, include.dirs = FALSE, no.. = FALSE)

#For Intigeo geos: read the data into a dataframe called d.lux
d.lux <- readMTlux("T:/McDuffie/SOSA_GEO/SOSA_R/SOSA_589/Z589_13Jun17_021544.lux") 

# Save only the columns named Date and Light
d.lux<-d.lux[,c("Date","Light")]

# Show first few lines of the data 
head(d.lux)

#Twilight adjustments--------------------------------------------------------------------------------------------------

# Set the capture coordinates for each bird #
CapLocs <- c(-149.750179,  61.284619)
cap.lat <-61.28
cap.long <- -149.75

# Set the seed
# The time when we know the geolocator recorded no light.
# You don't need that many replicates of the dates - I removed the ones that weren't needed
seed <- as.POSIXct(c("2016-12-01 04:00:00", 
                     "2017-01-01 04:00:00", 
                     "2017-02-01 04:00:00"), format = "%Y-%m-%d %H:%M:%S", tz = "GMT")

twl  <- findTwilights(tagdata = d.lux, 
                      threshold = 1.5, 
                      include = seed,
                      dark.min = 0) # minimum dark period in minutes

##### Plot the light data and the twilights ######
lightImage(tagdata = d.lux,
           offset = 18)
tsimagePoints(twl$Twilight, 
              offset = 19, 
              pch = 16, cex = 0.5, 
              col = ifelse(twl$Rise, "blue", "red"))


#Edit the transitions
#Moves or deletes outliers
twl <- twilightEdit(twilights = twl, 
                    window = 4,           # two days before and two days after
                    outlier.mins = 45,    # difference in mins
                    stationary.mins = 25, # are the other surrounding twilights within 25 mins of one another
                    plot = TRUE)

#In order to take advantage of the edits - you can subset the data using when Deleted == FALSE this removes 
#any of the twilights which were deleted.
#Here we adjust the twilight time back to account for the max light level readings every # minutes
twl <- twilightAdjust(twilights = twl, 
                      interval = 300) # The unit here is seconds

twl <- subset(twl, Deleted == FALSE)

# Zenith Angle--------------------------------------------------------------------------------------------------

# Create a vector with the dates known to be at deployment

# These dates should be changed a bit - as of right now the whole year is used as calibration dates
# what you want are dates when you know that the bird is at the deployment site. This can be a bit 
# tricky for near arctic breeding birds when tags are deployed during the summer. I moved some things
# around. Here, I use the lightImage and the cal.data to determine when the bird arrives in the spring
# I use those data until incubation starts to calculate zenith angles. Just after deployment doesn't 
# work that well in this case because of the nearly constant light issues. 

# Create a vector of the dates the geolocator was recording data
tm <- seq(d.lux[1,1], d.lux[nrow(d.lux),1], by = "day")
head (tm)

# Create a vector of TRUE /FALSE to determine whether it's sunrise or sunset
rise <- rep(c(TRUE, FALSE), length(tm))

# making predicted twilight times given location and zenith #
cal.data<- data.frame(Twilight = twilight(rep(tm, each = 2),
                           lon = CapLocs[1], 
                           lat = CapLocs[2], 
                           rise = rise, zenith = 95),
                           Rise = rise) 

##### Plot the light data and sunrise suset times at capture location ######
lightImage(tagdata = d.lux,
           offset = 18,
           zlim = c(0,5))
tsimagePoints(cal.data$Twilight, 
              offset = 19, 
              pch = 16, cex = 0.5, 
              col = ifelse(twl$Rise, "red", "blue"))
tsimagePoints(twl$Twilight, 
              offset = 19, 
              pch = 16, cex = 0.5, 
              col = ifelse(twl$Rise, "blue", "red"))

#Using this next bit of code you can click on the lightImage plot and find dates 
# as.POSIXct(locator(n = 2)$x,origin = "1970-01-01",tz = "GMT")
# "2017-05-08 04:15:19 GMT" "2017-05-17 17:40:58 GMT"

# This plots a polygon around the area that I used to determine Zenith Angles #
polygon(x = as.POSIXct(c("2017-05-08","2017-05-17","2017-05-17","2017-05-08")),
        y = c(4+24,4+24,14+24,14+24), 
        col = rgb(100,100,100,100,maxColorValue = 255),
        border = rgb(100,100,100,100,maxColorValue = 255))


calib.dates <- as.POSIXct(c("2017-05-08","2017-05-17"), tz = "GMT")

# Extract twilight data during calibration period
cal.data<-subset(twl,twl$Twilight>=calib.dates[1] & twl$Twilight<=calib.dates[2])

# Determine the sun elevation angle - here called the Zenith angle #
# Calculate solar time from calibration data 
sun  <- solar(cal.data[,1])

# Adjust the solar zenith angle for atmospheric refraction
z<- refracted(zenith(sun = sun,
                          lon = CapLocs[1],
                          lat = CapLocs[2]))

twl_t <- twilight(cal.data[,1], CapLocs[1], CapLocs[2], rise = cal.data[,2], zenith = quantile(z,probs=0.5))

# Determine the difference in minutes from when the sun rose and the geolocator said it rose 
twl_dev <- ifelse(cal.data$Rise, as.numeric(difftime(cal.data[,1], twl_t, units = "mins")),
                  as.numeric(difftime(twl_t, cal.data[,1], units = "mins")))

# Throw out values less than 0 - These values are not valid 
twl_dev<-subset(twl_dev, subset=twl_dev>=0)

# Describe the distribution of the error 
fitml <- fitdistr(twl_dev, "log-Normal")
head (fitml)

# save the Twilight model parameters
alpha<- c(fitml$estimate[1], fitml$estimate[2]) 
head (alpha)

# Make some plots to visualize the data 
seq <- seq(0,60, length = 100)
hist(twl_dev, freq = F,
     yaxt="n",
     ylim = c(0, 0.15),
     xlim = c(0, 60),
     breaks=3,
     main = "",
     xlab = "Twilight error (mins)")
axis(2,las=2)
lines(seq, dlnorm(seq, alpha[1], alpha[2]), col = "firebrick", lwd = 3,lty = 2)
#look at plot to see fit of data matches red dotted line, if not adjust fitml = "log-normal" to either normal or gamma

# Store the zenith (sun-elevation angle)
zenith0 <-quantile(z,prob=0.5)
zenith1 <-quantile(z,prob=0.95)

#Zenith angle plot
par(bty="l")
plot(median(z,na.rm=TRUE),ylim = c(80,120),pch=19,cex=1.25,ylab="Zenith Angle")
segments(1,quantile(z,probs=0.025),1,quantile(z,probs=0.95))
points(1,max(z,na.rm=TRUE),col="red",pch=20)

# Initial Locations----------------------------------------------------------------------------------------------

# How do I set tolernace values around equinox for South America?  
   # see code below
# How do I remove erroneous migration paths such as those along the east coast of the United States during equinox?
   # see tol code below
# How do I remove the migration paths into the Gulf of Alaska? These may be due to light-level fluctuations during
# incubation in June 2017?
   # There are several ways to deal with this - 1) set fixed = TRUE for the month of May - June. 2) Not include those
   # dates because it's already back on the breeding grounds  where you caught it (presumably) and locations are 
   # compromised because of incubation. Here - I set the fixed = TRUE for those dates where the incubation seems to
   # compromise the location data. 


# Obtain the initial path based on a simple threshold analysis

# Up above - I changed the calibration dates May 08 2017 - near the end of the file. 
# I didn't change this bit of code - so there are only a few days until you capture it
# from May 08 2017. I removed the calib.dates subset. It should now use all the data
# until captured that aren't deleted. 

d.twl <-subset(twl,!Deleted)

path <- thresholdPath(twilight = d.twl$Twilight, 
                      rise = d.twl$Rise, 
                      zenith=zenith0,
                      tol= c(0.16,0.09)) # Here is where you can mess with the tol values
                                          # the first value is fall - the second is spring

# NOTE - setting the tol values is a bit trial and error - I settled one these two values because they seem 
# to fit the general trajectory of the migration. For example, the spring and fall tracks over the Andes appears
# to be roughly in the same place. Larger tol values in the fall have the bird going up through Florida (didn't)
# seem right given a few points were near/on the Yucatan. Also, one of the potential stop-overs appears to be 
# similar in the fall and spring on Ecuador/Peru border. Larger tol values smooth that out. Ultimately, these
# values need to be justified i.e., maximized the amount of time on land, matched similar route in the spring etc.
# However, in the end these values (lat,long) are the starting values for the MCMC chain but they could potentially
# be informative. 

layout(matrix(c(1,1,3,3,3,
                2,2,3,3,3),2,5,byrow = TRUE))
       
par(mar = c(2,4,1,1)+0.1)
plot(path$time, path$x[, 2], type = "b", pch = 16, cex = 0.5, ylab = "Lat", xlab = '')
abline(h = cap.lat)
abline(v = as.POSIXct("2016-09-23"),col="red",lty=2,lwd=1.5)
abline(v = as.POSIXct("2017-03-20"),col="red",lty=2,lwd=1.5)

plot(path$time, path$x[,1],type="b",pch=16,cex=0.5,ylab="Long",xlab="")
abline(h = cap.long)
abline(v = as.POSIXct("2016-09-23"),col="red",lty=2,lwd=1.5)
abline(v = as.POSIXct("2017-03-20"),col="red",lty=2,lwd=1.5)

plot(path$x, type = "n",ylim=c(-60,60),xlim=c(-180,-55))
plot(wrld_simpl, add = T, col = "grey95")
lines(path$x, col = "blue")
points(path$x, pch = 16, cex = 0.5, col = "blue")
points(CapLocs[1],CapLocs[2],pch = 19, col = "red", cex = 2)

#Take the location estimates created above
x0<- path$x

# the model also needs the mid-points - generate those here
z0<- trackMidpts(x0)

# Set flight parameters
beta <- c(0.7, 0.08)



# This function sets the known location of the bird - capture location and recapture location
# This is where you can add ancillary resighting events; known locations throughout annual cycle
# Here I added Fixed = TRUE to the end of the file  - from May 08 on.

# Find which row in the twl file is May 08 2017 # 
May08 <- max(which(strptime(twl$Twilight,format = "%Y-%m-%d") == "2017-05-08"))

fixedx<- rep(F, nrow(x0))
fixedx[1:5] <- TRUE
# Fixed from May - end of file #
fixedx[May08:nrow(x0)] <- TRUE

x0[fixedx, 1] <- cap.long
x0[fixedx, 2] <- cap.lat

# Set what xlim and ylim values you need to span the range of your dataset
xlim <- c(-180, -30)
ylim <- c(-80, 90)

# Spatial Mask (do not use for shorebirds)

#Estelle Model#########################################################################################################

#The threshold.model function requires the following:
# 1)the twilight times and whether they are sunrise or sunset
# 2)a modified Log Normal model for the twilight errors (with relaxed assumptions to tune the proposals)
# 3)the parameters of the distribution of twilight errors (alpha)
# 4)the parameters of the distribution of travel speed (beta)
# 5)optional: the spatial mask
# 6)the initial x and z, and
# 7)the zenith angle that defines twilight.

model<- thresholdModel(d.twl$Twilight,d.twl$Rise,
                       twilight.model="ModifiedLogNormal",
                       alpha=alpha,
                       beta=beta,
                       x0=x0,
                       z0=z0,
                       zenith=zenith0,
                       fixedx=fixedx)
# Here you need to set the first few locations as "known locations"-fixed locations. These are periods  
# when you know the bird was at a specific location - capture site - when first deployed and captured. 

model$fixedx<-c(model$fixedx,rep(FALSE,(dim(model$x0)[1]-length(model$fixedx))))

model$fixedx[(length(model$fixedx)-3):length(model$fixedx)]<-TRUE

# This defines the error distribution around each location 

proposal.x <- mvnorm(S=diag(c(0.0025,0.0025)),n=nlocation(x0))
proposal.z <- mvnorm(S=diag(c(0.0025,0.0025)),n=nlocation(z0))          

fit <- estelleMetropolis(model = model,
                         proposal.x = proposal.x,
                         proposal.z = proposal.z,
                         iters=2000, # This value sets the number of iterations to run
                         thin=5,
                         chains=3)

xsum <- locationSummary(fit$x,collapse = TRUE)
zsum <- locationSummary(fit$z,collapse = TRUE)

x0 = cbind(xsum$'Lon.50%',xsum$'Lat.50%')
z0 = cbind(zsum$'Lon.50%',zsum$'Lat.50%')

fit <- estelleMetropolis(model = model,
                              proposal.x = proposal.x,
                              proposal.z = proposal.z,
                              x0 = x0,
                              z0 = z0,
                              iters=2000, # This value sets the number of iterations to run
                              thin=5,
                              chains=3)

#Final Run

xsum<-locationSummary(fit$x,collapse = TRUE)
zsum<-locationSummary(fit$z,collapse = TRUE)

proposal.x <- mvnorm(chainCov(fit$x),s=0.1)
proposal.z <- mvnorm(chainCov(fit$z),s=0.1)

fit <- estelleMetropolis(model = model,
                         proposal.x = proposal.x,
                         proposal.z = proposal.z,
                         x0 = x0,
                         z0 = z0,
                         iters=5000, # This value sets the number of iterations to run
                         thin=10,
                         chains=3)

#Plot--------------------------------------------------------------------------------------------------------------

#This step makes an empty raster #
r <- raster::raster(xmn = xlim[1],
            	  xmx = xlim[2], 
           		  ymn = ylim[1], 
                    ymx = ylim[2],
                    res = c(0.25,0.25))

# set plot parameters 
xlim = c(-180, 0)
ylim = c(-60, 80)

### Here is another way to determine the locations and where the bird was ###

S <- SGAT::slices(type="intermediate",
                 breaks="day",
                 weights = rep(0.5,length(fit$model$time)),
                 mcmc=fit,
                 grid=r)

## Quick lightImage plot to estimate stationary dates ## 
lightImage(tagdata = d.lux,
           offset = 18,
           zlim = c(0,5))
tsimagePoints(cal.data$Twilight, 
              offset = 19, 
              pch = 16, cex = 0.5, 
              col = ifelse(twl$Rise, "red", "blue"))
tsimagePoints(twl$Twilight, 
              offset = 19, 
              pch = 16, cex = 0.5, 
              col = ifelse(twl$Rise, "blue", "red"))

# These stationary dates were estimated from the lightImage #
# you can add as many of these as you want - each stationary time
# should be a row 

stationary.periods <- data.frame(start = "2016-11-01",
                                 stop = "2017-03-01")

SOSA <- LLmig::MigSchedule(MCMC = S, 
                           mig.quantile = 0.95,
                           stationary.periods = stationary.periods,
                           stationary.duration = 1,
                           prob = 0.95,
                           rm.lat.equinox = FALSE,
                           days.omit = 0)

# See what the MigSchedule function produces - 
str(SOSA,1)

#Wish List-------------------------------------------------------------------------------------------------
#Plot maps for a single southbound and single north bound migration route.
#Determine and then Plot important stop-over locations.

# Here is how to go about doing that using the results from the MigSchedule function 
# Take a look at the migration schedule
SOSA$Schedule

FallMig <- raster::spLines(sp::SpatialPoints(cbind(c(CapLocs[1],SOSA$Schedule[1:7,3]),# Longitudes
                                                   c(CapLocs[2],SOSA$Schedule[1:7,6]))))


SpringMig <- raster::spLines(sp::SpatialPoints(SOSA$Schedule[8:15,c(3,6)]))

plot(FallMig)
for(i in 1:7){
plot(SOSA$movements[[i]],add = TRUE,legend = FALSE)
}
plot(wrld_simpl, add = TRUE)
points(CapLocs[1],CapLocs[2],cex = 2,pch = 19)

plot(SpringMig, add = TRUE, col = "red")
for(i in 8:15){
plot(SOSA$movements[[i]],
     col = rev(bpy.colors(25)),
     add = TRUE,
     legend = FALSE)
}

#Plot maps on satellite imagery with state and country borders.
# There are quite a few ways to do this - but here is one way. 
# Obtain a BlueMarble image manually then plot the results onto that 
# The link below is a recent one - not much now cover
# https://neo.sci.gsfc.nasa.gov/view.php?datasetId=BlueMarbleNG-TB&date=2004-07-01 #

#Download and open raster.
list.files(path = "T:/McDuffie/SOSA_GEO/SOSA_R/SOSA_589", pattern = NULL,
           full.names = TRUE)
BlueMarb <- raster::brick("T:/McDuffie/SOSA_GEO/SOSA_R/SOSA_589/BlueMarbleNG.tiff") 

# Plot fall migration to set plot region #
plot(FallMig)

# add map using real-world colors #
plotRGB(BlueMarb, colNA = "black", maxpixels = 6480000, add = TRUE)

# add the migration path back #
plot(FallMig, add = TRUE, col = "red")

# add country borders # 
plot(wrld_simpl, add = TRUE, border = "yellow")

# add one of the stop-over locations #
plot(SOSA$movements[[5]],add = TRUE, col = bpy.colors(20), alpha = 0.4, legend = FALSE)

##### Here is another way to get border shapefiles ###
# You can get the data directly from the raster package #

# raster::getData("countries") # this returns an error but the data are available here 
                               #'http://biogeo.ucdavis.edu/data/gadm2.8/countries_gadm28.rds'

#Canada <- raster::getData("GADM", country = "Canada", level = 1)
#States <- raster::getData("GADM", country = "United States", level = 1)

#plot(Canada)
#plot(States, add = TRUE)

# Here is a third way # 

plot(dismo::gmap(FallMig, type = "satellite", lonlat = TRUE))
plot(FallMig, add = TRUE, col = "red")

#Set tolerance values for South American equinox to remove erroneous migration paths.
# See above

#Questions-------------------------------------------------------------------------------------------------------
#What is the ideal number of iterations to have for the MCMC model?
# The way you have it set up now - 2000 iterations - then summarize, then 2000 then summarize then 5000 
# is how I've done it in the past - you can do more but as long as you summarize between runs you 
# are using all the data. A few thousand is probably enough - it looks good for this bird. 

#During incubation, light-levels get skewed. Can I remove data points during the incubation period?
# I set those locations to fixed  - another option is to just remove them. Both methods work and 
# acheive pretty much the same thing. 









opar <- par(mar=c(2,2,2,2)+0.1)

# summary location from the model - 
s <- locationSummary(fit$x,time=model$time,collapse=T)

# set fixed locations - 
fixedz<- fixedx[-length(fixedx)] > 0 & fixedx[-length(fixedx)]==fixedx[-1]

# if fixed - set value to 0 if not use model dt
dt <- ifelse(fixedz,0,model$dt)

# Create the location image using the intermediate location
im <- locationImage(fit$z,xlim=xlim,ylim=ylim,nx=4*diff(xlim),ny=4*diff(ylim),
                    weight=dt,collapse=TRUE)

# Plot
opar <- par(mar=c(2,2,2,2)+0.1)
plot(NA, xlim = range(im$x), ylim = range(im$y), xlab = "Longitude", ylab = "Latitude")
plot(wrld_simpl,col= "grey90",border="grey10", add = T)
image(im$x,im$y,im$W,xlab="",ylab="",cex.axis=0.7, add = T, col = c("transparent", rev(topo.colors(100))))

plot(wrld_simpl, add = T)

lines(s$Lon.mean, s$Lat.mean,col=rgb(t(col2rgb("firebrick"))/255,alpha=0.5))
points(s$Lon.mean, s$Lat.mean,pch=16,cex=0.5,col=rgb(t(col2rgb("firebrick"))/255,alpha=0.2))
points(cap.long, cap.lat, pch = 22, bg = "white")
par(opar)

## Let's view the data by month as 2D histograms.
## First create an empty background raster
r <- raster(nrows = 80, ncols = 100, xmn = xlim[1], xmx = xlim[2], ymn = ylim[1], ymx = ylim[2])
## The slices function in SGAT can cut your data up into days, weeks, months, quarters, years.
s <- slices(type = "intermediate", breaks = "month", mcmc = fit, grid = r)

opar <- par(mfrow = c(4,3), mar = c(0,0,0,0), oma = c(0.4,0.4,0.4,0.4))
for (k in sliceIndices(s)) {   #loop through the slices to plot the data
  tm <- sliceInterval(s, k)    
  sk <- slice(s, k)
  plot(NA, xlim = range(im$x), ylim = range(im$y), xaxt = "n", yaxt = "n", xlab = "", ylab  = "")
  plot(wrld_simpl, col = "grey95", add = T)
  plot(sk, add = T, col = rev(heat.colors(64, alpha = 0.7)), useRaster = F, legend = F)
  
  mtext(sprintf("%s - %s", as.Date(tm[1]), as.Date(tm[2])), 1, line = -2, cex = 0.8)
}
par(opar)


