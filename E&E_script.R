#-------------Ground-truthing reveals limited predictive power of species distribution models for lizards------------
#------Setup work space----
#Clear work space
rm(list=ls(all=TRUE))

#Set working directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) #Sets to the location that the script is contained. Need to install rstudioapi if not previously installed.

#Packages
library(spocc)
library(ENMeval)
library(maptools)
library(rgeos)
library(raster)
library(sp)
library(rgdal)
library(sf)
library(rasterVis)
library(dplyr)
library(car)
library(pROC)
library(terra)
library(ecospat)

##------------------------Load all data into work space------------------------
###----------------------------Herptofauna database------------------------
#import df
all_interest_trimmed<-read.csv("jittered_all_interest.csv")
all_interest_trimmed$latitude<-as.numeric(all_interest_trimmed$latitude)
all_interest_trimmed$longitude<-as.numeric(all_interest_trimmed$longitude)

#Creating the subsets to build models with, based on how I have named the db
#GRA only clade 5
GRA<-subset(all_interest_trimmed,all_interest_trimmed$ScientificName=="Oligosoma aff. polychroma Clade 5" )

#MAC including lineoocellatum
MACplus<-subset(all_interest_trimmed,all_interest_trimmed$ScientificName=="Oligosoma prasinum"|all_interest_trimmed$ScientificName=="SPOT")

#MCC
MCC<-subset(all_interest_trimmed,all_interest_trimmed$ScientificName=="Oligosoma maccanni")

#RMM southern 
RMM<-subset(all_interest_trimmed,all_interest_trimmed$ScientificName=="RMM") 

#SAG woodworthia SA only
SAG<-subset(all_interest_trimmed,all_interest_trimmed$ScientificName=="SAG" )

#SCR only waimatense
SCR<-subset(all_interest_trimmed,all_interest_trimmed$ScientificName=="Oligosoma waimatense")

#RMM+Long
RMMlong<-subset(all_interest_trimmed,all_interest_trimmed$ScientificName=="RMM"|all_interest_trimmed$ScientificName=="Oligosoma longipes") 

rm(all_interest_trimmed)

###----------------------------Independent testing data------------------------

#Read in independent dataset
A<-read.csv("jittered_independent.csv")
A$Type<-as.factor(A$Type)
#Removing non trapping data
PresenceTEST<-A[1:9]
PresenceTEST<-subset(PresenceTEST,PresenceTEST$Type=='Presence'|PresenceTEST$Type=='Absence')

PresenceTEST$GRA = ifelse(is.na(PresenceTEST$Species), 0,ifelse(PresenceTEST$Species == "GRA",1,0))
PresenceTEST$MACplus = ifelse(is.na(PresenceTEST$Species), 0,ifelse(PresenceTEST$Species == "MAC",1,0))
PresenceTEST$MCC = ifelse(is.na(PresenceTEST$Species), 0,ifelse(PresenceTEST$Species == "MCC",1,0))
PresenceTEST$RMM = ifelse(is.na(PresenceTEST$Species), 0,ifelse(PresenceTEST$Species == "RMM",1,0))
PresenceTEST$SAG = ifelse(is.na(PresenceTEST$Species), 0,ifelse(PresenceTEST$Species == "SAG",1,0))
PresenceTEST$SCR = ifelse(is.na(PresenceTEST$Species), 0,ifelse(PresenceTEST$Species == "SCR",1,0))

####Trans_test------------
xyTEST<-PresenceTEST[,c(4,3)]
SpatialTEST<-SpatialPointsDataFrame(coords=xyTEST,proj4string = CRS ("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"), data = PresenceTEST)

rm(A,xyTEST)
###------------------------------Input rasters------------------------
#Input rasters
Aspect_North<-raster("./Rasters/NZEnvDS_V1.1_finalNZTM/final_layers_nztm/aspect_northness.tif")
Slope_Deg<-raster("./Rasters/NZEnvDS_V1.1_finalNZTM/final_layers_nztm/slope_deg.tif")
Elevation<-raster("./Rasters/NZEnvDS_V1.1_finalNZTM/final_layers_nztm/elevation.tif")
Temp_Min<-raster("./Rasters/NZEnvDS_V1.1_finalNZTM/final_layers_nztm/temp_minColdMonth.tif")
Temp_Max<-raster("./Rasters/NZEnvDS_V1.1_finalNZTM/final_layers_nztm/temp_maxWarmMonth.tif")
Geomorphons<-raster("./Rasters/NZEnvDS_V1.1_finalNZTM/final_layers_nztm/topo_geomorphons.tif")
sol_Ann<-raster("./Rasters/NZEnvDS_V1.1_finalNZTM/final_layers_nztm/solRad_meanAnn.tif")
precip_Ann<-raster("./Rasters/NZEnvDS_V1.1_finalNZTM/final_layers_nztm/precip_ann.tif")
landuse<-raster("./Rasters/Land_use/LandUse.tif")

#create raster stack
stack<-stack(Aspect_North,Slope_Deg,Elevation,Temp_Max,Temp_Min,Geomorphons,sol_Ann,precip_Ann,landuse)

#declare categorical variables
stack$topo_geomorphons<-as.factor(stack$topo_geomorphons)
stack$LandUse<-as.factor(stack$LandUse)

rm(Aspect_North,Slope_Deg,Elevation,Temp_Max,Temp_Min,Geomorphons,sol_Ann,precip_Ann,landuse)

##------------------------Species specific model set up------------------------
###-----------------------------GRA--------------------
#Full annotation is provided for grass skinks (GRA), but removed for the subsequent species which use an identical procedure
#Reading two columns of my data into spatial df MUST specifiy what projection it is in before projecting it
GRA_sp<-SpatialPointsDataFrame(coords = GRA[,c(5,4)],
                               data=GRA,
                               proj4string = CRS("+init=EPSG:4326"))

#Transforming/projecting the data to same as raster stack
GRA_sp<-spTransform(GRA_sp,
                    crs(stack[[1]]))

#Returning just the coordinate columns to a df
GRA_df<-as.data.frame(GRA_sp)
GRA_chck<-GRA_df[,c(12,13)]

#TEST of the data and raster 
plot(stack$aspect_northness)
points(GRA_chck)

#creating a boundary box around my spatial points
GRA_bb<-bbox(GRA_sp)

#Extended buffer to create a study zone by 50 kilometres 
GRA_bb.buf <- extent(GRA_bb[1]-50000, GRA_bb[3]+50000, GRA_bb[2]-50000, GRA_bb[4]+50000)

#crop rasters to new boundary box - TAKES TIME
GRA_stack<-crop(stack,GRA_bb.buf)

#Creating a simple frame object that the next code can read (spatial deos not work)
GRA_sf <- sf::st_as_sf(GRA_sp, coords = c("longitude","latitude"), crs = raster::crs(GRA_stack))

#Creating a buffer around points that will become our crop and study area
GRA_buf <- sf::st_buffer(GRA_sf, dist = 50000) %>% 
  sf::st_union() %>% 
  sf::st_sf() %>%
  sf::st_transform(crs = raster::crs(GRA_stack))

#check with plot - # To add sf objects to a plot, use add = TRUE
plot(GRA_stack[[1]],
     main = names(GRA_stack)[1])
points(GRA_chck)
plot(GRA_buf, border = "blue", lwd = 2, add = TRUE)

# Crop environmental rasters to match the study extent
GRA_bg <- raster::crop(GRA_stack, GRA_buf)

# Next, mask the rasters to the shape of the buffers
GRA_bg <- raster::mask(GRA_bg, GRA_buf)

#Check mask
plot(GRA_bg[[1]], main = names(GRA_bg)[1])
points(GRA_chck)
plot(GRA_buf, border = "blue", lwd = 3, add = TRUE)

#set seed to generate background points
set.seed(1)

GRA_bg_points <- dismo::randomPoints(GRA_bg[[1]], n = 10000) %>% as.data.frame()

# Notice how we have pretty good coverage (every cell).
plot(GRA_bg[[1]])
points(GRA_bg_points, pch = 20, cex = 0.2)

#Make the columns names the same
colnames(GRA_chck) <- c("x", "y")


rm(GRA_df,GRA_sp,GRA_buf,GRA_sf,GRA_bb,GRA_bb.buf,GRA_stack)

###-----------------------MACplus----------------------
MACplus_sp<-SpatialPointsDataFrame(coords = MACplus[,c(5,4)],
                               data=MACplus,
                               proj4string = CRS("+init=EPSG:4326"))

MACplus_sp<-spTransform(MACplus_sp,
                    crs(stack[[1]]))
MACplus_df<-as.data.frame(MACplus_sp)
MACplus_chck<-MACplus_df[,c(12,13)]


plot(stack$aspect_northness)
points(MACplus_chck)

MACplus_bb<-bbox(MACplus_sp)
MACplus_bb.buf <- extent(MACplus_bb[1]-50000, MACplus_bb[3]+50000, MACplus_bb[2]-50000, MACplus_bb[4]+50000)
MACplus_stack<-crop(stack,MACplus_bb.buf)

MACplus_sf <- sf::st_as_sf(MACplus_sp, coords = c("longitude","latitude"), crs = raster::crs(MACplus_stack))
MACplus_buf <- sf::st_buffer(MACplus_sf, dist = 50000) %>% 
  sf::st_union() %>% 
  sf::st_sf() %>%
  sf::st_transform(crs = raster::crs(MACplus_stack))

plot(MACplus_stack[[1]],
     main = names(MACplus_stack)[1])
points(MACplus_chck)
plot(MACplus_buf, border = "blue", lwd = 2, add = TRUE)

MACplus_bg <- raster::crop(MACplus_stack, MACplus_buf)
MACplus_bg <- raster::mask(MACplus_bg, MACplus_buf)

plot(MACplus_bg[[1]], main = names(MACplus_bg)[1])
points(MACplus_chck)
plot(MACplus_buf, border = "blue", lwd = 3, add = TRUE)

set.seed(1)

MACplus_bg_points <- dismo::randomPoints(MACplus_bg[[1]], n = 10000) %>% as.data.frame()

plot(MACplus_bg[[1]])
points(MACplus_bg_points, pch = 20, cex = 0.2)

colnames(MACplus_chck) <- c("x", "y")

rm(MACplus_df,MACplus_sp,MACplus_buf,MACplus_sf,MACplus_bb,MACplus_bb.buf,MACplus_stack)

###-----------------------------MCC--------------------
MCC_sp<-SpatialPointsDataFrame(coords = MCC[,c(5,4)],
                                   data=MCC,
                                   proj4string = CRS("+init=EPSG:4326"))

MCC_sp<-spTransform(MCC_sp,
                        crs(stack[[1]]))
MCC_df<-as.data.frame(MCC_sp)
MCC_chck<-MCC_df[,c(12,13)]


plot(stack$aspect_northness)
points(MCC_chck)

MCC_bb<-bbox(MCC_sp)
MCC_bb.buf <- extent(MCC_bb[1]-50000, MCC_bb[3]+50000, MCC_bb[2]-50000, MCC_bb[4]+50000)
MCC_stack<-crop(stack,MCC_bb.buf)

MCC_sf <- sf::st_as_sf(MCC_sp, coords = c("longitude","latitude"), crs = raster::crs(MCC_stack))
MCC_buf <- sf::st_buffer(MCC_sf, dist = 50000) %>% 
  sf::st_union() %>% 
  sf::st_sf() %>%
  sf::st_transform(crs = raster::crs(MCC_stack))

plot(MCC_stack[[1]],
     main = names(MCC_stack)[1])
points(MCC_chck)
plot(MCC_buf, border = "blue", lwd = 2, add = TRUE)

MCC_bg <- raster::crop(MCC_stack, MCC_buf)
MCC_bg <- raster::mask(MCC_bg, MCC_buf)

plot(MCC_bg[[1]], main = names(MCC_bg)[1])
points(MCC_chck)
plot(MCC_buf, border = "blue", lwd = 3, add = TRUE)

set.seed(1)

MCC_bg_points <- dismo::randomPoints(MCC_bg[[1]], n = 10000) %>% as.data.frame()

plot(MCC_bg[[1]])
points(MCC_bg_points, pch = 20, cex = 0.2)

colnames(MCC_chck) <- c("x", "y")

rm(MCC_df,MCC_sp,MCC_buf,MCC_sf,MCC_bb,MCC_bb.buf,MCC_stack)

###-----------------------------RMM--------------------
RMM_sp<-SpatialPointsDataFrame(coords = RMM[,c(5,4)],
                               data=RMM,
                               proj4string = CRS("+init=EPSG:4326"))

RMM_sp<-spTransform(RMM_sp,
                    crs(stack[[1]]))
RMM_df<-as.data.frame(RMM_sp)
RMM_chck<-RMM_df[,c(12,13)]


plot(stack$aspect_northness)
points(RMM_chck)

RMM_bb<-bbox(RMM_sp)
RMM_bb.buf <- extent(RMM_bb[1]-50000, RMM_bb[3]+50000, RMM_bb[2]-50000, RMM_bb[4]+50000)
RMM_stack<-crop(stack,RMM_bb.buf)

RMM_sf <- sf::st_as_sf(RMM_sp, coords = c("longitude","latitude"), crs = raster::crs(RMM_stack))
RMM_buf <- sf::st_buffer(RMM_sf, dist = 50000) %>% 
  sf::st_union() %>% 
  sf::st_sf() %>%
  sf::st_transform(crs = raster::crs(RMM_stack))

plot(RMM_stack[[1]],
     main = names(RMM_stack)[1])
points(RMM_chck)
plot(RMM_buf, border = "blue", lwd = 2, add = TRUE)

RMM_bg <- raster::crop(RMM_stack, RMM_buf)
RMM_bg <- raster::mask(RMM_bg, RMM_buf)

plot(RMM_bg[[1]], main = names(RMM_bg)[1])
points(RMM_chck)
plot(RMM_buf, border = "blue", lwd = 3, add = TRUE)

set.seed(1)

RMM_bg_points <- dismo::randomPoints(RMM_bg[[1]], n = 10000) %>% as.data.frame()

plot(RMM_bg[[1]])
points(RMM_bg_points, pch = 20, cex = 0.2)

colnames(RMM_chck) <- c("x", "y")

rm(RMM_df,RMM_sp,RMM_buf,RMM_sf,RMM_bb,RMM_bb.buf,RMM_stack)


###-----------------------------SAG--------------------
SAG_sp<-SpatialPointsDataFrame(coords = SAG[,c(5,4)],
                               data=SAG,
                               proj4string = CRS("+init=EPSG:4326"))

SAG_sp<-spTransform(SAG_sp,
                    crs(stack[[1]]))
SAG_df<-as.data.frame(SAG_sp)
SAG_chck<-SAG_df[,c(12,13)]


plot(stack$aspect_northness)
points(SAG_chck)

SAG_bb<-bbox(SAG_sp)
SAG_bb.buf <- extent(SAG_bb[1]-50000, SAG_bb[3]+50000, SAG_bb[2]-50000, SAG_bb[4]+50000)
SAG_stack<-crop(stack,SAG_bb.buf)

SAG_sf <- sf::st_as_sf(SAG_sp, coords = c("longitude","latitude"), crs = raster::crs(SAG_stack))
SAG_buf <- sf::st_buffer(SAG_sf, dist = 50000) %>% 
  sf::st_union() %>% 
  sf::st_sf() %>%
  sf::st_transform(crs = raster::crs(SAG_stack))

plot(SAG_stack[[1]],
     main = names(SAG_stack)[1])
points(SAG_chck)
plot(SAG_buf, border = "blue", lwd = 2, add = TRUE)

SAG_bg <- raster::crop(SAG_stack, SAG_buf)
SAG_bg <- raster::mask(SAG_bg, SAG_buf)

plot(SAG_bg[[1]], main = names(SAG_bg)[1])
points(SAG_chck)
plot(SAG_buf, border = "blue", lwd = 3, add = TRUE)

set.seed(1)

SAG_bg_points <- dismo::randomPoints(SAG_bg[[1]], n = 10000) %>% as.data.frame()

plot(SAG_bg[[1]])
points(SAG_bg_points, pch = 20, cex = 0.2)

colnames(SAG_chck) <- c("x", "y")

rm(SAG_df,SAG_sp,SAG_buf,SAG_sf,SAG_bb,SAG_bb.buf,SAG_stack)

###-----------------------------SCR--------------------
SCR_sp<-SpatialPointsDataFrame(coords = SCR[,c(5,4)],
                               data=SCR,
                               proj4string = CRS("+init=EPSG:4326"))

SCR_sp<-spTransform(SCR_sp,
                    crs(stack[[1]]))
SCR_df<-as.data.frame(SCR_sp)
SCR_chck<-SCR_df[,c(12,13)]


plot(stack$aspect_northness)
points(SCR_chck)

SCR_bb<-bbox(SCR_sp)
SCR_bb.buf <- extent(SCR_bb[1]-50000, SCR_bb[3]+50000, SCR_bb[2]-50000, SCR_bb[4]+50000)
SCR_stack<-crop(stack,SCR_bb.buf)

SCR_sf <- sf::st_as_sf(SCR_sp, coords = c("longitude","latitude"), crs = raster::crs(SCR_stack))
SCR_buf <- sf::st_buffer(SCR_sf, dist = 50000) %>% 
  sf::st_union() %>% 
  sf::st_sf() %>%
  sf::st_transform(crs = raster::crs(SCR_stack))

plot(SCR_stack[[1]],
     main = names(SCR_stack)[1])
points(SCR_chck)
plot(SCR_buf, border = "blue", lwd = 2, add = TRUE)

SCR_bg <- raster::crop(SCR_stack, SCR_buf)
SCR_bg <- raster::mask(SCR_bg, SCR_buf)

plot(SCR_bg[[1]], main = names(SCR_bg)[1])
points(SCR_chck)
plot(SCR_buf, border = "blue", lwd = 3, add = TRUE)

set.seed(1)

SCR_bg_points <- dismo::randomPoints(SCR_bg[[1]], n = 10000) %>% as.data.frame()

plot(SCR_bg[[1]])
points(SCR_bg_points, pch = 20, cex = 0.2)

colnames(SCR_chck) <- c("x", "y")

rm(SCR_df,SCR_sp,SCR_buf,SCR_sf,SCR_bb,SCR_bb.buf,SCR_stack)

###-----------------------------RMM+LONG--------------------
RMMlong_sp<-SpatialPointsDataFrame(coords = RMMlong[,c(5,4)],
                               data=RMMlong,
                               proj4string = CRS("+init=EPSG:4326"))

RMMlong_sp<-spTransform(RMMlong_sp,crs(stack[[1]]))

RMMlong_df<-as.data.frame(RMMlong_sp)
RMMlong_chck<-RMMlong_df[,c(12,13)]


plot(stack$aspect_northness)
points(RMMlong_chck)

RMMlong_bb<-bbox(RMMlong_sp)
RMMlong_bb.buf <- extent(RMMlong_bb[1]-50000, RMMlong_bb[3]+50000, RMMlong_bb[2]-50000, RMMlong_bb[4]+50000)
RMMlong_stack<-crop(stack,RMMlong_bb.buf)

RMMlong_sf <- sf::st_as_sf(RMMlong_sp, coords = c("longitude","latitude"), crs = raster::crs(RMMlong_stack))
RMMlong_buf <- sf::st_buffer(RMMlong_sf, dist = 50000) %>% 
  sf::st_union() %>% 
  sf::st_sf() %>%
  sf::st_transform(crs = raster::crs(RMMlong_stack))

plot(RMMlong_stack[[1]],
     main = names(RMMlong_stack)[1])
points(RMMlong_chck)
plot(RMMlong_buf, border = "blue", lwd = 2, add = TRUE)

RMMlong_bg <- raster::crop(RMMlong_stack, RMMlong_buf)
RMMlong_bg <- raster::mask(RMMlong_bg, RMMlong_buf)

plot(RMMlong_bg[[1]], main = names(RMMlong_bg)[1])
points(RMMlong_chck)
plot(RMMlong_buf, border = "blue", lwd = 3, add = TRUE)

set.seed(1)

RMMlong_bg_points <- dismo::randomPoints(RMMlong_bg[[1]], n = 10000) %>% as.data.frame()

plot(RMMlong_bg[[1]])
points(RMMlong_bg_points, pch = 20, cex = 0.2)

colnames(RMMlong_chck) <- c("x", "y")

rm(RMMlong_df,RMMlong_sp,RMMlong_buf,RMMlong_sf,RMMlong_bb,RMMlong_bb.buf,RMMlong_stack)


#------------------------MaxENT modelling all species------------------------
##---------------------------Full models-------------------------------------

GRA_mod<- ENMevaluate(occs = GRA_chck, envs = GRA_bg, bg=GRA_bg_points, algorithm = 'maxnet', partitions = 'checkerboard1',categoricals = c("topo_geomorphons", "LandUse"),tune.args = list(fc=c("L", "Q","LQ"), rm = 1:5))

MACplus_mod<- ENMevaluate(occs = MACplus_chck, envs = MACplus_bg, bg=MACplus_bg_points, algorithm = 'maxnet', partitions = 'checkerboard1',categoricals = c("topo_geomorphons", "LandUse"),tune.args = list(fc=c("L","Q","LQ"), rm = 1:5))

MCC_mod<- ENMevaluate(occs = MCC_chck, envs = MCC_bg, bg=MCC_bg_points, algorithm = 'maxnet', partitions = 'checkerboard1',categoricals = c("topo_geomorphons", "LandUse"),tune.args = list(fc=c("L", "Q","LQ"), rm = 1:5))

RMM_mod<- ENMevaluate(occs = RMM_chck, envs = RMM_bg, bg=RMM_bg_points, algorithm = 'maxnet', partitions = 'checkerboard1',categoricals = c("topo_geomorphons", "LandUse"),tune.args = list(fc=c("L", "Q","LQ"), rm = 1:5))

SAG_mod<- ENMevaluate(occs = SAG_chck, envs = SAG_bg, bg=SAG_bg_points, algorithm = 'maxnet', partitions = 'checkerboard1',categoricals = c("topo_geomorphons", "LandUse"),tune.args = list(fc=c("L", "Q","LQ"), rm = 1:5))

SCR_mod<- ENMevaluate(occs = SCR_chck, envs = SCR_bg, bg=SCR_bg_points, algorithm = 'maxnet', partitions = 'checkerboard1',categoricals = c("topo_geomorphons", "LandUse"),tune.args = list(fc=c("L", "Q","LQ"), rm = 1:5))

RMMlong_mod<- ENMevaluate(occs = RMMlong_chck, envs = RMMlong_bg, bg=RMMlong_bg_points, algorithm = 'maxnet', partitions = 'checkerboard1',categoricals = c("topo_geomorphons", "LandUse"),tune.args = list(fc=c("L", "Q","LQ"), rm = 1:5))

##-----------Dropping categorical levels to test over fitting issues-------------
GRA_bg_DROPPED <- dropLayer(GRA_bg,c("LandUse","topo_geomorphons"))
GRA_mod_DROPPED<- ENMevaluate(occs = GRA_chck, envs = GRA_bg_DROPPED, bg=GRA_bg_points, algorithm = 'maxnet', partitions = 'checkerboard1',tune.args = list(fc=c("L","Q","LQ"), rm = 1:5))
rm(GRA_bg_DROPPED)

MACplus_bg_DROPPED <- dropLayer(MACplus_bg,c("LandUse","topo_geomorphons"))
MACplus_mod_DROPPED<- ENMevaluate(occs = MACplus_chck, envs = MACplus_bg_DROPPED, bg=MACplus_bg_points, algorithm = 'maxnet', partitions = 'checkerboard1',tune.args = list(fc=c("L","Q","LQ"), rm = 1:5))
rm(MACplus_bg_DROPPED)

MCC_bg_DROPPED <- dropLayer(MCC_bg,c("LandUse","topo_geomorphons"))
MCC_mod_DROPPED<- ENMevaluate(occs = MCC_chck, envs = MCC_bg_DROPPED, bg=MCC_bg_points, algorithm = 'maxnet', partitions = 'checkerboard1',tune.args = list(fc=c("L","Q","LQ"), rm = 1:5))
rm(MCC_bg_DROPPED)

RMM_bg_DROPPED <- dropLayer(RMM_bg,c("LandUse","topo_geomorphons"))
RMM_mod_DROPPED<- ENMevaluate(occs = RMM_chck, envs = RMM_bg_DROPPED, bg=RMM_bg_points, algorithm = 'maxnet', partitions = 'checkerboard1',tune.args = list(fc=c("L","Q","LQ"), rm = 1:5))
rm(RMM_bg_DROPPED)

SAG_bg_DROPPED <- dropLayer(SAG_bg,c("LandUse","topo_geomorphons"))
SAG_mod_DROPPED<- ENMevaluate(occs = SAG_chck, envs = SAG_bg_DROPPED, bg=SAG_bg_points, algorithm = 'maxnet', partitions = 'checkerboard1',tune.args = list(fc=c("L","Q","LQ"), rm = 1:5))
rm(SAG_bg_DROPPED)

SCR_bg_DROPPED <- dropLayer(SCR_bg,c("LandUse","topo_geomorphons"))
SCR_mod_DROPPED<- ENMevaluate(occs = SCR_chck, envs = SCR_bg_DROPPED, bg=SCR_bg_points, algorithm = 'maxnet', partitions = 'checkerboard1',tune.args = list(fc=c("L","Q","LQ"), rm = 1:5))
rm(SCR_bg_DROPPED)

#--------------------------------------Evaluation--------------------------------
##-----------------------------------------GRA-----------------------------------
###-----------------------------Dependent data and validation--------------------
#check model
GRA_mod
#GRA_mod_DROPPED

#check str
str(GRA_mod, max.level=2)
#str(GRA_mod_DROPPED, max.level=2)

# Access algorithm, tuning settings, and partition method information.
eval.algorithm(GRA_mod)
eval.tune.settings(GRA_mod) %>% head()
eval.partition.method(GRA_mod)

#eval.algorithm(GRA_mod_DROPPED)
#eval.tune.settings(GRA_mod_DROPPED) %>% head()
#eval.partition.method(GRA_mod_DROPPED)
# Results table with summary statistics for cross validation on test data.
eval.results(GRA_mod) %>% head()
#eval.results(GRA_mod_DROPPED) %>% head()
# Results table with cross validation statistics for each test partition.
eval.results.partitions(GRA_mod)
#eval.results.partitions(GRA_mod_DROPPED)

# List of models with names corresponding to tune.args column label.
eval.models(GRA_mod) %>% str(max.level = 1)
#eval.models(GRA_mod_DROPPED) %>% str(max.level = 1)

# to tune.args column label.
eval.predictions(GRA_mod)
#eval.predictions(GRA_mod_DROPPED)

# Original occurrence data coordinates with associated predictor variable values.
eval.occs(GRA_mod) %>% head()
#eval.occs(GRA_mod_DROPPED) %>% head()

# Background data coordinates with associated predictor variable values.
eval.bg(GRA_mod) %>% head()
#eval.bg(GRA_mod_DROPPED) %>% head()

# Partition group assignments for occurrence data.
eval.occs.grp(GRA_mod) %>% str()
#eval.occs.grp(GRA_mod_DROPPED) %>% str()

#evaluation plots
evalplot.stats(e = GRA_mod, stats = c("or.10p", "auc.val"), color = "fc", x.var = "rm", error.bars = FALSE)
#evalplot.stats(e = GRA_mod_DROPPED, stats = c("or.10p", "auc.val"), color = "fc", x.var = "rm", error.bars = FALSE)

# Overall results
res <- eval.results(GRA_mod)
#res_DROPPED <- eval.results(GRA_mod_DROPPED)

# This dplyr operation executes the sequential criteria explained above.
opt.seq <- res %>% 
  filter(or.10p.avg == min(or.10p.avg)) %>% 
  filter(auc.val.avg == max(auc.val.avg))
opt.seq
tables<-as.data.frame(opt.seq)
tables$species<-"GRA"

#opt.seq_D <- res_DROPPED %>% 
  #filter(or.10p.avg == min(or.10p.avg)) %>% 
  #filter(auc.val.avg == max(auc.val.avg))
#opt.seq_D
#tables[2,1:19]<-opt.seq_D
#tables[2,"species"]<-"GRA_D"

# We can select a single model from the ENMevaluation object using the tune.args of our
# optimal model.
mod.seq <- eval.models(GRA_mod)[[opt.seq$tune.args]]
#mod.seq_D <- eval.models(GRA_mod_DROPPED)[[opt.seq_D$tune.args]]

# Here are the non-zero coefficients in our model.
mod.seq$betas
#mod.seq_D$betas
# And these are the marginal response curves for the predictor variables with non-zero 
# coefficients in our model. We define the y-axis to be the cloglog transformation, which
# is an approximation of occurrence probability (with assumptions) bounded by 0 and 1
# (Phillips et al. 2017).

m1.mx <- eval.models(GRA_mod)[["fc.LQ_rm.4"]]
m1.mx$betas
GRAmm<-m1.mx
plot(m1.mx, type="cloglog")
# model object above.
predGRA <- eval.predictions(GRA_mod)[[opt.seq$tune.args]]
plot(predGRA)
writeRaster(predGRA, "E:/Chapter_3 Modelling/QGIS/predGRA.tif")
#predGRA_D <- eval.predictions(GRA_mod_DROPPED)[[opt.seq_D$tune.args]]
#plot(predGRA_D)

# give the testing data the CRS of the prediction
transTEST<-spTransform(SpatialTEST,CRSobj= proj4string(predGRA))

###---------------*Independent data validation------------
####----------------**Modelling predictions------------------

#Creating prediction
PresenceTEST$GRAprediction_FULL<-extract(predGRA,transTEST)
#PresenceTEST$GRAprediction_DROPPED<-extract(predGRA_D,transTEST)

#AUC
tables[1,c("Ind_auc_low","Ind_auc","Ind_auc_high")] <- ci.auc(roc(PresenceTEST,response=GRA,predictor=GRAprediction_FULL))
#tables[2,c("Ind_auc_low","Ind_auc","Ind_auc_high")] <- ci.auc(roc(PresenceTEST,response=GRA,predictor=GRAprediction_DROPPED))

#plain roc curve
roc_GRA<-roc(PresenceTEST,response=GRA,predictor=GRAprediction_FULL)
#roc_GRA_D<-roc(PresenceTEST,response=GRA,predictor=GRAprediction_DROPPED)

#Determine value for maximum accuraccy threshold
roc_GRAs<-as.data.frame(roc_GRA$specificities)
roc_GRAs$sensitivities<-roc_GRA$sensitivities
roc_GRAs$thresholds<-roc_GRA$thresholds
colnames(roc_GRAs)[c(1,2,3)] <- c("spec", "sens",'thres')
roc_GRAs['sum'] = roc_GRAs['sens']+roc_GRAs['spec']

PresenceTEST$predTSS<-ifelse(PresenceTEST$GRAprediction_FULL>roc_GRAs[which(roc_GRAs$sum == max(roc_GRAs$sum)),"thres"],1,0)
PresenceTEST$correctprediction<-ifelse(PresenceTEST$predTSS==PresenceTEST$GRA,1,0)

#roc_GRAs_D<-as.data.frame(roc_GRA_D$specificities)
#roc_GRAs_D$sensitivities<-roc_GRA_D$sensitivities
#roc_GRAs_D$thresholds<-roc_GRA_D$thresholds
#colnames(roc_GRAs_D)[c(1,2,3)] <- c("spec", "sens",'thres')
#roc_GRAs_D['sum'] = roc_GRAs_D['sens']+roc_GRAs_D['spec']

#PresenceTEST$predTSS_D<-ifelse(PresenceTEST$GRAprediction_DROPPED>roc_GRAs_D[which(roc_GRAs_D$sum == max(roc_GRAs_D$sum)),"thres"],1,0)
#PresenceTEST$correctprediction_D<-ifelse(PresenceTEST$predTSS_D==PresenceTEST$GRA,1,0)

presenceGRA<-subset(PresenceTEST,PresenceTEST$GRA==1)
absenceGRA<-subset(PresenceTEST,PresenceTEST$GRA==0)

tables[1,"TSS_threshold"]<- roc_GRAs[which(roc_GRAs$sum == max(roc_GRAs$sum)),'thres']
#tables[2,"TSS_threshold"]<- roc_GRAs[which(roc_GRAs_D$sum == max(roc_GRAs_D$sum)),'thres']

tables[1,"TSS_Ind"]<-(length(subset(presenceGRA$correctprediction,presenceGRA$correctprediction==1))/nrow(presenceGRA))+ (length(subset(absenceGRA$correctprediction,absenceGRA$correctprediction==1))/nrow(absenceGRA))-1
#tables[2,"TSS_Ind"]<-(length(subset(presenceGRA$correctprediction_D,presenceGRA$correctprediction_D==1))/nrow(presenceGRA))+(length(subset(absenceGRA$correctprediction_D,absenceGRA$correctprediction_D==1))/nrow(absenceGRA))-1



##-----------------------MAC------------------------
###---------Dependent data and validation-----------
#check model
MACplus_mod
#MACplus_mod_DROPPED

#check str
str(MACplus_mod, max.level=2)
#str(MACplus_mod_DROPPED, max.level=2)

# Access algorithm, tuning settings, and partition method information.
eval.algorithm(MACplus_mod)
eval.tune.settings(MACplus_mod) %>% head()
eval.partition.method(MACplus_mod)

#eval.algorithm(MACplus_mod_DROPPED)
#eval.tune.settings(MACplus_mod_DROPPED) %>% head()
#eval.partition.method(MACplus_mod_DROPPED)
# Results table with summary statistics for cross validation on test data.
eval.results(MACplus_mod) %>% head()
#eval.results(MACplus_mod_DROPPED) %>% head()
# Results table with cross validation statistics for each test partition.
eval.results.partitions(MACplus_mod)
#eval.results.partitions(MACplus_mod_DROPPED)

# List of models with names corresponding to tune.args column label.
eval.models(MACplus_mod) %>% str(max.level = 1)
#eval.models(MACplus_mod_DROPPED) %>% str(max.level = 1)

# to tune.args column label.
eval.predictions(MACplus_mod)
#eval.predictions(MACplus_mod_DROPPED)

# Original occurrence data coordinates with associated predictor variable values.
eval.occs(MACplus_mod) %>% head()
#eval.occs(MACplus_mod_DROPPED) %>% head()

# Background data coordinates with associated predictor variable values.
eval.bg(MACplus_mod) %>% head()
#eval.bg(MACplus_mod_DROPPED) %>% head()

# Partition group assignments for occurrence data.
eval.occs.grp(MACplus_mod) %>% str()
#eval.occs.grp(MACplus_mod_DROPPED) %>% str()

#evaluation plots
evalplot.stats(e = MACplus_mod, stats = c("or.10p", "auc.val"), color = "fc", x.var = "rm", error.bars = FALSE)
#evalplot.stats(e = MACplus_mod_DROPPED, stats = c("or.10p", "auc.val"), color = "fc", x.var = "rm", error.bars = FALSE)

# Overall results
res <- eval.results(MACplus_mod)
#res_DROPPED <- eval.results(MACplus_mod_DROPPED)

# This dplyr operation executes the sequential criteria explained above.
opt.seq <- res %>% 
  filter(or.10p.avg == min(or.10p.avg)) %>% 
  filter(auc.val.avg == max(auc.val.avg))
opt.seq
tables[3,1:19]<-opt.seq
tables[3,"species"]<-"MAC"

#opt.seq_D <- res_DROPPED %>% 
  #filter(or.10p.avg == min(or.10p.avg)) %>% 
  #filter(auc.val.avg == max(auc.val.avg))
#opt.seq_D
#tables[4,1:19]<-opt.seq_D
#tables[4,"species"]<-"MAC_D"

# We can select a single model from the ENMevaluation object using the tune.args of our
# optimal model.
mod.seq <- eval.models(MACplus_mod)[[opt.seq$tune.args]]
#mod.seq_D <- eval.models(MACplus_mod_DROPPED)[[opt.seq_D$tune.args]]

# Here are the non-zero coefficients in our model.
mod.seq$betas
#mod.seq_D$betas
# And these are the marginal response curves for the predictor variables with non-zero 
# coefficients in our model. We define the y-axis to be the cloglog transformation, which
# is an approximation of occurrence probability (with assumptions) bounded by 0 and 1
# (Phillips et al. 2017).

m1.mx <- eval.models(MACplus_mod)[["fc.L_rm.1"]]
m1.mx$betas
MACmm<-m1.mx
# model object above.
predMACplus <- eval.predictions(MACplus_mod)[[opt.seq$tune.args]]
plot(predMACplus)

writeRaster(predMACplus, "E:/Chapter_3 Modelling/QGIS/predMACplus.tif")

#predMACplus_D <- eval.predictions(MACplus_mod_DROPPED)[[opt.seq_D$tune.args]]
#plot(predMACplus_D)

# give the testing data the CRS of the prediction
transTEST<-spTransform(SpatialTEST,CRSobj= proj4string(predMACplus))

###---------------*Independent data validation------------
####----------------**Modelling predictions------------------

#Creating prediction
PresenceTEST$MACplusprediction_FULL<-extract(predMACplus,transTEST)
#PresenceTEST$MACplusprediction_DROPPED<-extract(predMACplus_D,transTEST)

#AUC
tables[3,c("Ind_auc_low","Ind_auc","Ind_auc_high")] <- ci.auc(roc(PresenceTEST,response=MACplus,predictor=MACplusprediction_FULL))
#tables[4,c("Ind_auc_low","Ind_auc","Ind_auc_high")] <- ci.auc(roc(PresenceTEST,response=MACplus,predictor=MACplusprediction_DROPPED))

#plain roc curve
roc_MACplus<-roc(PresenceTEST,response=MACplus,predictor=MACplusprediction_FULL)
#roc_MACplus_D<-roc(PresenceTEST,response=MACplus,predictor=MACplusprediction_DROPPED)

#Determine value for maximum accuraccy threshold
roc_MACpluss<-as.data.frame(roc_MACplus$specificities)
roc_MACpluss$sensitivities<-roc_MACplus$sensitivities
roc_MACpluss$thresholds<-roc_MACplus$thresholds
colnames(roc_MACpluss)[c(1,2,3)] <- c("spec", "sens",'thres')
roc_MACpluss['sum'] = roc_MACpluss['sens']+roc_MACpluss['spec']

PresenceTEST$predTSS<-ifelse(PresenceTEST$MACplusprediction_FULL>roc_MACpluss[which(roc_MACpluss$sum == max(roc_MACpluss$sum)),"thres"],1,0)
PresenceTEST$correctprediction<-ifelse(PresenceTEST$predTSS==PresenceTEST$MACplus,1,0)

#roc_MACpluss_D<-as.data.frame(roc_MACplus_D$specificities)
#roc_MACpluss_D$sensitivities<-roc_MACplus_D$sensitivities
#roc_MACpluss_D$thresholds<-roc_MACplus_D$thresholds
#colnames(roc_MACpluss_D)[c(1,2,3)] <- c("spec", "sens",'thres')
#roc_MACpluss_D['sum'] = roc_MACpluss_D['sens']+roc_MACpluss_D['spec']

#PresenceTEST$predTSS_D<-ifelse(PresenceTEST$MACplusprediction_DROPPED>roc_MACpluss_D[which(roc_MACpluss_D$sum == max(roc_MACpluss_D$sum)),"thres"],1,0)
#PresenceTEST$correctprediction_D<-ifelse(PresenceTEST$predTSS_D==PresenceTEST$MACplus,1,0)

presenceMACplus<-subset(PresenceTEST,PresenceTEST$MACplus==1)
absenceMACplus<-subset(PresenceTEST,PresenceTEST$MACplus==0)

tables[3,"TSS_threshold"]<- roc_MACpluss[which(roc_MACpluss$sum == max(roc_MACpluss$sum)),'thres']
#tables[4,"TSS_threshold"]<- roc_MACpluss[which(roc_MACpluss_D$sum == max(roc_MACpluss_D$sum)),'thres']

tables[3,"TSS_Ind"]<-(length(subset(presenceMACplus$correctprediction,presenceMACplus$correctprediction==1))/nrow(presenceMACplus))+ (length(subset(absenceMACplus$correctprediction,absenceMACplus$correctprediction==1))/nrow(absenceMACplus))-1
#tables[4,"TSS_Ind"]<-(length(subset(presenceMACplus$correctprediction_D,presenceMACplus$correctprediction_D==1))/nrow(presenceMACplus))+(length(subset(absenceMACplus$correctprediction_D,absenceMACplus$correctprediction_D==1))/nrow(absenceMACplus))-1


##-----------------------MCC------------------------
###---------Dependent data and validation-----------
#check model
MCC_mod
#MCC_mod_DROPPED

#check str
str(MCC_mod, max.level=2)
#str(MCC_mod_DROPPED, max.level=2)

# Access algorithm, tuning settings, and partition method information.
eval.algorithm(MCC_mod)
eval.tune.settings(MCC_mod) %>% head()
eval.partition.method(MCC_mod)

#eval.algorithm(MCC_mod_DROPPED)
#eval.tune.settings(MCC_mod_DROPPED) %>% head()
#eval.partition.method(MCC_mod_DROPPED)
# Results table with summary statistics for cross validation on test data.
eval.results(MCC_mod) %>% head()
#eval.results(MCC_mod_DROPPED) %>% head()
# Results table with cross validation statistics for each test partition.
eval.results.partitions(MCC_mod)
#eval.results.partitions(MCC_mod_DROPPED)

# List of models with names corresponding to tune.args column label.
eval.models(MCC_mod) %>% str(max.level = 1)
#eval.models(MCC_mod_DROPPED) %>% str(max.level = 1)

# to tune.args column label.
eval.predictions(MCC_mod)
#eval.predictions(MCC_mod_DROPPED)

# Original occurrence data coordinates with associated predictor variable values.
eval.occs(MCC_mod) %>% head()
#eval.occs(MCC_mod_DROPPED) %>% head()

# Background data coordinates with associated predictor variable values.
eval.bg(MCC_mod) %>% head()
#eval.bg(MCC_mod_DROPPED) %>% head()

# Partition group assignments for occurrence data.
eval.occs.grp(MCC_mod) %>% str()
#eval.occs.grp(MCC_mod_DROPPED) %>% str()

#evaluation plots
evalplot.stats(e = MCC_mod, stats = c("or.10p", "auc.val"), color = "fc", x.var = "rm", error.bars = FALSE)
#evalplot.stats(e = MCC_mod_DROPPED, stats = c("or.10p", "auc.val"), color = "fc", x.var = "rm", error.bars = FALSE)

# Overall results
res <- eval.results(MCC_mod)
#res_DROPPED <- eval.results(MCC_mod_DROPPED)

# This dplyr operation executes the sequential criteria explained above.
opt.seq <- res %>% 
  filter(or.10p.avg == min(or.10p.avg)) %>% 
  filter(auc.val.avg == max(auc.val.avg))
opt.seq
tables[5,1:19]<-opt.seq
tables[5,"species"]<-"MCC"

#opt.seq_D <- res_DROPPED %>% 
 #filter(or.10p.avg == min(or.10p.avg)) %>% 
 #filter(auc.val.avg == max(auc.val.avg))
#opt.seq_D
#tables[6,1:19]<-opt.seq_D
#tables[6,"species"]<-"MCC_D"

# We can select a single model from the ENMevaluation object using the tune.args of our
# optimal model.
mod.seq <- eval.models(MCC_mod)[[opt.seq$tune.args]]
#mod.seq_D <- eval.models(MCC_mod_DROPPED)[[opt.seq_D$tune.args]]

# Here are the non-zero coefficients in our model.
mod.seq$betas
#mod.seq_D$betas
# And these are the marginal response curves for the predictor variables with non-zero 
# coefficients in our model. We define the y-axis to be the cloglog transformation, which
# is an approximation of occurrence probability (with assumptions) bounded by 0 and 1
# (Phillips et al. 2017).

m1.mx <- eval.models(MCC_mod)[["fc.LQ_rm.4"]]
m1.mx$betas
MCCmm<-m1.mx
# model object above.
predMCC <- eval.predictions(MCC_mod)[[opt.seq$tune.args]]
plot(predMCC)

writeRaster(predMCC, "E:/Chapter_3 Modelling/QGIS/predMCC.tif")
#predMCC_D <- eval.predictions(MCC_mod_DROPPED)[[opt.seq_D$tune.args]]
#plot(predMCC_D)

# give the testing data the CRS of the prediction
transTEST<-spTransform(SpatialTEST,CRSobj= proj4string(predMCC))

###---------------*Independent data validation------------
####----------------**Modelling predictions------------------

#Creating prediction
PresenceTEST$MCCprediction_FULL<-extract(predMCC,transTEST)
#PresenceTEST$MCCprediction_DROPPED<-extract(predMCC_D,transTEST)

#AUC
tables[5,c("Ind_auc_low","Ind_auc","Ind_auc_high")] <- ci.auc(roc(PresenceTEST,response=MCC,predictor=MCCprediction_FULL))
#tables[6,c("Ind_auc_low","Ind_auc","Ind_auc_high")] <- ci.auc(roc(PresenceTEST,response=MCC,predictor=MCCprediction_DROPPED))

#plain roc curve
roc_MCC<-roc(PresenceTEST,response=MCC,predictor=MCCprediction_FULL)
#roc_MCC_D<-roc(PresenceTEST,response=MCC,predictor=MCCprediction_DROPPED)

#Determine value for maximum accuraccy threshold
roc_MCCs<-as.data.frame(roc_MCC$specificities)
roc_MCCs$sensitivities<-roc_MCC$sensitivities
roc_MCCs$thresholds<-roc_MCC$thresholds
colnames(roc_MCCs)[c(1,2,3)] <- c("spec", "sens",'thres')
roc_MCCs['sum'] = roc_MCCs['sens']+roc_MCCs['spec']

PresenceTEST$predTSS<-ifelse(PresenceTEST$MCCprediction_FULL>roc_MCCs[which(roc_MCCs$sum == max(roc_MCCs$sum)),"thres"],1,0)
PresenceTEST$correctprediction<-ifelse(PresenceTEST$predTSS==PresenceTEST$MCC,1,0)

#roc_MCCs_D<-as.data.frame(roc_MCC_D$specificities)
#roc_MCCs_D$sensitivities<-roc_MCC_D$sensitivities
#roc_MCCs_D$thresholds<-roc_MCC_D$thresholds
#colnames(roc_MCCs_D)[c(1,2,3)] <- c("spec", "sens",'thres')
#roc_MCCs_D['sum'] = roc_MCCs_D['sens']+roc_MCCs_D['spec']

#PresenceTEST$predTSS_D<-ifelse(PresenceTEST$MCCprediction_DROPPED>roc_MCCs_D[which(roc_MCCs_D$sum == max(roc_MCCs_D$sum)),"thres"],1,0)
#PresenceTEST$correctprediction_D<-ifelse(PresenceTEST$predTSS_D==PresenceTEST$MCC,1,0)

presenceMCC<-subset(PresenceTEST,PresenceTEST$MCC==1)
absenceMCC<-subset(PresenceTEST,PresenceTEST$MCC==0)

tables[5,"TSS_threshold"]<- roc_MCCs[which(roc_MCCs$sum == max(roc_MCCs$sum)),'thres']
#tables[6,"TSS_threshold"]<- roc_MCCs[which(roc_MCCs_D$sum == max(roc_MCCs_D$sum)),'thres']

tables[5,"TSS_Ind"]<-(length(subset(presenceMCC$correctprediction,presenceMCC$correctprediction==1))/nrow(presenceMCC))+ (length(subset(absenceMCC$correctprediction,absenceMCC$correctprediction==1))/nrow(absenceMCC))-1
#tables[6,"TSS_Ind"]<-(length(subset(presenceMCC$correctprediction_D,presenceMCC$correctprediction_D==1))/nrow(presenceMCC))+(length(subset(absenceMCC$correctprediction_D,absenceMCC$correctprediction_D==1))/nrow(absenceMCC))-1

##-----------------------RMM------------------------
###---------Dependent data and validation-----------
#check model
RMM_mod
#RMM_mod_DROPPED

#check str
str(RMM_mod, max.level=2)
#str(RMM_mod_DROPPED, max.level=2)

# Access algorithm, tuning settings, and partition method information.
eval.algorithm(RMM_mod)
eval.tune.settings(RMM_mod) %>% head()
eval.partition.method(RMM_mod)

#eval.algorithm(RMM_mod_DROPPED)
#eval.tune.settings(RMM_mod_DROPPED) %>% head()
#eval.partition.method(RMM_mod_DROPPED)
# Results table with summary statistics for cross validation on test data.
eval.results(RMM_mod) %>% head()
#eval.results(RMM_mod_DROPPED) %>% head()
# Results table with cross validation statistics for each test partition.
eval.results.partitions(RMM_mod)
#eval.results.partitions(RMM_mod_DROPPED)

# List of models with names corresponding to tune.args column label.
eval.models(RMM_mod) %>% str(max.level = 1)
#eval.models(RMM_mod_DROPPED) %>% str(max.level = 1)

# to tune.args column label.
eval.predictions(RMM_mod)
#eval.predictions(RMM_mod_DROPPED)

# Original occurrence data coordinates with associated predictor variable values.
eval.occs(RMM_mod) %>% head()
#eval.occs(RMM_mod_DROPPED) %>% head()

# Background data coordinates with associated predictor variable values.
eval.bg(RMM_mod) %>% head()
#eval.bg(RMM_mod_DROPPED) %>% head()

# Partition group assignments for occurrence data.
eval.occs.grp(RMM_mod) %>% str()
#eval.occs.grp(RMM_mod_DROPPED) %>% str()

#evaluation plots
evalplot.stats(e = RMM_mod, stats = c("or.10p", "auc.val"), color = "fc", x.var = "rm", error.bars = FALSE)
#evalplot.stats(e = RMM_mod_DROPPED, stats = c("or.10p", "auc.val"), color = "fc", x.var = "rm", error.bars = FALSE)

# Overall results
res <- eval.results(RMM_mod)
#res_DROPPED <- eval.results(RMM_mod_DROPPED)

# This dplyr operation executes the sequential criteria explained above.
opt.seq <- res %>% 
  filter(or.10p.avg == min(or.10p.avg)) %>% 
  filter(auc.val.avg == max(auc.val.avg))
opt.seq
tables[7,1:19]<-opt.seq
tables[7,"species"]<-"RMM"

#opt.seq_D <- res_DROPPED %>% 
  #filter(or.10p.avg == min(or.10p.avg)) %>% 
  #filter(auc.val.avg == max(auc.val.avg))
#opt.seq_D
#tables[8,1:19]<-opt.seq_D
#tables[8,"species"]<-"RMM_D"

# We can select a single model from the ENMevaluation object using the tune.args of our
# optimal model.
mod.seq <- eval.models(RMM_mod)[[opt.seq$tune.args]]
#mod.seq_D <- eval.models(RMM_mod_DROPPED)[[opt.seq_D$tune.args]]

# Here are the non-zero coefficients in our model.
mod.seq$betas
#mod.seq_D$betas
# And these are the marginal response curves for the predictor variables with non-zero 
# coefficients in our model. We define the y-axis to be the cloglog transformation, which
# is an approximation of occurrence probability (with assumptions) bounded by 0 and 1
# (Phillips et al. 2017).

# model object above.
predRMM <- eval.predictions(RMM_mod)[[opt.seq$tune.args]]
plot(predRMM)

#predRMM_D <- eval.predictions(RMM_mod_DROPPED)[[opt.seq_D$tune.args]]
#plot(predRMM_D)

# give the testing data the CRS of the prediction
transTEST<-spTransform(SpatialTEST,CRSobj= proj4string(predRMM))

###---------------*Independent data validation------------
####----------------**Modelling predictions------------------

#Creating prediction
PresenceTEST$RMMprediction_FULL<-extract(predRMM,transTEST)
#PresenceTEST$RMMprediction_DROPPED<-extract(predRMM_D,transTEST)

#AUC
tables[7,c("Ind_auc_low","Ind_auc","Ind_auc_high")] <- ci.auc(roc(PresenceTEST,response=RMM,predictor=RMMprediction_FULL))
#tables[8,c("Ind_auc_low","Ind_auc","Ind_auc_high")] <- ci.auc(roc(PresenceTEST,response=RMM,predictor=RMMprediction_DROPPED))

#plain roc curve
roc_RMM<-roc(PresenceTEST,response=RMM,predictor=RMMprediction_FULL)
#roc_RMM_D<-roc(PresenceTEST,response=RMM,predictor=RMMprediction_DROPPED)

#Determine value for maximum accuraccy threshold
roc_RMMs<-as.data.frame(roc_RMM$specificities)
roc_RMMs$sensitivities<-roc_RMM$sensitivities
roc_RMMs$thresholds<-roc_RMM$thresholds
colnames(roc_RMMs)[c(1,2,3)] <- c("spec", "sens",'thres')
roc_RMMs['sum'] = roc_RMMs['sens']+roc_RMMs['spec']

PresenceTEST$predTSS<-ifelse(PresenceTEST$RMMprediction_FULL>roc_RMMs[which(roc_RMMs$sum == max(roc_RMMs$sum)),"thres"],1,0)
PresenceTEST$correctprediction<-ifelse(PresenceTEST$predTSS==PresenceTEST$RMM,1,0)

#roc_RMMs_D<-as.data.frame(roc_RMM_D$specificities)
#roc_RMMs_D$sensitivities<-roc_RMM_D$sensitivities
#roc_RMMs_D$thresholds<-roc_RMM_D$thresholds
#colnames(roc_RMMs_D)[c(1,2,3)] <- c("spec", "sens",'thres')
#roc_RMMs_D['sum'] = roc_RMMs_D['sens']+roc_RMMs_D['spec']

#PresenceTEST$predTSS_D<-ifelse(PresenceTEST$RMMprediction_DROPPED>roc_RMMs_D[which(roc_RMMs_D$sum == max(roc_RMMs_D$sum)),"thres"],1,0)
#PresenceTEST$correctprediction_D<-ifelse(PresenceTEST$predTSS_D==PresenceTEST$RMM,1,0)

presenceRMM<-subset(PresenceTEST,PresenceTEST$RMM==1)
absenceRMM<-subset(PresenceTEST,PresenceTEST$RMM==0)

tables[7,"TSS_threshold"]<- roc_RMMs[which(roc_RMMs$sum == max(roc_RMMs$sum)),'thres']
#tables[8,"TSS_threshold"]<- roc_RMMs[which(roc_RMMs_D$sum == max(roc_RMMs_D$sum)),'thres']

tables[7,"TSS_Ind"]<-(length(subset(presenceRMM$correctprediction,presenceRMM$correctprediction==1))/nrow(presenceRMM))+ (length(subset(absenceRMM$correctprediction,absenceRMM$correctprediction==1))/nrow(absenceRMM))-1
#tables[8,"TSS_Ind"]<-(length(subset(presenceRMM$correctprediction_D,presenceRMM$correctprediction_D==1))/nrow(presenceRMM))+(length(subset(absenceRMM$correctprediction_D,absenceRMM$correctprediction_D==1))/nrow(absenceRMM))-1

##-----------------------SAG------------------------
###---------Dependent data and validation-----------
#check model
SAG_mod
#SAG_mod_DROPPED

#check str
str(SAG_mod, max.level=2)
#str(SAG_mod_DROPPED, max.level=2)

# Access algorithm, tuning settings, and partition method information.
eval.algorithm(SAG_mod)
eval.tune.settings(SAG_mod) %>% head()
eval.partition.method(SAG_mod)

#eval.algorithm(SAG_mod_DROPPED)
#eval.tune.settings(SAG_mod_DROPPED) %>% head()
#eval.partition.method(SAG_mod_DROPPED)
# Results table with summary statistics for cross validation on test data.
eval.results(SAG_mod) %>% head()
#eval.results(SAG_mod_DROPPED) %>% head()
# Results table with cross validation statistics for each test partition.
eval.results.partitions(SAG_mod)
#eval.results.partitions(SAG_mod_DROPPED)

# List of models with names corresponding to tune.args column label.
eval.models(SAG_mod) %>% str(max.level = 1)
#eval.models(SAG_mod_DROPPED) %>% str(max.level = 1)

# to tune.args column label.
eval.predictions(SAG_mod)
#eval.predictions(SAG_mod_DROPPED)

# Original occurrence data coordinates with associated predictor variable values.
eval.occs(SAG_mod) %>% head()
#eval.occs(SAG_mod_DROPPED) %>% head()

# Background data coordinates with associated predictor variable values.
eval.bg(SAG_mod) %>% head()
#eval.bg(SAG_mod_DROPPED) %>% head()

# Partition group assignments for occurrence data.
eval.occs.grp(SAG_mod) %>% str()
#eval.occs.grp(SAG_mod_DROPPED) %>% str()

#evaluation plots
evalplot.stats(e = SAG_mod, stats = c("or.10p", "auc.val"), color = "fc", x.var = "rm", error.bars = FALSE)
#evalplot.stats(e = SAG_mod_DROPPED, stats = c("or.10p", "auc.val"), color = "fc", x.var = "rm", error.bars = FALSE)

# Overall results
res <- eval.results(SAG_mod)
#res_DROPPED <- eval.results(SAG_mod_DROPPED)

# This dplyr operation executes the sequential criteria explained above.
opt.seq <- res %>% 
  filter(or.10p.avg == min(or.10p.avg)) %>% 
  filter(auc.val.avg == max(auc.val.avg))
opt.seq
tables[9,1:19]<-opt.seq
tables[9,"species"]<-"SAG"

#opt.seq_D <- res_DROPPED %>% 
  #filter(or.10p.avg == min(or.10p.avg)) %>% 
  #filter(auc.val.avg == max(auc.val.avg))
#opt.seq_D
#tables[10,1:19]<-opt.seq_D
#tables[10,"species"]<-"SAG_D"

# We can select a single model from the ENMevaluation object using the tune.args of our
# optimal model.
mod.seq <- eval.models(SAG_mod)[[opt.seq$tune.args]]
#mod.seq_D <- eval.models(SAG_mod_DROPPED)[[opt.seq_D$tune.args]]

# Here are the non-zero coefficients in our model.
mod.seq$betas
#mod.seq_D$betas
# And these are the marginal response curves for the predictor variables with non-zero 
# coefficients in our model. We define the y-axis to be the cloglog transformation, which
# is an approximation of occurrence probability (with assumptions) bounded by 0 and 1
# (Phillips et al. 2017).

m1.mx <- eval.models(SAG_mod)[["fc.Q_rm.4"]]
m1.mx$betas
SAGmm<-m1.mx
# model object above.
predSAG <- eval.predictions(SAG_mod)[[opt.seq$tune.args]]
plot(predSAG)

writeRaster(predSAG, "E:/Chapter_3 Modelling/QGIS/predSAG.tif")
#predSAG_D <- eval.predictions(SAG_mod_DROPPED)[[opt.seq_D$tune.args]]
#plot(predSAG_D)

# give the testing data the CRS of the prediction
transTEST<-spTransform(SpatialTEST,CRSobj= proj4string(predSAG))

###---------------*Independent data validation------------
####----------------**Modelling predictions------------------

#Creating prediction
PresenceTEST$SAGprediction_FULL<-extract(predSAG,transTEST)
#PresenceTEST$SAGprediction_DROPPED<-extract(predSAG_D,transTEST)

#AUC
tables[9,c("Ind_auc_low","Ind_auc","Ind_auc_high")] <- ci.auc(roc(PresenceTEST,response=SAG,predictor=SAGprediction_FULL))
#tables[10,c("Ind_auc_low","Ind_auc","Ind_auc_high")] <- ci.auc(roc(PresenceTEST,response=SAG,predictor=SAGprediction_DROPPED))

#plain roc curve
roc_SAG<-roc(PresenceTEST,response=SAG,predictor=SAGprediction_FULL)
#roc_SAG_D<-roc(PresenceTEST,response=SAG,predictor=SAGprediction_DROPPED)

#Determine value for maximum accuraccy threshold
roc_SAGs<-as.data.frame(roc_SAG$specificities)
roc_SAGs$sensitivities<-roc_SAG$sensitivities
roc_SAGs$thresholds<-roc_SAG$thresholds
colnames(roc_SAGs)[c(1,2,3)] <- c("spec", "sens",'thres')
roc_SAGs['sum'] = roc_SAGs['sens']+roc_SAGs['spec']

PresenceTEST$predTSS<-ifelse(PresenceTEST$SAGprediction_FULL>roc_SAGs[which(roc_SAGs$sum == max(roc_SAGs$sum)),"thres"],1,0)
PresenceTEST$correctprediction<-ifelse(PresenceTEST$predTSS==PresenceTEST$SAG,1,0)

#roc_SAGs_D<-as.data.frame(roc_SAG_D$specificities)
#roc_SAGs_D$sensitivities<-roc_SAG_D$sensitivities
#roc_SAGs_D$thresholds<-roc_SAG_D$thresholds
#colnames(roc_SAGs_D)[c(1,2,3)] <- c("spec", "sens",'thres')
#roc_SAGs_D['sum'] = roc_SAGs_D['sens']+roc_SAGs_D['spec']

#PresenceTEST$predTSS_D<-ifelse(PresenceTEST$SAGprediction_DROPPED>roc_SAGs_D[which(roc_SAGs_D$sum == max(roc_SAGs_D$sum)),"thres"],1,0)
#PresenceTEST$correctprediction_D<-ifelse(PresenceTEST$predTSS_D==PresenceTEST$SAG,1,0)

presenceSAG<-subset(PresenceTEST,PresenceTEST$SAG==1)
absenceSAG<-subset(PresenceTEST,PresenceTEST$SAG==0)

tables[9,"TSS_threshold"]<- roc_SAGs[which(roc_SAGs$sum == max(roc_SAGs$sum)),'thres']
#tables[10,"TSS_threshold"]<- roc_SAGs[which(roc_SAGs_D$sum == max(roc_SAGs_D$sum)),'thres']

tables[9,"TSS_Ind"]<-(length(subset(presenceSAG$correctprediction,presenceSAG$correctprediction==1))/nrow(presenceSAG))+ (length(subset(absenceSAG$correctprediction,absenceSAG$correctprediction==1))/nrow(absenceSAG))-1
#tables[10,"TSS_Ind"]<-(length(subset(presenceSAG$correctprediction_D,presenceSAG$correctprediction_D==1))/nrow(presenceSAG))+(length(subset(absenceSAG$correctprediction_D,absenceSAG$correctprediction_D==1))/nrow(absenceSAG))-1

##-----------------------SCR------------------------
###---------Dependent data and validation-----------
#check model
SCR_mod
#SCR_mod_DROPPED

#check str
str(SCR_mod, max.level=2)
#str(SCR_mod_DROPPED, max.level=2)

# Access algorithm, tuning settings, and partition method information.
eval.algorithm(SCR_mod)
eval.tune.settings(SCR_mod) %>% head()
eval.partition.method(SCR_mod)

#eval.algorithm(SCR_mod_DROPPED)
#eval.tune.settings(SCR_mod_DROPPED) %>% head()
#eval.partition.method(SCR_mod_DROPPED)
# Results table with summary statistics for cross validation on test data.
eval.results(SCR_mod) %>% head()
#eval.results(SCR_mod_DROPPED) %>% head()
# Results table with cross validation statistics for each test partition.
eval.results.partitions(SCR_mod)
#eval.results.partitions(SCR_mod_DROPPED)

# List of models with names corresponding to tune.args column label.
eval.models(SCR_mod) %>% str(max.level = 1)
#eval.models(SCR_mod_DROPPED) %>% str(max.level = 1)

# to tune.args column label.
eval.predictions(SCR_mod)
#eval.predictions(SCR_mod_DROPPED)

# Original occurrence data coordinates with associated predictor variable values.
eval.occs(SCR_mod) %>% head()
#eval.occs(SCR_mod_DROPPED) %>% head()

# Background data coordinates with associated predictor variable values.
eval.bg(SCR_mod) %>% head()
#eval.bg(SCR_mod_DROPPED) %>% head()

# Partition group assignments for occurrence data.
eval.occs.grp(SCR_mod) %>% str()
#eval.occs.grp(SCR_mod_DROPPED) %>% str()

#evaluation plots
evalplot.stats(e = SCR_mod, stats = c("or.10p", "auc.val"), color = "fc", x.var = "rm", error.bars = FALSE)
#evalplot.stats(e = SCR_mod_DROPPED, stats = c("or.10p", "auc.val"), color = "fc", x.var = "rm", error.bars = FALSE)

# Overall results
res <- eval.results(SCR_mod)
#res_DROPPED <- eval.results(SCR_mod_DROPPED)

# This dplyr operation executes the sequential criteria explained above.
opt.seq <- res %>% 
  filter(or.10p.avg == min(or.10p.avg)) %>% 
  filter(auc.val.avg == max(auc.val.avg))
opt.seq
tables[11,1:19]<-opt.seq
tables[11,"species"]<-"SCR"

#opt.seq_D <- res_DROPPED %>% 
  #filter(or.10p.avg == min(or.10p.avg)) %>% 
  #filter(auc.val.avg == max(auc.val.avg))
#opt.seq_D
#tables[12,1:19]<-opt.seq_D
#tables[12,"species"]<-"SCR_D"

# We can select a single model from the ENMevaluation object using the tune.args of our
# optimal model.
mod.seq <- eval.models(SCR_mod)[[opt.seq$tune.args]]
#mod.seq_D <- eval.models(SCR_mod_DROPPED)[[opt.seq_D$tune.args]]

# Here are the non-zero coefficients in our model.
mod.seq$betas
#mod.seq_D$betas
# And these are the marginal response curves for the predictor variables with non-zero 
# coefficients in our model. We define the y-axis to be the cloglog transformation, which
# is an approximation of occurrence probability (with assumptions) bounded by 0 and 1
# (Phillips et al. 2017).
m1.mx <- eval.models(SCR_mod)[["fc.Q_rm.1"]]
m1.mx$betas
SCRmm<-m1.mx
# model object above.
predSCR <- eval.predictions(SCR_mod)[[opt.seq$tune.args]]
plot(predSCR)

writeRaster(predSCR, "E:/Chapter_3 Modelling/QGIS/predSCR.tif")

#predSCR_D <- eval.predictions(SCR_mod_DROPPED)[[opt.seq_D$tune.args]]
#plot(predSCR_D)

# give the testing data the CRS of the prediction
transTEST<-spTransform(SpatialTEST,CRSobj= proj4string(predSCR))

###---------------*Independent data validation------------
####----------------**Modelling predictions------------------

#Creating prediction
PresenceTEST$SCRprediction_FULL<-extract(predSCR,transTEST)
#PresenceTEST$SCRprediction_DROPPED<-extract(predSCR_D,transTEST)

#AUC
tables[11,c("Ind_auc_low","Ind_auc","Ind_auc_high")] <- ci.auc(roc(PresenceTEST,response=SCR,predictor=SCRprediction_FULL))
#tables[12,c("Ind_auc_low","Ind_auc","Ind_auc_high")] <- ci.auc(roc(PresenceTEST,response=SCR,predictor=SCRprediction_DROPPED))

#plain roc curve
roc_SCR<-roc(PresenceTEST,response=SCR,predictor=SCRprediction_FULL)
#roc_SCR_D<-roc(PresenceTEST,response=SCR,predictor=SCRprediction_DROPPED)

#Determine value for maximum accuraccy threshold
roc_SCRs<-as.data.frame(roc_SCR$specificities)
roc_SCRs$sensitivities<-roc_SCR$sensitivities
roc_SCRs$thresholds<-roc_SCR$thresholds
colnames(roc_SCRs)[c(1,2,3)] <- c("spec", "sens",'thres')
roc_SCRs['sum'] = roc_SCRs['sens']+roc_SCRs['spec']

PresenceTEST$predTSS<-ifelse(PresenceTEST$SCRprediction_FULL>roc_SCRs[which(roc_SCRs$sum == max(roc_SCRs$sum)),"thres"],1,0)
PresenceTEST$correctprediction<-ifelse(PresenceTEST$predTSS==PresenceTEST$SCR,1,0)

#roc_SCRs_D<-as.data.frame(roc_SCR_D$specificities)
#roc_SCRs_D$sensitivities<-roc_SCR_D$sensitivities
#roc_SCRs_D$thresholds<-roc_SCR_D$thresholds
#colnames(roc_SCRs_D)[c(1,2,3)] <- c("spec", "sens",'thres')
#roc_SCRs_D['sum'] = roc_SCRs_D['sens']+roc_SCRs_D['spec']

#PresenceTEST$predTSS_D<-ifelse(PresenceTEST$SCRprediction_DROPPED>roc_SCRs_D[which(roc_SCRs_D$sum == max(roc_SCRs_D$sum)),"thres"],1,0)
#PresenceTEST$correctprediction_D<-ifelse(PresenceTEST$predTSS_D==PresenceTEST$SCR,1,0)

presenceSCR<-subset(PresenceTEST,PresenceTEST$SCR==1)
absenceSCR<-subset(PresenceTEST,PresenceTEST$SCR==0)

tables[11,"TSS_threshold"]<- roc_SCRs[which(roc_SCRs$sum == max(roc_SCRs$sum)),'thres']
#tables[12,"TSS_threshold"]<- roc_SCRs[which(roc_SCRs_D$sum == max(roc_SCRs_D$sum)),'thres']

tables[11,"TSS_Ind"]<-(length(subset(presenceSCR$correctprediction,presenceSCR$correctprediction==1))/nrow(presenceSCR))+ (length(subset(absenceSCR$correctprediction,absenceSCR$correctprediction==1))/nrow(absenceSCR))-1
#[12,"TSS_Ind"]<-(length(subset(presenceSCR$correctprediction_D,presenceSCR$correctprediction_D==1))/nrow(presenceSCR))+(length(subset(absenceSCR$correctprediction_D,absenceSCR$correctprediction_D==1))/nrow(absenceSCR))-1

##-----------------------RMMlong------------------------
###---------Dependent data and validation-----------
#check model
RMMlong_mod
#RMMlong_mod_DROPPED

#check str
str(RMMlong_mod, max.level=2)
#str(RMMlong_mod_DROPPED, max.level=2)

# Access algorithm, tuning settings, and partition method information.
eval.algorithm(RMMlong_mod)
eval.tune.settings(RMMlong_mod) %>% head()
eval.partition.method(RMMlong_mod)

#eval.algorithm(RMMlong_mod_DROPPED)
#eval.tune.settings(RMMlong_mod_DROPPED) %>% head()
#eval.partition.method(RMMlong_mod_DROPPED)
# Results table with summary statistics for cross validation on test data.
eval.results(RMMlong_mod) %>% head()
#eval.results(RMMlong_mod_DROPPED) %>% head()
# Results table with cross validation statistics for each test partition.
eval.results.partitions(RMMlong_mod)
#eval.results.partitions(RMMlong_mod_DROPPED)

# List of models with names corresponding to tune.args column label.
eval.models(RMMlong_mod) %>% str(max.level = 1)
#eval.models(RMMlong_mod_DROPPED) %>% str(max.level = 1)

# to tune.args column label.
eval.predictions(RMMlong_mod)
#eval.predictions(RMMlong_mod_DROPPED)

# Original occurrence data coordinates with associated predictor variable values.
eval.occs(RMMlong_mod) %>% head()
#eval.occs(RMMlong_mod_DROPPED) %>% head()

# Background data coordinates with associated predictor variable values.
eval.bg(RMMlong_mod) %>% head()
#eval.bg(RMMlong_mod_DROPPED) %>% head()

# Partition group assignments for occurrence data.
eval.occs.grp(RMMlong_mod) %>% str()
#eval.occs.grp(RMMlong_mod_DROPPED) %>% str()

#evaluation plots
evalplot.stats(e = RMMlong_mod, stats = c("or.10p", "auc.val"), color = "fc", x.var = "rm", error.bars = FALSE)
#evalplot.stats(e = RMMlong_mod_DROPPED, stats = c("or.10p", "auc.val"), color = "fc", x.var = "rm", error.bars = FALSE)

# Overall results
res <- eval.results(RMMlong_mod)
#res_DROPPED <- eval.results(RMMlong_mod_DROPPED)

# This dplyr operation executes the sequential criteria explained above.
opt.seq <- res %>% 
  filter(or.10p.avg == min(or.10p.avg)) %>% 
  filter(auc.val.avg == max(auc.val.avg))
opt.seq
tables[13,1:19]<-opt.seq
tables[13,"species"]<-"RMMlong"

#opt.seq_D <- res_DROPPED %>% 
  #filter(or.10p.avg == min(or.10p.avg)) %>% 
  #filter(auc.val.avg == max(auc.val.avg))
#opt.seq_D
#tables[14,1:19]<-opt.seq_D
#tables[14,"species"]<-"RMMlong_D"

# We can select a single model from the ENMevaluation object using the tune.args of our
# optimal model.
mod.seq <- eval.models(RMMlong_mod)[[opt.seq$tune.args]]
#mod.seq_D <- eval.models(RMMlong_mod_DROPPED)[[opt.seq_D$tune.args]]

# Here are the non-zero coefficients in our model.
mod.seq$betas
#mod.seq_D$betas
# And these are the marginal response curves for the predictor variables with non-zero 
# coefficients in our model. We define the y-axis to be the cloglog transformation, which
# is an approximation of occurrence probability (with assumptions) bounded by 0 and 1
# (Phillips et al. 2017).

m1.mx <- eval.models(RMMlong_mod)[["fc.L_rm.5"]]
m1.mx$betas
RMMlongmm<-m1.mx
# model object above.
predRMMlong <- eval.predictions(RMMlong_mod)[[opt.seq$tune.args]]
plot(predRMMlong)

writeRaster(predRMMlong, "E:/Chapter_3 Modelling/QGIS/predRMMlong.tif")

#predRMMlong_D <- eval.predictions(RMMlong_mod_DROPPED)[[opt.seq_D$tune.args]]
#plot(predRMMlong_D)

# give the testing data the CRS of the prediction
transTEST<-spTransform(SpatialTEST,CRSobj= proj4string(predRMMlong))

###---------------*Independent data validation------------
####----------------**Modelling predictions------------------

#Creating prediction
PresenceTEST$RMMlongprediction_FULL<-extract(predRMMlong,transTEST)
#PresenceTEST$RMMlongprediction_DROPPED<-extract(predRMMlong_D,transTEST)

#AUC
tables[13,c("Ind_auc_low","Ind_auc","Ind_auc_high")] <- ci.auc(roc(PresenceTEST,response=RMM,predictor=RMMlongprediction_FULL))
#tables[14,c("Ind_auc_low","Ind_auc","Ind_auc_high")] <- ci.auc(roc(PresenceTEST,response=RMMlong,predictor=RMMlongprediction_DROPPED))

#plain roc curve
roc_RMMlong<-roc(PresenceTEST,response=RMM,predictor=RMMlongprediction_FULL)
#roc_RMMlong_D<-roc(PresenceTEST,response=RMMlong,predictor=RMMlongprediction_DROPPED)

#Determine value for maximum accuraccy threshold
roc_RMMlongs<-as.data.frame(roc_RMMlong$specificities)
roc_RMMlongs$sensitivities<-roc_RMMlong$sensitivities
roc_RMMlongs$thresholds<-roc_RMMlong$thresholds
colnames(roc_RMMlongs)[c(1,2,3)] <- c("spec", "sens",'thres')
roc_RMMlongs['sum'] = roc_RMMlongs['sens']+roc_RMMlongs['spec']

PresenceTEST$predTSS<-ifelse(PresenceTEST$RMMlongprediction_FULL>roc_RMMlongs[which(roc_RMMlongs$sum == max(roc_RMMlongs$sum)),"thres"],1,0)
PresenceTEST$correctprediction<-ifelse(PresenceTEST$predTSS==PresenceTEST$RMMlong,1,0)

#roc_RMMlongs_D<-as.data.frame(roc_RMMlong_D$specificities)
#roc_RMMlongs_D$sensitivities<-roc_RMMlong_D$sensitivities
#roc_RMMlongs_D$thresholds<-roc_RMMlong_D$thresholds
#colnames(roc_RMMlongs_D)[c(1,2,3)] <- c("spec", "sens",'thres')
#roc_RMMlongs_D['sum'] = roc_RMMlongs_D['sens']+roc_RMMlongs_D['spec']

#PresenceTEST$predTSS_D<-ifelse(PresenceTEST$RMMlongprediction_DROPPED>roc_RMMlongs_D[which(roc_RMMlongs_D$sum == max(roc_RMMlongs_D$sum)),"thres"],1,0)
#PresenceTEST$correctprediction_D<-ifelse(PresenceTEST$predTSS_D==PresenceTEST$RMMlong,1,0)

presenceRMMlong<-subset(PresenceTEST,PresenceTEST$RMMlong==1)
absenceRMMlong<-subset(PresenceTEST,PresenceTEST$RMMlong==0)

tables[13,"TSS_threshold"]<- roc_RMMlongs[which(roc_RMMlongs$sum == max(roc_RMMlongs$sum)),'thres']
#tables[14,"TSS_threshold"]<- roc_RMMlongs[which(roc_RMMlongs_D$sum == max(roc_RMMlongs_D$sum)),'thres']

tables[13,"TSS_Ind"]<-(length(subset(presenceRMMlong$correctprediction,presenceRMMlong$correctprediction==1))/nrow(presenceRMMlong))+ (length(subset(absenceRMMlong$correctprediction,absenceRMMlong$correctprediction==1))/nrow(absenceRMMlong))-1
#tables[14,"TSS_Ind"]<-(length(subset(presenceRMMlong$correctprediction_D,presenceRMMlong$correctprediction_D==1))/nrow(presenceRMMlong))+(length(subset(absenceRMMlong$correctprediction_D,absenceRMMlong$correctprediction_D==1))/nrow(absenceRMMlong))-1

#---------------Results table---------------
Results_table<-tables[c(1,3,5,7,9,11,13),c(4,8,9,12,20:25)]

#Results_table_Dropped<-tables[c(2,4,6,8,10,12),c(4,8,12,20:25)]

