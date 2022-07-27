#####Creat Shapefile with random distributed sample points within each strata####

#working directory
setwd("C:/Users/lesley.atwood/Desktop/GitHub/NTRI")

#### LOAD LIBRARIES ####
library(rgdal)
library(raster) #reads raster shapefiles
library(sp)
library(sf)
library(pwr)
library(stars)
library(ggplot2)
library(stats)
library(dplyr)
library(foster)
library(terra)

# import raster data (clipped using QGIS) ##########
#https://soilgrids.org/

r <- raster("soc_t..ha_0..20cm_wgs84_clipped.tif")
clusters <- raster("clusters.tif")


plot(r)
plot(clusters)


#import power analysis results
power <- read.csv("NTRI_strata_results.csv")

#total number of samples calcuated from power analysis
big.sample.total <- round(sum(power$n.power_0.35tC.ha.yr))

#proportion of samples allocated to each strata based on power analysis
prop.sample.strata <- power$n.power_0.35tC.ha.yr/big.sample.total

#use proportion of samples from power analysis to determine % of 330 samples randomly distributed to each strata

#strata 1####
clusters_s1 <- clusters
clusters_s1[clusters_s1 > 1] <- NA
clusters1<- mask(r, clusters_s1)

#identify sampling points
n.s1 <- round(329*prop.sample.strata[1])
cluster1.sample <- sampleRandom(clusters1, size = n.s1, sp=TRUE )
plot(cluster1.sample)

summary(cluster1.sample)


#strata 2 ####
clusters_s2 <- clusters
clusters_s2[clusters_s2 != 2] <- NA
clusters2<- mask(r, clusters_s2)
plot(clusters2)


#identify sampling points
n.s2 <- round(329*prop.sample.strata[2])
cluster2.sample <- sampleRandom(clusters2, size = n.s2, sp=TRUE)
plot(cluster2.sample)

summary(cluster2.sample)


#strata 3 ####
clusters_s3 <- clusters
clusters_s3[clusters_s3 != 3] <- NA
clusters3<- mask(r, clusters_s3)

#identify sampling points
n.s3 <- round(329*prop.sample.strata[3])
cluster3.sample <- sampleRandom(clusters3, size = n.s3, sp=TRUE)
plot(cluster3.sample)

summary(cluster3.sample)


#strata 4 ####
clusters_s4 <- clusters
clusters_s4[clusters_s4 != 4] <- NA
clusters4<- mask(r, clusters_s4)

#identify sampling points
n.s4 <- round(329*prop.sample.strata[4])+2
cluster4.sample <- sampleRandom(clusters4, size = n.s4, sp=TRUE)
plot(cluster4.sample)

summary(cluster4.sample)

####sum sample numbers to double check 330 in total
sum(c(n.s1, n.s2, n.s3, n.s4))


#stitch sp files together into one shapefile###

points <- rbind(cluster1.sample, cluster2.sample, cluster3.sample, cluster4.sample)
plot(points)

writeOGR(points, dsn ="NTRI_soilsamplelocations", layer="NTRI_soilsamplelocations",
         driver = "ESRI Shapefile" )

#ensure shapefile contains all 330 points
samplepoints <- st_read("NTRI_soilsamplelocations/NTRI_soilsamplelocations.shp")
plot(samplepoints)
summary(samplepoints)
