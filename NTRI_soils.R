
#working directory
setwd("C:/Users/lesley.atwood/Desktop/GitHub/NTRI")

#What this script does
#read in raster files
#plot raster
#conduct unspervised kMeans classification using satelitte imagery
#see: https://gis.stackexchange.com/questions/345823/identify-spatially-contiguous-clusters-in-raster-data-using-kmeans
#see: https://www.gis-blog.com/unsupervised-kmeans-classification-of-satellite-imagery-using-r/
#power analysis to determine # of samples

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


# import raster data (clipped using QGIS)
#https://soilgrids.org/

r <- raster("soc_t..ha_0..20cm_wgs84_clipped.tif")
#extent(r) #check size

#If needed, shrink map size
#drawExtent() #click on map area for extent
#e <- extent(806194.4, 1047181, 9436244, 9735652) #coordinates from drawExtent
#crop raster by the extent
#r_crop <- crop(r, e)
#plot(r_crop)
#r <- r_crop
plot(r)

###Stratify map by SOC using kmeans cluster ######
#convert raster to vector/matrix
nr <- getValues(r)

#set random number
set.seed(99)

######start fx #####

wss<-numeric(10)
for(i in 1:10){wss[i]<-sum(kmeans(na.omit(nr),centers=i)$withinss)}#sum within sum of squares and plot this

bss<-numeric(10)
for(i in 1:10){bss[i]<-sum(kmeans(na.omit(nr),centers=i)$betweenss)}#sum between sum of squares

plot(wss,type="b")
lines(bss, type = "b")

#where do the curves flatten?
c <- 4 # clusters

####apply clusters to raster####

i <- !is.na(nr)

#perform kmeans on matrix and inspect output, 4 clusters identified from kmeans loop
kmncluster <- kmeans(nr[i], centers = c)
nr[i] <- kmncluster$cluster
clusters <- setValues(r, nr)

par(mfrow = c(1,2))
plot(r, main = 'NTRI - SOC')

# add a color map with 5 colors
col=terrain.colors(4)

# add breaks to the colormap (6 breaks = 5 segments)
brk <- c(0,1,2,3,4)
plot(clusters, col= col, breaks = brk, main = 'Soil Strata Clusters')

par(mfrow=c(1,1))


totalarea <- cellStats(area(clusters, units= "km", na.rm=TRUE), 'sum')


#create four rasters where cells not in cluster = 0

strata1area <- clusters %>% 
                clamp(1,1, useValues = FALSE) %>% 
                area(units= "km", na.rm=TRUE) %>% 
                cellStats('sum')

strata2area <- clusters %>% 
              clamp(2,2, useValues = FALSE) %>% 
              area(units= "km", na.rm=TRUE) %>% 
              cellStats('sum')

strata3area <- clusters %>% 
                clamp(3,3, useValues = FALSE) %>% 
                area(units= "km", na.rm=TRUE) %>% 
                cellStats('sum')

strata4area <- clusters %>% 
                clamp(4,4, useValues = FALSE) %>% 
                area(units= "km", na.rm=TRUE) %>% 
                cellStats('sum')


###How many soil samples needed to detect change?#########################################
##############################Power Analysis##############################################



#reduce raster to polygons
#import shape file
polygon <- st_read("NTRI_SOC_CCROS/NTRI_SOC_CCROs.shp") #needed for plot, but not analysis


#project coordinate system of polygon to match raster
soc <- read_stars("soc_t..ha_0..20cm_wgs84_clipped.tif")


#ggplot() +
 # geom_stars(data = soc) +
  #geom_sf(data = polygon, alpha = 0)

#create df with cluster assignments
cluster_df <- as.data.frame(clusters, xy=TRUE) #dataframe with cluster assignment


soc_df <- as.data.frame(soc, xy=TRUE)

soc_clusters <- left_join(cluster_df, soc_df, by=c("x","y")) #takes time


#rename columns for analyses, manual input
soc_clusters <-rename(soc_clusters, "clusters" = "soc_t..ha_0..20cm_wgs84_clipped")
soc_clusters <-rename(soc_clusters, "socvalue" = "soc_t..ha_0..20cm_wgs84_clipped.tif")


#calculate soc means, std devs, total area of map, area per strata, percent of total area for each strata
means <- soc_clusters %>%
          group_by(clusters) %>%
          summarise(avg = mean(socvalue))
                    
                    
std <- soc_clusters %>%
          group_by(clusters) %>%
          summarise(stdev = sd(socvalue))

#means and stdevs for each soil strata (ided through kmeans cluster analysis)

avg <- means$avg
avg <- avg[complete.cases(avg)]
sd <- std$stdev
sd <- sd[complete.cases(sd)]

#Power analysis for each strata
#use soils revealed for intervention potential of area
#increase of 7 and 4.26/20 years to get ton/ha/yr
biggain <- 7/20
smallgain <- 4.26/20

s1_big <- power.t.test(n = NULL, delta = ((avg[1]+biggain)-avg[1]),
               sd = sd[1],
               sig.level = 0.05,
               power = 0.80,
               type = "one.sample")$n

s2_big <- power.t.test(n = NULL, delta = ((avg[2]+biggain)-avg[2]),
                       sd = sd[2],
                       sig.level = 0.05,
                       power = 0.80,
                       type = "one.sample")$n
s3_big <- power.t.test(n = NULL, delta = ((avg[3]+biggain)-avg[3]),
                       sd = sd[3],
                       sig.level = 0.05,
                       power = 0.80,
                       type = "one.sample")$n
s4_big <- power.t.test(n = NULL, delta = ((avg[4]+biggain)-avg[4]),
                       sd = sd[4],
                       sig.level = 0.05,
                       power = 0.80,
                       type = "one.sample")$n
n_biggains <- c(s1_big, s2_big, s3_big, s4_big)


###################small gain

s1_small <- power.t.test(n = NULL, delta = ((avg[1]+smallgain)-avg[1]),
                       sd = sd[1],
                       sig.level = 0.05,
                       power = 0.80,
                       type = "one.sample")$n

s2_small <- power.t.test(n = NULL, delta = ((avg[2]+smallgain)-avg[2]),
                       sd = sd[2],
                       sig.level = 0.05,
                       power = 0.80,
                       type = "one.sample")$n
s3_small <- power.t.test(n = NULL, delta = ((avg[3]+smallgain)-avg[3]),
                       sd = sd[3],
                       sig.level = 0.05,
                       power = 0.80,
                       type = "one.sample")$n
s4_small <- power.t.test(n = NULL, delta = ((avg[4]+smallgain)-avg[4]),
                       sd = sd[4],
                       sig.level = 0.05,
                       power = 0.80,
                       type = "one.sample")$n
n_smallgains <- c(s1_small, s2_small, s3_small, s4_small)

sum(n_biggains)
sum(n_smallgains)

###Results Table####

stratanumber <- c("1", "2", "3","4")
stratasoc_mean_tC_ha <- avg
stratasoc_sd_tC_ha <- sd


stratasize_km2 <- c(strata1area,strata2area,strata3area,strata4area)
totalsize_km2 <- c(totalarea, totalarea, totalarea, totalarea)
strataproportion <- stratasize/totalsize



results.df <- (data.frame(stratanumber,stratasoc_mean_tC_ha,
                         stratasoc_sd_tC_ha,stratasize_km2, 
                         totalsize_km2, strataproportion, 
                         n_biggains, n_smallgains))

write.csv(results.df, "NTRI_strata_results.csv")


