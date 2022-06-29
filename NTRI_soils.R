
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
extent(r) #check size

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
knr <- getValues(r)

i <- !is.na(knr)

#perform kmeans on matrix and inspect output, 4 clusters identified from kmeans loop
kmncluster <- kmeans(na.omit(nr[i]), centers = c)
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

###How many soil samples needed to detect change?#########################################
##############################Power Analysis##############################################



#reduce raster to polygons
#import shape file
polygon <- st_read("NTRI_SOC_CCROS/NTRI_SOC_CCROs.shp") #needed for plot, but not analysis


#project coordinate system of polygon to match raster
soc <- read_stars("soc_t..ha_0..20cm_wgs84_clipped.tif")


####begin next function#######
ggplot() +
  geom_stars(data = soc) +
  geom_sf(data = polygon, alpha = 0)

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

      
#area of entire site (km2)
totalarea_km2 <- function(filename) {
  cell_size <- area(r, na.rm=TRUE, weights=FALSE)      
  #delete NAs from vector of all raster cells
  ##NAs lie outside of the rastered region, can thus be omitted
  cell_size<-cell_size[!is.na(cell_size)]
        #compute area [km2] of all cells in geo_raster
        raster_area<-length(cell_size)*median(cell_size)
        print(raster_area)
}

NTRI_km2 <-totalarea_km2(r)

#area of each strata (km2)

####rework here##############################
#strata area size
strata1 <- soc_clusters[clusters==1]

strataarea_km2 <- soc_clusters %>%
  group_by(clusters) %>%
  summarise(avg = mean(socvalue))

strataarea_km2 <- function(filename) {
  starta_area <- r %>% 
    group_by(clusters)
  cell_size <- area(strata_area, na.rm=TRUE, weights=FALSE)      
#delete NAs from vector of all raster cells
##NAs lie outside of the rastered region, can thus be omitted
  cell_size<-cell_size[!is.na(cell_size)]
  length(cell_size)*median(cell_size)
  print(strata_area)
}

strataarea_km2(soc_clusters)

area_strata_km2 <- function(filename) {
  soc_clusters %>%
    group_by(clusters) %>%
    area_km2()
}






NTRI_strata_mean_sd <- left_join(means, std)
write.csv(NTRI_strata_mean_sd,"NTRI_strata_mean_sd.csv", row.names = FALSE)

#means and stdevs for each soil strata (ided through kmeans cluster analysis)

#(second map, africa specific map)     
#cluster   avg stdev
#1         11.5   0.568
#2         16.7   1.61 
#3         9.54  0.843
#4         13.6   0.720
c1 <- 11.5
c2 <-16.7
c3 <- 9.54
c4 <- 13.6

delta <- c(c1, c2, c3, c4)

c1_sd <- 0.568
c2_sd <-1.61
c3_sd <- 0.843
c4_sd <- 0.72

sd <- c(c1_sd, c2_sd, c3_sd, c4_sd)

#first map (global map)
 #cluster   avg stdev
#1           1  44.3  2.24
#2           2  31.9  2.80
#3           3  55.8  6.32
#4           4  38.2  1.62

c1 <- 44.3
c2 <-31.9
c3 <- 55.8
c4 <- 38.2

delta <- c(c1, c2, c3, c4)

c1_sd <- 2.24
c2_sd <-2.8
c3_sd <- 6.32
c4_sd <- 1.32

sd <- c(c1_sd, c2_sd, c3_sd, c4_sd)

#Power analysis for each strata
#use soils revealed for intervention potential of area
#increase of 7 and 4.26/20 years to get ton/ha/yr
biggain <- 7/20
smallgain <- 4.26/20

power.t.test(n = NULL, delta = ((c1+biggain)-c1),
               sd = c1_sd,
               sig.level = 0.05,
               power = 0.80,
               type = "one.sample")
c1big_n = 22.66698 #323.4161

power.t.test(n = NULL, delta = ((c2+biggain)-c2),
             sd = c2_sd,
             sig.level = 0.05,
             power = 0.80,
             type = "one.sample")
c2big_n = 168.0131 #504.2524

power.t.test(n = NULL, delta = ((c3+biggain)-c3),
             sd = c3_sd,
             sig.level = 0.05,
             power = 0.80,
             type = "one.sample")
c3big_n = 47.48959 #2561.13

power.t.test(n = NULL, delta = ((c4+biggain)-c4),
             sd = c4_sd,
             sig.level = 0.05,
             power = 0.80,
             type = "one.sample")
c4big_n = 35.18424 #113.5756

biggain_n <- c(c1big_n, c2big_n, c3big_n, c4big_n)

###################small gain

power.t.test(n = NULL, delta = ((c1+smallgain)-c1),
             sd = c1_sd,
             sig.level = 0.05,
             power = 0.80,
             type = "one.sample")
c1small_n = 57.76444 #869.9722

power.t.test(n = NULL, delta = ((c2+smallgain)-c2),
             sd = c2_sd,
             sig.level = 0.05,
             power = 0.80,
             type = "one.sample")
c2small_n = 450.3602 #1358.249

power.t.test(n = NULL, delta = ((c3+smallgain)-c3),
             sd = c3_sd,
             sig.level = 0.05,
             power = 0.80,
             type = "one.sample")
c3small_n = 124.8773 #6911.994

power.t.test(n = NULL, delta = ((c4+smallgain)-c4),
             sd = c4_sd,
             sig.level = 0.05,
             power = 0.80,
             type = "one.sample")
c4small_n = 91.62298 #303.3632

smallgain_n <- c(c1small_n, c2small_n, c3small_n, c4small_n)

strata_sample <- data.frame(biggain_n, smallgain_n)
strata_sample

#Africa specific map
#biggain_n smallgain_n
#1  22.66698    57.76444
#2 168.01310   450.36020
#3  47.48959   124.87730
#4  35.18424    91.62298

 sum(biggain_n) 
 #273.3539 total samples
 sum(smallgain_n)
 #724.6249 total samples

#   biggain_n smallgain_n
#1  323.4161  869.9722
#2  504.2524 1358.2490
#3 2561.1300 6911.9940
#4  113.5756  303.3632



power.t.test(n=NULL, delta = ((k1_mean+biggain)-k1_mean), #what is the expected change in SOC? What are the units for the base map?
             sd = k1_sd, sig.level = 0.05, power = 0.8, 
             type = "one.sample")
#n = 323.4161



 
#means and sd across all NTRI sites
SOCmean <- cellStats(r, stat = "mean", na.rm=TRUE)
SOCsd <- cellStats(r, stat = "sd", na.rm=TRUE)
cbind(SOCmean, SOCsd)

#  SOCmean   SOCsd
# 36.72706 5.467776 

################Use above for power analysis###################


#parameters needed
#effect size (0.1 is generally a small effect size, but maybe use data to inform this value)
#significance level (0.05 ?)
#power of the test (.8)




#one sample t test
power.t.test(n=NULL, delta = ((SOCmean+biggain)-SOCmean), #what is the expected change in SOC? What are the units for the base map?
             sd = SOCsd, sig.level = 0.05, power = 0.8, 
             type = "one.sample")
#n = 1917

#one sample t test
power.t.test(n=NULL, delta = ((SOCmean+smallgain)-SOCmean), #what is the expected change in SOC? What are the units for the base map?
             sd = SOCsd, sig.level = 0.05, power = 0.8, 
             type = "one.sample")
#n = 5174




#power.t.test(n=NULL, delta = (40-SOCmean), #what is the expected change in SOC? What are the units for the base map?
#             sd = SOCsd, sig.level = 0.05, power = 0.8, 
#             type = "two.sample")
#n = 44.79181

#power.t.test(n=NULL, delta = (40-SOCmean), #what is the expected change in SOC? What are the units for the base map?
#             sd = SOCsd, sig.level = 0.05, power = 0.8, 
#             type = "paired")
#n = 23.89726


#########################################
pwr.anova.test(k = 4, f = .1, sig.level = 0.5, power = .8)
#Balanced one-way analysis of variance power calculation 

#k = 4
#n = 74.37826, # samples for each group
#f = 0.1
#sig.level = 0.5
#power = 0.8

#now calculate based on acutual std devations of SOC values using raster data
#https://www.sheffield.ac.uk/polopoly_fs/1.885243!/file/117_Power_Analysis.pdf


##power analysis based on r raster file

##effect size calculation
#(mean of treatment - mean of control) / stdev(control)

#"Potential of grassland rehabilitation through high density-short duration grazing to sequester atmospheric carbon"

trt_mean <- 12.36 #g C/m^2/y where treatment is high density SD
control_mean <- 3.94 #g C/m^2/y where control is no management

trt_se <- 2.12
control_se <- 0.75

control_sd <- control_se * sqrt(12)


e_size <- (trt_mean-control_mean)/control_sd

#pwr.anova.test(k = 4, f = e_size, sig.level = 0.5, power = .8)

pwr.t.test(n = NULL, delta = 1.37, sd = 2.598076, sig.level = 0.05, 
           power = 0.8, 
           type = c("two.sample"))


library(pwr2)
n <- seq(2, 30, by=4)
f <- seq(0.1, 1.0, length.out=10)
pwr.plot(n=n, k=5, f=f, alpha=0.05)
