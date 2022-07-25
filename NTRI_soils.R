
#working directory
setwd("C:/Users/lesley.atwood/Desktop/GitHub/NTRI")

#What this script does
#read in raster files
#plot raster
#conducts unspurvised kMeans classification using satellite imagery
#power analysis to determine # of samples
#CDM eq. II to determine # of samples


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

######kmeans plots to determine cluster number #####

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

#Sample size: for big and small gains based on soils revealed data#########
#increase of 7 and 4.26/20 years to get ton/ha/yr
biggain <- 7/20
smallgain <- 4.26/20

#####Big Gain SOC###############
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
numsamples_0.35t_ha_yr <- c(s1_big, s2_big, s3_big, s4_big)


###################small gain SOC#####

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
numsamples_0.213t_ha_yr <- c(s1_small, s2_small, s3_small, s4_small)

sum(numsamples_0.35t_ha_yr)
sum(numsamples_0.213t_ha_yr)

###Results Table####
results.df <- matrix(NA, 4,8)
results.df[,1] <- 1:4
results.df[,2] <- avg
results.df[,3] <- sd
results.df[,4] <- c(strata1area,strata2area,strata3area,strata4area)
results.df[,5] <- totalarea
results.df[,6] <- results.df[,4]/results.df[,5] 
results.df[,7] <- numsamples_0.35t_ha_yr
results.df[,8] <- numsamples_0.213t_ha_yr
results.df <- data.frame(results.df)

#rename columns#######


results.df <- results.df %>%
                rename( strata = X1,
                        mean_SOC.tC.ha = X2,
                        sd_SOC.tC.ha = X3,
                        area_km2 = X4,
                        area_total_km2 = X5,
                        area_proportion = X6,
                        n.power_0.35tC.ha.yr = X7,
                        n.power_0.213tC.ha.yr = X8)

write.csv(results.df, "NTRI_strata_results.csv", row.names = FALSE)



##################CDM calculations#####################




#randomly sample your raster stratified across the distribution say, in 0.10 quantiles
#use the subsample to calculate your CI’s
#Note that if your rasters distribution is Gaussian, there is no need to stratify as a random sample should capture the sample distribution.
#To ensure that there is nothing hinky going on you could Bootstrap this procedure and output the p-values along with a distribution of the upper and lower CI’s
#this way you can evaluate if there is a sampling condition that is increasing the variance of the CI’s.
#Although, this is not entirely necessary if you just check to see if the sample distribution of the subsample is the same as your population (raster).
#I can almost guarantee that, even with a small sample fraction and using an equivalency test (eg., Kolmogorov–Smirnov), the sample distribution will match.

#Randomly sample each strata in 0.10 quantiles

#check distribution of data
hist(r)
#Gaussian


####Calculate CI for entire area (r) using quantiles
    soc_total <- quantile(r, probs=seq(0,1,1/10), na.rm=TRUE)
    ###issue with quantile results and how they integrate into sample stratified, see size = .1 below
    
    #variables needed for CDM calculation#####
    mean.total <- mean(soc_total)
    sd.total <- sd(soc_total)
    data.total <- unname(soc_total)
    

#strata 1####
      #distribution
      clusters_s1 <- clusters
      clusters_s1[clusters_s1 > 1] <- NA
      clusters1<- mask(r, clusters_s1)
      hist(clusters1)
      #biomodal
      
      #subsamples
      soc_s1 <- quantile(clusters1, probs=seq(0,1,1/10), na.rm=TRUE)
      mean_s1 <- mean(n.sample_s1)
      sd_s1 <- sd(n.sample_s1)
      
      
#strata 2 ####
      #distribution
      clusters_s2 <- clusters
      clusters_s2[clusters_s2 != 2] <- NA
      clusters2<- mask(r, clusters_s2)
      hist(clusters2)
      #skew right
      
      #subsamples
      soc_s2 <- quantile(clusters2, probs=seq(0,1,1/10), na.rm=TRUE)
      mean_s2 <- mean(n.sample_s2)
      sd_s2 <- sd(n.sample_s2)


#strata 3 ####
      #distribution
      clusters_s3 <- clusters
      clusters_s3[clusters_s3 != 3] <- NA
      clusters3<- mask(r, clusters_s3)
      hist(clusters3)
      #bimodal
      
      #subsamples
      soc_s3 <- quantile(clusters3, probs=seq(0,1,1/10), na.rm=TRUE)
      mean_s3 <- mean(n.sample_s3)
      sd_s3 <- sd(n.sample_s3)

#strata 4 ####
      #distribution
      clusters_s4 <- clusters
      clusters_s4[clusters_s4 != 4] <- NA
      clusters4<- mask(r, clusters_s4)
      hist(clusters4)
      #skew left
      
      #subsamples
      soc_s4 <- quantile(clusters4, probs=seq(0,1,1/10), na.rm=TRUE)
      mean_s4 <- mean(n.sample_s4)
      sd_s4 <- sd(n.sample_s4)

      
#variables needed for CDM calculation#####
means <- c(mean_s1, mean_s2, mean_s3, mean_s4)
sds <- c(sd_s1, sd_s2, sd_s3, sd_s4)      
areas <- c(strata1area, strata2area, strata3area, strata4area)
area.prop <- areas/totalarea
data <- c(unname(soc_s1), unname(soc_s2), unname(soc_s3), unname(soc_s4))




####Estimation of Sample Plots using CDM equations#######

##CDM Eq. II - CALCULATION OF NUMBER OF SAMPLE PLOTS REQUIRED FOR ESTIMATION OF C STOCK WITHIN THE PROJECT BOUNDARY
#equation results in value < 1, use simplified equation below for projects with small sampling fraction (<5% project area)
#t value at 95% confidence interval (t) and a  
t = 1.960

sd_total <- sd(data)
n.total <- length(data)

#margin of error(p)
E = (t*(sd_total/sqrt(n.total)))


#equation II
n = (t^2* sum(area.prop*sds)^2) / 
 (E^2 + t^2 * sum(area.prop*sds^2))
#<1 results, use simplified equation


#simplified equation, small sampling fraction (<5% of project area)
  n = (t/E)^2*(sum(area.prop*sds))^2
  #4.45891 total samples


#########################