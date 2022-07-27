#CDM calculations

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


##################CDM calculations#####################

#randomly sample raster stratified across the distribution say, in 0.10 quantiles
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
#soc_total <- quantile(r, probs=seq(0,1,1/10), na.rm=TRUE)
###issue with quantile results and how they integrate into sample stratified, see size = .1 below
#randomly sample within each quantile?

#randomly sample within each strata


#need stratified map to run CDM calculations
totalarea <- cellStats(area(clusters, units= "km", na.rm=TRUE), 'sum')


#strata 1####
#area
strata1area <- clusters %>% 
  clamp(1,1, useValues = FALSE) %>% 
  area(units= "km", na.rm=TRUE) %>% 
  cellStats('sum')

#distribution
clusters_s1 <- clusters
clusters_s1[clusters_s1 > 1] <- NA
clusters1<- mask(r, clusters_s1)
clusters1.df <- data.frame(clust)
hist(clusters1)
#bimodal

cluster1.sample <- sampleRandom(clusters1, 47 )
plot(clusters1)


#randomly sample within each strata area use set number -> use #'s from power analysis
#use those values to calculate below
?sample
#subsamples
mean_s1 <- mean(cluster1.sample)
sd_s1 <- sd(cluster1.sample)


#strata 2 ####
#area
strata2area <- clusters %>% 
  clamp(2,2, useValues = FALSE) %>% 
  area(units= "km", na.rm=TRUE) %>% 
  cellStats('sum')

#distribution
clusters_s2 <- clusters
clusters_s2[clusters_s2 != 2] <- NA
clusters2<- mask(r, clusters_s2)
hist(clusters2)
#skew right
cluster2.sample <- sampleRandom(clusters2, 30 )

#subsamples
mean_s2 <- mean(cluster2.sample)
sd_s2 <- sd(cluster2.sample)


#strata 3 ####
#area
strata3area <- clusters %>% 
  clamp(3,3, useValues = FALSE) %>% 
  area(units= "km", na.rm=TRUE) %>% 
  cellStats('sum')

#distribution
clusters_s3 <- clusters
clusters_s3[clusters_s3 != 3] <- NA
clusters3<- mask(r, clusters_s3)
hist(clusters3)
#bimodal
summary(clusters3)
cluster3.sample <- sampleRandom(clusters3, 12 )
plot(clusters3)

#subsamples
mean_s3 <- mean(cluster3.sample)
sd_s3 <- sd(cluster3.sample)

#strata 4 ####
#area
strata4area <- clusters %>% 
  clamp(4,4, useValues = FALSE) %>% 
  area(units= "km", na.rm=TRUE) %>% 
  cellStats('sum')

#distribution
clusters_s4 <- clusters
clusters_s4[clusters_s4 != 4] <- NA
clusters4<- mask(r, clusters_s4)
hist(clusters4)
#skew left

cluster4.sample <- sampleRandom(clusters4, 168 )

#subsamples
mean_s4 <- mean(cluster4.sample)
sd_s4 <- sd(cluster4.sample)


#create variables for CDM calculation#####

means <- c(mean_s1, mean_s2, mean_s3, mean_s4)
sds <- c(sd_s1, sd_s2, sd_s3, sd_s4)      
areas <- c(strata1area, strata2area, strata3area, strata4area)
area.prop <- areas/totalarea
data <- c(cluster1.sample, cluster2.sample, cluster3.sample, cluster4.sample)

plot(r)





####Estimation of Sample Plots using CDM equations#######

##CDM Eq. II - CALCULATION OF NUMBER OF SAMPLE PLOTS REQUIRED FOR ESTIMATION OF C STOCK WITHIN THE PROJECT BOUNDARY
#equation results in value < 1, use simplified equation below for projects with small sampling fraction (<5% project area)
#t value at 95% confidence interval (t) and a  
t = 1.960 #z score 95% confidence level

#second iteration, where t-val = n-1 degrees of freedom
#t = 3.182

sd_total <- sd(data)
n.total <- length(data)

sd.raster <- cellStats(r, sd)
mean.raster <- cellStats(r, mean)


#margin of error(p)
E = (t*(sd_total/sqrt(n.total)))
E

#equation II
n = (t^2* sum(area.prop*sds)^2) / 
  (E^2 + t^2 * sum(area.prop*sds^2*3))
n
#<1 results, use simplified equation


#simplified equation, small sampling fraction (<5% of project area)
n = (t/E)^2*(sum(area.prop*sds))^2
n
#4.45891 total samples


#########################