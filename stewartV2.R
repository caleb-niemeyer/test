
install.packages(c("RSQLite","tidyr","rgdal","foreach","doParallel","parallel","reshape2","GISTools", "ggplot2", "stringr", "spdep", "lubridate","gstat","dplyr","plyr","sf","scales","nngeo","measurements", "stars","starsExtra","tidyverse","raster","ggspatial","FNN"))


library(lubridate)
library(gdalcubes)
library(ggthemes)
library(reshape2)
library(dplyr)
library(plyr)
library(spdep)
library(sf)
library(measurements)
library(nngeo)
library(scales)
library(stars)
library(starsExtra)
library(tidyverse)
library(gstat)
library(sp)
library(raster)
library(FNN)
library(broom)
library(ggpubr)
library(tidymodels)
library(exactextractr)
library(dplyr)
library(plyr)
library(reshape2)
library(lemon)
library(rstatix)
library(cowplot)
library(esquisse)




Mode <- function(x) {
  if ( length(x) <= 2 ) return(x[1])
  if ( anyNA(x) ) x = x[!is.na(x)]
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

scaleFUN <- function(x) sprintf("%.2f", x)

cbp3a <- c("#56B4E9", "#009E73","#E69F00") # color blind palette
cbp3b <- c("#1b9e77","#d95f02","#7570b3")


setwd("C:\\Users\\cniemeyer\\OneDrive - University of Guelph\\Documents\\R\\Stewart")


tmpdir<-paste("MAIN\\tmp")
DATA_OUT <- paste("MAIN\\OUT")

path.gpkg<-"C:\\Users\\cniemeyer\\OneDrive - University of Guelph\\Documents\\R\\Stewart\\"

st_layers(paste(path.gpkg,"stewart.gpkg",sep=""))#list layers

CI_2class<-raster("CI_2classp1.tif") #convergence index :  concave convex ; you don't need the entire path to the raster because we already set the working directory

names(CI_2class)<-"CI"
spplot(CI_2class)

CI_2class<-as.factor(CI_2class)

SoilType<-raster("stewartDC.tif")
plot(SoilType)

m10<-st_read(paste(path.gpkg,"stewart.gpkg",sep=""), layer ="NStrips10" ,fid_column_name="fid")  # read 10 meter strip polygons
m20<-st_read(paste(path.gpkg,"stewart.gpkg",sep=""), layer ="NStrips20" ,fid_column_name="fid")
m30<-st_read(paste(path.gpkg,"stewart.gpkg",sep=""), layer ="NStrips30" ,fid_column_name="fid")
m40<-st_read(paste(path.gpkg,"stewart.gpkg",sep=""), layer ="NStrips40" ,fid_column_name="fid")

ObjectsVector <- c("m10", "m20", "m30","m40") # group NStrip objects together
list.strips<- lapply(ObjectsVector, get) # make list of NStrip objects
names(list.strips)<-c("10", "20", "30","40")


yield<-st_read(paste(path.gpkg,"stewart.gpkg",sep=""), layer ="Yield" ) # read yield data
yield<-st_intersection(m40,yield) # intersect yield with strip polygons (clipping)


yield.sp<-as(yield,"Spatial") # change yield from sf (newer format) to sp (older format)
#alpha angle at which this variogram will run below (find this in qgis) need to run for my strips
vg<-variogram(Yield~1,yield.sp,width=5,cutoff=100, alpha = c(0, 45, 90, 135)) # you need the direction the strips are oriented

ggplot(vg,aes(x=dist,y=gamma))+geom_point()+facet_wrap(~dir.hor)+ 
  scale_x_continuous(breaks=seq(0,100,20))+
  ylab("semivariance")+
  xlab("distance (m)")

yield.sp<-as(yield,"Spatial")


BigSF <- NULL # create empty container for loop output

for (i in 1:length(list.strips))
{ i<-2
  plots.sp<-as(list.strips[[i]],"Spatial")
  yieldbu<-over(plots.sp, yield.sp[,"Yield"], fn = mean) # average yield into plots
  plots.sp$Yield<-yieldbu$Yield
  plots.sf<-st_as_sf(plots.sp)
  Soils <- exactextractr::exact_extract(SoilType, plots.sf, 'mode',append_cols = c("fid"))%>% 
    mutate(SoilType=recode(mode, 
                           '1'="Well",
                           '2'="Imperfect",
                           '3'="Poor"))%>%
    select(-mode) # remove the column named "mode"}
  
  CI<-exactextractr::exact_extract(CI_2class, plots.sf, 'mode',append_cols = c("fid"))%>% # we need to add the column name "fid" to the result so we can join it back to the plot polygons
    mutate(CI=recode(mode,  # the lookup table fits better here
                     '1'="Concave",
                     '2'="Convex")) %>%
    select(-mode) 
  
  plots.sf<-st_as_sf(plots.sp)
  #join by attributes
  plots.sf<-merge(plots.sf,CI,by.x="fid") # combine zonal stats output with vector file
  plots.sf<-merge(plots.sf,Soils,by.x="fid")
  
  plots.sf <- plots.sf %>% 
    mutate(fungicide=recode(fid_2, 
                          `1`="trt",
                          `2`="notrt")) %>% 
    mutate(treatment=fid_2)
  
  plots.trt<-plots.sf %>% filter(treatment ==1) # sf subset
  plots.notrt<-plots.sf %>% filter(treatment ==2) # sf subset
  
  st_agr(plots.trt) = "constant" # prevent a warning
  cent.trt<-st_centroid(plots.trt)
  
  st_agr(plots.notrt) = "constant" # prevent a warning
  cent.notrt<-st_centroid(plots.notrt)  
  
  
  st_agr(plots.sf) = "constant" # prevent a warning
  plot.cents<-st_centroid(plots.sf) #centroid of polygon 
  # join nearest neighbour
  p.join.NN<-st_join(cent.trt,cent.notrt[c("fid","SoilType","CI","Yield")],join = st_nearest_feature)  
  # add plot length as a factor level
  p.join.NN$PlotLength<-names(list.strips)[i] #add plot length as factor level
  
  p.join.NN$Aggregation<-"Adjacent"  #add curvature class as factor level

  uniCI<-unique(plot.cents$CI) # find unique values of factor levels (e.g. convex/concave, soil types)
  p.join.CI <- NULL

  for (j in 1:length(uniCI)){
    plots.match<-plot.cents %>% filter(CI== uniCI[j])
    plots.trt<-plots.match %>% filter(treatment ==1)
    plots.notrt<-plots.match %>% filter(treatment ==2) # sf subset
    result<-st_join(plots.trt,plots.notrt[c("fid","SoilType","CI","Yield")],join = st_nearest_feature)
    result$PlotLength<-names(list.strips)[i] #add plot length as factor level
    result$Aggregation<-"CI" #add CI as factor level
    p.join.CI=rbind(p.join.CI,result)
  }
  
  uniSoil<-unique(plot.cents$SoilType)
  p.join.Soil <- NULL
  
  
  for (k in 1:length(uniSoil)){
    plots.match<-plot.cents %>% filter(SoilType == uniSoil[k])
    plots.trt<-plots.match %>% filter(treatment ==1)
    plots.notrt<-plots.match %>% filter(treatment ==2) # sf subset
    result<-st_join(plots.trt,plots.notrt[c("fid","SoilType","CI","Yield")],join = st_nearest_feature)
    result$PlotLength<-names(list.strips)[i] #add plot length as factor level
    result$Aggregation<-"SoilType" #add soil type as factor level
    p.join.Soil=rbind(p.join.Soil,result)
  }
  
  PlotLengthIterate<- rbind(rbind(p.join.NN, p.join.CI),p.join.Soil)
  
  BigSF<-rbind(BigSF,PlotLengthIterate)
  
} 
  st_write(BigSF, "stewart_SFoutput.gpkg", driver="GPKG")
  
  st_layers(paste(path.gpkg,"stewart_SFoutput.gpkg",sep=""))
  
  yield<-st_read(paste(path.gpkg,"stewart.gpkg",sep=""), layer ="yield" )
  
  d <- st_read(paste(path.gpkg,"stewart_SFoutput.gpkg", sep=""), layer ="stewart_SFoutput")

  d <- d %>%
    filter(PlotLength == "20" & Aggregation == "SoilType") 
  
  d$d.Yield<-d$Yield.x-d$Yield.y
  st_write(d, "stewart_SFoutputDY.gpkg", driver="GPKG")
  
  esquisser(d)
  
  ggplot(d) +
    aes(x = CI.y, y = d.Yield) +
    geom_boxplot(fill = "#112446") +
    theme_minimal()
  
    ggplot(d) +
    aes(x = SoilType.x, y = d.Yield) +
    geom_boxplot(fill = "#112446") +
    theme_minimal()

# Anova
    oneway <- aov(d.Yield ~ SoilType.x, data = d)
    summary(oneway)

 # pairwise.test <- d %>% 
 # anova_test(d.Yield ~ SoilType)
 # adjust_pvalue() %>%
 # add_significance("p.adj")
 # pairwise.test