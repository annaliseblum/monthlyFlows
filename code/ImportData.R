##Data import: USGS daily flows from NWIS, climate from Maurer et al and basin characteristics from GAGESII
#http://www.engr.scu.edu/~emaurer/gridded_obs/index_gridded_obs.html
#New Everything
#May 5, 2016
##Annalise blum annaliseblum@gmail.com

library(dataRetrieval) #to pull data from USGS site
library(geoR); library(fossil) #to do nearest neighbor stuff

#### 1 - USGS daily flows ####
##GAGES II REf sites for North East
NEsitedf<-read.table("data/NErefGAGESIIsites.txt", colClasses="character")
NEsites<-NEsitedf$V1

#import flow data from USGS site, for NEsites from 1950-2010
parameterCd<-"00060" #set the parameter of interest at cfs
rawDailyData<-readNWISdv(NEsites,parameterCd, startDate = "1950-01-01", endDate = "2010-12-31")
NEdailyRAW<-rawDailyData

#pull just the relevant columns
NEdailyRAW<-data.frame(NEdailyRAW$site_no,NEdailyRAW$Date,NEdailyRAW$X_00060_00003)
names(NEdailyRAW)<-c("site_no","Date","Flow")
save(NEdailyRAW, file="data/NEdailyRAW.rdata") #save

#pull sites with at least 10 years of data
#create variables for day, month year
df<-NEdailyRAW #rename shorter
df$year<-as.numeric(format(df$Date, "%Y")) #extract year variable 
df$month<-as.numeric(format(df$Date, "%m")) #extract month variable 
df$day<-as.numeric(format(df$Date, "%d")) #extract day variable

##2.2 remove years without all 365 days##
YearTally <- aggregate(df$Flow,by=list(df$site_no, df$year),FUN=length) #need to drop years without 365 days of data
names(YearTally)<-c("site_no","year","daysperyr")
sum(YearTally$daysperyr<365)
df1<-data.frame(merge(df, YearTally, by = c('site_no','year'))) #merge with flow data
df2<-df1[df1$daysperyr>364,] #remove years with less than 365 days of data

#now need to check years of data for each site and remove sites with <10 years of data
RLcheck <- aggregate(df2$Flow,by=list(df2$site_no),FUN=length) #need to drop years without 365 days of data
names(RLcheck)<-c("site_no","totaldays")
sum(RLcheck$totaldays<10*365) #19
df3<-data.frame(merge(df2, RLcheck, by = 'site_no')) #merge with flow data
df4<-df3[df3$totaldays>=3650,]

Sites<-unique(df4$site_no)
#with new list of sites, pull lat and long
SiteInfo<-readNWISsite(Sites) #get site information for the sites
SiteInfoKeep<-SiteInfo[c("site_no","station_nm","dec_lat_va","dec_long_va","alt_va","drain_area_va","huc_cd")] #keep relevant vars
SiteInfoKeep$site_no<-as.character(SiteInfoKeep$site_no)

#### 2 - Weather from Maurer et al ####
#Monthly precip, max temp, min temp and wind speed
# The links below are to files containing daily precipitation (mm/day), maximum and minimum temperature (C), 
# and 10-m wind speed (m/s) for each 1/8-degree grid cell, grouped (using the UNIX tar command) by hydrologic area. 
# The individual data files, (which can be extracted using tar zxvf <tar filename>) indicate the location center of 
# the grid cell in the file name: data_<latitude>_<longitude>. Each daily data file contains columns of year, month, 
# day, daily total precipitation, maximum temperature, minimum temperature, and average 10-meter wind speed in ascii 
# format. If you need sub-daily interpolation of temperature, or radiative and humidity forcings as well as those listed 
# above, the current VIC model code can be run with the OUTPUT_FORCE flag in user_def.h set to TRUE for the cells of 
# interest, and recompiling the program.

#http://stackoverflow.com/questions/5359517/import-multiple-text-files-in-r-and-assign-them-names-from-a-predetermined-list
# read txt files
txt_files = list.files(path="east/",pattern = 'data*') #get Lat and Long of grids from text file names

#use list of lat and long to find which ones I need for my sites: txt_files
str(txt_files)

lat<-rep(NA,length(txt_files))
long<-rep(NA,length(txt_files))

for (i in 1:length(txt_files)){
        lat[i]<-as.numeric(substr(txt_files[i],6,12))
        long[i]<-as.numeric(substr(txt_files[i],14,21))
}

lat<-as.numeric(lat)
long<-as.numeric(long)

#combine:
W_NDX <-as.character(seq(1, length(lat), 1))
LatLonW<-cbind.data.frame(W_NDX,lat,long,rep(0,length(lat))) #lat long of Weather grids
names(LatLonW)<-c("NDX","lat","long","gage")

#add the LagLong of the gages
LatLonG<-cbind.data.frame(SiteInfoKeep$site_no,SiteInfoKeep$dec_lat_va,SiteInfoKeep$dec_long_va,rep(1,length(SiteInfoKeep$dec_long_va))) #lat long of grids
names(LatLonG)<-c("NDX","lat","long","gage")

#rbind the two sets
LatLongSites<-rbind.data.frame(LatLonG,LatLonW)
LatLon<-LatLongSites[,2:3] #grab just the lat and long coordinates

##Find Nearest site: first calculate euclidean distance between site and all other sites
dist<-earth.dist(LatLon, dist = T) #create a distance matrix (lower triangle) between each site and each other site
dist<-data.frame(dist) #save in matrix format
dist[dist==0] <- NA #replace 0s with NA

dist.data<-data.frame(Sites) #make site list into a dataframe
names(dist.data)="site_no"

for (j in 1:nrow(dist)) { #looping through sites
        LatLongSites$closest[j]<-which.min(dist[j,]) # NDX is the index number of the closest site to site j
}
#get an index on LatLongSites to figure out closest
LatLongSites$order<-seq(1, length(LatLongSites$NDX), 1)
#pull out just the weather grids
LatLongGrid<-LatLongSites[LatLongSites$gage==0,]
LatLongGrid$gage<-NULL
#drop "closest" and rename order "closest" to merge by it (also drop gage)
LatLongGrid$closest<-NULL
LatLongGrid$closest<-LatLongGrid$order
LatLongGrid$order<-NULL
names(LatLongGrid)<-c("GridNDX","G_lat","G_long","closest")

#grab just the gage sites
LatLongGage<-LatLongSites[LatLongSites$gage==1,]
LatLongGage$order<-NULL; LatLongGage$gage<-NULL #drop unecessary vars

##Now merge by closest - losing 50 sites here// because of sites with closest was another site???? Yep! 55 sites...
NN<-merge(LatLongGage,LatLongGrid, by = "closest") #Nearest Neighbor grid
NN$GridNDX<-NULL; NN$closest<-NULL #drop extra vars 
colnames(NN)[1]<-"site_no"

#pull in appropriate grid weather data
# read txt files into a list (assuming separator is a comma)
txt_files = list.files(pattern = 'data_4*') #get Lat and Long of grids from text file names
data_list = lapply(txt_files, read.table, sep = "") #missing the lat and long, need to paste that somewhere
#lat and long in same order...So need to extract the data sets and add lat and long to them

n<-nrow(data_list[[1]])

weather<-do.call(rbind, data_list)
names(weather)<-c("year","month","precipmm","maxT","minT","wind")

#make index to add to weather dataframe
NDX<-sort(rep(1:length(data_list),n))
weather$NDX<-NDX

latlong<-cbind.data.frame(lat,long,seq(1, length(lat), 1))
names(latlong)<-c("G_lat","G_long","NDX")

weath_d<-merge(weather,latlong,by="NDX")

#now merge into NN
NNweather<-merge(NN,weath_d,by=c("G_lat","G_long"))
save(NNweather,file="output/NNweather.rdata")

#delete the txt files that are super obviously too far away: all with latitudes below 40
#just move all the files into main file to start


#### 2 - GAGES II Basin Characteristics ####
#import the GAGESII climate,hydro data
#first have to copy into new file and save each as a csv file, command just reads sheet 1
GAGESIIclim<-read.csv("data/GAGESIIclimate.csv",colClasses=c("character",rep("numeric",49))); colnames(GAGESIIclim)[1]<-"site_no"
GAGESIIhydro<-read.csv("data/GAGESIIhydro.csv",colClasses=c("character",rep("numeric",33))); colnames(GAGESIIhydro)[1]<-"site_no"
GAGESIIsoil <- read.csv("data/GAGESIIsoils.csv",colClasses=c("character",rep("numeric",24))); colnames(GAGESIIsoil)[1]<-"site_no"
GAGESIItopo <- read.csv("data/GAGESIItopo.csv",colClasses=c("character",rep("numeric",12))); colnames(GAGESIItopo)[1]<-"site_no"
GAGESIIClass<-read.csv("data/BasinClassGAGESII.csv",colClasses="character"); colnames(GAGESIIClass)[1]<-"site_no"
GAGESIIClass$HYDRO_DISTURB_INDX<-as.numeric(GAGESIIClass$HYDRO_DISTURB_INDX)

#merge all the datasets into one with all basin characteristics
#only 566 sites now, lost 24 in merge
data_BC<-merge(SiteInfoKeep,GAGESIIclim,by="site_no") 
data_BC<-merge(data_BC,GAGESIIhydro,by="site_no")
data_BC<-merge(data_BC,GAGESIIsoil,by="site_no")
data_BC<-merge(data_BC,GAGESIItopo,by="site_no")
data_BC<-merge(data_BC,GAGESIIClass,by="site_no")
save(data_BC, file="data_BC.rdata")
