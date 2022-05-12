rm(list=ls())
### read the packages
#install.packages(c("readxl","xlsx","dplyr","plyr","tidyr","reshape2",
#                   "codyn","stats","vegan","bipartite","bootcluster",
#                   "betapart","gplots","ggplot2","indicspecies",
#                   "NbClust","fmsb","ncdf4","fields","raster",
#                   "ncdf.tools","lattice","colorRamps","marmap",
#                   "maptools","maps","rgdal","latticeExtra",
#                   "rerddap","Hmsc","snow","olsrr"))
library(readxl)
library(xlsx)
library(dplyr)
library(plyr)
library(tibble)
library(tidyr)
library(reshape2)
library(codyn)
library(stats)
library(vegan)
library(bipartite)
library(bootcluster)
library(betapart)
library(gplots)
library(ggplot2)
library(indicspecies)
library(NbClust)
library(fmsb)
library(ncdf4)
library(fields)
library(raster)
library(ncdf.tools)
library(lattice)
library(colorRamps)
library(marmap)
library(maptools)
library(maps)
library(rgdal)
library(latticeExtra)
library(rerddap)
library(Hmsc)
library(snow)
library(olsrr)
devtools::install_github("cran/JBTools")
devtools::install_github("cran/ncdf.tools")

#########################
### Data importation  ###
#########################
setwd('') # set the working path
sst30.file <- nc_open('Linetal_dataset_NOAA_DHW_monthly_1985-2019.nc') # import 1985-2019 SST data
taiwan <- readOGR('Linetal_dataset_TWmap.shp') # import the shapefile of Taiwan
taiwan_map <- fortify(taiwan)  # transfer .shp file to dataframe
coord <- read.csv('Linetal_dataset_coord.csv' , head = T , sep = ",")
data.fct <- read_excel("Linetal_dataset_benthicmatrix.xlsx", col_names = TRUE) # read benthic data matrix
benthic.fct <- data.fct[,c(-1,-2)] # remove the OTU & category columns
benthic.fct <- aggregate(.~fct, as.data.frame(benthic.fct), FUN=sum) # aggregate to functional level
row.names(benthic.fct) <- benthic.fct[ ,1] # assign row names to transects
benthic.fct[ ,1] <- NULL # empty functional group column
benthic.fct <- benthic.fct[-c(which(row.names(benthic.fct)=='sub_stable'), # remove the stable substrate
                              which(row.names(benthic.fct)=='sub_unstable'), # remove the unstable substrate
                              which(row.names(benthic.fct)=='unk')),] # remove the unknown actegories
benthic.fct.hell <- as.data.frame(decostand(t(benthic.fct), "hell")) # Hellinger transformation
Latitude <- read_excel("Linetal_dataset_latitude.xlsx", col_names = TRUE) # read the latitude information of transects
fish.data <- read.csv("Linetal_dataset_fishtrait.csv",  header = TRUE) # read the fish trait information
fish.data <- fish.data[-c(which(fish.data$Species=='unknown'), # remove the fish species cannot be identified to species level
                          which(fish.data$Species=='Thalassoma_sp'), # remove the fish species cannot be identified to species level
                          which(fish.data$Species=='Sufflamen_sp'), # remove the fish species cannot be identified to species level
                          which(fish.data$Species=='Serranidae_sp'), # remove the fish species cannot be identified to species level
                          which(fish.data$Species=='Scarus_sp'), # remove the fish species cannot be identified to species level
                          which(fish.data$Species=='Scaridae_sp'), # remove the fish species cannot be identified to species level
                          which(fish.data$Species=='Pseudanthias_sp'), # remove the fish species cannot be identified to species level
                          which(fish.data$Species=='Pomacentridae_sp'), # remove the fish species cannot be identified to species level
                          which(fish.data$Species=='Naso_sp'), # remove the fish species cannot be identified to species level
                          which(fish.data$Species=='Lutjanus_sp'), # remove the fish species cannot be identified to species level
                          which(fish.data$Species=='Labroides_sp'), # remove the fish species cannot be identified to species level
                          which(fish.data$Species=='Labridae_sp'), # remove the fish species cannot be identified to species level
                          which(fish.data$Species=='Halichoeres_sp'), # remove the fish species cannot be identified to species level
                          which(fish.data$Species=='Chromis_sp'), # remove the fish species cannot be identified to species level
                          which(fish.data$Species=='Chaetodon_sp'), # remove the fish species cannot be identified to species level
                          which(fish.data$Species=='Caesionidae_sp'), # remove the fish species cannot be identified to species level
                          which(fish.data$Species=='Balistidae_sp'), # remove the fish species cannot be identified to species level
                          which(fish.data$Species=='Anampses_sp'), # remove the fish species cannot be identified to species level
                          which(fish.data$Species=='Acanthurus_sp'), # remove the fish species cannot be identified to species level
                          which(fish.data$Species=='Acanthuridae_sp'), # remove the fish species cannot be identified to species level
                          which(fish.data$Species=='Sphyraena_flavicauda'),# remove the pelagic species
                          which(fish.data$Species=='Mugil_cf_cephalus'),# remove the pelagic species
                          which(fish.data$Species=='Cirrhilabrus_cf_rubrimarginatus'), # remove the benthic species
                          which(fish.data$Species=='Cirrhitichthys_cf_falco'), # remove the benthic species
                          which(fish.data$Species=='Dendrochirus_zebra'), # remove the benthic species
                          which(fish.data$Species=='Paracirrhites_cf_arcatus'), # remove the benthic species
                          which(fish.data$Species=='Paracirrhites_cf_forsteri'), # remove the benthic species
                          which(fish.data$Species=='Parajulis_cf_poecilepterus'), # remove the benthic species
                          which(fish.data$Species=='Pteragogus_enneacanthus')),] # remove the benthic species
fish.abu <- read.csv("Linetal_dataset_fishmatrix.csv", header = TRUE, row.names = 1) # read fish abundance matrix
Labridae <- read.csv("Linetal_dataset_Labridaematrix.csv", header = TRUE, row.names = 1) # read Labridae abundance matrix
Pomacentridae <- read.csv("Linetal_dataset_Pomacentridaematrix.csv", header = TRUE, row.names = 1) # read Pomacentridae abundance matrix
Chaetodontidae <- read.csv("Linetal_dataset_Chaetodontidaematrix.csv", header = TRUE, row.names = 1) # read Chaetodontidae abundance matrix

###########
### SST ###
###########
### sst in the recent 5 years (2015-2019)
sst <-'sea_surface_temperature' # extract SST
lon <- ncvar_get(sst30.file,'longitude'); nlon <- dim(lon) # extract longtitude
lat <- ncvar_get(sst30.file,'latitude'); nlat <- dim(lat) # extract latitude
tm <- ncvar_get(sst30.file,'time') # extract time
nt<-dim(tm) # extract the length of time vector
nm<-12 # number of months = 12 
ny<-nt/nm # number of years = the length of time vector / 12
date<-as.POSIXct(tm,origin='1970-01-01',tz='') # conver the date-time format
year<- as.numeric(unlist(strsplit(as.character(date),'-'))[seq(1,nt*3,by=3)]) # year vector
month<-rep(seq(1:12),ny) # month vector
sst_array<-ncvar_get(sst30.file,sst) # read sst data from the sst file
nc_close(sst30.file) # close the sst file

begyr <- 2015; endyr <- 2019; nyrs <- endyr - begyr +1 # set up the extracted period 
begobs <- ((begyr - year[1]) * nm)+1 # the year begin to observe
endobs <- ((endyr - year[1] +1) * nm) # the year end to observe
base_period <- paste(as.character(begyr)," - ", as.character(endyr), sep="") # observed period
print(c(begyr, endyr, begobs, endobs, base_period))

sst_array_base <- array(dim = c(nlon, nlat, nyrs * nm)) # creat an empty array 
sst_array_base <- sst_array[,, begobs:endobs] # fill in the sst data in the empty array

# Chaojing
mean(sst_array_base[48,18,1:12]) # 2015 monthly average
mean(sst_array_base[48,18,13:24]) # 2016 monthly average
mean(sst_array_base[48,18,25:36]) # 2017 monthly average
mean(sst_array_base[48,18,37:48]) # 2018 monthly average
mean(sst_array_base[48,18,49:60]) # 2019 monthly average

mean(sst_array_base[48,18,1:2])   # 2015 winter monthly average
mean(sst_array_base[48,18,13:14]) # 2016 winter monthly average
mean(sst_array_base[48,18,25:26]) # 2017 winter monthly average
mean(sst_array_base[48,18,37:38]) # 2018 winter monthly average
mean(sst_array_base[48,18,49:50]) # 2019 winter monthly average

mean(sst_array_base[48,18,7:8]) # 2015 summer monthly average
mean(sst_array_base[48,18,19:20]) # 2016 summer monthly average
mean(sst_array_base[48,18,31:32]) # 2017 summer monthly average
mean(sst_array_base[48,18,43:44]) # 2018 summer monthly average
mean(sst_array_base[48,18,55:56]) # 2019 summer monthly average

# Older Rocks
mean(sst_array_base[43,80,1:12]) # 2015 monthly average
mean(sst_array_base[43,80,13:24]) # 2016 monthly average
mean(sst_array_base[43,80,25:36]) # 2017 monthly average
mean(sst_array_base[43,80,37:48]) # 2018 monthly average
mean(sst_array_base[43,80,49:60]) # 2019 monthly average

mean(sst_array_base[43,80,1:2]) # 2015 winter monthly average
mean(sst_array_base[43,80,13:14]) # 2016 winter monthly average
mean(sst_array_base[43,80,25:26]) # 2017 winter monthly average
mean(sst_array_base[43,80,37:38]) # 2018 winter monthly average
mean(sst_array_base[43,80,49:50]) # 2019 winter monthly average

mean(sst_array_base[43,80,7:8]) # 2015 summer monthly average
mean(sst_array_base[43,80,19:20]) # 2016 summer monthly average
mean(sst_array_base[43,80,31:32]) # 2017 summer monthly average
mean(sst_array_base[43,80,43:44]) # 2018 summer monthly average
mean(sst_array_base[43,80,55:56]) # 2019 summer monthly average

# 5 years average sst
sst_lty_5yrs<- array(NA, dim = c(nlon, nlat, 1)) # create an empty array for loop

for (j in 1:nlon) {
  for (k in 1:nlat) {
    if (!is.na(sst_array_base[j, k, 1])) {
      sst_lty_5yrs[j, k, 1] <- mean(sst_array_base[j, k, ],na.rm=T) # calculate the mean SST 
    }}
}

# 5 years summer sst
sst_summer_5yrs <- array(NA, dim = c(nlon, nlat,2)) # create an empty array for loop

for (j in 1:nlon) {
  for (k in 1:nlat) {
    if (!is.na(sst_array_base[j, k, 1])) {
      for (m in 1:2)
        sst_summer_5yrs[j, k, m] <- mean(sst_array_base[j, k, seq(m+6, ((m+6) + nm*nyrs-1), by=nm)]) # calculate the mean SST in two summer months
    }
  }
}

sst_summer_base_5yrs <- sst_summer_5yrs
sst_summer_5yrs<- array(NA, dim = c(nlon, nlat, 1)) # create an empty array for loop

for (j in 1:nlon) {
  for (k in 1:nlat) {
    if (!is.na(sst_summer_base_5yrs[j, k, 1])) {
      sst_summer_5yrs[j, k, 1] <- mean(sst_summer_base_5yrs[j, k, ],na.rm=T) # calculate the summer mean SST 
    }}
}


# 5 years winter sst
sst_winter_5yrs <- array(NA, dim = c(nlon, nlat, 2)) # create an empty array for loop

for (j in 1:nlon) {
  for (k in 1:nlat) {
    if (!is.na(sst_array_base[j, k, 1])) {
      for (m in 1:2)
        sst_winter_5yrs[j, k, m] <- mean(sst_array_base[j, k, seq(m, (m+ nm*nyrs-1), by=nm)]) # calculate the mean SST in two winter months
    }
  }
}

sst_winter_base_5yrs<-sst_winter_5yrs
sst_winter_5yrs<- array(NA, dim = c(nlon, nlat, 1)) # create an empty array for loop

for (j in 1:nlon) {
  for (k in 1:nlat) {
    if (!is.na(sst_winter_base_5yrs[j, k, 1])) {
      sst_winter_5yrs[j, k, 1] <- mean(sst_winter_base_5yrs[j, k, ],na.rm=T) # calculate the winter mean SST 
    }}
}

### sst in the past 30 years (1985-2019)
begyr <- 1985; endyr <- 2019; nyrs <- endyr - begyr +1 # set up the extracted period 
begobs <- ((begyr - year[1]) * nm)+1 # the year begin to observe
endobs <- ((endyr - year[1] +1) * nm) # the year end to observe
base_period <- paste(as.character(begyr)," - ", as.character(endyr), sep="") # observed period
print(c(begyr, endyr, begobs, endobs, base_period))

sst_array_base <- array(dim = c(nlon, nlat, nyrs * nm)) # creat an empty array 
sst_array_base <- sst_array[,, begobs:endobs] # fill in the sst data in the empty array

# 35 years average sst
sst_lty_35yrs <- array(NA, dim = c(nlon, nlat, 1)) # create an empty array for loop

for (j in 1:nlon) {
  for (k in 1:nlat) {
    if (!is.na(sst_array_base[j, k, 1])) {
      sst_lty_35yrs[j, k, 1] <- mean(sst_array_base[j, k, ],na.rm=T) # calculate the mean sst in the past 35 years
    }}
}

# 35 years summer sst
sst_summer_35yrs <- array(NA, dim = c(nlon, nlat, 2)) # create an empty array for loop

for (j in 1:nlon) {
  for (k in 1:nlat) {
    if (!is.na(sst_array_base[j, k, 1])) {
      for (m in 1:2)
        sst_summer_35yrs[j, k, m] <- mean(sst_array_base[j, k, seq(m+6, ((m+6) + nm*nyrs-1), by=nm)]) # calculate the mean SST in summer two months of past 35 years
    }
  }
}

sst_summer_base_35yrs<-sst_summer_35yrs
sst_summer_35yrs<- array(NA, dim = c(nlon, nlat, 1)) # create an empty array for loop

for (j in 1:nlon) {
  for (k in 1:nlat) {
    if (!is.na(sst_summer_base_35yrs[j, k, 1])) {
      sst_summer_35yrs[j, k, 1] <- mean(sst_summer_base_35yrs[j, k, ],na.rm=T) # calculate the mean SST in summer of past 35 years
    }} 
}

# 35 years winter sst
sst_winter_35yrs <- array(NA, dim = c(nlon, nlat, 2)) # create an empty array for loop

for (j in 1:nlon) {
  for (k in 1:nlat) {
    if (!is.na(sst_array_base[j, k, 1])) {
      for (m in 1:2)
        sst_winter_35yrs[j, k, m] <- mean(sst_array_base[j, k, seq(m, (m + nm*nyrs-1), by=nm)]) # calculate the mean SST in winter two months of past 35 years
    }
  }
}

sst_winter_base_35yrs <- sst_winter_35yrs
sst_winter_35yrs <- array(NA, dim = c(nlon, nlat, 1)) # create an empty array for loop

for (j in 1:nlon) {
  for (k in 1:nlat) {
    if (!is.na(sst_winter_base_35yrs[j, k, 1])) {
      sst_winter_35yrs[j, k, 1] <- mean(sst_winter_base_35yrs[j, k, ],na.rm=T) # calculate the mean SST in winter of past 35 years
    }}
}

### Calculate anomaly
summer_anomaly<-sst_summer_5yrs[,,1]-sst_summer_35yrs[,,1] # SST anomaly
winter_anomaly<-sst_winter_5yrs[,,1]-sst_winter_35yrs[,,1] # winter SST anomaly
annual_anomaly<-sst_lty_5yrs[,,1]-sst_lty_35yrs[,,1] # summer SST anomaly

#############################################################
### Partition on benthic habitats & perception of fish ###
#############################################################
# function of fish specialization degree (H2') to different partitions of benthic habitats
h.fishBC <- function(cbind.mat, ncluster){ 
  kmeans <- kmeans(benthic.fct.hell, centers= ncluster, nstart=1000, iter.max=1000) # assign the number of benthic partition using k-means clustering on benthic Hellinger transformed matrix
  cluster <- kmeans$cluster # save the partition information
  cluster <- cluster[sort(names(cluster),decreasing= FALSE)] # sort the v information by alphabetical order
  d.mat <- cbind(cluster, cbind.mat) # combine the partition information with the fish matrix
  d.matrix <- aggregate(.~cluster, data = d.mat, FUN = sum) # sum the fish abundance by the partition information 
  d.matrix[,1] <- NULL # remove the column with the partition information
  d.habitat <-  networklevel(web = d.matrix , index = "H2") # calculate the fish specialization degree (H2') to the assigned benthic partition
  return(d.habitat)
}

# fish specialization degree (H2') to 2-15 partitions of benthic habitats
H2 <- matrix() 
for (i in 2:15){ 
  H2[i] <- h.fishBC(fish.abu,i)
}
H2
which.max(H2) # the most specialized partition of Labridae

# Labridae specialization degree (H2') to 2-15 partitions of benthic habitats
H2.lab.test <- matrix() 
for (i in 2:15){
  H2.lab.test[i] <- h.fishBC(Labridae, i)
}
H2.lab.test
which.max(H2.lab.test) # the most specialized partition of Labridae

# Pomacentridae specialization degree (H2') to 2-15 partitions of benthic habitats
H2.pom.test <- matrix() 
for (i in 2:15){
  H2.pom.test[i] <- h.fishBC(Pomacentridae, i)
}
H2.pom.test
which.max(H2.pom.test) # the most specialized partition of Labridae

# Chaetodontidae specialization degree (H2') to 2-15 partitions of benthic habitats
H2.cha.test <- matrix() 
for (i in 2:15){
  H2.cha.test[i] <- h.fishBC(Chaetodontidae, i)
}
H2.cha.test
which.max(H2.cha.test) # the most specialized partition of Labridae

H2.table <- rbind(round(H2,2),round(H2.lab.test,2),round(H2.pom.test,2),round(H2.cha.test,2)) # extract H2' value and make a table
rownames(H2.table) <- c('All fish species','Labridae','Pomacentridae','Chaetodontidae') # give row names to the table
colnames(H2.table) <- c(1:15) ; H2.table[,-1] # give column name to the table
H2.table <-  H2.table[,-1]
# write.csv(H2.table,'Sup5.H2.csv') # Appendix S3

### k-means cascade - use 23 indices to determine the relevent number of benthic partitions
set.seed(1)
NBkmeans.kl <- NbClust(benthic.fct.hell, diss=NULL, min.nc=2, max.nc=15, method = "kmeans", index = "kl") 
NBkmeans.ch <- NbClust(benthic.fct.hell, diss=NULL, min.nc=2, max.nc=15, method = "kmeans", index = "ch") 
NBkmeans.hartigan <- NbClust(benthic.fct.hell, diss=NULL, min.nc=2, max.nc=15, method = "kmeans", index = "hartigan") 
NBkmeans.cindex <- NbClust(benthic.fct.hell, diss=NULL, min.nc=2, max.nc=15, method = "kmeans", index = "cindex") 
NBkmeans.db <- NbClust(benthic.fct.hell, diss=NULL, min.nc=2, max.nc=15, method = "kmeans", index = "db") 
NBkmeans.silhouette <- NbClust(benthic.fct.hell, diss=NULL, min.nc=2, max.nc=15, method = "kmeans", index = "silhouette")
NBkmeans.duda <- NbClust(benthic.fct.hell, diss=NULL, min.nc=2, max.nc=15, method = "kmeans", index = "duda")
NBkmeans.pseudot2 <- NbClust(benthic.fct.hell, diss=NULL, min.nc=2, max.nc=15, method = "kmeans", index = "pseudot2") 
NBkmeans.beale <- NbClust(benthic.fct.hell, diss=NULL, min.nc=2, max.nc=15, method = "kmeans", index = "beale", alphaBeale = 0.05)
NBkmeans.ratkowsky <- NbClust(benthic.fct.hell, diss=NULL,  min.nc=2, max.nc=15, method = "kmeans", index = "ratkowsky")
NBkmeans.ball <- NbClust(benthic.fct.hell, diss=NULL, min.nc=2, max.nc=15, method = "kmeans",index = "ball")
NBkmeans.ptbiserial <- NbClust(benthic.fct.hell, diss=NULL, min.nc=2, max.nc=15, method = "kmeans", index = "ptbiserial") 
NBkmeans.gap <- NbClust(benthic.fct.hell, diss=NULL, min.nc=2, max.nc=15, method = "kmeans", index = "gap")
NBkmeans.frey <- NbClust(benthic.fct.hell, diss=NULL, min.nc=2, max.nc=15, method = "kmeans", index = "frey") 
NBkmeans.mcclain <- NbClust(benthic.fct.hell, diss=NULL,min.nc=2, max.nc=15, method = "kmeans", index = "mcclain")
NBkmeans.gamma <- NbClust(benthic.fct.hell, diss=NULL, min.nc=2, max.nc=15, method = "kmeans",index = "gamma")
NBkmeans.gplus <- NbClust(benthic.fct.hell, diss=NULL,min.nc= 2, max.nc=15, method = "kmeans",index = "gplus") 
NBkmeans.tau <- NbClust(benthic.fct.hell, diss=NULL,min.nc=2, max.nc=15, method = "kmeans", index = "tau")
NBkmeans.dunn <- NbClust(benthic.fct.hell, diss=NULL, min.nc=2, max.nc=15, method = "kmeans", index = "dunn")
NBkmeans.hubert <- NbClust(benthic.fct.hell, diss=NULL, min.nc=2, max.nc=15, method = "kmeans", index = "hubert")
NBkmeans.sdindex <- NbClust(benthic.fct.hell, diss=NULL, min.nc=2, max.nc=15, method = "kmeans", index = "sdindex") 
NBkmeans.dindex <- NbClust(benthic.fct.hell, diss=NULL, min.nc=2, max.nc=15, method = "kmeans", index = "dindex")
NBkmeans.sdbw <- NbClust(benthic.fct.hell, diss=NULL, min.nc=2, max.nc=15, method = "kmeans", index = "sdbw") 
best.nb.par <- c(NBkmeans.kl$Best.nc[1], NBkmeans.ch$Best.nc[1], NBkmeans.hartigan$Best.nc[1], # extract the best number of partitions for all indices
                 NBkmeans.cindex$Best.nc[1], NBkmeans.db$Best.nc[1], NBkmeans.silhouette$Best.nc[1],
                 NBkmeans.duda$Best.nc[1], NBkmeans.pseudot2$Best.nc[1], NBkmeans.beale$Best.nc[1],
                 NBkmeans.ratkowsky$Best.nc[1], NBkmeans.ball$Best.nc[1], NBkmeans.ptbiserial$Best.nc[1],
                 NBkmeans.gap$Best.nc[1], NBkmeans.frey$Best.nc[1], NBkmeans.mcclain$Best.nc[1],
                 NBkmeans.gamma$Best.nc[1], NBkmeans.gplus$Best.nc[1], NBkmeans.tau$Best.nc[1],
                 NBkmeans.dunn$Best.nc[1], NBkmeans.sdindex$Best.nc[1], NBkmeans.sdbw$Best.nc[1],
                 names(which.max(NBkmeans.dindex$All.index)), names(which.max(NBkmeans.hubert$All.index))) 
indexes <- c('kl','ch','hartigan','cindex','db','silhouette','duda','pseudot2','beale','ratkowsky','ball', # vector of all indices
             'ptbiserial','gap','frey','mcclain','gamma','gplus','tau','dunn','sdindex','sdbw','dindex','hubert')
index.support.nb <- data.frame(indexes, best.nb.par) # combine the indices with the best number of partitions they identified
best.nb.par.table <- colSums(table(index.support.nb)) # the best number of partitions is 2
# write.csv(index.support.nb,'Sup9.index.support.nb.csv', row.names =T) # Appendix S4

kmeans2 <- kmeans(benthic.fct.hell, centers= 2, nstart=1000, iter.max=1000) # classify 120 transects into two partitions
cluster2 <- kmeans2$cluster[sort(names(kmeans2$cluster),decreasing= FALSE)] # save the partition information
cluster2[cluster2=='1'] <- 'Subtropical' # assign the first partition to the subtropical habitat
cluster2[cluster2=='2'] <- 'Tropical' # assign the second partition to the tropical habitat
mrpp(benthic.fct.hell, kmeans2$cluster, permutations = 999) # permutation test on the significance of the partition

# test if latitudes among partitions are different 
test.lat.data <- as.data.frame(cbind(Latitude$Latitude,as.character(cluster2))) # combine the transects with their latitude information
colnames(test.lat.data) <- c("latitude","partition") # rename the column name of the latitude information matrix
mod.lat <- aov(latitude~partition, data = test.lat.data) # ANOVA test to check if latitudes among partitions are different 
summary(mod.lat)

# look for the transitional transects between partitions
set.seed(4)
benthic.fct.hell.sort <- benthic.fct.hell[sort(rownames(benthic.fct.hell)),] # sort the benthic matrix with the alphabetical order of transects
stability <- stability(benthic.fct.hell, 2, B = 20, r = 5, scheme_2 = F) # estimate k-means partition stability with bootstrapping
cluster.stability <- data.frame(rownames(benthic.fct.hell.sort), stability$obs_wise)  # save the stability information
colnames(cluster.stability) <- c('transect','stability') # rename the column name of the stability information matrix
cluster.stability[,2] <- round(cluster.stability[,2],2) # show the value of stability until two decimal places
cluster.stability[which(cluster.stability[,2]<0.9),] # show the transitional transects with stability less than 0.9
# write.csv(cluster.stability,'Sup10.cluster.stability.csv') # Appendix S5

### benthic cover in each community
# major category level
major <- aggregate(.~Major_category, as.data.frame(data.fct[,c(-1,-3)]), FUN=sum) # aggregate benthic matrix to major category level
row.names(major) <- major[ ,1] ; major <- major[ ,-1] # assign names to benthic major categories
major <- as.data.frame(t(major)) # transport the matrix
major <- major[sort(rownames(major)),] # sort the matrix with the row name
major$cluster <- cluster2 # add the partition information
major$sub <- NULL ; major$unk <- NULL # remove the substrate and unknown categories
# functional group level
fct <- as.data.frame(t(benthic.fct)[sort(rownames(t(benthic.fct))),]) # sort the benthic matrix with the rowname
fct$cluster <- cluster2 # add the partition information

# subtropical community
# major category level
sub.major <- subset(major,cluster=='Subtropical') # subset the transects belonging to subtropical community
sub.major <- sub.major[,-15] # remove the partition information
sub.major.per <- sub.major/rowSums(sub.major[,1:14]) # calculate the percentage cover of benthic major groups in each transect
mean.sub <- apply(sub.major.per,2, mean) # calculate the mean cover of benthic major groups
sd.sub <- apply(sub.major.per,2, sd) # calculate the standard deviation of benthic major group cover
# functional group level
sub.fct <- subset(fct,cluster=='Subtropical') # subset the transects belonging to subtropical community
sub.fct <- sub.fct[,-43] # remove the partition information
sub.fct.per <- sub.fct/rowSums(sub.fct[,1:42]) # calculate the percentage cover of functional groups in each transect
sub.fct.mean <- colMeans(sub.fct.per) # calculate the mean cover of functional groups
sub.fct.sd <- apply(sub.fct.per,2, sd) # calculate the standard deviation  of functional group cover

# tropical community
# major category level
tro.major <- subset(major,cluster=='Tropical') # subset the transects belonging to tropical community
tro.major <- tro.major[,-15] # remove the partition information
tro.major.per <- tro.major/rowSums(tro.major[,1:14])  # calculate the percentage cover of benthic major groups in each transect
mean.tro <- apply(tro.major.per,2, mean) # calculate the mean cover of benthic major groups
sd.tro <- apply(tro.major.per,2, sd) # calculate the standard deviation of benthic major group cover
# functional group level
tro.fct <- subset(fct,cluster=='Tropical') # subset the transects belonging to tropical community
tro.fct <- tro.fct[,-43] # remove the partition information
tro.fct.per <- tro.fct/rowSums(tro.fct[,1:42]) # calculate the percentage cover of functional groups in each transect
tro.fct.mean <- colMeans(tro.fct.per) # calculate the mean cover of functional groups
tro.fct.sd <- apply(tro.fct.per,2, sd) # calculate the standard deviation of functional groups cover

##############################################################
### Identifying fish latitudinal specialists & generalists ###
##############################################################
set.seed(7)
fish.ind <-  multipatt(fish.abu, cluster2, restcomb = c(1,2,3), control = how(nperm=999)) # calculate indval values of fish species

# potential fish generalists
fish.ge <- fish.ind$sign[which(fish.ind$sign[1]== 1 & fish.ind$sign[2]== 1),] # fish generalist matrix
fish.ge$A <- fish.ind$A[which(fish.ind$sign[1]== 1 & fish.ind$sign[2]== 1),][,3] # component A: must be 1 for every species
fish.ge$B <- fish.ind$B[which(fish.ind$sign[1]== 1 & fish.ind$sign[2]== 1),][,3] # component B: species occurrence frequency at transects
fish.ge

# fish specialists
fish.ind.s <- fish.ind$sign[fish.ind$sign[,1]==1&fish.ind$sign[,5]<=0.05,] # extract the indval of subtropical specialists
fish.ind.s <- fish.ind.s[complete.cases(fish.ind.s), ] # remove the empty rows
fish.ind.s.A <- fish.ind$A[fish.ind$sign[,1]==1&fish.ind$sign[,5]<=0.05,] # extract the component A of subtropical specialists
fish.ind.s$A <- fish.ind.s.A[complete.cases(fish.ind.s.A), ][,1] # remove the empty rows
fish.ind.s.B <- fish.ind$B[fish.ind$sign[,1]==1&fish.ind$sign[,5]<=0.05,] # extract the component B of subtropical specialists
fish.ind.s$B <- fish.ind.s.B[complete.cases(fish.ind.s.B), ][,1] # remove the empty rows
fish.ind.t <- fish.ind$sign[fish.ind$sign[,2]==1&fish.ind$sign[,5]<=0.05,] # extract the indval of tropical specialists
fish.ind.t <- fish.ind.t[complete.cases(fish.ind.t), ] # remove the empty rows
fish.ind.t.A <- fish.ind$A[fish.ind$sign[,2]==1&fish.ind$sign[,5]<=0.05,] # extract the component A of tropical specialists
fish.ind.t$A <- fish.ind.t.A[complete.cases(fish.ind.t.A), ][,2] # remove the empty rows
fish.ind.t.B <- fish.ind$B[fish.ind$sign[,2]==1&fish.ind$sign[,5]<=0.05,] # extract the component B of tropical specialists
fish.ind.t$B <- fish.ind.t.B[complete.cases(fish.ind.t.B), ][,2] # remove the empty rows
fish.ind.data <- round(rbind(fish.ind.s,fish.ind.t,fish.ge)[,3:7],2) # combine the information of specialists and generalists
fish.ind.data$Habitat[fish.ind.data$index==1] <- 'Subtropical' # label the subtropical specialists
fish.ind.data$Habitat[fish.ind.data$index==2] <- 'Tropical' # label the tropical specialists
fish.ind.data$Habitat[fish.ind.data$index==3] <- 'Both' # label the potential generalists
fish.ind.data$index <- NULL # empty the index column
colnames(fish.ind.data) <- c('IndVal','p value','A','B','Habitat') # change the column name
# write.csv(fish.ind.data,'Sup11.fish.indicator.csv') # Appendix S6

# evaluate if fish???s IndVal in one single partition (specialization) varied with their latitudinal or temperature preferences
fish.ind.tro <- fish.ind$str[,2] # extract the indvals in tropical region
fish.ind.sub <- fish.ind$str[,1] # extract the indvals in subropical region

fish.ind.la <- cbind(Latitude[,2],fish.abu) # combine the latitude and fish abundance information
fish.ind.la <- aggregate(.~Latitude, fish.ind.la, sum) # aggregate fish abundance by latitude (location)
fish.ind.la[,2:185] [fish.ind.la[,2:185] >0] <-1 # change abundance to presence-absence matrix
latitude <- fish.ind.la[,1] # save the latitude information
fish.ind.sub.la <- t(fish.ind.sub * t(fish.ind.la[,2:185])) # calculate the subtropical indvals in each location 
fish.ind.sub.la <- rowMeans(fish.ind.sub.la) # calculate the mean indvals in each location 
fish.ind.sub.la <- as.data.frame(cbind(latitude ,fish.ind.sub.la)) # add the latitude information
colnames(fish.ind.sub.la) <- c('Latitude','subtropical') # add the columnname
lm.sub.la <- lm(subtropical~Latitude, fish.ind.sub.la)  # linear regression between indvals in subtropical habitat and latitude
summary(lm.sub.la) # check the result

fish.ind.t <- cbind(Latitude[,3],fish.abu) # combine the temperature and fish abundance information
fish.ind.t <- aggregate(.~avg_sst, fish.ind.t, sum) # aggregate fish abundance by 2019 monthly average SST
fish.ind.t[,2:185] [fish.ind.t[,2:185] >0] <-1 # change abundance to presence-absence matrix
avg_sst <- fish.ind.t[,1] # save the 2019 monthly average SST information
fish.ind.sub.t <- t(fish.ind.sub * t(fish.ind.t[,2:185])) # calculate the subtropical indvals in the location with the same SST
fish.ind.sub.t <- rowMeans(fish.ind.sub.t) # calculate the mean indvals in the location with the same SST
fish.ind.sub.t <- as.data.frame(cbind(avg_sst,fish.ind.sub.t)) # add the SST information
colnames(fish.ind.sub.t) <- c('avg_sst','subtropical') # add the column name
lm.sub <- lm(subtropical~avg_sst, fish.ind.sub.t) # linear regression between indvals in subtropical habitat and the 2019 monthly average SST
summary(lm.sub) # check the result

fish.ind.tro.la <- t(fish.ind.tro * t(fish.ind.la[,2:185])) # calculate the tropical indvals in each location 
fish.ind.tro.la  <- rowMeans(fish.ind.tro.la) # calculate the mean indvals in each location 
fish.ind.tro.la <- as.data.frame(cbind(latitude ,fish.ind.tro.la)) # add the latitude information
colnames(fish.ind.tro.la) <- c('Latitude','tropical') # add the column name
lm.tro.la <- lm(tropical~Latitude, fish.ind.tro.la) # linear regression between indvals in tropical habitat and latitude
summary(lm.tro.la) # check the result

fish.ind.tro.t <- t(fish.ind.tro * t(fish.ind.t[,2:185]) )# calculate the tropical indvals in the location with the same SST
fish.ind.tro.t  <- rowMeans(fish.ind.tro.t) # calculate the mean indvals in the location with the same SST
fish.ind.tro.t <- as.data.frame(cbind(avg_sst,fish.ind.tro.t)) # add the SST information
colnames(fish.ind.tro.t) <- c('avg_sst','tropical') # add the column name
lm.tro <- lm(tropical~avg_sst, fish.ind.tro.t) # linear regression between indvals in tropical habitat and the 2019 monthly average SST
summary(lm.tro) # check the result

################################################################
### Responses of fish specialists to SST and benthic habitat ###
################################################################
# environmental covariate - 5 yrs ave SST + benthic habitat
env <- data.frame(Latitude,cluster2)
env$avg_sst <- as.numeric(env$avg_sst)
env$cluster2 <- as.factor(env$cluster2)
X <- as.data.frame(env[,c('avg_sst','cluster2')])
XFormula <- ~ avg_sst + cluster2

# study design - transect + depth + site 
design_df <- matrix(unlist(strsplit(env$Transect,'_')), ncol = 4, byrow = T)
design_df_site <- design_df[,2]
design_df_depth <- design_df[,3]

design <- data.frame(as.factor(design_df_site),as.factor(design_df_depth),as.factor(env$Transect))
colnames(design) <- c('site','depth', 'transect')
site <- HmscRandomLevel(units = design$site)
depth <- HmscRandomLevel(units = design$depth)
transect <- HmscRandomLevel(units = design$transect)
ranlevels = list(site = site, depth = depth, transect = transect) 

# fish species matrix - specialists only
fish.ind.data$sp <- rownames(fish.ind.data)
fish.abu.hmsc <- as.data.frame(matrix(unlist(t(fish.abu)), ncol = 120, byrow = F))
colnames(fish.abu.hmsc) <- rownames(fish.abu)
rownames(fish.abu.hmsc) <- colnames(fish.abu)
fish.abu.hmsc$sp <- rownames(fish.abu.hmsc)
Y <- full_join(x = fish.ind.data, y = fish.abu.hmsc, by = 'sp')
Y_sp <- Y[Y$Habitat=='Tropical'|Y$Habitat=='Subtropical',]
Y_sp <- Y_sp[!is.na(Y_sp$IndVal),]
order_2 <- c('Subtropical','Tropical')
Y_sp <- Y_sp %>% 
  arrange(factor(Habitat, levels = order_2))
sp_sp <- Y_sp$sp
Y_sp_hmsc <- as.data.frame(Y_sp[,7:ncol(Y_sp)] )
rownames(Y_sp_hmsc) <- sp_sp
Y_sp_hmsc <- as.data.frame(t(Y_sp_hmsc))

# create model
simul_sp <- Hmsc(Y=Y_sp_hmsc, XData = X,
                 XFormula = XFormula,
                 studyDesign = design,
                 ranLevels  = ranlevels,
                 distr = "lognormal") 

# run model
thin = 100
samples = 10000
nChains = 4
set.seed(1)
mod_HMSC_sp = sampleMcmc(simul_sp,
                         samples = samples,
                         thin = thin,
                         transient = ceiling(0.5*samples*thin),
                         nChains = nChains, 
                         nParallel = nChains)
# it takes a while to run this model (~days), you can save it in .rds after
# finish running in case you gonna run it again in the future
# saveRDS(mod_HMSC_sp, "mod_HMSC_sp.rds")

# Convergence tests
mcoda_sp <- convertToCodaObject(mod_HMSC_sp)

# Visual check for different coefficients of interest 
plot(mcoda_sp$Beta[,1:6])
plot(mcoda_sp$Beta[,7:12])

# Gelman's diagnosis, which should be at most close to 1.0 for good convergence.
es.beta.sp <- effectiveSize(mcoda_sp$Beta)
ge.beta.sp <- gelman.diag(mcoda_sp$Beta,multivariate=F)$psrf
es.V.sp <- effectiveSize(mcoda_sp$V)

hist(ge.beta.sp[,1] , main ='MCMC convergence diagnostic',
     xlab = 'Potential scale reduction factors (parameter ??)') # Appendix S2
hist(es.beta.sp, main ='MCMC convergence diagnostic',
     xlab = 'Effective sample sizes (parameter ??)') # Appendix S2

## variance explained by the model
partition_sp <- createPartition(hM = mod_HMSC_sp, nfolds=10, column= "transect")
pred_Y_sp <- computePredictedValues(mod_HMSC_sp, expected = F)
MF_sp <- evaluateModelFit(hM = mod_HMSC_sp, predY = pred_Y_sp)
round(mean(MF_sp$AUC),2)

# variance partition
VP_sp <- computeVariancePartitioning(mod_HMSC_sp)

###################################################
### Fish characteristics among benthic habitats ###
###################################################
# fish evenness 
fish.even <- fish.abu # creat fish evenness matrix 
fish.even$transect <- rownames(fish.abu) # add the transect information to the matrix
fish.even <- gather(fish.even, key = "Fish", value = "Number",-transect) # transform the matrix with transect, fish, and number of fish as columns
fish.evenness <- community_structure(fish.even, replicate.var = 'transect', abundance.var = 'Number', metric ='Evar') # calculate the fish evenness and richness at transect level
rownames(fish.evenness) <- fish.evenness[,1] ;fish.evenness[,1] <- NULL # add the transect column to rownames and delete the column
fish.evenness <- as.data.frame(cbind(cluster2,fish.evenness)) # add the habitats that each transect belonging to  
colnames(fish.evenness) <- c('Habitat','Richness','Evar') # change the columnname
shapiro.test(fish.evenness$Evar) # test if evenness is normally distributed
kruskal.test(Evar ~ factor(Habitat), fish.evenness) # test if evenness is different among habitats

# fish richness 
shapiro.test(fish.evenness$Richness) # test if richness is normally distributed
kruskal.test(Richness ~ factor(Habitat), fish.evenness) # test if richness is different among habitats

# fish abundance
fish.abudance <- rowSums(fish.abu) # calculate the fish abundance in each transect
fish.abudance <- as.data.frame(cbind(cluster2 ,as.numeric(fish.abudance))) # add the habitats that each transect belonging to  
colnames(fish.abudance) <- c('Habitat','Abudance') # change the columnname
shapiro.test(as.numeric(fish.abudance$Abudance)) # test if abundance is normally distributed
kruskal.test(Abudance ~ Habitat, fish.abudance) # test if abundance is different among habitats

# fish biomass
fish.data$Tr<- paste(fish.data$Region, fish.data$Area, fish.data$Depth, fish.data$Transect, sep='_') # add the transect information to the matrix
fish.data$biomass <- round((fish.data$a * fish.data$Length ^ fish.data$b) * fish.data$Number) # calculate biomass for each fish individual
fish.bio <- dcast(fish.data, Tr~Species, value.var = 'biomass', sum) # fish biomass by transect
rownames(fish.bio) <- fish.bio[,1] ; fish.bio[,1] <- NULL # add transect names and remove the transect column
fish.bio <- rbind(fish.bio, 'E_JH_15_T1' = seq(0, by=0 ,length= nrow(fish.bio))) # add the transects without any fish
fish.bio <- rbind(fish.bio, 'E_JH_15_T3' = seq(0, by=0 ,length= nrow(fish.bio)))
fish.bio <- rbind(fish.bio, 'E_SY_5_T5' = seq(0, by=0 ,length= nrow(fish.bio)))
fish.bio <- fish.bio[order(row.names(fish.bio)), ] # reorder the matrix by fish name

fish.biomass <- rowSums(fish.bio, na.rm = T) # calculate the fish biomass in each transect
fish.biomass <- as.data.frame(cbind(cluster2 ,as.numeric(fish.biomass))) # add the habitats that each transect belonging to  
colnames(fish.biomass) <- c('Habitat','Biomass') # add transect names
shapiro.test(as.numeric(fish.biomass$Biomass)) # test if abundance is normally distributed
kruskal.test(Biomass ~ factor(Habitat), fish.biomass) # test if abundance is different among habitats

# beta diversity among habitats
cluster2.abu <- data.frame(cluster2,fish.abu) # add the habitats that each transect belonging to
cluster2.abu <- aggregate(.~cluster2, cluster2.abu, FUN=sum)  # sum the fish abundance by habitats
cluster2.abu <- cluster2.abu[,-1] # remove the habitat column
cluster2.pa <- replace(cluster2.abu,cluster2.abu>0,1) # change the abundance matrix to presence-absence matrix
cluster2.core <- betapart.core(cluster2.pa) # get betapart objects
cluster2.multi <-  beta.multi(cluster2.core, index.family="jaccard") # multiple habitat measures

# fish trophic structure among habitats
diet.richness <- data.frame(fish.data$Species,fish.data$diet) # extract the fish species and diet information
diet.richness <- t(sort(unique(diet.richness))) # assign the major diet to the species
colnames(diet.richness) <- diet.richness[1,] # add fish species as the column name
diet.richness <- diet.richness[-1,] # remove the row with fish species
diet.pa <- data.frame(diet.richness, t(cluster2.pa)) # combine the fish diet information with their occurrence
colnames(diet.pa) <- c('Diet', 'Subtropical','Tropical') # change the column name
diet.pa <- aggregate(.~Diet, diet.pa, FUN=sum) # calculate the number of species in 6 diet groups among two habitats
rownames(diet.pa) <- diet.pa[,1] # add fish species as the row name
diet.pa <- t(diet.pa[,-1]) # remove the column with fish species

diet.abu <- data.frame(diet.richness, t(cluster2.abu)) # combine the fish diet information with their abundance
colnames(diet.abu) <- c('Diet', 'Subtropical','Tropical') # change the column name
diet.abu <- aggregate(.~Diet, diet.abu, FUN=sum) # calculate the fish abundance in 6 diet groups among two habitats
rownames(diet.abu) <- diet.abu[,1] # add diet group as the row name
diet.abu <- t(diet.abu[,-1]) # remove the column with diet group
diet.abu.per <- diet.abu/rowSums(diet.abu) *100 # percentage of the diet group in each habitat

############################
### Figures & appendixes ###
############################
## Figure 1
createScaleBar <- function(lon,lat,distanceLon,distanceLat,distanceLegend, dist.units = "km"){
  # First rectangle
  bottomRight <- gcDestination(lon = lon, lat = lat, bearing = 90, dist = distanceLon, dist.units = dist.units, model = "WGS84")
  
  topLeft <- gcDestination(lon = lon, lat = lat, bearing = 0, dist = distanceLat, dist.units = dist.units, model = "WGS84")
  rectangle <- cbind(lon=c(lon, lon, bottomRight[1,"long"], bottomRight[1,"long"], lon),
                     lat = c(lat, topLeft[1,"lat"], topLeft[1,"lat"],lat, lat))
  rectangle <- data.frame(rectangle, stringsAsFactors = FALSE)
  
  # Second rectangle t right of the first rectangle
  bottomRight2 <- gcDestination(lon = lon, lat = lat, bearing = 90, dist = distanceLon*2, dist.units = dist.units, model = "WGS84")
  rectangle2 <- cbind(lon = c(bottomRight[1,"long"], bottomRight[1,"long"], bottomRight2[1,"long"], bottomRight2[1,"long"], bottomRight[1,"long"]),
                      lat=c(lat, topLeft[1,"lat"], topLeft[1,"lat"], lat, lat))
  rectangle2 <- data.frame(rectangle2, stringsAsFactors = FALSE)
  
  # Now let's deal with the text
  onTop <- gcDestination(lon = lon, lat = lat, bearing = 0, dist = distanceLegend, dist.units = dist.units, model = "WGS84")
  onTop2 <- onTop3 <- onTop
  onTop2[1,"long"] <- bottomRight[1,"long"]
  onTop3[1,"long"] <- bottomRight2[1,"long"]
  legend <- rbind(onTop, onTop2, onTop3)
  legend <- data.frame(cbind(legend, text = c(0, distanceLon, distanceLon*2)), stringsAsFactors = FALSE, row.names = NULL)
  return(list(rectangle = rectangle, rectangle2 = rectangle2, legend = legend))
}
scaleBar <- function(lon, lat, distanceLon, distanceLat, distanceLegend, dist.unit = "km", rec.fill = "white", rec.colour = "black", rec2.fill = "black", rec2.colour = "black", legend.colour = "black", legend.size = 3, orientation = TRUE, arrow.length = 500, arrow.distance = 300, arrow.North.size = 6){
  laScaleBar <- createScaleBar(lon = lon, lat = lat, distanceLon = distanceLon, distanceLat = distanceLat, distanceLegend = distanceLegend, dist.unit = dist.unit)
  # First rectangle
  rectangle1 <- geom_polygon(data = laScaleBar$rectangle, aes(x = lon, y = lat), fill = rec.fill, colour = rec.colour)
  
  # Second rectangle
  rectangle2 <- geom_polygon(data = laScaleBar$rectangle2, aes(x = lon, y = lat), fill = rec2.fill, colour = rec2.colour)
  
  # Legend
  scaleBarLegend <- annotate("text", label = paste(laScaleBar$legend[,"text"], dist.unit, sep=""), x = laScaleBar$legend[,"long"], y = laScaleBar$legend[,"lat"], size = legend.size, colour = legend.colour)
  
  res <- list(rectangle1, rectangle2, scaleBarLegend)
  
  if(orientation){# Add an arrow pointing North
    coordsArrow <- createOrientationArrow(scaleBar = laScaleBar, length = arrow.length, distance = arrow.distance, dist.unit = dist.unit)
    arrow <- list(geom_segment(data = coordsArrow$res, aes(x = x, y = y, xend = xend, yend = yend)), annotate("text", label = "N", x = coordsArrow$coordsN[1,"x"], y = coordsArrow$coordsN[1,"y"], size = arrow.North.size, colour = "black"))
    res <- c(res, arrow)
  }
  return(res)
}

ggplot() +
  geom_polygon(data= taiwan_map, aes(x= long, y= lat, group= group), fill='light blue') +
  geom_point(data= coord, aes(x = Longitude, y = Latitude), color="red", size = 1.5)+
  scale_x_continuous(limits = c(119, 123)) +                                                          # set x and y limt
  scale_y_continuous(limits = c(21.4, 25.5)) +
  annotate("text", x = 120.8, y = 23.6, label = "Taiwan", size=8) +
  theme_minimal() +                                                                                      # background                                                                      # delete legend
  theme(legend.position="none")+
  coord_map(projection ="mollweide")+ 
  scaleBar(lon = 122, lat = 21.4, distanceLon = 40, distanceLat = 4, 
           distanceLegend = 20, dist.unit = "km", orientation = F, legend.size = 3)

## Figure 2
# a. Yearly average of monthly SST in 2015-2019
sst_slice <- sst_lty_5yrs[,, 1]
grid <- expand.grid(lon=lon, lat=lat)
cutpts<-seq(13,31,by=0.5)
levelplot(sst_slice ~ lon * lat, data=grid, at=cutpts, cuts=36, pretty=T, 
          col.regions=(matlab.like2(480)),xlim = c(119.8,122.5),ylim = c(21.5,26),
          xlab = 'Longitude',ylab='Latitude')+
  layer(sp.lines(taiwan,fill='gray80'))

# b. Yearly average of summer monthly SST in 2015-2019
sst_slice <- sst_summer_5yrs[,, 1]
grid <- expand.grid(lon=lon, lat=lat)
cutpts<-seq(13,31,by=0.5)
levelplot(sst_slice ~ lon * lat, data=grid, at=cutpts, cuts=36, pretty=T, 
          col.regions=(matlab.like2(480)),xlim = c(119.8,122.5),ylim = c(21.5,26),
          xlab = 'Longitude',ylab='Latitude')+
  layer(sp.lines(taiwan,fill='gray80'))

# c. Yearly average of winter monthly SST in 2015-2019
sst_slice <- sst_winter_5yrs[,, 1]
grid <- expand.grid(lon=lon, lat=lat)
cutpts<-seq(13,31,by=0.5)
levelplot(sst_slice ~ lon * lat, data=grid, at=cutpts, cuts=36, pretty=T, 
          col.regions=(matlab.like2(540)),xlim = c(119.8,122.5),ylim = c(21.5,26),
          xlab = 'Longitude',ylab='Latitude')+
  layer(sp.lines(taiwan,fill='gray80'))

# d. SST anomaly comparing the differences between the recent 5 years and the past 30 years
bwr <- colorRampPalette(colors = c("blue", "white", "red"))
grid <-expand.grid(lon=lon, lat=lat)
cutpts <- seq(-2,2,by=0.1)
levelplot(summer_anomaly ~ lon * lat, data=grid, at=cutpts, cuts=40, pretty=T, 
          col.regions=(bwr(400)),xlim = c(119.8,122.5),ylim = c(21.5,26),
          xlab = 'Longitude',ylab='Latitude')+
  layer(sp.lines(taiwan,fill='gray80'))

# e. Winter SST anomaly comparing the differences between the recent 5 years and the past 30 years
grid <-expand.grid(lon=lon, lat=lat)
cutpts <- seq(-2,2,by=0.1)
levelplot(winter_anomaly ~ lon * lat, data=grid, at=cutpts, cuts=40, pretty=T, 
          col.regions=(bwr(400)),xlim = c(119.8,122.5),ylim = c(21.5,26),
          xlab = 'Longitude',ylab='Latitude')+
  layer(sp.lines(taiwan,fill='gray80'))

# f. Summer SST anomaly comparing the differences between the recent 5 years and the past 30 years
grid <-expand.grid(lon=lon, lat=lat)
cutpts <- seq(-2,2,by=0.1)
levelplot(annual_anomaly ~ lon * lat, data=grid, at=cutpts, cuts=40, pretty=T, margin = F,
          col.regions=(bwr(400)),xlim = c(119.8,122.5),ylim = c(21.5,26),
          xlab = 'Longitude',ylab='Latitude')+
  layer(sp.lines(taiwan,fill='gray80'))

# Fig 3 - bipartite plot visiaulizing fish perception to 2 benthic partitions
col <- c("firebrick1","DeepSkyBlue")
cluster2.abu <- data.frame(cluster2,fish.abu)
cluster2.web <- aggregate(.~cluster2, as.data.frame(cluster2.abu), FUN=sum)
cluster2.web <- t(cluster2.web[,-1])
colnames(cluster2.web) <- c('Subtropical','Tropical')
plotweb(sortweb(t(cluster2.web), sort.order="dec"),
        text.rot=90,
        high.spacing = 0.01,
        low.spacing = .1,
        ybig=1.2,
        high.lab.dis = 0,
        method = 'normal',
        col.low =rev(col)) 

# Labridae
cluster2.lab <- data.frame(cluster2,Labridae)
cluster2.lab <- aggregate(.~cluster2, as.data.frame(cluster2.lab), FUN=sum)
cluster2.lab <- cluster2.lab[,-1]
rownames(cluster2.lab) <- c('Subtropical','Tropical')
plotweb(sortweb(cluster2.lab, sort.order="dec"),
        text.rot=90,
        high.spacing = 0.01,
        low.spacing = .1,
        ybig=1.2,
        high.lab.dis = 0,
        method = 'normal',
        col.low = col) 

# Pomacentridae
cluster2.pom <- data.frame(cluster2,Pomacentridae)
cluster2.pom <- aggregate(.~cluster2, as.data.frame(cluster2.pom), FUN=sum)
cluster2.pom <- cluster2.pom[,-1]
rownames(cluster2.pom) <- c('Subtropical','Tropical')
plotweb(sortweb(cluster2.pom, sort.order="dec"),
        text.rot=90,
        high.spacing = 0.01,
        low.spacing = .1,
        ybig=1.2,
        high.lab.dis = 0,
        method = 'normal',
        col.low = rev(col)) 
dev.off()

# Chaetodontidae
cluster2.cha <- data.frame(cluster2,Chaetodontidae)
cluster2.cha <- aggregate(.~cluster2, as.data.frame(cluster2.cha), FUN=sum)
cluster2.cha <- cluster2.cha[,-1]
rownames(cluster2.cha) <- c('Subtropical','Tropical')
plotweb(sortweb(cluster2.cha, sort.order="dec"),
        text.rot=90,
        high.spacing = 0.01,
        low.spacing = .1,
        ybig=1.2,
        high.lab.dis = 0,
        method = 'normal',
        col.low = col) 

# Fig 4 - fish specialists
# a - IndVals of 58 fish specialists 
fish.ind.inf <- fish.ind.data
fish.ind.inf <- fish.ind.inf[fish.ind.inf$Habitat!='Both',]
fish.ind.inf$`p value` <- NULL
fish.ind.inf$sp <- rownames(fish.ind.inf)
fish.ind.inf <- melt(fish.ind.inf,id.vars = c( "sp",'Habitat'), variable.name='IndVal', value.name = "value")

specialists <- ggplot(fish.ind.inf ,aes(x = IndVal, y =  sp, fill = Habitat))+
  geom_point(aes(size = value), alpha = 0.5, shape = 21) +
  labs( x= "", y = "")  +
  scale_y_discrete(limits = rev(fish.ind.inf$sp[1:58]))+
  theme_bw()
specialists + theme(axis.text.y = element_text(face = "italic")) 
dev.off()

# b - responses of 58 fish specialists to SST and benthic habitat. 
par(mar=c(7,15,2.5,0))
plotBeta(mod_HMSC_sp,
         post = postBeta_sp, 
         plotTree = F,
         covNamesNumbers = c(T,F),
         spNamesNumbers = c(T,F))
dev.off()
# in this plot, the colors (red & blue) show at least 95% 
# posterior probability of being positive or negative

# Appendix S7
fish.sp.tro <- fish.abu[,c(4,5,10,11,17,24,30,34,35,40,41,46,48,58,64,67,69,76,77,89,99,
                           104,109,116,122,123,124,125,128,130,134,139,140,141,147,148,
                           149,151,152,153,155,163,167,169,172,175,177,180,181,182,184)]
fish.sp.sub <- fish.abu[,c(12,65,80,146,158,168,170)]

fish.la.dist.sub <- cbind(Latitude[,2],fish.sp.sub)
fish.la.dist.sub <- aggregate(.~Latitude, as.data.frame(fish.la.dist.sub), FUN=sum)
fish.la.dist.sub <- gather(fish.la.dist.sub,'Species','abundance',-Latitude)
fish.la.dist.sub <- uncount(fish.la.dist.sub,abundance)
Habitat <- rep('Subtropical',nrow(fish.la.dist.sub))
fish.la.dist.sub <- cbind(fish.la.dist.sub,Habitat)

fish.la.dist.tro <- cbind(Latitude[,2],fish.sp.tro)
fish.la.dist.tro <- aggregate(.~Latitude, as.data.frame(fish.la.dist.tro), FUN=sum)
fish.la.dist.tro <- gather(fish.la.dist.tro,'Species','abundance',-Latitude)
fish.la.dist.tro <- uncount(fish.la.dist.tro,abundance)
Habitat <- rep('Tropical',nrow(fish.la.dist.tro))
fish.la.dist.tro <- cbind(fish.la.dist.tro,Habitat)
fish.la.dist <- rbind(fish.la.dist.sub, fish.la.dist.tro)

ggplot(fish.la.dist,aes(x = reorder(Species,Latitude), y = Latitude, fill = Habitat))+
  geom_boxplot(outlier.size=1)+
  scale_fill_manual(values = rev(col)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  scale_y_continuous(name="Latitude") +
  theme(axis.text.x = element_text(face = "italic"))+
  labs(x= 'Species')

# Appendix S8
ggplot(fish.ind.sub.la, aes(x = Latitude, y = subtropical)) +
  geom_point(color= "DeepSkyBlue") +
  geom_smooth(method = "lm")+
  scale_y_continuous(name="Mean IndVals in the subtropical habitat") +
  scale_x_continuous(name ="Latitude")+
  theme(axis.text=element_text(size=12),text = element_text(size = 15))

ggplot(fish.ind.sub.t, aes(x = avg_sst, y = subtropical)) + 
  geom_point(color= "DeepSkyBlue") +
  geom_smooth(method = "lm")  +
  scale_y_continuous(name="Mean IndVals in the subtropical habitat") +
  scale_x_continuous(name ="2019 monthly average SST")+
  theme(axis.text=element_text(size=12),text = element_text(size = 15))

ggplot(fish.ind.tro.la , aes(x = Latitude, y = tropical)) + 
  geom_point(color="firebrick1") +
  geom_smooth(method = "lm")+
  scale_y_continuous(name="Mean IndVals in the tropical habitat") +
  scale_x_continuous(name ="Latitude")+
  theme(axis.text=element_text(size=12),text = element_text(size = 15))

ggplot(fish.ind.tro.t, aes(x = avg_sst, y = tropical)) +
  geom_point(color="firebrick1") +
  geom_smooth(method = "lm")+
  scale_y_continuous(name="Mean IndVals in the tropical habitat") +
  scale_x_continuous(name ="2019 monthly average SST")+
  theme(axis.text=element_text(size=12),text = element_text(size = 15))

# Appendix S9
par(mar=c(15,5,5,0))
barplot(VP_sp$vals,las = 2, horiz=F, col = c('red','coral','orange','yellow','gold'))
plotVariancePartitioning(mod_HMSC_sp, VP = VP_sp, las = 2, horiz=F)
dev.off()

# Appendix S10
ggplot(fish.evenness, aes( x= Habitat, y = Richness,fill=Habitat))+ 
  scale_fill_manual(values = rev(col))+
  geom_boxplot(outlier.size=1) + 
  theme(text = element_text(size = 18))

ggplot(fish.evenness, aes( x= Habitat, y = Evar,fill=Habitat))+ 
  scale_fill_manual(values = rev(col))+
  geom_boxplot(outlier.size=1) + 
  theme(text = element_text(size = 18))

ggplot(fish.abudance, aes( x= Habitat, y = as.numeric(Abudance), fill=Habitat))+ 
  scale_fill_manual(values = rev(col))+
  geom_boxplot(outlier.size=1)+ 
  scale_y_log10(name ="Abundance (ind.)")+ 
  theme(text = element_text(size = 18))

ggplot(fish.biomass, aes( x= Habitat, y = as.numeric(Biomass),fill=Habitat))+ 
  scale_fill_manual(values = rev(col))+
  geom_boxplot(outlier.size=1) + 
  scale_y_log10(name = "Biomass (g/100m^2)")+ 
  theme(text = element_text(size = 18))

# Appendix S11 - trophic structure of fish among habitats
max_min <- data.frame(HD = c(5, 0), MH = c(5, 0), MI = c(5, 0),
                      O = c(5, 0), PI = c(5, 0), PL = c(5, 0), SI = c(5, 0))
rownames(max_min) <- c("Max", "Min")
diet.abu.s.plot <- rbind(max_min, log10(diet.abu[1,]))
diet.abu.t.plot <- rbind(max_min, log10(diet.abu[2,]))
par(mfrow = c(1,2))
radarchart(diet.abu.s.plot, axistype=1, pcol= rgb(0,0.75,1,0.5) , pfcol=  rgb(0,0.75,1,0.5), plwd=2 , seg=5,
           vlcex= 1.5, cglcol="grey", cglty=1, axislabcol="grey", caxislabels=c(0,10,100,1000,10000,100000))

radarchart(diet.abu.t.plot, axistype=1, pcol= rgb(1,0.19,0.19,0.5) , pfcol=  rgb(1,0.19,0.19,0.5), plwd=2 , seg=5,
           vlcex= 1.5, cglcol="grey", cglty=1, axislabcol="grey", caxislabels=c(0,10,100,1000,10000,100000))
dev.off()
