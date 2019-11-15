

# Load libraries ----------------------------------------------------------
library(pacman)
pacman::p_load(tidyverse, raster, rgdal, cclust, dismo, gtools, sp, rgeos, FactoMineR, pROC, randomForest, Hmisc, velox, rgeos)
source("FunctionsRFclustering.R")
source("rfClust.R")

# Load data ---------------------------------------------------------------
nic <- raster::getData("GADM", country = "NIC", level = 1)
pnt <- shapefile("../data/occ/presences.shp")
plot(nic)
plot(pnt, add = TRUE, cex = 2, col = "red", pch = 16)
myproj <- CRS('+proj=longlat +datum=WGS84')

# Climate data
fls <- list.files("../data/climate/current", full.names = TRUE, pattern = ".asc$")
prec <- grep('prec', fls, value = TRUE) %>% stack()
tmax <- grep('tmax', fls, value = TRUE) %>% stack()
tmea <- grep('tmean', fls, value = TRUE) %>% stack()
tmin <- grep('tmin', fls, value = TRUE) %>% stack()
plot(prec[[1]])

crn <- biovars(prec, tmin, tmax)
plot(crn[[1]])

# Extraccion values to points
occ <- raster::extract(crn, pnt) %>% 
  cbind(coordinates(pnt), .) %>% 
  as.data.frame %>% 
  setNames(c('x', 'y', paste0('bio', 1:19)))

# Cluster analysis
env_values <- as.matrix(occ[,3:ncol(occ)])
datRF <- as.data.frame(occ[,3:ncol(occ)])
rfClust <- rf.clust(occ = occ, nforest = 25, ntrees = 100, nVars = 8, nclasses = 3)
labelRF <- rfClust[[1]]
clusterdata <- rfClust[[2]]
classdata <- cbind(pb = as.factor(labelRF), occ[,3:ncol(occ)])
clusteredpresdata <- cbind(occ, cluster = labelRF) %>% na.omit() %>% tbl_df()
no.clusters <- 3
dir.create('../rData/run1', recursive = TRUE)
save(datRF, file = paste0('../rData/run1/datRF.rData'))
save(clusterdata, file = paste0('../rData/run1/clusterdata.rData'))
save(occ, clusteredpresdata, no.clusters, labelRF, file = paste0('../rData/run1/clustereddata.rData'))

# Backgroun generate
SPspecies <- SpatialPoints(occ[,1:2]) 
crs(SPspecies) <- myproj
back_raster <- crn[[1]] * 0
speciescell <- raster::extract(crn[[1]] * 0, SPspecies, cellnumber = TRUE)
back_raster[speciescell[,1]]  <- NA
samplesize <- round(min(summary(as.factor(clusteredpresdata$cluster))) / 2, 0) 
NumberOfClusters <- max(clusteredpresdata$cluster) 
ratio <- NumberOfClusters/1
numberofpresences <- nrow(clusteredpresdata) 
crs(back_raster) <- myproj
back <- randomPoints(back_raster, 1*numberofpresences) %>%
  as_data_frame()
coordinates(back) <- ~ x + y
back_swd  <- raster::extract(crn, back) %>% 
  cbind(coordinates(back), .)
nrow(back_swd) == nrow(back_swd[complete.cases(back_swd),])
dir.create('../rf/input/points/run1', recursive = TRUE)
write.csv(back_swd, '../rf/input/points/run1/back_swd.csv', row.names = FALSE)
write.csv(occ, '../rf/input/points/run1/occ_swd.csv', row.names = FALSE)

# Cluster analysis to pseudabsences
# bckclust <- rf.clust(occ = back_swd[,3:nrow(back_swd)], nforest = 50, ntrees = 500, nVars = 8, nclasses = 2)
datRF <- as.data.frame(back_swd[,3:ncol(back_swd)])
attach(datRF)
no.forests <- 50#raw = 25
no.trees <- 500
distRF <- RFdist(datRF, mtry1 = 8, no.trees, no.forests, addcl1 = T, addcl2 = F, imp =T, oob.prox1 = T)
no.absenceclasses <- 2
labelRF <- pamNew(distRF$cl1, no.absenceclasses)
detach(datRF)
classdata <- cbind(pb = as.factor(labelRF), back_swd[,3:ncol(back_swd)])
presvalue_swd  <- clusteredpresdata[,3:ncol(clusteredpresdata)] %>%
  cbind(pb = (clusteredpresdata$cluster + no.absenceclasses), .) %>%
  na.omit() %>%
  as.data.frame() %>%
  mutate(cluster = cluster + no.absenceclasses)
presvalue_swd <- dplyr::select(presvalue_swd, pb, bio_1:bio?19)
presvalue_swd <- mutate(presvalue_swd, pb = as.factor(pb))
classdata_2 <- cbind(pb = as.data.frame(classdata)$pb, classdata[,2:ncol(classdata)]) # Background
presvalue_swd <- presvalue_swd %>% dplyr::select(-cluster)
dim(classdata_2); dim(presvalue_swd)
head(classdata_2)
head(presvalue_swd)

allclasses_swd <- rbind(classdata_2, presvalue_swd[,1:ncol(classdata_2)])
unique(allclasses_swd$pb)
write.csv(allclasses_swd, '../rf/input/points/run1/all_classes_swd.csv', row.names = FALSE)

# To make the random forest analysis --------------------------------------
vrs <- colnames(classdata) %>% grep('bio', ., value = TRUE)
vrs <- gsub('.asc', '', vrs) 
vrs <- gsub('\\$', '', vrs)
model1 <- as.formula(paste('factor(pb) ~', paste(paste(vrs), collapse = '+', sep =' ')))
rflist <- vector('list', 50) 
auc <- vector('list', 50)

for(repe in 1:50){ # 50 bosques
  
  print(repe)
  pressample <- list()
  
  for (i in 1:(NumberOfClusters+no.absenceclasses)){
    
    if(any(i==c(1:no.absenceclasses))) { 
      
      rows <- sample(rownames(allclasses_swd[allclasses_swd$pb==i,]), 
                     size = samplesize*NumberOfClusters/2/no.absenceclasses)
    } else {
      rows <- sample(rownames(allclasses_swd[allclasses_swd$pb==i,]), size=samplesize)
    }
    pressample[[i]] <- allclasses_swd[rows,] 
  }
  
  species <- na.omit(do.call(rbind, pressample)) 
  head(species)
  Samplesplit <- sample(rownames(species)) 
  
  envtrain <- species[Samplesplit[1:(0.8*nrow(species))],] 
  envtest <- species[Samplesplit[(0.8*nrow(species)):nrow(species)],] 
  
  rfmodel <- randomForest(model1, data = envtrain, ntree = 500, na.action = na.omit, nodesize = 2) 
  
  save(rfmodel, file = paste('../rf/output/run1/models/', NumberOfClusters, 'Prob_' , 'rep_' ,repe, '.rdata' ,sep=''))
  rflist[[repe]] <- rfmodel
  
  # AUC 
  predicted <- as.numeric(predict(rfmodel, envtest))
  observed <- as.vector(envtest[,'pb'])
  auc[[repe]] <- auc(observed, predicted) 
  rm(rfmodel)
  
  cat(auc[[repe]] ,'\n')
  
}

auc <- unlist(auc)
rff <- do.call(randomForest::combine, rflist)
importance <- as.data.frame(rff$importance)

run <- 'run1'
save(rflist, file = paste('../rData/', run, '/rflist_', NumberOfClusters, '.rdata', sep = ''))
save(importance, file = paste0('../rData/', run, '/importanceRF.rData'))
save(auc, file = paste0('../rData/', run, '/aucRF_dist.rData'))
save(rff, file = paste0('../rData/', run, '/rff_dist.rData'))

# Predict modell
climatevalues  <- data.frame(getValues(crn))
NumberOfClusters <- 3

rasterProbs <- predict(rff, climatevalues, type = 'prob') # proximity = T
rasterProbs_na <- na.omit(rasterProbs)
sum_rasterProbs_na <- apply(rasterProbs_na, 1, sum)

rasterRF <- rowSums(rasterProbs[,c(3:(NumberOfClusters+2))])
uncertainty <- apply(rasterProbs, 1, max)  

rasterRFprob <- crn[[1]]
values(rasterRFprob) <- rasterRF 

rasterRFuncertainty <- crn[[1]]
values(rasterRFuncertainty) <- uncertainty 

rasterRF <- max.col(rasterProbs, 'first')
rasterRFclass <- crn[[1]]
values(rasterRFclass) <- rasterRF

plot(rasterRFprob)

writeRaster(rasterRFclass, paste0('../rf/output/run1/results/raw/RF_5Clust_current.asc'), format = 'ascii', overwrite = T)
writeRaster(rasterRFprob, paste0('../rf/output/run1/results/raw/RF_5Prob_current.asc'), format = 'ascii', overwrite = T)
writeRaster(rasterRFuncertainty, paste0('../rf/output/run1/results/raw/RF_5Unc_current.asc'), format = 'ascii', overwrite = T)


fls <- list.files('../data/climate/future/2040_2069', full.names = TRUE, pattern = '.asc$')
prec <- grep('prec', fls, value = TRUE) %>% stack()
tmax <- grep('tmax', fls, value = TRUE) %>% stack()
tmea <- grep('tmean', fls, value = TRUE) %>% stack()
tmin <- grep('tmin', fls, value = TRUE) %>% stack()
ftr <- biovars(prec, tmin, tmax)

climatevalues <- data.frame(getValues(ftr))
rasterProbs <- predict(rff, climatevalues, type = 'prob') # proximity = T
rasterRF <- rowSums(rasterProbs[,c(3:(NumberOfClusters+2))])
uncertainty <- apply(rasterProbs, 1, max)  
rasterRFprob <- ftr[[1]]
values(rasterRFprob) <- rasterRF 
rasterRFuncertainty <- ftr[[1]]
values(rasterRFuncertainty) <- uncertainty 
rasterRF <- max.col(rasterProbs, 'first')
rasterRFclass <- ftr[[1]]
values(rasterRFclass) <- rasterRF
dir.create('../rf/output/run1/results/raw/2050s')
writeRaster(rasterRFclass, paste0('../rf/output/run1/results/raw/2050s/',  '/RF_', NumberOfClusters, 'Clust_', '2050s.asc', sep=''),  format = 'ascii', overwrite = T)
writeRaster(rasterRFprob, paste0('../rf/output/run1/results/raw/2050s/', '/RF_', NumberOfClusters, 'Prob_',  '2050s.asc', sep=''),  format = 'ascii', overwrite = T)
writeRaster(rasterRFuncertainty, paste('../rf/output/run1/results/raw/2050s/', '/RF_', NumberOfClusters, 'Unc_',  '2050s.asc', sep=''),  format = 'ascii', overwrite = T)

