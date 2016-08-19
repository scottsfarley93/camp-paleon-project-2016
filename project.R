library(raster)
library(neotoma)
library(rgdal)
library(fields)

setwd("/Users/scottsfarley/documents/paleon/project")
#r <- raster("/Users/scottsfarley/documents/paleon/project/LC_5min_global_2010.tif")


wiN <- 47.080799
wiE <- -86.805183
wiS <- 42.49192
wiW <- -92.88959

datasets <- get_dataset(loc=c(wiW, wiS, wiE, wiN), datasettype="pollen surface sample")

downloads <- get_download(datasets)

datasetID = 3802 ##Devil's Lake Modern Pollen Surface Sample, Webb, T 1970
siteID = 666

tsDatasetID <- 684 ##Devil's Lake Core, Maher, Louis J., Jr. 1980
tsDataset <- get_download(tsDatasetID)$`684`
dataset <- get_download(3802)$`3802`
counts <- dataset$counts


# wiBounds <- readOGR(dsn=".", layer = "WI_state_outline")
# wiBounds <- as(wi, "SpatialPolygonsDataFrame")
e <- extent(-92.88959, -86.805, 42.49192, 47.080799)
wiBounds <- as(e, "SpatialPolygons")

broadleaf <- raster("na-latlong-broadleaf.grd")
needleleaf <-raster("na-latlong-needleleaf.grd")
evergreen <- raster("na-latlong-evergreen.grd")
# deciduous <- raster("na-latlong-deciduous.grd")
tree.cover <- raster("na-latlong-treecover.grd")
broadleaf <- crop(broadleaf, extent(wiBounds))
broadleaf <- mask(broadleaf, wiBounds)
# 
needleleaf <- crop(needleleaf, extent(wiBounds))
needleleaf <- mask(needleleaf, wiBounds)
#
evergreen <- crop(evergreen, extent(wiBounds))
evergreen <- mask(evergreen, wiBounds)

deciduous <- crop(deciduous, extent(wiBounds))
deciduous <- mask(deciduous, wiBounds)
#
tree.cover <- crop(tree.cover, extent(wiBounds))
tree.cover <- mask(tree.cover, wiBounds)

writeRaster(needleleaf, "wi-needleleaf.grd", overwrite=TRUE)
writeRaster(broadleaf, "wi-broadleaf.grd", overwrite=TRUE)
writeRaster(evergreen, 'wi-evergreen.grd', overwrite=TRUE)
writeRaster(deciduous, 'wi-deciduous.grd', overwrite=TRUE)
writeRaster(tree.cover, 'wi-treecover.grd', overwrite=TRUE)



composition <- stack(needleleaf, broadleaf, evergreen, tree.cover)

tree.cover[tree.cover == 254] = 0
tree.cover[tree.cover == 255] = 0

tree.cover.points <- as(tree.cover, "SpatialPointsDataFrame")
needlelaf.points <- as(needleleaf, "SpatialPointsDataFrame")
broadleaf.points <- as(broadleaf, "SpatialPointsDataFrame")
evergreen.points <- as(evergreen, "SpatialPointsDataFrame")
deciduous.points <- as(deciduous, "SpatialPointsDataFrame")


siteLat <- dataset['lat']
siteLng <- dataset['long']

modernCompiled <- data.frame(compile_taxa(dataset, "P25")$counts)
tsCompiled <- data.frame(compile_taxa(tsDataset, "P25")$counts)

modernCompiled$total <- rowSums(modernCompiled)
tsCompiled$total <- rowSums(tsCompiled)

modernCompiled <- modernCompiled/modernCompiled$total
tsCompiled <- tsCompiled/tsCompiled$total

aggregateToTypes <- function(df){
  df$broadleaf <- 0
  df$needleleaf <- 0
  df$nonArboreal <- 0
  dfNames <- names(df)
  ##Williams and Jackson, 2003
  broadleafTypes = c("Acer", "Alnus", "Aquifoliaceae", "Betula", "Carya",
                     "Castanea", "Ceanothus", "Celtis", "Cephalanthus",
                     "Cercocarpus", "Clethra", "Corylus", "Empetrum",
                     "Ericaceae", "Fagus", "Fraxinus", "Juglans",
                     "Liquidambar", "Mimosa", "Myricaceae", "Nyssa",
                     "Ostrya", "Carpinus", "Ostrya.Carpinus", "Ostrya", "Carpinus", "Platanus", "Populus", "Quercus",
                     "Salix", "Shepherdia", "Tilia", "Ulmus")
  needleleafTypes <- c("Abies", "Cupressaceae.Taxaceae",  "Cupressaceae", "Taxaceae", "Picea", "Pinus",
                       "Taxodium", "Thuja", "Tsuga","Larix.Pseudotsuga", "Larix", "Pseudotsuga")
  nonArborealTypes <- c("Ambrosia", "Apiaceae", "Artemisia", "Asteraceae", 
                        "Brassicaceae", "Caryophyllaceae",
                        "Chenopodiaceae.Amaranthaceae", "Chenopodiaceae", "Amaranthaceae", "Cyperaceae",
                        "Dryas", "Ephedra", "Eriogonum", "Euphorbiaceae",
                        "Oxyria", "Poaceae", "Ranunculaceae", "Rubiaceae",
                        "Sarcobatus", "Saxifragaceae", "Sphaeralcea", "Prairie.Forbs")
  for (type in broadleafTypes){
    if(type %in% dfNames){
      df['broadleaf'] = df['broadleaf'] + df[type] 
      df[type] <- NULL
    }
  }
  for (type in needleleafTypes){
    if(type %in% dfNames){
      df['needleleaf'] = df['needleleaf'] + df[type]
      df[type] <- NULL
    }
  }
  for(type in nonArborealTypes){
    if(type %in% dfNames){
      df['nonArboreal'] = df['nonArboreal'] + df[type]
      df[type] <- NULL
    }
  }
  df$Other <- df$total - (df$nonArboreal + df$needleleaf + df$broadleaf)
  df$tree.sum <- df$needleleaf + df$broadleaf
  df$total <- NULL
  return(df)
}
# 
# tsAggregated <- aggregateToTypes(tsCompiled)
# modernAggregated <- aggregateToTypes(modernCompiled)
# 
# tsAges <- data.frame(tsDataset$sample.meta)$age
# tsAggregated$ages <- tsAges
# 
# plot(tsAggregated$nonArboreal~ tsAggregated$ages, type='l', col='green', ylim=c(0,1))
# lines(tsAggregated$broadleaf~ tsAggregated$ages, type='l', col='red', ylim=c(0,1))
# lines(tsAggregated$needleleaf~  tsAggregated$ages, type='l', col='purple', ylim=c(0,1))
# lines(tsAggregated$Other~ tsAggregated$ages, type='l', col='pink', ylim=c(0,1))
# legend("top", "", c("NonArboreal", "Broadleaf", "Needleleaf", "Other"), fill=c("green", "red",
#                                                                                       "purple", "pink"))
# 
# abline(modernAggregated$Other, 0, col='pink')
# abline(modernAggregated$needleleaf, 0, col='purple')
# abline(modernAggregated$broadleaf, 0, col='red')
# abline(modernAggregated$nonArboreal, 0, col='green')
# 
tree.cover.points.2 <- data.frame(rasterToPoints(tree.cover))
broadleaf.points.2 <- data.frame(rasterToPoints(broadleaf))
needleleaf.points.2 <- data.frame(rasterToPoints(needleleaf))

# tree.cover.locRound <- tree.cover.points.2
# broadleaf.locRound <- broadleaf.points.2
# needleleaf.locRound <- needleleaf.points.2
# 
# decPlaces = 1
# 
# tree.cover.locRound$x <- round(tree.cover.locRound$x, decPlaces)
# tree.cover.locRound$y <- round(tree.cover.locRound$y, decPlaces)
# broadleaf.locRound$x <- round(broadleaf.locRound$x, decPlaces)
# broadleaf.locRound$y <- round(broadleaf.locRound$y, decPlaces)
# needleleaf.locRound$x <- round(needleleaf.locRound$x, decPlaces)
# needleleaf.locRound$y <- round(needleleaf.locRound$y, decPlaces)
# 
# theCounts <- list()
# theCompiled <- list()
# theAggregated <- list()
# modernComp <- list()
# broadleaf_out <- data.frame(composition=vector(length=length(downloads)), pollen=vector(length=length(downloads)), numGridcells = vector(length=length(downloads)), siteName = vector(length=length(downloads)))
# tree.cover_out <- data.frame(composition=vector(length=length(downloads)), pollen=vector(length=length(downloads)), numGridcells = vector(length=length(downloads)), siteName = vector(length=length(downloads)))
# needleleaf_out <- data.frame(composition=vector(length=length(downloads)), pollen=vector(length=length(downloads)), numGridcells = vector(length=length(downloads)), siteName = vector(length=length(downloads)))
# for (ind in 1:length(downloads{1:10})){
#   download <- downloads[[ind]]
#   counts <- download$counts
#   theCounts[[ind]] <- counts
#   compiled <- data.frame(compile_taxa(counts, "P25"))
#   compiled$total <- rowSums(compiled)
#   compiled <- compiled/compiled$total
#   theCompiled[[ind]] <- compiled
#   aggregated <- aggregateToTypes(compiled)
#   theAggregated[[ind]] <- aggregated
#   pollenTreeCover <- aggregated$tree.sum
#   pollenBroadleaf <- aggregated$broadleaf
#   pollenNeedleleaf <- aggregated$needleleaf
#   latRound <- round(download$dataset$site.data$lat, decPlaces)
#   lngRound <- round(download$dataset$site.data$long, decPlaces)
#   lat <- download$dataset$site.data$lat
#   lng <- download$dataset$site.data$long
#   point <- as.matrix(data.frame(x=lng, y=lat))
#   gridpoints <- as.matrix(tree.cover.points.2[c("x", "y")])
#   near.tree.cover <- tree.cover.locRound[which(tree.cover.locRound$x == lngRound),]
#   near.tree.cover <- near.tree.cover[which(near.tree.cover$y == latRound),]
#   near.needleleaf <- needleleaf.locRound[which(needleleaf.locRound$x == lngRound),]
#   near.needleleaf <- near.needleleaf[which(near.needleleaf$y == latRound), ]
#   near.broadleaf <- broadleaf.locRound[which(broadleaf.locRound$x == lngRound),]
#   near.broadleaf <- near.broadleaf[which(near.broadleaf$y == latRound),]
#   avgTreeCover<-mean(near.tree.cover$na.latlong.treecover)/100
#   avgNeedleleaf <- mean(near.needleleaf$na.latlong.needleleaf)/100
#   avgBroadleaf <- mean(near.broadleaf$na.latlong.broadleaf)/100
#   broadleaf_out$pollen[ind] <- pollenBroadleaf
#   broadleaf_out$composition[ind] <- avgBroadleaf
#   needleleaf_out$pollen[ind] <- pollenNeedleleaf
#   needleleaf_out$composition[ind] <- avgNeedleleaf
#   tree.cover_out$pollen[ind] <- pollenTreeCover
#   tree.cover_out$composition[ind] <- avgTreeCover
#   numGridpoints <- nrow(near.needleleaf)
#   tree.cover_out$numGridcells[ind] <- numGridpoints
#   needleleaf_out$numGridcells[ind] <- numGridpoints
#   broadleaf_out$numGridcells[ind] <- numGridpoints
#   
#   tree.cover_out$siteName[ind] <- download$dataset$site.data$site.name
#   needleleaf_out$siteName[ind] <- download$dataset$site.data$site.name
#   broadleaf_out$siteName[ind] <- download$dataset$site.data$site.name
# }
# 
# 
# par(mfrow=c(3, 1))
# plot(tree.cover_out$pollen ~ tree.cover_out$composition, xlab=paste("AVHRR"), ylab="Pollen", main="Tree Cover", xlim=c(0,1), ylim=c(0,1))
# abline(lm(tree.cover_out$pollen~tree.cover_out$composition))
# plot(needleleaf_out$pollen ~ needleleaf_out$composition, xlab=paste("AVHRR"), ylab="Pollen", main="Needleleaf", xlim=c(0,1), ylim=c(0,1))
# abline(lm(needleleaf_out$pollen~needleleaf_out$composition))
# plot(broadleaf_out$pollen ~ broadleaf_out$composition, xlab=paste("AVHRR"), ylab="Pollen", main="Broadleaf", xlim=c(0,1), ylim=c(0,1))
# abline(lm(broadleaf_out$pollen~broadleaf_out$composition))
# hist(broadleaf_out$numGridcells)
# 

sourceAreas <- c(100)
broadleaf_out <- data.frame(idw=vector(length=length(downloads)),idw2=vector(length=length(downloads)),simple=vector(length=length(downloads)), pollen=vector(length=length(downloads)), numGridcells = vector(length=length(downloads)), siteName = vector(length=length(downloads)),  SA=vector(length=length(downloads)))
tree.cover_out <- data.frame(idw=vector(length=length(downloads)), idw2=vector(length=length(downloads)),simple=vector(length=length(downloads)),pollen=vector(length=length(downloads)), numGridcells = vector(length=length(downloads)), siteName = vector(length=length(downloads)), SA=vector(length=length(downloads)))
needleleaf_out <- data.frame(idw=vector(length=length(downloads)), idw2=vector(length=length(downloads)),simple=vector(length=length(downloads)),pollen=vector(length=length(downloads)), numGridcells = vector(length=length(downloads)), siteName = vector(length=length(downloads)), SA=vector(length=length(downloads)))

overAllInd <- 0
SAInd <- 0
for (area in sourceAreas){
  SAInd <- SAInd + 1
  for (ind in 1:length(downloads)){
    dataset <- downloads[[ind]]
    lat <- dataset$dataset$site.data$lat
    lng <- dataset$dataset$site.data$long
    siteName <- dataset$dataset$site.data$site.name
    point <- as.matrix(data.frame(x=lng, y=lat))
    gridpoints <- as.matrix(tree.cover.points.2[c("x", "y")])
    counts <- dataset$counts
    theCounts[[ind]] <- counts
    compiled <- data.frame(compile_taxa(counts, "P25"))
    compiled$total <- rowSums(compiled)
    #compiled <- compiled/compiled$total
    theCompiled[[ind]] <- compiled
    aggregated <- aggregateToTypes(compiled)
    theAggregated[[ind]] <- aggregated
    pollenTreeCover <- aggregated$tree.sum
    pollenBroadleaf <- aggregated$broadleaf
    pollenNeedleleaf <- aggregated$needleleaf
    dist <- rdist.earth(point, gridpoints)
    matches <- which(dist <= area)
    distMatches <- dist[matches]
    tcMatches <- tree.cover.points.2[matches,]
    broadleafMatches <- broadleaf.points.2[matches,]
    needleleafMatches <- needleleaf.points.2[matches,]
    tree.cover.mesh <- data.frame(x=tcMatches$x, y=tcMatches$y, gridValue=tcMatches$na.latlong.treecover, dist=distMatches)
    broadleaf.mesh <- data.frame(x=tcMatches$x, y=tcMatches$y, gridValue=broadleafMatches$na.latlong.broadleaf, dist=distMatches)
    needleleaf.mesh <- data.frame(x=tcMatches$x, y=tcMatches$y, gridValue=needleleafMatches$na.latlong.needleleaf, dist=distMatches)
    ##do simple averaging
    tree.cover.simple <- mean(tree.cover.mesh$gridValue) / 100
    broadleaf.simple <- mean(broadleaf.mesh$gridValue) / 100
    needleleaf.simple <- mean(needleleaf.mesh$gridValue) / 100
    ##do IDW
    tree.cover.idw <- sum(tree.cover.mesh$gridValue/tree.cover.mesh$dist)/sum(1/tree.cover.mesh$dist)/100
    broadleaf.idw <- sum(broadleaf.mesh$gridValue/broadleaf.mesh$dist)/sum(1/broadleaf.mesh$dist)/100
    needleleaf.idw <- sum(needleleaf.mesh$gridValue/needleleaf.mesh$dist)/sum(1/needleleaf.mesh$dist)/100
    
    ##do IDW^2
    tree.cover.idw2 <- sum(tree.cover.mesh$gridValue/(tree.cover.mesh$dist^2))/sum((1/tree.cover.mesh$dist^2))/100
    broadleaf.idw2 <- sum(broadleaf.mesh$gridValue/(broadleaf.mesh$dist^2))/sum((1/broadleaf.mesh$dist^2))/100
    needleleaf.idw2 <- sum(needleleaf.mesh$gridValue/(needleleaf.mesh$dist^2))/sum((1/needleleaf.mesh$dist^2))/100
    
    
    numCells <- length(matches)
    # dfInd <- SAInd * ind
    # print(dfInd)
    overAllInd <- overAllInd + 1
    ## record the output
    treeV <- c(tree.cover.idw, tree.cover.idw2, tree.cover.simple, pollenTreeCover, numCells, siteName, area)
    broadleafV <- c(broadleaf.idw, broadleaf.idw2, broadleaf.simple, pollenBroadleaf, numCells, siteName, area)
    needleleafV <- c(needleleaf.idw, needleleaf.idw2, needleleaf.simple, pollenNeedleleaf, numCells, siteName, area)
    print(overAllInd)
    tree.cover_out[overAllInd, ] <- treeV
    broadleaf_out[overAllInd, ] <- broadleafV
    needleleaf_out[overAllInd, ] <- needleleafV
    }
}

tree.cover.25 <- tree.cover_out[which(tree.cover_out$SA == 25),]
tree.cover.50 <- tree.cover_out[which(tree.cover_out$SA == 50), ]
tree.cover.75 <- tree.cover_out[which(tree.cover_out$SA == 75), ]
tree.cover.100 <- tree.cover_out[which(tree.cover_out$SA == 100), ]
tree.cover.150 <- tree.cover_out[which(tree.cover_out$SA == 150), ]
tree.cover.200 <- tree.cover_out[which(tree.cover_out$SA == 200), ]

broadleaf.25 <- broadleaf_out[which(broadleaf_out$SA == 25), ]
broadleaf.50 <- broadleaf_out[which(broadleaf_out$SA == 50), ]
broadleaf.75 <- broadleaf_out[which(broadleaf_out$SA == 75), ]
broadleaf.100 <- broadleaf_out[which(broadleaf_out$SA == 100), ]
broadleaf.150 <- broadleaf_out[which(broadleaf_out$SA == 150), ]
broadleaf.200 <- broadleaf_out[which(broadleaf_out$SA == 200), ]

needleleaf.25 <- needleleaf_out[which(needleleaf_out$SA == 25), ]
needleleaf.50 <- needleleaf_out[which(needleleaf_out$SA == 50), ]
needleleaf.75 <- needleleaf_out[which(needleleaf_out$SA == 75), ]
needleleaf.100 <- needleleaf_out[which(needleleaf_out$SA == 100), ]
needleleaf.150 <- needleleaf_out[which(needleleaf_out$SA == 150), ]
needleleaf.200 <- needleleaf_out[which(needleleaf_out$SA == 200), ]

save(broadleaf_out, file="/Users/scottsfarley/documents/paleon/project/broadleaf.RData")
save(needleleaf_out, file="/Users/scottsfarley/documents/paleon/project/needleleaf.RData")
save(tree.cover_out, file="/Users/scottsfarley/documents/paleon/project/treecover.RData")

par(mfrow=c(3, 2))
plot(tree.cover.25$pollen, tree.cover.25$idw, col='red', main='Tree Cover: 25km', xlim=c(0,1), ylim=c(0, 1))
points(tree.cover.25$pollen, tree.cover.25$idw2, col='green')
points(tree.cover.25$pollen, tree.cover.25$simple, col='blue')

plot(tree.cover.50$pollen, tree.cover.50$idw, col='red', main='Tree Cover: 50km', xlim=c(0,1), ylim=c(0, 1))
points(tree.cover.50$pollen, tree.cover.50$idw2, col='green')
points(tree.cover.50$pollen, tree.cover.50$simple, col='blue')

plot(tree.cover.75$pollen, tree.cover.75$idw, col='red', main='Tree Cover: 75km', xlim=c(0,1), ylim=c(0, 1))
points(tree.cover.75$pollen, tree.cover.75$idw2, col='green')
points(tree.cover.75$pollen, tree.cover.75$simple, col='blue')

plot(tree.cover.75$pollen, tree.cover.100$idw, col='red', main='Tree Cover: 100km', xlim=c(0,1), ylim=c(0, 1))
points(tree.cover.75$pollen, tree.cover.100$idw2, col='green')
points(tree.cover.75$pollen, tree.cover.100$simple, col='blue')

plot(tree.cover.200$pollen, tree.cover.200$idw, col='red', main='Tree Cover: 200km', xlim=c(0,1), ylim=c(0, 1))
points(tree.cover.200$pollen, tree.cover.200$idw2, col='green')
points(tree.cover.200$pollen, tree.cover.200$simple, col='blue')

plot(tree.cover, main="Tree Cover")


par(mfrow=c(3, 2))
plot(broadleaf.25$pollen, broadleaf.25$idw, col='red', main='Broadleaf Cover: 25km', xlim=c(0,1), ylim=c(0, 1))
points(broadleaf.25$pollen, broadleaf.25$idw2, col='green')
points(broadleaf.25$pollen, broadleaf.25$simple, col='blue')

plot(broadleaf.50$pollen, broadleaf.50$idw, col='red', main='Broadleaf Cover: 50km', xlim=c(0,1), ylim=c(0, 1))
points(broadleaf.50$pollen, broadleaf.50$idw2, col='green')
points(broadleaf.50$pollen, broadleaf.50$simple, col='blue')

plot(broadleaf.75$pollen, broadleaf.75$idw, col='red', main='Broadleaf Cover: 75km', xlim=c(0,1), ylim=c(0, 1))
points(broadleaf.75$pollen, broadleaf.75$idw2, col='green')
points(broadleaf.75$pollen, broadleaf.75$simple, col='blue')

plot(broadleaf.100$pollen, broadleaf.100$idw, col='red', main='Broadleaf Cover: 100km', xlim=c(0,1), ylim=c(0, 1))
points(broadleaf.100$pollen, broadleaf.100$idw2, col='green')
points(broadleaf.100$pollen, broadleaf.100$simple, col='blue')

plot(broadleaf.200$pollen, broadleaf.200$idw, col='red', main='Broadleaf Cover: 200km', xlim=c(0,1), ylim=c(0, 1))
points(broadleaf.200$pollen, broadleaf.200$idw2, col='green')
points(broadleaf.200$pollen, broadleaf.200$simple, col='blue')
plot(broadleaf, main="Broadleaf")

par(mfrow=c(3, 2))
plot(needleleaf.25$pollen, needleleaf.25$idw, col='red', main='Needleleaf Cover: 25km', xlim=c(0,1), ylim=c(0, 1))
points(needleleaf.25$pollen, needleleaf.25$idw2, col='green')
points(needleleaf.25$pollen, needleleaf.25$simple, col='blue')

plot(needleleaf.50$pollen, needleleaf.50$idw, col='red', main='Needleleaf Cover: 50km', xlim=c(0,1), ylim=c(0, 1))
points(needleleaf.50$pollen, needleleaf.50$idw2, col='green')
points(needleleaf.50$pollen, needleleaf.50$simple, col='blue')

plot(needleleaf.75$pollen, needleleaf.75$idw, col='red', main='Needleleaf Cover: 75km', xlim=c(0,1), ylim=c(0, 1))
points(needleleaf.75$pollen, needleleaf.75$idw2, col='green')
points(needleleaf.75$pollen, needleleaf.75$simple, col='blue')

plot(needleleaf.100$pollen, needleleaf.100$idw, col='red', main='Needleleaf Cover: 100km', xlim=c(0,1), ylim=c(0, 1))
points(needleleaf.100$pollen, needleleaf.100$idw2, col='green')
points(needleleaf.100$pollen, needleleaf.100$simple, col='blue')

plot(needleleaf.200$pollen, needleleaf.200$idw, col='red', main='Needleleaf Cover: 200km', xlim=c(0,1), ylim=c(0, 1))
points(needleleaf.200$pollen, needleleaf.200$idw2, col='green')
points(needleleaf.200$pollen, needleleaf.200$simple, col='blue')
plot(needleleaf, main="Needleleaf")


sourceAreas <- c(100)
out_broadleaf_by_taxon <- list()
out_needleaf_by_taxon <- list()
out_treecover_by_taxon <- list()
for (taxon in names(compiled)){
  df <-  data.frame(idw=vector(length=length(downloads)), idw2=vector(length=length(downloads)),simple=vector(length=length(downloads)),pollen=vector(length=length(downloads)), numGridcells = vector(length=length(downloads)), siteName = vector(length=length(downloads)), SA=vector(length=length(downloads)))
  out_broadleaf_by_taxon[[taxon]] = df
  out_needleaf_by_taxon[[taxon]] = df
  out_treecover_by_taxon[[taxon]] = df
}
overAllInd <- 0
SAInd <- 0
for (area in sourceAreas){
  SAInd <- SAInd + 1
  for (ind in 1:length(downloads)){
    dataset <- downloads[[ind]]
    lat <- dataset$dataset$site.data$lat
    lng <- dataset$dataset$site.data$long
    siteName <- dataset$dataset$site.data$site.name
    point <- as.matrix(data.frame(x=lng, y=lat))
    gridpoints <- as.matrix(tree.cover.points.2[c("x", "y")])
    counts <- dataset$counts
    theCounts[[ind]] <- counts
    compiled <- data.frame(compile_taxa(counts, "P25"))
    compiled$total <- rowSums(compiled)
    compiled <- compiled/compiled$total
    theCompiled[[ind]] <- compiled
    aggregated <- aggregateToTypes(compiled)
    theAggregated[[ind]] <- aggregated
    pollenTreeCover <- aggregated$tree.sum
    pollenBroadleaf <- aggregated$broadleaf
    pollenNeedleleaf <- aggregated$needleleaf
    dist <- rdist.earth(point, gridpoints)
    matches <- which(dist <= area)
    distMatches <- dist[matches]
    tcMatches <- tree.cover.points.2[matches,]
    for (taxon in names(out_treecover_by_taxon)){
      taxonPollen <- compiled[[taxon]]
      if(is.null(taxonPollen)){
        taxonPollen <- NA
      }
      broadleafMatches <- broadleaf.points.2[matches,]
      mesh <- data.frame(x=tcMatches$x, y=tcMatches$y, gridValue=broadleafMatches$na.latlong.broadleaf, dist=distMatches)
      needleleafMatches <- needleleaf.points.2[matches,]
      tree.cover.mesh <- data.frame(x=tcMatches$x, y=tcMatches$y, gridValue=tcMatches$na.latlong.treecover, dist=distMatches)
      broadleaf.mesh <- data.frame(x=tcMatches$x, y=tcMatches$y, gridValue=broadleafMatches$na.latlong.broadleaf, dist=distMatches)
      needleleaf.mesh <- data.frame(x=tcMatches$x, y=tcMatches$y, gridValue=needleleafMatches$na.latlong.needleleaf, dist=distMatches)
      ##do simple averaging
      tree.cover.simple <- mean(tree.cover.mesh$gridValue) / 100
      broadleaf.simple <- mean(broadleaf.mesh$gridValue) / 100
      needleleaf.simple <- mean(needleleaf.mesh$gridValue) / 100
      ##do IDW
      tree.cover.idw <- sum(tree.cover.mesh$gridValue/tree.cover.mesh$dist)/sum(1/tree.cover.mesh$dist)/100
      broadleaf.idw <- sum(broadleaf.mesh$gridValue/broadleaf.mesh$dist)/sum(1/broadleaf.mesh$dist)/100
      needleleaf.idw <- sum(needleleaf.mesh$gridValue/needleleaf.mesh$dist)/sum(1/needleleaf.mesh$dist)/100

      ##do IDW^2
      tree.cover.idw2 <- sum(tree.cover.mesh$gridValue/(tree.cover.mesh$dist^2))/sum((1/tree.cover.mesh$dist^2))/100
      broadleaf.idw2 <- sum(broadleaf.mesh$gridValue/(broadleaf.mesh$dist^2))/sum((1/broadleaf.mesh$dist^2))/100
      needleleaf.idw2 <- sum(needleleaf.mesh$gridValue/(needleleaf.mesh$dist^2))/sum((1/needleleaf.mesh$dist^2))/100


      numCells <- length(matches)
      # dfInd <- SAInd * ind
      # print(dfInd)
      ## record the output
      treeV <- c(tree.cover.idw, tree.cover.idw2, tree.cover.simple, taxonPollen, numCells, siteName, area)
      broadleafV <- c(broadleaf.idw, broadleaf.idw2, broadleaf.simple, taxonPollen, numCells, siteName, area)
      needleleafV <- c(needleleaf.idw, needleleaf.idw2, needleleaf.simple, taxonPollen, numCells, siteName, area)
      print(treeV)
      out_treecover_by_taxon[[taxon]][overAllInd, ] <- treeV
      out_broadleaf_by_taxon[[taxon]][overAllInd, ] <- broadleafV
      out_needleaf_by_taxon[[taxon]][overAllInd, ] <- needleleafV
    }
    overAllInd <- overAllInd + 1
    # 
  }
}



save(out_broadleaf_by_taxon, file="/Users/scottsfarley/documents/paleon/project/broadleaf_full.RData")
save(out_needleaf_by_taxon, file="/Users/scottsfarley/documents/paleon/project/needleleaf_full.RData")
save(out_treecover_by_taxon, file="/Users/scottsfarley/documents/paleon/project/treecover_full.RData")



for (ind in 1:length(names(out_broadleaf_by_taxon))){
  taxon <- names(out_broadleaf_by_taxon)[[ind]]
  print(taxon)
  pdf(paste("BL", taxon, ".pdf", sep=""))
  par(mfrow=c(3, 2))
  allRows <-  out_broadleaf_by_taxon[[taxon]]
  for (area in sourceAreas){
    rows <- out_broadleaf_by_taxon[[taxon]][which(out_broadleaf_by_taxon[[taxon]]$SA == area), ]
    plot(rows$pollen ~ rows$idw, col='red', main=paste(taxon, ":", area, "km"), xlim=c(0,1), ylim=c(0, 1), xlab="Broadleaf AVHRR Comp.", ylab=paste(taxon, "Pollen"))
    points(rows$pollen~ rows$idw2, col='green')
    points(rows$pollen ~ rows$simple, col='blue')
  }
  dev.off()
  print(paste("Done with: ", taxon))
}

library(R2jags)

model <- function(){
  ## inputs
    ## dobs <- vector of distance observations
    ## pObs <- vector(?) of observed pollen
    ## lObs <- vector of observed landcover percentage
    ## numCells <- number of gridcells  
  
  beta1 ~ dunif(0, 100) ## landcover observation error
  beta2 ~ dunif(0, 100) ## spatial uncertainty in landcover locations
  beta3 ~ dunif(0,100) ## pollen dispersal uncertainty
  beta4 ~ dunif(0, 100) ## pollen observation error
  beta5 ~ dunif(0, 100)
  
  ## the model
  for (i in 1:numCells){
    d[i] ~ dnorm(dObs[i], beta2)
    L[i] ~ dnorm(lObs[i], beta1)
    p_cell[i] ~ dnorm(((1/d[i])*L[i])/(1/d[i]), beta5)
  }
  P_total <- sum(p_cell)
  Pobs ~ dnorm(P, beta4)
  P ~ dnorm(P_total, beta3)
}

ind <- 1
area <- 25
dataset <- downloads[[ind]]
lat <- dataset$dataset$site.data$lat
lng <- dataset$dataset$site.data$long
siteName <- dataset$dataset$site.data$site.name
point <- as.matrix(data.frame(x=lng, y=lat))
gridpoints <- as.matrix(tree.cover.points.2[c("x", "y")])
counts <- dataset$counts
compiled <- data.frame(compile_taxa(counts, "P25"))
compiled$total <- rowSums(compiled)
compiled <- compiled/compiled$total
aggregated <- aggregateToTypes(compiled)

dist <- rdist.earth(point, gridpoints)
matches <- which(dist <= area)
distMatches <- as.vector(dist[matches])
tcMatches <- as.data.frame(tree.cover.points.2[matches,])$na.latlong.treecover
pObs <- aggregated$tree.sum
numCells <- length(distMatches)


tree.cover.jags <- jags(data = list(dObs = distMatches, pObs = pObs, numCells = numCells, lObs = tcMatches), 
            parameters.to.save = c('P', 'beta1', 'beta2', 'beta3', 'beta4', 'beta5'), 
            n.chains = 5, n.iter = 10000, n.burnin = 2000, 
            model.file = model, DIC = FALSE)



total <- (as.integer(needleleaf.100$pollen) + as.numeric(broadleaf.100$pollen))

n <- needleleaf.100
n$idw <- as.numeric(n$idw) / (as.numeric(n$idw) + as.numeric(broadleaf.100$idw))
rStar <- as.numeric(n$idw)

nSites <- length(rStar)
y <- as.integer(needleleaf.100$pollen)



model <- function(){
  beta1 ~ dnorm(0, 0.0001)
  beta0 ~ dnorm(0, 0.0001)
  sigma ~ dunif(0.0001, 10000)
  sigma2Inv <- 1/(sigma*sigma)
  for (i in 1:nSites){
    y[i] ~ dbin(p[i], total[i])
    p[i] <- ilogit(r[i])
    r[i] ~ dnorm(rStar[i]*beta1 + beta0, sigma2Inv)
  }
}
comp.jags <- jags(data = list(rStar = rStar, nSites = nSites, y=y, total=total), 
                        parameters.to.save = c("beta1", "beta0", "p"), 
                        n.chains = 1, n.iter = 10000, n.burnin = 1000, 
                        model.file = model, DIC = FALSE)



c1 <- as.mcmc(comp.jags)[[1]]
c2 <- as.mcmc(comp.jags)[[2]]


p <- colMeans(as.data.frame(c1))
p <- p[3:length(p)]

names(p) <- sapply(names(p), function(n){
  print(n)
    new <- gsub("p\\[", "", n)
    new <- gsub("]", "", new)
    return(new)
  })


idx <- sort(as.numeric(names(p)), index.return=T)$ix
p <- p[idx]

plot(c1)
new_names <- character(length=length(names(p)))
for (i in names(p)){
  name <- names(p)[i]
  name <- re
}

rStar2 <- data.frame(as.numeric(broadleaf.100$idw), as.numeric(needleleaf.100$idw))
rStar2[,1] <- rStar2[,1] * mean(c1[, "betaB"])
rStar2[,2] <- rStar2[,2] * mean(c1[, "betaN"]) + mean(c1[, "betaC"])
print(mean(c1[, "betaB"]))
print(mean(c1[, "betaN"]))
title("Model #1: Corrected")
plot(rStar2[,1], broadleaf.100$pollen, xlim=c(0,1), ylim=c(0,1), 
       xlab="Corrected IDW AVHRR", ylab="Pollen", col='red')
points(rStar2[,2], needleleaf.100$pollen, xlim=c(0,1), ylim=c(0,1), 
       xlab="Corrected IDW AVHRR", ylab="Pollen", col='blue')
abline(0,1)
points(needleleaf.100$idw, needleleaf.100$pollen, col='lightblue')
points(broadleaf.100$idw, broadleaf.100$pollen, col='pink')
#legend("topright", "", c("Corrected Broadleaf", "Uncorrected Broadleaf", "Corrected Needleleaf", "Uncorrected Needleleaf"), fill=c("red", "pink", "blue", "lightblue"))


model2 <- function(){
  betaB ~ dnorm(0, 0.0001) ## broadleaf effect
  betaN ~ dnorm(0, 0.0001) ## needleleaf effect
  for (i in 1:nSites){
    betai[i] ~ dunif(0, 10000) ## site level random effect
    r[i, 1] <- rStar[i, 1]  * betaB ## give each site a broadleaf correction
    r[i, 2] <- rStar[i, 2] * betaN ## give each site a needleleaf correction
    ri[i, 1] ~ dlnorm(r[i, 1], 1/betai[i]) ## give the broadleafs a site-level random effect 
    ri[i, 2] ~ dlnorm(r[i, 2], 1/betai[i]) ## give the needleleafs a site-level random effect
    y[i,] ~ ddirch(ri[i,]) ## give the composition of the site a multi-nomial distribution that sums to one
  }
}

comp.jags.2 <- jags(data = list(rStar = rStar, nSites = nSites, y=y), 
                  parameters.to.save = c("betaB", "betaN", "betai"), 
                  n.chains = 1, n.iter = 10000, n.burnin = 2500, n.thin=10,
                  model.file = model2, DIC = FALSE)


c1 <- as.mcmc(comp.jags.2)[[1]]
c2 <- as.mcmc(comp.jags.2)[[2]]

rStar2 <- data.frame(as.numeric(broadleaf.100$idw), as.numeric(needleleaf.100$idw))
rStar2[,1] <- rStar2[,1] * mean(c1[, "betaB"])
rStar2[,2] <- rStar2[,2] * mean(c1[, "betaN"])
title("Model #2: Corrected")
plot(rStar2[,2], broadleaf.100$pollen, xlim=c(0,1), ylim=c(0,1), 
     xlab="Corrected IDW AVHRR", ylab="Pollen", col='red')
points(rStar2[,2], needleleaf.100$pollen, xlim=c(0,1), ylim=c(0,1), 
       xlab="Corrected IDW AVHRR", ylab="Pollen", col='blue')
# abline(0,1)
plot(needleleaf.100$idw, needleleaf.100$pollen, col='pink')
points(broadleaf.100$idw, broadleaf.100$pollen, col='lightblue')
legend("topright", "", c("Corrected Broadleaf", "Uncorrected Broadleaf", "Corrected Needleleaf", "Uncorrected Needleleaf"), fill=c("red", "pink", "blue", "lightblue"))


needleLocs <- data.frame(x=vector(), y=vector(), pollen=vector())
broadLocs <- data.frame(x=vector(), y=vector(), pollen=vector())
for (ind in 1:length(downloads)){
  agg <- theAggregated[[ind]]
  needles <- agg$needleleaf
  broads <- agg$broadleaf
  dataset <- downloads[[ind]]
  lat <- dataset$dataset$site.data$lat
  lng <- dataset$dataset$site.data$long
  print(lat)
  print(lng)
  print(needles)
  vN <- c(lng, lat, needles)
  vB <- c(lng, lat, broads)
  needleLocs[ind,] <- vN
  broadLocs[ind,] <- vB
}

plot(broadleaf)
points(broadLocs$x, broadLocs$y, col=terrain.colors(broadLocs$pollen), pch=15)
points(broadLocs$x, broadLocs$y, pch=15, cex=0.25)
title("Broadleaf")


plot(needleleaf)
points(needleLocs$x, needleLocs$y, col=terrain.colors(needleLocs$pollen), pch=15)
points(broadLocs$x, broadLocs$y, pch=15, cex=0.25)
title("Needleleaf")
