# Preprocessing of geofiles 

# This script computes the values of the index for the all mesoscale catchment 
# (the study area encompassing the headwater catchments)
# for different source 

## reading Facc and FLS rasters
facc.rst <- raster(paste0("raster/", facc.name, ".tif"))
fls.rst <- raster(paste0("raster/", fls.name, ".tif"))

# Convert study area to raster 
if(!file.exists(paste0(cat.dir,"/", study.area.filename,".tif"))){
  sa.shp <- st_read(paste0(shp.dir,"/",study.area.filename,".shp"))
  sa.rst <- fasterize(sa.shp, facc.rst)
  sa.rst <- trim(sa.rst)
  writeRaster(sa.rst, paste0(cat.dir,"/", study.area.filename,".tif"))
}else {
    sa.rst <- raster(paste0(cat.dir,"/", study.area.filename,".tif"))
  }

# Convert headwater catchments from shp to raster 
hw.shp <- st_read(paste0(shp.dir,"/", hw.filename,".shp"))
hw.shp <- hw.shp %>% filter(.data[[hw.id]] %in% hw.select) %>% arrange(.data[[hw.id]])
hw.number <- nrow(hw.shp)

for(hw.i in seq(1, hw.number)){
  if(!file.exists(paste0(cat.dir,"/cat_", hw.shp$Name[hw.i], ".tif"))){
     hw.rst <- fasterize(hw.shp[hw.i,], facc.rst) * sa.rst
     writeRaster(hw.rst, paste0(cat.dir,"/cat_", hw.shp$Name[hw.i], ".tif"))
   }
}

## Function to read source (LU) raster 
ReadSourceRaster <- function(sourceName){
  sourcePath <- paste0(work.dir,"/raster/source/",sourceName,".tif")
  sourceRaster <- raster(sourcePath)
  return(sourceRaster)
}

## Function to read catchment delineation raster 
ReadCatRaster <- function(catName){
  return(raster(paste0(cat.dir,"/cat_",catName, ".tif")))
}

## function to trim rasters with a mask (usually catchment delineation)
# x is any raster that we want to mask and trim according to maskRaster
TrimRaster <- function(x, maskRaster){
  maskRaster <- maskRaster %in% c(NA) #create raster mask, NA are not in the catchment
  x[is.na(x)] <- 0 # change NA to 0 in the source raster
  maskedRaster <- mask(x, maskRaster , maskvalue = 1) #use the mask
  trimmedRaster <- trim(maskedRaster) #crop the raster (less memory)
  return(trimmedRaster)
}

# compute raster index for each pixel -> returns a raster
ComputeRasterIndex <- function(faccRst, faccCoef, flsRst, flsCoef, sourceRst){
  return(faccRst^faccCoef * flsRst^(-flsCoef) * sourceRst)
}

# Compute numerator or denominator of the LCI -> returns a number
ComputeMeanIndex <- function(rasterIndex){
  rasterIndex[is.infinite(rasterIndex)] <- NA
  return(cellStats(rasterIndex, "sum", na.rm = TRUE) / sum(!is.na(values(rasterIndex))) )
}

# return total weight of top pc % 
# Exemple : TotalWeightFirstXPc(rasterIndex, 0.05) = 0619 : 5% of the most weighted pixels contribute to 61.9% of the total weight of the LCI
TotalWeightFirstXPc <- function(rasterIndex, percentile){
  rasterIndex <- rasterIndex[!is.na(rasterIndex)]
  maxCells <- as.integer(percentile * length(rasterIndex)) #fix the number of cells to avoid more than percentile cells in case of facc = 0
  percentileValue <- as.numeric(quantile(rasterIndex, probs = 1-percentile))
  maskPercent <- rasterIndex >= percentileValue
  rasterIndexSupPercentile <- rasterIndex[maskPercent]
  sumSupPercentile <- sum(rasterIndexSupPercentile[1:maxCells])
  sumTot <- sum(rasterIndex)
  return(sumSupPercentile / sumTot)
}

# Select cells within the study area (Yvel's catchment : sa.rst)
sa.rst <- trim(sa.rst)
facc.rst <- facc.rst * sa.rst
fls.rst <-  fls.rst * sa.rst
ditch.stream.mask <- ditch.stream.mask * sa.rst

# In Facc and FLS rasters, set ditchs and streams to NA and convert Facc to m2 and fls to m
facc.rst <- facc.rst * ditch.stream.mask * facc.to.m2 # Conversion to m2
fls.rst <- fls.rst * ditch.stream.mask * fls.to.m # Conversion to m

# Make the sources data frame: contains info about the source(s) processed
sources <- cbind.data.frame(name = sources.name, type = sources.type, nut = sources.nut, year = sources.year)

# Prepare table of preprocessing results
# Compute : 
#   mean.index: the denominator of the LCI over the study area 
#   pc.05: weight (% of total LCI) of the top 5% most weighted pixels
# Quite long : a bit more than 1sec for 1 (a,b) pair * 1 LU (2.5 GhZ)
# 1600s for a source layer, 3200s for the results of the paper
# Ideas to improve speed if needed : 
#   1) Subsample the rasters (1/100), numerical results almost identical
#   2) parallelize (with package)
#   3) use *apply family function instead of for loops

preprocess.results <- data.frame(facc.name = rep(NA, length(a.facc.values)*length(b.fls.values)*length(sources.name)),
                                 fls.name = NA,
                                 source.name = NA, source.nut = NA, source.type = NA, source.year = NA, 
                                 coef.facc = NA_real_ , coef.fls = NA_real_,
                                 mean.index = NA_real_, pc.05 = NA_real_)

i <- 1
s <- 1
tic(msg = "Computing denominator for source(s) layer(s) for all (a,b)")
for(source.name in sources.name){
  print(paste("Processing", source.name))
  tic(msg = "Time for a source layer for all (a,b)")
  source.rst <- raster(paste0("raster/source/", source.name, ".tif"))
  source.rst[is.na(source.rst)] <- 0 # 0 codes no source while NA codes out of the zone and ditches/streams
  source.rst <- source.rst * ditch.stream.mask # burn streams/ditchs with NAs

  
  for(coef.facc in a.facc.values){
    tic(msg = paste(i, "/", nrow(preprocess.results), "iterations"))
    for(coef.fls in b.fls.values){
      preprocess.results$facc.name[i] <- facc.name
      preprocess.results$fls.name[i] <- fls.name
      preprocess.results$source.name[i] <- sources.name[s]
      preprocess.results$source.nut[i] <- sources.nut[s]
      preprocess.results$source.type[i] <- sources.type[s]
      preprocess.results$source.year[i] <- sources.year[s]
      preprocess.results$coef.facc[i] <- coef.facc
      preprocess.results$coef.fls[i] <- coef.fls
      
      raster.index <- ComputeRasterIndex(faccRst = facc.rst, flsRst = fls.rst,
                                         faccCoef = coef.facc, flsCoef = coef.fls,
                                         sourceRst = source.rst)
      
      preprocess.results$mean.index[i] <- ComputeMeanIndex(raster.index)
      preprocess.results$pc.05[i] <- TotalWeightFirstXPc(rasterIndex = raster.index, percentile = 0.05)
      
      i <- i + 1
    }# end loop on coef.fls
    toc() 
  }# end loop on coef.facc
  toc() 
  s <- s + 1 # s is the index for source number
} # end loop on sources
toc()

# Lets say that the sum of the 5% most weighted pixels cannot weight more than 95% of LCI
acceptable.coefs <- preprocess.results %>% 
  filter(pc.05<0.95) %>% 
  filter(source.name == "all1_ditchstream0")

# # plot percentage of weights explained by top 5%
# ggplot(data = preprocess.results, aes(x = coef.fls, y = coef.facc, fill = pc.05))+
#   geom_tile() +
#   scale_fill_viridis_b(breaks = c(0.25,0.5,0.75,0.9,0.95,0.99), option = "inferno", trans = "exp") +
#   scale_x_continuous(limits=c(-0.05, max(b.fls.values)), expand = c(0, -0.05)) + # 0.05 is half a tile (0.1)
#   scale_y_continuous(limits=c(-0.05, max(a.facc.values)), expand = c(0, -0.05)) +
#   coord_equal()+
#   theme_bw()
# ggsave(file = paste0(results.dir,"/pc05_allcoefs.png"))
# 
# # plot percentage of weights explained by top 5%, without elevated values
# ggplot(data = acceptable.coefs, aes(x = coef.fls, y = coef.facc, fill = pc.05))+
#   geom_tile() +
#   scale_fill_viridis_b(breaks = c(0.2,0.4,0.6,0.8,0.9), option = "inferno") +
#   scale_x_continuous(limits=c(-0.05, max(acceptable.coefs$coef.fls)), expand = c(0, -0.05)) +
#   scale_y_continuous(limits=c(-0.05, max(acceptable.coefs$coef.facc)), expand = c(0, -0.05)) +
#   coord_equal() +
#   theme_bw()
# ggsave(file = paste0(results.dir,"/pc05_acceptablecoefs.png"))


save.image(file = paste0(results.dir,"/1_preprocessdata.RData"))
