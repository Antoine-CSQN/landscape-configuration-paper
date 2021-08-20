# This script computes the values of the index for the all mesoscale catchment.


# Select mesoscale catchment / area of study
catchement <- c("YV1")

## List of coeficient of the power law we want to test
coef.list.facc = seq(0.00, 2.5, 0.1)
coef.list.fls = seq(0.00, 4, 0.1)

#Select source 
sources.name <- c("all_2010_2019",
                  "bh_2010_2019",
                  "ca_2010_2019",
                  "prairie_dont_bh_2010_2019",
                  "prairie_sans_bh_2010_2019",
                  "sau_sans_bh_2010_2019",
                  "sau_dont_bh_2010_2019"
) 
sources.type <- rep("long-term", 7)
sources.year <- rep(2010, 7)
sources.nut <- c("all",
                 "bh",
                 "ca",
                 "prairie_dont_bh",
                 "prairie_sans_bh",
                 "sau_sans_bh",
                 "sau_dont_bh") 
sources <- cbind.data.frame(name = sources.name, type = sources.type, nut = sources.nut, year = sources.year)

# Names of facc and fls topo layer
facc.name <- "faccMFDburnRS"
fls.name <- "FLenRSd8_taudem_corrected"

## reading flow acc and downstream slope raster
## ex : topo1 = facc et topo2 =fles ou 
facc <- raster(paste0("raster/", facc.name, ".tif"))
fls <- raster(paste0("raster/", fls.name, ".tif"))

# Set roads and streams to NA
road.streams.na <- raster("raster/roadstrmNAother1.tif") 
facc <- facc * road.streams.na * 100 # Conversion to m?
fls <- fls * road.streams.na

## Function to read source raster 
ReadSourceRaster <- function(sourceName){
  sourcePath <- paste0(work.dir,"/raster/source/",sourceName,".tif")
  sourceRaster <- raster(sourcePath)
  return(sourceRaster)
}

## Function to read catchement delienation raster 
ReadCatRaster <- function(catName){
  return(raster(paste0(wat.folder,"/wat_",catName, ".tif")))
}

## function to trim rasters with a mask (usually catchement delianation)
# x is any raster that we want to mask and trim according to maskRaster
TrimRaster <- function(x, maskRaster){
  maskRaster <- maskRaster %in% c(NA) #create raster mask, NA are not in the chatchment
  x[is.na(x)] <- 0 # change NA in 0 in the source raster
  maskedRaster <- mask(x, maskRaster , maskvalue = 1) #use the mask
  trimmedRaster <- trim(maskedRaster) #crop the raster (less memory)
  return(trimmedRaster)
}


TrimRasterYV1 <- function(x, maskRaster){
  maskRaster <- maskRaster %in% c(0) #create raster mask (only YV1 have 0 for masks)
  maskedRaster <- mask(x, maskRaster , maskvalue = 1) #use the mask
  trimmedRaster <- trim(maskedRaster) #crop the raster (less memory)
  return(trimmedRaster)
}

# compute raster index -> returns a raster
ComputeRasterIndex <- function(faccRst, faccCoef, flsRst, flsCoef, sourceRst){
  return((faccRst^faccCoef) * (flsRst^(-flsCoef)) * sourceRst)
}

# Compute numerator or denominator of the index -> returns a number
ComputeMeanIndex <- function(rasterIndex){
  return( cellStats(rasterIndex, "sum", na.rm = TRUE) / sum(!is.na(values(rasterIndex))) )
}

# return total weight of top pc % 
# Exemple : TotalWeightFirstXPc(rasterIndex, 0.05) = 0619 : 5% des pixels les plus élevés contribuent à 62% du total des poids
TotalWeightFirstXPc <- function(rasterIndex, percentile){
  rasterIndex <- rasterIndex[!is.na(rasterIndex)]
  maxCells <- as.integer(percentile * length(rasterIndex)) #fix the number of cells to avoid more than percentile cells in case od facc = 0
  percentileValue <- as.numeric(quantile(rasterIndex, probs = 1-percentile))
  maskPercent <- rasterIndex >= percentileValue
  rasterIndexSupPercentile <- rasterIndex[maskPercent]
  sumSupPercentile <- sum(rasterIndexSupPercentile[1:maxCells])
  sumTot <- sum(rasterIndex)
  return(sumSupPercentile / sumTot)
}

# Select YV1 cells
facc.rst <- TrimRasterYV1(facc, ReadCatRaster("YV1"))
fls.rst <- TrimRasterYV1(fls, ReadCatRaster("YV1"))


# Prepare table for preprocessing results
preprocess.results <- data.frame(facc.name = rep(NA, length(coef.list.facc)*length(coef.list.fls)*length(sources.name)),
                                 fls.name = NA,
                                 source.name = NA, source.nut = NA, source.type = NA, source.year = NA, 
                                 coef.facc = NA_real_ , coef.fls = NA_real_,
                                 mean.index = NA_real_, pc.05 = NA_real_, pc.01 = NA_real_)

i <- 1
s <- 1
tic(msg = "Total time")
for(source.name in sources.name){
  tic(msg = "Time for a lande use")
  source.rst <- raster(paste0("raster/source/", source.name, ".tif"))
  source.rst[is.na(source.rst)] <- 0 # 0 codes no source while NA codes out of the zone and roads/streams
  source.rst <- source.rst * road.streams.na # burn stream/ rods with NAs
  source.rst <- TrimRasterYV1(source.rst, ReadCatRaster("YV1")) # Cut YV1

  for(coef.facc in coef.list.facc){
    print(coef.facc)
    tic(msg = "Loop on fls")
    for(coef.fls in coef.list.fls){
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
      preprocess.results$pc.01[i] <- TotalWeightFirstXPc(rasterIndex = raster.index, percentile = 0.01)

      i <- i + 1
    }# end loop on coef.fls
    toc() 
  }# end loop on coef.facc
  toc() 
  s <- s + 1 # s is the index for source number
} # end loop on sources
toc()

# Lets say that the 5% heighest weigthed pixels cannot weight more than 95%
acceptable.coefs <- preprocess.results %>% filter(pc.05<0.95) %>% filter(source.name == "all_2010_2019")

# plot percentage of weights explained by top 5%
ggplot(data = preprocess.results, aes(x = coef.fls, y = coef.facc, fill = pc.05))+
  geom_tile() + 
  scale_fill_viridis_b(breaks = c(0.25,0.5,0.75,0.9,0.95,0.99), option = "inferno", trans = "exp") +
  scale_x_continuous(limits=c(-0.05, max(coef.list.fls)), expand = c(0, -0.05)) + # 0.05 is half a tile (0.1)
  scale_y_continuous(limits=c(-0.05, max(coef.list.facc)), expand = c(0, -0.05)) +
  coord_equal()+
  theme_bw()
ggsave(file = paste0(results.folder,"/pc05_allcoefs.png"))

# plot percentage of weights explained by top 5%, without elevated values
ggplot(data = acceptable.coefs, aes(x = coef.fls, y = coef.facc, fill = pc.05))+
  geom_tile() + 
  scale_fill_viridis_b(breaks = c(0.2,0.4,0.6,0.8,0.9), option = "inferno") +
  scale_x_continuous(limits=c(-0.05, max(acceptable.coefs$coef.fls)), expand = c(0, -0.05)) +
  scale_y_continuous(limits=c(-0.05, max(acceptable.coefs$coef.facc)), expand = c(0, -0.05)) +
  coord_equal() +
  theme_bw()
ggsave(file = paste0(results.folder,"/pc05_acceptablecoefs.png"))


save.image(file = paste0(results.folder,"/1_YV1preprocessdata.RData"))

