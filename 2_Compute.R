# First run 1_Preprocess

#Select headwater catchments : order<=3, no WWTP
list.select <- c("P01", "P02", "P03", "P04", "P05",
                 "P07", "P12", "P13", "P16", "P17",
                 "P18", "P19", "P20", "P22", "P23",
                 "P24", "P25", "P26", "P29")

## Read source rasters & pre_process
tic(name="Preprocessing of subcatchement sources")
subcat.sources <- list()
for(source.name in sources$name){
  source.path <- paste0(work.dir,"/raster/source/",source.name,".tif")
  source.raster <- raster(source.path)
  print(paste0("===== ", source.path, " processsing ======="))
  source.type <- sources$type[which(sources$name == source.name)]
  print(paste0("===== source type: ",source.type, " ======="))
  tic(name = source.name)
  
  for(subcat.name in list.select){
    trimmed.source.rst <- TrimRaster(source.raster, ReadCatRaster(subcat.name))
    source.rst <- trimmed.source.rst * road.streams.na # burn stream/ rods with NAs
    name.tr <- paste0(subcat.name,"_",source.name)
    subcat.sources[[name.tr]] <- trimmed.source.rst
  }# end loop on subcatchements"name.
  toc()
  
}# end loop on sources
toc() # approx 45s par source


## Preprocess fls & facc 
fls.subcat <- list()
facc.subcat <- list()
tic()
for(subcat.name in list.select){
  print(paste0("===== Trimming fls", subcat.name, " ======="))
  wat.fls <- TrimRaster(fls, ReadCatRaster(subcat.name))
  
  print(paste0("===== Trimming facc for ", subcat.name, " ======="))
  wat.facc <- TrimRaster(facc, ReadCatRaster(subcat.name))

  # Force roads and streams = NA
  wat.facc[wat.fls == 0] <- NA 
  wat.fls[wat.fls == 0] <- NA 
  fls.subcat[[subcat.name]] <- wat.fls
  facc.subcat[[subcat.name]] <- wat.facc
}
toc() #53 sec


#Preprocessing of all subcatchment raster index
raster.indexes <- list() 
tic("All raster indexes")

for(subcat.name in list.select){
  tic(paste("Point :", subcat.name))
  for(i in seq(1, nrow(acceptable.coefs))){
    a <- acceptable.coefs$coef.facc[i]
    b <- acceptable.coefs$coef.fls[i]
    
    raster.indexes[[paste0(subcat.name,"_", a, "_", b)]] <-
      facc.subcat[[subcat.name]]^a * fls.subcat[[subcat.name]]^(-b)
    
    # FIX FOR NA^0 = 1
    if((a == 0 & b == 0)){
      temp <- facc.subcat[[subcat.name]]
      temp[!is.na(temp)] <- 1
      raster.indexes[[paste0(subcat.name,"_", 0, "_", 0)]] <- temp
    }
  
  }
  
  toc()
}
toc() # 101 sec




MakeEmptyResDF <- function(acceptableCoefs, subcatName, sourceName, sourceType, sourceNut, sourceYear){
  meanForSource <- preprocess.results %>% filter(source.name == sourceName)
  
  acceptableCoefs <- left_join(acceptableCoefs %>% select(-mean.index),
                                meanForSource %>% select(coef.facc, coef.fls, mean.index),
                               by = c("coef.fls", "coef.facc")) 
  
  emptyResDF <- data.frame(facc.name = acceptableCoefs$facc.name,
                           fls.name = acceptableCoefs$fls.name,
                           subcat.name = rep(subcatName, nrow(acceptableCoefs)),
                           source.name = sourceName,
                           source.type = sourceType,
                           source.nut = sourceNut,
                           source.year = sourceYear,
                           coef.facc = acceptableCoefs$coef.facc,
                           coef.fls = acceptableCoefs$coef.fls,
                           numerator = NA,
                           denominator = as.numeric(acceptableCoefs$mean.index),
                           index.value = NA,
                           pc.05 = acceptableCoefs$pc.05)
  return(emptyResDF)
}


## Compute sensibility,
## factors that vary : sub.name and coef.facc, coef.fls
empty.res.df <- MakeEmptyResDF(acceptable.coefs, "P01", "sau_dont_bh_2010_2019", "long-term", "bh", 2010)
res.df = empty.res.df[0,]


for(n.source in seq(1, length(sources$name))){
  tic("Compute index for sources")
  source.name <- sources$name[n.source]
  source.type <- sources$type[n.source]
  source.nut <- sources$nut[n.source]
  source.year <- sources$year[n.source]
  
  for (n.subcat in seq(1, length(list.select))){ 
    tic(paste(n.subcat, "subcat"))
    subcat.name <- list.select[n.subcat]
    print(paste("Watershed n", as.character(n.subcat), ":", subcat.name, sep = " "))
    
    point.res.df <- MakeEmptyResDF(acceptable.coefs, subcat.name, source.name, source.type, source.nut,  source.year)
    
    for(i in seq(1, nrow(point.res.df))){
      a <- point.res.df$coef.facc[i]
      b <- point.res.df$coef.fls[i]
      
      temp.ras <- raster.indexes[[paste0(subcat.name,"_", a, "_", b)]] *
        subcat.sources[[paste0(subcat.name, "_", source.name)]]
      
      point.res.df$numerator[i] <- cellStats(temp.ras, "sum", na.rm = TRUE) / sum(!is.na(values(temp.ras)))
      point.res.df$index.value[i] <- point.res.df$numerator[i] / point.res.df$denominator[i]
    } #end loop on acceptable coefs
    toc()
    res.df <- rbind(res.df, point.res.df)
  }
  toc() #
}

#Plot Sensibiliy index relative to YV1 sensibility index
ggplot(data = res.df%>%filter(source.nut == "all"), aes(x = coef.fls, y = coef.facc, fill = index.value))+
  geom_tile() + 
  scale_fill_fermenter(palette = "PiYG", breaks = c(0.25,0.5,0.75,0.95,1.05,1.25,1.5,2), direction = -1) + 
  scale_x_continuous(limits=c(0, 4), expand = c(0, -0.05)) +
  scale_y_continuous(limits=c(0, 2.4), expand = c(0, -0.05)) +
  coord_equal() + 
  facet_wrap(~subcat.name) +
  theme_bw() + 
  theme(legend.key.size =  unit(1, "cm"))

ggsave(paste0(results.folder,"/Relative_Sensitivity_YV1.png"))


PlotMetricSource <- function(allRes, sourceName){
  dataPlot <- allRes %>% filter(source.name == sourceName)
  
  metricSourcePlot <- 
    ggplot(data = dataPlot, aes(x = coef.fls, y = coef.facc, fill = index.value))+
    geom_tile() + 
    scale_fill_fermenter(palette = "PiYG",
                         breaks = c(0.25,0.5,0.75,0.95,1.05,1.25,1.5,2),
                         limits = c(0,1e9)) + 
    scale_x_continuous(limits=c(-0.05, 4), expand = c(0, -0.05)) +
    scale_y_continuous(limits=c(-0.05, 2), expand = c(0, -0.05)) +
    coord_equal() + 
    facet_wrap(~subcat.name) +
    theme_bw() + 
    theme(legend.key.size =  unit(1, "cm"))
  
  ggsave(paste0(results.folder,"/", sourceName, "_indexValue", ".png"))
  
}

for(source.name in sources.name){
  PlotMetricSource(res.df,  source.name)
}

remove(raster.indexes)
gc()
print("Saving workspace")
save.image(file = paste0(results.folder,"/2_Compute.RData"))



