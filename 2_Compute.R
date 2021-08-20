# In this script, we compute the value of the LCI for
# all acceptable (a,b) (see limit of acceptability top5% < 95% total LCI)
# for all LU (sources) layers in the sampled subcatchments

## Extract Facc and FLS for all subcatchments
print("Preprocessing of subcatchment sources")
tic(name="Preprocessing of subcatchment sources: done")
subcat.sources <- list()
for(source.name in sources$name){
  source.path <- paste0(work.dir,"/raster/source/",source.name,".tif")
  source.raster <- raster(source.path)
  source.raster <- source.raster * ditch.stream.mask # burn stream/ rods with NAs
  print(paste0("===== ", source.path, " processsing ======="))
  source.type <- sources$type[which(sources$name == source.name)]
  print(paste0("===== source type: ",source.type, " ======="))
  tic(name = source.name)
  for(subcat.name in hw.select){
    trimmed.source.rst <- TrimRaster(source.raster, ReadCatRaster(subcat.name))
    name.tr <- paste0(subcat.name,"_",source.name)
    subcat.sources[[name.tr]] <- trimmed.source.rst * ditch.stream.mask # burn stream/ rods with NAs (0 kept for no source)
  }# end loop on subcatchments name.
  toc()
  
}# end loop on sources
toc() # approx 30s par source
# Random check that ditches/streams = NA, out of subcat = NA and sources is 0 or 1 
# plot(subcat.sources$P05_agri1_other0, colNA = "blue")

## Extract Facc and FLS for all subcatchments
fls.subcat <- list()
facc.subcat <- list()
print("Extracting Facc and FLS for all subcatchments")
tic("Extracting Facc and FLS for all subcatchments: done")
for(subcat.name in hw.select){
  print(paste0("===== Trimming Facc and FLS for ", subcat.name, " ======="))
  cat.fls <- TrimRaster(fls.rst, ReadCatRaster(subcat.name))
  cat.facc <- TrimRaster(facc.rst, ReadCatRaster(subcat.name))
  
  # Force ditches and streams = NA
  cat.facc[cat.fls == 0] <- NA 
  cat.fls[cat.fls == 0] <- NA 
  fls.subcat[[subcat.name]] <- cat.fls
  facc.subcat[[subcat.name]] <- cat.facc
}
toc() #60 sec
# Check patterns of FLS and FAcc on a random subcat
# plot(fls.subcat$P24, colNA = "blue")
# plot(facc.subcat$P24, colNA = "blue")



# Compute subcatchments facc^a * fls^-b rasters for acceptable (a,b)
# Pre processing step before final LCI computation
# Very memory intensive but speeds a lot the computation when 
# multiple nutrient sources layers (LU) are considered
# If a memory problem happens: compute facc^a * fls^-b rasters directly in the LCI computing loop as temp variables
raster.indexes <- list() 
print("Computing subcatchments facc^a * fls^-b rasters for acceptable (a,b)")
tic("Computing subcatchments facc^a * fls^-b rasters for acceptable (a,b): done")
for(subcat.name in hw.select){
  print(paste("Point :", subcat.name))
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
  } # End loop on (a,b)
} # End loop on subcat
toc() # 101 sec

#Function to generate an empty df in which we stock LCI results and parameters
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


# Initialize the results df 
empty.lci.df <- MakeEmptyResDF(acceptable.coefs,
                               hw.select[1],
                               sources.name[1],
                               sources.type[1],
                               sources.nut[1],
                               sources.year[1])
lci.df = empty.lci.df[0,]

# Loop on LU layers (source.name) and acceptable (a,b)
# LCI value is found in the index.value column
print("Computing subcatchments LCI values for acceptable (a,b)")
tic("Computing subcatchments LCI values for acceptable (a,b)")
for(n.source in seq(1, length(sources$name))){
  
  source.name <- sources$name[n.source]
  source.type <- sources$type[n.source]
  source.nut <- sources$nut[n.source]
  source.year <- sources$year[n.source]
  
  print(paste())
  for (n.subcat in seq(1, length(hw.select))){ 
    tic(paste(n.subcat, "subcat"))
    subcat.name <- hw.select[n.subcat]
    print(paste("Watershed n", as.character(n.subcat), ":", subcat.name, sep = " "))
    # point.lci.df is a results df (LCI values for 1 subcat and 1 LU layer)
    point.lci.df <- MakeEmptyResDF(acceptable.coefs, subcat.name, source.name, source.type, source.nut,  source.year)
    
    for(i in seq(1, nrow(point.lci.df))){
      a <- point.lci.df$coef.facc[i]
      b <- point.lci.df$coef.fls[i]
      temp.ras <- raster.indexes[[paste0(subcat.name,"_", a, "_", b)]] *
        subcat.sources[[paste0(subcat.name, "_", source.name)]]
      point.lci.df$numerator[i] <- cellStats(temp.ras, "sum", na.rm = TRUE) / sum(!is.na(values(temp.ras)))
      point.lci.df$index.value[i] <- point.lci.df$numerator[i] / point.lci.df$denominator[i]
    } #end loop on acceptable (a,b)
    
    lci.df <- rbind(lci.df, point.lci.df)
  } #end loop on subcat
  
} #end loop on sources layers (LU)
toc()

remove(raster.indexes)
gc()
print("Saving workspace")
save.image(file = paste0(results.dir,"/2_Compute.RData"))


# 
# PlotMetricSource <- function(allRes, sourceName){
#   dataPlot <- allRes %>% filter(source.name == sourceName)
#   
#   metricSourcePlot <- 
#     ggplot(data = dataPlot, aes(x = coef.fls, y = coef.facc, fill = index.value))+
#     geom_tile() + 
#     scale_fill_fermenter(palette = "PiYG",
#                          breaks = c(0.25,0.5,0.75,0.95,1.05,1.25,1.5,2),
#                          limits = c(0,1e9)) + 
#     scale_x_continuous(limits=c(-0.05, 4), expand = c(0, -0.05)) +
#     scale_y_continuous(limits=c(-0.05, 2), expand = c(0, -0.05)) +
#     coord_equal() + 
#     facet_wrap(~subcat.name) +
#     theme_bw() + 
#     theme(legend.key.size =  unit(1, "cm"))
#   
#   ggsave(paste0(results.dir,"/", sourceName, "_indexValue", ".png"))
#   
# }
# 
# for(source.name in sources.name){
#   PlotMetricSource(lci.df,  source.name)
# }
