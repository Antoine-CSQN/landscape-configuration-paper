# devtools::install_github("giswqs/whiteboxR")
library(whitebox) # Terrain processing
library(raster) # Raster manipulation
library(tidyverse) # Data wrangling
library(reshape2) # Data wrangling
library(sf) # Vector spatial data manipulation
library(patchwork) # Multi panels Figures
library(tictoc) # Time of execution
library(fasterize) # fast rasterization
library(stars) # Fast vectorization
library(ggspatial) # Plot RGB raster
library(ggrepel) # Repelling labels
library(scales) # Format numbers 
library(broom) # Organize lm results

raw.dem <-"F:/doctorat/R/IDW3/raster/DEM_raw.tif"
wbt_feature_preserving_smoothing(raw.dem, "raster/whitebox/smoothed.tif", filter=9, verbose_mode = TRUE)
wbt_breach_depressions("raster/whitebox/smoothed.tif", "raster/whitebox/breached.tif")

#Create isobasins given size in cells
wbt_isobasins(dem = "raster/whitebox/breached.tif", output = "raster/whitebox/isobas_1km.tif", size = 1e4, wd = NULL, verbose_mode = TRUE)
isobas.1km <- raster("raster/whitebox/isobas_1km.tif") 

raw.dem.rst <- raster(raw.dem)
dem.rst.50m <- aggregate(raw.dem.rst , fact = 5, fun=mean )
writeRaster(dem.rst.50m, "raster/whitebox/DEM_50m.tif")
dem.50m <-"F:/doctorat/R/IDW3/raster/whitebox/DEM_50m.tif"
wbt_feature_preserving_smoothing(dem.50m, "raster/whitebox/smoothed50m.tif", filter=9, verbose_mode = TRUE)
wbt_breach_depressions("raster/whitebox/smoothed50m.tif", "raster/whitebox/breached50m.tif")
#EAch pixel is 2500 m2 -- -> 1e4 = 25km2
wbt_isobasins(dem = "raster/whitebox/breached50m.tif", output = "raster/whitebox/isobas_25km.tif", size = 1e4, wd = NULL, verbose_mode = TRUE)
isobas.25km <- raster("raster/whitebox/isobas_25km.tif")
isobas.25km <- raster::disaggregate(isobas.25km, fact = 5, fun=majority)

#intersect isobas with YV1 (mesoscale catchement) - 1km2
isobas.1km <- isobas.1km * ReadCatRaster("YV1")
freq.1km <- as.data.frame(freq(isobas.1km))
freq.1km$reclass <- ifelse(freq.1km$count<1000,NA,freq.1km$value) #delete subcatment < 0.1 km
isobas.1km <- reclassify(isobas.1km, rcl = cbind(as.numeric(freq.1km$value), as.numeric(freq.1km$reclass)))
freq.1km <- as.data.frame(freq(isobas.1km)) %>% filter(!value == 0, !is.na(value))
plot(isobas.1km)

#intersect isobas with YV1 (mesoscale catchement) - 25km2
isobas.25km <- isobas.25km* ReadCatRaster("YV1")
freq.25km <- as.data.frame(freq(isobas.25km))
freq.25km$reclass <- ifelse(freq.25km$count<25000,NA,freq.25km$value) #delete subcatment < 2.5 km
isobas.25km <- reclassify(isobas.25km, rcl = cbind(as.numeric(freq.25km$value), as.numeric(freq.25km$reclass)))
freq.25km <- as.data.frame(freq(isobas.25km))  %>% filter(!value == 0, !is.na(value))
plot(isobas.25km)


# Return raster os the isobas n
SelectSubcat <- function(isoBas, n){
  maskRaster <- isoBas == c(n) #create raster mask, NA are not in the chatchment
  maskedRaster <- mask(maskRaster, maskRaster , maskvalue = 0) #use the mask
  trimmedRaster <- trim(maskedRaster) #crop the raster (less memory)
  return(trimmedRaster)
}

# plot(SelectSubcat(isobas.1km, 177)) #example


#Compute LCI-TP for each subcat
source.nut.rst <- raster("raster/source/sau_dont_bh_2010_2019.tif")
a <- 1.4
b <- 2.2
YV1.rst <- ReadCatRaster("YV1")
mean.index.YV1 <- cellStats(YV1.rst * source.nut.rst * facc.rst^a / fls.rst^b, stat ="mean", na.rm = TRUE)

#For 1km2 subassins
index.1km <- data.frame(id = as.numeric(freq.1km$value), value = NA)
for(n.subcat in index.1km$id){
  print(n.subcat)
  index.1km$value[which(index.1km$id == n.subcat)] <- cellStats(SelectSubcat(isobas.1km, n.subcat) *
                                                                  source.nut.rst * 
                                                                  facc.rst^a / fls.rst^b, stat ="mean", na.rm = TRUE) /
    mean.index.YV1
}

# For 25km2 subbasins
index.25km <- data.frame(id = as.numeric(freq.25km$value), value = NA)
for(n.subcat in index.25km$id){
  print(n.subcat)
  index.25km$value[which(index.25km$id == n.subcat)] <- cellStats(SelectSubcat(isobas.25km, n.subcat) *
                                                                    source.nut.rst * 
                                                                    facc.rst^a / fls.rst^b, stat ="mean", na.rm = TRUE) /
    mean.index.YV1
}

# vectorise with stars
isobas.1km.shp <- st_as_stars(isobas.1km)
isobas.1km.shp <- st_as_sf(isobas.1km.shp , as_points = FALSE, merge = TRUE)
isobas.25km.shp <- st_as_stars(isobas.25km)
isobas.25km.shp <- st_as_sf(isobas.25km.shp , as_points = FALSE, merge = TRUE)

## 1km2 aggregation map- Figure 5d
isobas.1km.sf <- st_as_sf(isobas.1km.shp)
isobas.1km.sf$area <- as.numeric(st_area(isobas.1km.sf))
#Join LCI-TP 1km2 values
isobas.1km.sf <- isobas.1km.sf %>%
  filter(layer %in% freq.1km$value) %>%
  filter(area > 50000) %>%
  left_join(index.1km, by=c("layer" = "id")) 

## Fig 5d
index.1km.map <- ggplot(isobas.1km.sf) + 
  geom_sf(aes(fill = value), lwd = 0.1) + 
  scale_fill_fermenter(palette = "PRGn", name = "LCI-TP",
                       breaks = c(0.25,0.5,0.75,0.95,1.05,1.25,2,4),
                       direction = -1) + 
  coord_sf(datum=st_crs(2154)) + 
  scale_x_continuous(breaks = c(300000,310000)) + 
  scale_y_continuous(breaks = c(6780000, 6790000, 6800000)) + 
  theme_bw() + 
  theme(plot.title = element_text(size = 10, face = "plain", colour = "black"),
        legend.title = element_text(size = 10, face = "plain", colour = "black"),
        legend.text = element_text(size = 8, face = "plain", colour = "black"),
        axis.title.x = element_text(size = 10, face = "plain", colour = "black"),    
        axis.title.y = element_text(size = 10, face = "plain", colour = "black"),    
        axis.text.x = element_text(size = 10, face = "plain", colour = "black"), 
        axis.text.y = element_text(size = 10, face = "plain", colour = "black"),
        strip.text.x = element_text(size = 10, face = "plain", colour = "black" ),
        strip.text.y = element_text(size = 10, face = "plain", colour = "black")) + 
  theme(axis.text.y = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
  theme(legend.position = c(0.8, 0.8),
        legend.background = element_rect(fill="transparent"),
        legend.key.height = unit(5, "mm"))
index.1km.map 

## 1km2 aggregation histogram - Figure 5f
hist.1km.map <- ggplot(isobas.1km.sf, aes(x = value)) +
  geom_histogram(binwidth = 0.1, colour = "black", lwd = 0.1, aes(fill = ..x..)) + 
  scale_x_continuous(breaks = seq(0,6,1), expand = c(0,0), limits = c(0,6)) +
  scale_y_continuous(expand = c(0,0)) +
  scale_fill_fermenter(palette = "PRGn",
                       breaks = c(0.25,0.5,0.75,0.95,1.05,1.25,2,4),
                       direction = -1) + 
  geom_vline(xintercept = 1, color = "red", linetype = "dashed", lwd = 0.2) + 
  geom_text(aes(x = 3.5, y = 30, label = paste0( "min - max = ", paste0(round( min(isobas.1km.sf$value),2 ),
                                                                 " - ", round( max(isobas.1km.sf$value),2 )) )), size = 3.5) + 
  geom_text(aes(x = 3.5, y = 25, label = paste0( "Median = ", round( median(isobas.1km.sf$value),2 ) )), size = 3.5) + 
  geom_text(aes(x = 3.5, y = 20, label = paste0( "Mean = ", round( mean(isobas.1km.sf$value),2 ) )), size = 3.5) + 
  geom_text(aes(x = 3.5, y = 15, label = paste0( "SD = ", round( sd(isobas.1km.sf$value),2 ) )), size = 3.5) + 
  labs(x=NULL, y=NULL) +
  theme_bw() + 
  theme_bw() + 
  theme(axis.text.x = element_text(size = 10, face = "plain", colour = "black"), 
        axis.text.y = element_text(size = 10, face = "plain", colour = "black")) + 
  theme(legend.position = "none")
hist.1km.map

# 25km2 aggregation
# Keep YV1 values and join them to index values
isobas.25km.sf <- st_as_sf(isobas.25km.shp)
isobas.25km.sf$area <- as.numeric(st_area(isobas.25km.sf))
isobas.25km.sf <- isobas.25km.sf %>%
  filter(layer %in% freq.25km$value) %>%
  filter(area > 50000) %>%
  left_join(index.25km, by=c("layer" = "id")) 

## 25 km2 aggregation map - Figure 5e
index.25km.map <- ggplot(isobas.25km.sf) + 
  geom_sf(aes(fill = value), lwd = 0.1) + 
  scale_fill_fermenter(palette = "PRGn", name = "LCI-TP",
                       breaks = c(0.25,0.5,0.75,0.95,1.05,1.25,2,4),
                       direction = -1, 
                       limits = c(0,8)) + 
  geom_rect(aes(xmin = 303500, xmax = 310000, ymin = 6790500, ymax = 6792500),
            fill = NA, colour = "black", lwd = 1) + #Exerpt of risk map and LCI-TP[field]
  coord_sf(datum=st_crs(2154)) + 
  scale_x_continuous(breaks = c(300000,310000)) + 
  scale_y_continuous(breaks = c(6780000, 6790000, 6800000)) + 
  theme_bw() + 
  theme(plot.title = element_text(size = 10, face = "plain", colour = "black"),
        legend.title = element_text(size = 10, face = "plain", colour = "black"),
        legend.text = element_text(size = 8, face = "plain", colour = "black"),
        axis.title.x = element_text(size = 10, face = "plain", colour = "black"),    
        axis.title.y = element_text(size = 10, face = "plain", colour = "black"),    
        axis.text.x = element_text(size = 10, face = "plain", colour = "black"), 
        axis.text.y = element_text(size = 10, face = "plain", colour = "black"),
        strip.text.x = element_text(size = 10, face = "plain", colour = "black" ),
        strip.text.y = element_text(size = 10, face = "plain", colour = "black")) + 
  theme(axis.text.y = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
  theme(legend.position = c(0.8, 0.8),
        legend.background = element_rect(fill="transparent"),
        legend.key.height = unit(5, "mm"))
index.25km.map 

## 25 km2 aggregation hist - Figure 5g
hist.25km.map <- ggplot(isobas.25km.sf,aes(x = value)) +
  geom_histogram(binwidth = 0.1, colour = "black", lwd = 0.1, aes(fill = ..x..)) + 
  scale_x_continuous(breaks = seq(0, 6, 1), expand = c(0,0), limits = c(0,6)) +
  scale_y_continuous(expand = c(0,0), limits = c(0,4.5)) +
  scale_fill_fermenter(palette = "PRGn",
                       breaks = c(0.25,0.5,0.75,0.95,1.05,1.25,2,4),
                       direction = -1) + 
  geom_vline(xintercept = 1, color = "red", linetype = "dashed", lwd = 0.2) + 
  geom_text(aes(x = 3.5, y = 4, label = paste0( "min - max = ", paste0(round( min(isobas.25km.sf$value),2 ),
                                                                       " - ", round( max(isobas.25km.sf$value),2 )) )), size = 3.5) + 
  geom_text(aes(x = 3.5, y = 3.5, label = paste0( "Median = ", round( median(isobas.25km.sf$value),2 ) )), size = 3.5) + 
  geom_text(aes(x = 3.5, y = 3, label = paste0( "Mean = ", round( mean(isobas.25km.sf$value),2 ) )), size = 3.5) + 
  geom_text(aes(x = 3.5, y = 2.5, label = paste0( "SD = ", round( sd(isobas.25km.sf$value),2 ) )), size = 3.5) + 
  labs(x=NULL, y=NULL) +
  theme_bw() + 
  theme(axis.text.x = element_text(size = 10, face = "plain", colour = "black"), 
        axis.text.y = element_text(size = 10, face = "plain", colour = "black")) + 
  theme(legend.position = "none")
hist.25km.map

# Low part of the Fig 5 panel
fig5.low <- (index.1km.map | index.25km.map) | (hist.1km.map / hist.25km.map) + 
  plot_layout(tag_level = 'new')
fig5.low  
ggsave("figure/paper/Fig5_low.tiff", fig5.low, width = 174, units = "mm", dpi = 700)
ggsave("figure/paper/Fig5_low.pdf", fig5.low, width = 174, units = "mm", dpi = 600)


## High part (A, B & C) of Figure 5 panel
# Valeur moyenne par parcelle RPG 
rpg.shp <- readOGR("shapefile/RPG_Yvel_2018.shp")
source.nut.rst <- raster("raster/source/sau_dont_bh_2010_2019.tif")
index.rst <- source.nut.rst * facc.rst^a / fls.rst^b

# Extract raster values to list object
tic()
r.vals <- raster :: extract(index.rst, rpg.shp, fun = mean, na.rm = TRUE)
toc() # about 20 min

# Join mean values to polygon data
rpg.shp@data <- data.frame(rpg.shp@data, mean.index = r.vals)
rpg.sf <- st_as_sf(rpg.shp)

rpg.values <- data.frame(mean.index = rpg.sf$mean.index[!is.na(rpg.sf$mean.index)]) / mean.index.YV1# 22.69 is the mean value of ses
ggplot(rpg.values, aes(x = mean.index)) + 
  geom_histogram() + 
  scale_x_log10(breaks = c(1e-4, 1e-3, 1e-2, 1e-1, 1, 1e1, 1e2, 1e3), limits = c(1e-5, 1e4)) + 
  scale_y_continuous(expand = c(0, 0)) + 
  geom_vline(xintercept = 1, color = "red", linetype = "dashed") + 
  theme_minimal() 

rpg.sf.YV1 <- rpg.sf %>%
  filter(mean.index > 0) %>%
  mutate(mean.index = mean.index / mean.index.YV1)

# Find how much plots are responsible for 95% of the sum index.
rpg.sf.YV1$sum_index <- rpg.sf.YV1$mean.index * rpg.sf.YV1$SURF_PARC
sum_index_YY1 <- sum(rpg.sf.YV1$sum_index)
rpg.sf.YV1 <- rpg.sf.YV1 %>% arrange(-sum_index) %>%
  mutate(cum_sum = cumsum(sum_index/sum_index_YY1),
         n_plot = seq(1, nrow(rpg.sf.YV1)))

rpg.sf.YV1 %>% filter(cum_sum > 0.95) %>% slice(1) %>% pull(n_plot) %>% as.numeric()

# How much weights for the top 20 % ? 1965 = 9826/5
rpg.sf.YV1[1965,]

# How much are less than one ?
nb_plot_inf1 <- rpg.sf.YV1 %>% 
  arrange(mean.index) %>% 
  mutate(n_plot = seq(1, nrow(rpg.sf.YV1))) %>%
  filter(mean.index < 1) %>% 
  arrange(-mean.index) %>% 
  slice(1) %>% pull(n_plot)
nb_plot_inf1/nrow(rpg.sf.YV1)


hydro <- st_read("shapefile/hydro_yvel_dissolve.shp")
roads <- st_read("shapefile/road_yvel_dissolve.shp")
ortho <- brick("F:/doctorat/bdd/geo/LacAuDuc/Images/OrthoYvel2012_2014_2M5.tif")
ortho_trimmed <- crop(ortho, extent(303500, 310000, 6790500, 6792500))

#Make the "top risk raster" based on top 5% and top 1% of pixels according to LCI-TP
optim.index.rst <- facc.rst^1.4 * fls.rst^(-2.2) * source.nut.rst * YV1.rst
thresholds <- quantile(optim.index.rst, probs = c(0.95,0.99))
optim.index.top01 <- optim.index.rst > thresholds[2]
optim.index.top05 <- optim.index.rst > thresholds[1]
optim.all <-  sum(optim.index.top01, optim.index.top05, na.rm = TRUE)
#Make the map for all YV1
optim.YV1 <- optim.all * source.nut.rst + source.nut.rst
optim.YV1 <- optim.YV1 * ReadCatRaster("YV1")
optim.all.sau.ll <- optim.YV1
optim.all.sau.ll[optim.all.sau.ll == 0 |optim.all.sau.ll == 1] <- NA
risk_trimmed <- crop(optim.all.sau.ll, extent(303500, 310000, 6790500, 6792500))

#Figure 5a : ortophoto + top risk pixels
risk.zoom <- ggplot()+ 
  layer_spatial(ortho_trimmed) +
  layer_spatial(risk_trimmed) +
  scale_fill_binned(na.value = NA, low="yellow", high="red", breaks = 2.5) + 
  geom_sf(data = roads, colour = "#000000", lwd = 0.5) +
  geom_sf(data = hydro, colour = "#5588EE", lwd = 0.5) + 
  coord_sf(xlim = c(303500, 310000), ylim = c(6790500, 6792500), datum=st_crs(2154)) + 
  scale_x_continuous(breaks = seq(304000, 309000, 1000), expand = c(0,0)) + 
  scale_y_continuous(breaks = c(6791000, 6792000), expand = c(0,0)) + 
  theme_bw() + 
  theme(legend.title = element_text(size = 10, face = "plain", colour = "black"),
        legend.text = element_text(size = 10, face = "plain", colour = "black"),
        axis.text.x = element_text(size = 10, face = "plain", colour = "black"), 
        axis.text.y = element_text(size = 10, face = "plain", colour = "black")) + 
  theme(axis.text.y = element_text(angle = 90, vjust = 0.5, hjust = 0.5)) + 
  theme(legend.position = "none")
risk.zoom 

mean.index.by.plot.zoom <- ggplot() + 
  geom_sf(data = rpg.sf.YV1, aes(fill = mean.index), colour = NA) + 
  scale_fill_fermenter(breaks = c(0.025, 0.1, 0.5, 0.75, 0.95, 1.05, 1.5, 5, 10, 40), limits = c(5e-4, 50),
                       trans = "log", palette = "PRGn", direction = -1, name = "LCI-TPfield") +
  geom_sf(data = roads, colour = "#000000", lwd = 0.5) + 
  geom_sf(data = hydro, colour = "#5588EE", lwd = 0.5) + 
  coord_sf(xlim = c(303500, 310000), ylim = c(6790500, 6792500), datum=st_crs(2154)) + 
  scale_x_continuous(breaks = seq(304000, 309000, 1000), expand = c(0,0)) + 
  scale_y_continuous(breaks = c(6791000, 6792000), expand = c(0,0)) + 
  theme_bw() + 
  theme(legend.title = element_text(size = 10, face = "plain", colour = "black"),
        legend.text = element_text(size = 8, face = "plain", colour = "black"),
        axis.text.x = element_text(size = 10, face = "plain", colour = "black"), 
        axis.text.y = element_text(size = 10, face = "plain", colour = "black")) + 
  theme(axis.text.y = element_text(angle = 90, vjust = 0.5, hjust = 0.5)) + 
  theme(legend.position = c(0.90,0.5),
        legend.key.height = unit(7, "mm")) 
mean.index.by.plot.zoom

hist.plot.map <- ggplot(rpg.sf.YV1, aes(x = mean.index)) +
  geom_histogram(bins = 21, colour = "black", lwd = 0.1, aes(fill = 10^(..x..) )) + 
  geom_vline(xintercept = mean(rpg.sf.YV1$mean.index), color = "red", linetype = "dashed", lwd = 0.5) + 
  geom_vline(xintercept = median(rpg.sf.YV1$mean.index), color = "blue", linetype = "dashed", lwd = 0.5) + 
  scale_x_log10(breaks = c(1e-3,1e-2,1e-1,1e0,1e1,1e2), expand = c(0,0), limits = c(1e-4,5e2)) +
  scale_y_continuous(expand = c(0,0)) +
  scale_fill_fermenter(palette = "PRGn",
                       breaks = c(0.025, 0.1, 0.5, 0.75, 0.95, 1.05, 1.5, 5, 10, 40),
                       direction = -1) + 
  geom_text(aes(x = 1e-3, y = 1250, label = paste0( "q5 - q95 = ", paste0(formatC( quantile(mean.index, probs = 0.05), 2, format = "e" ), " - ",
                                                                          formatC( quantile(mean.index, probs = 0.95), 2, format = "e" ) ) )),
            size = 3.0) +
  geom_text(aes(x = 1e-3, y = 1050, label = paste0( "q25 - q75 = ", paste0(formatC( quantile(mean.index, probs = 0.25), 2, format = "e" ), " - ",
                                                                           formatC( quantile(mean.index, probs = 0.75), 2, format = "e" ) ) )),
            size = 3.0) +
  geom_text(aes(x = 1e-3, y = 850, label = paste0( "Median = ", round( median(mean.index),2 ) )),
            size = 3.0, color = "blue") + 
  geom_text(aes(x = 1e-3, y = 650, label = paste0( "Mean = ", round( mean(mean.index),2 ) )),
            size = 3.0, color = "red") + 
  geom_text(aes(x = 1e-3, y = 450, label = paste0( "SD = ", round( sd(mean.index),2 ) )),
            size = 3.0) + 
  labs(x = NULL , y = NULL) +
  theme_bw() + 
  theme(legend.title = element_text(size = 10, face = "plain", colour = "black"),
        legend.text = element_text(size = 10, face = "plain", colour = "black"),
        axis.text.x = element_text(size = 10, face = "plain", colour = "black"), 
        axis.text.y = element_text(size = 10, face = "plain", colour = "black")) + 
  theme(axis.text.y = element_text(angle = 90, vjust = 0.5, hjust = 0.5)) + 
  theme(legend.position =  "none")
hist.plot.map

# fig.5 <- risk.zoom + mean.index.by.plot.zoom + hist.plot.map + fig5.low +
#   plot_layout(ncol = 1, heights = c(1,1,0.5,1.5), widths = c(1,1,1,1)) +
#   plot_annotation(tag_levels = "A")

fig.5.high <- risk.zoom + mean.index.by.plot.zoom + hist.plot.map +
  plot_layout(ncol = 1, heights = c(1,1,0.5)) 

ggsave("figure/paper/Fig5_high.tiff", fig.5.high, width = 174, units = "mm", dpi = 700)
ggsave("figure/paper/Fig5_high.pdf", fig.5.high, width = 174, units = "mm", dpi = 600)


## Figure 6 - Ratio LCI(0,0)/LCI-TP
random.cats <- st_read("shapefile/random_subcat_500.shp")
random.cats <- random.cats %>% rename(id = NumBV) %>%
  mutate(area.km2 = SurfaceHa/ 100) %>%
  mutate(compo = NA,
         distrib = NA)

source.nut.rst <- raster("raster/source/sau_dont_bh_2010_2019.tif")
index.rst <- source.nut.rst * facc.rst^a / fls.rst^b

mean.compo.YV1 <- as.numeric(preprocess.results %>%
                               filter(source.nut == "sau_dont_bh",
                                      coef.facc == 0, coef.fls == 0) %>% 
                               select(mean.index) )

tic()
for(ncat in seq(1, nrow(random.cats))){
  print(ncat)
  ncat.rst <- fasterize(random.cats[ncat,], index.rst, field="id", fun="last")
  ncat.rst <- trim(ncat.rst)
  ncat.rst[!is.na(ncat.rst)] <- 1
  source.temp <- source.nut.rst * ncat.rst
  index.temp <- index.rst * ncat.rst
  random.cats$compo[ncat] <- cellStats(source.temp, "mean", na.rm = TRUE) 
  random.cats$distrib[ncat] <- cellStats(index.temp, "mean", na.rm = TRUE) / mean.index.YV1
}
toc() # 315 sec
# # Write results
# st_write(random.cats, "shapefile/random_subcat_500_LCI.shp")


random.cats <- random.cats %>% filter(distrib > 0)

random.cats$compo2 <-random.cats$compo / 
  as.numeric(random.cats$compo[random.cats$area.km2 == max(random.cats$area.km2)]) # normalized by the mean compo of the biggest subcatchment
random.cats$compo.norm <- (random.cats$compo - min(random.cats$compo)) / (max(random.cats$compo) - min(random.cats$compo))
random.cats$distrib.norm <- (random.cats$distrib - min(random.cats$distrib)) / (max(random.cats$distrib) - min(random.cats$distrib))
random.cats$compo.ratio <- random.cats$compo2 / random.cats$distrib
random.cats$compo.ratio2 <- random.cats$distrib / random.cats$compo2

# Plot
fig6  <- ggplot(random.cats, aes(area.km2, compo.ratio2)) + 
  scale_x_sqrt(breaks = c(1, 5, 10, 25, 50, 100, 250), minor_breaks = NULL, expand = c(0.01,0.01)) + 
  geom_hline(aes(yintercept = 1), lwd = 0.25, lty = 2) + 
  geom_point(aes(color = compo*100), size = 1.5, alpha = 0.5) +
  scale_colour_viridis_b(breaks = 100*c(0.25,0.5,0.66,0.75), limits = c(20,80), name = "% agricultural surface") + 
  labs(x = "Area (km²)",
       y = "LCI-TP / LCI(0, 0)",
       title = NULL) + 
  theme_bw() + 
  theme(plot.title = element_text(size = 10, face = "plain", colour = "black"),
        legend.title = element_text(size = 10, face = "plain", colour = "black"),
        legend.text = element_text(size = 10, face = "plain", colour = "black"),
        axis.title.x = element_text(size = 10, face = "plain", colour = "black"),    
        axis.title.y = element_text(size = 10, face = "plain", colour = "black"),    
        axis.text.x = element_text(size = 10, face = "plain", colour = "black"), 
        axis.text.y = element_text(size = 10, face = "plain", colour = "black"),
        strip.text.x = element_text(size = 10, face = "plain", colour = "black" ),
        strip.text.y = element_text(size = 10, face = "plain", colour = "black")) + 
  theme(legend.position = c(0.7,0.85), 
        legend.direction = "horizontal",
        legend.key.width = unit (7,"mm"))
fig6
ggsave("figure/paper/fig6.tiff", plot = fig6, width = 174, height = 80, units = "mm", dpi = 700)
ggsave("figure/paper/fig6.pdf", plot = fig6, width = 174, height = 80, units = "mm", dpi = 600)

## Figure 1S
random.cats$area.cat <- cut(random.cats$area.km2, breaks = c(0,50, 400))
random.cats <- random.cats %>% mutate(area.cat.label = if_else(area.cat == levels(area.cat)[1],
                                                               "Area <= 50 km²", "Area > 50 km²"))

#garder les BV dont la compo est similaire à celle des BV de plus de 50 km2
min.compo.50 <- min(random.cats %>% filter(area.km2 > 50) %>% select(compo2) %>% st_drop_geometry()) 
max.compo.50 <- max(random.cats %>% filter(area.km2 > 50) %>% select(compo2) %>% st_drop_geometry()) 
filtered.cats <- random.cats %>% filter(compo2 >= min.compo.50, compo2 <= max.compo.50, distrib != 0)

sd.distrib.sup <- filtered.cats %>% filter(area.km2 > 50) %>% select(distrib) %>% st_drop_geometry()
sd.distrib.sup <- round(sd(unlist(sd.distrib.sup)), 3)
sd.distrib.inf <- filtered.cats %>% filter(area.km2 < 50) %>% select(distrib) %>% st_drop_geometry()
sd.distrib.inf <- round(sd(unlist(sd.distrib.inf)), 3)

sd.compo.sup <- filtered.cats %>% filter(area.km2 > 50) %>% select(compo2) %>% st_drop_geometry()
sd.compo.sup <- round(sd(unlist(sd.compo.sup)), 3)
sd.compo.inf <- filtered.cats %>% filter(area.km2 < 50) %>% select(compo2) %>% st_drop_geometry()
sd.compo.inf <- round(sd(unlist(sd.compo.inf)), 3)

mean.distrib.sup <- filtered.cats %>% filter(area.km2 > 50) %>% select(distrib) %>% st_drop_geometry()
mean.distrib.sup <- round(mean(unlist(mean.distrib.sup)), 3)
mean.distrib.inf <- filtered.cats %>% filter(area.km2 < 50) %>% select(distrib) %>% st_drop_geometry()
mean.distrib.inf <- round(mean(unlist(mean.distrib.inf)), 3)

mean.compo.sup <- filtered.cats %>% filter(area.km2 > 50) %>% select(compo2) %>% st_drop_geometry()
mean.compo.sup <- round(mean(unlist(mean.compo.sup)), 3)
mean.compo.inf <- filtered.cats %>% filter(area.km2 < 50) %>% select(compo2) %>% st_drop_geometry()
mean.compo.inf <- round(mean(unlist(mean.compo.inf)), 3)

## Fig 1S High

ggplot(random.cats, aes(x = compo2, y= distrib, color = area.cat.label, size = area.cat.label)) +
  geom_point(pch = 1) + 
  geom_abline(slope = 1.00, intercept = 0, linetype = "dashed", size = 0.75, alpha = 0.5) +
  geom_rect(xmin = min.compo.50, xmax = max.compo.50, ymin = 0.2, ymax = 2.5, colour="grey50", fill=NA, size = 0.25) + 
  scale_x_continuous(expand = c(0,0), limits = c(0,1.55)) + 
  scale_y_continuous(expand = c(0,0.05)) + 
  guides(size=guide_legend(title=NULL), color=guide_legend(title=NULL)) + 
  scale_color_manual(values = c("#5bc0de", "#337ab7")) + 
  labs(title = NULL,
       subtitle = NULL,
       x = "Landscape Composition: LCI(0, 0)",
       y = "Landscape Configuration: LCI-TP") + 
  theme_bw() + 
  theme(strip.background = element_rect(color="black", fill="white", size=1.5, linetype="solid")) + 
  theme(plot.title = element_text(size = 10, face = "plain", colour = "black"),
        legend.title = element_text(size = 10, face = "plain", colour = "black"),
        legend.text = element_text(size = 10, face = "plain", colour = "black"),
        axis.title.x = element_text(size = 10, face = "plain", colour = "black"),    
        axis.title.y = element_text(size = 10, face = "plain", colour = "black"),    
        axis.text.x = element_text(size = 10, face = "plain", colour = "black"), 
        axis.text.y = element_text(size = 10, face = "plain", colour = "black"),
        strip.text.x = element_text(size = 10, face = "plain", colour = "black" ),
        strip.text.y = element_text(size = 10, face = "plain", colour = "black")) + 
  theme(legend.position = c(0.2,0.8)) 
ggsave(filename = "figure/paper/Fig1S_high.tiff", height = 75, width = 174 , units = "mm", dpi = 600)
ggsave(filename = "figure/paper/Fig1S_high.pdf", height = 75, width = 174 , units = "mm", dpi = 600)


# non parametric tests
cats.inf <- filtered.cats%>% filter(area.cat.label == "Area <= 50 km²")
spear.cor.inf <- cor.test(cats.inf$compo, cats.inf$distrib, method = "spearman")
cats.sup <- filtered.cats%>% filter(area.cat.label == "Area > 50 km²")
spear.cor.sup <- cor.test(cats.sup$compo, cats.sup$distrib, method = "spearman")

lm.inf <- lm(distrib~compo, filtered.cats%>% filter(area.cat.label == "Area <= 50 km²"))
tidy.lm.inf <- tidy(lm.inf)
glance.lm.inf <- glance(lm.inf)
label.lm.inf <- paste0("R² = ", round(glance.lm.inf$r.squared, 3), ", p-val = ", sprintf('%.3f', glance.lm.inf$p.value))

lm.sup <- lm(distrib~compo, cats.sup50)
tidy.lm.sup <- tidy(lm.sup)
glance.lm.sup <- glance(lm.sup)
label.lm.sup <- paste0("R² = ", round(glance.lm.sup$r.squared, 3), ", p-val = ", sprintf('%.1g', glance.lm.sup$p.value))

label.compo.inf <- paste0("x: ", "mean", " = ", mean.compo.inf, ", ", "SD", " = ", sd.compo.inf)
label.distrib.inf <- paste0("y: ", "mean", " = ", mean.distrib.inf, ", ", "SD", " = ", sprintf('%.3f',sd.distrib.inf))

label.compo.sup <- paste0("x: ", "mean", " = ", mean.compo.sup, ", ", "SD", " = ", sd.compo.sup)
label.distrib.sup <- paste0("y: ", "mean", " = ", mean.distrib.sup, ", ", "SD", " = ", sd.distrib.sup)

dat_text <- data.frame(
  label = c(label.compo.inf, label.distrib.inf, label.lm.inf, label.compo.sup, label.distrib.sup, label.lm.sup),
  area.cat = levels(filtered.cats$area.cat)[c(1,1,1,2,2,2)],
  area.cat.label = levels(as.factor(filtered.cats$area.cat.label))[c(1,1,1,2,2,2)],
  x     = c(1.05, 1.05, 1.05, 1.05, 1.05, 1.05),
  y     = c(2.25, 2.10, 1.95, 2.25, 2.10, 1.95)
)


## Fig 1S Low
ggplot(filtered.cats,
       aes(x = compo2, y = distrib, color = area.cat )) + 
  geom_point(size = 3, pch = 1) +
  geom_smooth(method = "lm", linetype = "dashed", se = FALSE, size = 1.25, alpha = 0.5) + 
  geom_vline(xintercept = 1, size = 0.75, lty = "dashed", alpha = 0.5) + 
  geom_text(data = dat_text, aes(x = x, y = y, label=label), size = 3.5) +
  
  facet_wrap(~area.cat.label) +
  scale_color_manual(values = c("#5bc0de", "#337ab7")) + 
  labs(x = "Landscape Composition: LCI(0, 0)",
       y = "Landscape Configuration: LCI-TP") + 

  theme_bw() + 
  theme(strip.background = element_rect(color="black", fill="white", size=1.5, linetype="solid")) + 
  theme(plot.title = element_text(size = 10, face = "plain", colour = "black"),
        legend.title = element_text(size = 10, face = "plain", colour = "black"),
        legend.text = element_text(size = 10, face = "plain", colour = "black"),
        axis.title.x = element_text(size = 10, face = "plain", colour = "black"),    
        axis.title.y = element_text(size = 10, face = "plain", colour = "black"),    
        axis.text.x = element_text(size = 10, face = "plain", colour = "black"), 
        axis.text.y = element_text(size = 10, face = "plain", colour = "black"),
        strip.text.x = element_text(size = 10, face = "plain", colour = "black" ),
        strip.text.y = element_text(size = 10, face = "plain", colour = "black")) + 
  theme(legend.position = "none") 
ggsave(filename = "figure/paper/Fig1S_low.tiff", height = 75, width = 174 , units = "mm", dpi = 600)
ggsave(filename = "figure/paper/Fig1S_low.pdf", height = 75, width = 174 , units = "mm", dpi = 600)


