library(raster) # Raster operations
library(fasterize) # For fast rasterization
library(sf) # Shapefile operations 
library(tidyverse) # Data wrangling and plotting
library(reshape2) # Reshaping data frames
library(tictoc) # Timing the execution
library(lubridate) # Date operations

# Define working directory
work.dir <- "F:/doctorat/R/IDW_clean"
setwd(work.dir) 
# Subdirectories paths
shp.dir <- paste0(work.dir,"/shapefile")
results.dir <- paste0(work.dir,"/result")
cat.dir <- paste0(work.dir, "/raster/cat")
source.dir <- paste0(work.dir, "/raster/source")
fig.dir <- paste0(work.dir, "/figure")

# Names of flow accumulation and flow length to stream rasters
facc.name <- "Facc" 
fls.name <- "FLS" 
# Name of the mesoscale catchment shapefile (study area, Yvel's catchment in the paper)
study.area.filename <- "mesoscale_cat"
# Name of the headwater catchments shapefile 
hw.filename <- "hw_cat"
# Raster mask for streams and ditches
# raster with the same spatial resolution and extent as facc.rst and fls.rst
# NA is streams and ditches, 1 is everything else
ditch.stream.mask <- raster(paste0(source.dir, "/all1_ditchstreamNA.tif"))

# Column containing the name of the hw in the headwater catchments shapefile
hw.id <- "Name"
# Select headwater catchments : order<=3, no WWWTP, at least 20% agricultural 
hw.select <- c("P01", "P02", "P03", "P04", "P05",
                 "P07", "P12", "P13", "P16", "P17",
                 "P18", "P19", "P20", "P22", "P23",
                 "P24", "P25", "P26", "P29")

# Coefficient applied to convert Facc raster values to m2 and FLS to m
facc.to.m2 <- 100 
fls.to.m <- 1

## Range and increment of a(coef.facc) and b(coef.fls) values we want to explore
a.facc.values = seq(0.00, 2.5, 0.1) # Range from 0 to 2.5 with a 0.1 increment for parameter a
b.fls.values = seq(0.00, 4, 0.1) # Range from 0 to 4 with a 0.1 increment for parameter b

# In the paper, there is only one tested source (LU) which is agricultural areas (including riparian buffer strips), 
# taking the value of 1 (presence) or 0 (absence). 
# However the script is built to allow batch processing of different source layers
# It is also possible to use a weighted source layer with continuous values
# Here we have two sources : 
#     all (everything is a source, used to define the domain of exploration of parameters a and b)
#     agri1_other0 (described above)
sources.name <- c("all1_ditchstream0", "agri1_other0") # Name of the source raster (LU term of the LCI)
sources.type <- rep("long-term", 2) # For comparison of short-term (current) and long-term (legacy) sources
sources.year <- rep(2018, 2) # Year of the source
sources.nut <- c( "all", "agri") # what is considered source of nutrients ?

# Chosen parameters
water.params = c("TP","NO3")

# Read water chemistry data
chem.data <- read.csv2("data/StreamWaterChemistry.csv", stringsAsFactors = FALSE)
chem.data <- chem.data %>% 
  filter(p.sampling %in% hw.select) %>%
  rename(NO3 = nitrate) %>%
  mutate(date.sampling = as.Date(date.sampling, format="%d/%m/%Y"))
date.list <- unique(chem.data$date.sampling)

# compute dates with less than 2 HW not flowing
chem.datecount <- chem.data %>% 
  group_by(date.sampling) %>%
  mutate(n = sum(!is.na(TP))) %>%
  dplyr :: select(date.sampling,n) %>%
  distinct() %>%
  ungroup()

dates.flowing <- filter(chem.datecount, n>=(max(chem.datecount$n))-2) %>%
  dplyr :: select(date.sampling) 

# chem.median is a df containing median value for each water quality parameter 
# When at most 2 hw are not flowing
chem.median <- chem.data %>%
  filter(date.sampling %in% dates.flowing$date.sampling) %>%
  pivot_longer(cols = colnames(chem.data[,7:25]), names_to = "param") %>% #Columns 7:25 contain water quality data (one column per parameter)
  group_by(p.sampling, param) %>%
  summarise(median = median(value, na.rm = TRUE)) %>%
  ungroup() %>% 
  pivot_wider(names_from = param, values_from = median) 


source('1_PreProcess.R')
source('2_Compute.R')
source('3_Optimize.R')

