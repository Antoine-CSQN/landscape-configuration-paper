library(raster) # For raster import
library(tidyverse) # For data wrangling and plotting
library(reshape2) # Melting data
library(tictoc)
library(cowplot)
library(lubridate)



work.dir <- "F:/doctorat/R/IDW3"
setwd(work.dir) 
results.folder <- paste0(work.dir,"/results") 
# folder with catchments delineation rasters 
wat.folder <- paste0(work.dir, "/raster/wat") 

source('1_PreProcess.R')
source('2_Compute.R')
source('3_Optimize.R')
