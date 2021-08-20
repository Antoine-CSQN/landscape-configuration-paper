### Script to find the optimized values for facc (a) and fls (b) coef
# Run 1_Preprocess and 2_Compute first

lci.df$source.nut <- as.factor(lci.df$source.nut)
lci.df$source.name <- as.factor(lci.df$source.name)


#General function for computing the correlations
ComputeCorrelation <- function(geoChemData, allRes, paraM, datE, sourceNut, sourceType, verbose = FALSE){
  chemData <- geoChemData %>% select(subcat.name = p.sampling, date.sampling, !!paraM) %>%
    filter(date.sampling == datE)
  
  if(sourceType == "instant"){
    yearChem <- year(datE)
    indexData <- allRes %>% filter(source.nut == sourceNut,
                                   source.type == sourceType,
                                   source.year == yearChem)
    sourceName <- paste0(sourceNut,"_",yearChem)
  }
  if(sourceType == "long-term"){
    indexData <- allRes %>% filter(source.nut == sourceNut,
                                   source.type == sourceType)
    sourceName <- paste0(sourceNut,"_2010_2019")
  }
  
  corRes <- acceptable.coefs
  corRes$param <- paraM
  corRes$date <- datE
  corRes$source.name <- sourceName
  corRes$source.nut <- sourceNut
  corRes$source.type <- sourceType
  corRes$R <- NA
  corRes$p.value <- NA
  
  for(i in seq(1, nrow(corRes))){
    coefFacc <- as.numeric(corRes$coef.facc[i])
    coefFls <- as.numeric(corRes$coef.fls[i])
    
    indexDataCoef <- indexData %>% filter(coef.facc == coefFacc,
                                          coef.fls == coefFls)
    if(verbose == TRUE){ # debug
      data.corr <- as.data.frame(x = as.numeric(chemData[[paraM]]), y =  as.numeric(indexDataCoef$index.value))
      str(data.corr)
      str(chemData[[paraM]])
      str(indexDataCoef$index.value)
      str(cor.test(x = as.numeric(chemData[[paraM]]),
                   y = as.numeric(indexDataCoef$index.value),
                   method = 'spearman'))
    }
    
    corr <- cor.test(x = as.numeric(chemData[[paraM]]),
                     y = as.numeric(indexDataCoef$index.value),
                     method = 'spearman')
    
    corRes$R[i] <- corr$estimate
    corRes$p.value[i] <- corr$p.value
    

  }
  return(corRes)
}


####################
### Compute the correlations for acceptable (a,b) coefs
### for chosen parameters and all dates of sampling
###################


# Init table of correlations
all.correl <- ComputeCorrelation(chem.data, lci.df, water.params[1], date.list[1], sources.nut[1], sources.type[1])
all.correl <- all.correl[0,]


# Compute all correlations for each date
print("Computing all correlations for acceptable (a,b)")
tic("Compute correlations for param X source.nut X source.type X date")
for(param in water.params){
  print(paste("Parameter:", param))
  tic(param)
  
  for(s in seq(1, nrow(sources))){
    source.nut <- sources$nut[s]
    source.type <- sources$type[s]
    source.name <- sources$name[s]
    print(paste("Source:", source.name))
    tic(source.name)
    for(date in date.list){
      
      corrCombin <- ComputeCorrelation(chem.data, lci.df, param, as.Date(date, origin = '1970-01-01'), source.nut, source.type)
      all.correl <- rbind(all.correl, corrCombin)
      
    }# end loop on dates
    toc()
  }# end loop on source (approx 30 sec for 1 source* 1 param*30 dates*595 acceptable coefs)
  toc() 
}# end loop on water quality params
toc() 

all.correl <- all.correl %>% distinct()
all.correl$metric <- "l_config_idx"

#Keep max R for all param*metric*source.nut*source.type
best.correl <- all.correl %>% group_by(param, date, source.nut, source.type) %>%
  slice(which.max(R)) %>% ungroup()

#Get l_compo_idx value for all param*metric*source.nut*source.type
l_compo_idxCor <- all.correl  %>%
  filter(coef.facc == 0, coef.fls == 0) %>%
  mutate(metric = "l_compo_idx") %>%
  distinct()
# add it
all.correl <- rbind.data.frame(all.correl, l_compo_idxCor)
best.correl <- rbind.data.frame(best.correl, l_compo_idxCor)

# relevel metric
best.correl$metric <- factor(best.correl$metric,levels = c("l_compo_idx", "l_config_idx"))


####################
### Compute the correlations for acceptable (a,b) coefs
### for chosen parameters and median value
###################


# Initialisation du tableau
median.correl <- ComputeCorrelation(chem.data, lci.df, water.params[1], date.list[1], sources.nut[1], sources.type[1])
median.correl <- median.correl[0,]

#Compute all correlations for each date
print("Computing all correlations for acceptable (a,b)")
tic("Compute correlations for param X source.nut X source.type X date")
for(param in water.params){
  print(paste("Parameter:", param))
  tic(param)
  
  for(s in seq(1, nrow(sources))){
    source.type <- sources$type[s]
    source.nut <- sources$nut[s]
    source.name <- sources$name[s]
    print(paste("Source:", source.name))
    tic(source.name)
    
    chem.median$date.sampling <- "median"
    
    corrCombin <- ComputeCorrelation(chem.median, lci.df, param, "median" , source.nut, source.type)
    median.correl <- rbind(median.correl, corrCombin)
    
    toc()
  }# end loop on source 
  toc() 
}# end loop on water quality params
toc() 
median.correl$metric <- "l_config_idx"

# Save Workspace
save.image(file = paste0(results.dir,"/3_Optimized_.RData"))
