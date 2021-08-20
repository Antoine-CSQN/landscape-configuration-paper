### Script to find the optimized values for facc and fls coef
# Run 1_Preprocess and 2_Compute first

# Chosen parameters
var.names = c("TP","PP","TDP","SRP","NO3","DOC")
all.res <- res.df 

# Check all.res
all.res$source.nut <- as.factor(all.res$source.nut)
all.res$source.name <- as.factor(all.res$source.name)
summary(all.res$source.name)
summary(all.res$source.nut)

# get df with water chem data, geo data and discharge
chem.data <- read.csv2("data/StreamWaterChemistry.csv", stringsAsFactors = FALSE)
geo.data <- read.csv2("C:/doctorat/bdd/geo/LacAuDuc/zonalStats/watershed_properties_interregP0.csv")
geochem.data <- full_join(chem.data, geo.data, by = c("p.sampling"="NAME"))

geochem.data <- geochem.data %>% filter(p.sampling %in% list.select)#  %>% filter(!p.sampling %in% c("P12","P13",'P16','P18'))
geochem.data$PP <-  geochem.data$TP - geochem.data$TDP
geochem.data$NO3 <-  geochem.data$nitrate
geochem.data$date.sampling <- as.Date(geochem.data$date.sampling, format="%d/%m/%Y")
#geochem.data <-  geochem.data %>%
#filter(date.sampling > "2018-03-31") %>% # incomplete sampling
#filter(!date.sampling == '2018-10-02')  # not enough streams flowing

# Join discharge
debit <- read.csv2("C:/doctorat/bdd/eau/données/6.debit/debit1968_2018.csv")
debit$date <- as.Date(as.character(debit$Date),"%Y%m%d")
debit.tercile <- quantile(debit$Qls, probs = c(0.3333,0.6666))
debit <- debit %>% mutate(Qclass = ifelse(Qls <= debit.tercile[1], "Q<=Q33",
                                          ifelse(Qls <= debit.tercile[2], "Q33<Q<=Q67", "Q>Q67")))
debit$Qclass <- as.factor(debit$Qclass)
debit$Qclass <- relevel(debit$Qclass, "Q33<Q<=Q67")
debit$Qclass <- relevel(debit$Qclass, "Q<=Q33")

geochem.data <- left_join(geochem.data, debit[,c(2,3,4)], by = c("date.sampling" = "date"))
geochem.data$grouped <- FALSE

# compute mean by flow conditions and bind
geochem.grouped <- geochem.data %>% group_by(Qclass, p.sampling) %>%
  mutate(TP = mean(TP, na.rm = TRUE)) %>%
  mutate(TDP = mean(TDP, na.rm = TRUE)) %>%
  mutate(SRP = mean(SRP, na.rm = TRUE)) %>%
  mutate(DOC = mean(DOC, na.rm = TRUE)) %>%
  mutate(NO3 = mean(NO3, na.rm = TRUE)) %>%
  mutate(Qls = median(Qls, na.rm =TRUE)) %>%
  ungroup()
geochem.grouped <- distinct(geochem.grouped, Qclass, p.sampling, .keep_all = TRUE)
geochem.grouped$grouped <- TRUE
geochem.grouped <- geochem.grouped %>% mutate( date.sampling =
                                                 ifelse(Qclass == "Q<=Q33",'2001-01-01',ifelse(
                                                   Qclass == 'Q33<Q<=Q67', '2001-01-02','2001-01-03')) )
geochem.grouped$date.sampling <- as.Date(geochem.grouped$date.sampling)
geochem.data <- rbind.data.frame(geochem.data, geochem.grouped)

date.list <- unique(geochem.data$date.sampling)

# compute dates with less than 2 HW not flowing, and compute median
geochem.datecount <- geochem.data %>% 
  filter(grouped == FALSE) %>%
  group_by(date.sampling) %>%
  mutate(n = sum(!is.na(TP))) %>%
  dplyr :: select(date.sampling,n) %>%
  distinct() %>%
  ungroup()

dates.flowing <- filter(geochem.datecount, n>=(max(geochem.datecount$n))-2) %>%
  dplyr :: select(date.sampling) 

geochem.median <- geochem.data %>%
  filter(date.sampling %in% dates.flowing$date.sampling) %>%
  group_by(p.sampling) %>%
  mutate(TP = median(TP, na.rm = TRUE)) %>%
  mutate(TDP = median(TDP, na.rm = TRUE)) %>%
  mutate(SRP = median(SRP, na.rm = TRUE)) %>%
  mutate(DOC = median(DOC, na.rm = TRUE)) %>%
  mutate(NO3 = median(NO3, na.rm = TRUE)) %>%
  mutate(Qls = median(Qls, na.rm = TRUE)) %>%
  ungroup() %>% filter (date.sampling == date.list[1])

geochem.median$grouped <- TRUE
geochem.median$date.sampling <- '2001-01-04'
geochem.data <- rbind(geochem.data, geochem.median)
date.list <- unique(geochem.data$date.sampling)


ComputeCorrelation <- function(geoChemData, allRes, paraM, datE, sourceNut, sourceType, verbose = FALSE){
  chemData <- geoChemData %>% select(subcat.name = p.sampling, date.sampling, !!paraM, grouped) %>%
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
  corRes$grouped <- chemData$grouped[1]
  
  for(i in seq(1, nrow(corRes))){
    coefFacc <- as.numeric(corRes$coef.facc[i])
    coefFls <- as.numeric(corRes$coef.fls[i])
    
    indexDataCoef <- indexData %>% filter(coef.facc == coefFacc,
                                          coef.fls == coefFls)
    
    corr <- cor.test(x = as.numeric(chemData[,paraM]),
                     y = as.numeric(indexDataCoef$index.value),
                     method = 'spearman')
    
    corRes$R[i] <- corr$estimate
    corRes$p.value[i] <- corr$p.value
    
    if(verbose == TRUE){
      # data.corr <- as.data.frame(x = as.numeric(chemData[,paraM]), y =  as.numeric(indexDataCoef$index.value))
      # str(data.corr)
      str(chemData[,paraM])
      str(indexDataCoef$index.value)
      str(cor.test(x = as.numeric(chemData[,paraM]),
               y = as.numeric(indexDataCoef$index.value),
               method = 'spearman'))
    }
  }
  return(corRes)
}

# Initialisation du tableau
allCorr <- ComputeCorrelation(geochem.data, all.res, var.names[1], date.list[1], sources.nut[1], sources.type[1])
allCorr <- allCorr[0,]

tic("Compute correlations for param X source.nut X source.type X date")
for(param in var.names){
  tic(param)
  source.type <- "long-term"
  for(source.nut in sources.nut){
    tic(source.nut)
    #for(source.type in sources.type){
    for(date in date.list){

      corrCombin <- ComputeCorrelation(geochem.data, all.res, param, as.Date(date, origin = '1970-01-01'), source.nut, source.type)
      allCorr <- rbind(allCorr, corrCombin)
    }# end loop on dates
    # }# end loop on source type
    toc()
  }# end loop on source.nut
  toc() # 80 sec par indice pour un param
  
}# end loop on water quality param
toc() # 1400 sec

allCorr <- allCorr %>% distinct()

#join discharge on allCor and bestCor
allCorr$date <- as.Date(allCorr$date, origin = '1970-01-01')
allCorr <- inner_join(allCorr, debit %>% select(date, Qls))
allCorr <- inner_join(allCorr, geochem.data %>% select(date.sampling, Qclass), by = c("date" = "date.sampling"))
allCorr <- allCorr %>% distinct()                                     
allCorr$metric <- "mixed"

#Keep max R for all param*metric*source.nut*source.type
bestCor <- allCorr %>% group_by(param, date, source.nut, source.type) %>%
  slice(which.max(R)) %>% ungroup()

#Get lumped value for all param*metric*source.nut*source.type
lumpedCor <- allCorr  %>%
  filter(coef.facc == 0, coef.fls == 0) %>%
  mutate(metric = "lumped") %>%
  distinct()
# add it
allCorr <- rbind.data.frame(allCorr, lumpedCor)
bestCor <- rbind.data.frame(bestCor, lumpedCor)

# relevel metric
bestCor$metric <- factor(bestCor$metric,levels = c("lumped", "mixed"))


# Save Workspace
save.image(file = paste0(results.folder,"/3_Optimized_.RData"))



