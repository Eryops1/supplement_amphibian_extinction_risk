# Title: Evaluating the predicted extinction risk of living amphibian species with the fossil  record - Data processing
# Author: Melanie Tietje
  

library(ggplot2)
library(gridExtra)
library(grid)
library(scales)
library(gsubfn)
library(reshape)
library(fields)
library(rgdal)

theme=theme_set(theme_minimal()) 

# Intro
# This document contains the data processing and construction of the two datasets on extinct and living species.
# 
# This document is writen with the knitr package for the R statistical environment. The data processing can be recreated by running knitr::knit() on the underlying source code .Rmd file. The analysis and data files can be accessed in the Git repository: https://github.com/Eryops1/....
# 
# We used the following packages:
  
loadedNamespaces() # to make sure the list only contains packages used here, you might want to restart your R session


# Extinct species data

comb2 <- read.csv("comb2.csv")

# Data info: 
# mni_assumed in the comb2.csv dataset is coded as follows:
## y = info on specimens is missing, and the mni results in 1
## n = there is info on the specimens types and it`s either sure or likely (vertebrae) that it`s just one i

  
## Abundance
  
species <- unique(comb2$species)
mni_max <- c()
spec_max <- c()

for(i in 1:length(species)){
  temp <- comb2[comb2$species==species[i],]
  ifelse(all(is.na(temp$mni)), 
         mni_max <- c(mni_max, NA), mni_max <- c(mni_max, max(temp$mni, na.rm=TRUE)))
  ifelse(all(is.na(temp$specimen_no)),  
         spec_max <- c(spec_max, NA), spec_max <- c(spec_max, max(temp$specimen_no, na.rm=TRUE)))
}
abu <- data.frame(species, mni_max, spec_max)



### Abundance clustering
 # Abundance in the living dataset comes in categories, which is why we categorize our extinct data too.
extinct <- abu
extinct$complete_case <- complete.cases(extinct)


k2 <- kmeans(extinct[extinct$complete_case==TRUE,c("mni_max", "spec_max")], 2)
k3 <- kmeans(extinct[extinct$complete_case==TRUE,c("mni_max", "spec_max")], 3)
k4 <- kmeans(extinct[extinct$complete_case==TRUE,c("mni_max", "spec_max")], 4)
k5 <- kmeans(extinct[extinct$complete_case==TRUE,c("mni_max", "spec_max")], 5)
k6 <- kmeans(extinct[extinct$complete_case==TRUE,c("mni_max", "spec_max")], 6)

plot(c(2,3,4,5,6), 
     c(k2$betweenss/k2$totss, k3$betweenss/k3$totss, k4$betweenss/k4$totss, 
       k5$betweenss/k5$totss, k6$betweenss/k6$totss), type="b")

# Clustering the data (using both specimens and mni values) shows that between sum of squares divided 
# by total sum of squares (measure of the goodness of the classification k-means) there is much improvement 
# increasing the number of clusters from 2 to 3 and from 3 to 4, but not much from 4 clusters on. 
# Using 4 abundance categories is reasonable. 

set.seed(20)
mni_cluster <- kmeans(extinct[extinct$complete_case==TRUE,]$mni_max, 4)
abu_cat <- mni_cluster$cluster

# rename clusters to be intuitive
abu_cat[mni_cluster$cluster==order(mni_cluster$centers)[1]] <- 1
abu_cat[mni_cluster$cluster==order(mni_cluster$centers)[2]] <- 2
abu_cat[mni_cluster$cluster==order(mni_cluster$centers)[3]] <- 3
abu_cat[mni_cluster$cluster==order(mni_cluster$centers)[4]] <- 4

m <- cbind(extinct[extinct$complete_case==TRUE,], abu_cat) # 
abu <- merge(abu, m, all.x=TRUE)


## Geographic range size

dat <- comb2

## The maximum latitudinal range size of a species in a given stage
my.diff <- function(x){diff(range(x))}
sp <- c()
arr <- c()
species <- as.character(unique(dat$species))
for(i in 1:length(species)){
  temp <- dat[dat$species==species[i],]
  res <- tapply(temp$gp_lat, temp$ma_mid, my.diff)
  arr <- c(arr, res)
  sp <- c(sp, rep(species[i], length(res)))
}

arr <- as.data.frame(as.matrix(arr))
arr$age_mean <- rownames(arr)
rownames(arr) <- NULL
lat.range <- cbind(sp, arr)
lat.range <- as.data.frame(as.matrix(lat.range))
names(lat.range) <- c("species", "maxLATrange", "age_mean")


## The maximum GCD of a species in a given stage
species <- unique(dat$species)
maxGCD <- c()
species.name <- c()
slice <- c()
for(i in 1:length(species)){
  temp <- dat[dat$species==species[i],]
  slices <- unique(temp$slice)
  for(j in slices){
    temp2 <- temp[temp$slice==j,]
    temp3 <- rdist.earth(matrix(c(temp2$gp_lng, temp2$gp_lat), ncol=2), miles=FALSE)
    maxGCD <- c(maxGCD, max(temp3))
    species.name <- c(species.name, as.character(species[i]))
    slice <- c(slice, j)
  }
}
res <- data.frame(species.name, maxGCD, slice)
names(res)[1] <- "species"
range <- merge(lat.range, res)



## Get maximum ranges per species
lat_range <- aggregate(as.numeric(as.character(range$maxLATrange)), list(range$species), max)
names(lat_range)[2] <- "lat_range"
GCD <- aggregate(range$maxGCD, list(range$species), max)
names(GCD)[2] <- "gcd"

maxrange <- merge(lat_range, GCD)
names(maxrange)[1] <-"species"

trait <- merge(abu, maxrange, all.x=TRUE)


## Body size

duration <- unique(comb2[,c("species", "ma_range")])
trait <- merge(trait, duration, all.x=FALSE) # add duration

bodysize <- read.csv("bodysize.csv") # loads body size data

give.n <- function(x){
  return(c(y = mean(x)+50, label = length(x)))
}

ptl <- ggplot(bodysize, aes(y=tl_max, x=order))+
  geom_boxplot()+theme(axis.text.x = element_text(angle = 45, hjust = 1))+stat_summary(fun.data=give.n, geom = "text")
psvl <- ggplot(bodysize, aes(y=svl_max, x=order))+
  geom_boxplot()+theme(axis.text.x = element_text(angle = 45, hjust = 1))+stat_summary(fun.data=give.n, geom = "text")
pasl <- ggplot(bodysize, aes(y=asl, x=order))+
  geom_boxplot()+theme(axis.text.x = element_text(angle = 45, hjust = 1))+stat_summary(fun.data=give.n, geom = "text")
grid.arrange(ptl, psvl, pasl, ncol=3)


# Relation of absolute skull length (asl) to total length (tl) and snout-vent-length (svl). Grey shadows are 95% confidence intervals.

grid.arrange(
  ggplot(bodysize, aes(x=asl, y=svl_max))+geom_point(alpha=.3)+xlim(0,400)+geom_smooth(method="lm"),
  #ggplot(bodysize, aes(x=svl_max, y=tl_max))+geom_jitter(alpha=.3),
  ggplot(bodysize, aes(x=asl, y=tl_max))+geom_point(alpha=.3)+xlim(0,600)+geom_smooth(method="lm"),
  ncol=2)

lm1 <- lm(asl~0+svl_max, bodysize) # svl_max = predictor, asl = response
lm2 <- lm(asl~0+tl_max, bodysize) # same thing
# make sure the intercept is 0:0 to not allow for negative bodysizes
summary(lm1)
summary(lm2)

# Predicting ASL from TL and SVL
new <- data.frame(svl_max = bodysize$svl_max)
asl_pred <- predict.lm(lm1, newdata=new) # skull length data predicted using svl

new <- data.frame(tl_max = bodysize$tl_max)
asl_pred2 <- predict.lm(lm2, newdata=new) # skull length data predicted using tl

temp <- data.frame(asl_pred, asl_pred2, bodysize$asl)
cor.test(temp[,2], temp[,3])  # correlaton of asl values predicted by tl and svl

# Combine to one asl_pred
temp$asl_fin <- temp$bodysize.asl
for(i in 1:nrow(temp)){
  if(is.na(temp$bodysize.asl[i])){
    # check if others are not emtpy:
    if(!is.na(temp$asl_pred2[i])){temp$asl_fin[i] <- temp$asl_pred2[i]}
    if(!is.na(temp$asl_pred[i])){temp$asl_fin[i] <- temp$asl_pred[i]}
  }
}
bodysize$asl_pred <- temp$asl_fin

# Predicting SVL from ASL (TL is not possible as no intersection here)
lm1.1 <- lm(svl_max~asl_pred, bodysize) # svl_max = predictor, asl = response
new2 <- data.frame(asl_pred = bodysize$asl_pred) # asl_pred = original data + predicted from tl
svl_pred <- predict.lm(lm1.1, newdata=new2) # svl length data predicted using svl


temp <- data.frame(asl_pred, asl_pred2, bodysize$asl, bodysize$svl_max, svl_pred)
cor.test(temp[,2], temp[,3])  # correlaton of asl values with real values
cor.test(temp[,4], temp[,5])  # correlaton of svl values with real values

# add the predicted data
temp$svl_fin <- temp$bodysize.svl_max
for(i in 1:nrow(temp)){
  if(is.na(temp$bodysize.svl_max[i])){
    # check if others are not emtpy:
    if(!is.na(temp$svl_pred[i])){temp$svl_fin[i] <- temp$svl_pred[i]}
  }
}
bodysize$svl_pred <- temp$svl_fin



cor.test(bodysize$asl_pred, bodysize$svl_pred)


hist(bodysize$svl_pred, breaks=20)
#bodysize$asl_pred <- rescale(bodysize$asl_pred)
ggplot(bodysize, aes(x=asl_pred, fill=order))+geom_histogram()
ggplot(bodysize, aes(x=svl_pred, fill=order))+geom_histogram()


trait <- merge(trait, bodysize[c("species", "asl_pred", "svl_pred", "order")])

# Correlation is strong between predicted asl values and tl.
# By combining those original and predicted values, I get body size data for `r table(is.na(trait$asl_pred))[[1]]` species (SVL and ASL).




## Climatic preference (minimum and mean latitude of species occurrences)

# Mean latitudinal range
temp <- as.data.frame(tapply(comb2$gp_lat, comb2$species, mean, na.rm=TRUE))
temp$species <- rownames(temp)
names(temp)[1] <- "mean_lat"
trait <-merge(trait, temp, all.x=TRUE)

# Minimum latitude
temp <- as.data.frame(tapply(comb2$gp_lat, comb2$species, min, na.rm=TRUE))
temp$species <- rownames(temp)
names(temp)[1] <- "min_lat"
trait <-merge(trait, temp, all.x=TRUE)

# lat range
# my_range <- function(x)(max(x)-min(x))
# temp <- as.data.frame(tapply(comb2$gp_lat, comb2$species, my_range))
# temp$species <- rownames(temp)
# names(temp)[1] <- "lat_range"
# dat <-merge(dat, temp, all.x=TRUE)


## Add taxonomy

tax <- read.csv("taxonomy.csv")
trait <- merge(trait, unique(tax[,c("species", "order", "family", "genus")]), all.x=TRUE)






## Write file for modelling extinct species duration
# This is the file that`s being read into *Tietje_Roedel_2017_model_building.Rmd*
names(trait)[names(trait)=="svl_pred"] <- "svl"
trait <- droplevels(trait)
write.csv(trait, file="model_data.csv", row.names = FALSE)




















# Living species data

## Bodysize

mar <- read.csv("Hirschfeld_Roedel_bodysize_data.csv")
ead <- read.csv("European_Amphibians_Database3.csv")
rul <- read.csv("Ruland_Jeschke_bodysize_data.csv")

cau <- read.csv("bodysize_caudata.csv")

mar$species_binom <- paste(mar$Genus, mar$Species)
rul$species_binom <- paste(rul$genus, rul$species)

rul <- rul[rul$SVL!=99999,] # exclude NAs from Ruland & Jeschke data

species <- unique(c(rul$species_binom, mar$species_binom, as.character(ead$Species)))
dat <- data.frame(species, rep(NA, length(species)))
names(dat) <- c("species", "SVL")

dat <- merge(dat, mar[,c("species_binom", "SVL_f")],
             by.x="species", by.y="species_binom", all.x=TRUE)

dat <- merge(dat, rul[,c("species_binom", "SVL")],
             by.x="species", by.y="species_binom", all.x=TRUE) 

dat <- merge(dat, ead[,c("Species", "Adult.snout.to.vent.length")],
             by.x="species", by.y="Species", all.x=TRUE) 


names(dat) <- c("species", "SVL", "SVL_f_Mareike", "SVL_rul", "SVL_ead")
cor.test(dat$SVL_f_Mareike, dat$SVL_rul)
cor.test(dat$SVL_ead, dat$SVL_rul)

ggplot(dat, aes(x=SVL_rul, y=SVL_f_Mareike))+
  geom_point()+
  geom_point(aes(y=SVL_ead), col="blue")


dat$SVL <- dat$SVL_f_Mareike
for(i in 1:nrow(dat)){
  if(is.na(dat$SVL[i]))dat$SVL[i] <- dat$SVL_ead[i]
  if(is.na(dat$SVL[i]))dat$SVL[i] <- dat$SVL_rul[i]  
}

keeps <- c("species", "SVL")
dat <- dat[keeps]

# add caudata
cau <- na.omit(cau[,c("species", "SVL")])
any(dat$species %in% cau$species) # duplicates?
dat <- dat[-which(dat$species %in% cau$species),] # remove duplicates from the literature data and use amphibiaweb data
dat <- rbind(dat, cau)

# check for duplicates again
any(duplicated(dat$species))





## Abundance

# Collection of abundance data is done via text mining on the IUCN webpage. This takes quite long, therefore this script is set to read in the RData file from previous mining. 
# 
# IUCN red list population descriptions get scraped for keywords. All keywords were saved for each species. Keywords were categorized and species assigned according to their keywords. In cases of more than one keyword, I follow the conservative approach and use the lesser category. Example: common can be included as well in the description "the species is not common". Problems only occur with "common" and "abundant". This category gets assigned first. To account for the problem of descriptions sometimes including things like "the species is uncommon in xyz, but has large popuations in suitable habitats", categorys are assigned in a increasing order. Example: "very common" causes the species to end up in category 4, besides other metions of other keywords.

# Extract the species ID:
iucn <- read.csv("iucn_export-amphibia-03feb2017.csv")
iucn$species_binom <- paste(iucn$Genus, iucn$Species)
# reduce to species with body size data
id <- iucn$Species.ID[iucn$species_binom %in% dat$species]

keywords <- c("small", "small range", "large", "large range", "abundant", 
              "very abundant", "fairly abundant", "common", "fairly common",
              "rare", "locally common", "uncommon", "moderately abundant",
              "very common", "not common", "fewer than", "not abundant", "large population", "small population")



#library(XML)
# abu <- list()
# for(i in 1:length(id)){
#    # extract info from HTML
#    url <- paste0("http://www.iucnredlist.org/details/", id[i], "/0")
#    doc <- htmlParse(url)
#    doct <- xpathSApply(doc, "//*[@id='data_factsheet']/table[5]/tr[1]/td[2]", xmlValue)
#    # extract the xpath value, replace the double quotes with single quotes, and remove the tbody tag(!) from the path
#    abus <- c()
#    for(j in 1:length(keywords)){
#      if(grepl(keywords[j], doct)){
#        abus <- c(abus, keywords[j])
#      }
#    }
#    if(length(abus)==0){abu[[i]] <- "no keywords"}else{abu[[i]] <- abus}
#    print(i)
# }
# 
# # convert to dataframe
# #plyr::ldply(abu, rbind)
# mat <- matrix(ncol=length(keywords), nrow=length(abu))
# for(i in 1:length(abu)){
#   fill <- which(keywords %in% abu[[i]])
#   mat[i,fill] <- 1
# }
# mat <- as.data.frame(mat)
# names(mat) <- keywords
# mat <- cbind(mat, id)

#saveRDS(mat, file = "iucn_abundance_scraping.RData")
mat <- readRDS("iucn_abundance_scraping.RData")

# assign category for the species. cat3 needs conditioning
cat1 <- c("not common", "fewer than", "rare", "uncommon", "small", "small population", "not abundant")
cat2 <- c("fairly abundant", "fairly common", "moderately abundant")
cat3 <- c("abundant", "common", "large", "large populations")
cat4 <- c("very abundant", "very common")

abu_cat <- rep(NA, nrow(mat))
for(i in 1:nrow(mat)){
  temp <- mat[i,c(1:19)]
  if(all(is.na(temp))){abu_cat[i] <- NA}else{
    temp2 <- which(!is.na(temp))
    if(any(names(temp)[temp2] %in% cat3)){abu_cat[i] <- 3}
    # important: cat 3 has to be filled in first because it needs to be possible to replace this by other categories
    if(any(names(temp)[temp2] %in% cat1)){abu_cat[i] <- 1}
    if(any(names(temp)[temp2] %in% cat2)){abu_cat[i] <- 2}
    if(any(names(temp)[temp2] %in% cat4)){abu_cat[i] <- 4}
  }
  
}
iucn_abundance <- data.frame(abu_cat, mat$id)
iucn_abundance <- merge(iucn_abundance, iucn[,c("Species.ID", "species_binom")],
                        by.x = "mat.id", by.y = "Species.ID")

dat <- merge(dat, iucn_abundance[,c("abu_cat", "species_binom")], 
             by.x="species", by.y = "species_binom", all.x=TRUE)





## Geo range

# The following uncommented part of the script reads in the original shapefiles
# which can be obtained from the IUCN Red List website. The shapefiles are >100MB
# which is why I use the processed datafile generated by the following script 
# instead of the  actual shapefile. If you require to, the shapefiles can be obtained
# from http://www.iucnredlist.org/technical-documents/spatial-data (requires to generate
# a free account).

# geo <- readOGR("shapefiles_IUCN/AMPHIBIANS.shp")
#
# iucn.species <- sort(unique(geo$binomial))
# iucn.species <- iucn.species[iucn.species %in% dat$species]#
#
# shape.area <- c()
# threat.stat <- c()
# for(i in 1:length(iucn.species)){
#   temp <- geo[geo$binomial==iucn.species[i],]
#   shape.area <- c(shape.area, sum(temp$shape_Area)) # this combines all shapes, if multiple for one species
#   threat.stat <- c(threat.stat, as.character(unique(temp$rlcategory)[[1]]))
#   if(i %in% seq(10,2000,100)) print(i)
# }
# shape.area.extant <- data.frame(iucn.species, shape.area, threat.stat)
# shape.area.extant <- readRDS("shape.area.extant.RData")
#
#ggplot(shape.area.extant, aes(factor(threat.stat), shape.area))+
#  geom_boxplot(varwidth = TRUE)+
#  scale_y_log10()+
#  scale_x_discrete("IUCN Red List category")


## reduce the extraction to the needed species to reduce calculation time (still takes 10 minutes!) 
## extract all coordinates from the polygons

# system.time(
# for(i in 1:length(iucn.species)){ #
#  chapi <- geo[geo$binomial==iucn.species[i],] # polygons for single species
#  for(j in 1:length(chapi)){ # loop over all polygons from that species
#    coords <- unique(round(chapi@polygons[[j]]@Polygons[[1]]@coords, 2)) # round that values to speed it up!
#    if(j==1){mat <- coords}else{mat <- rbind(mat, coords)} # safe coordinates into matrix
#  }
#  if(i==1){fin <- data.frame(rep(iucn.species[i]), mat)}else{fin <- rbind(fin, data.frame(rep(iucn.species[i]), mat))} # # safe matrix for each species
#  if(i %in% seq(10,2000,10)) print(i)
# }
# )
# library(beepr)
#beep("treasure")
#saveRDS(fin, file="coords_from_shapefiles.RData")

fin <- readRDS("coords_from_shapefiles.RData")
names(fin) <- c("species", "lon", "lat")

# Lat range
iucnspecies <- unique(fin$species)
lrange <- matrix(ncol=2, nrow=length(iucnspecies), data=NA)
lrange <- as.data.frame(lrange)
for(i in 1:length(iucnspecies)){
  temp <- fin[fin$species==iucnspecies[i],]
  lrange[i,1] <- abs(diff(c(max(temp$lat), min(temp$lat)))) # get latitudinal range
  lrange[i,2] <- as.character(iucnspecies[i])
}
names(lrange) <- c("lat.range", "species")
dat <- merge(dat, lrange, all.x=TRUE)


# GCD
maxGCD_extant <- c()
system.time(
  for(i in 1:length(iucnspecies)){
    temp <- fin[fin$species==iucnspecies[i],]
    temp_reduced <- matrix(c(temp$lon, temp$lat), ncol=2)
    temp_reduced <- unique(round(temp_reduced)) # further reduction of precision due to limited memory allocation
    gcd <- rdist.earth(temp_reduced, miles=FALSE)
    maxGCD_extant <- c(maxGCD_extant, max(gcd)) 
    if(i %in% seq(10,2000,100)) print(i)
  }
)
res <- data.frame(iucnspecies, maxGCD_extant)
dat <- merge(dat, res, by.x="species", by.y="iucnspecies", all.x=TRUE)

# Mean and min lat
mean_lat <- c()
min_lat <- c()
for(i in 1:length(iucnspecies)){
  temp <- fin[fin$species==iucnspecies[i],]
  mean_lat <- c(mean_lat, mean(temp$lat))
  min_lat <- c(min_lat, min(temp$lat))
  if(i %in% seq(10,2000,100)) print(i)
}
mm_lat_range <- data.frame(mean_lat, min_lat, iucnspecies)

# Merge into dat
dat <- merge(dat, mm_lat_range, by.x="species", by.y="iucnspecies", all.x=TRUE)

# saveRDS(dat, file="dat_after_adding_shapefile_data.RData")
# dat <- readRDS("dat_after_adding_shapefile_data.RData")


## Extinction risk and taxonomy

iucn$species_binom <- paste(iucn$Genus, iucn$Species)
dat <- merge(dat, iucn[,c("species_binom", "Red.List.status", "Order", "Family", "Genus")],
             by.x="species", by.y="species_binom", all.x=TRUE)

write.csv(dat, file="living_data.csv")




# Data summary & Stratigraphy
## Temporal resolution
stages <- unique(comb2$stage)
stage_length <- numeric(length=length(stages))
no_occurrences <- numeric(length=length(stages))
for(i in 1:length(stages)){
  # get stage lengths
  temp <- comb2[comb2$stage==stages[i],]
  res <- unique(temp$ma_max - temp$ma_min)
  stage_length[i] <- res
  # get number of occurrences per stage
  no_occurrences[i] <- nrow(temp)
}
strat <- data.frame(stage_length, no_occurrences)

psl <- ggplot(strat, aes(x=stage_length))+
  geom_histogram(binwidth = 2)+
  scale_x_continuous("stage length in million years")+
  ggtitle("a")

comb2$stage_length <- comb2$ma_max - comb2$ma_min
pslno <- ggplot(strat, aes(x=stage_length, y=no_occurrences))+
  geom_point()+
  scale_y_continuous("number of occurrences")+
  scale_x_continuous("stage length in million years")
pslno <- pslno+ggtitle("b")

strat.plot <- grid.arrange(psl, pslno, ncol=2)
#save(psl, pslno, file="stratplot.RData")
ggsave("stratplot.jpg", strat.plot, dpi = 300, height=3, width=7)

median(stage_length)

## Data summary
### Occurrence map (contemporary coordinates)
library(maps)
#png("occurrence_map.png", height=5, width=8, res = 300)
map(fill = TRUE, col=c("grey90"), mar=c(0,0,0,0))
points(comb2$lngdec, comb2$latdec, pch=20, col="red")
#dev.off()

# creates a dataframe for the supplement Data summary
ds <- as.data.frame(table(comb2$species))
psych::describe(ds$Freq)


######## END ########################################################