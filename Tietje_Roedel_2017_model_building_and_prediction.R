# Title: Evaluating the predicted extinction risk of living amphibian species 
# with the fossil record SOM - Model building
# Author: Melanie Tietje
#
# Running this file for the first time might take about 30 minutes as cross-
 # validation and bootstrapping take place for several models. The process 
 # gets faster the second time as the most time-consuming computations will 
# be stored as RData files and then simply reloaded.

rm(list = ls()) # clear your workspace!

library(ggplot2)
library(gridExtra)
library(gsubfn)
library(randomForest)
library(mice)
library(reshape)
library(mgcv)
library(caret)
library(gbm)
library(simpleboot)


theme=theme_set(theme_minimal()+
                  theme(axis.text=element_text(size=7),
                        axis.title=element_text(size=9))) 

# Intro
# This document contains the model building, comparison, and prediction on living species.
# 
# This document is writen with the knitr package for the R statistical environment. The data processing can be recreated by running knitr::knit() on the underlying source code .Rmd file. The analysis and data files can be accessed in the Git repository: https://github.com/Eryops1/....
# 
# We used the following packages:

loadedNamespaces() # to make sure the list only contains packages used here, you might want to restart your R session


extinct.raw <- read.csv("model_data.csv")

extinct <- extinct.raw[,c("ma_range", "lat_range", "gcd", "svl", "mean_lat",
                 "min_lat", "abu_cat")]
extant <- read.csv("living_data.csv")
extant <- extant[,c("species", "SVL", "abu_cat", "lat.range", "maxGCD_extant", "mean_lat", "min_lat",
                    "Red.List.status", "Order", "Family", "Genus")]

# rename to match extinct dataset variable names
names(extant) <- c("species", "svl", "abu_cat", "lat_range", "gcd", "mean_lat", "min_lat",
                                 "Red.List.status", "order", "family", "genus")
extant <- na.omit(extant)



# Data imputation

# Data was imputed for missing cases of bodysize (svl) and abundance categories (abu_cat).Imputation was 
# performed using the mice package using the average results of 50 imputation repetitions. Further description
# of the procedure and  imputation validation can be found in the supplement file
# *Tietje_Roedel_2017_supplement_figures_tables_and_modelling_output.Rmd*. 

## load or create imputation
extinct.imp <- extinct.raw[,c("ma_range", "lat_range", "gcd", "svl", "mean_lat",
                              "min_lat", "abu_cat")]
extinct.imp$abu_cat <- as.factor(extinct.imp$abu_cat) # set abundance as factor

if(file.exists("imputed_data_extinct.RData")==TRUE){
  load("imputed_data_extinct.RData")}else{
    imputation_rf_svl <- matrix(ncol=50, nrow=nrow(extinct.imp))
    imputation_rf_abu <- matrix(ncol=50, nrow=nrow(extinct.imp))
    
    # Using randomForest algorithm
    for(i in 1:50){
      imputed <- rfImpute(extinct.imp, extinct$ma_range)
      imputation_rf_svl[,i] <- imputed$svl
      imputation_rf_abu[,i] <- imputed$abu_cat
    }
    
    # Using multivariate imputation by chained equations
    imputed <- mice(extinct.imp, maxit=20, m=50, print=FALSE, method = c("", "", "", "pmm", "", "", "polr")) # polr
    # use default imputation method for bodysize and polr instead of default polyreg (abu is an ordered factor)
    plot(imputed, c("svl", "abu_cat")) # check convergence of the predictor algorithm
    
    # combine the multiple iterations into one mean imputed value
    temp <- complete(imputed, action="broad", include=TRUE)
    imputation_mice_svl <- temp[,grepl(pattern = "svl.[1-9]", x = names(temp))] # grab all imputed svl columns
    imputation_mice_abu <- temp[,grepl(pattern = "abu_cat.[1-9]", x = names(temp))] # grab all imputed abu columns
    
    save(imputation_rf_abu, imputation_rf_svl,
         imputation_mice_abu, imputation_mice_svl,
         imputed,
         file="imputed_data_extinct.RData")
    }


# Add the imputed data
 ## the imputed data from mice needs not sorting in imputed and unimputed: the original values unchanged and are included
extinct.imp$svl <- rowMeans(imputation_mice_svl)
extinct.imp$abu_cat  <- as.factor(apply(imputation_mice_abu, 1, 
                                        function(x) names(sort(table(x), decreasing = TRUE))[1]))
extinct$abu_cat <- as.factor(extinct$abu_cat) # to keep the datasets comparable, transform extinct$abu_cat to factor too

extinct.imp <- cbind(extinct.imp, extinct.raw[,c("species", "order", "family", "genus", "complete_case")])
imp <- extinct.imp[,c("ma_range", "lat_range", "gcd", "svl", "mean_lat",
                 "min_lat", "abu_cat")]









# Generalized additive model

# GAM was fitted to the data to allow for non-linear relationships between the predictor variables and the species duration. 
# Normality (qq plot) and homogeneity are violated (residuals vs. fitted values (linear predictor)).
# 
# method and select of the mgcv::gam() are being chosen using the caret::train() function.

set.seed(432)
mySeeds <- sapply(simplify = FALSE, 1:31, function(u) sample(10^4, 5))
fitControl <- trainControl(method = "repeatedcv", number = 10,
                           repeats = 3, seeds = mySeeds,
                           returnResamp = "final",
                           savePredictions = "final") ## 3 times 10-fold CV
imp_for_gam <- imp
imp_for_gam$abu_cat <- as.numeric(as.character(imp_for_gam$abu_cat))
gamTrain <- train(ma_range ~ ., 
#                  knots=list(c(10,10,10,10,10,4)), 
                  data = imp_for_gam, method = "gam", 
                  trControl = fitControl, 
                  verbose = FALSE,
                  na.action = na.pass,
                  importance=TRUE)
gamTrain
summary(gamTrain$finalModel)
ggplot(gamTrain)
predict.gam <- predict.gam(gamTrain$finalModel, newdata = extant)
res <- extant
res$Red.List.status <- ordered(res$Red.List.status, levels = c("DD", "LC", "NT", "VU", "EN", "CR", "EW", "EX"))
res$predict.gam <- predict.gam
kruskal.test(res$predict.gam, res$Red.List.status)
pairwise.wilcox.test(res$predict.gam, res$Red.List.status, p.adjust.method = "fdr")

(gam.fin <- ggplot(res, aes(x=Red.List.status, y=log(predict.gam)))+
    geom_boxplot(outlier.colour=NULL, colour="#FF9999", fill="#FF9999")+
    scale_y_continuous("predicted duration (log ma)")+
    stat_summary(geom = "crossbar", 
                 width=0.65, fatten=.5, 
                 color="black", 
                 fun.data = function(x){return(c(y=median(x), ymin=median(x), ymax=median(x)))})
)  

res2 <- melt(res[,c("Red.List.status", "predict.gam")])
kruskal.test(res$predict.gam, res$Red.List.status)
pairwise.wilcox.test(res$predict.gam, res$Red.List.status, p.adjust.method = "fdr")











## Random Forest
# Set seeds and parameters for cross validation for parameter tuning with caret::train()
set.seed(432)
mySeeds <- sapply(simplify = FALSE, 1:31, function(u) sample(10^4, 5))

fitControl <- trainControl(method = "repeatedcv", number = 10,
                           repeats = 3, seeds = mySeeds,
                           returnResamp = "final",
                           savePredictions = "final") ## 3 times 10-fold CV
rfFit1 <- train(ma_range ~ ., data = imp, method = "rf", 
                trControl = fitControl, verbose = FALSE,
                na.action = na.pass,
                tuneGrid=data.frame(mtry=c(2,3,4,5,6)),
                importance=TRUE)
rfFit1
rfFit1$finalModel
ggplot(rfFit1)

plot(rfFit1$finalModel)
varImpPlot(rfFit1$finalModel)

### PREDICTION ####
if(class(extant$abu_cat)=="integer"){extant$abu_cat <- as.factor(extant$abu_cat)}
predict.rf <- predict(rfFit1, newdata = extant)
res$predict.rf <- predict.rf

res2 <- melt(res[,c("Red.List.status", "predict.gam", "predict.rf")])

ggplot(res2, aes(x=Red.List.status, y=log(value)))+
 geom_boxplot(aes(fill=variable))
#   geom_boxplot(outlier.colour=NULL, aes(fill=variable, colour=variable))+
#  scale_y_continuous("predicted duration (log ma)")+
#  stat_summary(geom = "crossbar", 
#               width=0.65, fatten=1, 
#               color="black", 
#               fun.data = function(x){return(c(y=median(x), ymin=median(x), ymax=median(x)))})

kruskal.test(res$predict.rf, res$Red.List.status)
pairwise.wilcox.test(res$predict.rf, res$Red.List.status, p.adjust.method = "fdr")








## Generalized boosted model
# Generalized boosted models are regular regression models that are being boosted, means that they 
# are being repeated according to the error parameters of the previous model, again and again. The 
# cross validation can be performed by classic division of data into training and evaluation set. 
# These models do not care about missing data, that`s a plus. They make minimal assumptions about 
# the parameters that are being entered, however collinearity can be a problem, GBMs proof to be 
# relatively tolerant to multicolinearity for prediction.

# parameter tuning using caret
gbmControl <- trainControl(method = "repeatedcv", number = 10,
                           repeats = 3, savePredictions = "final",
                           returnResamp ="final")

gbmGrid <-  expand.grid(interaction.depth = c(1, 2, 3), 
                        n.trees = (1:10)*50, 
                        shrinkage = c(0.1, 0.01, 0.001),
                        n.minobsinnode = c(5, 10, 15))
set.seed(825)
gbmFit1 <- train(ma_range ~ ., data = imp, method = "gbm", 
                 trControl = gbmControl, verbose=FALSE,
                 tuneGrid=gbmGrid)
gbmFit_no_imputation <- train(ma_range ~ ., data = na.omit(extinct), method = "gbm", 
                              trControl = gbmControl, verbose=FALSE,
                              tuneGrid=gbmGrid)

for(i in 1:8){
  temp <- plot.gbm(gbmFit_no_imputation$finalModel, i, return.grid = TRUE)
  if(i==1){marg <- temp}else{marg <- cbind(marg, temp)}
}
names(marg)[c(2,4,6,8,10,12,14,16)] <- c("y_lat", "y_gcd", "y_size", "y_latitude",
                                         "y_min_lat", "y_abu2", "y_abu3", "y_abu4")

a=ggplot(marg, aes(x=lat_range, y=y_lat))+
  geom_line()+xlab("latitudinal range")+ylab(" ")
b=ggplot(marg, aes(x=gcd, y=y_gcd))+
  geom_line()+xlab("geographic range size")+ylab(" ")+
  scale_x_continuous(breaks = c(4000,8000,12000))
c=ggplot(marg, aes(x=svl, y=y_size))+
  geom_line()+xlab("body size")+ylab(" ")
d=ggplot(marg, aes(x=mean_lat, y=y_latitude))+
  geom_line()+xlab("mean latitude")+ylab(" ")
e=ggplot(marg, aes(x=min_lat, y=y_min_lat))+
  geom_line()+xlab("minimum latitude")+ylab(" ")
f=ggplot(marg, aes(x=abu_cat2, y=y_abu2))+
  geom_line()+xlab("abundance 2")+ylab(" ")
g=ggplot(marg, aes(x=abu_cat3, y=y_abu3))+
  geom_line()+xlab("abundance 3")+ylab(" ")
h=ggplot(marg, aes(x=abu_cat4, y=y_abu4))+
  geom_line()+xlab("abundance 4")+ylab(" ")
grid.arrange(a,b,c,d,e,f,g,h, ncol=4)

summary(gbmFit_no_imputation)

gbmFit1
summary(gbmFit1); summary(gbmFit_no_imputation)
gbm.perf(gbmFit1$finalModel)
ggplot(gbmFit1, nameInStrip = TRUE)
ggplot(as.data.frame(summary(gbmFit1)),
                  aes(x=var, y=rel.inf))+
  geom_col()

## Prediction
predict.gbm1 <- predict(gbmFit1, newdata = extant, n.trees = gbmFit1$bestTune$n.trees)
res$predict.gbm1 <- predict.gbm1

### GBM without imputation
predict.gbm_no_imp <- predict(gbmFit_no_imputation, newdata = extant,
                              n.trees = gbmFit_no_imputation$bestTune$n.trees)
res$predict.gbm_no_imp <- predict.gbm_no_imp

### GBM on species not classified under a range size criterium
iucn <- read.csv("iucn_export-amphibia-03feb2017.csv")
### Remove all species listed under B and D2
crit_b <- grep("B", iucn$Red.List.criteria)
crit_d2 <- grep("D2", iucn$Red.List.criteria)
crit_sum <- c(crit_b, crit_d2)
no_range <- iucn[-c(crit_b, crit_d2),]
no_range$species <- paste(no_range$Genus, no_range$Species)
table(no_range$Red.List.criteria=="")
extant_no_range <- extant[which(extant$species %in% no_range$species),]
extant_no_range$Red.List.status <- ordered(extant_no_range$Red.List.status,
                                           levels = c("DD", "LC", "NT", "VU", "EN", "CR", "EW", "EX"))

predict.gbm_no_range <- predict(gbmFit1, newdata = extant_no_range,
                                n.trees = gbmFit1$bestTune$n.trees)
extant_no_range$predict.gbm_no_range <- predict.gbm_no_range

ggplot(extant_no_range, aes(x=Red.List.status, y=predict.gbm_no_range))+
  geom_boxplot(outlier.colour=NULL, colour="grey50", fill="grey50",  outlier.alpha = .5)+
  scale_y_continuous("predicted duration (ma)")+
  scale_x_discrete("Red List category")+
  stat_summary(geom = "crossbar", 
               width=0.65, fatten=.5, 
               color="black", 
               fun.data = function(x){return(c(y=median(x), ymin=median(x), ymax=median(x)))})

####
ggplot(res, aes(x=Red.List.status, y=predict.gbm1))+
    geom_boxplot(outlier.colour=NULL, colour="grey50", fill="grey50",  outlier.alpha = .5)+
    scale_y_continuous("predicted duration (ma)")+
    scale_x_discrete("Red List category")+
    stat_summary(geom = "crossbar", 
                 width=0.65, fatten=.5, 
                 color="black", 
                 fun.data = function(x){return(c(y=median(x), ymin=median(x), ymax=median(x)))})
#   geom_boxplot()+
#   scale_y_log10("predicted duration (log ma)")

#### add to other model results
res2 <- melt(res[,c("Red.List.status", "predict.gam", "predict.rf", 
                    "predict.gbm1", "predict.gbm_no_imp")])
ggplot(res2, aes(x=Red.List.status, y=value))+
  geom_boxplot(aes(fill=variable))+
  scale_y_log10("predicted duration (log ma)")#+
#  facet_wrap(~variable)

kruskal.test(res$predict.gbm1, res$Red.List.status)
pairwise.wilcox.test(res$predict.gbm1, res$Red.List.status, p.adjust.method = "fdr")
tapply(res$predict.gbm1, res$Red.List.status, psych::describe)


# More cross validation, and taxonomic-group level mean bias
# Assessing the prediction quality is done by using only part of the data to 
# build the model and evaluate it on the remaining data.

## Taxonomic level bias bootstrap
size_modelset <- 1/2
fin <- c()
n <- 500
set.seed(346)
interval.range <- 10
interval.size <- 0.5

for(i in 1:n){
  sub <- extinct.imp[sample(seq(1:nrow(extinct.imp)), nrow(extinct.imp)*size_modelset, replace = FALSE),]
  eval <- extinct.imp[-as.numeric(row.names(sub)),]
  gbm.sub <- gbm(ma_range ~ lat_range+gcd+svl+mean_lat+min_lat+abu_cat
                 , data = sub,
                 n.trees=gbmFit1$finalModel$n.trees, 
                 shrinkage = gbmFit1$finalModel$shrinkage,
                 interaction.depth=gbmFit1$finalModel$interaction.depth,
                 distribution="gaussian")
  best.iter.sub <- gbm.perf(gbm.sub, method="OOB", plot.it=FALSE)
  pre <- predict.gbm(gbm.sub, eval, best.iter.sub)
    # manually assign every prediction smaller than 0 to be 0 for technical reasons
  pre[pre<0] <- 0
  
  # group the predicted durations into bins
  g <- findInterval(pre, seq(0, interval.range, interval.size))
  cc <- data.frame(prediction=pre, interval=g, species=eval$species, ma_range=eval$ma_range, order=eval$order)
  
  # now calculate what was the actual observed duration of the species in these bins
  medians_predicted <- tapply(cc$prediction, cc$interval, median)
  m.df <- data.frame(medians_predicted, interval=as.numeric(names(medians_predicted)))
  medians_observed <- tapply(cc$ma_range, cc$interval, median)
  m.df <- data.frame(m.df, medians_observed)
  temp <- merge(cc, m.df, by="interval", all.x=TRUE)
  temp <- data.frame(temp, run=rep(i, nrow(temp)))
  # add interval name
  temp$interval2 <- seq(0, interval.range, interval.size)[temp$interval]
  # save stuff 
  if(i==1){fin <- temp}else{fin <- rbind(fin, temp)}
  # print progress
   if(i %% 10==0) { cat(paste0("iteration: ", i, "\n")) }
} # the bootstrap

# Add the frequency per interval for weighting the regression:
weights <- data.frame(table(fin$interval2))
fin <- merge(fin, weights, by.x="interval2", by.y="Var1", all.x=TRUE)

ggplot(fin[!is.na(fin$order),], aes(x=interval2, y=medians_observed, weight=Freq))+
  geom_point(alpha=1/50)+
#  geom_bin2d()+
  scale_x_continuous("gbm predicted duration (ma)")+
  scale_y_continuous("median observed duration per interval (ma)")+
  geom_smooth()+
  geom_abline(col="darkgrey", lty=2)+
  facet_wrap(~order, scales="free")

# ggplot(fin[!is.na(fin$order),], aes(x=prediction, y=ma_range))+
#   geom_point(alpha=1/50)+
# #  geom_bin2d()+
#   scale_x_continuous("gbm predicted duration (ma)")+
#   scale_y_continuous("median observed duration per interval (ma)")+
#   geom_smooth()+
#   facet_wrap(~order, scales="free")


# Get the taxonomic group level bias for each prediction bin. Take the predicted duration value for a living species from the gbm.predict and extract the empirically observed duration from the gam fit
formula = y ~ s(x, bs = "cs")
gam.boot.salientia <- gam(data=fin[fin$order=="Salientia",], medians_observed~s(interval2, bs="cs"), weights=Freq)
gam.boot.urodela <- gam(data=fin[fin$order=="Urodela",], medians_observed~s(interval2, bs="cs"), weights=Freq)

res$predict.gbm1.anura_correction <-  predict(gam.boot.salientia, data.frame(interval2=res$predict.gbm1))
res$predict.gbm1.caudata_correction <-  predict(gam.boot.urodela, data.frame(interval2=res$predict.gbm1))

# plot the effect
par(mfrow=c(1,2))
plot(res$predict.gbm1.anura_correction[res$order=="ANURA"]~res$predict.gbm1[res$order=="ANURA"],
     xlim=c(0,15), main="Anura", xlab="gbm duration (ma)", ylab="corrected gbm duration (ma)")
abline(a=0, b=1, lty=2)
plot(res$predict.gbm1.caudata_correction[res$order=="CAUDATA"]~res$predict.gbm1[res$order=="CAUDATA"],
     xlim=c(0,15), main="Caudata", xlab="gbm duration (ma)", ylab="corrected gbm duration (ma)")
abline(a=0, b=1, lty=2)
par(mfrow=c(1,1))

# mean difference anura:
mean(res$predict.gbm1[res$order=="ANURA"]-res$predict.gbm1.anura_correction[res$order=="ANURA"])

# Combine corrected predictions
predict.gbm1.comb <- c()
for(i in 1:nrow(res)){
  if(res$order[i]=="ANURA"){predict.gbm1.comb <- c(predict.gbm1.comb, res$predict.gbm1.anura_correction[i])}else{
    predict.gbm1.comb <- c(predict.gbm1.comb, res$predict.gbm1.caudata_correction[i])
  }
}
res$predict.gbm1.comb <- predict.gbm1.comb
ggplot(res, aes(x=Red.List.status, y=predict.gbm1.comb))+
  geom_boxplot(varwidth=TRUE)+
  scale_y_log10()

# add to other model results
res2 <- melt(res[,c("Red.List.status", "predict.gam", "predict.rf",
                    "predict.gbm1", "predict.gbm_no_imp", "predict.gbm1.comb")])
ggplot(res2, aes(x=Red.List.status, y=value))+
  geom_boxplot(aes(fill=variable))+
  scale_fill_discrete(name = "model types")+
  scale_y_log10("predicted duration (log ma)")+
  facet_wrap(~variable)

kruskal.test(res$predict.gbm1.comb, res$Red.List.status)
pairwise.wilcox.test(res$predict.gbm1.comb, res$Red.List.status, p.adjust.method = "fdr")


### Subset models (Lissamphbia, No-singletons)

#### Lissamphibia

liss.raw <- extinct.imp[extinct.imp$order %in% c("Urodela", "Salientia", "Parabatrachia"),]
nrow(liss.raw)
liss <- liss.raw[,c("ma_range", "lat_range", "gcd", "svl", "mean_lat",
                 "min_lat", "abu_cat")]

set.seed(825)
gbmFit_liss <- train(ma_range ~ ., data = liss, method = "gbm", 
                 trControl = gbmControl, verbose=FALSE,
                 tuneGrid=gbmGrid)

gbmFit_liss
gbm.perf(gbmFit_liss$finalModel)
ggplot(gbmFit_liss, nameInStrip = TRUE)
ggplot(as.data.frame(summary(gbmFit_liss)),
                  aes(x=var, y=rel.inf))+
  geom_col()

## Prediction
predict.gbm.liss <- predict(gbmFit_liss, newdata = extant, n.trees = gbmFit_liss$bestTune$n.trees)
res$predict.gbm.liss <- predict.gbm.liss

res2 <- melt(res[,c("Red.List.status", "predict.gam", "predict.rf",
                    "predict.gbm1", "predict.gbm_no_imp", "predict.gbm1.comb", "predict.gbm.liss")])
ggplot(res2, aes(x=Red.List.status, y=value))+
  geom_boxplot(aes(fill=variable))+
  scale_fill_discrete(name = "model types")+
  scale_y_log10("predicted duration (log ma)")+
  facet_wrap(~variable)

kruskal.test(res$predict.gbm.liss, res$Red.List.status)
pairwise.wilcox.test(res$predict.gbm.liss, res$Red.List.status, p.adjust.method = "fdr")



#### No single-interval species
imp_nosingles <- imp[imp$ma_range>0,]
nrow(imp_nosingles)
set.seed(825)
gbmFit_nosingles <- train(ma_range ~ ., data = imp_nosingles, method = "gbm", 
                 trControl = gbmControl, verbose=FALSE,
                 tuneGrid=gbmGrid)
gbmFit_nosingles
predict.gbm.nosingles <- predict(gbmFit_nosingles, newdata = extant, 
                                 n.trees = gbmFit_nosingles$bestTune$n.trees)
res$predict.gbm.nosingles <- predict.gbm.nosingles
ggplot(res, aes(x=Red.List.status, y=predict.gbm.nosingles))+
  geom_boxplot()+
  scale_y_log10("predicted duration (log ma)")
kruskal.test(res$predict.gbm.nosingles, res$Red.List.status)
pairwise.wilcox.test(res$predict.gbm.nosingles, res$Red.List.status, p.adjust.method = "fdr")

# res2 <- melt(res[,c("Red.List.status", "predict.gam.train", "predict.rf",
#                     "predict.gbm1", "predict.gbm1.comb", "predict.gbm.liss",
#                     "predict.gbm.nosingles")])
# ggplot(res2, aes(x=Red.List.status, y=value))+
#   geom_boxplot(aes(fill=variable))+
#   scale_fill_discrete(name = "model types")+
#   scale_y_log10("predicted duration (log ma)")+
#   facet_wrap(~variable)









# Null model

# Our null hypothesis is that we expect no connection between traits and extinction risk. 
# A null model should represent a fit to a random dataset. However, the dataset does not 
# have to be entirely random as this would create species trait combinations that do not 
# make sense (being very rare but extreme widely distributed), just random for the response 
# variable "duration". 

null.data <- imp

set.seed(2)
null.data$ma_range <- sample(null.data$ma_range, replace=FALSE)

set.seed(825)
gbmFit.null <- train(ma_range ~ ., data = null.data, method = "gbm", 
                 trControl = gbmControl, verbose=FALSE,
                 tuneGrid=gbmGrid)
gbmFit.null
gbm.perf(gbmFit.null$finalModel)
ggplot(gbmFit.null, nameInStrip = TRUE)
ggplot(as.data.frame(summary(gbmFit.null)),
                  aes(x=var, y=rel.inf))+
  geom_col()

## Prediction
predict.gbm.null <- predict(gbmFit.null, newdata = extant, n.trees = gbmFit1$bestTune$n.trees)
res$predict.gbm.null <- predict.gbm.null

ggplot(res, aes(x=Red.List.status, y=predict.gbm.null))+
  geom_boxplot()+
  scale_y_log10("predicted duration (log ma)")

kruskal.test(res$predict.gbm.null, res$Red.List.status)
#pairwise.wilcox.test(res$predict.gbm.null, res$Red.List.status, p.adjust.method = "fdr")


## BOOTSTRAP
iter <- 50
if(file.exists("null-boot.RData")==TRUE){
  load("null-boot.RData")}else{
    # Bootstrap - 50 iterations take up to 20 minutes
    store_kw <- c()
    store_rmse <- c()
    store_r2 <- c()
    store_median <- matrix(ncol=8, nrow=iter)
    for(i in 1:iter){
      null.data <- imp
      null.data$ma_range <- sample(null.data$ma_range, replace=FALSE)
      gbm.null.boot <- train(ma_range ~ ., data = null.data, method = "gbm",
                             trControl = gbmControl, verbose=FALSE,
                             tuneGrid=gbmGrid)
      predict.gbm.null.boot <- predict(gbm.null.boot, newdata = extant, n.trees=100)#
      temp <- data.frame(prediction=predict.gbm.null.boot, group=extant$Red.List.status)
      store_median[i,] <- tapply(temp$prediction, temp$group, median)
      store_kw <- c(store_kw, kruskal.test(predict.gbm.null.boot, extant$Red.List.status)$p.value)
      store_rmse <- c(store_rmse, min(gbm.null.boot$results$RMSE))
      store_r2 <- c(store_r2, max(gbm.null.boot$results$Rsquared))
      # print progress
      if(i %% 10==0) { cat(paste0("iteration: ", i, "\n")) }
    }
    save(store_kw, store_rmse, store_r2, store_median, file="null-boot.RData")
  }

par(mfrow=c(2,2))
hist(store_kw, breaks=40, xlab="Kruskal-Wallis p-value")
hist(store_rmse, breaks=40, xlab="Minimum RMSE of gbm")
hist(store_r2, breaks=40, xlab="Maximum Rsquared of gbm")
boxplot(store_median, names=c("DD", "LC", "NT", "VU", "EN", "CR", "EW", "EX"),
        main=paste0("bootstrap null model prediction, n=", iter),
        ylab="Median duration predicted (log)", xlab="IUCN Red List categories")
par(mfrow=c(1,1))
boot.k <- melt(store_median)
kruskal.test(boot.k$value, boot.k$X2)

# How do model performance metrics compare with the null model values (CIs etc)
shapiro.test(store_rmse) # normal distribution?
shapiro.test(store_r2) # normal distribution?

## collect the model performance metrics
mods <- c("rfFit1", "gbmFit1", "gbmFit_liss", "gbmFit.null")
rmse.mins <- c()
r2.maxs <- c()
for(i in 1:length(mods)){
  rmse.mins <- c(rmse.mins, min(eval(parse(text=paste0(mods[i], "$results$RMSE")))))
  r2.maxs <- c(r2.maxs, max(eval(parse(text=paste0(mods[i], "$results$Rsquared")))))
}

## check standard deviations
rmse.mins %in% range(c(mean(store_rmse)-2*sd(store_rmse), mean(store_rmse)+2*sd(store_rmse)))
r2.maxs %in% range(c(median(store_r2)-2*mad(store_r2), median(store_r2)+2*mad(store_r2)))

# Construct a 95% CI for the Null model RMSE
ci.upper <- mean(store_rmse)+2*(sd(store_rmse)/sqrt(length(store_rmse)))
ci.lower <- mean(store_rmse)-2*(sd(store_rmse)/sqrt(length(store_rmse)))
rmse.mins < ci.lower# & rmse.mins < ci.upper

# Construct a parametric 95% CI for the Null model r2
shapiro.test(log(store_r2))
ci.upper <- mean(log(store_r2))+2*(sd(log(store_r2))/sqrt(length(store_r2)))
ci.lower <- mean(log(store_r2))-2*(sd(log(store_r2))/sqrt(length(store_r2)))
r2.maxs > exp(ci.upper)# & r2.maxs < exp(ci.upper)

## alternative:
# 20% trimmed mean bootstrap
set.seed(1234)
b1 <- one.boot(store_r2, mean, R=2000, tr=.1)
boot.ci(b1, type=c("perc", "bca"))
r2.maxs %in% boot.ci(b1, type=c("perc", "bca"))$bca

# None of the final model RMSEs except the null model is within or below the 95%CI 
# of the bootstrapped null model RMSE distribution. Rsquared metrics are as well 
# way above the 95% nonparametric CI of the bootstrapped null model R squared distribution. 










# Comparing models
res2 <- melt(res[,c("Red.List.status", "predict.gbm1", "predict.gbm.null")]) 
res_all <- melt(res[,c("Red.List.status", "predict.gbm1", "predict.gbm_no_imp", "predict.gbm.null",
                       "predict.gam", "predict.rf", "predict.gbm.liss",
                       "predict.gbm.nosingles", "predict.gbm1.comb")])
labels <- c(predict.gbm1="GBM", predict.gbm.null="GBM-NULL")
labels2 <- c(predict.gbm1="GBM", predict.gbm_no_imp="GBM-NO_IMP", predict.gbm.null="GBM-NULL",
             predict.gam="GAM", predict.rf="RF", predict.gbm.liss="GBM-LISS",
             predict.gbm.nosingles="GBM-NO_SINGLES", predict.gbm1.comb="GBM-CORR")


my_col <- rgb(39, 123, 183, 254, maxColorValue=255) #poster?
# figure3
png("figure3.png", res = 600, unit="mm", width=80, height=60)
(figure3 <- ggplot(res2, aes(x=Red.List.status, y=value))+
    geom_boxplot(colour=my_col, fill=my_col,  outlier.alpha = .5, outlier.stroke = 0)+
    scale_y_continuous("predicted duration (ma)")+
    scale_x_discrete("Red List category")+
    stat_summary(geom = "crossbar", 
                 width=0.65, fatten=.5, 
                 color="black", 
                 fun.data = function(x){return(c(y=median(x), ymin=median(x), ymax=median(x)))})+
  facet_wrap(~variable, labeller=labeller(variable = labels))+
  theme(legend.position="none",
  axis.text.x=element_text(size=6))
)
dev.off()
ggsave("figure3.pdf", width=80, height=60, units="mm")
ggsave("figure3_wider.pdf", width=110, height=60, units="mm")


res_all_plot <- res
res_all_plot$predict.gam <- sqrt(res_all_plot$predict.gam)
res_all_plot <- melt(res_all_plot[,c("Red.List.status", "predict.gbm1", "predict.gbm_no_imp", 
                                     "predict.gbm.null",
                       "predict.gam", "predict.rf", "predict.gbm.liss",
                       "predict.gbm.nosingles", "predict.gbm1.comb")])
ggplot(res_all_plot, aes(x=Red.List.status, y=value))+
    geom_boxplot(colour=my_col, fill=my_col,  outlier.alpha = .5, outlier.stroke = 0)+
    scale_y_continuous("predicted duration (ma)")+
    scale_x_discrete("Red List category")+
    stat_summary(geom = "crossbar", 
                 width=0.65, fatten=.5, 
                 color="black", 
                 fun.data = function(x){return(c(y=median(x), ymin=median(x), ymax=median(x)))})+
  facet_wrap(~variable, labeller=labeller(variable = labels2))+
  theme(legend.position="none")
#ggsave("prediction_results_all.pdf", width = 6, height = 7)

tapply(res$predict.gbm1, res$Red.List.status, psych::describe)


resamps <- resamples(list(GBM = gbmFit1,
                          GBM_NO_IMP = gbmFit_no_imputation,
                          GBM_LISS = gbmFit_liss,
                          RF = rfFit1,
                          GBM_NULL = gbmFit.null,
                          GAM = gamTrain,
                          GBM_NO_SINGLES = gbmFit_nosingles))
# gets all the 30 resamples from cross validation from the final model
summary(resamps)
rmse.cols <- resamps$values[,c(1, which(grepl("RMSE", names(resamps$values))))]
r2.cols <- resamps$values[,c(1, which(grepl("Rsquared", names(resamps$values))))]
rmse.cols <- melt(rmse.cols)
r2.cols <- melt(r2.cols)
rmse.cols$variable <- gsub("~RMSE", "", rmse.cols$variable)
r2.cols$variable <- gsub("~Rsquared", "", r2.cols$variable)
names(rmse.cols)[3] <- "RMSE"
names(r2.cols)[3] <- "Rsquared"

resamps.gg <- merge(rmse.cols, r2.cols)
resamps.gg <- melt(resamps.gg, variable_name="var2")

ggplot(resamps.gg, aes(x=variable, y=value))+
  geom_boxplot()+
  scale_y_log10()+scale_x_discrete("model")+
  facet_wrap(~var2, scales = "free")+
  coord_flip()

dat.temp <- subset(resamps.gg, variable=="GAM")
tapply(dat.temp$value, dat.temp$var2, summary)

# The plot shows the performance metric distributions of the final model 
# (the one with the best parameters) for each of the 30 resamples performed 
# during cross validation. In model$resample the RMSE and Rsquared for each 
# of the gbms (based on the optimal number of trees) is stored. The mean of 
# these 30 RMSEs and Rsquares is used to chose the final model (=the best tuning parameters).


# run a t-test on the metrics
difValues <- diff(resamps, adjustment = "fdr")
difValues

summary(difValues, round=2)
bwplot(difValues, layout = c(2, 1))

# this function collects the RMSE, Rsquared, RMSE SD and Rsquared SD from 
# the final model, representing the mean values achieved by the cross validation 
# of this model. (note: the values could be obtained as well by calculating the 
# means etc from the resamps.gg dataframe)
collect.metrics <- function(x){ # provide x as character
  if(is.character(x)==FALSE) stop("Model names have to be given as characters")
    for(i in 1:length(x)){
      y <- get(x[i])
      temp.first <- y$results[which.min(y$results$RMSE),]#[,c((ncol(y$results)-3):ncol(y$results))]
      temp <- temp.first[c("RMSE", "RMSESD", "Rsquared", "RsquaredSD")]
      row.names(temp) <- x[i]
      if(i==1){
        temp2 <- temp
      }else{temp2 <- rbind(temp2, temp)}
    }
    return(temp2)
  }

model_metrics <- collect.metrics(c("gbmFit1", "gbmFit_no_imputation", "gbmFit_liss", "gbmFit.null",
                  "rfFit1", "gbmFit_nosingles", "gamTrain"))
round(model_metrics[], 2)








# Potential misclassifications
names(res)
plot(res$Red.List.status, res$predict.gbm1)
#mod <- res[,c("Red.List.status","predict.gbm1")]
#mod <- dcast(mod, predict.gbm1 ~ Red.List.status)

identify.out <- function(x){
  upper.limit <- quantile(x)[[4]] + 1.5*IQR(x)
  upper.out <- which(x > upper.limit)
  lower.limit <- quantile(x)[[2]] - 1.5*IQR(x)
  lower.out <- which(x < lower.limit)
  res <- list(lower.out=lower.out, upper.out=upper.out)
  return(res)
}

out <- tapply(res$predict.gbm1, res$Red.List.status, identify.out)


## Identify the species
species.longer <- c(
as.character(res$species[res$Red.List.status=="DD"][out$DD$upper.out]),
as.character(res$species[res$Red.List.status=="VU"][out$VU$upper.out]),
as.character(res$species[res$Red.List.status=="EN"][out$EN$upper.out]),
as.character(res$species[res$Red.List.status=="CR"][out$CR$upper.out]))

species.shorter <- c(
  as.character(res$species[res$Red.List.status=="DD"][out$DD$lower.out]),
  as.character(res$species[res$Red.List.status=="EN"][out$EN$lower.out]))

res$outlier <- rep(NA, nrow(res))
res$outlier[which(res$species %in% species.longer)] <- "longer"
res$outlier[which(res$species %in% species.shorter)] <- "shorter"

ggplot(res, aes(x=Red.List.status, fill=outlier))+
  geom_bar(position=position_fill())

number.out <- table(res$outlier, res$Red.List.status)
number.out <- as.data.frame.matrix(number.out)
number.out <- rbind(number.out, as.numeric(table(res$Red.List.status)))
row.names(number.out) <- c("Longer", "Shorter", "Total")
number.out

##################

#### Data summary 
type <- c(rep("paleo", nrow(extinct.raw)), rep("living", nrow(extant)))
mean_latitude <- c(extinct.raw$mean_lat, extant$mean_lat)
latitude_comp_living_fossil <- ggplot(data.frame(type,mean_latitude), aes(x=type,y=mean_latitude))+
  geom_boxplot(varwidth = TRUE)+
  scale_y_continuous("mean latitude")

#ggplot(extinct, aes(x=ma_range))+
#  geom_histogram()+
#  scale_y_log10()

labels_importance_plot <- c("abundance 2", "abundance 3", "abundance 4", "great circle distance",                            "latitudinal range", "mean latitude", "minimum latitude", "body size")
ggplot(as.data.frame(summary(gbmFit1$finalModel, plotit=FALSE)), aes(x=var, y=rel.inf))+
  geom_col()+
  xlab("variable")+ylab("relative influence")+
  scale_x_discrete(labels=labels_importance_plot)+
 theme(axis.text.x = element_text(angle = 45, hjust = 1),
       axis.text = element_text(size = 9))
kable(as.data.frame(summary(gbmFit1$finalModel, plotit=FALSE)))

#################


## Save workspace
save.image("Tietje_Roedel_2017_model_prediction_workspace.RData")


########## END ############


