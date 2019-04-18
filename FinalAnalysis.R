#################################################
####   Script for analysing body-mass data   ####
#################################################

##### goals #####

# 1. analyse body mass ~ forewing length + family + (1|species)
# 2. build predictive model using training dataset and body mass ~ forewing length + family.crude + (1|species)
# 3. analyse sample-level biomass ~ abundance, richness, diversity
# 4. use data from RIS Heslington trap to analyse annual biomass ~ abundance, richness, diversity, and trends in all four


### Clear the current workspace (DOESN'T CLEAR LOADED LIBRARIES)
rm(list=ls())

### install if necessary and then load the libraries you need

j <- c("rstudioapi","plyr","vegan","lme4","ggplot2","RColorBrewer","car","gridExtra","MuMIn","RVAideMemoire","emmeans","reshape2","moments","dplyr","ggmap","blighty","grid","egg","raster","MuMIn","nlme","maps","lmodel2")

new.packages <- j[!(j %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

lapply(j, require, character.only = TRUE)  # loads up any libraries that aren't already loaded


### set working directory to script's saved location

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# source a function
source("CheckResidsFunction.R")


## import the data

# first the individual-level moths we sampled, with their dry mass and forewing length
DF1 <- read.csv("Data/moth_data.csv", header=TRUE)
summary(DF1)


# and then the species-level data extracted from the field guides (min, max and median forewing length)
DF2 <-read.csv("Data/df2.csv", header=TRUE)
summary(DF2)


# try merging the two dataframes
DF.comb <- merge(DF1, DF2, by = 'BINOMIAL')    # this results in a frame that is 13 individuals shorter - which is the 13 unidentified micros
head(DF.comb)

DF.comb$SPECIES <- DF.comb$SPECIES.x
DF.comb$FAMILY <- DF.comb$FAMILY.x

DF.comb <- DF.comb[,-c(5,6,9,10)]



## before doing anything else, let's pick out training and testing datasets

# create frequency frame
frequencies <- ddply(DF.comb, .(SPECIES,BINOMIAL,FAMILY), summarise,
                     COUNT = length(DRY_MASS))

summary(frequencies)

write.csv(frequencies, file = "Data/Frequencies.csv", row.names = F)


# use those lists to pick out data
DF.comb.noct <- DF.comb[which(DF.comb$FAMILY=="Noctuidae"),]
DF.comb.good <- DF.comb[which(DF.comb$DRY_MASS>0),]

### a slight diversion to pull out data by site and species abundance

site.freq <- ddply (DF.comb, .(SPECIES,BINOMIAL,FAMILY,SITE), summarise,
                    COUNT = length(DRY_MASS))

summary(site.freq)

write.csv(site.freq, file = "Data/SiteFrequencies.csv", row.names = F)

## now we can move on to the actual analysis

# first we need to figure out which family of GLMM to use
hist(DF.comb.good$DRY_MASS)   # this looks neg binomial, but may be skewed by inclusion of some v. large moths of large families
hist(DF.comb.noct$DRY_MASS)   # still looks very skewed
hist(log(DF.comb.noct$DRY_MASS))  # fab - so our data are log-normal (that might have been expected!)



### now we can test objective 1 - what is the relationship between body mass, forewing length and family?

# fit the model, with both family and forewing length (i.e. different families can have different intercepts AND slopes)
mod1 <- lmer(log(DRY_MASS) ~ FOREWING_LENGTH * FAMILY + (1|SPECIES),
             data=DF.comb.good)

summary(mod1)
drop1(mod1, test = "Chi")
Anova(mod1, test = "Chisq")

chkres(mod1, DF.comb.good$FOREWING_LENGTH, DF.comb.good$FAMILY)  # these residuals are fine - actually pretty great!

r.squaredGLMM(mod1)

# try building an alternative, non-linear mixed effects model
hist(log(DF.comb.good$DRY_MASS))

mod1nl <- lmer(log(DRY_MASS) ~ log(FOREWING_LENGTH) * FAMILY + (1|SPECIES),
             data=DF.comb.good)

summary(mod1nl)
drop1(mod1nl, test = "Chi")
Anova(mod1nl, test = "Chisq")

chkres(mod1nl, log(DF.comb.good$FOREWING_LENGTH), DF.comb.good$FAMILY)

r.squaredGLMM(mod1nl)

## could also try a broken-stick model and see whether it is a better fit

# first set up a vector of positions that we allow the break to take up
# from visual inspection of the data it's likely to fall somewhere between 10 mm and 30 mm forewing length

breaks <- c(10:30)

mse <- numeric(length(breaks))

# loop over all possible break positions and pick out the mse of each

for (i in 1:length(breaks)){
  suppressMessages(piecewise1 <- lmer(log(DRY_MASS) ~ FAMILY * FOREWING_LENGTH*(FOREWING_LENGTH < breaks[i])
                                      + FAMILY * FOREWING_LENGTH*(FOREWING_LENGTH >= breaks[i])
                                      + (1|SPECIES),
                                      data = DF.comb.good))
  mse[i] <- summary(piecewise1)[6]
}

# pick out the model with the best mse
mse <- as.numeric(mse)
bestbreak <- breaks[which(mse==min(mse))]

# set up that model
suppressMessages(piecewise2 <- lmer(log(DRY_MASS) ~ FAMILY * FOREWING_LENGTH*(FOREWING_LENGTH < bestbreak)
                                    + FAMILY * FOREWING_LENGTH*(FOREWING_LENGTH >= bestbreak)
                                    + (1|SPECIES),
                                    data = DF.comb.good))

# pull out the BICs for comparison
BIC(piecewise2,mod1,mod1nl)
AIC(piecewise2,mod1,mod1nl)

# the lowest BIC is the best-fitting model after penalisation
# which is the non-linear model, by some distance! So let's use that one going forwards

### we also want to calculate within-groups and between-groups variance in measured body mass - 
# i.e. the variance between individuals of the same species, and the variance between species

# within-groups function

test.lme <- lme(DRY_MASS ~ 1, random = ~1|SPECIES, data = DF.comb.good)
VarCorr(test.lme)






### now start building our predictive model (towards objective 2)
# first pick out the families with at least 5 species included in the training dataset

fam.freqs <- ddply(frequencies, .(FAMILY), summarise,
                   NO.SPECIES = length(COUNT),
                   NO.INDIVS = sum(COUNT))
fam.freqs

# create a crudified family variable (grouping everything with low species richness under 'others')
families.to.keep <- levels(droplevels(fam.freqs[which(fam.freqs$NO.SPECIES >= 5),1]))

DF.comb.good$FAM.CRUDE <- as.factor(ifelse(DF.comb.good$FAMILY %in% families.to.keep,
                                           as.character(DF.comb.good$FAMILY),
                                           "Other"))


# make 'other' the intercept level
DF.comb.good$FAM.CRUDE <- relevel(DF.comb.good$FAM.CRUDE, "Other")
summary(DF.comb.good$FAM.CRUDE)

### now we can train the model!
mod2 <- lmer(log(DRY_MASS) ~ log(FOREWING_LENGTH) * FAM.CRUDE + (1|SPECIES),
             data=DF.comb.good)

summary(mod2)
drop1(mod2, test = "Chi")

chkres(mod2, log(DF.comb.good$FOREWING_LENGTH), DF.comb.good$FAM.CRUDE)  # again these are fine

r.squaredGLMM(mod2)
VarCorr(mod2)

# for the testing phase, we might also need a generic model (relationship only, no family adjustment at all!)
mod3 <- lmer(log(DRY_MASS) ~ log(FOREWING_LENGTH) + (1|SPECIES),
             data = DF.comb.good)

summary(mod3)
drop1(mod3, test = "Chi")

chkres(mod3, log(DF.comb.good$FOREWING_LENGTH), DF.comb.good$FAM.CRUDE)

r.squaredGLMM(mod3)

BIC(mod2,mod3)


# predict from these models
# family-adjustment
summary(mod2)


newdat1 <- DF2
newdat1$FAM.CRUDE <- as.factor(ifelse(newdat1$FAMILY %in% families.to.keep, as.character(newdat1$FAMILY), "Other"))
newdat1$FAM.CRUDE <- relevel(newdat1$FAM.CRUDE, "Other")
newdat1$FOREWING_LENGTH <- newdat1$Med_fw
newdat1$DRY_MASS <- 1


mm1 <- model.matrix(terms(mod2),newdat1)
newdat1$log_PRED_DRY_MASS = mm1 %*% fixef(mod2)
pvar1 <- diag(mm1 %*% tcrossprod(vcov(mod2),mm1))
newdat1$SE <- sqrt(pvar1)

newdat1$PRED_DRY_MASS <- exp(newdat1$log_PRED_DRY_MASS)

DF2 <- merge(DF2, newdat1)

# and generic
summary(mod3)

newdat2 <- DF2
newdat2$FOREWING_LENGTH <- newdat2$Med_fw
newdat2$DRY_MASS <- 1

mm2 <- model.matrix(terms(mod3),newdat2)
newdat2$log_PRED_DRY_MASS_GEN = mm2 %*% fixef(mod3)
pvar2 <- diag(mm2 %*% tcrossprod(vcov(mod3),mm2))
newdat2$SE_GEN <- sqrt(pvar2)

newdat2$PRED_DRY_MASS_GEN <- exp(newdat2$log_PRED_DRY_MASS_GEN)

DF2 <- merge(DF2, newdat2)


## 
### before moving on any further....
### plot a figure showing forewing length:body mass:family relationship

#p= graph using individuals as points
p1 <- ggplot(DF.comb.good, aes(x=FOREWING_LENGTH, y=DRY_MASS, colour=FAMILY, shape=FAMILY))+
  geom_point() +
  xlab("Forewing length (mm)")+
  ylab("Dry mass (mg)")+
  theme_light() +
  theme(panel.grid = element_blank())+
  scale_shape_manual(name = "Family", values = c(15,15,16,16,17,17,15,15,16,16,17))+
  scale_colour_manual(name = "Family", values = c("tan1","steelblue1", "palevioletred3", "skyblue4", "firebrick3", "seagreen3", "darkorchid3", "darkslategray", "orangered", "deepskyblue", "maroon3"))+
  scale_y_log10(limits = c(1,1000))+
  scale_x_log10(limits = c(5,50))+
  ggtitle("\nIndividuals")+
  theme(plot.title = element_text(hjust = 0.5))

p1 




p2 <- ggplot(DF.comb.good, aes(x=FOREWING_LENGTH, y=DRY_MASS, colour=FAM.CRUDE, shape=FAM.CRUDE))+
  geom_point()+
  xlab("Forewing length (mm)")+
  ylab("Dry mass (mg)")+
  theme_light()+
  theme(panel.grid = element_blank())+
  scale_shape_manual(name = "Family", values = c(15,15,15,16,17))+
  scale_colour_manual(name = "Family", values = c("grey80", "tan1", "steelblue1", "palevioletred3", "seagreen3"))+
  scale_y_log10(limits = c(1,1000))+
  scale_x_log10(limits = c(5,50))+
  ggtitle("\nIndividuals")+
  theme(plot.title = element_text(hjust = 0.5))


p2


## also want to do this with species averages as points (with error bars)

species <- ddply(DF1, .(SPECIES, FAMILY), summarise,
                 COUNT = length(FOREWING_LENGTH),
                 Mean.fw =mean(FOREWING_LENGTH),
                 Mean.bm =mean(DRY_MASS),
                 se.fw = sd(FOREWING_LENGTH)/sqrt(length(FOREWING_LENGTH)),
                 se.bm = sd(DRY_MASS)/sqrt(length(DRY_MASS)))

species$se.fw <- ifelse(species$se.fw == "Inf", 0, species$se.fw)

species$lci.fw <- species$Mean.fw - (1.96 * species$se.fw)

species$uci.fw <- species$Mean.fw + (1.96 * species$se.fw)

species$lci.bm <- species$Mean.bm - (1.96 * species$se.bm)

species$lci.bm <- ifelse(species$lci.bm < 0, 1, species$lci.bm)

species$uci.bm <- species$Mean.bm + (1.96 * species$se.bm)

species <- species[complete.cases(species), ]

crude <- merge(species, DF.comb.good)
str(crude)

#q= graph using species as points

q1 <- ggplot(species, aes(x=Mean.fw, y=Mean.bm, colour=FAMILY, shape=FAMILY, fill=FAMILY))+
  geom_point(size = 2.5)+
  geom_errorbar(aes(ymin=lci.bm, ymax=uci.bm))+
  geom_errorbarh(aes(xmin=lci.fw, xmax=uci.fw))+
  xlab("Mean forewing length (mm)")+
  ylab("Mean dry mass (mg)")+
  theme_light()  +
  theme(panel.grid = element_blank())+
  scale_shape_manual(name = "Family", values = c(15,15,16,16,25,24,15,15,16,16,15))+
  scale_colour_manual(name = "Family", values = c("tan1","steelblue1", "palevioletred3", "skyblue4", "firebrick3", "seagreen3", "darkorchid3", "darkslategray", "orangered", "deepskyblue", "maroon3")) +
  scale_fill_manual(name = "Family", values = c("tan1","steelblue1", "palevioletred3", "skyblue4", "firebrick3", "seagreen3", "darkorchid3", "darkslategray", "orangered", "deepskyblue", "maroon3")) +
  scale_y_log10(limits = c(0.5,1000))+
  scale_x_log10(limits = c(5,50))+
  ggtitle("\nSpecies means")+
  theme(plot.title = element_text(hjust = 0.5))


q1                                   


q2 <- ggplot(crude, aes(x=Mean.fw, y=Mean.bm, colour=FAM.CRUDE, shape=FAM.CRUDE))+
  geom_point(size = 2.5)+
  geom_errorbar(aes(ymin=lci.bm, ymax=uci.bm))+
  geom_errorbarh(aes(xmin=lci.fw, xmax=uci.fw))+
  xlab("Mean forewing length(mm)")+
  ylab("log(Mean dry mass (mg)")+
  theme_light()+ 
  theme(panel.grid = element_blank())+
  scale_y_log10(limits = c(0.5,1000))+
  scale_x_log10(limits = c(5,50))+
  scale_shape_manual(name = "Family", values = c(15,15,15,16,17))+
  scale_colour_manual(name = "Family", values = c("grey80", "tan1", "steelblue1", "palevioletred3", "seagreen3"))+
  ggtitle("\nSpecies means")+
  theme(plot.title = element_text(hjust = 0.5))

q2 


m1 <- grid.arrange(q1, p2, ncol = 1)


ggsave("Plots/Fig1.svg", m1, device = svg, width = 16, height = 24, units = "cm")


# and thirdly plot a figure showing the relationship between the expected fw lengths from the book and the actual measured forewing lengths

fwmod <- lm(data=DF.comb.good, Med_fw ~ FOREWING_LENGTH)
summary(fwmod)
summary(fwmod)$coefficients
drop1(fwmod, test = "F")

r.squaredGLMM(fwmod)

# actually a type II regression might be more theoretically appropriate here, so let's try one:

fwmodII <- lmodel2(FOREWING_LENGTH ~ Med_fw, data = DF.comb.good, nperm = 100)

fwmodII



fw <- ggplot(DF.comb.good, aes(x=Med_fw, y=FOREWING_LENGTH))+
  geom_point()+
  geom_abline(intercept = fwmodII$regression.results[2,2], slope = fwmodII$regression.results[2,3])+
  geom_abline(intercept = 0, slope = 1, colour = "royalblue")+
  xlim(0, 45)+ylim(0, 40)+
  xlab("Expected forewing length")+
  ylab("Measured forewing length")+
  theme_light()+
  theme(panel.grid = element_blank())
  
fw

ggsave("Plots/FigS3.svg", fw, device = svg, width = 12, height = 12, units = "cm")


# repeat this at species level

forewings.test <- ddply(DF.comb.good, .(BINOMIAL), summarise,
                        MEASURED = mean(FOREWING_LENGTH),
                        PREDICTED = mean(Med_fw))


# and thirdly plot a figure showing the relationship between the expected fw lengths from the book and the actual measured forewing lengths

fwmod.spec <- lm(data=forewings.test, PREDICTED ~ MEASURED)
summary(fwmod.spec)
summary(fwmod.spec)$coefficients
drop1(fwmod.spec, test = "F")

r.squaredGLMM(fwmod.spec)

# actually a type II regression might be more theoretically appropriate here, so let's try one:

fwmodII.spec <- lmodel2(MEASURED ~ PREDICTED, data = forewings.test, nperm = 100)

fwmodII.spec



fw.spec <- ggplot(forewings.test, aes(x=PREDICTED, y=MEASURED))+
  geom_point()+
  geom_abline(intercept = fwmodII.spec$regression.results[2,2], slope = fwmodII.spec$regression.results[2,3])+
  geom_abline(intercept = 0, slope = 1, colour = "royalblue")+
  xlim(0, 45)+ylim(0, 40)+
  xlab("Expected forewing length")+
  ylab("Measured forewing length")+
  theme_light()+
  theme(panel.grid = element_blank())

fw.spec

ggsave("Plots/FigS3a.svg", fw.spec, device = svg, width = 12, height = 12, units = "cm")





# also checking the percentage of measured forewing lengths that are within the range of min/max
DF.comb.good$FW_test <- ifelse(DF.comb.good$FOREWING_LENGTH >= DF.comb.good$Min_fw & DF.comb.good$FOREWING_LENGTH <= DF.comb.good$Max_fw, 1, 0)

summary(DF.comb.good$FW_test)


# look at some histograms of fw lengths for the four most abundant species (all noctuids), plus the most abundant erebid and lasiocampid
# - are they normally-distributed?

# list of species to test - these are the 15 most-abundant species in the dataset
normal.to.test <- c("Middle-barred Minor","Uncertain","Dark Arches","Double Square-spot","Heart & Dart",
                    "Smoky Wainscot","Common Footman","Common Rustic","Dun-bar","Large Yellow Underwing",
                    "Drinker","Marbled Minor","Dingy Footman","Agriphila straminella","Ingrailed Clay")

normal.test <- data.frame(FAMILY = factor(),
                          SPECIES = factor(),
                          W = numeric(),
                          P = numeric(),
                          SKEWNESS = numeric())


for (x in normal.to.test){
  spec <- DF.comb.good[which(DF.comb.good$SPECIES == x), ]
  
  FAMILY <- as.character(spec[1,'FAMILY'])
  SPECIES <- as.character(spec[1,'SPECIES'])
  
  W <- shapiro.test(spec$FOREWING_LENGTH)[[1]]
  P <- shapiro.test(spec$FOREWING_LENGTH)[[2]]
  
  SKEWNESS <- skewness(spec$FOREWING_LENGTH)
  
  out <- cbind(FAMILY,SPECIES,W,P,SKEWNESS)
  
  normal.test <- rbind(normal.test, out)
  
}

normal.test$P <- as.numeric(as.character(normal.test$P))

normal.test$P.adj <- p.adjust(normal.test$P, method = "BH")


## they are not normally distributed in most cases
# - which might explain why the midpoint of the range of forewing lengths tends to be slightly lower than the average measured length,
# except that skewness goes off in both directions!



### now do some testing - how well do sample-level biomass estimates (as would be predicted from the model)
# correlate with measured sample biomasses?

# first get we have to estimate the body mass of every individual... 

#### extract universal body masses ####

# as an added feature to the paper, we want to extract the estimated body mass of all UK macro-moths 
# plus the Crambidae (because they are a refined-estimate family) and the Pyralidae (because why not),
# not just those we need for this piece of analysis

# we've prepared all their forewing lengths in a .csv
forewings <- read.csv("../Agassiz/Species_macro.csv", header = T)

# rename a column to match up
forewings$BINOMIAL <- forewings$RIS_BINOMIAL

# first extract the median forewing length for each species
forewings$FOREWING_MED <- (forewings$FOREWING_LB + forewings$FOREWING_UB)/2


# now copy the formula down from above to estimate the LB, MED and UB of body mass for every species
summary(mod2)


newdat <- forewings
newdat$FAM.CRUDE <- as.factor(ifelse(newdat$FAMILY %in% families.to.keep, as.character(newdat$FAMILY), "Other"))
newdat$FAM.CRUDE <- relevel(newdat$FAM.CRUDE, "Other")
newdat$FOREWING_LENGTH <- newdat$FOREWING_MED
newdat$DRY_MASS <- 1

summary(newdat)

mm <- model.matrix(terms(mod2),newdat)
newdat$log_DRY_MASS = mm %*% fixef(mod2)
pvar <- diag(mm %*% tcrossprod(vcov(mod2),mm))
newdat$SE <- sqrt(pvar)

newdat$PRED_DRY_MASS <- exp(newdat$log_DRY_MASS)




# and now export this data out, first jiggling the order of columns

newdat <- newdat[,c(1:5,11,7:8,12,9,18,16:17,10)]

newdat_out <- newdat[,-6]

newdat_out <- newdat_out[which(newdat_out$EXCLUDE_FROM_FINAL == "N"), ]
newdat_out <- newdat_out[,-13]

write.csv(newdat_out, "Data/TableS1.csv", row.names = F)


# plot a figure with all species in 

newdat_out$BODY_MASS_LB <- newdat_out$PRED_DRY_MASS - (1.96 * newdat_out$SE)
newdat_out$BODY_MASS_UB <- newdat_out$PRED_DRY_MASS + (1.96 * newdat_out$SE)

p3 <- ggplot(newdat_out, aes(x=FOREWING_MED, y=PRED_DRY_MASS, colour=FAMILY, shape=FAMILY, fill=FAMILY))+
  geom_point(size = 2.5)+
  geom_errorbar(aes(ymin=BODY_MASS_LB, ymax=BODY_MASS_UB))+
  geom_errorbarh(aes(xmin=FOREWING_LB, xmax=FOREWING_UB))+
  xlim(0, 60) + ylim(0, 1100) +
  xlab("Forewing length (mm)")+
  ylab("Estimated dry mass (mg)")+
  theme_light()  +
  theme(panel.grid = element_blank())+
  scale_shape_manual(name = "Family", values = c(15,17,15,25,15,15,16,16,17,25,17,25,15,25,15,16,16))+
  scale_colour_manual(name = "Family", values = c("black","maroon3","tan1","darkgreen","tan1","steelblue1", "palevioletred3", "skyblue4", "firebrick3", "orangered4", "seagreen3", "grey70", "darkorchid3", "darkslategray", "gold2", "orangered", "deepskyblue")) +
  scale_fill_manual(name = "Family", values = c("black","maroon3","tan1","darkgreen","tan1","steelblue1", "palevioletred3", "skyblue4", "firebrick3", "orangered4", "seagreen3", "grey70", "darkorchid3", "darkslategray", "gold2", "orangered", "deepskyblue")) +
  scale_y_log10()+
  ggtitle("Species means")+
  theme(plot.title = element_text(hjust = 0.5))


p3     


# it looks quite messy so let's compare to a figure without error bars

p3a <- ggplot(newdat_out, aes(x=FOREWING_MED, y=PRED_DRY_MASS, colour=FAMILY, shape=FAMILY, fill=FAMILY))+
  geom_point(size = 2.5)+
  xlim(0, 60) + ylim(0, 1100) +
  xlab("Forewing length (mm)")+
  ylab("Estimated dry mass (mg)")+
  theme_light()  +
  theme(panel.grid = element_blank())+
  scale_shape_manual(name = "Family", values = c(15,15,17,25,15,16,16,25,17,17,25,15,15,25,16,16,15))+
  scale_colour_manual(name = "Family", values = c("grey80","tan1","darkgreen","darkblue","steelblue1", "palevioletred3", "skyblue4", "firebrick3", "orangered4", "seagreen3", "grey70", "darkorchid3", "darkslategray", "gold2", "deepskyblue", "orangered","maroon3")) +
  scale_fill_manual(name = "Family", values = c("grey80","tan1","darkgreen","darkblue","steelblue1", "palevioletred3", "skyblue4", "firebrick3", "orangered4", "seagreen3", "grey70", "darkorchid3", "darkslategray", "gold2", "deepskyblue","orangered","maroon3")) +
  scale_y_log10()+
  scale_x_log10()+
  ggtitle("Species means")+
  theme(plot.title = element_text(hjust = 0.5))


p3a
q1

ggsave("Plots/FigS6.svg", p3a, device = svg, width = 20, height = 15, units = "cm")






## now return to the analysis

# get rid of subspecies (or rather, smooth them out at species-level) as they confuse the next step
newdat_on <- ddply(newdat, .(FAMILY,BINOMIAL), numcolwise(mean))



## merge these predicted biomasses into the main dataset
DF1.pred <- merge(DF1,newdat_on, by = c('BINOMIAL','FAMILY'))


# this disposes of the 13 unidentified individuals, which is fine (for comparison's sake) because we can't predict their mass

## sum up measured and predicted moth biomasses at sample level

# summarise the data down into samples first (1 site on 1 night, which might include 1 or 2 traps' catch) 
# with a total count and a summed mass of individuals per species per sample


daily <- ddply(DF1.pred, .(DATE, SPECIES, SITE), summarise,
               COUNT = length(DRY_MASS),
               BIOMASS = sum(DRY_MASS),
               PRED_BIOMASS = sum(PRED_DRY_MASS))


# and now generate sample-level summary stats, including a range of diversity measures using vegan

srichnew <- ddply(daily, .(DATE, SITE), summarise,
                  Abundance = sum(COUNT),
                  Richness = length(COUNT),
                  Biomass = sum(BIOMASS),
                  Pred_Biomass = sum(PRED_BIOMASS),
                  Shannon = diversity(COUNT, index = "shannon"),
                  Simpson = diversity(COUNT, index = "simpson"),
                  invSimpson = diversity(COUNT, index = "invsimpson"),
                  Evenness = diversity(COUNT, index = "shannon")/log(specnumber(COUNT)))

# make site a factor and date behave as a date (rather than a factor)
srichnew$fSITE <- as.factor(srichnew$SITE)                  
srichnew$DATE <- as.Date(srichnew$DATE, format="%d/%m/%Y")



# now we can start analysing relationships between these variables -
# first, biomass vs pred_biomass
hist(srichnew$Pred_Biomass)
hist(srichnew$Biomass)

hist(log(srichnew$Pred_Biomass))
hist(log(srichnew$Biomass))

mod4 <- lmer(log(Biomass) ~ log(Pred_Biomass) + (1|fSITE),
             data=srichnew)

summary(mod4)
drop1(mod4, test="Chi")

chkres(mod4, srichnew$Pred_Biomass, srichnew$fSITE)

r.squaredGLMM(mod4)

## do this as a model II regression
mod4II <- lmodel2(log(Pred_Biomass) ~ log(Biomass),
             data=srichnew, nperm = 100)

mod4II



## now let's plot out a figure showing predicted vs measured biomass at a sample level
# and also look a bit at the scale and direction of errors

### examine disparity (direction and spread of error) both for all estimates and for misses only

# first check whether estimates tend to be larger or smaller than measured body mass (i.e. is there a systematic error?)

srichnew$ERROR_DIRECTION <- as.factor(ifelse(srichnew$Pred_Biomass > srichnew$Biomass, "Overestimate","Underestimate"))

summary(srichnew$ERROR_DIRECTION)

# and by how much?
srichnew$DIFF <- srichnew$Pred_Biomass - srichnew$Biomass

# overestimates...
srichnew.over <- srichnew[which(srichnew$ERROR_DIRECTION == "Overestimate"), ]
summary(srichnew.over$DIFF)
sd(abs(srichnew.over$DIFF))/sqrt(length(abs(srichnew.over$DIFF)))

# underestimates...
srichnew.under <- srichnew[which(srichnew$ERROR_DIRECTION == "Underestimate"), ]
summary(srichnew.under$DIFF)
sd(abs(srichnew.under$DIFF))/sqrt(length(abs(srichnew.under$DIFF)))


# more overestimates but larger underestimates BUT this is a product of the logging

# combined effect
summary(srichnew$DIFF)
sd(abs(srichnew$DIFF))/sqrt(length(abs(srichnew$DIFF)))

t.test(srichnew$DIFF)




### plot
# first, order all misses from biggest underestimate to biggest overestimate
error.sort <- srichnew[order(srichnew$DIFF), ]
error.sort$ORDER <- c(nrow(error.sort):1)



errors.plot <- ggplot(error.sort, aes(x = ORDER, y = abs(DIFF), colour = ERROR_DIRECTION))+
  geom_point()+
  scale_colour_manual(name = "Direction of error", values = c("goldenrod","royalblue"),labels = c("Overestimate","Underestimate"))+

  theme_light()+
  theme(panel.grid = element_blank())+
  xlab("Samples of moths ordered by margin and direction of error in biomass estimation")+
  ylab("Margin of error (mg)")+
  theme(legend.position = "bottom")


errors.plot


# now plot biomass vs predicted directly, colouring by sample size

pred.sample <- ggplot(srichnew, aes(x=Biomass,y=Pred_Biomass,fill=Abundance))+
  geom_abline(intercept = 0, slope = 1, colour = "royalblue")+
  geom_point(shape = 21, cex = 2)+
  scale_fill_distiller(palette = "YlOrRd", direction = 2)+
  theme_light()+
  theme(panel.grid = element_blank())+
  scale_x_log10(limits = c(10,2000))+
  scale_y_log10(limits = c(10,2000))+
  xlab("Measured dry mass (mg) of samples of moths\n")+
  ylab("Predicted dry mass\n (mg) of samples of moths")


pred.sample


# and now plot what's effectively the residuals of the model against sample abundance - 

resids.sample <- ggplot(srichnew, aes(x=Abundance, y = DIFF, fill = Abundance))+
  geom_hline(aes(yintercept = 0), size = 1, colour = "grey50")+
  geom_point(shape = 21, cex = 2)+
  scale_fill_distiller(palette = "YlOrRd", direction = 2)+
  theme_light()+
  theme(panel.grid = element_blank())+
  xlab("Abundance of moths in samples\n")+
  ylab("Residual of predicted dry mass\n (mg) of samples of moths")

resids.sample



# function to extract a legend
g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}


leg.misses <- g_legend(pred.sample)

m.misses <- grid.arrange(arrangeGrob(pred.sample + theme(legend.position="none"),
                                     resids.sample + theme(legend.position="none"),
                                     ncol = 1),
                         leg.misses, ncol = 2, widths = c(10,1))



ggsave("Plots/FigS4.svg", plot = m.misses, device = svg, width = 160, height = 240, units = "mm", limitsize = T)





## do some of the same plots and tests at the level of individual moths

# first, biomass vs pred_biomass
hist(DF1.pred$PRED_DRY_MASS)
hist(DF1.pred$DRY_MASS)

hist(log(DF1.pred$PRED_DRY_MASS))
hist(log(DF1.pred$DRY_MASS))

# we'll need to drop out the Narycia duplicella record (with 0 dry mass) to allow us to log
DF1.pred.nd <- DF1.pred[which(DF1.pred$BINOMIAL != "Narycia duplicella"), ]


mod4.in <- lmer(log(DRY_MASS) ~ log(PRED_DRY_MASS) + (1|SPECIES),
             data=DF1.pred.nd)

summary(mod4.in)
drop1(mod4.in, test="Chi")

chkres(mod4.in, DF1.pred.nd$PRED_DRY_MASS, DF1.pred.nd$SPECIES)

r.squaredGLMM(mod4.in)


## now let's plot out a figure showing predicted vs measured biomass at a sample level
# and also look a bit at the scale and direction of errors

### examine disparity (direction and spread of error) both for all estimates and for misses only

# first check whether estimates tend to be larger or smaller than measured body mass (i.e. is there a systematic error?)

DF1.pred.nd$ERROR_DIRECTION <- as.factor(ifelse(DF1.pred.nd$PRED_DRY_MASS > DF1.pred.nd$DRY_MASS, "Overestimate","Underestimate"))

summary(DF1.pred.nd$ERROR_DIRECTION)

# and by how much?
DF1.pred.nd$DIFF <- DF1.pred.nd$PRED_DRY_MASS - DF1.pred.nd$DRY_MASS

# overestimates...
DF1.pred.nd.over <- DF1.pred.nd[which(DF1.pred.nd$ERROR_DIRECTION == "Overestimate"), ]
summary(DF1.pred.nd.over$DIFF)
sd(abs(DF1.pred.nd.over$DIFF))/sqrt(length(abs(DF1.pred.nd.over$DIFF)))

# underestimates...
DF1.pred.nd.under <- DF1.pred.nd[which(DF1.pred.nd$ERROR_DIRECTION == "Underestimate"), ]
summary(DF1.pred.nd.under$DIFF)
sd(abs(DF1.pred.nd.under$DIFF))/sqrt(length(abs(DF1.pred.nd.under$DIFF)))


# more overestimates but larger underestimates BUT this is a product of the logging

# combined effect
summary(DF1.pred.nd$DIFF)
sd(abs(DF1.pred.nd$DIFF))/sqrt(length(abs(DF1.pred.nd$DIFF)))

t.test(DF1.pred.nd$DIFF)




### plot
# first, order all misses from biggest underestimate to biggest overestimate
error.sort.in <- DF1.pred.nd[order(DF1.pred.nd$DIFF), ]
error.sort.in$ORDER <- c(nrow(error.sort.in):1)



errors.plot.in <- ggplot(error.sort.in, aes(x = ORDER, y = abs(DIFF), colour = ERROR_DIRECTION))+
  geom_point()+
  scale_colour_manual(name = "Direction of error", values = c("goldenrod","royalblue"),labels = c("Overestimate","Underestimate"))+
  
  theme_light()+
  theme(panel.grid = element_blank())+
  xlab("Individual moths ordered by margin and direction of error in biomass estimation")+
  ylab("Margin of error (mg)")+
  theme(legend.position = "bottom")


errors.plot.in


# now plot biomass vs predicted directly, again colouring by direction of error

pred.in <- ggplot(DF1.pred.nd, aes(x=DRY_MASS,y=PRED_DRY_MASS))+
  geom_point()+
  geom_abline(intercept = 0, slope = 1, colour = "royalblue")+
  theme_light()+
  theme(panel.grid = element_blank())+
  scale_y_log10(limits = c(1,1000))+
  scale_x_log10(limits = c(1,1000))+
  xlab("Measured dry mass (mg) of samples of moths")+
  ylab("Predicted dry mass (mg) of samples of moths")


pred.in


m.misses.in <- grid.arrange(arrangeGrob(errors.plot + theme(legend.position="none"),
                                     errors.plot.in + theme(legend.position="none"),
                                     ncol = 2),
                         leg.misses, ncol = 2, widths=c(10,2))



ggsave("Plots/FigS5.svg", plot = m.misses.in, device = svg, width = 300, height = 200, units = "mm", limitsize = T)




## and do some of the same plots and tests at the level of moth species 

# starting with DF1.pred.nd, take the mean (and s.e.) of measured body mass for each species

spec.pred <- ddply(DF1.pred.nd, .(BINOMIAL,FAMILY,SPECIES), summarise,
                   MEAN_DRY_MASS = mean(DRY_MASS),
                   SD_DRY_MASS = sd(DRY_MASS),
                   SE_DRY_MASS = sd(DRY_MASS)/sqrt(length(DRY_MASS)),
                   PRED_DRY_MASS = mean(PRED_DRY_MASS),
                   SE_PRED = mean(SE),
                   N = length(DRY_MASS))

summary(spec.pred)

spec.pred$SD_DRY_MASS <- ifelse(is.na(spec.pred$SD_DRY_MASS), 0, spec.pred$SD_DRY_MASS)
spec.pred$SE_DRY_MASS <- ifelse(is.na(spec.pred$SE_DRY_MASS), 0, spec.pred$SE_DRY_MASS)


mod4.sp <- lmer(log(PRED_DRY_MASS) ~ log(MEAN_DRY_MASS) + (1|FAMILY),
                data=spec.pred)

summary(mod4.sp)
drop1(mod4.sp, test="Chi")

chkres(mod4.sp, spec.pred$PRED_DRY_MASS, spec.pred$FAMILY)

r.squaredGLMM(mod4.sp)


# do this as a model II regression
mod4.spII <- lmodel2(log(PRED_DRY_MASS) ~ log(MEAN_DRY_MASS),
                  data=spec.pred, nperm = 100)

mod4.spII


# try redoing this regression with the 10 smallest species removed
spec.pred.res <- spec.pred[which(spec.pred$MEAN_DRY_MASS >= 15), ]

mod4.spII.res <- lmodel2(log(PRED_DRY_MASS) ~ log(MEAN_DRY_MASS),
                     data=spec.pred.res, nperm = 100)

mod4.spII.res



## now let's plot out a figure showing predicted vs measured biomass at a sample level
# and also look a bit at the scale and direction of errors

### examine disparity (direction and spread of error) both for all estimates and for misses only

# first check whether estimates tend to be larger or smaller than measured body mass (i.e. is there a systematic error?)

spec.pred$ERROR_DIRECTION <- as.factor(ifelse(spec.pred$PRED_DRY_MASS > spec.pred$MEAN_DRY_MASS, "Overestimate","Underestimate"))

summary(spec.pred$ERROR_DIRECTION)

# and by how much?
spec.pred$DIFF <- spec.pred$PRED_DRY_MASS - spec.pred$MEAN_DRY_MASS

# overestimates...
spec.pred.over <- spec.pred[which(spec.pred$ERROR_DIRECTION == "Overestimate"), ]
summary(spec.pred.over$DIFF)
sd(abs(spec.pred.over$DIFF))/sqrt(length(abs(spec.pred.over$DIFF)))

# underestimates...
spec.pred.under <- spec.pred[which(spec.pred$ERROR_DIRECTION == "Underestimate"), ]
summary(spec.pred.under$DIFF)
sd(abs(spec.pred.under$DIFF))/sqrt(length(abs(spec.pred.under$DIFF)))


# more overestimates but larger underestimates BUT this is a product of the logging

# combined effect
summary(spec.pred$DIFF)
sd(abs(spec.pred$DIFF))/sqrt(length(abs(spec.pred$DIFF)))

t.test(spec.pred$DIFF)




### plot
# first, order all misses from biggest underestimate to biggest overestimate
error.sort.sp <- spec.pred[order(spec.pred$DIFF), ]
error.sort.sp$ORDER <- c(nrow(error.sort.sp):1)



errors.plot.sp <- ggplot(error.sort.sp, aes(x = ORDER, y = abs(DIFF), colour = ERROR_DIRECTION))+
  geom_point()+
  scale_colour_manual(name = "Direction of error", values = c("goldenrod","royalblue"),labels = c("Overestimate","Underestimate"))+
  
  theme_light()+
  theme(panel.grid = element_blank())+
  xlab("Moth species ordered by margin and direction of error in biomass estimation")+
  ylab("Margin of error (mg)")+
  theme(legend.position = "bottom")


errors.plot.sp


# now plot biomass vs predicted directly

pred.sp <- ggplot(spec.pred, aes(x=MEAN_DRY_MASS,y=PRED_DRY_MASS,fill=N))+
  geom_abline(intercept = 0, slope = 1, colour = "royalblue")+
  geom_point(shape = 21, cex = 2)+
  scale_fill_distiller(palette = "YlOrRd", direction = 2, name = "Abundance")+
  theme_light()+
  theme(panel.grid = element_blank())+
  scale_y_log10(limits = c(1,1000))+
  scale_x_log10(limits = c(1,1000))+
  xlab("Mean measured dry mass (mg) of moth species\n")+
  ylab("Predicted dry mass\n (mg) of moth species")


pred.sp


# and now plot what's effectively the residuals of the model against sample abundance - 

resids.sp <- ggplot(spec.pred, aes(x=N, y = DIFF, fill = N))+
  geom_hline(aes(yintercept = 0), size = 1, colour = "grey50")+
  geom_point(shape = 21, cex = 2)+
  scale_fill_distiller(palette = "YlOrRd", direction = 2)+
  theme_light()+
  theme(panel.grid = element_blank())+
 xlab("Sampled individuals of species\n")+
  ylab("Residual of predicted dry mass\n (mg) of moth species")

resids.sp




# pull together a few of these for a final figure

leg.misses.sp <- g_legend(pred.sp)

m.misses.sp <- grid.arrange(arrangeGrob(pred.sp + theme(legend.position="none"),
                                     resids.sp + theme(legend.position="none"),
                                     ncol = 1),
                            leg.misses.sp, ncol = 2, widths = c(10,1))



ggsave("Plots/FigS4a.svg", plot = m.misses.sp, device = svg, width = 160, height = 240, units = "mm", limitsize = T)



m.misses.all <- grid.arrange(arrangeGrob(pred.sp + theme(legend.position="none"),
                                         resids.sp + theme(legend.position="none"),
                                         pred.sample + theme(legend.position="none"),
                                         resids.sample + theme(legend.position="none"),
                                         ncol = 2),
                             leg.misses.sp, ncol = 2, widths = c(10,1.5))


ggsave("Plots/Fig2.svg", plot = m.misses.all, device = svg, width = 240, height = 160, units = "mm", limitsize = T)




### let's now turn to objective 3 - looking at the relationship between biomass and other variables

# first, biomass vs abundance

hist(srichnew$Biomass)
hist(log(srichnew$Biomass))


mod5 <- lmer(Biomass ~ Abundance + (1|fSITE),
             data=srichnew)


summary(mod5)
drop1(mod5, test="Chi")

chkres(mod5, srichnew$Abundance, srichnew$fSITE)  # model residuals are fine

# extract model R^2
r.squaredGLMM(mod5)

# retest using model II
mod5II <- lmodel2(Biomass ~ Abundance,
                  data = srichnew, nperm = 100)

mod5II


## plot it out

mod5plot <- ggplot(srichnew, aes(x=Abundance, y=Biomass))+
  geom_point()+
  geom_abline(intercept = summary(mod5)$coefficients[1,1], slope = summary(mod5)$coefficients[2,1])+
  xlim(0, 40) + ylim(20, 2000)+
  xlab("Abundance\n")+
  ylab("Biomass (mg)")+
  theme_light()+
  theme(panel.grid = element_blank())+
  ggtitle("Sampled moths")+
  theme(plot.title = element_text(hjust = 0.5))
mod5plot


# with predicted biomass

hist(srichnew$Pred_Biomass)
hist(log(srichnew$Pred_Biomass))


mod5p <- lmer(Pred_Biomass ~ Abundance + (1|fSITE),
             data=srichnew)


summary(mod5p)
drop1(mod5p, test="Chi")

chkres(mod5p, srichnew$Abundance, srichnew$fSITE)  # model residuals are fine

# extract model R^2
r.squaredGLMM(mod5p)

# retest using model II
mod5pII <- lmodel2(Pred_Biomass ~ Abundance,
                  data = srichnew, nperm = 100)

mod5pII


## plot it out

mod5plot.p <- ggplot(srichnew, aes(x=Abundance, y=Pred_Biomass))+
  geom_point()+
  geom_abline(intercept = summary(mod5p)$coefficients[1,1], slope = summary(mod5p)$coefficients[2,1])+
  xlim(0, 40) + ylim(20, 2000)+
  xlab("Abundance\n")+
  ylab("Predicted biomass (mg)")+
  theme_light()+
  theme(panel.grid = element_blank())+
  ggtitle("Sampled moths")+
  theme(plot.title = element_text(hjust = 0.5))
mod5plot.p




# next, biomass vs richness

mod6 <- lmer(Biomass ~ Richness + (1|fSITE),
             data=srichnew)

summary(mod6)
drop1(mod6, test="Chi")

chkres(mod6, srichnew$Abundance, srichnew$fSITE)

r.squaredGLMM(mod6)

# retest using model II
mod6II <- lmodel2(Biomass ~ Richness,
                   data = srichnew, nperm = 100)

mod6II


## plot it

mod6plot <- ggplot(srichnew, aes(x=Richness, y=Biomass))+
  geom_point()+
  geom_abline(intercept = summary(mod6)$coefficients[1,1], slope = summary(mod6)$coefficients[2,1])+
  xlim(0, 20)+ylim(20,2000)+
  xlab("Species richness\n")+
  ylab("Biomass(mg)")+
  theme_light()+
  theme(panel.grid = element_blank())
  

mod6plot


# predicted

mod6p <- lmer(Pred_Biomass ~ Richness + (1|fSITE),
             data=srichnew)

summary(mod6p)
drop1(mod6p, test="Chi")

chkres(mod6p, srichnew$Abundance, srichnew$fSITE)

r.squaredGLMM(mod6p)

# retest using model II
mod6pII <- lmodel2(Pred_Biomass ~ Richness,
                  data = srichnew, nperm = 100)

mod6pII


## plot it

mod6plot.p <- ggplot(srichnew, aes(x=Richness, y=Pred_Biomass))+
  geom_point()+
  geom_abline(intercept = summary(mod6p)$coefficients[1,1], slope = summary(mod6p)$coefficients[2,1])+
  xlim(0, 20)+ylim(20,2000)+
  xlab("Species richness\n")+
  ylab("Predicted biomass(mg)")+
  theme_light()+
  theme(panel.grid = element_blank())

mod6plot.p


# biomass vs diversity

mod7 <- lmer(Biomass ~ Shannon + (1|fSITE),
             data=srichnew)

summary(mod7)
drop1(mod7, test="Chi")

chkres(mod7, srichnew$Abundance, srichnew$fSITE)

r.squaredGLMM(mod7)

# retest using model II
mod7II <- lmodel2(Biomass ~ Shannon,
                   data = srichnew, nperm = 100)

mod7II


mod7plot <- ggplot(srichnew, aes(x=Shannon, y=Biomass))+
  geom_point()+
  geom_abline(intercept = summary(mod7)$coefficients[1,1], slope = summary(mod7)$coefficients[2,1])+
  xlim(0, 3) + ylim(20, 2000)+
  xlab("Shannon's H\n")+
  ylab("Biomass (mg)")+
  theme_light()+
  theme(panel.grid = element_blank())

mod7plot

# predicted

mod7p <- lmer(Pred_Biomass ~ Shannon + (1|fSITE),
             data=srichnew)

summary(mod7p)
drop1(mod7p, test="Chi")

chkres(mod7p, srichnew$Abundance, srichnew$fSITE)

r.squaredGLMM(mod7p)

# retest using model II
mod7pII <- lmodel2(Pred_Biomass ~ Shannon,
                  data = srichnew, nperm = 100)

mod7pII


mod7plot.p <- ggplot(srichnew, aes(x=Shannon, y=Pred_Biomass))+
  geom_point()+
  geom_abline(intercept = summary(mod7p)$coefficients[1,1], slope = summary(mod7p)$coefficients[2,1])+
  xlim(0, 3) + ylim(20, 2000)+
  xlab("Shannon's H\n")+
  ylab("Predicted biomass (mg)")+
  theme_light()+
  theme(panel.grid = element_blank())

mod7plot.p




# panel of graphs
m2a <- grid.arrange(mod5plot, mod6plot, mod7plot, ncol=1)


ggsave("Plots/Fig3a.svg", m2a, device = svg, width = 16, height = 24, units = "cm")


m2b <- grid.arrange(mod5plot.p, mod6plot.p, mod7plot.p, ncol=1)


ggsave("Plots/Fig3b.svg", m2b, device = svg, width = 16, height = 24, units = "cm")


m2c <- grid.arrange(mod5plot + ggtitle("Measured biomass"), mod5plot.p + ggtitle("Predicted biomass") + ylab(" "), 
                    mod6plot, mod6plot.p + ylab(" "), 
                    mod7plot, mod7plot.p + ylab(" "), ncol=2)


ggsave("Plots/Fig3c.svg", m2c, device = svg, width = 32, height = 24, units = "cm")







### finally, let's introduce objective 4

# we have the daily data for a Rothamsted trap run on the roof of the Biology Department
# between 1991-2006, by Terry Crawford

# we want to use this as an example dataset to demonstrate the use of our rules predicting body mass from forewing length
# specifically, we want to test whether there was any change in biomass over the 15 years for which the trap was run
# and whether this change is also visible in abundance, richness and diversity


### in a separate script, I've pulled out the records from the Heslington trap from the full Rothamsted dataset
# let's now read these data in

hes <- read.csv("Data/Heslington.csv", header = T)

# have a look at the data
summary(hes)

# we can see that there are three columns for species - giving the binomial name, the common name, and a unique code for each species
# there's also a column for date, and separately, year, and a column for abundance
# finally there's a column for trap name, but this is a legacy of the main dataset, as it gives "Heslington" in every row!


#### now let's play around, extracting plotting annual abundance, species-richness and diversity and looking for change over the 15-year period

# first we need yearly species totals (summing across dates)
# this will also usefully get rid of some redundant columns
yearly <- ddply(hes, .(YEAR,BINOMIAL,RIS_SPECIES), summarise,
                COUNT = sum(DAILY_COUNT))



### now, we have the formula above with which to calculate body mass for each of these moth species based on forewing length
# we'll want to read in the forewing lengths for each of these species 
# then use these to estimate body masses for each species
# and subsequently merge these back into the 'hes' dataframe
# and multiply by the number of individuals




## now merge the universal estimates into the yearly dataset

hes.yearly.mass <- merge(yearly,newdat_on)

# this has cut out a small number of things that were species aggs from the Heslington dataset - we can't accurately predict their wing-lengths
# and the only micro-moth recorded in this dataset (Plutella xylostella)

# we want to predict the body mass of each moth in this dataset using a normal distribution defined by the estimated dry mass and the standard error of that estimate
# to achieve this we need to melt down the dataset so that there is one row per individual moth
# (rather than one row per species, with a column for abundance)

hes.individuals <- data.frame(INDIVIDUAL = numeric(),
                              BINOMIAL = factor(),
                              YEAR = numeric(),
                              RIS_SPECIES = factor(),
                              FAMILY = factor(),
                              NAME = factor(),
                              FOREWING_LB = numeric(),
                              FOREWING_MED = numeric(),
                              FOREWING_UB = numeric(),
                              DRY_MASS = numeric(),
                              log_DRY_MASS = numeric(),
                              SE = numeric(),
                              COUNT = numeric())
                              

for (n in 1:nrow(hes.yearly.mass)){
  row <- hes.yearly.mass[n,]
  no.indivs <- row$COUNT[[1]]
  indivs <- data.frame(INDIVIDUAL = numeric(no.indivs))
  indivs$INDIVIDUAL <- c(1:no.indivs)
  indivs$BINOMIAL <- row$BINOMIAL[[1]]
  indivs$YEAR <- row$YEAR[[1]]
  indivs$RIS_SPECIES <- row$RIS_SPECIES[[1]]
  indivs$FAMILY <- row$FAMILY[[1]]
  indivs$NAME <- row$NAME[[1]]
  indivs$FOREWING_LB <- row$FOREWING_LB[[1]]
  indivs$FOREWING_MED <- row$FOREWING_MED[[1]]
  indivs$FOREWING_UB <- row$FOREWING_UB[[1]]
  indivs$DRY_MASS <- row$DRY_MASS[[1]]
  indivs$log_DRY_MASS <- row$log_DRY_MASS[[1]]
  indivs$SE <- row$SE[[1]]
  indivs$COUNT <- 1
  
  hes.individuals <- rbind(hes.individuals, indivs)
}

# to garner an estimate of the error around our estimates of sample-level biomass, we want to resample the *sample-level biomass* of each community multiple times
# that means, for each year's worth of data, repeatedly estimating the mass of each individual moth and summing the annual estimate
# then storing all these annual estimates, and measuring their mean and error

yearly.biomass <- data.frame(YEAR = numeric(),
                             EST.BIOMASS = numeric(),
                             EST.BIOMASS.SE = numeric())


for (year in 1991:2006){
  print(year)
  
  YEAR <- year
  
  year.indivs <- hes.individuals[which(hes.individuals$YEAR == year), ]     # extract that year's individuals

  year.biomass <- NULL     # seed an output vector
  
  # now generate 1000 independent estimates of the yearly biomass - each based on separate estimates of each moth's body mass

  for (x in 1:1000){
  
    if(round(x, -2) == x){
      print(x)
    }
    
    
    for (n in 1:nrow(year.indivs)){
      year.indivs$PRED_DRY_MASS[n] <- exp(rnorm(1, mean = year.indivs$log_DRY_MASS[n], sd = year.indivs$SE[n]))
    }
    biomass.estimate <- sum(year.indivs$PRED_DRY_MASS)
    
    year.biomass <- append(year.biomass, biomass.estimate)
  }
  
  # now take the mean and standard error of those estimates
  
  EST.BIOMASS <- mean(year.biomass)
  EST.BIOMASS.SE <- sd(year.biomass)/sqrt(length(year.biomass))
  
  out <- data.frame(YEAR,EST.BIOMASS,EST.BIOMASS.SE)
  yearly.biomass <- rbind(yearly.biomass, out)
}

summary(yearly.biomass)

# save and re-import this data as it took a while to generate
write.csv(yearly.biomass, "Data/EstBiomass.csv", row.names = F)

yearly.biomass <- read.csv("Data/EstBiomass.csv", header = T)




# then we can use ddply to calculate annual values for abundance, species richness and a range of diversity indices

srich <- ddply(hes.yearly.mass, .(YEAR), summarise,
               Abundance = sum(COUNT),
               Richness = length(COUNT),
               Richness1 = specnumber(COUNT),
               Shannon = diversity(COUNT, index = "shannon"),
               Simpson = diversity(COUNT, index = "simpson"),
               invSimpson = diversity(COUNT, index = "invsimpson"),
               Evenness = diversity(COUNT, index = "shannon")/log(specnumber(COUNT)))

# and merge in the biomass estimates

annual <- merge(yearly.biomass,srich)

## plot change in each variable over time

plot(annual$Abundance ~ annual$YEAR,
     ylim= range(0,500))
abline(lm(annual$Abundance ~ annual$YEAR))


plot(annual$Richness ~ annual$YEAR,
     ylim = range(0,100))
abline(lm(annual$Richness ~ annual$YEAR))


plot(annual$Shannon ~ annual$YEAR,
     ylim = range(0,5))
abline(lm(annual$Shannon ~ annual$YEAR))


plot(annual$Simpson ~ annual$YEAR,
     ylim = range(0,1))
abline(lm(annual$Simpson ~ annual$YEAR))


plot(annual$invSimpson ~ annual$YEAR,
     ylim = range(0,30))
abline(lm(annual$invSimpson ~ annual$YEAR))


plot(annual$Evenness ~ annual$YEAR,
     ylim = range(0,1))
abline(lm(annual$Evenness ~ annual$YEAR))


plot(annual$EST.BIOMASS~annual$YEAR,
     ylim = range(0,20000))
abline(lm(annual$EST.BIOMASS~annual$YEAR))


plot(annual$EST.BIOMASS~annual$Abundance,
     ylim = range(0,20000))
abline(lm(annual$EST.BIOMASS~annual$Abundance))



###HESLINGTON MODELS/GRAPHS- VARIABLES PLOTTED AGAINST EACH OTHER AND AGAINST YEAR###

# first, biomass vs abundance

hist(annual$EST.BIOMASS)
hist(log(annual$EST.BIOMASS))


mod8 <- lm(EST.BIOMASS ~ Abundance, data=annual)


summary(mod8)
drop1(mod8, test="F")

chkres(mod8, annual$Abundance, annual$YEAR)  # model residuals are mainly ok
#
#
#
#
#

# extract model R^2
r.squaredGLMM(mod8)

# retest using model II
mod8II <- lmodel2(EST.BIOMASS ~ Abundance,
                   data = annual, nperm = 100)

mod8II


## plot it out

mod8plot <- ggplot(annual, aes(x=Abundance, y=EST.BIOMASS))+
  geom_point()+
  geom_abline(intercept = summary(mod8)$coefficients[1,1], slope = summary(mod8)$coefficients[2,1])+
  xlim(0, 500) + ylim(0, 20000)+
  xlab("Abundance\n")+
  ylab("Biomass (mg)")+
  theme_light()+
  ggtitle("Historic dataset")+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(panel.grid = element_blank())
mod8plot



# next, biomass vs richness

mod9 <- lm(EST.BIOMASS ~ Richness,
           data=annual)

summary(mod9)
drop1(mod9, test="F")

chkres(mod9, annual$Richness, annual$YEAR)
#
#
#
#
#

r.squaredGLMM(mod9)

# retest using model II
mod9II <- lmodel2(EST.BIOMASS ~ Richness,
                  data = annual, nperm = 100)

mod9II



## plot it

mod9plot <- ggplot(annual, aes(x=Richness, y=EST.BIOMASS))+
  geom_point()+
  geom_abline(intercept = summary(mod9)$coefficients[1,1], slope = summary(mod9)$coefficients[2,1])+
  xlim(0, 100)+ylim(00,20000)+
  xlab("Species richness\n")+
  ylab("Biomass(mg)")+
  theme_light()+
  theme(panel.grid = element_blank())

mod9plot


# biomass vs diversity

mod10 <- lm(EST.BIOMASS ~ Shannon,
            data=annual)

summary(mod10)
drop1(mod10, test="F")

chkres(mod10, annual$Shannon, annual$YEAR)
#
#
#
#
#
#
r.squaredGLMM(mod10)

# retest using model II
mod10II <- lmodel2(EST.BIOMASS ~ Shannon,
                  data = annual, nperm = 100)

mod10II


mod10plot <- ggplot(annual, aes(x=Shannon, y=EST.BIOMASS))+
  geom_point()+
  xlim(0, 4) + ylim(0, 20000)+
  xlab("Shannon's H")+
  ylab("Biomass (mg)")+
  theme_light()+
  theme(panel.grid = element_blank())

mod10plot


# panel of graphs
m2d <- grid.arrange(mod8plot, mod9plot, mod10plot, ncol=1)


ggsave("Plots/Fig3b.svg", m2d, device = svg, width = 16, height = 24, units = "cm")

m2 <- grid.arrange(mod5plot + ggtitle("Sampled moths, \nmeasured biomass")  + xlab(" \n"),
                   mod5plot.p + ggtitle("Sampled moths, \npredicted biomass") + ylab(" "),
                   mod8plot + ggtitle("Historical dataset, \npredicted biomass") + xlab(" \n")+ ylab(" "), 
                   mod6plot + xlab(" \n"), mod6plot.p + ylab(" "), mod9plot + xlab(" \n") + ylab(" "), 
                   ncol=3)


ggsave("Plots/Fig3.svg", m2, device = svg, width = 28, height = 17, units = "cm")


## and now relate each thing to year

# abundance

mod11 <- lm(data=annual, Abundance~YEAR)

summary(mod11)
drop1(mod11, test="F")

chkres(mod11, annual$Abundance, annual$YEAR)
#
#
#
#
#

r.squaredGLMM(mod11)

# also do it with a rank correlation to check whether the overall decline is just a step function

cor.test(annual$YEAR, annual$Abundance,
         method = "spearman")

# and lastly establish the % change between 1991 and 2006
(1-(annual[16,'Abundance']/annual[1,'Abundance']))*100



mod11plot<- ggplot(annual, aes(x=YEAR, y=Abundance))+
  geom_point()+
  geom_abline(intercept=summary(mod11)$coefficients[1,1], slope = summary(mod11)$coefficients[2,1])+
  ylim(0, 500)+
  xlab("Year")+
  ylab("Abundance")+
  theme_light()+
  theme(panel.grid = element_blank())

mod11plot


# richness

mod12 <- lm(data = annual, Richness ~ YEAR)

summary(mod12)

drop1(mod12, test="F")

chkres(mod12, annual$Richness, annual$YEAR)
#
#
#
#
#

r.squaredGLMM(mod12)

# also do it with a rank correlation to check whether the overall decline is just a step function

cor.test(annual$YEAR, annual$Richness,
         method = "spearman")

# and lastly establish the % change between 1991 and 2006
(1-(annual[16,'Richness']/annual[1,'Richness']))*100




mod12plot <- ggplot(annual, aes(x=YEAR, y=Richness))+
  geom_point()+
  geom_abline(intercept=summary(mod12)$coefficients[1,1], slope = summary(mod12)$coefficients[2,1])+
  ylim(0, 100)+
  xlab("Year")+
  ylab("Species richness")+
  theme_light()+
  theme(panel.grid = element_blank())

mod12plot


# Shannon diversity

mod13 <- lm(data=annual, Shannon ~ YEAR)
summary(mod13)

drop1(mod13, test="F")

chkres(mod13, annual$Shannon, annual$YEAR)
#
#
#
#
#

r.squaredGLMM(mod13)


# also do it with a rank correlation to check whether the overall decline is just a step function

cor.test(annual$YEAR, annual$Shannon,
         method = "spearman")

# and lastly establish the % change between 1991 and 2006
(1-(annual[16,'Shannon']/annual[1,'Shannon']))*100


mod13plot <- ggplot(annual, aes(x=YEAR, y=Shannon))+
  geom_point()+
  ylim(0,4)+
  xlab("Year")+
  ylab("Shannon's H")+
  theme_light()+
  theme(panel.grid = element_blank())

mod13plot


# and biomass!

mod14 <- lm(data=annual, EST.BIOMASS ~ YEAR)
summary(mod14)

drop1(mod14, test="F")

chkres(mod14, annual$EST.BIOMASS, annual$YEAR)
#
#
#
#
#

r.squaredGLMM(mod14)

# also do it with a rank correlation to check whether the overall decline is just a step function

cor.test(annual$YEAR, annual$EST.BIOMASS,
         method = "spearman")

# and lastly establish the % change between 1991 and 2006
(1-(annual[16,'EST.BIOMASS']/annual[1,'EST.BIOMASS']))*100



mod14plot <- ggplot(annual, aes(x=YEAR, y=EST.BIOMASS))+
  geom_point()+
  geom_abline(intercept=summary(mod14)$coefficients[1,1], slope = summary(mod14)$coefficients[2,1])+
  ylim(0,20000)+
  xlab("Year")+
  ylab("Total biomass (mg)")+
  theme_light()+
  theme(panel.grid = element_blank())

mod14plot



m3 <- grid.arrange(mod11plot + xlab(" "), mod12plot + xlab(" "), mod14plot, ncol=1)

ggsave("Plots/Fig4.svg", m3, device = svg, width = 12, height = 24, units = "cm")



### drying curve ####
# now we want to plot out the loss in mass over time in a subset of ~100 moths
# the objective of this is to visually confirm that 7 days' drying is sufficient to accurately measure dry mass of our moths


# first read in the relevant data
dry <- read.csv("Data/drying_curve.csv", header=T)

summary(dry)

# there's a lot of missing data on the last day of measurements, so get rid of it
dry <- dry[,-length(dry)]

# wrangle it into a better shape
dry1 <- melt(dry, id.vars=c("ID"), variable.name = "Day", value.name = "DryMass")

summary(dry1)


# day is currently being treated as a factor so let's fix it to a numeric
# create relevant vectors
Day <- c("X0","X1","X2","X3","X4","X5","X6","X7","X8")
Day.num <- c(0:8)


# bind them together
days.frame <- data.frame(cbind(Day,Day.num))
days.frame

# merge it into the data
dry1 <- merge(dry1,days.frame)
dry1$Day.num <- as.numeric(as.character(dry1$Day.num))


# now, we have lots of NAs in this dataset, which are days on which given individuals were not weighed (because it was the weekend!)
dry1 <- dry1[complete.cases(dry1$DryMass), ]

summary(dry1)


# now calculate the mean daily mass of all individuals
mean.all <- ddply(dry1, .(Day.num), summarise,
                  DailyMean = mean(DryMass),
                  DailySE = sd(DryMass)/sqrt(length(DryMass)))


# calculate CIs
mean.all$UCI <- mean.all$DailyMean + (1.96*mean.all$DailySE)
mean.all$LCI <- mean.all$DailyMean - (1.96*mean.all$DailySE)


# separate out the top 20% of individuals by initial (wet) mass, excluding those not measured on day 8
cutoff <- round(nrow(dry)/5,0)

dry.good <- dry[which(!is.na(dry$X8)), ]

top20 <- top_n(dry.good, n = cutoff, wt = dry.good$X0)


# once again wrangle this into shape as before
top20a <- melt(top20, id.vars=c("ID"), variable.name = "Day", value.name = "DryMass")

summary(top20a)

top20a <- merge(top20a,days.frame)
top20a$Day.num <- as.numeric(as.character(top20a$Day.num))

top20a <- top20a[complete.cases(top20a$DryMass), ]

summary(top20a)



# take the mean of those top 20% as above
mean.top20 <- ddply(top20a, .(Day.num), summarise,
                    DailyMean = mean(DryMass),
                    DailySE = sd(DryMass)/sqrt(length(DryMass)))

mean.top20$UCI <- mean.top20$DailyMean + (1.96*mean.top20$DailySE)
mean.top20$LCI <- mean.top20$DailyMean - (1.96*mean.top20$DailySE)



# now we have all the means to graph it out

drying <- ggplot(mean.all, aes(x = Day.num, y = DailyMean)) +
  geom_point(data = top20a, aes(x = Day.num, y = DryMass, group = ID),
             shape = 19, colour = "grey70") +
  geom_line(data = top20a, aes(x = Day.num, y = DryMass, group = ID), colour = "grey70") +
  geom_point(shape = 19) + geom_line(linetype="dashed") +
  geom_errorbar(aes(ymin = LCI, ymax = UCI), width=0.5, size=0.5) +
  geom_point(data = mean.top20, aes(x = Day.num, y = DailyMean), shape = 19) +
  geom_line(data = mean.top20, aes(x = Day.num, y = DailyMean), size=0.75) +
  geom_errorbar(data = mean.top20, aes(ymin = LCI, ymax = UCI), width=0.5, size=0.5) +
  ylim(0,900)+
  scale_x_continuous(breaks = seq(0,10,2)) +
  xlab("Days since drying commenced")+
  ylab("Mass (mg)")+
  theme_light()+
  theme(panel.grid = element_blank())


drying

ggsave("Plots/FigS2.svg", drying, device = svg, width = 12, height = 12, units = "cm")






##### supplementary tests ####

### check weighing error
# to quantify the precision of Becci's weighing and measuring, we have done some re-measurement of a few individuals - 
# in total, 5 individuals per species of 2 species, weighing and measuring every individual 3 times

# read these measurements in
precis <- read.csv("Data/Reweighing.csv")

summary(precis)

# now collapse this down to get a mean and a standard error per individual, per dimension

precis.coll <- ddply(precis, .(Species, Individual, Bag.code, Dimension), summarise,
                     Mean = mean(Measurement),
                     SD = sd(Measurement),
                     SE = sd(Measurement)/sqrt(length(Measurement)),
                     CV = cv(Measurement))

# now take the means of these values for each species
precis.spec <- ddply(precis.coll, .(Species, Dimension), summarise,
                     Mean.SE = mean(SE),
                     SE.SE = sd(SE)/sqrt(length(SE)),
                     Mean.CV = mean(CV),
                     SE.CV = sd(CV)/sqrt(length(CV)))

precis.overall <- ddply(precis.coll, .(Dimension), summarise,
                        Mean.SE = mean(SE),
                        SE.SE = sd(SE)/sqrt(length(SE)),
                        Mean.CV = mean(CV),
                        SE.CV = sd(CV)/sqrt(length(CV)))

precis.overall


