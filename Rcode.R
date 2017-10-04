## R script for analysis in the article "Migratory animals feel the cost of infection: a meta-analysis across species"
#Alice Risely, Marcel Klaassen, & Bethany Hoye, J. Animal Ecology 2015
# Also available at github.com/Riselya

##R version 3.4.1
##contact riselya@gmail.com


library(metafor)
library(ggplot2)
library(tidyr)
library(dplyr)
library(glmulti)

#for phylogenetic tree
library(rotl)
library(ape)


##import data file "Supplementary File 3.csv" as 'all_data'

all_data <- read.csv("Supplementary File 3.csv")

#######################  DATA

### "all_data" is all data combined for infection status and intensity

###definitions for all variable

names(all_data)

#[1] "id"                   ##study ID                 
#[2] "rowid"                ##row ID (not continuous)             
#[3] "authors"            
#[4] "title"                
#[5] "journal"              
#[6] "year"                 
#[7] "abstract"             
#[8] "type"                 ##study design - observational or experimental           
#[9] "trait"                ##performance trait            
#[10] "measure"             ##how trait was actually measured           
#[11] "exp.methods"         ##the experimental methods used where relevant    
#[12] "species"             ##host species   
#[13] "species_latin"       ##host latin name
#[14] "taxa"                ##host kingdon   
#[15] "Order"               ##host order
#[16] "Family"              ##host family   
#[17] "strain"              ##parasite strain   
#[18] "parasite.taxa"       ##parasite type   
#[19] "life.history.measured"  ##list history of host when sampled
#[20] "setting"             ## host sampled in field or lab 
#[21] "subset"              ##if analyses were split by age/sex, which grouping was used 
#[22] "Infection.type"      ## single or multiple infection measured 
#[23] "migratory.leg2"      ##estimation of whether a migratory leg occured before infection and sampling  (Y/N)
#[24] "ss"                  ##sample size 
#[25] "ss_infected"         ## sample size of infected group 
#[26] "ss_healthy"          ##sample size of uninfected group 
#[27] "slope"               ##slope of relationship if presented 
#[28] "effect"              ##effect size presented 
#[29] "stat"                ##statistic used 
#[30] "presented.p"         ##p value presented 
#[31] "adj.p"               ## adjusted p - exact p value used in this study based on presented.p
#[32] "sig0.05"             ##whether study was significat to 0.05
#[33] "infection.measure"   ## whether observation on infection status or infection intensity 
#[34] "z"                   ##Fisher's z 
#[35] "var.z"               ## variation in Fisher's Z 
#[36] "g"                   ## Hedges' g 
#[37] "var.g"               ## variation in Hedges' g 
#[38] "var.g1"              ##capped var.g at 0.01 
#[39] "function."           ## the compute.es function used to calculate Fisher's z, Hedges' g and their variances. For info only
#[40] "Effect.direction"    ##effect direction (negative/positive)

unique(all_data$id) ##41 studies
str(all_data)

##recatgorise some variables
all_data$id<-factor(all_data$id)
all_data$rowid<-factor(all_data$rowid)
all_data$ss<-as.character(all_data$ss) #two step conversion
all_data$ss<-as.numeric(all_data$ss)
all_data$ss_infected<-as.character(all_data$ss_infected) #two step conversion
all_data$ss_infected<-as.numeric(all_data$ss_infected)
all_data$ss_healthy<-as.character(all_data$ss_healthy)
all_data$ss_healthy<-as.numeric(all_data$ss_healthy)
all_data$sig0.05<-factor(all_data$sig0.05)

##################### get phylogenetic tree for all host species to add OTT_ID (Tree of life) to dataframe. Need for later analyses.

############### supplementary material

## create seperate phylo trees for each status and intesity dataset

species<-unique(all_data$species_latin)

taxa<-tnrs_match_names(names= c("Acrocephalus arundinaceus",
                                "Anas crecca",
                                "Anas platyrhynchos",
                                "Anguilla anguilla",
                                "Arenaria interpres",
                                "Chen caerulescens",
                                "Clupea harengus",
                                "Cygnus columbianus",
                                "Danaus plexippus",
                                "Delichon urbica",
                                "Falco sparverius",
                                "Ficedula hypoleuca",
                                "Hirundo pyrrhonota",
                                "Hirundo rustica",
                                "Luscinia svecica",
                                "Merops apiaster",
                                "Oncorhynchus gorbuscha",
                                "Oncorhynchus mykiss",
                                "Oncorhynchus nerka",
                                "Progne subis",
                                "Salmo salar",
                                "Setophaga coronata",
                                "Setophaga magnolia",
                                "Sylvia atricapilla",
                                "Sylvia borin"))

taxa #all species assigned correctly

#make phylo tree

tree <- tol_induced_subtree(ott_ids = ott_id(taxa))
str(tree)

#plot
plot(tree, cex = .8, label.offset = .1, no.margin = TRUE)


##add species ott_ids to dataframe so that models can match correlation matrix phylo IDs to species in dataframe

tree$tip.label

all_data$species_ott[all_data$species_latin=="Acrocephalus arundinaceus"]<-"Acrocephalus_arundinaceus_ott800677"
all_data$species_ott[all_data$species_latin=="Anas crecca"]<-"Anas_crecca_ott656799"
all_data$species_ott[all_data$species_latin=="Anas platyrhynchos"]<-"Anas_platyrhynchos_ott765167"
all_data$species_ott[all_data$species_latin=="Anguilla anguilla"]<-"Anguilla_anguilla_ott854201"
all_data$species_ott[all_data$species_latin=="Arenaria interpres"]<-"Arenaria_interpres_ott821753"
all_data$species_ott[all_data$species_latin=="Chen caerulescens"]<-"Anser_caerulescens_ott190878"
all_data$species_ott[all_data$species_latin=="Clupea harengus"]<-"Clupea_harengus_ott1005932"
all_data$species_ott[all_data$species_latin=="Cygnus columbianus"]<-"Cygnus_columbianus_ott207360"
all_data$species_ott[all_data$species_latin=="Danaus plexippus"]<-"Danaus_plexippus_ott190091"
all_data$species_ott[all_data$species_latin=="Delichon urbica"]<-"Delichon_urbicum_ott36109"
all_data$species_ott[all_data$species_latin=="Falco sparverius"]<-"Falco_sparverius_ott964519"
all_data$species_ott[all_data$species_latin=="Ficedula hypoleuca"]<-"Ficedula_hypoleuca_ott22300"
all_data$species_ott[all_data$species_latin=="Hirundo pyrrhonota"]<-"Petrochelidon_pyrrhonota_ott302798"
all_data$species_ott[all_data$species_latin=="Hirundo rustica"]<-"Hirundo_rustica_ott1040135"
all_data$species_ott[all_data$species_latin=="Luscinia svecica"]<-"Luscinia_svecica_ott274225"
all_data$species_ott[all_data$species_latin=="Merops apiaster"]<-"Merops_apiaster_ott755107"
all_data$species_ott[all_data$species_latin=="Oncorhynchus gorbuscha"]<-"Oncorhynchus_gorbuscha_ott739927"
all_data$species_ott[all_data$species_latin=="Oncorhynchus mykiss"]<-"Oncorhynchus_mykiss_ott165368"
all_data$species_ott[all_data$species_latin=="Oncorhynchus nerka"]<-"Oncorhynchus_nerka_ott165375"
all_data$species_ott[all_data$species_latin=="Progne subis"]<-"Progne_subis_ott621799"
all_data$species_ott[all_data$species_latin=="Salmo salar"]<-"Salmo_salar_ott688328"
all_data$species_ott[all_data$species_latin=="Setophaga coronata"]<-"Setophaga_coronata_ott451168"
all_data$species_ott[all_data$species_latin=="Setophaga magnolia"]<-"Setophaga_magnolia_ott532751"
all_data$species_ott[all_data$species_latin=="Sylvia atricapilla"]<-"Sylvia_atricapilla_ott726312"
all_data$species_ott[all_data$species_latin=="Sylvia borin"]<-"Sylvia_borin_ott261410"




###################    METAANALYSIS ON EFFECT OF INFECTION STATUS ON PEFORMANCE

##subset observations on infection status

status<-subset(all_data, infection.measure=="Infection status") ##subset observations on infection status

unique(status$id) ##35 studies included
unique(status$measure) ##20 measures
table(status$Effect.direction)
table(status$sig0.05)

###order categories

status$trait<-factor(status$trait, levels=c("Body stores","Refuelling","Movement","Phenology","Survival"))
status$life.history.measured<-factor(status$life.history.measured, levels=c("Migration","Breeding","Non-breeding","Lab"))
status$parasite.taxa<-factor(status$parasite.taxa, levels=c("Protozoa","Virus","Mites","Helminth","Multiple"))
status$Order<-factor(status$Order, levels=c("Passeriformes","Anseriformes","Salmoniformes","Lepidoptera","Falconiformes",
                                            "Anguilliformes","Clupeiformes","Coraciiformes","Charadriiformes"))
status$species_latin<-factor(status$species_latin)
status<-status[order(status$trait, status$g),]

###forest plot

forest.default(status$g, status$var.g1, slab = status$id)

##change names of hedges' g and variance to standardized names

status$yi<-status$g
status$vi<-status$var.g ##variance not capped
status$vi1<-status$var.g1 ##capped to 0.01


####exclude data points that use same animals and same trait as non-independent
##Rowid's that use same animals and same traits

# 9 + 10 + 11 + 12 van Dijk et al. OIKOS
# 26 + 27 +28 Lopez et al. PLOS ONE
# 50 + 51 Arizaga et al. ARDEOLA
# 60 + 61 Latorre-Margalef et al. PROC B
# 63 + 64 Silvertsgard et al. HYDROBIOGIA
# 65 + 66 Silvertsgard et al. HYDROBIOGIA
# 83 + 84 Bradley et al ECOLOGY 
# 105 + 106 Ratii OECOLOGIA
# 155 + 156 Hoye et al. INTEGRATIVE & COMP BIOLOGY
# 157 + 158 Hoye et al. INTEGRATIVE & COMP BIOLOGY
# 159 + 160 Hoye et al. INTEGRATIVE & COMP BIOLOGY


status_independent<-subset(status, rowid!=27 & rowid!= 28  & 
                             rowid!= 10 & rowid!= 11 & rowid!= 12 & rowid!= 61 & 
                             rowid!= 84 & rowid!= 106 & rowid!= 64 & rowid!= 66 & 
                             rowid!= 160 & rowid!= 158 & rowid!= 156 & rowid!= 51 )

##repeat with different set of points

#status_independent1<-subset(status, rowid!=26 & rowid!= 28 & 
#                            rowid!= 9 & rowid!= 11 & rowid!= 12 & rowid!= 60 & 
#                        rowid!= 83 & rowid!= 105 & rowid!= 63 & rowid!= 65 & 
#                        rowid!= 159 & rowid!= 157 & rowid!= 155 & rowid!= 50 )

## run analysis excluding these points

status<-status_independent


##forest plots by trait

forest.default(status$g, status$vi1,  cex=1, xlab="Hedges g")


##Kendall rank correlation test

cor.test(status$yi, status$vi, method="kendall", alternative="less")
coef(lm(yi ~ vi, data = status))

##supplementary figure S7

ggplot(data=status, aes(x=vi, y = yi))+geom_point(size = 3)+ labs(x="Variance in Hedges' g", y = "Hedges' g")+
  geom_abline(intercept= -0.0399, slope = -2.61)  +
  theme_bw() + 
  theme(text = element_text(size=25),plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank() )+
  theme(panel.border= element_blank())+
  theme(axis.text=element_text(size = 12))+
  theme(axis.line.x = element_line(color="black", size = 0.5),
        axis.line.y = element_line(color="black", size = 0.5))

#############create phylogentic correlation matrix to control for phylogeny

unique(status$species_ott)
unique(status$species_latin)

taxa<-tnrs_match_names(names= c("Danaus plexippus",
                                "Sylvia atricapilla",
                                "Acrocephalus arundinaceus",
                                "Anas platyrhynchos",
                                "Hirundo pyrrhonota",
                                "Progne subis",
                                "Falco sparverius",
                                "Setophaga magnolia",
                                "Setophaga coronata",
                                "Salmo salar",
                                "Anguilla anguilla",
                                "Anas crecca",
                                "Merops apiaster",
                                "Oncorhynchus mykiss",
                                "Cygnus columbianus",
                                "Clupea harengus",
                                "Sylvia borin",
                                "Delichon urbica",
                                "Arenaria interpres",
                                "Hirundo rustica",
                                "Ficedula hypoleuca",
                                "Chen caerulescens",
                                "Luscinia svecica"))

tree <- tol_induced_subtree(ott_ids = ott_id(taxa))
str(tree)


##estimate branch lengths using APE


tree1<-compute.brlen(tree, method = "Grafen", power = 1)

plot(tree1, cex = .8, label.offset = .1, no.margin = TRUE)
str(tree1)

##create correlation matrix

corr_matrix<-vcv(tree1, corr=TRUE)
#corr_matrix


############################################ MODEL SELECTION ##################################################################################


# see Methods section
# explore which variables have an impact on effect size by comparing all models 
# use capped variance (vi1)

# use gmlulti package for model selction
# see http://www.metafor-project.org/doku.php/tips:model_selection_with_glmulti for coding methods

rma.glmulti <- function(formula, data, ...) {
  rma.mv(as.formula(paste(deparse(formula))), vi1, data=data, method="ML", ...)
}

### global model with 4 ecologically relevant variables. Max number of variables allowed in model is 2. Study ID and species (phylogeny) as random effects

res <- glmulti(yi ~ trait+life.history.measured+type+parasite.taxa,random = list(~1 | id, ~ 1 | species_ott), R = list(species_ott = corr_matrix), data=status,
               level=1, fitfunction=rma.glmulti, crit="aicc",  maxsize= 2)


tmp <- weightable(res)
tmp


##variable importance
plot(res, type="s")
plot(res, type="w")
print(res)

## lets check out the best models

summary(res@objects[[1]])
funnel(res@objects[[1]], yaxis="wi")

summary(res@objects[[2]])
funnel(res@objects[[2]], yaxis="wi")

summary(res@objects[[3]])
funnel(res@objects[[3]], yaxis="wi")


##################### Fit null and best models ###########################################


## first, null model (in text)

null.model<-rma.mv(yi, vi1,random = list(~1 | id, ~ 1 | species_ott), R = list(species_ott = corr_matrix), data=status, method="REML")
summary(null.model)
funnel(null.model, xlab = "Hedges' g", cex.lab = 1.5, cex.axis=1.5)

## heterogeneity (see http://www.metafor-project.org/doku.php/tips:i2_multilevel_multivariate for methods)

W <- diag(1/status$vi1)
X <- model.matrix(null.model)
P <- W - W %*% X %*% solve(t(X) %*% W %*% X) %*% t(X) %*% W
100 * sum(null.model$sigma2) / (sum(null.model$sigma2) + (null.model$k-null.model$p)/sum(diag(P)))
# = 56.8% #Total variance

100 * null.model$sigma2 / (sum(null.model$sigma2) + (null.model$k-null.model$p)/sum(diag(P)))
# sigma1 = 28.6%   sigma2 = 27.5 # all variance attributed to phylogeny


#with trait as moderator (best model)

best.model1<-rma.mv(yi, vi1, mods=~ trait, random = list(~1 | id, ~ 1 | species_ott), R = list(species_ott = corr_matrix), data=status)

#Table 3a
summary(best.model1)

funnel(best.model1)
funnel(best.model1, yaxis="wi") ##weighting of points

## heterogeneity

W <- diag(1/status$vi1)
X <- model.matrix(best.model1)
P <- W - W %*% X %*% solve(t(X) %*% W %*% X) %*% t(X) %*% W
100 * sum(best.model1$sigma2) / (sum(best.model1$sigma2) + (best.model1$k-best.model1$p)/sum(diag(P)))
# = 17.8% #Total variance

100 * best.model1$sigma2 / (sum(best.model1$sigma2) + (best.model1$k-best.model1$p)/sum(diag(P)))
# sigma1 = 0%   sigma2 = 17.8 # all variance attributed to phylogeny

## use uncapped variance to compare (table S5)

best.model2<-rma.mv(yi, vi, mods=~ trait ,random = list(~1 | id, ~ 1 | species_ott), R = list(species_ott = corr_matrix), data=status)
summary(best.model2)
funnel(best.model2, yaxis="wi") ##weighting of points


#heterogeneity

W <- diag(1/status$vi)
X <- model.matrix(best.model2)
P <- W - W %*% X %*% solve(t(X) %*% W %*% X) %*% t(X) %*% W
100 * sum(best.model2$sigma2) / (sum(best.model2$sigma2) + (best.model2$k-best.model2$p)/sum(diag(P)))
# = 89.0% #Total variance

100 * best.model2$sigma2 / (sum(best.model2$sigma2) + (best.model2$k-best.model2$p)/sum(diag(P)))
# sigma1 = 42.5%   sigma2 = 46.5 # all variance attributed to phylogeny

################### back to capped variaces

##make table of best model (table 3a) with capped variances

best.model.table<-data.frame(best.model1$b)
best.model.table$se<-best.model1$se
best.model.table$zval<-best.model1$zval
best.model.table$p<-best.model1$pval
best.model.table$l.ci<-best.model1$ci.lb
best.model.table$u.ci<-best.model1$ci.ub

#write.csv(best.model.table, "best.model.csv")


##plot resids

F1 <- fitted(best.model1)
E1 <- resid(best.model1)
plot(x = F1, 
     y = E1, 
     xlab = "Fitted values",
     ylab = "Residuals", 
     cex.lab = 1.5)
abline(h = 0, lty = 2)


##some outliers. 
hist(E1, breaks=10)

##save fitted and residuals and hat values

status$resids<-resid(best.model1)
status$predicted<-fitted(best.model1)
status$leverage<-hatvalues(best.model1)

mean(status$leverage)
#plot(status$leverage)

#subset(status, resids < -1.0)

##three outliers in residuals are from the van gils paper on avian influenza on phenology, movement and refuelling
##exclude these and rerun from line 205 - doesnt make much difference although reduces model residual heterogeneity:
#status<-subset(status, resids > -1) ##rerun without outliers

##PLOT ESTIMATES (Figure 3)


best.model1<-rma.mv(yi, vi1, mods=~ trait ,random = list(~1 | id, ~ 1 | species_ott), R = list(species_ott = corr_matrix), data=status)
summary(best.model1)

# get predicted values for each trait to plot
# see http://www.metafor-project.org/doku.php/tips:testing_factors_lincoms

trait.pred<-predict(best.model1, newmods=rbind(c(0,0,0,0),c(1,0,0,0),c(0,1,0,0),
                                               c(0,0,1,0),c(0,0,0,1)))

trait<-data.frame(trait.pred$pred)
trait$se<-trait.pred$se
trait$l.ci<-trait.pred$ci.lb
trait$u.ci<-trait.pred$ci.ub


trait$level<-factor(c("Body stores","Refuelling","Movement","Phenology","Survival"))
trait$level<-factor(trait$level, levels=c("Body stores","Refuelling","Movement","Phenology","Survival"))
colnames(trait)[1]<-"predicted"

predictedvalues<-trait

##create SE values

predictedvalues$l.se<-predictedvalues$predicted - predictedvalues$se
predictedvalues$u.se<-predictedvalues$predicted + predictedvalues$se


##boxplot

##overlayed with raw data with proportional weights

##Fig 3a)

ggplot(data=predictedvalues, aes(x=level, y = predicted, fill = Order))+
  geom_boxplot(aes(ymin=l.ci, lower=l.se, 
                   middle = predicted, upper = u.se, ymax=u.ci, width = 0.6), stat="identity", fill="lightgray", lwd=0.7)+
  geom_jitter(data = status, aes(x= trait, y = yi, size = 1/vi1, col = Order),width=0.15, alpha=1/2)+
  labs(y="Predicted Hedges' g", x ="")+geom_hline(yintercept=0, linetype="longdash")+
  theme_bw() +   theme(text = element_text(size=22),axis.text.x = element_text(angle=45, vjust=0.5),
                       plot.background = element_blank(),
                       panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank() )+
  theme(panel.border= element_blank())+
  theme(axis.line.x = element_line(color="black", size = 0.4),
        axis.line.y = element_line(color="black", size = 0.4))+
  scale_y_continuous(breaks=seq(-2,0.5, 0.5))+scale_color_brewer(type = "qual", palette = "Set1")+
  guides(colour = guide_legend(override.aes = list(size=4)))


### Fig S8 

##order parasite strains by type (Haemoparasites in birds, ectoparasites in birds, etc)

status$strain<-factor(status$strain, levels=c("Haemoparasites",
                                              "Leucocytozoon",
                                              "Trypanosoma",
                                              
                                              "Feather mites",
                                              "Ectoparasites",
                                              
                                              "Cestodes (general)",
                                              
                                              "LPAIV",
                                              
                                              "Anguillicola crassus",
                                              "Cloacotaenia megalops",
                                              "Ichthyophonus",
                                              "Ichthyophonus hoferi",
                                              "L. salmonis",
                                              
                                              "Ophryocystis elektroscirrha",
                                              
                                              "Multiple"))

ggplot(data=predictedvalues, aes(x=level, y = predicted, fill = strain))+
  geom_boxplot(aes(ymin=l.ci, lower=l.se, 
                   middle = predicted, upper = u.se, ymax=u.ci, width = 0.7), stat="identity", fill="lightgray", lwd=0.7)+
  geom_jitter(data = status, aes(x= trait, y = yi, size = 1/vi1, col = strain),width=0.3, alpha=0.7)+
  labs(y="Predicted Hedges' g", x ="")+geom_hline(yintercept=0, linetype="longdash")+
  theme_bw() +   theme(text = element_text(size=22),axis.text.x = element_text(angle=45, vjust=0.5),
                       plot.background = element_blank(),
                       panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank() )+
  theme(panel.border= element_blank())+
  theme(axis.line.x = element_line(color="black", size = 0.4),
        axis.line.y = element_line(color="black", size = 0.4))+
  scale_y_continuous(breaks=seq(-2,0.5, 0.5))+scale_color_manual(values= c("red","red","red","blue","blue",
                                                                           "magenta","forestgreen","goldenrod1","goldenrod1",
                                                                           "goldenrod1","goldenrod1","goldenrod1","gray33","cyan2"))+
  guides(colour = guide_legend(override.aes = list(size=4)))


#############################################################################################################


#########################modelling traits separately


survival<-subset(status, trait=="Survival")
Bodystores<-subset(status, trait=="Body stores")
phenology<-subset(status, trait=="Phenology")
movement<-subset(status, trait=="Movement")
refuelling<-subset(status, trait=="Refuelling")
#######################

null1<-rma.mv(yi, vi1, random = list(~1 | id, ~ 1 | species_ott), R = list(species_ott = corr_matrix), data=Bodystores) ##null model
summary(null1)

W <- diag(1/Bodystores$vi1)
X <- model.matrix(null1)
P <- W - W %*% X %*% solve(t(X) %*% W %*% X) %*% t(X) %*% W
100 * sum(null1$sigma2) / (sum(null1$sigma2) + (null1$k-null1$p)/sum(diag(P)))
# = 31.9% #Total variance

100 * null1$sigma2 / (sum(null1$sigma2) + (null1$k-null1$p)/sum(diag(P)))
# sigma1 = 1.9%   sigma2 = 30.02 # all variance attributed to phylogeny
#################

null2<-rma.mv(yi, vi1, random = list(~1 | id, ~ 1 | species_ott), R = list(species_ott = corr_matrix), data=refuelling) ##null model
summary(null2)

W <- diag(1/refuelling$vi1)
X <- model.matrix(null2)
P <- W - W %*% X %*% solve(t(X) %*% W %*% X) %*% t(X) %*% W
100 * sum(null2$sigma2) / (sum(null2$sigma2) + (null2$k-null2$p)/sum(diag(P)))
# = 0% #Total variance

100 * null2$sigma2 / (sum(null2$sigma2) + (null2$k-null2$p)/sum(diag(P)))
# sigma1 = 0%   sigma2 = 0%

###################

null3<-rma.mv(yi, vi1, random = list(~1 | id, ~ 1 | species_ott), R = list(species_ott = corr_matrix), data=movement) ##null model
summary(null3)

W <- diag(1/movement$vi1)
X <- model.matrix(null3)
P <- W - W %*% X %*% solve(t(X) %*% W %*% X) %*% t(X) %*% W
100 * sum(null3$sigma2) / (sum(null3$sigma2) + (null3$k-null3$p)/sum(diag(P)))
# = 67.7% #Total variance

100 * null3$sigma2 / (sum(null3$sigma2) + (null3$k-null3$p)/sum(diag(P)))
# sigma1 = 67.7%   sigma2 = 0%

#################################

null4<-rma.mv(yi, vi1, random = list(~1 | id, ~ 1 | species_ott), R = list(species_ott = corr_matrix), data=phenology) ##null model
summary(null4)

W <- diag(1/phenology$vi1)
X <- model.matrix(null4)
P <- W - W %*% X %*% solve(t(X) %*% W %*% X) %*% t(X) %*% W
100 * sum(null4$sigma2) / (sum(null4$sigma2) + (null4$k-null4$p)/sum(diag(P)))
# = 40.9% #Total variance

100 * null4$sigma2 / (sum(null4$sigma2) + (null4$k-null4$p)/sum(diag(P)))
# sigma1 = 40.9%   sigma2 = 0%

##########################

null5<-rma.mv(yi, vi1,random = list(~1 | id, ~ 1 | species_ott), R = list(species_ott = corr_matrix), data=survival) ##null model
summary(null5)

W <- diag(1/survival$vi1)
X <- model.matrix(null5)
P <- W - W %*% X %*% solve(t(X) %*% W %*% X) %*% t(X) %*% W
100 * sum(null5$sigma2) / (sum(null5$sigma2) + (null5$k-null5$p)/sum(diag(P)))
# = 40.9% #Total variance

100 * null5$sigma2 / (sum(null5$sigma2) + (null5$k-null5$p)/sum(diag(P)))
# sigma1 = 40.9%   sigma2 = 0%

################################### PLOT Fig 3b

pred1<-predict(null1)
pred2<-predict(null2)
pred3<-predict(null3)
pred4<-predict(null4)
pred5<-predict(null5)

trait<-data.frame(rbind(pred1,pred2,pred3,pred4,pred5))

trait$pred<-as.numeric(trait$pred)
trait$se<-as.numeric(trait$se)
trait$ci.lb<-as.numeric(trait$ci.lb)
trait$ci.ub<-as.numeric(trait$ci.ub)

trait<-trait[,1:4]



trait$level<-factor(c("Body stores","Refuelling","Movement","Phenology","Survival"))
trait$level<-factor(trait$level, levels=c("Body stores","Refuelling","Movement","Phenology","Survival"))
colnames(trait)[1]<-"predicted"
colnames(trait)[3]<-"l.ci"
colnames(trait)[4]<-"u.ci"

predictedvalues<-trait
str(predictedvalues)

predictedvalues$l.se<-predictedvalues$predicted - predictedvalues$se
predictedvalues$u.se<-predictedvalues$predicted + predictedvalues$se

##Figure 3b

ggplot(data=predictedvalues, aes(x=level, y = predicted, fill = parasite.taxa))+
  geom_boxplot(aes(ymin=l.ci, lower=l.se, 
                   middle = predicted, upper = u.se, ymax=u.ci, width = 0.6), stat="identity", fill="lightgray", lwd=0.7)+
  geom_jitter(data = status, aes(x= trait, y = yi, size = 1/vi1, col = parasite.taxa),width=0.15, alpha=0.4)+
  labs(y="Predicted Hedges' g", x ="")+geom_hline(yintercept=0, linetype="longdash")+
  theme_bw() +   theme(text = element_text(size=22),axis.text.x = element_text(angle=45, vjust=0.5),
                       plot.background = element_blank(),
                       panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank() )+
  theme(panel.border= element_blank())+
  theme(axis.line.x = element_line(color="black", size = 0.4),
        axis.line.y = element_line(color="black", size = 0.4))+
  scale_y_continuous(breaks=seq(-2,0.5, 0.5))+scale_color_manual(values=c("red","blue","cyan3", "gold","magenta"))+
  guides(colour = guide_legend(override.aes = list(size=4)))




###########################################SENSITIVITY ANALYSES##################################

###use uncapped variance instead of capped variance (replace vi1 with vi in all models) and rerun 
#best.model1<-rma.mv(yi, vi, mods=~ trait + Order, random = ~ 1 | id, data=status)
#summary(best.model1)


##exclude insects in case this biases models and rerun
#status<-subset(status, taxa!=”Insect”)

##exclude outliers for best model and rerun
#status<-subset(status, resids > -1)

###Does not change models overall

##################################################################################################################################



####################################           INTENSITY  META-ANALYSIS              ######################################
###same methods as with status

intensity<-subset(all_data, infection.measure=="Intensity")

###exclude points that use the same animals
#intensity<-subset(intensity, rowid != 4 & rowid!=5 & rowid!=134 & rowid!=22)

intensity<-intensity[order(intensity$trait, intensity$z),]
forest.default(intensity$z, intensity$var.z)


##order categories

intensity$trait<-factor(intensity$trait, levels=c("Body stores","Movement","Phenology","Survival"))

intensity$parasite.taxa<-factor(intensity$parasite.taxa, levels=c("Protozoa","Virus","Mites","Helminth","Myxospora"))
intensity$Order<-factor(intensity$Order, levels=c("Passeriformes","Anseriformes","Salmoniformes","Lepidoptera","Falconiformes",
                                                  "Anguilliformes","Clupeiformes","Coraciiformes","Charadriiformes"))

##make phylo tree

intensity$species_latin<-factor(intensity$species_latin)

unique(intensity$species_latin)
taxa<-tnrs_match_names(names= c("Acrocephalus arundinaceus",
                                "Anas platyrhynchos",
                                "Anguilla anguilla",
                                "Danaus plexippus",
                                "Falco sparverius",
                                "Hirundo rustica",
                                "Oncorhynchus gorbuscha",
                                "Oncorhynchus nerka",
                                "Setophaga coronata"))

tree <- tol_induced_subtree(ott_ids = ott_id(taxa))
str(tree)

##estimate branch lengths

tree1<-compute.brlen(tree, method = "Grafen", power = 1)

#plot
plot(tree1, cex = .8, label.offset = .1, no.margin = TRUE)
str(tree1)

##create correlation matrix

corr_matrix<-vcv(tree1, corr=TRUE)
str(corr_matrix)


###univariate models with each of the six explanatory variables

model.intensity<-rma.mv(z, var.z, mods = ~ trait, random = list(~1 | id, ~ 1 | species_ott), R = list(species_ott = corr_matrix),  data=intensity)
model.intensity2<-rma.mv(z, var.z, mods = ~ parasite.taxa, random = list(~1 | id, ~ 1 | species_ott), R = list(species_ott = corr_matrix),  data=intensity)
model.intensity3<-rma.mv(z, var.z, mods = ~ life.history.measured, random = list(~1 | id, ~ 1 | species_ott), R = list(species_ott = corr_matrix),  data=intensity)
model.intensity4<-rma.mv(z, var.z, mods = ~ type, random = list(~1 | id, ~ 1 | species_ott), R = list(species_ott = corr_matrix),  data=intensity)
model.intensity5<-rma.mv(z, var.z, random = list(~1 | id, ~ 1 | species_ott), R = list(species_ott = corr_matrix),  data=intensity) #null model


aicc(model.intensity) #trait
aicc(model.intensity2) #parasite type
aicc(model.intensity3) # Life history
aicc(model.intensity4) #study design
aicc(model.intensity5) #null


##best model includes just trait (AICC = 45.1)

summary(model.intensity) #best model
summary(model.intensity5) #null model

#I2 null model
W <- diag(1/intensity$var.z)
X <- model.matrix(model.intensity5)
P <- W - W %*% X %*% solve(t(X) %*% W %*% X) %*% t(X) %*% W
100 * sum(model.intensity5$sigma2) / (sum(model.intensity5$sigma2) + (model.intensity5$k-model.intensity5$p)/sum(diag(P)))
# = 66.9% #Total variance

100 * model.intensity5$sigma2 / (sum(model.intensity5$sigma2) + (model.intensity5$k-model.intensity5$p)/sum(diag(P)))
# sigma1 = 22.5 % sigma2 = 44.5%

##I2 best model

W <- diag(1/intensity$var.z)
X <- model.matrix(model.intensity)
P <- W - W %*% X %*% solve(t(X) %*% W %*% X) %*% t(X) %*% W
100 * sum(model.intensity$sigma2) / (sum(model.intensity$sigma2) + (model.intensity$k-model.intensity$p)/sum(diag(P)))
# = 73.08% #Total variance

100 * model.intensity$sigma2 / (sum(model.intensity$sigma2) + (model.intensity$k-model.intensity$p)/sum(diag(P)))
# sigma1 = 39.4%   sigma2 = 33.7%


##table - Table 4a

best.model.table<-data.frame(model.intensity$b)
best.model.table$se<-model.intensity$se
best.model.table$zval<-model.intensity$zval
best.model.table$p<-model.intensity$pval
best.model.table$l.ci<-model.intensity$ci.lb
best.model.table$u.ci<-model.intensity$ci.ub

#write.csv(best.model.table, "best.model.intensity.csv")

##plot model estimates (Fig 4)

trait.pred<-predict(model.intensity, newmods=rbind(c(0,0,0),c(1,0,0),c(0,1,0),
                                                   c(0,0,1)))

trait<-data.frame(trait.pred$pred)
trait$se<-trait.pred$se
trait$l.ci<-trait.pred$ci.lb
trait$u.ci<-trait.pred$ci.ub


trait$level<-factor(c("Body stores","Movement","Phenology","Survival"))
trait$level<-factor(trait$level, levels=c("Body stores","Movement","Phenology","Survival"))
colnames(trait)[1]<-"predicted"

predictedvalues<-trait

predictedvalues$l.se<-predictedvalues$predicted - predictedvalues$se
predictedvalues$u.se<-predictedvalues$predicted + predictedvalues$se


##boxplot

##overlayed with raw data with proportional weights

##Figure 4b

ggplot(data=predictedvalues, aes(x=level, y = predicted, fill = parasite.taxa))+
  geom_boxplot(aes(ymin=l.ci, lower=l.se, 
                   middle = predicted, upper = u.se, ymax=u.ci, width = 0.5), stat="identity", fill="lightgray", lwd=0.7)+
  geom_jitter(data = intensity, aes(x= trait, y = z, size = 1/var.z, col = parasite.taxa),width=0.15, alpha=1/2)+
  labs(y="Fisher's z", x = "")+geom_hline(yintercept=0, linetype="longdash")+
  theme_bw() +   theme(text = element_text(size=22),axis.text.x = element_text(angle=45, vjust=0.5),
                       plot.background = element_blank(),
                       panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank() )+
  theme(panel.border= element_blank())+
  theme(axis.line.x = element_line(color="black", size = 0.4),
        axis.line.y = element_line(color="black", size = 0.4))+
  scale_color_manual(values=c("red3","royalblue3","cyan3", "orange","magenta"))+
  guides(colour = guide_legend(override.aes = list(size=4)))


######################model intensity separately

survival<-subset(intensity, trait=="Survival")
Bodystores<-subset(intensity, trait=="Body stores")
phenology<-subset(intensity, trait=="Phenology")
movement<-subset(intensity, trait=="Movement")

########################################

null1<-rma.mv(z, var.z, random = list(~1 | id, ~ 1 | species_ott), R = list(species_ott = corr_matrix), data=Bodystores) ##null model
summary(null1)

##I2

W <- diag(1/Bodystores$var.z)
X <- model.matrix(null1)
P <- W - W %*% X %*% solve(t(X) %*% W %*% X) %*% t(X) %*% W
100 * sum(null1$sigma2) / (sum(null1$sigma2) + (null1$k-null1$p)/sum(diag(P)))
# = 17.34% #Total variance

100 * null1$sigma2 / (sum(null1$sigma2) + (null1$k-null1$p)/sum(diag(P)))
# sigma1 = 0%   sigma2 = 17.34%

#######################################

#Does not converge with phylogeny as random effect, so removed.

null3<-rma.mv(z, var.z, random = ~1 | id, data=movement) ##null model
summary(null3)

#i2
W <- diag(1/movement$var.z)
X <- model.matrix(null3)
P <- W - W %*% X %*% solve(t(X) %*% W %*% X) %*% t(X) %*% W
100 * sum(null3$sigma2) / (sum(null3$sigma2) + (null3$k-null3$p)/sum(diag(P)))
# = 0% #Total variance

100 * null3$sigma2 / (sum(null3$sigma2) + (null3$k-null3$p)/sum(diag(P)))
# sigma1 = 0%   # no sigma two

###########################################

null4<-rma.mv(z, var.z, random = list(~1 | id, ~ 1 | species_ott), R = list(species_ott = corr_matrix), data=phenology) ##null model
summary(null4)

#i2
W <- diag(1/phenology$var.z)
X <- model.matrix(null4)
P <- W - W %*% X %*% solve(t(X) %*% W %*% X) %*% t(X) %*% W
100 * sum(null4$sigma2) / (sum(null4$sigma2) + (null4$k-null4$p)/sum(diag(P)))
# = 88.8% #Total variance

100 * null4$sigma2 / (sum(null4$sigma2) + (null4$k-null4$p)/sum(diag(P)))
# sigma1 = 40.9%   # sigma2 = 47.9%

########################################

null5<-rma.mv(z, var.z, random = list(~1 | id, ~ 1 | species_ott), R = list(species_ott = corr_matrix), data=survival) ##null model
summary(null5)

#i2
W <- diag(1/survival$var.z)
X <- model.matrix(null5)
P <- W - W %*% X %*% solve(t(X) %*% W %*% X) %*% t(X) %*% W
100 * sum(null5$sigma2) / (sum(null5$sigma2) + (null5$k-null5$p)/sum(diag(P)))
# = 20.1% #Total variance

100 * null5$sigma2 / (sum(null5$sigma2) + (null5$k-null5$p)/sum(diag(P)))
# sigma1 = 0%   # sigma2 = 20.1%

################ plot Fig 4b

pred1<-predict(null1)
pred3<-predict(null3)
pred4<-predict(null4)
pred5<-predict(null5)

trait<-data.frame(rbind(pred1,pred3,pred4,pred5))

trait$pred<-as.numeric(trait$pred)
trait$se<-as.numeric(trait$se)
trait$ci.lb<-as.numeric(trait$ci.lb)
trait$ci.ub<-as.numeric(trait$ci.ub)

trait<-trait[,1:4]

trait$level<-factor(c("Body stores","Movement","Phenology","Survival"))
trait$level<-factor(trait$level, levels=c("Body stores","Movement","Phenology","Survival"))
colnames(trait)[1]<-"predicted"
colnames(trait)[3]<-"l.ci"
colnames(trait)[4]<-"u.ci"

predictedvalues<-trait
str(predictedvalues)

predictedvalues$l.se<-predictedvalues$predicted - predictedvalues$se
predictedvalues$u.se<-predictedvalues$predicted + predictedvalues$se

##Fig 4b

ggplot(data=predictedvalues, aes(x=level, y = predicted))+
  geom_boxplot(aes(ymin=l.ci, lower=l.se, 
                   middle = predicted, upper = u.se, ymax=u.ci, width = 0.5), stat="identity", fill="lightgray", lwd=0.7)+
  geom_jitter(data = intensity, aes(x= trait, y = z, size = 1/var.z, col = id),width=0.15, alpha=2/3)+
  labs(y="Fisher's z", x = "")+geom_hline(yintercept=0, linetype="longdash")+
  theme_bw() +   theme(text = element_text(size=22),axis.text.x = element_text(angle=45, vjust=0.5),
                       plot.background = element_blank(),
                       panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank() )+
  theme(panel.border= element_blank())+
  theme(axis.line.x = element_line(color="black", size = 0.4),
        axis.line.y = element_line(color="black", size = 0.4))+
  scale_color_manual(values=c("green","royalblue3","cyan", "red","magenta","gray40","black","sienna","purple","pink","green4","blue","gold"))+
  guides(colour = guide_legend(override.aes = list(size=4)))

######### END ####################################
