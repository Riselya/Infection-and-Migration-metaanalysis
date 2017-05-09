

## R script for analysis in the article "How does infection alter animal migration? A meta-analysis across experimental and observational studies"

##Have tried to annotate it do it makes sense


### Updated 26/4/17


library(metafor)
library(ggplot2)
library(tidyr)
library(dplyr)
library(glmulti)


##import data file "Data_table2_data.csv" (Data table 2) as 'basic'

basic <- read.csv("C:/Users/arisely/Dropbox/PhD/Lit Review/R/LitReviewR/Meta-analysis/Working datasheets/Risely_etal_data.csv")

####################### METAANALYSIS FOR INFECTION INTENSITY

### "Basic" is all data combined for infection status and intensity

###explanations for all columns

names(basic)

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
#[13] "latin"               ##host latin name
#[14] "taxa"                ##host kingdon   
#[15] "Order"               ##host order
#[16] "Family"              ##host family   
#[17] "strain"              ##parasite strain   
#[18] "parasite.taxa"       ##parasite type   
#[19] "life.history.measured"  ##list history of host when sampled
#[20] "setting"             ## host sampled in field or lab 
#[21] "subset"              ##if analyses were split by age/sex, which grouping was used 
#[22] "Infection.type"      ## single or multiple infection measured 
#[23] "migratory.leg2"      ##whether a migratory leg occured before infection and sampling 
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
#[39] "function."           ## the compute.es function used to calculate Fisher's z, Hedges' g and their variances 
#[40] "Effect.direction"    ##effect direction (negative/positive)

unique(basic$id) ##44 studies
str(basic)

##recatgorise some variables
basic$id<-factor(basic$id)
basic$rowid<-factor(basic$rowid)
basic$ss<-as.character(basic$ss)
basic$ss<-as.numeric(basic$ss)
basic$ss_infected<-as.character(basic$ss_infected)
basic$ss_infected<-as.numeric(basic$ss_infected)
basic$ss_healthy<-as.character(basic$ss_healthy)
basic$ss_healthy<-as.numeric(basic$ss_healthy)
basic$sig0.05<-factor(basic$sig0.05)

##########################METAANALYSIS ON EFFECT OF INFECTION STATUS ON PEFORMANCE

##subset observations on infection status

effect<-subset(basic, infection.measure=="Infection status") ##subset observations on infection status

unique(effect$id) ##35 studies included
unique(effect$measure) ##20 measures
table(effect$Effect.direction)
table(effect$sig0.05)

###order categories

effect$trait<-factor(effect$trait, levels=c("Body stores","Refuelling","Movement","Phenology","Survival"))
effect$life.history.measured<-factor(effect$life.history.measured, levels=c("Migration","Breeding","Non-breeding","Lab"))
effect$parasite.taxa<-factor(effect$parasite.taxa, levels=c("Protozoa","Virus","Mites","Helminth","Multiple"))
effect$Order<-factor(effect$Order, levels=c("Passeriformes","Anseriformes","Salmoniformes","Lepidoptera","Falconiformes",
                                            "Anguilliformes","Clupeiformes","Coraciiformes","Charadriiformes"))

effect<-effect[order(effect$trait, effect$g),]



forest.default(effect$g, effect$var.g1, slab = effect$id)


##change names of cohens d and variance to standardized names 

effect$yi<-effect$g
effect$vi<-effect$var.g ##variance not capped
effect$vi1<-effect$var.g1 ##capped to 0.01


####exclude data points that use same animals and same trait as non-independent

effect_independent<-subset(effect, rowid!=27 & rowid!= 28  & 
                             rowid!= 10 & rowid!= 11 & rowid!= 12 & rowid!= 61 & 
                             rowid!= 84 & rowid!= 106 & rowid!= 64 & rowid!= 66 & 
                             rowid!= 68 & rowid!= 160 & rowid!= 158 & rowid!= 156 & rowid!= 51 )

##repeat
#effect_independent1<-subset(effect, rowid!=26 & rowid!= 28 & 
#                            rowid!= 9 & rowid!= 11 & rowid!= 12 & rowid!= 60 & 
#                        rowid!= 83 & rowid!= 105 & rowid!= 63 & rowid!= 65 & 
#                        rowid!= 67 & rowid!= 159 & rowid!= 157 & rowid!= 155 & rowid!= 50 )


effect<-effect_independent


##forest plots by trait


forest.default(effect$g, effect$vi1,  cex=1, xlab="Hedges g", col=effect$trait)
forest.default(effect$g, effect$vi1,  cex=1, xlab="Hedges g", col=effect$sig0.05)
forest.default(effect$g, effect$vi1)



##Kendall rank correlation test

cor.test(effect$yi, effect$vi, method="kendall", alternative="less")
coef(lm(yi ~ vi, data = effect))


##supplementary fig

ggplot(data=effect, aes(x=vi, y = yi))+geom_point(size = 3)+ labs(x="Variance in Hedges' g", y = "Hedges' g")+
  geom_abline(intercept= -0.0562, slope = -2.58)  +
  theme_bw() + 
  theme(text = element_text(size=25),plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank() )+
  theme(panel.border= element_blank())+
  theme(axis.text=element_text(size = 14))+
  theme(axis.line.x = element_line(color="black", size = 0.5),
        axis.line.y = element_line(color="black", size = 0.5))



############################################ MODEL FITTING AND SELECTION ##################################################################################



##explore which variables have an impact on effect size by comparing all models 


rma.glmulti <- function(formula, data, ...) {
  rma.mv(as.formula(paste(deparse(formula))), vi1, data=data, method="ML", ...)
}


res <- glmulti(yi ~ trait+ parasite.taxa+ Order +life.history.measured+type+migratory.leg2,random = ~ 1 | id, data=effect,
               level=1, fitfunction=rma.glmulti, crit="aicc",  maxsize= 2)


tmp <- weightable(res)
tmp

#write.csv(tmp, "aic.csv")
##Phylogeny important


##variable importance
plot(res, type="s")
plot(res, type="w")
print(res)

##best model
summary(res@objects[[1]])
funnel(res@objects[[1]], yaxis="wi")

summary(res@objects[[2]])
funnel(res@objects[[2]], yaxis="wi")

summary(res@objects[[3]])
funnel(res@objects[[3]], yaxis="wi")




#####################best model####################################################################



##null model

null.model<-rma.mv(yi, vi1,random = ~ 1 | id, data=effect, method="REML")
summary(null.model)
funnel(null.model, xlab = "Hedges' g", cex.lab = 1.5, cex.axis=1.5)


#I2 = 100 * (Q - df)/Q
#==52%

#######best model - include trait with host order as a covariable (see Methods section)
##capped variance (vi1)

best.model1<-rma.mv(yi, vi1, mods=~ trait + Order ,random = ~ 1 | id, data=effect)
summary(best.model1)

#I2 = 40%

funnel(best.model1)
funnel(best.model1, yaxis="wi") ##weighting of points


##make table of best model (table 3a)

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
     ylab = "Pearson residuals", 
     cex.lab = 1.5)
abline(h = 0, lty = 2)


##some outliers. 
hist(E1, breaks=10)

##save fitted and residuals and hat values

effect$resids<-resid(best.model1, type ="pearson")
effect$predicted<-fitted(best.model1)
effect$leverage<-hatvalues(best.model1)

mean(effect$leverage)
plot(effect$leverage)

subset(effect, resids < -1)

##three outliers in residuals are from the van gils paper on avian influenza on phenology, movement and refuelling
##exclude these and rerun from line 205 - doesnt make much difference although reduces model residual heterogeneity
#effect<-subset(effect, resids > -1)

##PLOT ESTIMATES (Figure 3)

best.model1<-rma.mv(yi, vi1, mods=~ trait+Order,random = ~ 1 | id, data=effect)
summary(best.model1)
trait.pred<-predict(best.model1, newmods=rbind(c(0,0,0,0,0,0,0,0,0,0,0,0),c(1,0,0,0,0,0,0,0,0,0,0,0),c(0,1,0,0,0,0,0,0,0,0,0,0),
                                               c(0,0,1,0,0,0,0,0,0,0,0,0),c(0,0,0,1,0,0,0,0,0,0,0,0)))


trait<-data.frame(trait.pred$pred)
trait$se<-trait.pred$se
trait$l.ci<-trait.pred$ci.lb
trait$u.ci<-trait.pred$ci.ub


trait$level<-factor(c("Body stores","Refuelling","Movement","Phenology","Survival"))
trait$level<-factor(trait$level, levels=c("Body stores","Refuelling","Movement","Phenology","Survival"))
colnames(trait)[1]<-"predicted"

predictedvalues<-trait

predictedvalues$l.se<-predictedvalues$predicted - predictedvalues$se
predictedvalues$u.se<-predictedvalues$predicted + predictedvalues$se


##boxplot

##overlayed with raw data with proportional weights

###################fig 3a 

ggplot(data=predictedvalues, aes(x=level, y = predicted, fill = Order))+
  geom_boxplot(aes(ymin=l.ci, lower=l.se, 
                   middle = predicted, upper = u.se, ymax=u.ci, width = 0.5), stat="identity", fill="lightgray", lwd=0.7)+
  geom_jitter(data = effect, aes(x= trait, y = yi, size = 1/vi1, col = Order),width=0.15, alpha=1/2)+
  labs(y="Predicted Hedges' g", x ="")+geom_hline(yintercept=0, linetype="longdash")+
  theme_bw() +   theme(text = element_text(size=22),axis.text.x = element_text(angle=45, vjust=0.5),
                       plot.background = element_blank(),
                       panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank() )+
  theme(panel.border= element_blank())+
  theme(axis.line.x = element_line(color="black", size = 0.4),
        axis.line.y = element_line(color="black", size = 0.4))+
  scale_y_continuous(breaks=seq(-2,0.5, 0.5))+scale_color_brewer(type = "qual", palette = "Set1")+
  guides(colour = guide_legend(override.aes = list(size=6)))


####################fig 3b

ggplot(data=predictedvalues, aes(x=level, y = predicted))+
  geom_boxplot(aes(ymin=l.ci, lower=l.se, 
                   middle = predicted, upper = u.se, ymax=u.ci, width = 0.5), stat="identity", fill="lightgray", lwd=0.7)+
  geom_jitter(data = effect, aes(x= trait, y = yi, size = 1/vi1, col = parasite.taxa),width=0.1, alpha=1/2)+
  labs(y="Hedges' g", x = "")+geom_hline(yintercept=0, linetype="longdash")+
  theme_bw() +   theme(text = element_text(size=22),axis.text.x = element_text(angle=45, vjust=0.5),
                       plot.background = element_blank(),
                       panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank() )+
  theme(panel.border= element_blank())+
  theme(axis.line.x = element_line(color="black", size = 0.4),
        axis.line.y = element_line(color="black", size = 0.4))+
  scale_y_continuous(breaks=seq(-2,0.5, 0.5))+scale_color_manual(values=c("red3","royalblue3","darkgreen", "orange","magenta","cyan3"))+
  guides(colour = guide_legend(override.aes = list(size=6)))


#############################################################################################################


#########################modelling traits seperately


survival<-subset(effect, trait=="Survival")
Bodystores<-subset(effect, trait=="Body stores")
phenology<-subset(effect, trait=="Phenology")
movement<-subset(effect, trait=="Movement")
refuelling<-subset(effect, trait=="Refuelling")


null1<-rma.mv(yi, vi1, random = ~ 1 | id, data=Bodystores) ##null model
summary(null1)

#I2 = 100 * (28.7-16)/28.7

#I2 = 38

null2<-rma.mv(yi, vi1, random = ~ 1 | id, data=refuelling) ##null model
summary(null2)

#I2 = 44%

null3<-rma.mv(yi, vi1, random = ~ 1 | id, data=movement) ##null model
summary(null3)

#I2 = 65%


null4<-rma.mv(yi, vi1, random = ~ 1 | id, data=phenology) ##null model
summary(null4)

#I2 = 58%

null5<-rma.mv(yi, vi1, random = ~ 1 | id, data=survival) ##null model
summary(null5)


# I2 = 18%

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


ggplot(data=predictedvalues, aes(x=level, y = predicted, fill = parasite.taxa))+
  geom_boxplot(aes(ymin=l.ci, lower=l.se, 
                   middle = predicted, upper = u.se, ymax=u.ci, width = 0.5), stat="identity", fill="lightgray", lwd=0.7)+
  geom_jitter(data = effect, aes(x= trait, y = yi, size = 1/vi1, col = parasite.taxa),width=0.15, alpha=1/2)+
  labs(y="Predicted Hedges' g", x ="")+geom_hline(yintercept=0, linetype="longdash")+
  theme_bw() +   theme(text = element_text(size=22),axis.text.x = element_text(angle=45, vjust=0.5),
                       plot.background = element_blank(),
                       panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank() )+
  theme(panel.border= element_blank())+
  theme(axis.line.x = element_line(color="black", size = 0.4),
        axis.line.y = element_line(color="black", size = 0.4))+
  scale_y_continuous(breaks=seq(-2,0.5, 0.5))+scale_color_manual(values=c("red","blue","cyan3", "gold","magenta"))+
  guides(colour = guide_legend(override.aes = list(size=6)))




#model with just trait

best.model2<-rma.mv(yi, vi1, mods=~ trait,random = ~ 1 | id, data=effect)
summary(best.model2)
trait.pred<-predict(best.model2, newmods=rbind(c(0,0,0,0),c(1,0,0,0),c(0,1,0,0),
                                               c(0,0,1,0),c(0,0,0,1)))

trait<-data.frame(trait.pred$pred)
trait$se<-trait.pred$se
trait$l.ci<-trait.pred$ci.lb
trait$u.ci<-trait.pred$ci.ub


trait$level<-factor(c("Body stores","Refuelling","Movement","Phenology","Survival"))
trait$level<-factor(trait$level, levels=c("Body stores","Refuelling","Movement","Phenology","Survival"))
colnames(trait)[1]<-"predicted"

predictedvalues<-trait

predictedvalues$l.se<-predictedvalues$predicted - predictedvalues$se
predictedvalues$u.se<-predictedvalues$predicted + predictedvalues$se


##boxplot

##overlayed with raw data with proportional weights

ggplot(data=predictedvalues, aes(x=level, y = predicted, fill = parasite.taxa))+
  geom_boxplot(aes(ymin=l.ci, lower=l.se, 
                   middle = predicted, upper = u.se, ymax=u.ci, width = 0.5), stat="identity", fill="lightgray", lwd=0.7)+
  geom_jitter(data = effect, aes(x= trait, y = yi, size = 1/vi1, col = parasite.taxa),width=0.15, alpha=1/2)+
  labs(y="Predicted Hedges' g", x = "")+geom_hline(yintercept=0, linetype="longdash")+
  theme_bw() +   theme(text = element_text(size=22),axis.text.x = element_text(angle=45, vjust=0.5),
                       plot.background = element_blank(),
                       panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank() )+
  theme(panel.border= element_blank())+
  theme(axis.line.x = element_line(color="black", size = 0.4),
        axis.line.y = element_line(color="black", size = 0.4))+
  scale_y_continuous(breaks=seq(-2,0.5, 0.5))+scale_color_manual(values=c("red3","royalblue3","darkgreen", "orange","yellow"))+
  guides(colour = guide_legend(override.aes = list(size=6)))

##very similar to that of plot where all traits are modelled seperately.



###########################################SENSITIVITY ANALYSES##################################

###use uncapped varience instead of capped varience (replace vi1 with vi in all models) and rerun

##exlude insects in case this biases biases models and rerun

##exlude outliers and rerun

###Does not change models overall



##################################################################################################################################



####################################           INTENSITY  METANALYSIS              ######################################


intensity<-subset(basic, infection.measure=="Intensity")


table(intensity$sig0.05)

###exclude points that use the same animals
#intensity<-subset(intensity, rowid != 4 & rowid!=5 & rowid!=134 & rowid!=22)

intensity<-intensity[order(intensity$trait, intensity$z),]
forest.default(intensity$z, intensity$var.z)
forest.default(intensity$z, intensity$var.z, col = intensity$sig0.05)
forest.default(intensity$z, intensity$var.z, col = intensity$trait)


##exclude moller as as 9 datapoints so might bias results because so many 
#intensity<-subset(intensity, id !="62")

##order categories

intensity$trait<-factor(intensity$trait, levels=c("Body stores","Movement","Phenology","Survival"))

intensity$parasite.taxa<-factor(intensity$parasite.taxa, levels=c("Protozoa","Virus","Mites","Helminth","Myxospora"))
intensity$Order<-factor(intensity$Order, levels=c("Passeriformes","Anseriformes","Salmoniformes","Lepidoptera","Falconiformes",
                                                  "Anguilliformes","Clupeiformes","Coraciiformes","Charadriiformes"))


###univariate models with each of the six explanatory variables

model.intensity<-rma.mv(z, var.z, mods = ~ trait, random = ~1|id,  data=intensity)
model.intensity1<-rma.mv(z, var.z, mods = ~ Order, random = ~1|id,  data=intensity)
model.intensity2<-rma.mv(z, var.z, mods = ~ parasite.taxa, random = ~1|id,  data=intensity)
model.intensity3<-rma.mv(z, var.z, mods = ~ life.history.measured, random = ~1|id,  data=intensity)
model.intensity4<-rma.mv(z, var.z, mods = ~ migratory.leg2, random = ~1|id,  data=intensity)
model.intensity5<-rma.mv(z, var.z, mods = ~ type, random = ~1|id,  data=intensity)
model.intensity6<-rma.mv(z, var.z, random = ~1|id,  data=intensity)

aicc(model.intensity)
aicc(model.intensity1)
aicc(model.intensity2)
aicc(model.intensity3)
aicc(model.intensity4)
aicc(model.intensity5)
aicc(model.intensity6)
funnel(model.intensity6)
summary(model.intensity6)

##best model includes just trait (AICC = 42.4)

summary(model.intensity)

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
  guides(colour = guide_legend(override.aes = list(size=6)))




######################model intensity seperately

survival<-subset(intensity, trait=="Survival")
Bodystores<-subset(intensity, trait=="Body stores")
phenology<-subset(intensity, trait=="Phenology")
movement<-subset(intensity, trait=="Movement")



null1<-rma.mv(z, var.z, random = ~ 1 | id, data=Bodystores) ##null model
summary(null1)

#I2 = 100 * (28.7-16)/28.7

#I2 = 38

#I2 = 44%

null3<-rma.mv(z, var.z, random = ~ 1 | id, data=movement) ##null model
summary(null3)

#I2 = 65%


null4<-rma.mv(z, var.z, random = ~ 1 | id, data=phenology) ##null model
summary(null4)

#I2 = 58%

null5<-rma.mv(z, var.z, random = ~ 1 | id, data=survival) ##null model
summary(null5)


# I2 = 18%

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
  guides(colour = guide_legend(override.aes = list(size=6)))



##coloured by host phylogeny


ggplot(data=predictedvalues, aes(x=level, y = predicted))+
  geom_boxplot(aes(ymin=l.ci, lower=l.se, 
                   middle = predicted, upper = u.se, ymax=u.ci, width = 0.5), stat="identity", fill="lightgray", lwd=0.7)+
  geom_jitter(data = intensity, aes(x= trait, y = z, size = 1/var.z, col = Order),width=0.15, alpha=1/2)+
  labs(y="Fisher's z", x = "")+geom_hline(yintercept=0, linetype="longdash")+
  theme_bw() +   theme(text = element_text(size=22),axis.text.x = element_text(angle=45, vjust=0.5),
                       plot.background = element_blank(),
                       panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank() )+
  theme(panel.border= element_blank())+
  theme(axis.line.x = element_line(color="black", size = 0.4),
        axis.line.y = element_line(color="black", size = 0.4))+
  scale_color_brewer(type = "qual", palette = "Set1")+
  guides(colour = guide_legend(override.aes = list(size=6)))



