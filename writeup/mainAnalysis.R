# Load libraries

library(xtable)
library(sciplot)
library(languageR)
library(gdata)
library(ez)

# Load data

mainpath = "~/Dropbox/Current Work/RussianExtrapRCs/writeup"
setwd(mainpath)
files = dir(paste(mainpath,"/Results",sep=""), pattern = "dat", full.names=TRUE, recursive=FALSE)
source('SPRfnx.R')

data = {}

for (f in files){
	print(f)
	data = rbind(data,read.table(f,header=F,fileEncoding="latin1"))
}

# Label data, remove fillers, convert subjs, items to factors

colnames(data) = c("subj","exp","item","cond","pos","word","region","value")
data$subj = as.factor(data$subj)
data$item = as.factor(data$item)
data = subset(data,exp=="ExtrapRus")

# Create factors for structure and locality 

condlabels = strsplit(as.character(data$cond),split="_")
data$structure = as.factor(unlist(condlabels)[2*(1:length(data$cond))-1])
data$locality = as.factor(unlist(condlabels)[2*(1:length(data$cond))])

## Accuracy analysis ##
#######################

acc.data = drop.levels(subset(data,pos=="?" & exp=='ExtrapRus'),reorder=FALSE)
contrasts(data$structure) = contr.sum(levels(data$structure))/2
contrasts(data$locality) = contr.sum(levels(data$locality))/2
acc.data$region = as.integer(as.character(acc.data$region))
acc.model = glmer(region~structure*locality+(1+structure+locality|subj)+(1+ structure+locality |item),data=acc.data,control=glmerControl(optimizer="bobyqa"),family="binomial")

## RT analysis below here ##
############################

data = drop.levels(subset(data,pos != "?"),reorder=FALSE)    ## Why does drop.levels reset contrasts associated with factors? 
contrasts(data$structure) = contr.sum(levels(data$structure))/2
contrasts(data$locality) = contr.sum(levels(data$locality))/2
data$structureN = ifelse(data$structure %in% "NP",-0.5,0.5)
data$localityN = ifelse(data$locality %in% "loc",-0.5,0.5)
data$pos = as.integer(as.character(data$pos))
data$region = as.factor(data$pos)

# check Latin Square is correct; check (8/1)
xtabs(~subj+cond,data=subset(data,pos==0))

# Outlier removal
# Removed observations greater than 5000, less than 100; adopting thresholds from JML russian paper

N = nrow(data)
data = subset(data,value < 5000 & value > 100)
data = zscoreRegion(data,cutoff=3)
print(nrow(data))
trim.percentage = round((1-nrow(data)/N)*100,digits=1)
print(trim.percentage)

data$pos = data$pos + 1

# LME modeling; maximal effects structured adopted throughout
# Note, LMEs were used because of significant imbalances introduced in # of data points per subj
# given the outlier removal procedure. 

lmeModels = list()

for (curPos in 1:9) {
	lmeModels[[curPos]] = lmer(value~structureN*localityN+(1+ structureN*localityN|subj)+(1+ structureN*localityN|item),control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=25000)),data=subset(data,pos==curPos))
}


## Resolve critical interaction in region 6 using nested contrasts.

nestedModel = lmer(value~structure/localityN+(1+ structureN*localityN|subj)+(1+ structureN*localityN|item),control=lmerControl(optimizer="bobyqa"),data=subset(data,pos==6))

## Regions 6, 7 and 9 have significant or marginal effects
## Here summaries of fixed effects are created for easy referencing in the TeX document

Vcov = vcov(lmeModels[[6]], useScale = FALSE)
betas <- round(fixef(lmeModels[[6]]),digits = 0)
se <-  round(sqrt(diag(Vcov)),digits=0)
tval = betas / se
pvals = 2 * pnorm(abs(tval), lower.tail = FALSE)
region6.summary = cbind(betas,se,tval,pvals)

Vcov <- vcov(lmeModels[[7]], useScale = FALSE)
betas <- round(fixef(lmeModels[[7]]),digits = 0)
se <-  round(sqrt(diag(Vcov)),digits=0)
tval <- betas / se
pvals = 2 * pnorm(abs(tval), lower.tail = FALSE)
region7.summary = cbind(betas,se,tval,pvals)

Vcov <- vcov(lmeModels[[9]], useScale = FALSE)
betas <- round(fixef(lmeModels[[9]]),digits = 0)
se <-  round(sqrt(diag(Vcov)),digits=0)
tval <- betas / se
pvals = 2 * pnorm(abs(tval), lower.tail = FALSE)
region9.summary = cbind(betas,se,tval,pvals)

Vcov <- vcov(nestedModel, useScale = FALSE)
betas <- round(fixef(nestedModel),digits = 0)
se <-  round(sqrt(diag(Vcov)),digits=0)
tval <- betas / se
pvals = 2 * pnorm(abs(tval), lower.tail = FALSE)
nestedModel.summary = cbind(betas,se,tval,pvals)