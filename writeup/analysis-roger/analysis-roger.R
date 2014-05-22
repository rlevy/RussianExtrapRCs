library(languageR)
library(cplDataAnalysis)

# Load dat

mainpath = "~/ling/RussianExtrapRCs/writeup"
setwd(mainpath)
files = dir(paste(mainpath,"/Results",sep=""), pattern = "dat", full.names=TRUE, recursive=FALSE)
source('SPRfnx.R')

dat = {}

for (f in files){
        print(f)
        dat = rbind(dat,read.table(f,header=F,fileEncoding="MACCYRILLIC"))
}



# Label dat, remove fillers, convert subjs, items to factors

colnames(dat) = c("subj","expt","item","cond","pos","word","accuracy","rt")
dat$subj <- as.factor(paste("S",dat$subj))
dat$item <- as.factor(paste("I",dat$item))
dat$filler <- dat$expt != "ExtrapRus"
#dat <- droplevels(subset(dat, ! subj %in% "S 19", drop=T))

dat$structure <- factor(sapply(dat$cond, function(x) {
	if(grepl("^NP.*",as.character(x))) {
		"NP"
	} else if (grepl("^VP.*",as.character(x))){
		"VP"
	} else {
		NA	
	}
      }))

dat$locality <- factor(sapply(dat$cond, function(x) {
	if(grepl("^.*_loc$",as.character(x))) {
		"loc"
	} else if (grepl("^.*_nonloc$",as.character(x))){
		"nonloc"
	} else {
		NA	
	}
      }))

acc <- with(subset(dat,pos %in% "?"),tapply(accuracy,list(subj,expt,item),function(x) mean(as.numeric(as.character(x)))))
dat$correct <- with(dat,acc[cbind(subj,expt,item)])

dat <- subset(dat, ! (pos %in% "?" | expt %in% "practice"))
dat$wordlen <- nchar(as.character(dat$word))
dat$region <- factor(paste("R",dat$pos,sep=""))
dat$structureN <- ifelse(dat$structure %in% "NP",-0.5,0.5)
dat$localityN <- ifelse(dat$locality %in% "loc",-0.5,0.5)

regions.to.analyze <- paste("R",0:8,sep="")

## show items
words <- with(droplevels(subset(dat,!filler)),tapply(as.character(word),list(item,cond,pos),function(x) x[1]))
sentences <- apply(words,c(1,2),function(x) {
  tmp <- x[! is.na(x)]
  paste(tmp,collapse=" ")
})

words2 <- with(droplevels(subset(dat,!filler)),tapply(as.character(word),list(factor(paste(item,cond)),pos),function(x) x[1]))
sentences2 <- apply(words2,c(1),function(x) {
  tmp <- x[! is.na(x)]
  paste(tmp,collapse=" ")
})


pdf("rogers-results-condition-specific-outlier-removal.pdf",height=4,width=8)
dat.inliers <- analyze.spr(dat,factors=c("structure","locality"),region.list=regions.to.analyze,use.res=F,analyze.correct=F,col=c("green","magenta","green","magenta"),pch=c(21,21,23,23),lty=c(1,1,2,2),outlier.factors=c("structure","locality"))
dev.off()
system.time(m <- lmer(rt ~ structureN * localityN + (structureN*localityN|subj) + (structureN*localityN|item),subset(dat.inliers$data,region %in% "R4")))
system.time(m <- lmer(rt ~ structureN * localityN + (structureN*localityN|subj)+ (structureN+localityN|item) + (0+structureN:localityN|item),subset(dat.inliers$data,region %in% "R4")))

pdf("rogers-results-single-criterion-outlier-removal.pdf",height=4,width=8)
dat.inliers <- analyze.spr(dat,factors=c("structure","locality"),region.list=regions.to.analyze,use.res=F,analyze.correct=F,col=c("green","magenta","green","magenta"),pch=c(21,21,23,23),lty=c(1,1,2,2),outlier.factors=c())
dev.off()

system.time(m <- lmer(rt ~ structureN * localityN + (structureN*localityN|subj) + (structureN*localityN|item),subset(dat.inliers$data,region %in% "R4")))
system.time(m <- lmer(rt ~ structureN * localityN + (structureN+localityN|subj) + (0+structureN:localityN|subj)+ (structureN*localityN|item),subset(dat.inliers$data,region %in% "R5")))

## check that R5
head(subset(dat,region %in% "R5" & ! filler))
head(subset(dat,region %in% "R4" & ! filler))