library(lme4)
library(reshape)
library(languageR)

getSpill <- function(data,nspills=1) {
	spills <- matrix(NA,ncol=nspills,nrow=nrow(data))
	for (j in 1:nspills) {
		for (i in 1:nrow(data)) {
			if (data$pos[i] > j-1) {
				spills[i,j] <- data$value[i-j]
			} 
		}
	}
	data$spill1 = spills[,1]
	return(data)
}

zscore <- function(data,cutoff=3) {
	means <- tapply(data$value, list(data$region, data$cond), mean)
	sds <- tapply(data$value, list(data$region, data$cond), sd)
	vector.means <- {}
	vector.sds <- {}
	for (i in 1:nrow(data)) {
		vector.means <- c(vector.means,means[data$region[i],data$cond[i]])
		vector.sds <- c(vector.sds,sds[data$region[i],data$cond[i]])
	}
	zscore <- (data$value-vector.means)/vector.sds
	data <- cbind(data,zscore)
	data <- subset(data, abs(zscore) < cutoff) 
	return(data)
}


zscoreRegion <- function(data,cutoff=3) {
	means <- tapply(data$value, list(data$region), mean)
	sds <- tapply(data$value, list(data$region), sd)
	zscore <- (data$value-means[data$region])/sds[data$region]
	data <- cbind(data,zscore)
	data <- subset(data, abs(zscore) < cutoff )
	return(data)
}


fitResids <- function(data,spill=T) {
	rownames(data) = c(1:nrow(data))
	if (spill==F) {
		model <- lm(value~subj+pos+wordLength,data=data)
		data$resid <- residuals(model)
	} else {
		data$resid = rep(NA,nrow(data))
		for (reg in levels(data$pos)) {
			 if (reg != '0') {
						 curSet <- subset(data,pos==reg)
						 model <- lm(value~subj+wordLength+spill1,data=curSet)
						 data$resid[as.double(rownames(curSet))] <- residuals(model)
				}
		}
	}
	return(data)
}

p.print = function(pval) {
	if (pval < 0.001) {return(0.001)}
	else if (pval < 0.01) {return(0.01)}
	else if (pval < 0.05) {return(0.05)}
	else if (pval < 0.1) {return(0.1)}
	else return(pval)
}


