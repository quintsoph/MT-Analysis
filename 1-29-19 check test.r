#fix Pvalue issue
setwd("/home/quints1/Deseqstuff")
genesraw<- read.table("normalizedcnt.txt", header=TRUE, row.names=1)
newmatrix=NULL
for (i in colnames(genesraw))
{
	summation<-sum(genesraw[,i])
	newmatrix<- cbind(newmatrix, summation)
	
	}
colnames(newmatrix)=colnames(genesraw)
log2sub<- log2(newmatrix)
tlog2sub<- t(log2sub)
#vstsum<- vst(newmatrix, blind=FALSE)

#adrenal test ACTG1
library("stringr")
	setwd("/home/quints1/Deseqstuff")
	totaltissue<- read.table("allL1HSmatrix.txt", header=TRUE, row.names=1, sep="\t")
	matchcut_new<- readRDS("patientgene.rds")
	ttissue<- t(totaltissue)
	library("stringr")
	colnames(ttissue)=str_replace_all(colnames(ttissue), "-", ".")
	t_TEcntmatrix = t(ttissue)
	 t_cntmatrix = t(matchcut_new)
	m <- merge(t_cntmatrix, t_TEcntmatrix, by.x="row.names", by.y="row.names")
	new_patients = m[,1]
	new_genenames = colnames(m[,-1])
	new_cntmatrix = as.matrix(t(m[,-1]))
	finalpatlist<- data.frame(new_patients)
	rownames(finalpatlist)=finalpatlist[,1]
	setwd("/home/quints1/Deseqstuff")
	L1HSraw<- read.table("L1HS.VST.cnts.txt", header=TRUE, row.names=1)
	genesraw<- read.table("VSTcnt.txt", header=TRUE, row.names=1)
	tgenesraw<- t(genesraw)
	setwd("/home/quints1/Correlations/L1HSpatientedits")
	adrenalpats<- read.table("Adrenaledits.txt", header=TRUE, sep="\t")
	rownames(adrenalpats)=str_replace_all(adrenalpats[,1], "-", ".")
	adrenalgenes<- merge(adrenalpats, tgenesraw, by.x="row.names", by.y="row.names", drop=FALSE)
	adrenalgenesfinal<- adrenalgenes[,c(4:56207)]
	rownames(adrenalgenesfinal)=adrenalgenes[,1]
		#read in the covariate
		setwd("/home/quints1/Correlations/adrenal")
		adrenalcov<- read.table("adrenaleditcov.txt", header=TRUE, row.names=1, sep="\t")
	adrenalall<- merge(adrenalgenesfinal, adrenalcov, by.x="row.names", by.y="row.names")
	finaladrenalcov<- data.frame(adrenalall[,c(56206:56223)])	
	rownames(finaladrenalcov)=adrenalall[,1]
	finaladrenalL1HS<- data.frame(adrenalall[,c(56205)])
	rownames(finaladrenalL1HS)= adrenalall[,1]
	colnames(finaladrenalL1HS)= c("L1HS")
	finaladrenalgenes<- data.frame(adrenalall[,c(3:56204)])
	rownames(finaladrenalgenes)=adrenalall[,1]
adrenallogpat<-  merge(finaladrenalgenes, tlog2sub, by.x="row.names", by.y="row.names")
rownames(adrenallogpat)=adrenallogpat[,1]
finallog<- data.frame(adrenallogpat[,56204])
rownames(finallog)= adrenallogpat[,1]
adrenalACTG1<- data.frame(adrenallogpat[,c("ACTG1")])
rownames(adrenalACTG1)=adrenallogpat[,1]
library("ggplot2")
	#graphforadrenal
	mydata<- data.frame(merge(finallog,adrenalACTG1, by.x="row.names", by.y="row.names"))
	mydata<- mydata[-109,]
	myplotbottom<- ggplot(mydata, aes(mydata[,3], mydata[,2])) + geom_point(color="black") + geom_smooth(method="lm", se=FALSE, color="dark red") + theme_classic() +
					labs(x="ACTG1 Gene Expression", y="logsum Expression", title="Adrenal Test")
	#graph with L1HS and genes
	mydata<- data.frame(merge(finaladrenalL1HS,adrenalACTG1, by.x="row.names", by.y="row.names"))
	mydata<- mydata[-109,]
	myplot<- ggplot(mydata, aes(mydata[,3], mydata[,2])) + geom_point(color="black") + geom_smooth(method="lm", se=FALSE, color="dark red") + theme_classic() +
					labs(x="ACTG1 Gene Expression", y="L1HS Expression", title="Adrenal Test2")
	#graph with L1HS and sum
	mydata<- data.frame(merge(finallog, finaladrenalL1HS, by.x="row.names", by.y="row.names"))
	mydata<- mydata[-109,]
	myplot<- ggplot(mydata, aes(mydata[,3], mydata[,2])) + geom_point(color="black") + geom_smooth(method="lm", se=FALSE, color="dark red") + theme_classic() +
					labs(x="ACTG1 Gene Expression", y="L1HS Expression", title="Adrenal Test3")	
	#Breast
setwd("/home/quints1/Correlations/L1HSpatientedits")
	breastpats<- read.table("Breastedits.txt", header=TRUE, sep="\t")
	rownames(breastpats)=str_replace_all(breastpats[,1], "-", ".")
	breastgenes<- merge(breastpats, tgenesraw, by.x="row.names", by.y="row.names", drop=FALSE)
	breastgenesfinal<- breastgenes[,c(4:56207)]
	rownames(breastgenesfinal)=breastgenes[,1]
		#read in the covariate
		setwd("/home/quints1/Correlations/breast")
		breastcov<- read.table("breasteditcov.txt", header=TRUE, row.names=1, sep="\t")
	breastall<- merge(breastgenesfinal, breastcov, by.x="row.names", by.y="row.names")
	finalbreastcov<- data.frame(breastall[,c(56206:56238)])	
	rownames(finalbreastcov)=breastall[,1]
	finalbreastL1HS<- data.frame(breastall[,c(56205)])
	rownames(finalbreastL1HS)= breastall[,1]
	colnames(finalbreastL1HS)= c("L1HS")
	finalbreastgenes<- data.frame(breastall[,c(3:56204)])
	rownames(finalbreastgenes)=breastall[,1]
breastlogpat<-  merge(finalbreastgenes, tlog2sub, by.x="row.names", by.y="row.names")
rownames(breastlogpat)=breastlogpat[,1]
finallog<- data.frame(breastlogpat[,56204])
rownames(finallog)= breastlogpat[,1]
breastACTG1<- data.frame(breastlogpat[,c("ACTG1")])
rownames(breastACTG1)=breastlogpat[,1]
library("ggplot2")
	#graphforadrenal
	mydata<- data.frame(merge(finallog,breastACTG1, by.x="row.names", by.y="row.names"))
	mydata<- mydata[-109,]
	myplotbottom<- ggplot(mydata, aes(mydata[,3], mydata[,2])) + geom_point(color="black") + geom_smooth(method="lm", se=FALSE, color="dark red") + theme_classic() +
					labs(x="ACTG1 Gene Expression", y="logsum Expression", title="Breast Test")
	#graph with L1HS and genes
	mydata<- data.frame(merge(finalbreastL1HS,breastACTG1, by.x="row.names", by.y="row.names"))
	mydata<- mydata[-109,]
	myplot<- ggplot(mydata, aes(mydata[,3], mydata[,2])) + geom_point(color="black") + geom_smooth(method="lm", se=FALSE, color="dark red") + theme_classic() +
					labs(x="ACTG1 Gene Expression", y="L1HS Expression", title="Breast Test2")
	#graph with L1HS and sum
	mydata<- data.frame(merge(finallog, finalbreastL1HS, by.x="row.names", by.y="row.names"))
	mydata<- mydata[-109,]
	myplot<- ggplot(mydata, aes(mydata[,3], mydata[,2])) + geom_point(color="black") + geom_smooth(method="lm", se=FALSE, color="dark red") + theme_classic() +
					labs(x="L1HS Expression", y="log2sum", title="Breast Test3")	
####################################
################################
##########################
##################
#Start to rerun adrenal and breast with total sum

		#breast
	setwd("/home/quints1/Correlations/L1HSpatientedits")
	breastpats<- read.table("Breastedits.txt", header=TRUE, sep="\t")
	rownames(breastpats)=str_replace_all(breastpats[,1], "-", ".")
	breastgenes<- merge(breastpats, tgenesraw, by.x="row.names", by.y="row.names", drop=FALSE)
	breastgenesfinal<- breastgenes[,c(4:56207)]
	rownames(breastgenesfinal)=breastgenes[,1]
		#read in the covariate
		setwd("/home/quints1/Correlations/breast")
		breastcov<- read.table("breasteditcov.txt", header=TRUE, row.names=1, sep="\t")
	breastall<- merge(breastgenesfinal, breastcov, by.x="row.names", by.y="row.names")
	finalbreastcov<- data.frame(breastall[,c(56206:56238)])	
	rownames(finalbreastcov)=breastall[,1]
	finalbreastL1HS<- data.frame(breastall[,c(56205)])
	rownames(finalbreastL1HS)= breastall[,1]
	colnames(finalbreastL1HS)= c("L1HS")
	finalbreastgenes<- data.frame(breastall[,c(3:56204)])
	rownames(finalbreastgenes)=breastall[,1]
	breastlogpat<-  merge(finalbreastgenes, tlog2sub, by.x="row.names", by.y="row.names")
	log2sum<- data.frame(breastlogpat[,56204])
	colnames(log2sum)=c("log2sum")
	rownames(log2sum)=breastlogpat[,1]
	
	library(nlme)
	library(MuMIn)
	library(lmSupport)
	library(qvalue)
	library(matrixStats)
	
	#AICc
	finalAICc=NULL
		for (i in colnames(finalbreastgenes))
		{
		valuepc0 <-lm(finalbreastL1HS[,1] ~ finalbreastgenes[,c(i)] + finalbreastcov[,1] + log2sum[,1])
		valuepc1 <-lm(finalbreastL1HS[,1] ~ finalbreastgenes[,c(i)] + finalbreastcov[,2] + log2sum[,1])
		valuepc2 <-lm(finalbreastL1HS[,1] ~ finalbreastgenes[,c(i)] + finalbreastcov[,3] + log2sum[,1])
		valuein1 <-lm(finalbreastL1HS[,1] ~ finalbreastgenes[,c(i)] + finalbreastcov[,4] + log2sum[,1])
		valuein2 <-lm(finalbreastL1HS[,1] ~ finalbreastgenes[,c(i)] + finalbreastcov[,5] + log2sum[,1])
		valuein3 <-lm(finalbreastL1HS[,1] ~ finalbreastgenes[,c(i)] + finalbreastcov[,6] + log2sum[,1])
		valuein4 <-lm(finalbreastL1HS[,1] ~ finalbreastgenes[,c(i)] + finalbreastcov[,7] + log2sum[,1])
		valuein5 <-lm(finalbreastL1HS[,1] ~ finalbreastgenes[,c(i)] + finalbreastcov[,8] + log2sum[,1])
		valuein6 <-lm(finalbreastL1HS[,1] ~ finalbreastgenes[,c(i)] + finalbreastcov[,9] + log2sum[,1])
		valuein7 <-lm(finalbreastL1HS[,1] ~ finalbreastgenes[,c(i)] + finalbreastcov[,10] + log2sum[,1])
		valuein8 <-lm(finalbreastL1HS[,1] ~ finalbreastgenes[,c(i)] + finalbreastcov[,11] + log2sum[,1])
		valuein9 <-lm(finalbreastL1HS[,1] ~ finalbreastgenes[,c(i)] + finalbreastcov[,12] + log2sum[,1])
		valuein10 <-lm(finalbreastL1HS[,1] ~ finalbreastgenes[,c(i)] + finalbreastcov[,13] + log2sum[,1])
		valuein11 <-lm(finalbreastL1HS[,1] ~ finalbreastgenes[,c(i)] + finalbreastcov[,14] + log2sum[,1])
		valuein12 <-lm(finalbreastL1HS[,1] ~ finalbreastgenes[,c(i)] + finalbreastcov[,15] + log2sum[,1])
		valuein13 <-lm(finalbreastL1HS[,1] ~ finalbreastgenes[,c(i)] + finalbreastcov[,16] + log2sum[,1])
		valuein14 <-lm(finalbreastL1HS[,1] ~ finalbreastgenes[,c(i)] + finalbreastcov[,17] + log2sum[,1])
		valuein15 <-lm(finalbreastL1HS[,1] ~ finalbreastgenes[,c(i)] + finalbreastcov[,18] + log2sum[,1])
		valuein16 <-lm(finalbreastL1HS[,1] ~ finalbreastgenes[,c(i)] + finalbreastcov[,19] + log2sum[,1])
		valuein17 <-lm(finalbreastL1HS[,1] ~ finalbreastgenes[,c(i)] + finalbreastcov[,20] + log2sum[,1])
		valuein18 <-lm(finalbreastL1HS[,1] ~ finalbreastgenes[,c(i)] + finalbreastcov[,21] + log2sum[,1])
		valuein19 <-lm(finalbreastL1HS[,1] ~ finalbreastgenes[,c(i)] + finalbreastcov[,22] + log2sum[,1])
		valuein20 <-lm(finalbreastL1HS[,1] ~ finalbreastgenes[,c(i)] + finalbreastcov[,23] + log2sum[,1])
		valuein21 <-lm(finalbreastL1HS[,1] ~ finalbreastgenes[,c(i)] + finalbreastcov[,24] + log2sum[,1])
		valuein22 <-lm(finalbreastL1HS[,1] ~ finalbreastgenes[,c(i)] + finalbreastcov[,25] + log2sum[,1])
		valuein23 <-lm(finalbreastL1HS[,1] ~ finalbreastgenes[,c(i)] + finalbreastcov[,26] + log2sum[,1])
		valuein24 <-lm(finalbreastL1HS[,1] ~ finalbreastgenes[,c(i)] + finalbreastcov[,27] + log2sum[,1])
		valuein25 <-lm(finalbreastL1HS[,1] ~ finalbreastgenes[,c(i)] + finalbreastcov[,28] + log2sum[,1])
		valuein26 <-lm(finalbreastL1HS[,1] ~ finalbreastgenes[,c(i)] + finalbreastcov[,29] + log2sum[,1])
		valuein27 <-lm(finalbreastL1HS[,1] ~ finalbreastgenes[,c(i)] + finalbreastcov[,30] + log2sum[,1])
		valuein28 <-lm(finalbreastL1HS[,1] ~ finalbreastgenes[,c(i)] + finalbreastcov[,31] + log2sum[,1])
		valuein29 <-lm(finalbreastL1HS[,1] ~ finalbreastgenes[,c(i)] + finalbreastcov[,32] + log2sum[,1])
		valuein30 <-lm(finalbreastL1HS[,1] ~ finalbreastgenes[,c(i)] + finalbreastcov[,33] + log2sum[,1])
		valueall <-lm(finalbreastL1HS[,1] ~ finalbreastgenes[,c(i)] + finalbreastcov[,1] + finalbreastcov[,2] + 
					finalbreastcov[,3] + finalbreastcov[,4] + finalbreastcov[,5] + finalbreastcov[,6] + finalbreastcov[,7] + 
					finalbreastcov[,8] + finalbreastcov[,9] + finalbreastcov[,10] + finalbreastcov[,11] + finalbreastcov[,12] +
					finalbreastcov[,13] + finalbreastcov[,14] + finalbreastcov[,15] + finalbreastcov[,16] + finalbreastcov[,17] + finalbreastcov[,18]+
					finalbreastcov[,19] + finalbreastcov[,20] + finalbreastcov[,21] + finalbreastcov[,22] + finalbreastcov[,23] + finalbreastcov[,24] +
					finalbreastcov[,25] + finalbreastcov[,26]+ finalbreastcov[,27] + finalbreastcov[,28] + finalbreastcov[,29] + finalbreastcov[,30] + 
					finalbreastcov[,31] + finalbreastcov[,32] + finalbreastcov[,33] + log2sum[,1])
		AICc<- model.sel(valuepc0, valuepc1, valuepc2, valuein1, valuein2, valuein3, valuein4, valuein5, valuein6, valuein7, valuein8,
			valuein9, valuein10, valuein11, valuein12, valuein13, valuein14, valuein15, valuein16, valuein17, valuein18, valuein19,
			valuein20, valuein21, valuein22, valuein23, valuein24, valuein25, valuein26, valuein27, valuein28, valuein29, valuein30,
			valueall)
		evaluate<- data.frame(AICc[1,])
		result<- cbind(evaluate, rownames(evaluate))
		name<- paste(i, "topmodel", sep="")
		rownames(result)= name
		finalAICc<- rbind(finalAICc, result)
	}
	
	setwd("/home/quints1/Correlations/breast")
	write.table(finalAICc, "All AICc results breast redo.txt", sep="\t", quote=FALSE)
	
	#final analysis
	pvalues=NULL
	for (i in colnames(finalbreastgenes))
		{
		if(! all(is.na(finalbreastgenes[i]))){
		valuepc0 <-lm(finalbreastL1HS[,1] ~ finalbreastgenes[,c(i)] + finalbreastcov[,1] + log2sum[,1])
		valuepc1 <-lm(finalbreastL1HS[,1] ~ finalbreastgenes[,c(i)] + finalbreastcov[,2] + log2sum[,1])
		valuepc2 <-lm(finalbreastL1HS[,1] ~ finalbreastgenes[,c(i)] + finalbreastcov[,3] + log2sum[,1])
		valuein1 <-lm(finalbreastL1HS[,1] ~ finalbreastgenes[,c(i)] + finalbreastcov[,4] + log2sum[,1])
		valuein2 <-lm(finalbreastL1HS[,1] ~ finalbreastgenes[,c(i)] + finalbreastcov[,5] + log2sum[,1])
		valuein3 <-lm(finalbreastL1HS[,1] ~ finalbreastgenes[,c(i)] + finalbreastcov[,6] + log2sum[,1])
		valuein4 <-lm(finalbreastL1HS[,1] ~ finalbreastgenes[,c(i)] + finalbreastcov[,7] + log2sum[,1])
		valuein5 <-lm(finalbreastL1HS[,1] ~ finalbreastgenes[,c(i)] + finalbreastcov[,8] + log2sum[,1])
		valuein6 <-lm(finalbreastL1HS[,1] ~ finalbreastgenes[,c(i)] + finalbreastcov[,9] + log2sum[,1])
		valuein7 <-lm(finalbreastL1HS[,1] ~ finalbreastgenes[,c(i)] + finalbreastcov[,10] + log2sum[,1])
		valuein8 <-lm(finalbreastL1HS[,1] ~ finalbreastgenes[,c(i)] + finalbreastcov[,11] + log2sum[,1])
		valuein9 <-lm(finalbreastL1HS[,1] ~ finalbreastgenes[,c(i)] + finalbreastcov[,12] + log2sum[,1])
		valuein10 <-lm(finalbreastL1HS[,1] ~ finalbreastgenes[,c(i)] + finalbreastcov[,13] + log2sum[,1])
		valuein11 <-lm(finalbreastL1HS[,1] ~ finalbreastgenes[,c(i)] + finalbreastcov[,14] + log2sum[,1])
		valuein12 <-lm(finalbreastL1HS[,1] ~ finalbreastgenes[,c(i)] + finalbreastcov[,15] + log2sum[,1])
		valuein13 <-lm(finalbreastL1HS[,1] ~ finalbreastgenes[,c(i)] + finalbreastcov[,16] + log2sum[,1])
		valuein14 <-lm(finalbreastL1HS[,1] ~ finalbreastgenes[,c(i)] + finalbreastcov[,17] + log2sum[,1])
		valuein15 <-lm(finalbreastL1HS[,1] ~ finalbreastgenes[,c(i)] + finalbreastcov[,18] + log2sum[,1])
		valuein16 <-lm(finalbreastL1HS[,1] ~ finalbreastgenes[,c(i)] + finalbreastcov[,19] + log2sum[,1])
		valuein17 <-lm(finalbreastL1HS[,1] ~ finalbreastgenes[,c(i)] + finalbreastcov[,20] + log2sum[,1])
		valuein18 <-lm(finalbreastL1HS[,1] ~ finalbreastgenes[,c(i)] + finalbreastcov[,21] + log2sum[,1])
		valuein19 <-lm(finalbreastL1HS[,1] ~ finalbreastgenes[,c(i)] + finalbreastcov[,22] + log2sum[,1])
		valuein20 <-lm(finalbreastL1HS[,1] ~ finalbreastgenes[,c(i)] + finalbreastcov[,23] + log2sum[,1])
		valuein21 <-lm(finalbreastL1HS[,1] ~ finalbreastgenes[,c(i)] + finalbreastcov[,24] + log2sum[,1])
		valuein22 <-lm(finalbreastL1HS[,1] ~ finalbreastgenes[,c(i)] + finalbreastcov[,25] + log2sum[,1])
		valuein23 <-lm(finalbreastL1HS[,1] ~ finalbreastgenes[,c(i)] + finalbreastcov[,26] + log2sum[,1])
		valuein24 <-lm(finalbreastL1HS[,1] ~ finalbreastgenes[,c(i)] + finalbreastcov[,27] + log2sum[,1])
		valuein25 <-lm(finalbreastL1HS[,1] ~ finalbreastgenes[,c(i)] + finalbreastcov[,28] + log2sum[,1])
		valuein26 <-lm(finalbreastL1HS[,1] ~ finalbreastgenes[,c(i)] + finalbreastcov[,29] + log2sum[,1])
		valuein27 <-lm(finalbreastL1HS[,1] ~ finalbreastgenes[,c(i)] + finalbreastcov[,30] + log2sum[,1])
		valuein28 <-lm(finalbreastL1HS[,1] ~ finalbreastgenes[,c(i)] + finalbreastcov[,31] + log2sum[,1])
		valuein29 <-lm(finalbreastL1HS[,1] ~ finalbreastgenes[,c(i)] + finalbreastcov[,32] + log2sum[,1])
		valuein30 <-lm(finalbreastL1HS[,1] ~ finalbreastgenes[,c(i)] + finalbreastcov[,33] + log2sum[,1])
		valueall <-lm(finalbreastL1HS[,1] ~ finalbreastgenes[,c(i)] + finalbreastcov[,1] + finalbreastcov[,2] + 
					finalbreastcov[,3] + finalbreastcov[,4] + finalbreastcov[,5] + finalbreastcov[,6] + finalbreastcov[,7] + 
					finalbreastcov[,8] + finalbreastcov[,9] + finalbreastcov[,10] + finalbreastcov[,11] + finalbreastcov[,12] +
					finalbreastcov[,13] + finalbreastcov[,14] + finalbreastcov[,15] + finalbreastcov[,16] + finalbreastcov[,17] + finalbreastcov[,18]+
					finalbreastcov[,19] + finalbreastcov[,20] + finalbreastcov[,21] + finalbreastcov[,22] + finalbreastcov[,23] + finalbreastcov[,24] +
					finalbreastcov[,25] + finalbreastcov[,26]+ finalbreastcov[,27] + finalbreastcov[,28] + finalbreastcov[,29] + finalbreastcov[,30] + 
					finalbreastcov[,31] + finalbreastcov[,32] + finalbreastcov[,33] + log2sum[,1])
		AICc<- model.sel(valuepc0, valuepc1, valuepc2, valuein1, valuein2, valuein3, valuein4, valuein5, valuein6, valuein7, valuein8,
			valuein9, valuein10, valuein11, valuein12, valuein13, valuein14, valuein15, valuein16, valuein17, valuein18, valuein19,
			valuein20, valuein21, valuein22, valuein23, valuein24, valuein25, valuein26, valuein27, valuein28, valuein29, valuein30,
			valueall)
		bestmodel <- eval(getCall(AICc, 1))
		finalresult<- summary(bestmodel)
		pval<- finalresult$coefficients[1,"Pr(>|t|)"]
		if (finaladrenalgenes[,c(i)] != "L1HS" && !any(finalresult$aliased) && deviance(bestmodel) > 1.0e-08) {
        effectsize <- modelEffectSizes(bestmodel)
        partialeta2 <- effectsize$Effects[1,"pEta-sqr"]
      } else {
        partialeta2 <- NA
		}
		combo<- rbind(pval, partialeta2)}
			pvalues <- cbind(pvalues, combo) 
			}

		pvalues<- data.frame(pvalues)
		tpvalues<- t(pvalues)
		skip=NULL
		
		for (i in colnames(finalbreastgenes))
		{
			skipped=NULL
			if(!all(is.na(finalbreastgenes[i]))){
			skipped<- print(paste(i))
			}
			skip<- cbind(skip, skipped)
			}

		skip<- t(skip)
		rownames(tpvalues)= skip[,1]
	
	#rvalue
	rvalues=NULL
	for (i in colnames(finalbreastgenes))
		{
		if(! all(is.na(finalbreastgenes[i]))){
		valuepc0 <-lm(finalbreastL1HS[,1] ~ finalbreastgenes[,c(i)] + finalbreastcov[,1] + log2sum[,1])
		valuepc1 <-lm(finalbreastL1HS[,1] ~ finalbreastgenes[,c(i)] + finalbreastcov[,2] + log2sum[,1])
		valuepc2 <-lm(finalbreastL1HS[,1] ~ finalbreastgenes[,c(i)] + finalbreastcov[,3] + log2sum[,1])
		valuein1 <-lm(finalbreastL1HS[,1] ~ finalbreastgenes[,c(i)] + finalbreastcov[,4] + log2sum[,1])
		valuein2 <-lm(finalbreastL1HS[,1] ~ finalbreastgenes[,c(i)] + finalbreastcov[,5] + log2sum[,1])
		valuein3 <-lm(finalbreastL1HS[,1] ~ finalbreastgenes[,c(i)] + finalbreastcov[,6] + log2sum[,1])
		valuein4 <-lm(finalbreastL1HS[,1] ~ finalbreastgenes[,c(i)] + finalbreastcov[,7] + log2sum[,1])
		valuein5 <-lm(finalbreastL1HS[,1] ~ finalbreastgenes[,c(i)] + finalbreastcov[,8] + log2sum[,1])
		valuein6 <-lm(finalbreastL1HS[,1] ~ finalbreastgenes[,c(i)] + finalbreastcov[,9] + log2sum[,1])
		valuein7 <-lm(finalbreastL1HS[,1] ~ finalbreastgenes[,c(i)] + finalbreastcov[,10] + log2sum[,1])
		valuein8 <-lm(finalbreastL1HS[,1] ~ finalbreastgenes[,c(i)] + finalbreastcov[,11] + log2sum[,1])
		valuein9 <-lm(finalbreastL1HS[,1] ~ finalbreastgenes[,c(i)] + finalbreastcov[,12] + log2sum[,1])
		valuein10 <-lm(finalbreastL1HS[,1] ~ finalbreastgenes[,c(i)] + finalbreastcov[,13] + log2sum[,1])
		valuein11 <-lm(finalbreastL1HS[,1] ~ finalbreastgenes[,c(i)] + finalbreastcov[,14] + log2sum[,1])
		valuein12 <-lm(finalbreastL1HS[,1] ~ finalbreastgenes[,c(i)] + finalbreastcov[,15] + log2sum[,1])
		valuein13 <-lm(finalbreastL1HS[,1] ~ finalbreastgenes[,c(i)] + finalbreastcov[,16] + log2sum[,1])
		valuein14 <-lm(finalbreastL1HS[,1] ~ finalbreastgenes[,c(i)] + finalbreastcov[,17] + log2sum[,1])
		valuein15 <-lm(finalbreastL1HS[,1] ~ finalbreastgenes[,c(i)] + finalbreastcov[,18] + log2sum[,1])
		valuein16 <-lm(finalbreastL1HS[,1] ~ finalbreastgenes[,c(i)] + finalbreastcov[,19] + log2sum[,1])
		valuein17 <-lm(finalbreastL1HS[,1] ~ finalbreastgenes[,c(i)] + finalbreastcov[,20] + log2sum[,1])
		valuein18 <-lm(finalbreastL1HS[,1] ~ finalbreastgenes[,c(i)] + finalbreastcov[,21] + log2sum[,1])
		valuein19 <-lm(finalbreastL1HS[,1] ~ finalbreastgenes[,c(i)] + finalbreastcov[,22] + log2sum[,1])
		valuein20 <-lm(finalbreastL1HS[,1] ~ finalbreastgenes[,c(i)] + finalbreastcov[,23] + log2sum[,1])
		valuein21 <-lm(finalbreastL1HS[,1] ~ finalbreastgenes[,c(i)] + finalbreastcov[,24] + log2sum[,1])
		valuein22 <-lm(finalbreastL1HS[,1] ~ finalbreastgenes[,c(i)] + finalbreastcov[,25] + log2sum[,1])
		valuein23 <-lm(finalbreastL1HS[,1] ~ finalbreastgenes[,c(i)] + finalbreastcov[,26] + log2sum[,1])
		valuein24 <-lm(finalbreastL1HS[,1] ~ finalbreastgenes[,c(i)] + finalbreastcov[,27] + log2sum[,1])
		valuein25 <-lm(finalbreastL1HS[,1] ~ finalbreastgenes[,c(i)] + finalbreastcov[,28] + log2sum[,1])
		valuein26 <-lm(finalbreastL1HS[,1] ~ finalbreastgenes[,c(i)] + finalbreastcov[,29] + log2sum[,1])
		valuein27 <-lm(finalbreastL1HS[,1] ~ finalbreastgenes[,c(i)] + finalbreastcov[,30] + log2sum[,1])
		valuein28 <-lm(finalbreastL1HS[,1] ~ finalbreastgenes[,c(i)] + finalbreastcov[,31] + log2sum[,1])
		valuein29 <-lm(finalbreastL1HS[,1] ~ finalbreastgenes[,c(i)] + finalbreastcov[,32] + log2sum[,1])
		valuein30 <-lm(finalbreastL1HS[,1] ~ finalbreastgenes[,c(i)] + finalbreastcov[,33] + log2sum[,1])
		valueall <-lm(finalbreastL1HS[,1] ~ finalbreastgenes[,c(i)] + finalbreastcov[,1] + finalbreastcov[,2] + 
					finalbreastcov[,3] + finalbreastcov[,4] + finalbreastcov[,5] + finalbreastcov[,6] + finalbreastcov[,7] + 
					finalbreastcov[,8] + finalbreastcov[,9] + finalbreastcov[,10] + finalbreastcov[,11] + finalbreastcov[,12] +
					finalbreastcov[,13] + finalbreastcov[,14] + finalbreastcov[,15] + finalbreastcov[,16] + finalbreastcov[,17] + finalbreastcov[,18]+
					finalbreastcov[,19] + finalbreastcov[,20] + finalbreastcov[,21] + finalbreastcov[,22] + finalbreastcov[,23] + finalbreastcov[,24] +
					finalbreastcov[,25] + finalbreastcov[,26]+ finalbreastcov[,27] + finalbreastcov[,28] + finalbreastcov[,29] + finalbreastcov[,30] + 
					finalbreastcov[,31] + finalbreastcov[,32] + finalbreastcov[,33] + log2sum[,1])
		AICc<- model.sel(valuepc0, valuepc1, valuepc2, valuein1, valuein2, valuein3, valuein4, valuein5, valuein6, valuein7, valuein8,
			valuein9, valuein10, valuein11, valuein12, valuein13, valuein14, valuein15, valuein16, valuein17, valuein18, valuein19,
			valuein20, valuein21, valuein22, valuein23, valuein24, valuein25, valuein26, valuein27, valuein28, valuein29, valuein30,
			valueall)
		bestmodel <- eval(getCall(AICc, 1))
		finalresult<- summary(bestmodel)
		rval<- finalresult$adj.r.squared 
		
		}
			rvalues <- cbind(rvalues, rval) 
			}
		rvalues<- data.frame(rvalues)
		trvalues<- t(rvalues)
		
		skip=NULL
		
		
		for (i in colnames(finalbreastgenes))
		{
			skipped=NULL
			if(!all(is.na(finalbreastgenes[i]))){
			skipped<- print(paste(i))
			}
			skip<- cbind(skip, skipped)
			}	
			skip<- t(skip)
			rownames(trvalues)= skip[,1]
			
		#qvalues
		testforq<- tpvalues[,1]
		testforq[!is.finite(testforq)] <- NA
		testforq<- data.frame(na.omit(testforq))
		qvalues<- data.frame(qvalue(testforq[,1], pi0=1)$qvalues)
		rownames(qvalues)=rownames(tpvalues)
		results<- merge(tpvalues, trvalues, by.x="row.names", by.y="row.names")	
		rownames(results)=results[,1]
		results<- results[,2:4]
		colnames(results)=c("pvalues", "Partialeta", "rsqvalues")
		finalresults<- merge(results, qvalues, by.x="row.names", by.y="row.names")
		colnames(finalresults)=c("rownames", "pvalues", "Partialeta", "rsqvalues", "qvalues")
		setwd("/home/quints1/Correlations/breast")
		write.table(finalresults, "P and Q Breast redo.txt", sep="\t", quote=FALSE)
	
		rvalueresults=NULL	
		for (i in colnames(finalbreastgenes))
		{
		if(! all(is.na(finalbreastgenes[i]))){
			pearson<- cor.test(finalbreastgenes[,c(i)], finalbreastL1HS[,1], method="pearson")$estimate
			spearman<- cor.test(finalbreastgenes[,c(i)], finalbreastL1HS[,1], method="spearman")$estimate
			combinedr<- rbind(pearson, spearman)
			}
			else{
				print(paste("Skipping", i))
			}
			rvalueresults<- cbind(rvalueresults, combinedr)
			}
		rvalueresults<- data.frame(rvalueresults)
		trvalueresults<- t(rvalueresults)
		
		skip=NULL	
		for (i in colnames(finalbreastgenes))
		{
			skipped=NULL
			if(!all(is.na(finalbreastgenes[i]))){
			skipped<- print(paste(i))
			}
			skip<- cbind(skip, skipped)
			}	
			skip<- t(skip)
			rownames(trvalueresults)= skip[,1]
				
		results<- merge(tpvalues, trvalueresults, by.x="row.names", by.y="row.names")	
		rownames(results)=results[,1]
		results<- results[,2:5]
		colnames(results)=c("pvalues", "partialeta2", "pearsonr", "spearmanr")
		finalresults<- merge(results, qvalues, by.x="row.names", by.y="row.names")
		colnames(finalresults)=c("rownames", "pvalues", "partialeta2", "pearsonr", "spearmanr", "qvalues")
		setwd("/home/quints1/Correlations/breast")
		write.table(finalresults, "PRandQBreast redo.txt", sep="\t", quote=FALSE, row.names=FALSE,col.names=TRUE)

#breast
library("stringi")
setwd("/home/quints1/Correlations/breast")
breast<- read.table("PRandQBreast redo.txt", header=TRUE, row.names=1)
breastAIC<- read.table("All AICc results breast redo.txt", header=TRUE)
breastslope<- cbind(rownames(breastAIC), breastAIC[,3])
genes<-stri_sub(breastslope[,1], 0, -9)
rownames(breastslope)=genes
breast2<- merge(breast, breastslope, by.x="row.names", by.y="row.names", drop=FALSE)
rownames(breast2)=breast2[,1]
breast2<- breast2[,c(2,3,4,5,6,8)]
breastRsq<- read.table("P and Q Breast redo.txt", header=TRUE)
sqvals<- breastRsq[,c(1,4)]
rownames(sqvals)=sqvals[,1]
breast3<- merge(breast2, sqvals, by.x="row.names", by.y="row.names", drop=FALSE)
finalbreast<- breast3[,c(1,2,3,4,5,6,7,9)]
colnames(finalbreast)=c("rownames", "pvalues", "partialeta2", "pearsonr", "spearmanr", "qvalues", "AICcoeff", "rsqvalues")
test<- ifelse(as.numeric(as.character(finalbreast[,7])) < 0, TRUE, FALSE)
result<- cbind(finalbreast, test)
cutbreastneg<- subset(result, result[,9] == "TRUE", select=c(1,2,3,4,5,6,7,8))
negbreast<- cutbreastneg[order(-as.numeric(as.character(cutbreastneg[,7]))),]
breastcutoff<- subset(negbreast, negbreast[,8] > 0.5, select=c(1,2,3,4,5,6,7,8))
breastgenes<- data.frame(breastcutoff[c(1)])	
breastpvalcut<- subset(breastcutoff, breastcutoff[,2] < 0.005, select=c(1,2,3,4,5,6,7,8))
testbreast<- breastcutoff[order(-as.numeric(as.character(breastcutoff[,3]))),]
testbreast2<- subset(testbreast, testbreast[,3] > 0.3, select=c(1,2,3,4,5,6,7,8))
genelist<- testbreast2[,1]
write.table(genelist, "etatop22.txt", sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)

library("ggplot2")
	#ACTG1
	mydata<- data.frame(cbind(finalbreastgenes[,c("ACTG1")],finalbreastL1HS[,1]))
	rownames(mydata)=rownames(finalbreastgenes)
	myplotbottom<- ggplot(mydata, aes(mydata[,1], mydata[,2])) + geom_point(color="black") + geom_smooth(method="lm", se=FALSE, color="dark red") + theme_classic() +
					labs(x="ACTG1 Gene Expression", y="L1HS Expression", title="Breast Post Test")
	#PYGO2
	mydata<- data.frame(cbind(finalbreastgenes[,c("PYGO2")],finalbreastL1HS[,1]))
	rownames(mydata)=rownames(finalbreastgenes)
	myplotbottom<- ggplot(mydata, aes(mydata[,1], mydata[,2])) + geom_point(color="black") + geom_smooth(method="lm", se=FALSE, color="dark red") + theme_classic() +
					labs(x="PYGO2 Gene Expression", y="L1HS Expression", title="Breast Post Test PYGO2")
	#TUBB8P2 
	mydata<- data.frame(cbind(finalbreastgenes[,c("TUBB8P2")],finalbreastL1HS[,1]))
	rownames(mydata)=rownames(finalbreastgenes)
	myplotbottom<- ggplot(mydata, aes(mydata[,1], mydata[,2])) + geom_point(color="black") + geom_smooth(method="lm", se=FALSE, color="dark red") + theme_classic() +
					labs(x="TUBB8P2 Gene Expression", y="L1HS Expression", title="Breast Post Test TUBB8P2")
	#MT.ND5
	mydata<- data.frame(cbind(finalbreastgenes[,c("MT.ND5")],finalbreastL1HS[,1]))
	rownames(mydata)=rownames(finalbreastgenes)
	myplotbottom<- ggplot(mydata, aes(mydata[,1], mydata[,2])) + geom_point(color="black") + geom_smooth(method="lm", se=FALSE, color="dark red") + theme_classic() +
					labs(x="MT.ND5 Gene Expression", y="L1HS Expression", title="Breast Post Test MT.ND5")
test <- subset(breastcutoff, breastcutoff[,1] == "MT.ND5", select=c(1,2,3,4,5,6,7,8))

library("stringi")
setwd("/home/quints1/Correlations/breast")
breast<- read.table("PRandQBreast redo.txt", header=TRUE, row.names=1)
breastAIC<- read.table("All AICc results breast redo.txt", header=TRUE)
breastslope<- cbind(rownames(breastAIC), breastAIC[,3])
genes<-stri_sub(breastslope[,1], 0, -9)
rownames(breastslope)=genes
breast2<- merge(breast, breastslope, by.x="row.names", by.y="row.names", drop=FALSE)
rownames(breast2)=breast2[,1]
breast2<- breast2[,c(2,3,4,5,6,8)]
breastRsq<- read.table("P and Q Breast redo.txt", header=TRUE)
sqvals<- breastRsq[,c(1,4)]
rownames(sqvals)=sqvals[,1]
breast3<- merge(breast2, sqvals, by.x="row.names", by.y="row.names", drop=FALSE)
finalbreast<- breast3[,c(1,2,3,4,5,6,7,9)]
colnames(finalbreast)=c("rownames", "pvalues", "partialeta2", "pearsonr", "spearmanr", "qvalues", "AICcoeff", "rsqvalues")
testbreast2<- subset(finalbreast, finalbreast[,3] > 0.3, select=c(1,2,3,4,5,6,7,8))
genelist<- testbreast2[,1]
write.table(genelist, "etatop22.txt", sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
	#MRPS18B
	mydata<- data.frame(cbind(finalbreastgenes[,c("MRPS18B")],finalbreastL1HS[,1]))
	rownames(mydata)=rownames(finalbreastgenes)
	myplotbottom<- ggplot(mydata, aes(mydata[,1], mydata[,2])) + geom_point(color="black") + geom_smooth(method="lm", se=FALSE, color="dark red") + theme_classic() +
					labs(x="MRPS18B Gene Expression", y="L1HS Expression", title="Breast Post Test MRPS18B")
`#ZDHHC16
	mydata<- data.frame(cbind(finalbreastgenes[,c("ZDHHC16")],finalbreastL1HS[,1]))
	rownames(mydata)=rownames(finalbreastgenes)
	myplotbottom<- ggplot(mydata, aes(mydata[,1], mydata[,2])) + geom_point(color="black") + geom_smooth(method="lm", se=FALSE, color="dark red") + theme_classic() +
					labs(x="ZDHHC16 Gene Expression", y="L1HS Expression", title="Breast Post Test ZDHHC16")
