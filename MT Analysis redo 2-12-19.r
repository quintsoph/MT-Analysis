#Adrenal
library("stringi")
setwd("/home/quints1/Correlations/adrenal")
adrenal<- read.table("PRandQAdrenal redo.txt", header=TRUE, row.names=1)
adrenalAIC<- read.table("All AICc results redo.txt", header=TRUE)
adrenalslope<- cbind(rownames(adrenalAIC), adrenalAIC[,3])
genes<-stri_sub(adrenalslope[,1], 0, -9)
rownames(adrenalslope)=genes
adrenal2<- merge(adrenal, adrenalslope, by.x="row.names", by.y="row.names", drop=FALSE)
rownames(adrenal2)=adrenal2[,1]
adrenal2<- adrenal2[,c(2,3,4,5,6,8)]
adrenalRsq<- read.table("P and Q Adrenal redo.txt", header=TRUE)
sqvals<- adrenalRsq[,c(1,4)]
rownames(sqvals)=sqvals[,1]
adrenal3<- merge(adrenal2, sqvals, by.x="row.names", by.y="row.names", drop=FALSE)
finaladrenal<- adrenal3[,c(1,2,3,4,5,6,7,9)]
colnames(finaladrenal)=c("rownames", "pvalues", "partialeta2", "pearsonr", "spearmanr", "qvalues", "AICcoeff", "rsqvalues")
testadrenal2<- subset(finaladrenal, finaladrenal[,5] < 0, select=c(1,2,3,4,5,6,7,8))
testadrenal3<- subset(testadrenal2, testadrenal2[,6] < 0.1 & testadrenal2[,3] > 0.1, select=c(1,2,3,4,5,6,7,8))
adrenaleta<- testadrenal3[order(-as.numeric(as.character(testadrenal3[,3]))),]
adrenalcut<- adrenaleta
adrenalgenes<- data.frame(adrenalcut[c(1)])	
colnames(adrenalgenes)=c("adrenal")
rownames(adrenalgenes)= adrenalgenes[,1]
#0
#Breast
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
testbreast2<- subset(finalbreast, finalbreast[,5] < 0, select=c(1,2,3,4,5,6,7,8))
testbreast3<- subset(testbreast2, testbreast2[,6] < 0.1 & testbreast2[,3] > 0.1, select=c(1,2,3,4,5,6,7,8))
breasteta<- testbreast3[order(-as.numeric(as.character(testbreast3[,3]))),]
breastcut<- breasteta
breastgenes<- data.frame(breastcut[c(1)])	
colnames(breastgenes)=c("breast")
rownames(breastgenes)= breastgenes[,1]
#45
#colonsig
library("stringi")
setwd("/home/quints1/Correlations/colonsig")
colonsig<- read.table("PRandQColonSig redo.txt", header=TRUE, row.names=1)
colonsigAIC<- read.table("AllAICcresultscolonsig redo.txt", header=TRUE)
colonsigslope<- cbind(rownames(colonsigAIC), colonsigAIC[,3])
genes<-stri_sub(colonsigslope[,1], 0, -9)
rownames(colonsigslope)=genes
colonsig2<- merge(colonsig, colonsigslope, by.x="row.names", by.y="row.names", drop=FALSE)
rownames(colonsig2)=colonsig2[,1]
colonsig2<- colonsig2[,c(2,3,4,5,6,8)]
colonsigRsq<- read.table("P and Q colonsig redo.txt", header=TRUE)
sqvals<- colonsigRsq[,c(1,4)]
rownames(sqvals)=sqvals[,1]
colonsig3<- merge(colonsig2, sqvals, by.x="row.names", by.y="row.names", drop=FALSE)
finalcolonsig<- colonsig3[,c(1,2,3,4,5,6,7,9)]
colnames(finalcolonsig)=c("rownames", "pvalues", "partialeta2", "pearsonr", "spearmanr", "qvalues", "AICcoeff", "rsqvalues")
testcolonsig2<- subset(finalcolonsig, finalcolonsig[,5] < 0, select=c(1,2,3,4,5,6,7,8))
testcolonsig3<- subset(testcolonsig2, testcolonsig2[,6] < 0.1 & testcolonsig2[,3] > 0.1, select=c(1,2,3,4,5,6,7,8))
colonsigeta<- testcolonsig3[order(-as.numeric(as.character(testcolonsig3[,3]))),]
colonsigcut<- colonsigeta
colonsiggenes<- data.frame(colonsigcut[c(1)])	
colnames(colonsiggenes)=c("colonsig")
rownames(colonsiggenes)= colonsiggenes[,1]
#0

#colontrans
library("stringi")
setwd("/home/quints1/Correlations/colontrans")
colontrans<- read.table("PRandQColonTrans redo.txt", header=TRUE, row.names=1)
colontransAIC<- read.table("AllAICcresultscolontrans redo.txt", header=TRUE)
colontransslope<- cbind(rownames(colontransAIC), colontransAIC[,3])
genes<-stri_sub(colontransslope[,1], 0, -9)
rownames(colontransslope)=genes
colontrans2<- merge(colontrans, colontransslope, by.x="row.names", by.y="row.names", drop=FALSE)
rownames(colontrans2)=colontrans2[,1]
colontrans2<- colontrans2[,c(2,3,4,5,6,8)]
colontransRsq<- read.table("P and Q colontrans redo.txt", header=TRUE)
sqvals<- colontransRsq[,c(1,4)]
rownames(sqvals)=sqvals[,1]
colontrans3<- merge(colontrans2, sqvals, by.x="row.names", by.y="row.names", drop=FALSE)
finalcolontrans<- colontrans3[,c(1,2,3,4,5,6,7,9)]
colnames(finalcolontrans)=c("rownames", "pvalues", "partialeta2", "pearsonr", "spearmanr", "qvalues", "AICcoeff", "rsqvalues")
testcolontrans2<- subset(finalcolontrans, finalcolontrans[,5] < 0, select=c(1,2,3,4,5,6,7,8))
testcolontrans3<- subset(testcolontrans2, testcolontrans2[,6] < 0.1 & testcolontrans2[,3] > 0.1, select=c(1,2,3,4,5,6,7,8))
colontranseta<- testcolontrans3[order(-as.numeric(as.character(testcolontrans3[,3]))),]
colontranscut<- colontranseta
colontransgenes<- data.frame(colontranscut[c(1)])	
colnames(colontransgenes)=c("colontrans")
rownames(colontransgenes)= colontransgenes[,1]
#0

#esophgastro
library("stringi")
setwd("/home/quints1/Correlations/esophagusgastro")
esophgastro<- read.table("PRandQEsophGastro redo.txt", header=TRUE, row.names=1)
esophgastroAIC<- read.table("AllAICcresultsesophgastro redo.txt", header=TRUE)
esophgastroslope<- cbind(rownames(esophgastroAIC), esophgastroAIC[,3])
genes<-stri_sub(esophgastroslope[,1], 0, -9)
rownames(esophgastroslope)=genes
esophgastro2<- merge(esophgastro, esophgastroslope, by.x="row.names", by.y="row.names", drop=FALSE)
rownames(esophgastro2)=esophgastro2[,1]
esophgastro2<- esophgastro2[,c(2,3,4,5,6,8)]
esophgastroRsq<- read.table("P and Q esophgastro redo.txt", header=TRUE)
sqvals<- esophgastroRsq[,c(1,4)]
rownames(sqvals)=sqvals[,1]
esophgastro3<- merge(esophgastro2, sqvals, by.x="row.names", by.y="row.names", drop=FALSE)
finalesophgastro<- esophgastro3[,c(1,2,3,4,5,6,7,9)]
colnames(finalesophgastro)=c("rownames", "pvalues", "partialeta2", "pearsonr", "spearmanr", "qvalues", "AICcoeff", "rsqvalues")
testesophgastro2<- subset(finalesophgastro, finalesophgastro[,5] < 0, select=c(1,2,3,4,5,6,7,8))
testesophgastro3<- subset(testesophgastro2, testesophgastro2[,6] < 0.1 & testesophgastro2[,3] > 0.1, select=c(1,2,3,4,5,6,7,8))
esophgastroeta<- testesophgastro3[order(-as.numeric(as.character(testesophgastro3[,3]))),]
esophgastrocut<- esophgastroeta
esophgastrogenes<- data.frame(esophgastrocut[c(1)])	
colnames(esophgastrogenes)=c("esophgastro")
rownames(esophgastrogenes)= esophgastrogenes[,1]
#0

#esophmuscularis
library("stringi")
setwd("/home/quints1/Correlations/esophagusmuscularis")
esophmuscularis<- read.table("PRandQEsophMuscularis redo.txt", header=TRUE, row.names=1)
esophmuscularisAIC<- read.table("AllAICcresultsesophmuscularis redo.txt", header=TRUE)
esophmuscularisslope<- cbind(rownames(esophmuscularisAIC), esophmuscularisAIC[,3])
genes<-stri_sub(esophmuscularisslope[,1], 0, -9)
rownames(esophmuscularisslope)=genes
esophmuscularis2<- merge(esophmuscularis, esophmuscularisslope, by.x="row.names", by.y="row.names", drop=FALSE)
rownames(esophmuscularis2)=esophmuscularis2[,1]
esophmuscularis2<- esophmuscularis2[,c(2,3,4,5,6,8)]
esophmuscularisRsq<- read.table("P and Q esophmuscularis redo.txt", header=TRUE)
sqvals<- esophmuscularisRsq[,c(1,4)]
rownames(sqvals)=sqvals[,1]
esophmuscularis3<- merge(esophmuscularis2, sqvals, by.x="row.names", by.y="row.names", drop=FALSE)
finalesophmuscularis<- esophmuscularis3[,c(1,2,3,4,5,6,7,9)]
colnames(finalesophmuscularis)=c("rownames", "pvalues", "partialeta2", "pearsonr", "spearmanr", "qvalues", "AICcoeff", "rsqvalues")
testesophmuscularis2<- subset(finalesophmuscularis, finalesophmuscularis[,5] < 0, select=c(1,2,3,4,5,6,7,8))
testesophmuscularis3<- subset(testesophmuscularis2, testesophmuscularis2[,6] < 0.1 & testesophmuscularis2[,3] > 0.1, select=c(1,2,3,4,5,6,7,8))
esophmusculariseta<- testesophmuscularis3[order(-as.numeric(as.character(testesophmuscularis3[,3]))),]
esophmusculariscut<- esophmusculariseta
esophmuscularisgenes<- data.frame(esophmusculariscut[c(1)])	
colnames(esophmuscularisgenes)=c("esophmuscularis")
rownames(esophmuscularisgenes)= esophmuscularisgenes[,1]
#0

#esophmucosa
library("stringi")
setwd("/home/quints1/Correlations/esophagusmucosa")
esophmucosa<- read.table("PRandQEsophMucosa redo.txt", header=TRUE, row.names=1)
esophmucosaAIC<- read.table("AllAICcresultsesophmucosa redo.txt", header=TRUE)
esophmucosaslope<- cbind(rownames(esophmucosaAIC), esophmucosaAIC[,3])
genes<-stri_sub(esophmucosaslope[,1], 0, -9)
rownames(esophmucosaslope)=genes
esophmucosa2<- merge(esophmucosa, esophmucosaslope, by.x="row.names", by.y="row.names", drop=FALSE)
rownames(esophmucosa2)=esophmucosa2[,1]
esophmucosa2<- esophmucosa2[,c(2,3,4,5,6,8)]
esophmucosaRsq<- read.table("P and Q esophmucosa redo.txt", header=TRUE)
sqvals<- esophmucosaRsq[,c(1,4)]
rownames(sqvals)=sqvals[,1]
esophmucosa3<- merge(esophmucosa2, sqvals, by.x="row.names", by.y="row.names", drop=FALSE)
finalesophmucosa<- esophmucosa3[,c(1,2,3,4,5,6,7,9)]
colnames(finalesophmucosa)=c("rownames", "pvalues", "partialeta2", "pearsonr", "spearmanr", "qvalues", "AICcoeff", "rsqvalues")
testesophmucosa2<- subset(finalesophmucosa, finalesophmucosa[,5] < 0, select=c(1,2,3,4,5,6,7,8))
testesophmucosa3<- subset(testesophmucosa2, testesophmucosa2[,6] < 0.1 & testesophmucosa2[,3] > 0.1, select=c(1,2,3,4,5,6,7,8))
esophmucosaeta<- testesophmucosa3[order(-as.numeric(as.character(testesophmucosa3[,3]))),]
esophmucosacut<- esophmucosaeta
esophmucosagenes<- data.frame(esophmucosacut[c(1)])	
colnames(esophmucosagenes)=c("esophmucosa")
rownames(esophmucosagenes)= esophmucosagenes[,1]
#147

#kidney
library("stringi")
setwd("/home/quints1/Correlations/kidney")
kidney<- read.table("PRandQKidney redo.txt", header=TRUE, row.names=1)
kidneyAIC<- read.table("AllAICcresultskidney redo.txt", header=TRUE)
kidneyslope<- cbind(rownames(kidneyAIC), kidneyAIC[,3])
genes<-stri_sub(kidneyslope[,1], 0, -9)
rownames(kidneyslope)=genes
kidney2<- merge(kidney, kidneyslope, by.x="row.names", by.y="row.names", drop=FALSE)
rownames(kidney2)=kidney2[,1]
kidney2<- kidney2[,c(2,3,4,5,6,8)]
kidneyRsq<- read.table("P and Q kidney redo.txt", header=TRUE)
sqvals<- kidneyRsq[,c(1,4)]
rownames(sqvals)=sqvals[,1]
kidney3<- merge(kidney2, sqvals, by.x="row.names", by.y="row.names", drop=FALSE)
finalkidney<- kidney3[,c(1,2,3,4,5,6,7,9)]
colnames(finalkidney)=c("rownames", "pvalues", "partialeta2", "pearsonr", "spearmanr", "qvalues", "AICcoeff", "rsqvalues")
testkidney2<- subset(finalkidney, finalkidney[,5] < 0, select=c(1,2,3,4,5,6,7,8))
testkidney3<- subset(testkidney2, testkidney2[,6] < 0.1 & testkidney2[,3] > 0.1, select=c(1,2,3,4,5,6,7,8))
kidneyeta<- testkidney3[order(-as.numeric(as.character(testkidney3[,3]))),]
kidneycut<- kidneyeta
kidneygenes<- data.frame(kidneycut[c(1)])	
colnames(kidneygenes)=c("kidney")
rownames(kidneygenes)= kidneygenes[,1]
#142

#liver
library("stringi")
setwd("/home/quints1/Correlations/liver")
liver<- read.table("PRandQLiver redo.txt", header=TRUE, row.names=1)
liverAIC<- read.table("AllAICcresultsliver redo.txt", header=TRUE)
liverslope<- cbind(rownames(liverAIC), liverAIC[,3])
genes<-stri_sub(liverslope[,1], 0, -9)
rownames(liverslope)=genes
liver2<- merge(liver, liverslope, by.x="row.names", by.y="row.names", drop=FALSE)
rownames(liver2)=liver2[,1]
liver2<- liver2[,c(2,3,4,5,6,8)]
liverRsq<- read.table("P and Q liver redo.txt", header=TRUE)
sqvals<- liverRsq[,c(1,4)]
rownames(sqvals)=sqvals[,1]
liver3<- merge(liver2, sqvals, by.x="row.names", by.y="row.names", drop=FALSE)
finalliver<- liver3[,c(1,2,3,4,5,6,7,9)]
colnames(finalliver)=c("rownames", "pvalues", "partialeta2", "pearsonr", "spearmanr", "qvalues", "AICcoeff", "rsqvalues")
testliver2<- subset(finalliver, finalliver[,5] < 0, select=c(1,2,3,4,5,6,7,8))
testliver3<- subset(testliver2, testliver2[,6] < 0.1 & testliver2[,3] > 0.1, select=c(1,2,3,4,5,6,7,8))
livereta<- testliver3[order(-as.numeric(as.character(testliver3[,3]))),]
livercut<- livereta
livergenes<- data.frame(livercut[c(1)])	
colnames(livergenes)=c("liver")
rownames(livergenes)= livergenes[,1]
#0

#pancreas
library("stringi")
setwd("/home/quints1/Correlations/pancreas")
pancreas<- read.table("PRandQPancreas redo.txt", header=TRUE, row.names=1)
pancreasAIC<- read.table("AllAICcresultspancreas redo.txt", header=TRUE)
pancreasslope<- cbind(rownames(pancreasAIC), pancreasAIC[,3])
genes<-stri_sub(pancreasslope[,1], 0, -9)
rownames(pancreasslope)=genes
pancreas2<- merge(pancreas, pancreasslope, by.x="row.names", by.y="row.names", drop=FALSE)
rownames(pancreas2)=pancreas2[,1]
pancreas2<- pancreas2[,c(2,3,4,5,6,8)]
pancreasRsq<- read.table("P and Q pancreas redo.txt", header=TRUE)
sqvals<- pancreasRsq[,c(1,4)]
rownames(sqvals)=sqvals[,1]
pancreas3<- merge(pancreas2, sqvals, by.x="row.names", by.y="row.names", drop=FALSE)
finalpancreas<- pancreas3[,c(1,2,3,4,5,6,7,9)]
colnames(finalpancreas)=c("rownames", "pvalues", "partialeta2", "pearsonr", "spearmanr", "qvalues", "AICcoeff", "rsqvalues")
testpancreas2<- subset(finalpancreas, finalpancreas[,5] < 0, select=c(1,2,3,4,5,6,7,8))
testpancreas3<- subset(testpancreas2, testpancreas2[,6] < 0.1 & testpancreas2[,3] > 0.1, select=c(1,2,3,4,5,6,7,8))
pancreaseta<- testpancreas3[order(-as.numeric(as.character(testpancreas3[,3]))),]
pancreascut<- pancreaseta
pancreasgenes<- data.frame(pancreascut[c(1)])	
colnames(pancreasgenes)=c("pancreas")
rownames(pancreasgenes)= pancreasgenes[,1]
#0

#prostate
library("stringi")
setwd("/home/quints1/Correlations/prostate")
prostate<- read.table("PRandQProstate redo.txt", header=TRUE, row.names=1)
prostateAIC<- read.table("AllAICcresultsprostate redo.txt", header=TRUE)
prostateslope<- cbind(rownames(prostateAIC), prostateAIC[,3])
genes<-stri_sub(prostateslope[,1], 0, -9)
rownames(prostateslope)=genes
prostate2<- merge(prostate, prostateslope, by.x="row.names", by.y="row.names", drop=FALSE)
rownames(prostate2)=prostate2[,1]
prostate2<- prostate2[,c(2,3,4,5,6,8)]
prostateRsq<- read.table("P and Q prostate redo.txt", header=TRUE)
sqvals<- prostateRsq[,c(1,4)]
rownames(sqvals)=sqvals[,1]
prostate3<- merge(prostate2, sqvals, by.x="row.names", by.y="row.names", drop=FALSE)
finalprostate<- prostate3[,c(1,2,3,4,5,6,7,9)]
colnames(finalprostate)=c("rownames", "pvalues", "partialeta2", "pearsonr", "spearmanr", "qvalues", "AICcoeff", "rsqvalues")
testprostate2<- subset(finalprostate, finalprostate[,5] < 0, select=c(1,2,3,4,5,6,7,8))
testprostate3<- subset(testprostate2, testprostate2[,6] < 0.1 & testprostate2[,3] > 0.1, select=c(1,2,3,4,5,6,7,8))
prostateeta<- testprostate2[order(-as.numeric(as.character(testprostate3[,3]))),]
prostatecut<- prostateeta
prostategenes<- data.frame(prostatecut[c(1)])	
colnames(prostategenes)=c("prostate")
rownames(prostategenes)= prostategenes[,1]
#88

#stomach
library("stringi")
setwd("/home/quints1/Correlations/stomach")
stomach<- read.table("PRandQStomach redo.txt", header=TRUE, row.names=1)
stomachAIC<- read.table("AllAICcresultsstomach redo.txt", header=TRUE)
stomachslope<- cbind(rownames(stomachAIC), stomachAIC[,3])
genes<-stri_sub(stomachslope[,1], 0, -9)
rownames(stomachslope)=genes
stomach2<- merge(stomach, stomachslope, by.x="row.names", by.y="row.names", drop=FALSE)
rownames(stomach2)=stomach2[,1]
stomach2<- stomach2[,c(2,3,4,5,6,8)]
stomachRsq<- read.table("P and Q stomach redo.txt", header=TRUE)
sqvals<- stomachRsq[,c(1,4)]
rownames(sqvals)=sqvals[,1]
stomach3<- merge(stomach2, sqvals, by.x="row.names", by.y="row.names", drop=FALSE)
finalstomach<- stomach3[,c(1,2,3,4,5,6,7,9)]
colnames(finalstomach)=c("rownames", "pvalues", "partialeta2", "pearsonr", "spearmanr", "qvalues", "AICcoeff", "rsqvalues")
teststomach2<- subset(finalstomach, finalstomach[,5] < 0, select=c(1,2,3,4,5,6,7,8))
teststomach3<- subset(teststomach2, teststomach2[,6] < 0.1 & teststomach2[,3] > 0.1, select=c(1,2,3,4,5,6,7,8))
stomacheta<- teststomach3[order(-as.numeric(as.character(teststomach3[,3]))),]
stomachcut<- stomacheta
stomachgenes<- data.frame(stomachcut[c(1)])	
colnames(stomachgenes)=c("stomach")
rownames(stomachgenes)= stomachgenes[,1]
#0

#thyroid
library("stringi")
setwd("/home/quints1/Correlations/thyroid")
thyroid<- read.table("PRandQThyroid redo.txt", header=TRUE, row.names=1)
thyroidAIC<- read.table("AllAICcresultsthyroid redo.txt", header=TRUE)
thyroidslope<- cbind(rownames(thyroidAIC), thyroidAIC[,3])
genes<-stri_sub(thyroidslope[,1], 0, -9)
rownames(thyroidslope)=genes
thyroid2<- merge(thyroid, thyroidslope, by.x="row.names", by.y="row.names", drop=FALSE)
rownames(thyroid2)=thyroid2[,1]
thyroid2<- thyroid2[,c(2,3,4,5,6,8)]
thyroidRsq<- read.table("P and Q thyroid redo.txt", header=TRUE)
sqvals<- thyroidRsq[,c(1,4)]
rownames(sqvals)=sqvals[,1]
thyroid3<- merge(thyroid2, sqvals, by.x="row.names", by.y="row.names", drop=FALSE)
finalthyroid<- thyroid3[,c(1,2,3,4,5,6,7,9)]
colnames(finalthyroid)=c("rownames", "pvalues", "partialeta2", "pearsonr", "spearmanr", "qvalues", "AICcoeff", "rsqvalues")
testthyroid2<- subset(finalthyroid, finalthyroid[,5] < 0, select=c(1,2,3,4,5,6,7,8))
testthyroid3<- subset(testthyroid2, testthyroid2[,6] < 0.1 & testthyroid2[,3] > 0.1, select=c(1,2,3,4,5,6,7,8))
thyroideta<- testthyroid2[order(-as.numeric(as.character(testthyroid3[,3]))),]
thyroidcut<- thyroideta
thyroidgenes<- data.frame(thyroidcut[c(1)])	
colnames(thyroidgenes)=c("thyroid")
rownames(thyroidgenes)= thyroidgenes[,1]
#0

#merge all the tissues
merge1<- merge(breastgenes, adrenalgenes, by="row.names", all=TRUE)
rownames(merge1)=merge1[,1]
finalmerge1<- merge1[,c(2,3)]
merge2<- merge(colontransgenes, colonsiggenes, by="row.names", all=TRUE)
rownames(merge2)=merge2[,1]
finalmerge2<- merge2[,c(2,3)]
merge1_2<- merge(finalmerge2, finalmerge1, by="row.names", all=TRUE)
rownames(merge1_2)=merge1_2[,1]
finalmerge1_2<- merge1_2[,c(2,3,4,5)]
merge3<- merge(esophgastrogenes, esophmucosagenes, by="row.names", all=TRUE)
rownames(merge3)=merge3[,1]
finalmerge3<- merge3[,c(2,3)]
merge4<- merge(esophmuscularisgenes, kidneygenes, by="row.names", all=TRUE)
rownames(merge4)=merge4[,1]
finalmerge4<- merge4[,c(2,3)]
merge3_4<- merge(finalmerge3, finalmerge4, by="row.names", all=TRUE)
rownames(merge3_4)=merge3_4[,1]
finalmerge3_4<- merge3_4[,c(2,3,4,5)]
merge5<- merge(livergenes, pancreasgenes, by="row.names", all=TRUE)
rownames(merge5)=merge5[,1]
finalmerge5<- merge5[,c(2,3)]
merge6<- merge(prostategenes, stomachgenes, by="row.names", all=TRUE)
rownames(merge6)=merge6[,1]
finalmerge6<- merge6[,c(2,3)]
merge5_6<- merge(finalmerge5, finalmerge6, by="row.names", all=TRUE)
rownames(merge5_6)=merge5_6[,1]
finalmerge5_6<- merge5_6[,c(2,3,4,5)]
merge1234<- merge(finalmerge1_2, finalmerge3_4, by="row.names", all=TRUE)
rownames(merge1234)=merge1234[,1]
finalmerge1234<- merge1234[,c(2,3,4,5,6,7,8,9)]
finalmerge123456<- merge(finalmerge1234, finalmerge5_6, by="row.names", all=TRUE)
rownames(finalmerge123456)=finalmerge123456[,1]
result<- finalmerge123456[,c(2,3,4,5,6,7,8,9,10,11,12,13)]
finalallmerge<- merge(result, thyroidgenes, by="row.names", all=TRUE)
rownames(finalallmerge)=finalallmerge[,1]
final<- finalallmerge[,c(2,3,4,5,6,7,8,9,10,11,12,13,14)]

###Tally up
mylist=NULL
for (i in rownames(final))
	{
	test<- ifelse(final[c(i),] == "NA", 0, 1)
	summation<-sum(test, na.rm=TRUE)
	result<- c(i, summation)
	mylist<- rbind(mylist, result)
	}
rownames(mylist)=mylist[,1]
mylist2<- data.frame(as.numeric(as.character(mylist[,2])))
rownames(mylist2)=rownames(mylist)

#find the gene in the most tissues
test<- ifelse(mylist2[,1] > 1, TRUE, FALSE)
combined<- cbind(mylist2, test)
subbed<- subset(combined, combined[,2] == "TRUE", select=c(1,2))
saver<- rownames(subbed)
setwd("/home/quints1/Correlations/")
write.table(saver, "genesin2tissues redo2.txt", quote=FALSE, col.names=FALSE, row.names=FALSE)

#get from which tissues they appear
#tissueoutput=NULL
#for (i in rownames(subbed)){
#	y<- final[c(i),]
#	tissueoutput<- rbind(tissueoutput, y)
#	}
#write.table(tissueoutput, "genesin2finaloutput.txt", quote=FALSE, row.names=TRUE, col.names=TRUE, sep="\t")	
	
#generate graphs
	#upload all tissue files
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
	#adrenal
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
	rownames(adrenalall)=adrenalall[,1]
	adrenalall<- adrenalall[,-c(1,2)]
	adrenalallcombo<-  merge(adrenalall, tlog2sub, by.x="row.names", by.y="row.names")
	rownames(adrenalallcombo)=adrenalallcombo[,1]
	finaladrenalcov<- data.frame(adrenalallcombo[,c(56205:56222)])	
	rownames(finaladrenalcov)=adrenalallcombo[,1]
	finaladrenalL1HS<- data.frame(adrenalallcombo[,c(56204)])
	rownames(finaladrenalL1HS)= adrenalallcombo[,1]
	colnames(finaladrenalL1HS)= c("L1HS")
	finaladrenalgenes<- data.frame(adrenalallcombo[,c(2:56203)])
	log2sum<- data.frame(adrenalallcombo[,56223])
	colnames(log2sum)=c("log2sum")
	rownames(log2sum)=rownames(adrenalallcombo)
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
	rownames(breastall)=breastall[,1]
	breastall<- breastall[,-c(1,2)]
	breastallcombo<-  merge(breastall, tlog2sub, by.x="row.names", by.y="row.names")
	rownames(breastallcombo)=breastallcombo[,1]
	finalbreastcov<- data.frame(breastallcombo[,c(56205:56237)])	
	rownames(finalbreastcov)=breastallcombo[,1]
	finalbreastL1HS<- data.frame(breastallcombo[,c(56204)])
	rownames(finalbreastL1HS)= breastallcombo[,1]
	colnames(finalbreastL1HS)= c("L1HS")
	finalbreastgenes<- data.frame(breastallcombo[,c(2:56203)])
	log2sum<- data.frame(breastallcombo[,56238])
	colnames(log2sum)=c("log2sum")
	rownames(log2sum)=rownames(breastallcombo)
	#colontrans
	setwd("/home/quints1/Correlations/L1HSpatientedits")
	colonpats<- read.table("Colonedits.txt", header=TRUE, sep="\t")
	colontranspats<- subset(colonpats, colonpats[,3] == "Colon - Transverse", select=c(1,2,3))
	rownames(colontranspats)=str_replace_all(colontranspats[,1], "-", ".")
	colontransgenes<- merge(colontranspats, tgenesraw, by.x="row.names", by.y="row.names", drop=FALSE)
	colontransgenesfinal<- colontransgenes[,c(4:56207)]
	rownames(colontransgenesfinal)=colontransgenes[,1]
		#read in the covariate
		setwd("/home/quints1/Correlations/colontrans")
		colontranscov<- read.table("colontranseditcov.txt", header=TRUE, row.names=1, sep="\t")
	colontransall<- merge(colontransgenesfinal, colontranscov, by.x="row.names", by.y="row.names")
	rownames(colontransall)=colontransall[,1]
	colontransall<- colontransall[,c(3:56236)]
	colontransallcombo<-  merge(colontransall, tlog2sub, by.x="row.names", by.y="row.names")
	rownames(colontransallcombo)=colontransallcombo[,1]
	finalcolontranscov<- data.frame(colontransallcombo[,c(56205:56235)])	
	rownames(finalcolontranscov)=colontransallcombo[,1]
	finalcolontransL1HS<- data.frame(colontransallcombo[,c(56204)])
	rownames(finalcolontransL1HS)= colontransallcombo[,1]
	colnames(finalcolontransL1HS)= c("L1HS")
	finalcolontransgenes<- data.frame(colontransallcombo[,c(2:56203)])
	log2sum<- data.frame(colontransallcombo[,56236])
	colnames(log2sum)=c("log2sum")
	rownames(log2sum)=rownames(colontransallcombo)
	#colonsig
	setwd("/home/quints1/Correlations/L1HSpatientedits")
	colonpats<- read.table("Colonedits.txt", header=TRUE, sep="\t")
	colonsigpats<- subset(colonpats, colonpats[,3] == "Colon - Sigmoid", select=c(1,2,3))
	rownames(colonsigpats)=str_replace_all(colonsigpats[,1], "-", ".")
	colonsiggenes<- merge(colonsigpats, tgenesraw, by.x="row.names", by.y="row.names", drop=FALSE)
	colonsiggenesfinal<- colonsiggenes[,c(4:56207)]
	rownames(colonsiggenesfinal)=colonsiggenes[,1]
		#read in the covariate
		setwd("/home/quints1/Correlations/colonsig")
		colonsigcov<- read.table("colonsigeditcov.txt", header=TRUE, row.names=1, sep="\t")
	colonsigall<- merge(colonsiggenesfinal, colonsigcov, by.x="row.names", by.y="row.names")
	rownames(colonsigall)=colonsigall[,1]
	colonsigall<- colonsigall[,-c(1,2)]
	colonsigallcombo<-  merge(colonsigall, tlog2sub, by.x="row.names", by.y="row.names")
	rownames(colonsigallcombo)=colonsigallcombo[,1]
	finalcolonsigcov<- data.frame(colonsigallcombo[,c(56205:56237)])	
	rownames(finalcolonsigcov)=colonsigallcombo[,1]
	finalcolonsigL1HS<- data.frame(colonsigallcombo[,c(56204)])
	rownames(finalcolonsigL1HS)= colonsigallcombo[,1]
	colnames(finalcolonsigL1HS)= c("L1HS")
	finalcolonsiggenes<- data.frame(colonsigallcombo[,c(2:56203)])
	log2sum<- data.frame(colonsigallcombo[,56238])
	colnames(log2sum)=c("log2sum")
	rownames(log2sum)=rownames(colonsigallcombo)
	#esophagusmucosa
	setwd("/home/quints1/Correlations/L1HSpatientedits")
	esophpats<- read.table("Esophagusedits.txt", header=TRUE, sep="\t")
	esophpats<- subset(esophpats, esophpats[,3] == "Esophagus - Mucosa", select=c(1,2,3))
	rownames(esophpats)=str_replace_all(esophpats[,1], "-", ".")
	esophgenes<- merge(esophpats, tgenesraw, by.x="row.names", by.y="row.names", drop=FALSE)
	esophgenesfinal<- esophgenes[,c(4:56207)]
	rownames(esophgenesfinal)=esophgenes[,1]
		#read in the covariate
		setwd("/home/quints1/Correlations/esophagusmucosa")
		esophcov<- read.table("esophagusmucosaeditcov.txt", header=TRUE, row.names=1, sep="\t")
	esophall<- merge(esophgenesfinal, esophcov, by.x="row.names", by.y="row.names")
	rownames(esophall)=esophall[,1]
	esophall<- esophall[,-c(1,2)]
	esophallcombo<-  merge(esophall, tlog2sub, by.x="row.names", by.y="row.names")
	rownames(esophallcombo)=esophallcombo[,1]
	finalesophmucosacov<- data.frame(esophallcombo[,c(56205:56252)])	
	rownames(finalesophmucosacov)=esophallcombo[,1]
	finalesophmucosaL1HS<- data.frame(esophallcombo[,c(56204)])
	rownames(finalesophmucosaL1HS)= esophallcombo[,1]
	colnames(finalesophmucosaL1HS)= c("L1HS")
	finalesophmucosagenes<- data.frame(esophallcombo[,c(2:56203)])
	log2sum<- data.frame(esophallcombo[,56253])
	colnames(log2sum)=c("log2sum")
	rownames(log2sum)=rownames(esophallcombo)
	#esophagus - muscularis
		#Use the old L1HS patients to get correct patient output
	setwd("/home/quints1/Correlations/L1HSpatientedits")
	esophpats<- read.table("Esophagusedits.txt", header=TRUE, sep="\t")
	esophpats<- subset(esophpats, esophpats[,3] == "Esophagus - Muscularis", select=c(1,2,3))
	rownames(esophpats)=str_replace_all(esophpats[,1], "-", ".")
	esophgenes<- merge(esophpats, tgenesraw, by.x="row.names", by.y="row.names", drop=FALSE)
	esophgenesfinal<- esophgenes[,c(4:56207)]
	rownames(esophgenesfinal)=esophgenes[,1]
		#read in the covariate
		setwd("/home/quints1/Correlations/esophagusmuscularis")
		esophcov<- read.table("esophagusmusculariseditcov.txt", header=TRUE, row.names=1, sep="\t")
	esophall<- merge(esophgenesfinal, esophcov, by.x="row.names", by.y="row.names")
	rownames(esophall)=esophall[,1]
	esophall<- esophall[,-c(1,2)]
	esophallcombo<-  merge(esophall, tlog2sub, by.x="row.names", by.y="row.names")
	rownames(esophallcombo)=esophallcombo[,1]
	finalesophmusculariscov<- data.frame(esophallcombo[,c(56205:56252)])	
	rownames(finalesophmusculariscov)=esophallcombo[,1]
	finalesophmuscularisL1HS<- data.frame(esophallcombo[,c(56204)])
	rownames(finalesophmuscularisL1HS)= esophallcombo[,1]
	colnames(finalesophmuscularisL1HS)= c("L1HS")
	finalesophmuscularisgenes<- data.frame(esophallcombo[,c(2:56203)])
	log2sum<- data.frame(esophallcombo[,56253])
	colnames(log2sum)=c("log2sum")
	rownames(log2sum)=rownames(esophallcombo)
	#esophagus - gastroesophageal junction
		#Use the old L1HS patients to get correct patient output
	setwd("/home/quints1/Correlations/L1HSpatientedits")
	esophpats<- read.table("Esophagusedits.txt", header=TRUE, sep="\t")
	esophpats<- subset(esophpats, esophpats[,3] == "Esophagus - Gastroesophageal Junction", select=c(1,2,3))
	rownames(esophpats)=str_replace_all(esophpats[,1], "-", ".")
	esophgenes<- merge(esophpats, tgenesraw, by.x="row.names", by.y="row.names", drop=FALSE)
	esophgenesfinal<- esophgenes[,c(4:56207)]
	rownames(esophgenesfinal)=esophgenes[,1]
		#read in the covariate
		setwd("/home/quints1/Correlations/esophagusgastro")
		esophcov<- read.table("esophagusgastroeditcov.txt", header=TRUE, row.names=1, sep="\t")
	esophall<- merge(esophgenesfinal, esophcov, by.x="row.names", by.y="row.names")
	rownames(esophall)=esophall[,1]
	esophall<- esophall[,-c(1,2)]
	esophallcombo<-  merge(esophall, tlog2sub, by.x="row.names", by.y="row.names")
	rownames(esophallcombo)=esophallcombo[,1]
	finalesophgastrocov<- data.frame(esophallcombo[,c(56205:56222)])	
	rownames(finalesophgastrocov)=esophallcombo[,1]
	finalesophgastroL1HS<- data.frame(esophallcombo[,c(56204)])
	rownames(finalesophgastroL1HS)= esophallcombo[,1]
	colnames(finalesophgastroL1HS)= c("L1HS")
	finalesophgastrogenes<- data.frame(esophallcombo[,c(2:56203)])

	#kidney
		#Use the old L1HS patients to get correct patient output
	setwd("/home/quints1/Correlations/L1HSpatientedits")
	kidpats<- read.table("Kidneyedits.txt", header=TRUE, sep="\t")
	rownames(kidpats)=str_replace_all(kidpats[,1], "-", ".")
	kidgenes<- merge(kidpats, tgenesraw, by.x="row.names", by.y="row.names", drop=FALSE)
	kidgenesfinal<- kidgenes[,c(4:56207)]
	rownames(kidgenesfinal)=kidgenes[,1]
		#read in the covariate
		setwd("/home/quints1/Correlations/kidney")
		kidneycov<- read.table("kidneyeditcov.txt", header=TRUE, row.names=1, sep="\t")
	kidall<- merge(kidgenesfinal, kidneycov, by.x="row.names", by.y="row.names")
	rownames(kidall)=kidall[,1]
	kidall<- kidall[,-c(1,2)]
	kidallcombo<-  merge(kidall, tlog2sub, by.x="row.names", by.y="row.names")
	rownames(kidallcombo)=kidallcombo[,1]
	finalkidneycov<- data.frame(kidallcombo[,c(56205:56222)])	
	rownames(finalkidneycov)=kidallcombo[,1]
	finalkidneyL1HS<- data.frame(kidallcombo[,c(56204)])
	rownames(finalkidneyL1HS)= kidallcombo[,1]
	colnames(finalkidneyL1HS)= c("L1HS")
	finalkidneygenes<- data.frame(kidallcombo[,c(2:56203)])
	log2sum<- data.frame(kidallcombo[,56223])
	colnames(log2sum)=c("log2sum")
	rownames(log2sum)=rownames(kidallcombo)
	#liver
		#Use the old L1HS patients to get correct patient output
	setwd("/home/quints1/Correlations/L1HSpatientedits")
	livpats<- read.table("Liveredits.txt", header=TRUE, sep="\t")
	rownames(livpats)=str_replace_all(livpats[,1], "-", ".")
	livgenes<- merge(livpats, tgenesraw, by.x="row.names", by.y="row.names", drop=FALSE)
	livgenesfinal<- livgenes[,c(4:56207)]
	rownames(livgenesfinal)=livgenes[,1]
		#read in the covariate
		setwd("/home/quints1/Correlations/liver")
		livcov<- read.table("livereditcov.txt", header=TRUE, row.names=1, sep="\t")
	livall<- merge(livgenesfinal, livcov, by.x="row.names", by.y="row.names")
	rownames(livall)=livall[,1]
	livall<- livall[,-c(1,2)]
	livallcombo<-  merge(livall, tlog2sub, by.x="row.names", by.y="row.names")
	rownames(livallcombo)=livallcombo[,1]
	finallivercov<- data.frame(livallcombo[,c(56205:56222)])	
	rownames(finallivercov)=livallcombo[,1]
	finalliverL1HS<- data.frame(livallcombo[,c(56204)])
	rownames(finalliverL1HS)= livallcombo[,1]
	colnames(finalliverL1HS)= c("L1HS")
	finallivergenes<- data.frame(livallcombo[,c(2:56203)])
	log2sum<- data.frame(livallcombo[,56223])
	colnames(log2sum)=c("log2sum")
	rownames(log2sum)=rownames(livallcombo)
	#pancreas
		#Use the old L1HS patients to get correct patient output
	setwd("/home/quints1/Correlations/L1HSpatientedits")
	panpats<- read.table("Pancreasedits.txt", header=TRUE, sep="\t")
	rownames(panpats)=str_replace_all(panpats[,1], "-", ".")
	pangenes<- merge(panpats, tgenesraw, by.x="row.names", by.y="row.names", drop=FALSE)
	pangenesfinal<- pangenes[,c(4:56207)]
	rownames(pangenesfinal)=pangenes[,1]
		#read in the covariate
		setwd("/home/quints1/Correlations/pancreas")
		pancov<- read.table("pancreaseditcov.txt", header=TRUE, row.names=1, sep="\t")
	panall<- merge(pangenesfinal, pancov, by.x="row.names", by.y="row.names")
	rownames(panall)=panall[,1]
	panall<- panall[,-c(1,2)]
	panallcombo<-  merge(panall, tlog2sub, by.x="row.names", by.y="row.names")
	rownames(panallcombo)=panallcombo[,1]
	finalpancreascov<- data.frame(panallcombo[,c(56205:56237)])	
	rownames(finalpancreascov)=panallcombo[,1]
	finalpancreasL1HS<- data.frame(panallcombo[,c(56204)])
	rownames(finalpancreasL1HS)= panallcombo[,1]
	colnames(finalpancreasL1HS)= c("L1HS")
	finalpancreasgenes<- data.frame(panallcombo[,c(2:56203)])
	log2sum<- data.frame(panallcombo[,56238])
	colnames(log2sum)=c("log2sum")
	rownames(log2sum)=rownames(panallcombo)
	#prostate
	#Use the old L1HS patients to get correct patient output
	setwd("/home/quints1/Correlations/L1HSpatientedits")
	propats<- read.table("Prostateedits.txt", header=TRUE, sep="\t")
	rownames(propats)=str_replace_all(propats[,1], "-", ".")
	progenes<- merge(propats, tgenesraw, by.x="row.names", by.y="row.names", drop=FALSE)
	progenesfinal<- progenes[,c(4:56207)]
	rownames(progenesfinal)=progenes[,1]
		#read in the covariate
		setwd("/home/quints1/Correlations/prostate")
		procov<- read.table("prostateeditcov.txt", header=TRUE, row.names=1, sep="\t")
	proall<- merge(progenesfinal, procov, by.x="row.names", by.y="row.names")
	rownames(proall)=proall[,1]
	proall<- proall[,-c(1,2)]
	proallcombo<-  merge(proall, tlog2sub, by.x="row.names", by.y="row.names")
	rownames(proallcombo)=proallcombo[,1]
	finalprostatecov<- data.frame(proallcombo[,c(56205:56222)])	
	rownames(finalprostatecov)=proallcombo[,1]
	finalprostateL1HS<- data.frame(proallcombo[,c(56204)])
	rownames(finalprostateL1HS)= proallcombo[,1]
	colnames(finalprostateL1HS)= c("L1HS")
	finalprostategenes<- data.frame(proallcombo[,c(2:56203)])
	log2sum<- data.frame(proallcombo[,56223])
	colnames(log2sum)=c("log2sum")
	rownames(log2sum)=rownames(proallcombo)
	#stomach
		#Use the old L1HS patients to get correct patient output
	setwd("/home/quints1/Correlations/L1HSpatientedits")
	stopats<- read.table("Stomachedits.txt", header=TRUE, sep="\t")
	rownames(stopats)=str_replace_all(stopats[,1], "-", ".")
	stogenes<- merge(stopats, tgenesraw, by.x="row.names", by.y="row.names", drop=FALSE)
	stogenesfinal<- stogenes[,c(4:56207)]
	rownames(stogenesfinal)=stogenes[,1]
		#read in the covariate
		setwd("/home/quints1/Correlations/stomach")
		stocov<- read.table("stomacheditcov.txt", header=TRUE, row.names=1, sep="\t")
	stoall<- merge(stogenesfinal, stocov, by.x="row.names", by.y="row.names")
	rownames(stoall)=stoall[,1]
	stoall<- stoall[,-c(1,2)]
	stoallcombo<-  merge(stoall, tlog2sub, by.x="row.names", by.y="row.names")
	rownames(stoallcombo)=stoallcombo[,1]
	finalstomachcov<- data.frame(stoallcombo[,c(56205:56237)])	
	rownames(finalstomachcov)=stoallcombo[,1]
	finalstomachL1HS<- data.frame(stoallcombo[,c(56204)])
	rownames(finalstomachL1HS)= stoallcombo[,1]
	colnames(finalstomachL1HS)= c("L1HS")
	finalstomachgenes<- data.frame(stoallcombo[,c(2:56203)])
	log2sum<- data.frame(stoallcombo[,56238])
	colnames(log2sum)=c("log2sum")
	rownames(log2sum)=rownames(stoallcombo)
	#thyroid
		#Use the old L1HS patients to get correct patient output
	setwd("/home/quints1/Correlations/L1HSpatientedits")
	thypats<- read.table("Thyroidedits.txt", header=TRUE, sep="\t")
	rownames(thypats)=str_replace_all(thypats[,1], "-", ".")
	thygenes<- merge(thypats, tgenesraw, by.x="row.names", by.y="row.names", drop=FALSE)
	thygenesfinal<- thygenes[,c(4:56207)]
	rownames(thygenesfinal)=thygenes[,1]
		#read in the covariate
		setwd("/home/quints1/Correlations/thyroid")
		thycov<- read.table("thyroideditcov.txt", header=TRUE, row.names=1, sep="\t")
	thyall<- merge(thygenesfinal, thycov, by.x="row.names", by.y="row.names")
	rownames(thyall)=thyall[,1]
	thyall<- thyall[,-c(1,2)]
	thyallcombo<-  merge(thyall, tlog2sub, by.x="row.names", by.y="row.names")
	rownames(thyallcombo)=thyallcombo[,1]
	finalthyroidcov<- data.frame(thyallcombo[,c(56205:56252)])	
	rownames(finalthyroidcov)=thyallcombo[,1]
	finalthyroidL1HS<- data.frame(thyallcombo[,c(56204)])
	rownames(finalthyroidL1HS)= thyallcombo[,1]
	colnames(finalthyroidL1HS)= c("L1HS")
	finalthyroidgenes<- data.frame(thyallcombo[,c(2:56203)])
	log2sum<- data.frame(thyallcombo[,56253])
	colnames(log2sum)=c("log2sum")
	rownames(log2sum)=rownames(thyallcombo)
	
#iterate across all the tissues where that specific gene is found and plot	
library("ggplot2")
setwd("/home/quints1/Correlations/graphsforoutput/graphs2")
y=NULL
tissues=NULL
myplot=NULL
for (i in rownames(subbed)){
	y<- final[c(i),]
	output<- y[, colSums(is.na(y)) != nrow(y)]
	tissues<- data.frame(colnames(output))
	rownames(tissues)=tissues[,1]
	for (tiss in rownames(tissues)){
		selectedtissue<- tiss
		genes<- get(paste("final", tiss, "genes", sep=""))
		L1HS<- get(paste("final", tiss, "L1HS", sep=""))
		mydata<- data.frame(cbind(genes[,c(i)], L1HS[,1]))
		rownames(mydata)=rownames(genes)
		myplot<- ggplot(mydata, aes(mydata[,1], mydata[,2])) + geom_point(color="black") + geom_smooth(method="lm", se=FALSE, color="dark red") + theme_classic() +
		labs(x=paste(i, "Gene Expression"), y="L1HS Expression", title= paste(tiss, "final result", i))
		ggsave(paste(i, tiss, ".png", sep=""))
		}
	}
