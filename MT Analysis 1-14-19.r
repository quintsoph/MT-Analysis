###STEP 3: Sort for MT genes
#adrenal
library("stringi")
setwd("/home/quints1/Correlations/adrenal")
adrenal<- read.table("PRandQAdrenal.txt", header=TRUE, row.names=1)
adrenalAIC<- read.table("All AICc results.txt", header=TRUE)
adrenalslope<- cbind(rownames(adrenalAIC), adrenalAIC[,3])
genes<-stri_sub(adrenalslope[,1], 0, -9)
rownames(adrenalslope)=genes
adrenal2<- merge(adrenal, adrenalslope, by.x="row.names", by.y="row.names", drop=FALSE)
rownames(adrenal2)=adrenal2[,1]
adrenal2<- adrenal2[,c(2,3,4,5,7)]
adrenalRsq<- read.table("P and Q Adrenal.txt", header=TRUE)
sqvals<- adrenalRsq[,c(1,3)]
rownames(sqvals)=sqvals[,1]
adrenal3<- merge(adrenal2, sqvals, by.x="row.names", by.y="row.names", drop=FALSE)
finaladrenal<- adrenal3[,c(1,2,3,4,5,6,8)]
colnames(finaladrenal)=c("rownames", "pvalues", "pearsonr", "spearmanr", "qvalues", "AICcoeff", "rsqvalues")
test<- ifelse(as.numeric(as.character(finaladrenal[,6])) < 0, TRUE, FALSE)
result<- cbind(finaladrenal, test)
cutadrenalneg<- subset(result, result[,8] == "TRUE", select=c(1,2,3,4,5,6,7,8))
negadrenal<- cutadrenalneg[order(-as.numeric(as.character(cutadrenalneg[,7]))),]
adrenalcutoff<- subset(negadrenal, negadrenal[,7] > 0.5, select=c(1,2,3,4,5,6,7))
adrenalcut2<- adrenalcutoff[order(as.numeric(as.character(adrenalcutoff[,6]))),]
adrenalgenes<- data.frame(adrenalcutoff[c(1)])
colnames(adrenalgenes)=c("adrenal")
rownames(adrenalgenes)= adrenalgenes[,1]
test<- adrenalcutoff[order(-as.numeric(as.character(adrenalcutoff[,2]))),]

#breast
library("stringi")
setwd("/home/quints1/Correlations/breast")
breast<- read.table("PRandQBreast.txt", header=TRUE, row.names=1)
breastAIC<- read.table("All AICc results breast.txt", header=TRUE)
breastslope<- cbind(rownames(breastAIC), breastAIC[,3])
genes<-stri_sub(breastslope[,1], 0, -9)
rownames(breastslope)=genes
breast2<- merge(breast, breastslope, by.x="row.names", by.y="row.names", drop=FALSE)
rownames(breast2)=breast2[,1]
breast2<- breast2[,c(2,3,4,5,7)]
breastRsq<- read.table("P and Q Breast.txt", header=TRUE)
sqvals<- breastRsq[,c(1,3)]
rownames(sqvals)=sqvals[,1]
breast3<- merge(breast2, sqvals, by.x="row.names", by.y="row.names", drop=FALSE)
finalbreast<- breast3[,c(1,2,3,4,5,6,8)]
colnames(finalbreast)=c("rownames", "pvalues", "pearsonr", "spearmanr", "qvalues", "AICcoeff", "rsqvalues")
test<- ifelse(as.numeric(as.character(finalbreast[,6])) < 0, TRUE, FALSE)
result<- cbind(finalbreast, test)
cutbreastneg<- subset(result, result[,8] == "TRUE", select=c(1,2,3,4,5,6,7,8))
negbreast<- cutbreastneg[order(-as.numeric(as.character(cutbreastneg[,7]))),]
breastcutoff<- subset(negbreast, negbreast[,7] > 0.5, select=c(1,2,3,4,5,6,7))
breastgenes<- data.frame(breastcutoff[c(1)])	
colnames(breastgenes)=c("breast")
rownames(breastgenes)= breastgenes[,1]

#colonsig
setwd("/home/quints1/Correlations/colonsig")
colonsig<- read.table("PRandQColonSig.txt", header=TRUE, row.names=1)
colonsigAIC<- read.table("AllAICcresultscolonsig.txt", header=TRUE)
colonsigslope<- cbind(rownames(colonsigAIC), colonsigAIC[,3])
genes<-stri_sub(colonsigslope[,1], 0, -9)
rownames(colonsigslope)=genes
colonsig2<- merge(colonsig, colonsigslope, by.x="row.names", by.y="row.names", drop=FALSE)
rownames(colonsig2)=colonsig2[,1]
colonsig2<- colonsig2[,c(2,3,4,5,7)]
colonsigRsq<- read.table("P and Q colonsig.txt", header=TRUE)
sqvals<- colonsigRsq[,c(1,3)]
rownames(sqvals)=sqvals[,1]
colonsig3<- merge(colonsig2, sqvals, by.x="row.names", by.y="row.names", drop=FALSE)
finalcolonsig<- colonsig3[,c(1,2,3,4,5,6,8)]
colnames(finalcolonsig)=c("rownames", "pvalues", "pearsonr", "spearmanr", "qvalues", "AICcoeff", "rsqvalues")
test<- ifelse(as.numeric(as.character(finalcolonsig[,6])) < 0, TRUE, FALSE)
result<- cbind(finalcolonsig, test)
cutcolonsigneg<- subset(result, result[,8] == "TRUE", select=c(1,2,3,4,5,6,7,8))
negcolonsig<- cutcolonsigneg[order(-as.numeric(as.character(cutcolonsigneg[,7]))),]
colonsigcutoff<- subset(negcolonsig, negcolonsig[,7] > 0.5, select=c(1,2,3,4,5,6,7))
colonsiggenes<- data.frame(colonsigcutoff[,1])	
colnames(colonsiggenes)=c("colonsig")
rownames(colonsiggenes)= colonsiggenes[,1]

#colontrans
setwd("/home/quints1/Correlations/colontrans")
colontrans<- read.table("PRandQColonTrans.txt", header=TRUE, row.names=1)
colontransAIC<- read.table("AllAICcresultscolontrans.txt", header=TRUE)
colontransslope<- cbind(rownames(colontransAIC), colontransAIC[,3])
genes<-stri_sub(colontransslope[,1], 0, -9)
rownames(colontransslope)=genes
colontrans2<- merge(colontrans, colontransslope, by.x="row.names", by.y="row.names", drop=FALSE)
rownames(colontrans2)=colontrans2[,1]
colontrans2<- colontrans2[,c(2,3,4,5,7)]
colontransRsq<- read.table("P and Q colontrans.txt", header=TRUE)
sqvals<- colontransRsq[,c(1,3)]
rownames(sqvals)=sqvals[,1]
colontrans3<- merge(colontrans2, sqvals, by.x="row.names", by.y="row.names", drop=FALSE)
finalcolontrans<- colontrans3[,c(1,2,3,4,5,6,8)]
colnames(finalcolontrans)=c("rownames", "pvalues", "pearsonr", "spearmanr", "qvalues", "AICcoeff", "rsqvalues")
test<- ifelse(as.numeric(as.character(finalcolontrans[,6])) < 0, TRUE, FALSE)
result<- cbind(finalcolontrans, test)
cutcolontransneg<- subset(result, result[,8] == "TRUE", select=c(1,2,3,4,5,6,7,8))
negcolontrans<- cutcolontransneg[order(-as.numeric(as.character(cutcolontransneg[,7]))),]
colontranscutoff<- subset(negcolontrans, negcolontrans[,7] > 0.5, select=c(1,2,3,4,5,6,7))
colontransgenes<- data.frame(colontranscutoff[,1])	
colnames(colontransgenes)=c("colontrans")
rownames(colontransgenes)= colontransgenes[,1]

#esophgastro
setwd("/home/quints1/Correlations/esophagusgastro")
esophgastro<- read.table("PRandQEsophGastro.txt", header=TRUE, row.names=1)
esophgastroAIC<- read.table("AllAICcresultsesophgastro.txt", header=TRUE)
esophgastroslope<- cbind(rownames(esophgastroAIC), esophgastroAIC[,3])
genes<-stri_sub(esophgastroslope[,1], 0, -9)
rownames(esophgastroslope)=genes
esophgastro2<- merge(esophgastro, esophgastroslope, by.x="row.names", by.y="row.names", drop=FALSE)
rownames(esophgastro2)=esophgastro2[,1]
esophgastro2<- esophgastro2[,c(2,3,4,5,7)]
esophgastroRsq<- read.table("P and Q esophgastro.txt", header=TRUE)
sqvals<- esophgastroRsq[,c(1,3)]
rownames(sqvals)=sqvals[,1]
esophgastro3<- merge(esophgastro2, sqvals, by.x="row.names", by.y="row.names", drop=FALSE)
finalesophgastro<- esophgastro3[,c(1,2,3,4,5,6,8)]
colnames(finalesophgastro)=c("rownames", "pvalues", "pearsonr", "spearmanr", "qvalues", "AICcoeff", "rsqvalues")
test<- ifelse(as.numeric(as.character(finalesophgastro[,6])) < 0, TRUE, FALSE)
result<- cbind(finalesophgastro, test)
cutesophgastroneg<- subset(result, result[,8] == "TRUE", select=c(1,2,3,4,5,6,7,8))
negesophgastro<- cutesophgastroneg[order(-as.numeric(as.character(cutesophgastroneg[,7]))),]
esophgastrocutoff<- subset(negesophgastro, negesophgastro[,7] > 0.5, select=c(1,2,3,4,5,6,7))
esophgastrogenes<- data.frame(esophgastrocutoff[,1])	
colnames(esophgastrogenes)=c("esophgastro")
rownames(esophgastrogenes)= esophgastrogenes[,1]

#esophmucosa
setwd("/home/quints1/Correlations/esophagusmucosa")
esophmucosa<- read.table("PRandQEsophMucosa.txt", header=TRUE, row.names=1)
esophmucosaAIC<- read.table("AllAICcresultsesophmucosa.txt", header=TRUE)
esophmucosaslope<- cbind(rownames(esophmucosaAIC), esophmucosaAIC[,3])
genes<-stri_sub(esophmucosaslope[,1], 0, -9)
rownames(esophmucosaslope)=genes
esophmucosa2<- merge(esophmucosa, esophmucosaslope, by.x="row.names", by.y="row.names", drop=FALSE)
rownames(esophmucosa2)=esophmucosa2[,1]
esophmucosa2<- esophmucosa2[,c(2,3,4,5,7)]
esophmucosaRsq<- read.table("P and Q esophmucosa.txt", header=TRUE)
sqvals<- esophmucosaRsq[,c(1,3)]
rownames(sqvals)=sqvals[,1]
esophmucosa3<- merge(esophmucosa2, sqvals, by.x="row.names", by.y="row.names", drop=FALSE)
finalesophmucosa<- esophmucosa3[,c(1,2,3,4,5,6,8)]
colnames(finalesophmucosa)=c("rownames", "pvalues", "pearsonr", "spearmanr", "qvalues", "AICcoeff", "rsqvalues")
test<- ifelse(as.numeric(as.character(finalesophmucosa[,6])) < 0, TRUE, FALSE)
result<- cbind(finalesophmucosa, test)
cutesophmucosaneg<- subset(result, result[,8] == "TRUE", select=c(1,2,3,4,5,6,7,8))
negesophmucosa<- cutesophmucosaneg[order(-as.numeric(as.character(cutesophmucosaneg[,7]))),]
esophmucosacutoff<- subset(negesophmucosa, negesophmucosa[,7] > 0.5, select=c(1,2,3,4,5,6,7))
esophmucosagenes<- data.frame(esophmucosacutoff[,1])	
colnames(esophmucosagenes)=c("esophmucosa")
rownames(esophmucosagenes)= esophmucosagenes[,1]

#esophmuscularis
setwd("/home/quints1/Correlations/esophagusmuscularis")
esophmuscularis<- read.table("PRandQEsophMuscularis.txt", header=TRUE, row.names=1)
esophmuscularisAIC<- read.table("AllAICcresultsesophmuscularis.txt", header=TRUE)
esophmuscularisslope<- cbind(rownames(esophmuscularisAIC), esophmuscularisAIC[,3])
genes<-stri_sub(esophmuscularisslope[,1], 0, -9)
rownames(esophmuscularisslope)=genes
esophmuscularis2<- merge(esophmuscularis, esophmuscularisslope, by.x="row.names", by.y="row.names", drop=FALSE)
rownames(esophmuscularis2)=esophmuscularis2[,1]
esophmuscularis2<- esophmuscularis2[,c(2,3,4,5,7)]
esophmuscularisRsq<- read.table("P and Q esophmuscularis.txt", header=TRUE)
sqvals<- esophmuscularisRsq[,c(1,3)]
rownames(sqvals)=sqvals[,1]
esophmuscularis3<- merge(esophmuscularis2, sqvals, by.x="row.names", by.y="row.names", drop=FALSE)
finalesophmuscularis<- esophmuscularis3[,c(1,2,3,4,5,6,8)]
colnames(finalesophmuscularis)=c("rownames", "pvalues", "pearsonr", "spearmanr", "qvalues", "AICcoeff", "rsqvalues")
test<- ifelse(as.numeric(as.character(finalesophmuscularis[,6])) < 0, TRUE, FALSE)
result<- cbind(finalesophmuscularis, test)
cutesophmuscularisneg<- subset(result, result[,8] == "TRUE", select=c(1,2,3,4,5,6,7,8))
negesophmuscularis<- cutesophmuscularisneg[order(-as.numeric(as.character(cutesophmuscularisneg[,7]))),]
esophmusculariscutoff<- subset(negesophmuscularis, negesophmuscularis[,7] > 0.5, select=c(1,2,3,4,5,6,7))
esophmuscularisgenes<- data.frame(esophmusculariscutoff[,1])	
colnames(esophmuscularisgenes)=c("esophmuscularis")
rownames(esophmuscularisgenes)= esophmuscularisgenes[,1]

#kidney
setwd("/home/quints1/Correlations/kidney")
kidney<- read.table("PRandQKidney.txt", header=TRUE, row.names=1)
kidneyAIC<- read.table("AllAICcresultskidney.txt", header=TRUE)
kidneyslope<- cbind(rownames(kidneyAIC), kidneyAIC[,3])
genes<-stri_sub(kidneyslope[,1], 0, -9)
rownames(kidneyslope)=genes
kidney2<- merge(kidney, kidneyslope, by.x="row.names", by.y="row.names", drop=FALSE)
rownames(kidney2)=kidney2[,1]
kidney2<- kidney2[,c(2,3,4,5,7)]
kidneyRsq<- read.table("P and Q kidney.txt", header=TRUE)
sqvals<- kidneyRsq[,c(1,3)]
rownames(sqvals)=sqvals[,1]
kidney3<- merge(kidney2, sqvals, by.x="row.names", by.y="row.names", drop=FALSE)
finalkidney<- kidney3[,c(1,2,3,4,5,6,8)]
colnames(finalkidney)=c("rownames", "pvalues", "pearsonr", "spearmanr", "qvalues", "AICcoeff", "rsqvalues")
test<- ifelse(as.numeric(as.character(finalkidney[,6])) < 0, TRUE, FALSE)
result<- cbind(finalkidney, test)
cutkidneyneg<- subset(result, result[,8] == "TRUE", select=c(1,2,3,4,5,6,7,8))
negkidney<- cutkidneyneg[order(-as.numeric(as.character(cutkidneyneg[,7]))),]
kidneycutoff<- subset(negkidney, negkidney[,7] > 0.5, select=c(1,2,3,4,5,6,7))
kidneygenes<- data.frame(kidneycutoff[,1])	
colnames(kidneygenes)=c("kidney")
rownames(kidneygenes)= kidneygenes[,1]

#liver
setwd("/home/quints1/Correlations/liver")
liver<- read.table("PRandQLiver.txt", header=TRUE, row.names=1)
liverAIC<- read.table("AllAICcresultsliver.txt", header=TRUE)
liverslope<- cbind(rownames(liverAIC), liverAIC[,3])
genes<-stri_sub(liverslope[,1], 0, -9)
rownames(liverslope)=genes
liver2<- merge(liver, liverslope, by.x="row.names", by.y="row.names", drop=FALSE)
rownames(liver2)=liver2[,1]
liver2<- liver2[,c(2,3,4,5,7)]
liverRsq<- read.table("P and Q liver.txt", header=TRUE)
sqvals<- liverRsq[,c(1,3)]
rownames(sqvals)=sqvals[,1]
liver3<- merge(liver2, sqvals, by.x="row.names", by.y="row.names", drop=FALSE)
finalliver<- liver3[,c(1,2,3,4,5,6,8)]
colnames(finalliver)=c("rownames", "pvalues", "pearsonr", "spearmanr", "qvalues", "AICcoeff", "rsqvalues")
test<- ifelse(as.numeric(as.character(finalliver[,6])) < 0, TRUE, FALSE)
result<- cbind(finalliver, test)
cutliverneg<- subset(result, result[,8] == "TRUE", select=c(1,2,3,4,5,6,7,8))
negliver<- cutliverneg[order(-as.numeric(as.character(cutliverneg[,7]))),]
livercutoff<- subset(negliver, negliver[,7] > 0.5, select=c(1,2,3,4,5,6,7))
livergenes<- data.frame(livercutoff[,1])	
colnames(livergenes)=c("liver")
rownames(livergenes)= livergenes[,1]

#pancreas
setwd("/home/quints1/Correlations/pancreas")
pancreas<- read.table("PRandQPancreas.txt", header=TRUE, row.names=1)
pancreasAIC<- read.table("AllAICcresultspancreas.txt", header=TRUE)
pancreasslope<- cbind(rownames(pancreasAIC), pancreasAIC[,3])
genes<-stri_sub(pancreasslope[,1], 0, -9)
rownames(pancreasslope)=genes
pancreas2<- merge(pancreas, pancreasslope, by.x="row.names", by.y="row.names", drop=FALSE)
rownames(pancreas2)=pancreas2[,1]
pancreas2<- pancreas2[,c(2,3,4,5,7)]
pancreasRsq<- read.table("P and Q pancreas.txt", header=TRUE)
sqvals<- pancreasRsq[,c(1,3)]
rownames(sqvals)=sqvals[,1]
pancreas3<- merge(pancreas2, sqvals, by.x="row.names", by.y="row.names", drop=FALSE)
finalpancreas<- pancreas3[,c(1,2,3,4,5,6,8)]
colnames(finalpancreas)=c("rownames", "pvalues", "pearsonr", "spearmanr", "qvalues", "AICcoeff", "rsqvalues")
test<- ifelse(as.numeric(as.character(finalpancreas[,6])) < 0, TRUE, FALSE)
result<- cbind(finalpancreas, test)
cutpancreasneg<- subset(result, result[,8] == "TRUE", select=c(1,2,3,4,5,6,7,8))
negpancreas<- cutpancreasneg[order(-as.numeric(as.character(cutpancreasneg[,7]))),]
pancreascutoff<- subset(negpancreas, negpancreas[,7] > 0.5, select=c(1,2,3,4,5,6,7))
pancreasgenes<- data.frame(pancreascutoff[,1])	
colnames(pancreasgenes)=c("pancreas")
rownames(pancreasgenes)= pancreasgenes[,1]

#prostate
setwd("/home/quints1/Correlations/prostate")
prostate<- read.table("PRandQProstate.txt", header=TRUE, row.names=1)
prostateAIC<- read.table("AllAICcresultprostate.txt", header=TRUE)
prostateslope<- cbind(rownames(prostateAIC), prostateAIC[,3])
genes<-stri_sub(prostateslope[,1], 0, -9)
rownames(prostateslope)=genes
prostate2<- merge(prostate, prostateslope, by.x="row.names", by.y="row.names", drop=FALSE)
rownames(prostate2)=prostate2[,1]
prostate2<- prostate2[,c(2,3,4,5,7)]
prostateRsq<- read.table("P and Q prostate.txt", header=TRUE)
sqvals<- prostateRsq[,c(1,3)]
rownames(sqvals)=sqvals[,1]
prostate3<- merge(prostate2, sqvals, by.x="row.names", by.y="row.names", drop=FALSE)
finalprostate<- prostate3[,c(1,2,3,4,5,6,8)]
colnames(finalprostate)=c("rownames", "pvalues", "pearsonr", "spearmanr", "qvalues", "AICcoeff", "rsqvalues")
test<- ifelse(as.numeric(as.character(finalprostate[,6])) < 0, TRUE, FALSE)
result<- cbind(finalprostate, test)
cutprostateneg<- subset(result, result[,8] == "TRUE", select=c(1,2,3,4,5,6,7,8))
negprostate<- cutprostateneg[order(-as.numeric(as.character(cutprostateneg[,7]))),]
prostatecutoff<- subset(negprostate, negprostate[,7] > 0.5, select=c(1,2,3,4,5,6,7))
prostategenes<- data.frame(prostatecutoff[,1])	
colnames(prostategenes)=c("prostate")
rownames(prostategenes)= prostategenes[,1]

#stomach
setwd("/home/quints1/Correlations/stomach")
stomach<- read.table("PRandQStomach.txt", header=TRUE, row.names=1)
stomachAIC<- read.table("AllAICcresultsstomach.txt", header=TRUE)
stomachslope<- cbind(rownames(stomachAIC), stomachAIC[,3])
genes<-stri_sub(stomachslope[,1], 0, -9)
rownames(stomachslope)=genes
stomach2<- merge(stomach, stomachslope, by.x="row.names", by.y="row.names", drop=FALSE)
rownames(stomach2)=stomach2[,1]
stomach2<- stomach2[,c(2,3,4,5,7)]
stomachRsq<- read.table("P and Q stomach.txt", header=TRUE)
sqvals<- stomachRsq[,c(1,3)]
rownames(sqvals)=sqvals[,1]
stomach3<- merge(stomach2, sqvals, by.x="row.names", by.y="row.names", drop=FALSE)
finalstomach<- stomach3[,c(1,2,3,4,5,6,8)]
colnames(finalstomach)=c("rownames", "pvalues", "pearsonr", "spearmanr", "qvalues", "AICcoeff", "rsqvalues")
test<- ifelse(as.numeric(as.character(finalstomach[,6])) < 0, TRUE, FALSE)
result<- cbind(finalstomach, test)
cutstomachneg<- subset(result, result[,8] == "TRUE", select=c(1,2,3,4,5,6,7,8))
negstomach<- cutstomachneg[order(-as.numeric(as.character(cutstomachneg[,7]))),]
stomachcutoff<- subset(negstomach, negstomach[,7] > 0.5, select=c(1,2,3,4,5,6,7))
stomachgenes<- data.frame(stomachcutoff[,1])
colnames(stomachgenes)=c("stomach")
rownames(stomachgenes)= stomachgenes[,1]

#thyroid
setwd("/home/quints1/Correlations/thyroid")
thyroid<- read.table("PRandQThyroid.txt", header=TRUE, row.names=1)
thyroidAIC<- read.table("AllAICcresultsthyroid.txt", header=TRUE)
thyroidslope<- cbind(rownames(thyroidAIC), thyroidAIC[,3])
genes<-stri_sub(thyroidslope[,1], 0, -9)
rownames(thyroidslope)=genes
thyroid2<- merge(thyroid, thyroidslope, by.x="row.names", by.y="row.names", drop=FALSE)
rownames(thyroid2)=thyroid2[,1]
thyroid2<- thyroid2[,c(2,3,4,5,7)]
thyroidRsq<- read.table("P and Q thyroid.txt", header=TRUE)
sqvals<- thyroidRsq[,c(1,3)]
rownames(sqvals)=sqvals[,1]
thyroid3<- merge(thyroid2, sqvals, by.x="row.names", by.y="row.names", drop=FALSE)
finalthyroid<- thyroid3[,c(1,2,3,4,5,6,8)]
colnames(finalthyroid)=c("rownames", "pvalues", "pearsonr", "spearmanr", "qvalues", "AICcoeff", "rsqvalues")
test<- ifelse(as.numeric(as.character(finalthyroid[,6])) < 0, TRUE, FALSE)
result<- cbind(finalthyroid, test)
cutthyroidneg<- subset(result, result[,8] == "TRUE", select=c(1,2,3,4,5,6,7,8))
negthyroid<- cutthyroidneg[order(-as.numeric(as.character(cutthyroidneg[,7]))),]
thyroidcutoff<- subset(negthyroid, negthyroid[,7] > 0.5, select=c(1,2,3,4,5,6,7))
thyroidgenes<- data.frame(thyroidcutoff[,1])
colnames(thyroidgenes)=c("thyroid")
rownames(thyroidgenes)= thyroidgenes[,1]

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
test<- ifelse(mylist2[,1] > 2, TRUE, FALSE)
combined<- cbind(mylist2, test)
subbed<- subset(combined, combined[,2] == "TRUE", select=c(1,2))
saver<- rownames(subbed)
setwd("/home/quints1/Correlations/")
write.table(saver, "genesin2tissues redo.txt", quote=FALSE, col.names=FALSE, row.names=FALSE, sep="\t")

#make graphs of the adrenal and breast
topadrenal<- adrenalcutoff[1,]
bestmodel<- adrenalAIC[c("PHBP8"),]
#bestmodel=valueall

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
	library("stringi")
	
