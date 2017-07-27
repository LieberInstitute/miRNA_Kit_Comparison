########################Reproducibility analysis of miRNA expression of 4 Kits#####################
library(data.table)

###Accuracyfullset###
Accur_mir_RPM<-read.table("/home/carrie/miRge/Hiseq_newrun_Accuracy/miR.RPM.csv", header = T, sep = ",")
rownames(Accur_mir_RPM)<- Accur_mir_RPM$miRNA
Accur_mir_RPM<-Accur_mir_RPM[,2:length(colnames(Accur_mir_RPM))]
threekitsAccur<-Accur_mir_RPM
###Reprofullset###
Repro_mir_RPM<-read.table("/home/carrie/miRge/Hiseq_newrun_Repro/miR.RPM.csv", header = T, sep = ",")
rownames(Repro_mir_RPM)<- Repro_mir_RPM$miRNA
Repro_mir_RPM<-Repro_mir_RPM[,2:length(colnames(Repro_mir_RPM))]
threekitsRepro<-Repro_mir_RPM
#ClontechAccur
Accur_mir_RPM<-read.table("/home/carrie/miRge/Hiseq_newrun_Accuracy/clontech/miR.RPM.csv", header = T, sep = ",")
rownames(Accur_mir_RPM)<- Accur_mir_RPM$miRNA
Accur_mir_RPM<-Accur_mir_RPM[,2:length(colnames(Accur_mir_RPM))]
clontechAccur<-Accur_mir_RPM
#ClontechRepro
Repro_mir_RPM<-read.table("/home/carrie/miRge/Hiseq_newrun_Repro/clontech/miR.RPM.csv", header = T, sep = ",")
rownames(Repro_mir_RPM)<- Repro_mir_RPM$miRNA
Repro_mir_RPM<-Repro_mir_RPM[,2:length(colnames(Repro_mir_RPM))]
clontechRepro<-Repro_mir_RPM

#remove synthetic for later tests
Nth.delete<-function(dataframe, n)dataframe[-(seq(n,to=ncol(dataframe),by=n))]
onemgAcc<- Nth.delete(threekitsAccur, 4)
#onemgAcc_clontech<- Nth.delete(clontechAccur, 4)

Repro_mir_RPM<- data.frame(clontechRepro, threekitsRepro)


###Reproducibility within experiment
Accur_mir_RPM<-onemgAcc 
clontechA<-clontechAccur#clontech was run in miRge separately (needed to be redone)... it has about 30 more miRNAs

save(Repro_mir_RPM, Accur_mir_RPM, clontechA, file = "/home/carrie/miRge/Reproducibility4kits.rda")
#Filtering Each Kit for Expression  with RPM>0
clontechA<- clontechAccur[1:3]
row_sub = which(rowMeans(clontechA)>0)
clontechA<-clontechA[row_sub,]
IllA<- Accur_mir_RPM[1:3]
row_sub = which(rowMeans(IllA)>0)
IllA<-IllA[row_sub,]
NEBA<- Accur_mir_RPM[4:6]
row_sub = which(rowMeans(NEBA)>0)
NEBA<-NEBA[row_sub,]
NEXTA<- Accur_mir_RPM[7:9]
row_sub = which(rowMeans(NEXTA)>0)
NEXTA<-NEXTA[row_sub,]

hist(rowMeans(clontechA))
hist(rowMeans(IllA))
hist(rowMeans(NEBA))
hist(rowMeans(NEXTA))

#Filtering Each Kit for Expression  with RPM>10
clontechA<- clontechAccur[1:3]
row_sub = which(rowMeans(clontechA)>10)
clontechA<-clontechA[row_sub,]
IllA<- Accur_mir_RPM[1:3]
row_sub = which(rowMeans(IllA)>10)
IllA<-IllA[row_sub,]
NEBA<- Accur_mir_RPM[4:6]
row_sub = which(rowMeans(NEBA)>10)
NEBA<-NEBA[row_sub,]
NEXTA<- Accur_mir_RPM[7:9]
row_sub = which(rowMeans(NEXTA)>10)
NEXTA<-NEXTA[row_sub,]

hist(rowMeans(clontechA))
hist(rowMeans(IllA))
hist(rowMeans(NEBA))
hist(rowMeans(NEXTA))


###log 2 normalization
clontechA<- log2(clontechA +1)
IllA<-log2(IllA +1)
NEBA<-log2(NEBA +1)
NEXTA<-log2(NEXTA +1)

hist(rowMeans(clontechA))
hist(rowMeans(IllA))
hist(rowMeans(NEBA))
hist(rowMeans(NEXTA))


##number of detected miRNAs
dim(clontechA)
dim(IllA)
dim(NEBA)
dim(NEXTA)


###Test the correlation between each triplicate for each kit
Tests <- list()
Tests[[1]]<-cor.test(clontechA$Clontech_acc_trimmed.1_R1.fq, clontechA$Clontech_acc_trimmed.2_R1.fq)
Tests[[2]]<-cor.test(clontechA$Clontech_acc_trimmed.1_R1.fq, clontechA$Clontech_acc_trimmed.3_R1.fq)
Tests[[3]]<-cor.test(clontechA$Clontech_acc_trimmed.3_R1.fq, clontechA$Clontech_acc_trimmed.2_R1.fq)

Tests[[4]]<-cor.test(IllA$Illumina_acc_trimmed.1_R1.fq, IllA$Illumina_acc_trimmed.2_R1.fq)
Tests[[5]]<-cor.test(IllA$Illumina_acc_trimmed.1_R1.fq, IllA$Illumina_acc_trimmed.3_R1.fq)
Tests[[6]]<-cor.test(IllA$Illumina_acc_trimmed.3_R1.fq, IllA$Illumina_acc_trimmed.2_R1.fq)

Tests[[7]]<-cor.test(NEBA$NEB_acc_trimmed.1_R1.fq, NEBA$NEB_acc_trimmed.2_R1.fq)
Tests[[8]]<-cor.test(NEBA$NEB_acc_trimmed.1_R1.fq, NEBA$NEB_acc_trimmed.3_R1.fq)
Tests[[9]]<-cor.test(NEBA$NEB_acc_trimmed.3_R1.fq, NEBA$NEB_acc_trimmed.1_R1.fq)

Tests[[10]]<-cor.test(NEXTA$NEXT_acc_trimmed.1_R1.fq, NEXTA$NEXT_acc_trimmed.2_R1.fq)
Tests[[11]]<-cor.test(NEXTA$NEXT_acc_trimmed.1_R1.fq, NEXTA$NEXT_acc_trimmed.3_R1.fq)
Tests[[12]]<-cor.test(NEXTA$NEXT_acc_trimmed.3_R1.fq, NEXTA$NEXT_acc_trimmed.2_R1.fq)


# extract your values using `sapply`
Acc_CorrStats<-sapply(Tests, function(x) {
  c(format(x$estimate,digits = 3),
    ci = format(x$conf.int, digits = 3),
    p.value = format(x$p.value, scientific = TRUE))
})

colnames(Acc_CorrStats)<-sapply(Tests, function(x) {
  c(test = x$data.name)
})


boxplot(IllA, lwd = 2, ylab = 'counts across triplicates(log2(RPM +1))', notch = TRUE, outline = FALSE, at =1:3, xlim=c(1,12), ylim=c(0,16), col ='turquoise1', xaxt='n') #xaxt= 'n' gets rid of x labels
axis(1, at =2, labels = c("Illumina"))
boxplot(NEBA, lwd = 2, notch = TRUE, outline = FALSE, add = TRUE, at =4:6, col = 'springgreen1', xaxt='n')
axis(1, at =5, labels = c("NEB"))
boxplot(NEXTA, lwd = 2, notch = TRUE, outline = FALSE, add =TRUE, at =7:9, col = "maroon1", xaxt='n')
axis(1, at =8, labels = c("NEXT"))
boxplot(clontechA, lwd = 2, xlab = 'Kits', notch = TRUE, outline = FALSE, add = TRUE, at = 10:12, col = "yellow", xaxt='n')
axis(1, at =11, labels = c("Clontech"))

t(Acc_CorrStats)
test<-as.data.frame(t(Acc_CorrStats))

test$kit<-(c(rep("clontech",3), rep("Ill", 3), rep("NEB",3), rep("NEXT", 3)))
test$trip<-rep(c(1:3),4)
#test$kit<-paste0(rep(c(1:3),4), test$kit)
test$kit=factor(test$kit)
Accur_Corr<-data.frame(test$kit, test$cor) ####displays the data nicely
boxplot(as.numeric(levels(test$cor)[test$cor])~test$kit, ylab = "Correlation")

###ttests of correlation values
test$cor<-as.numeric(levels(test$cor)[test$cor])
ttests<-list()
ttests[[1]]<- t.test(test$cor[grep("Ill",test$kit)],test$cor[grep("NEXT",test$kit)])
ttests[[2]]<- t.test(test$cor[grep("Ill",test$kit)],test$cor[grep("NEB",test$kit)])
ttests[[3]]<- t.test(test$cor[grep("Ill",test$kit)],test$cor[grep("clontech",test$kit)])
ttests[[4]]<- t.test(test$cor[grep("NEXT",test$kit)],test$cor[grep("NEB",test$kit)])
ttests[[5]]<- t.test(test$cor[grep("NEXT",test$kit)],test$cor[grep("clontech",test$kit)])
ttests[[6]]<- t.test(test$cor[grep("NEB",test$kit)],test$cor[grep("clontech",test$kit)])

# extract your values using `sapply`
Acc_CorrStatsttest<-sapply(ttests, function(x) {
  c(format(x$estimate,digits = 3),
    t = format(x$statistic, digits = 3),
    p.value = format(x$p.value, scientific = FALSE, digits = 2),
    bonferroni = format(.05/6, digits = 2),
    sig = ifelse(x$p.value<(.05/6), "yes", "no"))
})

colnames(Acc_CorrStatsttest)<-sapply(ttests, function(x) {
  c(test = x$data.name)
})



IllA<-as.matrix(IllA)
Ill1<-abs(IllA[,1]-IllA[,2])
Ill2<-abs(IllA[,1]-IllA[,3])
Ill3<-abs(IllA[,2]-IllA[,3])
Illresid<-data.frame(Ill1, Ill2, Ill3)
boxplot(Illresid, lwd = 2, ylab = 'absolute value of difference b/w counts across triplicates', notch = TRUE, outline = FALSE, at =1:3, xlim=c(0,12), ylim=c(0,1.3), col ='turquoise1', xaxt ='n')
axis(1, at =2, labels = c("Illumina"))

NEBA<-as.matrix(NEBA)
NEB1<-abs(NEBA[,1]-NEBA[,2])
NEB2<-abs(NEBA[,1]-NEBA[,3])
NEB3<-abs(NEBA[,2]-NEBA[,3])
NEBresid<-data.frame(NEB1, NEB2, NEB3)
boxplot(NEBresid, lwd = 2, notch = TRUE, outline = FALSE, add = TRUE, at =4:6, col = 'springgreen1', xaxt ="n")
axis(1, at = 5, labels = c("NEB"))

NEXTA<-as.matrix(NEXTA)
NEXT1<-abs(NEXTA[,1]-NEXTA[,2])
NEXT2<-abs(NEXTA[,1]-NEXTA[,3])
NEXT3<-abs(NEXTA[,2]-NEXTA[,3])
NEXTresid<-data.frame(NEXT1, NEXT2, NEXT3)
boxplot(NEXTresid, lwd = 2, notch = TRUE, outline = FALSE, add =TRUE, at =7:9, col = "maroon1", xaxt ="n")
axis(1, at = 8, labels = c("NEXT"))

clontechA<-as.matrix(clontechA)
clontechResid1<-data.frame
clontech1<-abs(clontechA[,1]-clontechA[,2])
clontech2<-abs(clontechA[,1]-clontechA[,3])
clontech3<-abs(clontechA[,2]-clontechA[,3])
clontechresid<-data.frame(clontech1, clontech2, clontech3)
boxplot(clontechresid, lwd = 2, notch = TRUE, outline = FALSE, add = TRUE, at = 10:12, col = "yellow", xaxt ="n")
axis(1, at = 11, labels = c("Clontech"))

#resid_accur<-data.frame(clontechresid, Illresid, NEBresid, NEXTresid)

###########in the Repro dataset##############
Repro_mir_RPM<- data.frame(clontechRepro, threekitsRepro)
#Repro_mir_RPM <- log2(Repro_mir_RPM + 1)
temp<-colnames(Repro_mir_RPM)
temp<-gsub("_R1.fq", ".fq", temp)
colnames(Repro_mir_RPM)<-temp

#500ngdata####################################################################
Repro_mir_RPM<- data.frame(clontechRepro, threekitsRepro)
clontechR<- Repro_mir_RPM[16:18]#7:9
row_sub = which(rowMeans(clontechR)>10)
clontechR<-clontechR[row_sub,]
NEBR<- Repro_mir_RPM[37:39]#7:9
row_sub = which(rowMeans(NEBR)>10)
NEBR<-NEBR[row_sub,]
NEXTR<- Repro_mir_RPM[55:57]#7:9
row_sub = which(rowMeans(NEXTR)>10)
NEXTR<-NEXTR[row_sub ,]

clontechR<- log2(clontechR +1)
NEBR<-log2(NEBR +1)
NEXTR<-log2(NEXTR +1)

dim(clontechR)
dim(NEBR)
dim(NEXTR)

Tests <- list()
Tests[[1]]<-cor.test(clontechR$Clontech_trimmed.7_R1.fq, clontechR$Clontech_trimmed.8_R1.fq)
Tests[[2]]<-cor.test(clontechR$Clontech_trimmed.7_R1.fq, clontechR$Clontech_trimmed.9_R1.fq)
Tests[[3]]<-cor.test(clontechR$Clontech_trimmed.9_R1.fq, clontechR$Clontech_trimmed.8_R1.fq)

Tests[[4]]<-cor.test(NEBR$NEB_trimmed.7_R1.fq, NEBR$NEB_trimmed.8_R1.fq)
Tests[[5]]<-cor.test(NEBR$NEB_trimmed.7_R1.fq, NEBR$NEB_trimmed.9_R1.fq)
Tests[[6]]<-cor.test(NEBR$NEB_trimmed.9_R1.fq, NEBR$NEB_trimmed.8_R1.fq)

Tests[[7]]<-cor.test(NEXTR$NEXT_trimmed.7_R1.fq, NEXTR$NEXT_trimmed.8_R1.fq)
Tests[[8]]<-cor.test(NEXTR$NEXT_trimmed.7_R1.fq, NEXTR$NEXT_trimmed.9_R1.fq)
Tests[[9]]<-cor.test(NEXTR$NEXT_trimmed.9_R1.fq, NEXTR$NEXT_trimmed.8_R1.fq)

Repro_CorrStats<-sapply(Tests, function(x) {
  c(format(x$estimate[1], digits = 2),
    ci.lower = format(x$conf.int[1], digits = 2),
    ci.upper = format(x$conf.int[2], digits = 2),
    p.value = format(x$p.value, scientific = TRUE, digits = 2))
})
colnames(Repro_CorrStats)<-sapply(Tests, function(x) {
  c(test = x$data.name)
})


t(Repro_CorrStats)

boxplot(NEBR, lwd = 2, notch = TRUE, ylab = 'counts across triplicates(log2(RPM +1))', xlim =c(1, 9), outline = FALSE, at =1:3, col = 'springgreen1', xaxt='n')
axis(1, at =2, labels = c("NEB"))
boxplot(NEXTR, lwd = 2, notch = TRUE, outline = FALSE, add =TRUE, at =4:6, col = "maroon1", xaxt='n')
axis(1, at =5, labels = c("NEXT"))
boxplot(clontechR, lwd = 2, xlab = 'Kits', notch = TRUE, outline = FALSE, add = TRUE, at = 7:9, col = "yellow", xaxt='n')
axis(1, at =8, labels = c("Clontech"))

test<-as.data.frame(t(Repro_CorrStats))
test$kit<-(c(rep("clontech", 3), rep("NEB",3), rep("NEXT", 3)))
test$trip<-rep(c(1:3),3)
#test$kit<-paste0(rep(c(1:3),4), test$kit)
test$kit=factor(test$kit)
Accur_Corr<-data.frame(test$kit, test$cor)
boxplot(as.numeric(levels(test$cor)[test$cor])~test$kit, ylab = "Correlation")
test$cor<-as.numeric(levels(test$cor)[test$cor])
ttests<-list()
ttests[[1]]<-t.test(test$cor[grep("NEXT",test$kit)],test$cor[grep("NEB",test$kit)])
ttests[[2]]<-t.test(test$cor[grep("NEXT",test$kit)],test$cor[grep("clontech",test$kit)])
ttests[[3]]<-t.test(test$cor[grep("NEB",test$kit)],test$cor[grep("clontech",test$kit)])
# extract your values using `sapply`
Repro_CorrStatsttest<-sapply(ttests, function(x) {
  c(format(x$estimate,digits = 3),
    t = format(x$statistic, digits = 3),
    p.value = format(x$p.value, scientific = FALSE, digits = 2),
    bonferroni = format(.05/6, digits = 2),
    sig = ifelse(x$p.value<(.05/6), "yes", "no"))
})

colnames(Repro_CorrStatsttest)<-sapply(ttests, function(x) {
  c(test = x$data.name)
})

t(Repro_CorrStatsttest)


NEBR<-as.matrix(NEBR)
NEB1<-abs(NEBR[,1]-NEBR[,2])
NEB2<-abs(NEBR[,1]-NEBR[,3])
NEB3<-abs(NEBR[,2]-NEBR[,3])
NEBresid<-data.frame(NEB1, NEB2, NEB3)
boxplot(NEBresid, lwd = 2, ylab = 'absolute value of difference b/w counts across triplicates', notch = TRUE, outline = FALSE, at =1:3, xlim=c(0,10), ylim = c(0,1.8), col = 'springgreen1', xaxt ="n")
axis(1, at = 2, labels = c("NEB"))

NEXTR<-as.matrix(NEXTR)
NEXT1<-abs(NEXTR[,1]-NEXTR[,2])
NEXT2<-abs(NEXTR[,1]-NEXTR[,3])
NEXT3<-abs(NEXTR[,2]-NEXTR[,3])
NEXTresid<-data.frame(NEXT1, NEXT2, NEXT3)
boxplot(NEXTresid, lwd = 2, notch = TRUE, outline = FALSE, add =TRUE, at =4:6, col = "maroon1", xaxt ="n")
axis(1, at = 5, labels = c("NEXT"))

clontechA<-as.matrix(clontechA)
clontechResid1<-data.frame
clontech1<-abs(clontechR[,1]-clontechR[,2])
clontech2<-abs(clontechR[,1]-clontechR[,3])
clontech3<-abs(clontechR[,2]-clontechR[,3])
clontechresid<-data.frame(clontech1, clontech2, clontech3)
boxplot(clontechresid, lwd = 2, notch = TRUE, outline = FALSE, add = TRUE, at = 7:9, col = "yellow", xaxt ="n")
axis(1, at = 8, labels = c("Clontech"))


#1 microgram data####################################################################
Repro_mir_RPM<- data.frame(clontechRepro, threekitsRepro)
clontechR<- Repro_mir_RPM[1:3]#10:12
row_sub = which(rowMeans(clontechR)>10)
clontechR<-clontechR[row_sub,]
IllR<- Repro_mir_RPM[19:21]#1:3
row_sub = which(rowMeans(IllR)>10)
IllR<-IllR[row_sub,]
NEBR<- Repro_mir_RPM[28:30]#10:12
row_sub = which(rowMeans(NEBR)>10)
NEBR<-NEBR[row_sub,]
NEXTR<- Repro_mir_RPM[40:42]#10:12
row_sub = which(rowMeans(NEXTR)>10)
NEXTR<-NEXTR[row_sub ,]

clontechR<- log2(clontechR +1)
IllR<-log2(IllR +1)
NEBR<-log2(NEBR +1)
NEXTR<-log2(NEXTR +1)

dim(clontechR)
dim(IllR)
dim(NEBR)
dim(NEXTR)


Tests <- list()

Tests[[1]]<-cor.test(clontechR$Clontech_trimmed.10_R1.fq, clontechR$Clontech_trimmed.11_R1.fq)
Tests[[2]]<-cor.test(clontechR$Clontech_trimmed.10_R1.fq, clontechR$Clontech_trimmed.12_R1.fq)
Tests[[3]]<-cor.test(clontechR$Clontech_trimmed.12_R1.fq, clontechR$Clontech_trimmed.11_R1.fq)

Tests[[4]]<-cor.test(IllR$Illumina_trimmed.1_R1.fq, IllR$Illumina_trimmed.2_R1.fq)
Tests[[5]]<-cor.test(IllR$Illumina_trimmed.1_R1.fq, IllR$Illumina_trimmed.3_R1.fq)
Tests[[6]]<-cor.test(IllR$Illumina_trimmed.3_R1.fq, IllR$Illumina_trimmed.2_R1.fq)

Tests[[7]]<-cor.test(NEBR$NEB_trimmed.10_R1.fq, NEBR$NEB_trimmed.11_R1.fq)
Tests[[8]]<-cor.test(NEBR$NEB_trimmed.10_R1.fq, NEBR$NEB_trimmed.12_R1.fq)
Tests[[9]]<-cor.test(NEBR$NEB_trimmed.12_R1.fq , NEBR$NEB_trimmed.11_R1.fq)

Tests[[10]]<-cor.test(NEXTR$NEXT_trimmed.10_R1.fq, NEXTR$NEXT_trimmed.11_R1.fq)
Tests[[11]]<-cor.test(NEXTR$NEXT_trimmed.10_R1.fq, NEXTR$NEXT_trimmed.12_R1.fq)
Tests[[12]]<-cor.test(NEXTR$NEXT_trimmed.12_R1.fq, NEXTR$NEXT_trimmed.11_R1.fq)


Repro_CorrStats<-sapply(Tests, function(x) {
  c(x$estimate[1],
    ci.lower = x$conf.int[1],
    ci.upper = x$conf.int[2],
    p.value = format(x$p.value, scientific = TRUE))
})
t(Repro_CorrStats)

boxplot(IllR, lwd = 2, ylab = 'counts across triplicates(log2(RPM +1))', notch = TRUE, outline = FALSE, at =1:3, xlim=c(1,12), ylim=c(0,16), col ='turquoise1', xaxt='n') #xaxt= 'n' gets rid of x labels
axis(1, at =2, labels = c("Illumina"))
boxplot(NEBR, lwd = 2, notch = TRUE, outline = FALSE, add = TRUE, at =4:6, col = 'springgreen1', xaxt='n')
axis(1, at =5, labels = c("NEB"))
boxplot(NEXTR, lwd = 2, notch = TRUE, outline = FALSE, add =TRUE, at =7:9, col = "maroon1", xaxt='n')
axis(1, at =8, labels = c("NEXT"))
boxplot(clontechR, lwd = 2, xlab = 'Kits', notch = TRUE, outline = FALSE, add = TRUE, at = 10:12, col = "yellow", xaxt='n')
axis(1, at =11, labels = c("Clontech"))

test<-as.data.frame(t(Repro_CorrStats))

test$kit<-(c(rep("clontech",3), rep("Ill", 3), rep("NEB",3), rep("NEXT", 3)))
test$trip<-rep(c(1:3),4)
#test$kit<-paste0(rep(c(1:3),4), test$kit)
test$kit=factor(test$kit)
Accur_Corr<-data.frame(test$kit, test$cor)
boxplot(as.numeric(levels(test$cor)[test$cor])~test$kit, ylab = "Correlation")
test$cor<-as.numeric(levels(test$cor)[test$cor])
ttests<-list()
ttests[[1]]<- t.test(test$cor[grep("Ill",test$kit)],test$cor[grep("NEXT",test$kit)])
ttests[[2]]<- t.test(test$cor[grep("Ill",test$kit)],test$cor[grep("NEB",test$kit)])
ttests[[3]]<- t.test(test$cor[grep("Ill",test$kit)],test$cor[grep("clontech",test$kit)])
ttests[[4]]<- t.test(test$cor[grep("NEXT",test$kit)],test$cor[grep("NEB",test$kit)])
ttests[[5]]<- t.test(test$cor[grep("NEXT",test$kit)],test$cor[grep("clontech",test$kit)])
ttests[[6]]<- t.test(test$cor[grep("NEB",test$kit)],test$cor[grep("clontech",test$kit)])
Repro_CorrStatsttest<-sapply(ttests, function(x) {
  c(format(x$estimate,digits = 3),
    t = format(x$statistic, digits = 3),
    p.value = format(x$p.value, scientific = FALSE, digits = 2),
    bonferroni = format(.05/6, digits = 2),
    sig = ifelse(x$p.value<(.05/6), "yes", "no"))
})

colnames(Repro_CorrStatsttest)<-sapply(ttests, function(x) {
  c(test = x$data.name)
})

t(Repro_CorrStatsttest)

IllR<-as.matrix(IllR)
Ill1<-abs(IllR[,1]-IllR[,2])
Ill2<-abs(IllR[,1]-IllR[,3])
Ill3<-abs(IllR[,2]-IllR[,3])
Illresid<-data.frame(Ill1, Ill2, Ill3)
boxplot(Illresid, lwd = 2, ylab = 'absolute value of difference b/w counts across triplicates', notch = TRUE, outline = FALSE, at =1:3, xlim=c(0,12), ylim=c(0,1.3), col ='turquoise1', xaxt ='n')
axis(1, at =2, labels = c("Illumina"))

NEBR<-as.matrix(NEBR)
NEB1<-abs(NEBR[,1]-NEBR[,2])
NEB2<-abs(NEBR[,1]-NEBR[,3])
NEB3<-abs(NEBR[,2]-NEBR[,3])
NEBresid<-data.frame(NEB1, NEB2, NEB3)
boxplot(NEBresid, lwd = 2, notch = TRUE, outline = FALSE, add = TRUE, at =4:6, col = 'springgreen1', xaxt ="n")
axis(1, at = 5, labels = c("NEB"))

NEXTR<-as.matrix(NEXTR)
NEXT1<-abs(NEXTR[,1]-NEXTR[,2])
NEXT2<-abs(NEXTR[,1]-NEXTR[,3])
NEXT3<-abs(NEXTR[,2]-NEXTR[,3])
NEXTresid<-data.frame(NEXT1, NEXT2, NEXT3)
boxplot(NEXTresid, lwd = 2, notch = TRUE, outline = FALSE, add =TRUE, at =7:9, col = "maroon1", xaxt ="n")
axis(1, at = 8, labels = c("NEXT"))

clontechR<-as.matrix(clontechR)
clontechResid1<-data.frame
clontech1<-abs(clontechR[,1]-clontechR[,2])
clontech2<-abs(clontechR[,1]-clontechR[,3])
clontech3<-abs(clontechR[,2]-clontechR[,3])
clontechresid<-data.frame(clontech1, clontech2, clontech3)
boxplot(clontechresid, lwd = 2, notch = TRUE, outline = FALSE, add = TRUE, at = 10:12, col = "yellow", xaxt ="n")
axis(1, at = 11, labels = c("Clontech"))

#resid_accur<-data.frame(clontechresid, Illresid, NEBresid, NEXTresid)

###Reproducibility accross experiment##########################################################################################################
clontechA<-data.frame(clontechA)
clontechR<-data.frame(clontechR)
IllR<-data.frame(IllR)
NEBR<-data.frame(NEBR)
NEXTR<-data.frame(NEXTR)
IllA<-data.frame(IllA)
NEBA<-data.frame(NEBA)
NEXTA<-data.frame(NEXTA)

Tests <- list()
diffA<-(which(rownames(clontechA) %in% rownames(clontechR)))
diffB<-(which(rownames(clontechR) %in% rownames(clontechA)))
length(diffA)
Tests[[1]]<-cor.test(clontechA$Clontech_acc_trimmed.1_R1.fq[diffA], clontechR$Clontech_trimmed.10_R1.fq[diffB])
Tests[[2]]<-cor.test(clontechA$Clontech_acc_trimmed.1_R1.fq[diffA], clontechR$Clontech_trimmed.11_R1.fq[diffB])
Tests[[3]]<-cor.test(clontechA$Clontech_acc_trimmed.1_R1.fq[diffA], clontechR$Clontech_trimmed.12_R1.fq[diffB])
Tests[[4]]<-cor.test(clontechA$Clontech_acc_trimmed.2_R1.fq[diffA], clontechR$Clontech_trimmed.10_R1.fq[diffB])
Tests[[5]]<-cor.test(clontechA$Clontech_acc_trimmed.2_R1.fq[diffA], clontechR$Clontech_trimmed.11_R1.fq[diffB])
Tests[[6]]<-cor.test(clontechA$Clontech_acc_trimmed.2_R1.fq[diffA], clontechR$Clontech_trimmed.12_R1.fq[diffB])
Tests[[7]]<-cor.test(clontechA$Clontech_acc_trimmed.3_R1.fq[diffA], clontechR$Clontech_trimmed.10_R1.fq[diffB])
Tests[[8]]<-cor.test(clontechA$Clontech_acc_trimmed.3_R1.fq[diffA], clontechR$Clontech_trimmed.11_R1.fq[diffB])
Tests[[9]]<-cor.test(clontechA$Clontech_acc_trimmed.3_R1.fq[diffA], clontechR$Clontech_trimmed.12_R1.fq[diffB])

diffA<-(which(rownames(IllA) %in% rownames(IllR)))
diffB<-(which(rownames(IllR) %in% rownames(IllA)))
length(diffA)
Tests[[10]]<-cor.test(IllA$Illumina_acc_trimmed.1_R1.fq[diffA], IllR$Illumina_trimmed.1_R1.fq[diffB])
Tests[[11]]<-cor.test(IllA$Illumina_acc_trimmed.1_R1.fq[diffA], IllR$Illumina_trimmed.2_R1.fq[diffB])
Tests[[12]]<-cor.test(IllA$Illumina_acc_trimmed.1_R1.fq[diffA], IllR$Illumina_trimmed.3_R1.fq[diffB])
Tests[[13]]<-cor.test(IllA$Illumina_acc_trimmed.2_R1.fq[diffA], IllR$Illumina_trimmed.1_R1.fq[diffB])
Tests[[14]]<-cor.test(IllA$Illumina_acc_trimmed.2_R1.fq[diffA], IllR$Illumina_trimmed.2_R1.fq[diffB])
Tests[[15]]<-cor.test(IllA$Illumina_acc_trimmed.2_R1.fq[diffA], IllR$Illumina_trimmed.3_R1.fq[diffB])
Tests[[16]]<-cor.test(IllA$Illumina_acc_trimmed.3_R1.fq[diffA], IllR$Illumina_trimmed.1_R1.fq[diffB])
Tests[[17]]<-cor.test(IllA$Illumina_acc_trimmed.3_R1.fq[diffA], IllR$Illumina_trimmed.2_R1.fq[diffB])
Tests[[18]]<-cor.test(IllA$Illumina_acc_trimmed.3_R1.fq[diffA], IllR$Illumina_trimmed.3_R1.fq[diffB])

diffA<-(which(rownames(NEBA) %in% rownames(NEBR)))
diffB<-(which(rownames(NEBR) %in% rownames(NEBA)))
length(diffA)
Tests[[19]]<-cor.test(NEBA$NEB_acc_trimmed.1_R1.fq[diffA], NEBR$NEB_trimmed.10_R1.fq[diffB])
Tests[[20]]<-cor.test(NEBA$NEB_acc_trimmed.1_R1.fq[diffA], NEBR$NEB_trimmed.11_R1.fq[diffB])
Tests[[21]]<-cor.test(NEBA$NEB_acc_trimmed.1_R1.fq[diffA], NEBR$NEB_trimmed.12_R1.fq[diffB])
Tests[[22]]<-cor.test(NEBA$NEB_acc_trimmed.2_R1.fq[diffA], NEBR$NEB_trimmed.10_R1.fq[diffB])
Tests[[23]]<-cor.test(NEBA$NEB_acc_trimmed.2_R1.fq[diffA], NEBR$NEB_trimmed.11_R1.fq[diffB])
Tests[[24]]<-cor.test(NEBA$NEB_acc_trimmed.2_R1.fq[diffA], NEBR$NEB_trimmed.12_R1.fq[diffB])
Tests[[25]]<-cor.test(NEBA$NEB_acc_trimmed.3_R1.fq[diffA], NEBR$NEB_trimmed.10_R1.fq[diffB])
Tests[[26]]<-cor.test(NEBA$NEB_acc_trimmed.3_R1.fq[diffA], NEBR$NEB_trimmed.11_R1.fq[diffB])
Tests[[27]]<-cor.test(NEBA$NEB_acc_trimmed.3_R1.fq[diffA], NEBR$NEB_trimmed.12_R1.fq[diffB])

diffA<-(which(rownames(NEXTA) %in% rownames(NEXTR)))
diffB<-(which(rownames(NEXTR) %in% rownames(NEXTA)))
length(diffA)
Tests[[28]]<-cor.test(NEXTA$NEXT_acc_trimmed.1_R1.fq[diffA], NEXTR$NEXT_trimmed.10_R1.fq[diffB])
Tests[[29]]<-cor.test(NEXTA$NEXT_acc_trimmed.1_R1.fq[diffA], NEXTR$NEXT_trimmed.11_R1.fq[diffB])
Tests[[30]]<-cor.test(NEXTA$NEXT_acc_trimmed.1_R1.fq[diffA], NEXTR$NEXT_trimmed.12_R1.fq[diffB])
Tests[[31]]<-cor.test(NEXTA$NEXT_acc_trimmed.2_R1.fq[diffA], NEXTR$NEXT_trimmed.10_R1.fq[diffB])
Tests[[32]]<-cor.test(NEXTA$NEXT_acc_trimmed.2_R1.fq[diffA], NEXTR$NEXT_trimmed.11_R1.fq[diffB])
Tests[[33]]<-cor.test(NEXTA$NEXT_acc_trimmed.2_R1.fq[diffA], NEXTR$NEXT_trimmed.12_R1.fq[diffB])
Tests[[34]]<-cor.test(NEXTA$NEXT_acc_trimmed.3_R1.fq[diffA], NEXTR$NEXT_trimmed.10_R1.fq[diffB])
Tests[[35]]<-cor.test(NEXTA$NEXT_acc_trimmed.3_R1.fq[diffA], NEXTR$NEXT_trimmed.11_R1.fq[diffB])
Tests[[36]]<-cor.test(NEXTA$NEXT_acc_trimmed.3_R1.fq[diffA], NEXTR$NEXT_trimmed.12_R1.fq[diffB])


Repro_accross_exp_CorrStats<-sapply(Tests, function(x) {
  c(x$estimate[1],
    ci.lower = x$conf.int[1],
    ci.upper = x$conf.int[2],
    p.value = format(x$p.value, scientific = TRUE))
})
colnames(Repro_accross_exp_CorrStats)<-sapply(Tests, function(x) {
  c(test = x$data.name)
})
t(Repro_accross_exp_CorrStats)


boxplot(IllA, lwd = 2, ylab = 'counts across triplicates(log2(RPM +1))', notch = TRUE, outline = FALSE, at =1:3, xlim=c(1,24), ylim=c(0,16) , col ='turquoise1', xaxt='n') #xaxt= 'n' gets rid of x labels
boxplot(IllR, lwd = 2, notch = TRUE, outline = FALSE, at =4:6, add = TRUE, col ='turquoise1', xaxt='n') #xaxt= 'n' gets rid of x labels
axis(1, at =3, labels = c("Illumina"))
boxplot(NEBA, lwd = 2, notch = TRUE, outline = FALSE, add = TRUE, at =7:9, col = 'springgreen1', xaxt='n')
boxplot(NEBR, lwd = 2, notch = TRUE, outline = FALSE, add = TRUE, at =10:12, col = 'springgreen1', xaxt='n')
axis(1, at =9, labels = c("NEB"))
boxplot(NEXTA, lwd = 2, notch = TRUE, outline = FALSE, add =TRUE, at =13:15, col = "maroon1", xaxt='n')
boxplot(NEXTR, lwd = 2, notch = TRUE, outline = FALSE, add =TRUE, at =16:18, col = "maroon1", xaxt='n')
axis(1, at =15, labels = c("NEXT"))
boxplot(clontechA, lwd = 2, xlab = 'Kits', notch = TRUE, outline = FALSE, add = TRUE, at = 19:21, col = "yellow", xaxt='n')
boxplot(clontechR, lwd = 2, xlab = 'Kits', notch = TRUE, outline = FALSE, add = TRUE, at = 22:24, col = "yellow", xaxt='n')
axis(1, at =21, labels = c("Clontech"))

test<-as.data.frame(t(Repro_accross_exp_CorrStats))

test$kit<-(c(rep("clontech",9), rep("Ill", 9), rep("NEB",9), rep("NEXT", 9)))
#test$trip<-rep(c(1:3),4)
#test$kit<-paste0(rep(c(1:3),4), test$kit)
test$kit=factor(test$kit)
Accur_Corr<-data.frame(test$kit, test$cor)
boxplot(as.numeric(levels(test$cor)[test$cor])~test$kit, ylab = "Correlation")
test$cor<-as.numeric(levels(test$cor)[test$cor])

ttests<-list()
ttests[[1]]<-t.test(test$cor[grep("Ill",test$kit)],test$cor[grep("NEXT",test$kit)])
ttests[[2]]<-t.test(test$cor[grep("Ill",test$kit)],test$cor[grep("NEB",test$kit)])
ttests[[3]]<-t.test(test$cor[grep("Ill",test$kit)],test$cor[grep("clontech",test$kit)])
ttests[[4]]<-t.test(test$cor[grep("NEXT",test$kit)],test$cor[grep("NEB",test$kit)])
ttests[[5]]<-t.test(test$cor[grep("NEXT",test$kit)],test$cor[grep("clontech",test$kit)])
ttests[[6]]<-t.test(test$cor[grep("NEB",test$kit)],test$cor[grep("clontech",test$kit)])
Repro_CorrStatsttest<-sapply(ttests, function(x) {
  c(format(x$estimate,digits = 3),
    t = format(x$statistic, digits = 3),
    p.value = format(x$p.value, scientific = FALSE, digits = 2),
    bonferroni = format(.05/6, digits = 2),
    sig = ifelse(x$p.value<(.05/6), "yes", "no"))
})

colnames(Repro_CorrStatsttest)<-sapply(ttests, function(x) {
  c(test = x$data.name)
})
t(Repro_CorrStatsttest)

IllA<-as.matrix(IllA)
IllR<-as.matrix(IllR)
diffA<-(which(rownames(IllA) %in% rownames(IllR)))
diffB<-(which(rownames(IllR) %in% rownames(IllA)))
Ill1<-abs(IllA[,1][diffA]-IllR[,1][diffB])
Ill2<-abs(IllA[,1][diffA]-IllR[,2][diffB])
Ill3<-abs(IllA[,1][diffA]-IllR[,3][diffB])
Ill4<-abs(IllA[,2][diffA]-IllR[,1][diffB])
Ill5<-abs(IllA[,2][diffA]-IllR[,2][diffB])
Ill6<-abs(IllA[,2][diffA]-IllR[,3][diffB])
Ill7<-abs(IllA[,3][diffA]-IllR[,1][diffB])
Ill8<-abs(IllA[,3][diffA]-IllR[,2][diffB])
Ill9<-abs(IllA[,3][diffA]-IllR[,3][diffB])
Illresid<-data.frame(Ill1, Ill2, Ill3, Ill4, Ill5, Ill6, Ill7, Ill8, Ill9)
boxplot(Illresid, lwd = 2, ylab = 'absolute value of difference b/w counts across studies', notch = TRUE, outline = FALSE, at =1:9, xlim=c(0,36), ylim=c(0,2.2), col ='turquoise1', xaxt ='n')
axis(1, at =5, labels = c("Illumina"))

NEBA<-as.matrix(NEBA)
NEBR<-as.matrix(NEBR)
diffA<-(which(rownames(NEBA) %in% rownames(NEBR)))
diffB<-(which(rownames(NEBR) %in% rownames(NEBA)))
NEB1<-abs(NEBA[,1][diffA]-NEBR[,1][diffB])
NEB2<-abs(NEBA[,1][diffA]-NEBR[,2][diffB])
NEB3<-abs(NEBA[,1][diffA]-NEBR[,3][diffB])
NEB4<-abs(NEBA[,2][diffA]-NEBR[,1][diffB])
NEB5<-abs(NEBA[,2][diffA]-NEBR[,2][diffB])
NEB6<-abs(NEBA[,2][diffA]-NEBR[,3][diffB])
NEB7<-abs(NEBA[,3][diffA]-NEBR[,1][diffB])
NEB8<-abs(NEBA[,3][diffA]-NEBR[,2][diffB])
NEB9<-abs(NEBA[,3][diffA]-NEBR[,3][diffB])
NEBresid<-data.frame(NEB1, NEB2, NEB3, NEB4, NEB5, NEB6, NEB7, NEB8, NEB9)
boxplot(NEBresid, lwd = 2, notch = TRUE, outline = FALSE, add = TRUE, at =10:18, col = 'springgreen1', xaxt ="n")
axis(1, at = 14, labels = c("NEB"))

NEXTA<-as.matrix(NEXTA)
NEXTR<-as.matrix(NEXTR)
diffA<-(which(rownames(NEXTA) %in% rownames(NEXTR)))
diffB<-(which(rownames(NEXTR) %in% rownames(NEXTA)))
NEXT1<-abs(NEXTA[,1][diffA]-NEXTR[,1][diffB])
NEXT2<-abs(NEXTA[,1][diffA]-NEXTR[,2][diffB])
NEXT3<-abs(NEXTA[,1][diffA]-NEXTR[,3][diffB])
NEXT4<-abs(NEXTA[,2][diffA]-NEXTR[,1][diffB])
NEXT5<-abs(NEXTA[,2][diffA]-NEXTR[,2][diffB])
NEXT6<-abs(NEXTA[,2][diffA]-NEXTR[,3][diffB])
NEXT7<-abs(NEXTA[,3][diffA]-NEXTR[,1][diffB])
NEXT8<-abs(NEXTA[,3][diffA]-NEXTR[,2][diffB])
NEXT9<-abs(NEXTA[,3][diffA]-NEXTR[,3][diffB])
NEXTresid<-data.frame(NEXT1, NEXT2, NEXT3, NEXT4, NEXT5, NEXT6, NEXT7, NEXT8, NEXT9)
boxplot(NEXTresid, lwd = 2, notch = TRUE, outline = FALSE, add =TRUE, at =19:27, col = "maroon1", xaxt ="n")
axis(1, at = 23, labels = c("NEXT"))

clontechA<-as.matrix(clontechA)
clontechR<-as.matrix(clontechR)
diffA<-(which(rownames(clontechA) %in% rownames(clontechR)))
diffB<-(which(rownames(clontechR) %in% rownames(clontechA)))
clontech1<-abs(clontechA[,1][diffA]-clontechR[,1][diffB])
clontech2<-abs(clontechA[,1][diffA]-clontechR[,2][diffB])
clontech3<-abs(clontechA[,1][diffA]-clontechR[,3][diffB])
clontech4<-abs(clontechA[,2][diffA]-clontechR[,1][diffB])
clontech5<-abs(clontechA[,2][diffA]-clontechR[,2][diffB])
clontech6<-abs(clontechA[,2][diffA]-clontechR[,3][diffB])
clontech7<-abs(clontechA[,3][diffA]-clontechR[,1][diffB])
clontech8<-abs(clontechA[,3][diffA]-clontechR[,2][diffB])
clontech9<-abs(clontechA[,3][diffA]-clontechR[,3][diffB])
clontechresid<-data.frame(clontech1, clontech2, clontech3, clontech4, clontech5, clontech6, clontech7, clontech8, clontech9)
boxplot(clontechresid, lwd = 2, notch = TRUE, outline = FALSE, add = TRUE, at = 28:36, col = "yellow", xaxt ="n")
axis(1, at = 32, labels = c("Clontech"))


####################################venndiagram##################################
library(VennDiagram)
library(venneuler)
###synth miRNA detected

#doesn't include overlaps and greater than 4 expression value
clontech<-as.vector(union(rownames(clontechA), rownames(clontechR)))
Illumina<-as.vector(union(rownames(IllA), rownames(IllR)))
NEB<-as.vector(union(rownames(NEBA), rownames(NEBR)))
NEXT<-as.vector(union(rownames(NEXTA), rownames(NEXTR)))

miRNA_det<-list("clontech" = clontech, "Illumina" = Illumina, "NEB" = NEB, "NextFlex" =NEXT)

synthoverlap<- calculate.overlap(x = miRNA_det)
# 
# test<- venneuler(testdata)
#vennsynthmiRNA<-venn.diagram(Synth_det, filename = NULL)
# grid.draw(vennsynthmiRNA)
MyVenn <- venneuler(c(Clontech=length(clontech),Illumina=length(Illumina),NEB=length(NEB),NextFlex=length(NEXT),
                      "Clontech&Illumina"=length(intersect(clontech, Illumina)),
                      "Clontech&NEB"=length(intersect(clontech, NEB)), 
                      "Clontech&NextFlex"=length(intersect(clontech, NEXT)), 
                      "Clontech&NEB&NextFlex"=length(intersect(intersect(clontech,NEB), NEXT)),
                      "Clontech&Illumina&NEB"=length(intersect(intersect(clontech, Illumina), NEB)), 
                      "Clontech&Illumina&NextFlex"=length(intersect(intersect(clontech, Illumina), NEXT)),
                      "Clontech&Illumina&NextFlex&NEB"=length(intersect(intersect(clontech, Illumina), (intersect(NEB, NEXT)))),
                      "Illumina&NEB"=length(intersect(Illumina, NEB)), 
                      "Illumina&NextFlex"= length(intersect(Illumina, NEXT)),
                      "NEB&NextFlex"=length(intersect(NEB,NEXT)),
                      "NEB&NextFlex&Illumina"=length(intersect(intersect(NEB, NEXT), Illumina))))
plot(MyVenn)

# MyVenn <- venneuler(c(Clontech=21,Illumina=0,NEB=2,NextFlex=2,"Clontech&Illumina"=4,
#                       "Clontech&NEB"=7, "Clontech&NextFlex"=9, "Clontech&NEB&NextFlex"=25, "Clontech&Illumina&NEB"=3, "Clontech&Illumina&NextFlex"=32, 
#                       "Illumina&NEB"=0, "Illumina&NextFlex"=1,"NEB&NextFlex&Illumina"=0,
#                       "NEB&NextFlex"=0))
# plot(MyVenn)
extractDiff <- function(P){
  sampleset = sample(nrow(P), 15, replace=FALSE) #select the first 15 rows, note replace=FALSE
  subA <- P[sampleset, ] # takes the 15 selected rows
  subB <- P[-sampleset, ] # takes the remaining rows in the set
  meanA <- mean(subA$val)
  meanB <- mean(subB$val)
  diff <- abs(meanA-meanB)
  outdf <- c(mA = meanA, mB= meanB, diffAB = diff)
  return(outdf)
}

#####radial plots############################
library(plotrix)
###accuracy

load("/home/carrie/synthetic_miRNA/accuracy_human_Final_stats.rda")
human<-RMSE
load("/home/carrie/synthetic_miRNA/accuracy_RNO_Final_stats.rda")
rat<-RMSE
load("/home/carrie/synthetic_miRNA/accuracy_MMU_Final_stats.rda")
mouse<-RMSE
load("/home/carrie/synthetic_miRNA/accuracy_Virus_Final_stats.rda")
virus<-RMSE


Accuracytestatributes<-data.frame(t(human), t(rat), t(mouse), t(virus))

#diamondplot(Accuracytestatributes)

kit<-c(1:4)
kit.names <-names(kit)<-c("Illumina", "NEB", "NEXT", "Clon")
kit<-human[1, 1:4]
radial.plot(kit, labels = kit.names, rp.type = "p", line.col = "blue", radial.lim =c(0,11))

load("/home/carrie/synthetic_miRNA/accuracy_human_Final_stats.rda")
human<-Var
load("/home/carrie/synthetic_miRNA/accuracy_RNO_Final_stats.rda")
rat<-Var
load("/home/carrie/synthetic_miRNA/accuracy_MMU_Final_stats.rda")
mouse<-Var
load("/home/carrie/synthetic_miRNA/accuracy_Virus_Final_stats.rda")
virus<-Var


Accuracytestatributes<-data.frame(t(human), t(rat), t(mouse), t(virus))

#diamondplot(Accuracytestatributes)
#radial.grid(Accuracytestatributes, grid.pos = 1)

kit<-c(1:4)
kit.names <-names(kit)<-c("Illumina", "NEB", "NEXT", "Clon")
kit<-human[1, 1:4]
radial.plot(kit, labels = kit.names, rp.type = "p", line.col = "red", radial.lim = c(0,11), add = TRUE)

###################################radialplot with each species together
load("/home/carrie/synthetic_miRNA/accuracy_human_Final_stats.rda")
human<-RMSE
load("/home/carrie/synthetic_miRNA/accuracy_RNO_Final_stats.rda")
rat<-RMSE
load("/home/carrie/synthetic_miRNA/accuracy_MMU_Final_stats.rda")
mouse<-RMSE
load("/home/carrie/synthetic_miRNA/accuracy_Virus_Final_stats.rda")
virus<-RMSE

kit<-c(1:4)
kit.names <-names(kit)<-c("Illumina", "NEB", "NEXT", "Clon")
kit<-human[1, 1:4]
radial.plot(kit, labels = kit.names, rp.type = "p", line.col = "blue", radial.lim =c(0,11))
kit<-c(1:4)
kit.names <-names(kit)<-c("Illumina", "NEB", "NEXT", "Clon")
kit<-rat[1, 1:4]
radial.plot(kit, labels = kit.names, rp.type = "p", line.col = "red", radial.lim =c(0,11), add =TRUE)
kit<-c(1:4)
kit.names <-names(kit)<-c("Illumina", "NEB", "NEXT", "Clon")
kit<-mouse[1, 1:4]
radial.plot(kit, labels = kit.names, rp.type = "p", line.col = "green", radial.lim =c(0,11), add =TRUE)
kit<-c(1:4)
kit.names <-names(kit)<-c("Illumina", "NEB", "NEXT", "Clon")
kit<-virus[1, 1:4]
radial.plot(kit, labels = kit.names, rp.type = "p", line.col = "purple", radial.lim =c(0,11), add =TRUE)





##################################heatmap
# ReproCorr<-as.data.frame(Repro_CorrStats)# this gives the confidence intervals
# #colnames(ReproCorr)<-colnames(Repro_mir_RPM)
# Next<-ReproCorr[,40:57]
# Clontec<-ReproCorr[,1:18]
# Illumina<-ReproCorr[,19:27]
# NEB<-ReproCorr[,28:39]
# 
# Nextcor<-t(Next[1,])
# Clontecor<-t(Clontec[1,])
# Illuminacor<-t(Illumina[1,])
# NEBcor<-t(NEB[1,])
# 
# Nextcor<-matrix(as.numeric(Nextcor), nrow = 3, ncol =6)
# Clontecor<-matrix(as.numeric(Clontecor), nrow =3, ncol =6)
# Illuminacor<-matrix(as.numeric(Illuminacor), nrow = 3, ncol =3)
# NEBcor<-matrix(as.numeric(NEBcor), nrow = 3, ncol = 4)
# 
# amounts<-c("100ng","250ng","500ng","1000ng","1500ng","2000ng")
# colnames(Nextcor)<-amounts
# colnames(Clontecor)<-amounts
# colnames(Illuminacor)<-amounts[4:6]
# colnames(NEBcor)<-amounts[1:4]
# 
# 
# Next_heatmap <- heatmap(Nextcor, Rowv=NA, Colv=NA, col = cm.colors(256), scale="column", margins=c(5,10))
# Next_heatmap <- heatmap(Nextcor, Rowv=NA, Colv=NA, col = heat.colors(256), scale="column", margins=c(5,10))
# 
# library(gridExtra)
# library(grid)
# library(ggplot2)
# Next <-melt(Nextcor)
# X <- ggplot(data=Next, aes(x=Var1, y=Var2, fill=value)) + geom_tile()
# X <- X+xlab("Sample Triplicate") +ylab("Starting Amount") + ggtitle("NextFlex")
# X<-X + theme(axis.text.y = element_text(size =10), axis.text.x = element_text(size = 25, colour = "black"), axis.title.y = element_text(size = 20, vjust = 2), title = element_text(size = 25))
# X<-X + scale_fill_continuous(guide = guide_legend(title = "R^2"))
# X<-X+theme_bw()
# #X<-X + coord_fixed(ratio = 10)
# 
# NEB <-melt(NEBcor)
# N<- ggplot(data=NEB, aes(x=Var1, y=Var2, fill=value)) + geom_tile()
# N <- N+xlab("Sample Triplicate") +ylab("Starting Amount") + ggtitle("NEB")
# N <-N + theme(axis.text.y = element_text(size =10), axis.text.x = element_text(size = 25, colour = "black"), axis.title.y = element_text(size = 20, vjust = 2), title = element_text(size = 25))
# N<-N + scale_fill_continuous(guide = guide_legend(title = "R^2"))
# #N<-N + coord_fixed(ratio = 10)
# 
# Illumina <-melt(Illuminacor)
# I <- ggplot(data=Illumina, aes(x=Var1, y=Var2, fill=value)) + geom_tile()
# I <- I+xlab("Sample Triplicate") +ylab("Starting Amount") + ggtitle("Illumina")
# I<-I + theme(axis.text.y = element_text(size =10), axis.text.x = element_text(size = 25, colour = "black"), axis.title.y = element_text(size = 20, vjust = 2), title = element_text(size = 25))
# I<-I + scale_fill_continuous(guide = guide_legend(title = "R^2"))
# #I<-I + coord_fixed(ratio = 10)
# 
# Clontech <-melt(Clontecor)
# C<- ggplot(data=Next, aes(x=Var1, y=Var2, fill=value)) + geom_tile()
# C <- C+xlab("Sample Triplicate") +ylab("Starting Amount") + ggtitle("Clontech")
# C<-C + theme(axis.text.y = element_text(size =10), axis.text.x = element_text(size = 25, colour = "black"), axis.title.y = element_text(size = 20, vjust = 2), title = element_text(size = 25))
# C<-C + scale_fill_continuous(guide = guide_legend(title = "R^2"))
# #C<-C + coord_fixed(ratio = 10)
# 
# grid.arrange(X, N, I, C, ncol =4)
# 
# 
# maxWidth = grid::unit.pmax(gp1$widths[2:3],gp2$widths[2:3], gp3$widths[2:3], gp4$widths[2:3], gpt$widths[2:3])
# gp1$widths[2:3] <- as.list(maxWidth)
# gp2$widths[2:3] <- as.list(maxWidth)
# gp3$widths[2:3] <- as.list(maxWidth)
# gp4$widths[2:3] <- as.list(maxWidth)
# gpt$widths[2:3] <- as.list(maxWidth)
# #With functions grid.arrange() and arrangeGrob() arrange both plots and legend in one plot.
# 
# #grid.arrange(arrangeGrob(arrangeGrob(X,N, I,C,heights=c(3/4,1/4),ncol=4),widths=c(7/8,1/8),ncol=2))
# 
# gp1 <- ggplot_gtable(ggplot_build(X))
# gp2 <- ggplot_gtable(ggplot_build(N))
# gp3 <- ggplot_gtable(ggplot_build(I))
# gp4 <- ggplot_gtable(ggplot_build(C))
# 
# # maxWidth = grid::unit.pmax(gp1$widths[2:3],gp2$widths[2:3], gp3$widths[2:3], gp4$widths[2:3], gpt$widths[2:3])
# # gp1$widths[2:3] <- as.list(maxWidth)
# # gp2$widths[2:3] <- as.list(maxWidth)
# # gp3$widths[2:3] <- as.list(maxWidth)
# # gp4$widths[2:3] <- as.list(maxWidth)
# # gpt$widths[2:3] <- as.list(maxWidth)
# #With functions grid.arrange() and arrangeGrob() arrange both plots and legend in one plot.
# #grid.arrange(arrangeGrob(arrangeGrob(gp1,gp2,gp3,gp4,gpt,heights=c(3/4,1/4),ncol=4),widths=c(7/8,1/8),nrow =2))
# grid.arrange(gp1,gp2,gp3,gp4, layout_matrix = rbind(c(1,2), c(3,4)))


