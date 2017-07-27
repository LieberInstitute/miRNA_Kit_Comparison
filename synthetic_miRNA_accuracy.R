#The synthetic miRNA were mapped using bowtie...directly to the sequence that was known. The synth sequences were equimolar.
#this is only for human synth sequences
setwd("~/synthetic_miRNA/counts/")#set directory where all of the files are
#Clontech<-read.table("~/synthetic_miRNA/Clontech.txt", header =  T) #name, length, #mapped reads, #unmapped reads from samtools idxstat

all.the.files <-list.files()#make a list of all the files in the directory
all.the.data <- lapply(all.the.files, read.table, header = F)#read the data from each file into a list of list
names(all.the.data)<- c(all.the.files)#name each of the lists the name of the file
#miRNA<-lapply(all.the.data, function(x) x[1])
synthCounts<-data.frame(all.the.data)
synthCounts<-synthCounts[1:(length(rownames(synthCounts))-1),]#last row needs to be removed

#### How many of the synthetic sequencess are not detected for each kit?
which(synthCounts$Clontech.txt.V3 == 0) # 555 out of 556 are detected
which(synthCounts$Illumina.txt.V3 == 0) # all 556 are detected
which(synthCounts$NEB.txt.V3 == 0)# all 556 are detected
which(synthCounts$NEXT.txt.V3 == 0)# all 556 are detected

#### How many of the synthetic sequences are not mapped for each kit?
which(synthCounts$Clontech.txt.V4 != 0) # no unmapped reads
which(synthCounts$Illumina.txt.V4 != 0) # no unmapped reads
which(synthCounts$NEB.txt.V4 != 0)# no unmapped reads
which(synthCounts$NEXT.txt.V4 != 0)# no unmapped reads

#cleaning the data
index<-c(1, seq(from =3, by =4, to= length(colnames(synthCounts))))#index to just get the columns we want--- the count columns...as opposed to 1)geneID, 2)sequence length, 4)unmapped reads
trimmedCounts<- synthCounts[,index]

#mapped for all sequences
Mapped<-data.frame(Clontech = 69275085, Illumina = 39625101, NEB = 19793767, NEXT = 31817768)# number of mapped reads from bowtie output
#this is Total mapped to all sequences... or should I do mapped to human etc?.... I guess all instead of species specific makes sense

#justmapped for human sequences
totalmapped<- colSums(trimmedCounts[2:5])

#RPM Normalization
Clontech <- trimmedCounts$Clontech.txt.V3/(Mapped$Clontech/1e6)
Illumina<- trimmedCounts$Illumina.txt.V3/(Mapped$Illumina/1e6)
NEB <- trimmedCounts$NEB.txt.V3/(Mapped$NEB/1e6)
NEXT <- trimmedCounts$NEXT.txt.V3/(Mapped$NEXT/1e6)
RPM<-data.frame(Clontech, Illumina, NEB, NEXT)
rownames(RPM)<- trimmedCounts$Clontech.txt.V1#give rowname sequence identity

save(synthCounts, trimmedCounts, Mapped, totalmapped, Clontech, Illumina, NEB, NEXT, RPM, file ="/home/carrie/synthetic_miRNA/synthmiRNAdata.rda")

RPMfilt<-RPM[which(rowMeans(RPM)>1),] ##if we want to threshold for expression...not sure if we should... and if so how much... often RPM >10 for biological studies...but this is a technical study

#log2 normalization
Counts<-log2(RPMfilt +1)
#Counts<-RPM # if we dont want to do log2 normalization


#makes a plot of the data
boxplot(Counts, lwd = 2, ylab = 'log2(Counts(RPM)+1)', notch = TRUE, outline = FALSE)
stripchart(Counts, vertical = TRUE, method = "jitter", add = TRUE, pch = 20, col = 'blue')


###Now calculate the average for each kit... todetermine the expected value... if there are miRNAs not detected this will pull the average down...so I did not include the one zero value for Clontech

#expectedCount<-(sum(Counts$Illumina.txt.V3, Counts$NEB.txt.V3, Counts$NNEXT.txt.V3))/(3*length(Counts$NEB.txt.V3)) # I did not include clontec here... but maybe it would be the average for each kit??/
expCountIll<-(sum(Counts$Illumina))/ (length(Counts$Illumina))
expCountNEB<-(sum(Counts$NEB))/ (length(Counts$NEB))
expCountNEXT<-(sum(Counts$NEX))/ (length(Counts$NEXT))
expCountClon<-(sum(Counts$Clontech[which(Counts$Clontech>0)]))/ (length(Counts$Clontech -1))#missing one sequence that was not detected - should I account for it?
expCountClon<-(sum(Counts$Clontech))/ (length(Counts$Clontech))#missing one sequence that was not detected - shouls I not account for it?

###Now calculate the difference from the mean

diff_Ill<- abs(Counts$Illumina - expCountIll)
diff_NEB<- abs(Counts$NEB - expCountNEB)
diff_NEXT <- abs(Counts$NEXT - expCountNEXT)
diff_Clon <- abs(Counts$Clontech[which(Counts$Clontech>0)]) - expCountClon # maybe should only do those that are detected?
diff_Clon <- abs(Counts$Clontech - expCountClon)

diff<-data.frame(diff_Ill, diff_NEB, diff_NEXT, diff_Clon)


###ttest of difference from the mean for between each kit
ttests <- list()
ttests[[1]]<-t.test(diff_Ill, diff_NEXT)
ttests[[2]]<-t.test(diff_Ill, diff_NEB)
ttests[[3]]<-t.test(diff_Ill, diff_Clon)
ttests[[4]]<-t.test(diff_NEXT, diff_NEB)
ttests[[5]]<-t.test(diff_NEXT, diff_Clon)
ttests[[6]]<-t.test(diff_NEB, diff_Clon)

# extract your values using `sapply`
Acc_ttestStats<-sapply(ttests, function(x) {
  c(t =format(x$statistic, digits = 2),
    df = format(x$parameter, digits = 0),
    p.value = format(x$p.value, scientific = FALSE, digits = 2),
    bonferroni = format(.05/6, digits = 2),
    sig = ifelse(x$p.value<(.05/6), "yes", "no"))
})
colnames(Acc_ttestStats)<-sapply(ttests, function(x) {
  c(test = x$data.name)
})

#plot of diff
boxplot(diff, lwd = 2, ylab = 'Absolute value of difference from expected mean', notch = TRUE, outline = FALSE, ylim =c(0,9), col = c('turquoise1', 'springgreen1', 'maroon1', 'yellow'))
#abline(a=0, b=0)
#stripchart(diff, vertical = TRUE, method = "jitter", add = TRUE, pch = 20, col = c('lightblue', 'green', 'red', 'khaki1'))
text(x=1, y=3, "RMSE=3.05")
text(x=2, y=2.8, "RMSE=2.89")
text(x=3, y=2.2, "RMSE=2.35")
text(x=4, y=2, "RMSE=2.14")
##Root mean squared error for each sequence
RMSE_Ill <- sqrt(mean((diff_Ill)^2))
RMSE_NEB <- sqrt(mean((diff_NEB)^2))
RMSE_NEXT <- sqrt(mean((diff_NEXT)^2))
RMSE_Clon <- sqrt(mean((diff_Clon)^2))

#RMSE
RMSE_Ill
RMSE_NEB
RMSE_NEXT
RMSE_Clon
RMSE<-data.frame(RMSE_Ill, RMSE_NEB, RMSE_NEXT, RMSE_Clon)
test<-(melt(RMSE))
#barplot(RMSE)
barplot(test$value, names.arg = test$variable, xlab = test$variable)

#Ranges
max(Counts$Clontech)
max(Counts$Illumina)
max(Counts$NEB)
max(Counts$NEXT)
min(Counts$Clontech[which(Counts$Clontech>0)])
min(Counts$Clontech)
min(Counts$Illumina)
min(Counts$NEB)
min(Counts$NEXT)


#sumsquareddifference
SSD_Ill <- sum(diff_Ill^2)
SSD_NEB <- sum(diff_NEB^2)
SSD_NEXT <- sum(diff_NEXT^2)
SSD_Clon <- sum(diff_Clon^2)

SSD_Ill
SSD_NEB
SSD_NEXT
SSD_Clon
SSD<-(data.frame(SSD_Ill, SSD_NEB, SSD_NEXT, SSD_Clon))
test<-(melt(SSD))
barplot(test$value, names.arg = test$Var1)

IllVar<-var(Counts$Illumina)
NEBvar<-var(Counts$NEB)
NEXTvar<-var(Counts$NEXT)
Clonvar<-var(Counts$Clontech)
Var<-data.frame(IllVar, NEBvar, NEXTvar, Clonvar)
save(RMSE, diff, SSD, Var, Acc_ttestStats, file = "/home/carrie/synthetic_miRNA/accuracy_human_Final_stats.rda")

#################################################################Rat

setwd("~/synthetic_miRNA/nonHumanAligned/")#set directory where all of the files are
#Clontech<-read.table("~/synthetic_miRNA/Clontech.txt", header =  T) #name, length, #mapped reads, #unmapped reads from samtools idxstat

all.the.files <-list.files()#make a list of all the files in the directory
all.the.data <- lapply(all.the.files, read.table, header = F)#read the data from each file into a list of list
names(all.the.data)<- c(all.the.files)#name each of the lists the name of the file

RNOcounts <- data.frame(all.the.data[grep("RNO",names(all.the.data))])
RNOcounts<- RNOcounts[1:(length(rownames(RNOcounts))-1),]#last row needs to be removed

which(RNOcounts$Clontech_acc_trimmed.4_R1_synRNO.idxstats.V3 == 0) # 12 out of 390 are not detected(378)####ok so these are duplicates....only 379 unqiue sequences... and last one is nothing so...378
which(RNOcounts$Illumina_acc_trimmed.4_R1_synRNO.idxstats.V3 == 0) # 12 out of 390 are not detected
which(RNOcounts$NEB_acc_trimmed.4_R1_synRNO.idxstats.V3 == 0)# 12 out of 390 are not detected
which(RNOcounts$NEXT_acc_trimmed.4_R1_synRNO.idxstats.V3 == 0)# 12 out of 390 are not detected

which(RNOcounts$Clontech_acc_trimmed.4_R1_synRNO.idxstats.V4 != 0) # no unmapped reads
which(RNOcounts$Illumina_acc_trimmed.4_R1_synRNO.idxstats.V4 != 0) # no unmapped reads
which(RNOcounts$NEB_acc_trimmed.4_R1_synRNO.idxstats.V4 != 0)# no unmapped reads
which(RNOcounts$NEXT_acc_trimmed.4_R1_synRNO.idxstats.V4!= 0)# no unmapped reads


Mapped<-data.frame(Clontech = 69275085, Illumina = 39625101, NEB = 19793767, NEXT = 31817768)# number of mapped reads from bowtie output

index<-c(1, seq(from =3, by =4, to= length(colnames(RNOcounts))))#index to just get the columns we want--- the count columns...as opposed to 1)geneID, 2)sequence length, 4)unmapped reads
trimmedCounts<-RNOcounts[,index]

#RPM Normalization
Clontech <- trimmedCounts$Clontech_acc_trimmed.4_R1_synRNO.idxstats.V3/(Mapped$Clontech/1e6)
Illumina<- trimmedCounts$Illumina_acc_trimmed.4_R1_synRNO.idxstats.V3/(Mapped$Illumina/1e6)
NEB <- trimmedCounts$NEB_acc_trimmed.4_R1_synRNO.idxstats.V3/(Mapped$NEB/1e6)
NEXT <- trimmedCounts$NEXT_acc_trimmed.4_R1_synRNO.idxstats.V3/(Mapped$NEXT/1e6)

RPM<-data.frame(Clontech, Illumina, NEB, NEXT)
RPM<-RPM[which(rowMeans(RPM)>0),]#to get rid of the empty rows for the duplicates
rownames(RPM)<- unique(trimmedCounts$Clontech_acc_trimmed.4_R1_synRNO.idxstats.V1)#give rowname sequence identity

#RPM<-RPM[which(rowMeans(RPM)>10),]#little bit of noise... down from 378 to 372
#log2 normalization

Counts<-log2(RPM +1)
#Counts<-RPM

boxplot(Counts, lwd = 2, ylab = 'log2(Counts(RPM)+1)', notch = TRUE, outline = FALSE)
stripchart(Counts, vertical = TRUE, method = "jitter", add = TRUE, pch = 20, col = 'blue')


###Now calculate the average for each kit... todetermine the expected value... if there are miRNAs not detected this will pull the average down...so I did not include the one zero value for Clontech

#expectedCount<-(sum(Counts$Illumina.txt.V3, Counts$NEB.txt.V3, Counts$NNEXT.txt.V3))/(3*length(Counts$NEB.txt.V3)) # I did not include clontec here... but maybe it would be the average for each kit??/
expCountIll<-(sum(Counts$Illumina))/ (length(Counts$Illumina))
expCountNEB<-(sum(Counts$NEB))/ (length(Counts$Illumina))
expCountNEXT<-(sum(Counts$NEX))/ (length(Counts$Illumina))
expCountClon<-(sum(Counts$Clontech[which(Counts$Clontech>0)]))/ (length(Counts$Illumina)-1)#missing one sequence that was not detected
expCountClon<-(sum(Counts$Clontech[which(Counts$Clontech>0)]))/ (length(Counts$Illumina))#missing one sequence that was not detected

###Now calculate the difference from the mean

diff_Ill<- abs(Counts$Illumina - expCountIll)
diff_NEB<- abs(Counts$NEB - expCountNEB)
diff_NEXT <- abs(Counts$NEXT - expCountNEXT)
diff_Clon <- abs(Counts$Clontech) - expCountClon # maybe should only do those that are detected?
diff_Clon <- abs(Counts$Clontech - expCountClon)
diff<-data.frame(diff_Ill, diff_NEB, diff_NEXT, diff_Clon)
###ttest of difference from the mean for between each kit
ttests <- list()
ttests[[1]]<-t.test(diff_Ill, diff_NEXT)
ttests[[2]]<-t.test(diff_Ill, diff_NEB)
ttests[[3]]<-t.test(diff_Ill, diff_Clon)
ttests[[4]]<-t.test(diff_NEXT, diff_NEB)
ttests[[5]]<-t.test(diff_NEXT, diff_Clon)
ttests[[6]]<-t.test(diff_NEB, diff_Clon)

# extract your values using `sapply`
Acc_ttestStats<-sapply(ttests, function(x) {
  c(t =format(x$statistic, digits = 2),
    df = format(x$parameter, digits = 0),
    p.value = format(x$p.value, scientific = FALSE, digits = 2),
    bonferroni = format(.05/6, digits = 2),
    sig = ifelse(x$p.value<(.05/6), "yes", "no"))
})
colnames(Acc_ttestStats)<-sapply(ttests, function(x) {
  c(test = x$data.name)
})


boxplot(diff, lwd = 2, ylab = 'A bsolute value of difference from expected mean', notch = TRUE, outline = FALSE, col = c('turquoise1', 'springgreen1', 'maroon1', 'yellow'))
#abline(a=0, b=0)
stripchart(diff, vertical = TRUE, method = "jitter", add = TRUE, pch = 20, col = c('lightblue', 'green', 'red', 'khaki1'))
text(x=1, y=0.4, "RMSE=3.11")
text(x=2, y=0.5, "RMSE=2.96")
text(x=3, y=1, "RMSE=2.49")
text(x=4, y=0.9, "RMSE=2.22")
##Root mean squared error for each sequence
RMSE_Ill <- sqrt(mean((diff_Ill)^2))
RMSE_NEB <- sqrt(mean((diff_NEB)^2))
RMSE_NEXT <- sqrt(mean((diff_NEXT)^2))
RMSE_Clon <- sqrt(mean((diff_Clon)^2))

#RMSE
RMSE_Ill
RMSE_NEB
RMSE_NEXT
RMSE_Clon
RMSE<-data.frame(RMSE_Ill, RMSE_NEB, RMSE_NEXT, RMSE_Clon)
test<-(melt(RMSE))
barplot(RMSE)
barplot(test$value, names.arg = test$variable, xlab = test$variable)

#Ranges
max(Counts$Clontech[which(Counts$Clontech>0)])
max(Counts$Illumina)
max(Counts$NEB)
max(Counts$NEXT)
min(Counts$Clontech[which(Counts$Clontech>0)])
min(Counts$Illumina)
min(Counts$NEB)
min(Counts$NEXT)


#sumsquareddifference
SSD_Ill <- sum(diff_Ill^2)
SSD_NEB <- sum(diff_NEB^2)
SSD_NEXT <- sum(diff_NEXT^2)
SSD_Clon <- sum(diff_Clon^2)

SSD_Ill
SSD_NEB
SSD_NEXT
SSD_Clon
SSD<-(data.frame(SSD_Ill, SSD_NEB, SSD_NEXT, SSD_Clon))
test<-(melt(SSD))
barplot(test$value, names.arg = test$Var1)


IllVar<-var(Counts$Illumina)
NEBvar<-var(Counts$NEB)
NEXTvar<-var(Counts$NEXT)
Clonvar<-var(Counts$Clontech)
Var<-data.frame(IllVar, NEBvar, NEXTvar, Clonvar)
save(RMSE, diff, SSD, Var, Acc_ttestStats, file = "/home/carrie/synthetic_miRNA/accuracy_RNO_Final_stats.rda")

#################################################################################################Mouse########################################################

setwd("~/synthetic_miRNA/nonHumanAligned/")#set directory where all of the files are
#Clontech<-read.table("~/synthetic_miRNA/Clontech.txt", header =  T) #name, length, #mapped reads, #unmapped reads from samtools idxstat

all.the.files <-list.files()#make a list of all the files in the directory
all.the.data <- lapply(all.the.files, read.table, header = F)#read the data from each file into a list of list
names(all.the.data)<- c(all.the.files)#name each of the lists the name of the file

RNOcounts <- data.frame(all.the.data[grep("MMU",names(all.the.data))])
RNOcounts<- RNOcounts[1:(length(rownames(RNOcounts))-1),]#last row needs to be removed

which(RNOcounts$Clontech_acc_trimmed.4_R1_synMMU.idxstats.V3 == 0) # 16 out of 499are not detected(378)####ok so these are duplicates....so 486 for the 13
which(RNOcounts$Illumina_acc_trimmed.4_R1_synMMU.idxstats.V3 == 0) # 13 out of 499 are not detected
which(RNOcounts$NEB_acc_trimmed.4_R1_synMMU.idxstats.V3 == 0)# 13 out of 499 are not detected
which(RNOcounts$NEXT_acc_trimmed.4_R1_synMMU.idxstats.V3 == 0)# 13 out of 499 are not detected

which(RNOcounts$Clontech_acc_trimmed.4_R1_synMMU.idxstats.V4 != 0) # no unmapped reads
which(RNOcounts$Illumina_acc_trimmed.4_R1_synMMU.idxstats.V4!= 0) # no unmapped reads
which(RNOcounts$NEB_acc_trimmed.4_R1_synMMU.idxstats.V4 != 0)# no unmapped reads
which(RNOcounts$NEXT_acc_trimmed.4_R1_synMMU.idxstats.V4!= 0)# no unmapped reads


Mapped<-data.frame(Clontech = 69275085, Illumina = 39625101, NEB = 19793767, NEXT = 31817768)# number of mapped reads from bowtie output

index<-c(1, seq(from =3, by =4, to= length(colnames(RNOcounts))))#index to just get the columns we want--- the count columns...as opposed to 1)geneID, 2)sequence length, 4)unmapped reads
trimmedCounts<-RNOcounts[,index]

#RPM Normalization
Clontech <- trimmedCounts$Clontech_acc_trimmed.4_R1_synMMU.idxstats.V3/(Mapped$Clontech/1e6)
Illumina<- trimmedCounts$Illumina_acc_trimmed.4_R1_synMMU.idxstats.V3/(Mapped$Illumina/1e6)
NEB <- trimmedCounts$NEB_acc_trimmed.4_R1_synMMU.idxstats.V3/(Mapped$NEB/1e6)
NEXT <- trimmedCounts$NEXT_acc_trimmed.4_R1_synMMU.idxstats.V3/(Mapped$NEXT/1e6)

RPM<-data.frame(Clontech, Illumina, NEB, NEXT)
RPM<-RPM[which(rowMeans(RPM)>0),]#to get rid of the empty rows for the duplicates
rownames(RPM)<- unique(trimmedCounts$Clontech_acc_trimmed.4_R1_synRNO.idxstats.V1)#give rowname sequence identity

#RPM<-RPM[which(rowMeans(RPM)>10),]#little bit of noise... down from 378 to 372
#log2 normalization

Counts<-log2(RPM +1)
#Counts<-RPM

boxplot(Counts, lwd = 2, ylab = 'log2(Counts(RPM)+1)', notch = TRUE, outline = FALSE)
stripchart(Counts, vertical = TRUE, method = "jitter", add = TRUE, pch = 20, col = 'blue')


###Now calculate the average for each kit... todetermine the expected value... if there are miRNAs not detected this will pull the average down...so I did not include the one zero value for Clontech

#expectedCount<-(sum(Counts$Illumina.txt.V3, Counts$NEB.txt.V3, Counts$NNEXT.txt.V3))/(3*length(Counts$NEB.txt.V3)) # I did not include clontec here... but maybe it would be the average for each kit??/
expCountIll<-(sum(Counts$Illumina))/ (length(Counts$Illumina))
expCountNEB<-(sum(Counts$NEB))/ (length(Counts$Illumina))
expCountNEXT<-(sum(Counts$NEX))/ (length(Counts$Illumina))
expCountClon<-(sum(Counts$Clontech[which(Counts$Clontech>0)]))/ (length(Counts$Illumina)-1)#missing one sequence that was not detected
expCountClon<-(sum(Counts$Clontech[which(Counts$Clontech>0)]))/ (length(Counts$Illumina))#missing one sequence that was not detected

###Now calculate the difference from the mean

diff_Ill<- abs(Counts$Illumina - expCountIll)
diff_NEB<- abs(Counts$NEB - expCountNEB)
diff_NEXT <- abs(Counts$NEXT - expCountNEXT)
diff_Clon <- abs(Counts$Clontech) - expCountClon # maybe should only do those that are detected?
diff_Clon <- abs(Counts$Clontech - expCountClon)
diff<-data.frame(diff_Ill, diff_NEB, diff_NEXT, diff_Clon)

###ttest of difference from the mean for between each kit
ttests <- list()
ttests[[1]]<-t.test(diff_Ill, diff_NEXT)
ttests[[2]]<-t.test(diff_Ill, diff_NEB)
ttests[[3]]<-t.test(diff_Ill, diff_Clon)
ttests[[4]]<-t.test(diff_NEXT, diff_NEB)
ttests[[5]]<-t.test(diff_NEXT, diff_Clon)
ttests[[6]]<-t.test(diff_NEB, diff_Clon)

# extract your values using `sapply`
Acc_ttestStats<-sapply(ttests, function(x) {
  c(t =format(x$statistic, digits = 2),
    df = format(x$parameter, digits = 0),
    p.value = format(x$p.value, scientific = FALSE, digits = 2),
    bonferroni = format(.05/6, digits = 2),
    sig = ifelse(x$p.value<(.05/6), "yes", "no"))
})
colnames(Acc_ttestStats)<-sapply(ttests, function(x) {
  c(test = x$data.name)
})


boxplot(diff, lwd = 2, ylab = 'A bsolute value of difference from expected mean', notch = TRUE, outline = FALSE, col = c('turquoise1', 'springgreen1', 'maroon1', 'yellow'))
#abline(a=0, b=0)
stripchart(diff, vertical = TRUE, method = "jitter", add = TRUE, pch = 20, col = c('lightblue', 'green', 'red', 'khaki1'))
text(x=1, y=0.4, "RMSE=3.11")
text(x=2, y=0.5, "RMSE=2.96")
text(x=3, y=1, "RMSE=2.49")
text(x=4, y=0.9, "RMSE=2.22")
##Root mean squared error for each sequence
RMSE_Ill <- sqrt(mean((diff_Ill)^2))
RMSE_NEB <- sqrt(mean((diff_NEB)^2))
RMSE_NEXT <- sqrt(mean((diff_NEXT)^2))
RMSE_Clon <- sqrt(mean((diff_Clon)^2))

#RMSE
RMSE_Ill
RMSE_NEB
RMSE_NEXT
RMSE_Clon
RMSE<-data.frame(RMSE_Ill, RMSE_NEB, RMSE_NEXT, RMSE_Clon)
test<-(melt(RMSE))
barplot(RMSE)
barplot(test$value, names.arg = test$variable, xlab = test$variable)

#Ranges
max(Counts$Clontech[which(Counts$Clontech>0)])
max(Counts$Illumina)
max(Counts$NEB)
max(Counts$NEXT)
min(Counts$Clontech[which(Counts$Clontech>0)])
min(Counts$Illumina)
min(Counts$NEB)
min(Counts$NEXT)


#sumsquareddifference
SSD_Ill <- sum(diff_Ill^2)
SSD_NEB <- sum(diff_NEB^2)
SSD_NEXT <- sum(diff_NEXT^2)
SSD_Clon <- sum(diff_Clon^2)

SSD_Ill
SSD_NEB
SSD_NEXT
SSD_Clon
SSD<-(data.frame(SSD_Ill, SSD_NEB, SSD_NEXT, SSD_Clon))
test<-(melt(SSD))
barplot(test$value, names.arg = test$Var1)


IllVar<-var(Counts$Illumina)
NEBvar<-var(Counts$NEB)
NEXTvar<-var(Counts$NEXT)
Clonvar<-var(Counts$Clontech)
Var<-data.frame(IllVar, NEBvar, NEXTvar, Clonvar)
save(RMSE, diff, SSD, Var, Acc_ttestStats, file = "/home/carrie/synthetic_miRNA/accuracy_MMU_Final_stats.rda")

#################################################################################################Virus########################################################

setwd("~/synthetic_miRNA/nonHumanAligned/")#set directory where all of the files are
#Clontech<-read.table("~/synthetic_miRNA/Clontech.txt", header =  T) #name, length, #mapped reads, #unmapped reads from samtools idxstat

all.the.files <-list.files()#make a list of all the files in the directory
all.the.data <- lapply(all.the.files, read.table, header = F)#read the data from each file into a list of list
names(all.the.data)<- c(all.the.files)#name each of the lists the name of the file

RNOcounts <- data.frame(all.the.data[grep("Virus",names(all.the.data))])
RNOcounts<- RNOcounts[1:(length(rownames(RNOcounts))-1),]#last row needs to be removed

which(RNOcounts$Clontech_acc_trimmed.4_R1_synAllVirus.idxstats.V3 == 0) # all103 detected
which(RNOcounts$Illumina_acc_trimmed.4_R1_synAllVirus.idxstats.V3 == 0) # all103 detected
which(RNOcounts$NEB_acc_trimmed.4_R1_synAllVirus.idxstats.V3 == 0)# 1 out of 103 are not detected
which(RNOcounts$NEXT_acc_trimmed.4_R1_synAllVirus.idxstats.V3 == 0)# all103 detected

which(RNOcounts$Clontech_acc_trimmed.4_R1_synAllVirus.idxstats.V4 != 0) # no unmapped reads
which(RNOcounts$Illumina_acc_trimmed.4_R1_synAllVirus.idxstats.V4 != 0) # no unmapped reads
which(RNOcounts$NEB_acc_trimmed.4_R1_synAllVirus.idxstats.V4 != 0)# no unmapped reads
which(RNOcounts$NEXT_acc_trimmed.4_R1_synAllVirus.idxstats.V4 != 0)# no unmapped reads


Mapped<-data.frame(Clontech = 69275085, Illumina = 39625101, NEB = 19793767, NEXT = 31817768)# number of mapped reads from bowtie output

index<-c(1, seq(from =3, by =4, to= length(colnames(RNOcounts))))#index to just get the columns we want--- the count columns...as opposed to 1)geneID, 2)sequence length, 4)unmapped reads
trimmedCounts<-RNOcounts[,index]

#RPM Normalization
Clontech <- trimmedCounts$Clontech_acc_trimmed.4_R1_synAllVirus.idxstats.V3/(Mapped$Clontech/1e6)
Illumina<- trimmedCounts$Illumina_acc_trimmed.4_R1_synAllVirus.idxstats.V3/(Mapped$Illumina/1e6)
NEB <- trimmedCounts$NEB_acc_trimmed.4_R1_synAllVirus.idxstats.V3/(Mapped$NEB/1e6)
NEXT <- trimmedCounts$NEXT_acc_trimmed.4_R1_synAllVirus.idxstats.V3/(Mapped$NEXT/1e6)

RPM<-data.frame(Clontech, Illumina, NEB, NEXT)
RPM<-RPM[which(rowMeans(RPM)>0),]#to get rid of the empty rows for the duplicates
rownames(RPM)<- unique(trimmedCounts$Clontech_acc_trimmed.4_R1_synRNO.idxstats.V1)#give rowname sequence identity

#RPM<-RPM[which(rowMeans(RPM)>10),]#little bit of noise... down from 378 to 372
#log2 normalization

Counts<-log2(RPM +1)
#Counts<-RPM

boxplot(Counts, lwd = 2, ylab = 'log2(Counts(RPM)+1)', notch = TRUE, outline = FALSE)
stripchart(Counts, vertical = TRUE, method = "jitter", add = TRUE, pch = 20, col = 'blue')


###Now calculate the average for each kit... todetermine the expected value... if there are miRNAs not detected this will pull the average down...so I did not include the one zero value for Clontech

#expectedCount<-(sum(Counts$Illumina.txt.V3, Counts$NEB.txt.V3, Counts$NNEXT.txt.V3))/(3*length(Counts$NEB.txt.V3)) # I did not include clontec here... but maybe it would be the average for each kit??/
expCountIll<-(sum(Counts$Illumina))/ (length(Counts$Illumina))
expCountNEB<-(sum(Counts$NEB))/ (length(Counts$Illumina))
expCountNEXT<-(sum(Counts$NEX))/ (length(Counts$Illumina))
expCountClon<-(sum(Counts$Clontech[which(Counts$Clontech>0)]))/ (length(Counts$Illumina)-1)#missing one sequence that was not detected
expCountClon<-(sum(Counts$Clontech[which(Counts$Clontech>0)]))/ (length(Counts$Illumina))#missing one sequence that was not detected

###Now calculate the difference from the mean

diff_Ill<- abs(Counts$Illumina - expCountIll)
diff_NEB<- abs(Counts$NEB - expCountNEB)
diff_NEXT <- abs(Counts$NEXT - expCountNEXT)
diff_Clon <- abs(Counts$Clontech) - expCountClon # maybe should only do those that are detected?
diff_Clon <- abs(Counts$Clontech - expCountClon)
diff<-data.frame(diff_Ill, diff_NEB, diff_NEXT, diff_Clon)

###ttest of difference from the mean for between each kit
ttests <- list()
ttests[[1]]<-t.test(diff_Ill, diff_NEXT)
ttests[[2]]<-t.test(diff_Ill, diff_NEB)
ttests[[3]]<-t.test(diff_Ill, diff_Clon)
ttests[[4]]<-t.test(diff_NEXT, diff_NEB)
ttests[[5]]<-t.test(diff_NEXT, diff_Clon)
ttests[[6]]<-t.test(diff_NEB, diff_Clon)

# extract your values using `sapply`
Acc_ttestStats<-sapply(ttests, function(x) {
  c(t =format(x$statistic, digits = 2),
    df = format(x$parameter, digits = 0),
    p.value = format(x$p.value, scientific = FALSE, digits = 2),
    bonferroni = format(.05/6, digits = 2),
    sig = ifelse(x$p.value<(.05/6), "yes", "no"))
})
colnames(Acc_ttestStats)<-sapply(ttests, function(x) {
  c(test = x$data.name)
})


boxplot(diff, lwd = 2, ylab = 'A bsolute value of difference from expected mean', notch = TRUE, outline = FALSE, col = c('turquoise1', 'springgreen1', 'maroon1', 'yellow'))
#abline(a=0, b=0)
stripchart(diff, vertical = TRUE, method = "jitter", add = TRUE, pch = 20, col = c('lightblue', 'green', 'red', 'khaki1'))
text(x=1, y=0.4, "RMSE=3.11")
text(x=2, y=0.5, "RMSE=2.96")
text(x=3, y=1, "RMSE=2.49")
text(x=4, y=0.9, "RMSE=2.22")
##Root mean squared error for each sequence
RMSE_Ill <- sqrt(mean((diff_Ill)^2))
RMSE_NEB <- sqrt(mean((diff_NEB)^2))
RMSE_NEXT <- sqrt(mean((diff_NEXT)^2))
RMSE_Clon <- sqrt(mean((diff_Clon)^2))

#RMSE
RMSE_Ill
RMSE_NEB
RMSE_NEXT
RMSE_Clon
RMSE<-data.frame(RMSE_Ill, RMSE_NEB, RMSE_NEXT, RMSE_Clon)
test<-(melt(RMSE))
barplot(RMSE)
barplot(test$value, names.arg = test$variable, xlab = test$variable)

#Ranges
max(Counts$Clontech[which(Counts$Clontech>0)])
max(Counts$Illumina)
max(Counts$NEB)
max(Counts$NEXT)
min(Counts$Clontech[which(Counts$Clontech>0)])
min(Counts$Illumina)
min(Counts$NEB)
min(Counts$NEXT)


#sumsquareddifference
SSD_Ill <- sum(diff_Ill^2)
SSD_NEB <- sum(diff_NEB^2)
SSD_NEXT <- sum(diff_NEXT^2)
SSD_Clon <- sum(diff_Clon^2)

SSD_Ill
SSD_NEB
SSD_NEXT
SSD_Clon
SSD<-(data.frame(SSD_Ill, SSD_NEB, SSD_NEXT, SSD_Clon))
test<-(melt(SSD))
barplot(test$value, names.arg = test$Var1)


IllVar<-var(Counts$Illumina)
NEBvar<-var(Counts$NEB)
NEXTvar<-var(Counts$NEXT)
Clonvar<-var(Counts$Clontech)
Var<-data.frame(IllVar, NEBvar, NEXTvar, Clonvar)
save(RMSE, diff, SSD, Var, Acc_ttestStats, file = "/home/carrie/synthetic_miRNA/accuracy_Virus_Final_stats.rda")

