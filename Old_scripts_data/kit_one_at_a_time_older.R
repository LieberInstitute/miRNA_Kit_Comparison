###Get The DATA##################
miR_counts<-read.table("/home/carrie/UMI/after_UMI_tools_Nextflex_kitcomparison/UMI_NEXTflex_and_all_kits_1000ng/All_combined_1000ng/miR.Counts.csv", header = TRUE, sep = ",")
rownames(miR_counts)<- miR_counts$miRNA#make miRNA rownames
miR_counts<-miR_counts[,2:length(colnames(miR_counts))]#remove miRNA col
miR_counts<-miR_counts[-1,]#remove miRNA col

miR_RPMs<-read.table("/home/carrie/UMI/after_UMI_tools_Nextflex_kitcomparison/UMI_NEXTflex_and_all_kits_1000ng/All_combined_1000ng/miR.RPM.csv", header = TRUE, sep = ",")
rownames(miR_RPMs)<- miR_RPMs$miRNA#make miRNA rownames
miR_RPMs<-miR_RPMs[,2:length(colnames(miR_RPMs))]#remove miRNA col
Pheno<- read.table("/home/carrie/UMI/Pheno_1000ng_full", header = T)
pheno_all <- Pheno[-grep('deduped', Pheno$File), ]
pheno_all$Kit <- relevel(droplevels(pheno_all$Kit), ref = "NextFlex")
load("/home/carrie/UMI/All_intersect_over10RPM.rda")

yGene<-as.matrix(miR_RPMs)

library(rafalib)
library(dendextend)


pdD<-Pheno
miRNAc <- hclust(dist(t(yGene)))
# dendA<-as.dendrogram(miRNAc)
# labels(dendA)<-pdD$TriplicateGroup
# labels_colors(dendA)<-c(colors)
#par(mar=c(bottom,left,top,right)) oma instea of mar is for outer margin
par(mar=c(10,2,3,0))
myplclust(miRNAc, labels=pdD$TriplicateGroup,lab.col=(as.numeric(pdD$Kit)))
newcol<-recode(colors,"1= '5'")
myplclust(miRNAc, labels=pdD$TriplicateGroup,lab.col=(as.numeric(newcol)))
colors = as.numeric(pdD$Kit)
newcol<-replace(colors, 4, '5')
library(car) 

# SchoolData$Grade<-recode(SchoolData$Grade,"5=6;6=7")
# myplclust(miRNAc, labels=pdD$TriplicateGroup,lab.col=(colors))
# 
# myplclust(miRNAc, labels=pdD$TriplicateGroup,lab.col=as.numeric(pdD$Kit, colors = "red","yellow", "pink", "green"))
# miRNAc <- hclust(dist(t(yGene[,which(pd$Dataset =="Repro")])))
# pdrepro<-pd[which(pd$Dataset=="Repro"),]
# 
# myplclust(miRNAc, labels=pdrepro$TriplicateGroup,lab.col=as.numeric(pdrepro$Kit))
# miRNAc <- hclust(dist(t(yGene[,which(pd$Dataset =="Accur")])))
# 
# pdaccur<-pd[which(pd$Dataset=="Accur"),]
# myplclust(miRNAc, labels=pdaccur$TriplicateGroup,lab.col=as.numeric(pdaccur$Kit))

###An MA-plot is a plot of log-intensity ratios (M-values) versus log-intensity averages (A-values)####

###need to add more samples...

MA <- new("MAList")
counts<-miR_RPMs[grep('NEXT|deduped', colnames(miR_RPMs))]
### check for zeros
row_sub = apply(counts, 1, function(row) all(row !=0 ))
counts<-counts[row_sub,]
sample10<-(counts[grep("10",colnames(counts))])+1
MA$M<-log2(sample10$deduped_10_bam2fastq.fq/sample10$mm0_NEXT_trimmed.10_R1.fq)
MA$A<-log2(sqrt(sample10$deduped_10_bam2fastq.fq*sample10$mm0_NEXT_trimmed.10_R1.fq))
limma::plotMA(MA,main="sample10", cex = 2)
ma.plot( rowMeans(log2(sample10)), log2(sample10[, 1])-log2(sample10[, 2]), cex=1, main = "deduped vs raw")
deduped2<-counts[1:2]+1
ma.plot( rowMeans(log2(deduped2)), log2(deduped2[, 1])-log2(deduped2[, 2]), cex=1, main = "two deduped")
raw2<- counts[4:5]+1
ma.plot( rowMeans(log2(raw2)), log2(raw2[, 1])-log2(raw2[, 2]), cex=1, main = "two raw")

plot(log2(raw2[, 1])-log2(raw2[, 2]), log2(deduped2[, 1])-log2(deduped2[, 2]))

plot(log2(raw2[,1])-log2(deduped2[,1]))

MA$A <-log2(sqrt(miR_counts$mm0_NEXT_trimmed.10_R1.fq* miR_counts$deduped_10_bam2fastq.fq)+1)
MA$M <- log2((miR_counts$mm0_NEXT_trimmed.10_R1.fq/miR_counts$deduped_10_bam2fastq.fq)+1)
limma::plotMA(sample10,main="sample10", cex = 2)
limma::plotMA(sample10,main="sample10", cex = 2)
MA$A <-rowMeans(data.frame(miR_RPMs$mm0_NEXT_trimmed.11_R1.fq, miR_RPMs$deduped_11_bam2fastq.fq))
MA$M <- miR_RPMs$mm0_NEXT_trimmed.11_R1.fq - miR_RPMs$deduped_11_bam2fastq.fq
limma::plotMA(MA,main="sample11")
MA$A <-rowMeans(data.frame(miR_RPMs$mm0_NEXT_trimmed.12_R1.fq, miR_RPMs$deduped_12_bam2fastq.fq))
MA$M <- miR_RPMs$mm0_NEXT_trimmed.12_R1.fq - miR_RPMs$deduped_12_bam2fastq.fq
limma::plotMA(MA,main="sample12")
MA$A <-rowMeans(data.frame(miR_RPMs$deduped_12_bam2fastq.fq, miR_RPMs$mm0_NEXT_trimmed.12_R1.fq))
MA$M <- miR_RPMs$deduped_12_bam2fastq.fq - miR_RPMs$mm0_NEXT_trimmed.12_R1.fq
limma::plotMA(MA,main="sample12")

single <- function(kit, voom = FALSE) {
    dat <- RPM_all[, grep(kit, tolower(colnames(RPM_all)))]
    p <- pheno_all[grep(kit, tolower(pheno_all[, 1])), ]
    design <- model.matrix(~p$Batch)
    if(!voom) {
        dat <- log2(dat + 1)
    } else {
        dge <- DGEList(counts = dat)
        dge <- calcNormFactors(dge)
        dat <- voom(dge, design, plot = TRUE)
        Sys.sleep(3)
    }
    f <- lmFit(dat, design)
    f <- eBayes(f)
}

kits <- c('illumina', 'next', 'neb', 'clon')

####Choose the data###############
### For starting with RPM
RPM_all<-miR_RPMs[which(rownames(miR_RPMs) %in% All_10), -grep('deduped|dedupped', colnames(miR_RPMs))]
fits <- lapply(kits, single)

## For starting with counts
RPM_all<-miR_counts[which(rownames(miR_counts) %in% All_10), -grep('deduped|dedupped', colnames(miR_counts))]
fits <- lapply(kits, single, voom = TRUE)
##########pvalues etc.#####
pval_df_single <- as.data.frame(do.call(cbind, lapply(fits, function(f) { f$p.value[, 2] })))
colnames(pval_df_single) <- c('Illumina', 'NEXTflex', 'NEB', 'Clontech')
qval_df_single <- as.data.frame(sapply(pval_df_single, p.adjust, method = 'fdr'))
sapply(qval_df_single, function(x) { table(x > 0.05 )}) # not sig accross batch is true
round(sapply(qval_df_single, function(x) { table(x > 0.05 )}) / nrow(qval_df_single) * 100, 2)
qval_df <- -log10(qval_df_single)
logging<- function(kit){
  dat <- rowMeans(log2(RPM_all[, grep(kit, tolower(colnames(RPM_all)))]+1))}
exp <- lapply(kits, logging)
make_plot <- function(All_test_single, coef = FALSE) {
    plot1001<-ggplot(data = All_test_single, aes(x = variable, y = value, color= variable))+geom_boxplot(aes(fill = variable, alpha =.7))+
      ggtitle(label = "Inconsistency Across Batch")+ 
      labs(y = ifelse(coef, "Coef for batch effect \n on indv. miRNAs", "-log10(q) for batch effect \n on indv. miRNAs"))+
      theme(axis.title.x = element_text(size =0), 
            plot.title = element_text(size = 60, face = "bold", hjust = 0.5), 
            axis.text.x = element_text(size = 40),
            axis.text.y = element_text(size = 40), 
            axis.title.y = element_text(size =35),
            legend.position= "none") 
    plot1001 +scale_fill_manual(values=c("firebrick3", "blue","green3", "black"))+ scale_color_manual(values=rep("black", 4))
}

make_plot(melt(-log10(qval_df_single)))

####relationship of expression and sig batch effect##########
Expr_plot_Data <- function(All_test_single, exp) {
  All_test_single <- -log10(All_test_single)
  Alldata<-data.frame(All_test_single,exp)
  colnames(Alldata)<-c("Illumina_q", "NEXTflex_q", "NEB_q", "Clontech_q","Illumina_exp", "NEXTflex_exp", "NEB_exp", "Clontech_exp")
  print(Alldata)
}
expdata<- Expr_plot_Data(qval_df, exp)
ggplot(data = expdata, aes(x = expdata$Illumina_q, y = expdata$Illumina_exp))+ 
geom_smooth(method = "lm", se=FALSE, color="red", formula = y ~ x)  +geom_point(aes(expdata$Illumina_q,expdata$Illumina_exp))

plot(-log10(qval_df$Illumina)~ rowMeans(log2(RPM_all[grep("Illumina", colnames(RPM_all))]+1)))
abline(lm(-log10(qval_df_single$Illumina)~ rowMeans(log2(RPM_all[grep("Illumina", colnames(RPM_all))]+1))), col="red")
plot(-log10(qval_df$NEB)~ rowMeans(log2(RPM_all[grep("NEB", colnames(RPM_all))]+1)))
abline(lm(-log10(qval_df$NEB)~ rowMeans(log2(RPM_all[grep("NEB", colnames(RPM_all))]+1))), col="red")
plot(-log10(qval_df_single$Clontech) ~ rowMeans(log2(RPM_all[grep("Clontech", colnames(RPM_all))]+1)))
abline(lm(-log10(qval_df_single$Clontech)~ rowMeans(log2(RPM_all[grep("Clontech", colnames(RPM_all))]+1))), col="red")
plot(-log10(qval_df_single$NEXTflex) ~ rowMeans(log2(RPM_all[grep("NEXT", colnames(RPM_all))]+1)))
abline(lm(-log10(qval_df_single$NEXTflex)~ rowMeans(log2(RPM_all[grep("NEXT", colnames(RPM_all))]+1))), col="red")

summary(lm(-log10(qval_df_single$Illumina)~ rowMeans(log2(RPM_all[grep("Illumina", colnames(RPM_all))]+1))))
summary(lm(-log10(qval_df_single$NEB)~ rowMeans(log2(RPM_all[grep("NEB", colnames(RPM_all))]+1))))
summary(lm(-log10(qval_df_single$Clontech)~ rowMeans(log2(RPM_all[grep("Clontech", colnames(RPM_all))]+1))))#####signifcant relationship
summary(lm(-log10(qval_df_single$NEXT)~ rowMeans(log2(RPM_all[grep("NEXT", colnames(RPM_all))]+1))))#####signifcant relationship
##makes sense... if batch effect is somewhat due to adapter ligation efficency then this might happen to all miRNAs in general irregardless of expression for NEB and NEXTflex??? maybe??

###ttests_single####
ttests_single <-list()
ttests_single[[1]]<-t.test(-log10(qval_df_single$NEXTflex), -log10(qval_df_single$Illumina), paired = TRUE)
ttests_single[[2]]<-t.test(-log10(qval_df_single$NEXTflex), -log10(qval_df_single$NEB), paired = TRUE)
ttests_single[[3]]<-t.test(-log10(qval_df_single$NEXTflex), -log10(qval_df_single$Clontech), paired = TRUE)
ttests_single[[4]]<-t.test(-log10(qval_df_single$Illumina), -log10(qval_df_single$NEB), paired = TRUE)
ttests_single[[5]]<-t.test(-log10(qval_df_single$Clontech), -log10(qval_df_single$NEB), paired = TRUE)
ttests_single[[6]]<-t.test(-log10(qval_df_single$Illumina), -log10(qval_df_single$Clontech), paired = TRUE)
# extract your values using `sapply`
Acc_ttests_singletats<-sapply(ttests_single, function(x) {
  c(t =format(x$statistic, digits = 2),
    df = format(x$parameter, digits = 0),
    p.value = format(x$p.value, scientific = FALSE, digits = 2),
    bonferroni = format(.05/6, digits = 2),
    sig = ifelse(x$p.value<(.05/6), "yes", "no"))
})
colnames(Acc_ttests_singletats)<-sapply(ttests_single, function(x) {
  c(test = x$data.name)
})

Acc_ttests_singletats



sig_Table_single<-round(sapply(qval_df_single, function(x) { table(x < 0.05 )}) / nrow(qval_df_single) * 100, 2)
rownames(sig_Table_single)<-c("No Batch Effect", "Sig Batch Effect")
plot(sig_Table_single)
DF1_single <- melt(sig_Table_single)

ggplot(DF1_single, aes(x = Var2, y = value, fill = Var1)) +
  geom_bar(stat = "identity") +scale_fill_manual(values=c("green3", "black"))+
  guides(fill=guide_legend(title=""))+
  ggtitle(label = "Percentage of miRNAs with significant batch effect")+ 
  ylab("Percentage of miRNAs")+
  xlab("Kit")+
  theme(axis.title.x = element_text(size =0), 
        plot.title = element_text(size = 50, face = "bold", hjust = 0.5), 
        axis.text.x = element_text(size = 40),
        axis.text.y = element_text(size = 40), 
        axis.title.y = element_text(size =35),
        legend.text=element_text(size=20))



#To make coeffecient plot####
# coef_df_single <- as.data.frame(do.call(cbind, lapply(fits, function(f) { f$coefficients[, 2] })))
# colnames(coef_df_single) <- c('Illumina', 'NEXTflex', 'NEB', 'Clontech')
# make_plot(melt(abs(coef_df_single)), coef = TRUE)


####Pairwise comparisons
### For starting with RPM
#####################################loop

Filenames<-unique(sub("\\_.*","",pheno_all$File))
Kitpairs<-(combn(Filenames, 2, FUN = NULL))
combName = function(pairs) {
  paste0(pairs[1], "|", pairs[2])
}
Kitpairname <- apply(Kitpairs, 2, FUN =combName)

for (pair in Kitpairname){
  x<-data.frame(name=pair)
  int<-gsub(pattern="[[:punct:]]", x$name, replacement="_")
  y = paste0(int, "_10")
  b = paste0(paste(strsplit(int, '_', fixed=T)[[1]][c(2,1)], collapse='_'), "_10")
  if (exists(as.character(substitute(y)))) {
    testdata <- get(y)
  } else {testdata <-get(b)}
  RPM_all<- miR_RPMs[which(rownames(miR_RPMs) %in% testdata), grep(x$name, colnames(miR_RPMs))]
  kits2 <- kits[grep(tolower(x$name), kits)]
  #Sys.sleep(50)
  fits <- lapply(kits2, single)
  #Sys.sleep(50)
  pval_df_single <- as.data.frame(do.call(cbind, lapply(fits, function(f) { f$p.value[, 2] })))
  #Sys.sleep(50)
  colnames(pval_df_single)<-c(strsplit(int, '_', fixed=T)[[1]][1], strsplit(int, '_', fixed=T)[[1]][2])
  qval_df_single <- as.data.frame(sapply(pval_df_single, p.adjust, method = 'fdr'))
  sapply(qval_df_single, function(x) { table(x > 0.05 )}) # not sig accross batch is true
  round(sapply(qval_df_single, function(x) { table(x > 0.05 )}) / nrow(qval_df_single) * 100, 2)
  qval_df <- -log10(qval_df_single)
  make_plot <- function(All_test_single, coef = FALSE) {
    plot1001<-ggplot(data = All_test_single, aes(x = variable, y = value, color= variable))+geom_boxplot(aes(fill = variable, alpha =.7))+
      ggtitle(label = "Inconsistency Across Batch")+ 
      labs(y = ifelse(coef, "Coef for batch effect \n on indv. miRNAs", "-log10(q) for batch effect \n on indv. miRNAs"))+
      theme(axis.title.x = element_text(size =0), 
            plot.title = element_text(size = 60, face = "bold", hjust = 0.5), 
            axis.text.x = element_text(size = 40),
            axis.text.y = element_text(size = 40), 
            axis.title.y = element_text(size =35),
            legend.position= "none") 
    plot1001 +scale_fill_manual(values=c("firebrick3", "blue"))+ scale_color_manual(values=rep("black", 4))
  }
  
  make_plot(melt(-log10(qval_df_single)))
}


t.test(-log10(qval_df_single$NEXTflex), -log10(qval_df_single$Illumina), paired = TRUE)
#########################################################Pairwise####log2#####


kits <- c( 'next', 'neb','illumina', 'clontech')

x = data.frame(name ='NEXT|Illumina')
x = data.frame(name ='NEXT|NEB')
x = data.frame(name ='NEXT|Clontech')
x = data.frame(name ='Illumina|Clontech')
x = data.frame(name ='NEB|Clontech')
x = data.frame(name ='NEB|Illumina')

x = data.frame(name ='deduped|Illumina')
x = data.frame(name ='deduped|NEB')
x = data.frame(name ='deduped|Clontech')

Pheno<- read.table("/home/carrie/UMI/Pheno_1000ng_full", header = T)
pheno_all <- Pheno[-grep('NEXT', Pheno$File), ]
pheno_all$Kit <- relevel(droplevels(pheno_all$Kit), ref = "NextFlex_deduped")
kits <- c('deduped', 'neb', 'illumina','clontech')

load("/home/carrie/UMI/All_intersect_over10RPM_deduped.rda")
miR_RPMs<-read.table("/home/carrie/UMI/after_UMI_tools_Nextflex_kitcomparison/UMI_NEXTflex_and_all_kits_1000ng/All_combined_1000ng/miR.RPM.csv", header = TRUE, sep = ",")
rownames(miR_RPMs)<- miR_RPMs$miRNA#make miRNA rownames
miR_RPMs<-miR_RPMs[,2:length(colnames(miR_RPMs))]#remove miRNA col
colnames(miR_RPMs)<-gsub("dedupped", "deduped", colnames(miR_RPMs))

int<-gsub(pattern="[[:punct:]]", x$name, replacement="_")
y = paste0(int, "_10")
b = paste0(paste(strsplit(int, '_', fixed=T)[[1]][c(2,1)], collapse='_'), "_10")
if (exists(as.character(substitute(y)))) {
  testdata <- get(y)
} else {testdata <-get(b)}
RPM_all<- miR_RPMs[which(rownames(miR_RPMs) %in% testdata), grep(x$name, colnames(miR_RPMs))]
RPM_all<- RPM_all[-(which(rownames(RPM_all) %in% All_10)),]
#RPM_all<- RPM_all[-(which(rownames(RPM_all) %in% All_10_deduped)),]
kits2 <- kits[grep(tolower(x$name), kits)]####need to make just first part of name elements
#kits2 <- kits[grep("neb|clon", kits)]

fits <- lapply(kits2, single)
pval_df_single <- as.data.frame(do.call(cbind, lapply(fits, function(f) { f$p.value[, 2] })))
kits2
colnames(pval_df_single)<-c(strsplit(int, '_', fixed=T)[[1]][1], strsplit(int, '_', fixed=T)[[1]][2])
colnames(pval_df_single)
qval_df_single <- as.data.frame(sapply(pval_df_single, p.adjust, method = 'fdr'))
sapply(qval_df_single, function(x) { table(x > 0.05 )}) # not sig accross batch is true
round(sapply(qval_df_single, function(x) { table(x > 0.05 )}) / nrow(qval_df_single) * 100, 2)
qval_df <- -log10(qval_df_single)
make_plot <- function(All_test_single, coef = FALSE) {
  plot1001<-ggplot(data = All_test_single, aes(x = variable, y = value, color= variable))+geom_boxplot(aes(fill = variable, alpha =.7))+
    ggtitle(label = "Inconsistency Across Batch")+ 
    labs(y = ifelse(coef, "Coef for batch effect \n on indv. miRNAs", "-log10(q) for batch effect \n on indv. miRNAs"))+
    theme(axis.title.x = element_text(size =0), 
          plot.title = element_text(size = 60, face = "bold", hjust = 0.5), 
          axis.text.x = element_text(size = 40),
          axis.text.y = element_text(size = 40), 
          axis.title.y = element_text(size =35),
          legend.position= "none") 
  plot1001 +scale_fill_manual(values=c("firebrick3", "blue"))+ scale_color_manual(values=rep("black", 4))
}

make_plot(melt(-log10(qval_df_single)))
t.test(-log10(qval_df_single[,1]), -log10(qval_df_single[,2]), paired = TRUE)
#########################################################Pairwise####voom#####
colnames(miR_counts)<-gsub("dedupped", "deduped", colnames(miR_counts))
Pheno<- read.table("/home/carrie/UMI/Pheno_1000ng_full", header = T)
pheno_all <- Pheno[-grep('deduped|dedupped', Pheno$File), ]
pheno_all$Kit <- relevel(droplevels(pheno_all$Kit), ref = "NextFlex")


kits <- c( 'next', 'neb','illumina', 'clontech')

x = data.frame(name ='NEXT|Illumina')
x = data.frame(name ='NEXT|NEB')
x = data.frame(name ='NEXT|Clontech')
x = data.frame(name ='Illumina|Clontech')
x = data.frame(name ='NEB|Clontech')
x = data.frame(name ='NEB|Illumina')


# kits <- c('deduped', 'neb', 'illumina','clontech')
# load("/home/carrie/UMI/All_intersect_over10RPM_deduped.rda")
# colnames(miR_RPMs)<-gsub("dedupped", "deduped", colnames(miR_RPMs))
#x = data.frame(name ='deduped|Illumina')
#x = data.frame(name ='deduped|NEB')
#x = data.frame(name ='deduped|Clontech')

# Pheno<- read.table("/home/carrie/UMI/Pheno_1000ng_full", header = T)
# pheno_all <- Pheno[-grep('NEXT', Pheno$File), ]
# pheno_all$Kit <- relevel(droplevels(pheno_all$Kit), ref = "NextFlex_deduped")
# 
# kits <- c('illumina', 'deduped', 'neb', 'clontech')
# load("/home/carrie/UMI/All_intersect_over10RPM_deduped.rda")
# colnames(miR_RPMs)<-gsub("dedupped", "deduped", colnames(miR_RPMs))
int<-gsub(pattern="[[:punct:]]", x$name, replacement="_")
y = paste0(int, "_10")
b = paste0(paste(strsplit(int, '_', fixed=T)[[1]][c(2,1)], collapse='_'), "_10")
if (exists(as.character(substitute(y)))) {
  testdata <- get(y)
} else {testdata <-get(b)}
#testdata<- NEB_Illumina_10
RPM_all<- miR_counts[which(rownames(miR_counts) %in% testdata), grep(x$name, colnames(miR_counts))]
#RPM_all<- RPM_all[-(which(rownames(RPM_all) %in% All_10)),]
#RPM_all<- RPM_all[-(which(rownames(RPM_all) %in% All_10_deduped)),]
kits2 <- kits[grep(tolower(x$name), kits)]####need to make just first part of name elements
#kits2 <- kits[grep("neb|clon", kits)]

fits <- lapply(kits2, single, voom = TRUE)
pval_df_single <- as.data.frame(do.call(cbind, lapply(fits, function(f) { f$p.value[, 2] })))
kits2
colnames(pval_df_single)<-c(strsplit(int, '_', fixed=T)[[1]][1], strsplit(int, '_', fixed=T)[[1]][2])
colnames(pval_df_single)
qval_df_single <- as.data.frame(sapply(pval_df_single, p.adjust, method = 'fdr'))
sapply(qval_df_single, function(x) { table(x > 0.05 )}) # not sig accross batch is true
round(sapply(qval_df_single, function(x) { table(x > 0.05 )}) / nrow(qval_df_single) * 100, 2)
qval_df <- -log10(qval_df_single)
make_plot <- function(All_test_single, coef = FALSE) {
  plot1001<-ggplot(data = All_test_single, aes(x = variable, y = value, color= variable))+geom_boxplot(aes(fill = variable, alpha =.7))+
    ggtitle(label = "Inconsistency Across Batch")+ 
    labs(y = ifelse(coef, "Coef for batch effect \n on indv. miRNAs", "-log10(q) for batch effect \n on indv. miRNAs"))+
    theme(axis.title.x = element_text(size =0), 
          plot.title = element_text(size = 60, face = "bold", hjust = 0.5), 
          axis.text.x = element_text(size = 40),
          axis.text.y = element_text(size = 40), 
          axis.title.y = element_text(size =35),
          legend.position= "none") 
  plot1001 +scale_fill_manual(values=c("firebrick3", "blue"))+ scale_color_manual(values=rep("black", 4))
}

make_plot(melt(-log10(qval_df_single)))
ttests_single<-list()
ttests_single[[1]]<-t.test(-log10(qval_df_single[,1]), -log10(qval_df_single[,2]), paired = TRUE)


Acc_ttests_singletats<-sapply(ttests_single, function(x) {
  c(t =format(x$statistic, digits = 2),
    df = format(x$parameter, digits = 0),
    p.value = format(x$p.value, scientific = FALSE, digits = 2),
    bonferroni = format(.05/6, digits = 2),
    sig = ifelse(x$p.value<(.05/6), "yes", "no"))
})
colnames(Acc_ttests_singletats)<-sapply(ttests_single, function(x) {
  c(test = x$data.name)
})

Acc_ttests_singletats