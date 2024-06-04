library(tidyverse)
library(dbplyr)
library(vroom)
setwd("D:/课题/project10")
rm(list=ls())
getRawRes <- function(x){
  a <- filter(x,as.numeric(nRegion)>=3) %>% 
    arrange(as.numeric(P.value)) %>% 
    group_by(CANCER_TYPE,AGE_GROUP) %>%
    mutate(FDR=p.adjust(as.numeric(P.value),method = "BH")) %>% ungroup
  
  b <- filter(x,as.numeric(nRegion)<3) %>%
    mutate(P.value=1,FDR=1)
  c <- rbind(a,b)
  # rerun(c)
}
hg19_gtf <- read.table("D:/数据/基因组/GENCODE/gencode.v19.gtf.basic.txt",header=T,stringsAsFactors=F,sep='\t')
hg19_gtf <- filter(hg19_gtf,Gene_Type=="protein_coding",Type=="gene") 
gencode.hg19.gene <- unique(hg19_gtf$Hugo_Symbol)
raw <- read_tsv("./result/TCGA/1.突变优选/1.TCGA_raw_domain_singleCancer.tsv") %>% getRawRes %>% 
  mutate(Region=paste(Transcript_ID, from, to, `hmm acc`, `hmm name`, clan,sep=";")) %>%
  filter(as.numeric(nRegion)>=3,as.numeric(FDR)<0.05)

#YA
raw_ya <- filter(raw ,AGE_GROUP=="YA")
cgc <- read_tsv("D:/数据/CGC数据/CGC.tsv")
cgcSymbol <- unique(cgc$`Gene Symbol`)

both_symbol <- intersect(raw_ya$Hugo_Symbol,cgcSymbol)

overlap100000 <- c()
for (i in 1:100000) {
  if(i%%10000==0)print(i)
  sampleOpt <- unique(gencode.hg19.gene)[sample(1:length(unique(gencode.hg19.gene)),length(both_symbol))]
  sampleOptCGC <- intersect(sampleOpt,cgcSymbol)
  overlap100000 <- rbind(overlap100000,data.frame(ID=i,overlap=length(sampleOptCGC)))
}
options(scipen=200)
EP <- length(which(overlap100000$overlap>length(both_symbol)))/100000
O <- length(both_symbol)
M1 <- length(unique(raw_ya$Hugo_Symbol))
M2 <- length(cgcSymbol)
N <- length(unique(gencode.hg19.gene))
OE <- O/((M1*M2)/N)
histData <- overlap100000$overlap
pdf("./result/TCGA/2.优选基因的重要性评估/1.1YA_domain优选基因与CGC.pdf",width = 8,height = 5)
h <-hist(histData, breaks=30, col="#C89F62", xlab="Numer of overlapping genes",
         main=paste("p is",EP,";intersect is",length(both_symbol),";IDR genes:",M1,";CGC genes:",M2,";OE:",OE ,sep=' '),)
xfit<-seq(min(histData),max(histData),length=100)
yfit<-dnorm(xfit,mean=mean(histData),sd=sd(histData))
yfit <- yfit*diff(h$mids[1:2])*length(histData)
lines(xfit, yfit, col="red", lwd=2) 
dev.off()



cm <- read_tsv("D:/数据/CancerMiner数据/cancermine_collated.tsv")
cmSymbol <- unique(cm$gene_normalized)

both_symbol <- intersect(raw_ya$Hugo_Symbol,cmSymbol)
overlap100000 <- c()
for (i in 1:100000) {
  if(i%%10000==0)print(i)
  sampleOpt <- unique(gencode.hg19.gene)[sample(1:length(unique(gencode.hg19.gene)),length(both_symbol))]
  sampleOptCGC <- intersect(sampleOpt,cmSymbol)
  overlap100000 <- rbind(overlap100000,data.frame(ID=i,overlap=length(sampleOptCGC)))
}
options(scipen=300)
EP <- length(which(overlap100000$overlap>length(both_symbol)))/100000
O <- length(both_symbol)
M1 <- length(unique(raw_ya$Hugo_Symbol))
M2 <- length(cmSymbol)
N <- length(unique(gencode.hg19.gene))
OE <- O/((M1*M2)/N)
histData <- overlap100000$overlap
pdf("./result/TCGA/2.优选基因的重要性评估/1.1YA_domain优选基因与CancerMine.pdf",width = 8,height = 5)
h <-hist(histData, breaks=30, col="#C89F62", xlab="Numer of overlapping genes",
         main=paste("p is",EP,";intersect is",length(both_symbol),";IDR genes:",M1,";CancerMine genes:",M2,";OE:",OE ,sep=' '),)
xfit<-seq(min(histData),max(histData),length=100)
yfit<-dnorm(xfit,mean=mean(histData),sd=sd(histData))
yfit <- yfit*diff(h$mids[1:2])*length(histData)
lines(xfit, yfit, col="red", lwd=2) 
dev.off()


#OA
raw_oa <- filter(raw ,AGE_GROUP=="OA")
cgc <- read_tsv("D:/数据/CGC数据/CGC.tsv")
cgcSymbol <- unique(cgc$`Gene Symbol`)

both_symbol <- intersect(raw_oa$Hugo_Symbol,cgcSymbol)

overlap100000 <- c()
for (i in 1:100000) {
  if(i%%10000==0)print(i)
  sampleOpt <- unique(gencode.hg19.gene)[sample(1:length(unique(gencode.hg19.gene)),length(both_symbol))]
  sampleOptCGC <- intersect(sampleOpt,cgcSymbol)
  overlap100000 <- rbind(overlap100000,data.frame(ID=i,overlap=length(sampleOptCGC)))
}
options(scipen=200)
EP <- length(which(overlap100000$overlap>length(both_symbol)))/100000
O <- length(both_symbol)
M1 <- length(unique(raw_oa$Hugo_Symbol))
M2 <- length(cgcSymbol)
N <- length(unique(gencode.hg19.gene))
OE <- O/((M1*M2)/N)
histData <- overlap100000$overlap
pdf("./result/TCGA/2.优选基因的重要性评估/1.1OA_domain优选基因与CGC.pdf",width = 8,height = 5)
h <-hist(histData, breaks=30, col="#C89F62", xlab="Numer of overlapping genes",
         main=paste("p is",EP,";intersect is",length(both_symbol),";IDR genes:",M1,";CGC genes:",M2,";OE:",OE ,sep=' '),)
xfit<-seq(min(histData),max(histData),length=100)
yfit<-dnorm(xfit,mean=mean(histData),sd=sd(histData))
yfit <- yfit*diff(h$mids[1:2])*length(histData)
lines(xfit, yfit, col="red", lwd=2) 
dev.off()



cm <- read_tsv("D:/数据/CancerMiner数据/cancermine_collated.tsv")
cmSymbol <- unique(cm$gene_normalized)

both_symbol <- intersect(raw_oa$Hugo_Symbol,cmSymbol)
overlap100000 <- c()
for (i in 1:100000) {
  if(i%%10000==0)print(i)
  sampleOpt <- unique(gencode.hg19.gene)[sample(1:length(unique(gencode.hg19.gene)),length(both_symbol))]
  sampleOptCGC <- intersect(sampleOpt,cmSymbol)
  overlap100000 <- rbind(overlap100000,data.frame(ID=i,overlap=length(sampleOptCGC)))
}
options(scipen=300)
EP <- length(which(overlap100000$overlap>length(both_symbol)))/100000
O <- length(both_symbol)
M1 <- length(unique(raw_oa$Hugo_Symbol))
M2 <- length(cmSymbol)
N <- length(unique(gencode.hg19.gene))
OE <- O/((M1*M2)/N)
histData <- overlap100000$overlap
pdf("./result/TCGA/2.优选基因的重要性评估/1.1OA_domain优选基因与CancerMine.pdf",width = 8,height = 5)
h <-hist(histData, breaks=30, col="#C89F62", xlab="Numer of overlapping genes",
         main=paste("p is",EP,";intersect is",length(both_symbol),";IDR genes:",M1,";CancerMine genes:",M2,";OE:",OE ,sep=' '),)
xfit<-seq(min(histData),max(histData),length=100)
yfit<-dnorm(xfit,mean=mean(histData),sd=sd(histData))
yfit <- yfit*diff(h$mids[1:2])*length(histData)
lines(xfit, yfit, col="red", lwd=2) 
dev.off()
