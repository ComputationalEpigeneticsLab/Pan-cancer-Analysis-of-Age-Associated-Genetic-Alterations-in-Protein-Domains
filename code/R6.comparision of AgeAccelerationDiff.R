library(tidyverse)
library(dbplyr)
library(vroom)
library(survival)
library(survminer)
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
}

error_proc <- function(x){
  if(inherits(x, "try-error")) {
    x <- NA  # 为结果赋一个新的值
    
  }
  return(x)
}
error_proc_2 <- function(x){
  if(inherits(x, "try-error")) {
    cat(paste("An error occurred in ",cancer,",nrow is ",ii,"\n"))
    next  # 为结果赋一个新的值
    
  }
}

tot_mut_sample <- read_tsv("./result/TCGA/0.1癌症突变_样本统计/3.各个癌症中基因的突变样本.tsv") %>% mutate(submitter_id.samples=substr(Tumor_Sample_Barcode,1,16))
dnam_age <- read_tsv("D:/参考文献/年龄相关癌症基因组/34879281-泛癌解释与年龄相关的分子模式/Epigenetic_age_TCGA.tsv") %>% mutate(SampleID=gsub("[.]","-",SampleID),submitter_id.samples=substr(SampleID,1,16))

mut_sample <- read_tsv("./result/TCGA/1.突变优选/1.位于domain区域内的样本.tsv")
raw <- read_tsv("./result/TCGA/1.突变优选/1.TCGA_raw_domain_singleCancer.tsv") %>% getRawRes %>%
  mutate(Region=paste(Transcript_ID, from, to, `hmm acc`, `hmm name`, clan,sep=";")) %>%
  filter(as.numeric(nRegion)>=3,as.numeric(FDR)<0.05)

library(patchwork)
age_group <- c(40,50,60,70)
p5 <- list()
p5_1 <- list()
for (ii in age_group) {
  dnam_age <- read_tsv("F:/参考文献/年龄相关癌症基因组/34879281-泛癌解释与年龄相关的分子模式/Epigenetic_age_TCGA.tsv") %>%
    mutate(SampleID=gsub("[.]","-",SampleID),submitter_id.samples=substr(SampleID,1,16),
           SAMPLE_TYPE=ifelse(submitter_id.samples %in% substr(mut_sample$Tumor_Sample_Barcode,1,16),"InRegion","Others"),
           AGE_GROUP=ifelse(as.numeric(Age)<=ii,"YA",ifelse(as.numeric(Age)>ii,"OA","None")),
           SAMPLE_TYPE_1 = paste(AGE_GROUP,SAMPLE_TYPE,sep = "-")) 
  
  dnam_age$SAMPLE_TYPE_1 <- factor(dnam_age$SAMPLE_TYPE_1,levels = c("YA-InRegion", "YA-Others", "OA-InRegion", "OA-Others"))
  my_comparisons <- list(c("OA-InRegion", "OA-Others"),c("YA-InRegion", "YA-Others"),c("YA-InRegion","OA-InRegion"))
  p1 <- ggboxplot(dnam_age ,x="SAMPLE_TYPE_1",y="DNAmAge",fill="SAMPLE_TYPE_1",bxp.errorbar = T,bxp.errorbar.width = 0.2,outlier.shape = NA)+
     ylab("DNAmAge")+xlab("")+
    theme(plot.title = element_text(size = 14, face =  "bold"),
          text = element_text(size = 12),
          axis.title = element_text(face="bold"),
          axis.text.x=element_text(size = 11))+
    stat_compare_means(comparisons = my_comparisons, label = "p.format")+
    scale_fill_manual(values=c("YA-InRegion"="#0b6ba7","YA-Others"="#C6C5C1","OA-InRegion"="#b23b28","OA-Others"="#808081"))+theme_classic()
  
  p1_1 <- ggplot(dnam_age,aes(DNAmAge,fill=SAMPLE_TYPE_1, fill=SAMPLE_TYPE_1,color=NULL)) + 
    xlab("DNAmAge") + #xlim(0,0.03) +
    geom_density() +  #alpha = 0.6
    theme_bw()+
    scale_fill_manual(values=c("YA-InRegion"="#0b6ba7","YA-Others"="#C6C5C1","OA-InRegion"="#b23b28","OA-Others"="#808081"))+theme_classic()

  p1+p1_1
   p2 <- ggboxplot(dnam_age ,x="SAMPLE_TYPE_1",y="AgeAccelerationDiff",fill="SAMPLE_TYPE_1",bxp.errorbar = T,bxp.errorbar.width = 0.2,outlier.shape = NA)+
    ylab("AgeAccelerationDiff")+xlab("")+
    theme(plot.title = element_text(size = 14, face =  "bold"),
          text = element_text(size = 12),
          axis.title = element_text(face="bold"),
          axis.text.x=element_text(size = 11))+
    stat_compare_means(comparisons = my_comparisons, label = "p.format")+
    scale_fill_manual(values=c("YA-InRegion"="#0b6ba7","YA-Others"="#C6C5C1","OA-InRegion"="#b23b28","OA-Others"="#808081"))+theme_classic()
  ggsave("./result/TCGA/10.表观年龄加速/4.3区域内与其他样本的DNAmAge差异_40y-AgeAccelerationDiff.pdf",p2,width = 5,height = 4)
  p2_1 <- ggplot(dnam_age,aes(AgeAccelerationDiff,fill=SAMPLE_TYPE_1, fill=SAMPLE_TYPE_1,color=NULL)) + 
    xlab("AgeAccelerationDiff") + #xlim(0,0.03) +
    geom_density() +  #alpha = 0.6
    theme_bw()+
    scale_fill_manual(values=c("YA-InRegion"="#0b6ba7","YA-Others"="#C6C5C1","OA-InRegion"="#b23b28","OA-Others"="#808081"))+theme_classic()  + coord_flip()
  
  p2+p2_1
  ggsave("./result/TCGA/10.表观年龄加速/4.3分布_区域内与其他样本的DNAmAge差异_40y-AgeAccelerationDiff.pdf",p2_1,width = 4,height = 4)
  
  p3 <- ggboxplot(dnam_age ,x="SAMPLE_TYPE_1",y="Age",fill="SAMPLE_TYPE_1",bxp.errorbar = T,bxp.errorbar.width = 0.2,outlier.shape = NA)+
    ylab("Age")+xlab("")+
    theme(plot.title = element_text(size = 14, face =  "bold"),
          text = element_text(size = 12),
          axis.title = element_text(face="bold"),
          axis.text.x=element_text(size = 11))+
    stat_compare_means(comparisons = my_comparisons, label = "p.format")+
    scale_fill_manual(values=c("YA-InRegion"="#0b6ba7","YA-Others"="#C6C5C1","OA-InRegion"="#b23b28","OA-Others"="#808081"))+theme_classic()

  p3_1 <- ggplot(dnam_age,aes(Age,fill=SAMPLE_TYPE_1, fill=SAMPLE_TYPE_1,color=NULL)) + 
    xlab("Age") + #xlim(0,0.03) +
    geom_density() +  #alpha = 0.6
    theme_bw()+
    scale_fill_manual(values=c("YA-InRegion"="#0b6ba7","YA-Others"="#C6C5C1","OA-InRegion"="#b23b28","OA-Others"="#808081"))+theme_classic()
  
  p3+p3_1
  
  p4 <- p1/p2/p3
  p4_1 <- p1_1/p2_1/p3_1
  
  p5[[paste0(ii,"y")]] <- p4
  p5_1[[paste0(ii,"y")]] <- p4_1
}
p6 <- p5[["40y"]]|p5[["50y"]]|p5[["60y"]]|p5[["70y"]]
p6_1 <- p5_1[["40y"]]|p5_1[["50y"]]|p5_1[["60y"]]|p5_1[["70y"]]
ggsave("./result/TCGA/10.表观年龄加速/4.1区域内与其他样本的DNAmAge差异_40y_50y_60y_70y.pdf",p6,width =20,height = 10)
ggsave("./result/TCGA/10.表观年龄加速/4.2分布_区域内与其他样本的DNAmAge差异_40y_50y_60y_70y.pdf",p6_1,width = 20,height = 5)