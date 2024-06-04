library(dbplyr)
library(tidyverse)
library(ggrepel)
library(ggBubbles)
library(ggpubr)
library(patchwork)
rm(list=ls())
setwd("D:/课题/project10")

ref_df_SAMPLE_TYPE <- read_tsv("./result/TCGA/0.1癌症突变_样本统计/调试13.4转移状态分组的样本统计.tsv")  %>% rename(SAMPLE_TYPE_1=SAMPLE_TYPE)%>% mutate(SAMPLE_TYPE=ifelse(SAMPLE_TYPE_1=="Primary","Primary","MR"))
ref_df_SAMPLE_TYPE <- ref_df_SAMPLE_TYPE %>% group_by(CANCER_TYPE, order, AGE_GROUP,  ID, SAMPLE_TYPE,nMut) %>% summarise(SAMPLE_TYPE_Sample=sum(SAMPLE_TYPE_Sample)) %>% arrange(order) %>% ungroup %>%
  mutate(percent_Sample=SAMPLE_TYPE_Sample/nMut)


ya_sample <- ref_df_SAMPLE_TYPE %>% group_by(CANCER_TYPE,AGE_GROUP) %>% summarise(Sample_N=sum(SAMPLE_TYPE_Sample)) %>%filter(AGE_GROUP=="YA")
oa_sample <- ref_df_SAMPLE_TYPE %>% group_by(CANCER_TYPE,AGE_GROUP) %>% summarise(Sample_N=sum(SAMPLE_TYPE_Sample)) %>%filter(AGE_GROUP=="OA")
ya_oa_sample <- left_join(ya_sample,oa_sample,by="CANCER_TYPE") %>% mutate(Sample_Ratio=Sample_N.x/Sample_N.y)

ref_df_SAMPLE_TYPE <- left_join(ref_df_SAMPLE_TYPE,ya_oa_sample[,c("CANCER_TYPE","Sample_Ratio")],by="CANCER_TYPE")

cancer <- unique(ref_df_SAMPLE_TYPE$CANCER_TYPE)
ii="PanCancer"
SAMPLE_TYPE_res <-c()
for (ii in cancer) {
  dat_1 <- filter(ref_df_SAMPLE_TYPE,CANCER_TYPE==ii)
  ya_Primary <- dat_1$SAMPLE_TYPE_Sample[which(dat_1$AGE_GROUP=="YA"&dat_1$SAMPLE_TYPE=="Primary")]
  if(length(ya_Primary)==0){ya_Primary <- 0}
  ya_MR <- dat_1$SAMPLE_TYPE_Sample[which(dat_1$AGE_GROUP=="YA"&dat_1$SAMPLE_TYPE=="MR")]
  if(length(ya_MR)==0){ya_MR <- 0}
  
  oa_Primary <- dat_1$SAMPLE_TYPE_Sample[which(dat_1$AGE_GROUP=="OA"&dat_1$SAMPLE_TYPE=="Primary")]
  if(length(oa_Primary)==0){oa_Primary <- 0}
  oa_MR <- dat_1$SAMPLE_TYPE_Sample[which(dat_1$AGE_GROUP=="OA"&dat_1$SAMPLE_TYPE=="MR")]
  if(length(oa_MR)==0){oa_MR <- 0}
  
  dat_mtx <- matrix(c(ya_Primary,ya_MR,oa_Primary,oa_MR),nrow=2)
  colnames(dat_mtx) <- c("YA", "OA") 
  rownames(dat_mtx) <- c("Primary", "MR")
  dat_mtx
  p.val <- fisher.test(x = dat_mtx, alternative = "two.sided")$p.value
  odds_ratio <- fisher.test(x = dat_mtx, alternative = "two.sided")$estimate
  
  SAMPLE_TYPE_res <- rbind(SAMPLE_TYPE_res,data.frame(CANCER_TYPE=ii,P.value=p.val,OR=odds_ratio))
}

ref_df_SAMPLE_TYPE_fisher <- left_join(ref_df_SAMPLE_TYPE,SAMPLE_TYPE_res,by="CANCER_TYPE")

Primary_dat <- ref_df_SAMPLE_TYPE_fisher[which(ref_df_SAMPLE_TYPE_fisher$SAMPLE_TYPE=="Primary"),] %>% reshape2::dcast(CANCER_TYPE+P.value+OR+Sample_Ratio~AGE_GROUP,value.var = "percent_Sample") %>% 
  mutate(Sig=ifelse(as.numeric(P.value)<0.05,"Sig","NonSig"),Label=ifelse(as.numeric(P.value)<0.05,CANCER_TYPE,""))

p1 <- ggplot(data=Primary_dat,aes(x=OA,y=YA, color = Sig)) +
  geom_point(position = position_surround(),aes(size=10))+theme_bw()+
  scale_color_manual(values = c("Sig"="orange","NonSig"="#C9C9CA"))+
  geom_abline(slope = 1, intercept = 0) +
  xlim(c(0, 1)) + ylim(c(0, 1)) + ylab("YA Primary ratio") + xlab("OA Primary ratio")+
  geom_text_repel(aes(label = Label,color="black"),
                  xlim=c(NA, 0),#控制标签在左边
                  direction = 'y',#nudge_x = 200,hjust=1,
                  fontface='italic',max.iter = 3e3,force=1,max.overlaps =200,
                  size = 5)  + theme_bw() 
p1 


MR_dat <- ref_df_SAMPLE_TYPE_fisher[which(ref_df_SAMPLE_TYPE_fisher$SAMPLE_TYPE=="MR"),] %>% reshape2::dcast(CANCER_TYPE+P.value+OR+Sample_Ratio~AGE_GROUP,value.var = "percent_Sample") %>% 
  mutate(Sig=ifelse(as.numeric(P.value)<0.05,"Sig","NonSig"),Label=ifelse(as.numeric(P.value)<0.05,CANCER_TYPE,""))

p2 <- ggplot(data=MR_dat,aes(x=OA,y=YA, color = Sig)) +
  geom_point(position = position_surround(),aes(size=10))+theme_bw()+
  scale_color_manual(values = c("Sig"="orange","NonSig"="#C9C9CA"))+
  geom_abline(slope = 1, intercept = 0) +
  xlim(c(0, 1)) + ylim(c(0, 1))+ ylab("YA MR ratio") + xlab("OA MR ratio")+
  geom_text_repel(aes(label = Label,color="black"),
                  xlim=c(NA, 0),#控制标签在左边
                  direction = 'y',#nudge_x = 200,hjust=1,
                  fontface='italic',max.iter = 3e3,force=1,max.overlaps =200,
                  size = 5)  + theme_bw() 
p2 

p3 <- p1+p2
ggsave("./result/TCGA/0.1癌症突变_样本统计/调试13.2年龄分组间样本类别关联.pdf",p3,width = 10,height = 4)

