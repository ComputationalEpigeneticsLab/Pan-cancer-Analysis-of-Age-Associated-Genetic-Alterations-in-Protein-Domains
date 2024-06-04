library(dbplyr)
library(tidyverse)
library(ggrepel)
library(ggBubbles)
library(ggpubr)
library(patchwork)
rm(list=ls())
setwd("D:/课题/project10")

ref_df_gender <- read_tsv("./result/TCGA/0.1癌症突变_样本统计/调试13.3年龄分组的性别统计.tsv")
ya_sample <- ref_df_gender %>% group_by(CANCER_TYPE,AGE_GROUP) %>% summarise(Sample_N=sum(Gender_Sample)) %>%filter(AGE_GROUP=="YA")
oa_sample <- ref_df_gender %>% group_by(CANCER_TYPE,AGE_GROUP) %>% summarise(Sample_N=sum(Gender_Sample)) %>%filter(AGE_GROUP=="OA")
ya_oa_sample <- left_join(ya_sample,oa_sample,by="CANCER_TYPE") %>% mutate(Sample_Ratio=Sample_N.x/Sample_N.y)

ref_df_gender <- left_join(ref_df_gender,ya_oa_sample[,c("CANCER_TYPE","Sample_Ratio")],by="CANCER_TYPE")

cancer <- unique(ref_df_gender$CANCER_TYPE)
ii="PanCancer"
gender_res <-c()
for (ii in cancer) {
  dat_1 <- filter(ref_df_gender,CANCER_TYPE==ii)
  ya_male <- dat_1$Gender_Sample[which(dat_1$AGE_GROUP=="YA"&dat_1$Gender=="male")]
  if(length(ya_male)==0){ya_male <- 0}
  ya_female <- dat_1$Gender_Sample[which(dat_1$AGE_GROUP=="YA"&dat_1$Gender=="female")]
  if(length(ya_female)==0){ya_female <- 0}
  
  oa_male <- dat_1$Gender_Sample[which(dat_1$AGE_GROUP=="OA"&dat_1$Gender=="male")]
  if(length(oa_male)==0){oa_male <- 0}
  oa_female <- dat_1$Gender_Sample[which(dat_1$AGE_GROUP=="OA"&dat_1$Gender=="female")]
  if(length(oa_female)==0){oa_female <- 0}
  
  dat_mtx <- matrix(c(ya_male,ya_female,oa_male,oa_female),nrow=2)
  colnames(dat_mtx) <- c("YA", "OA") 
  rownames(dat_mtx) <- c("male", "female")
  dat_mtx
  p.val <- fisher.test(x = dat_mtx, alternative = "two.sided")$p.value
  odds_ratio <- fisher.test(x = dat_mtx, alternative = "two.sided")$estimate
  
  gender_res <- rbind(gender_res,data.frame(CANCER_TYPE=ii,P.value=p.val,OR=odds_ratio))
}

ref_df_gender_fisher <- left_join(ref_df_gender,gender_res,by="CANCER_TYPE")

#male，缺少："CESC" "OV"   "UCEC" "UCS" 癌症
male_dat <- ref_df_gender_fisher[which(ref_df_gender_fisher$Gender=="male"),] %>% reshape2::dcast(CANCER_TYPE+P.value+OR+Sample_Ratio~AGE_GROUP,value.var = "percent_Sample") %>% 
  mutate(Sig=ifelse(as.numeric(P.value)<0.05,"Sig","NonSig"),Label=ifelse(as.numeric(P.value)<0.05,CANCER_TYPE,""))

p1 <- ggplot(data=male_dat,aes(x=OA,y=YA, color = Sig)) +
  geom_point(position = position_surround(),aes(size=10))+theme_bw()+
  scale_color_manual(values = c("Sig"="orange","NonSig"="#C9C9CA"))+
  geom_abline(slope = 1, intercept = 0) +
  xlim(c(0, 1)) + ylim(c(0, 1)) + ylab("YA male ratio") + xlab("OA male ratio")+
  geom_text_repel(aes(label = Label,color="black"),
                  xlim=c(NA, 0),#控制标签在左边
                  direction = 'y',#nudge_x = 200,hjust=1,
                  fontface='italic',max.iter = 3e3,force=1,max.overlaps =200,
                  size = 5)  + theme_bw() 
p1 


#female，缺少"PRAD" "TGCT"
female_dat <- ref_df_gender_fisher[which(ref_df_gender_fisher$Gender=="female"),] %>% reshape2::dcast(CANCER_TYPE+P.value+OR+Sample_Ratio~AGE_GROUP,value.var = "percent_Sample") %>% 
  mutate(Sig=ifelse(as.numeric(P.value)<0.05,"Sig","NonSig"),Label=ifelse(as.numeric(P.value)<0.05,CANCER_TYPE,""))

p2 <- ggplot(data=female_dat,aes(x=OA,y=YA, color = Sig)) +
  geom_point(position = position_surround(),aes(size=10))+theme_bw()+
  scale_color_manual(values = c("Sig"="orange","NonSig"="#C9C9CA"))+
  geom_abline(slope = 1, intercept = 0) +
  xlim(c(0, 1)) + ylim(c(0, 1))+ ylab("YA female ratio") + xlab("OA female ratio")+
  geom_text_repel(aes(label = Label,color="black"),
                  xlim=c(NA, 0),#控制标签在左边
                  direction = 'y',#nudge_x = 200,hjust=1,
                  fontface='italic',max.iter = 3e3,force=1,max.overlaps =200,
                  size = 5)  + theme_bw() 
p2 

p3 <- p1+p2
ggsave("./result/TCGA/0.1癌症突变_样本统计/调试13.2年龄分组间性别关联.pdf",p3,width = 10,height = 4)

