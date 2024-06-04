library(dbplyr)
library(tidyverse)
library(ggrepel)
library(ggBubbles)
library(ggpubr)
library(patchwork)
rm(list=ls())
setwd("D:/课题/project10")
ref_df <- data.frame(CANCER_TYPE=c("PanCancer","GBM","LGG","UVM","HNSC","ACC","THCA","PCPG","DLBC","THYM","LUAD","LUSC","MESO","BRCA","ESCA","STAD","LIHC","CHOL","PAAD","READ","COAD","KICH","KIRC","KIRP","BLCA","PRAD","TGCT","CESC","OV","UCEC","UCS","SKCM","LAML","SARC"),
                     order=1:34)
ref_df_race <- read_tsv("./result/TCGA/0.1癌症突变_样本统计/调试13.5RACE状态分组的样本统计.tsv")
ya_sample <- ref_df_race %>% group_by(CANCER_TYPE,AGE_GROUP) %>% summarise(Sample_N=sum(RACE_Sample)) %>%filter(AGE_GROUP=="YA")
oa_sample <- ref_df_race %>% group_by(CANCER_TYPE,AGE_GROUP) %>% summarise(Sample_N=sum(RACE_Sample)) %>%filter(AGE_GROUP=="OA")
ya_oa_sample <- left_join(ya_sample,oa_sample,by="CANCER_TYPE") %>% mutate(Sample_Ratio=Sample_N.x/Sample_N.y)

ref_df_race <- left_join(ref_df_race,ya_oa_sample[,c("CANCER_TYPE","Sample_Ratio")],by="CANCER_TYPE")

cancer <- unique(ref_df_race$CANCER_TYPE)
ii="LIHC"
race_res <-c()
for (ii in cancer) {
  dat_1 <- filter(ref_df_race,CANCER_TYPE==ii)
  races <- unique(dat_1$RACE)
  for (jj in races) {
    dat_2 <- dat_1 %>% mutate(RACE=ifelse(RACE==jj,RACE,"others")) %>% group_by(CANCER_TYPE,order,AGE_GROUP,RACE) %>% summarise(RACE_Sample=sum(RACE_Sample))
    ya_sele_race <- dat_2$RACE_Sample[which(dat_2$AGE_GROUP=="YA"&dat_2$RACE==jj)]
    if(length(ya_sele_race)==0){ya_sele_race <- 0}
    ya_other_race <- dat_2$RACE_Sample[which(dat_2$AGE_GROUP=="YA"&dat_2$RACE!=jj)]
    if(length(ya_other_race)==0){ya_other_race <- 0}
    
    oa_sele_race <- dat_2$RACE_Sample[which(dat_2$AGE_GROUP=="OA"&dat_2$RACE==jj)]
    if(length(oa_sele_race)==0){oa_sele_race <- 0}
    oa_other_race <- dat_2$RACE_Sample[which(dat_2$AGE_GROUP=="OA"&dat_2$RACE!=jj)]
    if(length(oa_other_race)==0){oa_other_race <- 0}
    
    dat_mtx <- matrix(c(ya_sele_race,ya_other_race,oa_sele_race,oa_other_race),nrow=2)
    colnames(dat_mtx) <- c("YA", "OA") 
    rownames(dat_mtx) <- c(jj, "others")
    dat_mtx
    p.val <- fisher.test(x = dat_mtx, alternative = "greater")$p.value
    odds_ratio <- fisher.test(x = dat_mtx, alternative = "greater")$estimate
    
    race_res <- rbind(race_res,data.frame(CANCER_TYPE=ii,RACE=jj,P.value=p.val,OR=odds_ratio))
  }
}

sele_race_OR_mtx <- race_res %>% reshape2::dcast(CANCER_TYPE~RACE,value.var = "OR") 
sele_race_OR_mtx <- left_join(ref_df,sele_race_OR_mtx,by="CANCER_TYPE") %>% column_to_rownames(var = "CANCER_TYPE") %>% select(!order) %>% as.matrix() %>% t()
v1 <- c(sele_race_OR_mtx)
v1 <- v1[which(!is.na(v1)&v1!=Inf)]
max_v1 <- max(v1)
sele_race_OR_mtx[which(sele_race_OR_mtx==Inf)] <- max_v1+5

sele_race_P_mtx <- race_res %>% reshape2::dcast(CANCER_TYPE~RACE,value.var = "P.value")  #%>% as.matrix()
sele_race_P_mtx <- left_join(ref_df,sele_race_P_mtx,by="CANCER_TYPE") %>% column_to_rownames(var = "CANCER_TYPE") %>% select(!order) %>% as.matrix() %>% t()

sele_race_P_mtx[which(is.na(sele_race_P_mtx))] <- 1

bk <- c(seq(0,max(v1),by=0.01),seq(0.01,5,by=0.01))
bk <-  c(seq(0,1,length.out=4999),seq(1.000001,max_v1,length.out=4999),seq(max_v1+0.00001,max_v1+5,length.out=2))
pdf("./result/TCGA/0.1癌症突变_样本统计/调试13.2年龄分组间RACE关联.pdf",width = 8,height = 4)
pheatmap::pheatmap(sele_race_OR_mtx,
         scale = "none",
         c(colorRampPalette(c("navy","white","firebrick3"))(10000)),
         breaks = bk,
         fontsize_row=10,fontsize_col=10,cellwidth=15,cellheight=15,cluster_rows =FALSE,cluster_cols=FALSE,
         display_numbers =matrix(ifelse(sele_race_P_mtx<0.05, "**", ""), nrow(sele_race_P_mtx)),
         number_color="black",
         angle_col = 45)
dev.off()

p1 <- ggplot(data=sele_race_dat,aes(x=OA,y=YA, color = Sig)) +
  geom_point(position = position_surround(),aes(size=10))+theme_bw()+
  scale_color_manual(values = c("Sig"="orange","NonSig"="#C9C9CA"))+
  geom_abline(slope = 1, intercept = 0) +
  xlim(c(0, 1)) + ylim(c(0, 1)) + ylab("YA sele_race ratio") + xlab("OA sele_race ratio")+
  geom_text_repel(aes(label = Label,color="black"),
                  xlim=c(NA, 0),#控制标签在左边
                  direction = 'y',#nudge_x = 200,hjust=1,
                  fontface='italic',max.iter = 3e3,force=1,max.overlaps =200,
                  size = 5)  + theme_bw() 
p1 



other_race_dat <- ref_df_race_fisher[which(ref_df_race_fisher$RACE=="other_race"),] %>% reshape2::dcast(CANCER_TYPE+P.value+OR+Sample_Ratio~AGE_GROUP,value.var = "percent_Sample") %>% 
  mutate(Sig=ifelse(as.numeric(P.value)<0.05,"Sig","NonSig"),Label=ifelse(as.numeric(P.value)<0.05,CANCER_TYPE,""))

p2 <- ggplot(data=other_race_dat,aes(x=OA,y=YA, color = Sig)) +
  geom_point(position = position_surround(),aes(size=10))+theme_bw()+
  scale_color_manual(values = c("Sig"="orange","NonSig"="#C9C9CA"))+
  geom_abline(slope = 1, intercept = 0) +
  xlim(c(0, 1)) + ylim(c(0, 1))+ ylab("YA other_race ratio") + xlab("OA other_race ratio")+
  geom_text_repel(aes(label = Label,color="black"),
                  xlim=c(NA, 0),#控制标签在左边
                  direction = 'y',#nudge_x = 200,hjust=1,
                  fontface='italic',max.iter = 3e3,force=1,max.overlaps =200,
                  size = 5)  + theme_bw() 
p2 

p3 <- p1+p2
ggsave("./result/TCGA/0.1癌症突变_样本统计/调试13.2年龄分组间性别关联.pdf",p3,width = 10,height = 4)

