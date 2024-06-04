#以这个为主
library(dbplyr)
library(tidyverse)
rm(list=ls())
setwd("D:/课题/project10")
get_vcColors <- function (alpha = 1, websafe = FALSE, named = TRUE) {
  if (websafe) {
    col = c("#F44336", "#E91E63", "#9C27B0", 
            "#673AB7", "#3F51B5", "#2196F3", 
            "#03A9F4", "#00BCD4", "#009688", 
            "#4CAF50", "#8BC34A", "#CDDC39", 
            "#FFEB3B", "#FFC107", "#FF9800", 
            "#FF5722", "#795548", "#9E9E9E", 
            "#607D8B")
  }
  else {
    col = c(RColorBrewer::brewer.pal(11, name = "Paired"), 
            RColorBrewer::brewer.pal(11, name = "Spectral")[1:3], 
            "black", "violet", "royalblue", 
            "#7b7060", "#535c68")
    col = grDevices::adjustcolor(col = col, alpha.f = alpha)
  }
  if (named) {
    names(col) = names = c("Nonstop_Mutation", "Frame_Shift_Del", 
                           "IGR", "Missense_Mutation", "Silent", 
                           "Nonsense_Mutation", "RNA", "Splice_Site", 
                           "Intron", "Frame_Shift_Ins", "In_Frame_Del", 
                           "ITD", "In_Frame_Ins", "Translation_Start_Site", 
                           "Multi_Hit", "Amp", "Del", "Complex_Event", 
                           "pathway")
  }
  col
}
cancer <- c('ACC','BLCA','BRCA','CESC','CHOL','COAD','DLBC','ESCA','GBM','HNSC','KICH',
            'KIRC','KIRP','LAML','LGG','LIHC','LUAD','LUSC','MESO','OV','PAAD','PCPG','PRAD','READ','SARC','SKCM','STAD','TGCT','THCA','THYM','UCEC','UCS','UVM')


muts <- c()
for (ii in cancer) {
  print(ii)
  surv <- read_tsv(paste0("D:/数据/表达谱数据/",ii,"/TCGA-",ii,".survival.tsv")) %>% dplyr::select(sample,OS,OS.time)
  
  pheno <- read_tsv(paste0("D:/数据/表达谱数据/",ii,"/TCGA-",ii,".GDC_phenotype.tsv")) %>%
    select(submitter_id.samples,age_at_initial_pathologic_diagnosis,gender.demographic,race.demographic,sample_type.samples) %>%
    mutate(CANCER_TYPE=ii,sampleID=substr(submitter_id.samples,1,15)) %>% 
    left_join(.,surv,by=c("submitter_id.samples"="sample")) %>% select(!sampleID) %>%
    rename(Gender=gender.demographic,RACE=race.demographic,SAMPLE_TYPE=sample_type.samples)
  mut <- read_tsv(paste0("D:/数据/TCGA突变数据/TCGA_mut/",ii,"_mut.txt")) %>% 
    mutate(MUT_ID=paste(ii,1:nrow(.),sep = "_"),submitter_id.samples=substr(Tumor_Sample_Barcode,1,16),CANCER_TYPE=ii) %>%
    select(MUT_ID,CANCER_TYPE,Variant_Classification,Variant_Type,Tumor_Sample_Barcode,submitter_id.samples)
  nrow(mut)
  mut_pheno <- left_join(mut,pheno,by=c("CANCER_TYPE","submitter_id.samples")) %>%
    filter(!is.na(age_at_initial_pathologic_diagnosis))  %>% 
    mutate(AGE_GROUP=ifelse(as.numeric(age_at_initial_pathologic_diagnosis)<=60,"YA",ifelse(as.numeric(age_at_initial_pathologic_diagnosis)>60,"OA","None"))) %>% 
    unique()
  mut_pheno <-  mut_pheno %>% mutate(patient=substr(Tumor_Sample_Barcode,1,12))
  
  muts <- rbind(muts,mut_pheno)
  
}
muts_phenotype <- muts %>% unique()
write_tsv(muts_phenotype,"./result/TCGA/0.1癌症突变_样本统计/调试13.1泛癌的临床信息.tsv")

ref_df <- data.frame(CANCER_TYPE=c("PanCancer",'ACC','BLCA','BRCA','CESC','CHOL','COAD','DLBC','ESCA','GBM','HNSC','KICH','KIRC',
                                   'KIRP','LAML','LGG','LIHC','LUAD','LUSC','MESO','OV','PAAD','PCPG','PRAD','READ',
                                   'SARC','SKCM','STAD','TGCT','THCA','THYM','UCEC','UCS','UVM')
)
ref_df <- data.frame(CANCER_TYPE=c("PanCancer","GBM","LGG","UVM","HNSC","ACC","THCA","PCPG","DLBC","THYM","LUAD","LUSC","MESO","BRCA","ESCA","STAD","LIHC","CHOL","PAAD","READ","COAD","KICH","KIRC","KIRP","BLCA","PRAD","TGCT","CESC","OV","UCEC","UCS","SKCM","LAML","SARC"),
                         order=1:34)
#过滤出SNP
dat_1 <- data.table::fread("./result/TCGA/0.1癌症突变_样本统计/调试13.1泛癌的临床信息.tsv") %>% 
  mutate(Variant_GROUP=paste(Variant_Classification,Variant_Type,sep = "-"))
dat_2 <- data.table::fread("./result/TCGA/0.1癌症突变_样本统计/调试13.1泛癌的临床信息.tsv") %>% 
  mutate(Variant_GROUP=paste(Variant_Classification,Variant_Type,sep = "-"),CANCER_TYPE="PanCancer")

tot_mut <- dat_3 %>% group_by(CANCER_TYPE,AGE_GROUP) %>% summarise(nMut=n()) %>% ungroup() 
tot_sample <- dat_3 %>% select(CANCER_TYPE,Tumor_Sample_Barcode,AGE_GROUP) %>% unique() %>% group_by(CANCER_TYPE,AGE_GROUP) %>% summarise(nMut=n()) %>% ungroup() 

col = get_vcColors(alpha = 0.7, named = TRUE)
nonSilent = c("Frame_Shift_Del", "Frame_Shift_Ins",
              "Splice_Site", "Translation_Start_Site",
              "Nonsense_Mutation", "Nonstop_Mutation",
              "In_Frame_Del", "In_Frame_Ins", "Missense_Mutation")

GVariant_Classification_nMut <- dat_3 %>% 
  filter(Variant_Classification %in% nonSilent) %>%
  group_by(CANCER_TYPE,AGE_GROUP,Variant_Classification) %>% summarise(Variant_Classification_Mut=n()) %>% ungroup() %>%
  left_join(.,tot_mut,by=c("CANCER_TYPE","AGE_GROUP")) %>% mutate(ID=paste(CANCER_TYPE, AGE_GROUP,sep = "-")) %>% mutate(percent_Sample=Variant_Classification_Mut/nMut)
ref_df_vc <- left_join(ref_df,GVariant_Classification_nMut,by="CANCER_TYPE") %>% #%>% mutate(Order=paste(CANCER_TYPE, AGE_GROUP, Gender))
  arrange(order,desc(AGE_GROUP))
ref_df_vc$ID <- factor(ref_df_vc$ID,levels = rev(unique(ref_df_vc$ID)))


write_tsv(ref_df_vc ,"./result/TCGA/0.1癌症突变_样本统计/调试13.2年龄分组的突变类型统计.tsv")

p0 <- ggplot()+
  geom_bar(data=ref_df_vc,
           aes(x=percent_Sample,y=ID,
               fill=Variant_Classification),
           stat="identity",width = 0.8) +
  scale_fill_manual(values = col) +
  theme_classic()

ggsave("./result/TCGA/0.1癌症突变_样本统计/调试13.2年龄分组的突变类型统计.pdf",p0,height=8,width=5)
paste(names(col),as.character(col),sep = "'='") %>% paste(.,collapse = "','")


#各个癌症类型中的性别样本分布
Gender_nSample <- dat_3 %>% select(CANCER_TYPE,Tumor_Sample_Barcode,AGE_GROUP,Gender) %>% unique() %>% group_by(CANCER_TYPE,AGE_GROUP,Gender) %>% summarise(Gender_Sample=n()) %>% ungroup()%>% 
  mutate(ID=paste(CANCER_TYPE, AGE_GROUP,sep = "-")) %>% left_join(.,tot_sample,by=c("CANCER_TYPE","AGE_GROUP")) %>% mutate(percent_Sample=Gender_Sample/nMut)
ref_df_gender <- left_join(ref_df,Gender_nSample,by="CANCER_TYPE") %>% #%>% mutate(Order=paste(CANCER_TYPE, AGE_GROUP, Gender))
  arrange(order,desc(AGE_GROUP))
ref_df_gender$ID <- factor(ref_df_gender$ID,levels = rev(unique(ref_df_gender$ID)))
write_tsv(ref_df_gender,"./result/TCGA/0.1癌症突变_样本统计/调试13.3年龄分组的性别统计.tsv")
p1 <- ggplot()+
  geom_bar(data=ref_df_gender,
           aes(x=percent_Sample,y=ID,
               fill=Gender),
           stat="identity",width = 0.8) +
  scale_fill_manual(values = c("male"="#6B889A","female"="#F0CA80")) +
  theme_classic()

ggsave("./result/TCGA/0.1癌症突变_样本统计/调试13.3年龄分组的性别统计.pdf",p1,height=8,width=3)

SAMPLE_TYPE_nSample <- dat_3 %>% select(CANCER_TYPE,Tumor_Sample_Barcode,AGE_GROUP,SAMPLE_TYPE) %>% unique() %>% 
  mutate(SAMPLE_TYPE=ifelse(grepl("Primary",SAMPLE_TYPE),"Primary",ifelse(grepl("Metastatic",SAMPLE_TYPE),"Metastatic","Recurrent")))%>% 
  group_by(CANCER_TYPE,AGE_GROUP,SAMPLE_TYPE) %>% summarise(SAMPLE_TYPE_Sample=n()) %>% ungroup()%>% 
  mutate(ID=paste(CANCER_TYPE, AGE_GROUP,sep = "-")) %>% left_join(.,tot_sample,by=c("CANCER_TYPE","AGE_GROUP")) %>% mutate(percent_Sample=SAMPLE_TYPE_Sample/nMut)
ref_df_SAMPLE_TYPE <- left_join(ref_df,SAMPLE_TYPE_nSample,by="CANCER_TYPE") %>% #%>% mutate(Order=paste(CANCER_TYPE, AGE_GROUP, Gender))
  arrange(order,desc(AGE_GROUP))
ref_df_SAMPLE_TYPE$ID <- factor(ref_df_SAMPLE_TYPE$ID,levels = rev(unique(ref_df_SAMPLE_TYPE$ID)))
ref_df_SAMPLE_TYPE$SAMPLE_TYPE <- factor(ref_df_SAMPLE_TYPE$SAMPLE_TYPE,levels = c("Metastatic","Recurrent","Primary"))
write_tsv(ref_df_SAMPLE_TYPE,"./result/TCGA/0.1癌症突变_样本统计/调试13.4转移状态分组的样本统计.tsv")

p2 <- ggplot()+
  geom_bar(data=ref_df_SAMPLE_TYPE,
           aes(x=percent_Sample,y=ID,
               fill=SAMPLE_TYPE),
           stat="identity",width = 0.8) +
  scale_fill_manual(values = c("Primary"="#7FBADC","Metastatic"="#9832CA","Recurrent"="#FD7F4E")) +
  theme_classic()

ggsave("./result/TCGA/0.1癌症突变_样本统计/调试13.4转移状态分组的样本统计.pdf",p2,height=8,width=3)

tot_sample_1 <- dat_3 %>%  filter(RACE!="not reported") %>% select(CANCER_TYPE,Tumor_Sample_Barcode,AGE_GROUP) %>% unique() %>% group_by(CANCER_TYPE,AGE_GROUP) %>% summarise(nMut=n()) %>% ungroup() 
Race_nSample <- dat_3 %>% filter(RACE!="not reported") %>% select(CANCER_TYPE,Tumor_Sample_Barcode,AGE_GROUP,RACE) %>% unique() %>% group_by(CANCER_TYPE,AGE_GROUP,RACE) %>% summarise(RACE_Sample=n()) %>% ungroup()%>% 
  mutate(ID=paste(CANCER_TYPE, AGE_GROUP,sep = "-")) %>% left_join(.,tot_sample_1 ,by=c("CANCER_TYPE","AGE_GROUP")) %>% mutate(percent_Sample=RACE_Sample/nMut)# 
ref_df_race <- left_join(ref_df,Race_nSample,by="CANCER_TYPE") %>% #%>% mutate(Order=paste(CANCER_TYPE, AGE_GROUP, Race))
  arrange(order,desc(AGE_GROUP))
ref_df_race$ID <- factor(ref_df_race$ID,levels = rev(unique(ref_df_race$ID)))
write_tsv(ref_df_race,"./result/TCGA/0.1癌症突变_样本统计/调试13.5RACE状态分组的样本统计.tsv")

p3 <- ggplot()+
  geom_bar(data=ref_df_race,
           aes(x=percent_Sample,y=ID,
               fill=RACE),
           stat="identity",width = 0.8) +
  theme_classic()
p3
ggsave("./result/TCGA/0.1癌症突变_样本统计/调试13.5RACE状态分组的样本统计.pdf",p3,height=8,width=5)

#