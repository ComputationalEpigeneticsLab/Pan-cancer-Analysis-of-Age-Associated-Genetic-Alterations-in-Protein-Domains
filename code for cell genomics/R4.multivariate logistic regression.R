library(tidyverse)
library(dbplyr)
library(vroom)
setwd("D:/课题/project10")
rm(list=ls())
getROI <- function(x){
  if(is.na(as.numeric(x["from"]))){
    x["type"] <- "None"
    a <- x
  }else if(as.numeric(x["Protein_position"])>=as.numeric(x["from"])&as.numeric(x["Protein_position"])<=as.numeric(x["to"])){
    a <- x
  }else{
    x["type"] <- "nonDomain"
    a <- x
  }
  return(a)
}
getRawRes <- function(x){
  a <- filter(x,as.numeric(nRegion)>=3) %>% 
    arrange(as.numeric(P.value)) %>% 
    group_by(CANCER_TYPE,AGE_GROUP) %>%
    mutate(FDR=p.adjust(as.numeric(P.value),method = "BH")) %>% ungroup
  
  b <- filter(x,as.numeric(nRegion)<3) %>%
    mutate(P.value=1,FDR=1)
  c <- rbind(a,b)
}
raw <- read_tsv("./result/TCGA/1.突变优选/1.TCGA_raw_domain_singleCancer.tsv") %>% getRawRes %>% 
  mutate(Region=paste(Transcript_ID, from, to, `hmm acc`, `hmm name`, clan,sep=";")) %>%
  filter(as.numeric(nRegion)>=3,as.numeric(FDR)<0.05)

raw %>% group_by(CANCER_TYPE,Hugo_Symbol,   Transcript_ID, from , to, `hmm acc`, `hmm name`,   clan) %>% summarise(n=n()) %>% arrange(desc(n)) %>% ungroup
raw_stat <- raw %>% group_by(CANCER_TYPE,AGE_GROUP) %>% summarise(nCase=n()) %>% reshape2::dcast(.,CANCER_TYPE~AGE_GROUP,value.var = "nCase",fill = 0) %>% ungroup
tcga_sub <- read_tsv("D:/数据/TCGA/癌症亚型/TCGASubtype.20170308.tsv") %>% select(sampleID,Subtype_Selected) %>% rename(Subtype=Subtype_Selected)

mut_sample <- read_tsv("./result/TCGA/1.突变优选/1.位于domain区域内的样本.tsv") %>% mutate(patient=substr(Tumor_Sample_Barcode,1,12))
cancer <- unique(raw$CANCER_TYPE)
for (ii in cancer) {
  print(ii)
  pheno <- read_tsv(paste0("D:/数据/表达谱数据/",ii,"/TCGA-",ii,".GDC_phenotype.tsv")) %>%
    select(submitter_id.samples,age_at_initial_pathologic_diagnosis,gender.demographic,race.demographic,sample_type.samples) %>% mutate(CANCER_TYPE=ii,sampleID=substr(submitter_id.samples,1,15)) %>% 
    left_join(.,tcga_sub,by="sampleID") %>% mutate(Subtype=ifelse(grepl(".NA",Subtype),NA,Subtype)) %>% select(!sampleID) %>% rename(Gender=gender.demographic,RACE=race.demographic,SAMPLE_TYPE=sample_type.samples)
  mut <- read_tsv(paste0("D:/数据/TCGA突变数据/TCGA_mut/",ii,"_mut.txt")) %>% 
    filter(Variant_Classification=="Missense_Mutation",Variant_Type=="SNP") %>% 
    mutate(submitter_id.samples=substr(Tumor_Sample_Barcode,1,16),CANCER_TYPE=ii,MUT_ID=paste(ii,1:nrow(.),sep = "_"))%>%
    select(MUT_ID,CANCER_TYPE,Tumor_Sample_Barcode,submitter_id.samples)
  nrow(mut)
  mut_pheno <- left_join(mut,pheno,by=c("CANCER_TYPE","submitter_id.samples")) %>%
    filter(!is.na(age_at_initial_pathologic_diagnosis))  %>% 
    mutate(AGE_GROUP=ifelse(as.numeric(age_at_initial_pathologic_diagnosis)<=60,"YA",ifelse(as.numeric(age_at_initial_pathologic_diagnosis)>60,"OA","None"))) %>% 
    unique()
  muts <-  mut_pheno %>% mutate(patient=substr(Tumor_Sample_Barcode,1,12))
  muts_phenotype <- select(muts,CANCER_TYPE,AGE_GROUP,patient,colnames( pheno)) %>% unique()
  raw_sig <- filter(raw,nRegion>=3,FDR<0.05,CANCER_TYPE==ii) %>% ungroup()#%>% group_by(Hugo_Symbol) %>% summarise(n=n())
  regions <- raw_sig  %>% 
    select(Hugo_Symbol,   Transcript_ID, from , to, `hmm acc`, `hmm name`,   clan) %>%
    unique()
  res <- c()
  for (jj in 1:nrow(regions)) {
    if(jj%%100==0)print(jj)
    region <- regions[jj,]
    region_dat <- filter(raw_sig,Hugo_Symbol==as.character(region$Hugo_Symbol),Transcript_ID==as.character(region$Transcript_ID),
                         from==as.character(region$from),to==as.character(region$to),
                         `hmm acc`==as.character(region$`hmm acc`), `hmm name`==as.character(region$`hmm name`),clan==as.character(region$clan))# %>% select(!AGE_GROUP)
    m_sample <- left_join(region, mut_sample,by=intersect(colnames(region),colnames(mut_sample))) %>% filter(CANCER_TYPE==ii)
    sample_mut <- data.frame(CANCER_TYPE=ii,region,Tumor_Sample_Barcode=m_sample$Tumor_Sample_Barcode,patient=m_sample$patient)
    region_dat_2_wt <- filter(muts,!patient %in% m_sample$patient) %>% unique
    sample_wt <- data.frame(CANCER_TYPE=ii,region,Tumor_Sample_Barcode=region_dat_2_wt$Tumor_Sample_Barcode,patient=region_dat_2_wt$patient)
    
    muts_phenotype_1 <- filter(muts_phenotype,CANCER_TYPE==ii) %>%
      mutate(MUT_STATUS="None")
    muts_phenotype_1$MUT_STATUS[which(muts_phenotype_1$patient %in% m_sample$patient)] <- 1
    muts_phenotype_1$MUT_STATUS[which(muts_phenotype_1$patient %in% region_dat_2_wt$patient)] <- 0
    muts_phenotype_1$MUT_STATUS <- as.numeric(muts_phenotype_1$MUT_STATUS)
    
    muts_phenotype_1$Subtype  <-  as.numeric(factor(muts_phenotype_1$Subtype))-1
    muts_phenotype_1$SAMPLE_TYPE  <-  as.numeric(factor(muts_phenotype_1$SAMPLE_TYPE))-1
    muts_phenotype_1$Gender  <-  as.numeric(factor(muts_phenotype_1$Gender))-1
    muts_phenotype_1$RACE  <-  as.numeric(factor(muts_phenotype_1$RACE))-1
    
    muts_phenotype_1$AGE_STATUS <- ifelse(muts_phenotype_1$AGE_GROUP=="YA",1,0)
    
    muts_phenotype_1$AGE_STATUS <- factor(muts_phenotype_1$AGE_STATUS,levels = c(0,1),labels = c("Others","YA"))
    muts_phenotype_1$MUT_STATUS <- factor(muts_phenotype_1$MUT_STATUS,levels = c(0,1),labels = c("WT","MUT"))
    
    Subtype_stat <- filter(muts_phenotype_1,!is.na(Subtype))
    if(nrow(Subtype_stat)==0){
      fit1 <- glm(MUT_STATUS~AGE_STATUS+SAMPLE_TYPE+Gender+RACE,data = muts_phenotype_1,family = binomial())
      summary(fit1)
      coef(fit1)#取系数
      fit_res_1 <- broom::tidy(fit1,exponentiate = TRUE) %>% as.data.frame() %>% mutate(AGE_GROUP="YA") %>% select(AGE_GROUP,everything())
      fit_res_2 <- mutate( fit_res_1,CANCER_TYPE=ii,region,Sig=AGE_GROUP,coef=log(estimate)) %>%
        select(CANCER_TYPE,colnames(regions),everything())
    }else{
      fit1 <- glm(MUT_STATUS~AGE_STATUS+Subtype+SAMPLE_TYPE+Gender+RACE,data = muts_phenotype_1,family = binomial())
      summary(fit1)
      coef(fit1)#取系数
      fit_res_1 <- broom::tidy(fit1,exponentiate = TRUE) %>% as.data.frame() %>% mutate(AGE_GROUP="YA") %>% select(AGE_GROUP,everything())
      fit_res_2 <- mutate( fit_res_1,CANCER_TYPE=ii,region,Sig=AGE_GROUP,coef=log(estimate)) %>%
        select(CANCER_TYPE,colnames(regions),everything())
    }
    res <- rbind(res,fit_res_2)
  }
  res <- unique(res) %>% rename(Logistic.Pval=p.value)
  write_tsv(res,paste0("./result/TCGA/1.突变优选/1.5年龄偏倚/test_",ii,"_raw_年龄偏倚的优选事件_优化写法.tsv"))
  res_1 <- filter(res,grepl("AGE_STATUS",term))
  res_2 <- filter(res_1,Logistic.Pval<0.05)
  write_tsv(res_2,paste0("./result/TCGA/1.突变优选/1.5年龄偏倚/test_",ii,"_年龄偏倚的优选事件_优化写法.tsv"))
}





cancer <- c( 'ACC',  'BLCA', 'BRCA','CESC','CHOL','COAD','DLBC','ESCA','GBM','HNSC','KICH',
             'KIRC','KIRP','LAML','LGG','LIHC','LUAD','LUSC','MESO','OV','PAAD','PCPG','PRAD','READ','SARC','SKCM','STAD','TGCT','THCA','THYM','UCEC','UCS','UVM')

tot_dat <- c()
for (ii in cancer) {
  dat_1 <- read_tsv(paste0("./result/TCGA/1.突变优选/1.5年龄偏倚/test_",ii,"_raw_年龄偏倚的优选事件_优化写法.tsv")) %>%
    filter(grepl("AGE_STATUS",term)) %>% mutate(coef=log(estimate))
  if(nrow(dat_1)!=0){
    dat_1$Logistic.FDR <- p.adjust(dat_1$Logistic.Pval)
    
    tot_dat <- rbind(tot_dat,dat_1)
  }
}
tot_dat <- mutate(tot_dat,Region=paste0(Hugo_Symbol,":", from, "-",to, ":",`hmm name`))
tot_dat_sig <- filter(tot_dat,Logistic.Pval<0.05)


write_tsv(tot_dat,"./result/TCGA/1.突变优选/1.5年龄偏倚/tot_青年群体偏倚事件.tsv")

