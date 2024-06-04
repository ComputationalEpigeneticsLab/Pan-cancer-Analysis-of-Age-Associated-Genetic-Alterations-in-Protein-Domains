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
domain <- read_tsv("D:/数据/domain区域/2023预测/2.结果整合及过滤/2.gencode.hg19_domain.tsv",col_names = T) %>%
  separate(`seq.id`,into = c("Transcript_ID","Gene_ID","AA_length"),sep = "[|]") %>% 
  dplyr::rename(from=`envelope start`,to=`envelope end`) %>%
  mutate(Transcript_ID=substr(Transcript_ID,1,15),Gene_ID=substr(Gene_ID,1,15),LRegion=as.numeric(to)-as.numeric(from)+1) %>% 
  filter(`E-value`<0.0001,type=="Domain") %>% #
  dplyr::select(Transcript_ID,AA_length,from,to,LRegion,`hmm acc`,`hmm name`,clan) %>%
  as.data.frame() 

cancer <- c('ACC','BLCA','BRCA','CESC','CHOL','COAD','DLBC','ESCA','GBM','HNSC','KICH',
            'KIRC','KIRP','LAML','LGG','LIHC','LUAD','LUSC','MESO','OV','PAAD','PCPG','PRAD','READ','SARC','SKCM','STAD','TGCT','THCA','THYM','UCEC','UCS','UVM')
nSampleInRegion <- c()
nSampleInSymbol <- c()
sampleInRegion <- c()
raw_res <- c()
tot_opt <- c()
for (ii in cancer) {
  print(ii)
  pheno <- read_tsv(paste0("D:/数据/表达谱数据/",ii,"/TCGA-",ii,".GDC_phenotype.tsv")) %>%
    select(submitter_id.samples,age_at_initial_pathologic_diagnosis) %>% mutate(CANCER_TYPE=ii)
  mut <- read_tsv(paste0("D:/数据/TCGA突变数据/TCGA_mut/",ii,"_mut.txt")) %>% 
    filter(Variant_Classification=="Missense_Mutation",Variant_Type=="SNP") %>% 
    mutate(submitter_id.samples=substr(Tumor_Sample_Barcode,1,16),CANCER_TYPE=ii)
  nrow(mut)
  mut_pheno <- left_join(mut,pheno,by=c("CANCER_TYPE","submitter_id.samples")) %>%
    filter(!is.na(age_at_initial_pathologic_diagnosis))  %>% 
    mutate(AGE_GROUP=ifelse(as.numeric(age_at_initial_pathologic_diagnosis)<=60,"YA",ifelse(as.numeric(age_at_initial_pathologic_diagnosis)>60,"OA","None"))) %>% 
    unique()
  select(mut_pheno,CANCER_TYPE,AGE_GROUP,submitter_id.samples) %>% unique() %>%group_by(CANCER_TYPE,AGE_GROUP) %>%
    summarise(nMisSample=n()) %>% ungroup
  select(mut_pheno,CANCER_TYPE,AGE_GROUP,submitter_id.samples) %>%group_by(CANCER_TYPE,AGE_GROUP) %>%
    summarise(nMisMuts=n()) %>% ungroup
  
  muts <- mut_pheno
  nTot <- muts %>% group_by(CANCER_TYPE,AGE_GROUP,Hugo_Symbol,Transcript_ID) %>% summarise(nTot=n()) %>% ungroup
  AGE_GROUP <- unique(muts$AGE_GROUP)
 
  for (jj in AGE_GROUP) {
    mut <- filter(muts,AGE_GROUP==jj) %>% filter(!is.na(Protein_position)&Protein_position!=".")
    mut_in_domain <- left_join(mut,domain,by=c("Transcript_ID")) %>% filter(!is.na(`hmm acc`))
    if(nrow(mut_in_domain)==0){
      next
    }
    mut_in_domain$type <- "domain"
    mut_in_domain_ROI <- apply(mut_in_domain, 1, getROI) %>%
      t() %>% 
      as.data.frame() %>% 
      filter(type=="domain") %>% 
      mutate(from=as.numeric(from),to=as.numeric(to))
    
    nSampleInRegion <- mut_in_domain_ROI %>%
      select(CANCER_TYPE,AGE_GROUP,Hugo_Symbol,Transcript_ID,from,to,`hmm acc`,`hmm name`,clan,Tumor_Sample_Barcode) %>% 
      unique() %>% 
      group_by(CANCER_TYPE,AGE_GROUP,Hugo_Symbol,Transcript_ID,from,to,`hmm acc`,`hmm name`,clan) %>% 
      summarise(nSampleInRegion=n()) %>% rbind(nSampleInRegion,.) %>% ungroup
    
    sampleInRegion_1 <- mut_in_domain_ROI %>%
      select(CANCER_TYPE,AGE_GROUP,Hugo_Symbol,Transcript_ID,from,to,`hmm acc`,`hmm name`,clan,Tumor_Sample_Barcode) %>% 
      unique()
    
    sampleInRegion <- rbind(sampleInRegion ,sampleInRegion_1)
    
    nSampleInSymbol <- mut_in_domain_ROI %>%
      select(CANCER_TYPE,AGE_GROUP,Hugo_Symbol,Tumor_Sample_Barcode) %>% 
      unique() %>% 
      group_by(CANCER_TYPE,AGE_GROUP,Hugo_Symbol) %>% 
      summarise(nSampleInSymbol=n()) %>% rbind(nSampleInSymbol,.) %>% ungroup
    
    nRegion <- mut_in_domain_ROI %>% unique() %>% 
      group_by(CANCER_TYPE,AGE_GROUP,Hugo_Symbol,Transcript_ID,from,to,`hmm acc`,`hmm name`,clan) %>% 
      summarise(nRegion=n()) %>% 
      filter(nRegion>=1) %>% mutate(from=as.numeric(from),to=as.numeric(to)) %>% ungroup
     nMut_inform <- left_join(nRegion,nTot,by=c("CANCER_TYPE","AGE_GROUP","Hugo_Symbol","Transcript_ID")) %>%
      as.data.frame()
    
    mut_in_domain_ROI_inform <- left_join(mut_in_domain_ROI,nMut_inform,by=intersect(colnames(mut_in_domain_ROI),colnames(nMut_inform)))
    
      dim(mut_in_domain_ROI_inform);length(unique(mut_in_domain_ROI_inform$Hugo_Symbol))
    
    nMut_inform__ROI <- left_join(nMut_inform,mut_in_domain_ROI,by=intersect(colnames(mut_in_domain_ROI),colnames(nMut_inform))) %>% 
      select(colnames(nMut_inform),LRegion,AA_length) %>% unique()
    
    nMut_inform__ROI <- nMut_inform__ROI %>% mutate(E.value=(as.numeric(nRegion)/as.numeric(nTot))/(as.numeric(LRegion)/as.numeric(AA_length)),
                                                    P.value=(1-pbinom(as.numeric(nRegion),as.numeric(nTot),as.numeric(LRegion)/as.numeric(AA_length))))
    raw_res <- rbind(raw_res,data.frame(nMut_inform__ROI,check.names = F))
    
    
    nMut_inform__ROI <- filter(nMut_inform__ROI,nRegion>=3)
    if(nrow(nMut_inform__ROI)==0){
      res_1 <- c(CANCER_TYPE=ii,rep("-",length(colnames(nMut_inform__ROI))))  %>% t()%>% as.data.frame()
      colnames(res_1) <- c("CANCER_TYPE","AGE_GROUP","Hugo_Symbol","Transcript_ID","from","to","hmm acc","hmm name","clan","nRegion","nTot","LRegion","AA_length","E.value","P.value","FDR")
      ncase=0
      ngene=0
      nmut=0
    
      tot_opt <- rbind(tot_opt,res_1)
      next
    }
    nMut_inform__ROI$FDR <- p.adjust(nMut_inform__ROI$P.value,method = "BH")
    
    nMut_inform__ROI_nRegion3 <- nMut_inform__ROI %>% filter(FDR<0.05) %>% select(!CANCER_TYPE) # E.value>2,
     if(nrow(nMut_inform__ROI_nRegion3)!=0){
      res_1 <- data.frame(CANCER_TYPE=ii,nMut_inform__ROI_nRegion3,check.names = F)
      ncase=nrow(res_1)
      ngene=length(unique(res_1$Hugo_Symbol))
      nmut=sum(res_1$nRegion)
    }else{
      res_1 <- c(CANCER_TYPE=ii,rep("-",length(colnames(nMut_inform__ROI_nRegion3))))  %>% t()%>% as.data.frame()
      colnames(res_1) <- c("CANCER_TYPE","AGE_GROUP","Hugo_Symbol","Transcript_ID","from","to","hmm acc","hmm name","clan","nRegion","nTot","LRegion","AA_length","E.value","P.value","FDR")
      ncase=0
      ngene=0
      nmut=0
    }
    tot_opt <- rbind(tot_opt,res_1)
    
  }
}
write_tsv(sampleInRegion,"./result/TCGA/1.突变优选/1.位于domain区域内的样本.tsv")
write_tsv(raw_res,"./result/TCGA/1.突变优选/1.TCGA_raw_domain_singleCancer.tsv")
write_tsv(tot_opt,"./result/TCGA/1.突变优选/3.TCGA_domain_singleCancer.tsv")

