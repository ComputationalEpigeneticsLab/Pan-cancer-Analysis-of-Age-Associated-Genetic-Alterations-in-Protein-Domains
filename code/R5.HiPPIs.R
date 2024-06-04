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
}

raw <- read_tsv("./result/TCGA/1.突变优选/1.TCGA_raw_domain_singleCancer.tsv") %>% getRawRes %>% 
  mutate(Region=paste(Transcript_ID, from, to, `hmm acc`, `hmm name`, clan,sep=";")) %>%
  filter(as.numeric(nRegion)>=3,as.numeric(FDR)<0.05)
both_gene <- read_tsv("./result/TCGA/3.优选基因的交集/3.1各癌症中交集的优选基因.tsv") %>% separate_rows(Hugo_Symbol,sep = ";")

domain_int <- read_tsv("D:/数据/domain区域/2023预测/4.domain参与的互作/6.gencode.hg19_domain_INSIDER.tsv") %>% 
  filter(P1_IntFace2Domain1!="nonOverlap"|P2_IntFace2Domain2!="nonOverlap") %>% 
  mutate(Region=paste(Transcript_ID, from, to, hmm.acc, hmm.name, clan,sep=";")) %>% select(Region,INSIDER_ID,Trans_Status)
int <- read_tsv("D:/数据/互作接口数据/INSIDER/H_sapiens_interfacesHQ_anno_PDB.txt") %>%
  mutate(P1_Trans=substr(P1_Trans,1,15),P2_Trans=substr(P2_Trans,1,15),
  ) %>%
  filter(!is.na(P1_Trans)&!is.na(P2_Trans)) %>% select(INSIDER_ID,PDB_ID,P1,P2,P1_Trans,P1_Gene,P2_Trans,P2_Gene) %>% unique()
gtf <- read_tsv("D:/数据/基因组/GENCODE/gencode.v19.gtf.basic.txt") %>% select(Hugo_Symbol,Transcript_ID,Gene_ID) %>% unique()
colnames(gtf) <- c("gtf_Hugo_Symbol","gtf_Transcript_ID","gtf_Gene_ID")

int_gtf <- left_join(int ,gtf,by=c("P1_Trans"="gtf_Transcript_ID")) %>% rename(P1_gtf_Hugo_Symbol=gtf_Hugo_Symbol,P1_gtf_Gene_ID=gtf_Gene_ID) %>%
  left_join(.,gtf,by=c("P2_Trans"="gtf_Transcript_ID")) %>% rename(P2_gtf_Hugo_Symbol=gtf_Hugo_Symbol,P2_gtf_Gene_ID=gtf_Gene_ID) %>%  
  filter(P1_gtf_Hugo_Symbol!=P2_gtf_Hugo_Symbol) %>%
  apply(.,1,function(x){
    if(as.character(x['P1_gtf_Hugo_Symbol'])>as.character(x['P2_gtf_Hugo_Symbol'])){
      x['P1'] <- as.character(x['P1_gtf_Hugo_Symbol'])
      x['P2'] <- as.character(x['P2_gtf_Hugo_Symbol'])
    }else{
      x['P1'] <- as.character(x['P2_gtf_Hugo_Symbol'])
      x['P2'] <- as.character(x['P1_gtf_Hugo_Symbol'])
    }
    return(x)
  }) %>% t() %>%as.data.frame() %>% mutate(IntAct=paste(P1,P2,sep = ";"))
tot_int <- unique(int_gtf$IntAct)

raw_domain_int <- left_join(raw,domain_int,by="Region") %>%
  left_join(.,int,by="INSIDER_ID")
raw_domain_int_gtf <- left_join(raw_domain_int ,gtf,by=c("P1_Trans"="gtf_Transcript_ID")) %>% rename(P1_gtf_Hugo_Symbol=gtf_Hugo_Symbol,P1_gtf_Gene_ID=gtf_Gene_ID) %>%
  left_join(.,gtf,by=c("P2_Trans"="gtf_Transcript_ID")) %>% rename(P2_gtf_Hugo_Symbol=gtf_Hugo_Symbol,P2_gtf_Gene_ID=gtf_Gene_ID)


raw_domain_int_gtf_symbol <- select(raw_domain_int_gtf,CANCER_TYPE,AGE_GROUP,Hugo_Symbol,Transcript_ID,from,to,`hmm acc`,`hmm name`,clan,nRegion,nTot,LRegion,AA_length,E.value,P.value,FDR,
                                    Region,INSIDER_ID,Trans_Status,P1_gtf_Hugo_Symbol,P1_gtf_Gene_ID,P2_gtf_Hugo_Symbol,P2_gtf_Gene_ID) %>%
  filter(P1_gtf_Hugo_Symbol!=P2_gtf_Hugo_Symbol) %>% unique()%>%  apply(.,1,function(x){
                                      if(as.character(x['P1_gtf_Hugo_Symbol'])>as.character(x['P2_gtf_Hugo_Symbol'])){
                                        x['P1'] <- as.character(x['P1_gtf_Hugo_Symbol'])
                                        x['P2'] <- as.character(x['P2_gtf_Hugo_Symbol'])
                                      }else{
                                        x['P1'] <- as.character(x['P2_gtf_Hugo_Symbol'])
                                        x['P2'] <- as.character(x['P1_gtf_Hugo_Symbol'])
                                      }
                                      return(x)
                                    }) %>% t() %>%as.data.frame() %>% mutate(IntAct=paste(P1,P2,sep = ";"))
write_tsv(unique(raw_domain_int_gtf_symbol),"./result/TCGA/4.优选事件参与的互作/4.0优选事件参与的互作.tsv")