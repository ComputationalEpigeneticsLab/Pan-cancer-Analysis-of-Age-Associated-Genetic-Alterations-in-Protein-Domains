library(tidyverse)
library(dbplyr)
library(vroom)
rm(list=ls())
gc()
setwd("D:/课题/project10")
getRawRes <- function(x){
  # x=dom_ya
  # a <- filter(x,as.numeric(nRegion)>=3) %>% 
  #   arrange(as.numeric(P.value)) %>% 
  #   mutate(FDR=p.adjust(as.numeric(P.value),method = "BH"))
  a <- filter(x,as.numeric(nRegion)>=3) %>% 
    arrange(as.numeric(P.value)) %>% 
    group_by(CANCER_TYPE) %>%
    mutate(FDR=p.adjust(as.numeric(P.value),method = "BH"))
  
  b <- filter(x,as.numeric(nRegion)<3) %>%
    mutate(P.value=1,FDR=1)
  c <- rbind(a,b)
  # rerun(c)
}
getROI <- function(x){
  # x=mut_in_domain[1,]
  if(as.numeric(x["Protein_position"])>=as.numeric(x["from"])&as.numeric(x["Protein_position"])<=as.numeric(x["to"])){
    a <- x
  }else{
    x["type"] <- "nonDomain"
    a <- x
  }
  return(a)
}
domain <- read_tsv("D:/数据/domain区域/2023预测/2.结果整合及过滤/2.gencode.hg19_domain.tsv",col_names = T) %>%
  separate(`seq.id`,into = c("Transcript_ID","Gene_ID","AA_length"),sep = "[|]") %>% 
  rename(from=`envelope start`,to=`envelope end`) %>%
  mutate(Transcript_ID=substr(Transcript_ID,1,15),Gene_ID=substr(Gene_ID,1,15),LRegion=as.numeric(to)-as.numeric(from)+1) %>% 
  filter(`E-value`<0.0001,type=="Domain") %>% #
  select(Transcript_ID,AA_length,from,to,LRegion,`hmm acc`,`hmm name`,clan) %>%
  as.data.frame() 

cancer <- c( 'ACC',  'BLCA', 'BRCA','CESC','CHOL','COAD','DLBC','ESCA','GBM','HNSC','KICH',
             'KIRC','KIRP','LAML','LGG','LIHC','LUAD','LUSC','MESO','OV','PAAD','PCPG','PRAD','READ','SARC','SKCM','STAD','TGCT','THCA','THYM','UCEC','UCS','UVM')
cancer_color <- data.frame(colors=c("#F2C1A6","#DAA682","#F1A776","#D3AB95","#f0cfb7","#ddac9e","#ef9767","#f1bda6","#BA5D2A","#bbddb4","#6eab59","#95c487","#499436","#afd9f0","#61c0e1","#92cde6","#3bbcdf","#f9dcdf","#ea8d8b","#e25367","#67c44d","#41960F","#4CC642","#8d78b5","#c8bddc","#50479a","#efe88a","#e5d839","#3e3a39","#8082a4","#9d7873","#b9a08f","#bcbbbb"),
                           cancer=c("LUAD",    "LUSC",    "MESO",   "BLCA",   "KICH",   "KIRC",   "KIRP",   "PRAD",   "TGCT",   "COAD",   "ESCA",   "READ",   "STAD",   "CESC",   "OV",     "UCEC",    "UCS",    "ACC",   "PCPG",   "THCA",   "CHOL",   "LIHC",   "PAAD",   "DLBC",   "LAML",   "THYM",   "GBM",    "LGG",    "SARC",   "SKCM",   "HNSC",   "UVM",   "BRCA"))
cancer_color <- c('LUAD'='#F2C1A6','LUSC'='#DAA682','MESO'='#F1A776','BLCA'='#D3AB95','KICH'='#f0cfb7','KIRC'='#ddac9e','KIRP'='#ef9767','PRAD'='#f1bda6','TGCT'='#BA5D2A','COAD'='#bbddb4','ESCA'='#6eab59','READ'='#95c487','STAD'='#499436','CESC'='#afd9f0','OV'='#61c0e1','UCEC'='#92cde6','UCS'='#3bbcdf','ACC'='#f9dcdf','PCPG'='#ea8d8b','THCA'='#e25367','CHOL'='#67c44d',
                  'LIHC'='#41960F','PAAD'='#4CC642','DLBC'='#8d78b5','LAML'='#c8bddc','THYM'='#50479a','GBM'='#efe88a','LGG'='#e5d839','SARC'='#3e3a39','SKCM'='#8082a4','HNSC'='#9d7873','UVM'='#b9a08f','BRCA'='#bcbbbb')
muts <- c()
muts_sample <- c()
for (ii in cancer) {
  print(ii)
  pheno <- read_tsv(paste0("D:/数据/表达谱数据/",ii,"/TCGA-",ii,".GDC_phenotype.tsv")) %>%
    select(submitter_id.samples,age_at_initial_pathologic_diagnosis,gender.demographic) %>% mutate(CANCER_TYPE=ii)
  mut <- read_tsv(paste0("D:/数据/TCGA突变数据/TCGA_mut/",ii,"_mut.txt")) %>% 
    filter(Variant_Classification=="Missense_Mutation",Variant_Type=="SNP",Protein_position!=".") %>% 
    mutate(submitter_id.samples=substr(Tumor_Sample_Barcode,1,16),CANCER_TYPE=ii,Mut_ID=paste0("MUTATION_",1:nrow(.)))
  
  nrow(mut)
  mut_pheno <- left_join(mut,pheno,by=c("CANCER_TYPE","submitter_id.samples")) %>%
    filter(!is.na(age_at_initial_pathologic_diagnosis))  %>% 
    mutate(AGE_GROUP=ifelse(as.numeric(age_at_initial_pathologic_diagnosis)<=60,"YA",ifelse(as.numeric(age_at_initial_pathologic_diagnosis)>60,"OA","None"))) %>% 
    unique() %>% select(Mut_ID,CANCER_TYPE,Hugo_Symbol,Protein_position, Transcript_ID,Tumor_Sample_Barcode,age_at_initial_pathologic_diagnosis,AGE_GROUP)
  muts <- rbind(muts,mut_pheno)
  
  muts_sample <- rbind(muts_sample,unique(mut_pheno[,c("CANCER_TYPE","Tumor_Sample_Barcode","age_at_initial_pathologic_diagnosis","AGE_GROUP")]))
}

#domain
muts_dom <- left_join(muts,domain,by=c("Transcript_ID"))  %>% filter(!is.na(`hmm acc`))
muts_dom$type <- "domain"
muts_dom_2 <- apply(muts_dom, 1, getROI) %>%
  t() %>% 
  as.data.frame() %>% 
  filter(type=="domain") %>% unique() %>% 
  mutate(ID=paste(CANCER_TYPE,Hugo_Symbol, Transcript_ID, from, to, `hmm acc`, `hmm name`, clan, LRegion, AA_length,sep = ";"),
         from=as.numeric(from),to=as.numeric(to),TYPE="MUT") %>% arrange(Protein_position) %>% unique()

range(as.numeric(muts_dom_2$age_at_initial_pathologic_diagnosis))
muts_dom_2$AGE_Fragment <- cut(as.numeric(muts_dom_2$age_at_initial_pathologic_diagnosis),breaks = c(18, 30, 40, 50, 60, 70,80,100), 
                               labels = c("<=30y", "30-40y", "40-50y", "50-60y", "60-70y", "70-80y", ">80y"), include.lowest = TRUE)


age_frag_inCancer <- lapply(cancer, function(x) {
  data.frame(CANCER_TYPE=x,AGE_Fragment=c("<=30y", "30-40y", "40-50y", "50-60y", "60-70y", "70-80y", ">80y"))
}) %>% do.call("rbind",.)


data_dom <- muts_dom_2 %>% 
  group_by(CANCER_TYPE,AGE_Fragment) %>%
  summarise(nMut=n()) %>%
  mutate(Type="Domain") %>%
  left_join(age_frag_inCancer,.,by=c("CANCER_TYPE","AGE_Fragment")) %>%
  mutate(nMut=ifelse(is.na(nMut),0,nMut),Type=ifelse(is.na(Type),"Domain",Type))

data_dis_dom <- data_dom

library(RColorBrewer)
mycolors <- brewer.pal(12,"Set3")
barplot(rep(1,9),col=mycolors) 
data_dis_dom$AGE_Fragment <- factor(data_dis_dom$AGE_Fragment,levels = c("<=30y", "30-40y", "40-50y", "50-60y", "60-70y", "70-80y", ">80y"))
p1 <- ggplot(data_dis_dom,aes(x=AGE_Fragment,y=nMut,fill=CANCER_TYPE)) +
  geom_bar(position = "stack",stat = "identity") +
  theme_bw() + 
  theme(axis.text.x = element_text(angle=70, vjust=0.5))+
  scale_fill_manual(values = cancer_color)
p1
ggsave("1.distribution_of_samples.pdf.pdf",p1,width=8,height=5)

#统计各癌症AGEE_Fragment对应的样本数
data_sample <- select(muts_dom_2,CANCER_TYPE,AGE_Fragment,Tumor_Sample_Barcode) %>% unique() %>% group_by(CANCER_TYPE,AGE_Fragment) %>%
  summarise(nSample=n()) %>% filter(!is.na(AGE_Fragment))
p2 <- ggplot(data_sample,aes(x=AGE_Fragment,y=nSample,fill=CANCER_TYPE)) +
  geom_bar(position = "stack",stat = "identity") +
  theme_bw() + 
  theme(axis.text.x = element_text(angle=70, vjust=0.5))+
  scale_fill_manual(values = cancer_color)
ggsave("1.distribution_of_samples.pdf",p2,width=8,height=5)
