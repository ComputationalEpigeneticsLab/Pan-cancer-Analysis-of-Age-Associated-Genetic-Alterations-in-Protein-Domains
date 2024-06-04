library(tidyverse)
library(dbplyr)
library(vroom)
library(VennDiagram)
library (venn)
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
  mutate(Region=paste(Transcript_ID, from, to, `hmm acc`, `hmm name`, clan,sep=";"),Score=E.value*(-log10(FDR+10e-16)),Score_2=-log10(FDR+10e-16)) %>%
  filter(as.numeric(nRegion)>=3,as.numeric(FDR)<0.05)
order_id <- openxlsx::read.xlsx("./result/TCGA/人体图-final3癌症归属系统及排序.xlsx")


raw %>% group_by(`hmm name`,AGE_GROUP) %>% summarise(n=n()) %>% arrange(desc(n))


group_by(raw,`hmm name`,AGE_GROUP) %>% summarise(n=n()) %>% arrange(desc(n))

domain_freq <- group_by(raw,`hmm name`,AGE_GROUP,CANCER_TYPE) %>% summarise(nCase=n()) %>% arrange(desc(nCase)) %>% ungroup()


color=c('LUAD'='#F2C1A6','LUSC'='#DAA682','MESO'='#F1A776','BLCA'='#D3AB95','KICH'='#f0cfb7','KIRC'='#ddac9e','KIRP'='#ef9767','PRAD'='#f1bda6','TGCT'='#BA5D2A',
        'COAD'='#bbddb4','ESCA'='#6eab59','READ'='#95c487','STAD'='#499436','CESC'='#afd9f0','OV'='#61c0e1','UCEC'='#92cde6','UCS'='#3bbcdf','ACC'='#f9dcdf','PCPG'='#ea8d8b','THCA'='#e25367',
        'CHOL'='#67c44d','LIHC'='#41960F','PAAD'='#4CC642',
        'DLBC'='#8d78b5','LAML'='#c8bddc','THYM'='#50479a','GBM'='#efe88a','LGG'='#e5d839','SARC'='#3e3a39','SKCM'='#8082a4','HNSC'='#9d7873','UVM'='#b9a08f','BRCA'='#bcbbbb')

top_domain <- raw %>% group_by(CANCER_TYPE,`hmm name`,AGE_GROUP) %>% summarise(n=n()) %>% arrange(desc(n)) %>% ungroup %>% group_by(CANCER_TYPE,AGE_GROUP) %>% top_n(n = 10,wt = "n") %>% filter(n>=15)

domain_freq_ya <- filter(domain_freq,AGE_GROUP=="YA")
domain_freq_oa <- filter(domain_freq,AGE_GROUP=="OA")

top_dat_ya <- filter(domain_freq_ya,`hmm name` %in% top_domain$`hmm name`)
top_dat_ya$`hmm name` <- factor(top_dat_ya$`hmm name`,levels = unique(top_domain$`hmm name`))
top_dat_ya$CANCER_TYPE <- factor(top_dat_ya$CANCER_TYPE,levels = order_id$CANCER_TYPE)
p1 <- ggplot(top_dat_ya,aes(x=`hmm name`,y=nCase,fill=CANCER_TYPE)) +
  geom_bar(position = "stack",stat = "identity")+ ylim(c(0,400))+
  scale_fill_manual(values = color[intersect(names(color),top_dat_ya$CANCER_TYPE)])+ theme_classic() +theme(axis.text.x = element_text(angle = 90, hjust = 1)) 


top_dat_oa <- filter(domain_freq_oa,`hmm name` %in% top_domain$`hmm name`)
top_dat_oa$`hmm name` <- factor(top_dat_oa$`hmm name`,levels = unique(top_domain$`hmm name`))
top_dat_oa$CANCER_TYPE <- factor(top_dat_oa$CANCER_TYPE,levels = order_id$CANCER_TYPE)
p2 <- ggplot(top_dat_oa,aes(x=`hmm name`,y=nCase,fill=CANCER_TYPE)) +
  geom_bar(position = "stack",stat = "identity")+ ylim(c(0,400))+
  scale_fill_manual(values = color[intersect(names(color),top_dat_oa$CANCER_TYPE)])+ theme_classic() +theme(axis.text.x = element_text(angle = 90, hjust = 1)) 

p3 <- p1/p2
p3
ggsave("./result/TCGA/1.突变优选/1.7优选domain的分析/1.柱状图-distribution_of_cancer_across_domain.pdf",p3,width = 8,height=6)

library(ggvenn)
p <- list(YA=raw$`hmm name`[which(raw$AGE_GROUP=="YA")],OA=raw$`hmm name`[which(raw$AGE_GROUP=="OA")])  %>%
  ggvenn(show_percentage = T,show_elements = F,label_sep = ",",
         digits = 1,stroke_color = "white",stroke_alpha=1,
         fill_color = c("#4C5AA7", "#EA5B5A"),
  )
p
ggsave("./result/TCGA/1.突变优选/1.7优选domain的分析/2.饼图-distribution_of_domain_across_age_group.pdf",p,width = 8,height=8)

top_domain_t2 <- raw %>% group_by(`hmm name`,AGE_GROUP) %>% summarise(n=n()) %>% arrange(desc(n)) %>% ungroup #%>% group_by(CANCER_TYPE,AGE_GROUP)

dom_in_ya <- setdiff(raw$`hmm name`[which(raw$AGE_GROUP=="YA")],raw$`hmm name`[which(raw$AGE_GROUP=="OA")])

raw_in_domYa <- filter(top_domain_t2,`hmm name` %in% dom_in_ya)
raw_in_domYa_top <- unique(raw_in_domYa$`hmm name`)[1:20] #YA中排序靠前的domain
dom_in_oa <- setdiff(raw$`hmm name`[which(raw$AGE_GROUP=="OA")],raw$`hmm name`[which(raw$AGE_GROUP=="YA")])
raw_in_domOa <- filter(top_domain_t2,`hmm name` %in% dom_in_oa)
raw_in_domOa_top <- unique(raw_in_domOa$`hmm name`)[1:20]#OA中排序靠前的domain

top_dat_ya <- filter(domain_freq,AGE_GROUP=="YA",`hmm name` %in% raw_in_domYa_top)
top_dat_ya$`hmm name` <- factor(top_dat_ya$`hmm name`,levels = raw_in_domYa_top)

top_dat_ya$CANCER_TYPE <- factor(top_dat_ya$CANCER_TYPE,levels = order_id$CANCER_TYPE)
p1 <- ggplot(top_dat_ya,aes(x=`hmm name`,y=nCase,fill=CANCER_TYPE)) +
  geom_bar(position = "stack",stat = "identity")+ #ylim(c(0,400))+
  scale_fill_manual(values = color[intersect(names(color),top_dat_ya$CANCER_TYPE)])+ theme_classic() +theme(axis.text.x = element_text(angle = 90, hjust = 1)) 


top_dat_oa <- filter(domain_freq,AGE_GROUP=="OA",`hmm name` %in% raw_in_domOa_top)
top_dat_oa$`hmm name` <- factor(top_dat_oa$`hmm name`,levels = raw_in_domOa_top)
top_dat_oa$CANCER_TYPE <- factor(top_dat_oa$CANCER_TYPE,levels = order_id$CANCER_TYPE)
p2 <- ggplot(top_dat_oa,aes(x=`hmm name`,y=nCase,fill=CANCER_TYPE)) +
  geom_bar(position = "stack",stat = "identity")+ #ylim(c(0,400))+
  scale_fill_manual(values = color[intersect(names(color),top_dat_oa$CANCER_TYPE)])+ theme_classic() +theme(axis.text.x = element_text(angle = 90, hjust = 1)) 

p3 <- p1/p2
ggsave("./result/TCGA/1.突变优选/1.7优选domain的分析/3.YA_OA独有的结构域top20.pdf",p3,width = 8,height=6)

library(RColorBrewer)
only_in_ya_raw <- filter(raw,`hmm name` %in% dom_in_ya) %>%
  dplyr::select(CANCER_TYPE,AGE_GROUP,`hmm name`) %>% 
  group_by(`hmm name`,AGE_GROUP) %>% 
  summarise(n=n())%>% ungroup %>% 
  arrange(desc(n)) %>% mutate(GROUP=n) %>% group_by(GROUP) %>% summarise(nGROUP=n()) %>%mutate(ratio=round(nGROUP/sum(.$nGROUP),3)*100,label=paste0(GROUP,"(",paste0(ratio,"%"),")"))
only_in_ya_raw$GROUP <- as.character(only_in_ya_raw$GROUP)
pdf("./result/TCGA/1.突变优选/1.7优选domain的分析/4.1YA中独有domain的癌症特异性.pdf",width = 5,height = 5)
pie(only_in_ya_raw$nGROUP, border="white", col=brewer.pal(n=10,name="Greys")[3:10] , label=only_in_ya_raw$label)
dev.off()
only_in_oa_raw <- filter(raw,`hmm name` %in% dom_in_oa) %>% 
  dplyr::select(CANCER_TYPE,AGE_GROUP,`hmm name`) %>% 
  group_by(`hmm name`,AGE_GROUP) %>% 
  summarise(n=n())%>% ungroup %>% 
  arrange(desc(n)) %>% mutate(GROUP=n) %>% group_by(GROUP) %>% summarise(nGROUP=n()) %>% mutate(ratio=round(nGROUP/sum(.$nGROUP),3)*100,label=paste0(GROUP,"(",paste0(ratio,"%"),")"))
only_in_oa_raw$GROUP <- as.character(only_in_oa_raw$GROUP)


library(ggpubr)
pdf("./result/TCGA/1.突变优选/1.7优选domain的分析/4.2OA中独有domain的癌症特异性.pdf",width = 5,height = 5)
pie(only_in_oa_raw$nGROUP, border="white", col=brewer.pal(n=10,name="Greys")[3:10], label=only_in_oa_raw$label)
dev.off()
