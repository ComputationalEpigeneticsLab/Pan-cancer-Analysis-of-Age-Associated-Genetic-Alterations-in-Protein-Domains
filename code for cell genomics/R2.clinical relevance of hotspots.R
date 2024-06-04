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
raw_nCase <- raw %>% group_by(CANCER_TYPE,AGE_GROUP) %>% summarise(nTotCase=n()) %>% ungroup
tot_dat <- read_tsv("./result/TCGA/5.优选事件及互作的表达和预后效应/4.3.4优选domain的表达和预后效应_改-优化-改3/5.4tot_Prognosis_Expr-改2.tsv") %>% rename(`hmm name`=hmm.name,`hmm acc`=hmm.acc) %>%
  select(colnames(raw),Surv.Pval,Case_nor_wt_p.val, Case_nor_inR_p.val, Case_nor_outR_p.val, Case_wt_inR_p.val, Case_wt_outR_p.val, Case_inR_outR_p.val) %>% unique()
drug_r <- read_tsv("./result/TCGA/1.突变优选/1.4药物反应/Logistic_Response_WT_vs_Mut_raw.tsv") #%>% filter(Response_Logistic.Pval<0.1)
tot_dat <- full_join(tot_dat,drug_r,by=intersect(colnames(tot_dat),colnames(drug_r)))
surv_sig_n <- tot_dat %>% filter(Surv.Pval<0.05) %>% group_by(CANCER_TYPE,AGE_GROUP) %>% summarise(nSigCase=n()) %>% ungroup
raw_nCase_surv <- left_join(raw_nCase,surv_sig_n,by=c("CANCER_TYPE","AGE_GROUP")) %>% mutate(nSigCase=ifelse(is.na(nSigCase),0,nSigCase),Ratio=as.numeric(nSigCase)/as.numeric(nTotCase),AGE_GROUP=paste0(AGE_GROUP,"_Surv"))

raw_nCase_surv_plot_mtx <- raw_nCase_surv %>% reshape2::dcast(CANCER_TYPE~AGE_GROUP,value.var = "nSigCase") %>% column_to_rownames(var = "CANCER_TYPE") %>% as.matrix()
raw_nCase_surv_plot_mtx[which(is.na(raw_nCase_surv_plot_mtx))] <- 0


expr_sig_n <- tot_dat %>% filter(Case_wt_inR_p.val<0.05|Case_inR_outR_p.val<0.05) %>% group_by(CANCER_TYPE,AGE_GROUP) %>% summarise(nSigCase=n()) %>% ungroup
raw_nCase_expr <- left_join(raw_nCase,expr_sig_n,by=c("CANCER_TYPE","AGE_GROUP")) %>% mutate(nSigCase=ifelse(is.na(nSigCase),0,nSigCase),Ratio=as.numeric(nSigCase)/as.numeric(nTotCase),AGE_GROUP=paste0(AGE_GROUP,"_Expr"))

raw_nCase_expr_plot_mtx <- raw_nCase_expr %>% reshape2::dcast(CANCER_TYPE~AGE_GROUP,value.var = "nSigCase") %>% column_to_rownames(var = "CANCER_TYPE") %>% as.matrix()
raw_nCase_expr_plot_mtx[which(is.na(raw_nCase_expr_plot_mtx))] <- 0
drug_sig_n <- tot_dat %>% filter(Response_Logistic.Pval<0.05) %>%  group_by(CANCER_TYPE,AGE_GROUP) %>% summarise(nSigCase=n()) %>% ungroup
raw_nCase_drug <- left_join(raw_nCase,drug_sig_n,by=c("CANCER_TYPE","AGE_GROUP")) %>% mutate(nSigCase=ifelse(is.na(nSigCase),0,nSigCase),Ratio=as.numeric(nSigCase)/as.numeric(nTotCase),AGE_GROUP=paste0(AGE_GROUP,"_DrugResponse"))

raw_nCase_drug_plot_mtx <- raw_nCase_drug %>% reshape2::dcast(CANCER_TYPE~AGE_GROUP,value.var = "nSigCase") %>% column_to_rownames(var = "CANCER_TYPE") %>% as.matrix()
raw_nCase_drug_plot_mtx[which(is.na(raw_nCase_drug_plot_mtx))] <- 0
surv_expr_sig_n <- tot_dat %>% filter(Surv.Pval<0.05,Case_wt_inR_p.val<0.05|Case_inR_outR_p.val<0.05) %>% group_by(CANCER_TYPE,AGE_GROUP) %>% summarise(nSigCase=n()) %>% ungroup
raw_nCase_surv_expr <- left_join(raw_nCase,surv_expr_sig_n,by=c("CANCER_TYPE","AGE_GROUP")) %>% mutate(nSigCase=ifelse(is.na(nSigCase),0,nSigCase),Ratio=as.numeric(nSigCase)/as.numeric(nTotCase),AGE_GROUP=paste0(AGE_GROUP,"_Surv_Expr"))

raw_nCase_surv_expr_plot_mtx <- raw_nCase_surv_expr %>% reshape2::dcast(CANCER_TYPE~AGE_GROUP,value.var = "nSigCase") %>% column_to_rownames(var = "CANCER_TYPE") %>% as.matrix()
raw_nCase_surv_expr_plot_mtx[which(is.na(raw_nCase_surv_expr_plot_mtx))] <- 0
tot_dat %>% filter(Surv.Pval<0.05,Case_wt_inR_p.val<0.05|Case_inR_outR_p.val<0.05,Response_Logistic.Pval<0.05) #%>% group_by(CANCER_TYPE,AGE_GROUP) %>% summarise(nSigCase=n()) %>% ungroup

plot_dat <- cbind(raw_nCase_surv_plot_mtx,
                  raw_nCase_expr_plot_mtx ,
                  raw_nCase_drug_plot_mtx ,
                  raw_nCase_surv_expr_plot_mtx) %>% t()
pheatmap::pheatmap(plot_dat,cluster_rows = F,cluster_cols = F) #,scale = "row"


order_id <- openxlsx::read.xlsx("./result/TCGA/人体图-final3癌症归属系统及排序.xlsx")
plot_dat <- plot_dat[,order_id$CANCER_TYPE]

plot(density(plot_dat))
pdf("./result/TCGA/5.优选事件及互作的表达和预后效应/调试19.2YA和OA之间表达和生存差异的事件占比_pval05.pdf",width = 8,height = 6)
pheatmap::pheatmap(plot_dat,fontsize = 6.5,
                   # breaks = bk,
                   col = c(colorRampPalette(c("#FEF5DC","#DB5E7B","#712787"))(10000)),
                   show_colnames = T,show_rownames = T,#cutree_cols = 2,
                   cluster_rows = F,cluster_cols = F,
                   cellwidth=18,cellheight=8,
                   display_numbers =matrix(ifelse(plot_dat!=0, round(plot_dat,3), ""), nrow(plot_dat)),
                   gaps_row = c(2,4,6)
                   # annotation_col=name_col
)#,gaps_row = c(23,24)
dev.off()


