library(tidyverse)
library(vroom)
library(ggrepel)
library(fgsea)
setwd("D:/课题/project10")
rm(list=ls())
least_phase_1 <- read_tsv("D:/数据/药物数据/CMAP/CLUE/repurposing_samples_20180907.txt",skip = 9)
file_lists <- list.files("./result/TCGA/11.domain突变与药效/Cmap_CLUE/v1/",pattern = "_1")

ii=file_lists[5]
dat_top20 <- c()
tot_dat <- c()
moa_list <- c()
gsea_cdk_list <- c()
for (ii in file_lists) {
  print(ii)
  dat_1 <- vroom(paste0("./result/TCGA/11.domain突变与药效/Cmap_CLUE/v1/",ii,"/query_result.gct"),skip = 2) %>%
    filter(id!="desc",pert_type=="trt_cp",) %>% #,fdr_q_nlog10>(-log10(0.05))
    arrange(norm_cs) 
  moa <- dat_1 %>% select(pert_iname,moa) %>% unique() %>% arrange(pert_iname,desc(moa)) %>% distinct(pert_iname,.keep_all = T)
  moa_list <- rbind(moa_list,moa)
  dupli_drug <- dat_1$pert_iname[which(duplicated(dat_1$pert_iname))]

  uni_cs_res <- dat_1 %>%
    mutate(norm_cs=as.numeric(norm_cs)) %>% 
    group_by(pert_iname) %>% summarise(
      hi_67 = quantile(as.numeric(norm_cs), probs = 0.67),
      lo_33 = quantile(as.numeric(norm_cs), probs = 0.33)
    ) %>%
    mutate(cell_norm_cs = ifelse(abs(hi_67) >= abs(lo_33), hi_67, lo_33)) %>% ungroup %>%
    filter(pert_iname %in% least_phase_1$pert_iname)%>% 
    left_join(.,moa,by="pert_iname")
  tot_dat <- rbind(tot_dat,data.frame(File=ii,uni_cs_res))
  
  top_20_neg <- uni_cs_res %>% 
    top_n(n=-20,wt = as.numeric(cell_norm_cs)) %>% mutate(type="negative")
  top_20_pos <- uni_cs_res %>%
    top_n(n=20,wt = as.numeric(cell_norm_cs)) %>% mutate(type="positive")
  
  dat_3 <- rbind(top_20_neg,top_20_pos) %>% arrange(cell_norm_cs) # %>% mutate(moa=ifelse(is.na(moa),"-666",moa))
  dat_top20 <- rbind(dat_top20 ,data.frame(File=ii,dat_3))
  
  dat_4 <- rbind(top_20_neg,data.frame(pert_iname="test",hi_67=0,lo_33=0,cell_norm_cs=0,moa="test",type="test"),top_20_pos) %>% arrange(cell_norm_cs)
  plot_dat_1 <- dat_4 %>% 
    mutate(virtual_col=cell_norm_cs) %>% 
    select(pert_iname,cell_norm_cs) %>% 
    column_to_rownames(var="pert_iname") #%>% 
  
  col<-c(colorRampPalette(c('#3D53A3','#FFFFFF','#E94B49'))(1000))
  bk <- unique(c(seq(-2.158828,2.069100,length=1000)))
  pdf(paste0("./result/TCGA/11.domain突变与药效/Cmap_CLUE/v1/",ii,"/",ii,"_cs_pheatmap.pdf"),width = 5,height = 6)
  pheatmap::pheatmap(as.matrix(plot_dat_1),#scale="column",
                     color = col,
                     breaks = bk,
                     cluster_rows = F,show_rownames = T,show_colnames = T,cluster_cols = F,
                     cellwidth=15,cellheight=8,
                     display_numbers =matrix(dat_4$moa, ncol=1)
  )
  dev.off()

  cdk_list  <- filter(moa,grepl("CDK",moa)) %>% mutate(Term="CDK inhibitor") %>% select(!moa) %>% rename(`CDK Inhibitor`=pert_iname)

  
  uni_cs_res_2 <- uni_cs_res %>% filter(cell_norm_cs<0) %>% arrange(desc(cell_norm_cs))
  uni_cs_res_v <- uni_cs_res_2$cell_norm_cs
  names(uni_cs_res_v) <- uni_cs_res_2$pert_iname
  
  gsea_cdk <-  fgsea(pathways =  cdk_list, 
        stats = uni_cs_res_v,
        minSize=5,
        maxSize=500,
        nperm=10000)
  gsea_cdk_list <- rbind(gsea_cdk_list, data.frame(File=ii,gsea_cdk))

}

write_tsv(tot_dat,"./result/TCGA/11.domain突变与药效/Cmap_CLUE/v1/3.0normalized_connectivity_score_least-phase-I.tsv")
write_tsv(dat_top20,"./result/TCGA/11.domain突变与药效/Cmap_CLUE/v1/3.1top20_normalized_connectivity_score_least-phase-I.tsv")

write_tsv(gsea_cdk_list,"./result/TCGA/11.domain突变与药效/Cmap_CLUE/v1/4.1CDK_Inhibitor_GSEA.tsv")

tot_dat <- read_tsv("./result/TCGA/11.domain突变与药效/Cmap_CLUE/v1/3.0normalized_connectivity_score_least-phase-I.tsv")
gsea_cdk_list <- read_tsv("./result/TCGA/11.domain突变与药效/Cmap_CLUE/v1/4.1CDK_Inhibitor_GSEA.tsv")

mean_cdk <- tot_dat %>% filter(grepl("CDK",moa)) %>% group_by(File) %>% summarise(mean_cdk=mean(cell_norm_cs))

plot_gsea_cdk <- left_join(gsea_cdk_list,mean_cdk,by=c("File"))  %>% arrange(desc(NES)) %>% mutate(Label = ifelse(padj<0.05,paste(File,"*",sep = "_"),paste(File,"#",sep = "_")))
plot_gsea_cdk$Label <- factor(plot_gsea_cdk$Label,levels = unique(plot_gsea_cdk$Label))
p1 <- ggplot(
  plot_gsea_cdk ,
  aes(Label, NES)) + #reorder(Label, NES)
  geom_col(aes(fill= mean_cdk)) +
  coord_flip() +
  scale_fill_gradient2(low = 'blue',mid = 'white',high = 'red') +
  labs(x="CDK Inhibitor", y="Normalized Enrichment Score",title="KEGG Gene Sets NES from GSEA") + theme_bw()
ggsave("./result/TCGA/11.domain突变与药效/Cmap_CLUE/v1/4.1CDK_Inhibitor_GSEA.pdf",p1,width = 5,height = 6)