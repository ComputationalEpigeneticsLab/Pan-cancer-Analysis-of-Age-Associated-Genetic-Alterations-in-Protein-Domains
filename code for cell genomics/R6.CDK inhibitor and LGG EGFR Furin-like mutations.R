library(tidyverse)
library(survival)
library(survminer)
library(patchwork)
setwd("D:/课题/project10")
rm(list=ls())
mut_sample <- read_tsv("./result/TCGA/1.突变优选/1.位于domain区域内的样本.tsv") %>%
  mutate(ID=paste(CANCER_TYPE,AGE_GROUP,Hugo_Symbol,Transcript_ID, from, to, `hmm acc`, `hmm name`, clan,sep = ";")) %>% select(ID,Tumor_Sample_Barcode) %>% unique()
input_dat <- read_tsv("./result/TCGA/11.domain突变与药效/Cmap_CLUE/v1/input_data.tsv") %>% select(ID) %>% 
  unique %>%
  separate(ID,into = c("CANCER_TYPE","AGE_GROUP","Hugo_Symbol","Transcript_ID","from","to","hmm.acc","hmm.name","clan"),sep = ";", remove = F) %>%
  mutate(File=paste(CANCER_TYPE,Hugo_Symbol,1,sep = "_")) %>% select(ID,File,CANCER_TYPE,AGE_GROUP)

input_dat_sample <- left_join(input_dat,mut_sample,by="ID") %>% unique()
input_dat_sample %>% filter(File %in% input_dat$File)  %>% group_by(ID,File,CANCER_TYPE,AGE_GROUP) %>% summarise(n=n()) %>% arrange(desc(n))
tot_cs <- read_tsv("./result/TCGA/11.domain突变与药效/Cmap_CLUE/v1/3.0normalized_connectivity_score_least-phase-I.tsv") %>%
  filter(grepl("CDK",moa)) %>% mutate(Term="CDK inhibitor") %>% select(!moa) %>% rename(`CDK Inhibitor`=pert_iname)

file_lists <- list.files("./result/TCGA/11.domain突变与药效/Cmap_CLUE/v1/",pattern = "_1$")

raw_dat <- c()
for (ii in file_lists) {
  print(ii)
  dat_1 <- vroom::vroom(paste0("./result/TCGA/11.domain突变与药效/Cmap_CLUE/v1/",ii,"/query_result.gct"),skip = 2) %>%
    filter(id!="desc",pert_type=="trt_cp",) %>% #,fdr_q_nlog10>(-log10(0.05))
    arrange(norm_cs) %>%
    select(pert_id,pert_iname,moa,target_name) %>% unique()
  raw_dat <- rbind(raw_dat,data.frame(File=ii,dat_1))
}

raw_dat_1 <- filter(raw_dat,target_name!="-666")
tot_cs_target <- left_join(tot_cs,raw_dat,by=c("File","CDK Inhibitor"="pert_iname")) %>% select(File,`CDK Inhibitor`,target_name,cell_norm_cs) %>% unique() %>% separate_rows(target_name,sep = "[|]")

case_sele <- c("COAD_TP53_1","HNSC_TP53_1","BRCA_TP53_1","LGG_EGFR_1","LUAD_TP53_1","UCEC_CDK2_1","LIHC_TP53_1","SARC_TP53_1","UCEC_SMC1A_1") #这个列表是通过图6E中Reverser的事件选择（即gsea显著富集，且mean drug connectivity为负数的）
neg_cs_target_sample <- filter(tot_cs_target, cell_norm_cs<0) %>% left_join(.,input_dat_sample,by="File") %>% mutate(submitter_id.samples=substr(Tumor_Sample_Barcode,1,16))
cancers <- unique(neg_cs_target_sample$CANCER_TYPE)
ya_can_pheno <- c()
for (cancer in cancers) {
  print(cancer)
  pheno <- read_tsv(paste0("D:/数据/表达谱数据/",cancer,"/TCGA-",cancer,".GDC_phenotype.tsv")) %>%
    dplyr::select(submitter_id.samples,age_at_initial_pathologic_diagnosis) %>% 
    mutate(AGE_GROUP=ifelse(age_at_initial_pathologic_diagnosis<=60,"YA",ifelse(age_at_initial_pathologic_diagnosis>60,"OA","NA")))
  surv <- read_tsv(paste0("D:/数据/表达谱数据/",cancer,"/TCGA-",cancer,".survival.tsv")) %>% dplyr::select(sample,OS,OS.time)
  pheno <- left_join(pheno,surv,by=c("submitter_id.samples"="sample"))
  can_exprs <- read.table(paste0("D:/数据/表达谱数据/",cancer,"/cancer_exp.tsv"),header = T,stringsAsFactors = F,sep = '\t',row.names = 1,check.names = F) %>%
    rownames_to_column(var = "Hugo_Symbol") %>% 
    reshape2::melt() %>% 
    left_join(.,pheno,by=c("variable"="submitter_id.samples")) %>%
    dplyr::select(!age_at_initial_pathologic_diagnosis) %>%
    mutate(Type="Cancer")
  
  ya_can_exprs <- filter(can_exprs,AGE_GROUP=="YA",Hugo_Symbol %in% unique(neg_cs_target_sample$target_name))
  ya_can_pheno <- rbind(ya_can_pheno,data.frame(CANCER_TYPE=cancer,ya_can_exprs))
}




neg_cs_target_sample_id <- neg_cs_target_sample %>%
  mutate(ID2=paste(File,`CDK Inhibitor`,target_name,sep = ";")) %>% filter(target_name != "-666") 
microglia_df <- read_tsv("./result/TCGA/12.domain突变与免疫/LGG_Microglia_Signature.tsv")
timer <- read_csv("D:/数据/TIMER/infiltration_estimation_for_tcga.csv") %>% #reshape2::melt() %>% 
  dplyr::rename(submitter_id.samples=cell_type) %>% select(submitter_id.samples,colnames(.)[grepl("_XCELL",colnames(.))])


cases <- unique(neg_cs_target_sample_id$ID2)
tot_res <- c()
infi_res <- c()
for (id in cases) {
  print(id)
  dat_1 <- filter(neg_cs_target_sample_id , ID2==id)
  can_1 <- unique(dat_1$CANCER_TYPE)
  symbol_1 <- unique(dat_1$target_name)
  
  ya_can_pheno_1 <- filter(ya_can_pheno,CANCER_TYPE==can_1,Hugo_Symbol==symbol_1 ) %>% 
    mutate(MUT_STATUS=ifelse(variable %in% dat_1$submitter_id.samples,"InRegion","Others"),GROUP=ifelse(as.numeric(value)>median(as.numeric(value)),"High","Low"))
  select(ya_can_pheno_1,MUT_STATUS,GROUP)
  
  fisher_dat <- table(ya_can_pheno_1$GROUP,ya_can_pheno_1$MUT_STATUS) %>% as.matrix()
  f_t <- fisher.test(fisher_dat ,alternative = "greater")
  fhsier_p <- f_t$p.value
  fisher_or <- f_t$estimate
  
  deg_p <- wilcox.test(ya_can_pheno_1$value[which(ya_can_pheno_1$MUT_STATUS=="InRegion")],ya_can_pheno_1$value[which(ya_can_pheno_1$MUT_STATUS=="Others")])$p.value
  meanHi <- mean(ya_can_pheno_1$value[which(ya_can_pheno_1$MUT_STATUS=="InRegion")])
  meanOthers <- mean(ya_can_pheno_1$value[which(ya_can_pheno_1$MUT_STATUS=="Others")])
  if(meanHi>meanOthers){
    deg_type="highHi"
  }else if(meanHi<meanOthers){
    deg_type="lowHi"
  }else{
    deg_type="equal"
  }
  
  sdf  <- survdiff(Surv(time = OS.time,event = OS ) ~ GROUP, data =  ya_can_pheno_1)# %>% try(silent = TRUE) %>% error_proc_2
  p.val = 1 - pchisq(sdf$chisq, length(sdf$n) - 1)
  lo_order <- grep("GROUP=Low",names(sdf$n))
  hi_order <- grep("GROUP=High",names(sdf$n))
  HR = (sdf$obs[hi_order]/sdf$exp[hi_order])/(sdf$obs[lo_order]/sdf$exp[lo_order])#相较于out，In的生存状态
  if(length(HR)==0){HR=NAN}
  fit_single_cox <- survfit(Surv(time = OS.time,event = OS ) ~ GROUP, data =  ya_can_pheno_1)
  g1 <- ggsurvplot(fit_single_cox,
                   pval = T,risk.table = F,conf.int = T,#title=ii,#paste0("p = ",data$cox_p_values)
                   xlab = "Time(days)",
                   ncensor.plot = F,linetype = "strata", # Change line type by groups
                   #surv.median.line = "hv", # Specify median survival
                   legend.title = "Group",
                   font.x = c(14,"bold.italic"),
                   font.y = c(14,"bold.italic"),
                   font.legend = c(12,"bold.italic"),
                   font.tickslab = c(12,"bold.italic"),
                   legend.labs = c("High","Low"),
                   ggtheme = theme_bw(), # Change ggplot2 theme
                   palette = c("High"="#fe6a6a","Low"="#3577AD")
  )
  if(grepl("LGG_EGFR_1",id)&p.val<0.05){
    ggsave(paste0("./result/TCGA/11.domain突变与药效/Cmap_CLUE/v1/LGG_EGFR_1/药物靶基因分高低表达的生存差异/Surv.",gsub(";","__",id),".pdf"),g1[[1]],width = 5,height = 5)
  }
  
  if(grepl("LGG_",id)){
    ya_can_pheno_2 <- left_join(ya_can_pheno_1,microglia_df,by=c("variable"="submitter_id.samples")) %>% mutate(Microglia_Polarization=Microglia_Signature_M2-c(0.65*Microglia_Signature_M1))
    mt_deg_p <- wilcox.test(Microglia_Signature ~ GROUP, data = ya_can_pheno_2)$p.value
    mt_deg_p_m1 <- wilcox.test(Microglia_Signature_M1 ~ GROUP, data = ya_can_pheno_2)$p.value
    mt_deg_p_m2 <- wilcox.test(Microglia_Signature_M2 ~ GROUP, data = ya_can_pheno_2)$p.value
    mt_deg_p_polar <- wilcox.test(Microglia_Polarization ~ GROUP, data = ya_can_pheno_2)$p.value
    
    mean_group <- ya_can_pheno_2 %>% group_by(GROUP) %>% 
      summarise(mean_Microglia_Signature=mean(Microglia_Signature),mean_Microglia_Signature_M1=mean(Microglia_Signature_M1),mean_Microglia_Signature_M2=mean(Microglia_Signature_M2),
                mean_Microglia_Polarization=mean(Microglia_Polarization)) %>% 
      column_to_rownames(var = "GROUP")%>% 
      t() %>% as.data.frame() %>% mutate(TYPE=ifelse(High>Low,"High","Low"))
    mt_group_res <- data.frame(Microglia_Signature_GROUP.Pval=mt_deg_p,Microglia_Signature_GROUP_Hi_mean=mean_group["mean_Microglia_Signature","High"],Microglia_Signature_GROUP_Lo_mean=mean_group["mean_Microglia_Signature","Low"],
                               Microglia_Signature_M1_GROUP.Pval=mt_deg_p_m1,Microglia_Signature_M1_GROUP_Hi_mean=mean_group["mean_Microglia_Signature_M1","High"],Microglia_Signature_M1_GROUP_Lo_mean=mean_group["mean_Microglia_Signature_M1","Low"], #M1
                               Microglia_Signature_M2_GROUP.Pval=mt_deg_p_m2,Microglia_Signature_M2_GROUP_Hi_mean=mean_group["mean_Microglia_Signature_M2","High"],Microglia_Signature_M2_GROUP_Lo_mean=mean_group["mean_Microglia_Signature_M2","Low"], #M2
                               Microglia_Polarization_GROUP.Pval=mt_deg_p_polar,Microglia_Polarization_GROUP_Hi_mean=mean_group["mean_Microglia_Polarization","High"],Microglia_Polarization_GROUP_Lo_mean=mean_group["mean_Microglia_Polarization","Low"]  #极化
                               )
    

  }else{
    mt_group_res <- data.frame(rep(NA,12)) %>% t() %>% as.data.frame(row.names = 1)
    colnames(mt_group_res) <- c('Microglia_Signature_GROUP.Pval','Microglia_Signature_GROUP_Hi_mean','Microglia_Signature_GROUP_Lo_mean',
                                'Microglia_Signature_M1_GROUP.Pval','Microglia_Signature_M1_GROUP_Hi_mean','Microglia_Signature_M1_GROUP_Lo_mean',
                                'Microglia_Signature_M2_GROUP.Pval','Microglia_Signature_M2_GROUP_Hi_mean','Microglia_Signature_M2_GROUP_Lo_mean',
                                'Microglia_Polarization_GROUP.Pval','Microglia_Polarization_GROUP_Hi_mean','Microglia_Polarization_GROUP_Lo_mean')
  }
  if(grepl("LGG_",id)){
    ya_can_pheno_2 <- left_join(ya_can_pheno_1,microglia_df,by=c("variable"="submitter_id.samples"))%>% mutate(Microglia_Polarization=Microglia_Signature_M2-c(0.65*Microglia_Signature_M1))
    mt_deg_p <- wilcox.test(Microglia_Signature ~ MUT_STATUS, data = ya_can_pheno_2)$p.value
    mt_deg_p_m1 <- wilcox.test(Microglia_Signature_M1 ~ MUT_STATUS, data = ya_can_pheno_2)$p.value
    mt_deg_p_m2 <- wilcox.test(Microglia_Signature_M2 ~ MUT_STATUS, data = ya_can_pheno_2)$p.value
    mt_deg_p_polar <- wilcox.test(Microglia_Polarization ~ MUT_STATUS, data = ya_can_pheno_2)$p.value
    mean_region <- ya_can_pheno_2 %>% group_by(MUT_STATUS) %>% 
      summarise(mean_Microglia_Signature=mean(Microglia_Signature),mean_Microglia_Signature_M1=mean(Microglia_Signature_M1),mean_Microglia_Signature_M2=mean(Microglia_Signature_M2),
                mean_Microglia_Polarization=mean(Microglia_Polarization)) %>% 
      column_to_rownames(var = "MUT_STATUS")%>% 
      t() %>% as.data.frame() %>% mutate(TYPE=ifelse(InRegion>Others,"HighInR","LowInR"))
    mt_region_res <- data.frame(Microglia_Signature_MUT_STATUS.Pval=mt_deg_p,Microglia_Signature_MUT_STATUS_InR_mean=mean_region["mean_Microglia_Signature","InRegion"],Microglia_Signature_MUT_STATUS_Others_mean=mean_region["mean_Microglia_Signature","Others"], 
                               Microglia_Signature_M1_MUT_STATUS.Pval=mt_deg_p_m1,Microglia_Signature_M1_MUT_STATUS_InR_mean=mean_region["mean_Microglia_Signature_M1","InRegion"],Microglia_Signature_M1_MUT_STATUS_Others_mean=mean_region["mean_Microglia_Signature_M1","Others"], #M1
                               Microglia_Signature_M2_MUT_STATUS.Pval=mt_deg_p_m2,Microglia_Signature_M2_MUT_STATUS_InR_mean=mean_region["mean_Microglia_Signature_M2","InRegion"],Microglia_Signature_M2_MUT_STATUS_Others_mean=mean_region["mean_Microglia_Signature_M2","Others"], #M2
                               Microglia_Polarization_MUT_STATUS.Pval=mt_deg_p_polar,Microglia_Polarization_MUT_STATUS_InR_mean=mean_region["mean_Microglia_Polarization","InRegion"],Microglia_Polarization_MUT_STATUS_Others_mean=mean_region["mean_Microglia_Polarization","Others"]  #极化
                               )
    
    
  }else{
    mt_region_res <- data.frame(rep(NA,12)) %>% t() %>% as.data.frame(row.names = 1)
    colnames(mt_region_res) <- c('Microglia_Signature_MUT_STATUS.Pval','Microglia_Signature_MUT_STATUS_InR_mean','Microglia_Signature_MUT_STATUS_Others_mean',
                                 'Microglia_Signature_M1_MUT_STATUS.Pval','Microglia_Signature_M1_MUT_STATUS_InR_mean','Microglia_Signature_M1_MUT_STATUS_Others_mean',
                                 'Microglia_Signature_M2_MUT_STATUS.Pval','Microglia_Signature_M2_MUT_STATUS_InR_mean','Microglia_Signature_M2_MUT_STATUS_Others_mean',
                                 'Microglia_Polarization_MUT_STATUS.Pval','Microglia_Polarization_MUT_STATUS_InR_mean','Microglia_Polarization_MUT_STATUS_Others_mean'
                                )
  }
  can_estimate <- read_tsv(paste0("D:/数据/表达谱数据/",can_1,"/cancer_exp_estimate_score.gct"),skip = 2) %>%
    select(!Description) %>% column_to_rownames(var = "NAME") %>% t() %>% as.data.frame() %>% rownames_to_column(var = "variable") %>%
    mutate(variable=gsub("[.]","-",variable))
  ya_can_pheno_3 <- mutate(ya_can_pheno_1, estimate_sample=substr(variable,1,15)) %>%
    left_join(.,can_estimate,by=c("variable")) %>% left_join(.,timer,by=c("estimate_sample"="submitter_id.samples"))

  sele_cells <- c('StromalScore','ImmuneScore','ESTIMATEScore','TumorPurity',colnames(timer)[which(colnames(timer)!="submitter_id.samples")])

  
  for (cell in sele_cells) {
    ya_can_pheno_4 <- select(ya_can_pheno_3,variable,GROUP,cell) #%>% filter(!is.na(cell))
    ya_can_pheno_4 <- ya_can_pheno_4[!is.na(ya_can_pheno_4[,3]),]
    infi_p <- wilcox.test(ya_can_pheno_4[which(ya_can_pheno_4$GROUP=="High"),3],ya_can_pheno_4[which(ya_can_pheno_4$GROUP=="Low"),3])$p.value
    infi_meanHi <- mean(ya_can_pheno_4[which(ya_can_pheno_4$GROUP=="High"),3])
    infi_meanLo <- mean(ya_can_pheno_4[which(ya_can_pheno_4$GROUP=="Low"),3])
    if(infi_meanHi>infi_meanLo){
      infi_deg_type="High"
    }else if(infi_meanHi<infi_meanLo){
      infi_deg_type="Low"
    }else{
      infi_deg_type="equal"
    }
    infi_res <- rbind(infi_res,data.frame(unique(dat_1[,c('File','CDK Inhibitor','target_name','cell_norm_cs','ID','ID2')]),CELL_TYPE=cell,Infi_meanHi=infi_meanHi,Infi_meanLo=infi_meanLo,Infi.Pval=infi_p,Infi_TYPE=infi_deg_type))
  }
  ya_can_pheno_marc <- select(ya_can_pheno_3,variable,MUT_STATUS,GROUP,`Macrophage M1_XCELL`,`Macrophage M2_XCELL`) %>% mutate(Macrophage_Polarization=`Macrophage M2_XCELL`-c(0.65*`Macrophage M1_XCELL`))
  marc_polar_p <- wilcox.test(Macrophage_Polarization ~ MUT_STATUS,ya_can_pheno_marc)$p.value
  mean_region <- ya_can_pheno_marc %>% group_by(MUT_STATUS) %>% 
    summarise(mean_Macrophage_Polarization=mean(Macrophage_Polarization)) %>% 
    column_to_rownames(var = "MUT_STATUS")%>% 
    t() %>% as.data.frame() %>% mutate(TYPE=ifelse(InRegion>Others,"High","Low"))
  marc_polar_region <- data.frame(unique(dat_1[,c('File','CDK Inhibitor','target_name','cell_norm_cs','ID','ID2')]),CELL_TYPE="Macrophage_Polarization_Region",
                           Infi_meanHi=mean_region["mean_Macrophage_Polarization","InRegion"],Infi_meanLo=mean_region["mean_Macrophage_Polarization","Others"],
                           Infi.Pval=marc_polar_p,Infi_TYPE=ifelse(mean_region["mean_Macrophage_Polarization","InRegion"]>mean_region["mean_Macrophage_Polarization","Others"],"High","Low"))
  marc_polar_p <- wilcox.test(Macrophage_Polarization ~ GROUP,ya_can_pheno_marc)$p.value
  mean_group <- ya_can_pheno_marc %>% group_by(GROUP) %>% 
    summarise(mean_Macrophage_Polarization=mean(Macrophage_Polarization)) %>% 
    column_to_rownames(var = "GROUP")%>% 
    t() %>% as.data.frame() %>% mutate(TYPE=ifelse(High>Low,"High","Low"))
  marc_polar_group <- data.frame(unique(dat_1[,c('File','CDK Inhibitor','target_name','cell_norm_cs','ID','ID2')]),CELL_TYPE="Macrophage_Polarization_Group",
             Infi_meanHi=mean_group["mean_Macrophage_Polarization","High"],Infi_meanLo=mean_group["mean_Macrophage_Polarization","Low"],
             Infi.Pval=marc_polar_p,Infi_TYPE=ifelse(mean_group["mean_Macrophage_Polarization","High"]>mean_group["mean_Macrophage_Polarization","Low"],"High","Low"))
  
  
  
  infi_res <- rbind(infi_res,marc_polar_region,marc_polar_group)
  res <- data.frame(unique(dat_1[,c('File','CDK Inhibitor','target_name','cell_norm_cs','ID','ID2')]),fisher.OR=fisher_or,fisher.Pval=fhsier_p,meanHi,meanOthers,DEG_TYPE=deg_type,DEG.Pval=deg_p,Surv.Pval=p.val,Surv.HR=HR,
                    mt_group_res, #小胶质细胞的浸润比例：有microglia、其m1、m2；靶基因高低表达分组的差异
                    mt_region_res, #小胶质细胞的浸润比例：有microglia、其m1、m2；区域内外分组的差异
                    Others_Sample_Num=as.numeric(table(ya_can_pheno_1$MUT_STATUS)["Others"]),InRegion_Sample_Num=as.numeric(table(ya_can_pheno_1$MUT_STATUS)["InRegion"]),
                    Low_Sample_Num=as.numeric(table(ya_can_pheno_1$GROUP)["Low"]),High_Sample_Num=as.numeric(table(ya_can_pheno_1$GROUP)["High"]))
  
  tot_res <- rbind(tot_res,res)
}
write_tsv(tot_res,"./result/TCGA/11.domain突变与药效/Cmap_CLUE/v1/5.0tot_CLUE_phenotype_difference-调试40.5.tsv")
write_tsv(infi_res,"./result/TCGA/11.domain突变与药效/Cmap_CLUE/v1/5.0XCELL_Infiltration_phenotype_difference-调试40.5.tsv")

tot_res <- read_tsv("./result/TCGA/11.domain突变与药效/Cmap_CLUE/v1/5.0tot_CLUE_phenotype_difference-调试40.5.tsv")
aa <- filter(tot_res,as.numeric(fisher.Pval)<0.05,as.numeric(Surv.Pval)<0.05,as.numeric(DEG.Pval)<0.05,grepl("CDK",target_name))

tot_res_cdk_target <- filter(tot_res) %>% mutate(FISHER_STATUS=ifelse(as.numeric(fisher.Pval)<0.05,"Sig","non-Sig"),
                                                                          DEG_STATUS=ifelse(as.numeric(Surv.Pval)<0.05,"Sig","non-Sig"),
                                                                          SURV_STATUS=ifelse(as.numeric(DEG.Pval)<0.05,"Sig","non-Sig"),
                                                                          # MICROGLIA_TARGET_STATUS=ifelse(as.numeric(microglia.target.Pval)<0.05,"Sig","non-Sig"),
                                                                          # MICROGLIA_REGION_STATUS=ifelse(as.numeric(microglia.region.Pval)<0.05,"Sig","non-Sig")
)
write_tsv(tot_res_cdk_target,"./result/TCGA/11.domain突变与药效/Cmap_CLUE/v1/5.1CDK_inhibitor_phenotype_difference.tsv")

tot_res_cdk_target <- read_tsv("./result/TCGA/11.domain突变与药效/Cmap_CLUE/v1/5.1CDK_inhibitor_phenotype_difference.tsv")

lgg_dat <- filter(tot_res_cdk_target,grepl("LGG",File)) %>% mutate(mRNA_FC=as.numeric(meanHi)/as.numeric(meanOthers),
                                                                   microglia_GROUP_FC=as.numeric(Microglia_Signature_GROUP_Hi_mean)/as.numeric(Microglia_Signature_GROUP_Lo_mean),
                                                                   microglia_m1_GROUP_FC=as.numeric(Microglia_Signature_M1_GROUP_Hi_mean)/as.numeric(Microglia_Signature_M1_GROUP_Lo_mean),
                                                                   microglia_m2_GROUP_FC=as.numeric(Microglia_Signature_M2_GROUP_Hi_mean)/as.numeric(Microglia_Signature_M2_GROUP_Lo_mean),
                                                                   Microglia_Polarization_GROUP_FC=as.numeric(Microglia_Polarization_GROUP_Hi_mean)/as.numeric(Microglia_Polarization_GROUP_Lo_mean),
                                                                   
                                                                   microglia_MUT_STATUS_FC=as.numeric(Microglia_Signature_MUT_STATUS_InR_mean)/as.numeric(Microglia_Signature_MUT_STATUS_Others_mean),
                                                                   microglia_m1_MUT_STATUS_FC=as.numeric(Microglia_Signature_M1_MUT_STATUS_InR_mean)/as.numeric(Microglia_Signature_M1_MUT_STATUS_Others_mean),
                                                                   microglia_m2_MUT_STATUS_FC=as.numeric(Microglia_Signature_M2_MUT_STATUS_InR_mean)/as.numeric(Microglia_Signature_M2_MUT_STATUS_Others_mean),
                                                                   Microglia_Polarization_MUT_STATUS_FC=as.numeric(Microglia_Polarization_MUT_STATUS_InR_mean)/as.numeric(Microglia_Polarization_MUT_STATUS_Others_mean)
                                                                   )

lgg_dat$target_name <- factor(lgg_dat$target_name,levels = rev(sort(unique(lgg_dat$target_name))))
p1 <- ggplot(lgg_dat, aes(x = CDK.Inhibitor, y = target_name, fill = as.numeric(cell_norm_cs))) +
  geom_tile(
    color = "#FAFAFA", # 设置小格子边框颜色
    lwd = 0.8, # 设置小格子的边框宽度
    linetype = 1) + # 设置小格子的边框线型
  scale_fill_gradient2(low = '#619CFF',
                       mid = 'white',
                       high = '#F8766D',
                       midpoint = 0)+ #设置中间颜色对应的数值
  geom_text(aes(label = ifelse(as.numeric(fisher.Pval) < 0.05, "*","")), 
            size = 5, # 设置标签文本大小
            color = "black") + # 设置标签文本颜色
  theme_classic()+
  theme(panel.grid = element_blank(),# 去除网格线
        axis.text.x = element_text(angle = 45, hjust = 1))+ 
  coord_fixed() # 设置x轴与y轴的单位长度相等
p1

lgg_dat$ID3 <- 1
p2<- ggplot(lgg_dat, aes(x = ID3, y = target_name, fill = mRNA_FC)) +
  geom_tile(
    color = "#FAFAFA", # 设置小格子边框颜色
    lwd = 0.8, # 设置小格子的边框宽度
    linetype = 1) + # 设置小格子的边框线型
  scale_fill_gradient2(low = '#619CFF',
                       mid = 'white',
                       high = '#F8766D',
                       midpoint = 1)+ #设置中间颜色对应的数值
  geom_text(aes(label = ifelse(as.numeric(DEG.Pval) < 0.05, "*","")), 
            size = 5, # 设置标签文本大小
            color = "black") + # 设置标签文本颜色
  theme_minimal()+
  theme(panel.grid = element_blank(),# 去除网格线
        axis.text.x = element_text(angle = 45, hjust = 1))+ 
  xlab("")+ # 去除横轴标题
  ylab("")+ # 去除纵轴标题
  coord_fixed() # 设置x轴与y轴的单位长度相等

p3 <- ggplot(lgg_dat, aes(x = ID3, y = target_name, fill = Surv.HR)) +
  geom_tile(
    color = "#FAFAFA", # 设置小格子边框颜色
    lwd = 0.8, # 设置小格子的边框宽度
    linetype = 1) + # 设置小格子的边框线型
  scale_fill_gradient2(low = '#619CFF',
                       mid = 'white',
                       high = '#F8766D',
                       midpoint = 1)+ #设置中间颜色对应的数值
  geom_text(aes(label = ifelse(as.numeric(Surv.Pval) < 0.05, "*","")), 
            size = 5, # 设置标签文本大小
            color = "black") + # 设置标签文本颜色
  theme_minimal()+
  theme(panel.grid = element_blank(),# 去除网格线
        axis.text.x = element_text(angle = 45, hjust = 1))+ 
  xlab("")+ # 去除横轴标题
  ylab("")+ # 去除纵轴标题
  coord_fixed() # 设置x轴与y轴的单位长度相等

p4 <- p1+p2+p3 + plot_layout(guides = 'collect')
p4 
ggsave("./result/TCGA/11.domain突变与药效/Cmap_CLUE/v1/6.1.1LGG_EGFR_1-cs_deg_surv.pdf",p4, width = 10,height=8)

sele_example <- filter(tot_res_cdk_target,  as.numeric(fisher.Pval)<0.05,as.numeric(DEG.Pval)<0.05,as.numeric(Surv.Pval)<0.05,DEG_TYPE=="highHi",as.numeric(Surv.HR)>1) %>% as.data.frame()
sele_example$File %>% unique()

lgg_egfr_sample <- filter(input_dat_sample, File=="LGG_EGFR_1")
lgg_mut <- read_tsv("D:/数据/TCGA突变数据/TCGA_mut/LGG_mut.txt") %>% mutate(submitter_id.samples=substr(Tumor_Sample_Barcode,1,16),CANCER_TYPE="LGG")
lgg_pheno <- read_tsv(paste0("D:/数据/表达谱数据/LGG/TCGA-LGG.GDC_phenotype.tsv")) %>%
  select(submitter_id.samples,age_at_initial_pathologic_diagnosis) %>%
  mutate(CANCER_TYPE="LGG",AGE_GROUP=ifelse(as.numeric(age_at_initial_pathologic_diagnosis)<=60,"YA",ifelse(as.numeric(age_at_initial_pathologic_diagnosis)>60,"OA","None")))
lgg_mut <- filter(lgg_mut, submitter_id.samples %in% lgg_pheno$submitter_id.samples[which(lgg_pheno$AGE_GROUP=="YA")])

lgg_mut_sample <- select(lgg_mut ,Tumor_Sample_Barcode) %>% unique()
lgg_comb <- expand.grid(Tumor_Sample_Barcode = unique(lgg_mut$Tumor_Sample_Barcode), Hugo_Symbol = unique(lgg_dat$target_name))
lgg_mut_cdk <- filter(lgg_mut,Hugo_Symbol %in% unique(lgg_dat$target_name)) %>% select(Hugo_Symbol,Variant_Classification,Tumor_Sample_Barcode) #%>% mutate(TYPE="MUT")
lgg_mut_cdk_tot <- left_join(lgg_comb,lgg_mut_cdk,by=c("Tumor_Sample_Barcode","Hugo_Symbol")) %>% mutate(Variant_Classification=ifelse(is.na(Variant_Classification),"WT",Variant_Classification))
CDK_mut_sample <- lgg_mut_cdk_tot$Tumor_Sample_Barcode[which(lgg_mut_cdk_tot$Variant_Classification!="WT")]

lgg_mut_cdk_tot_sele <- filter(lgg_mut_cdk_tot, Tumor_Sample_Barcode %in% c(as.character(CDK_mut_sample), #
                                                                            unique(lgg_egfr_sample$Tumor_Sample_Barcode)))
lgg_mut_cdk_tot_sele$Hugo_Symbol <- factor(lgg_mut_cdk_tot_sele$Hugo_Symbol,levels = rev(sort(unique(lgg_mut_cdk_tot_sele$Hugo_Symbol))))
lgg_mut_cdk_tot_sele$Tumor_Sample_Barcode <- factor(lgg_mut_cdk_tot_sele$Tumor_Sample_Barcode,levels = c(unique(as.character(CDK_mut_sample)),#CDK基因突变的样本
                                                                                                         unique(lgg_egfr_sample$Tumor_Sample_Barcode))) #EGFR突变的样本
lgg_mut_cdk_tot_sele$Variant_Classification[which(lgg_mut_cdk_tot_sele$Tumor_Sample_Barcode %in% unique(lgg_egfr_sample$Tumor_Sample_Barcode))] <- "InRegion" #
cols=c(
  "WT"="#BEBDBD","Nonsense_Mutation"="#B2CB73","Missense_Mutation"="#DC5957",
  "3'UTR"="#8BCFD6","Silent"="#83852B", "InRegion"="#EEEFEF"
)
p5 <- ggplot(lgg_mut_cdk_tot_sele,aes(x=Tumor_Sample_Barcode,y=Hugo_Symbol))+
  geom_tile(aes(fill=Variant_Classification),color="white",size=1)+ #color和size分别指定方块边线的颜色和粗细
  scale_x_discrete("",expand = c(0,0))+ #不显示横纵轴的label文本；画板不延长
  scale_y_discrete("",expand = c(0,0))+
  scale_fill_manual(values = cols)+ #指定自定义的颜色
  theme(
    axis.text.x.bottom = element_text(size=10),axis.text.y.left = element_text(size = 12), #修改坐标轴文本大小
    axis.ticks = element_blank(), #不显示坐标轴刻度
    legend.title = element_blank() #不显示图例title
    
  )

ab <- colnames(lgg_dat)[grep("Microglia",colnames(lgg_dat))]
ac <- ab[grep("GROUP",ab)]
ad <- ac[-grep("Pval",ab)]
range(as.matrix(lgg_dat[,c('Microglia_Signature_GROUP_Hi_mean','Microglia_Signature_GROUP_Lo_mean',
                           'Microglia_Signature_M1_GROUP_Hi_mean','Microglia_Signature_M1_GROUP_Lo_mean',
                           'Microglia_Signature_M2_GROUP_Hi_mean','Microglia_Signature_M2_GROUP_Lo_mean',
                           'Microglia_Polarization_GROUP_Hi_mean','Microglia_Polarization_GROUP_Lo_mean')]))


infi_res <- read_tsv("./result/TCGA/11.domain突变与药效/Cmap_CLUE/v1/5.0XCELL_Infiltration_phenotype_difference-调试40.5.tsv")
lgg_infi_res <- filter(infi_res, grepl("LGG_",File),grepl("CDK",target_name)) %>% select(!CDK.Inhibitor&!cell_norm_cs&!ID&!ID2) %>% unique()
lgg_infi_res$target_name <- factor(lgg_infi_res$target_name,levels = rev(sort(unique(lgg_infi_res$target_name))))
lgg_macrophage <- lgg_infi_res %>% filter(grepl("Macrophage_XCELL",CELL_TYPE)|grepl("TumorPurity",CELL_TYPE)) #只用巨噬细胞

microglia_res <- lgg_dat %>% select(File,target_name,microglia_GROUP_FC,Microglia_Signature_GROUP.Pval) %>% unique() %>% mutate(CELL_TYPE="Microglia") %>% 
  rename(Infi.FC=microglia_GROUP_FC,Infi.Pval=Microglia_Signature_GROUP.Pval) %>% 
  mutate(log2FC=log2(Infi.FC))
macrophage_res <- lgg_macrophage %>% mutate(Infi.FC=Infi_meanHi/Infi_meanLo,log2FC=log2(Infi.FC)) %>% select(colnames(microglia_res))
brain_immune_infi <- rbind(microglia_res,macrophage_res)
p12 <- ggplot(brain_immune_infi, aes(x = CELL_TYPE, y = target_name, fill = Infi.FC)) +
  geom_tile(
    color = "#FAFAFA", # 设置小格子边框颜色
    lwd = 0.8, # 设置小格子的边框宽度
    linetype = 1) + # 设置小格子的边框线型
  scale_fill_gradient2(low = '#619CFF',
                       mid = 'white',
                       high = '#F8766D',
                       midpoint = 1)+ #设置中间颜色对应的数值
  geom_text(aes(label = ifelse(as.numeric(Infi.Pval) < 0.05, "*","")),
            size = 5, # 设置标签文本大小
            color = "black") + # 设置标签文本颜色
  theme_classic()+
  theme(panel.grid = element_blank(),# 去除网格线
        axis.text.x = element_text(angle = 45, hjust = 1))+
  coord_fixed() # 设置x轴与y轴的单位长度相等

p6 <- p1+p2+p3+p12+p5+ plot_layout(guides = 'collect',ncol=5)
p6
ggsave("./result/TCGA/11.domain突变与药效/Cmap_CLUE/v1/6.1.2LGG_EGFR_1-cs_deg_surv-microglia_marcophage-target_mut.pdf",p6, width = 20,height=10)