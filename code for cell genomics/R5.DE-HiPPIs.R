library(tidyverse)
library(dbplyr)
library(vroom)
library(survival)
library(survminer)
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

error_proc <- function(x){
  if(inherits(x, "try-error")) {
    x <- NA
    
  }
  return(x)
}
error_proc_2 <- function(x){
  if(inherits(x, "try-error")) {
    cat(paste("An error occurred in ",cancer,",nrow is ",ii,"\n"))
    next  # 为结果赋一个新的值
    
  }
}

tot_mut_sample <- read_tsv("./result/TCGA/0.1癌症突变_样本统计/3.各个癌症中基因的突变样本.tsv") %>% mutate(submitter_id.samples=substr(Tumor_Sample_Barcode,1,16))
timer <- read_csv("D:/数据/TIMER/infiltration_estimation_for_tcga.csv") %>% reshape2::melt() %>% dplyr::rename(submitter_id.samples=cell_type) %>% filter(grepl("_TIMER",variable))

mut_sample <- read_tsv("./result/TCGA/1.突变优选/1.位于domain区域内的样本.tsv")
raw <- read_tsv("./result/TCGA/1.突变优选/1.TCGA_raw_domain_singleCancer.tsv") %>% getRawRes %>%
  mutate(Region=paste(Transcript_ID, from, to, `hmm acc`, `hmm name`, clan,sep=";")) %>%
  filter(as.numeric(nRegion)>=3,as.numeric(FDR)<0.05)

int_inform <- read_tsv("./result/TCGA/4.优选事件参与的互作/4.0优选事件参与的互作.tsv") %>% 
  dplyr::select(Region, INSIDER_ID, Trans_Status, P1_gtf_Hugo_Symbol, P1_gtf_Gene_ID, P2_gtf_Hugo_Symbol, P2_gtf_Gene_ID, P1, P2, IntAct) %>% unique()

cancers <- unique(raw$CANCER_TYPE)

cancers <-sort(cancers)

prognosis_error <- c()
for (cancer in cancers) {
  prog_res_out <- c()
  expr_deg_df_res <- c()
  cell_res_out <- c()
  
  print(cancer)
  
  raw_1 <- filter(raw,CANCER_TYPE==cancer)
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
  file_1 <- list.files(paste0("D:/数据/表达谱数据/",cancer))
  if("normal_exp.tsv" %in% file_1){
    nor_exprs <- read.table(paste0("D:/数据/表达谱数据/",cancer,"/normal_exp.tsv"),header = T,stringsAsFactors = F,sep = '\t',row.names = 1,check.names = F) %>% 
      rownames_to_column(var = "Hugo_Symbol") %>% 
      reshape2::melt() %>% 
      left_join(.,pheno,by=c("variable"="submitter_id.samples")) %>%
      dplyr::select(!age_at_initial_pathologic_diagnosis) %>%
      mutate(Type="Normal")
    
  }else{
    nor_exprs <- NULL
  }
  tot_expr <- rbind(can_exprs,nor_exprs)# %>% filter(Hugo_Symbol=="HDAC4")
  ii=which(raw_1$Hugo_Symbol=="NOD1")
  for (ii in 1:nrow(raw_1)) {
    if(ii %% 100==0){print(ii)}
    raw_2 <- raw_1[ii,]
    age_g <- raw_1$AGE_GROUP[ii]
    gene <- raw_1$Hugo_Symbol[ii]
    mut_s <- raw_2 %>% 
      left_join(.,mut_sample,by=intersect(colnames(.),colnames(mut_sample))) %>%
      mutate(submitter_id.samples=substr(Tumor_Sample_Barcode,1,16))
    
    tot_mut_s_outR <- filter(tot_mut_sample,AGE_GROUP==raw_2$AGE_GROUP,CANCER_TYPE==raw_2$CANCER_TYPE,Hugo_Symbol==raw_2$Hugo_Symbol,Transcript_ID==raw_2$Transcript_ID,!Tumor_Sample_Barcode %in% mut_s$Tumor_Sample_Barcode)
    pheno_type <- filter(pheno,AGE_GROUP==age_g)%>% mutate(MUT_STATUS=ifelse(submitter_id.samples %in% mut_s$submitter_id.samples,"InRegion",ifelse(!submitter_id.samples %in% nor_exprs$variable,"WT","Normal")))
    pheno_type_1 <- filter(pheno_type,MUT_STATUS=="WT") %>% mutate(MUT_STATUS=ifelse(submitter_id.samples %in% tot_mut_s_outR$submitter_id.samples,"OutRegion",MUT_STATUS))
    length(intersect(pheno_type_1$submitter_id.samples[which(pheno_type_1$MUT_STATUS=="OutRegion")],tot_mut_s_outR$submitter_id.samples))
    pheno_type <- filter(pheno_type,MUT_STATUS!="WT") %>% rbind(.,pheno_type_1)
    
    
    
    pheno_type$MUT_STATUS <- factor( pheno_type$MUT_STATUS,levels = c("OutRegion","InRegion","WT","Normal"))
    a1 <- pheno_type[which(pheno_type$MUT_STATUS!="Normal"),] %>% filter(!is.na(OS))
    a1_freq <- table(a1$MUT_STATUS)[which(table(a1$MUT_STATUS)!=0)]
    raw_2_int <- left_join(raw_2,int_inform,by="Region") %>% filter(!is.na(INSIDER_ID))
    gene_lists <- unique(c(raw_2_int$P1_gtf_Hugo_Symbol,raw_2_int$P2_gtf_Hugo_Symbol,raw_2$Hugo_Symbol))
    tot_expr_select <- tot_expr %>% filter(Hugo_Symbol %in% gene_lists) %>% left_join(.,pheno_type,by=c("variable"="submitter_id.samples","AGE_GROUP","OS","OS.time"))
    tot_g_num <- length(gene_lists)
    p<-c()
    for (jj in gene_lists) {
      inR_expr <- filter(tot_expr_select,Hugo_Symbol==jj ,MUT_STATUS=="InRegion",Type=="Cancer",AGE_GROUP==age_g)
      if(nrow(inR_expr)==0){next}
      outR_expr <- filter(tot_expr_select,Hugo_Symbol==jj ,MUT_STATUS=="OutRegion",Type=="Cancer",AGE_GROUP==age_g)
      wt_expr <- filter(tot_expr_select,Hugo_Symbol==jj ,MUT_STATUS=="WT",Type=="Cancer",AGE_GROUP==age_g)
      nor_expr <- filter(tot_expr_select,Hugo_Symbol==jj ,MUT_STATUS=="Normal",Type=="Normal",AGE_GROUP==age_g)
      wt_outR_expr <- rbind(outR_expr,wt_expr)
    
      mean_inR <- mean(inR_expr$value)
      mean_outR <- mean(outR_expr$value)
      mean_wt <- mean(wt_expr$value)
      mean_nor <- mean(nor_expr$value)
      
      mean_wt_outR <- mean(wt_outR_expr$value)
     
      nor_wt_p.val <- wilcox.test(nor_expr$value,wt_expr$value)$p.value %>% try(silent = TRUE) %>% error_proc
      nor_inR_p.val <- wilcox.test(nor_expr$value,inR_expr$value)$p.value %>% try(silent = TRUE) %>% error_proc
      nor_outR_p.val <- wilcox.test(nor_expr$value,outR_expr$value)$p.value %>% try(silent = TRUE) %>% error_proc
      inR_outR_p.val <- wilcox.test(inR_expr$value,outR_expr$value)$p.value %>% try(silent = TRUE) %>% error_proc
      wt_inR_p.val <- wilcox.test(wt_expr$value,inR_expr$value)$p.value %>% try(silent = TRUE) %>% error_proc
      wt_outR_p.val <- wilcox.test(wt_expr$value,outR_expr$value)$p.value %>% try(silent = TRUE) %>% error_proc
      
      nor_wt_outR_p.val <- wilcox.test(nor_expr$value,wt_outR_expr$value)$p.value %>% try(silent = TRUE) %>% error_proc
      
      wt_outR_inR_p.val <- wilcox.test(wt_outR_expr$value,inR_expr$value)$p.value %>% try(silent = TRUE) %>% error_proc
      
      expr_deg_df <- data.frame(Symbol=jj,mean_nor,mean_wt,mean_inR,mean_outR,mean_wt_outR,
                                nor_wt_p.val,nor_inR_p.val,nor_outR_p.val,wt_inR_p.val,wt_outR_inR_p.val,wt_outR_p.val,inR_outR_p.val)
      if(jj==gene){
        expr_deg_df$GROUP="Case"
      }else{
        expr_deg_df$GROUP="Partner"
      }
      expr_deg_df_res <- rbind(expr_deg_df_res,data.frame(raw_2,expr_deg_df))

    }
    
    if(min(a1_freq)<2|length(na.omit(a1_freq["InRegion"]))==0){
      print(paste("Prognosis error in ",ii))
      prognosis_error <- rbind(prognosis_error, raw_2)
      next
    }
    sdf  <- survdiff(Surv(time = OS.time,event = OS ) ~ MUT_STATUS, data =  pheno_type[which(pheno_type$MUT_STATUS!="Normal"),])# %>% try(silent = TRUE) %>% error_proc_2
    p.val = 1 - pchisq(sdf$chisq, length(sdf$n) - 1)
    inR_order <- grep("MUT_STATUS=InRegion",names(sdf$n))
    outR_order <- grep("MUT_STATUS=OutRegion",names(sdf$n))
    wt_order <- grep("MUT_STATUS=WT",names(sdf$n))

    HR1 = (sdf$obs[inR_order]/sdf$exp[inR_order])/(sdf$obs[outR_order]/sdf$exp[outR_order])#相较于out，In的生存状态
    if(length(HR1)==0){HR1=NA}
    
    HR2 = (sdf$obs[inR_order]/sdf$exp[inR_order])/(sdf$obs[wt_order]/sdf$exp[wt_order])#相较于WT，In的生存状态
    if(length(HR2)==0){HR2=NA}
    
    HR3 = (sdf$obs[outR_order]/sdf$exp[outR_order])/(sdf$obs[wt_order]/sdf$exp[wt_order])#相较于WT，Out的生存状态
    if(length(HR3)==0){HR3=NA}
    
    prog_res <- data.frame(Symbol=gene,Surv.Pval=p.val,Surv.HR.between.Out.In=HR1,Surv.HR.between.WT.In=HR2,Surv.HR.between.WT.Out=HR3,
                           Surv.STATUS.between.Out.In=ifelse(HR1>1,"InRegion with poor Survival than OutRegion",ifelse(HR1<1,"InRegion with Better Survival than OutRegion","InRegion with Stable Survival than OutRegion")),
                           Surv.STATUS.between.WT.In=ifelse(HR2>1,"InRegion with poor Survival than WT",ifelse(HR2<1,"InRegion with Better Survival than WT","InRegion with Stable Survival than WT")),
                           Surv.STATUS.between.WT.Out=ifelse(HR3>1,"OutRegion with poor Survival than WT",ifelse(HR3<1,"OutRegion with Better Survival than WT","OutRegion with Stable Survival than WT")))
    prog_res_out <- rbind(prog_res_out,data.frame(raw_2,prog_res[,colnames(prog_res)!="Symbol"]))
    if(p.val<0.1){
      file_n <- raw_2 %>% mutate(ID=paste(CANCER_TYPE, AGE_GROUP, Hugo_Symbol, Transcript_ID, from, to, `hmm acc`, `hmm name`,sep = "_"))
      fit_single_cox <- survfit(Surv(OS.time,OS) ~ MUT_STATUS, pheno_type[which(pheno_type$MUT_STATUS!="Normal"),])
      p1 <- ggsurvplot(fit_single_cox,
                       pval = T,risk.table = F,conf.int = T,#title=ii,#paste0("p = ",data$cox_p_values)
                       xlab = "Time(days)",
                       ncensor.plot = F,linetype = "strata", # Change line type by groups
                       #surv.median.line = "hv", # Specify median survival
                       legend.title = "Group",
                       font.x = c(14,"bold.italic"),
                       font.y = c(14,"bold.italic"),
                       font.legend = c(12,"bold.italic"),
                       font.tickslab = c(12,"bold.italic"),
                       # legend.labs = c("High","Low"),
                       ggtheme = theme_bw(), # Change ggplot2 theme
                       # palette = c("#cd2626", "#0a0a85")
      )
      ggsave(paste0("./result/TCGA/5.优选事件及互作的表达和预后效应/4.3.4优选domain的表达和预后效应_改-优化-改3/1.优选基因的预后/",file_n$ID,".pdf"),p1[[1]],width=5,height=5)
    }
  }
  
  write_tsv(prog_res_out,paste0("./result/TCGA/5.优选事件及互作的表达和预后效应/4.3.4优选domain的表达和预后效应_改-优化-改3/5.1",cancer,"_prognosis.tsv"))
  write_tsv(expr_deg_df_res,paste0("./result/TCGA/5.优选事件及互作的表达和预后效应/4.3.4优选domain的表达和预后效应_改-优化-改3/5.2",cancer,"_Expr_DEG.tsv"))
}

cancers <- c("GBM","LGG","UVM","HNSC","ACC","THCA","PCPG","DLBC","THYM","LUAD","LUSC","MESO","BRCA","ESCA","STAD","LIHC","CHOL","PAAD","READ","COAD","KICH","KIRC","KIRP","BLCA","PRAD","TGCT","CESC","OV","UCEC","UCS","SKCM","LAML","SARC")
tot_dat <- c()
for (cancer in cancers) {
  prog_res_out <- data.table::fread(paste0("./result/TCGA/5.优选事件及互作的表达和预后效应/4.3.4优选domain的表达和预后效应_改-优化-改3/5.1",cancer,"_prognosis.tsv")) %>% unique()
  expr_deg_df_res <- data.table::fread(paste0("./result/TCGA/5.优选事件及互作的表达和预后效应/4.3.4优选domain的表达和预后效应_改-优化-改3/5.2",cancer,"_Expr_DEG.tsv")) %>% select(!E.value&!P.value&!FDR) %>% unique()
  expr_case <- filter(expr_deg_df_res,GROUP=="Case") %>% select(!GROUP&!Symbol)
  colnames(expr_case)[grep("mean",colnames(expr_case))] <- paste0("Case_",colnames(expr_case)[grep("mean",colnames(expr_case))])
  colnames(expr_case)[grep("_p.val",colnames(expr_case))] <- paste0("Case_",colnames(expr_case)[grep("_p.val",colnames(expr_case))])
  
  expr_partner <- filter(expr_deg_df_res,GROUP=="Partner") %>% select(!GROUP)
  colnames(expr_partner)[grep("mean",colnames(expr_partner))] <- paste0("Partner_",colnames(expr_partner)[grep("mean",colnames(expr_partner))])
  colnames(expr_partner)[grep("_p.val",colnames(expr_partner))] <- paste0("Partner_",colnames(expr_partner)[grep("_p.val",colnames(expr_partner))])
  if(nrow(expr_case)+nrow(expr_partner)!=nrow(expr_deg_df_res)){print(paste("Expr Error!!"))}
  
  expr_deg_df_res_2 <- left_join(expr_case,expr_partner,by=intersect(colnames(expr_case),colnames(expr_partner)))

  dat_1 <- full_join(prog_res_out,expr_deg_df_res_2,by=intersect(colnames(prog_res_out),colnames(expr_deg_df_res_2))) #%>%
  tot_dat <- rbind(tot_dat,dat_1)
}
write_tsv(tot_dat,"./result/TCGA/5.优选事件及互作的表达和预后效应/4.3.4优选domain的表达和预后效应_改-优化-改3/5.4tot_Prognosis_Expr-改2.tsv")