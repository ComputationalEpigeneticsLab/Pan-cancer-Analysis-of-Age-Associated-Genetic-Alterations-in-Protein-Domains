library(tidyverse)
library(dbplyr)
library(vroom)
library(VennDiagram)
library (venn)
library(clusterProfiler)
library(org.Hs.eg.db)
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

print(cancer)
YA <- filter(raw,AGE_GROUP=="YA")
OA <- filter(raw,AGE_GROUP=="OA")
#YA
ya_res <- enrichGO(gene = unique(YA$Hugo_Symbol),
                   OrgDb = "org.Hs.eg.db",
                   keyType = "SYMBOL",#这里指定ID类型
                   ont = "ALL", # "BP", "BP", "CC" 
                   pvalueCutoff = 0.05,
                   pAdjustMethod = "BH",
                   qvalueCutoff = 0.05,
                   minGSSize = 10,# 最少的基因数量
                   maxGSSize = 500, # 最大的基因数量
                   readable = T # 把ENTREZID转换为SYMBOL
)
class(ya_res)
write_tsv(as.data.frame(ya_res),"./result/TCGA/1.突变优选/1.2功能富集/simplifyEnrichment/GO_ALL_YA-功能富集.tsv")
go_id_bp <- ya_res[ya_res$ONTOLOGY == "BP", "ID"]

length(go_id_bp);length(go_id_cc);length(go_id_bp)

library(simplifyEnrichment)
mat <- GO_similarity(go_id_bp, ont = "BP", db="org.Hs.eg.db")

pdf("./result/TCGA/1.突变优选/1.2功能富集/simplifyEnrichment/GO_BP_YA-功能富集.pdf",width = 8,height = 5)
df <- simplifyGO(mat, plot = T)
dev.off()



#OA
oa_res <- enrichGO(gene = unique(OA$Hugo_Symbol),
                   OrgDb = "org.Hs.eg.db",
                   keyType = "SYMBOL",#这里指定ID类型
                   ont = "ALL", # "BP", "BP", "CC" 
                   pvalueCutoff = 0.05,
                   pAdjustMethod = "BH",
                   qvalueCutoff = 0.05,
                   minGSSize = 10,# 最少的基因数量
                   maxGSSize = 500, # 最大的基因数量
                   readable = T # 把ENTREZID转换为SYMBOL
)
class(oa_res)
write_tsv(as.data.frame(oa_res),"./result/TCGA/1.突变优选/1.2功能富集/simplifyEnrichment/GO_ALL_OA-功能富集.tsv")

go_id_bp <- oa_res[oa_res$ONTOLOGY == "BP", "ID"]

length(go_id_bp);length(go_id_cc);length(go_id_bp)

library(simplifyEnrichment)
mat <- GO_similarity(go_id_bp, ont = "BP", db="org.Hs.eg.db")

# 聚类并画图
pdf("./result/TCGA/1.突变优选/1.2功能富集/simplifyEnrichment/GO_BP_OA-功能富集.pdf",width = 8,height = 5)
df <- simplifyGO(mat, plot = T)
dev.off()

