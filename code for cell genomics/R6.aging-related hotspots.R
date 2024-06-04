library(tidyverse)
library(vroom)
library(venn)        
library(VennDiagram) 
library(cols4all)
setwd("D:/课题/project10")
rm(list=ls())
gtf <- vroom("D:/数据/基因组/GENCODE/gencode.v19.gtf.basic.txt") %>% filter(Gene_Type=="protein_coding",Type=="gene")
int_res <-  read_tsv("./result/TCGA/5.优选事件及互作的表达和预后效应/4.3.4优选domain的表达和预后效应_改-优化-改3/5.4tot_Prognosis_Sig_Expr-改2.tsv") %>% mutate(ID=paste(CANCER_TYPE,AGE_GROUP,Hugo_Symbol,Region,sep = ";")) #%>%
pathways <- clusterProfiler::read.gmt("D:/数据/mSigDB/msigdb.v2022.1.Hs.symbols.gmt") %>% filter(grepl("SENESCENCE",term))
backgroud <- unique(gtf$Hugo_Symbol)

id <- unique(int_res$ID)
name <- unique(pathways$term)

enrich_result <- c()
for (ii in id) {
  for (jj in name) {
    ints <- filter(int_res, ID==ii)
    int_lists <- unique(c(ints$Hugo_Symbol,ints$Symbol))
    path_1 <- filter(pathways,term==jj)
     both_gene <- intersect(int_lists,path_1$gene)
    if(length(both_gene)==0){
      dat <- data.frame(Pathway_Name=jj,unique(ints[,c("CANCER_TYPE","AGE_GROUP","Hugo_Symbol","Transcript_ID","from","to","hmm.acc","hmm.name","clan","nRegion","nTot","LRegion","AA_length","Region","ID")]),
                        count=0,rich_score=NA,Pvalue=1)
      enrich_result <- rbind(enrich_result,dat)
      next
    }
    back_in_gene <- intersect(backgroud,path_1$gene)
    rich_score <- length(both_gene)/length(path_1$gene)
    pvalue=1-phyper(length(both_gene),# 差异基因中，位于通路中基因数量
                    length(int_lists), # 差异基因的数量
                    length(backgroud)-length(int_lists), # 全部基因的数量 - 差异数量
                    length(back_in_gene))  # 全部基因中，位于通路中基因数量
    dat <- data.frame(Pathway_Name=jj,unique(ints[,c("CANCER_TYPE","AGE_GROUP","Hugo_Symbol","Transcript_ID","from","to","hmm.acc","hmm.name","clan","nRegion","nTot","LRegion","AA_length","Region","ID")]),
                      count=length(both_gene),rich_score=rich_score,Pvalue=pvalue)
    enrich_result <- rbind(enrich_result,dat)
  }
}
write_tsv(enrich_result,"./result/TCGA/10.表观年龄加速/10.3年龄表达差异网络的功能富集-SENESCENCE.tsv")
enrich_result <- read_tsv("./result/TCGA/10.表观年龄加速/10.3年龄表达差异网络的功能富集-SENESCENCE.tsv")
sig_enrich_ya <- enrich_result  %>% filter(count>=2,AGE_GROUP=="YA") %>% group_by(CANCER_TYPE) %>% mutate(FDR=p.adjust(as.numeric(Pvalue),method = "BH"))  %>% mutate(logFDR=-log10(FDR+2.819966e-15))  %>% ungroup()
range(sig_enrich_ya$FDR)

unique(sig_enrich_ya$ID)
int_res_sig_enrich_ya <- filter(int_res, ID %in% sig_enrich_ya$ID)
sig_enrich_oa <- enrich_result  %>% filter(count>=2,AGE_GROUP=="OA") %>% group_by(CANCER_TYPE) %>% mutate(FDR=p.adjust(as.numeric(Pvalue),method = "BH"))  %>% mutate(logFDR=-log10(FDR+2.819966e-15)) %>% ungroup()
unique(sig_enrich_oa$ID)

int_res_sig_enrich_oa <- filter(int_res, ID %in% sig_enrich_oa$ID)

colors=c('LUAD'='#F2C1A6','LUSC'='#DAA682','MESO'='#F1A776','BLCA'='#D3AB95','KICH'='#f0cfb7','KIRC'='#ddac9e','KIRP'='#ef9767','PRAD'='#f1bda6','TGCT'='#BA5D2A',
        'COAD'='#bbddb4','ESCA'='#6eab59','READ'='#95c487','STAD'='#499436','CESC'='#afd9f0','OV'='#61c0e1','UCEC'='#92cde6','UCS'='#3bbcdf','ACC'='#f9dcdf','PCPG'='#ea8d8b','THCA'='#e25367',
        'CHOL'='#67c44d','LIHC'='#41960F','PAAD'='#4CC642',
        'DLBC'='#8d78b5','LAML'='#c8bddc','THYM'='#50479a','GBM'='#efe88a','LGG'='#e5d839','SARC'='#3e3a39','SKCM'='#8082a4','HNSC'='#9d7873','UVM'='#b9a08f','BRCA'='#bcbbbb')
colors_sele <- colors[unique(sig_enrich_ya$CANCER_TYPE)]
#仅针对YA绘图
sig_enrich_ya$Pathway_Name <- factor(sig_enrich_ya$Pathway_Name,levels = sort(unique(sig_enrich_ya$Pathway_Name)))
p1 <- ggplot(sig_enrich_ya, aes(x = Pathway_Name, y = logFDR,color=CANCER_TYPE)) +
  geom_jitter(width = 0.4, height = 0)  + scale_x_discrete(labels= gsub("_","\n",sort(unique(sig_enrich_ya$Pathway_Name)))) +
  geom_hline(yintercept=-log10(0.05+2.819966e-15)) +  
  geom_vline(xintercept = c(1:length(unique(sig_enrich_ya$Pathway_Name))) - 0.5, color = "black") +
  scale_color_manual(values=colors_sele ) + theme_classic()#+ theme(axis.text.x = element_text(angle = 45, hjust = 1))
p1 
ggsave("./result/TCGA/10.表观年龄加速/6.1对照-年龄表达差异网络的功能富集-SENESCENCE.pdf",p1,width = 18,height = 4)
p2 <- ggplot(sig_enrich_ya, aes(x = Pathway_Name, y = logFDR,color=CANCER_TYPE)) +
  geom_jitter(width = 0.4, height = 0)  +
  geom_hline(yintercept=-log10(0.05+2.819966e-15)) +  
  geom_vline(xintercept = c(1:length(unique(sig_enrich_ya$Pathway_Name))) - 0.5, color = "black") +
  scale_color_manual(values=colors_sele ) + theme_classic()+ theme(axis.text.x = element_text(angle = 45, hjust = 1))
p2 
ggsave("./result/TCGA/10.表观年龄加速/6.1留横坐标-年龄表达差异网络的功能富集-SENESCENCE.pdf",p2,width = 10,height = 4)

select(sig_enrich_ya,Pathway_Name,Hugo_Symbol) %>% unique() %>% group_by(Pathway_Name) %>% summarise(Symbols=paste0(Hugo_Symbol,collapse = ","))


library(ggsankey)
library(cols4all)
sankey <- sig_enrich_ya %>% select(Pathway_Name,CANCER_TYPE,Hugo_Symbol) %>% unique() %>% arrange(Pathway_Name)  %>% mutate(id=1:nrow(.)) %>% select(id,everything())
axis_df <- as.data.frame(table(sankey$Pathway_Name)) %>% 
  mutate(xmax= cumsum(Freq)) %>%
  mutate(xmin = xmax -Freq) %>%
  mutate(label = (xmin + xmax)/2)
sig_enrich_ya$Pathway_Name <- factor(sig_enrich_ya$Pathway_Name,levels = sort(unique(sig_enrich_ya$Pathway_Name)))
plot_dat <- sig_enrich_ya %>% left_join(.,axis_df,by=c("Pathway_Name"="Var1"))
p1 <- ggplot(plot_dat, aes(x = label, y = logFDR,color=CANCER_TYPE)) + #xlim(-3,160) + 
  geom_jitter(width = 1, height = 0,size=3)  +  scale_x_discrete(labels= unique(sort(plot_dat$label))) +#scale_x_discrete(labels= gsub("_","\n",sort(unique(sig_enrich_ya$Pathway_Name)))) +
  geom_hline(yintercept=-log10(0.05+2.819966e-15)) +  
  scale_color_manual(values=colors_sele ) + theme_classic()#+ theme(axis.text.x = element_text(angle = 45, hjust = 1))
p1 
ggsave("./result/TCGA/10.表观年龄加速/6.1年龄表达差异网络的功能富集-SENESCENCE.pdf",p1,width = 18,height = 4)



order_id <- openxlsx::read.xlsx("./result/TCGA/人体图-final3癌症归属系统及排序.xlsx")
order_id$CANCER_TYPE[order_id$CANCER_TYPE%in%sankey$CANCER_TYPE]

df_1 <- sankey %>%
  make_long(Hugo_Symbol,CANCER_TYPE, Pathway_Name) 
head(df_1)
df_1$node <- factor(df_1$node,levels = c(
  sankey$Pathway_Name %>% unique()%>% as.character %>% rev(),
  order_id$CANCER_TYPE[order_id$CANCER_TYPE%in%sankey$CANCER_TYPE] %>% rev(),
  sankey$Hugo_Symbol %>% unique() %>% rev())
)

colors=c('LUAD'='#F2C1A6','LUSC'='#DAA682','MESO'='#F1A776','BLCA'='#D3AB95','KICH'='#f0cfb7','KIRC'='#ddac9e','KIRP'='#ef9767','PRAD'='#f1bda6','TGCT'='#BA5D2A',
         'COAD'='#bbddb4','ESCA'='#6eab59','READ'='#95c487','STAD'='#499436','CESC'='#afd9f0','OV'='#61c0e1','UCEC'='#92cde6','UCS'='#3bbcdf','ACC'='#f9dcdf','PCPG'='#ea8d8b','THCA'='#e25367',
         'CHOL'='#67c44d','LIHC'='#41960F','PAAD'='#4CC642',
         'DLBC'='#8d78b5','LAML'='#c8bddc','THYM'='#50479a','GBM'='#efe88a','LGG'='#e5d839','SARC'='#3e3a39','SKCM'='#8082a4','HNSC'='#9d7873','UVM'='#b9a08f','BRCA'='#bcbbbb')
colors_df <- as.data.frame(colors)%>% rownames_to_column(var="CANCER_TYPE")
mycol_df <- data.frame(Node=c(unique(sankey$Hugo_Symbol),rev(order_id$CANCER_TYPE[order_id$CANCER_TYPE%in%sankey$CANCER_TYPE]),as.character(unique(sankey$Pathway_Name))),
                       Type=c(rep("Hugo_Symbol",length(unique(sankey$Hugo_Symbol))),rep("CANCER_TYPE",length(unique(sankey$CANCER_TYPE))),rep("Pathway_Name",length(unique(sankey$Pathway_Name)))),
                       colors=NA)

mycol_df$colors[na.omit(match(colors_df$CANCER_TYPE,mycol_df$Node))] <- colors_df$colors[which(!is.na(match(colors_df$CANCER_TYPE,mycol_df$Node)))]
mycol_df$colors[na.omit(match(colors_df$CANCER_TYPE,mycol_df$Node))]
mycol_df$colors[which(mycol_df$Type=="Hugo_Symbol")] <- "gray"
c4acol <- c4a('rainbow_wh_rd',length(as.character(unique(sankey$Pathway_Name))))
mycol_df$colors[which(mycol_df$Type=="Pathway_Name")] <- c4acol
mycol <- mycol_df$colors
p4 <- ggplot(df_1, aes(x = x,
                       next_x = next_x,
                       node = node,
                       next_node = next_node,
                       fill = node,
                       label = node)) +
  geom_sankey(flow.alpha = 0.5,
              flow.fill = 'grey',
              flow.color = 'grey80', #条带描边色
              node.fill = mycol, #节点填充色
              smooth = 8,
              width = 0.08) +
  geom_sankey_text(size = 3.2,
                   color = "black")+
  theme_void() +
  theme(legend.position = 'none')
ggsave("./result/TCGA/10.表观年龄加速/7.1桑基图-Hugo_Symbol-CANCER_TYPE-Pathway_Name.pdf",p4,width = 8,height=10)
