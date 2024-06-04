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
raw <- read_tsv("D:/课题/project10/result/TCGA/1.突变优选/1.TCGA_raw_domain_singleCancer.tsv") %>% getRawRes %>% 
  mutate(Region=paste(Transcript_ID, from, to, `hmm acc`, `hmm name`, clan,sep=";")) %>%
  filter(as.numeric(nRegion)>=3,as.numeric(FDR)<0.05)

gtf <- read_tsv("D:/数据/基因组/GENCODE/gencode.v19.gtf.basic.txt") %>% filter(Type=="gene",Gene_Type=="protein_coding") %>% select(Gene_ID,Hugo_Symbol) %>% unique()

Net <- read_tsv("D:/数据/蛋白互作/其它/HPRD_all.txt",col_names = F) %>% dplyr::rename(P1=X1,P2=X2)
nrow(Net[,1:2])
# [1] 36867
length(unique(c(Net$P1,Net$P2))) #[1] 9453
length(unique(c(Net$X3,Net$X4))) #[1] 9453
library(igraph)
humanP_network<-make_graph(t(Net[,c("P1","P2")]),directed = F)


all_humanP_degree<-degree(humanP_network,mode = "total")
all_humanP_degree<-as.data.frame(all_humanP_degree) %>% rownames_to_column(var = "Hugo_Symbol")

YA_degree <-all_humanP_degree[which(all_humanP_degree$Hugo_Symbol %in% raw$Hugo_Symbol[which(raw$AGE_GROUP=="YA")]),]
YA_degree$group<-c("YA")

OA_degree <-all_humanP_degree[which(all_humanP_degree$Hugo_Symbol %in% raw$Hugo_Symbol[which(raw$AGE_GROUP=="OA")]),]
OA_degree$group<-c("OA")
others1 <- all_humanP_degree[-which(all_humanP_degree$Hugo_Symbol %in% unique(c(raw$Hugo_Symbol[which(raw$AGE_GROUP=="YA")],raw$Hugo_Symbol[which(raw$AGE_GROUP=="OA")]))),]
others1$group<-c("others")

main_degree_data <-rbind(YA_degree,OA_degree,others1) %>% rename(degree=all_humanP_degree)
main_degree_data$log<-log2(main_degree_data$degree)

library(ggplot2)
library(ggthemes)
library(ggprism)
library(ggpubr)
my_theme2<-theme_prism()+ #背景变为白色
  theme(axis.text.x=element_text(colour="black",family="Times",size=10), #设置x轴刻度标签的字体显示倾斜角度为15度，并向下调整1(hjust = 1)，字体簇为Times大小为20
        axis.title.x=element_text(family="Times",size = 13),
        axis.text.y=element_text(family="Times",size=10,face="plain"), #设置y轴刻度标签的字体簇，字体大小，字体样式为plain
        axis.title.y=element_text(family="Times",size = 13,), #设置y轴标题的字体属性
        panel.border = element_blank(),axis.line = element_line(colour = "black",size=1), #去除默认填充的灰色，并将x=0轴和y=0轴加粗显示(size=1)
        legend.text=element_text( family="Times", colour="black",  #设置图例的子标题的字体属性
                                  size=12),
        legend.title=element_text( family="Times", colour="black", #设置图例的总标题的字体属性
                                   size=14),
        panel.grid.major = element_blank(),   #不显示网格线
        panel.grid.minor = element_blank())


my_comparisons<-list(c("YA","OA"),c("OA","others"),c("YA","others"))
main_degree_data$group <- factor(main_degree_data$group, levels = c("YA","OA","others"))
p1 <- ggplot(main_degree_data,
             aes(x=group,y=log))+
  geom_violin(aes(fill=group),trim = F)+
  geom_boxplot(width=0.1,fill="white",alpha=1,outlier.shape = NA)+
  scale_fill_manual(values = c("YA"="#4A579C","OA"="#EA5B5A","Others"="#E5E4E2"))+
  stat_compare_means(comparisons=my_comparisons,label.y = c(14,15,16))+
  stat_compare_means(label.y = 16)+
  labs(x='gene type',y='log10(degree)')+
  my_theme2
p1 
ggsave("./result/TCGA/1.突变优选/调试36.1优选事件的网络拓扑属性/HPRD_degree.pdf",p1,width = 5, height = 4)




all_humanP_closeness<-closeness(humanP_network,mode="in")###in表示无向图
all_humanP_closeness<-as.data.frame(all_humanP_closeness) %>% rownames_to_column(var = "Hugo_Symbol")


YA_closeness<-all_humanP_closeness[which(all_humanP_closeness$Hugo_Symbol %in% raw$Hugo_Symbol[which(raw$AGE_GROUP=="YA")]),]
YA_closeness$group<-c("YA")

OA_closeness<-all_humanP_closeness[which(all_humanP_closeness$Hugo_Symbol %in% raw$Hugo_Symbol[which(raw$AGE_GROUP=="OA")]),]
OA_closeness$group<-c("OA")

others1<-all_humanP_closeness[-which(all_humanP_closeness$Hugo_Symbol %in%  unique(c(raw$Hugo_Symbol[which(raw$AGE_GROUP=="YA")],raw$Hugo_Symbol[which(raw$AGE_GROUP=="OA")]))),]
others1$group<-c("others")

mian_closeness_data<-rbind(YA_closeness,OA_closeness,others1) %>% rename(closeness=all_humanP_closeness)
mian_closeness_data$log<-log10(mian_closeness_data$closeness)



compare_means(log~group,data = mian_closeness_data)
my_comparisons<-list(c("YA","OA"),c("OA","others"),c("YA","others"))
mian_closeness_data$group <- factor(mian_closeness_data$group, levels = c("YA","OA","others"))

p2 <- ggplot(mian_closeness_data,
             aes(x=group,y=log))+
  geom_violin(aes(fill=group),trim = F)+
  geom_boxplot(width=0.1,fill="white",alpha=1,outlier.shape = NA)+
  scale_fill_manual(values =  c("YA"="#4A579C","OA"="#EA5B5A","Others"="#E5E4E2"))+
  stat_compare_means(comparisons=my_comparisons,label.y = c(-6.376,-6.376,-6.376))+
  stat_compare_means(label.y = -6.375)+
  ylim(c(-6.39,-6.37)) +
  labs(x='gene type',y='log10(closeness)')+
  my_theme2
p2
ggsave("./result/TCGA/1.突变优选/调试36.1优选事件的网络拓扑属性/HPRD_closeness.pdf",p2,width = 5, height = 4)


all_humanP_betweenness<-betweenness(humanP_network,normalized = T)
all_humanP_betweenness<-as.data.frame(all_humanP_betweenness)%>% rownames_to_column(var = "Hugo_Symbol")

YA_betweenness<-all_humanP_betweenness[which(all_humanP_betweenness$Hugo_Symbol %in% raw$Hugo_Symbol[which(raw$AGE_GROUP=="YA")]),]
YA_betweenness$group<-c("YA")

OA_betweenness<-all_humanP_betweenness[which(all_humanP_betweenness$Hugo_Symbol %in% raw$Hugo_Symbol[which(raw$AGE_GROUP=="OA")]),]
OA_betweenness$group<-c("OA")

others1<-all_humanP_betweenness[-which(all_humanP_betweenness$Hugo_Symbol %in% unique(c(raw$Hugo_Symbol[which(raw$AGE_GROUP=="YA")],raw$Hugo_Symbol[which(raw$AGE_GROUP=="OA")]))),]
others1$group<-c("others")

main_betweenness_data<-rbind(YA_betweenness,OA_betweenness,others1) %>% rename(betweenness=all_humanP_betweenness)
main_betweenness_data$log<-log10(main_betweenness_data$betweenness)

compare_means(log~group,data = main_betweenness_data)
my_comparisons<-list(c("YA","OA"),c("OA","others"),c("YA","others"))

main_betweenness_data$group <- factor(main_betweenness_data$group, levels = c("YA","OA","others"))

p3 <- ggplot(main_betweenness_data,
             aes(x=group,y=log))+
  geom_violin(aes(fill=group),trim = F)+
  geom_boxplot(width=0.1,fill="white",alpha=1,outlier.shape = NA)+
  scale_fill_manual(values = c("YA"="#4A579C","OA"="#EA5B5A","Others"="#E5E4E2"))+
  stat_compare_means(comparisons=my_comparisons,label.y = c(0,1,1))+
  stat_compare_means(label.y = -2.5)+
  ylim(c(-8,2)) +
  labs(x='gene type',y='log10(betweenness)')+
  my_theme2
p3
ggsave("./result/TCGA/1.突变优选/调试36.1优选事件的网络拓扑属性/HPRD_betweenness.pdf",p3,width = 5, height = 4)

