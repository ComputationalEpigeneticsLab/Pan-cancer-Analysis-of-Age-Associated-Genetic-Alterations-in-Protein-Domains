library(ggpubr)
library(ggpointdensity)
library(tidyverse)
library(ggplot2)
setwd("D:/课题/project10")
rm(list=ls())
mut_sample <- read_tsv("./result/TCGA/1.突变优选/1.位于domain区域内的样本.tsv")
dnam_age <- read_tsv("F:/参考文献/年龄相关癌症基因组/34879281-泛癌解释与年龄相关的分子模式/Epigenetic_age_TCGA.tsv") %>% 
  mutate(SampleID=gsub("[.]","-",SampleID),submitter_id.samples=substr(SampleID,1,16),
         MUT_STATUS=ifelse(submitter_id.samples %in% substr(mut_sample$Tumor_Sample_Barcode,1,16),"InRegion","Others"),
         MUT_STATUS_int=ifelse(submitter_id.samples %in% substr(mut_sample$Tumor_Sample_Barcode,1,16),1,0))

g1 <- ggscatter(dnam_age, x = "Age", y = "DNAmAge",size = 1,
               add = "reg.line", conf.int = F,
               cor.coef = TRUE, cor.method = "pearson",#color = "density",
               xlab = " Age", ylab = " DNAmAge",add.params = list(color="black"),cor.coef.size = 6) +
  geom_pointdensity(aes(x=Age,y=DNAmAge),size=2) +
  scale_colour_gradientn(colours = terrain.colors(10))+
  theme_classic()
g1

ggscatter(dnam_age[which(dnam_age$MUT_STATUS=="InRegion"),], x = "Age", y = "DNAmAge",size = 1,
          add = "reg.line", conf.int = F,
          cor.coef = TRUE, cor.method = "pearson",#color = "density",
          xlab = " Age", ylab = " DNAmAge",add.params = list(color="black"),cor.coef.size = 6) +
  geom_pointdensity(aes(x=Age,y=DNAmAge),size=3)+
  geom_smooth(aes(x=Age,y=DNAmAge),color = "black" ,method = "lm", se = FALSE) +
  scale_colour_gradientn(colours = terrain.colors(10))+
  theme_classic() #+theme(axis.text=element_text(size=20))

ggscatter(dnam_age[which(dnam_age$MUT_STATUS=="Others"),], x = "Age", y = "DNAmAge",size = 1,
          add = "reg.line", conf.int = F,
          cor.coef = TRUE, cor.method = "pearson",#color = "density",
          xlab = " Age", ylab = " DNAmAge",add.params = list(color="black"),cor.coef.size = 6) +
  geom_pointdensity(aes(x=Age,y=DNAmAge),size=3)+
  geom_smooth(aes(x=Age,y=DNAmAge),color = "black" ,method = "lm", se = FALSE) +
  scale_colour_gradientn(colours = terrain.colors(10))+
  theme_classic() #+theme(axis.text=element_text(size=20))

ggscatter(dnam_age, x = "Age", y = "AgeAccelerationDiff",size = 1,
          # add = "reg.line", conf.int = F,
          # cor.coef = TRUE, cor.method = "pearson",#color = "density",
          xlab = " Age", ylab = " AgeAccelerationDiff",add.params = list(color="black"),cor.coef.size = 6) +
  geom_pointdensity(aes(x=Age,y=AgeAccelerationDiff),size=3)+
  geom_smooth(aes(x=Age,y=AgeAccelerationDiff,group=MUT_STATUS),color = "black" ,method = "lm", se = FALSE) +
  scale_colour_gradientn(colours = terrain.colors(10))+
  theme_classic() #+theme(axis.text=element_text(size=20))



g3 <- ggscatter(dnam_age, x = "Age", y = "AgeAccelerationDiff", 
          color = "MUT_STATUS", # 根据"group"列来分组
          palette = c("InRegion"="#C12324","Others"="#7A9DAE"), # 为不同组设置颜色
          alpha = 0.3,
          add = "reg.line", # 添加回归线
          conf.int = TRUE, # 添加回归线的置信区间
          cor.coef = TRUE, cor.method = "pearson",
          legend.title = "Group")
g3

g3_1 <- ggscatter(dnam_age[which(dnam_age$MUT_STATUS=="InRegion"),], x = "Age", y = "AgeAccelerationDiff", 
          color = "MUT_STATUS", # 根据"group"列来分组
          palette = c("InRegion"="#C12324"), # 为不同组设置颜色
          alpha = 0.2,
          add = "reg.line", # 添加回归线
          conf.int = TRUE, # 添加回归线的置信区间
          cor.coef = TRUE, cor.method = "pearson",
          legend.title = "Group")

g3_2 <- ggscatter(dnam_age[which(dnam_age$MUT_STATUS=="Others"),], x = "Age", y = "AgeAccelerationDiff", 
          color = "MUT_STATUS", # 根据"group"列来分组
          palette = c("Others"="#7A9DAE"), # 为不同组设置颜色
          alpha = 0.1,
          add = "reg.line", # 添加回归线
          conf.int = TRUE, # 添加回归线的置信区间
          cor.coef = TRUE, cor.method = "pearson",
          legend.title = "Group")
g3/g3_1/g3_2
g4 <- g3/g3_1/g3_2
ggsave("./result/TCGA/10.表观年龄加速/5.1表观修饰年龄指标与实际年龄的相关性-AgeAccelerationDiff.pdf",g4,width = 6,height = 15)
