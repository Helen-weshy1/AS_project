# Fig1
library(Seurat)
library(ggplot2)
library(data.table)
library(dplyr)

# prepare data
law_sc_gene_all

##### panel A #####
DimPlot(object = law_sc_gene_all, reduction = "tsne", pt.size = 2, label = T,
        cols = c( 
          "#1F77B4FF", "#FF7F0EFF", "#bebebe", "#2CA02CFF" ,"#D62728FF", "#9467BDFF",
          "#AEC7E8FF", "#FFBB78FF", "#98DF8AFF", "#FF9896FF",
          "#8C564BFF" ,"#E377C2FF", "#C49C94FF", "#BCBD22FF", "#17BECFFF"))
DimPlot(object = law_sc_gene_all, reduction = "tsne", pt.size = 2, label = T,
        split.by = 'Status',
        cols = c( 
          "#1F77B4FF", "#FF7F0EFF", "#bebebe", "#2CA02CFF" ,"#D62728FF", "#9467BDFF",
          "#AEC7E8FF", "#FFBB78FF", "#98DF8AFF", "#FF9896FF",
          "#8C564BFF" ,"#E377C2FF", "#C49C94FF", "#BCBD22FF", "#17BECFFF"))

##### panel B #####
#### celltype barplot 
ct_sta_all=as.data.frame(table(law_sc_gene_all@meta.data$celltype,law_sc_gene_all@meta.data$Status))
ggplot(ct_sta_all,aes(x=Var2,y=Freq,fill=Var1))+  
  geom_bar( stat="identity", position="fill",width = 0.7,size=1)+
  scale_fill_manual(values=c("#1F77B4FF", "#FF7F0EFF", "#2CA02CFF","#D62728FF",  "#a020ef",
                             "#8C564BFF" ,"#E377C2FF","#C49C94FF", "#BCBD22FF"))+
  theme_bw()+
  ylab('Percentage')

##### panel C #####
DimPlot(object = law_na100, reduction = "tsne", pt.size = 5, label = T,
        cols = c(
          "#00688b", "#ff3030", "#00ffff","#4876ff", 
          "#ff1494","#00ff7f", "#ce960d",
          "#000080", "#a020ef",
          "#64b8fe","#BFEFFF", "#8DB6CD",
          "#ffc0cb", 
          "#99fa99", "#cebe71","#ff6ab4"))
ggsave('all_as_tsne.pdf', width = 7.5, height = 6.5)

##### panel D #####
DimPlot(object = law_na100, reduction = "tsne", pt.size = 5, label = T,
        group.by = 'celltype',
        cols = c("#8C564BFF" ,"#1F77B4FF", "#FF7F0EFF", "#D62728FF",
                 "#E377C2FF","#9467BDFF" , "#2CA02CFF", "#a020ef", "#BCBD22FF" ))
ggsave('all_as_ct_tsne.pdf', width = 8, height = 6.5)
