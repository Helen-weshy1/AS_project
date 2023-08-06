rm(list=ls())
setwd('/media/user/sdh/cmq_wj/others')
.libPaths("/media/user/sdf/R.lib/library/wj/sc")
library(Seurat)
library(ggplot2)
AS1_4<-readRDS('/media/user/sdh/cmq_wj/others/c1_4_cytotrace_as.rds')

pdf('mature.pdf',width = 10,height = 40)
VlnPlot(AS1_4,features=c('INS','SLC2A2','MAFA','RFX6','PDX1','CHL1'),group.by = 'new_as_cluster')
dev.off()


pdf('AS.pdf',width = 10,height = 40)
VlnPlot(AS1_4,features=c('FXR1','FMR1','HNRNPH2','HNRNPK'),group.by = 'new_as_cluster')
dev.off()

