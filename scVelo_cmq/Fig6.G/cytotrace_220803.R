.libPaths("/media/user/sdf/R.lib/library/wj/sc")
library(Seurat)
library(ggplot2)
library(RColorBrewer)
ss_obj<-readRDS('/media/user/sdd/zhaolabguest/wsy/sc_as/lawlor.1050/diff/endo_diff_cluster/data/law_c1to4.rds')
name<-as.data.frame( ss_obj@assays$RNA@counts)
write.table(name,'/media/user/sdh/cmq_wj/others/c1_4_ND_T2D.counts2cytotarce.txt')



cyto<-read.csv('/media/user/sdh/cmq_wj/others/CytoTRACE_results_c1_4.csv')
cyto_1<-as.data.frame(cyto[,2]) 
rownames(cyto_1) <-cyto[,1]
colnames(cyto_1)<-'CytoTRACE'
ss_obj<-AddMetaData(ss_obj, cyto_1, col.name = NULL)
# Idents(ss_obj)<-'new_as_cluster'
# Idents(ss_obj)<-'celltype'
FeaturePlot(ss_obj,features = 'MAFA',label = T,reduction ='new_as_cluster')

# saveRDS(ss_obj,'/media/user/sdh/cmq_wj/others/CYTOTRACE_C1_4.RDS')

AS<- readRDS('/media/user/sdh/cmq_wj/others/as_ND_T2D.RDS') 
cell = AS@meta.data[AS@meta.data$new_as_cluster %in% c(1:4),]
AS1_4<-subset(AS,cells= rownames(cell) )
AS1_4<-AddMetaData(AS1_4, cyto_1, col.name = NULL)
Idents(AS1_4)<-'new_as_cluster'
FeaturePlot(AS1_4,features = 'CytoTRACE',label = T,reduction ='as',
            cols=colorRampPalette(brewer.pal(11, "RdYlBu")[c(11:1)])(100))



ggsave('/media/user/sdh/cmq_wj/others/c1_4_cytotrace_as.pdf',width = 3.6,height = 3)

saveRDS(AS1_4,'/media/user/sdh/cmq_wj/others/c1_4_cytotrace_as.rds')
