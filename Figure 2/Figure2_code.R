# Fig2
library(Seurat)
library(ggplot2)
library(data.table)
library(dplyr)
library(reticulate)


# Endocrine cell AS seurat object from the Lawlor dataset -----------------------------------------

table(law_AS_sce_all$celltype)
Idents(law_AS_sce_all)='celltype'
endo_as_sce=subset(law_AS_sce_all,
                   cells = WhichCells(law_AS_sce_all,idents = c('Alpha',
                                                           'Beta',
                                                           'Delta',
                                                           'PP',
                                                           'Multiple')))
dim(endo_as_sce)
endo_as_sce=FindVariableFeatures(endo_as_sce,selection.method = "vst")
endo_as_sce=RunPCA(endo_as_sce)
endo_as_sce <- JackStraw(endo_as_sce,num.replicate = 100)
endo_as_sce<- ScoreJackStraw(endo_as_sce, dims = 1:20)
ElbowPlot(endo_as_sce)
endo_as_sce <- FindNeighbors(endo_as_sce, reduction = 'pca',dims = 1:10) 
endo_as_sce <- FindClusters(endo_as_sce, resolution =2)
endo_as_sce=RunTSNE(endo_as_sce,reduction = 'pca',dims = 1:5)
new.cluster.ids <- c('2','3','5','7','6','1','9','4','10','8')
names(new.cluster.ids) <- levels(endo_as_sce)
endo_as_sce <- RenameIdents(endo_as_sce, new.cluster.ids)
mylevel=c('1','2','3','4','5','6','7','8','9','10')
Idents(endo_as_sce)=factor(Idents(endo_as_sce),levels = mylevel)
endo_as_sce$new_as_cluster=endo_as_sce@active.ident
endogroup=endo_as_sce@meta.data


# * panel A ---------------------------------------------------------------
DimPlot(object = endo_as_sce, reduction = "tsne", pt.size = 3, label = T,
        cols = c("#8C564BFF" ,"#E377C2FF", "#C49C94FF", "#BCBD22FF", "#17BECFFF",
                 "#1F77B4FF", "#FF7F0EFF", "#2CA02CFF" ,"#D62728FF", "#9467BDFF",
                 "#AEC7E8FF", "#FFBB78FF", "#98DF8AFF", "#FF9896FF" ))
ggsave('endo_as_tsne.pdf', width = 8, height = 6.5)


# * panel B ---------------------------------------------------------------
DimPlot(object = endo_as_sce, reduction = "tsne", pt.size = 3, label = T,
        group.by = 'celltype',
        cols = c("#1F77B4FF", "#FF7F0EFF", "#D62728FF",
                 "#2CA02CFF", "#a020ef" ))
ggsave('endo_as_ct_tsne.pdf', width = 8, height = 6.5)


# * panel C ---------------------------------------------------------------
overlap1=as.data.frame(table(endogroup$new_as_cluster,
                             endogroup$celltype))
colnames(overlap1)=c('cluster','celltype','Freq')
x=vector(mode = 'numeric')
for (i in 1:10){
  j=unique(overlap1$cluster)[i]
  a=subset(overlap1,cluster==j)
  prop.table(a$Freq)->sumnum
  x=append(x,sumnum)
}
overlapdata=data.frame(celltype=c('Alpha','Beta', 'Delta','Multiple', 'PP'),
                       v1=x[1:5],
                       v2=x[6:10],
                       v3=x[11:15],
                       v4=x[16:20],
                       v5=x[21:25],
                       v6=x[26:30],
                       v7=x[31:35],
                       v8=x[36:40],
                       v9=x[41:45],
                       v10=x[46:50])
rownames(overlapdata)=overlapdata$celltype
overlapdata=overlapdata[,-1]
colnames(overlapdata)=1:10
apply(overlapdata,2,function(x) sum(x))
overlapdata=t(overlapdata)%>%as.data.frame()
library(pheatmap)
pheatmap(overlapdata)
overlapnew=data.frame(Beta=overlapdata$Beta,
                      Multiple=overlapdata$Multiple,
                      Alpha=overlapdata$Alpha,
                      Delta=overlapdata$Delta,
                      PP=overlapdata$PP)
pdf('endo_ct_cluster_overlap_heatmap.pdf',width = 5,height = 6)
pheatmap(overlapnew, border=F,
         cluster_rows = F,cluster_cols=F,
         color = colorRampPalette(c("white",'#fbcdd5', 'firebrick3','#960001'))(100),
         display_numbers = F,border_color = "white",
         angle_col = '45')
dev.off()


# * panel D ---------------------------------------------------------------
# Using QUANTAS detects differential spliced exons
up_PSI=character()
for (i in 1:10){
  cldata=read.table(paste0('diff_cluster/endo_c',i,'_dI0.2_ratio0.3.diff.txt'),
                    sep = ',',header = T)
  cldata=cldata[order(cldata$dI_g1_vs_g2,decreasing = T),]
  up_PSI=c(up_PSI,as.character(cldata$name[c(1:10)]))
  }
up_PSI_d=as.data.frame(unique(up_PSI))
colnames(up_PSI_d)='name'
up_PSI_d=left_join(up_PSI_d,law_endo_cluster_psi[,c(1:7,11:20)])
up_PSI_d_p=reshape2::melt(up_PSI_d[,c(1,2,8:17)],id=c('name','gene'))
colnames(up_PSI_d_p)[4]='PSI'

up_FDR=left_join(up_PSI_d,law_endo_cluster_FDR[,c(5,15:24)],by='name')
up_FDR_p=reshape2::melt(up_FDR[,c(1,2,18:27)],id=c('name','gene'))
up_PSI_d_p$FDR=up_FDR_p$value
up_PSI_d_p_f=up_PSI_d_p
up_PSI_d_p_f_FDR=ifelse(up_PSI_d_p$FDR>0.01,NA,up_PSI_d_p$FDR)
up_PSI_d_p_f_PSI=ifelse(up_PSI_d_p$FDR>0.01,NA,up_PSI_d_p$PSI)

up_PSI_d_p_f$FDR=up_PSI_d_p_f_FDR
up_PSI_d_p_f$PSI=up_PSI_d_p_f_PSI
ggplot(up_PSI_d_p_f,aes(x=variable,y=factor(up_PSI_d_p$name,levels =rev(unique(up_PSI)) )))+
  geom_point(aes(size=PSI,color=FDR))+
  theme_bw()+scale_size_continuous(range = c(0,7))+
  theme(panel.grid = element_blank(),
        axis.text.x=element_text(angle=90,hjust = 1,vjust=0.5))+
  #scale_color_gradient(low="#F08080",high="#FF0000")+
  scale_color_gradientn(colors=c(
    #mid = "#E61C15",
    '#FF8247',
    #'orange',
    "white"))+
  labs(x=NULL,y=NULL)
ggsave('endo_up_marker10_retio0.3_dotplot.pdf',width = 7,height = 18)


# * panel E ---------------------------------------------------------------
law_sc_gene_endo=subset(law_sc_gene_all,
                   cells = WhichCells(law_sc_gene_all,idents = c('Alpha',
                                                                'Beta',
                                                                'Delta',
                                                                'PP',
                                                                'Multiple')))
VlnPlot(law_sc_gene_endo,group.by = 'celltype',features =c('SEC13','C7odf44'))

# * panel F ---------------------------------------------------------------
VlnPlot(law_sc_gene_endo,group.by = 'celltype',features =c('KARS','SBDSP1'))

