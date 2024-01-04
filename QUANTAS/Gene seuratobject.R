
library(Seurat)
library(ggplot2)
# library(clusterProfiler)
library(data.table)
library(dplyr)


# meta info
meta_info=read.table('SraRunTable.txt',sep = ',',header = T)
status_info=readxl::read_xlsx('patient_info.xlsx')
repeat_info=read.table('83repeat.csv',sep = ',',header = T)
colnames(repeat_info)
repeat_info=repeat_info[,1:5]
head(meta_info)
colnames(meta_info)
meta_info_s=meta_info[,c(1,2,27,29,30)]
colnames(meta_info_s)
colnames(repeat_info)
colnames(repeat_info)[1]="Sample.Name"
meta_info_s=left_join(meta_info_s,repeat_info[,c(1,5)],by='Sample.Name')

tmp=as.data.frame(unique(status_info[,c('Age','Status')]))
tmp
meta_info_s$Status=ifelse(meta_info_s$AGE==42|
                            meta_info_s$AGE==51|
                            meta_info_s$AGE==55,'T2D','ND')
table(meta_info_s$AGE,meta_info_s$Status)
table(meta_info_s$Status)
write.csv(meta_info_s,file = 'meta_dataset/1157_meta.csv',row.names = F,quote = F)

# gene seurat object
meta_info_s=read.csv('meta_dataset/1157_meta.csv',header = T)
View(meta_info_s)
table(meta_info_s$repeat.sample.)
lawlor_meta=na.omit(meta_info_s[meta_info_s$repeat.sample.=='NO',])
table(lawlor_meta$Status)
gene_sc_data=read.table('raw_dataset/gene_sc_data.csv',sep = ',',header = T)
head(sort(table(gene_sc_data$gene_symbol),decreasing =T ))
grep('^MT-',gene_sc_data$gene_symbol)#none
grep('No match found',gene_sc_data$gene_symbol)
grep('DUX4',gene_sc_data$gene_symbol)#remove duplicate
uniq=gene_sc_data[-c(17401,21757,22616,21850),]
rownames(uniq)=uniq$gene_symbol
uniq=uniq[,colnames(uniq)%in%lawlor_meta$Run]
genecellinfo=subset(lawlor_meta,lawlor_meta$Run%in%colnames(uniq))
View(genecellinfo)
rownames(genecellinfo)=genecellinfo$Run

law_sc_gene_all=CreateSeuratObject(counts = uniq,
                                   meta.data =genecellinfo ,
                                   project = 'lawlor',
                                   min.cells = 0,
                                   min.features = 0)
VlnPlot(object = law_sc_gene_all,features = c("nFeature_RNA", "nCount_RNA"), 
        pt.size = 0.5,group.by ='orig.ident')
ggsave(filename = "gene_sce/01_A_FeatureScatter.pdf",width = 5,height = 4.15)
law_sc_gene_all
law_sc_gene_all <- subset(law_sc_gene_all, subset = nFeature_RNA > 2500&nFeature_RNA<10000 )


law_sc_gene_all<-NormalizeData(law_sc_gene_all, normalization.method = "LogNormalize", scale.factor = 10000)
law_sc_gene_all <- FindVariableFeatures(law_sc_gene_all, selection.method = "vst", nfeatures = 2500)
top20 <- head(VariableFeatures(law_sc_gene_all), 20)
plot1 <- VariableFeaturePlot(law_sc_gene_all)
plot2 <- LabelPoints(plot = plot1, points = top20, repel = TRUE)
CombinePlots(plots = list(plot1, plot2))
ggsave(filename = "gene_sce/03_A_FeatureSelection.pdf",width = 15,height = 8)
all.genes<-rownames(law_sc_gene_all)
law_sc_gene_all<-ScaleData(law_sc_gene_all,features = all.genes)
law_sc_gene_all<-RunPCA(law_sc_gene_all, features = VariableFeatures(object = law_sc_gene_all))
law_sc_gene_all <- JackStraw(law_sc_gene_all, num.replicate = 100)
law_sc_gene_all <- ScoreJackStraw(law_sc_gene_all, dims = 1:20)
JackStrawPlot(law_sc_gene_all, dims = 1:15)
ggsave(filename = "gene_sce/07_A_JackStrawPlot.pdf",width = 10.15,height = 5)
ElbowPlot(law_sc_gene_all)
ggsave(filename = "gene_sce/08_A_ElbowPlot.pdf",width = 10.15,height = 6.29)

law_sc_gene_all <- FindNeighbors(law_sc_gene_all, dims = 1:10)
law_sc_gene_all <- FindClusters(law_sc_gene_all, resolution = 3)
library(clustree)
clustree(law_sc_gene_all@meta.data,prefix = 'RNA_snn_res.')
ggsave(filename = "clustree.pdf",width = 10.15,height = 6.29)

Idents(law_sc_gene_all)='RNA_snn_res.3'
#TSNE
for(i in 2:15){
  pdf.path <- paste0("gene_sce/tsne/A_TSNE_1_dim", i, ".pdf")
  law_sc_gene_all <- RunTSNE(object = law_sc_gene_all, dims = 1:i)
  DimPlot(object = law_sc_gene_all, reduction = "tsne", pt.size = 4, label = T,
          cols = c( 
            "#bebebe", "#ce960d", "#ff3030", "#00ffff", 
            "#4876ff", "#99fa99", "#cebe71", "#ff6ab4", 
            "#64b8fe","#BFEFFF", "#8DB6CD", "#000080", "#a020ef", 
            "#ffc0cb", "#00ff7f", "#ff1494", "#00688b"))
  ggsave(filename = pdf.path, width = 7.74, height = 7)
}

law_sc_gene_all <- RunTSNE(object = law_sc_gene_all, dims = 1:15, do.fast = TRUE,do.label=T)
DimPlot(object = law_sc_gene_all, reduction = "tsne", pt.size = 4, label = T,group.by = 'seurat_clusters',
        # cols = c( 
        #   "#1F77B4FF", "#FF7F0EFF", "#2CA02CFF" ,"#D62728FF", "#9467BDFF",
        #   "#AEC7E8FF", "#FFBB78FF", "#98DF8AFF", "#FF9896FF",
        #   "#8C564BFF" ,"#E377C2FF", "#C49C94FF", "#BCBD22FF", "#17BECFFF","#bebebe")
)
ggsave(filename = 'gene_sce/all_cluster_TSNE_dim15.pdf', width = 7.74, height = 7)

FeaturePlot(law_sc_gene_all, reduction = 'tsne', pt.size = 2,
            features = c( "GCG","INS", "SST",'PPY', "KRT19", "PRSS1", "COL1A1",'GHRL','CD86','TPSAB1','COL1A2','PLVAP'))
ggsave(filename = "gene_sce/all_FeaturePlot_tsne.pdf",width = 30,height = 20)
new.cluster.ids <- c("Endocrine",
                     "Endocrine", "Endocrine", "Ductal", "Endocrine","Endocrine",
                     "Endocrine", "Endocrine",'Endocrine','Endocrine','Stellate',
                     "Endocrine", "Endocrine", "Acinar", "Endocrine", "Acinar", 
                     "Endothelial", "Endocrine")
names(new.cluster.ids) <- levels(law_sc_gene_all)
law_sc_gene_all <- RenameIdents(law_sc_gene_all, new.cluster.ids)
law_sc_gene_all$seurat_clusters=law_sc_gene_all$RNA_snn_res.1
table(Idents(law_sc_gene_all))
law_sc_gene_all$celltype=Idents(law_sc_gene_all)
DimPlot(object = law_sc_gene_all, reduction = "tsne", pt.size = 4, label = T,
        cols = c( 
          "#1F77B4FF", "#FF7F0EFF", "#2CA02CFF" ,"#D62728FF", "#9467BDFF",
          "#AEC7E8FF", "#FFBB78FF", "#98DF8AFF", "#FF9896FF",
          "#8C564BFF" ,"#E377C2FF", "#C49C94FF", "#BCBD22FF", "#17BECFFF"))
ggsave(filename = 'gene_sce/all_celltype_TSNE_dim15.pdf', width = 7.74, height = 7)
saveRDS(law_sc_gene_all, file = "gene_sce/law_sc_gene_all.rds")


# endocrine
law_sc_gene_endo=subset(law_sc_gene_all,cells=WhichCells(law_sc_gene_all,idents = 'Endocrine'))
law_sc_gene_endo=FindVariableFeatures(law_sc_gene_endo,selection.method = "vst", nfeatures = 2500)
law_sc_gene_endo=RunPCA(law_sc_gene_endo,features = VariableFeatures(object = law_sc_gene_endo))
law_sc_gene_endo <- FindNeighbors(law_sc_gene_endo, dims = 1:15) 
law_sc_gene_endo <- FindClusters(law_sc_gene_endo, resolution = 2)
Idents(law_sc_gene_endo)='RNA_snn_res.2'
set.seed(1234)
for(i in 4:20){
  pdf.path <- paste0("gene_sce/gene_endo_sce/tsne/A_TSNE_endocrine_dim", i, ".pdf")
  law_sc_gene_endo <- RunTSNE(object = law_sc_gene_endo, dims = 1:i)
  DimPlot(object = law_sc_gene_endo, reduction = "tsne", pt.size = 2, label = T,
          cols = c("#1F77B4FF", "#FF7F0EFF", "#2CA02CFF" ,"#D62728FF", "#9467BDFF",
                   "#AEC7E8FF", "#FFBB78FF", "#98DF8AFF", "#FF9896FF",
                   "#8C564BFF" ,"#E377C2FF", "#C49C94FF", "#BCBD22FF", "#17BECFFF"))
  ggsave(filename = pdf.path, width = 7.74, height = 7)
}
law_sc_gene_endo=RunTSNE(object = law_sc_gene_endo,dims = 1:14)
DimPlot(object = law_sc_gene_endo, reduction = "tsne", pt.size = 4, label = T,
        #group.by = 'seurat_clusters',
        # cols = c(
        #   "#1F77B4FF", "#FF7F0EFF", "#2CA02CFF" ,"#D62728FF", "#9467BDFF",
        #   "#AEC7E8FF", "#FFBB78FF", "#98DF8AFF", "#FF9896FF",
        #   "#8C564BFF" ,"#E377C2FF", "#C49C94FF", "#BCBD22FF","#bebebe" ,"#17BECFFF")
)
ggsave(filename = 'gene_sce/gene_endo_sce/endo_cluster_TSNE_dim14_res2.pdf', width = 7.74, height = 7)

FeaturePlot(law_sc_gene_endo, reduction = 'tsne', pt.size = 4,
            features = c("GCG", "ARX",
                         "INS",'NPTX2',
                         'PPY',
                         "SST"))
ggsave(filename = "gene_sce/gene_endo_sce/endo_FeaturePlot_tsne.pdf",width = 15,height = 15)
table(Idents(law_sc_gene_endo))
new.cluster.ids2=c('Beta',
                   'Alpha','Alpha','Multiple','Beta','Beta',
                   'Alpha','Beta','PP','Delta','Alpha')
names(new.cluster.ids2) <- levels(law_sc_gene_endo)
law_sc_gene_endo <- RenameIdents(law_sc_gene_endo, new.cluster.ids2)
mylevel=c('Alpha','Beta','Multiple','Delta','PP')
Idents(law_sc_gene_endo)=factor(Idents(law_sc_gene_endo),levels = mylevel)
law_sc_gene_endo$celltype=law_sc_gene_endo@active.ident
law_sc_gene_endo$seurat_clusters=law_sc_gene_endo$RNA_snn_res.2
saveRDS(law_sc_gene_endo, file = "gene_sce/gene_endo_sce/law_sc_gene_endo.rds")

Idents(law_sc_gene_endo)='RNA_snn_res.2'
table(Idents(law_sc_gene_endo))
law_sc_gene_c7=subset(law_sc_gene_endo,cells=WhichCells(law_sc_gene_endo,idents = '7'))
law_sc_gene_c7=FindVariableFeatures(law_sc_gene_c7,selection.method = "vst", nfeatures = 2500)
law_sc_gene_c7=RunPCA(law_sc_gene_c7,features = VariableFeatures(object = law_sc_gene_c7))
law_sc_gene_c7 <- FindNeighbors(law_sc_gene_c7, dims = 1:10) 
law_sc_gene_c7 <- FindClusters(law_sc_gene_c7, resolution =1)
Idents(law_sc_gene_c7)='RNA_snn_res.1'
FeaturePlot(law_sc_gene_c7, reduction = 'tsne', pt.size = 4,
            features = c("GCG","INS",'PPY',"SST"))
ggsave(filename = "gene_sce/gene_endo_sce/endo_c7_FeaturePlot_tsne.pdf",width = 11,height = 10)
DimPlot(object = law_sc_gene_c7, reduction = "tsne", pt.size = 4, label = T,
        #group.by = 'seurat_clusters',
        cols = c( 
          "#1F77B4FF", "#FF7F0EFF", "#bebebe", "#2CA02CFF" ,"#D62728FF", "#9467BDFF",
          "#AEC7E8FF", "#FFBB78FF", "#98DF8AFF", "#FF9896FF",
          "#8C564BFF" ,"#E377C2FF", "#C49C94FF", "#BCBD22FF", "#17BECFFF"))
new.cluster.ids3=c('Beta','Alpha')
names(new.cluster.ids3) <- levels(law_sc_gene_c7)
law_sc_gene_c7 <- RenameIdents(law_sc_gene_c7, new.cluster.ids3)
mylevel=c('Alpha','Beta')
Idents(law_sc_gene_c7)=factor(Idents(law_sc_gene_c7),levels = mylevel)
law_sc_gene_c7$celltype=law_sc_gene_c7@active.ident

table(law_sc_gene_c7$celltype)
table(law_sc_gene_endo$celltype)
law_sc_gene_c7$celltype=as.character(law_sc_gene_c7$celltype)
law_sc_gene_endo$celltype=as.character(law_sc_gene_endo$celltype)
law_sc_gene_endo$celltype[rownames(law_sc_gene_c7@meta.data)]<-law_sc_gene_c7@meta.data$celltype
table(law_sc_gene_endo$celltype)

mylevel=c('Alpha','Beta','Multiple','Delta','PP')
Idents(law_sc_gene_endo)='celltype'
Idents(law_sc_gene_endo)=factor(Idents(law_sc_gene_endo),levels = mylevel)
law_sc_gene_endo$celltype=Idents(law_sc_gene_endo)
DimPlot(object = law_sc_gene_endo, reduction = "tsne", pt.size = 4, label = T,
        #group.by = 'seurat_clusters',
        cols = c( 
          "#1F77B4FF", "#FF7F0EFF", "#bebebe", "#2CA02CFF" ,"#D62728FF", "#9467BDFF",
          "#AEC7E8FF", "#FFBB78FF", "#98DF8AFF", "#FF9896FF",
          "#8C564BFF" ,"#E377C2FF", "#C49C94FF", "#BCBD22FF", "#17BECFFF"))
ggsave(filename = 'gene_sce/gene_endo_sce/endo_celltype_TSNE_dim14_res2.pdf', width = 7.74, height = 7)
saveRDS(law_sc_gene_endo, file = "gene_sce/gene_endo_sce/law_sc_gene_endo.rds")


##
table(Idents(law_sc_gene_endo))
law_gene_beta=subset(law_sc_gene_endo,cells = WhichCells(law_sc_gene_endo,idents = 'Beta'))
Idents(law_gene_beta)='Status'
table(Idents(law_gene_beta))
VlnPlot(object = law_gene_beta,features=c('STX1A','DLK1'))
ggsave(filename = 'gene_sce/gene_endo_sce/STX1A_DLK1_vln.pdf', width = 7.74, height = 7)
Idents(law_sc_gene_endo)='celltype'
mylevel=c('Alpha','Beta','Multiple','Delta','PP')
Idents(law_sc_gene_endo)=factor(Idents(law_sc_gene_endo),levels = mylevel)
VlnPlot(object = law_sc_gene_endo,features=c('STX1A','DLK1',
                                             'GDA', 'CD36', 'RCOR1', 'LAPTM4B'),split.by = 'Status')
ggsave(filename = 'gene_sce/gene_endo_sce/STX1A_DLK1_vln.pdf', width = 7.74, height = 7)


# endo info mapping return to all cell
table(law_sc_gene_all$celltype)
table(law_sc_gene_endo$celltype)
law_sc_gene_all$celltype=as.character(law_sc_gene_all$celltype)
law_sc_gene_endo$celltype=as.character(law_sc_gene_endo$celltype)
law_sc_gene_all$celltype[rownames(law_sc_gene_endo@meta.data)]<-law_sc_gene_endo@meta.data$celltype
table(law_sc_gene_all$celltype)
mylevel=c('Alpha','Beta','Multiple','Delta','PP','Acinar','Ductal','Stellate','Endothelial')
Idents(law_sc_gene_all)=factor(law_sc_gene_all$celltype,levels = mylevel)
table(Idents(law_sc_gene_all))
DimPlot(law_sc_gene_all,pt.size = 2,label = T,repel=T,
        cols = c( 
          "#1F77B4FF", "#FF7F0EFF",  "#2CA02CFF","#D62728FF", "#a020ef",
          "#8C564BFF" ,"#E377C2FF", "#BCBD22FF" ,"#9467BDFF", "#bebebe" ))+labs(title = 'T2D and ND')
ggsave(filename = "gene_sce/all_detail_celltype_tsne.pdf",width = 8, height =6.5)
DimPlot(law_sc_gene_all,pt.size = 2,label = T,repel=T,split.by = 'Status',
        cols = c( 
          "#1F77B4FF", "#FF7F0EFF",  "#2CA02CFF","#D62728FF", "#a020ef",
          "#8C564BFF" ,"#E377C2FF", "#BCBD22FF" ,"#9467BDFF", "#bebebe" ))
ggsave(filename = "gene_sce/all_detail_celltype_tsne_split.pdf",width = 10, height =5)

saveRDS(law_sc_gene_all, file = "gene_sce/law_sc_gene_all.rds")

