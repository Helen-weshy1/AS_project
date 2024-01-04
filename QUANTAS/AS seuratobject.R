library(Seurat)
library(ggplot2)
library(data.table)
library(dplyr)
library(reticulate)
use_python("/python",required = T)
py_module_available("umap")

law_psi_junc=read.csv('/sc_as/lawlor.1050/raw_dataset/psi_junc_sc_data.csv')
dim(law_psi_junc)
a=read.table('/sc-AS/diab_pancr/Lawlor.Genom.Res.2016/count_AS_m4/SRR3617192/cass.count.txt', 
             header = F,sep = '\t')
rownames(law_psi_junc)=a$V4
law_psi_junc[1:4,1:4]

# splicing seuratobject
fil_law_psi=law_psi_junc[apply(law_psi_junc,1,function(x) (1133-sum(is.na(x)))>1133*0.3),]
fil_law_psi=as.matrix(fil_law_psi)
na100_law_psi=ifelse(is.na(fil_law_psi),100,fil_law_psi)
na100_law_psi[1:4,1:4]

law_sc_gene_all=readRDS('/sc_as/lawlor.1050/gene_sce/law_sc_gene_all.rds')
head(colnames(law_sc_gene_all))
na100_law_psi=as.data.frame(na100_law_psi)
count=na100_law_psi[,colnames(na100_law_psi)%in%colnames(law_sc_gene_all)]
count[1:4,1:4]
dim(law_sc_gene_all@meta.data)

law_na100 <- CreateSeuratObject(counts = count, 
                                meta.data =law_sc_gene_all@meta.data,
                                min.cells = 0, 
                                min.features =0,
                                project = "splicing")
law_na100
law_na100=FindVariableFeatures(law_na100,selection.method = "vst",nfeatures = 4616)
####### don't scale
count=as.matrix(law_na100@assays$RNA@counts)
law_na100@assays$RNA@scale.data=count

law_na100=RunPCA(law_na100)
law_na100 <- JackStraw(law_na100,num.replicate = 100)
law_na100<- ScoreJackStraw(law_na100, dims = 1:20)
JackStrawPlot(law_na100, dims = 1:15)
ElbowPlot(law_na100)
ggsave(filename = "spl_all_sce/08_A_ElbowPlot.pdf",width = 10.15,height = 6.29)

law_na100 <- FindNeighbors(law_na100, reduction = 'pca',dims = 1:5) 
law_na100 <- FindClusters(law_na100, resolution = seq(0,2,.1))
table(Idents(law_na100)) 

Idents(law_na100)='RNA_snn_res.2'
set.seed(1234)
for(i in 2:20){
  pdf.path <- paste0("spl_all_sce/umap/A_UMAP_1.1_dim", i, ".pdf")
  law_na100 <- RunUMAP(object = law_na100, dims = 1:i)
  DimPlot(object = law_na100, reduction = "umap", pt.size = 4, label = T,
          cols = c(
            "#bebebe", "#ce960d", "#ff3030", "#00ffff",
            "#4876ff", "#99fa99", "#cebe71", "#ff6ab4",
            "#64b8fe","#BFEFFF", "#8DB6CD", "#000080", "#a020ef",
            "#ffc0cb", "#00ff7f", "#ff1494", "#00688b"))
  ggsave(filename = pdf.path, width = 7.74, height = 7)
}
for(i in 2:20){
  pdf.path <- paste0("spl_all_sce/tsne/A_TSNE_1.1_dim", i, ".pdf")
  law_na100 <- RunTSNE(object = law_na100, dims = 1:i)
  DimPlot(object = law_na100, reduction = "tsne", pt.size = 4, label = T,
          cols = c(
            "#bebebe", "#ce960d", "#ff3030", "#00ffff",
            "#4876ff", "#99fa99", "#cebe71", "#ff6ab4",
            "#64b8fe","#BFEFFF", "#8DB6CD", "#000080", "#a020ef",
            "#ffc0cb", "#00ff7f", "#ff1494", "#00688b"))
  ggsave(filename = pdf.path, width = 7.74, height = 7)
}

law_na100=RunTSNE(law_na100,reduction = 'pca',dims = 1:4)
DimPlot(object = law_na100, reduction = "tsne", pt.size = 4, label = T,
        #group.by = 'celltype',
        cols = c(
          "#ce960d", "#ff3030", "#00ffff",
          "#4876ff", "#99fa99", "#cebe71", "#ff6ab4",
          "#64b8fe","#BFEFFF", "#8DB6CD", "#000080", "#a020ef",
          "#ffc0cb", "#00ff7f", "#ff1494", "#00688b"))
DimPlot(object = law_na100, reduction = "tsne", pt.size = 5, label = T,
        group.by = 'celltype',
        cols = c("#8C564BFF" ,"#1F77B4FF", "#FF7F0EFF", "#D62728FF",
                 "#E377C2FF","#9467BDFF" , "#2CA02CFF", "#a020ef", "#BCBD22FF" ))
ggsave('all_as_ct_tsne.pdf', width = 8, height = 6.5)
DimPlot(object = law_na100, reduction = "tsne", pt.size = 5, label = T,split.by = 'Status',
        group.by = 'celltype',
        cols = c("#8C564BFF" ,"#1F77B4FF", "#FF7F0EFF", "#D62728FF",
                 "#E377C2FF","#9467BDFF" , "#2CA02CFF", "#a020ef", "#BCBD22FF" ))
ggsave('all_as_ct_tsne_split.pdf', width = 10.5, height = 5)

DimPlot(object = law_na100, reduction = "tsne", pt.size = 5, label = T,
        cols = c(
          "#00688b", "#ff3030", "#00ffff","#4876ff", 
          "#ff1494","#00ff7f", "#ce960d",
          "#000080", "#a020ef",
          "#64b8fe","#BFEFFF", "#8DB6CD",
          "#ffc0cb", 
          "#99fa99", "#cebe71","#ff6ab4"))
ggsave('all_as_tsne.pdf', width = 7.5, height = 6.5)
DimPlot(object = law_na100, reduction = "tsne", pt.size = 5, label = T,split.by = 'Status',
        cols = c(
          "#00688b", "#ff3030", "#00ffff","#4876ff", 
          "#ff1494","#00ff7f", "#ce960d",
          "#000080", "#a020ef",
          "#64b8fe","#BFEFFF", "#8DB6CD",
          "#ffc0cb", 
          "#99fa99", "#cebe71","#ff6ab4"))
ggsave('all_as_tsne_split.pdf', width = 10, height = 5)
saveRDS(law_na100, file = "law_allcell_splicing_na100.rds")
