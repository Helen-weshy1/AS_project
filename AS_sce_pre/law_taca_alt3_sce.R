library(Seurat)
library(ggplot2)
# library(clusterProfiler)
library(data.table)
library(dplyr)
library(reticulate)

law_psi=read.table('/media/user/sdg/WShi/AS/lawlor/merge/all_as_PSI.txt',
                        sep = '\t',header = T)
dim(law_psi)
law_psi[1:4,1:10]
rownames(law_psi)=law_psi$V4
law_psi[1:4,1:4]

################################# taca #########################################
# splicing seuratobject
law_psi_taca=law_psi[law_psi$V7=='taca',]
law_psi_taca=law_psi_taca[,-c(1:7)]


law_sc_gene_all=readRDS('/media/user/sdd/zhaolabguest/wsy/sc_as/lawlor.1050/gene_sce/law_sc_gene_all.rds')
law_sc_gene_all

law_psi_taca=law_psi_taca[,colnames(law_psi_taca)%in%colnames(law_sc_gene_all)]
fil_law_psi=law_psi_taca[apply(law_psi_taca,1,function(x) (972-sum(is.na(x)))>972*0.3),]
fil_law_psi=as.matrix(fil_law_psi)
na100_law_psi=ifelse(is.na(fil_law_psi),100,fil_law_psi)
na100_law_psi[1:4,1:4]

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
law_na100=FindVariableFeatures(law_na100,selection.method = "vst",nfeatures = 3465)
####### don't scale
count=as.matrix(law_na100@assays$RNA@counts)
law_na100@assays$RNA@scale.data=count

law_na100=RunPCA(law_na100)
law_na100 <- JackStraw(law_na100,num.replicate = 100)
law_na100<- ScoreJackStraw(law_na100, dims = 1:20)
JackStrawPlot(law_na100, dims = 1:15)
ElbowPlot(law_na100)

law_na100 <- FindNeighbors(law_na100, reduction = 'pca',dims = 1:4) 
law_na100 <- FindClusters(law_na100, resolution = 1)
table(Idents(law_na100)) 

set.seed(1234)

law_na100=RunTSNE(law_na100,reduction = 'pca',dims = 1:4)
p1=DimPlot(object = law_na100, reduction = "tsne", pt.size = 1, label = T,
        #group.by = 'celltype',
        cols = c(
          "#ce960d", "#ff3030", "#00ffff",
          "#4876ff", "#99fa99", "#cebe71", "#ff6ab4",
          "#64b8fe","#BFEFFF", "#8DB6CD", "#000080", "#a020ef",
          "#ffc0cb", "#00ff7f", "#ff1494", "#00688b"))
p2=DimPlot(object = law_na100, reduction = "tsne", pt.size = 1, label = T,
        group.by = 'celltype',
        cols = c("#8C564BFF" ,"#1F77B4FF", "#FF7F0EFF", "#D62728FF",
                 "#E377C2FF","#9467BDFF" , "#2CA02CFF", "#a020ef", "#BCBD22FF" ))
p3=DimPlot(object = law_na100, reduction = "tsne", pt.size = 1, label = T,split.by = 'Status',
        group.by = 'celltype',
        cols = c("#8C564BFF" ,"#1F77B4FF", "#FF7F0EFF", "#D62728FF",
                 "#E377C2FF","#9467BDFF" , "#2CA02CFF", "#a020ef", "#BCBD22FF" ))

p1+p2+p3+plot_layout(widths = c(1,1,2))
ggsave('/media/user/sdd/zhaolabguest/wsy/sc_as/lawlor.1050/other_AS/law_taca_1.pdf',
       width = 16,height = 4)

#saveRDS(law_na100,file='/media/user/sdd/zhaolabguest/wsy/sc_as/lawlor.1050/other_AS/law_taca_as_sce.rds')


################################# alt3 #########################################
# splicing seuratobject
law_psi_alt3=law_psi[law_psi$V7=='alt3',]
law_psi_alt3=law_psi_alt3[,-c(1:7)]

law_psi_alt3=law_psi_alt3[,colnames(law_psi_alt3)%in%colnames(law_sc_gene_all)]
fil_law_psi=law_psi_alt3[apply(law_psi_alt3,1,function(x) (972-sum(is.na(x)))>972*0.3),]
fil_law_psi=as.matrix(fil_law_psi)
na100_law_psi=ifelse(is.na(fil_law_psi),100,fil_law_psi)
na100_law_psi[1:4,1:4]

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
law_na100=FindVariableFeatures(law_na100,selection.method = "vst",nfeatures = 2187)
####### don't scale
count=as.matrix(law_na100@assays$RNA@counts)
law_na100@assays$RNA@scale.data=count

law_na100=RunPCA(law_na100)
law_na100 <- JackStraw(law_na100,num.replicate = 100)
law_na100<- ScoreJackStraw(law_na100, dims = 1:20)
JackStrawPlot(law_na100, dims = 1:15)
ElbowPlot(law_na100)

law_na100 <- FindNeighbors(law_na100, reduction = 'pca',dims = 1:3) 
law_na100 <- FindClusters(law_na100, resolution = 1)
table(Idents(law_na100)) 

law_na100=RunTSNE(law_na100,reduction = 'pca',dims = 1:3)
p1=DimPlot(object = law_na100, reduction = "tsne", pt.size = 1, label = T,
           #group.by = 'celltype',
           cols = c(
             "#ce960d", "#ff3030", "#00ffff",
             "#4876ff", "#99fa99", "#cebe71", "#ff6ab4",
             "#64b8fe","#BFEFFF", "#8DB6CD", "#000080", "#a020ef",
             "#ffc0cb", "#00ff7f", "#ff1494", "#00688b"))
p2=DimPlot(object = law_na100, reduction = "tsne", pt.size = 1, label = T,
           group.by = 'celltype',
           cols = c("#8C564BFF" ,"#1F77B4FF", "#FF7F0EFF", "#D62728FF",
                    "#E377C2FF","#9467BDFF" , "#2CA02CFF", "#a020ef", "#BCBD22FF" ))
p3=DimPlot(object = law_na100, reduction = "tsne", pt.size = 1, label = T,split.by = 'Status',
           group.by = 'celltype',
           cols = c("#8C564BFF" ,"#1F77B4FF", "#FF7F0EFF", "#D62728FF",
                    "#E377C2FF","#9467BDFF" , "#2CA02CFF", "#a020ef", "#BCBD22FF" ))

p1+p2+p3+plot_layout(widths = c(1,1,2))
ggsave('/media/user/sdd/zhaolabguest/wsy/sc_as/lawlor.1050/other_AS/law_alt3_1.pdf',
       width = 16,height = 4)

#saveRDS(law_na100,file='/media/user/sdd/zhaolabguest/wsy/sc_as/lawlor.1050/other_AS/law_alt3_as_sce.rds')
