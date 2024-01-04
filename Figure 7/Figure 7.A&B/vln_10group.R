library(Seurat)
library(ggplot2)

colors<-c("#8C564BFF" ,"#8C564BFF" ,"#9467BDFF","#9467BDFF",
          "#E377C2FF", "#E377C2FF", 
          "#C49C94FF", "#C49C94FF",
          "#BCBD22FF","#BCBD22FF", 
          "#17BECFFF","#17BECFFF",
          "#1F77B4FF","#1F77B4FF",
          "#FF7F0EFF", "#FF7F0EFF",
          "#2CA02CFF" ,"#2CA02CFF" ,
          "#D62728FF","#D62728FF")
colors_2<-c("#8C564BFF" ,"#E377C2FF", "#C49C94FF", "#BCBD22FF", "#17BECFFF",
            "#1F77B4FF", "#FF7F0EFF", "#2CA02CFF" ,"#D62728FF", "#9467BDFF",
            "#8C564BFF" ,"#E377C2FF", "#C49C94FF", "#BCBD22FF", "#17BECFFF",
            "#1F77B4FF", "#FF7F0EFF", "#2CA02CFF" ,"#D62728FF", "#9467BDFF")
ss_obj<-readRDS('/as_ND_T2D.RDS')

ND_OBJ <-readRDS("/law_sc_gene_endo_nd.rds")
T2D_OBJ<-readRDS("/law_sc_gene_endo_t2d.rds")

Idents(ND_OBJ)<-'new_as_cluster'
Idents(T2D_OBJ)<-'new_as_cluster'
aa<-as.data.frame( T2D_OBJ@active.ident)
aa[,1]<- (as.numeric( aa$`T2D_OBJ@active.ident`) -1)
# aa<-factor(aa)

aa1<-as.data.frame( ND_OBJ@active.ident)
aa1[,1]<- (as.numeric( aa1$`ND_OBJ@active.ident`) -1)
# aa1<-as.factor(aa1)
ss_obj<-AddMetaData(ss_obj,c(paste0(aa$`T2D_OBJ@active.ident`,'_T2D'),
                             paste0(aa1$`ND_OBJ@active.ident`,'_ND')  ), col.name = 'vlnplot_group') #添加分组为metadata
Idents(ss_obj)<-'vlnplot_group'


ss_obj<-AddMetaData(ss_obj,c(paste0(T2D_OBJ@active.ident,'_T2D'),
                             paste0(ND_OBJ@active.ident,'_ND')  ), col.name = 'vlnplot_group') #添加分组为metadata
Idents(ss_obj)<-'vlnplot_group'

ss_obj<-AddMetaData(ss_obj,c(paste0('T2D_',T2D_OBJ@active.ident),
                             paste0('ND_',ND_OBJ@active.ident)  ), col.name = 'vlnplot_group2') #添加分组为metadata
Idents(ss_obj)<-'vlnplot_group2'
# saveRDS(ss_obj,'20GROUP_ND_T2D.rds')
ss_obj<-readRDS('20GROUP_ND_T2D.rds')





pdf('dediff_ma.pdf',width = 20,height = 15)
VlnPlot(ss_obj,features = c('INS','SLC2A2','MAFA','RFX6','PDX1','CHL1','GCK',
                            'PPARGC1A','ESRRG','MDH1','NEUROD1','CREB1','G6PC2',
                            'PFKFB2','PFKM','SIX2','SIX3','ENTPD3','GPD2','DNMT3A',
                            'SLC25A1','MTOR'),
        group.by = 'vlnplot_group',cols=rep(c("#add8e6","#ee5e49"),10),
        pt.size=0.1,ncol=4 )
dev.off()


pdf('transdiff_a.pdf',width = 15,height = 6)
VlnPlot(ss_obj,features = c('RFX6','GCG','TM4SF4',
                            'RGS4','ARX','SMARCA1'),
        group.by = 'vlnplot_group',cols=rep(c("#add8e6","#ee5e49"),10),
        pt.size=0.1,ncol=3 )
dev.off()

pdf('score_3.pdf',width = 10,height = 3)
VlnPlot(ss_obj,features = c("MURARO_PANCREAS_ALPHA_CELL.v7.5.11",
                            'MURARO_PANCREAS_BETA_CELL.v7.5.11'),
        group.by = 'vlnplot_group',cols=rep(c("#add8e6","#ee5e49"),10),
        pt.size=0,ncol=2 )+ 
 geom_boxplot(width=.2,col="black",fill="white")+  
 NoLegend()
dev.off()


######### 
 aa1<-c('TM4SF4','ARX','SMARCA1','MAFA','RFX6','CHL1','ENTPD3')
r=1
pdf('vlnplot_0806_1-4mean.pdf',width = 5,height = 3)
VlnPlot(ss_obj,features =aa1[r],
        group.by = 'vlnplot_group',cols=rep(c("#add8e6","#ee5e49"),10),
        pt.size=0.1,ncol=1 )+
  stat_summary(fun= mean, geom = "point",
               shape = 1, size = 2, color = "black")+NoLegend()
r=r+1
dev.off()

################
pdf('vlnplot_0806_5-9mean.pdf',width = 5,height = 3)
aa<-c('RFX6',
      'PFKFB2','CREB1','ESRRG')
i=1
VlnPlot(ss_obj,features =aa[i],
        group.by = 'vlnplot_group',cols=rep(c("#add8e6","#ee5e49"),10),
        pt.size=0.1,ncol=1 )+
  stat_summary(fun= mean, geom = "point",
               shape = 1, size = 2, color = "black")+NoLegend()
i=i+1
dev.off()


######
pdf('dediff_imma-0806median.pdf',width = 5,height = 3)
VlnPlot(ss_obj,features ='HES1' ,
        group.by = 'vlnplot_group',cols=rep(c("#add8e6","#ee5e49"),10),
        pt.size=0.1)+
  stat_summary(fun= median, geom = "point",
               shape = 1, size = 3, color = "black")+NoLegend()

VlnPlot(ss_obj,features ='CD81' ,
          group.by = 'vlnplot_group',cols=rep(c("#add8e6","#ee5e49"),10),
          pt.size=0.1)+
    stat_summary(fun= median, geom = "point",
                 shape = 1, size = 2, color = "black")+NoLegend()

  
dev.off()

pdf('dediff_imma-0806mean.pdf',width = 5,height = 3)
VlnPlot(ss_obj,features ='HES1' ,
        group.by = 'vlnplot_group',cols=rep(c("#add8e6","#ee5e49"),10),
        pt.size=0.1)+
  stat_summary(fun= mean, geom = "point",
               shape = 1, size = 3, color = "black")+NoLegend()
VlnPlot(ss_obj,features ='CD81' ,
        group.by = 'vlnplot_group',cols=rep(c("#add8e6","#ee5e49"),10),
        pt.size=0.1)+
  stat_summary(fun= mean, geom = "point",
               shape = 1, size = 2, color = "black")+NoLegend()

dev.off()







#############              score
# rm(list=ls())
# setwd('/media/user/sdh/cmq_wj/others')
source('/MAP.R')
library(Seurat)
library(ggplot2)
ss_obj
######333333   go_geneset
for (group in c('ND','T2D')){
  
  name1 <- c(rownames(ss_obj@meta.data)[grep(group,ss_obj$Status)])
  ss_obj_SUB <- subset(ss_obj,cells = name1)
  
tsne_result.ch<-data.frame(ss_obj_SUB@reductions[["as"]]@cell.embeddings); tsne_result.ch[1:8,]
path<-dir('/GO_GENE')
for (i1 in path){
  he<-read.table(paste0 ('/GO_GENE/',i1))
  
  gene <- list(he[1,2:length(he[1,])])
  ss_obj_SUB <- AddModuleScore(  #对每个模块的基因集打分 
    object = ss_obj_SUB,
    features = gene,
    ctrl = 100, #默认值是100
    name = strsplit(i1,'.gmt')[[1]][1])
}

if (F)
  {for (i in path){

  name<-paste0(strsplit(i,'.gmt')[[1]][1],'1')
  expr.ch<-data.frame(ss_obj_SUB@meta.data[[name]]);colnames(expr.ch)<-name
  rownames(expr.ch)<-rownames(tsne_result.ch)
  expr.ch<-as.data.frame(t(expr.ch))
  #embedding <- reducedDim(tsne.ch, "tSNE")
  
  
  pdf(paste0(group,"_IsoheightMap_",name,".pdf"), width=8, height=6)
  p <- list()
  for (genes in name) {
    p[[genes]] <- IsoheightMap(
      tsne_result.ch,
      genes, 
      expr.ch)
  }
  do.call(grid.arrange, c(p, ncol=1)) #(3,3) #(15,3*)
  dev.off()
  }}

}
ss_obj_T2D<-ss_obj_SUB
ss_obj_ND<-ss_obj_SUB

i=1
ND_SCORE<-list()
T2D_SCORE<-list()

  for( name in c("MURARO_PANCREAS_ALPHA_CELL.v7.5.11",
                 'MURARO_PANCREAS_BETA_CELL.v7.5.11')){
    for (nn in c(1:10)){ 

  for (group in c('ND','T2D')){

    if(group=='ND'){
      
      ND.name <- c(rownames(ss_obj_ND@meta.data)[grep(nn,ss_obj_ND$new_as_cluster )])
      ND.umap <- subset(ss_obj_ND,cells = ND.name)

      ND_RNA1 <-ND.umap@meta.data[[name]]
    }
    
    if(group=='T2D'){
      T2D.name <- c(rownames(ss_obj_T2D@meta.data)[grep(nn,ss_obj_T2D$new_as_cluster )])
      T2D.umap <- subset(ss_obj_T2D,cells = T2D.name)

      T2D_RNA1 <-T2D.umap@meta.data[[name]]
    }

      }     
      test<- ks.test(x=ND_RNA1 ,
                     y=T2D_RNA1 )
      test.matrix <- as.data.frame( do.call(rbind, test))
      colnames(test.matrix)<-paste0('group_',nn,name)
      
      if(i==1){p_value_dataframe<-test.matrix}
      if(i>1){
        p_value_dataframe<-cbind(p_value_dataframe,test.matrix)
       
  }
    i=i+1
}  
 
  }
write.csv(p_value_dataframe,'p_value_score.csv')
View(p_value_dataframe)
class(ss_obj_T2D@meta.data[[name]])
##############################
library(ggpubr)

library(ggsci)


 name<-as.data.frame( ss_obj@assays$RNA@counts)
kk_cluster <-(as.matrix( ss_obj@meta.data$vlnplot_group)) # kk_cluster是细胞分组信息

source('/violin_gene_exp_function_0.2_jitter.R')

pdf('average_group_vio.pdf',width = 12,height = 5)
for (i in 27:28){
  
  kk<-as.data.frame(t(ss_obj@meta.data[i]) )   # kk是表达矩阵
  colnames(kk)<-colnames(name)
  rownames(kk)<-colnames(ss_obj@meta.data)[i]
  violin_gene_exp(colnames(ss_obj@meta.data)[i],kk,kk_cluster)
}

dev.off()


