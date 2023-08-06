###################################### IsoheightMap #######################################
######################################## comb data ########################################
rm(list=ls())
setwd('/media/user/sdh/cmq_wj/others')
.libPaths("/media/user/sdf/R.lib/library/wj/sc")
library(Seurat)
library(ggplot2)
AS1_4<-readRDS('/media/user/sdh/cmq_wj/others/c1_4_cytotrace_as.rds')

aa_3<-list(c('INS','SLC2A2','MAFA','RFX6','PDX1','CHL1',
           'GCK','PPARGC1A','MDH1','NEUROD1',
           'CREB1','G6PC2','PFKFB2','PFKM','SIX2','SIX3',
           'ENTPD3','GPD2','DNMT3A','MTOR'))

AS1_4 <- AddModuleScore(  #对每个模块的基因集打分
  object = AS1_4,
  features = aa_3,
  ctrl = 100, #默认值是100
  name = 'aa_3')

#########
# comb
# isoheight map
# setwd("/media/user/sdb/lys.data/wnt10b_msc_sc_smart/SalmonAH/analysis/2.output/new.plot")
# tsne.ch<-readRDS("../A.Files.Seurat/Salmon_counts_nfea2000_Int_tsne.ch5.rds"); head(Idents(tsne.ch))

tsne_result.ch<-data.frame(AS1_4@reductions[["as"]]@cell.embeddings); tsne_result.ch[1:8,]
expr.ch<-data.frame(AS1_4@meta.data[["aa_31"]]);colnames(expr.ch)<-"mature_score"

rownames(expr.ch)<-rownames(tsne_result.ch)
expr.ch<-as.data.frame(t(expr.ch))
#embedding <- reducedDim(tsne.ch, "tSNE")

source('/media/user/sdh/CMQ_data/MAP.R') #等高线图代码

pdf("IsoheightMap_tsne.comb_df.aa3_1.5_12r.pdf", width=3.3, height=3.3)
p <- list()
for (genes in rownames(expr.ch) ) {
  p[[genes]] <- IsoheightMap(
    tsne_result.ch,
    genes, 
    expr.ch)
}
do.call(grid.arrange, c(p, ncol=1)) #(3,3) #(15,3*)
dev.off()

