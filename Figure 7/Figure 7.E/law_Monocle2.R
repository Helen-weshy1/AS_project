devtools::load_all('/library/monocle/')
library(data.table)
library(Seurat)
library(dplyr)
library(stringr)
library(harmony)
library(ggplot2)
library(ggsci)

law_sce=readRDS('/law_sc_gene_endocrine.rds')
law_sce
################################### Monocle ####################################
table(Idents(law_sce))
law_beta_sce=law_sce
Idents(law_beta_sce)=law_beta_sce$Status
table(Idents(law_beta_sce))

# * ND --------------------------------------------------------------------
law_beta_sce=subset(law_beta_sce,idents='ND')
Idents(law_beta_sce)=law_beta_sce$new_as_cluster
law_beta_sce=subset(law_beta_sce,idents=c(1,2,3,4))
table(Idents(law_beta_sce))

cds=as.CellDataSet(law_beta_sce)
cds=estimateSizeFactors(cds)
cds=estimateDispersions(cds)

marker.var=VariableFeatures(law_beta_sce)
cds=setOrderingFilter(cds,marker.var)
plot_ordering_genes(cds)
cds=reduceDimension(cds,reduction_method = 'DDRTree')
cds=orderCells(cds)
plot_cell_trajectory(cds,color_by = 'new_as_cluster'
                     ,cell_size =1)+
  scale_color_manual( values=c("#8C564BFF" ,"#E377C2FF", "#C49C94FF", "#BCBD22FF", "#17BECFFF",
                               "#1F77B4FF", "#FF7F0EFF", "#2CA02CFF" ,"#D62728FF", "#9467BDFF")) + 
  theme(legend.position = "right")
plot_complex_cell_trajectory(cds,x=1,y=2,
                             color_by = 'new_as_cluster',
                             cell_size =1)+
  scale_color_manual( values=c("#8C564BFF" ,"#E377C2FF", "#C49C94FF", "#BCBD22FF", "#17BECFFF",
                               "#1F77B4FF", "#FF7F0EFF", "#2CA02CFF" ,"#D62728FF", "#9467BDFF"))+
  facet_wrap('~new_as_cluster',nrow = 1)+ 
  theme(legend.title = element_blank())

# * T2D --------------------------------------------------------------------
law_beta_sce=subset(law_beta_sce,idents='T2D')
Idents(law_beta_sce)=law_beta_sce$new_as_cluster
law_beta_sce=subset(law_beta_sce,idents=c(1,2,3,4))
table(Idents(law_beta_sce))

cds=as.CellDataSet(law_beta_sce)
cds=estimateSizeFactors(cds)
cds=estimateDispersions(cds)

marker.var=VariableFeatures(law_beta_sce)
cds=setOrderingFilter(cds,marker.var)
plot_ordering_genes(cds)
cds=reduceDimension(cds,reduction_method = 'DDRTree')
cds=orderCells(cds)
plot_cell_trajectory(cds,color_by = 'new_as_cluster'
                     ,cell_size =1)+
  scale_color_manual( values=c("#8C564BFF" ,"#E377C2FF", "#C49C94FF", "#BCBD22FF", "#17BECFFF",
                               "#1F77B4FF", "#FF7F0EFF", "#2CA02CFF" ,"#D62728FF", "#9467BDFF")) + 
  theme(legend.position = "right")
plot_complex_cell_trajectory(cds,x=1,y=2,
                             color_by = 'new_as_cluster',
                             cell_size =1)+
  scale_color_manual( values=c("#8C564BFF" ,"#E377C2FF", "#C49C94FF", "#BCBD22FF", "#17BECFFF",
                               "#1F77B4FF", "#FF7F0EFF", "#2CA02CFF" ,"#D62728FF", "#9467BDFF"))+
  facet_wrap('~new_as_cluster',nrow = 1)+ 
  theme(legend.title = element_blank())

