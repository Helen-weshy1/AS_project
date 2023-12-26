library(Seurat)
library(ggplot2)
library(data.table)
library(dplyr)


# Prepare endocrine celltype differential splicing exons ------------------

# * calculate endocrine celltype DE ---------------------------------------
# law_Alpha_ct_as
# law_Beta_ct_as
# law_Multiple_ct_as
# law_Delta_ct_as
# law_PP_ct_as

# * DE filter matrix ------------------------------------------------------
# law_endo_filter
#law_sc_gene_endo
#law_endo_meta
#law_endo_filter
law_psi_junc=read.csv('/media/user/sdd/zhaolabguest/wsy/sc_as/lawlor.1050/raw_dataset/psi_junc_sc_data.csv')
dim(law_psi_junc)
a=read.table('/media/user/sdf/temp_wj/sc-AS/diab_pancr/Lawlor.Genom.Res.2016/count_AS_m4/SRR3617192/cass.count.txt', 
             header = F,sep = '\t')
rownames(law_psi_junc)=a$V4
law_psi_junc[1:4,1:4]
dim(law_psi_junc)
law_sc_gene_endo=readRDS('/media/user/sdd/zhaolabguest/wsy/sc_as/lawlor.1050/gene_sce/gene_endo_sce/law_sc_gene_endo.rds')
head(colnames(law_sc_gene_endo))
law_psi_junc_endo=law_psi_junc[,colnames(law_psi_junc)%in%colnames(law_sc_gene_endo)]
write.csv(law_psi_junc_endo,file = 'law_psi_junc_endo.csv',quote = F)
hist(apply(law_psi_junc_endo,1,function(x) ncol(law_psi_junc_endo)-sum(is.na(x)) ),
     breaks=100)

law_endo_meta=law_sc_gene_endo@meta.data
law_endo_filter=read.table('/media/user/sdd/zhaolabguest/wsy/sc_as/lawlor.1050/diff/diff_result/BetaTvN.diff.txt',
                           sep = '\t',header = T)[,1:7]
table(law_endo_meta$celltype)

# beta
beta=subset(law_endo_meta,law_endo_meta$celltype=='Beta')
beta_psi_junc=law_psi_junc_endo[,colnames(law_psi_junc_endo)%in%rownames(beta)]
beta_ratio=apply(beta_psi_junc,1,function(x) (ncol(beta_psi_junc)-sum(is.na(x)))/ncol(beta_psi_junc) )
law_endo_filter$beta_ratio=beta_ratio

nbeta=subset(law_endo_meta,law_endo_meta$celltype!='Beta')
nbeta_psi_junc=law_psi_junc_endo[,colnames(law_psi_junc_endo)%in%rownames(nbeta)]
nbeta_ratio=apply(nbeta_psi_junc,1,function(x) (ncol(nbeta_psi_junc)-sum(is.na(x)))/ncol(nbeta_psi_junc) )
law_endo_filter$nbeta_ratio=nbeta_ratio

# alpha
alpha=subset(law_endo_meta,law_endo_meta$celltype=='Alpha')
alpha_psi_junc=law_psi_junc_endo[,colnames(law_psi_junc_endo)%in%rownames(alpha)]
alpha_ratio=apply(alpha_psi_junc,1,function(x) (ncol(alpha_psi_junc)-sum(is.na(x)))/ncol(alpha_psi_junc) )
law_endo_filter$alpha_ratio=alpha_ratio

nalpha=subset(law_endo_meta,law_endo_meta$celltype!='Alpha')
nalpha_psi_junc=law_psi_junc_endo[,colnames(law_psi_junc_endo)%in%rownames(nalpha)]
nalpha_ratio=apply(nalpha_psi_junc,1,function(x) (ncol(nalpha_psi_junc)-sum(is.na(x)))/ncol(nalpha_psi_junc) )
law_endo_filter$nalpha_ratio=nalpha_ratio

# multiple
multiple=subset(law_endo_meta,law_endo_meta$celltype=='Multiple')
multiple_psi_junc=law_psi_junc_endo[,colnames(law_psi_junc_endo)%in%rownames(multiple)]
multiple_ratio=apply(multiple_psi_junc,1,function(x) (ncol(multiple_psi_junc)-sum(is.na(x)))/ncol(multiple_psi_junc) )
law_endo_filter$multiple_ratio=multiple_ratio

nmultiple=subset(law_endo_meta,law_endo_meta$celltype!='Multiple')
nmultiple_psi_junc=law_psi_junc_endo[,colnames(law_psi_junc_endo)%in%rownames(nmultiple)]
nmultiple_ratio=apply(nmultiple_psi_junc,1,function(x) (ncol(nmultiple_psi_junc)-sum(is.na(x)))/ncol(nmultiple_psi_junc) )
law_endo_filter$nmultiple_ratio=nmultiple_ratio

# delta
delta=subset(law_endo_meta,law_endo_meta$celltype=='Delta')
delta_psi_junc=law_psi_junc_endo[,colnames(law_psi_junc_endo)%in%rownames(delta)]
delta_ratio=apply(delta_psi_junc,1,function(x) (ncol(delta_psi_junc)-sum(is.na(x)))/ncol(delta_psi_junc) )
law_endo_filter$delta_ratio=delta_ratio

ndelta=subset(law_endo_meta,law_endo_meta$celltype!='Delta')
ndelta_psi_junc=law_psi_junc_endo[,colnames(law_psi_junc_endo)%in%rownames(ndelta)]
ndelta_ratio=apply(ndelta_psi_junc,1,function(x) (ncol(ndelta_psi_junc)-sum(is.na(x)))/ncol(ndelta_psi_junc) )
law_endo_filter$ndelta_ratio=ndelta_ratio

# pp
pp=subset(law_endo_meta,law_endo_meta$celltype=='PP')
pp_psi_junc=law_psi_junc_endo[,colnames(law_psi_junc_endo)%in%rownames(pp)]
pp_ratio=apply(pp_psi_junc,1,function(x) (ncol(pp_psi_junc)-sum(is.na(x)))/ncol(pp_psi_junc) )
law_endo_filter$pp_ratio=pp_ratio

npp=subset(law_endo_meta,law_endo_meta$celltype!='PP')
npp_psi_junc=law_psi_junc_endo[,colnames(law_psi_junc_endo)%in%rownames(npp)]
npp_ratio=apply(npp_psi_junc,1,function(x) (ncol(npp_psi_junc)-sum(is.na(x)))/ncol(npp_psi_junc) )
law_endo_filter$npp_ratio=npp_ratio

fm_law_Alpha_ct_as=law_Alpha_ct_as[law_Alpha_ct_as$coverage>=20&
                                     law_Alpha_ct_as$FDR<=0.05&
                                     abs(law_Alpha_ct_as$dI_g1_vs_g2)>=0.1&
                                     law_Alpha_ct_as$name%in%law_endo_filter[law_endo_filter$alpha_ratio>=0.1|
                                                                               law_endo_filter$nalpha_ratio>=0.1,]$name,]
fm_law_Beta_ct_as=law_Beta_ct_as[law_Beta_ct_as$coverage>=20&
                                   law_Beta_ct_as$FDR<=0.05&
                                   abs(law_Beta_ct_as$dI_g1_vs_g2)>=0.1&
                                   law_Beta_ct_as$name%in%law_endo_filter[law_endo_filter$beta_ratio>=0.1|
                                                                            law_endo_filter$nbeta_ratio>=0.1,]$name,]
fm_law_Multiple_ct_as=law_Multiple_ct_as[law_Multiple_ct_as$coverage>=20&
                                           law_Multiple_ct_as$FDR<=0.05&
                                           abs(law_Multiple_ct_as$dI_g1_vs_g2)>=0.1&
                                           law_Multiple_ct_as$name%in%law_endo_filter[law_endo_filter$multiple_ratio>=0.1|
                                                                                        law_endo_filter$nmultiple_ratio>=0.1,]$name,]
fm_law_Delta_ct_as=law_Delta_ct_as[law_Delta_ct_as$coverage>=20&
                                     law_Delta_ct_as$FDR<=0.05&
                                     abs(law_Delta_ct_as$dI_g1_vs_g2)>=0.1&
                                     law_Delta_ct_as$name%in%law_endo_filter[law_endo_filter$delta_ratio>=0.1|
                                                                               law_endo_filter$ndelta_ratio>=0.1,]$name,]
fm_law_PP_ct_as=law_PP_ct_as[law_PP_ct_as$coverage>=20&
                               law_PP_ct_as$FDR<=0.05&
                               abs(law_PP_ct_as$dI_g1_vs_g2)>=0.1&
                               law_PP_ct_as$name%in%law_endo_filter[law_endo_filter$pp_ratio>=0.1|
                                                                      law_endo_filter$npp_ratio>=0.1,]$name,]


# panel A -----------------------------------------------------------------
# combine diff exon
law_endo_ct_as=cbind(law_Beta_ct_as[,c(1:5,13)],
                     law_Multiple_ct_as$dI_g1_vs_g2,
                     law_Alpha_ct_as$dI_g1_vs_g2,
                     law_Delta_ct_as$dI_g1_vs_g2,
                     law_PP_ct_as$dI_g1_vs_g2)
colnames(law_endo_ct_as)[6:10]=paste0(c('β','β/α','α','δ','PP'),'.','dI')
fm_law_Alpha_ct_as$celltype='α'
fm_law_Beta_ct_as$celltype='β'
fm_law_Multiple_ct_as$celltype='β/α'
fm_law_Delta_ct_as$celltype='δ'
fm_law_PP_ct_as$celltype='PP'
# differential inclusive splicing exons
fm_law_endo_ct_as_in=rbind(fm_law_Beta_ct_as[fm_law_Beta_ct_as$dI_g1_vs_g2>0,c(1:5,28)],
                           fm_law_Multiple_ct_as[fm_law_Multiple_ct_as$dI_g1_vs_g2>0,c(1:5,28)],
                           fm_law_Alpha_ct_as[fm_law_Alpha_ct_as$dI_g1_vs_g2>0,c(1:5,28)],
                           fm_law_Delta_ct_as[fm_law_Delta_ct_as$dI_g1_vs_g2>0,c(1:5,28)],
                           fm_law_PP_ct_as[fm_law_PP_ct_as$dI_g1_vs_g2>0,c(1:5,28)])
fm_law_endo_ct_as_in=left_join(fm_law_endo_ct_as_in,law_endo_ct_as[,5:10],by='name')

library(pheatmap)
count=na.omit(fm_law_endo_ct_as_in[,6:11])
count_s=count[,-1]
n=scale(count_s,  scale = TRUE,center = T)
n[n>2]=2
n[n< -2]= -2
pheatmap(n,show_colnames =F,show_rownames = F,cluster_rows = F,cluster_cols=F,
         color = colorRampPalette(c( "white",'#efdadc',"#b71819"))(100)) 
# differential exclusive splicing exons
fm_law_endo_ct_as_ex=rbind(fm_law_Beta_ct_as[fm_law_Beta_ct_as$dI_g1_vs_g2<0,c(1:5,28)],
                           fm_law_Multiple_ct_as[fm_law_Multiple_ct_as$dI_g1_vs_g2<0,c(1:5,28)],
                           fm_law_Alpha_ct_as[fm_law_Alpha_ct_as$dI_g1_vs_g2<0,c(1:5,28)],
                           fm_law_Delta_ct_as[fm_law_Delta_ct_as$dI_g1_vs_g2<0,c(1:5,28)],
                           fm_law_PP_ct_as[fm_law_PP_ct_as$dI_g1_vs_g2<0,c(1:5,28)])
fm_law_endo_ct_as_ex=left_join(fm_law_endo_ct_as_ex,law_endo_ct_as[,5:10],by='name')
count=na.omit(fm_law_endo_ct_as_ex[,6:11])
count_s=count[,-1]
n=scale(count_s,  scale = TRUE,center = T)
n[n>2]=2
n[n< -2]= -2
pheatmap(n,show_colnames =F,show_rownames = F,cluster_rows = F,cluster_cols=F,
         color = colorRampPalette(c( '#161af9','#c2c3ef',"white"))(100)) 


# panel B -----------------------------------------------------------------
# significant GO term intersect
# Alpha_ct_ego,Beta_ct_ego,Multiple_ct_ego,Delta_ct_ego,PP_ct_ego
library(pheatmap)
library(VennDiagram)
library(stringr)
venn.plot=venn.diagram(list(Alpha_ct_ego=Alpha_ct_ego@result$Description,
                            Beta_ct_ego=Beta_ct_ego@result$Description,
                            Multiple_ct_ego=Multiple_ct_ego@result$Description,
                            Delta_ct_ego=Delta_ct_ego@result$Description,
                            PP_ct_ego=PP_ct_ego@result$Description),
                       filename = NULL)
grid.draw(venn.plot)

inter <- get.venn.partitions(list(Alpha_ct_ego=Alpha_ct_ego@result$Description,
                                  Beta_ct_ego=Beta_ct_ego@result$Description,
                                  Multiple_ct_ego=Multiple_ct_ego@result$Description,
                                  Delta_ct_ego=Delta_ct_ego@result$Description,
                                  PP_ct_ego=PP_ct_ego@result$Description))
for (i in 1:nrow(inter)) inter[i,'values'] <- paste(inter[[i,'..values..']], collapse = ', ')
list1=inter$values[1]
list2=as.character(str_split_fixed(list1,', ',20))

alpha_m=Alpha_ct_ego@result[Alpha_ct_ego$Description%in%list2,]
beta_m=Beta_ct_ego@result[Beta_ct_ego$Description%in%list2,]
multi_m=Multiple_ct_ego@result[Multiple_ct_ego$Description%in%list2,]
delta_m=Delta_ct_ego@result[Delta_ct_ego$Description%in%list2,]
pp_m=PP_ct_ego@result[PP_ct_ego$Description%in%list2,]
alpha_m=alpha_m[,c(3,11)]
beta_m=beta_m[,c(3,11)]
multi_m=multi_m[,c(3,11)]
delta_m=delta_m[,c(3,11)]
pp_m=pp_m[,c(3,11)]
inter_heatmap=left_join(beta_m,
                        multi_m,
                        by='Description')
inter_heatmap=left_join(inter_heatmap,
                        alpha_m,
                        by='Description')
inter_heatmap=left_join(inter_heatmap,
                        delta_m,
                        by='Description')
inter_heatmap=left_join(inter_heatmap,
                        pp_m,
                        by='Description')
rownames(inter_heatmap)=inter_heatmap$Description
inter_heatmap=inter_heatmap[,-1]
colnames(inter_heatmap)[1:5]=paste0(c('β','β/α','α','δ','pp'),'_enrichment_fold')
n=scale(inter_heatmap,  scale = TRUE,center = T)
n[n>2]=2
n[n< -1]= -1
pheatmap(n,show_colnames =T,show_rownames = T,cluster_rows = F,cluster_cols=F,
         color = colorRampPalette(c( "white",'#f6cfa1','#eea146',"#ea8209"))(100),
         border_color = 'white')

# all GO term
#endo_ct_go_all.Rdata
alpha_m_all=Alpha_ct_ego_all@result
beta_m_all=Beta_ct_ego_all@result
multi_m_all=Multiple_ct_ego_all@result
delta_m_all=Delta_ct_ego_all@result
pp_m_all=PP_ct_ego_all@result
alpha_m_all=alpha_m_all[,c(3,11)]
beta_m_all=beta_m_all[,c(3,11)]
multi_m_all=multi_m_all[,c(3,11)]
delta_m_all=delta_m_all[,c(3,11)]
pp_m_all=pp_m_all[,c(3,11)]
inter_heatmap=left_join(beta_m_all,
                        multi_m_all,
                        by='Description')
inter_heatmap=left_join(inter_heatmap,
                        alpha_m_all,
                        by='Description')
inter_heatmap=left_join(inter_heatmap,
                        delta_m_all,
                        by='Description')
inter_heatmap=left_join(inter_heatmap,
                        pp_m_all,
                        by='Description')
go_list=rbind(Alpha_ct_ego@result,
              Beta_ct_ego@result,
              Multiple_ct_ego@result,
              Delta_ct_ego@result,
              PP_ct_ego@result)
go_list_name=go_list[,c('Description','ID')]%>%unique()
inter_heatmap_all=left_join(go_list_name,inter_heatmap,by='Description')
rownames(inter_heatmap_all)=inter_heatmap_all$Description
colnames(inter_heatmap_all)[3:7]=paste0(c('β','β/α','α','δ','pp'),'_enrichment_fold')
n=scale(inter_heatmap_all[,3:7],  scale = T,center = T)
n[n>2]=2
n[n< -1]= -1
pheatmap(n,show_colnames =T,show_rownames = T,cluster_rows = F,cluster_cols=F,
         color = colorRampPalette(c( "white",'#f6cfa1','#eea146',"#ea8209"))(100),
         border_color = 'white')


# panel C -----------------------------------------------------------------
# endo gene marker heatmap
endo_g_marker
top200 <- endo_g_marker %>% group_by(cluster) %>% top_n(n = 200, wt = avg_logFC)
fil_feature=rbind(top200[top200$cluster==c('Beta'),],
                  top200[top200$cluster==c('Multiple'),],
                  top200[top200$cluster==c('Alpha'),],
                  top200[top200$cluster==c('Delta'),],
                  top200[top200$cluster==c('PP'),])
DoHeatmap(law_sc_gene_endo, features = fil_feature$gene,group.by = 'celltype')
ggsave(filename = "endo_g_top200markers_Heatmap.pdf",width = 15,height = 15)


# panel D -----------------------------------------------------------------
# intersection of DE and marker genes
## venn
marker_gene
alpha_venn=venn.diagram(list(exon=unique(fm_law_Alpha_ct_as$gene),
                             gene_expr=marker_gene[marker_gene$cluster=='Alpha',]$gene),
                        fill=c("#1F77B4FF", "#FF7F0EFF"),
                        cat.cex=1,
                        cex=2,
                        height = 1500, width =1500,
                        cat.col = c("#1F77B4FF", "#FF7F0EFF"),
                        cat.fontfamily = "serif",
                        cat.fontface = "bold",
                        margin = 0.05,
                        filename = NULL)
grid.draw(alpha_venn)
beta_venn=venn.diagram(list(exon=unique(fm_law_Beta_ct_as$gene),
                            gene_expr=marker_gene[marker_gene$cluster=='Beta',]$gene),
                       fill=c("#1F77B4FF", "#FF7F0EFF"),
                       cat.cex=1,
                       cex=2,
                       height = 1500, width =1500,
                       cat.col = c("#1F77B4FF", "#FF7F0EFF"),
                       cat.fontfamily = "serif",
                       cat.fontface = "bold",
                       margin = 0.05,
                       filename = NULL)
grid.draw(beta_venn)
multi_venn=venn.diagram(list(exon=unique(fm_law_Multiple_ct_as$gene),
                             gene_expr=marker_gene[marker_gene$cluster=='Multiple',]$gene),
                        fill=c("#1F77B4FF", "#FF7F0EFF"),
                        cat.cex=1,
                        cex=2,
                        height = 1500, width =1500,
                        cat.col = c("#1F77B4FF", "#FF7F0EFF"),
                        cat.fontfamily = "serif",
                        cat.fontface = "bold",
                        margin = 0.05,
                        filename = NULL)
grid.draw(multi_venn)
delta_venn=venn.diagram(list(exon=unique(fm_law_Delta_ct_as$gene),
                             gene_expr=marker_gene[marker_gene$cluster=='Delta',]$gene),
                        fill=c("#1F77B4FF", "#FF7F0EFF"),
                        cat.cex=1,
                        cex=2,
                        height = 1500, width =1500,
                        cat.col = c("#1F77B4FF", "#FF7F0EFF"),
                        cat.fontfamily = "serif",
                        cat.fontface = "bold",
                        margin = 0.05,
                        filename = NULL)
grid.draw(delta_venn)
PP_venn=venn.diagram(list(exon=unique(fm_law_PP_ct_as$gene),
                          gene_expr=marker_gene[marker_gene$cluster=='PP',]$gene),
                     fill=c("#1F77B4FF", "#FF7F0EFF"),
                     cat.cex=1,
                     cex=2,
                     height = 1500, width =1500,
                     cat.col = c("#1F77B4FF", "#FF7F0EFF"),
                     cat.fontfamily = "serif",
                     cat.fontface = "bold",
                     margin = 0.05,
                     filename = NULL)
grid.draw(PP_venn)


# panel E -----------------------------------------------------------------
## ct dI and gene expression cor
law_endo_allg_fc=FindAllMarkers(law_sc_gene_endo,logfc.threshold = 0,return.thresh = 1)
Beta_DE_FC=left_join(fm_law_Beta_ct_as[,c(1,13)],
                     law_endo_allg_fc[law_endo_allg_fc$cluster=='Beta' 
                                      & law_endo_allg_fc$p_val<=0.05,
                                      c(2,7)],
                     by='gene')
ggplot(Beta_DE_FC,aes(x=avg_logFC,y=dI_g1_vs_g2))+
  geom_point(size=1,color='#162B87')+
  ggtitle( 'Beta')+
  theme_test()
Multiple_DE_FC=left_join(fm_law_Multiple_ct_as[,c(1,13)],
                         law_endo_allg_fc[law_endo_allg_fc$cluster=='Multiple' 
                                          & law_endo_allg_fc$p_val<=0.05,
                                          c(2,7)],
                         by='gene')
ggplot(Multiple_DE_FC,aes(x=avg_logFC,y=dI_g1_vs_g2))+
  geom_point(size=1,color='#162B87')+
  ggtitle( 'Multiple')+
  theme_test()
Alpha_DE_FC=left_join(fm_law_Alpha_ct_as[,c(1,13)],
                      law_endo_allg_fc[law_endo_allg_fc$cluster=='Alpha' 
                                       & law_endo_allg_fc$p_val<=0.05,
                                       c(2,7)],
                      by='gene')
ggplot(Alpha_DE_FC,aes(x=avg_logFC,y=dI_g1_vs_g2))+
  geom_point(size=1,color='#162B87')+
  ggtitle( 'Alpha')+
  theme_test()
Delta_DE_FC=left_join(fm_law_Delta_ct_as[,c(1,13)],
                      law_endo_allg_fc[law_endo_allg_fc$cluster=='Delta' 
                                       & law_endo_allg_fc$p_val<=0.05,
                                       c(2,7)],
                      by='gene')
ggplot(Delta_DE_FC,aes(x=avg_logFC,y=dI_g1_vs_g2))+
  geom_point(size=1,color='#162B87')+
  ggtitle( 'Delta')+
  theme_test()
PP_DE_FC=left_join(fm_law_PP_ct_as[,c(1,13)],
                   law_endo_allg_fc[law_endo_allg_fc$cluster=='PP' 
                                    & law_endo_allg_fc$p_val<=0.05,
                                    c(2,7)],
                   by='gene')
ggplot(PP_DE_FC,aes(x=avg_logFC,y=dI_g1_vs_g2))+
  geom_point(size=1,color='#162B87')+
  ggtitle( 'PP')+
  theme_test()


# panel F, H --------------------------------------------------------------
# ct DE example--using VALERIE
# chr9:139744465:139744589:+@chr9:139744958:139745012:+@chr9:139745207:139745474
PlotPSI(
  tran_id="chr4:84378062:84378111:+@chr4:84379499:84379582:+@chr4:84380893:84380950",
  event.type="SE",
  strand="positive",
  Bam="/media/user/sdg/WShi/AS/lawlor/marvel/bam/",
  BamPheno=BamPheno,
  cell.types=c( 'Alpha','Beta','Multiple','Delta','PP'),
  min.coverage=10,
  cons.exon.cutoff=100,
  method="ks",
  method.adj="bonferroni",
  cell.types.colors="ggplot.default",
  plot.title="MRPS18C",
  plot.width=5,
  plot.height=8,
  plot.out="/media/user/sdg/WShi/AS/lawlor/VALERIE/figure/MRPS18C_ct.pdf"
)

# panel G, I --------------------------------------------------------------
## gene expression
law_sc_gene_endo
VlnPlot(law_sc_gene_endo,features = 'SLC30A8',pt.size = 0.2,group.by = 'celltype')+
  scale_fill_manual(values=c( "#FF7F0EFF",  "#2CA02CFF", "#1F77B4FF","#D62728FF", "#a020ef"))
ggsave(file='alpha_SLC30A8_gene_vln.pdf',width=4,height=3)
VlnPlot(law_sc_gene_endo,features = 'MRPS18C',pt.size = 0.2,group.by = 'celltype')+
  scale_fill_manual(values=c("#FF7F0EFF",  "#2CA02CFF", "#1F77B4FF", "#D62728FF", "#a020ef"))
ggsave(file='alpha_MRPS18C_gene_vln.pdf',width=4,height=3)
