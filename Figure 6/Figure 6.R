library(SingleCellExperiment)
library(miloR)
library(Seurat)
library(statmod)
library(scater)
library(dplyr)
library(patchwork)


# panel A -----------------------------------------------------------------
DimPlot(object = endo_as_sce, reduction = "tsne", pt.size = 3, label = T,
        group.by = 'celltype',
        cols = c("#1F77B4FF", "#FF7F0EFF", "#D62728FF",
                 "#2CA02CFF", "#a020ef" ))
ggsave('spl_endo_sce/endo_as_ct_tsne.pdf', width = 8, height = 6.5)

# panel B -----------------------------------------------------------------
DimPlot(object = endo_as_sce, reduction = "tsne", pt.size = 3, label = T,split.by = 'Status',
        group.by = 'celltype',
        cols = c("#1F77B4FF", "#FF7F0EFF", "#D62728FF",
                 "#2CA02CFF", "#a020ef" ))
ggsave('spl_endo_sce/endo_as_ct_tsne_split.pdf', width = 10.5, height = 5.5)


# panel C -----------------------------------------------------------------
# endo_as_sce  endogroup
table(endo_as_sce$new_as_cluster)
endo_as_sce$milo_Sample=paste0(endo_as_sce$Status,'_',endo_as_sce$new_as_cluster)

endo_as_sce1=DietSeurat(endo_as_sce,dimreducs = 'tsne')
endo_as_sce1 <- as.SingleCellExperiment(endo_as_sce1)
traj_milo <- Milo(endo_as_sce1)
reducedDim(traj_milo, "TSNE") <- reducedDim(endo_as_sce1, "TSNE")

traj_milo <- buildGraph(traj_milo, k = 20, d = 20)
traj_milo <- makeNhoods(traj_milo, prop = 0.1, k = 20, d=20, refined = TRUE,
                        reduced_dims = 'TSNE')
plotNhoodSizeHist(traj_milo)

colnames(endo_as_sce@meta.data)
head(endo_as_sce@meta.data)
traj_milo <- countCells(traj_milo, meta.data = data.frame(colData(traj_milo)), 
                        samples="milo_Sample")
head(nhoodCounts(traj_milo))

traj_design <- data.frame(colData(traj_milo))[,c("milo_Sample", "Status")]
traj_design <- distinct(traj_design)
rownames(traj_design) <- traj_design$milo_Sample
## Reorder rownames to match columns of nhoodCounts(milo)
traj_design <- traj_design[colnames(nhoodCounts(traj_milo)), , drop=FALSE]

traj_design

traj_milo <- calcNhoodDistance(traj_milo, d=20)

rownames(traj_design) <- traj_design$milo_Sample
da_results <- testNhoods(traj_milo, design = ~ Status, design.df = traj_design)

write.csv(da_results,file = '/lawlor/miloR/lawlor_endo_as_miloR_k20d20.csv',
          quote = F)

da_results %>%
  arrange(- SpatialFDR) %>%
  head() 

traj_milo <- buildNhoodGraph(traj_milo)

plotNhoodGraphDA(traj_milo, da_results, 
                 alpha=1,layout ='TSNE') +
  plot_layout(guides="collect")
ggsave(file='/lawlor/miloR/lawlor_endo_as_miloR_FDR1_d20.pdf',
       width=3,height=1.6)



# panel D -----------------------------------------------------------------
law_beta_1to4=readRDS('/media/user/sdd/zhaolabguest/wsy/sc_as/lawlor.1050/diff/endo_diff_cluster/endoc1to4_betaonly/data/law_sc_gene_beta_cl1to4.rds')

gene_marker=FindAllMarkers(law_beta_1to4,group.by = 'new_as_cluster')
ENTREZID=bitr(gene_marker$gene,fromType = 'SYMBOL',
              toType = ('ENTREZID'),OrgDb = 'org.Hs.eg.db')
colnames(ENTREZID)[1]='gene'
gene_marker=left_join(gene_marker,ENTREZID)
write.table(gene_marker, 
            "/media/user/sdd/zhaolabguest/wsy/sc_as/lawlor.1050/supp_data/R3_Q12/beta_c1to4_gene_marker.txt",            
            row.names=FALSE,col.names=TRUE,sep="\t")

top5 <- gene_marker %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
law_beta_1to4$new_as_cluster=factor(law_beta_1to4$new_as_cluster,
                                    levels = c(4,3,2,1))
DotPlot(law_beta_1to4,features =top5$gene ,
        cols = c("white", 
                 '#ff181e'),
        group.by = 'new_as_cluster')+
  theme(axis.text.x = element_text (angle = 45,vjust = 0.5#,hjust = 0.7
  ))
ggsave(filename = "/media/user/sdd/zhaolabguest/wsy/sc_as/lawlor.1050/supp_data/R3_Q12/all_top5markers_dotplot.pdf",
       width =6,height = 3)


# panel E -----------------------------------------------------------------
ego <- enrichGO(gene          =gene_marker[gene_marker$cluster=='4',]$ENTREZID,
                OrgDb         = org.Hs.eg.db,
                pAdjustMethod = "BH",
                ont           = 'ALL' ,
                pvalueCutoff  = 0.1,
                #qvalueCutoff  = 0.99,
                readable      = TRUE)
dotplot(ego)


# panel F -----------------------------------------------------------------
VlnPlot(law_beta_1to4,features =c('INS','CHL1','SLC2A2','MAFA','RFX6','PDX1'),
        group.by = 'new_as_cluster',
        cols = c("#8C564BFF" ,"#E377C2FF", "#C49C94FF", "#BCBD22FF"))


# panel G -----------------------------------------------------------------
library(ggplot2)
library(grid)
library(gridExtra)
library(RColorBrewer)
# IsoheightMap

IsoheightMap <- function(tsne_result, gene, gene_expression){
  colnames(tsne_result)<- c("tSNE_1", "tSNE_2")
  warm <- colorRampPalette(c('#ffffb2','#fecc5c','#e31a1c','#800026'))
  mypalette <- c(warm(20))
  
  cold <- colorRampPalette(c('#f7fcf0','#41b6c4','lightgrey'))
  warm <- colorRampPalette(c('#ffffb2','#fecc5c','#e31a1c'))
  mypalette <- c(rev(cold(11)), warm(10))
  mypalette <-colorRampPalette(brewer.pal(11, "RdYlBu")[c(11:1)])(1000)
  rpkm <- gene_expression
  gene_exp <- rpkm[gene,,drop=FALSE]
  # gene_exp <- log(gene_exp+1)
  density_exp <- tsne_result[rownames(tsne_result) %in% names(gene_exp[,gene_exp>rowMeans(gene_exp)*1.2,drop=FALSE]),]
  print(density_exp)
  
  title <- grobTree(textGrob(gene, x=0.05, y=0.93, hjust = 0, gp=gpar(fontsize=10, fontface="bold.italic")))
  
  p <- ggplot(tsne_result, aes(tSNE_1, tSNE_2)) +
    geom_point(shape = 21, stroke=0.25, aes(colour=as.numeric(rpkm[gene,]) , fill=as.numeric(rpkm[gene,])), size = 2) +
    geom_density_2d(data=density_exp, aes(x=tSNE_1, y=tSNE_2), bins = 12, # 等高线的数量
                    colour="black",size=0.5) +
    # scale_fill_gradient2(high="darkred", low="yellow")+
    scale_fill_gradientn(colours = mypalette)+
    scale_colour_gradientn(colours = mypalette)+
    # scale_fill_viridis_c()+
    theme_bw() +
    # ggtitle(gene) +
    annotation_custom(title) +
    xlab("t-SNE 1") +
    ylab("t-SNE 2") +
    theme(
      plot.title = element_text(size=28, face="bold.italic", hjust = 0),
      # axis.text=element_text(size=12),
      # axis.title=element_text(size=16),
      axis.title=element_blank(),
      axis.text=element_blank(),
      legend.text = element_text(size =12),
      legend.title=element_blank(),
      aspect.ratio=1, # landscape
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position="none",
      axis.ticks = element_blank()
    )
  # print(p)
}
aa_3<-list(c('INS','SLC2A2','MAFA','RFX6','PDX1','CHL1',
             'GCK','PPARGC1A','MDH1','NEUROD1',
             'CREB1','G6PC2','PFKFB2','PFKM','SIX2','SIX3',
             'ENTPD3','GPD2','DNMT3A','MTOR'))

law_beta_1to4 <- AddModuleScore(  #对每个模块的基因集打分
  object = law_beta_1to4,
  features = aa_3,
  ctrl = 100, #默认值是100
  name = 'aa_3')


tsne_result.ch<-data.frame(law_beta_1to4@reductions[["as"]]@cell.embeddings); tsne_result.ch[1:8,]
expr.ch<-data.frame(law_beta_1to4@meta.data[["aa_31"]]);colnames(expr.ch)<-"mature_score"

rownames(expr.ch)<-rownames(tsne_result.ch)
expr.ch<-as.data.frame(t(expr.ch))
#embedding <- reducedDim(tsne.ch, "tSNE")

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



# panel H -----------------------------------------------------------------
library(Seurat)
library(ggplot2)
library(RColorBrewer)
ss_obj<-readRDS('/media/user/sdd/zhaolabguest/wsy/sc_as/lawlor.1050/diff/endo_diff_cluster/data/law_c1to4.rds')
name<-as.data.frame( ss_obj@assays$RNA@counts)
write.table(name,'/media/user/sdh/cmq_wj/others/c1_4_ND_T2D.counts2cytotarce.txt')
# web
cyto<-read.csv('/media/user/sdh/cmq_wj/others/CytoTRACE_results_c1_4.csv')
cyto_1<-as.data.frame(cyto[,2]) 
rownames(cyto_1) <-cyto[,1]
colnames(cyto_1)<-'CytoTRACE'
ss_obj<-AddMetaData(ss_obj, cyto_1, col.name = NULL)
FeaturePlot(ss_obj,features = 'MAFA',label = T,reduction ='new_as_cluster')

AS<- readRDS('/media/user/sdh/cmq_wj/others/as_ND_T2D.RDS') 
cell = AS@meta.data[AS@meta.data$new_as_cluster %in% c(1:4),]
law_beta_1to4<-subset(AS,cells= rownames(cell) )
law_beta_1to4<-AddMetaData(law_beta_1to4, cyto_1, col.name = NULL)
Idents(law_beta_1to4)<-'new_as_cluster'
FeaturePlot(law_beta_1to4,features = 'CytoTRACE',label = T,reduction ='as',
            cols=colorRampPalette(brewer.pal(11, "RdYlBu")[c(11:1)])(100))
ggsave('/media/user/sdh/cmq_wj/others/c1_4_cytotrace_as.pdf',width = 3.6,height = 3)

saveRDS(law_beta_1to4,'/media/user/sdh/cmq_wj/others/c1_4_cytotrace_as.rds')


# panel I -----------------------------------------------------------------
# lawlor beta cell monocle 
law_beta_sce

cds=as.CellDataSet(law_beta_sce)
cds=estimateSizeFactors(cds)
cds=estimateDispersions(cds)
marker.var=VariableFeatures(law_beta_sce)
cds=setOrderingFilter(cds,marker.var)
plot_ordering_genes(cds)

cds=reduceDimension(cds,reduction_method = 'DDRTree')
cds=orderCells(cds)
p1=plot_cell_trajectory(cds,color_by = 'new_as_cluster'
                        ,cell_size =1)+
  scale_color_manual( values=c("#8C564BFF" ,"#E377C2FF", "#C49C94FF", "#BCBD22FF", "#17BECFFF",
                               "#1F77B4FF", "#FF7F0EFF", "#2CA02CFF" ,"#D62728FF", "#9467BDFF")) + 
  theme(legend.position = "right")
p2=plot_complex_cell_trajectory(cds,x=1,y=2,
                                color_by = 'new_as_cluster',
                                cell_size =1)+
  scale_color_manual( values=c("#8C564BFF" ,"#E377C2FF", "#C49C94FF", "#BCBD22FF", "#17BECFFF",
                               "#1F77B4FF", "#FF7F0EFF", "#2CA02CFF" ,"#D62728FF", "#9467BDFF"))+
  facet_wrap('~new_as_cluster',nrow = 1)+ 
  theme(legend.title = element_blank())
p1|p2


# panel K -----------------------------------------------------------------
# beta cl4 AS marker--prepare for QUANTAS
table(law_sc_gene_beta_cl1to4$new_as_cluster)
beta_ascl=data.frame(Run=law_sc_gene_beta_cl1to4$Run,
                     new_as_cluster=law_sc_gene_beta_cl1to4$new_as_cluster)
head(beta_ascl)
for (i in 1:4){
  sub_group=beta_ascl
  sub_group$group=ifelse(sub_group$new_as_cluster==i,'group1','group2')
  sub_group=sub_group[,c('Run','group')]
  sub_group=sub_group[order(sub_group$group),]
  write.table(sub_group,file =paste0('/media/user/sdd/zhaolabguest/wsy/sc_as/lawlor.1050/diff/diff_conf/endoc1to4_betaonly/onlybeta_c',
                                     i,'.dataset.group.conf') ,
              sep='\t',row.names = F,col.names=F,quote = F)}

## volcano
law_endo_filter
law_beta_cl4_DE=read.table('/media/user/sdd/zhaolabguest/wsy/sc_as/lawlor.1050/diff/endo_diff_cluster/endoc1to4_betaonly/data/endo_onlybeta_c4_DE_dI0.1_ratio0.1.csv',
                           sep = ',',header = T)
id=str_split_fixed(law_beta_cl4_DE$gene,'//',2)
colnames(id)=c('ensemble','symbol')
law_beta_cl4_DE=cbind(id,law_beta_cl4_DE)
law_fmcov_obc4=law_beta_cl4_DE[law_beta_cl4_DE$coverage>=20,c(1:17)]
colnames(law_endo_filter)
fmdata=left_join(law_fmcov_obc4,law_endo_filter[,c(5,56,57)],by='name')
fmdata=subset(fmdata,fmdata[,18]>=0.1|
                fmdata[,19]>=0.1)
law_fmcov_obc4=subset(fmdata,fmdata[,18]>=0.1|
                        fmdata[,19]>=0.1)
law_fmcov_obc4$change = ifelse(law_fmcov_obc4$FDR < 0.05 & abs(law_fmcov_obc4$dI_g1_vs_g2) >= 0.1, 
                               ifelse(law_fmcov_obc4$dI_g1_vs_g2> 0.1 ,'Inclusion','exclusion'),
                               'Stable')
p <- ggplot(
  # 数据、映射、颜色
  law_fmcov_obc4, aes(x = dI_g1_vs_g2, y = -log10(FDR), colour=change)) +
  geom_point( alpha=0.8,size=2) +
  scale_color_manual(values=c(exclusion="#546de5", Stable="#d2dae2",Inclusion="#ff4757"))+
  # 辅助线
  geom_vline(xintercept=c(-0.1,0.1),lty=1,col="black",lwd=0.8) +
  geom_hline(yintercept = -log10(0.05),lty=1,col="black",lwd=0.8) +
  # 坐标轴
  labs(x="delta PSI",
       y="-log10 FDR",
       title='law.beta.ascluster4.DE')+
  xlim(-1,1)+
  theme_bw()+theme(panel.grid=element_blank())+
  # 图例
  theme(plot.title = element_text(hjust = 0.5), 
        legend.position="right", 
        legend.title = element_blank(),
        axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15))+
  geom_text_repel(data =law_fmcov_obc4[law_fmcov_obc4$symbol%in%c('SEC31A')|
                                         (law_fmcov_obc4$FDR<=10^(-75)&
                                            law_fmcov_obc4$coverage>=20&
                                            abs(law_fmcov_obc4$dI_g1_vs_g2)>=0.1),] ,
                  # law_fmcov_obc4[(law_fmcov_obc4$coverage>=20&
                  #                            abs(law_fmcov_obc4$dI_g1_vs_g2)>=0.4&
                  #                            law_fmcov_obc4$FDR<=0.05)|
                  #                           (law_fmcov_obc4$FDR<=10^(-75)&
                  #                              law_fmcov_obc4$coverage>=20&
                  #                              abs(law_fmcov_obc4$dI_g1_vs_g2)>=0.1),],
                  size = 3,color='black',
                  #direction = 'y',
                  box.padding = unit(0.5, "lines"),aes(label = symbol) ,
                  point.padding = unit(0.2, "lines"), segment.color = "black", show.legend = FALSE )
ggsave(p,filename = '/media/user/sdd/zhaolabguest/wsy/sc_as/lawlor.1050/diff/endo_diff_cluster/endoc1to4_betaonly/vol_law_as_c4_beta.pdf',
       width=8,height=7)


# panel L -----------------------------------------------------------------
## beta c4 up AS GO 
dim(law_beta_cl4_DE)
library(clusterProfiler)
law_beta_cl4_as_id=bitr(law_beta_cl4_DE$symbol,fromType = 'SYMBOL',
                        toType = c('ENSEMBL','ENTREZID'),
                        OrgDb = 'org.Hs.eg.db')
law_beta_cl4_ego <- enrichGO(gene      =law_beta_cl4_as_id$ENTREZID ,
                             OrgDb         = org.Hs.eg.db,
                             pAdjustMethod = "BH",
                             ont           = 'ALL' ,
                             readable      = TRUE)
View(law_beta_cl4_ego@result)
write.table(law_beta_cl4_ego@result,file = '/media/user/sdd/zhaolabguest/wsy/sc_as/lawlor.1050/diff/endo_diff_cluster/endoc1to4_betaonly/data/endo_onlybeta_c4_up_DE_GO_result.txt',
            quote = F,row.names = F,sep = '\t')

new_matrix=law_beta_cl4_ego@result[c(1,12,13,15,16,18,22,23,25,45),]
new_beta_cl4_ego=law_beta_cl4_ego
new_beta_cl4_ego@result=new_matrix
dotplot(new_beta_cl4_ego,showCategory=10)

dotplot(law_beta_cl4_ego,showCategory = c(11:25))+
  scale_y_discrete(labels=function(x) str_wrap(x, width=60))+
  guides(colour=guide_colorbar(reverse = T))
ggsave('/media/user/sdd/zhaolabguest/wsy/sc_as/lawlor.1050/diff/endo_diff_cluster/endoc1to4_betaonly/betaonly_c4_DE_up_dotplot.pdf',width = 5.8,height = 3)


# panel M -----------------------------------------------------------------

for (i in 1:4){
  law_onlybeta_c1_tvn_as=read.table(paste0('/media/user/sdd/zhaolabguest/wsy/sc_as/lawlor.1050/diff/diff_result/endoc1to4_betaonly/onlybeta_c',
                                           i,'.diff.txt'),
                                    sep = '\t',header = T)
  id=str_split_fixed(law_onlybeta_c1_tvn_as$gene,'//',2)
  colnames(id)=c('ensemble','symbol')
  law_onlybeta_c1_tvn_as=cbind(id,law_onlybeta_c1_tvn_as)
  law_fm_obc1=law_onlybeta_c1_tvn_as[law_onlybeta_c1_tvn_as$coverage>=20&
                                       abs(law_onlybeta_c1_tvn_as$dI_g1_vs_g2)>=0.1&
                                       law_onlybeta_c1_tvn_as$FDR<=0.05,c(1:17)]
  fmdata=left_join(law_fm_obc1,law_endo_filter[,c(5,(50+2*(i-1)),(51+2*(i-1)))],by='name')
  colnames(fmdata)[18:19]=paste0('law.',colnames(fmdata)[18:19])
  fmdata=subset(fmdata,fmdata[,18]>=0.1|
                  fmdata[,19]>=0.1)
  write.csv(fmdata,file = paste0('/media/user/sdd/zhaolabguest/wsy/sc_as/lawlor.1050/diff/endo_diff_cluster/endoc1to4_betaonly/data/endo_onlybeta_c',
                                 i,'_DE_dI0.1_ratio0.1.csv'),quote = F,row.names = F)
  fmdata_up=fmdata[fmdata$dI_g1_vs_g2>=0.1,]
  write.csv(fmdata_up,file = paste0('/media/user/sdd/zhaolabguest/wsy/sc_as/lawlor.1050/diff/endo_diff_cluster/endoc1to4_betaonly/data/endo_onlybeta_c',
                                    i,'_upDE_dI0.1_ratio0.1.csv'),quote = F,row.names = F)
  fmdata_down=fmdata[fmdata$dI_g1_vs_g2<=-0.1,]
  write.csv(fmdata_down,file = paste0('/media/user/sdd/zhaolabguest/wsy/sc_as/lawlor.1050/diff/endo_diff_cluster/endoc1to4_betaonly/data/endo_onlybeta_c',
                                      i,'_downDE_dI0.1_ratio0.1.csv'),quote = F,row.names = F)
}

## prepare for rMAPS
site_func=function(multiple_DE,all_exon_id){
  can_exon_name=str_split_fixed(multiple_DE$name,'-',6)%>%as.data.frame()
  colnames(can_exon_name)=paste0('c',1:6)
  str_name=str_split_fixed(can_exon_name$c6,'[INC]',2)%>%as.data.frame()
  can_exon_name$c7=can_exon_name$c6
  can_exon_name$c6=str_name$V1
  
  can_exon_name$exon_l=paste0('EX-',can_exon_name$c2,'-',can_exon_name$c3,'[')
  can_exon_name$exon_m=paste0('EX-',can_exon_name$c2,'-',can_exon_name$c4,'-',can_exon_name$c5,'[')
  can_exon_name$exon_r=paste0('EX-',can_exon_name$c2,'-',can_exon_name$c6)
  
  exon_l_site=left_join(can_exon_name,all_exon_id[,c(1:3,11)],by='exon_l')
  exon_m_site=left_join(can_exon_name,all_exon_id[,c(1:3,12)],by='exon_m')
  exon_r_site=left_join(can_exon_name,all_exon_id[,c(6,1:3,13)],by='exon_r')
  exon_site=cbind(exon_r_site,exon_m_site[,12:13],exon_l_site[,12:13])%>%na.omit()
  colnames(exon_site)[c(2,13:18)]=c('ensemble','firstFlankingExonStart','firstFlankingExonEnd','exonStart','exonEnd','secondFlankingExonStart','secondFlankingExonEnd')
  multiple_DE_exon_site=exon_site[,c(12,11,15,16,13,14,17,18)]
  return(multiple_DE_exon_site)
}
exon_matrix=read.table('/media/user/sda/MyDoc/Quantas/index_annot/hg19/annotation/hg19.exon.trio.hmr.nr.bed',
                       sep = '\t')
all_exon_id=str_split_fixed(exon_matrix$V4,'[.]',2)%>%as.data.frame()
all_exon_id=str_split_fixed(all_exon_id$V1,'-',4)%>%as.data.frame()
all_exon_id=cbind(exon_matrix[,1:6],all_exon_id)
colnames(all_exon_id)=c('chr','start','end','name','p','strand','type','gene','exonstart','exonend')
all_exon_id=all_exon_id[all_exon_id$type=='EX',]
all_exon_id$exon_l=paste0('EX-',all_exon_id$gene,'-',all_exon_id$exonend)
all_exon_id$exon_m=paste0('EX-',all_exon_id$gene,'-',all_exon_id$exonstart,'-',all_exon_id$exonend)
all_exon_id$exon_r=paste0('EX-',all_exon_id$gene,'-',all_exon_id$exonstart,'[')
head(all_exon_id)
for (i in 1:4){
  endoc1_obup_DE=read.csv(paste0('/media/user/sdd/zhaolabguest/wsy/sc_as/lawlor.1050/diff/endo_diff_cluster/endoc1to4_betaonly/data/endo_onlybeta_c',
                                 i,'_upDE_dI0.1_ratio0.1.csv'),
                          header = T)
  endoc1_obup_DE=site_func(endoc1_obup_DE,all_exon_id)
  write.table(endoc1_obup_DE,file = paste0('/media/user/sdd/zhaolabguest/wsy/sc_as/lawlor.1050/diff/endo_diff_cluster/endoc1to4_betaonly/data/for_rMAPS/endo_onlybeta_c',
                                           i,'_up_DE_exon_site_for_rMAPS.txt'),
              sep = '\t',col.names = F,row.names = F,quote = F)
  endoc1_obdown_DE=read.csv(paste0('/media/user/sdd/zhaolabguest/wsy/sc_as/lawlor.1050/diff/endo_diff_cluster/endoc1to4_betaonly/data/endo_onlybeta_c',
                                   i,'_downDE_dI0.1_ratio0.1.csv'),
                            header = T)
  endoc1_obdown_DE=site_func(endoc1_obdown_DE,all_exon_id)
  write.table(endoc1_obdown_DE,file = paste0('/media/user/sdd/zhaolabguest/wsy/sc_as/lawlor.1050/diff/endo_diff_cluster/endoc1to4_betaonly/data/for_rMAPS/endo_onlybeta_c',
                                             i,'_down_DE_exon_site_for_rMAPS.txt'),
              sep = '\t',col.names = F,row.names = F,quote = F)
}

# beta as cluster background exon
law_onlybeta_c4_tvn_as=read.table(paste0('/media/user/sdd/zhaolabguest/wsy/sc_as/lawlor.1050/diff/diff_result/endoc1to4_betaonly/onlybeta_c4.diff.txt'),
                                  sep = '\t',header = T)
id=str_split_fixed(law_onlybeta_c4_tvn_as$gene,'//',2)
colnames(id)=c('ensemble','symbol')
law_onlybeta_c4_tvn_as=cbind(id,law_onlybeta_c4_tvn_as)
law_backg_beta_as_cl=law_onlybeta_c4_tvn_as[law_onlybeta_c4_tvn_as$coverage>=20&
                                              abs(law_onlybeta_c4_tvn_as$dI_g1_vs_g2)<0.1
                                            ,c(1:17)]
write.csv(law_backg_beta_as_cl,file = '/media/user/sdd/zhaolabguest/wsy/sc_as/lawlor.1050/diff/endo_diff_cluster/endoc1to4_betaonly/data/endo_onlybeta_background.csv',
          quote = F,row.names = F)

backg_exon=read.csv('/media/user/sdd/zhaolabguest/wsy/sc_as/lawlor.1050/diff/endo_diff_cluster/endoc1to4_betaonly/data/endo_onlybeta_background.csv',
                    header = T)
backg_exon_site=site_func(law_backg_beta_as_cl,all_exon_id)
write.table(backg_exon_site,file = '/media/user/sdd/zhaolabguest/wsy/sc_as/lawlor.1050/diff/endo_diff_cluster/endoc1to4_betaonly/data/for_rMAPS/endo_onlybeta_background_exon_site_for_rMAPS.txt',
            sep = '\t',col.names = F,row.names = F,quote = F)


# panel N -----------------------------------------------------------------
VlnPlot(law_beta_1to4,features =c('FXR1','FMR1','HNRNPH2','HNRNPK'),
        group.by = 'new_as_cluster',
        cols = c("#8C564BFF" ,"#E377C2FF", "#C49C94FF", "#BCBD22FF"))
