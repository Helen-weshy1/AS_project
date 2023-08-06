library(data.table)
library(Seurat)
library(dplyr)
library(stringr)
library(ggplot2)
library(ggsci)

xin_sce=readRDS('E:/windows/AS/20230212Fig/data/Fig.S2/A_cluster_sc.rds')
xin_sce
################################### Monocle ####################################
.libPaths('F:/Rprogram/R-4.0.3/library/')
devtools::load_all('F:/Rprogram/R-4.0.3/library/monocle/')
table(Idents(xin_sce))
xin_beta_sce=subset(xin_sce,idents='Beta')
table(Idents(xin_beta_sce))

colnames(xin_beta_sce@meta.data)
cellinfo=subset(xin_beta_sce@meta.data,select = c("orig.ident",
                                                   'celltype',
                                                   'condition'))
xin_beta_sce=CreateSeuratObject(xin_beta_sce@assays$RNA@counts,
                                 meta.data = cellinfo)
xin_beta_sce=SCTransform(xin_beta_sce)
xin_beta_sce=RunPCA(xin_beta_sce)

ElbowPlot(xin_beta_sce,ndims = 50)

xin_beta_sce <- RunUMAP(object = xin_beta_sce, dims = 1:10)
p1=DimPlot(xin_beta_sce)+scale_color_npg()
p2=DimPlot(xin_beta_sce,split.by = 'condition')+scale_color_npg()
pdf('E:/windows/AS/article/dataset/xin/sce_split_status_dimplot.pdf',
    width=12,height = 4)
p1+p2+plot_layout(widths=c(1,2))
dev.off()
saveRDS(xin_beta_sce,file = 'E:/windows/AS/article/dataset/xin/xin_beta_sce.rds')

table(Idents(xin_beta_sce))
colnames(xin_beta_sce@meta.data)
cds=as.CellDataSet(xin_beta_sce)
cds=estimateSizeFactors(cds)
cds=estimateDispersions(cds)
###
# marker.cluster=FindAllMarkers(xin_beta_sce)
# cds=setOrderingFilter(cds,marker.cluster$gene)
# plot_ordering_genes(cds)
# ###
marker.var=VariableFeatures(xin_beta_sce)
cds=setOrderingFilter(cds,marker.var)
plot_ordering_genes(cds)
# ###
# marker.mono=dispersionTable(cds)
# colnames(marker.mono)
# marker.mono=subset(marker.mono,mean_expression>=0.1 & 
#                      dispersion_empirical>=1*dispersion_fit)$gene_id
# cds=setOrderingFilter(cds,marker.mono)
# plot_ordering_genes(cds)
# ###
# expressed_genes=VariableFeatures(xin_beta_sce)
# diff_test_res <- differentialGeneTest(cds[expressed_genes,],
#                                       fullModelFormulaStr = "~new_as_cluster")
# ordering_genes <- row.names (subset(diff_test_res, qval < 0.01)) ## 不要也写0.1 ，而是要写0.01。
# 
# cds <- setOrderingFilter(cds, ordering_genes)
# plot_ordering_genes(cds)



cds=reduceDimension(cds,reduction_method = 'DDRTree')
cds=orderCells(cds)
#cds=orderCells(cds,root_state = 4)

p1=DimPlot(xin_beta_sce)
p2=plot_cell_trajectory(cds,cell_size =1)
p3=plot_cell_trajectory(cds,color_by = 'Pseudotime',cell_size =1)
p4=plot_cell_trajectory(cds,color_by = 'SCT_snn_res.1'
                        ,cell_size =1,
                        # cols = c("#8C564BFF" ,"#E377C2FF", "#C49C94FF", "#BCBD22FF", "#17BECFFF",
                        #          "#1F77B4FF", "#FF7F0EFF", "#2CA02CFF" ,"#D62728FF", "#9467BDFF"
                        #          #"#AEC7E8FF", "#FFBB78FF", "#98DF8AFF", "#FF9896FF"
                        #          )
)+scale_color_npg()
p5=plot_cell_trajectory(cds,color_by = 'condition',
                        cols = c('#abd4de','#d95b48'),
                        cell_size =1)
pdf('E:/windows/AS/article/dataset/1.analysis/figure/xin_beta_monocle_1_variable.pdf',
    width=20,height = 5)
p1|p2|p3|p4|p5
dev.off()

pdf('E:/windows/AS/article/dataset/1.analysis/figure/xin_beta_monocle_2_subcluster_faced.pdf',
    width=15,height = 5)
plot_cell_trajectory(cds,color_by = 'SCT_snn_res.1',
                     cell_size =1,
                     # cols = c("#8C564BFF" ,"#E377C2FF", "#C49C94FF", "#BCBD22FF", "#17BECFFF",
                     #          "#1F77B4FF", "#FF7F0EFF", "#2CA02CFF" ,"#D62728FF", "#9467BDFF",
                     #          #"#AEC7E8FF", "#FFBB78FF", "#98DF8AFF", "#FF9896FF"
                     # )
)+scale_color_npg()+
  facet_wrap('~SCT_snn_res.1',nrow = 1)
dev.off()


p6=plot_complex_cell_trajectory(cds,x=1,y=2,
                                color_by = 'SCT_snn_res.1',
                                cell_size =1)+
  scale_color_npg()+
  theme(legend.title = element_blank())
p7=plot_complex_cell_trajectory(cds,x=1,y=2,
                                color_by = 'condition',
                                cell_size =1,
                                cols = c('#abd4de','#d95b48'))+
  theme(legend.title = element_blank())
pdf('E:/windows/AS/article/dataset/1.analysis/figure/xin_beta_monocle_3_subcluster_tree.pdf',
    width=5,height = 4)
p6|p7
dev.off()

saveRDS(cds,file = 'E:/windows/AS/article/dataset/xin/xin_cds.rds')



################################# mature score #################################
library(ggplot2)
library(grid)
library(gridExtra)
library(RColorBrewer)



aa_3<-list(c('INS','SLC2A2','MAFA','RFX6','PDX1','CHL1',
             'GCK','PPARGC1A','MDH1','NEUROD1',
             'CREB1','G6PC2','PFKFB2','PFKM','SIX2','SIX3',
             'ENTPD3','GPD2','DNMT3A','MTOR'))

xin_beta_sce <- AddModuleScore(  #对每个模块的基因集打分
  object = xin_beta_sce,
  features = aa_3,
  ctrl = 100, #默认值是100
  name = 'mature_score')

tsne_result.ch<-data.frame(xin_beta_sce@reductions[["umap"]]@cell.embeddings); tsne_result.ch[1:8,]
expr.ch<-data.frame(xin_beta_sce@meta.data[["mature_score1"]]);colnames(expr.ch)<-"mature_score"

rownames(expr.ch)<-rownames(tsne_result.ch)
expr.ch<-as.data.frame(t(expr.ch))


IsoheightMap <- function(tsne_result, gene, gene_expression){
  colnames(tsne_result)<- c("tSNE_1", "tSNE_2")
  warm <- colorRampPalette(c('#ffffb2','#fecc5c','#e31a1c','#800026'))
  mypalette <- c(warm(20))
  
  
  # pheatmap(aa,
  #          cluster_cols = F,
  #          cluster_rows = F,
  #          scale = "row",
  #          cellheight = 10,
  #          cellwidth = 16.18,
  #          border_color = "white",
  #          colorRampPalette(brewer.pal(11, "RdYlBu")[c(11:1)])(1000))
  # 
  # cold <- colorRampPalette(c('#f7fcf0','#41b6c4','#253494')) # blue
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
                    colour="black",size=0.7) +
    # scale_fill_gradient2(high="darkred", low="yellow")+
    scale_fill_gradientn(colours = mypalette)+
    scale_colour_gradientn(colours = mypalette)+
    # scale_fill_viridis_c()+
    theme_bw() +
    # ggtitle(gene) +
    annotation_custom(title) +
    xlab("UMAP_1") +
    ylab("UMAP_2") +
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

pdf('E:/windows/AS/article/dataset/1.analysis/figure/xin_beta_mature_score_IsoheightMap.pdf', width=5, height=5)
p <- list()
for (genes in rownames(expr.ch) ) {
  p[[genes]] <- IsoheightMap(
    tsne_result.ch,
    genes, 
    expr.ch)
}
do.call(grid.arrange, c(p, ncol=1)) #(3,3) #(15,3*)
dev.off()



################################## CytoTRACE ###################################

library(CytoTRACE)

dir.create('E:/windows/AS/article/dataset/xin/CytoTRACE/')
# T cell
####提取表型文件
phe <- Idents(xin_beta_sce)
phe = as.character(phe)
names(phe) <- rownames(xin_beta_sce@meta.data)
####提取表达矩阵
counts <- as.matrix(xin_beta_sce@assays$RNA@counts)
counts[1:4,1:4]
results <- CytoTRACE(mat = counts)
plotCytoGenes(results, numOfGenes = 10,outputDir = "E:/windows/AS/article/dataset/xin/CytoTRACE/"
)
plotCytoTRACE(results, phenotype = phe,gene = "INS",outputDir = "E:/windows/AS/article/dataset/xin/CytoTRACE/")


cyto<-read.table('E:/windows/AS/article/dataset/xin/CytoTRACE/CytoTRACE_plot_table.txt',
                 sep = '\t',header = T)
cyto_1<-as.data.frame(cyto[,1]) 
rownames(cyto_1) <-rownames(cyto)
colnames(cyto_1)<-'CytoTRACE'
xin_beta_sce<-AddMetaData(xin_beta_sce, cyto_1, col.name = NULL)
FeaturePlot(xin_beta_sce,features = 'CytoTRACE',label = T,
            reduction ='umap',pt.size = 1,
            cols=colorRampPalette(brewer.pal(11, "RdYlBu")[c(11:1)])(100))
ggsave('E:/windows/AS/article/dataset/xin/CytoTRACE/CytoTRACE_plot.pdf',width = 5,height = 4.5)


