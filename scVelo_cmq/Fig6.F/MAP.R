.libPaths("/media/user/sdf/R.lib/library/wj/reg")

library(ggplot2)
library(grid)
library(gridExtra)
library(RColorBrewer)
# IsoheightMap

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