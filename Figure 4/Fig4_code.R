library(data.table)
library(dplyr)
library(stringr)


# panel A -----------------------------------------------------------------
# volcano
## lawlor
law_beta_tvn_as
fm_law_beta_tvn_as=law_beta_tvn_as[law_beta_tvn_as$beta_T2D_ratio>0.1|
                                     law_beta_tvn_as$beta_ND_ratio>0.1,]
fm_law_beta_tvn_as$mark=ifelse(fm_law_beta_tvn_as$coverage>=20&
                                 abs(fm_law_beta_tvn_as$dI_g1_vs_g2)>=0.1&
                                 fm_law_beta_tvn_as$FDR<=0.05,
                               ifelse(fm_law_beta_tvn_as$dI_g1_vs_g2>0.1,'Up','Down'),
                               'None')
p <- ggplot(
  fm_law_beta_tvn_as, aes(x = dI_g1_vs_g2, y = -log10(FDR), colour=mark)) +
  geom_point( alpha=0.6,size=3.5) +
  scale_color_manual(values=c("#546de5", "#d2dae2","#ff4757"))+
  geom_vline(xintercept=c(-0.1,0.1),lty=1,col="black",lwd=0.8) +
  geom_hline(yintercept = -log10(0.05),lty=1,col="black",lwd=0.8) +
  labs(x="delta PSI",
       y="-log10 FDR",
       title='law.beta.DE')+
  xlim(-1,1)+
  theme_bw()+theme(panel.grid=element_blank())+
  theme(plot.title = element_text(hjust = 0.5), 
        legend.position="right", 
        legend.title = element_blank(),
        axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15))+
  geom_text_repel(data = fm_law_beta_tvn_as[(fm_law_beta_tvn_as$coverage>=20&
                                               abs(fm_law_beta_tvn_as$dI_g1_vs_g2)>=0.7&
                                               fm_law_beta_tvn_as$FDR<=0.05)|
                                              (fm_law_beta_tvn_as$FDR<=10^(-120)&
                                                 fm_law_beta_tvn_as$coverage>=20&
                                                 abs(fm_law_beta_tvn_as$dI_g1_vs_g2)>=0.1),],
                  size = 3,color='black',
                  box.padding = unit(0.5, "lines"),aes(label = symbol) ,
                  point.padding = unit(0.2, "lines"), segment.color = "black", show.legend = FALSE )
ggsave(p,filename = 'vol_law_beta.pdf',width=5,height=4)


# panel B -----------------------------------------------------------------
# GO
ensemble.id=fm_law_beta_tvn_as[fm_law_beta_tvn_as$mark=='Up'|
                                 fm_law_beta_tvn_as$mark=='Down',]$symbol
law_beta_id=bitr(ensemble.id,fromType = 'SYMBOL',
                 toType = c('ENSEMBL','ENTREZID'),
                 OrgDb = 'org.Hs.eg.db')
ego <- enrichGO(gene          =unique(law_beta_id$ENTREZID) ,
                OrgDb         = org.Hs.eg.db,
                pAdjustMethod = "BH",
                ont           = 'ALL' ,
                pvalueCutoff  = 0.05,
                #qvalueCutoff  = 0.99,
                readable      = TRUE)
law_beta_ego=ego
write.table(law_beta_ego@result, 'law_beta_GO_result.csv', row.names = FALSE, sep = '\t', quote = FALSE)

dat=law_beta_ego@result[order(law_beta_ego@result$p.adjust,decreasing = F),]
dat=dat[1:20,]
dat=dat[order(dat$p.adjust,decreasing = T),]
dat$Description=factor(dat$Description,levels = dat$Description)
dat=dat[c(1:5,9,10,13,20),]
ggplot(dat, aes(x=Description, y=-log10(p.adjust), fill=Count)) + 
  geom_bar(stat="identity")+
  scale_fill_gradientn(colours=c("#9098A0","#284048"))+
  coord_flip()+
  theme_bw()+theme(panel.grid=element_blank())+
  xlab("GO term") +
  theme(axis.text = element_text(size = 30))+
  ggtitle("lawlor beta DE GO")
ggsave('law_beta_GO_bar.pdf',width = 11,height = 5)


# panel C -----------------------------------------------------------------
law_endo_filter
xin_filter
per_law_fm_beta=left_join(law_fm_beta,law_endo_filter[,c(5,18,19)],by='name')
per_xin_fm_beta=left_join(xin_fm_beta,xin_filter[,c(5,40,41)],by='name')
per_law_fm_beta=subset(per_law_fm_beta,per_law_fm_beta$beta_T2D_ratio>0.1|
                         per_law_fm_beta$beta_ND_ratio>0.1)
per_xin_fm_beta=subset(per_xin_fm_beta,per_xin_fm_beta$beta_T2D_ratio>0.1|
                         per_xin_fm_beta$beta_ND_ratio>0.1)
venn.diagram(list(law.beta.DE=unique(per_law_fm_beta[per_law_fm_beta$dI_g1_vs_g2>0,]$name),
                  xin.beta.DE=unique(per_xin_fm_beta[per_xin_fm_beta$dI_g1_vs_g2>0,]$name)),
             resolution = 300, imagetype = "png", 
             alpha=c(0.5,0.5),
             #cat.fontface=4,fontfamily=3,
             fill=c("#F8B098",'#F0C050'), 
             cat.cex=1,
             cex=2,
             height = 1500, width =1500,
             cat.col = c("#F8B098",'#F0C050'),
             cat.fontfamily = "serif",
             cat.fontface = "bold",
             margin = 0.05,
             cat.dist = c(0.09,0.09),
             filename = "Venn_beta_xin_up_overlap2.png")
venn.diagram(list(law.beta.DE=unique(per_law_fm_beta[per_law_fm_beta$dI_g1_vs_g2<0,]$name),
                  xin.beta.DE=unique(per_xin_fm_beta[per_xin_fm_beta$dI_g1_vs_g2<0,]$name)),
             resolution = 300, imagetype = "png", 
             alpha=c(0.5,0.5),
             #cat.fontface=4,fontfamily=3,
             fill=c("#F8B098",'#F0C050'), 
             cat.cex=1,
             cex=2,
             height = 1500, width =1500,
             cat.col = c("#F8B098",'#F0C050'),
             cat.fontfamily = "serif",
             cat.fontface = "bold",
             margin = 0.05,
             cat.dist = c(0.09,0.09),
             filename = "Venn_beta_xin_down_overlap2.png")


# panel F -----------------------------------------------------------------
# APTX  chr9:32985969:32986028:-@chr9:32984629:32984855:-@chr9:32974455:32974559 strand="negative"
PlotPSI(
  tran_id="chr9:32985969:32986028:-@chr9:32984629:32984855:-@chr9:32974455:32974559",
  event.type="SE",
  strand="negative",
  Bam="/lawlor/marvel/bam/",
  BamPheno=BamPheno,
  cell.types=c( "T2D","ND"),
  min.coverage=5,
  cons.exon.cutoff=100,
  method="ks",
  method.adj="bonferroni",
  cell.types.colors="ggplot.default",
  plot.title="APTX",
  plot.width=5,
  plot.height=8,
  plot.out="APTX.pdf"
)


# panel G -----------------------------------------------------------------
# EXOSC3  chr9:37783911:37784060:-@chr9:37781983:37782134:-@chr9:37779711:37780877 strand="negative"
PlotPSI(
  tran_id="chr9:37783911:37784060:-@chr9:37781983:37782134:-@chr9:37779711:37780877",
  event.type="SE",
  strand="negative",
  Bam="/lawlor/marvel/bam/",
  BamPheno=BamPheno,
  cell.types=c( "T2D","ND"),
  min.coverage=5,
  cons.exon.cutoff=100,
  method="ks",
  method.adj="bonferroni",
  cell.types.colors="ggplot.default",
  plot.title="EXOSC3",
  plot.width=5,
  plot.height=8,
  plot.out="EXOSC3.pdf"
)


# panel H -----------------------------------------------------------------
Beta_tvn_DE_FC=left_join(fm_law_Beta_tvn_as[,c(1,13)],
                     law_Beta_tvng_fc[law_Beta_tvng_fc$p_val<=0.05,
                                      c(2,7)],
                     by='gene')
ggplot(Beta_tvn_DE_FC,aes(x=avg_logFC,y=dI_g1_vs_g2))+
  geom_point(size=1,color='#162B87')+
  ggtitle( 'Beta')+
  theme_test()
