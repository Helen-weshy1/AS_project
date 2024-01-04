library(ggplot2)
library(dplyr)
library(Seurat)



# panel E, F --------------------------------------------------------------
law_beta_gene <- readRDS("/law_beta_gene.rds")
xin_beta_gene<- readRDS("/xin_beta_gene.rds")


plot.features=c('HNRNPH2',
                'FXR1',
                'FMR1')
plots=list()
for (i in seq_along(plot.features)){
  plots[[i]]=VlnPlot(law_sc_gene_endo,pt.size = 0.1,
                     cols = c('#abd4de','#d95b48'),group.by = 'Status',
                     features = plot.features[i])+
    stat_summary(fun.y="mean")+
    NoLegend()
}
violin=wrap_plots(plots=plots,nrow = 1)
ggsave('/law_beta_tvn_RBP_vlnplot.pdf',
       plot = violin,width=7,height = 3)


# panel H, J --------------------------------------------------------------
ST7_PSI=data.frame(PSI=c(0.303370787,0.803837953,
                         0.654676259,0.76969697),
                   dataset=c('Lawlor','Lawlor','Xin','Xin'),
                   Status=c('T2D','ND','T2D','ND'))
ggplot(ST7_PSI, aes(x=dataset, y=PSI, fill=Status)) + 
  geom_bar(stat="identity", color="black", 
           position=position_dodge(),width = 0.5) +
  # geom_errorbar(aes(ymin=value-sd, ymax=value+sd), width=.2,
  #               position=position_dodge(0.9)) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45,hjust=1, vjust=1) )+
  scale_fill_manual(values=c('#a8ced9','#ce5a46'))
ggsave('/law_xin_ST7_bar.pdf',
       width=4,height=3)

TMEM222_PSI=data.frame(PSI=c(0.01734104,0.205741627,
                             0.008547009,0.204395604),
                       dataset=c('Lawlor','Lawlor','Xin','Xin'),
                       Status=c('T2D','ND','T2D','ND'))
ggplot(TMEM222_PSI, aes(x=dataset, y=PSI, fill=Status)) + 
  geom_bar(stat="identity", color="black", 
           position=position_dodge(),width = 0.5) +
  # geom_errorbar(aes(ymin=value-sd, ymax=value+sd), width=.2,
  #               position=position_dodge(0.9)) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45,hjust=1, vjust=1) )+
  scale_fill_manual(values=c('#a8ced9','#ce5a46'))
ggsave('/law_xin_TMEM222_bar.pdf',
       width=4,height=3)
