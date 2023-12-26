violin_gene_exp<-function(gene, rpkm, conditions=conditions,  test=TRUE){
  # exp<-as.numeric(log(rpkm[gene,]+1))
  exp<-as.numeric(rpkm[gene,])
  gene_exp<-data.frame(
    cell=colnames(rpkm),
    clusters=conditions,
    gene=exp
  )
  
  if(test==TRUE){
    #Perform pairwise comparisons
    test<-compare_means(gene~clusters, data=gene_exp, method="wilcox.test")
    #print(test)
    #print(compare_means(gene ~ clusters, data = gene_exp, , method = "wilcox.test", ref.group = ".all."))
    
    p<-ggboxplot(gene_exp, x="clusters", y="gene", color="white")+
      # geom_violin(scale="width", width=0.7, adjust=.5,aes(fill=clusters)) +
      stat_summary(fun.y=mean, geom="point", shape=21, size=3,  stroke=1, fill="white")+
      geom_boxplot(fill=rep(c("#add8e6","#ee5e49"),2), outlier.shape = NA, width = 0.5)+ #箱子的颜色
      
      geom_jitter(size=0.3)+  	  
      # geom_hline(yintercept=mean(gene_exp$gene), linetype=2)+
      # geom_hline(yintercept=mean(gene_exp[494:542,]$gene), linetype=2)+  #虚线设为 T-HSC的中值
      
      # scale_fill_aaas()+
      # scale_fill_manual(values = c('#bb0021','#008280','#008b45','#631879','#5f559b','#a20056','#3b4992','#ee0000'))+
      scale_fill_manual(values =rep(c("#add8e6","#ee5e49"),2))+
      theme_bw()+
      ggtitle(gene)+
      expand_limits(y=c(0, max(gene_exp$gene)*1.1))+
      # scale_x_discrete(labels = c('T-HSC','LT-HSC','CLP','CMP','Tcm(LN)','Tcm(SP)','CD4T','CD8T')) +
      stat_compare_means(method = "wilcox.test", label.y = max(gene_exp$gene)+1,label.x = c(1,3,5,7))+      # Add global p-value
      stat_compare_means(
        label="p.signif",  #p.format
        method="wilcox.test",
        ref.group="G7",    #  .all.
        label.y=max(gene_exp$gene)*1.05,
        size=6
      )+#Pairwise comparison against all
      theme(
        plot.title=element_text(size=18, face="bold.italic", hjust=0.5),
        axis.text=element_text(size=6),
        #axis.title=element_text(size=16),
        axis.title=element_blank(),
        legend.text=element_text(size=16),
        legend.title=element_blank(),
        aspect.ratio=0.5,
        legend.position="none",
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank()
      )
    
  }else{
    p<-ggplot(gene_exp, aes(clusters, gene))+
      geom_violin(scale="width", width=0.7, adjust=.5,aes(fill=clusters)) +
      stat_summary(fun.y=mean, geom="point", shape=21, size=3,  stroke=1, fill="white")+
      #geom_jitter(size=0.3)+ 
      scale_fill_manual(
        values=colours
      ) +
      theme_bw()+
      ggtitle(gene)+
      theme(
        plot.title=element_text(size=18, face="bold.italic", hjust=0.5),
        axis.text=element_text(size=16),
        #axis.title=element_text(size=16),
        axis.title=element_blank(),
        legend.text=element_text(size=16),
        legend.title=element_blank(),
        aspect.ratio=0.5,
        legend.position="none"
        #panel.grid.major = element_blank(),
        #panel.grid.minor = element_blank()
      )
  }
  
  
  print(p)
}

