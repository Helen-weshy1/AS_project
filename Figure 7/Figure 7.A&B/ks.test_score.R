.libPaths("/media/user/sdf/R.lib/library/wj/sc")
library(Seurat)
library(ggplot2)
setwd('/media/user/sdh/cmq_wj/others/vlnplot')

ss_obj<-readRDS('20GROUP_ND_T2D.rds')
i=1
ND_SCORE<-list()
T2D_SCORE<-list()

for( name in c("MURARO_PANCREAS_ALPHA_CELL.v7.5.11",
               'MURARO_PANCREAS_BETA_CELL.v7.5.11')){
  for (nn in c(1:10)){ 
    
    for (group in c('ND','T2D')){
      
      if(group=='ND'){
        name1 <- c(rownames(ss_obj@meta.data)[grep(group,ss_obj$Status)])
        ss_obj_ND <- subset(ss_obj,cells = name1)
        
        ND.name <- c(rownames(ss_obj_ND@meta.data)[grep(nn,ss_obj_ND$new_as_cluster )])
        ND.umap <- subset(ss_obj_ND,cells = ND.name)
        
        ND_RNA1 <-ND.umap@meta.data[[name]]
      }
      
      if(group=='T2D'){
        name1 <- c(rownames(ss_obj@meta.data)[grep(group,ss_obj$Status)])
        ss_obj_T2D <- subset(ss_obj,cells = name1)
        
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
View(p_value_dataframe)

write.csv(p_value_dataframe,'p_value_score.csv')