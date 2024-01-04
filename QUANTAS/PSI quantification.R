getwd()
rm(list = ls()) 
options(stringsAsFactors = F) 
library(reticulate)
use_python("/python",required = T)
py_config()
py_module_available("umap")
library(BiocManager)
library(tidygraph)
library(clusterProfiler)
library(data.table)
library(reshape2)
library(dplyr)
library(Seurat)
library(ggplot2)

##### integrate data
#sed header
a=do.call(rbind,lapply(list.files('/sc_as/Lawlor.Genom.Res.2016/mis4_as/count_as_mis4_law/'),
                       function(x){ 
                         read.table(file.path('/sc_as/Lawlor.Genom.Res.2016/mis4_as/count_as_mis4_law/',x),
                                    sep = '\t')}))
law_mix=a
saveRDS(law_mix,file = 'law_mis4_mix_data.rds')
a<- readRDS("/sc_as/Lawlor.Genom.Res.2016/law_raw_mix_data.rds")
head(a)
##### calculate psi/junction reads/filter criteria
a$psiiso=a$V9/(a$V9+a$V10)
a$psijunc=(a$V12+a$V13)/(a$V12+a$V13+2*a$V14)
a$junc_cov=a$V12+a$V13+a$V14
a$junc1.2=a$V12/a$V13

head(a)
psiiso.data=a$psiiso
psijunc.data=a$psijunc
junc_cov.data=a$junc_cov
junc1.2.data=a$junc1.2
setwd('/sc_as/Lawlor.Genom.Res.2016/mis4_as/')
write.csv(psiiso.data,file='psi_iso.csv',sep = ',',quote = F,row.names = F)
write.csv(psijunc.data,file = 'psi_junc.csv',sep=',',quote = F,row.names = F)
write.csv(junc_cov.data,file = 'junc_cov.csv',sep=',',quote = F,row.names = F)
write.csv(junc1.2.data,file = 'junc1.2.data.csv',sep=',',quote = F,row.names = F)

##### make matrix
# pay attention to sed
file_split <- function(filename,eachfile_lines_num){ 
  c <- file(filename,"r") 
  varnames <- paste("splitfile", 1:1157, sep = "_")   
  i <- 1
  while(TRUE){  
    assign(varnames[i],value = readLines(c,n = eachfile_lines_num))  
    write.csv(get(varnames[i]),paste(varnames[i],".csv",sep = ""),row.names = F,
              quote = F) 
    if (length(get(varnames[i])) < eachfile_lines_num) break             
    else i <- i + 1 
  }  
  return(i)
}  
setwd('/sc_as/Lawlor.Genom.Res.2016/mis4_as/split_iso_csv/')
file_split('/sc_as/Lawlor.Genom.Res.2016/mis4_as/psi_iso.csv',42761)
setwd('/sc_as/Lawlor.Genom.Res.2016/mis4_as/split_junc_csv/')
file_split('/sc_as/Lawlor.Genom.Res.2016/mis4_as/psi_junc.csv',42761)
setwd('/sc_as/Lawlor.Genom.Res.2016/mis4_as/split_cov_csv/')
file_split('/sc_as/Lawlor.Genom.Res.2016/mis4_as/junc_cov.csv',42761)
setwd('/sc_as/Lawlor.Genom.Res.2016/mis4_as/split_junc1.2_csv/')
file_split('/sc_as/Lawlor.Genom.Res.2016/mis4_as/junc1.2.data.csv',42761)

#pay attention to sed
psiisodataset=do.call(cbind,lapply(list.files('/sc_as/Lawlor.Genom.Res.2016/mis4_as/split_iso_csv/'),
                                   function(x){ 
                                     read.table(file.path('/sc_as/Lawlor.Genom.Res.2016/mis4_as/split_iso_csv/',x),
                                                sep = ',')[,1]})) 
psijuncdataset=do.call(cbind,lapply(list.files('/sc_as/Lawlor.Genom.Res.2016/mis4_as/split_junc_csv/'),
                                    function(x){ 
                                      read.table(file.path('/sc_as/Lawlor.Genom.Res.2016/mis4_as/split_junc_csv/',x),
                                                 sep = ',')[,1]}))
junc_cov.data=do.call(cbind,lapply(list.files('/sc_as/Lawlor.Genom.Res.2016/mis4_as/split_cov_csv/'),
                                   function(x){ 
                                     read.table(file.path('/sc_as/Lawlor.Genom.Res.2016/mis4_as/split_cov_csv/',x),
                                                sep = ',' )[,1]}))
junc1.2.data=do.call(cbind,lapply(list.files('/sc_as/Lawlor.Genom.Res.2016/mis4_as/split_junc1.2_csv/'),
                                  function(x){ 
                                    read.table(file.path('/sc_as/Lawlor.Genom.Res.2016/mis4_as/split_junc1.2_csv/',x),
                                               sep = ',' )[,1]}))
#gene expression
lawgeneexpr=do.call(cbind,lapply(list.files('/sc_as/Lawlor.Genom.Res.2016/mis4_as/gene_mis4_law/'),function(x){ 
  read.table(file.path('/sc_as/Lawlor.Genom.Res.2016/mis4_as/gene_mis4_law/',x),
             header = F,sep = '\t')[,3]}))
#each junc reads
junc_tag1=do.call(cbind,lapply(list.files('/sc_as/Lawlor.Genom.Res.2016/mis4_as/count_as_mis4_law/'),
                               function(x){ 
                                 read.table(file.path('/sc_as/Lawlor.Genom.Res.2016/mis4_as/count_as_mis4_law/',x),
                                            sep = '\t')[,12]}))
junc_tag2=do.call(cbind,lapply(list.files('/sc_as/Lawlor.Genom.Res.2016/mis4_as/count_as_mis4_law/'),
                               function(x){ 
                                 read.table(file.path('/sc_as/Lawlor.Genom.Res.2016/mis4_as/count_as_mis4_law/',x),
                                            sep = '\t')[,13]}))
skip_junc=do.call(cbind,lapply(list.files('/sc_as/Lawlor.Genom.Res.2016/mis4_as/count_as_mis4_law/'),
                               function(x){ 
                                 read.table(file.path('/sc_as/Lawlor.Genom.Res.2016/mis4_as/count_as_mis4_law/',x),
                                            sep = '\t')[,14]}))
setwd('/sc_as/Lawlor.Genom.Res.2016/mis4_as/matrix/')
write.csv(junc_tag1,file = "junc_tag1.csv",sep=',',quote = F,row.names = F)
write.csv(junc_tag2,file = "junc_tag2.csv",sep=',',quote = F,row.names = F)
write.csv(skip_junc,file = "skip_junc.csv",sep=',',quote = F,row.names = F)
write.csv(junc1.2.data,file = "junc1.2.data.csv",sep=',',quote = F,row.names = F)

#colname
fs=list.files('/sc_as/Lawlor.Genom.Res.2016/mis4_as/gene_mis4_law/')
colnames(psiisodataset)=gsub('.txt','',
                             substring(fs,1,nchar("SRR3617192")))
psiisodataset[1:4,1:4]
colnames(lawgeneexpr)=gsub('.txt','',
                           substring(fs,1,nchar("SRR3617192")))
lawgeneexpr[1:4,1:4] 
colnames(psijuncdataset)=gsub('.txt','',
                              substring(fs,1,nchar("SRR3617192")))
psijuncdataset[1:4,1:2] 
colnames(junc_cov.data)=gsub('.txt','',
                             substring(fs,1,nchar("SRR3617192")))
junc_cov.data[1:4,1:4]

#rowname
a=read.table('/sc-AS/diab_pancr/Lawlor.Genom.Res.2016/count_AS_m4/SRR3617192/cass.count.txt', 
             header = F,sep = '\t')
b=read.table('/sc-AS/diab_pancr/Lawlor.Genom.Res.2016/gene_exp_m4/SRR3617192.txt',
             header=F,sep = '\t')
lawgeneexpr=as.data.frame(lawgeneexpr)
psiisodataset=as.data.frame(psiisodataset)
psijuncdataset=as.data.frame(psijuncdataset)
junc_cov.data=as.data.frame(junc_cov.data)

psi_iso_bulk_data=psiisodataset[,2:25]
psi_iso_bulk_data[1:4,1:4]
psi_iso_sc_data=psiisodataset[,c(1,26:1157)]
psi_junc_bulk_data=psijuncdataset[,2:25]
psi_junc_sc_data=psijuncdataset[,c(1,26:1157)]
junc_cov_bulk_data=junc_cov.data[,2:25]
junc_cov_sc_data=junc_cov.data[,c(1,26:1157)]

cass_info=a[,c(1:4,6:8)]
head(cass_info)
library(dplyr)
colnames(cass_info)=c('chrom','chromStart','chromEnd','name','strand','type','isoformIDs')
psi_iso_bulk_data=bind_cols(cass_info,psi_iso_bulk_data)
psi_iso_sc_data=bind_cols(cass_info,psi_iso_sc_data)
psi_junc_bulk_data=bind_cols(cass_info,psi_junc_bulk_data)
psi_junc_sc_data=bind_cols(cass_info,psi_junc_sc_data)

gene_bulk_data=lawgeneexpr[,2:25]
gene_sc_data=lawgeneexpr[,c(1,26:1157)]
gene_bulk_data$gene_symbol=b$V2
gene_sc_data$gene_symbol=b$V2

setwd('/sc_as/Lawlor.Genom.Res.2016/mis4_as/matrix/')
write.csv(psi_iso_bulk_data,file = "psi_iso_bulk_data.csv",sep=',',quote = F,row.names = F)
write.csv(psi_iso_sc_data,file = "psi_iso_sc_data.csv",sep=',',quote = F,row.names = F)
write.csv(psi_junc_bulk_data,file = "psi_junc_bulk_data.csv",sep=',',quote = F,row.names = F)
write.csv(psi_junc_sc_data,file = "psi_junc_sc_data.csv",sep=',',quote = F,row.names = F)
write.csv(gene_bulk_data,file = "gene_bulk_data.csv",sep=',',quote = F,row.names = F)
write.csv(gene_sc_data,file = "gene_sc_data.csv",sep=',',quote = F,row.names = F)
write.csv(junc_cov_sc_data,file = "junc_cov_sc_data.csv",sep=',',quote = F,row.names = F)
write.csv(junc_cov_bulk_data,file = "junc_cov_bulk_data.csv",sep=',',quote = F,row.names = F)

junc_cov_sc_data$sum=apply(junc_cov_sc_data,1,sum)
rownames(junc_cov_sc_data)=cellinfo$name
new_junc_cov=subset(junc_cov_sc_data,junc_cov_sc_data$sum>19)

junc_cov.data=read.csv('law_data_matrix/junc_cov_sc_data.csv',header = T,row.names = 1)

table(apply(junc_cov.data,1,function(x) sum(x>19)>225))
new_junc_cov$name=rownames(new_junc_cov)
new_psi_iso_sc=subset(psi_iso_sc_data,psi_iso_sc_data$name%in%new_junc_cov$name)
rownames(new_psi_iso_sc)=new_psi_iso_sc$name
# filter junction coverage
col1=colnames(new_junc_cov)
col2=colnames(new_psi_iso_sc)
col3=intersect(col1,col2)
new_psi_iso_sc=subset(new_psi_iso_sc,select = col3)
getwd()
setwd('/sc_as/Lawlor.Genom.Res.2016/law_data_matrix/')
write.csv(new_psi_iso_sc,file='fil_psi_iso_data.csv')
write.csv(junc_cov_bulk_data,file='junc_cov_bulk_data.csv')
write.csv(junc_cov_sc_data,file='junc_cov_sc_data.csv')
