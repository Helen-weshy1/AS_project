exon_matrix=read.table('/Quantas/index_annot/hg19/annotation/hg19.exon.trio.hmr.nr.bed',
                       sep = '\t')


# alpha
alpha_DE=read.csv('/alpha_diff/alpha_allDE_dI0.1_ratio0.1.csv',
                  header = T)
head(exon_matrix)

can_exon_name=str_split_fixed(alpha_DE$name,'-',6)
head(can_exon_name)
class(can_exon_name)
can_exon_name=as.data.frame(can_exon_name)
colnames(can_exon_name)=paste0('c',1:6)
head(can_exon_name)
str_name=str_split_fixed(can_exon_name$c6,'[INC]',2)
head(str_name)
str_name=as.data.frame(str_name)
can_exon_name$c7=can_exon_name$c6
can_exon_name$c6=str_name$V1
head(can_exon_name)

all_exon_id=str_split_fixed(exon_matrix$V4,'[.]',2)%>%as.data.frame()
head(all_exon_id)
all_exon_id=str_split_fixed(all_exon_id$V1,'-',4)%>%as.data.frame()
head(all_exon_id)
all_exon_id=cbind(exon_matrix[,1:6],all_exon_id)
colnames(all_exon_id)=c('chr','start','end','name','p','strand','type','gene','exonstart','exonend')
all_exon_id=all_exon_id[all_exon_id$type=='EX',]
head(all_exon_id)
all_exon_id$exon_l=paste0('EX-',all_exon_id$gene,'-',all_exon_id$exonend)
all_exon_id$exon_m=paste0('EX-',all_exon_id$gene,'-',all_exon_id$exonstart,'-',all_exon_id$exonend)
all_exon_id$exon_r=paste0('EX-',all_exon_id$gene,'-',all_exon_id$exonstart,'[')
head(all_exon_id)

head(can_exon_name)
can_exon_name$exon_l=paste0('EX-',can_exon_name$c2,'-',can_exon_name$c3,'[')
can_exon_name$exon_m=paste0('EX-',can_exon_name$c2,'-',can_exon_name$c4,'-',can_exon_name$c5,'[')
can_exon_name$exon_r=paste0('EX-',can_exon_name$c2,'-',can_exon_name$c6)
head(can_exon_name)

exon_l_site=left_join(can_exon_name,all_exon_id[,c(1:3,11)],by='exon_l')
exon_m_site=left_join(can_exon_name,all_exon_id[,c(1:3,12)],by='exon_m')
exon_r_site=left_join(can_exon_name,all_exon_id[,c(6,1:3,13)],by='exon_r')
head(exon_l_site)
head(exon_m_site)
head(exon_r_site)
exon_site=cbind(exon_r_site,exon_m_site[,12:13],exon_l_site[,12:13])
head(exon_site)
exon_site=na.omit(exon_site)
colnames(exon_site)[c(2,13:18)]=c('ensemble','firstFlankingExonStart','firstFlankingExonEnd','exonStart','exonEnd','secondFlankingExonStart','secondFlankingExonEnd')
alpha_DE_exon_site=exon_site[,c(12,11,15,16,13,14,17,18)]
head(alpha_DE_exon_site)
write.table(alpha_DE_exon_site,file = '/alpha_diff/law_alpha_DE_exon_site_for_rMAPS.txt',
            sep = '\t',col.names = F,row.names = F,quote = F)

# beta
beta_DE=read.csv('/beta_diff/beta_allDE_dI0.1_ratio0.1.csv',
                 header = T)

can_exon_name=str_split_fixed(beta_DE$name,'-',6)
head(can_exon_name)
class(can_exon_name)
can_exon_name=as.data.frame(can_exon_name)
colnames(can_exon_name)=paste0('c',1:6)
head(can_exon_name)
str_name=str_split_fixed(can_exon_name$c6,'[INC]',2)
head(str_name)
str_name=as.data.frame(str_name)
can_exon_name$c7=can_exon_name$c6
can_exon_name$c6=str_name$V1
head(can_exon_name)

can_exon_name$exon_l=paste0('EX-',can_exon_name$c2,'-',can_exon_name$c3,'[')
can_exon_name$exon_m=paste0('EX-',can_exon_name$c2,'-',can_exon_name$c4,'-',can_exon_name$c5,'[')
can_exon_name$exon_r=paste0('EX-',can_exon_name$c2,'-',can_exon_name$c6)
head(can_exon_name)

exon_l_site=left_join(can_exon_name,all_exon_id[,c(1:3,11)],by='exon_l')
exon_m_site=left_join(can_exon_name,all_exon_id[,c(1:3,12)],by='exon_m')
exon_r_site=left_join(can_exon_name,all_exon_id[,c(6,1:3,13)],by='exon_r')
head(exon_l_site)
head(exon_m_site)
head(exon_r_site)
exon_site=cbind(exon_r_site,exon_m_site[,12:13],exon_l_site[,12:13])
head(exon_site)
exon_site=na.omit(exon_site)
colnames(exon_site)[c(2,13:18)]=c('ensemble','firstFlankingExonStart','firstFlankingExonEnd','exonStart','exonEnd','secondFlankingExonStart','secondFlankingExonEnd')
beta_DE_exon_site=exon_site[,c(12,11,15,16,13,14,17,18)]
head(beta_DE_exon_site)
write.table(beta_DE_exon_site,file = '/beta_diff/law_beta_DE_exon_site_for_rMAPS.txt',
            sep = '\t',col.names = F,row.names = F,quote = F)

# beta overlap
beta_DE=read.csv('/media/user/sdd/zhaolabguest/wsy/sc_as/lawlor.1050/diff/beta_diff/overlap_beta_allDE_dI0.1_ratio0.1.csv',
                 header = T)

can_exon_name=str_split_fixed(beta_DE$name,'-',6)
head(can_exon_name)
class(can_exon_name)
can_exon_name=as.data.frame(can_exon_name)
colnames(can_exon_name)=paste0('c',1:6)
head(can_exon_name)
str_name=str_split_fixed(can_exon_name$c6,'[INC]',2)
head(str_name)
str_name=as.data.frame(str_name)
can_exon_name$c7=can_exon_name$c6
can_exon_name$c6=str_name$V1
head(can_exon_name)

can_exon_name$exon_l=paste0('EX-',can_exon_name$c2,'-',can_exon_name$c3,'[')
can_exon_name$exon_m=paste0('EX-',can_exon_name$c2,'-',can_exon_name$c4,'-',can_exon_name$c5,'[')
can_exon_name$exon_r=paste0('EX-',can_exon_name$c2,'-',can_exon_name$c6)
head(can_exon_name)

exon_l_site=left_join(can_exon_name,all_exon_id[,c(1:3,11)],by='exon_l')
exon_m_site=left_join(can_exon_name,all_exon_id[,c(1:3,12)],by='exon_m')
exon_r_site=left_join(can_exon_name,all_exon_id[,c(6,1:3,13)],by='exon_r')
head(exon_l_site)
head(exon_m_site)
head(exon_r_site)
exon_site=cbind(exon_r_site,exon_m_site[,12:13],exon_l_site[,12:13])
head(exon_site)
exon_site=na.omit(exon_site)
colnames(exon_site)[c(2,13:18)]=c('ensemble','firstFlankingExonStart','firstFlankingExonEnd','exonStart','exonEnd','secondFlankingExonStart','secondFlankingExonEnd')
beta_DE_exon_site=exon_site[,c(12,11,15,16,13,14,17,18)]
head(beta_DE_exon_site)
write.table(beta_DE_exon_site,file = '/media/user/sdd/zhaolabguest/wsy/sc_as/lawlor.1050/diff/beta_diff/overlap_beta_DE_exon_site_for_rMAPS.txt',
            sep = '\t',col.names = F,row.names = F,quote = F)


############## site_func
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

alpha_up_DE=read.csv('/media/user/sdd/zhaolabguest/wsy/sc_as/lawlor.1050/diff/alpha_diff/alpha_upDE_dI0.1_ratio0.1.csv',
                     header = T)
alpha_down_DE=read.csv('/media/user/sdd/zhaolabguest/wsy/sc_as/lawlor.1050/diff/alpha_diff/alpha_downDE_dI0.1_ratio0.1.csv',
                       header = T)
beta_up_DE=read.csv('/media/user/sdd/zhaolabguest/wsy/sc_as/lawlor.1050/diff/beta_diff/beta_upDE_dI0.1_ratio0.1.csv',
                    header = T)
beta_down_DE=read.csv('/media/user/sdd/zhaolabguest/wsy/sc_as/lawlor.1050/diff/beta_diff/beta_downDE_dI0.1_ratio0.1.csv',
                      header = T)
alpha_up_DE_site=site_func(alpha_up_DE,all_exon_id)
alpha_down_DE_site=site_func(alpha_down_DE,all_exon_id)
beta_up_DE_site=site_func(beta_up_DE,all_exon_id)
beta_down_DE_site=site_func(beta_down_DE,all_exon_id)
write.table(alpha_up_DE_site,file = '/media/user/sdd/zhaolabguest/wsy/sc_as/lawlor.1050/diff/alpha_diff/law_alpha_up_DE_exon_site_for_rMAPS.txt',
            sep = '\t',col.names = F,row.names = F,quote = F)
write.table(alpha_down_DE_site,file = '/media/user/sdd/zhaolabguest/wsy/sc_as/lawlor.1050/diff/alpha_diff/law_alpha_down_DE_exon_site_for_rMAPS.txt',
            sep = '\t',col.names = F,row.names = F,quote = F)
write.table(beta_up_DE_site,file = '/media/user/sdd/zhaolabguest/wsy/sc_as/lawlor.1050/diff/beta_diff/law_beta_up_DE_exon_site_for_rMAPS.txt',
            sep = '\t',col.names = F,row.names = F,quote = F)
write.table(beta_down_DE_site,file = '/media/user/sdd/zhaolabguest/wsy/sc_as/lawlor.1050/diff/beta_diff/law_beta_down_DE_exon_site_for_rMAPS.txt',
            sep = '\t',col.names = F,row.names = F,quote = F)

backg_exon=read.csv('/media/user/sdd/zhaolabguest/wsy/sc_as/lawlor.1050/diff/beta_diff/law_backg_beta.csv',
                    header = T)
backg_exon_site=site_func(backg_exon,all_exon_id)
write.table(backg_exon_site,file = '/media/user/sdd/zhaolabguest/wsy/sc_as/lawlor.1050/diff/beta_diff/law_beta_background_exon_site_for_rMAPS.txt',
            sep = '\t',col.names = F,row.names = F,quote = F)

### all exon site
law_endo_c1_as=read.table(paste0('/media/user/sdd/zhaolabguest/wsy/sc_as/lawlor.1050/diff/diff_result/endo_as_cluster/endo_c',1,'.diff.txt'),
                          sep = '\t',header = T)
dim(law_endo_c1_as)
can_exon_name=str_split_fixed(law_endo_c1_as$name,'-',6)%>%as.data.frame()
colnames(can_exon_name)=paste0('c',1:6)
str_name=str_split_fixed(can_exon_name$c6,'[INC]',2)%>%as.data.frame()
can_exon_name$c7=can_exon_name$c6
can_exon_name$c6=str_name$V1
can_exon_name$name=law_endo_c1_as$name
can_exon_name$exon_l=paste0('EX-',can_exon_name$c2,'-',can_exon_name$c3,'[')
can_exon_name$exon_m=paste0('EX-',can_exon_name$c2,'-',can_exon_name$c4,'-',can_exon_name$c5,'[')
can_exon_name$exon_r=paste0('EX-',can_exon_name$c2,'-',can_exon_name$c6)

exon_l_site=left_join(can_exon_name,all_exon_id[,c(1:3,11)],by='exon_l')
exon_m_site=left_join(can_exon_name,all_exon_id[,c(1:3,12)],by='exon_m')
exon_r_site=left_join(can_exon_name,all_exon_id[,c(6,1:3,13)],by='exon_r')
exon_site=cbind(exon_r_site,exon_m_site[,13:14],exon_l_site[,13:14])%>%na.omit()
colnames(exon_site)
head(exon_site)
colnames(exon_site)[c(2,14:19)]=c('ensemble','firstFlankingExonStart','firstFlankingExonEnd','exonStart','exonEnd','secondFlankingExonStart','secondFlankingExonEnd')
all_exon_site=exon_site[,c(2,8,13,12,16,17,14,15,18,19)]
head(all_exon_site)
write.csv(all_exon_site,file = '/media/user/sdd/zhaolabguest/wsy/sc_as/lawlor.1050/diff/endo_diff_cluster/data/all_exon_site.csv',
          quote = F,row.names = F)

can_exon_name=str_split_fixed(law_endo_c1_as$name,'-',6)%>%as.data.frame()
colnames(can_exon_name)=paste0('c',1:6)
str_name=str_split_fixed(can_exon_name$c6,'[INC]',2)%>%as.data.frame()
can_exon_name$c7=can_exon_name$c6
can_exon_name$c6=str_name$V1
can_exon_name$name=law_endo_c1_as$name
can_exon_name$exon_m=paste0('EX-',can_exon_name$c2,'-',can_exon_name$c4,'-',can_exon_name$c5,'[')

exon_m_site=left_join(can_exon_name,all_exon_id[,c(1:3,12)],by='exon_m')%>%na.omit()
head(exon_m_site)
colnames(exon_m_site)[c(2,13,14)]=c('ensemble','exonStart_alone','exonEnd_alone')
all_exon_site=exon_m_site[,c(2,8,12,13,14)]
head(all_exon_site)
write.csv(all_exon_site,file = '/media/user/sdd/zhaolabguest/wsy/sc_as/lawlor.1050/diff/endo_diff_cluster/data/all_exon_site_single_exon.csv',
          quote = F,row.names = F)
