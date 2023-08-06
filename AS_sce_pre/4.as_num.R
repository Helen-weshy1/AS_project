
library(data.table)

path <- "/media/user/sdg/WShi/AS/STAR_map/bam2/media/user/sdf/temp_wj/sc-AS/diab_pancr/Lawlor.Genom.Res.2016/trim/trimmed.fq/"
file <- list.files("/media/user/sdf/temp_wj/sc-AS/diab_pancr/Lawlor.Genom.Res.2016/count_iret2/")
dir = paste(path,
            file[1:3],
            'SJ.out.tab',
            sep="") 

library(tidyverse)
library(parallel)
sample_meta = data.table::fread("./SJ_phenoData.txt",data.table=F)
samples=list.files("/media/user/sdf/temp_wj/sc-AS/diab_pancr/Lawlor.Genom.Res.2016/count_iret2/")
psi.list = mclapply(samples, function(sample){   
  psi_sp = data.table::fread(paste0("/media/user/sdg/WShi/AS/lawlor/iret2/psi_txt/",sample,".txt"), data.table=F)
  if(nrow(psi_sp)==0) return(NULL)
  psi_sp = psi_sp %>% 
    dplyr::select(V1,V2,V3,V4,V5,V6, V7,PSI)
  colnames(psi_sp)[8] = sample
  return(psi_sp)
},mc.cores=10)
names(psi.list) = samples
psi.list = psi.list[!unlist(lapply(psi.list, is.null))]

psi_res = psi.list[[1]]
for(x in psi.list[-1]){
  psi_res = dplyr::full_join(psi_res, x)
}

for(sample in setdiff(samples,names(psi.list))){
  psi_res[,sample]=NA
}
dim(psi_res)
write.table(psi_res, file="/media/user/sdg/WShi/AS/lawlor/merge/iret_PSI.txt",
            quote=F, row.names=F, sep="\t")


# txt=list.files("/media/user/sdg/WShi/AS/lawlor/taca/psi_txt/")
# txt=str_split_fixed(txt,'[.]',2)[,1]
# setdiff(list.files("/media/user/sdf/temp_wj/sc-AS/diab_pancr/Lawlor.Genom.Res.2016/count_iret2/"),
#         txt)
# 
# samples=samples[1:10]
# samples=list.files("/media/user/sdf/temp_wj/sc-AS/diab_pancr/Lawlor.Genom.Res.2016/count_iret2/")
# psi.list = mclapply(samples, function(sample){   
#   psi_sp = data.table::fread(paste0("/media/user/sdg/WShi/AS/lawlor/alt3/psi_txt/",sample,".txt"), data.table=F)
#   if(nrow(psi_sp)==0) return(NULL)
#   psi_sp=psi_sp[!duplicated(psi_sp$V4),]
#   psi_sp = psi_sp %>% 
#     dplyr::select(V1,V2,V3,V4,V5,V6, V7,PSI)
#   colnames(psi_sp)[8] = sample
#   return(psi_sp)
# },mc.cores=10)
# names(psi.list) = samples
# psi.list = psi.list[!unlist(lapply(psi.list, is.null))]
# 
# psi_res = psi.list[[1]]
# for(x in psi.list[-1]){
#   psi_res = dplyr::full_join(psi_res, x)
# }
# 
# for(sample in setdiff(samples,names(psi.list))){
#   psi_res[,sample]=NA
# }
# dim(psi_res)
# write.table(psi_res, file="/media/user/sdg/WShi/AS/lawlor/merge/alt3_PSI.txt",
#             quote=F, row.names=F, sep="\t")


