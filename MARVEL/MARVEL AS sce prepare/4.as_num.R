
library(data.table)

path <- "/sc-AS/diab_pancr/Lawlor.Genom.Res.2016/trim/trimmed.fq/"
file <- list.files("/sc-AS/diab_pancr/Lawlor.Genom.Res.2016/count_iret2/")
dir = paste(path,
            file[1:3],
            'SJ.out.tab',
            sep="") 

library(tidyverse)
library(parallel)
sample_meta = data.table::fread("./SJ_phenoData.txt",data.table=F)
samples=list.files("/sc-AS/diab_pancr/Lawlor.Genom.Res.2016/count_iret2/")
psi.list = mclapply(samples, function(sample){   
  psi_sp = data.table::fread(paste0("/lawlor/iret2/psi_txt/",sample,".txt"), data.table=F)
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
write.table(psi_res, file="/lawlor/merge/iret_PSI.txt",
            quote=F, row.names=F, sep="\t")


