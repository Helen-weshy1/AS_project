
path <- "/sc-AS/diab_pancr/Lawlor.Genom.Res.2016/count_iret/"
file <- list.files("/sc-AS/diab_pancr/Lawlor.Genom.Res.2016/count_iret/")
# sj <- as.data.frame(fread(paste(path, file, sep=""), sep="\t", header=TRUE, stringsAsFactors=FALSE, na.strings="NA"))
# 
# sj[!is.na(sj[,2]), ][1:5,1:5]

dir = paste("/sc-AS/diab_pancr/Lawlor.Genom.Res.2016/count_iret/",
            file,
            '/iret.count.txt',
            sep="")                 
n = length(dir)                                       
n
head(dir)


for (i in 1:1157){
  tryCatch({new.data = read.table(file =  dir[i],header=F, sep="\t")
  new.data$PSI=new.data$V9/(new.data$V9+new.data$V10)
  write.table(new.data,file = paste0("/media/user/sdg/WShi/AS/lawlor/iret/psi_txt/",
                                     file[i],'.txt'),row.names=FALSE,quote = F)},
  error=function(e){print('error')})
  
}

dir = paste("/sc-AS/diab_pancr/Lawlor.Genom.Res.2016/count_alt3/",
            file,
            '/alt3.count.txt',
            sep="")                 
n = length(dir)                                       
n
head(dir)


for (i in 1:1157){
  tryCatch({new.data = read.table(file =  dir[i],header=F, sep="\t")
  new.data$PSI=new.data$V9/(new.data$V9+new.data$V10)
  write.table(new.data,file = paste0("/media/user/sdg/WShi/AS/lawlor/alt3/psi_txt/",
                                     file[i],'.txt'),row.names=FALSE,quote = F)},
  error=function(e){print('error')})
  
}


dir = paste("/sc-AS/diab_pancr/Lawlor.Genom.Res.2016/count_alt5/",
            file,
            '/alt5.count.txt',
            sep="")                 
n = length(dir)                                       
n
head(dir)


for (i in 1:1157){
  tryCatch({new.data = read.table(file =  dir[i],header=F, sep="\t")
  new.data$PSI=new.data$V9/(new.data$V9+new.data$V10)
  write.table(new.data,file = paste0("/media/user/sdg/WShi/AS/lawlor/alt5/psi_txt/",
                                     file[i],'.txt'),row.names=FALSE,quote = F)},
  error=function(e){print('error')})
  
}



### iret
path <- "/sc-AS/diab_pancr/Lawlor.Genom.Res.2016/count_iret/"
file <- list.files("/sc-AS/diab_pancr/Lawlor.Genom.Res.2016/count_iret/")
# sj <- as.data.frame(fread(paste(path, file, sep=""), sep="\t", header=TRUE, stringsAsFactors=FALSE, na.strings="NA"))
# 
# sj[!is.na(sj[,2]), ][1:5,1:5]

dir = paste("/sc-AS/diab_pancr/Lawlor.Genom.Res.2016/count_iret/",
            file,
            '/iret.count.txt',
            sep="")                 
n = length(dir)                                       
n
head(dir)

for (i in 1:1157){
  tryCatch({new.data = read.table(file =  dir[i],header=F, sep="\t")
  new.data$PSI=new.data$V9/(new.data$V9+new.data$V10)
  write.table(new.data,file = paste0("/media/user/sdg/WShi/AS/lawlor/iret/psi_txt/",
                                     file[i],'.txt'),row.names=FALSE,quote = F)},
  error=function(e){print('error')})
  
}

dir = paste("/sc-AS/diab_pancr/Lawlor.Genom.Res.2016/count_taca/",
            file,
            '/taca.count.txt',
            sep="")                 
n = length(dir)                                       
n
head(dir)
for (i in 1:1157){
  tryCatch({new.data = read.table(file =  dir[i],header=F, sep="\t")
  new.data$PSI=new.data$V9/(new.data$V9+new.data$V10)
  write.table(new.data,file = paste0("/media/user/sdg/WShi/AS/lawlor/taca/psi_txt/",
                                     file[i],'.txt'),row.names=FALSE,quote = F)},
  error=function(e){print('error')})
  
}
