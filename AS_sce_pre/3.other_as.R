setwd('/media/user/sdd/zhaolabguest/wsy/sc_as/lawlor.1050/other_AS/')

path <- "/media/user/sdf/temp_wj/sc-AS/diab_pancr/Lawlor.Genom.Res.2016/count_iret/"
file <- list.files("/media/user/sdf/temp_wj/sc-AS/diab_pancr/Lawlor.Genom.Res.2016/count_iret/")
# sj <- as.data.frame(fread(paste(path, file, sep=""), sep="\t", header=TRUE, stringsAsFactors=FALSE, na.strings="NA"))
# 
# sj[!is.na(sj[,2]), ][1:5,1:5]

#命令构建路径变量dir（方便更改），也可以不构建，后面示例                                        
dir = paste("/media/user/sdf/temp_wj/sc-AS/diab_pancr/Lawlor.Genom.Res.2016/count_iret/",
            file,
            '/iret.count.txt',
            sep="")                 
#读取dir长度，也就是文件夹下的文件个数
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

dir = paste("/media/user/sdf/temp_wj/sc-AS/diab_pancr/Lawlor.Genom.Res.2016/count_alt3/",
            file,
            '/alt3.count.txt',
            sep="")                 
#读取dir长度，也就是文件夹下的文件个数
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


dir = paste("/media/user/sdf/temp_wj/sc-AS/diab_pancr/Lawlor.Genom.Res.2016/count_alt5/",
            file,
            '/alt5.count.txt',
            sep="")                 
#读取dir长度，也就是文件夹下的文件个数
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



### iret2
path <- "/media/user/sdf/temp_wj/sc-AS/diab_pancr/Lawlor.Genom.Res.2016/count_iret2/"
file <- list.files("/media/user/sdf/temp_wj/sc-AS/diab_pancr/Lawlor.Genom.Res.2016/count_iret2/")
# sj <- as.data.frame(fread(paste(path, file, sep=""), sep="\t", header=TRUE, stringsAsFactors=FALSE, na.strings="NA"))
# 
# sj[!is.na(sj[,2]), ][1:5,1:5]

#命令构建路径变量dir（方便更改），也可以不构建，后面示例                                        
dir = paste("/media/user/sdf/temp_wj/sc-AS/diab_pancr/Lawlor.Genom.Res.2016/count_iret2/",
            file,
            '/iret.count.txt',
            sep="")                 
#读取dir长度，也就是文件夹下的文件个数
n = length(dir)                                       
n
head(dir)

for (i in 1:1157){
  tryCatch({new.data = read.table(file =  dir[i],header=F, sep="\t")
  new.data$PSI=new.data$V9/(new.data$V9+new.data$V10)
  write.table(new.data,file = paste0("/media/user/sdg/WShi/AS/lawlor/iret2/psi_txt/",
                                     file[i],'.txt'),row.names=FALSE,quote = F)},
  error=function(e){print('error')})
  
}

dir = paste("/media/user/sdf/temp_wj/sc-AS/diab_pancr/Lawlor.Genom.Res.2016/count_taca/",
            file,
            '/taca.count.txt',
            sep="")                 
#读取dir长度，也就是文件夹下的文件个数
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
