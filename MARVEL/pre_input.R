# Load packages
library(MARVEL)
library(data.table)
library(parallel)

# Read GTF file
# This should be the same file and version as that used in your rMAST step
gtf <- as.data.frame(data.table::fread('/index/gtf/gencode.v38lift37.annotation.gtf', 
                                       sep="\t", header=FALSE, stringsAsFactors=FALSE, quote=""))

samples=list.files("/AS/lawlor/marvel/rMATS/ASanno/")

df.list = mclapply(samples, function(sample){   
  df_sp = read.table(paste0("/AS/lawlor/marvel/rMATS/ASanno/",sample,'/fromGTF.RI.txt'), 
                     sep="\t", header=TRUE, stringsAsFactors=FALSE)
  if(nrow(df_sp)==0) return(NULL)
  df_sp = Preprocess_rMATS(file=df_sp, GTF=gtf, EventType="RI")
  return(df_sp)
},mc.cores=10)
names(df.list) = samples
df.list = df.list[!unlist(lapply(df.list, is.null))]

df_res = df.list[[1]]
for(x in df.list[-1]){
  df_res = dplyr::full_join(df_res, x)
}

write.table(df_res, file="/AS/lawlor/marvel/input_file/RI_res.txt",
            quote=F, row.names=F, sep="\t")

# A5SS
df.list = mclapply(samples, function(sample){   
  df_sp = read.table(paste0("/AS/lawlor/marvel/rMATS/ASanno/",sample,'/fromGTF.A5SS.txt'), 
                     sep="\t", header=TRUE, stringsAsFactors=FALSE)
  if(nrow(df_sp)==0) return(NULL)
  df_sp = Preprocess_rMATS(file=df_sp, GTF=gtf, EventType="A5SS")
  return(df_sp)
},mc.cores=10)
names(df.list) = samples
df.list = df.list[!unlist(lapply(df.list, is.null))]

df_res = df.list[[1]]
for(x in df.list[-1]){
  df_res = dplyr::full_join(df_res, x)
}

write.table(df_res, file="/AS/lawlor/marvel/input_file/A5SS_res.txt",
            quote=F, row.names=F, sep="\t")

# A3SS
df.list = mclapply(samples, function(sample){   
  df_sp = read.table(paste0("/AS/lawlor/marvel/rMATS/ASanno/",sample,'/fromGTF.A3SS.txt'), 
                     sep="\t", header=TRUE, stringsAsFactors=FALSE)
  if(nrow(df_sp)==0) return(NULL)
  df_sp = Preprocess_rMATS(file=df_sp, GTF=gtf, EventType="A3SS")
  return(df_sp)
},mc.cores=10)
names(df.list) = samples
df.list = df.list[!unlist(lapply(df.list, is.null))]

df_res = df.list[[1]]
for(x in df.list[-1]){
  df_res = dplyr::full_join(df_res, x)
}

write.table(df_res, file="/AS/lawlor/marvel/input_file/A3SS_res.txt",
            quote=F, row.names=F, sep="\t")
