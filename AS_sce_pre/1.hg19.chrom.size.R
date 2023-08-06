chr_size = data.table::fread(paste0("/media/user/sdg/WShi/STAR/STAR_index/hg19/chrNameLength.txt"),header=F,data.table=F)
chr_order = data.table::fread(paste0("/media/user/sdg/WShi/AS/lawlor/marvel/bam_sort/", k, "/sorted_chr_in_bam.txt"),header=F,data.table=F)
#chr_out = chr_size[match(unique(intron_coord$chr), chr_size$V1),]
chr_out = dplyr::inner_join(chr_order, chr_size)
#chr_out = chr_out[grep("_", chr_out$V1,fixed=T, invert=T),]
#chr_out$V1 = factor(chr_out$V1, levels=c(paste0("chr",1:22),"chrX","chrY", "chrM"), order=T)
chr_out=chr_out[order(chr_out$V1),]

write.table(chr_out, file=paste0("/media/user/sdg/WShi/AS/lawlor/marvel/intron/hg19.chrom.sizes.txt"), row.names=F, col.names=F,
            quote=F, sep="\t")