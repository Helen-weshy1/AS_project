# https://github.com/lishensuo/utils/blob/main/marvel/rscript_bedtools_input.R
args = commandArgs(T)
k = args[1]

introns = data.table::fread(paste0("/marvel/rMATS/ASanno/",k,"/fromGTF.RI.txt"))
intron_coord = introns[,c(4, 9, 10)]
intron_coord$chr = factor(intron_coord$chr, levels=c(paste0("chr",c(11:19,1,20:22,2:9)), "chrM","chrX","chrY"), order=T)
intron_coord = intron_coord[order(intron_coord$chr, 
                                  intron_coord$upstreamEE),]
intron_coord$chr=as.character(intron_coord$chr)
intron_coord = na.omit(intron_coord)

write.table(intron_coord, file=paste0("/marvel/intron/", k, "/RI_Coordinates_sorted.bed"), row.names=F, col.names=F,
            quote=F, sep="\t")


chr_size = data.table::fread(paste0("/STAR_index/hg19/chrNameLength.txt"),header=F,data.table=F)
chr_order = data.table::fread(paste0("/marvel/bam_sort/", k, "/sorted_chr_in_bam.txt"),header=F,data.table=F)
chr_out = dplyr::inner_join(chr_order, chr_size)

write.table(chr_out, file=paste0("/marvel/intron/", k, "/hg19.chrom.sizes.txt"), row.names=F, col.names=F,
            quote=F, sep="\t")