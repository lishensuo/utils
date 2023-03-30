args = commandArgs(T)
id = args[1]
# Dir_pare = args[2]
# id = "ERR1562273"
# Dir_pare = "/home/lishensuo/STUDY/marvel/practice2"

chr_size = data.table::fread(paste0("./refdata/hg38.chrom.sizes"),header=F,data.table=F)

chr_order = data.table::fread(paste0("./work/star/", id, "/sorted_chr_in_bam.txt"),header=F,data.table=F)
# table(chr_order[,1] %in% chr_size[,1])
# chr_out = dplyr::inner_join(chr_order, chr_size)
# chr_out = chr_out[grep("_", chr_out$V1,fixed=T, invert=T),]

introns = data.table::fread(paste0("./work/rmats/",id,"/fromGTF.RI.txt"))
intron_coord = introns[,c(4, 9, 10)]
intron_coord$chr = factor(intron_coord$chr, levels=c(paste0("chr",1:22),"chrX","chrY", "chrM"), order=T)
intron_coord = intron_coord[order(intron_coord$chr, intron_coord$upstreamEE),]
intron_coord$chr=as.character(intron_coord$chr)
intron_coord = na.omit(intron_coord)

chr_size = data.table::fread(paste0("./refdata/hg38.chrom.sizes"),header=F,data.table=F)
chr_order = data.table::fread(paste0("./work/star/", id, "/sorted_chr_in_bam.txt"),header=F,data.table=F)
# chr_out = chr_size[match(unique(intron_coord$chr), chr_size$V1),]
chr_out = dplyr::inner_join(chr_order, chr_size)
chr_out = chr_out[grep("_", chr_out$V1,fixed=T, invert=T),]

write.table(intron_coord, file=paste0("./work/intron/", id, "/RI_Coordinates_sorted.bed"), row.names=F, col.names=F,
     quote=F, sep="\t")

write.table(chr_out, file=paste0("./work/intron/", id, "/hg38.chrom.sizes.txt"), row.names=F, col.names=F,
     quote=F, sep="\t")
