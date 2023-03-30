library(tidyverse)
library(parallel)
sample_meta = data.table::fread("./SJ_phenoData.txt",data.table=F)
samples=sample_meta$sample.id

# step1: 汇总剪切位点矩阵

# samples = list.files(paste0("./work/star/"))
# sample="SRR3934349"
sj.list = mclapply(samples, function(sample){   
     sj_sp = data.table::fread(paste0("./work/star/",sample,"/",sample,".SJ.out.tab"), data.table=F)
     if(nrow(sj_sp)==0) return(NULL)
     sj_sp = sj_sp %>% 
       dplyr::mutate(`coord.intron`=paste0(V1,":",V2,":",V3)) %>% 
       dplyr::select(`coord.intron`, V7)
     colnames(sj_sp)[2] = sample
     return(sj_sp)
},mc.cores=10)
names(sj.list) = samples
sj.list = sj.list[!unlist(lapply(sj.list, is.null))]

sj_res = sj.list[[1]]
for(x in sj.list[-1]){
     sj_res = dplyr::full_join(sj_res, x)
}

for(sample in setdiff(samples,names(sj.list))){
     sj_res[,sample]=NA
}
dim(sj_res)

write.table(sj_res, file=paste0("./work/merge/SJ.txt"),
     quote=F, row.names=F, sep="\t")


## step2:整理AS注释结果
source(paste0("./refdata/rscript_Anno_SJ_rMATs.R"))

AS_types = c("SE", "RI", "MXE", "A3SS", "A5SS")

for(AS_type in AS_types){
     # AS_type = AS_types[1]
     print(AS_type)
     fls_AS = list.files(paste0("./work/rmats"), pattern=paste0("fromGTF.",AS_type,".txt"),
          recursive=T, full.name=T)

     if (AS_type=="SE"){
          func = se_func
     } else if (AS_type=="RI"){
          func = ri_func 
     } else if (AS_type=="MXE"){
          func = mxe_func 
     } else if (AS_type=="A3SS"){
          func = a3ss_func 
     } else if (AS_type=="A5SS"){
          func = a5ss_func
     }

     df_AS = mclapply(fls_AS, function(fls){
          rmats_AS = data.table::fread(fls, data.table=F)
          rmats_AS = rmats_AS %>%
               dplyr::mutate(tran_id = func(.)) %>% 
               dplyr::rename(gene_id=GeneID,gene_short_name=geneSymbol) %>% 
               dplyr::select(tran_id, gene_id, gene_short_name)
          return(rmats_AS)
     },mc.cores=10) %>% do.call(rbind, .) %>% 
          dplyr::distinct()
     write.table(df_AS, row.names=F, sep="\t",quote=F,
          file=paste0("./work/merge/as_",AS_type,"_featureData.txt"))
}





## step3: 汇总intron矩阵
# samples = list.files(paste0("./work/intron/"))[100]
int_mt.list = mclapply(samples, function(sample){
     # sample = samples[1]
     int_mt = data.table::fread(paste0("./work/intron/",sample,"/intron_count.txt")) 
     int_mt = int_mt %>% 
          dplyr::mutate(`coord.intron`=paste(V1,V2,V3, sep=":")) %>% 
          dplyr::group_by(coord.intron) %>% 
          dplyr::summarize(count=sum(V5)) %>% 
          dplyr::select(coord.intron, count) %>% as.data.frame()
     colnames(int_mt)[2] = sample
     return(int_mt)
},mc.cores=10)

intron_res = int_mt.list[[1]]
for(x in int_mt.list[-1]){
     intron_res = dplyr::full_join(intron_res, x)
}
dim(intron_res)

write.table(intron_res, file=paste0("./work/merge/intron_count_by_region.txt"),
     quote=F, row.names=F, sep="\t")




## step4: 汇总gene矩阵
# samples = list.files(paste0("./work/rsem/"))
gene_mt.list = mclapply(samples, function(sample){
     fls=paste0("./work/rsem/",sample,"/",sample,".genes.results")
     gene_mt = data.table::fread(fls, data.table=F) %>% 
          dplyr::select(gene_id, TPM)
     colnames(gene_mt)[2] = sample
     return(gene_mt)
},mc.cores=10)
gene_res = gene_mt.list[[1]]
for(x in gene_mt.list[-1]){
     gene_res = dplyr::full_join(gene_res, x)
}
dim(gene_res)
write.table(gene_res, file=paste0("./work/merge/rsem_tpm.txt"),,
     quote=F, row.names=F, sep="\t")



## step5: 初步创建MARVEL对象
library(tidyverse)
library(MARVEL)

# sample_meta = data.table::fread("./work/SraRunTable.txt",data.table=F)
# samples=sample_meta$Run[1:500]

# (1) 细胞表型分组
df.pheno = data.table::fread("./SJ_phenoData.txt",data.table=F)

# (2) Gene文件
df.tpm <- read.table("./work/merge/rsem_tpm.txt", sep="\t", header=TRUE)
dim(df.tpm)


# (3) GTF文件
gtf <- data.table::fread("../basic/gtf/gencode.v31.annotation.gtf", sep="\t", header=FALSE, na.strings="NA",quote="\"") %>% 
     as.data.frame()

# (4) Gene类型
gene.feature = subset(gtf, V3=="gene") 
df.tpm.feature <- str_split(gene.feature$V9, ";", simplify=T) %>% 
     as.data.frame() %>%
     dplyr::mutate(gene_id=str_match(V1, '"(.*)"')[,2]) %>% 
     dplyr::mutate(gene_short_name=str_match(V3, '"(.*)"')[,2]) %>%      
     dplyr::mutate(gene_type=str_match(V2, '"(.*)"')[,2]) %>% 
     dplyr::select(gene_id, gene_short_name, gene_type) %>% 
     dplyr::filter(gene_id %in% df.tpm$gene_id)
df.tpm.feature = df.tpm.feature[match(df.tpm$gene_id, df.tpm.feature$gene_id),]
rownames(df.tpm.feature) = seq(nrow(df.tpm.feature))


# (5) SJ文件
sj = data.table::fread("./work/merge/SJ.txt", sep="\t", header=TRUE, na.strings="NA") %>% 
     as.data.frame()
dim(sj)



# (6) AS文件
df.feature.list = lapply(c("SE", "MXE", "RI", "A5SS", "A3SS"), function(x){
     # x = "SE"
     df.feature = read.table(paste0("./work/merge/as_",x,"_featureData.txt"), 
          sep="\t", header=TRUE, na.strings="NA") %>% 
          dplyr::distinct(tran_id, .keep_all=TRUE) %>% 
          dplyr::left_join(df.tpm.feature)
     return(df.feature)
})
names(df.feature.list) <- c("SE", "MXE", "RI", "A5SS", "A3SS")



# (7) Intron文件
df.intron.counts <- data.table::fread("./work/merge/intron_count_by_region.txt", 
     sep="\t", header=TRUE, na.strings="NA") %>% 
     as.data.frame()
dim(df.intron.counts)




## MARVEL object

marvel <- CreateMarvelObject(SpliceJunction=sj,
                             SplicePheno=df.pheno,
                             SpliceFeature=df.feature.list,
                             IntronCounts=df.intron.counts,
                             GeneFeature=df.tpm.feature,
                             Exp=df.tpm,
                             GTF=gtf
                             )
# 预处理
marvel <- TransformExpValues(MarvelObject=marvel,
                             offset=1,
                             transformation="log2",
                             threshold.lower=1
                            )

marvel <- DetectEvents(MarvelObject=marvel,
                       min.cells=20, #细胞群的表达百分比
                       min.expr=1,   #认为细胞表达该基因的阈值
                       track.progress=TRUE,
                       EventType="AFE"
                       )
marvel <- DetectEvents(MarvelObject=marvel,
                       min.cells=20,
                       min.expr=1,
                       track.progress=FALSE,
                       EventType="ALE"
                       )
length(marvel$SpliceFeature)

# 计算PSI
marvel <- ComputePSI(MarvelObject=marvel,
                     CoverageThreshold=10,
                     UnevenCoverageMultiplier=10,
                     EventType="SE"
                     )

marvel <- ComputePSI(MarvelObject=marvel,
                     CoverageThreshold=10,
                     UnevenCoverageMultiplier=10,
                     EventType="MXE"
                     )

marvel <- ComputePSI(MarvelObject=marvel,
                     CoverageThreshold=5,
                     EventType="RI",
                     thread=4  # only support RI
                     )

marvel <- ComputePSI(MarvelObject=marvel,
                     CoverageThreshold=10,
                     EventType="A5SS"
                     )

marvel <- ComputePSI(MarvelObject=marvel,
                     CoverageThreshold=10,
                     EventType="A3SS"
                     )
 
marvel <- ComputePSI(MarvelObject=marvel,
                     CoverageThreshold=10,
                     EventType="AFE"
                     )
    
marvel <- ComputePSI(MarvelObject=marvel,
                     CoverageThreshold=10,
                     EventType="ALE"
                     )



marvel <- CheckAlignment(MarvelObject=marvel, level="SJ")
marvel <- CheckAlignment(MarvelObject=marvel, level="splicing")
marvel <- CheckAlignment(MarvelObject=marvel, level="gene")
marvel <- CheckAlignment(MarvelObject=marvel, level="splicing and gene")

save(marvel, file="./work/MARVEL.RData")






# ## intron coord
# bed = data.table::fread("/home/lishensuo/STUDY/marvel/practice/work/intron/SRR3934349/RI_Coordinates_sorted.bed")
# table(bed$V2 < bed$V3)

# bed_ref = data.table::fread("/home/lishensuo/STUDY/marvel/plate_test/Data/MARVEL/PSI/RI/Counts_by_Region.txt")
# bed_ref = str_split(bed_ref$coord.intron, ":", simplify=T) %>% as.data.frame()
# table(bed_ref$V2 < bed_ref$V3)



# ri_ref = data.table::fread("/home/lishensuo/STUDY/marvel/plate_test/Data/rMATS/RI/RI_featureData.txt")

# ri_ref = str_split(gsub("+@","",ri_ref$tran_id), ":", simplify=TRUE) %>% 
#      as.data.frame()
# summary((as.numeric(ri_ref$V5)-as.numeric(ri_ref$V3)))
