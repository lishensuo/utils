# conda activate marvel_plate
# mkdir -p work/fq work/star work/rmats work/intron work/merge  work/rsem 


id=ERR1562083
# id=$1
## step-1 下载数据
mkdir -p ./work/fq/${id}
ascp -QT -l 300m -P33001  \
-i ~/miniconda3/envs/marvel_plate/etc/asperaweb_id_dsa.openssh   \
era-fasp@fasp.sra.ebi.ac.uk:/vol1/fastq/${id:0:6}/00${id:0-1}/${id}/${id}_1.fastq.gz  ./work/fq/${id}

ascp -QT -l 300m -P33001  \
-i ~/miniconda3/envs/marvel_plate/etc/asperaweb_id_dsa.openssh   \
era-fasp@fasp.sra.ebi.ac.uk:/vol1/fastq/${id:0:6}/00${id:0-1}/${id}/${id}_2.fastq.gz  ./work/fq/${id}


## step-2 star比对
mkdir -p ./work/star/${id}
ref_idx_star=../basic/gtf/star_index
ref_gtf=../basic/gtf/gencode.v31.annotation.gtf
echo ${id} is runing......
#1st pass mode
STAR --runThreadN 16 \
--genomeDir ${ref_idx_star} \
--readFilesCommand zcat \
--readFilesIn ./work/fq/${id}/${id}_1.fastq.gz ./work/fq/${id}/${id}_2.fastq.gz \
--outFileNamePrefix ./work/star/${id}/${id}. \
--outSAMtype None
#2nd pass mode
STAR --runThreadN 16 \
--genomeDir ${ref_idx_star} \
--readFilesCommand zcat \
--readFilesIn ./work/fq/${id}/${id}_1.fastq.gz ./work/fq/${id}/${id}_2.fastq.gz \
--outFileNamePrefix ./work/star/${id}/${id}. \
--sjdbFileChrStartEnd ./work/star/${id}/*SJ.out.tab \
--sjdbGTFfile ${ref_gtf} \
--outSAMtype BAM SortedByCoordinate \
--outSAMattributes NH HI AS nM XS \
--quantMode TranscriptomeSAM


## step-3 rmats注释
mkdir -p ./work/rmats/${id}

echo ${id} is running......
echo "./work/star/${id}/${id}.Aligned.sortedByCoord.out.bam" > ./work/star/${id}/BAM_fls.txt

rmats.py \
--b1 ./work/star/${id}//BAM_fls.txt \
--gtf ${ref_gtf} \
--od ./work/rmats/${id} \
--tmp ./work/rmats/${id}/tmp \
-t paired \
--readLength 125 \
--variable-read-length \
--nthread 16 \
--statoff
# readLength参考如下
# zcat ./work/fq/${id}_1.fastq.gz | awk '{if(NR%4==2) print NR"\t"$0"\t"length($0)}' | less


## step-4 内含子定量
mkdir -p ./work/intron/${id}
echo ${id} is running......
samtools view -H ./work/star/${id}/${id}.Aligned.sortedByCoord.out.bam | \
  grep SQ | cut -f 2 | awk '{ sub(/^SN:/, ""); print;}' > ./work/star/${id}/sorted_chr_in_bam.txt
# 预定义R脚本，生成bedtools计算所需的两个文件：染色体大小以及内含子坐标
Rscript ../basic/rscript_bedtools_input.R ${id}

bedtools coverage \
-g ./work/intron/${id}/hg38.chrom.sizes.txt \
-split \
-sorted \
-a ./work/intron/${id}/RI_Coordinates_sorted.bed \
-b ./work/star/${id}/${id}.Aligned.sortedByCoord.out.bam > \
./work/intron/${id}/intron_count.txt \
-d


## step-5 基因表达定量
mkdir -p ./work/rsem/${id}
rsem-calculate-expression \
--bam \
--paired-end \
-p 16 \
./work/star/${id}/${id}.Aligned.toTranscriptome.out.bam \
../basic/gtf/rsem_index/hg38 \
./work/rsem/${id}/${id}


## 删除耗费内存文件
rm ./work/fq/${id}/${id}*
rm  ./work/star/${id}/*bam
rm -rf ./work/rmats/${id}/tmp/*
rm ./work/rsem/${id}/${id}.transcript.bam
rm ./work/rsem/${id}/${id}.isoforms.results



# chmod u+x ../basic/single_sample_upstream.sh
# ../basic/single_sample_upstream.sh ERR1562273

# cat ./work/SraRunTable.txt | grep -v "Run," | cut -d , -f 1 | head 100



# nohup parallel -j 5 ../basic/single_sample_upstream.sh {} ::: \
#   $(grep -v "sample" ./SJ_phenoData.txt | awk  '{print $1}') \
#   1> ./work/single_sample_upstream.log 2>&1 &
# # 138323


