bam=/picb/qilab/stu.04.huzhuocheng/NGS/PON/DLBCL/DLBCL_RA202406120161/03MarkBQSR/BQSR_DLBCL_RA202406120161.bam
dir=/picb/qilab/stu.04.huzhuocheng/code/NGS/CNV_CNVkit/codeByStep0516/TestMQ
bed=/picb/qilab/stu.04.huzhuocheng/NGS/CNVkit/0625Test/Bed/BQSRmergedTarget.bed

# 过滤出 MQ 为 0 的 reads 并输出为新的 BAM 文件
samtools view -b -q 0 -e 'mapq == 0' $bam > $dir/mq0.bam

# 将过滤后的 BAM 文件转换为 BED 文件
bedtools bamtobed -i $dir/mq0.bam > $dir/mq0.bed

# 计算 MQ0 的每个位置的深度
bedtools coverage -d -a $bed -b $dir/mq0.bed > $dir/mq0_coverage.txt

# 计算 MQ20 的每个位置的深度
bedtools coverage -d -a $bed -b $bam > $dir/mq_coverage.txt

