##output (if any reads in region)
##
## SAMPLE_unknown_R1.fastq
## SAMPLE_unknown_R2.fastq
## SAMPLE_V1V2_R1.fastq
## SAMPLE_V1V2_R2.fastq
## SAMPLE_V2V3_R1.fastq
## SAMPLE_V2V3_R2.fastq
## SAMPLE_V3V4_R1.fastq
## SAMPLE_V3V4_R2.fastq
## SAMPLE_V4V5_R1.fastq
## SAMPLE_V4V5_R2.fastq
## SAMPLE_V5V7_R1.fastq
## SAMPLE_V5V7_R2.fastq
## SAMPLE_V5V9_R1.fastq
## SAMPLE_V7V9_R2.fastq
## SAMPLE_ITS_R1.fastq
## SAMPLE_ITS_R2.fastq


cutadapt -g file:qiaseq_primers_fwd.fa -G file:qiaseq_primers_rev.fa --pair-adapters --no-indels -e 0.1 --suffix ':region={name}' -o [SAMPLE]_{name}_R1.fastq -p [SAMPLE]_{name}_R2.fastq data/[SAMPLE]_R1.fastq data/[SAMPLE]_R2.fastq > [SAMPLE]_qiaseq_demultiplex.log


# python qiaseq_region_summary.py *_qiaseq_demultiplex.log > qiaseq_regions_mqc.yaml
