#!/bin/bash
#program: bowtie2-2.3.4.1

#Directories
input_d=raw
result_d=bowtie2
log_d=logs/bowtie2

#Options
index=/home/refseq/bowtie2_indexes/hg38
multimapping=4
threads=20

mkdir ${result_d}
mkdir ${log_d}

#bowtie2-align-s v2.3.4.1
for sample in $(cat list/list_sample);do
  fq1=${input_d}/${sample}_1.fastq.gz
  fq2=${input_d}/${sample}_2.fastq.gz
  bowtie2 -k $multimapping -X2000 --end-to-end --mm -x $index -t -p $threads -q -1 $fq1 -2 $fq2 \
    2> ${log_d}/${sample}.bowtie2.log | samtools view -@ 2 -bS - > ${result_d}/${sample}.bam
done
