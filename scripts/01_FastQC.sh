#!/bin/bash
#program: FastQC-v0.11.5

#Directories
input_d="raw"
result_d="qc/fastqc"
log_d="logs/fastqc"

#Options
threads=10

mkdir -p ${result_d} ${log_d}

for sample in $(cat list/list_sample);do
  fastqc ${input_d}/${sample}_1.fastq -t ${threads} -o ${result_d} &> ${log_d}/${sample}_1.log &
  fastqc ${input_d}/${sample}_2.fastq -t ${threads} -o ${result_d} &> ${log_d}/${sample}_2.log
done
