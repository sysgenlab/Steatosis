#!/bin/bash
#program: MACS2-2.1.1, UCSCtoolkit, bedtools-2.27.1

blacklist=/home/refseq/Blacklist/hg38-blacklist.v2.bed.gz

#Directories
input_d=browser
result_d=macs2
log_d=logs/macs2

#Options
##https://deeptools.readthedocs.io/en/develop/content/feature/effectiveGenomeSize.html
eff_gsize=hs
q_thresh=0.05

mkdir -p ${result_d} ${log_d}

#Peak Calling
##macs2-2.1.1
##shiftsize=smooth_window/2, smooth_window=median_insert_size from InsertSizeMetircs(picard)
for sample in $(cat list/list_sample);do
  macs2 callpeak -t ${input_d}/${sample}.filt.nmsrt.nodup.tn5.bed.gz -f BED \
    -n ${result_d}/${sample}.tn5 -g ${eff_gsize} -q ${q_thresh} --nomodel --shift -100 --extsize 200 \
    -B --SPMR --keep-dup all --tempdir tmp 2> ${log_d}/macs2_callpeak_${sample}.log

  sort -k 8gr,8gr ${result_d}/${sample}.tn5_peaks.narrowPeak | awk 'BEGIN{OFS="\t"}{$4="Peak_"NR;print $0}' | gzip -nc > ${result_d}/${sample}.tn5.narrowPeak.gz
  #rm -f ${result_d}/${sample}_peaks.narrowPeak
  #rm -f ${result_d}/${sample}_peaks.xls
  #rm -f ${result_d}/${sample}_summits.bed

  #Remove blacklist regions
  bedtools intersect -v -a ${result_d}/${sample}.tn5.narrowPeak.gz -b ${blacklist} \
    | grep -P 'chr[\dXY]+[\t]' | gzip -nc > ${result_d}/${sample}.tn5.narrowPeak.filt.gz

  #Create browser bed file
  zcat ${result_d}/${sample}.tn5.narrowPeak.filt.gz | sortBed -i stdin > ${result_d}/${sample}.tn5_peaks.bed
done
