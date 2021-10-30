#!/bin/bash
#program: samtools-1.7, picard-tools-2.18.1

#Directories
input_d=bowtie2
result_d=${input_d}/Filtering
log_d=logs/MarkDuplicates
qc=${input_d}/qc

#Options
threads=20
picard=/data/program/picard-2.18.1/picard.jar

mkdir ${result_d} ${log_d} ${qc}

for sample in $(cat list/list_sample);do
  #Remove unmapped, unpaired mates, reads failing platform
  ##nmsrt: read name sorted; srt: position sorted
  samtools view -@ ${threads} -F 524 -u ${input_d}/${sample}.bam | samtools sort -@ ${threads} -n -o tmp.${sample}.filt.srt.nmsrt.multimap.bam -
  samtools view -@ ${threads} -h tmp.${sample}.filt.srt.nmsrt.multimap.bam | \
    python ./scripts/assign_multimappers.py -k 4 --paired-end | \
    samtools view -@ ${threads} -bS - > tmp.${sample}.filt.srt.nmsrt.bam
  rm -f tmp.${sample}.filt.srt.nmsrt.multimap.bam

  #Remove orphan reads and read pairs mapping to different chromosomes
  samtools fixmate -@ ${threads} -r tmp.${sample}.filt.srt.nmsrt.bam tmp.${sample}.fixmate
  samtools view -@ ${threads} -F 1804 -f 2 -u tmp.${sample}.fixmate | samtools sort -@ ${threads} -o ${result_d}/${sample}.filt.srt.bam -
  samtools flagstat -@ ${threads} ${result_d}/${sample}.filt.srt.bam > ${qc}/${sample}.filt.srt.flagstat.qc
  rm tmp.${sample}.filt.srt.nmsrt.bam
  rm tmp.${sample}.fixmate

  #Mark duplicate by Picard
  ##picard-tools-2.18.1
  java -Xmx16g -Djava.io.tmpdir=./tmp -jar $picard MarkDuplicates \
      INPUT=${result_d}/${sample}.filt.srt.bam OUTPUT=tmp.${sample}.filt.srt.dupmark.bam \
      METRICS_FILE=${qc}/${sample}.filt.srt.dup.qc VALIDATION_STRINGENCY=LENIENT \
      ASSUME_SORTED=true REMOVE_DUPLICATES=false \
      &> ${log_d}/${sample}.MarkDup.log
  mv tmp.${sample}.filt.srt.dupmark.bam ${result_d}/${sample}.filt.srt.bam

  #Remove duplicates
  samtools view -@ ${threads} -F 1804 -f 2 -b ${result_d}/${sample}.filt.srt.bam > ${result_d}/${sample}.filt.srt.nodup.bam
  samtools index -@ ${threads} ${result_d}/${sample}.filt.srt.nodup.bam ${result_d}/${sample}.filt.srt.nodup.bai
  samtools flagstat -@ ${threads} ${result_d}/${sample}.filt.srt.nodup.bam > ${qc}/${sample}.filt.srt.nodup.flagstat.qc

  #Sort reads by name > used for creating BEDPE and tagAlign file
  samtools sort -@ ${threads} -n -o ${result_d}/${sample}.filt.nmsrt.nodup.bam ${result_d}/${sample}.filt.srt.nodup.bam

  #Split BAM into chrM and non-chrM
  nonMitoChromosomes=$(samtools view -H ${result_d}/${sample}.filt.srt.nodup.bam | grep chr | cut -f2 | sed 's/SN://g' | grep -v chrM)
  samtools view -b ${result_d}/${sample}.filt.srt.nodup.bam ${nonMitoChromosomes} > ${result_d}/${sample}.filt.srt.nodup.nonchrM.bam
  samtools flagstat -@ ${threads} ${result_d}/${sample}.filt.srt.nodup.nonchrM.bam > ${qc}/${sample}.filt.srt.nodup.nonchrM.flagstat.qc
  samtools view -b ${result_d}/${sample}.filt.srt.nodup.bam chrM > ${result_d}/${sample}.filt.srt.nodup.chrM.bam
done
