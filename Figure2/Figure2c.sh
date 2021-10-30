#!/bin/bash
#program: picard-tools-2.18.1

#Directories
input_d=../bowtie2/Filtering
metrics_d=picard
picard=/data/program/picard-2.18.1/picard.jar

mkdir -p $metrics_d

#Measure metrics for validating library construction
##picard-2.18.1
for sample in $(cat list/list_sample);do
  java -Xmx16g -Djava.io.tmpdir=./tmp -jar $picard \
      CollectInsertSizeMetrics INPUT=${input_d}/${sample}.filt.srt.nodup.bam \
      OUTPUT=${metrics_d}/${sample}.inserts.hist_data.log \
      H=${metrics_d}/${sample}.inserts.hist_graph.pdf VERBOSITY=ERROR QUIET=TRUE \
      W=1000 STOP_AFTER=5000000
done
