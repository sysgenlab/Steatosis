#!/bin/bash
#program: bedtools-2.26.0, deepTools-2.3.5, picard-tools-2.18.1

#Directories
input_d=bowtie2/Filtering
qc_d=bowtie2/qc

summary_d=deeptools/summary
heatmap_d=deeptools/heatmap
metrics_d=picard/CollectInsertSizeMetrics

#Options
threads=20
picard=/data/program/picard-2.18.1/picard.jar

mkdir -p ${summary_d} ${heatmap_d} ${metrics_d}

#Compute library complexity
##bedtools-2.26.0
###TotalReadPairs [tab] DistinctReadPairs [tab] OneReadPair [tab] TwoReadPairs [tab] 
###NRF=Distinct/Total [tab] PBC1=OnePair/Distinct [tab] PBC2=OnePair/TwoPair
for sample in $(cat list/list_sample);do
  samtools sort -@ ${threads} -n -o tmp.${sample}.filt.srt.bam ${input_d}/${sample}.filt.srt.bam
  bedtools bamtobed -bedpe -i tmp.${sample}.filt.srt.bam | \
      awk 'BEGIN{OFS="\t"}{print $1,$2,$4,$6,$9,$10}' | grep -v 'chrM' | \
      sort | uniq -c | awk 'BEGIN{mt=0;m0=0;m1=0;m2=0} ($1==1){m1=m1+1} \
      ($1==2){m2=m2+1} {m0=m0+1} {mt=mt+$1} END{printf \
      "%d\t%d\t%d\t%d\t%f\t%f\t%f\n",mt,m0,m1,m2,m0/mt,m1/m0,m1/m2}' \
      > ${qc_d}/${sample}.filt.srt.nodup.pbc.qc
  rm tmp.${sample}.filt.srt.bam
done

#Computes the read coverages for genomic regions
##deepTools-2.3.5
for sample in $(cat list/list_sample);do
  samtools index -@ ${threads} ${input_d}/${sample}.filt.srt.nodup.nonchrM.bam
done

BAM=$(ls -d ${input_d}/* | grep nonchrM | grep -v bai)
multiBamSummary bins --bamfiles ${BAM} -out ${summary_d}/pilot.npz \
    --outRawCounts ${summary_d}/pilot.tab -p ${threads} --ignoreDuplicates

#Generate a heatmap
plotCorrelation -in ${summary_d}/pilot.npz --corMethod spearman --skipZeros \
    --removeOutliers --plotTitle "Spearman Correlation" --whatToPlot heatmap \
    --colorMap coolwarm --plotNumbers -o ${heatmap_d}/heatmap.SpearmanCorr.pilot.svg \
    --outFileCorMatrix ${heatmap_d}/SpearmanCorr.pilot.tab

#Measure metrics for validating library construction
##picard-2.18.1
for sample in $(cat list/list_sample);do
  java -Xmx16g -Djava.io.tmpdir=./tmp -jar $picard \
      CollectInsertSizeMetrics INPUT=${input_d}/${sample}.filt.srt.nodup.bam \
      OUTPUT=${metrics_d}/${sample}.inserts.hist_data.log \
      H=${metrics_d}/${sample}.inserts.hist_graph.pdf VERBOSITY=ERROR QUIET=TRUE \
      W=1000 STOP_AFTER=5000000
done
