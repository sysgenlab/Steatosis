#!/bin/bash
#program: bedtools-2.27.1, UCSCtoolkit

genome=/home/refseq/fetch/human.GRCh38.chr.genome

#Directories
input_d=bowtie2/Filtering
result_d=browser
bedgraph_d=${result_d}/bedgraph
bigwig_d=${result_d}/bigwig

mkdir ${result_d} ${bedgraph_d} ${bigwig_d}
cp ${genome} ${result_d}

#Create bedgraph and bigwig files
##bedtools-2.27.1
##scale = 1,000,000 / # of De-duplicated pairs aligned
##genome files are ignored when the input file is bam format
for sample in $(cat list/list_sample);do
  DEDUP=`head -1 bowtie2/qc/${sample}.filt.srt.nodup.flagstat.qc | sed 's/ .*/\/2/g' | bc`
  scale=`echo "1000000/$DEDUP" | bc -l | sed 's/^/0/g'`

  genomeCoverageBed -bg -ibam ${input_d}/${sample}.filt.srt.nodup.bam -scale ${scale} > ${bedgraph_d}/${sample}.filt.srt.nodup.norm.bedgraph
  bedSort ${bedgraph_d}/${sample}.filt.srt.nodup.norm.bedgraph ${bedgraph_d}/${sample}.filt.srt.nodup.norm.bedgraph
  bedGraphToBigWig ${bedgraph_d}/${sample}.filt.srt.nodup.norm.bedgraph ${genome} ${bigwig_d}/${sample}.filt.srt.nodup.norm.bw
done

for sample in $(cat list/list_sample);do
  #Convert BAM to BED and adjust ends by Tn5 insertion position
  ##Tn5 adjustment: add 4 to '+' strand and subtract 5 from '-' strand
  bamToBed -i ${input_d}/${sample}.filt.nmsrt.nodup.bam | awk -F $'\t' 'BEGIN {OFS = FS}{ if ($6 == "+") {$2 = $2 + 4} else if ($6 == "-") {$3 = $3 - 5} print $0}' | \
      gzip -c > ${result_d}/${sample}.filt.nmsrt.nodup.tn5.bed.gz
done
