mkdir -p DiffBind/sheet

sheet=DiffBind/sheet/samplesheet.DiffBind.FRiP20.csv
echo SampleID,Tissue,Condition,Replicate,bamReads,Peaks,PeakCaller > $sheet

N_N=0
N_S=0
N_F=0

for sample in $(cat list/list_sample_FRiP20);do
  Tissue=Liver
  Condition=$(echo $sample | sed 's/^.*._N/N/' | sed 's/^.*._S/S/' | sed 's/^.*._f/f/' | sed 's/_.*//')
  bamReads=bowtie2/Filtering/${sample}.filt.srt.nodup.bam
  Peaks=macs2/${sample}.tn5_peaks.bed
  PeakCaller=narrow

  if [ $Condition == "Normal" ];then
    let "N_N=N_N+1"
    Rep=$N_N
  elif [ $Condition == "Steatosis" ];then
    let "N_S=N_S+1"
    Rep=$N_S
  elif [ $Condition == "fNASH" ];then
    let "N_F=N_F+1"
    Rep=$N_F
  fi

  echo $sample,$Tissue,$Condition,$Rep,$bamReads,$Peaks,$PeakCaller
done >> $sheet
