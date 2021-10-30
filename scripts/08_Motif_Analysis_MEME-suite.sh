FASTA=/home/refseq/fasta/hg38/hg38_dna_UCSC/hg38.fa
jaspar=/data/program/meme-5.0.5/database/motif_databases/JASPAR/JASPAR2020_CORE_vertebrates_non-redundant_pfms.meme

mkdir -p MEME-suite/DREME MEME-suite/TOMTOM MEME-suite/AME MEME-suite/fasta
mkdir -p logs/DREME logs/TOMTOM logs/AME

for sample in FvsN FvsS SvsN;do
  for reg in downregulated upregulated;do
    BED_s=DiffBind/Result.DiffBind.$sample.$reg.bed
    FASTA_s=MEME-suite/fasta/Result.DiffBind.$sample.$reg.fasta
    
    bedtools getfasta -fi $FASTA -bed $BED_s -name > $FASTA_s

    dreme -oc MEME-suite/DREME/$sample.$reg -p $FASTA_s -dna -m 20 -maxk 16 -png &> logs/DREME/$sample.$reg.DREME.log
    tomtom -oc MEME-suite/TOMTOM/$sample.$reg MEME-suite/DREME/$sample.$reg/dreme.txt $jaspar \
      -min-overlap 5 -dist pearson -thresh 0.05 -no-ssc &> logs/TOMTOM/$sample.$reg.TOMTOM.log 
    ame -oc MEME-suite/AME/$sample.$reg -control --shuffle-- -scoring avg -method fisher $FASTA_s $jaspar &> logs/AME/$sample.$reg.AME.log 
  done
done
