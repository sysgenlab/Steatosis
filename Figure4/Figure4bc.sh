up_bed=../data/DiffBind/Result.DiffBind.FvsS.upregulated.bed
down_bed=../data/DiffBind/Result.DiffBind.FvsS.downregulated.bed

up_C=$(wc -l ${up_bed} | sed 's/ .*//g' | numfmt --g)
down_C=$(wc -l ${down_bed} | sed 's/ .*//g' | numfmt --g)

#Draw Heatmap
computeMatrix reference-point \
  -S $(sed 's/^/browser\/bigwig\//' list/list_sample_FRiP20 | sed 's/$/.filt.srt.nodup.norm.bw/') \
  -R $up_bed $down_bed \
  --referencePoint center -p 40 -a 2000 -b 2000 \
  --skipZeros --missingDataAsZero \
  --samplesLabel $(sed -r 's/_Normal.*|_Steatosis.*|_fNASH.*//' list/list_sample_FRiP20) \
  -out ../data/computeMatrix/matrix.ATAC.macs2.DiffBind.FvsS.mat.gz

plotHeatmap -m ../data/computeMatrix/matrix.ATAC.macs2.DiffBind.FvsS.mat.gz \
  --colorMap $(echo Greens$_{1..1} Blues$_{1..3} Reds$_{1..10}) --refPointLabel "" --whatToShow "heatmap and colorbar" \
  --regionsLabel "Up-regulated regions (n=$up_C)" "Down-regulated regions (n=$down_C)" \
  --zMin 0 --zMax 0.5 \
  -out image/Figure4b.pdf

plotProfile -m ../data/computeMatrix/matrix.ATAC.macs2.DiffBind.FvsS.mat.gz \
  --colors $(echo green$_{1..1} blue$_{1..3} red$_{1..10}) --refPointLabel "" \
  --perGroup --plotWidth 8 --yMax 0.35 \
  --regionsLabel "Up-regulated regions (n=$up_C)" "Down-regulated regions (n=$down_C)" \
  -out image/Figure4c.pdf
