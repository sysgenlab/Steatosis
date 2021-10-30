options(stringsAsFactors=FALSE)
library(ggplot2)

RPKM <- read.delim(paste0("../data/RPKM.FRiP20.tsv"))
myRPKM <- as.matrix(RPKM[,-c(1,2,3)])
FRiP <- read.delim("Summary/Summary_FRiP", row.names=1)
pca <- prcomp(t(log10(myRPKM+0.01)))

df <- as.data.frame(pca$x)
df$FRiP <- as.numeric(FRiP[3,rownames(df)])
type <- sub("X.*_S", "S", sub("X.*_N", "N", sub("X.*_f", "f", rownames(df))))
df$type <- factor(type, levels=c("Normal_1", "Steatosis_1", "Steatosis_2", "fNASH_1", "fNASH_2"))

xlab <- paste0("PC 1 (", round(100*summary(pca)$importance[2,1]), "%)")
ylab <- paste0("PC 2 (", round(100*summary(pca)$importance[2,2]), "%)")

p <- ggplot(df, aes(x=PC1, y=PC2, color=type))
p <- p + geom_point(size=5)
p <- p + theme_bw()
p <- p + labs(x=xlab, y=ylab)
p <- p + guides(color=guide_legend(title=""))
p <- p + scale_color_manual(values=c("gray", "blue", "navy", "red", "brown"))
p <- p + theme(axis.title.x = element_text(size=14, face="bold"),
                           axis.title.y = element_text(size=14, face="bold"),
                           axis.text.x = element_text(size=11),
                           axis.text.y = element_text(size=11),
                           legend.text = element_text(size=11))

pdf("imageA/Figure3a.pdf", width=7, height=5)
print(p)
dev.off()
