library(gprofiler2)
library(sessioninfo)
library(dplyr)
library(readr)
library(TxDb.Mmusculus.UCSC.mm10.ensGene)
library(tximport)
library("DESeq2")
library(EnhancedVolcano)
library("pheatmap")
library(ggplot2)
library(eulerr)


## list all installed packages
pkgs = installed.packages()[,'Package']

## read description files & parse using sessioninfo
desc = lapply(pkgs, utils::packageDescription)
source = vapply(desc, sessioninfo:::pkg_source, character(1))

## combine and filter for BioC only
bioc_pkgs = data_frame(pkgs, source) %>%
  filter(source == "Bioconductor")
bioc_pkgs

samples = read.csv("reads_treatments.csv", header=TRUE)
files = file.path("quants_prebuilt", samples$Run, "quant.sf")
txdb = TxDb.Mmusculus.UCSC.mm10.ensGene
k = keys(txdb, keytype = "TXNAME")
tx2gene = select(txdb, k, "GENEID", "TXNAME")
txi = tximport(files, type = "salmon", tx2gene = tx2gene, ignoreTxVersion=TRUE)
ddsTxi = DESeqDataSetFromTximport(txi, colData = samples, design = ~ Treatment)
txsByGene = transcriptsBy(txdb,"gene")
lengthData = median(width(txsByGene))

#PF4_EGFplus vs. NaCl_EGFplus
keep = rowSums(counts(ddsTxi)) >= 10
dds = ddsTxi[keep,]
dds$Treatment = relevel(dds$Treatment, ref = "NaCl_EGFplus")
dds = DESeq(dds)
res = results(dds, contrast=c("Treatment","PF4_EGFplus","NaCl_EGFplus"))
summary(res)
resOrdered = res[order(res$pvalue),]
write.csv(as.data.frame(resOrdered), 
          file="paper_analysis/PF4_EGFplus_vs_NaCl_EGFplus.csv")

# Volcano plot
resLFC = lfcShrink(dds, coef=4, type="ashr")
keyvals = ifelse(
  resLFC$log2FoldChange < -2 & resLFC$padj<0.05 , 'royalblue',
  ifelse(resLFC$log2FoldChange > 2 & resLFC$padj<0.05, 'red',
         'grey'))
keyvals[is.na(keyvals)] <- 'grey'
names(keyvals)[keyvals == 'royalblue'] <- 'LFC<-2 & Padj<0.05'
names(keyvals)[keyvals == 'grey'] <- 'Other'
names(keyvals)[keyvals == 'red'] <- 'LFC>2 & Padj<0.05'

png(file="paper_analysis/PF4_EGFplus vs. NaCl_EGFplus volcano plot.png",
    width=20, height=30, res=600, units="cm")
EnhancedVolcano(resLFC,
                lab = rownames(resLFC),
                x = 'log2FoldChange',
                y = 'padj',
                title = 'PF4_EGFplus vs. NaCl_EGFplus volcano plot',
                pCutoff = 0.05,
                FCcutoff = 2,
                gridlines.major = FALSE,
                gridlines.minor = FALSE,
                labSize=0,
                pointSize=4,
                colCustom = keyvals,
                caption = 'LFC cutoff, 2; p-adj cutoff, 0.05',
                )
dev.off()

#MA plot
png(file="paper_analysis/PF4_EGFplus vs. NaCl_EGFplus MA plot.png",
    width=15, height=10, res=600, units="cm")
pma = plotMA(resLFC, ylim=c(-28,28), returnData = TRUE)
pma$condition = with(pma, ifelse(lfc < 0 & isDE==TRUE, "-1", 
                                 ifelse(lfc>0 & isDE==TRUE, "1", "0")))
pmagg = ggplot(pma, aes(x=mean, y=lfc))+geom_point(size=1, aes(color=factor(condition)))+scale_x_continuous(trans='log10')
pmagg+ theme_classic() +scale_color_manual(values = c("royalblue", "grey", "red")) + labs(x="Mean abundance", y="Log fold change") + theme(legend.position = "none")
dev.off()

#Heatmap
res05 = results(dds, contrast=c("Treatment","PF4_EGFplus","NaCl_EGFplus"), alpha=0.05)
res05 = res05[!is.na(res05$padj), ]
res05 = res05[res05$padj<0.05, ]
summary(res05)

ntd = normTransform(dds)
sig_ntd=ntd[rownames(res05),]
genes = order(res05$log2FoldChange, decreasing=TRUE)[1:nrow(res05)]

png(file="paper_analysis/PF4_EGFplus vs. NaCl_EGFplus heatmap.png",
    width=10, height=16, res=600, units="cm")
pheatmap(assay(sig_ntd)[genes,c(7,8,9,10,11,12,1,3,4,5,2,6)], 
         cluster_rows=FALSE,show_rownames=FALSE, show_colnames = FALSE,
         cluster_cols=FALSE, fontsize = 7, 
         color=colorRampPalette(c("royalblue", "white", "red"))(50))
dev.off()


#Enrichment analysis

#Up-regulated genes (LFC>2)
upreg = as.data.frame(res05[res05$log2FoldChange > 2, ])
plus_upreg_genes = row.names(upreg)
write_lines(plus_upreg_genes, 
            "paper_analysis/PF4_EGFplus_v_NaCl_EGFplus_upreg.txt")
gostres_up = gost(query = plus_upreg_genes, organism = "mmusculus", 
                  significant=TRUE, user_threshold=0.05)
gostres_up_plot = gostplot(gostres_up, capped = FALSE)
gostres_up_plot

#Down-regulated genes (LFC<-2)
downreg = as.data.frame(res05[res05$log2FoldChange < -2, ])
plus_downreg_genes = row.names(downreg)
write_lines(plus_downreg_genes, 
            "paper_analysis/PF4_EGFplus_v_NaCl_EGFplus_downreg.txt")
gostres_down = gost(query = plus_downreg_genes, organism = "mmusculus", 
                    significant=TRUE, user_threshold=0.05)
gostres_down_plot = gostplot(gostres_down, capped = FALSE)
gostres_down_plot

#PF4_EGFminus vs. NaCl_EGFminus
keep = rowSums(counts(ddsTxi)) >= 10
dds = ddsTxi[keep,]
dds$Treatment = relevel(dds$Treatment, ref = "NaCl_EGFminus")
dds = DESeq(dds)
res = results(dds, contrast=c("Treatment","PF4_EGFminus","NaCl_EGFminus"))
summary(res)
resOrdered = res[order(res$pvalue),]
write.csv(as.data.frame(resOrdered), 
          file="paper_analysis/PF4_EGFminus_vs_NaCl_EGFminus.csv")

# Volcano plot
resLFC = lfcShrink(dds, coef=3, type="ashr")
keyvals = ifelse(
  resLFC$log2FoldChange < -2 & resLFC$padj<0.05 , 'royalblue',
  ifelse(resLFC$log2FoldChange > 2 & resLFC$padj<0.05, 'red',
         'grey'))
keyvals[is.na(keyvals)] <- 'grey'
names(keyvals)[keyvals == 'royalblue'] <- 'LFC<-2 & Padj<0.05'
names(keyvals)[keyvals == 'grey'] <- 'Other'
names(keyvals)[keyvals == 'red'] <- 'LFC>2 & Padj<0.05'

png(file="paper_analysis/PF4_EGFminus vs. NaCl_EGFminus volcano plot.png",
    width=20, height=30, res=600, units="cm")           
EnhancedVolcano(resLFC,
               lab = rownames(resLFC),
               x = 'log2FoldChange',
               y = 'padj',
               pCutoff = 0.05,
               FCcutoff = 2,
               gridlines.major = FALSE,
               gridlines.minor = FALSE,
               labSize=0,
               pointSize=4,
               colCustom = keyvals,
               caption = 'LFC cutoff, 2; p-adj cutoff, 0.05',
                )
dev.off()

#MA plot
png(file="paper_analysis/PF4_EGFminus vs. NaCl_EGFminus MA plot.png",
    width=15, height=10, res=600, units="cm")
pma = plotMA(resLFC, ylim=c(-28,28), returnData = TRUE)
pma$condition = with(pma, ifelse(lfc < 0 & isDE==TRUE, "-1", 
                                 ifelse(lfc>0 & isDE==TRUE, "1", "0")))
pmagg = ggplot(pma, aes(x=mean, y=lfc))+geom_point(size=1, aes(color=factor(condition)))+scale_x_continuous(trans='log10')
pmagg+ theme_classic() +scale_color_manual(values = c("royalblue", "grey", "red")) + labs(x="Mean abundance", y="Log fold change") + theme(legend.position = "none")
dev.off()

#Heatmap
res05 = results(dds, contrast=c("Treatment","PF4_EGFminus","NaCl_EGFminus"), 
                alpha=0.05)
res05 = res05[!is.na(res05$padj), ]
res05 = res05[res05$padj<0.05, ]
summary(res05)

ntd = normTransform(dds)
sig_ntd=ntd[rownames(res05),]
genes = order(res05$log2FoldChange, decreasing=TRUE)[1:nrow(res05)]

png(file="paper_analysis/PF4_EGFminus vs. NaCl_EGFminus heatmap.png",
    width=10, height=16, res=600, units="cm")
pheatmap(assay(sig_ntd)[genes,c(19,20,21,22,23,24,13,14,15,16,17,18)], 
         cluster_rows=FALSE,show_rownames=FALSE, show_colnames = FALSE,
         cluster_cols=FALSE, fontsize = 7, 
         color=colorRampPalette(c("royalblue", "white", "red"))(50))
dev.off()

#PCA plot
png(file="paper_analysis/PCA.png",
    width=10, height=8, res=600, units="cm")
z = plotPCA(ntd, intgroup=c("Treatment"), ntop = 500, returnData = TRUE)
y = ggplot(z, aes(x= PC1, y = PC2))+geom_point(size= 2, aes(colour=Treatment), alpha=0.8, stroke = 0.4) 
y + theme_classic() + theme(text = element_text(size = 8)) + xlim(-150,150) +ylim(-100,100) +scale_color_manual(values = c("lightgrey", "lightblue","grey27","royalblue"))
dev.off()

#Enrichment analysis

#Up-regulated genes (LFC>2)
upreg = as.data.frame(res05[res05$log2FoldChange > 2, ])
minus_upreg_genes = row.names(upreg)
write_lines(minus_upreg_genes, 
            "paper_analysis/PF4_EGFminus_v_NaCl_EGFminus_upreg.txt")
gostres_up = gost(query = minus_upreg_genes, organism = "mmusculus", 
                  significant=TRUE, user_threshold=0.05)
gostres_up_plot = gostplot(gostres_up, capped = FALSE)
gostres_up_plot

#Down-regulated genes (LFC<-2)
downreg = as.data.frame(res05[res05$log2FoldChange < -2, ])
minus_downreg_genes = row.names(downreg)
write_lines(minus_downreg_genes, 
            "paper_analysis/PF4_EGFminus_v_NaCl_EGFminus_downreg.txt")
gostres_down = gost(query = minus_downreg_genes, organism = "mmusculus", 
                    significant=TRUE, user_threshold=0.05)

#Intersection plus - minus - upreg
intersect_upreg = intersect(plus_upreg_genes, minus_upreg_genes)
intersect = length(intersect_upreg)
plus_upreg = length(plus_upreg_genes)-length(intersect_upreg)
minus_upreg = length(minus_upreg_genes)-length(intersect_upreg)
fit = euler(c("PF4-EGFplus up"=plus_upreg, 
              "PF4-EGFminus up"=minus_upreg,
              "PF4-EGFplus up&PF4-EGFminus up"=intersect
                        ))
png(file="paper_analysis/euler_upreg.png",
    width=16, height=13, res=600, units="cm")
plot(fit,shape = "ellipse",
     quantities = F,
     fill = "transparent",
     lty = 1:3,
     labels = "")
dev.off()

#Intersection plus - minus - downreg
intersect_downreg = intersect(plus_downreg_genes, minus_downreg_genes)
intersect = length(intersect_downreg)
plus_downreg = length(plus_downreg_genes)-length(intersect_downreg)
minus_downreg = length(minus_downreg_genes)-length(intersect_downreg)
fit = euler(c("PF4-EGFplus down"=plus_downreg, 
              "PF4-EGFminus down"=minus_downreg,
              "PF4-EGFplus down&PF4-EGFminus down"=intersect
))
png(file="paper_analysis/euler_downreg.png",
    width=18, height=14, res=600, units="cm")
plot(fit,shape = "ellipse",
     quantities = F,
     fill = "transparent",
     lty = 1:3,
     labels = "")
dev.off()
