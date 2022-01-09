# Script to process BOLT GWAS output

# Dependencies
library(data.table)
library(qqman)
library(R.utils)

# Read in data
gwas_filtered <- fread(file='lvmi_seg_adjusted_bolt.imputed.filtered.tsv')

# Format
gwas_filtered <- gwas_filtered[order(gwas_filtered$P_BOLT_LMM)]

results <- data.table(SNP = gwas_filtered$SNP, CHROM = gwas_filtered$CHR, 
                      POS = gwas_filtered$BP, P = gwas_filtered$P_BOLT_LMM)

## Manhattan
png("manhattan_plot_reg.png",
    width=3000, height=2250, res=300)
par(mar=c(3,3,1,1),oma=c(3,3,1,1))
print(manhattan(results,
                chr="CHROM",bp="POS",p="P",snp="SNP",
                suggestiveline = -log10(1e-6),
                genomewideline = -log10(5e-8),
                #col=reds.c,
                #cex=1.1,
                ylim = c(1.5, 11),
                annotatePval = -log10(5e-8),
                annotateTop = TRUE,
                main = "ML4H Reg Manhattan Plot"))
dev.off()

## QQ
png("seg_QQ_lvmi.png")
qq(results$P, main = "Q-Q plot of ML4H Seg LVM GWAS p-values : log")
dev.off()