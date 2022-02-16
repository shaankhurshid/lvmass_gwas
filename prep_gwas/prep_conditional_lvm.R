# Script to format SNP list to run a conditional analysis

# Dependencies
library(data.table)
library(stringr)

# Load SNPs
lead_snps <- fread(file='lvm_seg_adjusted_leadsnps_filtered.txt')

# Filter
lead_snps <- lead_snps[A1FREQ >= 0.01 & A1FREQ <= 0.99]
lead_snps <- lead_snps[str_detect(SNP,'rs')]

# Format for tool
lead_snps <- lead_snps[,c("SNP","CHR")]

# Write out
write.table(lead_snps,file='snps.list.lvm',sep='\t',row.names=F,col.names=F,quote=F)

# To bgen2bigquery 
#./bgen2bigquery.osx -snps snps.list.lvm -assembly "grch37" -bgen-template "gs://fc-7d5088b4-7673-45b5-95c2-17ae00a04183/imputed/ukb_imp_chr%s_v3.bgen" -sample "gs://ukbb_v2/projects/jamesp/projects/gwas/lv20k/data/ukb7089_imp_chr15_v3_s487395.sample" > long_out.txt

# Now load long output
long_out <- fread(file='long_out_lvm.txt')

# Load GWAS sample
lvm_gwas <- fread(file='lvm_seg_adjusted.tsv')

# Limit long output to individuals in GWAS sample
long_out <- long_out[sample_id %in% lvm_gwas$FID]
long_out <- long_out[,c('sample_id','rsid','alt_allele_dosage')]

# Widen
wide_out <- dcast(long_out,formula = sample_id ~ rsid,value.var='alt_allele_dosage')

# Merge with pheno file
setkey(lvm_gwas,FID); setkey(wide_out,sample_id)
lvm_gwas_conditional <- lvm_gwas[wide_out,nomatch=0]

# Write out
write.table(lvm_gwas_conditional,file='lvm_seg_adjusted_conditional.tsv',
          sep='\t',row.names=F,quote=F)
