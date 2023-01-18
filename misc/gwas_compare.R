# Depends
library(data.table)
library(plyr)
library(stringr)

# Load GWAS files
seg <- fread(file='v27_lvmi_seg_adjusted_std_bolt.imputed.tsv')
ts <- fread(file='v28_lvmi_seg_adjusted_27_std_bolt.imputed.tsv')
int <- fread(file='v22_lvmi_seg_adjusted_int_bolt.imputed.filtered.tsv')
lean <- fread(file='v29_lvmi_nfmass_std_bolt.imputed.tsv')

# Load final SNP list
snps <- fread(file='seg_snp_list.csv')

# Reduce
seg <- seg[SNP %in% snps$SNP]
ts <- ts[SNP %in% snps$SNP]
int <- int[SNP %in% snps$SNP]
lean <- lean[SNP %in% snps$SNP]

# Set keys
setkey(seg,SNP); setkey(ts,SNP); setkey(int,SNP); setkey(lean,SNP); setkey(snps,SNP)

# Merge
seg[int,':='(int_beta = i.BETA, int_se = i.SE, int_p = i.P_BOLT_LMM)]
seg[ts,':='(ts_beta = i.BETA, ts_se = i.SE, ts_p = i.P_BOLT_LMM)]
seg[lean,':='(lean_beta = i.BETA, lean_se = i.SE, lean_p = i.P_BOLT_LMM)]

# Save out
write.csv(seg,file='gwas_compare.csv',row.names=F)
