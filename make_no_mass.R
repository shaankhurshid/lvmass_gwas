# Script to create cohort of interest in UKBB WITHOUT CMR

# Load CMR LV mass list
all_lvm <- fread(file='lvmi_seg_adjusted.tsv')

# Load all UKBB
all_ukbb <- fread(file='censor_202106.csv')

# Remove exclusions
exclusions <- fread(file='w7089_20210809.csv')
all_ukbb <- all_ukbb[!(sample_id %in% exclusions$V1)]

# Remove LV mass
all_ukbb <- all_ukbb[!(sample_id %in% all_lvm$FID)]

# Remove failure of QC
qc <-fread("ukb_sqc_v2_7089.tsv",header=T) 
setkey(all_ukbb,'sample_id'); setkey(qc,'eid')
all_ukbb <- all_ukbb[qc,nomatch=0]
all_ukbb[,':='(array_UKBB = ifelse(genotyping_array=='UKBB',1,0))]
print(paste0('N after merge with sample QC file:',nrow(all_ukbb)))

##remove poor quality
all_ukbb[,':='(ex_poor = ifelse(het_missing_outliers==1 | putative_sex_chromosome_aneuploidy==1,1,0),
               ex_sex = ifelse(Submitted_Gender==Inferred_Gender,0,1),
               ex_misKin = ifelse(excluded_from_kinship_inference==1,1,0))]

#high quality data
all_ukbb <- all_ukbb[all_ukbb]
print(paste0('N after removal of sex mismatch:',nrow(all_ukbb))) 
all_ukbb <- all_ukbb[all_ukbb$ex_poor==0]
print(paste0('N after removal of poor:',nrow(all_ukbb)))

# remove individual level missingness filter
missing <- fread(file='v2_gt_missing_indiv10.csv')
all_ukbb <- all_ukbb[!(sample_id %in% missing$IID)]

# Save out
## Seed file
write.table(all_ukbb,file='no_mass_seed.tsv',quote=F,row.names=F,
            sep='\t')

## List of IDs for plink
no_mass_ids <- all_ukbb[,'sample_id']
setnames(no_mass_ids,'sample_id','FID')
no_mass_ids[,':='(IID = FID)]
write.table(no_mass_ids,file='no_mass_list.tsv',quote=F,row.names=F,
            sep='\t')