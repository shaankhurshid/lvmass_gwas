# Depends
library(data.table)

# Load phenotypes
seg_lvm <- fread(file='seg_lvm.csv')
censor <- fread(file='censor_202106.csv')
race <- fread(file='race_0.csv')
mi <- fread(file='mi_icd_202106.tsv')
hf <- fread(file='hf_inclusive_202106.tsv')
htn <- fread(file='htn_202106.csv')
dm <- fread(file='any_dm_202106.csv')
af <- fread(file='af_202106.csv')
instance2 <- fread(file='instance2_dates.csv')
sbp <- fread(file='bp_combined_instance02_202006.csv')
dbp <- fread(file='bp_combined_instance02_202006.csv')

# MRI date per person
instance2[,value := as.Date(value,format='%Y-%m-%d')]

# Join 
setkey(seg_lvm,sample_id); setkey(instance2,sample_id); setkey(censor,sample_id); setkey(dm,sample_id)
setkey(mi,sample_id); setkey(htn,sample_id); setkey(hf,sample_id); setkey(af,sample_id)
setkey(sbp,sample_id); setkey(dbp,sample_id)

seg_lvm[instance2,mri_date := i.value]
seg_lvm[censor,':='(enroll_date = i.enroll_date,birthdate = i.birthdate,enroll_age = i.enroll_age)]

seg_lvm[mi,':=' (updated_mi = i.has_disease,updated_mi_date = i.censor_date)]
seg_lvm[hf,':=' (updated_hf = i.has_disease,updated_hf_date = i.censor_date)]
seg_lvm[htn,':=' (updated_htn = i.has_disease,updated_htn_date = i.censor_date)]
seg_lvm[dm,':=' (updated_dm = i.has_disease,updated_dm_date = i.censor_date)]
seg_lvm[af,':=' (updated_af = i.has_disease,updated_af_date = i.censor_date)]
seg_lvm[sbp,':=' (sbp = i.sbp_combined)]
seg_lvm[dbp,':=' (dbp = i.dbp_combined)]

# Date casting
for (j in c('mri_date','enroll_date',
            'updated_mi_date','updated_hf_date','updated_htn_date',
            'birthdate','updated_dm_date','updated_af_date'))
{set(seg_lvm,j=j,value=as.Date(seg_lvm[[j]],format='%Y-%m-%d'))}

seg_lvm[,':='(updated_htn_age = as.numeric(updated_htn_date - birthdate)/365.25,
                    updated_hf_age = as.numeric(updated_hf_date - birthdate)/365.25,
                    updated_mi_age = as.numeric(updated_mi_date - birthdate)/365.25,
                    updated_af_age = as.numeric(updated_af_date - birthdate)/365.25,
                    updated_dm_age = as.numeric(updated_dm_date - birthdate)/365.25)]

# Age at MRI
seg_lvm[,age_at_mri := as.numeric(enroll_age + (mri_date - enroll_date)/365.25)]

# Prevalent defs
seg_lvm[c((updated_mi==1) & (updated_mi_age <= age_at_mri)),prev_mi_mri := 1]
seg_lvm[c((updated_hf==1) & (updated_hf_age <= age_at_mri)),prev_hf_mri := 1]
seg_lvm[c((updated_htn==1) & (updated_htn_age <= age_at_mri)),prev_htn_mri := 1]
seg_lvm[c((updated_dm==1) & (updated_dm_age <= age_at_mri)),prev_dm_mri := 1]
seg_lvm[c((updated_af==1) & (updated_af_age <= age_at_mri)),prev_af_mri := 1]

# Granular race
race[,race_categories := ifelse(value %in% c(1,1001,1002,1003),'white',
                                ifelse(value %in% c(2,2001,2002,2003,2004),'mixed',
                                       ifelse(value %in% c(3,3001,3002,3003,3004,5),'asian',
                                              ifelse(value %in% c(4,4001,4002,4003),'black',
                                                     ifelse(value %in% 6,'other',NA)))))]
setkey(race,sample_id)
seg_lvm[race,race_category := i.race_categories]

## Ensure exclusions removed
withdrawals <- fread('w7089_20210809.csv') # UKBB withdrawals
seg_lvm <- seg_lvm[!(sample_id %in% withdrawals$V1)]

# Save
save(seg_lvm,file='seg_lvm_pheno_111921.RData')

############## GWAS set
# Load GWAS set
gwas_set <- fread(file='lvmi_seg_adjusted.tsv')

# Reduce phenotype set to included in GWAS
seg_lvm_gwas <- seg_lvm[sample_id %in% gwas_set$FID]

# Save out
save(seg_lvm_gwas,file='seg_lvm_gwas_pheno.RData')
