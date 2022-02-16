# Depends
library(data.table)

# Load phenotypes
seg_prs <- fread(file='prs_no_mass_processed.csv')
censor <- fread(file='censor_202106.csv')
race <- fread(file='race_0.csv')
mi <- fread(file='mi_icd_202106.tsv')
hf <- fread(file='hf_inclusive_202106.tsv')
htn <- fread(file='htn_202106.csv')
dm <- fread(file='any_dm_202106.csv')
af <- fread(file='af_202106.csv')
instance2 <- fread(file='instance2_dates.csv')
sbp <- fread(file='sbp_combined_instance02_202006.csv')
dbp <- fread(file='dbp_combined_instance02_202006.csv')

# MRI date per person
instance2[,value := as.Date(value,format='%Y-%m-%d')]

# Join 
setkey(seg_prs,sample_id); setkey(instance2,sample_id); setkey(censor,sample_id); setkey(dm,sample_id)
setkey(mi,sample_id); setkey(htn,sample_id); setkey(hf,sample_id); setkey(af,sample_id)
setkey(sbp,sample_id); setkey(dbp,sample_id)

seg_prs[instance2,mri_date := i.value]
seg_prs[censor,':='(enroll_date = i.enroll_date,birthdate = i.birthdate,
                    enroll_age = i.enroll_age,
                    sex = i.sex)]

seg_prs[mi,':=' (prev_mi = i.prevalent_disease)]
seg_prs[hf,':=' (prev_hf = i.prevalent_disease)]
seg_prs[htn,':=' (prev_htn = i.prevalent_disease)]
seg_prs[dm,':=' (prev_dm = i.prevalent_disease)]

# Granular race
race[,race_categories := ifelse(value %in% c(1,1001,1002,1003),'white',
                                ifelse(value %in% c(2,2001,2002,2003,2004),'mixed',
                                       ifelse(value %in% c(3,3001,3002,3003,3004,5),'asian',
                                              ifelse(value %in% c(4,4001,4002,4003),'black',
                                                     ifelse(value %in% 6,'other',NA)))))]
setkey(race,sample_id)
seg_prs[race,race_category := i.race_categories]

# Save
save(seg_prs,file='prs_no_mass_pheno.RData')
