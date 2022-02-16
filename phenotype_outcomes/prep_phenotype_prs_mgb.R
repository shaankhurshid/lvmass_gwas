# Script to assess association between LVM PRS and outcomes in MGB

# Depends
library(data.table)

# Load PRS
prs <- fread(file='lvm_prs_all_com.txt')

# Define outcomes
## Load comorb definitions
comorbs <- fread(file='comorbidities.tsv')

## Isolate comorbs to diseases of interest
af_dia <- comorbs[cov=='af_aflt' & cov.code.type %in% c('ICD9','ICD10')]$cov.code
af_prc <- comorbs[cov=='af_aflt' & cov.code.type %in% c('CPT')]$cov.code
mi_dia <- comorbs[cov=='mi' & cov.code.type %in% c('ICD9','ICD10')]$cov.code
hf_dia <- comorbs[cov=='heartFailure' & cov.code.type %in% c('ICD9','ICD10')]$cov.code
vt_dia <- comorbs[cov=='vt' & cov.code.type %in% c('ICD9','ICD10')]$cov.code
icd_prc <- comorbs[cov=='defib' & cov.code.type %in% c('CPT')]$cov.code
dcm_dia <- c('I42.0')
hcm_dia <- c('I42.1','I42.2','425.11','425.18')

## Load long dia file
load(file='mgb_dia_genetic.RData')

## Load long prc file
load(file='mgb_prc_genetic.RData')

# Format Date_Age column
dia_prs[,Date_Age := as.numeric(str_extract(Date_Age,'\\d+'))]
prc_prs[,Date_Age := as.numeric(str_extract(Date_Age,'\\d+'))]

## Define lookup function
icd_lookup <- function(query,dia_list,id_list,disease_name){
  setkey(id_list,linker_id)
  all_match <- dia_list[Code %in% query]
  first_match <- all_match[,.SD[which.min(Date_Age)],by='linker_id']
  setkey(first_match,linker_id)
  id_list[first_match,paste0(disease_name) := i.Date_Age]
  return(id_list)
}

# Run function
prs <- icd_lookup(query=af_dia,dia_list=dia_prs,id_list=prs,disease_name='af')
prs <- icd_lookup(query=mi_dia,dia_list=dia_prs,id_list=prs,disease_name='mi')
prs <- icd_lookup(query=vt_dia,dia_list=dia_prs,id_list=prs,disease_name='vt')
prs <- icd_lookup(query=hf_dia,dia_list=dia_prs,id_list=prs,disease_name='hf')
prs <- icd_lookup(query=dcm_dia,dia_list=dia_prs,id_list=prs,disease_name='dcm')
prs <- icd_lookup(query=hcm_dia,dia_list=dia_prs,id_list=prs,disease_name='hcm')
prs <- icd_lookup(query=af_prc,dia_list=prc_prs,id_list=prs,disease_name='af_prc')
prs <- icd_lookup(query=icd_prc,dia_list=prc_prs,id_list=prs,disease_name='icd')

# Unify AF (earlier of dia or prc)
prs[,af_combined := pmin(af,af_prc,na.rm=T)]

# Days to years
for (j in c('af_combined','mi','hf','vt','dcm','hcm','icd')){
  set(prs,j=j,value=prs[[j]]/365.25)
}

# Death and last enc
## Can steal this from accel prs
accel <- fread(file='accel_prs_mgb.csv')
setkey(accel,linker_id); setkey(prs,linker_id)
prs[accel,':='(last_enc = i.last_enc, censor_age = i.censor_age)]

# Scope columns
prs <- prs[,.SD,.SDcols=names(prs)[!(names(prs) %in% c('af','af_prc'))]]
setnames(prs,'af_combined','af')

# Save out
write.csv(prs,file='lvm_prs_phenos.csv',row.names=F)

