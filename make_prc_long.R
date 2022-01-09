# Script to generate giant long file of prc codes for phecode processing

# Depends
library(data.table)

# Load file lists
prs <- fread('accel_prs_mgb.csv')
prs_eur <- fread('accel_prs_mgb_eur.csv')

# Read in prc flat files
out <- list()
for (i in 1:6){
  current_path <- paste0("phb_flat_files/phb",i,'/')
  prc <- fread(paste0(current_path,list.files(current_path)[str_detect(list.files(current_path),'Prc')]),
               sep='|',quote="")
  out[[i]] <- prc[Code_Type=='CPT',c('linker_id','Code','Code_Type','Date_Age')]
  print(paste0('I just finished loading prc file ',i,' out of 6!'))
}
all_prc <- do.call(rbind,out)
setDT(all_prc)

# Create sets corresponding to those with a) genetic data, and b) genetic data + EUR
## All
save(all_prc,file='mgb_all_prc.RData')

## PRS
prc_prs <- all_prc[linker_id %in% prs$linker_id]
save(prc_prs,file='mgb_prc_genetic.RData')

## PRS-EUR
prc_prs_eur <- all_prc[linker_id %in% prs_eur$linker_id]
save(prc_prs_eur,file='mgb_prc_genetic_eur.RData')
