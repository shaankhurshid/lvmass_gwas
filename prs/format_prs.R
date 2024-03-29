# Process PRS results

# Dependencies
library(data.table)
library(stringr)

# Read in results
path_to_prs <- '/prs_no_mass/'
score_files <- list.files(path_to_prs)[str_detect(list.files(path_to_prs),'sscore')]

score_combined <- list()
for (file in score_files){
  name <- str_extract(file,'chr\\d+')
  score_combined[[name]] <- fread(paste0(path_to_prs,file))}

# Merge into a unified wide file
current <- score_combined[[1]]; setkey(current,'#IID')
setnames(current,'NAMED_ALLELE_DOSAGE_SUM','dose_1'); setnames(current,'SCORE1_SUM','score_1')

for (num in 2:length(score_combined)){
  next_table <- score_combined[[num]]
  setkey(next_table,'#IID')
  current[next_table,c(paste0('dose_',num),paste0('score_',num)) := list(i.NAMED_ALLELE_DOSAGE_SUM,i.SCORE1_SUM)]
}

# Now add to get overall score
score_names <- names(current)[str_detect(names(current),"score")]
current[,sum_score := apply(.SD, FUN=sum, MARGIN=1),.SDcols=score_names]

# Format and save out
prs <- current[,c('#IID','sum_score')]
setnames(prs,'#IID','sample_id')
prs[,':='(prs_std = (sum_score - mean(sum_score)) / sd(sum_score),
          prs_sum = sum_score)]

write.csv(prs,file=paste0(path_to_prs,'prs_no_mass_processed.csv'),row.names=F)