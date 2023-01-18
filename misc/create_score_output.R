# Script to create score output file for PRS transfer

# Depends
library(data.table)

# Create function to provide standardized score output for PRS replication
create_score <- function(gwas,clump,output){
  
  # Load filtered GWAS
  gwas <- fread(file=paste0(gwas))
  
  # Load score
  score_id <- fread(file=paste0(clump),header=F)
  
  # Merge
  score <- gwas[SNP %in% score_id$V1]
  
  # Isolate duplicates
  duplicated <- score[SNP %in% score[duplicated(SNP)]$SNP]
  duplicated[,BETA := lapply(.SD,mean),by='SNP',.SDcols='BETA']
  special <- unique(duplicated,by='SNP')
  
  # Add back to master dataset
  score <- score[!(SNP %in% score[duplicated(SNP)]$SNP)]
  score <- rbind(score,special)
  
  # Select relevant columns
  score_output <- score[,c('SNP','CHR','BP','ALLELE1','ALLELE0','A1FREQ','BETA')]
  
  # Write out
  write.table(score_output,file=paste0(output),sep='\t',quote=F,row.names=F)
}

create_score(gwas='lvmi_seg_adjusted_bolt.imputed.filtered.tsv',
             clump='bolt_clumped.tsv',
             output='prs_output.tsv')