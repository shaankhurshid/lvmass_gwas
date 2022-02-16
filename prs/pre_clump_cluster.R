# Script to prep datasets for PRS building using PLINK

# Dependencies
library(data.table)
library(stringr)

# This version depends on the presence of a 1000G EUR LD file that has already been pre-processed

pre_clump <- function(results_dir,input_gwas,MAF_min){
  # STEP 1: QC the GWAS
  raw_gwas <- fread(paste0(results_dir,input_gwas,'.tsv'))
  
  # Apply MAF filter
  gwas_filtered <- raw_gwas[c((A1FREQ >= MAF_min) & (A1FREQ <= (1-MAF_min)))]
  
  # Save out
  write.table(gwas_filtered,file=paste0(results_dir,input_gwas,'.filtered.tsv'),sep='\t',row.names=F,quote=F)
  
  # Remove rows with ambiguous variants
  gwas_unambiguous <- gwas_filtered[!c((ALLELE1=='A' & ALLELE0 == "T") | 
                                       (ALLELE1=='T' & ALLELE0 == "A") |
                                       (ALLELE1=='C' & ALLELE0 == "G") |
                                       (ALLELE1=='G' & ALLELE0 == "C"))]
  
  # Remove indels
  gwas_unambiguous <- gwas_unambiguous[c((ALLELE1 %in% c('A','T','C','G')) & (ALLELE0 %in% c('A','T','C','G')))]
  
  # Save out
  write.table(gwas_unambiguous,file=paste0(results_dir,input_gwas,'.filtered.unambiguous.noindel.tsv'),sep='\t',row.names=F,quote=F)
 
  # STEP 2: Script to prep datasets for PRS building using PLINK
  gwas_unambiguous[,':='(ID = paste0(CHR,':',BP,':',ALLELE0,':',ALLELE1))]
  
  formatted_output <- gwas_unambiguous[,c('ID','P_BOLT_LMM')]
  
  write.table(formatted_output,file=paste0(results_dir,input_gwas,'.filtered.unambiguous.noindel.forplinkprs.tsv'),sep='\t',row.names=F,quote=F)
}

pre_clump(results_dir='lvm/',
          input_gwas='lvm_seg_adjusted_bolt.imputed',MAF_min=0.01)