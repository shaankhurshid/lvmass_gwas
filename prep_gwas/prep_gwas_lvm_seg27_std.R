  # Script for formatting acceleration data for GWAS
  
  # Depends
  library(data.table)
  library(sktools)
  
  # Load ECG segmentation data
  seg_lvm <- fread(file='seg_lvm_std.csv')
  
  # Load censor file (for age)
  censor <- fread(file='censor_202106.csv')
  
  # Load instance2 data (for age at MRI)
  instance2 <- fread(file='instance2_dates.csv')
  
  # Join to get MRI date
  setkey(instance2,sample_id);setkey(seg_lvm,sample_id)
  seg_lvm[instance2,':='(mri_date = i.value)]
  
  # Join to get birthdate
  setkey(censor,sample_id)
  seg_lvm[censor,':='(birthdate = i.birthdate)]
  
  # Remove withdrawals
  withdrawals <- fread('w7089_20220222.csv') # UKBB withdrawals
  
  seg_lvm <- seg_lvm[!(sample_id %in% withdrawals$V1)]
  
  # Remove in plink but not imputed
  not_imputed <- fread('bolt.in_plink_but_not_imputed.FID_IID.968.txt') # UKBB withdrawals
  seg_lvm <- seg_lvm[!(sample_id %in% not_imputed$V1)]
  
  # Keep only in PLINK
  gt <- fread(file='in_v2_gt.csv')
  seg_lvm <- seg_lvm[(sample_id %in% gt$FID)]
  
  # Remove individual missing GT (>0.1)
  missing_indiv <- fread(file='v2_gt_missing_indiv10.csv')
  seg_lvm <- seg_lvm[!(sample_id %in% missing_indiv$IID)]
  
  # Age at MRI
  format_date(seg_lvm,cols=c('mri_date','birthdate'))
  seg_lvm[,age_at_mri := as.numeric(mri_date - birthdate)/365.25]
  
  ## PLINK output function
  create<-function(trait,exclude_all_both=NULL,exclude_all_cases=NULL,exclude_all_controls=NULL,
                   exclude_incident_cases=NULL,exclude_flexible=NULL,data){
    ##phenotype file
    a <- data
    a<-a[!is.na(get(trait))]
    print(paste0('Total N:',nrow(a)))
    
    ##sample qc file
    b<-fread("ukb_sqc_v2_7089.tsv",header=T) 
    setkey(a,'sample_id'); setkey(b,'eid')
    ab <- a[b,nomatch=0]
    ab[,':='(array_UKBB = ifelse(genotyping_array=='UKBB',1,0))]
    print(paste0('N after merge with sample QC file:',nrow(ab)))
    
    ##remove poor quality
    ab[,':='(ex_poor = ifelse(het_missing_outliers==1 | putative_sex_chromosome_aneuploidy==1,1,0),
             ex_sex = ifelse(Submitted_Gender==Inferred_Gender,0,1),
             ex_misKin = ifelse(ab$excluded_from_kinship_inference==1,1,0))]
    
    #high quality data
    ab <- ab[ab$ex_sex==0]
    print(paste0('N after removal of sex mismatch:',nrow(ab)))
    ab <- ab[ab$ex_poor==0]
    print(paste0('N after removal of poor:',nrow(ab)))
    #ab <- ab[ab$ex_misKin==0]
    #print(paste0('N after removal of missing kinship inference:',nrow(ab)))
    
    # Loop over "exclude all both" phenotypes - all individuals with exclusion phenotype at any time removed for both cases/controls
    if (length(exclude_all_both)!=0){
      for (i in exclude_all_both){
        exclude <- fread(paste0("/phenotypes/",i,'.tab.tsv'),header=T)
        setkey(ab,sample_id); setkey(exclude,sample_id)
        ab[exclude,':='(exclude_prev=i.prevalent_disease,exclude_incd=i.incident_disease,exclude_censor = i.censor_age)]
        ab[,exclude := ifelse(c(c(!is.na(exclude_prev) & exclude_prev==1) | c(!is.na(exclude_incd) & exclude_incd == 1)),1,0)]
        print(paste0('I am going to exclude ',sum(ab$exclude),' individuals for diagnosis: ',i,' occurring at any time for cases/controls'))
        ab <- ab[exclude==0]
        ab <- ab[,!(c('exclude_prev','exclude_incd','exclude_censor','exclude'))]
      }}
    
    # Loop over "exclude all cases" phenotypes - all individuals with exclusion phenotype at any time removed for cases
    if (length(exclude_all_cases)!=0){
      for (i in exclude_all_cases){
        exclude <- fread(paste0("/phenotypes/",i,'.tab.tsv'),header=T)
        setkey(ab,sample_id); setkey(exclude,sample_id)
        ab[exclude,':='(exclude_prev=i.prevalent_disease,exclude_incd=i.incident_disease,exclude_censor = i.censor_age)]
        ab[,exclude := ifelse(c(c(!is.na(get(trait)) & get(trait)==1) & 
                                  c(c(!is.na(exclude_prev) & exclude_prev==1) | c(!is.na(exclude_incd) & exclude_incd==1))),1,0)]
        print(paste0('I am going to exclude ',sum(ab$exclude),' individuals for diagnosis: ',i,' occurring at any time for cases only'))
        ab <- ab[exclude==0]
        ab <- ab[,!(c('exclude_prev','exclude_incd','exclude_censor','exclude'))]
      }}
    
    # Loop over "exclude all controls" phenotypes - all individuals with exclusion phenotype at any time removed for controls
    if (length(exclude_all_controls)!=0){
      for (i in exclude_all_controls){
        exclude <- fread(paste0("/phenotypes/",i,'.tab.tsv'),header=T)
        setkey(ab,sample_id); setkey(exclude,sample_id)
        ab[exclude,':='(exclude_prev=i.prevalent_disease,exclude_incd=i.incident_disease,exclude_censor = i.censor_age)]
        ab[,exclude := ifelse(c(c(get(trait)==0 | is.na(get(trait))) & 
                                  c(c(!is.na(exclude_prev) & exclude_prev==1) | c(!is.na(exclude_incd) & exclude_incd==1))),1,0)]
        print(paste0('I am going to exclude ',sum(ab$exclude),' individuals for diagnosis: ',i,' occurring at any time for controls only'))
        ab <- ab[exclude==0]
        ab <- ab[,!(c('exclude_prev','exclude_incd','exclude_censor','exclude'))]
      }}
    
    # Loop over "exclude incident" - only cases with exclusion phenotype before disease removed
    if (length(exclude_incident_cases)!=0){
      for (i in exclude_incident_cases){
        exclude <- fread(paste0("/phenotypes/",i,'.tab.tsv'),header=T)
        setkey(ab,sample_id); setkey(exclude,sample_id)
        ab[exclude,':='(exclude_disease = i.has_disease, exclude_prev = i.prevalent_dsease, exclude_censor = i.censor_date)]
        ab[,exclude := ifelse(c(c(!is.na(get(trait)) & get(trait)==1) & 
                                  c(c(!is.na(exclude_disease) & (exclude_censor <= censor_date)) |
                                      c(!is.na(exclude_disease) & exclude_prev==1))),1,0)]
        print(paste0('I am going to exclude ',sum(ab$exclude),' individuals for diagnosis: ',i,' occurring before case diagnosis'))
        ab <- ab[exclude==0]
        ab <- ab[,!(c('exclude_prev','exclude_censor','exclude_disease','exclude'))]
      }}
    
    # Loop over "exclude flexible" - excludes any instance of exclusion phenotype among controls, and only exclusion phenotype prior to disease for cases
    if (length(exclude_flexible)!=0){
      for (i in exclude_flexible){
        exclude <- fread(paste0("/phenotypes/",i,'.tab.tsv'),header=T)
        setkey(ab,sample_id); setkey(exclude,sample_id)
        ab[exclude,':='(exclude_incd = i.incident_disease, exclude_disease = i.has_disease, exclude_prev = i.prevalent_disease, exclude_censor = i.censor_date)]
        ab[,exclude := ifelse(c(!is.na(get(trait)) & get(trait)==1),
                              ifelse(c(!is.na(exclude_prev) & exclude_prev==1),1,
                                     ifelse(c(!is.na(exclude_incd) & (exclude_incd==1) & (exclude_censor <= censor_date)),1,0)),
                              ifelse(c(!is.na(exclude_disease) & exclude_disease==1),1,0))]
        print(paste0('I am going to exclude ',sum(ab$exclude),' individuals for diagnosis: ',i,' occurring before case diagnosis or at any time for controls'))
        ab <- ab[exclude==0]
        ab <- ab[,!(c('exclude_disease','exclude_prev','exclude_incd','exclude_censor','exclude'))]
      }}
    
    ab[,':='(male = ifelse(Inferred_Gender=='M',1,0))]
    
    #######
    ###test AF related PCs in each cleaned dataset
    #######
    
    #dim(subset(ab, used_in_pca_calculation==1 & ex_sex==0 & ex_poor==0 & ex_misKin==0 & white==1))
    #######
    ##all white, no relatives
    #######
    form1<-formula(paste0(trait,"~age_at_mri + ",paste0("PC",1:40,collapse="+"),"+ array_UKBB + male",collapse="+"))
    s1<-summary(lm(form1,data=ab))$coefficients
    s1<-s1[substring(rownames(s1),1,2)=="PC",]
    ab.1<-ab[,.SD,.SDcols=c("age_at_mri",trait)]
    allN_1<-nrow(ab.1)
    male_1<-nrow(ab[ab[["male"]]==1 & !is.na(ab.1[[trait]])])
    pcs1<-paste(rownames(s1)[1:5],collapse=",") # First 5 PCs
    
    #######
    ##create summary file
    #######
    t1<-c("all",allN_1,"mean_age",round(mean(ab[!is.na(trait)]$age_at_mri),2),"sd_age",round(sd(ab[!is.na(trait)]$age_at_mri),2),"male_N",male_1,"male%",round(mean(ab[!is.na(trait)]$male)*100,2),"related-PCs",ifelse((pcs1!=""),pcs1,"None"))
    write.table(t1,file=paste0("/summary_",trait,".txt"),row.names=F,quote=F,sep="\t")
    
    #######
    ##create phenotype file
    #######
    ## Choose columns
    pheno<-ab[,c("sample_id",trait,"age_at_mri",rownames(s1)[1:5],"array_UKBB","male"),with=F]
    ## Format for PLINK
    setnames(pheno,"sample_id","FID")
    pheno[,':='(IID = FID)]
    setcolorder(pheno,c('FID','IID'))
    print(paste0('Final phenotype N: ',nrow(pheno)))
    write.table(pheno,file=paste0('/lvm_gwas/',trait,".tsv"),sep="\t",col.names =T,row.names = F,quote = F)
  }
  
  create(trait="lvmi_seg_adjusted_27",data=seg_lvm)
  
  # Exclusions list for BOLT
  ## Read processed phenotype
  seg_lvm <- fread(file='lvmi_seg_adjusted_27.tsv')
  ## Load list of all individuals in UKBB (includes people in plink GT set that are not in censor files)
  ukbb_all <- fread(file='ukbb_all.csv')
  exclusions <- ukbb_all[!(FID %in% seg_lvm$FID)]
  ## Add individuals that are in the plink dataset but not the censor file
  write.table(exclusions,file='lvmi_seg_exclusion_27.tsv',sep="\t",col.names =T,row.names = F,quote = F)
  
