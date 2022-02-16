## Tally individual level missingness for UKBB GT data
smiss_sum<-data.table()
for (i in 1:22){
  smiss <- fread(file=paste0('qc_gt/chr',i,'/plink.imiss'))
  smiss_sum<-rbind(smiss_sum,smiss)
}
smiss_sum <- smiss_sum[,c(1,4,5)]
smiss_sum <-aggregate(. ~ FID, smiss_sum, sum)
smiss_sum$mis<-smiss_sum$N_MISS/smiss_sum$N_GENO
summary(smiss_sum$mis)

## List of individuals with > 10% missing genotypes
missing_indiv <- smiss_sum[smiss_sum$mis > 0.1,]
missing_indiv <- data.table(IID=missing_indiv$FID, FID=missing_indiv$FID)

## Write out
write.csv(missing_indiv,file='v2_gt_missing_indiv10.csv',row.names=F)