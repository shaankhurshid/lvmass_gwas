# Script for formatting acceleration data for GWAS

# Depends
library(data.table)

# Load ECG segmentation data
seg_lvm <- fread(file='lvm_from_inlinevf_and_ml4h_segmentation.tsv')

# Load hand-labeled values for recalibration
labeled <- fread(file='lvm_labeled.csv')

# Load censor data
censor <- fread(file='censor_202006.csv')

# Perform calibration to sex-specific mean
## Get sex
setkey(censor,sample_id); setkey(labeled,sample_id); setkey(seg_lvm,sample_id)
labeled[censor,sex := i.sex]
seg_lvm[censor,sex := i.sex]

## One person does not have sex (N = 44537 - 1 = 44536)
seg_lvm <- seg_lvm[!is.na(sex)]

## Perform recalibration
correction_f <- mean(seg_lvm[sex=='Female']$ml4h_lvm,na.rm=T) - mean(labeled[sex=='Female']$LVM,na.rm=T)
correction_m <- mean(seg_lvm[sex=='Male']$ml4h_lvm,na.rm=T) -  mean(labeled[sex=='Male']$LVM,na.rm=T)
seg_lvm[sex=="Female",lvm_seg_adjusted := ml4h_lvm-(correction_f)]
seg_lvm[sex=="Male",lvm_seg_adjusted := ml4h_lvm-(correction_m)]

## Get BMI to get LVMI
## Ht/wt
### Height
ht0 <- fread(file='height0.csv')
ht1 <- fread(file='height1.csv')
ht2 <- fread(file='height2.csv')

# Join
setkey(ht0,sample_id); setkey(ht1,sample_id); setkey(ht2,sample_id)
ht0[ht1,ht1 := i.height]
ht0[ht2,ht2 := i.height]
ht0[,":="(ht_unified = ifelse(!is.na(ht2),ht2,
                              ifelse(!is.na(ht1),ht1,height)))]

### Weight
wt0 <- fread(file='weight0.csv')
wt1 <- fread(file='weight1.csv')
wt2 <- fread(file='weight2.csv')

# Join
setkey(wt0,sample_id); setkey(wt1,sample_id); setkey(wt2,sample_id)
wt0[wt1,wt1 := i.weight]
wt0[wt2,wt2 := i.weight]
wt0[,":="(wt_unified = ifelse(!is.na(wt2),wt2,
                              ifelse(!is.na(wt1),wt1,weight)))]

# Join to main
setkey(ht0,sample_id); setkey(wt0,sample_id)
seg_lvm[ht0,ht_unified := i.ht_unified]
seg_lvm[wt0,wt_unified := i.wt_unified]

# Index the LV mass
seg_lvm[c(!is.na(ht_unified) & !is.na(wt_unified)),bsa := (0.007184 * ht_unified^0.725 * wt_unified^0.425)]
seg_lvm[!is.na(bsa),':='(lvmi_seg_adjusted = lvm_seg_adjusted / bsa)]

## Remove no BMI (N = 44536 - 63 = 44473)
seg_lvm <- seg_lvm[!is.na(bsa)]

## Outlier cleanup 
## First, remove < 20 (N=44473 - 89 = 44384)
seg_lvm <- seg_lvm[lvm_seg_adjusted > 20]

## Second, remove outside 5xIQR (N=44384 - 2 = 44382)
iqr <- as.numeric(quantile(seg_lvm$lvm_seg_adjusted,0.75) - quantile(seg_lvm$lvm_seg_adjusted,0.25))
seg_lvm <- seg_lvm[c((lvm_seg_adjusted > (median(lvm_seg_adjusted) - 5*(iqr)))
                     & (lvm_seg_adjusted < (median(lvm_seg_adjusted) + 5*(iqr)))),]

seg_lvm[!is.na(bsa),':='(lvmi_seg_adjusted_int = qnorm((rank(lvmi_seg_adjusted,na.last="keep")-0.5)/sum(!is.na(lvmi_seg_adjusted))))]
seg_lvm[!is.na(bsa),':='(lvmi_seg_adjusted_27 = lvm_seg_adjusted / ((ht_unified/100)^2.7))]

## Remove withdrawals (N=44382 - 7 = 44375)
withdrawals <- fread('w7089_20220222.csv') # UKBB withdrawals
seg_lvm <- seg_lvm[!(sample_id %in% withdrawals$V1)]

## Scope columns
seg_lvm <- seg_lvm[,c('sample_id','sex','lvm_seg_adjusted','lvmi_seg_adjusted','lvmi_seg_adjusted_27','lvmi_seg_adjusted_int')]

# Save out
write.csv(seg_lvm,'seg_lvm.csv',row.names=F,quote=F)

