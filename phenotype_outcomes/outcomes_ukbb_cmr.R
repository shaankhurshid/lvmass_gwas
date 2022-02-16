# Script to assess outcomes in UKBB using CMR truth only

# Dependencies
library(data.table)
library(survival)
library(prodlim)
library(plyr)
library(Cairo)

# Load files
censor <- fread(file='censor_202106.csv')
af <- fread(file='af_202106.csv')
mi <- fread(file='mi_icd_202106.tsv')
hf <- fread(file='hf_inclusive_202106.tsv')
vt <- fread(file='vt_202106.csv')
hcm <- fread(file='hcm_202106.csv')
dcm <- fread(file='dcm_202106.csv')
icd <- fread(file='icd_202106.csv')
instance2 <- fread(file='instance2_dates.csv')

# Load MRI data
seg_lvm <- fread(file='seg_lvm.csv')

# Censor at death or last follow-up for non-events
af[,censor_date := as.Date(ifelse(has_disease==0 & has_died == 1,
                                  pmin(death_censor_date,censor_date),censor_date),origin='1970-01-01')]
mi[,censor_date := as.Date(ifelse(has_disease==0 & has_died == 1,
                                  pmin(death_censor_date,censor_date),censor_date),origin='1970-01-01')]
hf[,censor_date := as.Date(ifelse(has_disease==0 & has_died == 1,
                                  pmin(death_censor_date,censor_date),censor_date),origin='1970-01-01')]
vt[,censor_date := as.Date(ifelse(has_disease==0 & has_died == 1,
                                  pmin(death_censor_date,censor_date),censor_date),origin='1970-01-01')]
dcm[,censor_date := as.Date(ifelse(has_disease==0 & has_died == 1,
                                   pmin(death_censor_date,censor_date),censor_date),origin='1970-01-01')]
hcm[,censor_date := as.Date(ifelse(has_disease==0 & has_died == 1,
                                   pmin(death_censor_date,censor_date),censor_date),origin='1970-01-01')]
icd[,censor_date := as.Date(ifelse(has_disease==0 & has_died == 1,
                                   pmin(death_censor_date,censor_date),censor_date),origin='1970-01-01')]

# Generate phenotypes (based on instance 2)
## Clean up and prejoin
setkey(af,sample_id); setkey(mi,sample_id); setkey(seg_lvm,sample_id)
setkey(dcm,sample_id); setkey(hcm,sample_id); setkey(icd,sample_id)
setkey(hf,sample_id); setkey(censor,sample_id); setkey(vt,sample_id); setkey(instance2,sample_id)

setnames(af,'has_disease','has_af'); setnames(af,'incident_disease','incident_af'); 
setnames(af,'prevalent_disease','prevalent_af'); setnames(af,'censor_date','af_censor_date')

## Joins
af[mi,":="(has_mi = i.has_disease, prevalent_mi = i.prevalent_disease, 
           incident_mi = i.incident_disease, mi_censor_date = i.censor_date)]
af[hf,":="(has_hf = i.has_disease, prevalent_hf = i.prevalent_disease, 
           incident_hf = i.incident_disease, hf_censor_date = i.censor_date)]
af[vt,":="(has_vt = i.has_disease, prevalent_vt = i.prevalent_disease, 
           incident_vt = i.incident_disease, vt_censor_date = i.censor_date)]
af[icd,":="(has_icd = i.has_disease, prevalent_icd = i.prevalent_disease, 
           incident_icd = i.incident_disease, icd_censor_date = i.censor_date)]
af[dcm,":="(has_dcm = i.has_disease, prevalent_dcm = i.prevalent_disease, 
           incident_dcm = i.incident_disease, dcm_censor_date = i.censor_date)]
af[hcm,":="(has_hcm = i.has_disease, prevalent_hcm = i.prevalent_disease, 
           incident_hcm = i.incident_disease, hcm_censor_date = i.censor_date)]
af[censor,":="(enroll_date = i.enroll_date,death_date = i.death_date, death_age = i.death_age,
               death_censor_date = i.death_censor_date)]
af[instance2,":="(mri_date = i.value)]

## Ensure exclusions removed
withdrawals <- fread('w7089_20210809.csv') # UKBB withdrawals
all_data <- af[!(sample_id %in% withdrawals$V1)]

# Fix censoring dates
# Format dates
dates <- c('birthdate','af_censor_date','mi_censor_date','vt_censor_date',
           'dcm_censor_date','hcm_censor_date','icd_censor_date',
           'hf_censor_date','enroll_date','mri_date')
for (j in dates){set(all_data,j=j,value=as.Date(all_data[[j]],format='%Y-%m-%d'))}

# Load center categories
center <- fread(file='center0.csv')
center_lookup <- fread(file='enrollment_correspondences.csv')

# Add center value to dataset
setkey(all_data,sample_id); setkey(center,sample_id)
all_data[center,':='(center_code = i.value)]

setkey(all_data,center_code); setkey(center_lookup,Code)
all_data[center_lookup,':='(center_location = i.Region)]

# Now correct censor dates based on location
all_data[,':='(af_censor_date = as.Date(ifelse(c(has_af==1 | center_location=='England'),af_censor_date,
                                               ifelse(center_location=='Scotland',pmin(af_censor_date,as.Date('2021-03-31',format='%Y-%m-%d')),
                                                      pmin(af_censor_date,as.Date('2018-02-28',format='%Y-%m-%d')))),origin='1970-01-01'))]
all_data[,':='(mi_censor_date = as.Date(ifelse(c(has_mi==1 | center_location=='England'),mi_censor_date,
                                               ifelse(center_location=='Scotland',pmin(mi_censor_date,as.Date('2021-03-31',format='%Y-%m-%d')),
                                                      pmin(mi_censor_date,as.Date('2018-02-28',format='%Y-%m-%d')))),origin='1970-01-01'))]
all_data[,':='(hf_censor_date = as.Date(ifelse(c(has_hf==1 | center_location=='England'),hf_censor_date,
                                               ifelse(center_location=='Scotland',pmin(hf_censor_date,as.Date('2021-03-31',format='%Y-%m-%d')),
                                                      pmin(hf_censor_date,as.Date('2018-02-28',format='%Y-%m-%d')))),origin='1970-01-01'))]
all_data[,':='(vt_censor_date = as.Date(ifelse(c(has_vt==1 | center_location=='England'),vt_censor_date,
                                               ifelse(center_location=='Scotland',pmin(vt_censor_date,as.Date('2021-03-31',format='%Y-%m-%d')),
                                                      pmin(vt_censor_date,as.Date('2018-02-28',format='%Y-%m-%d')))),origin='1970-01-01'))]
all_data[,':='(icd_censor_date = as.Date(ifelse(c(has_icd==1 | center_location=='England'),icd_censor_date,
                                               ifelse(center_location=='Scotland',pmin(icd_censor_date,as.Date('2021-03-31',format='%Y-%m-%d')),
                                                      pmin(icd_censor_date,as.Date('2018-02-28',format='%Y-%m-%d')))),origin='1970-01-01'))]
all_data[,':='(dcm_censor_date = as.Date(ifelse(c(has_dcm==1 | center_location=='England'),dcm_censor_date,
                                               ifelse(center_location=='Scotland',pmin(dcm_censor_date,as.Date('2021-03-31',format='%Y-%m-%d')),
                                                      pmin(dcm_censor_date,as.Date('2018-02-28',format='%Y-%m-%d')))),origin='1970-01-01'))]
all_data[,':='(hcm_censor_date = as.Date(ifelse(c(has_hcm==1 | center_location=='England'),hcm_censor_date,
                                               ifelse(center_location=='Scotland',pmin(hcm_censor_date,as.Date('2021-03-31',format='%Y-%m-%d')),
                                                      pmin(hcm_censor_date,as.Date('2018-02-28',format='%Y-%m-%d')))),origin='1970-01-01'))]

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
all_data[ht0,ht_unified := i.ht_unified]
all_data[wt0,wt_unified := i.wt_unified]

## Create specific variables for analysis
# Disease
all_data[!is.na(mri_date),":="(incd_af_after_mri = ifelse(c(!is.na(incident_af) & (incident_af==1) 
                                                            & (af_censor_date > mri_date)),1,0),
                               incd_mi_after_mri = ifelse(c(!is.na(incident_mi) & (incident_mi==1) 
                                                            & (mi_censor_date > mri_date)),1,0),
                               incd_vt_after_mri = ifelse(c(!is.na(incident_vt) & (incident_vt==1) 
                                                            & (vt_censor_date > mri_date)),1,0),
                               incd_hf_after_mri = ifelse(c(!is.na(incident_hf) & (incident_hf==1) 
                                                            & (hf_censor_date > mri_date)),1,0),
                               incd_dcm_after_mri = ifelse(c(!is.na(incident_dcm) & (incident_dcm==1) 
                                                            & (dcm_censor_date > mri_date)),1,0),
                               incd_hcm_after_mri = ifelse(c(!is.na(incident_hcm) & (incident_hcm==1) 
                                                            & (hcm_censor_date > mri_date)),1,0),
                               incd_icd_after_mri = ifelse(c(!is.na(incident_icd) & (incident_icd==1) 
                                                            & (icd_censor_date > mri_date)),1,0))]

# Time
all_data[!is.na(mri_date),":="(time_mri_to_af = as.numeric(af_censor_date - mri_date)/365.25,
                               time_mri_to_mi = as.numeric(mi_censor_date - mri_date)/365.25,
                               time_mri_to_vt = as.numeric(vt_censor_date - mri_date)/365.25,
                               time_mri_to_hf = as.numeric(hf_censor_date - mri_date)/365.25,
                               time_mri_to_icd = as.numeric(icd_censor_date - mri_date)/365.25,
                               time_mri_to_hcm = as.numeric(hcm_censor_date - mri_date)/365.25,
                               time_mri_to_dcm = as.numeric(dcm_censor_date - mri_date)/365.25)]

# Age at MRI
all_data[,age_at_mri := as.numeric(mri_date - birthdate)/365.25]

# LVH variables
## Join the inference values
setkey(all_data,sample_id); setkey(seg_lvm,sample_id)
all_data[seg_lvm, ':='(lvm = i.lvm_seg_adjusted,
                             lvmi = i.lvmi_seg_adjusted,
                             sex = i.sex)]

## BSA
all_data[c(!is.na(ht_unified) & !is.na(wt_unified)),bsa := (0.007184 * ht_unified^0.725 * wt_unified^0.425)]

# LVH variables
all_data[!is.na(lvmi),lvmi_std := (lvmi-mean(lvmi))/sd(lvmi)]

all_data[c(!is.na(sex) & !is.na(lvmi)),
         ':='(lvh = ifelse(c((sex=='Female') & (lvmi > 55)),1,
                                     ifelse(lvmi > 72,1,0)))]

all_data[!is.na(lvmi),':='(lvh90 = ifelse(sex=='Female',
                                                         ifelse(lvmi > quantile(all_data[sex=='Female']$lvmi,probs=0.90,na.rm=T),1,0),
                                                         ifelse(lvmi > quantile(all_data[sex=='Male']$lvmi,probs=0.90,na.rm=T),1,0)))]

# Incident disease associations
incident_data <- all_data[!is.na(lvmi)]

# AF
mod_af <- coxph(Surv(time_mri_to_af,incd_af_after_mri) ~ lvmi_std + age_at_mri + sex,data=incident_data[time_mri_to_af>0])
mod_af_lvh <- coxph(Surv(time_mri_to_af,incd_af_after_mri) ~ lvh + age_at_mri + sex,data=incident_data[time_mri_to_af>0])
mod_af_lvh90 <- coxph(Surv(time_mri_to_af,incd_af_after_mri) ~ lvh90 + age_at_mri + sex,data=incident_data[time_mri_to_af>0])

# Graphical variables
incident_data[,inferred_lvh_graphical := factor(ifelse(lvh==0,"No inferred LVH","Inferred LVH"),levels=c('Inferred LVH','No inferred LVH'))]
incident_data[,inferred_lvh_graphical90 := factor(ifelse(lvh90==0,"No inferred LVH","Inferred LVH"),levels=c('Inferred LVH','No inferred LVH'))]

# AF KM
prodlim_af <- prodlim(Hist(time_mri_to_af,incd_af_after_mri)~inferred_lvh_graphical,data=incident_data[time_mri_to_af>0])

CairoPDF(file='km_af_ukbb.pdf',height=3,width=3.5,
         pointsize=4)
par(oma=c(3,1,1,1),mar=c(3,1,1,1))
plot(prodlim_af,"cuminc",ylim=c(0,0.04),xlim=c(0,4),
     axis2.at=seq(0,0.04,0.01),axis2.las=2,lwd=1.5,background=F,
     axis1.at=c(0,1,2,3,4),axis1.labels=as.character(0:4),
     atrisk.times=c(0,1,2,3,4,5),col=c("#fc4e2a",'#1c9099'),atrisk.col='black',confint=FALSE,
     legend.x=0,legend.y=0.04,axis1.cex.axis=2.5,axis2.cex.axis=2.5,axis1.padj=0.5,
     legend.cex=2.2,legend.legend=c("CMR-derived LVH","No CMR-derived LVH"),
     atrisk.title=("                     "),atrisk.pos=-0.5,atrisk.line=c(1.2,2.5),
     atrisk.cex=1.8,atrisk.interspace=1.4,xlab='',ylab='')
mtext("Cumulative risk (%)",side=2,line=-1.2,at=0.02,cex=2.5)
mtext("Years",side=1, line=-0.6,cex=2.5)
mtext('Stratum',side=1, line=-0.3,cex=1.8,at=-0.92)
dev.off()

# MI
mod_mi <- coxph(Surv(time_mri_to_mi,incd_mi_after_mri) ~ lvmi_std + age_at_mri + sex,data=incident_data[time_mri_to_mi>0])
mod_mi_lvh <- coxph(Surv(time_mri_to_mi,incd_mi_after_mri) ~ lvh + age_at_mri + sex,data=incident_data[time_mri_to_mi>0])
mod_mi_lvh90 <- coxph(Surv(time_mri_to_mi,incd_mi_after_mri) ~ lvh90 + age_at_mri + sex,data=incident_data[time_mri_to_mi>0])

# MI KM
prodlim_mi <- prodlim(Hist(time_mri_to_mi,incd_mi_after_mri)~inferred_lvh_graphical,data=incident_data[time_mri_to_mi>0])

CairoPDF(file='km_mi_ukbb.pdf',height=3,width=3.5,
         pointsize=4)
par(oma=c(3,1,1,1),mar=c(3,1,1,1))
plot(prodlim_mi,"cuminc",ylim=c(0,0.04),xlim=c(0,4),
     axis2.at=seq(0,0.04,0.01),axis2.las=2,lwd=1.5,background=F,
     axis1.at=c(0,1,2,3,4),axis1.labels=as.character(0:4),
     atrisk.times=c(0,1,2,3,4,5),col=c("#fc4e2a",'#1c9099'),atrisk.col='black',confint=FALSE,
     legend.x=0,legend.y=0.04,axis1.cex.axis=2.5,axis2.cex.axis=2.5,axis1.padj=0.5,
     legend.cex=2.2,legend.legend=c("CMR-derived LVH","No CMR-derived LVH"),
     atrisk.title=("                     "),atrisk.pos=-0.5,atrisk.line=c(1.2,2.5),
     atrisk.cex=1.8,atrisk.interspace=1.4,xlab='',ylab='')
mtext("Cumulative risk (%)",side=2,line=-1.2,at=0.02,cex=2.5)
mtext("Years",side=1, line=-0.6,cex=2.5)
mtext('Stratum',side=1, line=-0.3,cex=1.8,at=-0.92)
dev.off()

# HF
mod_hf <- coxph(Surv(time_mri_to_hf,incd_hf_after_mri) ~ lvmi_std + age_at_mri + sex,data=incident_data[time_mri_to_hf>0])
mod_hf_lvh <- coxph(Surv(time_mri_to_hf,incd_hf_after_mri) ~ lvh + age_at_mri + sex,data=incident_data[time_mri_to_hf>0])
mod_hf_lvh90 <- coxph(Surv(time_mri_to_hf,incd_hf_after_mri) ~ lvh90 + age_at_mri + sex,data=incident_data[time_mri_to_hf>0])

# HF KM
prodlim_hf <- prodlim(Hist(time_mri_to_hf,incd_hf_after_mri)~inferred_lvh_graphical,data=incident_data[time_mri_to_hf>0])

CairoPDF(file='km_hf_ukbb.pdf',height=3,width=3.5,
         pointsize=4)
par(oma=c(3,1,1,1),mar=c(3,1,1,1))
plot(prodlim_hf,"cuminc",ylim=c(0,0.04),xlim=c(0,4),
     axis2.at=seq(0,0.04,0.01),axis2.las=2,lwd=1.5,background=F,
     axis1.at=c(0,1,2,3,4),axis1.labels=as.character(0:4),
     atrisk.times=c(0,1,2,3,4,5),col=c("#fc4e2a",'#1c9099'),atrisk.col='black',confint=FALSE,
     legend.x=0,legend.y=0.04,axis1.cex.axis=2.5,axis2.cex.axis=2.5,axis1.padj=0.5,
     legend.cex=2.2,legend.legend=c("CMR-derived LVH","No CMR-derived LVH"),
     atrisk.title=("                     "),atrisk.pos=-0.5,atrisk.line=c(1.2,2.5),
     atrisk.cex=1.8,atrisk.interspace=1.4,xlab='',ylab='')
mtext("Cumulative risk (%)",side=2,line=-1.2,at=0.02,cex=2.5)
mtext("Years",side=1, line=-0.6,cex=2.5)
mtext('Stratum',side=1, line=-0.3,cex=1.8,at=-0.92)
dev.off()

# VT
mod_vt <- coxph(Surv(time_mri_to_vt,incd_vt_after_mri) ~ lvmi_std + age_at_mri + sex,data=incident_data[time_mri_to_vt>0])
mod_vt_lvh <- coxph(Surv(time_mri_to_vt,incd_vt_after_mri) ~ lvh + age_at_mri + sex,data=incident_data[time_mri_to_vt>0])
mod_vt_lvh90 <- coxph(Surv(time_mri_to_vt,incd_vt_after_mri) ~ lvh90 + age_at_mri + sex,data=incident_data[time_mri_to_vt>0])

# VT KM
prodlim_vt <- prodlim(Hist(time_mri_to_vt,incd_vt_after_mri)~inferred_lvh_graphical,data=incident_data[time_mri_to_vt>0])

CairoPDF(file='km_vt_ukbb.pdf',height=3,width=3.5,
         pointsize=4)
par(oma=c(3,1,1,1),mar=c(3,1,1,1))
plot(prodlim_vt,"cuminc",ylim=c(0,0.04),xlim=c(0,4),
     axis2.at=seq(0,0.04,0.01),axis2.las=2,lwd=1.5,background=F,
     axis1.at=c(0,1,2,3,4),axis1.labels=as.character(0:4),
     atrisk.times=c(0,1,2,3,4),col=c("#fc4e2a",'#1c9099'),atrisk.col='black',confint=FALSE,
     legend.x=0,legend.y=0.04,axis1.cex.axis=2.5,axis2.cex.axis=2.5,axis1.padj=0.5,
     legend.cex=2.2,legend.legend=c("CMR-derived LVH","No CMR-derived LVH"),
     atrisk.title=("                     "),atrisk.pos=-0.5,atrisk.line=c(1.2,2.5),
     atrisk.cex=1.8,atrisk.interspace=1.4,xlab='',ylab='')
mtext("Cumulative risk (%)",side=2,line=-1.2,at=0.02,cex=2.5)
mtext("Years",side=1, line=-0.6,cex=2.5)
mtext('Stratum',side=1, line=-0.3,cex=1.8,at=-0.92)
dev.off()

# DCM
mod_dcm <- coxph(Surv(time_mri_to_dcm,incd_dcm_after_mri) ~ lvmi_std + age_at_mri + sex,data=incident_data[time_mri_to_dcm>0])
mod_dcm_lvh <- coxph(Surv(time_mri_to_dcm,incd_dcm_after_mri) ~ lvh + age_at_mri + sex,data=incident_data[time_mri_to_dcm>0])
mod_dcm_lvh90 <- coxph(Surv(time_mri_to_dcm,incd_dcm_after_mri) ~ lvh90 + age_at_mri + sex,data=incident_data[time_mri_to_dcm>0])

# DCM KM
prodlim_dcm <- prodlim(Hist(time_mri_to_dcm,incd_dcm_after_mri)~inferred_lvh_graphical,data=incident_data[time_mri_to_dcm>0])

CairoPDF(file='km_dcm_ukbb.pdf',height=3,width=3.5,
         pointsize=4)
par(oma=c(3,1,1,1),mar=c(3,1,1,1))
plot(prodlim_dcm,"cuminc",ylim=c(0,0.04),xlim=c(0,4),
     axis2.at=seq(0,0.04,0.01),axis2.las=2,lwd=1.5,background=F,
     axis1.at=c(0,1,2,3,4),axis1.labels=as.character(0:4),
     atrisk.times=c(0,1,2,3,4),col=c("#fc4e2a",'#1c9099'),atrisk.col='black',confint=FALSE,
     legend.x=0,legend.y=0.04,axis1.cex.axis=2.5,axis2.cex.axis=2.5,axis1.padj=0.5,
     legend.cex=2.2,legend.legend=c("CMR-derived LVH","No CMR-derived LVH"),
     atrisk.title=("                     "),atrisk.pos=-0.5,atrisk.line=c(1.2,2.5),
     atrisk.cex=1.8,atrisk.interspace=1.4,xlab='',ylab='')
mtext("Cumulative risk (%)",side=2,line=-1.2,at=0.02,cex=2.5)
mtext("Years",side=1, line=-0.6,cex=2.5)
mtext('Stratum',side=1, line=-0.3,cex=1.8,at=-0.92)
dev.off()

# HCM
mod_hcm <- coxph(Surv(time_mri_to_hcm,incd_hcm_after_mri) ~ lvmi_std + age_at_mri + sex,data=incident_data[time_mri_to_hcm>0])
mod_hcm_lvh <- coxph(Surv(time_mri_to_hcm,incd_hcm_after_mri) ~ lvh + age_at_mri + sex,data=incident_data[time_mri_to_hcm>0])
mod_hcm_lvh90 <- coxph(Surv(time_mri_to_hcm,incd_hcm_after_mri) ~ lvh90 + age_at_mri + sex,data=incident_data[time_mri_to_hcm>0])

# HCM KM
prodlim_hcm <- prodlim(Hist(time_mri_to_hcm,incd_hcm_after_mri)~inferred_lvh_graphical,data=incident_data[time_mri_to_hcm>0])

CairoPDF(file='km_hcm_ukbb.pdf',height=3,width=3.5,
         pointsize=4)
par(oma=c(3,1,1,1),mar=c(3,1,1,1))
plot(prodlim_hcm,"cuminc",ylim=c(0,0.04),xlim=c(0,4),
     axis2.at=seq(0,0.04,0.01),axis2.las=2,lwd=1.5,background=F,
     axis1.at=c(0,1,2,3,4),axis1.labels=as.character(0:4),
     atrisk.times=c(0,1,2,3,4),col=c("#fc4e2a",'#1c9099'),atrisk.col='black',confint=FALSE,
     legend.x=0,legend.y=0.04,axis1.cex.axis=2.5,axis2.cex.axis=2.5,axis1.padj=0.5,
     legend.cex=2.2,legend.legend=c("CMR-derived LVH","No CMR-derived LVH"),
     atrisk.title=("                     "),atrisk.pos=-0.5,atrisk.line=c(1.2,2.5),
     atrisk.cex=1.8,atrisk.interspace=1.4,xlab='',ylab='')
mtext("Cumulative risk (%)",side=2,line=-1.2,at=0.02,cex=2.5)
mtext("Years",side=1, line=-0.6,cex=2.5)
mtext('Stratum',side=1, line=-0.3,cex=1.8,at=-0.92)
dev.off()

# Spline
hcm_spline <- coxph(Surv(time_mri_to_hcm,incd_hcm_after_mri) ~ pspline(lvmi_std,df=0,caic=T) + age_at_mri + sex,data=incident_data[time_mri_to_hcm>0])

# Plot spline
resid_hcm_spline <- termplot(hcm_spline,term=1,se=TRUE,plot=F)
value_term <- resid_hcm_spline$lvmi_std
center <- with(value_term, y[x==min(abs((x-0)))])
y <- value_term$y + outer(value_term$se,c(0,-1.96,1.96),'*')

CairoPDF(file='hcm_cmr_spline.pdf',height=3,width=3.5,
         pointsize=4)
par(oma=c(3,4,1,1),mar=c(3,4,1,1),xpd=FALSE)
matplot(value_term$x,exp(y-center),log='y',type='l',xaxt='n',xlab='',ylab='',bty='n',yaxt='n',
        xlim=c(-4,4),ylim=c(0.01,100),col=c('black','#d73027','#1a9850'))
axis(1,at=seq(-4,4,1),cex.axis=1.6)
axis(2,las=2,cex.axis=1.6,pos=-4.2,at=c(0.01,0.1,1,10,100),labels=c('0.01','0.1','1','10','100'))
legend(1,0.05,c('estimate','upper 95% CI','lower 95% CI'),bty='n',lty=1,
       col=c('black','#d73027','#1a9850'),cex=1.4)
mtext(expression(paste("Standardized LVMI (g/m"^"2",")")),1,line=3.5,cex=1.8)
mtext("Relative hazard for HCM compared to average LVMI",2,line=5.5,cex=1.8)
dev.off()

# ICD
mod_icd <- coxph(Surv(time_mri_to_icd,incd_icd_after_mri) ~ lvmi_std + age_at_mri + sex,data=incident_data[time_mri_to_icd>0])
mod_icd_lvh <- coxph(Surv(time_mri_to_icd,incd_icd_after_mri) ~ lvh + age_at_mri + sex,data=incident_data[time_mri_to_icd>0])
mod_icd_lvh90 <- coxph(Surv(time_mri_to_icd,incd_icd_after_mri) ~ lvh90 + age_at_mri + sex,data=incident_data[time_mri_to_icd>0])

# ICD KM
prodlim_icd <- prodlim(Hist(time_mri_to_icd,incd_icd_after_mri)~inferred_lvh_graphical,data=incident_data[time_mri_to_icd>0])

CairoPDF(file='km_icd_ukbb.pdf',height=3,width=3.5,
         pointsize=4)
par(oma=c(3,1,1,1),mar=c(3,1,1,1))
plot(prodlim_icd,"cuminc",ylim=c(0,0.04),xlim=c(0,4),
     axis2.at=seq(0,0.04,0.01),axis2.las=2,lwd=1.5,background=F,
     axis1.at=c(0,1,2,3,4),axis1.labels=as.character(0:4),
     atrisk.times=c(0,1,2,3,4),col=c("#fc4e2a",'#1c9099'),atrisk.col='black',confint=FALSE,
     legend.x=0,legend.y=0.04,axis1.cex.axis=2.5,axis2.cex.axis=2.5,axis1.padj=0.5,
     legend.cex=2.2,legend.legend=c("CMR-derived LVH","No CMR-derived LVH"),
     atrisk.title=("                     "),atrisk.pos=-0.5,atrisk.line=c(1.2,2.5),
     atrisk.cex=1.8,atrisk.interspace=1.4,xlab='',ylab='')
mtext("Cumulative risk (%)",side=2,line=-1.2,at=0.02,cex=2.5)
mtext("Years",side=1, line=-0.6,cex=2.5)
mtext('Stratum',side=1, line=-0.3,cex=1.8,at=-0.92)
dev.off()