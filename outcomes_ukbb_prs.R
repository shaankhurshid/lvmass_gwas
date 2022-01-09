# Script to assess outcomes in UKBB using CMR truth only

# Dependencies
library(data.table)
library(survival)
library(prodlim)
library(plyr)
library(sktools)
library(Cairo)
library(ggplot2)

# Load files
censor <- fread(file='censor_202106.csv')
af <- fread(file='af_202106.csv')
mi <- fread(file='mi_icd_202106.tsv')
hf <- fread(file='hf_inclusive_202106.tsv')
vt <- fread(file='vt_202106.csv')
hcm <- fread(file='hcm_202106.csv')
dcm <- fread(file='dcm_202106.csv')
icd <- fread(file='icd_202106.csv')
pcs <- fread(file='pcs15.csv')

# Load PRS data
prs <- fread(file='prs_no_mass_processed.csv')

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
setkey(hf,sample_id); setkey(censor,sample_id); setkey(vt,sample_id); setkey(prs,sample_id); 
setkey(hcm,sample_id); setkey(dcm,sample_id); setkey(icd,sample_id)

setnames(af,'has_disease','has_af'); setnames(af,'incident_disease','incident_af'); 
setnames(af,'prevalent_disease','prevalent_af'); setnames(af,'censor_date','af_censor_date')

## Joins
af[mi,":="(has_mi = i.has_disease, prevalent_mi = i.prevalent_disease, 
           incident_mi = i.incident_disease, mi_censor_date = i.censor_date)]
af[hf,":="(has_hf = i.has_disease, prevalent_hf = i.prevalent_disease, 
           incident_hf = i.incident_disease, hf_censor_date = i.censor_date)]
af[vt,":="(has_vt = i.has_disease, prevalent_vt = i.prevalent_disease, 
           incident_vt = i.incident_disease, vt_censor_date = i.censor_date)]
af[hcm,":="(has_hcm = i.has_disease, prevalent_hcm = i.prevalent_disease, 
           incident_hcm = i.incident_disease, hcm_censor_date = i.censor_date)]
af[dcm,":="(has_dcm = i.has_disease, prevalent_dcm = i.prevalent_disease, 
           incident_dcm = i.incident_disease, dcm_censor_date = i.censor_date)]
af[icd,":="(has_icd = i.has_disease, prevalent_icd = i.prevalent_disease, 
            incident_icd = i.incident_disease, icd_censor_date = i.censor_date)]
af[censor,":="(enroll_date = i.enroll_date,death_date = i.death_date, death_age = i.death_age,
               death_censor_date = i.death_censor_date, enroll_age = i.enroll_age,
               sex = i.sex)]

## Ensure exclusions removed
withdrawals <- fread('w7089_20210809.csv') # UKBB withdrawals
all_data <- af[!(sample_id %in% withdrawals$V1)]

# Fix censoring dates
# Format dates
dates <- c('birthdate','af_censor_date','mi_censor_date','vt_censor_date',
           'hf_censor_date','enroll_date','dcm_censor_date','hcm_censor_date','icd_censor_date')
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
all_data[,':='(dcm_censor_date = as.Date(ifelse(c(has_dcm==1 | center_location=='England'),dcm_censor_date,
                                               ifelse(center_location=='Scotland',pmin(dcm_censor_date,as.Date('2021-03-31',format='%Y-%m-%d')),
                                                      pmin(dcm_censor_date,as.Date('2018-02-28',format='%Y-%m-%d')))),origin='1970-01-01'))]
all_data[,':='(hcm_censor_date = as.Date(ifelse(c(has_hcm==1 | center_location=='England'),hcm_censor_date,
                                               ifelse(center_location=='Scotland',pmin(hcm_censor_date,as.Date('2021-03-31',format='%Y-%m-%d')),
                                                      pmin(hcm_censor_date,as.Date('2018-02-28',format='%Y-%m-%d')))),origin='1970-01-01'))]
all_data[,':='(icd_censor_date = as.Date(ifelse(c(has_icd==1 | center_location=='England'),icd_censor_date,
                                                ifelse(center_location=='Scotland',pmin(icd_censor_date,as.Date('2021-03-311',format='%Y-%m-%d')),
                                                       pmin(icd_censor_date,as.Date('2018-02-28',format='%Y-%m-%d')))),origin='1970-01-01'))]
## Create specific variables for analysis
# Time
all_data[!is.na(enroll_date),":="(time_to_af = as.numeric(af_censor_date - enroll_date)/365.25,
                               time_to_mi = as.numeric(mi_censor_date - enroll_date)/365.25,
                               time_to_vt = as.numeric(vt_censor_date - enroll_date)/365.25,
                               time_to_hcm = as.numeric(hcm_censor_date - enroll_date)/365.25,
                               time_to_dcm = as.numeric(dcm_censor_date - enroll_date)/365.25,
                               time_to_hf = as.numeric(hf_censor_date - enroll_date)/365.25,
                               time_to_icd = as.numeric(icd_censor_date - enroll_date)/365.25)]

# Incident disease associations
setkey(all_data,sample_id)
incident_data <- all_data[prs,nomatch=0]

# Add PCs
setkey(pcs,sample_id)
incident_data <- incident_data[pcs,nomatch=0]

# Plot PRS distribution
ggplot() + geom_density(data=incident_data,aes(x=incident_data$prs_std),fill="#f03b20",alpha=0.55) +
  scale_x_continuous(breaks=seq(-3,3,0.5),expand=c(0.01,0),limits=c(-3,3)) +
  scale_y_continuous(breaks=seq(0,0.5,0.1),expand=c(0,0),limits=c(0,0.5)) +
  theme(panel.background=element_blank(),axis.line=element_line(color='black'),legend.position='None',
        axis.text=element_text(size=20,color='black'),plot.margin=unit(c(0.5,0.6,0.5,0.5),'cm'),
        axis.title.y = element_text(size=20,margin = margin(t = 0, r = 10, b = 0, l = 0)),
        axis.title.x = element_text(size=20),legend.text=element_text(size=20)) +
  labs(x=expression(paste("LVMI Standardized PRS")),y='Density')
ggsave(filename='prs_std_density_outcomes.pdf',height=1.8,width=2.53,
       scale=4,device='pdf')

# Create analysis sets
## AF
af_set <- incident_data[time_to_af > 0 & prevalent_af == 0]
af_set[, prs_std := (prs_std - mean(prs_std))/sd(prs_std)]
af_set[,':='(prs_decile = quantilize(prs_std,10),
             prs_ventile = quantilize(prs_std,20))]
af_set[,':='(prs_95 = ifelse(prs_ventile==20,1,0),
             prs_90 = ifelse(prs_decile==10,1,0))]

## MI
mi_set <- incident_data[time_to_mi > 0 & prevalent_mi == 0]
mi_set[, prs_std := (prs_std - mean(prs_std))/sd(prs_std)]
mi_set[,':='(prs_decile = quantilize(prs_std,10),
             prs_ventile = quantilize(prs_std,20))]
mi_set[,':='(prs_95 = ifelse(prs_ventile==20,1,0),
             prs_90 = ifelse(prs_decile==10,1,0))]

## HF
hf_set <- incident_data[time_to_hf > 0 & prevalent_hf == 0]
hf_set[, prs_std := (prs_std - mean(prs_std))/sd(prs_std)]
hf_set[,':='(prs_decile = quantilize(prs_std,10),
             prs_ventile = quantilize(prs_std,20))]
hf_set[,':='(prs_95 = ifelse(prs_ventile==20,1,0),
             prs_90 = ifelse(prs_decile==10,1,0))]

## VT
vt_set <- incident_data[time_to_vt > 0 & prevalent_vt == 0]
vt_set[, prs_std := (prs_std - mean(prs_std))/sd(prs_std)]
vt_set[,':='(prs_decile = quantilize(prs_std,10),
             prs_ventile = quantilize(prs_std,20))]
vt_set[,':='(prs_95 = ifelse(prs_ventile==20,1,0),
             prs_90 = ifelse(prs_decile==10,1,0))]

## HCM
hcm_set <- incident_data[time_to_hcm > 0 & prevalent_hcm == 0]
hcm_set[, prs_std := (prs_std - mean(prs_std))/sd(prs_std)]
hcm_set[,':='(prs_decile = quantilize(prs_std,10),
             prs_ventile = quantilize(prs_std,20))]
hcm_set[,':='(prs_95 = ifelse(prs_ventile==20,1,0),
             prs_90 = ifelse(prs_decile==10,1,0))]

## DCM
dcm_set <- incident_data[time_to_dcm > 0 & prevalent_dcm == 0]
dcm_set[, prs_std := (prs_std - mean(prs_std))/sd(prs_std)]
dcm_set[,':='(prs_decile = quantilize(prs_std,10),
             prs_ventile = quantilize(prs_std,20))]
dcm_set[,':='(prs_95 = ifelse(prs_ventile==20,1,0),
             prs_90 = ifelse(prs_decile==10,1,0))]

## ICD
icd_set <- incident_data[time_to_icd > 0 & prevalent_icd == 0]
icd_set[, prs_std := (prs_std - mean(prs_std))/sd(prs_std)]
icd_set[,':='(prs_decile = quantilize(prs_std,10),
              prs_ventile = quantilize(prs_std,20))]
icd_set[,':='(prs_95 = ifelse(prs_ventile==20,1,0),
              prs_90 = ifelse(prs_decile==10,1,0))]

### AF ANALYSES
af_time <- survSplit(Surv(time_to_af,incident_af) ~ ., data=af_set, cut=10,episode='tgroup')

mod_af <- coxph(Surv(time_to_af,incident_af) ~ prs_std + enroll_age:strata(tgroup) + sex:strata(tgroup) + pc1 + pc2 + pc3 + pc4 + pc5,data=af_time)
mod_af90 <- coxph(Surv(time_to_af,incident_af) ~ factor(prs_decile) + enroll_age:strata(tgroup) + sex:strata(tgroup) + pc1 + pc2 + pc3 + pc4 + pc5,data=af_time)
mod_af95 <- coxph(Surv(time_to_af,incident_af) ~ factor(prs_ventile) + enroll_age:strata(tgroup) + sex:strata(tgroup) + pc1 + pc2 + pc3 + pc4 + pc5,data=af_time)
mod_af_top10 <- coxph(Surv(time_to_af,incident_af) ~ prs_90 + enroll_age:strata(tgroup) + sex:strata(tgroup) + pc1 + pc2 + pc3 + pc4 + pc5,data=af_time)
mod_af_top5 <- coxph(Surv(time_to_af,incident_af) ~ prs_95 + enroll_age:strata(tgroup) + sex:strata(tgroup) + pc1 + pc2 + pc3 + pc4 + pc5,data=af_time)

# Spline
af_spline <- coxph(Surv(time_to_af,incident_af) ~ pspline(prs_std,df=0) + enroll_age + sex + pc1 + pc2 + pc3 + pc4 + pc5,data=af_set)

# Plot spline
resid_af_spline <- termplot(af_spline,term=1,se=TRUE,plot=F)
value_term <- resid_af_spline$prs_std
center <- with(value_term, y[x==min(abs((x-0)))])
y <- value_term$y + outer(value_term$se,c(0,-1.96,1.96),'*')

CairoPDF(file='af_spline.pdf',height=3,width=3.5,
         pointsize=4)
par(oma=c(3,4,1,1),mar=c(3,4,1,1),xpd=FALSE)
matplot(value_term$x,exp(y-center),log='y',type='l',xaxt='n',xlab='',ylab='',bty='n',yaxt='n',
        xlim=c(-2.4,2.4),ylim=c(0.5,1.5),col=c('black','#d73027','#1a9850'))
axis(1,at=seq(-2.5,2.5,0.5),cex.axis=1.4)
axis(2,at=seq(0.5,1.5,0.25),las=2,cex.axis=1.4,pos=-2.6)
legend(-4,0.05,c('estimate','upper 95% CI','lower 95% CI'),bty='n',lty=1,
       col=c('black','#d73027','#1a9850'),cex=1.4)
mtext("Standardized PRS",1,line=3,cex=1.5)
mtext("Relative hazard for AF compared to PRS of zero",2,line=5.5,cex=1.5)
dev.off()

# Graphical variables
af_set[,prs_graphical90 := factor(ifelse(prs_ventile==20,"Above 95th percentile","Below 95th percentile"),levels=c('Above 95th percentile','Below 95th percentile'))]

# AF KM
prodlim_af <- prodlim(Hist(time_to_af,incident_af)~prs_graphical90,data=af_set)

CairoPDF(file='km_af_prs.pdf',height=3,width=3.5,
         pointsize=4)
par(oma=c(3,1,1,1),mar=c(3,1,1,1))
plot(prodlim_af,"cuminc",ylim=c(0,0.10),xlim=c(0,12),
     axis2.at=seq(0,0.10,0.01),axis2.las=2,lwd=1.5,background=F,
     axis1.at=seq(0,12,2),axis1.labels=as.character(seq(0,12,2)),
     atrisk.times=c(0,4,8,12),col=c("#fc4e2a",'#1c9099'),atrisk.col='black',confint=FALSE,
     legend.x=0,legend.y=0.10,axis1.cex.axis=2.5,axis2.cex.axis=2.5,axis1.padj=0.5,
     legend.cex=2.2,legend.legend=c("Above 95th percentile","Below 95th percentile"),
     atrisk.title='',atrisk.pos=-1,atrisk.line=c(1.2,2.5),
     atrisk.cex=1.8,atrisk.interspace=1.4,xlab='',ylab='')
mtext("Cumulative risk of AF (%)",side=2,line=-1.2,at=0.05,cex=2.5)
mtext("Years",side=1, line=-0.6,cex=2.5)
mtext('Stratum',side=1, line=-0.3,cex=1.8,at=-2.5)
dev.off()

# HR plot for AF
pdf('hr_plot_af.pdf',pointsize=6,
    height=4,width=8)
par(oma=c(4,5,2,2))
par(mar=c(4,5,2,2))
col <- '#88419d'

est <- c(1,summary(mod_af95)$conf.int[1:19,1])
lower <- c(1,summary(mod_af95)$conf.int[1:19,3])
upper <- c(1,summary(mod_af95)$conf.int[1:19,4])

# Plot
plot(x=1:20,y=est,pch=c(18,rep(19,19)),xlab='',ylab='',xaxt='n',yaxt='n',frame=F,
     xlim=c(1,20),ylim=c(0,1.5),cex=c(2.5,rep(1.5,19)),col=c('darkgray',rep(col,19)))

# CIs
for (i in 1:20){
  segments(i,lower[i],i,upper[i],col=col,lwd=1.7)
}

# Axes
axis(1,at=1:20,cex.axis=1.8,labels=as.character(1:20))
axis(2,at=seq(0,1.5,0.3),cex.axis=1.8,las=1)
mtext(side=1,'PRS demi-decile',cex=2.2,line=5)
mtext(side=2,'Hazard ratio for AF',cex=2.2,line=5.7)

segments(-1,1,20,1,lty=5)

par(xpd=NA)

dev.off()

### MI ANALYSES
mi_set[,sex_numeric := ifelse(sex=='Male',1,0)]
mi_time <- survSplit(Surv(time_to_mi,incident_mi) ~ ., data=mi_set, cut=9,episode='tgroup')

mod_mi <- coxph(Surv(time_to_mi,incident_mi) ~ prs_std + enroll_age:strata(tgroup) + sex_numeric:strata(tgroup) + pc1 + pc2 + pc3 + pc4:strata(tgroup) + pc5,data=mi_time)
mod_mi90 <- coxph(Surv(time_to_mi,incident_mi) ~ factor(prs_decile) + enroll_age:strata(tgroup) + sex_numeric:strata(tgroup) + pc1 + pc2 + pc3 + pc4:strata(tgroup) + pc5,data=mi_time)
mod_mi95 <- coxph(Surv(time_to_mi,incident_mi) ~ factor(prs_ventile) + enroll_age:strata(tgroup) + sex_numeric:strata(tgroup) + pc1 + pc2 + pc3 + pc4:strata(tgroup) + pc5,data=mi_time)
mod_mi_top10 <- coxph(Surv(time_to_mi,incident_mi) ~ prs_90 + enroll_age:strata(tgroup) + sex_numeric:strata(tgroup) + pc1 + pc2 + pc3 + pc4:strata(tgroup) + pc5,data=mi_time)
mod_mi_top5 <- coxph(Surv(time_to_mi,incident_mi) ~ prs_95 + enroll_age:strata(tgroup) + sex_numeric:strata(tgroup) + pc1 + pc2 + pc3 + pc4:strata(tgroup) + pc5,data=mi_time)

# Spline
mi_spline <- coxph(Surv(time_to_mi,incident_mi) ~ pspline(prs_std,df=0) + enroll_age + sex + pc1 + pc2 + pc3 + pc4 + pc5,data=mi_set)

# Plot spline
resid_mi_spline <- termplot(mi_spline,term=1,se=TRUE,plot=F)
value_term <- resid_mi_spline$prs_std
center <- with(value_term, y[x==min(abs((x-0)))])
y <- value_term$y + outer(value_term$se,c(0,-1.96,1.96),'*')

CairoPDF(file='mi_spline.pdf',height=3,width=3.5,
         pointsize=4)
par(oma=c(3,4,1,1),mar=c(3,4,1,1),xpd=FALSE)
matplot(value_term$x,exp(y-center),log='y',type='l',xaxt='n',xlab='',ylab='',bty='n',yaxt='n',
        xlim=c(-2.4,2.4),ylim=c(0.5,1.5),col=c('black','#d73027','#1a9850'))
axis(1,at=seq(-2.5,2.5,0.5),cex.axis=1.4)
axis(2,at=seq(0.5,1.5,0.25),las=2,cex.axis=1.4,pos=-2.6)
legend(-4,0.05,c('estimate','upper 95% CI','lower 95% CI'),bty='n',lty=1,
       col=c('black','#d73027','#1a9850'),cex=1.4)
mtext("Standardized PRS",1,line=3,cex=1.5)
mtext("Relative hazard for MI compared to PRS of zero",2,line=5.5,cex=1.5)
dev.off()

# HR plot for MI
pdf('hr_plot_mi.pdf',pointsize=6,
    height=4,width=8)
par(oma=c(4,5,2,2))
par(mar=c(4,5,2,2))
col <- '#78c679'

est <- c(1,summary(mod_mi95)$conf.int[1:19,1])
lower <- c(1,summary(mod_mi95)$conf.int[1:19,3])
upper <- c(1,summary(mod_mi95)$conf.int[1:19,4])

# Plot
plot(x=1:20,y=est,pch=c(18,rep(19,19)),xlab='',ylab='',xaxt='n',yaxt='n',frame=F,
     xlim=c(1,20),ylim=c(0,1.5),cex=c(2.5,rep(1.5,19)),col=c('darkgray',rep(col,19)))

# CIs
for (i in 1:20){
  segments(i,lower[i],i,upper[i],col=col,lwd=1.7)
}

# Axes
axis(1,at=1:20,cex.axis=1.8,labels=as.character(1:20))
axis(2,at=seq(0,1.5,0.3),cex.axis=1.8,las=1)
mtext(side=1,'PRS demi-decile',cex=2.2,line=5)
mtext(side=2,'Hazard ratio for MI',cex=2.2,line=5.7)

segments(-1,1,20,1,lty=5)

par(xpd=NA)

dev.off()

# MI KM
mi_set[,prs_graphical90 := factor(ifelse(prs_ventile==20,"Above 95th percentile","Below 95th percentile"),levels=c('Above 95th percentile','Below 95th percentile'))]
prodlim_mi <- prodlim(Hist(time_to_mi,incident_mi)~prs_graphical90,data=mi_set)

CairoPDF(file='km_mi_prs.pdf',height=3,width=3.5,
         pointsize=4)
par(oma=c(3,1,1,1),mar=c(3,1,1,1))
plot(prodlim_mi,"cuminc",ylim=c(0,0.05),xlim=c(0,12),
     axis2.at=seq(0,0.05,0.01),axis2.las=2,lwd=1.5,background=F,
     axis1.at=seq(0,12,2),axis1.labels=as.character(seq(0,12,2)),
     atrisk.times=c(0,4,8,12),col=c("#fc4e2a",'#1c9099'),atrisk.col='black',confint=FALSE,
     legend.x=0,legend.y=0.05,axis1.cex.axis=2.5,axis2.cex.axis=2.5,axis1.padj=0.5,
     legend.cex=2.2,legend.legend=c("Above 95th percentile","Below 95th percentile"),
     atrisk.title='',atrisk.pos=-1,atrisk.line=c(1.2,2.5),
     atrisk.cex=1.8,atrisk.interspace=1.4,xlab='',ylab='')
mtext("Cumulative risk of MI (%)",side=2,line=-1.2,at=0.025,cex=2.5)
mtext("Years",side=1, line=-0.6,cex=2.5)
mtext('Stratum',side=1, line=-0.3,cex=1.8,at=-2.5)
dev.off()

### HF ANALYSES
hf_set[,sex_numeric := ifelse(sex=='Male',1,0)]
hf_time <- survSplit(Surv(time_to_hf,incident_hf) ~ ., data=hf_set, cut=9,episode='tgroup')

mod_hf <- coxph(Surv(time_to_hf,incident_hf) ~ prs_std + enroll_age + sex_numeric:strata(tgroup) + pc1 + pc2 + pc3 + pc4 + pc5,data=hf_time)
mod_hf90 <- coxph(Surv(time_to_hf,incident_hf) ~ factor(prs_decile) + enroll_age + sex_numeric:strata(tgroup) + pc1 + pc2 + pc3 + pc4 + pc5,data=hf_time)
mod_hf95 <- coxph(Surv(time_to_hf,incident_hf) ~ factor(prs_ventile) + enroll_age + sex_numeric:strata(tgroup) + pc1 + pc2 + pc3 + pc4 + pc5,data=hf_time)
mod_hf_top10 <- coxph(Surv(time_to_hf,incident_hf) ~ prs_90 + enroll_age + sex_numeric:strata(tgroup) + pc1 + pc2 + pc3 + pc4 + pc5,data=hf_time)
mod_hf_top5 <- coxph(Surv(time_to_hf,incident_hf) ~ prs_95 + enroll_age + sex_numeric:strata(tgroup) + pc1 + pc2 + pc3 + pc4 + pc5,data=hf_time)

# Spline
hf_spline <- coxph(Surv(time_to_hf,incident_hf) ~ pspline(prs_std,df=0) + enroll_age + sex + pc1 + pc2 + pc3 + pc4 + pc5,data=hf_set)

# Plot spline
resid_hf_spline <- termplot(hf_spline,term=1,se=TRUE,plot=F)
value_term <- resid_hf_spline$prs_std
center <- with(value_term, y[x==min(abs((x-0)))])
y <- value_term$y + outer(value_term$se,c(0,-1.96,1.96),'*')

CairoPDF(file='hf_spline.pdf',height=3,width=3.5,
         pointsize=4)
par(oma=c(3,4,1,1),mar=c(3,4,1,1),xpd=FALSE)
matplot(value_term$x,exp(y-center),log='y',type='l',xaxt='n',xlab='',ylab='',bty='n',yaxt='n',
        xlim=c(-2.4,2.4),ylim=c(0.5,1.5),col=c('black','#d73027','#1a9850'))
axis(1,at=seq(-2.5,2.5,0.5),cex.axis=1.4)
axis(2,at=seq(0.5,1.5,0.25),las=2,cex.axis=1.4,pos=-2.6)
legend(-4,0.05,c('estimate','upper 95% CI','lower 95% CI'),bty='n',lty=1,
       col=c('black','#d73027','#1a9850'),cex=1.4)
mtext("Standardized PRS",1,line=3,cex=1.5)
mtext("Relative hazard for HF compared to PRS of zero",2,line=5.5,cex=1.5)
dev.off()

# HR plot for HF
pdf('hr_plot_hf.pdf',pointsize=6,
    height=4,width=8)
par(oma=c(4,5,2,2))
par(mar=c(4,5,2,2))
col <- '#f03b20'

est <- c(1,summary(mod_hf95)$conf.int[1:19,1])
lower <- c(1,summary(mod_hf95)$conf.int[1:19,3])
upper <- c(1,summary(mod_hf95)$conf.int[1:19,4])

# Plot
plot(x=1:20,y=est,pch=c(18,rep(19,19)),xlab='',ylab='',xaxt='n',yaxt='n',frame=F,
     xlim=c(1,20),ylim=c(0,1.5),cex=c(2.5,rep(1.5,19)),col=c('darkgray',rep(col,19)))

# CIs
for (i in 1:20){
  segments(i,lower[i],i,upper[i],col=col,lwd=1.7)
}

# Axes
axis(1,at=1:20,cex.axis=1.8,labels=as.character(1:20))
axis(2,at=seq(0,1.5,0.3),cex.axis=1.8,las=1)
mtext(side=1,'PRS demi-decile',cex=2.2,line=5)
mtext(side=2,'Hazard ratio for HF',cex=2.2,line=5.7)

segments(-1,1,20,1,lty=5)

par(xpd=NA)

dev.off()

# HF KM
hf_set[,prs_graphical90 := factor(ifelse(prs_ventile==20,"Above 95th percentile","Below 95th percentile"),levels=c('Above 95th percentile','Below 95th percentile'))]
prodlim_hf <- prodlim(Hist(time_to_hf,incident_hf)~prs_graphical90,data=hf_set)

CairoPDF(file='km_hf_prs.pdf',height=3,width=3.5,
         pointsize=4)
par(oma=c(3,1,1,1),mar=c(3,1,1,1))
plot(prodlim_hf,"cuminc",ylim=c(0,0.05),xlim=c(0,12),
     axis2.at=seq(0,0.05,0.01),axis2.las=2,lwd=1.5,background=F,
     axis1.at=seq(0,12,2),axis1.labels=as.character(seq(0,12,2)),
     atrisk.times=c(0,4,8,12),col=c("#fc4e2a",'#1c9099'),atrisk.col='black',confint=FALSE,
     legend.x=0,legend.y=0.05,axis1.cex.axis=2.5,axis2.cex.axis=2.5,axis1.padj=0.5,
     legend.cex=2.2,legend.legend=c("Above 95th percentile","Below 95th percentile"),
     atrisk.title='',atrisk.pos=-1,atrisk.line=c(1.2,2.5),
     atrisk.cex=1.8,atrisk.interspace=1.4,xlab='',ylab='')
mtext("Cumulative risk of HF (%)",side=2,line=-1.2,at=0.025,cex=2.5)
mtext("Years",side=1, line=-0.6,cex=2.5)
mtext('Stratum',side=1, line=-0.3,cex=1.8,at=-2.5)
dev.off()

### VT ANALYSES
mod_vt <- coxph(Surv(time_to_vt,incident_vt) ~ prs_std + enroll_age + sex + pc1 + pc2 + pc3 + pc4 + pc5,data=vt_set)
mod_vt90 <- coxph(Surv(time_to_vt,incident_vt) ~ factor(prs_decile) + enroll_age + sex + pc1 + pc2 + pc3 + pc4 + pc5,data=vt_set)
mod_vt95 <- coxph(Surv(time_to_vt,incident_vt) ~ factor(prs_ventile) + enroll_age + sex + pc1 + pc2 + pc3 + pc4 + pc5,data=vt_set)
mod_vt_top10 <- coxph(Surv(time_to_vt,incident_vt) ~ prs_90 + enroll_age + sex + pc1 + pc2 + pc3 + pc4 + pc5,data=vt_set)
mod_vt_top5 <- coxph(Surv(time_to_vt,incident_vt) ~ prs_95 + enroll_age + sex + pc1 + pc2 + pc3 + pc4 + pc5,data=vt_set)

# Spline
vt_spline <- coxph(Surv(time_to_vt,incident_vt) ~ pspline(prs_std,df=0) + enroll_age + sex + pc1 + pc2 + pc3 + pc4 + pc5,data=vt_set)

# Plot spline
resid_vt_spline <- termplot(vt_spline,term=1,se=TRUE,plot=F)
value_term <- resid_vt_spline$prs_std
center <- with(value_term, y[which.min(abs((x-0)))])
y <- value_term$y + outer(value_term$se,c(0,-1.96,1.96),'*')

CairoPDF(file='vt_spline.pdf',height=3,width=3.5,
         pointsize=4)
par(oma=c(3,4,1,1),mar=c(3,4,1,1),xpd=FALSE)
matplot(value_term$x,exp(y-center),log='y',type='l',xaxt='n',xlab='',ylab='',bty='n',yaxt='n',
        xlim=c(-2.4,2.4),ylim=c(0.5,1.5),col=c('black','#d73027','#1a9850'))
axis(1,at=seq(-2.5,2.5,0.5),cex.axis=1.4)
axis(2,at=seq(0.5,1.5,0.25),las=2,cex.axis=1.4,pos=-2.6)
legend(-4,0.05,c('estimate','upper 95% CI','lower 95% CI'),bty='n',lty=1,
       col=c('black','#d73027','#1a9850'),cex=1.4)
mtext("Standardized PRS",1,line=3,cex=1.5)
mtext("Relative hazard for VT compared to PRS of zero",2,line=5.5,cex=1.5)
dev.off()

# HR plot for VT
pdf('hr_plot_vt.pdf',pointsize=6,
    height=4,width=8)
par(oma=c(4,5,2,2))
par(mar=c(4,5,2,2))
col <- '#2b8cbe'

est <- c(1,summary(mod_vt95)$conf.int[1:19,1])
lower <- c(1,summary(mod_vt95)$conf.int[1:19,3])
upper <- c(1,summary(mod_vt95)$conf.int[1:19,4])

# Plot
plot(x=1:20,y=est,pch=c(18,rep(19,19)),xlab='',ylab='',xaxt='n',yaxt='n',frame=F,
     xlim=c(1,20),ylim=c(0,1.5),cex=c(2.5,rep(1.5,19)),col=c('darkgray',rep(col,19)))

# CIs
for (i in 1:20){
  segments(i,lower[i],i,upper[i],col=col,lwd=1.7)
}

# Axes
axis(1,at=1:20,cex.axis=1.8,labels=as.character(1:20))
axis(2,at=seq(0,1.5,0.3),cex.axis=1.8,las=1)
mtext(side=1,'PRS demi-decile',cex=2.2,line=5)
mtext(side=2,'Hazard ratio for VT',cex=2.2,line=5.7)

segments(-1,1,20,1,lty=5)

par(xpd=NA)

dev.off()

# VT KM
vt_set[,prs_graphical90 := factor(ifelse(prs_ventile==20,"Above 95th percentile","Below 95th percentile"),levels=c('Above 95th percentile','Below 95th percentile'))]
prodlim_vt <- prodlim(Hist(time_to_vt,incident_vt)~prs_graphical90,data=vt_set)

CairoPDF(file='km_vt_prs.pdf',height=3,width=3.5,
         pointsize=4)
par(oma=c(3,1,1,1),mar=c(3,1,1,1))
plot(prodlim_vt,"cuminc",ylim=c(0,0.02),xlim=c(0,12),
     axis2.at=seq(0,0.02,0.005),axis2.las=2,lwd=1.5,background=F,
     axis1.at=seq(0,12,2),axis1.labels=as.character(seq(0,12,2)),
     atrisk.times=c(0,4,8,12),col=c("#fc4e2a",'#1c9099'),atrisk.col='black',confint=FALSE,
     legend.x=0,legend.y=0.02,axis1.cex.axis=2.5,axis2.cex.axis=2.5,axis1.padj=0.5,
     legend.cex=2.2,legend.legend=c("Above 95th percentile","Below 95th percentile"),
     atrisk.title='',atrisk.pos=-1,atrisk.line=c(1.2,2.5),
     atrisk.cex=1.8,atrisk.interspace=1.4,xlab='',ylab='')
mtext("Cumulative risk of VT (%)",side=2,line=-1.2,at=0.01,cex=2.5)
mtext("Years",side=1, line=-0.6,cex=2.5)
mtext('Stratum',side=1, line=-0.3,cex=1.8,at=-2.5)
dev.off()

### HCM ANALYSES
mod_hcm <- coxph(Surv(time_to_hcm,incident_hcm) ~ prs_std + enroll_age + sex + pc1 + pc2 + pc3 + pc4 + pc5,data=hcm_set)
mod_hcm90 <- coxph(Surv(time_to_hcm,incident_hcm) ~ factor(prs_decile) + enroll_age + sex + pc1 + pc2 + pc3 + pc4 + pc5,data=hcm_set)
mod_hcm95 <- coxph(Surv(time_to_hcm,incident_hcm) ~ factor(prs_ventile) + enroll_age + sex + pc1 + pc2 + pc3 + pc4 + pc5,data=hcm_set)
mod_hcm_top10 <- coxph(Surv(time_to_hcm,incident_hcm) ~ prs_90 + enroll_age + sex + pc1 + pc2 + pc3 + pc4 + pc5,data=hcm_set)
mod_hcm_top5 <- coxph(Surv(time_to_hcm,incident_hcm) ~ prs_95 + enroll_age + sex + pc1 + pc2 + pc3 + pc4 + pc5,data=hcm_set)

# Spline
hcm_spline <- coxph(Surv(time_to_hcm,incident_hcm) ~ pspline(prs_std,df=0,caic=T) + enroll_age + sex + pc1 + pc2 + pc3 + pc4 + pc5,data=hcm_set)

# Plot spline
resid_hcm_spline <- termplot(hcm_spline,term=1,se=TRUE,plot=F)
value_term <- resid_hcm_spline$prs_std
center <- with(value_term, y[which.min(abs((x-0)))])
y <- value_term$y + outer(value_term$se,c(0,-1.96,1.96),'*')

CairoPDF(file='hcm_spline.pdf',height=3,width=3.5,
         pointsize=4)
par(oma=c(3,4,1,1),mar=c(3,4,1,1),xpd=FALSE)
matplot(value_term$x,exp(y-center),log='y',type='l',xaxt='n',xlab='',ylab='',bty='n',yaxt='n',
        xlim=c(-3.8,3.8),ylim=c(0.25,4),col=c('black','#d73027','#1a9850'))
axis(1,at=seq(-4,4,1),cex.axis=1.6)
axis(2,las=2,cex.axis=1.6,pos=-4.3,at=c(0.25,0.5,1,2,4),labels=c('0.25','0.5','1','2','4'))
legend(1.5,0.4,c('estimate','upper 95% CI','lower 95% CI'),bty='n',lty=1,
       col=c('black','#d73027','#1a9850'),cex=1.4)
mtext("Standardized PRS",1,line=3,cex=1.8)
mtext("Relative hazard for HCM compared to average PRS",2,line=5.5,cex=1.8)
dev.off()

# Graphical variables
hcm_set[,prs_graphical90 := factor(ifelse(prs_ventile==20,"Above 95th percentile","Below 95th percentile"),levels=c('Above 95th percentile','Below 95th percentile'))]

# HCM KM
prodlim_hcm <- prodlim(Hist(time_to_hcm,incident_hcm)~prs_graphical90,data=hcm_set)

CairoPDF(file='km_hcm_prs.pdf',height=3,width=3.5,
         pointsize=4)
par(oma=c(3,1,1,1),mar=c(3,1,1,1))
plot(prodlim_hcm,"cuminc",ylim=c(0,0.10),xlim=c(0,12),
     axis2.at=seq(0,0.10,0.01),axis2.las=2,lwd=1.5,background=F,
     axis1.at=seq(0,12,2),axis1.labels=as.character(seq(0,12,2)),
     atrisk.times=c(0,4,8,12),col=c("#fc4e2a",'#1c9099'),atrisk.col='black',confint=FALSE,
     legend.x=0,legend.y=0.10,axis1.cex.axis=2.5,axis2.cex.axis=2.5,axis1.padj=0.5,
     legend.cex=2.2,legend.legend=c("Above 95th percentile","Below 95th percentile"),
     atrisk.title='',atrisk.pos=-1,atrisk.line=c(1.2,2.5),
     atrisk.cex=1.8,atrisk.interspace=1.4,xlab='',ylab='')
mtext("Cumulative risk of hcm (%)",side=2,line=-1.2,at=0.05,cex=2.5)
mtext("Years",side=1, line=-0.6,cex=2.5)
mtext('Stratum',side=1, line=-0.3,cex=1.8,at=-2.5)
dev.off()

### DCM ANALYSES
mod_dcm <- coxph(Surv(time_to_dcm,incident_dcm) ~ prs_std + enroll_age + sex + pc1 + pc2 + pc3 + pc4 + pc5,data=dcm_set)
mod_dcm90 <- coxph(Surv(time_to_dcm,incident_dcm) ~ factor(prs_decile) + enroll_age + sex + pc1 + pc2 + pc3 + pc4 + pc5,data=dcm_set)
mod_dcm95 <- coxph(Surv(time_to_dcm,incident_dcm) ~ factor(prs_ventile) + enroll_age + sex + pc1 + pc2 + pc3 + pc4 + pc5,data=dcm_set)
mod_dcm_top10 <- coxph(Surv(time_to_dcm,incident_dcm) ~ prs_90 + enroll_age + sex + pc1 + pc2 + pc3 + pc4 + pc5,data=dcm_set)
mod_dcm_top5 <- coxph(Surv(time_to_dcm,incident_dcm) ~ prs_95 + enroll_age + sex + pc1 + pc2 + pc3 + pc4 + pc5,data=dcm_set)

# Spline
dcm_spline <- coxph(Surv(time_to_dcm,incident_dcm) ~ pspline(prs_std,df=0) + enroll_age + sex + pc1 + pc2 + pc3 + pc4 + pc5,data=dcm_set)

# Plot spline
resid_dcm_spline <- termplot(dcm_spline,term=1,se=TRUE,plot=F)
value_term <- resid_dcm_spline$prs_std
center <- with(value_term, y[x==min(abs((x-0)))])
y <- value_term$y + outer(value_term$se,c(0,-1.96,1.96),'*')

CairoPDF(file='dcm_spline.pdf',height=3,width=3.5,
         pointsize=4)
par(oma=c(3,4,1,1),mar=c(3,4,1,1),xpd=FALSE)
matplot(value_term$x,exp(y-center),log='y',type='l',xaxt='n',xlab='',ylab='',bty='n',yaxt='n',
        xlim=c(-2.4,2.4),ylim=c(0.5,1.5),col=c('black','#d73027','#1a9850'))
axis(1,at=seq(-2.5,2.5,0.5),cex.axis=1.4)
axis(2,at=seq(0.5,1.5,0.25),las=2,cex.axis=1.4,pos=-2.6)
legend(-4,0.05,c('estimate','upper 95% CI','lower 95% CI'),bty='n',lty=1,
       col=c('black','#d73027','#1a9850'),cex=1.4)
mtext("Standardized PRS",1,line=3,cex=1.5)
mtext("Relative hazard for dcm compared to PRS of zero",2,line=5.5,cex=1.5)
dev.off()

# Graphical variables
dcm_set[,prs_graphical90 := factor(ifelse(prs_ventile==20,"Above 95th percentile","Below 95th percentile"),levels=c('Above 95th percentile','Below 95th percentile'))]

# DCM KM
prodlim_dcm <- prodlim(Hist(time_to_dcm,incident_dcm)~prs_graphical90,data=dcm_set)

CairoPDF(file='km_dcm_prs.pdf',height=3,width=3.5,
         pointsize=4)
par(oma=c(3,1,1,1),mar=c(3,1,1,1))
plot(prodlim_dcm,"cuminc",ylim=c(0,0.10),xlim=c(0,12),
     axis2.at=seq(0,0.10,0.01),axis2.las=2,lwd=1.5,background=F,
     axis1.at=seq(0,12,2),axis1.labels=as.character(seq(0,12,2)),
     atrisk.times=c(0,4,8,12),col=c("#fc4e2a",'#1c9099'),atrisk.col='black',confint=FALSE,
     legend.x=0,legend.y=0.10,axis1.cex.axis=2.5,axis2.cex.axis=2.5,axis1.padj=0.5,
     legend.cex=2.2,legend.legend=c("Above 95th percentile","Below 95th percentile"),
     atrisk.title='',atrisk.pos=-1,atrisk.line=c(1.2,2.5),
     atrisk.cex=1.8,atrisk.interspace=1.4,xlab='',ylab='')
mtext("Cumulative risk of dcm (%)",side=2,line=-1.2,at=0.05,cex=2.5)
mtext("Years",side=1, line=-0.6,cex=2.5)
mtext('Stratum',side=1, line=-0.3,cex=1.8,at=-2.5)
dev.off()

### ICD ANALYSES
mod_icd <- coxph(Surv(time_to_icd,incident_icd) ~ prs_std + enroll_age + sex + pc1 + pc2 + pc3 + pc4 + pc5,data=icd_set)
mod_icd90 <- coxph(Surv(time_to_icd,incident_icd) ~ factor(prs_decile) + enroll_age + sex + pc1 + pc2 + pc3 + pc4 + pc5,data=icd_set)
mod_icd95 <- coxph(Surv(time_to_icd,incident_icd) ~ factor(prs_ventile) + enroll_age + sex + pc1 + pc2 + pc3 + pc4 + pc5,data=icd_set)
mod_icd_top10 <- coxph(Surv(time_to_icd,incident_icd) ~ prs_90 + enroll_age + sex + pc1 + pc2 + pc3 + pc4 + pc5,data=icd_set)
mod_icd_top5 <- coxph(Surv(time_to_icd,incident_icd) ~ prs_95 + enroll_age + sex + pc1 + pc2 + pc3 + pc4 + pc5,data=icd_set)

# Spline
icd_spline <- coxph(Surv(time_to_icd,incident_icd) ~ pspline(prs_std,df=0) + enroll_age + sex + pc1 + pc2 + pc3 + pc4 + pc5,data=icd_set)

# Plot spline
resid_icd_spline <- termplot(icd_spline,term=1,se=TRUE,plot=F)
value_term <- resid_icd_spline$prs_std
center <- with(value_term, y[x==min(abs((x-0)))])
y <- value_term$y + outer(value_term$se,c(0,-1.96,1.96),'*')

CairoPDF(file='icd_spline.pdf',height=3,width=3.5,
         pointsize=4)
par(oma=c(3,4,1,1),mar=c(3,4,1,1),xpd=FALSE)
matplot(value_term$x,exp(y-center),log='y',type='l',xaxt='n',xlab='',ylab='',bty='n',yaxt='n',
        xlim=c(-2.4,2.4),ylim=c(0.5,1.5),col=c('black','#d73027','#1a9850'))
axis(1,at=seq(-2.5,2.5,0.5),cex.axis=1.4)
axis(2,at=seq(0.5,1.5,0.25),las=2,cex.axis=1.4,pos=-2.6)
legend(-4,0.05,c('estimate','upper 95% CI','lower 95% CI'),bty='n',lty=1,
       col=c('black','#d73027','#1a9850'),cex=1.4)
mtext("Standardized PRS",1,line=3,cex=1.5)
mtext("Relative hazard for icd compared to PRS of zero",2,line=5.5,cex=1.5)
dev.off()

# Graphical variables
icd_set[,prs_graphical90 := factor(ifelse(prs_ventile==20,"Above 95th percentile","Below 95th percentile"),levels=c('Above 95th percentile','Below 95th percentile'))]

# icd KM
prodlim_icd <- prodlim(Hist(time_to_icd,incident_icd)~prs_graphical90,data=icd_set)

CairoPDF(file='km_icd_prs.pdf',height=3,width=3.5,
         pointsize=4)
par(oma=c(3,1,1,1),mar=c(3,1,1,1))
plot(prodlim_icd,"cuminc",ylim=c(0,0.10),xlim=c(0,12),
     axis2.at=seq(0,0.10,0.01),axis2.las=2,lwd=1.5,background=F,
     axis1.at=seq(0,12,2),axis1.labels=as.character(seq(0,12,2)),
     atrisk.times=c(0,4,8,12),col=c("#fc4e2a",'#1c9099'),atrisk.col='black',confint=FALSE,
     legend.x=0,legend.y=0.10,axis1.cex.axis=2.5,axis2.cex.axis=2.5,axis1.padj=0.5,
     legend.cex=2.2,legend.legend=c("Above 95th percentile","Below 95th percentile"),
     atrisk.title='',atrisk.pos=-1,atrisk.line=c(1.2,2.5),
     atrisk.cex=1.8,atrisk.interspace=1.4,xlab='',ylab='')
mtext("Cumulative risk of icd (%)",side=2,line=-1.2,at=0.05,cex=2.5)
mtext("Years",side=1, line=-0.6,cex=2.5)
mtext('Stratum',side=1, line=-0.3,cex=1.8,at=-2.5)
dev.off()