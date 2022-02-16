# Script to perform time to event analyses using LVM PRS in MGB

# Dependencies
library(data.table)
library(survival)
library(prodlim)
library(plyr)
library(Cairo)
library(sktools)

# Load phenotype data
prs <- fread(file='lvm_prs_phenos.csv')

## Step 1: Create cox variables
# Define prevalent disease
prs[,':='(prev_mi = ifelse(c(!is.na(mi) & mi <= collect_age),1,0),
          prev_af = ifelse(c(!is.na(af) & af <= collect_age),1,0),
          prev_vt = ifelse(c(!is.na(vt) & vt <= collect_age),1,0),
          prev_hf = ifelse(c(!is.na(hf) & hf <= collect_age),1,0),
          prev_dcm = ifelse(c(!is.na(dcm) & dcm <= collect_age),1,0),
          prev_hcm = ifelse(c(!is.na(hcm) & hcm <= collect_age),1,0),
          prev_icd = ifelse(c(!is.na(icd) & icd <= collect_age),1,0))]

# Define incident disease
prs[,':='(incd_mi = ifelse(c(!is.na(mi) & mi > collect_age),1,0),
          incd_af = ifelse(c(!is.na(af) & af > collect_age),1,0),
          incd_vt = ifelse(c(!is.na(vt) & vt > collect_age),1,0),
          incd_hf = ifelse(c(!is.na(hf) & hf > collect_age),1,0),
          incd_dcm = ifelse(c(!is.na(dcm) & dcm > collect_age),1,0),
          incd_hcm = ifelse(c(!is.na(hcm) & hcm > collect_age),1,0),
          incd_icd = ifelse(c(!is.na(icd) & icd > collect_age),1,0))]

# Define time vars
prs[,':='(time_to_mi = ifelse(c(prev_mi == 1 | incd_mi == 1),mi-collect_age,censor_age-collect_age),
          time_to_af = ifelse(c(prev_af == 1 | incd_af == 1),af-collect_age,censor_age-collect_age),
          time_to_vt = ifelse(c(prev_vt == 1 | incd_vt == 1),vt-collect_age,censor_age-collect_age),
          time_to_hf = ifelse(c(prev_hf == 1 | incd_hf == 1),hf-collect_age,censor_age-collect_age),
          time_to_dcm = ifelse(c(prev_dcm == 1 | incd_dcm == 1),dcm-collect_age,censor_age-collect_age),
          time_to_hcm = ifelse(c(prev_hcm == 1 | incd_hcm == 1),hcm-collect_age,censor_age-collect_age),
          time_to_icd = ifelse(c(prev_icd == 1 | incd_icd == 1),icd-collect_age,censor_age-collect_age))]

# Create analysis sets
## AF
af_set <- prs[time_to_af > 0 & prev_af == 0]
af_set[, lvm_prs_std := (lvm_prs_std - mean(lvm_prs_std))/sd(lvm_prs_std)]
af_set[,':='(prs_decile = quantilize(lvm_prs_std,10),
             prs_ventile = quantilize(lvm_prs_std,20))]
af_set[,':='(prs_95 = ifelse(prs_ventile==20,1,0),
             prs_90 = ifelse(prs_decile==10,1,0))]

## MI
mi_set <- prs[time_to_mi > 0 & prev_mi == 0]
mi_set[, lvm_prs_std := (lvm_prs_std - mean(lvm_prs_std))/sd(lvm_prs_std)]
mi_set[,':='(prs_decile = quantilize(lvm_prs_std,10),
             prs_ventile = quantilize(lvm_prs_std,20))]
mi_set[,':='(prs_95 = ifelse(prs_ventile==20,1,0),
             prs_90 = ifelse(prs_decile==10,1,0))]

## HF
hf_set <- prs[time_to_hf > 0 & prev_hf == 0]
hf_set[, lvm_prs_std := (lvm_prs_std - mean(lvm_prs_std))/sd(lvm_prs_std)]
hf_set[,':='(prs_decile = quantilize(lvm_prs_std,10),
             prs_ventile = quantilize(lvm_prs_std,20))]
hf_set[,':='(prs_95 = ifelse(prs_ventile==20,1,0),
             prs_90 = ifelse(prs_decile==10,1,0))]

## VT
vt_set <- prs[time_to_vt > 0 & prev_vt == 0]
vt_set[, lvm_prs_std := (lvm_prs_std - mean(lvm_prs_std))/sd(lvm_prs_std)]
vt_set[,':='(prs_decile = quantilize(lvm_prs_std,10),
             prs_ventile = quantilize(lvm_prs_std,20))]
vt_set[,':='(prs_95 = ifelse(prs_ventile==20,1,0),
             prs_90 = ifelse(prs_decile==10,1,0))]

## HCM
hcm_set <- prs[time_to_hcm > 0 & prev_hcm == 0]
hcm_set[, lvm_prs_std := (lvm_prs_std - mean(lvm_prs_std))/sd(lvm_prs_std)]
hcm_set[,':='(prs_decile = quantilize(lvm_prs_std,10),
              prs_ventile = quantilize(lvm_prs_std,20))]
hcm_set[,':='(prs_95 = ifelse(prs_ventile==20,1,0),
              prs_90 = ifelse(prs_decile==10,1,0))]

## DCM
dcm_set <- prs[time_to_dcm > 0 & prev_hcm == 0]
dcm_set[, lvm_prs_std := (lvm_prs_std - mean(lvm_prs_std))/sd(lvm_prs_std)]
dcm_set[,':='(prs_decile = quantilize(lvm_prs_std,10),
              prs_ventile = quantilize(lvm_prs_std,20))]
dcm_set[,':='(prs_95 = ifelse(prs_ventile==20,1,0),
              prs_90 = ifelse(prs_decile==10,1,0))]

## ICD
icd_set <- prs[time_to_icd > 0 & prev_hcm == 0]
icd_set[, lvm_prs_std := (lvm_prs_std - mean(lvm_prs_std))/sd(lvm_prs_std)]
icd_set[,':='(prs_decile = quantilize(lvm_prs_std,10),
              prs_ventile = quantilize(lvm_prs_std,20))]
icd_set[,':='(prs_95 = ifelse(prs_ventile==20,1,0),
              prs_90 = ifelse(prs_decile==10,1,0))]


### AF ANALYSES
mod_af <- coxph(Surv(time_to_af,incd_af) ~ lvm_prs_std + collect_age + sex + PC1 + PC2 + PC3 + PC4 + PC5,data=af_set)
mod_af90 <- coxph(Surv(time_to_af,incd_af) ~ factor(prs_decile) + collect_age + sex + PC1 + PC2 + PC3 + PC4 + PC5,data=af_set)
mod_af95 <- coxph(Surv(time_to_af,incd_af) ~ factor(prs_ventile) + collect_age + sex + PC1 + PC2 + PC3 + PC4 + PC5,data=af_set)
mod_af_top10 <- coxph(Surv(time_to_af,incd_af) ~ prs_90 + collect_age + sex + PC1 + PC2 + PC3 + PC4 + PC5,data=af_set)
mod_af_top5 <- coxph(Surv(time_to_af,incd_af) ~ prs_95 + collect_age + sex + PC1 + PC2 + PC3 + PC4 + PC5,data=af_set)

# Spline
af_spline <- coxph(Surv(time_to_af,incd_af) ~ pspline(lvm_prs_std,df=0) + collect_age + sex + PC1 + PC2 + PC3 + PC4 + PC5,data=af_set)

# Plot spline
resid_af_spline <- termplot(af_spline,term=1,se=TRUE,plot=F)
value_term <- resid_af_spline$lvm_prs_std
center <- with(value_term, y[which.min(abs((x-0)))])
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
prodlim_af <- prodlim(Hist(time_to_af,incd_af)~prs_graphical90,data=af_set)

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
mod_mi <- coxph(Surv(time_to_mi,incd_mi) ~ lvm_prs_std + collect_age + sex + PC1 + PC2 + PC3 + PC4 + PC5,data=mi_set)
mod_mi90 <- coxph(Surv(time_to_mi,incd_mi) ~ factor(prs_decile) + collect_age + sex + PC1 + PC2 + PC3 + PC4 + PC5,data=mi_set)
mod_mi95 <- coxph(Surv(time_to_mi,incd_mi) ~ factor(prs_ventile) + collect_age + sex + PC1 + PC2 + PC3 + PC4 + PC5,data=mi_set)
mod_mi_top10 <- coxph(Surv(time_to_mi,incd_mi) ~ prs_90 + collect_age + sex + PC1 + PC2 + PC3 + PC4 + PC5,data=mi_set)
mod_mi_top5 <- coxph(Surv(time_to_mi,incd_mi) ~ prs_95 + collect_age + sex + PC1 + PC2 + PC3 + PC4 + PC5,data=mi_set)

# Spline
mi_spline <- coxph(Surv(time_to_mi,incd_mi) ~ pspline(lvm_prs_std,df=0) + collect_age + sex + PC1 + PC2 + PC3 + PC4 + PC5,data=mi_set)

# Plot spline
resid_mi_spline <- termplot(mi_spline,term=1,se=TRUE,plot=F)
value_term <- resid_mi_spline$lvm_prs_std
center <- with(value_term, y[which.min(abs((x-0)))])
y <- value_term$y + outer(value_term$se,c(0,-1.96,1.96),'*')

CairoPDF(file='mi_spline.pdf',height=3,width=3.5,
         pointsize=4)
par(oma=c(3,4,1,1),mar=c(3,4,1,1),xpd=FALSE)
matplot(value_term$x,exp(y-center),log='y',type='l',xaxt='n',xlab='',ylab='',bty='n',yaxt='n',
        xlim=c(-3,3),ylim=c(0.5,1.5),col=c('black','#d73027','#1a9850'))
axis(1,at=seq(-3,3,1),cex.axis=1.4)
axis(2,at=seq(0.5,1.5,0.25),las=2,cex.axis=1.4,pos=-3.3)
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
prodlim_mi <- prodlim(Hist(time_to_mi,incd_mi)~prs_graphical90,data=mi_set)

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
mod_hf <- coxph(Surv(time_to_hf,incd_hf) ~ lvm_prs_std + collect_age + sex + PC1 + PC2 + PC3 + PC4 + PC5,data=hf_set)
mod_hf90 <- coxph(Surv(time_to_hf,incd_hf) ~ factor(prs_decile) + collect_age + sex + PC1 + PC2 + PC3 + PC4 + PC5,data=hf_set)
mod_hf95 <- coxph(Surv(time_to_hf,incd_hf) ~ factor(prs_ventile) + collect_age + sex + PC1 + PC2 + PC3 + PC4 + PC5,data=hf_set)
mod_hf_top10 <- coxph(Surv(time_to_hf,incd_hf) ~ prs_90 + collect_age + sex + PC1 + PC2 + PC3 + PC4 + PC5,data=hf_set)
mod_hf_top5 <- coxph(Surv(time_to_hf,incd_hf) ~ prs_95 + collect_age + sex + PC1 + PC2 + PC3 + PC4 + PC5,data=hf_set)

# Spline
hf_spline <- coxph(Surv(time_to_hf,incd_hf) ~ pspline(lvm_prs_std,df=0) + collect_age + sex + PC1 + PC2 + PC3 + PC4 + PC5,data=hf_set)

# Plot spline
resid_hf_spline <- termplot(hf_spline,term=1,se=TRUE,plot=F)
value_term <- resid_hf_spline$lvm_prs_std
center <- with(value_term, y[which.min(abs((x-0)))])
y <- value_term$y + outer(value_term$se,c(0,-1.96,1.96),'*')

CairoPDF(file='hf_spline.pdf',height=3,width=3.5,
         pointsize=4)
par(oma=c(3,4,1,1),mar=c(3,4,1,1),xpd=FALSE)
matplot(value_term$x,exp(y-center),log='y',type='l',xaxt='n',xlab='',ylab='',bty='n',yaxt='n',
        xlim=c(-3,3),ylim=c(0.5,1.5),col=c('black','#d73027','#1a9850'))
axis(1,at=seq(-3,3,1),cex.axis=1.4)
axis(2,at=seq(0.5,1.5,0.25),las=2,cex.axis=1.4,pos=-3.3)
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
prodlim_hf <- prodlim(Hist(time_to_hf,incd_hf)~prs_graphical90,data=hf_set)

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
mod_vt <- coxph(Surv(time_to_vt,incd_vt) ~ lvm_prs_std + collect_age + sex + PC1 + PC2 + PC3 + PC4 + PC5,data=vt_set)
mod_vt90 <- coxph(Surv(time_to_vt,incd_vt) ~ factor(prs_decile) + collect_age + sex + PC1 + PC2 + PC3 + PC4 + PC5,data=vt_set)
mod_vt95 <- coxph(Surv(time_to_vt,incd_vt) ~ factor(prs_ventile) + collect_age + sex + PC1 + PC2 + PC3 + PC4 + PC5,data=vt_set)
mod_vt_top10 <- coxph(Surv(time_to_vt,incd_vt) ~ prs_90 + collect_age + sex + PC1 + PC2 + PC3 + PC4 + PC5,data=vt_set)
mod_vt_top5 <- coxph(Surv(time_to_vt,incd_vt) ~ prs_95 + collect_age + sex + PC1 + PC2 + PC3 + PC4 + PC5,data=vt_set)

# Spline
vt_spline <- coxph(Surv(time_to_vt,incd_vt) ~ pspline(lvm_prs_std,df=0) + collect_age + sex + PC1 + PC2 + PC3 + PC4 + PC5,data=vt_set)

# Plot spline
resid_vt_spline <- termplot(vt_spline,term=1,se=TRUE,plot=F)
value_term <- resid_vt_spline$lvm_prs_std
center <- with(value_term, y[which.min(abs((x-0)))])
y <- value_term$y + outer(value_term$se,c(0,-1.96,1.96),'*')

CairoPDF(file='vt_spline.pdf',height=3,width=3.5,
         pointsize=4)
par(oma=c(3,4,1,1),mar=c(3,4,1,1),xpd=FALSE)
matplot(value_term$x,exp(y-center),log='y',type='l',xaxt='n',xlab='',ylab='',bty='n',yaxt='n',
        xlim=c(-3,3),ylim=c(0.5,1.5),col=c('black','#d73027','#1a9850'))
axis(1,at=seq(-3,3,1),cex.axis=1.4)
axis(2,at=seq(0.5,1.5,0.25),las=2,cex.axis=1.4,pos=-3.3)
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
prodlim_vt <- prodlim(Hist(time_to_vt,incd_vt)~prs_graphical90,data=vt_set)

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
mod_hcm <- coxph(Surv(time_to_hcm,incd_hcm) ~ lvm_prs_std + collect_age + sex + PC1 + PC2 + PC3 + PC4 + PC5,data=hcm_set)
mod_hcm90 <- coxph(Surv(time_to_hcm,incd_hcm) ~ factor(prs_decile) + collect_age + sex + PC1 + PC2 + PC3 + PC4 + PC5,data=hcm_set)
mod_hcm95 <- coxph(Surv(time_to_hcm,incd_hcm) ~ factor(prs_ventile) + collect_age + sex + PC1 + PC2 + PC3 + PC4 + PC5,data=hcm_set)
mod_hcm_top10 <- coxph(Surv(time_to_hcm,incd_hcm) ~ prs_90 + collect_age + sex + PC1 + PC2 + PC3 + PC4 + PC5,data=hcm_set)
mod_hcm_top5 <- coxph(Surv(time_to_hcm,incd_hcm) ~ prs_95 + collect_age + sex + PC1 + PC2 + PC3 + PC4 + PC5,data=hcm_set)

# Spline
hcm_spline <- coxph(Surv(time_to_hcm,incd_hcm) ~ pspline(lvm_prs_std,df=0) + collect_age + sex + PC1 + PC2 + PC3 + PC4 + PC5,data=hcm_set)

# Plot spline
resid_hcm_spline <- termplot(hcm_spline,term=1,se=TRUE,plot=F)
value_term <- resid_hcm_spline$lvm_prs_std
center <- with(value_term, y[x==min(abs((x-0)))])
y <- value_term$y + outer(value_term$se,c(0,-1.96,1.96),'*')

CairoPDF(file='hcm_spline.pdf',height=3,width=3.5,
         pointsize=4)
par(oma=c(3,4,1,1),mar=c(3,4,1,1),xpd=FALSE)
matplot(value_term$x,exp(y-center),log='y',type='l',xaxt='n',xlab='',ylab='',bty='n',yaxt='n',
        xlim=c(-2.4,2.4),ylim=c(0.5,1.5),col=c('black','#d73027','#1a9850'))
axis(1,at=seq(-2.5,2.5,0.5),cex.axis=1.4)
axis(2,at=seq(0.5,1.5,0.25),las=2,cex.axis=1.4,pos=-2.6)
legend(-4,0.05,c('estimate','upper 95% CI','lower 95% CI'),bty='n',lty=1,
       col=c('black','#d73027','#1a9850'),cex=1.4)
mtext("Standardized PRS",1,line=3,cex=1.5)
mtext("Relative hazard for hcm compared to PRS of zero",2,line=5.5,cex=1.5)
dev.off()

# Graphical variables
hcm_set[,prs_graphical90 := factor(ifelse(prs_ventile==20,"Above 95th percentile","Below 95th percentile"),levels=c('Above 95th percentile','Below 95th percentile'))]

# HCM KM
prodlim_hcm <- prodlim(Hist(time_to_hcm,incd_hcm)~prs_graphical90,data=hcm_set)

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
mod_dcm <- coxph(Surv(time_to_dcm,incd_dcm) ~ lvm_prs_std + collect_age + sex + PC1 + PC2 + PC3 + PC4 + PC5,data=dcm_set)
mod_dcm90 <- coxph(Surv(time_to_dcm,incd_dcm) ~ factor(prs_decile) + collect_age + sex + PC1 + PC2 + PC3 + PC4 + PC5,data=dcm_set)
mod_dcm95 <- coxph(Surv(time_to_dcm,incd_dcm) ~ factor(prs_ventile) + collect_age + sex + PC1 + PC2 + PC3 + PC4 + PC5,data=dcm_set)
mod_dcm_top10 <- coxph(Surv(time_to_dcm,incd_dcm) ~ prs_90 + collect_age + sex + PC1 + PC2 + PC3 + PC4 + PC5,data=dcm_set)
mod_dcm_top5 <- coxph(Surv(time_to_dcm,incd_dcm) ~ prs_95 + collect_age + sex + PC1 + PC2 + PC3 + PC4 + PC5,data=dcm_set)

# Spline
dcm_spline <- coxph(Surv(time_to_dcm,incd_dcm) ~ pspline(lvm_prs_std,df=0) + collect_age + sex + PC1 + PC2 + PC3 + PC4 + PC5,data=dcm_set)

# Plot spline
resid_dcm_spline <- termplot(dcm_spline,term=1,se=TRUE,plot=F)
value_term <- resid_dcm_spline$lvm_prs_std
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
prodlim_dcm <- prodlim(Hist(time_to_dcm,incd_dcm)~prs_graphical90,data=dcm_set)

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
mod_icd <- coxph(Surv(time_to_icd,incd_icd) ~ lvm_prs_std + collect_age + sex + PC1 + PC2 + PC3 + PC4 + PC5,data=icd_set)
mod_icd90 <- coxph(Surv(time_to_icd,incd_icd) ~ factor(prs_decile) + collect_age + sex + PC1 + PC2 + PC3 + PC4 + PC5,data=icd_set)
mod_icd95 <- coxph(Surv(time_to_icd,incd_icd) ~ factor(prs_ventile) + collect_age + sex + PC1 + PC2 + PC3 + PC4 + PC5,data=icd_set)
mod_icd_top10 <- coxph(Surv(time_to_icd,incd_icd) ~ prs_90 + collect_age + sex + PC1 + PC2 + PC3 + PC4 + PC5,data=icd_set)
mod_icd_top5 <- coxph(Surv(time_to_icd,incd_icd) ~ prs_95 + collect_age + sex + PC1 + PC2 + PC3 + PC4 + PC5,data=icd_set)

# Spline
icd_spline <- coxph(Surv(time_to_icd,incd_icd) ~ pspline(lvm_prs_std,df=0) + collect_age + sex + PC1 + PC2 + PC3 + PC4 + PC5,data=icd_set)

# Plot spline
resid_icd_spline <- termplot(icd_spline,term=1,se=TRUE,plot=F)
value_term <- resid_icd_spline$lvm_prs_std
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
mtext("Relative hazard for ICD compared to PRS of zero",2,line=5.5,cex=1.5)
dev.off()

# Graphical variables
icd_set[,prs_graphical90 := factor(ifelse(prs_ventile==20,"Above 95th percentile","Below 95th percentile"),levels=c('Above 95th percentile','Below 95th percentile'))]

# icd KM
prodlim_icd <- prodlim(Hist(time_to_icd,incd_icd)~prs_graphical90,data=icd_set)

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
mtext("Cumulative risk of ICD (%)",side=2,line=-1.2,at=0.05,cex=2.5)
mtext("Years",side=1, line=-0.6,cex=2.5)
mtext('Stratum',side=1, line=-0.3,cex=1.8,at=-2.5)
dev.off()