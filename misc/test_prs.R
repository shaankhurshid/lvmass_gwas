# Script to analyze PRS-pheWAS

# Depends 
library(data.table)
library(ggplot2)

# Load PRS
prs <- fread(file='prs_processed.csv')

# Load LV mass
seg <- fread(file='lvmi_seg_adjusted.tsv')

# Load processed LV mass
seg_lvm <- fread(file='seg_lvm.csv')
seg_lvm[c(!is.na(sex) & !is.na(lvmi_seg_adjusted)),
         ':='(lvh_ukbb = ifelse(c((sex=='Female') & (lvmi_seg_adjusted > 55)),1,
                           ifelse(lvmi_seg_adjusted > 72,1,0)))]

# Join
setkey(prs,sample_id); setkey(seg,IID); setkey(seg_lvm,sample_id)
seg[prs,':='(prs_std = i.prs_std)]
seg[seg_lvm,lvh_ukbb := i.lvh_ukbb]

# Remove non-white/no PRS
seg <- seg[!is.na(prs_std)]

# Check exposure ~ PRS
mod <- lm(lvmi_seg_adjusted ~ prs_std,data=seg)
mod_adj <- lm(lvmi_seg_adjusted ~ prs_std + male + age_at_mri + PC1 + PC2 + PC3 + PC4 + PC5,data=seg)

# Corr
corr <- cor.test(seg$prs_std,y=seg$lvmi_seg_adjusted)

# Plot
pdf(file='prs_lvmi_corr_ukbb.pdf',height=3,width=3,pointsize=5)
par(mar=c(3,3,1,1),oma=c(2,2,1,1))
plot(x=seg$prs_std,y=seg$lvmi_seg_adjusted,bty='n',xlab='',ylab='',xaxt='n',yaxt='n',pch=19,col='#2171b58C',
     xlim=c(-5,5),ylim=c(0,150))
axis(1,cex.axis=1,at=c(-4:4))
axis(2,cex.axis=1,at=seq(0,150,25),las=2)
mtext("Standardized PRS",1,line=3,cex=1.5)
mtext(expression(paste("Indexed LV mass (g/m"^2,")")),2,line=3,cex=1.5)
text(-2.5,145,labels="r=0.29, p<0.01")
text(-2.5,138,labels="95% CI 0.28-0.30")
dev.off()

# Plot distribution stratified by LVH
lvh <- seg[lvh_ukbb==1]
no_lvh <- seg[lvh_ukbb==0]
x <- list(v1=lvh$prs_std,v2=no_lvh$prs_std)
data <- melt(x)

ggplot() + geom_density(data=data,aes(x=value,fill=L1),alpha=0.55) +
  scale_x_continuous(breaks=seq(-3.5,3.5,0.5),expand=c(0.01,0),limits=c(-3.5,3.5)) +
  scale_y_continuous(breaks=seq(0,0.5,0.1),expand=c(0,0),limits=c(0,0.5)) +
  scale_fill_manual(values=c("#2b8cbe","#f03b20"),name='',labels=c('LVH','No LVH')) +
  theme(panel.background=element_blank(),axis.line=element_line(color='black'),legend.position=c(0.20,0.90),
        axis.text=element_text(size=20,color='black'),plot.margin=unit(c(0.5,0.6,0.5,0.5),'cm'),
        axis.title.y = element_text(size=20,margin = margin(t = 0, r = 10, b = 0, l = 0)),
        axis.title.x = element_text(size=20),legend.text=element_text(size=20)) +
  labs(x=expression(paste("LVMI Standardized PRS")),y='Density')
ggsave(filename='prs_std_density_lvh_strat_eur.pdf',height=1.8,width=2.53,
       scale=4,device='pdf')
