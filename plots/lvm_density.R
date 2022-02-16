# Script to visualize LV mass estimates

# Dependencies
library(data.table)
library(survival)
library(prodlim)
library(plyr)
library(ggplot2)

# Load MRI data
seg_lvm <- fread(file='seg_lvm.csv')

## Ensure exclusions removed
withdrawals <- fread('w7089_20210809.csv') # UKBB withdrawals
seg_lvm <- seg_lvm[!(sample_id %in% withdrawals$V1)]

# Density plot for unindexed
ggplot() + geom_density(data=seg_lvm,aes(x=seg_lvm$lvm_seg_adjusted),fill="#fb6a4a",alpha=0.55) +
  scale_x_continuous(breaks=seq(0,250,50),expand=c(0,0),limits=c(0,250)) +
  scale_y_continuous(breaks=seq(0,0.02,0.01),expand=c(0,0),limits=c(0,0.02)) +
  theme(panel.background=element_blank(),axis.line=element_line(color='black'),legend.position='None',
        axis.text=element_text(size=20,color='black'),plot.margin=unit(c(0.5,0.6,0.5,0.5),'cm'),
        axis.title.y = element_text(size=20,margin = margin(t = 0, r = 10, b = 0, l = 0)),
        axis.title.x = element_text(size=20),legend.text=element_text(size=20)) +
  labs(x="CMR-derived LV mass (g)",y='Density')
ggsave(filename='lvm_density_ukbb.pdf',height=1.8,width=2.53,
       scale=4,device='pdf')

# Density plot for indexed
ggplot() + geom_density(data=seg_lvm,aes(x=seg_lvm$lvmi_seg_adjusted),fill="#2b8cbe",alpha=0.55) +
  scale_x_continuous(breaks=seq(0,125,25),expand=c(0,0),limits=c(0,125)) +
  scale_y_continuous(breaks=seq(0,0.04,0.01),expand=c(0,0),limits=c(0,0.04)) +
  theme(panel.background=element_blank(),axis.line=element_line(color='black'),legend.position='None',
        axis.text=element_text(size=20,color='black'),plot.margin=unit(c(0.5,0.6,0.5,0.5),'cm'),
        axis.title.y = element_text(size=20,margin = margin(t = 0, r = 10, b = 0, l = 0)),
        axis.title.x = element_text(size=20),legend.text=element_text(size=20)) +
  labs(x=expression(paste("CMR-derived indexed LV mass (g/m"^"2",")")),y='Density')
ggsave(filename='lvmi_density_ukbb.pdf',height=1.8,width=2.53,
       scale=4,device='pdf')