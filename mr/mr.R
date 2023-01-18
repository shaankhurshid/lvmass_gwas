  library(data.table)
  library(MendelianRandomization)
  library(ggplot2)
  
  mr <- function(exposures, outcomes, output_dir) {
    exp <- NULL; out <- NULL
    for (i in exposures) {
      # Read in data
      if (i == "SBP") {
        exp <- fread("/Volumes/medpop_esp2/karagam/UK_Biobank/Kiran_B/Projects/LVM_MR/MR/data/SBP_exp_weights.csv")
      } else if (i == "DBP") {
        exp <- fread("/Volumes/medpop_esp2/karagam/UK_Biobank/Kiran_B/Projects/LVM_MR/MR/data/DBP_exp_weights.csv")
      } else if (i == "T2D") {
        exp <- fread("/Volumes/medpop_esp2/karagam/UK_Biobank/Kiran_B/Projects/LVM_MR/MR/data/T2D_exp_weights.csv")
      } else {
        print("Invalid exposure."); break
      }
      
      # Format data
      exp <- exp[, c("rsid", "Effect_allele", "Ref_allele", "Effect_beta", "SE")]
      colnames(exp)[4:5] <- c("beta_exp", "se_exp")
      
      for (j in outcomes) {
        # make directory for each outcome
        dir.create(file.path(output_dir, j), showWarnings = FALSE)
        
        # Read in data -- CHANGE PATHS HERE IF UPDATING DATA
        if (j == "LVMi") {
          out <- fread("/Volumes/medpop_afib/skhurshid/lvm_gwas/gwas/v23_seg_white/v23_lvmi_seg_white_bolt_imputed.filtered.tsv")
        } else if (j == "LVM") {
          out <- fread("/Volumes/medpop_afib/skhurshid/lvm_gwas/gwas/v44_lvm_white/v44_lvm_seg_adjusted_white_bolt_imputed.filtered.tsv")
        } else {
          print("Invalid outcome")
          break
        }
        
        # Format data
        out$beta_out <- out$BETA
        out$se_out <- out$SE
        
        # align data
        df <- merge(exp, out[,c("SNP", "ALLELE1", "beta_out", "se_out")], by.x = "rsid", by.y = "SNP")
        table(df$ALLELE1 != df$Effect_allele)
        df$beta_out[df$ALLELE1 != df$Effect_allele] <- -1 * df$beta_out[df$ALLELE1 != df$Effect_allele]
        df$ALT <- df$Effect_allele
        df$REF <- df$Ref_allele
        df$ALT[df$ALLELE1 != df$Effect_allele] <- df$REF[df$ALLELE1 != df$Effect_allele]
        df$REF[df$ALLELE1 != df$Effect_allele] <- df$Effect_allele[df$ALLELE1 != df$Effect_allele]
        df$ALLELE1 <- NULL
        
        # remove any SNPs with different allele patterns
        df <- df[!is.na(df$beta_out),]
        df <- df[!is.na(df$beta_exp),]
        
        # Mendelian randomization
        wm_beta_list <- c(); wm_se_list <- c(); wm_lower_95_list <- c(); wm_upper_95_list <- c(); wm_pval_list <- c()
        ivw_beta_list <- c(); ivw_se_list <- c(); ivw_lower_95_list <- c(); ivw_upper_95_list <- c(); ivw_pval_list <- c()
        egger_beta_list <- c(); egger_se_list <- c(); egger_lower_95_list <- c(); egger_upper_95_list <- c(); egger_pval_list <- c()
        egger_intercept_list <- c()
        
        # Input
        input <- mr_input(bx = df$beta_exp,bxse = df$se_exp, by = df$beta_out, 
                          byse = df$se_out); test <- mr_allmethods(input)
        # Weighted Median Results
        wm_beta_list <- c(wm_beta_list, test$Values[2,2]); wm_se_list <- c(wm_se_list, test$Values[2,3])
        wm_lower_95_list <- c(wm_lower_95_list, test$Values[2,4]); wm_upper_95_list <- c(wm_upper_95_list, test$Values[2,5])
        wm_pval_list <- c(wm_pval_list, test$Values[2,6])
        # IVW Results
        ivw_beta_list <- c(ivw_beta_list, test$Values[4,2]); ivw_se_list <- c(ivw_se_list, test$Values[4,3])
        ivw_lower_95_list <- c(ivw_lower_95_list, test$Values[4,4]); ivw_upper_95_list <- c(ivw_upper_95_list, test$Values[4,5])
        ivw_pval_list <- c(ivw_pval_list, test$Values[4,6])
        # MR-Egger Results
        egger_beta_list <- c(egger_beta_list, test$Values[8,2]); egger_se_list <- c(egger_se_list, test$Values[8,3])
        egger_lower_95_list <- c(egger_lower_95_list, test$Values[8,4]); egger_upper_95_list <- c(egger_upper_95_list, test$Values[8,5])
        egger_pval_list <- c(egger_pval_list, test$Values[8,6])
        # MR-Egger Intercept Results
        egger_intercept_list <- c(egger_intercept_list, paste0(round(test$Values[9,2],2), " (", round(test$Values[9,6],2), ")"))
        
        # Write data
        WeightedMedian <- cbind(wm_beta_list, wm_se_list, wm_lower_95_list, wm_upper_95_list, wm_pval_list)
        colnames(WeightedMedian) <- c("Beta", "SE", "95% CI Lower", "95% CI Upper", "P-value")
        rownames(WeightedMedian) <- i
        write.csv(WeightedMedian, paste0(output_dir, "/", j, "/", i, "-", j,"-WeightedMedian-2SMRResults.csv"))
        IVW <- cbind(ivw_beta_list, ivw_se_list, ivw_lower_95_list, ivw_upper_95_list, ivw_pval_list)
        colnames(IVW) <- c("Beta", "SE", "95% CI Lower", "95% CI Upper", "P-value")
        rownames(IVW) <- i
        print(IVW)
        print(paste0(output_dir, "/", j, "/", i, "-", j, "-IVW-2SMRResults.csv"))
        write.csv(IVW, paste0(output_dir, "/", j, "/", i, "-", j, "-IVW-2SMRResults.csv"))
        egger <- cbind(egger_beta_list, egger_se_list, egger_lower_95_list, egger_upper_95_list, egger_pval_list, egger_intercept_list)
        colnames(egger) <- c("Beta", "SE", "95% CI Lower", "95% CI Upper", "P-value", "MR-Egger Intercept (P-value)")
        rownames(egger) <- i
        write.csv(egger, paste0(output_dir, "/", j, "/", i, "-", j, "-Egger-2SMRResults.csv"))
        
        plot_path <- paste0(output_dir, "/", j, "/", i, "-", j, "-ivw_plot.jpg")
        mr_plot(input,line='ivw',interactive=FALSE,orientate=TRUE,error=FALSE)
        ggsave(plot_path,height=5,width=7,dpi=300,device='jpeg',units='in')
      }
    }
  }
  
  # Run with specified outcome directory, will create subdirectories for each pheno (MAKE SURE NO "/" AT END OF OUTPUT DIR)
  
  output_path <- "~/Documents/MGH Research/lvm_gwas/mr"
  
  exposures <- c("SBP","DBP","T2D")
  outcomes <- c("LVMi","LVM")
  out <- mr(exposures=exposures,outcomes=outcomes,output_dir=output_path)
