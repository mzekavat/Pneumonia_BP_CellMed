### Performing 2-sample MR in R-3.5 between SBP/DBP variants from the ICBP consortiom (exposure) and Pneumonia (from the UK Biobank Pneumonia GWAS): 
### Author: Seyedeh Maryam Zekavat 

R
library(MendelianRandomization)
library(data.table)

### Loading in the SBP and DBP instruments (which include SNP, Beta, SE, P for each of the 75 independent variants associated with SBP and DBP, respectively)
SBPinstr = fread("/medpop/esp2/mzekavat/Miscellaneous_Tables/GWAS_Results/ICBP_Instrument/ICBP.SBP_vars.forMR.txt", header=TRUE, sep="\t")
SBPinstr = data.frame(SBPinstr)

DBPinstr = fread("/medpop/esp2/mzekavat/Miscellaneous_Tables/GWAS_Results/ICBP_Instrument/ICBP.DBP_vars.forMR.txt", header=TRUE, sep="\t")
DBPinstr = data.frame(DBPinstr)

### Loading in the SBP and DBP respective variants from the pneumonia GWAS: 
SBPinstr_PneumOutc = fread("/medpop/esp2/mzekavat/UKBB/Pneumonia_HTN/gwas_SBPvar.logreg_wald_All_Pneumonia.tsv.gz", header=TRUE, sep="\t")
SBPinstr_PneumOutc = data.frame(SBPinstr_PneumOutc)

DBPinstr_PneumOutc = fread("/medpop/esp2/mzekavat/UKBB/Pneumonia_HTN/gwas_DBPvar.logreg_wald_All_Pneumonia.tsv.gz", header=TRUE, sep="\t")
DBPinstr_PneumOutc = data.frame(SBPinstr_PneumOutc)

colnames(SBPinstr_PneumOutc) = paste("Pneum_", colnames(SBPinstr_PneumOutc), sep="")
colnames(DBPinstr_PneumOutc) = paste("Pneum_", colnames(DBPinstr_PneumOutc), sep="")
colnames(SBPinstr) = paste("SBP_",colnames(SBPinstr),sep="")
colnames(DBPinstr) = paste("DBP_",colnames(DBPinstr),sep="")

merged_SBP_Pneum = merge(SBPinstr,SBPinstr_PneumOutc,by.x=1, by.y=3, all.x=TRUE)
merged_DBP_Pneum = merge(DBPinstr,DBPinstr_PneumOutc,by.x=1, by.y=3, all.x=TRUE)

###  setting up the betas such that they reflect the same effect allele for both BP and pneumonia:
merged_DBP_Pneum$new_Pneum_beta = ifelse(merged_DBP_Pneum$Pneum_alt ==merged_DBP_Pneum$DBP_Effect.allele, merged_DBP_Pneum$Pneum_beta,
										ifelse(merged_DBP_Pneum$Pneum_alt ==merged_DBP_Pneum$DBP_ICBP_A2, -1*merged_DBP_Pneum$Pneum_beta, NA))

merged_SBP_Pneum$new_Pneum_beta = ifelse(merged_SBP_Pneum$Pneum_alt ==merged_SBP_Pneum$SBP_Effect.allele, merged_SBP_Pneum$Pneum_beta,
										ifelse(merged_SBP_Pneum$Pneum_alt ==merged_SBP_Pneum$SBP_ICBP_A2, -1*merged_SBP_Pneum$Pneum_beta, NA))

# Running MR for SBP on Pneumonia:
mr_allmethods(mr_input(bx = merged_SBP_Pneum$SBP_Beta, bxse = merged_SBP_Pneum$SBP_ICBP_se_SBP,by = merged_SBP_Pneum$new_Pneum_beta,byse =  merged_SBP_Pneum$Pneum_standard_error),method = "all")

# Running MR for DBP on Pneumonia:
mr_allmethods(mr_input(bx = merged_DBP_Pneum$DBP_Beta, bxse = merged_DBP_Pneum$DBP_ICBP_se_DBP,by = merged_DBP_Pneum$new_Pneum_beta,byse =  merged_DBP_Pneum$Pneum_standard_error),method = "all")

### Sensitivity analysis using the TwoSampleMR package: 
library(TwoSampleMR)

####### DBP
outcome_dat <- read_outcome_data(
    filename = "/medpop/esp2/mzekavat/UKBB/Pneumonia_HTN/merged_DBP_Pneum.phenotypes.txt",
    sep = "\t",
    snp_col = "DBP_SNP",
    beta_col = "new_Pneum_beta",
    se_col = "Pneum_standard_error",
    effect_allele_col = "DBP_Effect.allele",
    other_allele_col = "DBP_ICBP_A2",
    eaf_col = "Pneum_AF",
    pval_col = "Pneum_p_value"

    )
#    eaf_col = "a1_freq",

exposure_dat <- read_exposure_data(
    filename = "/medpop/esp2/mzekavat/UKBB/Pneumonia_HTN/merged_DBP_Pneum.phenotypes.txt",
    sep = "\t",
    snp_col = "DBP_SNP",
    beta_col = "DBP_Beta",
    se_col = "DBP_ICBP_se_DBP",
    effect_allele_col = "DBP_Effect.allele",
    other_allele_col = "DBP_ICBP_A2",
    eaf_col = "Pneum_AF",
    pval_col = "DBP_ICBP_P_DBP"
    )

dat = harmonise_data(outcome_dat, exposure_dat)

### Performing MR using mr_raps and mr_ivw and making scatter plots:
res <- mr(dat, method_list=c("mr_raps", "mr_ivw"))
p1 <- mr_scatter_plot(res, dat)+theme(axis.text.y=element_text(size=18, hjust=1, color='black'),axis.text.x=element_text(size=18, hjust=1, color='black'),axis.title.x= element_text(size=18))+ xlab("SNP effect on SBP ")+ ylab("SNP effect on Pneumonia")

pdf(paste("/medpop/esp2/mzekavat/UKBB/Pneumonia_HTN/Pneumonia_SBP_MRReport/MRadjforest_DBP_Pneum.xyPLOT.pdf",sep=""), width = 4, height= 5)
p1[[1]]
dev.off()

#########SBP:
outcome_dat <- read_outcome_data(
    filename = "/medpop/esp2/mzekavat/UKBB/Pneumonia_HTN/merged_SBP_Pneum.phenotypes.txt",
    sep = "\t",
    snp_col = "SBP_SNP",
    beta_col = "new_Pneum_beta",
    se_col = "Pneum_standard_error",
    effect_allele_col = "SBP_Effect.allele",
    other_allele_col = "SBP_ICBP_A2",
    eaf_col = "Pneum_AF",
    pval_col = "Pneum_p_value"

    )
#    eaf_col = "a1_freq",

exposure_dat <- read_exposure_data(
    filename = "/medpop/esp2/mzekavat/UKBB/Pneumonia_HTN/merged_SBP_Pneum.phenotypes.txt",
    sep = "\t",
    snp_col = "SBP_SNP",
    beta_col = "SBP_Beta",
    se_col = "SBP_ICBP_se_SBP",
    effect_allele_col = "SBP_Effect.allele",
    other_allele_col = "SBP_ICBP_A2",
    eaf_col = "Pneum_AF",
    pval_col = "SBP_ICBP_P_SBP"
    )

dat = harmonise_data(outcome_dat, exposure_dat)

### Performing MR using mr_raps and mr_ivw and making scatter plots:
res <- mr(dat, method_list=c("mr_raps", "mr_ivw"))
p1 <- mr_scatter_plot(res, dat)+theme(axis.text.y=element_text(size=18, hjust=1, color='black'),axis.text.x=element_text(size=18, hjust=1, color='black'),axis.title.x= element_text(size=18))+ xlab("SNP effect on SBP ")+ ylab("SNP effect on Pneumonia")


pdf(paste("/medpop/esp2/mzekavat/UKBB/Pneumonia_HTN/Pneumonia_SBP_MRReport/MRadjforest_SBPactual_Pneum.xyPLOT.pdf",sep=""), width = 4, height= 5)
p1[[1]]
dev.off()
