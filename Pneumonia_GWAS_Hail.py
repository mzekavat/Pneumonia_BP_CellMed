### Pneumonia GWAS followed by 2-sample MR IN UKBB using Hail-0.2:
### Author: Seyedeh Maryam Zekavat 

### Run the line below to activate hail on terminal via: 
conda activate hail

### sample code for making a Hail cluster and open up a jupyter notebook to run Hail interactively: 
hailctl dataproc start mz02 --master-machine-type n1-highmem-16 --worker-machine-type n1-highmem-16 --worker-boot-disk-size 200 --num-workers 40 --num-preemptible-workers 30 --master-boot-disk-size 100 --region us-east1 --zone us-east1-d --requester-pays-allow-all --vep GRCh37 --properties "spark:spark.driver.memory=90G,spark:spark.driver.maxResultSize=50G,spark:spark.kryoserializer.buffer.max=1G,spark:spark.task.maxFailures=20,spark:spark.driver.extraJavaOptions=-Xss4M,spark:spark.executor.extraJavaOptions=-Xss4M,spark:spark.speculation=true"
hailctl dataproc connect mz02 notebook --zone us-east1-d --region us-east1

### Step 1: annotating UKBB variants with VEP
import hail as hl
from pprint import pprint
from bokeh.io import output_notebook,show,save
from bokeh.layouts import gridplot
from bokeh.models import Span
import hail.expr.aggregators as agg
from bokeh.plotting import figure, output_file
import numpy as np
​
from bokeh.io import show, output_notebook
from bokeh.layouts import gridplot
output_notebook()
​
hl.init(default_reference='GRCh37')
​
## Variant level annotations (VEP annotations; annotated separately)
mt5 = hl.read_table('gs://ukbb_v2/projects/mzekavat/ukbb_v3.AllAutosomalANDchrX.annotations.ht')
## UKBB imputed bgens:
ds = hl.import_bgen('gs://fc-7d5088b4-7673-45b5-95c2-17ae00a04183/imputed/ukb_imp_chr{1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22}_v3.bgen',entry_fields = ['GT'],sample_file='gs://ukbb_v2/data/ukb7089_imp_chr3_v3_s487395.sample')
## Phenotype file
phenos = hl.import_table('gs://ukbb_v2/projects/mzekavat/Pneumonia_GWAS/ukbb_PhenoFile.ALL_500k_incidPrevCases.plusRespPhenos.plusBPMeds.plusPFTs.plusCHIP.QCed.txt.gz',force_bgz=True,key = 'id',types={'id':hl.tstr},impute=True)

ds = ds.annotate_rows(**mt5.index(ds.row_key))
ds = ds.annotate_cols(pheno = phenos[ds.col_key])
ds = ds.annotate_cols(array = hl.if_else((ds.pheno.genotyping_array == "UKBB"), 1, 0))
ds = ds.filter_cols(hl.is_defined(ds.pheno.age), keep=True)

### variant qc 
mt = hl.variant_qc(ds,name='variant_qc')
mt = mt.filter_rows( ((mt.variant_qc.AF[1] > 0.0001) & (mt.variant_qc.AF[1] < 1) & (mt.info>0.4) & (mt.variant_qc.p_value_hwe >= 0.0000000001)),keep = True )
final= mt.annotate_rows(AF = mt.variant_qc.AF[1],AC = mt.variant_qc.AC[1],AN = mt.variant_qc.AN)
#final_annot = final.annotate_rows(HWE = final.variant_qc.p_value_hwe, callRate = final.variant_qc.call_rate)
#final_annot = final_annot.drop('variant_qc').rows()
### gwas logistic regression wald
gwas = hl.logistic_regression_rows(test='wald',\
									y=final.pheno.All_Pneumonia,\
									x=final.GT.n_alt_alleles(),\
									covariates=[1, final.pheno.age,final.pheno.age2, final.pheno.Sex_numeric, final.pheno.ever_smoked, final.pheno.PC1,final.pheno.PC2,final.pheno.PC3,final.pheno.PC4,final.pheno.PC5,final.pheno.PC6,final.pheno.PC7,final.pheno.PC8,final.pheno.PC9,final.pheno.PC10,final.array],
									pass_through=['rsid','Gene','Consequence','clin_sig', 'metasvm','LOF_LOFTEE','PolyPhen','SIFT','hgvsp','AF', 'AC', 'AN','info'])
​
### Writting out the annotated GWAS results:
gwas.flatten().export('gs://ukbb_v2/projects/mzekavat/Pneumonia_GWAS/logreg_wald_All_Pneumonia.tsv.bgz')
gwas.write('gs://ukbb_v2/projects/mzekavat/Pneumonia_GWAS/logreg_wald_All_Pneumonia.ht')
gwas = hl.read_table('gs://ukbb_v2/projects/mzekavat/Pneumonia_GWAS/logreg_wald_All_Pneumonia.ht')
gwas_v2 = gwas.filter(gwas.p_value<0.0001, keep=True)


### Filtering the pneumonia GWAS to just the SBP and DBP SNPs:
SBPsnps = hl.import_table('gs://ukbb_v2/projects/mzekavat/Pneumonia_GWAS/SBP_75SNP_instrument_hg37.txt', impute = True) 
DBPsnps = hl.import_table('gs://ukbb_v2/projects/mzekavat/Pneumonia_GWAS/DBP_75SNP_instrument_hg37.txt', impute = True) 


SBP_SNPs =  [row['SNP'] for row in SBPsnps.select(SBPsnps.SNP).collect()]

gwas_SBPvar =gwas.filter(hl.literal(SBP_SNPs).contains(gwas.rsid), keep=True)

DBP_SNPs =  [row['SNP'] for row in DBPsnps.select(DBPsnps.SNP).collect()]

gwas_DBPvar =gwas.filter(hl.literal(DBP_SNPs).contains(gwas.rsid), keep=True)


gwas_SBPvar.export('gs://ukbb_v2/projects/mzekavat/Pneumonia_GWAS/gwas_SBPvar.logreg_wald_All_Pneumonia.tsv.bgz')
gwas_DBPvar.export('gs://ukbb_v2/projects/mzekavat/Pneumonia_GWAS/gwas_DBPvar.logreg_wald_All_Pneumonia.tsv.bgz')


### Performing 2-sample MR in R-3.5: 

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

