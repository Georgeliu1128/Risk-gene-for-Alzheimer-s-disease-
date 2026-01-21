#!/usr/bin/env Rscript

.libPaths("/nfs/yaofss2/data/shared/Software/R_libs/4.4")
rm(list = ls())

library(ggplot2)
library(dplyr)
library(tidyr)
library(readr)
library(stringr)
library(forcats)
library(tibble)
library(purrr)
#source("/nfs/yangfss2/public/Useful_Tools/util_func.r")

# cd /nfs/yangfss2/data/commons/jyang51/ROSMAP_motorTrip_analysis/

out_dir = "/home/zliu9/zliu9/ROSMAP2025/proteimics_input_for_GEMMA/modified_cov/"

tissue = "Spinal_cord_protein"

protein_file = "/nfs/yangfss2/data/shared/AMP-AD/ROSMAP/Proteomics/2025-11-07_protein_500_spinal_cord/ROSMAP_500_spinal_cord.txt"

######### Load gene annotation #############
gene_anno <- read_tsv("/home/zliu9/yangfss2/public/Gene_Annotation/gene_annot.tsv")

############ Load phenotype data #############
pheno_data <- read_tsv("/home/zliu9/yangfss2/shared/AMP-AD/ROSMAP/Phenotype/Phenotypes_2022/processed/data_long_last_basic_combined_nov_2024.tsv")

######### load proteomics data ############
Protein_dat <- read_tsv(protein_file) %>% data.frame() # n=480; 9324 proteins update
dim(Protein_dat)

#convert "SampleName"  to "projid" match the phenotype
protein_projid <- Protein_dat %>%
  mutate(projid = str_remove(SampleName, "^X")) %>%   # Keep as character to preserve leading zeros
  select(projid)

dim(protein_projid)
head(protein_projid)
colnames(protein_projid) = "projid"

write_delim(protein_projid, "projid.txt", delim = "\t")

#Only keep the protein data and the na.number is less than half
Protein_dt <- Protein_dat[, -(1:8)] # should be 8
na_number <- apply(Protein_dt, 2, function(x){sum(is.na(x))})
Protein_dt <- Protein_dt[, na_number <= 240]  #half of samples
Protein_dt <- Protein_dt %>% mutate(across(everything(),
                ~ ifelse(is.na(.), mean(., na.rm = TRUE), .)))

protein_list <- colnames(Protein_dt) # 9870 proteins; gene_names
n_proteins <- length(protein_list)
n_proteins

#compute and select the top5 PCs
pca_protein <- prcomp(Protein_dt, center = TRUE, scale = TRUE)
top_5_pcs <- pca_protein$x[, 1:5]
Protein_dt <- cbind(Protein_dt, top_5_pcs)

#protein_projid <- read.table("/home/zliu9/zliu9/ROSMAP2025/Spinal_cord_projid.txt", header = TRUE, fill = TRUE)  
Protein_dt$projid <- unlist(protein_projid)

#### Subset phenotype data for the sequenced samples ###
pheno_data_sub <- pheno_data %>% filter(projid %in% Protein_dt$projid) %>% select(c(projid, age_death, msex, educ, pmi, bmi, gait_speed, motor_dexterity, motor_gait, motor_handstreng, motor10, bradysc, gaitsc, parksc, rigidsc, tremsc, cogdx, LewyBody, gpath, amyloid, tangles, caa_4gp, cogn_global))
dim(pheno_data_sub) # 

## Categorical variable: LewyBody, caa_4gp, ADD
pheno_data_sub$ADD <- cut(pheno_data_sub$cogdx, breaks = c(0, 3.1, 5.1, 7), labels = c("0", "1", NA))
pheno_data_sub$caa <- ifelse(pheno_data_sub$caa_4gp > 1, 1, 0)
pheno_data_sub$amyloid <- sqrt(pheno_data_sub$amyloid)
pheno_data_sub$tangles <- sqrt(pheno_data_sub$tangles)
pheno_data_sub$gpath <- sqrt(pheno_data_sub$gpath)
head(pheno_data_sub)


########### get data frame for DGE analysis:  ###############
dim(Protein_dt) # 480; 9870
DGE_lm_dt <- merge(data.frame(Protein_dt), pheno_data_sub, by = "projid", all = TRUE, sort = FALSE) # keep all, should keep protein data
dim(DGE_lm_dt)  


## annotation file
unique(protein_list) %>% length()
sum(protein_list %in% gene_anno$gene_name)
anno_dt <- filter(gene_anno, gene_name %in% protein_list) %>% select(c("gene_name", "start", "chr", "ENSE_ID")) %>% data.frame() %>% distinct(gene_name, .keep_all = TRUE)
row.names(anno_dt) <- anno_dt$gene_name
dim(anno_dt)
head(anno_dt)
write_tsv(anno_dt[, 1:3], file = paste0(out_dir, tissue, "_anno.tsv"), col_names = FALSE)

## Build a GEMMA‑style data frame
n_proteins = nrow(anno_dt)
Exp_dt_gemma <- data.frame(gene = anno_dt$gene_name, A = rep("A", n_proteins), T= rep("T", n_proteins), t(select(DGE_lm_dt,  anno_dt$gene_name)) )

head(Exp_dt_gemma[, 1:10])
dim(Exp_dt_gemma)

out_file_name = paste0(out_dir, tissue, "_dt.tsv")
write_tsv(Exp_dt_gemma, file = out_file_name, col_names = FALSE)
system(paste("gzip -f ", out_file_name))

## Covariate file
n_samples = nrow(DGE_lm_dt)
n_samples
# , "pmi" 
cov_list <- c("PC1", "PC2", "PC3", "PC4", "PC5", "age_death", "msex", "educ", "bmi")
cov_bim <-data.frame(intercept = rep(1, n_samples), select(DGE_lm_dt,  all_of(cov_list) ) )
write_tsv(cov_bim, file = paste0(out_dir, tissue, "_cov_bim.tsv"), col_names = FALSE)

select(DGE_lm_dt, c(projid, age_death)) 

## only consider continuous phenotypes; phenotype file
pheno_list <- c("gait_speed", "motor_dexterity", "motor_gait", "motor_handstreng", "motor10", "bradysc", "gaitsc", "parksc", "rigidsc", "tremsc", "LewyBody", "gpath", "amyloid", "tangles", "cogn_global", "ADD", "caa")
for(pheno in pheno_list){
    write_tsv(select(DGE_lm_dt, all_of(pheno)), file = paste0(out_dir, tissue, ".", pheno, ".tsv"), col_names = FALSE)
}



