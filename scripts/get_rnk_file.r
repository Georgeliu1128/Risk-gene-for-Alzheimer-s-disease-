#!/usr/bin/env Rscript
.libPaths("/nfs/yaofss2/data/shared/Software/R/x86_64-pc-linux-gnu-library/4.4")
rm(list = ls())

#load library

library(dplyr)
library(tidyr)
library(readr)
library(org.Hs.eg.db)
source("/home/zliu9/zliu9/Scripts/util_func.r")

# Set data paths
file_dir <- "/home/zliu9/zliu9/ROSMAP2025/GEMMA_Data_Spinal_cord_out/"
out_dir  <- "/home/zliu9/zliu9/ROSMAP2025/gsea_analysis_01082026/prernk_files/"

#if out_dir does not exist then creat
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE) 

# Sample lists for analysis
sig_list <- c(
  "DLPFC.gpath", "DLPFC.amyloid", "DLPFC_protein.gpath",
  "DLPFC_protein.amyloid", "DLPFC_protein.tangles",
  "DLPFC_protein.caa", "Vh_protein.LewyBody", "Vh.tremsc", "Quad.tremsc", "DLPFC.cogn_global",
  "DLPFC_protein.cogn_global", "DLPFC.ADD",
  "DLPFC_protein.tremsc", "Sma.gpath", "Sma.tangles", "Sma.cogn_global", "Sma.motor_dexterity",
  "Sma.ADD", "Vh.gpath", "Vh.tangles", "Vh.cogn_global", "Quad.gait_speed",
  "Quad.ADD","DLPFC.tangles", "DLPFC.parksc", "DLPFC.motor_dexterity",
   "DLPFC_protein.motor_gait", "DLPFC_protein.motor10",
   "Quad_protein.gait_speed", "Quad_protein.motor_handstreng",
   "Vh_protein.cogn_global"
)

rna_traits <- c(
  "DLPFC.gpath", "DLPFC.amyloid",
  "Vh.tremsc", "Quad.tremsc", "DLPFC.cogn_global",
  "DLPFC.ADD",
  "Sma.gpath", "Sma.tangles", "Sma.cogn_global",
  "Sma.motor_dexterity", "Sma.ADD", "Vh.gpath",
  "Vh.tangles", "Vh.cogn_global", "Quad.gait_speed",
  "Quad.ADD", "DLPFC.tangles", "DLPFC.parksc",
  "DLPFC.motor_dexterity"
)


for (c_name in sig_list)  {
  is_rna <- c_name %in% rna_traits

  # Read association results
  DGE_results <- read_tsv(paste0(file_dir, c_name, ".lmm.assoc.txt"))
  
  # GC correction
  gc <- round(GetGC(DGE_results$p_wald), 4)
  message(c_name, " | GC = ", gc)
  DGE_results$p_wald_correct <- CorrectGC(DGE_results$p_wald)
  
  # RNA branch: convert Ensembl → SYMBOL
if (is_rna) {
    gene_map <- AnnotationDbi::select(
      org.Hs.eg.db,
      keys = DGE_results$rs,
      keytype = "ENSEMBL",
      columns = "SYMBOL"
    )

    # Standardize column names
    colnames(gene_map) <- tolower(colnames(gene_map))

    # Detect correct columns
    ensembl_col <- grep("ensembl", colnames(gene_map), value = TRUE)
    symbol_col  <- grep("symbol",  colnames(gene_map), value = TRUE)

    # Rename safely
    gene_map <- dplyr::rename(
      gene_map,
      rs   = all_of(ensembl_col),
      gene = all_of(symbol_col)
    )

    # Select significant genes
    sig_selected <- DGE_results |>
      dplyr::filter(p_wald_correct < 0.05) |>
      dplyr::arrange(beta) |>
      dplyr::select(rs, beta) |>
      dplyr::left_join(gene_map, by = "rs") |>
      dplyr::filter(!is.na(gene)) |>
      dplyr::select(gene, beta) |>
      dplyr::group_by(gene) |>
      dplyr::slice_max(abs(beta), n = 1) |>
      dplyr::ungroup()
} else {
    # Non‑RNA branch: keep rs + beta
    sig_selected <- DGE_results |>
      dplyr::filter(p_wald_correct < 0.05) |>
      dplyr::arrange(beta) |>
      dplyr::select(rs, beta)

  }
  
  # Write .rnk file
  write.table(
    sig_selected,
    file = paste0(out_dir, c_name,".rnk"),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    col.names = FALSE
  )
}

