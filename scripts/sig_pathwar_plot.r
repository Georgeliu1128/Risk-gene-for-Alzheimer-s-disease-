#!/usr/bin/env Rscript
.libPaths("/nfs/yaofss2/data/shared/Software/R/x86_64-pc-linux-gnu-library/4.4")
rm(list = ls())

library(ggplot2)
library(dplyr)
library(tidyr)
library(readr)
library(stringr)
library(forcats)
library(tibble)
library(purrr)
library(fs)

base_dir <- "/home/zliu9/zliu9/ROSMAP2025/gsea_analysis_01082026/gsea_output/fourth_run"
out_dir_pdf <- "/home/zliu9/zliu9/ROSMAP2025/gsea_analysis_01082026/visualization/PDF/fourth_run"
out_dir_npg <- "/home/zliu9/zliu9/ROSMAP2025/gsea_analysis_01082026/visualization/npg/fourth_run"
dir.create(out_dir_pdf, showWarnings = FALSE, recursive = TRUE)
dir.create(out_dir_npg, showWarnings = FALSE, recursive = TRUE)

AD_trait <- c(
  "gpath",
  "amyloid",
  "tangles",
  "cogn_global",
  "ADD"
)

motor_trait <- c(
  "gait_speed",
  "motor_dexterity",
  "motor_gait",
  "motor_handstreng",
  "motor10",
  "bradysc",
  "gaitsc",
  "parksc",
  "rigidsc",
  "tremsc",
  "LewyBody"
)


# find all matching files recursively
files <- list.files(
  path = base_dir,
  pattern = "gsea_report_for_na_.*\\.tsv$",
  full.names = TRUE,
  recursive = TRUE
)

# read + annotate each file
results <- map_df(files, function(f) {

  df <- read_tsv(f, show_col_types = FALSE)
  
  num_cols <- c("SIZE", "ES", "NES", "NOM p-val", "FDR q-val", 
              "FWER p-val", "RANK AT MAX")

  df <- df %>%
    mutate(across(all_of(num_cols), ~ suppressWarnings(as.numeric(.x))))



  # extract relative directory path
  rel_path <- strsplit(sub(base_dir, "", dirname(f)), "/")[[1]]
  clean_path <- rel_path[rel_path != ""]

  lib_name <- clean_path[1]
  second_subdir <- if (length(clean_path) >= 2) clean_path[2] else NA

  parts <- unlist(strsplit(second_subdir, "\\."))
  tissue <- parts[1]
  trait  <- parts[2]

  # detect pos or neg from filename
  direction <- case_when(
    str_detect(basename(f), "pos") ~ "pos",
    str_detect(basename(f), "neg") ~ "neg",
    TRUE ~ "unknown"
  )

  df %>%
    mutate(
      file = f,
      library = lib_name,
      tissue = tissue,
      trait = trait,
      direction = direction
    )
})

# output file paths
out_results <- file.path(base_dir, "all_pathways.tsv")

# save full results
write_tsv(results, out_results)

cat("Saved all pathway results to:\n", out_results, "\n")

# find the significant pathways and summarize
sig_path <- results |>
  filter(`FDR q-val` < 0.05) |>
  group_by(library, tissue, trait, direction) |>
  summarise(n_sets = n(), .groups = "drop")

# output file paths
out_sig <- file.path(base_dir, "significant_pathway_summary_fourth.tsv")
write_tsv(sig_path, out_sig)

#shared AD_trait and motor_trait
df <- results |>
  mutate(
    pathology = case_when(
      trait %in% AD_trait ~ "AD_trait",
      trait %in% motor_trait ~ "motor_trait",
      TRUE ~ "other"
    )
  ) |>
  filter(`FDR q-val` < 0.05)

shared_pathways <- df |>
  filter(pathology != "other") |>
  group_by(library, tissue, NAME) |>
  summarise(n_path = n_distinct(pathology), .groups = "drop") |>
  filter(n_path >= 2)

df_shared <- df |>
  semi_join(shared_pathways, by = c("library", "tissue", "NAME"))

distinct_pathways <- df_shared |>
  distinct(library, tissue, NAME)  # keep only unique pathways per group

write_tsv(distinct_pathways, file.path(base_dir, "distinct_pathway_both2.tsv"))  
  
write_tsv(df_shared, file.path(base_dir, "pathology_pathway_both2.tsv"))

# list of tissues
tissues <- unique(df_shared$tissue)

# loop over tissues
for (t in tissues) {

  df_t <- df_shared %>%
    filter(tissue == t) %>%
    mutate(
      NAME_trait = paste0(NAME, " (", trait, ")"),
      NAME_trait = fct_reorder(NAME_trait, NES)
    )

  df_t$NAME_trait <- as.character(df_t$NAME_trait)

  bar_plot <- ggplot(df_t, aes(
    x = NES,
    y = NAME_trait,
    fill = trait
    )) +
    geom_col() +
    facet_grid(. ~ library, scales = "free_y", space = "free_y") +
    scale_fill_brewer(palette = "Set3") +
    labs(
      title = paste("Pathway NES for Tissue:", t),
      x = "NES",
      y = "Pathway (NAME × Trait)",
      fill = "trait"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 14, face = "bold"),
      axis.title.x = element_text(size = 12, face = "bold"),
      axis.title.y = element_text(size = 12, face = "bold"),
      axis.text.x  = element_text(size = 12, face = "bold"),
      axis.text.y  = element_text(size = 12, face = "bold"),
      strip.text = element_text(size = 12, face = "bold")
    )

  # dynamic sizing
  max_label_len <- max(nchar(df_t$NAME_trait))
  dynamic_width  <- 6 + 0.12 * max_label_len
  dynamic_height <- 0.4 * nrow(df_t)

  # filenames
  prefix <- paste0("NES_trait_", t)

  # save PDF
  ggsave(
    file.path(out_dir_pdf, paste0(prefix, ".pdf")),
    bar_plot,
    width = dynamic_width,
    height = dynamic_height,
    dpi = 600
  )

  # save PNG
  ggsave(
    file.path(out_dir_npg, paste0(prefix, ".png")),
    bar_plot,
    width = dynamic_width,
    height = dynamic_height,
    dpi = 600,
    bg = "white"
  )
}