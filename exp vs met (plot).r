library(dplyr)
library(tidyr)
library(ggplot2)

plot_expr_meth <- function(data, gene_name, project_name = "TCGA-BRCA") {
  
  # ------------------------------------------------------------
  # Harmonize to sample level
  # ------------------------------------------------------------
  gene_long <- data |>
    filter(gene == gene_name) |>
    mutate(sample_id = substr(barcode, 1, 16)) |>
    select(sample_id, omic, value)
  
  # ------------------------------------------------------------
  # Collapse duplicates per sample × omic
  # ------------------------------------------------------------
  gene_collapsed <- gene_long |>
    group_by(sample_id, omic) |>
    summarize(value = mean(value, na.rm = TRUE), .groups = "drop")
  
  # ------------------------------------------------------------
  # Pivot to wide format
  # ------------------------------------------------------------
  gene_wide <- gene_collapsed |>
    pivot_wider(names_from = omic, values_from = value)
  
  # ------------------------------------------------------------
  # Keep paired samples only
  # ------------------------------------------------------------
  gene_paired <- gene_wide |>
    filter(complete.cases(RNAseq, Methylation))
  
  # If too few paired samples, stop gracefully
  if (nrow(gene_paired) < 10) {
    return(list(
      gene = gene_name,
      n_paired = nrow(gene_paired),
      correlation = NA_real_,
      plot = NULL
    ))
  }
  
  # ------------------------------------------------------------
  # Correlation
  # ------------------------------------------------------------
  cor_val <- cor(
    gene_paired$RNAseq,
    gene_paired$Methylation,
    method = "spearman"
  )
  
  # ------------------------------------------------------------
  # Plot
  # ------------------------------------------------------------
  p <- ggplot(gene_paired, aes(x = Methylation, y = RNAseq)) +
    geom_point(alpha = 0.7) +
    geom_smooth(method = "lm", se = FALSE) +
    labs(
      title = paste0(gene_name, ": expression vs methylation"),
      subtitle = paste0(
        project_name,
        " | Spearman ρ = ",
        round(cor_val, 3),
        " | n = ",
        nrow(gene_paired)
      ),
      x = "methylation (mean beta)",
      y = "RNA-seq expression (log2 CPM)"
    ) +
    theme_bw()
  
  list(
    gene = gene_name,
    n_paired = nrow(gene_paired),
    correlation = cor_val,
    plot = p
  )
}

gene_results <- lapply(genes, function(g) {
  plot_expr_meth(final_results, g)
})

names(gene_results) <- genes

gene_results[["APC"]]$plot
gene_results[["YTHDF1"]]$plot
gene_results[["METTL14"]]$plot
gene_results[["LILRB4"]]$plot
