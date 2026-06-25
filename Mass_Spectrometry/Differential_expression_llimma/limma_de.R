#!/usr/bin/env Rscript

# Rscript limma_massspec.R /path/to/file.csv /path/to/outdir
suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
  library(limma)
  library(ggrepel)
})

# --------------------------
# Args
# --------------------------
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  cat("Usage:\n  Rscript limma_massspec.R <input_csv_or_tsv> <output_dir>\n\n")
  quit(status = 1)
}

infile <- args[1]
outdir <- args[2]

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# --------------------------
# Helpers
# --------------------------
sanitize_tag <- function(x) {
  x <- gsub("[^A-Za-z0-9._-]+", "_", x)
  x <- gsub("_+", "_", x)
  x <- gsub("^_|_$", "", x)
  x
}

save_base_plot <- function(filename, expr) {
  png(filename, width = 1800, height = 900, res = 200)
  on.exit(dev.off(), add = TRUE)
  expr
}

# --------------------------
# Read data
# --------------------------
dt <- fread(infile)

# Detect gene column (default to first column)
gene_col <- which(tolower(names(dt)) %in% c("genes", "gene", "symbol", "geneid", "gene_id"))
gene_col <- if (length(gene_col) > 0) gene_col[1] else 1

gene_col_name <- names(dt)[gene_col]
sample_cols <- setdiff(names(dt), gene_col_name)

if (length(sample_cols) < 2) {
  stop("Could not find >=2 sample columns. Your file should have 1 gene column + sample columns.")
}

# --------------------------
# De-duplicate genes
# --------------------------
dup_tbl <- dt %>%
  count(.data[[gene_col_name]], sort = TRUE) %>%
  filter(n > 1)

genes_removed <- dup_tbl[[gene_col_name]]

dup_report <- sprintf(
  "The following %d genes that are repeated are removed: (%s).\n",
  length(genes_removed),
  paste(genes_removed, collapse = ", ")
)
cat(dup_report)

dt_dedup <- dt %>% distinct(.data[[gene_col_name]], .keep_all = TRUE)

# # --------------------------
# # Remove rows with zero variance across sample columns
# # --------------------------
# vals <- as.matrix(sapply(dt_dedup[, ..sample_cols, drop = FALSE], as.numeric))
# 
# rv <- apply(vals, 1, var, na.rm = TRUE)            # row-wise variance  :contentReference[oaicite:0]{index=0}
# keep <- is.finite(rv) & rv > 0                     # drop constant rows (and all-NA rows)
# 
# cat(sprintf("Removed %d zero-variance rows; keeping %d\n", sum(!keep), sum(keep)))
# 
# dt_dedup <- dt_dedup[keep, ]

# idx <- 2:5
# idx <- setdiff(idx, gene_col)  # gene_col is the index of gene_col_name in dt/dt_dedup
# 
# dt_dedup[, idx := lapply(.SD, function(x) {
#   x <- as.numeric(x)
#   x[x < 0] <- 0
#   x
# }), .SDcols = idx]

# --------------------------
# Build matrix X
# --------------------------
genes <- dt_dedup[[gene_col_name]]

X_df <- dt_dedup[, ..sample_cols]
X <- as.matrix(sapply(X_df, as.numeric))
rownames(X) <- genes

# --------------------------
# Infer groups (control vs experiment)
# --------------------------
nm <- colnames(X)

# is_control <- grepl("control|ctrl|igg|input", nm, ignore.case = TRUE)
# is_experiment <- grepl("experiment|exp|treat|tert|case", nm, ignore.case = TRUE)
# 
# control_cols <- nm[is_control & !is_experiment]
# experiment_cols <- nm[is_experiment & !is_control]

# If regex inference fails, fall back to "first half control, second half experiment"
# if (length(control_cols) == 0 || length(experiment_cols) == 0) {
  n <- length(nm)
  split_idx <- floor(n / 2)
  control_cols <- nm[seq_len(split_idx)]
  experiment_cols <- nm[(split_idx + 1):n]
# }

# Final safety checks
if (length(control_cols) < 1 || length(experiment_cols) < 1) {
  stop("Could not assign control/experiment columns. Rename columns or ensure >=2 samples.")
}

group <- ifelse(nm %in% control_cols, "control", "experiment")
group <- factor(group, levels = c("control", "experiment"))

design <- model.matrix(~ 0 + group)
colnames(design) <- levels(group)

# Tag used in ALL filenames (as requested)
tag <- paste0(
  "experiment_vs_control__experiment_",
  paste(experiment_cols, collapse = "_"),
  "__control_",
  paste(control_cols, collapse = "_")
)
tag <- sanitize_tag(tag)

# Save grouping info
grouping_lines <- c(
  paste0("Input: ", infile),
  paste0("Output: ", outdir),
  paste0("Control columns (", length(control_cols), "): ", paste(control_cols, collapse = ", ")),
  paste0("Experiment columns (", length(experiment_cols), "): ", paste(experiment_cols, collapse = ", ")),
  paste0("Design matrix columns: ", paste(colnames(design), collapse = ", "))
)

cat(paste(grouping_lines, collapse = "\n"), "\n", sep = "")

# --------------------------
# QC plots (boxplots + densities) -- saved with column names
# --------------------------
X_raw <- X
X_norm <- sweep(X_raw, 2, apply(X_raw, 2, median, na.rm = TRUE), FUN = "-")

save_base_plot(file.path(outdir, paste0("qc_boxplot_median_centering_", tag, ".png")), {
  par(mfrow = c(1, 2))
  boxplot(X_raw, las = 2, main = "Before", ylab = "value")
  boxplot(X_norm, las = 2, main = "After hypothetical median-centering", ylab = "value")
  par(mfrow = c(1, 1))
})

save_base_plot(file.path(outdir, paste0("qc_densities_median_centering_", tag, ".png")), {
  par(mfrow = c(1, 2))
  limma::plotDensities(X_raw, main = "Densities: before")
  limma::plotDensities(X_norm, main = "Densities: after hypothetical median-centering")
  par(mfrow = c(1, 1))
})

# --------------------------
# LIMMA differential analysis
# --------------------------
fit <- lmFit(X, design)
contr <- makeContrasts(experiment_vs_control = experiment - control, levels = design)
fit2 <- contrasts.fit(fit, contr)
fit2 <- eBayes(fit2)

res <- topTable(fit2, coef = "experiment_vs_control", number = Inf, sort.by = "P")

# Ensure gene column exists + make volcano-ready columns
res <- res %>%
  mutate(
    Genes = rownames(.),
    adj.P.Val = pmax(adj.P.Val, .Machine$double.xmin),
    neglog10FDR = -log10(adj.P.Val),
    status = case_when(
      adj.P.Val <= 0.05 & logFC >=  1 ~ "Up",
      adj.P.Val <= 0.05 & logFC <= -1 ~ "Down",
      TRUE ~ "NS"
    )
  )

# Save results table
fwrite(
  as.data.table(res),
  file.path(outdir, paste0("limma_results_", tag, ".tsv")),
  sep = "\t"
)

res_filtered <- res %>% 
  filter(logFC >= 1 ) %>% 
  filter(P.Value <= 0.05) %>% 
  select(Genes,logFC,P.Value) %>% 
  arrange(desc(logFC))

fwrite(
  as.data.table(res_filtered),
  file.path(outdir, paste0("limma_results_", tag, ".significant.tsv")),
  sep = "\t"
)

# --------------------------
# Volcano plot (saved) + 3-line title
# --------------------------
top_labs <- res %>%
  filter(status != "NS") %>%
  arrange(adj.P.Val) %>%
  slice_head(n = 20)

vol_title <- paste0(
  "experiment vs control\n",
  "experiment: ", paste(experiment_cols, collapse = ", "), "\n",
  "control: ", paste(control_cols, collapse = ", ")
)

p_volcano <- ggplot(res, aes(x = logFC, y = neglog10FDR, color = status)) +
  geom_point(alpha = 0.6, size = 1.3) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  ggrepel::geom_text_repel(
    data = top_labs,
    aes(label = Genes),
    size = 3,
    max.overlaps = Inf
  ) +
  theme_bw() +
  labs(
    title = vol_title,
    x = "log2 Fold Change (logFC)",
    y = "-log10(adj.P.Val)",
    color = ""
  )

ggsave(
  filename = file.path(outdir, paste0("volcano_", tag, ".png")),
  plot = p_volcano,
  width = 9, height = 6, dpi = 300
)

cat("\nDone.\n")
cat("Results saved to:\n  ", outdir, "\n \n", sep = "")




















# dup_tbl <- Genewisecounts %>%
#   count(Genes, sort = TRUE) %>%
#   filter(n > 1)
# 
# genes_removed <- dup_tbl$Genes
# 
# cat(
#   sprintf(
#     "The following %d genes that are repeated are removed: (%s).\n",
#     length(genes_removed),
#     paste(genes_removed, collapse = ", ")
#   )
# )
# 
# Genewisecounts_dedup <- Genewisecounts 
#   %>% distinct(Genes, .keep_all = TRUE)
# 
# 
# # STARTING LIMMA
# genes <- Genewisecounts[[1]]
# X <- as.matrix(Genewisecounts[, 2:5, with = FALSE])
# X <- apply(X, 2, as.numeric)
# rownames(X) <- genes
# 
# # checking if normalization is needed 
# X_raw <- X
# X_norm <- sweep(X_raw, 2, apply(X_raw, 2, median, na.rm = TRUE), FUN = "-")
# 
# par(mfrow = c(1,2))
# boxplot(X_raw,  las=2, main="Before", ylab="value")
# boxplot(X_norm, las=2, main="After hypothetical median-centering", ylab="value")
# par(mfrow = c(1,1))
# 
# 
# par(mfrow = c(1,2))
# plotDensities(X_raw,  main="Densities: before",)
# plotDensities(X_norm, main="Densities: after hypothetical median-centering ")
# par(mfrow = c(1,1))
# 
# 
# group <- factor(c(rep("control", 2), rep("experiment", 2)),
#                 levels=c("control","experiment"))
# 
# fit <- lmFit(X, design)
# contr <- makeContrasts(experiment_vs_control = experiment - control, levels=design)
# fit2 <- contrasts.fit(fit, contr)
# fit2 <- eBayes(fit2)
# 
# res <- topTable(fit2, coef="experiment_vs_control", number=Inf, sort.by="P")
# 
# res <- res %>%
#   mutate(
#     Genes = if ("Genes" %in% names(.)) Genes else rownames(.),
#     # avoid Inf if adj.P.Val is 0
#     neglog10FDR = -log10(adj.P.Val),
#     status = case_when(
#       adj.P.Val <= 0.05 & logFC >=  1 ~ "Up",
#       adj.P.Val <= 0.05 & logFC <= -1 ~ "Down",
#       TRUE ~ "NS"
#     )
#   )
# 
# top_labs <- res %>%
#   filter(status != "NS") %>%
#   arrange(adj.P.Val) %>%
#   slice_head(n = 20)
# 
# ggplot(res, aes(x = logFC, y = neglog10FDR, color = status)) +
#   geom_point(alpha = 0.6, size = 1.3) +
#   geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
#   geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
#   ggrepel::geom_text_repel(
#     data = top_labs,
#     aes(label = Genes),
#     size = 3,
#     max.overlaps = Inf
#   ) +
#   theme_bw() +
#   labs(
#     x = "log2 Fold Change (logFC)",
#     y = "-log10(adj.P.Val)",
#     color = ""
#   )
