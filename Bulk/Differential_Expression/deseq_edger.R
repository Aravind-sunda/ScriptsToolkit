suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
  library(DESeq2)
  library(edgeR)
  library(gplots)
  library(RColorBrewer)
  library(EnhancedVolcano)
})

# --------------------------
# Args & Setup
# --------------------------
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  cat("Usage:\n  Rscript de_pipeline.R <counts.csv/tsv> <metadata.csv/tsv> <output_dir>\n\n")
  quit(status = 1)
}

counts_file <- args[1]
meta_file   <- args[2]
outdir      <- args[3]

outdir_deseq <- paste0(outdir,"/DESeq2")
outdir_edger <-  paste0(outdir,"/edgeR")

dir.create(outdir_deseq, recursive = TRUE, showWarnings = FALSE)
dir.create(outdir_edger, recursive = TRUE, showWarnings = FALSE)


# --------------------------
# Read Data
# --------------------------
cat("[1/6] Reading data...\n")
counts_df <- fread(counts_file)
meta_df   <- fread(meta_file)

# Assume col 1 of counts is GeneID, col 1 of meta is SampleID, col 2 of meta is Condition
gene_col  <- names(counts_df)[1]
genes     <- counts_df[[gene_col]]

counts_mat <- as.matrix(counts_df[, -1, with = FALSE])
rownames(counts_mat) <- genes

# Prepare metadata
sample_col <- names(meta_df)[1]
cond_col   <- names(meta_df)[2]

meta <- data.frame(
  Condition = as.factor(meta_df[[cond_col]]),
  row.names = meta_df[[sample_col]]
)

# Ensure matrices align
common_samples <- intersect(colnames(counts_mat), rownames(meta))
if (length(common_samples) < 2) stop("Not enough matching samples between counts and metadata.")

counts_mat <- counts_mat[, common_samples]
counts_mat[is.na(counts_mat)] <- 0
cat("NAs successfully replaced with zeros. Matrix is ready for DESeq2/edgeR.\n")
meta       <- meta[common_samples, , drop = FALSE]

# --------------------------
# Filter Lowly Expressed Genes
# --------------------------
# 1. Evaluate counts > 10, ignoring NAs so they don't break the sum
keep_genes <- rowSums(counts_mat > 10, na.rm = TRUE) >= 1

# 2. Safely convert any lingering NAs in the logical vector to FALSE. rows sums of NA will inject NAs into the whole row. 
# (Prevents R from injecting blank NA rows when we subset)
keep_genes[is.na(keep_genes)] <- FALSE

# 3. Apply the filter to the matrix
counts_mat <- counts_mat[keep_genes, ,drop = FALSE]

# 4. sum() will now work perfectly
cat(sprintf("Filtering of low expression and NA value rows complete: kept %d genes out of %d original genes.\n", 
            sum(keep_genes), length(keep_genes)))

# --------------------------
# Dynamic Setup (Pairs & Colors)
# --------------------------
conditions <- levels(meta$Condition)
if (length(conditions) < 2) stop("You need at least 2 distinct sample types for differential expression.")

# Generate all pairwise combinations
pairs <- combn(conditions, 2, simplify = FALSE)
cat(sprintf("Found %d conditions. Running %d pairwise comparisons.\n", length(conditions), length(pairs)))

# Setup dynamic colors for heatmaps based on the number of conditions
num_conds <- length(conditions)
pal_colors <- if(num_conds <= 9) brewer.pal(max(3, num_conds), "Set1") else colorRampPalette(brewer.pal(9, "Set1"))(num_conds)
cond_colors <- setNames(pal_colors[1:num_conds], conditions)
col_side_colors <- cond_colors[as.character(meta$Condition)]


# --------------------------
# Helpers for Plotting
# --------------------------
plot_volcano <- function(res_df, title, filename, lfc_col="log2FoldChange", pval_col="padj") {
  p <- EnhancedVolcano(res_df,
    lab = rownames(res_df),
    x = lfc_col,
    y = pval_col,
    title = title,
    pCutoff = 0.05,
    FCcutoff = 1,
    pointSize = 2.0,
    labSize = 4.0
  )
  ggsave(filename, plot = p, width = 8, height = 8, dpi = 300)
}

plot_heatmap <- function(mat, title, filename) {
  png(filename, width = 1000, height = 1000, res = 150)
  heatmap.2(mat,
            scale = "row",
            Rowv = TRUE, Colv = TRUE,
            trace = "none",
            dendrogram = "both",
            srtCol = 45,
            margins = c(12, 10),
            main = title,
            col = colorRampPalette(rev(brewer.pal(9, "RdBu")))(255),
            ColSideColors = col_side_colors)
  dev.off()
}


# --------------------------
# DESeq2 Pipeline
# --------------------------
cat("[2/6] Running DESeq2 overall model...\n")
dds <- DESeqDataSetFromMatrix(countData = counts_mat, colData = meta, design = ~ Condition)
dds <- DESeq(dds)

# Get Variance Stabilized counts for heatmaps
vsd <- DESeq2::vst(dds,blind = FALSE)
vsd_mat <- assay(vsd)

cat("[3/6] Extracting DESeq2 pairwise results...\n")
for (p in pairs) {
  grp1 <- p[1]
  grp2 <- p[2]
  comp_name <- paste0(grp1, "_vs_", grp2)
  cat("  ->", comp_name, "\n")

  # Note: contrast is c("factorName", "numerator", "denominator")
  res_deseq <- results(dds, contrast = c("Condition", grp1, grp2))
  res_df <- as.data.frame(res_deseq)

  # Save raw results
  fwrite(as.data.table(res_df, keep.rownames = "Gene"),
         file.path(outdir_deseq, paste0("DESeq2_", comp_name, ".csv")))

  # Volcano Plot
  plot_volcano(res_df, paste("DESeq2:", grp1, "vs", grp2),
               file.path(outdir_deseq, paste0("DESeq2_volcano_", comp_name, ".png")))

  # Heatmap of top 50 DE genes (by lowest padj)
  top_genes <- rownames(res_df[order(res_df$padj), ])[1:50]
  top_genes <- top_genes[!is.na(top_genes)]
  if(length(top_genes) > 1) {
    plot_heatmap(vsd_mat[top_genes, ], paste("DESeq2 Top 50:", comp_name),
                 file.path(outdir_deseq, paste0("DESeq2_heatmap_", comp_name, ".png")))
  }
}


# --------------------------
# edgeR Pipeline
# --------------------------
cat("[4/6] Running edgeR overall model...\n")
y <- DGEList(counts = counts_mat, group = meta$Condition)

# Filter lowly expressed genes
keep <- filterByExpr(y)
y <- y[keep, , keep.lib.sizes = FALSE]

# Normalize and estimate dispersion
y <- calcNormFactors(y)
design <- model.matrix(~ 0 + Condition, data = meta)
colnames(design) <- levels(meta$Condition)

y <- estimateDisp(y, design)
fit <- glmQLFit(y, design)

# Get logCPM values for heatmaps
logCPM_mat <- cpm(y, prior.count = 2, log = TRUE)

cat("[5/6] Extracting edgeR pairwise results...\n")
for (p in pairs) {
  grp1 <- p[1]
  grp2 <- p[2]
  comp_name <- paste0(grp1, "_vs_", grp2)
  cat("  ->", comp_name, "\n")

  # Make contrast array dynamically
  contr_vector <- rep(0, ncol(design))
  names(contr_vector) <- colnames(design)
  contr_vector[grp1] <- 1
  contr_vector[grp2] <- -1

  qlf <- glmQLFTest(fit, contrast = contr_vector)
  res_edger <- topTags(qlf, n = Inf)$table

  # Save raw results
  fwrite(as.data.table(res_edger, keep.rownames = "Gene"),
         file.path(outdir_edger, paste0("edgeR_", comp_name, ".csv")))

  # Volcano Plot (edgeR uses logFC and FDR)
  plot_volcano(res_edger, paste("edgeR:", grp1, "vs", grp2),
               file.path(outdir_edger, paste0("edgeR_volcano_", comp_name, ".png")),
               lfc_col = "logFC", pval_col = "FDR")

  # Heatmap of top 50 DE genes (by lowest FDR)
  top_genes_edger <- rownames(res_edger[order(res_edger$FDR), ])[1:50]
  top_genes_edger <- top_genes_edger[!is.na(top_genes_edger)]
  if(length(top_genes_edger) > 1) {
    plot_heatmap(logCPM_mat[top_genes_edger, ], paste("edgeR Top 50:", comp_name),
                 file.path(outdir_edger, paste0("edgeR_heatmap_", comp_name, ".png")))
  }
}

cat("[6/6] Pipeline complete! Results saved in:", outdir, "\n")


# --------------------------
# Compare DESeq2 and edgeR Overlaps
# --------------------------
library(ggVennDiagram)

cat("[7/7] Generating overlap Venn diagrams...\n")

for (p in pairs) {
  grp1 <- p[1]
  grp2 <- p[2]
  comp_name <- paste0(grp1, "_vs_", grp2)
  cat("  ->", comp_name, "\n")
  
  # Read the saved full results back in
  deseq_file <- file.path(outdir_deseq, paste0("DESeq2_", comp_name, ".csv"))
  edger_file <- file.path(outdir_edger, paste0("edgeR_", comp_name, ".csv"))
  
  if (file.exists(deseq_file) && file.exists(edger_file)) {
    deseq_res <- fread(deseq_file)
    edger_res <- fread(edger_file)
    
    # Filter for significant genes (Adjust thresholds here if needed)
    # Using padj/FDR < 0.05 and absolute logFC > 1
    sig_deseq <- deseq_res[padj < 0.05 & abs(log2FoldChange) > 1, Gene]
    sig_edger <- edger_res[FDR < 0.05 & abs(logFC) > 1, Gene]
    
    # Create a named list of the significant genes
    gene_lists <- list(
      DESeq2 = sig_deseq,
      edgeR  = sig_edger
    )
    
    # Generate the Venn diagram
    p_venn <- ggVennDiagram(gene_lists, label_alpha = 0) +
      ggtitle(paste("Significant DE Genes:", comp_name)) +
      theme(plot.title = element_text(hjust = 0.5, face = "bold")) +
      scale_fill_gradient(low = "white", high = "dodgerblue") +
    scale_x_continuous(expand = expansion(mult = 0.2)) + # Adds 20% padding left/right
      scale_y_continuous(expand = expansion(mult = 0.2))   # Adds 20% padding top/bottom
    
    # Save the plot in the main output directory
    ggsave(
      filename = file.path(outdir, paste0("Venn_", comp_name, ".png")), 
      plot = p_venn, 
      width = 8, height = 5, dpi = 300
    )
  } else {
    cat("     Skipping: Results files not found for", comp_name, "\n")
  }
}

cat("Overlap analysis complete!\n")

