#!/usr/bin/env Rscript
# =============================================================================
# go_enrichment.R -- general-purpose Gene Ontology enrichment
#
# Design rationale (see README.md for citations):
#   * ORA engine  = clusterProfiler::enrichGO  (Frontiers 2026 benchmark: one of
#     four tools -- with DAVID, Enrichr, topGO -- that consistently rank the
#     expected target ontologies near the top and return small, deep, i.e.
#     specific, GO terms; the only one of the four that is fully scriptable,
#     supports custom backgrounds, and is not a web service.)
#   * Optional topGO weight01 cross-check (most conservative scriptable tool;
#     prunes redundant parent terms using the GO graph).
#   * Optional FCS/GSEA via clusterProfiler::gseGO when a ranked statistic is
#     available (Wijesooriya/Ziemann: prefer FCS over ORA when you have ranks).
#   * Background is a first-class argument and the script refuses to guess
#     quietly -- inappropriate backgrounds are the single most common error in
#     the published ORA literature.
#   * Up- and down-regulated genes are tested separately by default.
#   * Full result tables are written unfiltered (raw + adjusted p) so nothing
#     is silently lost, plus annotation size and GO depth per term as
#     interpretability metrics.
#
# Usage examples at the bottom of this file (--help also prints them).
# =============================================================================

suppressWarnings(suppressPackageStartupMessages({
  library(optparse)
}))

# ------------------------------------------------------------------ CLI ------

option_list <- list(
  make_option("--genes", type = "character", default = NULL,
              help = "Query gene list file (one gene per line, or first column of a table). Comma-separate several files to compare lists Metascape-style. Mutually exclusive with --de-table."),
  make_option("--labels", type = "character", default = NULL,
              help = "Comma-separated names for the lists in --genes [default: file basenames]."),
  make_option("--background", type = "character", default = NULL,
              help = "File with the background/universe gene list. STRONGLY recommended: all genes detected/tested in the experiment, NOT the whole genome."),
  make_option("--de-table", type = "character", default = NULL,
              help = "Differential-expression table (csv/tsv). Query and background are derived from it automatically."),
  make_option("--gene-col", type = "character", default = NULL,
              help = "Gene ID column in --de-table [default: first column]."),
  make_option("--lfc-col", type = "character", default = "log2FoldChange",
              help = "Log fold-change column in --de-table [default %default]."),
  make_option("--padj-col", type = "character", default = "padj",
              help = "Adjusted p-value column in --de-table [default %default]."),
  make_option("--pval-col", type = "character", default = NULL,
              help = "Raw p-value column [default: auto-detect pvalue / PValue / P.Value]. Used by --rank-metric signed_logp."),
  make_option("--stat-col", type = "character", default = NULL,
              help = "Test-statistic column in --de-table used for GSEA ranking. Signed statistics (DESeq2 'stat', limma 't') are used as-is; unsigned ones (edgeR 'F' or 'LR') are converted to sign(logFC)*sqrt(stat)."),
  make_option("--lfc-cutoff", type = "double", default = 1,
              help = "|log2FC| threshold for calling a gene significant [default %default]."),
  make_option("--padj-cutoff", type = "double", default = 0.05,
              help = "Adjusted p threshold for calling a gene significant [default %default]."),
  make_option("--split", type = "character", default = "both",
              help = "Which DE genes to test: up | down | both (separately) | all (pooled) [default %default]."),

  make_option("--organism", type = "character", default = "human",
              help = "human | mouse | rat | zebrafish | fly | worm | yeast, or an OrgDb package name [default %default]."),
  make_option("--id-type", type = "character", default = "auto",
              help = "auto | SYMBOL | ENSEMBL | ENTREZID | REFSEQ | UNIPROT [default %default]."),
  make_option("--ont", type = "character", default = "BP",
              help = "GO domain: BP | MF | CC | ALL [default %default]."),

  make_option("--min-size", type = "integer", default = 10,
              help = "Minimum genes annotated to a GO term (in the background) [default %default]."),
  make_option("--max-size", type = "integer", default = 500,
              help = "Maximum genes annotated to a GO term. Caps uninformative root-level terms [default %default]."),
  make_option("--min-count", type = "integer", default = 3,
              help = "Minimum query genes overlapping a term for it to be reported as significant [default %default]."),
  make_option("--padj", type = "double", default = 0.05,
              help = "Adjusted p cutoff for the *filtered* result table [default %default]."),
  make_option("--p-method", type = "character", default = "BH",
              help = "Multiple-testing correction: BH | BY | holm | bonferroni [default %default]."),

  make_option("--simplify", action = "store_true", default = FALSE,
              help = "Also emit a redundancy-reduced table (clusterProfiler::simplify, Wang semantic similarity)."),
  make_option("--simplify-cutoff", type = "double", default = 0.7,
              help = "Semantic-similarity cutoff for --simplify [default %default]."),
  make_option("--topgo", action = "store_true", default = FALSE,
              help = "Run a topGO weight01 cross-check and merge its p-values into the output."),
  make_option("--gsea", action = "store_true", default = FALSE,
              help = "Also run GSEA (gseGO) on the ranked gene list. Requires --de-table."),
  make_option("--rank-metric", type = "character", default = "auto",
              help = paste("GSEA ranking: auto | stat | signed_sqrt | signed_logp | lfc.",
                           "'auto' prefers a signed moderated statistic (DESeq2 'stat', limma 't'),",
                           "then sign(logFC)*sqrt(F or LR) for edgeR, then signed -log10(p)",
                           "[default %default].")),

  make_option("--heatmap", action = "store_true", default = FALSE,
              help = "When >1 gene set is analysed, build a Metascape-style cross-list heatmap of non-redundant term clusters."),
  make_option("--kappa", type = "double", default = 0.3,
              help = "Kappa similarity cutoff for grouping redundant terms in the heatmap [default %default]."),
  make_option("--heatmap-top", type = "integer", default = 20,
              help = "Number of term clusters shown in the heatmap [default %default]."),

  make_option("--outdir", type = "character", default = "go_results",
              help = "Output directory [default %default]."),
  make_option("--prefix", type = "character", default = "GO",
              help = "Prefix for output files [default %default]."),
  make_option("--top-n", type = "integer", default = 25,
              help = "Number of terms shown in plots [default %default]."),
  make_option("--seed", type = "integer", default = 42,
              help = "RNG seed (GSEA permutations) [default %default]."),
  make_option("--allow-genome-background", action = "store_true", default = FALSE,
              help = "Permit running with the whole-annotation background. Off by default on purpose.")
)

opt <- parse_args(OptionParser(
  option_list = option_list,
  usage = "%prog [--genes LIST --background LIST | --de-table TABLE] [options]"
))

# optparse's naming of multi-word flags varies by version (gene-col / gene_col /
# gene.col), so look the key up tolerantly.
g <- function(k) {
  for (kk in unique(c(k, gsub("-", "_", k), gsub("-", ".", k))))
    if (kk %in% names(opt)) return(opt[[kk]])
  NULL
}

# ------------------------------------------------------- validate inputs -----

msg  <- function(...) cat(sprintf(...), "\n", sep = "")
fail <- function(...) { cat("ERROR: ", sprintf(...), "\n", sep = ""); quit(status = 1) }

if (is.null(g("genes")) && is.null(g("de-table")))
  fail("supply either --genes (plus --background) or --de-table.")
if (!is.null(g("genes")) && !is.null(g("de-table")))
  fail("--genes and --de-table are mutually exclusive.")
if (g("gsea") && is.null(g("de-table")))
  fail("--gsea needs a ranked list, so it requires --de-table.")
if (!g("ont") %in% c("BP", "MF", "CC", "ALL"))
  fail("--ont must be BP, MF, CC or ALL.")
if (!g("split") %in% c("up", "down", "both", "all"))
  fail("--split must be up, down, both or all.")

set.seed(g("seed"))
dir.create(g("outdir"), showWarnings = FALSE, recursive = TRUE)

# ------------------------------------------------------------- packages ------

orgdb_for <- c(
  human = "org.Hs.eg.db", mouse = "org.Mm.eg.db", rat = "org.Rn.eg.db",
  zebrafish = "org.Dr.eg.db", fly = "org.Dm.eg.db", worm = "org.Ce.eg.db",
  yeast = "org.Sc.sgd.db"
)
orgdb_name <- if (g("organism") %in% names(orgdb_for)) orgdb_for[[g("organism")]] else g("organism")

need <- c("clusterProfiler", "GO.db", "ggplot2", orgdb_name)
if (g("topgo")) need <- c(need, "topGO")
missing <- need[!vapply(need, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing))
  fail("missing R packages: %s\n  BiocManager::install(c(%s))",
       paste(missing, collapse = ", "),
       paste(sprintf('"%s"', missing), collapse = ", "))

suppressWarnings(suppressPackageStartupMessages({
  library(clusterProfiler); library(GO.db); library(ggplot2)
  library(orgdb_name, character.only = TRUE)
}))
OrgDb <- get(orgdb_name)

# ------------------------------------------------------------- helpers -------

read_gene_file <- function(path) {
  if (!file.exists(path)) fail("file not found: %s", path)
  x <- readLines(path, warn = FALSE)
  x <- x[nzchar(trimws(x))]
  # take the first field if the file is actually a table
  x <- vapply(strsplit(x, "[\t,]"), `[`, character(1), 1)
  x <- trimws(gsub('"', "", x))
  # drop a header-looking first line
  if (length(x) && grepl("^(gene|gene_?id|gene_?name|symbol|id)$", x[1], ignore.case = TRUE))
    x <- x[-1]
  unique(x[nzchar(x)])
}

read_table_auto <- function(path) {
  if (!file.exists(path)) fail("file not found: %s", path)
  sep <- if (grepl("\\.csv$", path, ignore.case = TRUE)) "," else "\t"
  df  <- read.delim(path, sep = sep, header = TRUE, check.names = FALSE,
                    stringsAsFactors = FALSE)
  if (ncol(df) == 1) df <- read.delim(path, sep = ",", header = TRUE,
                                      check.names = FALSE, stringsAsFactors = FALSE)
  df
}

strip_ens_version <- function(x) sub("\\.\\d+$", "", x)

guess_id_type <- function(x) {
  x <- x[!is.na(x)][1:min(200, length(x))]
  if (mean(grepl("^ENS[A-Z]*[GT]?\\d{6,}", x)) > 0.5) return("ENSEMBL")
  if (mean(grepl("^\\d+$", x))                > 0.8) return("ENTREZID")
  if (mean(grepl("^N[MR]_|^X[MR]_", x))       > 0.5) return("REFSEQ")
  "SYMBOL"
}

# Map arbitrary IDs -> ENTREZID, reporting the mapping rate (Ziemann et al.
# flag unreported ID loss as a reproducibility problem).
to_entrez <- function(ids, from, label) {
  if (from == "ENSEMBL") ids <- strip_ens_version(ids)
  ids <- unique(ids[!is.na(ids) & nzchar(ids)])
  if (from == "ENTREZID") {
    valid <- ids[ids %in% keys(OrgDb, keytype = "ENTREZID")]
    msg("  %-12s %d input -> %d mapped (%.1f%%)", label, length(ids), length(valid),
        100 * length(valid) / max(1, length(ids)))
    return(valid)
  }
  m <- suppressWarnings(suppressMessages(
    tryCatch(bitr(ids, fromType = from, toType = "ENTREZID", OrgDb = OrgDb),
             error = function(e) NULL)))
  if (is.null(m) || !nrow(m)) fail("could not map any %s IDs of type %s.", label, from)
  out <- unique(m$ENTREZID)
  msg("  %-12s %d input -> %d mapped (%.1f%%)", label, length(ids), length(out),
      100 * length(unique(m[[from]])) / max(1, length(ids)))
  out
}

# Longest path from a GO term to the ontology root. Deeper == more specific.
go_depth <- function(ids, ont) {
  PAR <- switch(ont, BP = GO.db::GOBPPARENTS, MF = GO.db::GOMFPARENTS,
                CC = GO.db::GOCCPARENTS, NULL)
  if (is.null(PAR)) return(rep(NA_integer_, length(ids)))
  parents <- as.list(PAR)
  memo <- new.env(hash = TRUE, parent = emptyenv())
  depth_of <- function(id) {
    if (exists(id, envir = memo, inherits = FALSE)) return(get(id, envir = memo))
    p <- parents[[id]]
    p <- if (is.null(p)) character(0) else setdiff(unname(p), "all")
    d <- if (!length(p)) 0L else max(vapply(p, depth_of, integer(1))) + 1L
    assign(id, d, envir = memo)
    d
  }
  vapply(ids, function(i) if (i %in% names(parents)) depth_of(i) else NA_integer_,
         integer(1), USE.NAMES = FALSE)
}

ratio_num <- function(x) as.numeric(sub("/.*", "", x))
ratio_den <- function(x) as.numeric(sub(".*/", "", x))

# --------------------------------------------------- assemble gene sets ------

msg("\n== Input ==")
id_type <- g("id-type")
sets <- list()          # name -> character vector of query genes (original IDs)
bg_raw <- NULL
de <- NULL

if (!is.null(g("de-table"))) {
  de <- read_table_auto(g("de-table"))
  gcol <- if (!is.null(g("gene-col"))) g("gene-col") else names(de)[1]
  for (cc in c(gcol, g("lfc-col"), g("padj-col")))
    if (!cc %in% names(de)) fail("column '%s' not in --de-table (have: %s)", cc,
                                 paste(names(de), collapse = ", "))
  de <- de[!is.na(de[[gcol]]) & nzchar(de[[gcol]]), , drop = FALSE]
  de <- de[!duplicated(de[[gcol]]), , drop = FALSE]

  bg_raw <- de[[gcol]]                                  # <- correct background
  lfc  <- suppressWarnings(as.numeric(de[[g("lfc-col")]]))
  padj <- suppressWarnings(as.numeric(de[[g("padj-col")]]))
  sig  <- !is.na(padj) & padj < g("padj-cutoff") &
          !is.na(lfc)  & abs(lfc) >= g("lfc-cutoff")

  if (g("split") == "all") {
    sets[["all"]] <- de[[gcol]][sig]
  } else {
    if (g("split") %in% c("up",   "both")) sets[["up"]]   <- de[[gcol]][sig & lfc > 0]
    if (g("split") %in% c("down", "both")) sets[["down"]] <- de[[gcol]][sig & lfc < 0]
  }
  msg("DE table: %d genes tested; %d significant (|log2FC| >= %g, %s < %g)",
      nrow(de), sum(sig), g("lfc-cutoff"), g("padj-col"), g("padj-cutoff"))
} else {
  files <- trimws(strsplit(g("genes"), ",")[[1]])
  labs  <- if (!is.null(g("labels"))) trimws(strsplit(g("labels"), ",")[[1]])
           else sub("\\.[^.]*$", "", basename(files))
  if (length(labs) != length(files))
    fail("--labels has %d entries but --genes has %d files.", length(labs), length(files))
  if (anyDuplicated(labs)) fail("--labels must be unique.")
  for (i in seq_along(files)) sets[[labs[i]]] <- read_gene_file(files[i])
  msg("Query lists: %s", paste(sprintf("%s (n=%d)", labs, lengths(sets)), collapse = ", "))
  if (!is.null(g("background"))) bg_raw <- read_gene_file(g("background"))
}

if (id_type == "auto") {
  id_type <- guess_id_type(unlist(sets, use.names = FALSE))
  msg("Detected gene ID type: %s", id_type)
}

if (is.null(bg_raw)) {
  if (!g("allow-genome-background"))
    fail(paste0("no --background supplied.\n",
      "  An ORA background must be the set of genes that could plausibly have been\n",
      "  detected in your experiment (e.g. all genes with non-zero expression, or all\n",
      "  genes tested for peaks) -- using the whole annotation inflates enrichment and\n",
      "  is the most common error in the published literature.\n",
      "  Pass --background FILE, or --allow-genome-background to override deliberately."))
  msg("WARNING: using the whole %s annotation as background (--allow-genome-background).", orgdb_name)
}

msg("ID mapping:")
bg_entrez <- if (!is.null(bg_raw)) to_entrez(bg_raw, id_type, "background") else NULL
sets_entrez <- lapply(names(sets), function(n) to_entrez(sets[[n]], id_type, n))
names(sets_entrez) <- names(sets)

if (!is.null(bg_entrez)) {
  for (n in names(sets_entrez)) {
    lost <- setdiff(sets_entrez[[n]], bg_entrez)
    if (length(lost))
      msg("  note: %d '%s' gene(s) not in background; adding them to the universe.", length(lost), n)
  }
  bg_entrez <- union(bg_entrez, unlist(sets_entrez, use.names = FALSE))
}

# --------------------------------------------------------------- ORA ---------

annotate_terms <- function(res_df) {
  if (!nrow(res_df)) return(res_df)
  ont_col <- if ("ONTOLOGY" %in% names(res_df)) as.character(res_df$ONTOLOGY)
             else rep(g("ont"), nrow(res_df))
  res_df$AnnotationSize <- ratio_num(res_df$BgRatio)          # term size in background
  res_df$FoldEnrichment <- (ratio_num(res_df$GeneRatio) / ratio_den(res_df$GeneRatio)) /
                           (ratio_num(res_df$BgRatio)   / ratio_den(res_df$BgRatio))
  res_df$Depth <- NA_integer_
  for (o in unique(ont_col)) {
    idx <- which(ont_col == o)
    res_df$Depth[idx] <- go_depth(res_df$ID[idx], o)
  }
  res_df
}

run_topgo <- function(sig_ids, bg_ids, ont) {
  gl <- factor(as.integer(bg_ids %in% sig_ids), levels = c(0, 1))
  names(gl) <- bg_ids
  GOdata <- suppressMessages(new("topGOdata", ontology = ont, allGenes = gl,
                                 nodeSize = g("min-size"), annot = topGO::annFUN.org,
                                 mapping = orgdb_name, ID = "entrez"))
  res <- suppressMessages(topGO::runTest(GOdata, algorithm = "weight01", statistic = "fisher"))
  n   <- length(topGO::usedGO(GOdata))
  tab <- topGO::GenTable(GOdata, weight01 = res, topNodes = n, numChar = 1000)
  p   <- suppressWarnings(as.numeric(gsub("[<>= ]", "", tab$weight01)))
  p[is.na(p)] <- 1e-30
  data.frame(ID = tab$GO.ID, topGO_weight01_p = p,
             topGO_weight01_padj = p.adjust(p, method = g("p-method")),
             stringsAsFactors = FALSE)
}

dot_plot <- function(df, title, file) {
  df <- head(df[order(df$p.adjust), ], g("top-n"))
  if (!nrow(df)) return(invisible(NULL))
  df$Description <- factor(df$Description, levels = rev(df$Description[order(df$FoldEnrichment)]))
  p <- ggplot(df, aes(x = FoldEnrichment, y = Description,
                      size = Count, colour = -log10(p.adjust))) +
    geom_point() +
    scale_colour_gradient(low = "#4575b4", high = "#d73027", name = "-log10 FDR") +
    scale_size_continuous(name = "genes") +
    labs(title = title, x = "Fold enrichment", y = NULL) +
    theme_bw(base_size = 11) +
    theme(panel.grid.minor = element_blank(),
          plot.title = element_text(size = 11, face = "bold"))
  h <- max(3, 0.32 * nrow(df) + 1.4)
  ggsave(paste0(file, ".png"), p, width = 9, height = h, dpi = 300)
  ggsave(paste0(file, ".pdf"), p, width = 9, height = h)   # editable in Illustrator
  invisible(p)
}

# ---- Metascape-style cross-list comparison -----------------------------------
# Zhou et al. 2019 (Nat Commun 10:1523) group enriched terms by pairwise Kappa
# similarity of their gene membership, cut the hierarchical tree at 0.3, and
# show the most significant term per cluster as one heatmap row. Reimplemented
# here on clusterProfiler output so the statistics stay the ones chosen above.

kappa_sim <- function(term_genes, universe) {
  M  <- vapply(term_genes, function(gs) as.integer(universe %in% gs),
               integer(length(universe)))
  if (is.null(dim(M))) M <- matrix(M, ncol = length(term_genes))
  n  <- nrow(M)
  A  <- crossprod(M)                      # a[i,j] = genes in both terms
  rs <- diag(A)
  Ri <- matrix(rs, nrow(A), ncol(A))       # rs[i]
  Rj <- t(Ri)                              # rs[j]
  B  <- Ri - A; C <- Rj - A; D <- n - A - B - C
  Po <- (A + D) / n
  Pe <- ((A + B) * (A + C) + (C + D) * (B + D)) / n^2
  K  <- (Po - Pe) / (1 - Pe)
  K[!is.finite(K)] <- 0
  diag(K) <- 1
  dimnames(K) <- list(names(term_genes), names(term_genes))
  K
}

cross_list_heatmap <- function(sig_store) {
  ids <- unique(unlist(lapply(sig_store, `[[`, "ID"), use.names = FALSE))
  if (length(ids) < 2) { msg("\n  <2 significant terms overall; skipping heatmap."); return(invisible(NULL)) }
  msg("\n== Cross-list comparison (%d sets, %d significant terms) ==",
      length(sig_store), length(ids))

  # gene membership per term, pooled over lists (readable symbols from enrichGO)
  term_genes <- lapply(ids, function(id) {
    unique(unlist(lapply(sig_store, function(d) {
      r <- d[d$ID == id, ]
      if (!nrow(r)) character(0) else strsplit(r$geneID[1], "/")[[1]]
    }), use.names = FALSE))
  })
  names(term_genes) <- ids
  universe <- unique(unlist(term_genes, use.names = FALSE))

  K  <- kappa_sim(term_genes, universe)
  hc <- hclust(as.dist(1 - K), method = "average")
  cl <- cutree(hc, h = 1 - g("kappa"))

  # best adjusted p per term across lists, and per cluster
  best_p <- vapply(ids, function(id)
    min(vapply(sig_store, function(d) {
      r <- d[d$ID == id, ]; if (nrow(r)) r$p.adjust[1] else 1
    }, numeric(1))), numeric(1))
  desc <- vapply(ids, function(id) {
    for (d in sig_store) { r <- d[d$ID == id, ]; if (nrow(r)) return(r$Description[1]) }
    id }, character(1))

  memb <- data.frame(cluster = cl[ids], ID = ids, Description = desc,
                     best_padj = best_p, stringsAsFactors = FALSE)
  memb <- memb[order(memb$cluster, memb$best_padj), ]
  reps <- memb[!duplicated(memb$cluster), ]
  reps <- head(reps[order(reps$best_padj), ], g("heatmap-top"))
  memb$is_representative <- memb$ID %in% reps$ID
  write.table(memb, file.path(g("outdir"), sprintf("%s_term_clusters.tsv", g("prefix"))),
              sep = "\t", row.names = FALSE, quote = FALSE)
  msg("  %d terms -> %d clusters at kappa >= %.2f; showing top %d",
      length(ids), length(unique(cl)), g("kappa"), nrow(reps))

  # matrix of -log10(adjusted p); NA where the term was not significant in a list
  mat <- do.call(cbind, lapply(sig_store, function(d)
    -log10(d$p.adjust[match(reps$ID, d$ID)])))
  colnames(mat) <- names(sig_store); rownames(mat) <- reps$Description
  write.table(data.frame(ID = reps$ID, Description = reps$Description, mat,
                         check.names = FALSE),
              file.path(g("outdir"), sprintf("%s_heatmap_matrix.tsv", g("prefix"))),
              sep = "\t", row.names = FALSE, quote = FALSE)

  long <- data.frame(
    term = factor(rep(rownames(mat), ncol(mat)), levels = rev(rownames(mat))),
    set  = factor(rep(colnames(mat), each = nrow(mat)), levels = colnames(mat)),
    val  = as.vector(mat), stringsAsFactors = FALSE)
  p <- ggplot(long, aes(x = set, y = term, fill = val)) +
    geom_tile(colour = "white", linewidth = 0.4) +
    scale_fill_gradient(low = "#fee8c8", high = "#b30000",
                        na.value = "grey90", name = "-log10 FDR") +
    labs(x = NULL, y = NULL,
         title = sprintf("Shared and list-specific GO %s clusters", g("ont")),
         caption = "grey = not significant in that list; one row per kappa cluster") +
    theme_minimal(base_size = 11) +
    theme(panel.grid = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1),
          plot.title = element_text(size = 11, face = "bold"))
  base <- file.path(g("outdir"), sprintf("%s_cross_list_heatmap", g("prefix")))
  h <- max(3, 0.28 * nrow(mat) + 1.6); w <- max(5, 1.1 * ncol(mat) + 5)
  ggsave(paste0(base, ".png"), p, width = w, height = h, dpi = 300)
  ggsave(paste0(base, ".pdf"), p, width = w, height = h)
  invisible(p)
}

summary_rows <- list()
sig_store    <- list()

for (n in names(sets_entrez)) {
  q <- sets_entrez[[n]]
  msg("\n== ORA: %s (%d genes, %d in universe) ==", n, length(q),
      if (is.null(bg_entrez)) NA_integer_ else length(bg_entrez))
  if (length(q) < 5) { msg("  too few mapped genes; skipping."); next }

  ego <- enrichGO(gene          = q,
                  universe      = bg_entrez,
                  OrgDb         = OrgDb,
                  keyType       = "ENTREZID",
                  ont           = g("ont"),
                  pAdjustMethod = g("p-method"),
                  pvalueCutoff  = 1,      # keep everything; filter downstream
                  qvalueCutoff  = 1,
                  minGSSize     = g("min-size"),
                  maxGSSize     = g("max-size"),
                  readable      = TRUE)

  if (is.null(ego) || !nrow(as.data.frame(ego))) { msg("  no GO terms tested."); next }
  res <- annotate_terms(as.data.frame(ego))

  if (g("topgo")) {
    onts <- if (g("ont") == "ALL") c("BP", "MF", "CC") else g("ont")
    tg <- do.call(rbind, lapply(onts, function(o) run_topgo(q, bg_entrez, o)))
    res <- merge(res, tg, by = "ID", all.x = TRUE)
  }

  res <- res[order(res$p.adjust, res$pvalue), ]
  keep <- res$p.adjust < g("padj") & res$Count >= g("min-count")
  sig  <- res[keep, , drop = FALSE]

  base <- file.path(g("outdir"), sprintf("%s_%s", g("prefix"), n))
  write.table(res, paste0(base, "_all_terms.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
  write.table(sig, paste0(base, "_significant.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
  msg("  %d terms tested, %d significant (%s < %g & Count >= %d)",
      nrow(res), nrow(sig), g("p-method"), g("padj"), g("min-count"))
  if (nrow(sig))
    msg("  median annotation size %.0f, median depth %.0f  (smaller/deeper = more specific)",
        median(sig$AnnotationSize, na.rm = TRUE), median(sig$Depth, na.rm = TRUE))

  if (nrow(sig)) dot_plot(sig, sprintf("GO %s enrichment - %s", g("ont"), n), paste0(base, "_dotplot"))

  if (g("simplify") && nrow(sig) && g("ont") != "ALL") {
    tryCatch({
      ego@result <- ego@result[ego@result$ID %in% sig$ID, , drop = FALSE]
      ego@pvalueCutoff <- 1; ego@qvalueCutoff <- 1
      s <- clusterProfiler::simplify(ego, cutoff = g("simplify-cutoff"),
                                     by = "p.adjust", select_fun = min,
                                     measure = "Wang")
      sdf <- annotate_terms(as.data.frame(s))
      write.table(sdf, paste0(base, "_significant_simplified.tsv"),
                  sep = "\t", row.names = FALSE, quote = FALSE)
      dot_plot(sdf, sprintf("GO %s (redundancy-reduced) - %s", g("ont"), n),
               paste0(base, "_dotplot_simplified"))
      msg("  simplify(): %d -> %d terms", nrow(sig), nrow(sdf))
    }, error = function(e) msg("  simplify() failed: %s", conditionMessage(e)))
  }

  if (nrow(sig)) sig_store[[n]] <- sig

  summary_rows[[n]] <- data.frame(
    set = n, genes_input = length(sets[[n]]), genes_mapped = length(q),
    terms_tested = nrow(res), terms_significant = nrow(sig),
    median_annotation_size = if (nrow(sig)) median(sig$AnnotationSize, na.rm = TRUE) else NA,
    median_depth = if (nrow(sig)) median(sig$Depth, na.rm = TRUE) else NA,
    stringsAsFactors = FALSE)
}

if (g("heatmap")) {
  if (length(sig_store) < 2)
    msg("\n  --heatmap needs >=2 sets with significant terms; skipping.")
  else cross_list_heatmap(sig_store)
}

# -------------------------------------------------------------- GSEA ---------

if (g("gsea")) {
  msg("\n== GSEA (gseGO) ==")
  gcol <- if (!is.null(g("gene-col"))) g("gene-col") else names(de)[1]
  lfc  <- suppressWarnings(as.numeric(de[[g("lfc-col")]]))
  rmet <- g("rank-metric")
  statcol <- g("stat-col")
  signed_cols   <- c("stat", "t", "statistic", "z")       # DESeq2, limma
  unsigned_cols <- c("F", "LR", "f", "lr", "Chisq")       # edgeR QLF / LRT
  if (is.null(statcol)) {
    hits <- c(intersect(signed_cols, names(de)), intersect(unsigned_cols, names(de)))
    statcol <- if (length(hits)) hits[1] else NULL
  }
  is_unsigned <- !is.null(statcol) && statcol %in% unsigned_cols

  pcol <- g("pval-col")
  if (is.null(pcol)) {
    hits <- intersect(c("pvalue", "PValue", "P.Value", "p_value", "pval", "p"), names(de))
    pcol <- if (length(hits)) hits[1] else g("padj-col")
  }

  if (rmet == "auto") {
    rmet <- if (is.null(statcol)) "signed_logp" else if (is_unsigned) "signed_sqrt" else "stat"
    msg("  --rank-metric auto -> %s%s", rmet,
        if (rmet %in% c("stat", "signed_sqrt")) sprintf(" (column '%s')", statcol)
        else sprintf(" (column '%s')", pcol))
  }
  if (rmet == "shrunk_lfc") rmet <- "lfc"

  metric <- switch(rmet,
    lfc = lfc,
    stat = { if (is.null(statcol) || !statcol %in% names(de))
               fail("--rank-metric stat requires a valid --stat-col.")
             suppressWarnings(as.numeric(de[[statcol]])) },
    # edgeR: QL F-test on a single contrast has df1 = 1, so sqrt(F) == |moderated t|
    # and the LRT statistic is chi-square_1, so sqrt(LR) == |z|. Re-attaching the
    # sign of logFC recovers the signed moderated statistic edgeR never prints.
    signed_sqrt = { if (is.null(statcol) || !statcol %in% names(de))
                      fail("--rank-metric signed_sqrt requires a valid --stat-col (e.g. F or LR).")
                    s <- suppressWarnings(as.numeric(de[[statcol]]))
                    sign(lfc) * sqrt(pmax(s, 0)) },
    signed_logp = {
      praw <- suppressWarnings(as.numeric(de[[pcol]]))
      sign(lfc) * -log10(pmax(praw, .Machine$double.xmin))
    },
    fail("--rank-metric must be auto, stat, signed_sqrt, signed_logp or lfc"))

  keepg <- !is.na(metric) & !is.na(de[[gcol]])
  ids <- de[[gcol]][keepg]; metric <- metric[keepg]
  if (id_type == "ENSEMBL") ids <- strip_ens_version(ids)
  if (id_type == "ENTREZID") {
    ent <- ids
  } else {
    m <- suppressWarnings(suppressMessages(bitr(ids, id_type, "ENTREZID", OrgDb)))
    idx <- match(m[[id_type]], ids); ent <- m$ENTREZID; metric <- metric[idx]
  }
  rk <- tapply(metric, ent, function(v) v[which.max(abs(v))])
  rk <- sort(rk, decreasing = TRUE)
  msg("  ranked list: %d genes (metric = %s)", length(rk), rmet)
  nties <- sum(duplicated(round(rk, 10)))
  if (nties > 0.02 * length(rk))
    msg(paste("  WARNING: %d/%d ranks are tied (%.1f%%). GSEA's running statistic",
              "is unstable with heavy ties -- prefer --rank-metric stat."),
        nties, length(rk), 100 * nties / length(rk))

  gs <- tryCatch(
    gseGO(geneList = rk, OrgDb = OrgDb, keyType = "ENTREZID", ont = g("ont"),
          minGSSize = g("min-size"), maxGSSize = g("max-size"),
          pvalueCutoff = 1, pAdjustMethod = g("p-method"),
          eps = 0, seed = TRUE, verbose = FALSE),
    error = function(e)            # older clusterProfiler has no `eps`
      gseGO(geneList = rk, OrgDb = OrgDb, keyType = "ENTREZID", ont = g("ont"),
            minGSSize = g("min-size"), maxGSSize = g("max-size"),
            pvalueCutoff = 1, pAdjustMethod = g("p-method"),
            seed = TRUE, verbose = FALSE))
  gdf <- as.data.frame(gs)
  if (nrow(gdf)) {
    gdf <- as.data.frame(setReadable(gs, OrgDb, keyType = "ENTREZID"))
    ont_col <- if ("ONTOLOGY" %in% names(gdf)) gdf$ONTOLOGY else rep(g("ont"), nrow(gdf))
    gdf$Depth <- if (length(unique(ont_col)) == 1) go_depth(gdf$ID, unique(ont_col)[1]) else NA
    gdf <- gdf[order(gdf$p.adjust), ]
    base <- file.path(g("outdir"), sprintf("%s_gsea", g("prefix")))
    write.table(gdf, paste0(base, "_all_terms.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
    sigg <- gdf[gdf$p.adjust < g("padj"), , drop = FALSE]
    write.table(sigg, paste0(base, "_significant.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
    msg("  %d gene sets tested, %d significant", nrow(gdf), nrow(sigg))
    if (nrow(sigg)) {
      d <- head(sigg[order(sigg$p.adjust), ], g("top-n"))
      d$Description <- factor(d$Description, levels = rev(d$Description[order(d$NES)]))
      p <- ggplot(d, aes(x = NES, y = Description, size = setSize, colour = -log10(p.adjust))) +
        geom_point() + geom_vline(xintercept = 0, linetype = 2, colour = "grey50") +
        scale_colour_gradient(low = "#4575b4", high = "#d73027", name = "-log10 FDR") +
        labs(title = sprintf("GSEA GO %s", g("ont")), x = "Normalised enrichment score", y = NULL) +
        theme_bw(base_size = 11) + theme(panel.grid.minor = element_blank())
      ggsave(paste0(base, "_dotplot.png"), p, width = 9,
             height = max(3, 0.32 * nrow(d) + 1.4), dpi = 300)
    }
  } else msg("  no gene sets returned.")
}

# ------------------------------------------------------ provenance log -------

if (length(summary_rows)) {
  s <- do.call(rbind, summary_rows)
  write.table(s, file.path(g("outdir"), sprintf("%s_summary.tsv", g("prefix"))),
              sep = "\t", row.names = FALSE, quote = FALSE)
  msg("\n== Summary ==")
  print(s, row.names = FALSE)
}

go_version <- tryCatch(as.character(packageVersion("GO.db")), error = function(e) NA)
log_path <- file.path(g("outdir"), sprintf("%s_run_log.txt", g("prefix")))
sink(log_path)
cat("go_enrichment.R run log\n")
cat("date:            ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), "\n")
cat("working dir:     ", getwd(), "\n")
cat("command:         ", paste(commandArgs(trailingOnly = FALSE), collapse = " "), "\n\n")
cat("parameters:\n")
for (k in names(opt)) if (k != "help") cat(sprintf("  %-24s %s\n", k,
    paste(as.character(opt[[k]]), collapse = ",")))
cat("\nannotation:\n")
cat("  OrgDb                 ", orgdb_name, as.character(packageVersion(orgdb_name)), "\n")
cat("  GO.db                 ", go_version, "\n")
cat("  clusterProfiler       ", as.character(packageVersion("clusterProfiler")), "\n")
if (g("topgo")) cat("  topGO                 ", as.character(packageVersion("topGO")), "\n")
cat("  resolved ID type      ", id_type, "\n")
cat("  background size       ", if (is.null(bg_entrez)) "whole annotation" else length(bg_entrez), "\n\n")
print(sessionInfo())
sink()

msg("\nOutputs written to %s/ (provenance in %s)", g("outdir"), basename(log_path))

# -----------------------------------------------------------------------------
# EXAMPLES
#
# 1) Plain gene list vs. expressed-gene background (the CLIP / peak-gene case):
#    Rscript go_enrichment.R \
#      --genes  idr_peak_genes.txt \
#      --background expressed_genes.txt \
#      --organism human --ont BP --simplify --topgo \
#      --outdir go/TERT_full --prefix TERT_full
#
# 2) DESeq2 / edgeR table, up and down tested separately, background = all
#    genes tested (derived automatically), plus GSEA on the ranked list:
#    Rscript go_enrichment.R \
#      --de-table deseq2_results.csv --gene-col gene_id \
#      --lfc-col log2FoldChange --padj-col padj \
#      --organism mouse --ont BP --split both --gsea --rank-metric stat --stat-col stat \
#      --outdir go/RPS2_STAMP --prefix RPS2
#
# 3) edgeR topTags table (columns logFC, logCPM, F, PValue, FDR). The gene IDs
#    live in the rownames, so export them first:
#      tt <- topTags(qlf, n = Inf)$table
#      tt$gene_id <- rownames(tt)
#      write.csv(tt[, c("gene_id", setdiff(names(tt), "gene_id"))],
#                "edger_TE_results.csv", row.names = FALSE)
#    Rscript go_enrichment.R \
#      --de-table edger_TE_results.csv --gene-col gene_id \
#      --lfc-col logFC --padj-col FDR --stat-col F \
#      --organism human --ont BP --split both --gsea \
#      --outdir go/TE --prefix TE
#    (--rank-metric auto sees the unsigned F and uses sign(logFC)*sqrt(F).)
#
# 4) Quick exploratory run, no background available (documented override):
#    Rscript go_enrichment.R --genes list.txt --allow-genome-background --outdir go/quick
# -----------------------------------------------------------------------------
