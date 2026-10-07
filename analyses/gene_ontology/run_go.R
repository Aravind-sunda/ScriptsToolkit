# =============================================================================
# run_go.R -- GO enrichment: edit the settings below, then run the whole file.
#
#   Rscript run_go.R              # from the shell
#   source("run_go.R")            # from R / RStudio
#
# Keep one copy of this file per project or per analysis. The engine lives in
# go_functions.R and does not need editing.
#
# Anything you do not want, set to NULL or FALSE.
# =============================================================================


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                               SETTINGS
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# --- the gene lists to test -------------------------------------------------
# Give each list a name. The value is a file (one gene per line, or a table
# whose first column is the gene) or a character vector you already have in R.
# One list is fine; several get compared against each other.

GENE_LISTS <- c(
  dTERT = "lists/dTERT_targets.txt",
  dRT   = "lists/dRT_targets.txt",
  dIDR  = "lists/dIDR_targets.txt"
)


# --- the background ---------------------------------------------------------
# Every gene that COULD have ended up in your lists: genes expressed in that
# cell type, genes with read coverage, genes tested by DESeq2 -- not all 20,000
# genes in the genome. This is the single most important setting in the file.
#
# Set to NULL only if you truly have no background, and then also set
# ALLOW_GENOME_BACKGROUND <- TRUE below to confirm you meant it.

BACKGROUND <- "lists/expressed_genes.txt"

ALLOW_GENOME_BACKGROUND <- FALSE


# --- organism and gene IDs --------------------------------------------------
ORGANISM <- "human"   # human, mouse, rat, zebrafish, fly, worm, yeast,
                      # or a package name like "org.Hs.eg.db"

ID_TYPE  <- NULL      # NULL = detect automatically.
                      # Or force: "SYMBOL", "ENSEMBL", "ENTREZID", "REFSEQ"


# --- which gene sets to test ------------------------------------------------
# Any combination of the three GO domains. Each one is run separately and gets
# its own output files.
#   BP = biological process,  MF = molecular function,  CC = cellular component

ONTOLOGIES <- c("BP", "MF")


# --- statistics -------------------------------------------------------------
MIN_TERM_SIZE <- 10     # ignore GO terms with fewer genes than this in the background
MAX_TERM_SIZE <- 500    # ignore huge, generic terms ("metabolic process")
MIN_GENES_HIT <- 3      # a term needs at least this many of your genes to be reported
PADJ_CUTOFF   <- 0.05   # adjusted p threshold for the significant table
P_ADJUST      <- "BH"   # BH, BY, holm or bonferroni


# --- redundancy and plots ---------------------------------------------------
SIMPLIFY          <- TRUE   # flag which significant terms are representative
SIMPLIFY_CUTOFF   <- 0.7    # semantic similarity above which terms are redundant

TOP_N             <- 20     # terms per dotplot
PLOT_RANK_BY      <- "p.adjust"   # "p.adjust" (most significant) or
                                  # "FoldEnrichment" (strongest effect)

# Cross-list comparison. Only does anything when you give more than one list.
COMPARE_LISTS     <- TRUE
KAPPA_CUTOFF      <- 0.3    # how aggressively to merge redundant terms; lower merges more


# --- output -----------------------------------------------------------------
OUTDIR <- "go_results"
PREFIX <- "GO"


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                    nothing below here needs editing
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Find go_functions.R, however this script was launched: next to it, one level
# up (if you keep configs in a subfolder), or in the working directory.
.this_dir <- tryCatch({
  args <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", args[grep("^--file=", args)])
  if (length(f) > 0) dirname(normalizePath(f)) else dirname(sys.frame(1)$ofile)
}, error = function(e) ".")
if (is.null(.this_dir) || is.na(.this_dir) || !nzchar(.this_dir)) .this_dir <- "."

.candidates <- file.path(c(.this_dir, file.path(.this_dir, ".."), "."),
                         "go_functions.R")
.engine <- .candidates[file.exists(.candidates)]
if (length(.engine) == 0) {
  stop("Cannot find go_functions.R. Looked in:\n  ",
       paste(normalizePath(.candidates, mustWork = FALSE), collapse = "\n  "),
       call. = FALSE)
}
source(.engine[1])


# --- 1. checks --------------------------------------------------------------

if (!all(ONTOLOGIES %in% c("BP", "MF", "CC")))
  stop("ONTOLOGIES must be any of \"BP\", \"MF\", \"CC\" (case-sensitive). Got: ",
       paste(ONTOLOGIES, collapse = ", "), call. = FALSE)

if (!PLOT_RANK_BY %in% c("p.adjust", "FoldEnrichment"))
  stop('PLOT_RANK_BY must be "p.adjust" or "FoldEnrichment". Got: ', PLOT_RANK_BY,
       call. = FALSE)

if (!P_ADJUST %in% p.adjust.methods)
  stop("P_ADJUST must be one of: ", paste(p.adjust.methods, collapse = ", "),
       call. = FALSE)

if (is.null(BACKGROUND) && !ALLOW_GENOME_BACKGROUND) {
  stop("No BACKGROUND set.\n",
       "  An ORA background must be the genes that could plausibly have appeared\n",
       "  in your lists -- expressed genes, tested genes, covered genes. Using the\n",
       "  whole genome instead inflates enrichment and is the most common error in\n",
       "  the published literature (95% of papers, Wijesooriya et al. 2022).\n",
       "  Set BACKGROUND, or set ALLOW_GENOME_BACKGROUND <- TRUE on purpose.",
       call. = FALSE)
}

dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)


# --- 2. load packages and read the genes ------------------------------------

org <- go_setup(ORGANISM)

message("\n== Input ==")
gene_lists <- read_all_gene_lists(GENE_LISTS)
message("Query lists: ",
        paste(sprintf("%s (n=%d)", names(gene_lists), lengths(gene_lists)),
              collapse = ", "))

id_type <- ID_TYPE
if (is.null(id_type)) {
  id_type <- guess_id_type(unlist(gene_lists, use.names = FALSE))
  message("Detected gene ID type: ", id_type)
}

message("ID mapping:")
if (!is.null(BACKGROUND)) {
  # Read it as its own step: passing read_gene_list() straight into
  # map_to_entrez() means a "file not found" surfaces from deep inside a lazily
  # evaluated argument, which buries the message.
  background_genes <- read_gene_list(BACKGROUND)
  background <- map_to_entrez(background_genes, id_type, org$db, "background")
} else {
  background <- NULL
  message("  WARNING: no background -- using the whole ", org$name, " annotation.")
}

lists_entrez <- lapply(names(gene_lists), function(nm)
  map_to_entrez(gene_lists[[nm]], id_type, org$db, nm))
names(lists_entrez) <- names(gene_lists)

# A query gene missing from the background would be tested against a universe
# it isn't in, which is incoherent. Add them and say so.
if (!is.null(background)) {
  for (nm in names(lists_entrez)) {
    outside <- setdiff(lists_entrez[[nm]], background)
    if (length(outside) > 0)
      message(sprintf("  note: %d '%s' gene(s) were not in the background; added.",
                      length(outside), nm))
  }
  background <- union(background, unlist(lists_entrez, use.names = FALSE))
}


# --- 3. run ORA for every list x every ontology -----------------------------

summary_rows <- list()

for (ont in ONTOLOGIES) {

  message("\n== GO:", ont, " ==")

  # The honest denominator for multiple-testing correction: every term that
  # could have been tested, not only the ones your genes happened to hit.
  # With no background of your own, the universe is the whole annotation.
  counting_universe <- if (is.null(background))
    AnnotationDbi::keys(org$db, keytype = "ENTREZID") else background
  n_testable <- count_testable_terms(counting_universe, ont, org$db,
                                     MIN_TERM_SIZE, MAX_TERM_SIZE)
  message("Testable GO:", ont, " terms in this background: ", n_testable)

  full_tables <- list()
  sig_tables  <- list()

  for (nm in names(lists_entrez)) {

    result <- run_ora(lists_entrez[[nm]], background, ont, org$db,
                      MIN_TERM_SIZE, MAX_TERM_SIZE, P_ADJUST)

    if (is.null(result)) {
      message(sprintf("  %-12s no terms could be tested", nm))
      next
    }

    terms <- as.data.frame(result)
    terms <- recorrect_pvalues(terms, n_testable, P_ADJUST)
    terms <- add_term_metrics(terms, ont)

    sig <- filter_significant(terms, PADJ_CUTOFF, MIN_GENES_HIT)

    # Which of the significant terms are representative rather than restatements
    # of a neighbour. Adds a column; removes nothing.
    terms$Representative <- NA
    if (SIMPLIFY && nrow(sig) > 1) {
      rep_ids <- mark_representative_terms(result, sig$ID, SIMPLIFY_CUTOFF)
      terms$Representative[terms$ID %in% sig$ID] <- terms$ID[terms$ID %in% sig$ID] %in% rep_ids
      sig$Representative <- sig$ID %in% rep_ids
    }

    stem <- file.path(OUTDIR, paste(PREFIX, safe_name(nm), ont, sep = "_"))
    write_tsv(terms, paste0(stem, "_all_terms.tsv"))
    write_tsv(sig,   paste0(stem, "_significant.tsv"))
    plot_dotplot(sig, sprintf("%s  --  GO:%s", nm, ont), paste0(stem, "_dotplot"),
                 top_n = TOP_N, rank_by = PLOT_RANK_BY)

    full_tables[[nm]] <- terms
    sig_tables[[nm]]  <- sig

    message(sprintf("  %-12s %5d terms tested, %4d significant (FDR < %g)",
                    nm, nrow(terms), nrow(sig), PADJ_CUTOFF))

    # Median annotation size and depth tell you whether your output is specific
    # or has drifted to "biological regulation" territory.
    summary_rows[[length(summary_rows) + 1]] <- data.frame(
      List           = nm,
      Ontology       = ont,
      GenesIn        = length(gene_lists[[nm]]),
      GenesMapped    = length(lists_entrez[[nm]]),
      TermsTested    = n_testable,
      TermsSignif    = nrow(sig),
      TermsRepresent = if (SIMPLIFY && "Representative" %in% names(sig))
                         sum(sig$Representative, na.rm = TRUE) else NA_integer_,
      MedianAnnotSize = if (nrow(sig)) median(sig$AnnotationSize) else NA_real_,
      MedianDepth     = if (nrow(sig)) median(sig$Depth, na.rm = TRUE) else NA_real_,
      MedianFoldEnr   = if (nrow(sig)) round(median(sig$FoldEnrichment), 2) else NA_real_,
      stringsAsFactors = FALSE
    )
  }

  # --- 4. compare the lists to each other -----------------------------------

  if (COMPARE_LISTS && length(sig_tables) > 1) {
    message("Comparing lists...")
    clusters <- compare_lists_heatmap(
      sig_tables, full_tables,
      file_stem    = file.path(OUTDIR, paste(PREFIX, "comparison", ont, sep = "_")),
      kappa_cutoff = KAPPA_CUTOFF, top_n = TOP_N, padj_cutoff = PADJ_CUTOFF
    )
    if (!is.null(clusters))
      write_tsv(clusters, file.path(OUTDIR,
                paste0(paste(PREFIX, "comparison", ont, sep = "_"), "_term_clusters.tsv")))
  }
}


# --- 5. summary and run log -------------------------------------------------

if (length(summary_rows) > 0) {
  summary_table <- do.call(rbind, summary_rows)
  write_tsv(summary_table, file.path(OUTDIR, paste0(PREFIX, "_summary.tsv")))
  message("\n== Summary ==")
  print(summary_table, row.names = FALSE)
}

write_run_log(
  file.path(OUTDIR, paste0(PREFIX, "_run_log.txt")),
  settings = list(
    gene_lists = paste(names(GENE_LISTS), unlist(GENE_LISTS), sep = "="),
    background = if (is.null(BACKGROUND)) "NONE (whole genome)" else BACKGROUND,
    background_size = length(background),
    organism = ORGANISM, orgdb = org$name, id_type = id_type,
    ontologies = ONTOLOGIES,
    min_term_size = MIN_TERM_SIZE, max_term_size = MAX_TERM_SIZE,
    min_genes_hit = MIN_GENES_HIT, padj_cutoff = PADJ_CUTOFF,
    p_adjust = P_ADJUST, simplify = SIMPLIFY, simplify_cutoff = SIMPLIFY_CUTOFF,
    kappa_cutoff = KAPPA_CUTOFF
  ),
  extra_lines = c(
    "p.adjust           = clusterProfiler's own correction, over terms it tested",
    "p.adjust.allterms  = correction over every testable term (what to report)",
    "enrichment_force_universe was set to FALSE: it is miscalibrated in",
    "  DOSE 4.4.0. See the comment in go_functions.R::go_setup()."
  )
)

message("\nDone. Results in: ", normalizePath(OUTDIR))
