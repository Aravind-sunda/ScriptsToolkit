# =============================================================================
# go_functions.R -- the engine for GO over-representation analysis (ORA)
#
# This file only DEFINES functions. It changes nothing when you source it.
# Set your options in run_go.R and run that instead.
#
# You can also source this file inside an .Rmd and call the functions directly.
#
# Reading order, top to bottom:
#   1. Setup ................ load the right annotation package
#   2. Input ................ read gene lists, work out the ID type, map to Entrez
#   3. ORA .................. run the hypergeometric test, one list x one ontology
#   4. Term metrics ......... fold enrichment, annotation size, GO depth
#   5. Redundancy ........... mark which terms are representative
#   6. Plots ................ dotplot per list
#   7. Compare lists ........ kappa clustering + cross-list heatmap
#   8. Output ............... write tables and a run log
#
# Why these choices are what they are, with citations, is in README.md.
# =============================================================================


# =============================================================================
# 1. SETUP
# =============================================================================

# Short organism names -> Bioconductor annotation package.
# Add a row here if you need another organism.
ORGDB_PACKAGES <- c(
  human     = "org.Hs.eg.db",
  mouse     = "org.Mm.eg.db",
  rat       = "org.Rn.eg.db",
  zebrafish = "org.Dr.eg.db",
  fly       = "org.Dm.eg.db",
  worm      = "org.Ce.eg.db",
  yeast     = "org.Sc.sgd.db"
)


#' Load clusterProfiler and the annotation package for one organism.
#'
#' @param organism  "human", "mouse", ... or a package name like "org.Hs.eg.db"
#' @return list(db = the OrgDb object, name = the package name)
go_setup <- function(organism) {

  pkg <- if (organism %in% names(ORGDB_PACKAGES)) ORGDB_PACKAGES[[organism]] else organism

  needed <- c("clusterProfiler", "GO.db", "ggplot2", pkg)
  missing <- needed[!vapply(needed, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing) > 0) {
    stop("Missing R packages: ", paste(missing, collapse = ", "), "\n",
         '  Install with: BiocManager::install(c("',
         paste(missing, collapse = '", "'), '"))', call. = FALSE)
  }

  suppressPackageStartupMessages({
    library(clusterProfiler)
    library(ggplot2)
    library(pkg, character.only = TRUE)
  })

  # Ziemann et al. (2024) point out that ORA tools quietly drop background genes
  # with no annotation, which shrinks the denominator. clusterProfiler has an
  # undocumented `enrichment_force_universe` option that is meant to fix this.
  #
  # DO NOT TURN IT ON. Tested against DOSE 4.4.0: it replaces the background
  # denominator with the full universe but leaves the query denominator at the
  # annotated count, so the two sides of the hypergeometric test no longer
  # match. A random 150-gene list drawn from a 6000-gene background returned
  # 383 significant BP terms with it on, and 0 with it off. The default
  # behaviour is internally consistent -- annotated query tested against
  # annotated background -- and is correctly calibrated.
  #
  # Set explicitly rather than left alone, in case it is on in your .Rprofile.
  options(enrichment_force_universe = FALSE)

  list(db = get(pkg), name = pkg)
}


# =============================================================================
# 2. INPUT
# =============================================================================

#' Read one gene list.
#'
#' Accepts either a file path or a character vector you already have in R.
#' Files may be one gene per line, or a table -- the first column is used.
#' A header line that looks like a column name is dropped.
read_gene_list <- function(x) {

  # Already a vector of genes? Just clean it up.
  if (length(x) > 1 || !file.exists(x[1])) {

    # ...unless it looks like a file path, in which case it is almost certainly
    # a typo. Without this check a misspelled path is silently accepted as a
    # gene list containing one gene named "lists/my_typo.txt", and you don't
    # find out until the mapping step fails for a confusing reason.
    if (length(x) == 1 &&
        (grepl("[/\\\\]", x) ||
         grepl("\\.(txt|csv|tsv|list|genes|dat)$", x, ignore.case = TRUE))) {
      stop("File not found: ", x,
           "\n  (it was read as a gene name because no such file exists)",
           call. = FALSE)
    }
    return(unique(trimws(x[nzchar(trimws(x))])))
  }

  lines <- readLines(x, warn = FALSE)
  lines <- lines[nzchar(trimws(lines))]

  # If it is a table, take the first field of each line.
  genes <- vapply(strsplit(lines, "[\t,]"), function(f) f[1], character(1))
  genes <- trimws(gsub('"', "", genes))

  # Drop a header row like "gene", "gene_id", "symbol".
  if (length(genes) > 0 &&
      grepl("^(gene|gene_?id|gene_?name|symbol|id)$", genes[1], ignore.case = TRUE)) {
    genes <- genes[-1]
  }

  unique(genes[nzchar(genes)])
}


#' Read all the query lists into a named list of character vectors.
#'
#' @param spec  named vector or named list; each element is a file path or a
#'              character vector of genes.
read_all_gene_lists <- function(spec) {

  if (is.null(names(spec)) || any(!nzchar(names(spec)))) {
    stop("GENE_LISTS must be named, e.g. c(treated = 'up.txt', control = 'down.txt')",
         call. = FALSE)
  }

  lists <- lapply(spec, read_gene_list)
  names(lists) <- names(spec)

  empty <- names(lists)[lengths(lists) == 0]
  if (length(empty) > 0) stop("These gene lists are empty: ",
                              paste(empty, collapse = ", "), call. = FALSE)
  lists
}


#' Guess whether genes are symbols, Ensembl IDs, Entrez IDs or RefSeq.
#'
#' Looks at the first 200 genes and picks whichever pattern most of them match.
guess_id_type <- function(genes) {
  sample_genes <- genes[!is.na(genes)][seq_len(min(200, length(genes)))]

  if (mean(grepl("^ENS[A-Z]*[GT]?\\d{6,}", sample_genes)) > 0.5) return("ENSEMBL")
  if (mean(grepl("^\\d+$",                 sample_genes)) > 0.8) return("ENTREZID")
  if (mean(grepl("^[NX][MR]_",             sample_genes)) > 0.5) return("REFSEQ")
  "SYMBOL"
}


#' Convert gene IDs to Entrez IDs, which is what clusterProfiler works in.
#'
#' Prints the mapping rate. Unreported ID loss is a known reproducibility
#' problem (Wijesooriya et al. 2022), so this is deliberately noisy.
#'
#' @return character vector of unique Entrez IDs
map_to_entrez <- function(genes, id_type, orgdb, label = "genes") {

  # Ensembl IDs often carry a version suffix (ENSG00000141510.16) that no
  # annotation package recognises.
  if (id_type == "ENSEMBL") genes <- sub("\\.\\d+$", "", genes)
  genes <- unique(genes[!is.na(genes) & nzchar(genes)])

  if (id_type == "ENTREZID") {
    mapped     <- genes[genes %in% AnnotationDbi::keys(orgdb, keytype = "ENTREZID")]
    n_recognised <- length(mapped)
  } else {
    hits <- suppressWarnings(suppressMessages(
      clusterProfiler::bitr(genes, fromType = id_type,
                            toType = "ENTREZID", OrgDb = orgdb)
    ))
    if (is.null(hits) || nrow(hits) == 0) {
      stop("Could not map any ", label, " from ID type ", id_type,
           ". Check ID_TYPE and ORGANISM.", call. = FALSE)
    }
    mapped <- unique(hits$ENTREZID)
    # Mapping is 1:many -- a few symbols hit two Entrez IDs -- so the success
    # rate is "how many of my genes were recognised", not "how many IDs came out".
    n_recognised <- length(unique(hits[[id_type]]))
  }

  message(sprintf("  %-22s %6d in -> %6d recognised (%.1f%%) -> %6d Entrez IDs",
                  label, length(genes), n_recognised,
                  100 * n_recognised / max(1, length(genes)), length(mapped)))
  mapped
}


# =============================================================================
# 3. ORA
# =============================================================================

#' Count every GO term that could have been tested.
#'
#' clusterProfiler only tests terms that contain at least one query gene, and
#' corrects p-values across just those. A term with zero overlap is still a test
#' that was performed, and leaving it out makes the FDR too small
#' (Ziemann et al. 2024). This counts the real number of tests: every GO term in
#' this ontology whose size *in your background* is within the size limits.
#'
#' @return integer, the correct denominator for multiple-testing correction
count_testable_terms <- function(background, ontology, orgdb,
                                 min_size, max_size) {

  # GOALL propagates a gene up to every ancestor term, which is what enrichGO
  # tests against -- so term sizes here match the ones enrichGO uses.
  ann <- suppressMessages(AnnotationDbi::select(
    orgdb, keys = background, keytype = "ENTREZID",
    columns = c("GOALL", "ONTOLOGYALL")
  ))
  ann <- ann[!is.na(ann$GOALL) & ann$ONTOLOGYALL == ontology, ]

  term_sizes <- tapply(ann$ENTREZID, ann$GOALL, function(g) length(unique(g)))
  sum(term_sizes >= min_size & term_sizes <= max_size)
}


#' Run ORA for one gene list against one ontology.
#'
#' Everything is returned, not just the significant terms, so nothing is lost
#' silently and you can see the shape of the whole result.
#'
#' @return the enrichGO result object, or NULL if nothing was testable
run_ora <- function(genes, background, ontology, orgdb,
                    min_size, max_size, p_method) {

  result <- clusterProfiler::enrichGO(
    gene          = genes,
    universe      = background,
    OrgDb         = orgdb,
    ont           = ontology,
    keyType       = "ENTREZID",
    minGSSize     = min_size,
    maxGSSize     = max_size,
    pAdjustMethod = p_method,
    pvalueCutoff  = 1,        # keep every term; we filter ourselves, later
    qvalueCutoff  = 1,
    readable      = TRUE      # geneID column comes back as gene symbols
  )

  if (is.null(result) || nrow(as.data.frame(result)) == 0) return(NULL)
  result
}


#' Redo the multiple-testing correction over ALL testable terms.
#'
#' p.adjust's `n` argument exists for exactly this: it lets you correct as if
#' you had run `n` tests when your vector only holds the non-zero ones.
#' Adds `p.adjust.allterms`, and keeps clusterProfiler's original `p.adjust`
#' beside it so you can see what the difference was.
recorrect_pvalues <- function(terms, n_testable, p_method) {

  n <- max(n_testable, nrow(terms))   # never correct over fewer tests than we have
  terms$p.adjust.allterms <- p.adjust(terms$pvalue, method = p_method, n = n)
  terms$n.terms.tested    <- n
  terms
}


# =============================================================================
# 4. TERM METRICS
# =============================================================================

# clusterProfiler stores ratios as text like "12/300". These pull the numbers out.
ratio_top    <- function(x) as.numeric(sub("/.*", "", x))
ratio_bottom <- function(x) as.numeric(sub(".*/", "", x))


#' How far a GO term sits below the root of its ontology.
#'
#' Deeper = more specific. "cell cycle" is shallow, "mitotic spindle midzone
#' assembly" is deep. The Frontiers 2026 benchmark uses depth and annotation
#' size to judge whether a tool returns informative terms or generic ones, so
#' we report the same two numbers for your own results.
go_depth <- function(go_ids, ontology) {

  parent_map <- switch(ontology,
                       BP = GO.db::GOBPPARENTS,
                       MF = GO.db::GOMFPARENTS,
                       CC = GO.db::GOCCPARENTS)
  if (is.null(parent_map)) return(rep(NA_integer_, length(go_ids)))

  parents <- as.list(parent_map)
  cache   <- new.env(hash = TRUE, parent = emptyenv())

  # Longest path to the root, walking up the graph. Cached because GO terms
  # share ancestors heavily and this would otherwise be very slow.
  depth_of <- function(id) {
    if (!is.null(cache[[id]])) return(cache[[id]])
    above <- setdiff(unname(parents[[id]]), "all")
    d <- if (length(above) == 0) 0L else max(vapply(above, depth_of, integer(1))) + 1L
    cache[[id]] <- d
    d
  }

  vapply(go_ids,
         function(id) if (id %in% names(parents)) depth_of(id) else NA_integer_,
         integer(1), USE.NAMES = FALSE)
}


#' Add effect size and specificity columns to a result table.
#'
#' FoldEnrichment  how much more often the term appears in your list than in
#'                 the background. Reporting an effect size, not just a p-value,
#'                 is one of the ten common mistakes in Bora et al. 2026.
#' AnnotationSize  genes annotated to the term in your background. Big = generic.
#' Depth           distance from the ontology root. Deep = specific.
add_term_metrics <- function(terms, ontology) {

  in_list <- ratio_top(terms$GeneRatio)  / ratio_bottom(terms$GeneRatio)
  in_bg   <- ratio_top(terms$BgRatio)    / ratio_bottom(terms$BgRatio)

  terms$FoldEnrichment <- in_list / in_bg
  terms$AnnotationSize <- ratio_top(terms$BgRatio)
  terms$Depth          <- go_depth(terms$ID, ontology)
  terms$Ontology       <- ontology
  terms
}


#' Keep the terms that pass the significance and size filters.
filter_significant <- function(terms, padj_cutoff, min_count,
                               padj_column = "p.adjust.allterms") {
  keep <- terms[[padj_column]] < padj_cutoff & terms$Count >= min_count
  terms[which(keep), , drop = FALSE]
}


# =============================================================================
# 5. REDUNDANCY
# =============================================================================

#' Flag which significant terms are representative rather than redundant.
#'
#' GO is a hierarchy, so a real signal shows up as a dozen nested restatements
#' of the same thing. clusterProfiler::simplify groups terms by Wang semantic
#' similarity and keeps the most significant one per group.
#'
#' Note it sometimes keeps the generic parent and drops the specific child, so
#' this ADDS a `Representative` column rather than deleting rows. Nothing is
#' thrown away; you decide.
mark_representative_terms <- function(ora_result, significant_ids, cutoff) {

  if (length(significant_ids) < 2) return(significant_ids)
  if (!requireNamespace("GOSemSim", quietly = TRUE)) {
    warning("GOSemSim is not installed, skipping redundancy reduction.", call. = FALSE)
    return(significant_ids)
  }

  # simplify() works on the result object, so hand it a copy holding only the
  # significant terms -- semantic similarity over thousands of terms is slow.
  sig_only <- ora_result
  sig_only@result <- ora_result@result[ora_result@result$ID %in% significant_ids, ,
                                       drop = FALSE]

  reduced <- try(clusterProfiler::simplify(sig_only, cutoff = cutoff,
                                           by = "p.adjust", select_fun = min),
                 silent = TRUE)
  if (inherits(reduced, "try-error")) {
    warning("simplify() failed, keeping all terms.", call. = FALSE)
    return(significant_ids)
  }
  as.data.frame(reduced)$ID
}


# =============================================================================
# 6. PLOTS
# =============================================================================

#' Dotplot of the top terms for one gene list.
#'
#' x = fold enrichment (the effect size), dot size = how many of your genes hit
#' the term, colour = significance. Bubble charts carrying both significance and
#' effect size are what Bora et al. (2026) ask for.
#'
#' @param rank_by  "p.adjust" to show the most significant terms, or
#'                 "FoldEnrichment" to show the strongest ones.
plot_dotplot <- function(terms, title, file_stem, top_n = 20,
                         rank_by = "p.adjust", padj_column = "p.adjust.allterms") {

  if (nrow(terms) == 0) return(invisible(NULL))

  ord <- if (rank_by == "FoldEnrichment") order(-terms$FoldEnrichment)
         else                             order(terms[[padj_column]])
  top <- terms[head(ord, top_n), , drop = FALSE]

  # Long GO names wreck the layout, so shorten them. GO has plenty of terms that
  # are identical for their first 52 characters ("regulation of transcription
  # from RNA polymerase II promoter in response to ..."), and duplicate factor
  # levels are an error -- make.unique keeps every row addressable.
  top$Label <- ifelse(nchar(top$Description) > 55,
                      paste0(substr(top$Description, 1, 52), "..."),
                      top$Description)
  top$Label <- make.unique(top$Label)
  top$Label <- factor(top$Label, levels = rev(top$Label[order(top$FoldEnrichment)]))

  p <- ggplot2::ggplot(top, ggplot2::aes(x = FoldEnrichment, y = Label)) +
    ggplot2::geom_point(ggplot2::aes(size = Count,
                                     colour = -log10(.data[[padj_column]]))) +
    ggplot2::scale_colour_gradient(low = "#4575b4", high = "#d73027") +
    ggplot2::labs(title = title, x = "Fold enrichment", y = NULL,
                  size = "Genes", colour = "-log10(FDR)") +
    ggplot2::theme_bw(base_size = 11) +
    ggplot2::theme(panel.grid.minor = ggplot2::element_blank())

  height <- max(3, 0.28 * nrow(top) + 1.5)
  ggplot2::ggsave(paste0(file_stem, ".png"), p, width = 8, height = height, dpi = 300)
  ggplot2::ggsave(paste0(file_stem, ".pdf"), p, width = 8, height = height)
  invisible(p)
}


# =============================================================================
# 7. COMPARE LISTS
# =============================================================================

#' Kappa similarity between every pair of terms, based on which genes they hit.
#'
#' Cohen's kappa on the two binary membership vectors, i.e. "do these two terms
#' pick out the same genes, more than you would expect by chance". This is the
#' measure Metascape (Zhou et al. 2019) uses to group redundant terms.
kappa_matrix <- function(term_genes) {

  all_genes <- sort(unique(unlist(term_genes)))
  # rows = terms, cols = genes, TRUE if the term contains that gene
  membership <- t(vapply(term_genes, function(g) all_genes %in% g,
                         logical(length(all_genes))))
  n_terms <- nrow(membership)
  n_genes <- length(all_genes)

  both    <- membership %*% t(membership)             # in A and in B
  size    <- diag(both)
  neither <- n_genes - outer(size, size, "+") + both  # in neither

  observed <- (both + neither) / n_genes
  # chance agreement from the marginal totals
  expected <- (outer(size, size) +
               outer(n_genes - size, n_genes - size)) / n_genes^2

  kappa <- (observed - expected) / (1 - expected)
  kappa[!is.finite(kappa)] <- 0
  diag(kappa) <- 1
  dimnames(kappa) <- list(names(term_genes), names(term_genes))
  kappa
}


#' Group redundant terms into clusters and pick one representative each.
#'
#' Average-linkage hierarchical clustering on (1 - kappa), cut at the chosen
#' similarity. The representative is the term with the best adjusted p across
#' all the lists.
#'
#' @return data.frame: ID, Description, Cluster, Representative (TRUE/FALSE)
cluster_terms <- function(term_genes, best_padj, descriptions, kappa_cutoff) {

  if (length(term_genes) < 2) {
    return(data.frame(ID = names(term_genes),
                      Description = descriptions[names(term_genes)],
                      Cluster = 1L, Representative = TRUE,
                      stringsAsFactors = FALSE))
  }

  k <- kappa_matrix(term_genes)
  tree <- hclust(as.dist(1 - k), method = "average")
  cluster <- cutree(tree, h = 1 - kappa_cutoff)

  out <- data.frame(
    ID          = names(cluster),
    Description = descriptions[names(cluster)],
    Cluster     = as.integer(cluster),
    BestPadj    = best_padj[names(cluster)],
    stringsAsFactors = FALSE
  )
  # Best term in each cluster represents it.
  out$Representative <- FALSE
  for (cl in unique(out$Cluster)) {
    rows <- which(out$Cluster == cl)
    out$Representative[rows[which.min(out$BestPadj[rows])]] <- TRUE
  }
  out[order(out$Cluster, out$BestPadj), ]
}


#' Cross-list comparison heatmap, Metascape style.
#'
#' Each list was enriched separately against the same background. Redundant
#' terms are collapsed into clusters so each row is a distinct theme rather than
#' the same GO branch restated ten times. Colour is -log10(FDR); grey means the
#' term was not significant in that list.
#'
#' @param sig_tables   named list of significant-term data.frames, one per list
#' @param full_tables  named list of the corresponding complete tables -- used to
#'                     look up a term's p-value in a list where it missed the cutoff
#' @return data.frame of the term clusters (also written to disk by the caller)
compare_lists_heatmap <- function(sig_tables, full_tables, file_stem,
                                  kappa_cutoff = 0.3, top_n = 20,
                                  padj_cutoff = 0.05,
                                  padj_column = "p.adjust.allterms") {

  sig_tables <- sig_tables[vapply(sig_tables, nrow, integer(1)) > 0]
  if (length(sig_tables) < 2) {
    message("  Fewer than 2 lists have significant terms; skipping heatmap.")
    return(NULL)
  }

  pooled <- do.call(rbind, lapply(names(sig_tables), function(nm) {
    d <- sig_tables[[nm]]
    data.frame(List = nm, ID = d$ID, Description = d$Description,
               padj = d[[padj_column]], geneID = d$geneID,
               stringsAsFactors = FALSE)
  }))

  # One entry per term: its genes, its name, and its best p across the lists.
  term_genes   <- tapply(pooled$geneID, pooled$ID,
                         function(x) unique(unlist(strsplit(x, "/"))))
  descriptions <- tapply(pooled$Description, pooled$ID, function(x) x[1])
  best_padj    <- tapply(pooled$padj, pooled$ID, min)

  clusters <- cluster_terms(term_genes, best_padj, descriptions, kappa_cutoff)

  # Show the best `top_n` clusters, one row each.
  reps <- clusters[clusters$Representative, ]
  reps <- head(reps[order(reps$BestPadj), ], top_n)

  # Look each representative up in every list's FULL table, so a term that just
  # missed the cutoff in one list still shows a value rather than a hole.
  cells <- expand.grid(ID = reps$ID, List = names(full_tables),
                       stringsAsFactors = FALSE)
  cells$padj <- mapply(function(id, lst) {
    row <- full_tables[[lst]][full_tables[[lst]]$ID == id, ]
    if (nrow(row) == 0) NA_real_ else row[[padj_column]][1]
  }, cells$ID, cells$List)

  cells$score <- ifelse(is.na(cells$padj) | cells$padj >= padj_cutoff,
                        NA_real_, -log10(cells$padj))

  # Shorten long GO names once, and use the SAME shortened text for the values
  # and for the row order -- otherwise a truncated label falls outside its own
  # factor levels and the row plots as NA.
  shorten <- function(x) ifelse(nchar(x) > 55, paste0(substr(x, 1, 52), "..."), x)

  # One label per representative GO ID, made unique so that two terms sharing
  # their first 52 characters stay two rows instead of silently merging into one.
  row_labels <- make.unique(shorten(as.character(descriptions[reps$ID])))
  names(row_labels) <- reps$ID

  # Rows in cluster order (best at the top), columns in the order given.
  cells$Description <- factor(row_labels[cells$ID], levels = rev(row_labels))
  cells$List <- factor(cells$List, levels = names(full_tables))

  p <- ggplot2::ggplot(cells, ggplot2::aes(x = List, y = Description, fill = score)) +
    ggplot2::geom_tile(colour = "white", linewidth = 0.4) +
    ggplot2::scale_fill_gradient(low = "#fee5d9", high = "#a50f15",
                                 na.value = "grey90") +
    ggplot2::labs(x = NULL, y = NULL, fill = "-log10(FDR)",
                  title = "Enriched GO terms across gene lists",
                  subtitle = sprintf("terms grouped at kappa >= %.2f; grey = n.s.",
                                     kappa_cutoff)) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(panel.grid = ggplot2::element_blank(),
                   axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))

  width  <- max(5, 0.7 * length(full_tables) + 5)
  height <- max(3, 0.32 * nrow(reps) + 2)
  ggplot2::ggsave(paste0(file_stem, ".png"), p, width = width, height = height, dpi = 300)
  ggplot2::ggsave(paste0(file_stem, ".pdf"), p, width = width, height = height)

  clusters
}


# =============================================================================
# 8. OUTPUT
# =============================================================================

write_tsv <- function(df, path) {
  write.table(df, path, sep = "\t", quote = FALSE, row.names = FALSE, na = "")
}


#' Make a list name safe to use inside a filename.
#'
#' A name like "up/down" is a perfectly reasonable thing to call a gene list and
#' a broken file path. Plots and tables keep the original name; only the
#' filename is sanitised.
safe_name <- function(x) gsub("[^A-Za-z0-9._-]+", "_", x)


#' Write down everything needed to reproduce the run.
#'
#' GO annotations change every month, so the same script on the same list will
#' not give the same answer next year. Keep this file with your figures.
write_run_log <- function(path, settings, extra_lines = character(0)) {

  con <- file(path, "w")
  on.exit(close(con))

  writeLines(c("GO enrichment run log",
               paste("date:", Sys.time()),
               "", "--- settings ---"), con)

  for (nm in names(settings)) {
    value <- settings[[nm]]
    value <- if (is.null(value)) "NULL" else paste(value, collapse = ", ")
    writeLines(sprintf("%-22s %s", nm, value), con)
  }

  writeLines(c("", "--- annotation versions ---",
               paste("GO.db          ", as.character(utils::packageVersion("GO.db"))),
               paste("clusterProfiler", as.character(utils::packageVersion("clusterProfiler")))), con)

  if (length(extra_lines) > 0) writeLines(c("", "--- notes ---", extra_lines), con)

  writeLines(c("", "--- sessionInfo ---"), con)
  writeLines(utils::capture.output(utils::sessionInfo()), con)
}
