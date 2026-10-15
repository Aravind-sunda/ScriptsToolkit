#!/usr/bin/env Rscript
# top_go_terms.R -- top N enriched GO terms per Category x GeneList from a Metascape run.
#
# Reads Enrichment_GO/GO_AllLists.csv (long format: one row per term per gene list),
# filters to significant terms, ranks them, keeps the top N per group, and writes a
# CSV table plus a PNG dot plot.
#
# Usage:
#   Rscript top_go_terms.R --input <.../Enrichment_GO/GO_AllLists.csv> --outdir <dir> [options]
#
# Options:
#   --input FILE         GO_AllLists.csv from Metascape (required)
#   --outdir DIR         output directory, created if needed         [default: .]
#   --prefix NAME        output file name prefix                     [default: top_go_terms]
#   --top-n N            terms to keep per Category x GeneList       [default: 5]
#   --rank-by METRIC     logq | logp | enrichment | z                [default: logq]
#                          logq        Log(q-value), BH-adjusted p (most negative first)
#                          logp        LogP, raw p-value
#                          enrichment  fold enrichment (largest first)
#                          z           Z-score (largest first)
#   --q-cutoff Q         keep terms with q-value below Q; 1 = no filter   [default: 0.05]
#   --min-hits N         keep terms with at least N genes from the list   [default: 1]
#   --distinct-parent    keep only the best term per PARENT_GO within each group
#   --no-plot            write the CSV only
#   --help               show this help
#
# Ties are broken by larger Enrichment (then by Log(q-value) if ranking by enrichment).

suppressPackageStartupMessages(library(tidyverse))

# Pin the dplyr verbs used below, so packages loaded earlier or later in the same R
# session (MASS, plyr, AnnotationDbi, ...) cannot mask them with their own select/rename/etc.
select  <- dplyr::select
rename  <- dplyr::rename
filter  <- dplyr::filter
mutate  <- dplyr::mutate
arrange <- dplyr::arrange
count   <- dplyr::count
desc    <- dplyr::desc
slice_head <- dplyr::slice_head


# ============================== OPTIONS ==============================
defaults <- list(
  input           = "",
  outdir          = ".",
  prefix          = "top_go_terms",
  top_n           = 5,
  rank_by         = "logq",
  q_cutoff        = 0.05,
  min_hits        = 1,
  distinct_parent = FALSE,
  no_plot         = FALSE,
  help            = FALSE
)

# "--top-n 10" -> opts$top_n = 10; a logical default is a switch that needs no value
parse_args <- function(args, defaults) {
  opts <- defaults
  i <- 1
  while (i <= length(args)) {
    key <- args[i] %>% str_remove("^--") %>% str_replace_all("-", "_")
    if (!key %in% names(defaults)) stop("Unknown option: ", args[i], call. = FALSE)

    if (is.logical(defaults[[key]])) {
      opts[[key]] <- TRUE
      i <- i + 1
    } else {
      value <- args[i + 1]
      opts[[key]] <- if (is.numeric(defaults[[key]])) as.numeric(value) else value
      i <- i + 2
    }
  }
  opts
}

# Print the comment header of this script as the help text
show_help <- function() {
  script_path <- commandArgs() %>% str_subset("^--file=") %>% str_remove("^--file=")

  readLines(script_path) %>%
    str_subset("^#") %>%
    str_remove("^# ?") %>%
    str_subset("^!", negate = TRUE) %>%
    walk(message)
}

check_options <- function(opts) {
  if (!file.exists(opts$input))
    stop("--input not found: '", opts$input, "'", call. = FALSE)
  if (!opts$rank_by %in% c("logq", "logp", "enrichment", "z"))
    stop("--rank-by must be one of: logq, logp, enrichment, z", call. = FALSE)
  invisible(opts)
}


# ============================== STEPS ==============================
# 1. Read the Metascape table and give the columns plain names
read_go_results <- function(path) {
  read_csv(path, show_col_types = FALSE) %>%
    rename(
      category    = Category,
      go          = GO,
      description = Description,
      parent_go   = PARENT_GO,
      log_p       = LogP,
      enrichment  = Enrichment,
      z_score     = `Z-score`,
      n_hits      = `#GeneInGOAndHitList`,
      log_q       = `Log(q-value)`,
      gene_list   = GeneList,
      hits        = Hits
    ) %>%
    select(category, gene_list, go, description, parent_go,
           log_p, log_q, enrichment, z_score, n_hits, hits)
}

# 2. Keep significant terms with enough hit genes (log_q is log10 of the q-value)
filter_terms <- function(terms, q_cutoff, min_hits) {
  terms %>%
    filter(log_q < log10(q_cutoff), n_hits >= min_hits)
}

# 3. Put the chosen ranking metric in `rank_value`, where bigger always means better
add_rank_value <- function(terms, rank_by) {
  metric <- c(logq = "log_q", logp = "log_p", enrichment = "enrichment", z = "z_score")[[rank_by]]
  smaller_is_better <- rank_by %in% c("logq", "logp")

  terms %>%
    mutate(rank_value = if (smaller_is_better) -.data[[metric]] else .data[[metric]])
}

# 4. Order terms best-first within each Category x GeneList; ties -> larger enrichment
order_terms <- function(terms) {
  terms %>%
    arrange(category, gene_list, desc(rank_value), desc(enrichment), log_q)
}

# 5. Optional: one term per PARENT_GO per group (the best-ranked one).
#    Terms with no parent (e.g. all Molecular Functions) are never collapsed.
keep_distinct_parent <- function(terms) {
  terms %>%
    mutate(parent_key = if_else(is.na(parent_go) | parent_go == "", go, parent_go)) %>%
    distinct(category, gene_list, parent_key, .keep_all = TRUE) %>%
    select(-parent_key)
}

# 6. Top N per Category x GeneList, with a rank column
take_top_n <- function(terms, top_n) {
  terms %>%
    group_by(category, gene_list) %>%
    slice_head(n = top_n) %>%
    mutate(rank = row_number()) %>%
    ungroup() %>%
    relocate(rank, .after = gene_list)
}

# 7. Dot plot: x = -log10(q), colour = enrichment, size = genes in term
plot_top_terms <- function(top_terms, top_n) {
  plot_data <- top_terms %>%
    # description is not unique across panels, so the y key carries list + category
    mutate(
      y_key = str_c(description, gene_list, category, sep = "___"),
      y_key = fct_rev(fct_inorder(y_key))
    )

  ggplot(plot_data, aes(x = -log_q, y = y_key, colour = enrichment, size = n_hits)) +
    geom_point() +
    facet_grid(gene_list ~ category, scales = "free_y", space = "free_y",
               labeller = labeller(gene_list = label_wrap_gen(18))) +
    scale_y_discrete(labels = function(x) x %>% str_remove("___.*") %>% str_wrap(width = 45)) +
    scale_colour_viridis_c(option = "plasma", end = 0.9) +
    labs(x = "-log10(q-value)", y = NULL, colour = "Fold\nenrichment", size = "Genes\nin term",
         title = str_glue("Top {top_n} GO terms per gene list")) +
    theme_bw(base_size = 11) +
    theme(strip.text.y = element_text(angle = 0), panel.grid.minor = element_blank())
}

# Show how many terms each group ended up with, and flag groups below top_n
report_groups <- function(top_terms, top_n) {
  top_terms %>%
    count(category, gene_list, name = "n_terms") %>%
    mutate(note = if_else(n_terms < top_n, "fewer than top_n pass the filters", "")) %>%
    print(n = Inf)
}


# ============================== MAIN ==============================
main <- function() {
  opts <- parse_args(commandArgs(trailingOnly = TRUE), defaults)
  if (opts$help) { show_help(); quit(status = 0) }
  check_options(opts)
  dir.create(opts$outdir, recursive = TRUE, showWarnings = FALSE)

  ranked_terms <- read_go_results(opts$input) %>%
    filter_terms(opts$q_cutoff, opts$min_hits) %>%
    add_rank_value(opts$rank_by) %>%
    order_terms()

  if (opts$distinct_parent) ranked_terms <- keep_distinct_parent(ranked_terms)

  top_terms <- take_top_n(ranked_terms, opts$top_n)

  csv_path <- file.path(opts$outdir, str_c(opts$prefix, ".csv"))
  write_csv(top_terms, csv_path)
  report_groups(top_terms, opts$top_n)
  message("Wrote: ", csv_path)

  if (!opts$no_plot) {
    png_path <- file.path(opts$outdir, str_c(opts$prefix, ".png"))
    plot_height <- 2 + 0.3 * n_distinct(top_terms$gene_list) * opts$top_n

    top_terms %>%
      plot_top_terms(opts$top_n) %>%
      ggsave(png_path, plot = ., width = 13, height = plot_height, dpi = 200)
    message("Wrote: ", png_path)
  }
}

main()
