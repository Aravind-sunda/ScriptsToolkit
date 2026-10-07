# =============================================================================
# make_test_datasets.R -- regenerate the synthetic GO test data.
#
# Run this only if you want to rebuild the .txt files. They are already in this
# directory; the seed is fixed so you will get the same ones back.
#
#   Rscript make_test_datasets.R
#
# The idea: build gene lists whose correct answer is known in advance, so you
# can tell whether an enrichment script is working rather than guessing.
#
#   POSITIVE control -- real genes from one GO term, diluted with unrelated
#                       genes. You know which term must come out on top.
#   NEGATIVE control -- genes drawn at random FROM THE BACKGROUND. This is the
#                       null hypothesis of the hypergeometric test made literal,
#                       so a correct method must return ~nothing. Anything it
#                       does return is a false positive.
#
# The negative controls are the ones that catch real bugs. They are how the
# miscalibration in clusterProfiler's `enrichment_force_universe` option was
# found (see ../README.md).
# =============================================================================

suppressPackageStartupMessages({
  library(org.Hs.eg.db)
  library(org.Mm.eg.db)
  library(AnnotationDbi)
})

set.seed(42)
outdir <- dirname(sub("^--file=", "", commandArgs(FALSE)[grep("^--file=", commandArgs(FALSE))]))
if (length(outdir) == 0 || !nzchar(outdir)) outdir <- "."


#' All gene symbols annotated to a GO term, including its descendants.
#' GOALL does the ancestor propagation, so this is every gene the term covers.
genes_in_term <- function(go_id, orgdb) {
  hits <- suppressMessages(AnnotationDbi::select(
    orgdb, keys = go_id, keytype = "GOALL", columns = "SYMBOL"))
  unique(na.omit(hits$SYMBOL))
}

save_list <- function(genes, filename) {
  writeLines(genes, file.path(outdir, filename))
  message(sprintf("  %-52s %4d genes", filename, length(genes)))
}


# --- human ------------------------------------------------------------------
# Three unrelated themes, so a list built from one of them cannot be confused
# with a list built from another.

splicing  <- genes_in_term("GO:0000398", org.Hs.eg.db)  # mRNA splicing, via spliceosome
ribosome  <- genes_in_term("GO:0042254", org.Hs.eg.db)  # ribosome biogenesis
cellcycle <- genes_in_term("GO:0000278", org.Hs.eg.db)  # mitotic cell cycle

all_human <- unique(na.omit(AnnotationDbi::keys(org.Hs.eg.db, keytype = "SYMBOL")))
unrelated <- setdiff(all_human, c(splicing, ribosome, cellcycle))

# The background stands in for "genes expressed in this experiment": ~6000 genes,
# containing a realistic slice of each theme plus a lot of unrelated genes.
themed_part <- unique(c(sample(splicing,  120),
                        sample(ribosome,  120),
                        sample(cellcycle, 150)))
background <- unique(c(themed_part,
                       sample(unrelated, 6000 - length(themed_part))))

message("\nHuman:")
save_list(background, "background_human_6000genes.txt")

# Positive controls. 60 genes of known biology diluted into 90 unrelated genes
# -- a 40% signal, which is roughly what a decent real gene list looks like.
save_list(unique(c(sample(intersect(splicing, background), 60),
                   sample(intersect(unrelated, background), 90))),
          "positive_splicing_60planted_of_150.txt")

save_list(unique(c(sample(intersect(ribosome, background), 60),
                   sample(intersect(unrelated, background), 90))),
          "positive_ribosome_60planted_of_150.txt")

# Two themes in one list, so the cross-list heatmap has something to separate.
save_list(unique(c(sample(intersect(splicing,  background), 30),
                   sample(intersect(cellcycle, background), 40),
                   sample(intersect(unrelated, background), 80))),
          "positive_mixed_splicing30_cellcycle40_of_150.txt")

# Negative controls: drawn straight from the background, no structure at all.
# Three replicates, because one draw is not evidence.
for (i in 1:3) {
  save_list(sample(background, 150),
            sprintf("negative_random_from_background_150_rep%d.txt", i))
}

# A second negative control, drawn from the background AFTER removing all three
# themes. The worry was that carving out three large GO branches would bias what
# remains and create artefactual enrichment. Tested: it does not -- this returns
# 0 terms too. Kept as a stricter null than the plain random draw.
save_list(sample(intersect(unrelated, background), 150),
          "negative_themes_removed_150.txt")

# A list whose correct answer is a GENERIC one. "developmental process" is an
# enormous, shallow GO term, so a list built from it should return terms with
# large AnnotationSize and small Depth. This is what uninformative output looks
# like, and it is the dataset to check the two specificity columns against --
# contrast its medians with the splicing list.
broad <- genes_in_term("GO:0032502", org.Hs.eg.db)   # developmental process
save_list(unique(c(sample(intersect(broad, background), 90),
                   sample(intersect(unrelated, background), 60))),
          "positive_broad_developmental_shallow_terms_150.txt")


# --- mouse ------------------------------------------------------------------
# Same construction, to check the organism switch works.

mm_ribosome  <- genes_in_term("GO:0042254", org.Mm.eg.db)
all_mouse    <- unique(na.omit(AnnotationDbi::keys(org.Mm.eg.db, keytype = "SYMBOL")))
mm_unrelated <- setdiff(all_mouse, mm_ribosome)

mm_background <- unique(c(sample(mm_ribosome, 100), sample(mm_unrelated, 5000)))

message("\nMouse:")
save_list(mm_background, "background_mouse_5100genes.txt")
save_list(unique(c(sample(intersect(mm_ribosome,  mm_background), 50),
                   sample(intersect(mm_unrelated, mm_background), 100))),
          "positive_mouse_ribosome_50planted_of_150.txt")

message("\nDone.\n")
