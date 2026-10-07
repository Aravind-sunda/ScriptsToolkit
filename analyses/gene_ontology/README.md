# GO enrichment

Two files:

| File | What it is | Do you edit it? |
|---|---|---|
| `run_go.R` | settings block at the top, then ~130 lines of driver | **yes** — this is the one you edit |
| `go_functions.R` | the engine: ~20 short functions, each doing one thing | no |

Keep a copy of `run_go.R` per project. `go_functions.R` stays put.

```bash
Rscript run_go.R          # from the shell
```
```r
source("run_go.R")        # from R / RStudio
```

You can also `source("go_functions.R")` inside an `.Rmd` and call the functions
yourself. Sourcing it has no side effects — it only defines functions.

## Install

```r
BiocManager::install(c("clusterProfiler", "GO.db", "GOSemSim", "org.Hs.eg.db"))
install.packages("ggplot2")
```

Swap in `org.Mm.eg.db` for mouse, etc. `GOSemSim` is only needed for
`SIMPLIFY <- TRUE`. Already installed on this machine (R 4.5.3,
clusterProfiler 4.18.4).

## The settings

Everything lives in one block at the top of `run_go.R`. Anything you don't want,
set to `NULL` or `FALSE`.

```r
GENE_LISTS <- c(
  dTERT = "lists/dTERT_targets.txt",
  dRT   = "lists/dRT_targets.txt",
  dIDR  = "lists/dIDR_targets.txt"
)
BACKGROUND <- "lists/expressed_genes.txt"
ORGANISM   <- "human"
ID_TYPE    <- NULL            # NULL = detect automatically
ONTOLOGIES <- c("BP", "MF")   # any of BP, MF, CC
```

**Gene lists.** One or many. Each value is either a file (one gene per line, or
a table whose first column is the gene — a header row is detected and dropped)
or a character vector you already have in R. These work:

```r
GENE_LISTS <- c(peaks = "idr_peak_genes.txt")                 # one file
GENE_LISTS <- list(up = up_genes, down = down_genes)          # vectors from R
GENE_LISTS <- list(clip = "targets.txt", de = res$gene[res$padj < 0.05])  # mixed
```

If you have a DE table, pull the up and down genes out yourself and pass them as
two vectors. Analysing directions separately is deliberate — pooling them loses
sensitivity (Bora et al. 2026, mistake 6).

**Background.** Every gene that *could* have ended up in your lists: genes
expressed in that cell type, genes with read coverage, every gene DESeq2 tested.
Not all 20,000 genes in the genome. The script refuses to run without one unless
you also set `ALLOW_GENOME_BACKGROUND <- TRUE`.

This is not pedantry. Wijesooriya et al. (2022) screened 186 papers: 95% of ORA
analyses either used an inappropriate background or never described it. In their
own RNA-seq example, swapping the detected-gene background for the whole genome
took the result overlap down to a Jaccard index of 0.41 and made 26 gene sets
come out simultaneously up- and down-regulated. On the test data here, the same
gene list has 1,130 testable GO:BP terms against a 6,000-gene expressed
background and 6,381 against the whole annotation.

**Gene IDs** are auto-detected (symbol / Ensembl / Entrez / RefSeq), Ensembl
version suffixes are stripped, and the mapping rate is printed and logged.

### The rest

| Setting | Default | Change it when |
|---|---|---|
| `MIN_TERM_SIZE` / `MAX_TERM_SIZE` | 10 / 500 | lower the max to 300 to push toward specificity; raise it only if you actually want broad terms |
| `MIN_GENES_HIT` | 3 | a term hit by 2 genes is rarely worth reporting |
| `PADJ_CUTOFF` | 0.05 | tighten to 0.01 if the output is huge — that is the benchmark's advice, rather than scrolling further |
| `P_ADJUST` | `"BH"` | — |
| `SIMPLIFY` | `TRUE` | costs a little time; adds a `Representative` column, deletes nothing |
| `SIMPLIFY_CUTOFF` | 0.7 | lower merges more aggressively |
| `TOP_N` | 20 | terms per plot |
| `PLOT_RANK_BY` | `"p.adjust"` | `"FoldEnrichment"` to show the strongest terms rather than the most significant ones — worth trying, ranking by p-value tends to surface shallow terms like "gene expression" |
| `COMPARE_LISTS` | `TRUE` | only does anything with >1 list |
| `KAPPA_CUTOFF` | 0.3 | Metascape's value; lower merges more |

## Output

```
PREFIX_<list>_<ONT>_all_terms.tsv     every term tested, raw + adjusted p
PREFIX_<list>_<ONT>_significant.tsv   filtered by PADJ_CUTOFF and MIN_GENES_HIT
PREFIX_<list>_<ONT>_dotplot.{png,pdf} fold enrichment, sized by gene count
PREFIX_comparison_<ONT>.{png,pdf}     cross-list heatmap (>1 list)
PREFIX_comparison_<ONT>_term_clusters.tsv  what got collapsed into what
PREFIX_summary.tsv                    per-list counts + specificity metrics
PREFIX_run_log.txt                    every setting, package versions, sessionInfo()
```

Four columns most tools don't give you:

- **FoldEnrichment** — the effect size. Reporting only p-values is mistake 3 in
  Bora et al. (2026).
- **AnnotationSize** — genes annotated to the term in your background. Large = generic.
- **Depth** — longest path to the ontology root. Deep = specific.
- **p.adjust.allterms** — see below. **This is the column to report**;
  `p.adjust` is kept beside it so you can see the difference.

`AnnotationSize` and `Depth` come from the Frontiers 2026 benchmark's
informativeness criterion, and `PREFIX_summary.tsv` reports the median of each
across your significant terms. They are a sanity check you can actually act on.
From the test data in this directory:

| list | significant BP terms | median annotation size | median depth |
|---|---|---|---|
| spliceosome genes | 40 | 105 | 5 |
| random genes | 0 | — | — |

When a list *does* return only generic output you see it immediately: median
depth 3 and median annotation size in the hundreds means you have
"biological regulation"-grade results, and should tighten `MAX_TERM_SIZE`.

## Why it's built this way

**de Oliveira, Gomes & Feltes (2026)**, *Front. Bioinform.* 6:1755664 —
the paper you pointed me at. They benchmarked 12 ORA tools (DAVID, PANTHER,
WebGestalt, Enrichr, ShinyGO, limma/goana, topGO, GOstats, clusterProfiler,
g:Profiler, ClueGO, BiNGO) on random negative controls, two positive-control
sets, and a real lung-cancer dataset, at 50/100/200/500 genes. Their headline is
that tools using the *same* method on the *same* database diverge badly, and the
divergence is about ranking and term specificity rather than which biology
appears somewhere in the output.

| Their finding | What it changed here |
|---|---|
| DAVID, Enrichr, clusterProfiler and topGO returned the smallest, deepest (most specific) terms and kept the known targets top-ranked at every list size | clusterProfiler is the engine — the only one of those four that is fully scriptable, takes a custom background, and isn't a web service |
| goana and g:Profiler returned the largest, shallowest terms | not used |
| ClueGO and goana returned significant terms from *randomised* lists even after correction | not used |
| Larger inputs disperse true targets down the ranking | the summary table reports median annotation size and depth so you can see this happening |
| They suggest filtering by annotation size and depth, and clustering for interpretability | `MIN/MAX_TERM_SIZE`, the `Depth` column, and the kappa-clustered heatmap |

**Wijesooriya, Jadaan, Perera, Kaur & Ziemann (2022)**, *PLoS Comput. Biol.*
18:e1009935 — the background requirement, and the fact that 43% of the papers
they screened didn't correct p-values at all.

**Ziemann, Schroeter & Bora (2024)**, *Bioinform. Adv.* 4:vbae159, "Two subtle
problems with overrepresentation analysis". Their second problem is the one
worth fixing: clusterProfiler only tests GO terms that contain at least one of
your genes, and corrects across just those. A term with zero overlap is still a
test that was performed. Leaving it out makes the FDR too small.

So the script counts every GO term whose size *in your background* falls within
the size limits, and re-corrects across that number using `p.adjust(n = ...)`.
That's the `p.adjust.allterms` column. On the test data it corrects across 1,130
BP terms where clusterProfiler used 588 — always more conservative, never less.
The count is verified against clusterProfiler's own internal GO data.

**Bora, McKenzie & Ziemann (2026)**, *PLoS Comput. Biol.* 22:e1014122, "Ten
common mistakes that could ruin your enrichment analysis" — the source for
reporting an effect size, ranking by it, separating directions, and logging
enough to reproduce the run.

**Zhou et al. (2019)**, *Nat. Commun.* 10:1523 (Metascape) — the cross-list
heatmap. Each list is enriched separately against the shared background; pairwise
Cohen's kappa is computed between the gene memberships of every significant term;
that matrix is hierarchically clustered and cut at 0.3; the most significant term
in each cluster becomes one heatmap row, coloured by -log10(FDR), grey where the
term wasn't significant in that list. The clustering is the part that matters —
without it a cross-list heatmap is 40 rows of the same GO branch restated. On the
test data, 139 significant terms collapse to 16 rows.

Two deliberate differences from Metascape: it enriches against the whole genome
by default and reports raw *p* < 0.01, whereas you keep your own background and
an adjusted cutoff; and cutting an average-linkage tree at kappa 0.3 approximates
their tree-trimming rather than reimplementing it. `*_term_clusters.tsv` shows
exactly what was collapsed into what.

Metascape's own docs are blunt about a caveat that applies here too: enrichment
on lists of several thousand genes washes out. Their server refuses lists over
3,000. If your lists are that big, threshold harder before comparing them.

### One thing to know about `enrichment_force_universe`

Ziemann et al. (2024)'s *first* problem is that ORA tools silently drop
background genes with no annotation, shrinking the denominator. They note
clusterProfiler added an undocumented `enrichment_force_universe = TRUE` option
to address it.

**Don't turn it on.** Tested here against DOSE 4.4.0: it replaces the background
denominator with the full universe but leaves the query denominator at the
annotated count, so the two sides of the hypergeometric test no longer match. A
random 150-gene list drawn from a 6,000-gene background returned **383
significant BP terms** with it on and **0** with it off. The script sets it to
`FALSE` explicitly, in case it's on in your `.Rprofile`. clusterProfiler's
default behaviour — annotated query against annotated background — is internally
consistent and correctly calibrated.

## How it was checked

`go_functions.R` and `run_go.R` were run on synthetic lists with known biology
planted in them (150 genes each: 60 real spliceosome or ribosome-biogenesis
genes plus 90 fillers, against a 6,000-gene background):

- the spliceosome list returns "mRNA splicing, via spliceosome" at the top, depth 12
- the ribosome list returns "ribosome biogenesis" at the top
- five independent random 150-gene lists return **0** significant terms each
- redundancy marking collapsed three nested splicing restatements to one
  representative and kept the *most specific* one
- BP, MF and CC all run; human and mouse both run; symbol, Ensembl (with and
  without version suffix), Entrez and RefSeq IDs are all detected correctly
- a missing background is refused; `ALLOW_GENOME_BACKGROUND` overrides it
- a list with no significant terms writes its tables and skips the plot without
  erroring

What has *not* been tested: rat, zebrafish, fly, worm and yeast (the code path is
identical, only the annotation package changes), and real data of your own.

## Caveats

- GO annotations change monthly. The same script on the same list will not give
  the same answer next year. `PREFIX_run_log.txt` records the GO.db and
  clusterProfiler versions — keep it with your figures.
- `SIMPLIFY` can keep a generic parent and drop the specific child. It adds a
  `Representative` column rather than deleting rows, so check what it flagged.
- ORA has lower sensitivity than functional class scoring. If you have a ranked
  statistic rather than a list, GSEA is the better tool — but note Geistlinger
  et al. (2021) found plain ORA was not significantly worse than PADOG on their
  benchmark, so a correctly-run ORA with the right background is not a
  second-class result. The background matters more than the method.
- The old `old_test/go_enrichment.R` also does GSEA and DE-table input via a
  command-line interface. This rewrite is deliberately ORA-and-gene-lists only.

## References

- de Oliveira FHS, Gomes FA, Feltes BC (2026) Benchmarking multiple gene ontology
  enrichment tools reveals high biological significance, ranking, and stringency
  heterogeneity among datasets. *Front. Bioinform.* 6:1755664.
  doi:10.3389/fbinf.2026.1755664
- Wijesooriya K, Jadaan SA, Perera KL, Kaur T, Ziemann M (2022) Urgent need for
  consistent standards in functional enrichment analysis. *PLoS Comput. Biol.*
  18:e1009935. doi:10.1371/journal.pcbi.1009935
- Ziemann M, Schroeter B, Bora A (2024) Two subtle problems with
  overrepresentation analysis. *Bioinform. Adv.* 4:vbae159.
  doi:10.1093/bioadv/vbae159
- Bora A, McKenzie M, Ziemann M (2026) Ten common mistakes that could ruin your
  enrichment analysis. *PLoS Comput. Biol.* 22:e1014122.
  doi:10.1371/journal.pcbi.1014122
- Zhou Y, Zhou B, Pache L, Chang M, Khodabakhshi AH, Tanaseichuk O, Benner C,
  Chanda SK (2019) Metascape provides a biologist-oriented resource for the
  analysis of systems-level datasets. *Nat. Commun.* 10:1523.
- Geistlinger L, Csaba G, Santarelli M, et al. (2021) Toward a gold standard for
  benchmarking gene set enrichment analysis. *Brief. Bioinform.* 22:545–556.
- Wu T, Hu E, Xu S, et al. (2021) clusterProfiler 4.0: A universal enrichment
  tool for interpreting omics data. *The Innovation* 2:100141.
