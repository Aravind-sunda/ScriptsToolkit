# Synthetic test datasets

Gene lists whose correct answer is known before you run anything. Use them to
check the script works, to see what good and bad output look like side by side,
or as a template if you want to build controls for your own organism.

Human gene symbols, one per line. Regenerate with `Rscript make_test_datasets.R`
(seed is fixed, you get the same files back).

## The files

### Background

| File | Genes | What it is |
|---|---|---|
| `background_human_6000genes.txt` | 6000 | Stands in for "genes expressed in this experiment". Contains a realistic slice of each planted theme plus ~5600 unrelated genes. **Every human list below is a subset of this** — use it as `BACKGROUND`. |
| `background_mouse_5100genes.txt` | 5100 | Same idea for mouse. |

### Positive controls — you know what must come out

| File | Genes | Planted signal | Expected top term |
|---|---|---|---|
| `positive_splicing_60planted_of_150.txt` | 150 | 60 from GO:0000398 *mRNA splicing, via spliceosome* + 90 unrelated | RNA splicing / spliceosome, **deep** |
| `positive_ribosome_60planted_of_150.txt` | 150 | 60 from GO:0042254 *ribosome biogenesis* + 90 unrelated | ribosome biogenesis |
| `positive_mixed_splicing30_cellcycle40_of_150.txt` | 150 | 30 splicing + 40 from GO:0000278 *mitotic cell cycle* + 80 unrelated | both themes — this is the one that gives the cross-list heatmap something to separate |
| `positive_broad_developmental_shallow_terms_150.txt` | 150 | 90 from GO:0032502 *developmental process*, an enormous shallow term | "developmental process", **shallow** — this is what uninformative output looks like |
| `positive_mouse_ribosome_50planted_of_150.txt` | 150 | 50 mouse ribosome-biogenesis genes + 100 unrelated | ribosome biogenesis (set `ORGANISM <- "mouse"`) |

40% signal is roughly what a decent real gene list looks like. The
`broad_developmental` one is deliberately 60%, because broad terms give weak
fold enrichment and need more genes to show up.

### Negative controls — nothing should come out

| File | Genes | What it is |
|---|---|---|
| `negative_random_from_background_150_rep1.txt` (also `_rep2`, `_rep3`) | 150 | Drawn at random straight from the background. This is the null hypothesis of the hypergeometric test made literal, so a correct method must return ~nothing. Three replicates because one draw is not evidence. |
| `negative_themes_removed_150.txt` | 150 | Drawn from the background *after* removing all three planted themes. A stricter null — it checks that carving out three large GO branches doesn't create artefactual enrichment in what remains. It doesn't. |

**The negative controls are the ones that catch real bugs.** They are how the
miscalibration in clusterProfiler's `enrichment_force_universe` option was found
(see `../README.md`). A positive control that works tells you much less — a
badly broken method will still return splicing terms for a list of splicing
genes.

## Verified results

`Rscript example_run_go.R` runs everything against the shared background.
Actual output, GO:BP, FDR < 0.05, 1168 testable terms:

| List | Sig. terms | Top term | Top-term depth | Median depth of top 20 | Median fold enr. |
|---|---|---|---|---|---|
| splicing | 45 | RNA splicing, via transesterification reactions | 10 | 8.0 | 6.5 |
| ribosome | 24 | ribosome biogenesis | 5 | 6.0 | — |
| mixed | 82 | mitotic cell cycle | 3 | 5.5 | — |
| broad_generic | 268 | developmental process | 1 | 3.0 | 2.7 |
| **neg_random** | **0** | — | — | — | — |
| **neg_depleted** | **0** | — | — | — | — |

Mouse, run separately against `background_mouse_5100genes.txt` (1071 testable
terms): 41 significant terms, top term *ribosome biogenesis*, depth 5, 9.2×
enriched.

Read the top and bottom rows together — that contrast is the whole point:

- **Both negative controls return exactly 0.** All three random replicates do.
- **`Depth` separates good output from generic output; `AnnotationSize` mostly
  doesn't.** Splicing has median depth 8.0 in its top 20, broad_generic has 3.0.
  But their median annotation sizes are 140 and 130 — nearly identical, because
  `MAX_TERM_SIZE = 500` has already truncated the genuinely huge terms. If you
  only check one specificity column, check **Depth**.
- **Generic input produces a flood.** 268 significant terms versus 45, at less
  than half the fold enrichment. A very long result list is itself a warning
  sign, not a richer result.

## Using them

```bash
Rscript example_run_go.R      # runs all six human lists, writes example_output/
```

`example_run_go.R` is just `run_go.R` with the settings block pointed at these
files. To test the mouse path, set:

```r
GENE_LISTS <- c(mouse_ribosome = "positive_mouse_ribosome_50planted_of_150.txt")
BACKGROUND <- "background_mouse_5100genes.txt"
ORGANISM   <- "mouse"
```

To see the `enrichment_force_universe` failure for yourself, add
`options(enrichment_force_universe = TRUE)` after the `source()` line in
`example_run_go.R` and re-run. `neg_random` goes from 0 significant terms to 168.
