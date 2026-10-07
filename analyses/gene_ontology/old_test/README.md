# `go_enrichment.R` — general-purpose GO enrichment

One script, any gene list. Defaults are set from the current benchmarking literature rather than
from what is convenient.

## Why these tools

**de Oliveira, Gomes & Feltes (2026), *Front. Bioinform.* 6:1755664** benchmarked 12 ORA-based GO
tools (DAVID, PANTHER, WebGestalt, Enrichr, ShinyGO, limma/goana, topGO, GOstats, clusterProfiler,
g:Profiler, ClueGO, BiNGO) on random negative controls, two positive-control gene sets, and a real
DE dataset, at 50/100/200/500 genes. Their headline finding is that tools using the *same* method
and the *same* database diverge badly, and that the divergence is mostly about **ranking and term
specificity**, not about which biology is present somewhere in the output.

Findings that drove the design here:

| Finding | Consequence for this script |
|---|---|
| DAVID, Enrichr, clusterProfiler and topGO returned the smallest, deepest (i.e. most specific) GO terms and were the only four tools that kept the known target ontologies near the top of the ranking at every list size | clusterProfiler is the ORA engine; topGO is the optional cross-check. DAVID and Enrichr are web-only, so they don't fit a scriptable pipeline |
| goana (limma) and g:Profiler returned the largest, shallowest terms — "metabolic process", "response to stimulus" | not used |
| ClueGO and goana returned significant terms from *randomised* gene lists even after correction | not used |
| ShinyGO and ClueGO had the highest false-positive rates (0.8–1.0 in places), partly because they return almost nothing that isn't significant | not used |
| Most tools disperse the true targets down the ranking as input size grows | the script reports median annotation size and median depth of your significant terms, so you can see when your output has drifted generic |
| Nominal p-values give hits on random gene lists; adjusted p-values mostly don't | correction is mandatory, BH by default, and both raw and adjusted p are written out |
| topGO's `weight01` prunes redundant parent terms using the GO graph, but identified only ~2/3 of the planted targets | offered as `--topgo`, as a stringency check, never as the sole result |

**Wijesooriya, Jadaan, Kaur, Perera & Ziemann (2022), *PLoS Comput. Biol.* 18:e1009935** screened
186 papers reporting enrichment results. 95% of over-representation analyses either failed to use an appropriate background gene list or failed to describe it, and 43% did not correct p-values for multiple testing. In their own RNA-seq example, swapping the detected-genes background for the whole genome dropped the overlap between the two result sets to a Jaccard index of 0.41, and made 26 gene sets come out simultaneously up- and down-regulated. Their other recommendations:
prefer functional class scoring over ORA where possible, since ORA has lower sensitivity and is rarely used correctly, and analyse up- and down-regulated genes separately.

Hence: `--background` is required unless you explicitly override it, `--split both` is the default
for DE tables, and `--gsea` is available whenever you have a ranked statistic.

**Geistlinger et al. (2021), *Brief. Bioinform.* 22:545** is worth knowing about as the
counterweight: on their benchmark, PADOG (a permutation-based FCS method) beat GSEA, and plain ORA
was not significantly worse than PADOG. So GSEA is not automatically better than a correctly-run
ORA — the background is what matters most.

## Install

```r
BiocManager::install(c("clusterProfiler", "org.Hs.eg.db", "GO.db", "topGO", "GOSemSim"))
install.packages(c("optparse", "ggplot2"))
```
Swap `org.Hs.eg.db` for `org.Mm.eg.db` etc. as needed. `topGO` and `GOSemSim` are only needed for
`--topgo` and `--simplify`.

## Usage

**A plain gene list** (CLIP target genes, peak-associated genes, a cluster from a heatmap):

```bash
Rscript go_enrichment.R \
  --genes idr_peak_genes.txt \
  --background expressed_genes.txt \
  --organism human --ont BP \
  --simplify --topgo \
  --outdir go/TERT_full --prefix TERT_full
```

The background here should be every gene that *could* have appeared in the query — for CLIP, the
genes expressed in that cell type or the genes with any read coverage, not all 20,000 protein-coding
genes.

**A DESeq2/edgeR table** — background is derived automatically from every gene tested, up and down
are analysed separately, and GSEA runs on the full ranked list:

```bash
Rscript go_enrichment.R \
  --de-table deseq2_results.csv --gene-col gene_id \
  --lfc-col log2FoldChange --padj-col padj --stat-col stat \
  --organism mouse --ont BP --split both \
  --gsea --rank-metric stat \
  --outdir go/RPS2_STAMP --prefix RPS2
```

Gene IDs are auto-detected (symbol / Ensembl / Entrez / RefSeq), Ensembl version suffixes are
stripped, and the mapping rate is printed and logged.

## Comparing several gene lists (the Metascape heatmap)

Yes — this fits the script well, and it's now built in. Pass several files to `--genes` with
`--heatmap`:

```bash
Rscript go_enrichment.R \
  --genes dTERT.txt,dRT.txt,dIDR.txt,HA-GFP.txt \
  --labels dTERT,dRT,dIDR,GFP \
  --background expressed_genes.txt \
  --organism human --ont BP --heatmap \
  --outdir go/domain_compare --prefix domains
```

It also fires automatically for `--split both` on a DE table, giving you an up-vs-down heatmap.

What it reproduces from Metascape (Zhou et al. 2019, *Nat. Commun.* 10:1523): each list is enriched
separately against a shared background; pairwise **Kappa** similarity is computed between the gene
memberships of every significant term; the similarity matrix is hierarchically clustered and cut at
0.3; the most significant term in each cluster represents that cluster as one heatmap row, coloured
by -log10(FDR) with grey for "not significant in this list". That clustering step is the part that
matters — without it a cross-list heatmap is 40 rows of the same GO branch restated.

What differs, deliberately:

- **Statistics stay the ones chosen for the rest of the script.** Metascape's public server enriches
  against the whole genome by default and reports terms at raw *p* < 0.01 with enrichment factor
  > 1.5; here you keep your own background and a BH-adjusted cutoff.
- **Cutting an average-linkage tree at kappa 0.3** is an approximation of Metascape's tree-trimming,
  not a reimplementation of it. Cluster membership is written to `*_term_clusters.tsv` so you can see
  exactly what was collapsed into what.
- No Circos/MCODE/PPI layers. If you want the enrichment network with pie-chart nodes, export
  `*_term_clusters.tsv` and the kappa matrix into Cytoscape.

One caution the Metascape docs are blunt about and which applies here too: enrichment on lists of
several thousand genes washes out. Their server refuses lists over 3000 and they suggest keeping to
a few hundred. If your CLIP target lists are that large, threshold harder before comparing them.

## Which GSEA ranking metric

Short answer: **the moderated test statistic** — DESeq2's `stat` (Wald), limma-voom's `t`. That is
now what `--rank-metric auto` picks, falling back to signed -log10(p) only when no such column
exists.

Reasoning:

- **Zyla et al. (2017), *BMC Bioinformatics* 18:256** tested 16 ranking metrics on 28 benchmark
  datasets and found the choice materially changes results — a default metric can give poor output.
  Their best performers were the moderated Welch test statistic and signal-to-noise, both of which
  shrink the variance estimate. DESeq2's `stat` and limma's `t` are the practical members of that
  family for RNA-seq. Note their top-ranked metrics were the *absolute* versions, which detect sets
  changing in either direction but throw away the sign — only useful if you aren't going to
  interpret NES direction, which usually you are.
- **Raw log2FC is the weakest common choice.** It ignores the variance entirely, so a low-count gene
  with a 5-fold swing outranks a well-measured 1.5-fold change. Use it only if the column already
  holds shrunken LFCs (apeglm/ashr), which restores most of what plain LFC loses.
- **Signed -log10(p) has a ties problem.** DESeq2 returns p-values that underflow to 0 for the
  strongest genes, so the top of the ranking becomes a block of ties at `Inf`, exactly where the
  running enrichment statistic is most sensitive. The script clamps to the smallest representable
  double and warns when more than 2% of ranks are tied.
- **edgeR** returns an unsigned statistic (`F` from `glmQLFTest`, `LR` from `glmLRT`), handled
  below.

### edgeR tables specifically

edgeR never prints a signed moderated statistic, but you can recover one exactly. A QL F-test on a
single coefficient or contrast has numerator df = 1, so `F` is the square of the moderated *t*, and
the LRT statistic is chi-square on 1 df, so `LR` is the square of a *z*. Re-attaching the sign of
`logFC` gives back the statistic you want:

```
rank = sign(logFC) * sqrt(F)      # glmQLFTest
rank = sign(logFC) * sqrt(LR)     # glmLRT
```

`--rank-metric auto` now detects an `F` or `LR` column and does this for you. Export the table first,
since edgeR keeps gene IDs in the rownames:

```r
tt <- topTags(qlf, n = Inf)$table
tt$gene_id <- rownames(tt)
write.csv(tt[, c("gene_id", setdiff(names(tt), "gene_id"))],
          "edger_TE_results.csv", row.names = FALSE)
```

```bash
Rscript go_enrichment.R \
  --de-table edger_TE_results.csv --gene-col gene_id \
  --lfc-col logFC --padj-col FDR --stat-col F \
  --organism human --ont BP --split both --gsea \
  --outdir go/TE --prefix TE
```

Note the column names: edgeR uses `logFC`, `PValue` and `FDR`, not DESeq2's `log2FoldChange`,
`pvalue` and `padj`, so `--lfc-col` and `--padj-col` have to be set explicitly. The raw p-value
column is auto-detected.

The sign trick is only valid when the test has one degree of freedom in the numerator. If you ran an
ANOVA-style test across several coefficients, `F` is unsigned for a real reason and `logFC` no longer
summarises the contrast — rank on a specific pairwise contrast instead.

**For an offset-based model** (translation efficiency from edit counts with `log(RNA counts)` as a
second offset, for instance) two things carry over into the enrichment step. The `logFC` is a TE fold
change, not an expression fold change, so NES sign means "translated more per unit mRNA" — worth
spelling out in the figure legend, because readers will assume expression. And the background must be
the genes that survived your coverage filter and actually entered the model, not all expressed genes;
with `--de-table` the script takes care of that automatically, since every row of the table is by
definition a gene that was tested.

Two further points worth knowing. Candia & Ferrucci (2024) found gene-set permutation gave the best
sensitivity/specificity trade-off among GSEA modalities (AUC ≈ 0.99), which is what `fgsea`/`gseGO`
does by default. And Stead et al. (2025) showed GSEA readily returns significant gene sets from data
with no differentially expressed genes at all — so if your DE analysis is empty, a full GSEA hit
list is a warning sign, not a rescue.

## Output

```
PREFIX_up_all_terms.tsv              every term tested, raw + adjusted p
PREFIX_up_significant.tsv            filtered by --padj and --min-count
PREFIX_up_significant_simplified.tsv redundancy-reduced (--simplify)
PREFIX_up_dotplot.{png,pdf}          fold enrichment, sized by gene count
PREFIX_gsea_*.tsv                    GSEA results (--gsea)
PREFIX_summary.tsv                   per-set counts + specificity metrics
PREFIX_run_log.txt                   every parameter, package version, sessionInfo()
```

Two columns worth reading that most tools don't give you:

- **AnnotationSize** — genes annotated to the term in your background. Large = generic.
- **Depth** — longest path to the ontology root. Deep = specific.

Both come straight from the benchmark's informativeness criterion. If your top 20 terms have a
median annotation size in the thousands and a depth of 3, you have "biological regulation"-grade
output and should tighten `--max-size` or use `--topgo`.

## Defaults and when to change them

| Flag | Default | Change it when |
|---|---|---|
| `--min-size` / `--max-size` | 10 / 500 | raise `--max-size` only if you specifically want broad terms; lower to 300 to push toward specificity |
| `--min-count` | 3 | a term hit by 2 genes is rarely worth reporting |
| `--padj` | 0.05 | the benchmark's advice for large outputs is a stricter cutoff, not more browsing |
| `--p-method` | BH | — |
| `--split` | both | `all` only if direction is meaningless for your list |
| `--heatmap` | off | on whenever you analyse >1 list; free for `--split both` |
| `--kappa` | 0.3 | Metascape's value; lower merges more aggressively |
| `--rank-metric` | auto | override only if you know your table has a better column |
| `--simplify` | off | on when the significant table is long and visibly redundant; check what it removed, since redundancy reduction can drop the specific term and keep the generic parent |
| `--topgo` | off | on when you want a conservative second opinion; treat terms significant in both as your confident set |

## Caveats

- I wrote this against the clusterProfiler/topGO APIs but had no R available to execute it, so run
  it once on a small list first. If anything breaks it will most likely be a package-version
  argument mismatch (`eps` in `gseGO` is already guarded).
- `--simplify` needs `GOSemSim` and only works with a single ontology (`--ont BP`, not `ALL`).
- topGO's `weight01` p-values are not independent across terms, so the FDR column for it is
  approximate — the benchmark applied `p.adjust` the same way, but treat it as a ranking aid.
- Record the `run_log.txt` alongside your figures. GO annotations change monthly, and the same
  script on the same list will not give the same answer next year.

## References

- de Oliveira FHS, Gomes FA, Feltes BC (2026) Benchmarking multiple gene ontology enrichment tools
  reveals high biological significance, ranking, and stringency heterogeneity among datasets.
  *Front. Bioinform.* 6:1755664. doi:10.3389/fbinf.2026.1755664
- Wijesooriya K, Jadaan SA, Perera KL, Kaur T, Ziemann M (2022) Urgent need for consistent standards
  in functional enrichment analysis. *PLoS Comput. Biol.* 18:e1009935.
- Geistlinger L et al. (2021) Toward a gold standard for benchmarking gene set enrichment analysis.
  *Brief. Bioinform.* 22:545–556.
- Candia J, Ferrucci L (2024) Assessment of Gene Set Enrichment Analysis using curated RNA-seq-based
  benchmarks. *PLoS ONE* 19:e0302696.
- Zhou Y, Zhou B, Pache L, Chang M, Khodabakhshi AH, Tanaseichuk O, Benner C, Chanda SK (2019)
  Metascape provides a biologist-oriented resource for the analysis of systems-level datasets.
  *Nat. Commun.* 10:1523.
- Zyla J, Marczyk M, Weiner J, Polanska J (2017) Ranking metrics in gene set enrichment analysis:
  do they matter? *BMC Bioinformatics* 18:256.
- Stead JDH et al. (2025) Gene Set Enrichment Analysis in zebrafish embryos is susceptible to
  false-positive results in the absence of differentially expressed genes.
  *Bioinform. Biol. Insights* 19.
- Klopfenstein DV et al. (2018) GOATOOLS: a Python library for Gene Ontology analyses.
  *Sci. Rep.* 8:10872. (Python alternative if you'd rather not use R.)
