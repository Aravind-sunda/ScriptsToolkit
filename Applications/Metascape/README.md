# Metascape MSBio runner (`run_metascape.sh`)

Runs Metascape's command-line tool (MSBio, `msbio2.sif`) on a multi-gene-list file
and writes the full Metascape report (GO, heatmap, PPI, QC, circos, Excel, HTML, PPTX).

Metascape docs: https://metascape.org/gp/index.html#/menu/msbio
(The page is JavaScript-only, so it can't be fetched from the terminal. Options below
were taken from `msbio2 --help`, the vendor `bin/sms.sh`, and the example files.)

## Quick start

```bash
cd /condo/brannanlab/tmhaxs421/scripts/ScriptsToolkit/Applications/Metascape

./run_metascape.sh -i /path/to/multi_gene_list.txt -o /path/to/output_dir
```

Results end up in `output_dir/results/<input file name>/`. Open `AnalysisReport.html` first.

Example (ACLY Venn regions):

```bash
./run_metascape.sh \
  -i /condo/brannanlab/tmhaxs421/Bulk_RNAseq/ACLY/analysis/ACLY_AC_vs_NC_CLIP_venn_regions_multi_gene_list.txt \
  -o /condo/brannanlab/tmhaxs421/Bulk_RNAseq/ACLY/metascape
```

## Input format

Tab-separated, one gene list per row: `Name<TAB>gene,gene,gene`. Lines starting with `#`
are ignored. `write_multi_gene_list()` in the analysis Rmds writes exactly this.
List names should avoid `/` (the Rmd replaces it with `_AND_`).

For a single gene list (not multi-list) use `-u`.

## Files and layout

```
Metascape/
├── run_metascape.sh          # this runner
├── README.md
├── instructions.sh           # how the image/zip was installed
├── msbio_v3.5.20260201.zip   # vendor download
└── metascape/                # unzipped vendor install (gitignored)
    ├── msbio2.sif            # the container image (3.1 GB)
    ├── license/license.txt   # Metascape license
    ├── bin/sms.sh            # vendor singularity launcher (what run_metascape.sh is modelled on)
    └── data/example/         # example inputs, test.job, option.json
```

## How the script works

One CONFIG block at the top, then small single-purpose functions, then `main`.

| Step | Function | What it does |
|---|---|---|
| 1 | CONFIG block | All defaults: paths, species, CPUs, `EXTRA_FLAGS`. Edit here for a new project, or override with flags |
| 2 | `parse_args` | Reads `-i -o -s -t -O -c -u -n -h` and overwrites the config values |
| 3 | `check_inputs` | Fails early if the input, sif, license or option file is missing |
| 4 | `stage_files` | Copies the input (and option.json) into `outdir/input/`, creates `results/` and `logs/` |
| 5 | `build_msbio_args` | Builds the msbio2 arguments using container paths (`/data/...`); optional flags are added only if set |
| 6 | `run_msbio` | Prints and runs the `apptainer run` command (same call as vendor `sms.sh`) |
| 7 | `main` | Runs the steps in order |

**Why the staging step exists:** the container only sees folders that are bound into it.
The script binds `outdir` as `/data` and the license folder as `/license`, so the input
file can live anywhere. Inside the container, paths are `/data/input/...` and
`/data/results/...`.

## Options

| Flag | Meaning | Default |
|---|---|---|
| `-i FILE` | input gene list file (TXT/CSV/XLSX) | required |
| `-o DIR` | output directory (created if needed) | required |
| `-s TAXID` | species of the input genes (9606 human, 10090 mouse, ...) | 9606 |
| `-t TYPE` | gene id type: Entrez, RefSeq, Symbol, dbxref | auto-detect |
| `-O FILE` | option.json to change analysis settings | Metascape defaults |
| `-c N` | CPUs | 4 |
| `-u` | input is a single list, not a multi-list file | off |
| `-n` | dry run: print the command only | off |
| `-h` | help | |

Config-only settings (edit the CONFIG block):
- `TARGET_TAX_ID`: analysis species if different from the source (maps genes via orthologs).
- `EXTRA_FLAGS`: passed straight to msbio2. Useful ones:

| Flag | Skips |
|---|---|
| `-P` | Cytoscape plotting (much faster) |
| `-G` | GO analysis |
| `-p` | PPI analysis |
| `-E` | GPEC prioritization |
| `-Q` | QC sets |
| `-r`, `-w` | list colours / weights |

Full list: run `msbio2 --help` inside the image (see Troubleshooting).

## Custom analysis settings (option.json)

Copy `metascape/data/example/option.json` and edit, then pass with `-O`. The file allows
`#` comments. Commonly changed:
- `S_go_category`: which databases to use (e.g. 19 GO BP, 23 Hallmark, 11 Canonical, 6 Reactome, 24 KEGG)
- `p_cutoff`, `min_overlap`, `min_enrichment`
- `l_PPI`, `l_QC`, `l_GO`: turn analyses on/off
- `S_ann_id`: which gene annotation columns appear in the Excel output

## Notes and gotchas

- **License:** `metascape/license/license.txt` is tied to aravind.sunda@gmail.com, expires
  2027-10-07, non-commercial. Replace the file when it expires.
- **Exit message:** `INFO: Cleanup error ... fuse-overlayfs exited` is printed at the end of every
  run and is harmless.
- **"Deleting ... .pdf" lines** in the log are Metascape cleaning intermediate PPI plots.
- **Species:** the default is human. Mouse data needs `-s 10090`.
- **Gene id auto-detect:** works by omitting `-t`. If a list is mostly Ensembl IDs or a
  mix, check the "unmapped" counts in the report before trusting the enrichment.
- **Small lists:** Metascape needs at least 3 genes per list for GO (`min_list_size_for_go`)
  and caps GO at 3000 genes (`max_list_size_for_go`).
- **Tested:** default multi-list run on `data/example/multiple_list_symbol.txt` completed with
  the full report. `-O`, `-t`, `-u` and batch scheduler submission were not tested.

## Troubleshooting

Print msbio2's own help:

```bash
cd metascape
apptainer exec --writable-tmpfs --bind $PWD/data:/data --bind $PWD/license:/license msbio2.sif \
  sh -c 'export HOME=/tmp MPLCONFIGDIR=/tmp; cd /msbio && ./mylib/ms/run_msbio2.py --help'
```

Check the image is intact: `apptainer sif list metascape/msbio2.sif` (should show a Squashfs
object of about 3.1 GB). `apptainer verify` will say "signature not found" because the
image is unsigned; that does not mean it's corrupt.

Interactive shell in the container: `metascape/bin/debug_sms.sh`.
