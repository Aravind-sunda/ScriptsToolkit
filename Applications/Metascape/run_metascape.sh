#!/usr/bin/env bash
# run_metascape.sh -- run Metascape MSBio (msbio2.sif) on a multi-gene-list file.
#
# Skeleton for reuse: edit the CONFIG block, or override anything from the command line.
#
# Usage:
#   run_metascape.sh -i multi_gene_list.txt -o /path/to/output_dir [options]
#
# Options (all override the CONFIG defaults below):
#   -i FILE   input gene list file (TXT/CSV/XLSX; multi-list = Name<TAB>Genes per line)
#   -o DIR    output directory (created if needed; results go to DIR/results/<input name>)
#   -s TAXID  source species tax id of the input genes (9606 human, 10090 mouse, ...)
#   -t TYPE   gene id type: Entrez,RefSeq,Symbol,dbxref. Default = auto-detect
#   -O FILE   option.json to change analysis settings (see data/example/option.json)
#   -c N      CPUs
#   -u        input is a single gene list (not a multi-list file)
#   -n        dry run: print the apptainer command, do not run
#   -h        show this help
#
# Notes:
#   * The container only sees directories bound to it. The input and option files are
#     copied into DIR/input/ and DIR is bound as /data, so any paths work for -i / -O.
#   * Metascape's docs: https://metascape.org/gp/index.html#/menu/msbio

set -euo pipefail

# ============================== CONFIG ==============================
METASCAPE_DIR="/condo/brannanlab/tmhaxs421/scripts/ScriptsToolkit/Applications/Metascape/metascape"
SIF="${METASCAPE_DIR}/msbio2.sif"
LICENSE_DIR="${METASCAPE_DIR}/license"       # must contain license.txt

INPUT=""                  # -i
OUTDIR=""                 # -o
SOURCE_TAX_ID="9606"      # -s   species of the input genes
TARGET_TAX_ID=""          # -T   leave empty = same as source (analysis species)
ID_TYPE=""                # -t   empty = auto-detect gene id type
OPTION_JSON=""            # -O   empty = Metascape defaults
CPU=4                     # -c
ONE_LIST=false            # -u

# Extra flags passed straight to msbio2 (see `msbio2 --help`), e.g.
#   -P skip Cytoscape plotting   -G skip GO   -p skip PPI   -E skip GPEC   -Q skip QC
EXTRA_FLAGS=()
DRY_RUN=false
# ====================================================================

usage() { sed -n '2,/^set -euo/p' "$0" | sed '$d' | sed 's/^# \{0,1\}//'; }

parse_args() {
  while getopts "i:o:s:t:O:c:unh" opt; do
    case $opt in
      i) INPUT=$OPTARG ;;
      o) OUTDIR=$OPTARG ;;
      s) SOURCE_TAX_ID=$OPTARG ;;
      t) ID_TYPE=$OPTARG ;;
      O) OPTION_JSON=$OPTARG ;;
      c) CPU=$OPTARG ;;
      u) ONE_LIST=true ;;
      n) DRY_RUN=true ;;
      h) usage; exit 0 ;;
      *) usage; exit 1 ;;
    esac
  done
}

check_inputs() {
  [[ -n $INPUT && -n $OUTDIR ]] || { echo "ERROR: -i and -o are required" >&2; usage; exit 1; }
  [[ -f $INPUT ]]               || { echo "ERROR: input not found: $INPUT" >&2; exit 1; }
  [[ -f $SIF ]]                 || { echo "ERROR: sif not found: $SIF" >&2; exit 1; }
  [[ -f $LICENSE_DIR/license.txt ]] || { echo "ERROR: no license.txt in $LICENSE_DIR" >&2; exit 1; }
  [[ -z $OPTION_JSON || -f $OPTION_JSON ]] || { echo "ERROR: option file not found: $OPTION_JSON" >&2; exit 1; }
}

# Copy input (and option file) into OUTDIR/input so they are visible inside the container as /data/input
stage_files() {
  mkdir -p "$OUTDIR/input" "$OUTDIR/results" "$OUTDIR/logs"
  chmod a+rwx "$OUTDIR" "$OUTDIR/input" "$OUTDIR/results" "$OUTDIR/logs"
  cp -f "$INPUT" "$OUTDIR/input/"
  [[ -n $OPTION_JSON ]] && cp -f "$OPTION_JSON" "$OUTDIR/input/"
  return 0
}

# Build the msbio2 argument list (paths as seen INSIDE the container)
build_msbio_args() {
  local name; name=$(basename "$INPUT")
  MSBIO_ARGS=("/data/input/$name" -o "/data/results/${name%.*}" -S "$SOURCE_TAX_ID" -c "$CPU")
  [[ -n $TARGET_TAX_ID ]] && MSBIO_ARGS+=(-T "$TARGET_TAX_ID")
  [[ -n $ID_TYPE ]]       && MSBIO_ARGS+=(-t "$ID_TYPE")
  [[ -n $OPTION_JSON ]]   && MSBIO_ARGS+=(--option "/data/input/$(basename "$OPTION_JSON")")
  $ONE_LIST               && MSBIO_ARGS+=(-u)
  MSBIO_ARGS+=("${EXTRA_FLAGS[@]}")
}

# Same invocation as the vendor bin/sms.sh, with our own data dir
run_msbio() {
  local cmd=(apptainer run --writable-tmpfs
             --bind "$OUTDIR:/data" --bind "$LICENSE_DIR:/license"
             "$SIF" /msbio/mylib/ms/smsbio2.sh "${MSBIO_ARGS[@]}")
  echo "Command: ${cmd[*]}"
  $DRY_RUN && return 0
  "${cmd[@]}"
}

main() {
  parse_args "$@"
  check_inputs
  OUTDIR=$(realpath -m "$OUTDIR")
  stage_files
  build_msbio_args
  run_msbio
  echo "Done. Results: $OUTDIR/results"
}

main "$@"
