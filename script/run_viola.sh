#!/bin/bash

# ================================
# VIOLA: Run the 3-step pipeline
# ================================

# ---- Default values ----
OUTPUT_DIR="results/viola_run"

# ---- Help message ----
usage() {
  echo "Usage: $0 -v VEP_INPUT -c CADD_INPUT -h HPO_TABLE -f UNIQUE_VARIANTS_VCF -r RESOURCES [-o OUTPUT_DIR]"
  echo
  echo "Arguments:"
  echo "  -v, --vep        Path to VEP input file"
  echo "  -c, --cadd       Path to CADD input file"
  echo "  -h, --hpo        Path to patient HPO table"
  echo "  -f, --vcf        Path to VCF file of unique rare variants"
  echo "  -r, --resources  Path to resources folder (contains ClinVar file)"
  echo "  -o, --out        Output directory (default: results/viola_run)"
  echo "  --help           Show this help message"
  echo
  exit 1
}

# ---- Parse arguments ----
# $# = number of arguments
# -gt 0 = greater than 0
while [[ "$#" -gt 0 ]]; do
    case $1 in
        -v|--vep) VEP_INPUT="$2"; shift ;;
        -c|--cadd) CADD_INPUT="$2"; shift ;;
        -h|--hpo) HPO_TABLE="$2"; shift ;;
        -f|--vcf) UNIQUE_VARIANTS_VCF="$2"; shift ;;
        -r|--resources) RESOURCES="$2"; shift ;;
        -o|--out) OUTPUT_DIR="$2"; shift ;;
        --help) usage ;;
        *) echo "Unknown parameter: $1"; usage ;;
    esac
    shift
done

# ---- Check required arguments ----
if [[ -z "$VEP_INPUT" || -z "$CADD_INPUT" || -z "$HPO_TABLE" || -z "$UNIQUE_VARIANTS_VCF" || -z "$RESOURCES" ]]; then
  echo "Error: Missing required arguments"
  usage
fi

# ---- Create output directory ----
mkdir -p "$OUTPUT_DIR"

# ---- Step 1: Merge annotations ----
echo "Running Step 1: merge annotations..."
# ClinVar file is assumed to be in resources
CLINVAR_INPUT="$RESOURCES/clinvar_210125_hg38_cleaned.tsv"

Rscript ./viola_step1_merge.R \
  -v "$VEP_INPUT" \
  -c "$CADD_INPUT" \
  -l "$CLINVAR_INPUT" \
  -o "$OUTPUT_DIR"

STEP1_OUT=$(ls "$OUTPUT_DIR"/*_rare_variant_cadd_input_viola.csv 2>/dev/null | head -n 1)
if [[ -z "$STEP1_OUT" ]]; then
  echo "Error: Step 1 output file not found!"
  exit 1
fi


# ---- Step 2: VAE + DBSCAN clustering ----
echo "Running Step 2: clustering..."
python ./viola_step2_cluster.py \
  -f "$STEP1_OUT" \
  -o "$OUTPUT_DIR"

STEP2_OUT=$(ls "$OUTPUT_DIR"/*_res_dbscan.csv 2>/dev/null | head -n 1)
if [[ -z "$STEP2_OUT" ]]; then
  echo "Error: Step 2 output file not found!"
  exit 1
fi

# ---- Step 3: Filtering + Ranking ----
echo "Running Step 3: filtering & ranking..."
Rscript ./viola_step3_rank.R \
  -f "$STEP2_OUT" \
  -o "$OUTPUT_DIR" \
  -t "$HPO_TABLE" \
  -p "$UNIQUE_VARIANTS_VCF" \
  -r "$RESOURCES"

echo "VIOLA pipeline completed! Results are in $OUTPUT_DIR"
