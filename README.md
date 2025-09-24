# VIOLA

VIOLA is a patient-specific pipeline for variant prioritization in mitochondrial diseases, 
integrating genomics, transcriptomics, and phenotype data with machine learning.

# VIOLA workflow

<img width="4588" height="3511" alt="viola_workflow" src="https://github.com/user-attachments/assets/b2560189-5baa-45c1-a19e-844d1c54eac8" />

## 📦 Requirements

- R (≥ 4.2)
- Python (≥ 3.9)

# Scripts

The VIOLA pipeline is composed of 3 scripts :
- `viola_step1_merge.R`
- `viola_step2_cluster.py`
- `viola_step3_rank.R`

## viola_step1_merge

This script merges annotation datasets (e.g. VEP and CADD annotations) into a unified input table for downstream analysis.

### R requirements

The following R libraries are required:
- dplyr
- tidyr

### Required files

#### VEP file

VIOLA requires variants annotated with VEP (Variant Effect Predictor).
Please run VEP with the following options:

```bash
vep -i input_file.vcf -o output_file.vcf --cache --offline --assembly GRCh38 --vcf \
--check_existing --af --af_1kg --af_gnomad --af_esp
```

Then, run filter_vep to select only rare variants with the following options:

```bash
filter_vep -i input_file.vcf -o output_file.vcf -filter "SYMBOL and ((AF < 0.01 or gnomAD_AF < 0.01) or (not AF and not gnomAD_AF and not EUR_AF))"
```

Finally, run bcftools with the plugin split-vep to obtain a tabulated file.

```bash
echo -e "CHROM\tPOS\tREF\tALT\tQUAL\tGT\tAD\tDP\t$(bcftools +split-vep -l input_file.vcf | cut -f 2 | tr '\n' '\t' | sed 's/\t$//')" > output_file.tsv ; \
bcftools +split-vep -f '%CHROM\t%POS\t%REF\t%ALT\t%QUAL\t[%GT]\t[%AD]\t[%DP]\t%CSQ\n' -d -A tab input_file.vcf >> output_file.tsv
```

The output file (.tsv) is used as input for viola_step1_merge.R.

#### CADD file

VIOLA requires CADD (Combined Annotation Dependent Depletion) scores.
Please run CADD with the following options:

```bash
CADD.sh -a -g GRCh38 -o output.tsv.gz input_file.vcf
```

The output file (.tsv) is used as input for viola_step1_merge.R.

#### ClinVar input file

A preprocessed ClinVar file is provided in the repository `resources`: `clinvar_210125_hg38_cleaned.tsv`.

### Usage

`Rscript viola_step1_merge.R -h` will give you the following help message:

```bash
Usage: viola_step1_merge.R [options]

Options:
	-v VEP_INPUT, --vep_input=VEP_INPUT
		Path to VEP input file

	-c CADD_INPUT, --cadd_input=CADD_INPUT
		Path to CADD input file

	-l CLINVAR_INPUT, --clinvar_input=CLINVAR_INPUT
		Path to ClinVar input file

	-o OUTPUT_PATH, --output_path=OUTPUT_PATH
		Path to output directory

	-h, --help
		Show this help message and exit
```

## viola_step2_cluster

This script runs the **Variational Autoencoder (VAE)** for dimensionality reduction and applies **DBSCAN** clustering to group outlier variants.

### Python requirements

The following Python libraries are required:
- tensorflow
- sklearn
- pandas
- numpy

### Required file

The input file to process is the output file of `viola_step1_merge.R` script.

### Usage

`python viola_step2_cluster.py -h` will give you the following help message:

```bash
usage: viola_step2_cluster.py [-h] -f FILE_PATH -o OUTPUT_FOLDER_PATH

optional arguments:
  -h, --help            show this help message and exit
  -f FILE_PATH, --file_path FILE_PATH
                        Path to the file to process
  -o OUTPUT_FOLDER_PATH, --output_folder_path OUTPUT_FOLDER_PATH
                        Path to the output folder
```


## viola_step3_rank

This script applies **filtering** (quality, biotype, and Variant Allele Frequency), integrates **HPO** terms, and generates the final variant **ranking**.

### R requirements

The following R libraries are required:
- dplyr
- tidyr
- stringr
- data.table
- ontologyIndex

### Required files

#### Input file

The input file is the one of the outputs of the script `viola_step2_cluster.py` and follows the pattern: "res_dbscan.csv".

#### HPO table

The HPO table must be provided by the user. This is a 2-column CSV file like:

| patient_id |    hpo_id  |
|------------|------------|
| patient1   | HP:0001250 |
| patient1   | HP:0000518 |
| patient2   | HP:0001638 |

#### VCF of unique rare variants

This file can be obtained by filtering the original VCF using bcftools.
If no cohort is available to determine whether a variant is unique, please provide a VCF file containing only rare variants.

#### Resources folder

Resources folder contains reference files for transcriptomic co-expression matrices and mitochondrial gene lists.
The provided repository already includes the necessary files in resources.

### Usage

`Rscript viola_step3_rank.R -h` will give you the following help message:

```bash
Usage: viola_step3_rank.R [options]


Options:
	-f INPUT, --input=INPUT
		Path to input file

	-o OUTPUT_PATH, --output_path=OUTPUT_PATH
		Path to output directory

	-t HPO_TABLE, --hpo_table=HPO_TABLE
		Path to patient HPO table

	-p PATH_UNIQUE_VAR, --path_unique_var=PATH_UNIQUE_VAR
		Path to VCF of unique rare variants

	-r RESOURCES, --resources=RESOURCES
		Path to resources folder containing all required files

	-h, --help
		Show this help message and exit
```


# Installation
- Clone this repository
- Install dependencies (Python, R, etc.)
- Prepare input VCF and configuration files
---


