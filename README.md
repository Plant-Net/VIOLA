# VIOLA

VIOLA is a patient-specific pipeline for variant prioritization in mitochondrial diseases, 
integrating genomics, transcriptomics, and phenotype data with machine learning.

# VIOLA workflow

<img src="viola_workflow.png" alt="VIOLA workflow" width="600"/>

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

#### R requirements

The following R libraries are required:
- dplyr
- tidyr


#### Usage

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

The ClinVar file can be found in `resources` (`clinvar_210125_hg38_cleaned.tsv`).


## viola_step2_cluster

This script runs the **Variational Autoencoder (VAE)** for dimensionality reduction and applies **DBSCAN** clustering to group outlier variants.

#### Python requirements

The following Python libraries are required:
- tensorflow
- sklearn
- pandas
- numpy


#### Usage

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

#### R requirements

The following R libraries are required:
- dplyr
- tidyr
- stringr
- data.table
- ontologyIndex


#### Usage

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


